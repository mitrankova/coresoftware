#include "SiSeedingOrchestrator.h"

#include "SiHitSeedData.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/TpcCrossingDecisionv1.h>
#include <tpctrackreco/IdealPadMap.h>  // complete type for TpcCrossingFinder's unique_ptr<IdealPadMap> (destructor)
#include <tpctrackreco/TpcCrossingFinder.h>
#include <tpctrackreco/Tpc_AssembledTrack.h>
#include <tpctrackreco/Tpc_AssembledTrackContainer.h>
#include <tpctrackreco/Tpc_CrossingTrajectory.h>
#include <tpctrackreco/Tpc_CrossingTrajectoryContainer.h>

#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertexMap_v1.h>
#include <globalvertex/SvtxVertex_v3.h>

#include <trackbase/InttDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

  // a pass-2 candidate chain of one assembled track
  struct Candidate
  {
    int chain = -1;       // index in SiHitSeedEvent::chains
    short crossing = 0;   // crossing of its trajectory
    double z = 0;         // Si z at the beam [cm]
    double vertexZ = 0;   // bootstrap vertex it joins [cm]
    double pull = 0;      // (z - vertexZ) / sigma
    double score = 0;
  };

  bool isSelectedStatus(unsigned char s)
  {
    using S = TpcCrossingStatus;
    return s == static_cast<unsigned char>(S::SelectedByVertex) || s == static_cast<unsigned char>(S::SelectedByVertexAmbiguous) ||
           s == static_cast<unsigned char>(S::SelectedByVertexLoose) || s == static_cast<unsigned char>(S::SelectedByContainment) ||
           s == static_cast<unsigned char>(S::SelectedByContainmentAmbiguous);
  }
}  // namespace

SiSeedingOrchestrator::SiSeedingOrchestrator(const std::string& name)
  : SubsysReco(name)
  , m_finder(std::make_unique<SiHitSeedReco>("SiHitSeedReco"))
  , m_resolver(std::make_unique<TpcCrossingResolver>("TpcCrossingResolver"))
  , m_fitter(std::make_unique<SiTrajectoryFitter>("SiTrajectoryFitter"))
  , m_tpcOnlyFinder(std::make_unique<TpcCrossingFinder>("TpcCrossingFinderTpcOnly"))
{
  // pass-2 defaults: two MVTX layers are a chain (the vertex constraint replaces the third point)
  m_pass2Cfg.minLayers = 2;
  m_pass2Cfg.minMvtxLayers = 2;
}

SiSeedingOrchestrator::~SiSeedingOrchestrator() = default;

int SiSeedingOrchestrator::InitRun(PHCompositeNode* topNode)
{
  // the workers read each other's nodes: keep the names consistent
  m_resolver->setSiSeedNodeName(m_finder->outputNodeName());
  if (m_finder->usesCrossingTrajectories())
  {
    m_resolver->setTrajectoryNodeName(m_finder->trajectoryNodeName());
  }
  m_fitter->setInputNodeName(m_finder->outputNodeName());
  // crossings of the tracks without Si chain: here (tpcOnlyFinder) instead of the resolver's 3-point match
  m_tpcOnlyActive = m_tpcOnly && !m_triggeredMode && m_runResolver;
  m_resolver->setTpcOnlyMatch(!m_tpcOnlyActive);

  int ret = m_finder->InitRun(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  if (m_runResolver)
  {
    ret = m_resolver->InitRun(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return ret;
    }
  }
  if (m_runFit)
  {
    ret = m_fitter->InitRun(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return ret;
    }
  }
  if (!m_triggeredMode)
  {
    // collision vertex map (created before tpcOnlyFinder's InitRun, which looks for it)
    PHNodeIterator iter(topNode);
    auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
    }
    PHNodeIterator dstIter(dstNode);
    auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
    if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      dstNode->addNode(svtxNode);
    }
    if (!findNode::getClass<SvtxVertexMap>(topNode, m_collisionMapName))
    {
      svtxNode->addNode(new PHIODataNode<PHObject>(new SvtxVertexMap_v1, m_collisionMapName, "PHObject"));
    }
  }
  if (m_tpcOnlyActive)
  {
    m_tpcOnlyFinder->setTriggeredMode(false);
    m_tpcOnlyFinder->setInputNodeName(m_assembledNodeName);
    m_tpcOnlyFinder->setOutputNodeName(m_tpcOnlyDecisionName);
    m_tpcOnlyFinder->setVertexMapNodeName(m_collisionMapName);
    ret = m_tpcOnlyFinder->InitRun(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return ret;
    }
  }

  auto opt = [](const auto& o) { return o ? std::to_string(*o) : std::string("as pass 1"); };
  std::cout << Name() << "::InitRun - " << (m_triggeredMode ? "triggered" : "streaming") << " mode; workers: SiHitSeedReco"
            << (m_runResolver ? " + TpcCrossingResolver" : "") << (m_tpcOnlyActive ? " + TpcCrossingFinder (tracks without Si chain)" : "")
            << (m_runFit ? " + SiTrajectoryFitter" : "") << std::endl
            << Name() << ":   good pass-1 chain: >= " << m_goodMinLayers << " layers, >= " << m_goodMinMvtx << " MVTX"
            << (m_goodRequireIntt ? ", INTT hit" : "") << ", chi2/ndf <= " << m_goodMaxChi2 << ", tier <= " << m_goodMaxTier << std::endl
            << Name() << ":   bootstrap vertices: good chains, window " << m_vertexWindow << " cm, >= " << m_vertexMinTracks
            << " tracks; single-track Si z sigma " << m_trackZSigma << " cm" << std::endl;
  if (m_pass2)
  {
    std::cout << Name() << ":   pass 2: tracks without chain" << (m_pass2RedoWeak ? " + weak chains" : "")
              << " on free blobs; min layers " << opt(m_pass2Cfg.minLayers) << ", min MVTX " << opt(m_pass2Cfg.minMvtxLayers)
              << ", anchor z " << opt(m_pass2Cfg.anchorZCm) << " cm, seed z " << opt(m_pass2Cfg.seedZHalfCm)
              << " cm, MVTX propagation z " << opt(m_pass2Cfg.propZHalfBins) << " bins, agreement " << opt(m_pass2Cfg.agreementMaxChi2)
              << "; selection chi2/ndf <= " << m_pass2MaxChi2 << ", vertex within " << m_pass2VertexNSigma << " sigma, ambiguity ratio "
              << m_pass2AmbiguityRatio << ", tier " << static_cast<int>(m_pass2Tier) << std::endl;
  }
  if (!m_triggeredMode)
  {
    std::cout << Name() << ":   collisions: Si-matched tracks of tier <= " << m_collisionMaxTier << ", >= " << m_collisionMinTracks
              << " per cluster -> " << m_collisionMapName << " (z = mean TPC z0)" << std::endl;
  }
  if (m_tpcOnlyActive)
  {
    std::cout << Name() << ":   tracks without Si chain: TpcCrossingFinder against the collisions -> " << m_tpcOnlyDecisionName
              << "; vertex-matched -> tier " << m_tpcOnlyVertexTier << (m_tpcOnlyRejectAmbiguous ? " (ambiguous -> 2)" : "")
              << ", else tier 2" << (m_tpcOnlyAddMissing ? "; tracks without resolver decision added" : "") << std::endl;
  }
  if (m_pass3)
  {
    std::cout << Name() << ":   pass 3: hit recovery" << (m_pass3Cfg.mvtx ? " MVTX" : "") << (m_pass3Cfg.intt ? " INTT" : "")
              << ", windows MVTX +-" << m_pass3Cfg.mvtxPhiHalfBins << " Uphi / +-" << m_pass3Cfg.mvtxZHalfBins << " z bins, INTT +-"
              << m_pass3Cfg.inttPhiHalfBins << " / +-" << m_pass3Cfg.inttZHalfBins << ", chi2/ndf <= " << m_pass3Cfg.maxChi2
              << " and <= before + " << m_pass3Cfg.maxChi2Increase << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiSeedingOrchestrator::process_event(PHCompositeNode* topNode)
{
  ++m_nEvents;

  // ================================ PASS 1 ================================
  int ret = m_finder->process_event(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  if (m_runResolver)
  {
    ret = m_resolver->process_event(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return ret;
    }
  }

  SiHitSeedEvent* ev = m_finder->eventData();
  if (!ev)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  auto* decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  auto* trajectories = m_finder->usesCrossingTrajectories()
                           ? findNode::getClass<Tpc_CrossingTrajectoryContainer>(topNode, m_finder->trajectoryNodeName())
                           : nullptr;

  // TPC z0 of a TPC input (the resolver's definition with crossing trajectories)
  auto tpcZ0 = [&](int index) -> double
  {
    if (trajectories && index >= 0)
    {
      const Tpc_CrossingTrajectory* t = trajectories->get_trajectory(static_cast<unsigned int>(index));
      if (t && t->get_fit_status() == 1 && std::isfinite(t->get_z()))
      {
        return t->get_z();
      }
    }
    return m_finder->tpcInputZ0(index);
  };
  // Si z at the beam = TPC z0 + Si - TPC z offset
  auto siZ = [&](const SiHitChain& ch) { return tpcZ0(ch.tpc_track_index) + ch.tpc_dz_offset; };
  auto tierOf = [&](unsigned int id) -> int
  {
    if (!decisions)
    {
      return 0;
    }
    const TpcCrossingDecision* d = decisions->get_decision(id);
    return d ? static_cast<int>(d->get_selected_tier()) : 255;
  };
  auto betterChain = [](const SiHitChain& a, const SiHitChain& b)
  {
    if (a.n_layers() != b.n_layers())
    {
      return a.n_layers() > b.n_layers();
    }
    const bool fa = std::isfinite(a.tpc_chi2), fb = std::isfinite(b.tpc_chi2);
    if (fa != fb)
    {
      return fa;
    }
    return fa && a.tpc_chi2 < b.tpc_chi2;
  };

  // ---- classification: one selected chain per assembled track (the best if there are several)
  std::set<unsigned int> allTracks;
  for (const auto& s : ev->tpc_seeds)
  {
    allTracks.insert(s.tpc_assembled_track_id);
  }
  std::map<unsigned int, int> selectedOf;  // assembled id -> chain index
  for (std::size_t c = 0; c < ev->chains.size(); ++c)
  {
    const SiHitChain& ch = ev->chains[c];
    if (!ch.selected)
    {
      continue;
    }
    const auto [it, inserted] = selectedOf.try_emplace(ch.tpc_assembled_track_id, static_cast<int>(c));
    if (!inserted && betterChain(ch, ev->chains[it->second]))
    {
      it->second = static_cast<int>(c);
    }
  }
  auto isGood = [&](const SiHitChain& ch)
  {
    return ch.n_layers() >= m_goodMinLayers && ch.n_mvtx >= m_goodMinMvtx && (!m_goodRequireIntt || ch.n_intt > 0) &&
           std::isfinite(ch.tpc_chi2) && ch.tpc_chi2 <= m_goodMaxChi2 && tierOf(ch.tpc_assembled_track_id) <= m_goodMaxTier;
  };
  std::vector<char> good(ev->chains.size(), 0);
  std::set<unsigned int> weakIds, noneIds;
  for (const unsigned int id : allTracks)
  {
    const auto it = selectedOf.find(id);
    if (it == selectedOf.end())
    {
      noneIds.insert(id);
    }
    else if (isGood(ev->chains[it->second]))
    {
      good[it->second] = 1;
    }
    else
    {
      weakIds.insert(id);
    }
  }
  const unsigned int nGood = static_cast<unsigned int>(std::count(good.begin(), good.end(), 1));
  m_nGood += nGood;
  m_nWeak += weakIds.size();
  m_nNone += noneIds.size();

  // ================================ BOOTSTRAP VERTICES ================================
  std::vector<ZEntry> goodZ;
  for (std::size_t c = 0; c < ev->chains.size(); ++c)
  {
    if (good[c])
    {
      const double z = siZ(ev->chains[c]);
      if (std::isfinite(z))
      {
        goodZ.push_back({ev->chains[c].tpc_crossing, z, static_cast<int>(c)});
      }
    }
  }
  std::vector<BootVertex> vertices;
  for (const auto& cl : clusterZ(std::move(goodZ), m_vertexMinTracks))
  {
    BootVertex v;
    v.crossing = cl.front().crossing;
    for (const auto& e : cl)
    {
      v.z += e.z;
      v.members.push_back(e.index);
    }
    v.n = static_cast<unsigned int>(cl.size());
    v.z /= v.n;
    vertices.push_back(std::move(v));
  }
  m_nVertices += vertices.size();
  for (const auto& v : vertices)
  {
    m_nVertexTracks += v.n;
  }

  // ---- INTT time buckets of the vertices (one entry per INTT layer of every member chain) and
  // vertices in neighbouring crossings: is the smaller one made of the bigger one's tracks, placed one
  // crossing off (the +-1 INTT crossing window lets INTT hits of crossing c time crossing c +- 1)?
  auto bucketsOf = [&](const BootVertex& v)
  {
    std::map<int, unsigned int> b;
    for (const int c : v.members)
    {
      const SiHitChain& ch = ev->chains[c];
      for (std::size_t i = 0; i < ch.layers.size() && i < ch.point_hit_ids.size(); ++i)
      {
        const int hid = ch.point_hit_ids[i];
        if (ch.layers[i] >= 3 && hid >= 0 && hid < static_cast<int>(ev->hits.size()))
        {
          ++b[static_cast<int>(InttDefs::getTimeBucketId(ev->hits[hid].hitsetkey))];
        }
      }
    }
    return b;
  };
  std::vector<std::string> adjacentLines;
  for (std::size_t i = 0; i < vertices.size(); ++i)
  {
    for (std::size_t j = 0; j < vertices.size(); ++j)
    {
      // ordered pair (big, small): |crossing difference| = 1, small = fewer tracks (ties: later crossing)
      const auto& big = vertices[i];
      const auto& small = vertices[j];
      if (std::abs(big.crossing - small.crossing) != 1 || big.n < small.n || (big.n == small.n && big.crossing > small.crossing))
      {
        continue;
      }
      ++m_nAdjacentPairs;
      unsigned int own = 0, neighbour = 0, other = 0;
      for (const auto& [b, n] : bucketsOf(small))
      {
        (b == small.crossing ? own : (b == big.crossing ? neighbour : other)) += n;
      }
      m_adjSmallOwn += own;
      m_adjSmallNeighbour += neighbour;
      m_adjSmallOther += other;
      m_adjDz.push_back(small.z - big.z);
      if (Verbosity() > 0)
      {
        char line[256];
        std::snprintf(line, sizeof(line),
                      "  neighbouring vertices: crossing %d (z %.2f cm, %u tracks) and %d (z %.2f cm, %u tracks), dz %.2f cm; INTT buckets of the "
                      "smaller one: own crossing %u, the other vertex's crossing %u, other %u",
                      big.crossing, big.z, big.n, small.crossing, small.z, small.n, small.z - big.z, own, neighbour, other);
        adjacentLines.emplace_back(line);
      }
    }
  }

  // ================================ PASS 2 ================================
  unsigned int nRescued = 0, nImproved = 0, nP2Chains = 0;
  if (m_pass2 && !vertices.empty())
  {
    std::set<unsigned int> targets = noneIds;
    if (m_pass2RedoWeak)
    {
      targets.insert(weakIds.begin(), weakIds.end());
    }
    if (!targets.empty())
    {
      m_nPass2Targets += targets.size();
      // free blobs: everything except the chains that are kept as they are
      std::vector<char> excluded(ev->blobs.size(), 0);
      for (const auto& ch : ev->chains)
      {
        const bool redo = m_pass2RedoWeak && weakIds.count(ch.tpc_assembled_track_id) != 0;
        if (!ch.selected || redo)
        {
          continue;
        }
        for (const int b : ch.blob_ids)
        {
          if (b >= 0 && b < static_cast<int>(excluded.size()))
          {
            excluded[b] = 1;
          }
        }
      }
      // TPC inputs to chain again: the target tracks' trajectories, without the ones that cannot change
      // (no usable TPC prediction does not depend on the free blobs) and, in streaming mode, only at
      // crossings that have a bootstrap vertex (the selection below needs one anyway)
      std::set<short> vertexCrossings;
      for (const auto& v : vertices)
      {
        vertexCrossings.insert(v.crossing);
      }
      std::set<int> inputs;
      for (const auto& info : ev->tpc_seeds)
      {
        if (targets.count(info.tpc_assembled_track_id) == 0)
        {
          continue;
        }
        if (m_pass2SkipUnpredictable && (info.status == SiTpcSeedInfo::TooFewClusters || info.status == SiTpcSeedInfo::FitFailed ||
                                         info.status == SiTpcSeedInfo::NoPrediction))
        {
          ++m_nPass2SkippedNoPred;
          continue;
        }
        if (m_pass2VertexCrossingsOnly && !m_triggeredMode && vertexCrossings.count(info.tpc_crossing) == 0)
        {
          ++m_nPass2SkippedNoVertex;
          continue;
        }
        inputs.insert(info.tpc_track_index);
      }
      m_nPass2Inputs += inputs.size();
      const std::vector<int> newChains = m_finder->rechainInputs(inputs, excluded, m_pass2Cfg, 2);
      nP2Chains = static_cast<unsigned int>(newChains.size());
      m_nPass2Chains += newChains.size();

      // ---- vertex-constrained candidates
      std::map<unsigned int, std::vector<Candidate>> cands;
      for (const int cid : newChains)
      {
        const SiHitChain& ch = ev->chains[cid];
        if (!std::isfinite(ch.tpc_chi2) || ch.tpc_chi2 > m_pass2MaxChi2)
        {
          continue;
        }
        const double z = siZ(ch);
        if (!std::isfinite(z))
        {
          continue;
        }
        const BootVertex* best = nullptr;
        double bestPull = 0;
        for (const auto& v : vertices)
        {
          if (v.crossing != ch.tpc_crossing)
          {
            continue;
          }
          const double sigma = m_trackZSigma * std::sqrt(1.0 + 1.0 / v.n);
          const double pull = (z - v.z) / sigma;
          if (!best || std::fabs(pull) < std::fabs(bestPull))
          {
            best = &v;
            bestPull = pull;
          }
        }
        if (!best || std::fabs(bestPull) > m_pass2VertexNSigma)
        {
          ++m_nPass2NoVertex;
          continue;
        }
        Candidate c;
        c.chain = cid;
        c.crossing = ch.tpc_crossing;
        c.z = z;
        c.vertexZ = best->z;
        c.pull = bestPull;
        c.score = std::min(1.0, ch.n_layers() / 7.0) * std::exp(-0.5 * ch.tpc_chi2) * std::exp(-0.5 * bestPull * bestPull) *
                  (ch.n_intt > 0 ? 1.0 : m_noInttFactor);
        cands[ch.tpc_assembled_track_id].push_back(c);
      }

      // ---- one crossing per track: the best candidate must beat every other crossing
      struct Choice
      {
        unsigned int id = 0;
        std::vector<Candidate> list;  // candidates of the chosen crossing, best first
      };
      std::vector<Choice> choices;
      for (auto& [id, list] : cands)
      {
        std::sort(list.begin(), list.end(), [](const Candidate& a, const Candidate& b) { return a.score > b.score; });
        const Candidate& best = list.front();
        double second = 0;
        for (const auto& c : list)
        {
          if (c.crossing != best.crossing)
          {
            second = std::max(second, c.score);
          }
        }
        if (second > 0 && best.score < m_pass2AmbiguityRatio * second)
        {
          ++m_nPass2Ambiguous;
          continue;
        }
        Choice ch;
        ch.id = id;
        for (const auto& c : list)
        {
          if (c.crossing == best.crossing)
          {
            ch.list.push_back(c);
          }
        }
        choices.push_back(std::move(ch));
      }
      std::sort(choices.begin(), choices.end(), [](const Choice& a, const Choice& b) { return a.list.front().score > b.list.front().score; });

      // ---- ownership among the pass-2 winners, best tracks first.  The weak chains own their blobs
      // until their track takes a better chain (a weak track that keeps its chain keeps its hits).
      std::unordered_map<int, unsigned int> owner;  // blob -> assembled id
      if (m_pass2RedoWeak)
      {
        for (const unsigned int id : weakIds)
        {
          for (const int b : ev->chains[selectedOf.at(id)].blob_ids)
          {
            owner.emplace(b, id);
          }
        }
      }
      for (const auto& choice : choices)
      {
        const Candidate* taken = nullptr;
        for (const auto& c : choice.list)
        {
          unsigned int shared = 0;
          for (const int b : ev->chains[c.chain].blob_ids)
          {
            const auto it = owner.find(b);
            shared += it != owner.end() && it->second != choice.id;
          }
          if (shared < m_pass2MinShared)
          {
            taken = &c;
            break;
          }
        }
        if (!taken)
        {
          ++m_nPass2Shared;
          continue;
        }
        SiHitChain& nc = ev->chains[taken->chain];
        const auto is = selectedOf.find(choice.id);
        const int old = is != selectedOf.end() ? is->second : -1;
        if (old >= 0)
        {
          // weak track: only a better chain replaces its pass-1 chain
          if (!betterChain(nc, ev->chains[old]))
          {
            continue;
          }
          ev->chains[old].selected = false;
          for (const int b : ev->chains[old].blob_ids)
          {
            const auto it = owner.find(b);
            if (it != owner.end() && it->second == choice.id)
            {
              owner.erase(it);
            }
          }
          for (auto& info : ev->tpc_seeds)
          {
            if (info.selected_chain == ev->chains[old].id)
            {
              info.selected_chain = -1;
            }
          }
          ++nImproved;
        }
        else
        {
          ++nRescued;
        }
        for (const int b : nc.blob_ids)
        {
          owner.emplace(b, choice.id);
        }
        nc.selected = true;
        selectedOf[choice.id] = taken->chain;
        for (auto& info : ev->tpc_seeds)
        {
          if (info.tpc_track_index == nc.tpc_track_index)
          {
            info.selected_chain = nc.id;
            info.status = SiTpcSeedInfo::Chained;
          }
        }

        // decision: the crossing of the vertex the chain joins (streaming only; triggered = crossing fixed)
        if (!m_triggeredMode && decisions)
        {
          TpcCrossingDecision* d = decisions->get_decision(choice.id);
          if (!d)
          {
            d = new TpcCrossingDecisionv1();
            d->set_assembled_track_id(choice.id);
            decisions->add_decision(d);
          }
          const bool sameCrossing = d->get_selected_tier() != std::numeric_limits<unsigned char>::max() &&
                                    d->get_selected_crossing() == taken->crossing;
          const unsigned char tier = sameCrossing ? std::min(d->get_selected_tier(), m_pass2Tier) : m_pass2Tier;
          d->set_selected_crossing(taken->crossing);
          d->set_selected_tier(tier);
          d->set_selected_score(static_cast<float>(taken->score));
          d->set_tpc_z0(static_cast<float>(tpcZ0(nc.tpc_track_index)));
          d->set_silicon_vertex_z(static_cast<float>(taken->vertexZ));
          d->set_delta_z(static_cast<float>(nc.tpc_dz_offset));
          d->set_best_abs_delta_z(static_cast<float>(std::fabs(nc.tpc_dz_offset)));
          d->set_status(TpcCrossingStatus::SelectedByVertexLoose);
        }
      }
      m_nPass2Rescued += nRescued;
      m_nPass2Improved += nImproved;
      if (m_pass2RedoWeak)
      {
        m_nPass2WeakKept += weakIds.size() - nImproved;
      }
    }
  }
  else if (m_pass2 && m_pass2RedoWeak)
  {
    m_nPass2WeakKept += weakIds.size();
  }

  // ================================ PASS 3 ================================
  std::array<unsigned int, 7> added{};
  if (m_pass3)
  {
    std::vector<char> excluded(ev->blobs.size(), 0);
    std::vector<int> ids;
    for (std::size_t c = 0; c < ev->chains.size(); ++c)
    {
      const SiHitChain& ch = ev->chains[c];
      if (!ch.selected)
      {
        continue;
      }
      for (const int b : ch.blob_ids)
      {
        if (b >= 0 && b < static_cast<int>(excluded.size()))
        {
          excluded[b] = 1;
        }
      }
      if (ch.n_layers() < 7)
      {
        ids.push_back(static_cast<int>(c));
      }
    }
    // best chains first: they get the free blobs if two chains want the same one
    std::sort(ids.begin(), ids.end(), [&](int a, int b) { return betterChain(ev->chains[a], ev->chains[b]); });
    m_nPass3Chains += ids.size();
    // bootstrap vertex z of every selected chain (its crossing, within the vertex window), else NaN
    std::vector<double> vertexZ(ev->chains.size(), std::numeric_limits<double>::quiet_NaN());
    for (const int c : ids)
    {
      const double z = siZ(ev->chains[c]);
      double best = m_vertexWindow;
      for (const auto& v : vertices)
      {
        if (v.crossing == ev->chains[c].tpc_crossing && std::isfinite(z) && std::fabs(z - v.z) <= best)
        {
          best = std::fabs(z - v.z);
          vertexZ[c] = v.z;
        }
      }
      m_nPass3WithVertex += std::isfinite(vertexZ[c]);
    }
    added = m_finder->recoverHits(ids, excluded, m_pass3Cfg, &vertexZ);
    for (int l = 0; l < 7; ++l)
    {
      m_pass3Added[l] += added[l];
    }
  }

  // ================================ COLLISIONS + TPC-ONLY ================================
  if (!m_triggeredMode && decisions)
  {
    makeCollisions(topNode, selectedOf);
    if (m_tpcOnlyActive)
    {
      ret = m_tpcOnlyFinder->process_event(topNode);
      if (ret != Fun4AllReturnCodes::EVENT_OK)
      {
        return ret;
      }
      decideTpcOnly(topNode, selectedOf);
    }
  }

  if (Verbosity() > 0)
  {
    unsigned int nAdded = 0;
    for (const auto n : added)
    {
      nAdded += n;
    }
    std::cout << Name() << ": event " << m_nEvents << " tracks " << allTracks.size() << " | pass 1: good " << nGood << ", weak "
              << weakIds.size() << ", no chain " << noneIds.size() << " | bootstrap vertices " << vertices.size()
              << " | pass 2: " << nP2Chains << " candidate chains, rescued " << nRescued << ", improved " << nImproved
              << " | pass 3: layers added " << nAdded << std::endl;
    if (Verbosity() > 1)
    {
      for (const auto& v : vertices)
      {
        std::cout << Name() << ":   bootstrap vertex crossing " << v.crossing << " z " << v.z << " cm, " << v.n << " good tracks, INTT time buckets:";
        for (const auto& [b, n] : bucketsOf(v))
        {
          std::cout << " " << b << "(" << n << ")";
        }
        std::cout << std::endl;
      }
    }
    for (const auto& l : adjacentLines)
    {
      std::cout << Name() << ":" << l << std::endl;
    }
  }

  // ================================ FIT ================================
  if (m_runFit)
  {
    ret = m_fitter->process_event(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return ret;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiSeedingOrchestrator::End(PHCompositeNode* topNode)
{
  m_finder->End(topNode);
  if (m_runResolver)
  {
    m_resolver->End(topNode);
  }
  if (m_runFit)
  {
    m_fitter->End(topNode);
  }

  const double ne = std::max<unsigned long>(1, m_nEvents);
  std::cout << Name() << "::End - " << m_nEvents << " events" << std::endl
            << Name() << ":   pass 1: good " << m_nGood << ", weak " << m_nWeak << ", no chain " << m_nNone << " assembled tracks" << std::endl
            << Name() << ":   bootstrap vertices " << m_nVertices << " (" << m_nVertices / ne << " per event, "
            << (m_nVertices > 0 ? static_cast<double>(m_nVertexTracks) / m_nVertices : 0.0) << " good tracks per vertex)" << std::endl;
  if (m_pass2)
  {
    std::cout << Name() << ":   pass 2: " << m_nPass2Targets << " target tracks, " << m_nPass2Inputs << " TPC inputs chained again (skipped: "
              << m_nPass2SkippedNoPred << " without usable TPC prediction, " << m_nPass2SkippedNoVertex << " at crossings without a bootstrap vertex), "
              << m_nPass2Chains << " candidate chains ("
              << m_nPass2NoVertex << " not at a bootstrap vertex of their crossing); rescued " << m_nPass2Rescued << ", improved weak "
              << m_nPass2Improved << ", weak kept " << m_nPass2WeakKept << ", ambiguous between crossings " << m_nPass2Ambiguous
              << ", all candidates taken by better tracks " << m_nPass2Shared << std::endl;
  }
  if (m_pass3)
  {
    unsigned long nAdded = 0;
    for (const auto n : m_pass3Added)
    {
      nAdded += n;
    }
    std::cout << Name() << ":   pass 3: " << m_nPass3Chains << " selected chains with missing layers (" << m_nPass3WithVertex
              << " at a bootstrap vertex), " << nAdded << " layers added (L0-6:";
    for (const auto n : m_pass3Added)
    {
      std::cout << " " << n;
    }
    std::cout << ")" << std::endl;
  }
  if (!m_triggeredMode)
  {
    std::cout << Name() << ":   collisions " << m_nCollisions << " (" << m_nCollisions / ne << " per event, "
              << (m_nCollisions > 0 ? static_cast<double>(m_nCollisionTracks) / m_nCollisions : 0.0) << " Si-matched tracks per collision), mean Si z - TPC z0 "
              << (m_nCollisionTracks > 0 ? m_sumSiMinusTpc / m_nCollisionTracks : kNaN) << " cm" << std::endl;
  }
  if (m_tpcOnlyActive)
  {
    std::cout << Name() << ":   tracks without Si chain: " << m_nTpcOnlyTracks << " decided by TpcCrossingFinder (no finder decision "
              << m_nTpcOnlyNoDecision << ", added " << m_nTpcOnlyAdded << ", crossing changed " << m_nTpcOnlyMoved << "); tier old -> new:";
    for (const auto& [k, n] : m_tpcOnlyTiers)
    {
      std::cout << " " << k.first << "->" << k.second << " " << n;
    }
    std::cout << std::endl;
  }
  {
    const unsigned long nb = m_adjSmallOwn + m_adjSmallNeighbour + m_adjSmallOther;
    std::sort(m_adjDz.begin(), m_adjDz.end());
    std::cout << Name() << ":   bootstrap vertices in neighbouring crossings: " << m_nAdjacentPairs << " pairs";
    if (nb > 0)
    {
      std::cout << "; INTT hits of the smaller vertex: at its own crossing " << m_adjSmallOwn << ", at the other vertex's crossing "
                << m_adjSmallNeighbour << ", elsewhere " << m_adjSmallOther
                << " (mostly 'other vertex' = the same collision split over two crossings)";
    }
    if (!m_adjDz.empty())
    {
      std::cout << "; z(small) - z(big) median " << m_adjDz[m_adjDz.size() / 2] << " cm";
    }
    std::cout << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________________________________
std::vector<std::vector<SiSeedingOrchestrator::ZEntry>> SiSeedingOrchestrator::clusterZ(std::vector<ZEntry> entries, unsigned int minTracks) const
{
  std::sort(entries.begin(), entries.end(), [](const ZEntry& a, const ZEntry& b)
            { return a.crossing != b.crossing ? a.crossing < b.crossing : a.z < b.z; });
  std::vector<std::vector<ZEntry>> out;
  std::vector<ZEntry> cur;
  double sum = 0;
  auto flush = [&]()
  {
    if (cur.size() >= minTracks && !cur.empty())
    {
      out.push_back(cur);
    }
    cur.clear();
    sum = 0;
  };
  for (const auto& e : entries)
  {
    if (!cur.empty() && (e.crossing != cur.front().crossing || std::fabs(e.z - sum / cur.size()) > m_vertexWindow))
    {
      flush();
    }
    cur.push_back(e);
    sum += e.z;
  }
  flush();
  return out;
}

//____________________________________________________________________________________________________
// final collisions of the event -> collision vertex map (z in the TPC z0 frame, see the class comment)
void SiSeedingOrchestrator::makeCollisions(PHCompositeNode* topNode, const std::map<unsigned int, int>& selectedOf)
{
  auto* map = findNode::getClass<SvtxVertexMap>(topNode, m_collisionMapName);
  const auto* decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  if (!map || !decisions)
  {
    return;
  }
  map->clear();
  std::vector<ZEntry> entries;
  std::vector<std::pair<double, unsigned int>> tpcZ;  // per entry: TPC z0, assembled id
  for (const auto& [id, chain] : selectedOf)
  {
    (void) chain;
    const TpcCrossingDecision* d = decisions->get_decision(id);
    if (!d || d->get_selected_crossing() == SHRT_MAX || static_cast<unsigned int>(d->get_selected_tier()) > m_collisionMaxTier)
    {
      continue;
    }
    const double z0 = d->get_tpc_z0(), dz = d->get_delta_z();
    if (!std::isfinite(z0) || !std::isfinite(dz) || !std::isfinite(d->get_silicon_vertex_z()))  // untimed lone: no collision
    {
      continue;
    }
    entries.push_back({d->get_selected_crossing(), z0 + dz, static_cast<int>(tpcZ.size())});
    tpcZ.emplace_back(z0, id);
  }
  for (const auto& cl : clusterZ(std::move(entries), m_collisionMinTracks))
  {
    double zSi = 0, zTpc = 0;
    auto v = std::make_unique<SvtxVertex_v3>();
    for (const auto& e : cl)
    {
      zSi += e.z;
      zTpc += tpcZ[e.index].first;
      v->insert_track(tpcZ[e.index].second);
    }
    const double n = static_cast<double>(cl.size());
    zSi /= n;
    zTpc /= n;
    v->set_x(0);
    v->set_y(0);
    v->set_z(static_cast<float>(zTpc));
    v->set_t(0);
    v->set_chisq(0);
    v->set_ndof(0);
    v->set_beam_crossing(cl.front().crossing);
    for (unsigned int a = 0; a < 3; ++a)
    {
      for (unsigned int b = 0; b < 3; ++b)
      {
        v->set_error(a, b, 0);
      }
    }
    v->set_error(2, 2, static_cast<float>(m_tpcZSigma * m_tpcZSigma / n));
    map->insert(v.release());
    ++m_nCollisions;
    m_nCollisionTracks += cl.size();
    m_sumSiMinusTpc += (zSi - zTpc) * n;
    if (Verbosity() > 1)
    {
      std::cout << Name() << ":   collision crossing " << cl.front().crossing << ": " << cl.size() << " Si-matched tracks, Si z " << zSi
                << " cm, TPC z0 " << zTpc << " cm" << std::endl;
    }
  }
}

//____________________________________________________________________________________________________
// tracks without a selected Si chain take the TpcCrossingFinder decision (see the class comment)
void SiSeedingOrchestrator::decideTpcOnly(PHCompositeNode* topNode, const std::map<unsigned int, int>& selectedOf)
{
  auto* decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  const auto* finder = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_tpcOnlyDecisionName);
  const auto* assembled = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_assembledNodeName);
  if (!decisions || !finder || !assembled)
  {
    return;
  }
  for (unsigned int i = 0; i < assembled->size(); ++i)
  {
    const Tpc_AssembledTrack* t = assembled->get_track(i);
    if (!t)
    {
      continue;
    }
    const unsigned int id = t->get_track_id();
    if (selectedOf.count(id) != 0)
    {
      continue;  // Si-matched: resolver / pass-2 decision stays
    }
    TpcCrossingDecision* d = decisions->get_decision(id);
    if (!d && !m_tpcOnlyAddMissing)
    {
      continue;
    }
    ++m_nTpcOnlyTracks;
    const TpcCrossingDecision* f = finder->get_decision(id);
    if (!f)
    {
      ++m_nTpcOnlyNoDecision;
      continue;
    }
    const int oldTier = d ? static_cast<int>(d->get_selected_tier()) : -1;
    const short oldCrossing = d ? d->get_selected_crossing() : static_cast<short>(SHRT_MAX);
    const bool valid = isSelectedStatus(f->get_status()) && f->get_selected_crossing() != SHRT_MAX;
    const bool ambiguous = f->get_status() == static_cast<unsigned char>(TpcCrossingStatus::SelectedByVertexAmbiguous);
    const unsigned char tier = (valid && f->get_selected_tier() == 0 && !(m_tpcOnlyRejectAmbiguous && ambiguous))
                                   ? static_cast<unsigned char>(m_tpcOnlyVertexTier)
                                   : static_cast<unsigned char>(2);
    if (!d)
    {
      if (!valid)
      {
        continue;
      }
      auto* nd = new TpcCrossingDecisionv1();
      nd->set_assembled_track_id(id);
      decisions->add_decision(nd);
      d = decisions->get_decision(id);
      if (!d)
      {
        continue;
      }
      ++m_nTpcOnlyAdded;
    }
    if (valid)
    {
      d->set_selected_crossing(f->get_selected_crossing());
      d->set_tpc_z0(f->get_tpc_z0());
      d->set_silicon_vertex_id(f->get_silicon_vertex_id());
      d->set_silicon_vertex_z(f->get_selected_tier() <= 1 ? f->get_silicon_vertex_z() : std::numeric_limits<float>::quiet_NaN());
      d->set_delta_z(f->get_delta_z());
      d->set_best_abs_delta_z(f->get_best_abs_delta_z());
      d->set_second_best_abs_delta_z(f->get_second_best_abs_delta_z());
      d->set_selected_score(f->get_selected_score());
      d->set_number_of_available_crossings(f->get_number_of_available_crossings());
      d->set_number_of_allowed_crossings(f->get_number_of_allowed_crossings());
      d->set_number_of_tested_crossings(f->get_number_of_tested_crossings());
      d->set_number_of_tpc_valid_crossings(f->get_number_of_tpc_valid_crossings());
      d->set_number_of_vertex_compatible_crossings(f->get_number_of_vertex_compatible_crossings());
      d->set_status(f->get_status());
      m_nTpcOnlyMoved += oldCrossing != SHRT_MAX && oldCrossing != f->get_selected_crossing();
    }
    d->set_selected_tier(tier);
    ++m_tpcOnlyTiers[{oldTier, static_cast<int>(tier)}];
    if (Verbosity() > 1)
    {
      std::cout << Name() << ":   assembled " << id << " (no Si chain): crossing " << oldCrossing << " tier " << oldTier << " -> "
                << (valid ? "crossing " + std::to_string(f->get_selected_crossing()) : std::string("no valid crossing")) << " tier "
                << static_cast<int>(tier) << " (finder tier " << static_cast<int>(f->get_selected_tier()) << ", TPC z0 " << f->get_tpc_z0()
                << ", collision z " << f->get_silicon_vertex_z() << ", dz " << f->get_delta_z() << ", 2nd best |dz| "
                << f->get_second_best_abs_delta_z() << ")" << std::endl;
    }
  }
}
