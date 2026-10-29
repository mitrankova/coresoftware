#include "SiSeedingOrchestrator.h"

#include "SiHitSeedData.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/TpcCrossingDecisionv1.h>
#include <tpctrackreco/Tpc_CrossingTrajectory.h>
#include <tpctrackreco/Tpc_CrossingTrajectoryContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <set>
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
}  // namespace

SiSeedingOrchestrator::SiSeedingOrchestrator(const std::string& name)
  : SubsysReco(name)
  , m_finder(std::make_unique<SiHitSeedReco>("SiHitSeedReco"))
  , m_resolver(std::make_unique<TpcCrossingResolver>("TpcCrossingResolver"))
  , m_fitter(std::make_unique<SiTrajectoryFitter>("SiTrajectoryFitter"))
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

  auto opt = [](const auto& o) { return o ? std::to_string(*o) : std::string("as pass 1"); };
  std::cout << Name() << "::InitRun - " << (m_triggeredMode ? "triggered" : "streaming") << " mode; workers: SiHitSeedReco"
            << (m_runResolver ? " + TpcCrossingResolver" : "") << (m_runFit ? " + SiTrajectoryFitter" : "") << std::endl
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
  std::map<short, std::vector<double>> zByCrossing;
  for (std::size_t c = 0; c < ev->chains.size(); ++c)
  {
    if (good[c])
    {
      const double z = siZ(ev->chains[c]);
      if (std::isfinite(z))
      {
        zByCrossing[ev->chains[c].tpc_crossing].push_back(z);
      }
    }
  }
  std::vector<BootVertex> vertices;
  for (auto& [crossing, zs] : zByCrossing)
  {
    std::sort(zs.begin(), zs.end());
    double sum = 0;
    unsigned int n = 0;
    auto flush = [&]()
    {
      if (n >= m_vertexMinTracks)
      {
        vertices.push_back({crossing, sum / n, n});
      }
      sum = 0;
      n = 0;
    };
    for (const double z : zs)
    {
      if (n > 0 && std::fabs(z - sum / n) > m_vertexWindow)
      {
        flush();
      }
      sum += z;
      ++n;
    }
    flush();
  }
  m_nVertices += vertices.size();
  for (const auto& v : vertices)
  {
    m_nVertexTracks += v.n;
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
      const std::vector<int> newChains = m_finder->rechain(targets, excluded, m_pass2Cfg, 2);
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
    added = m_finder->recoverHits(ids, excluded, m_pass3Cfg);
    for (int l = 0; l < 7; ++l)
    {
      m_pass3Added[l] += added[l];
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
        std::cout << Name() << ":   bootstrap vertex crossing " << v.crossing << " z " << v.z << " cm, " << v.n << " good tracks" << std::endl;
      }
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
    std::cout << Name() << ":   pass 2: " << m_nPass2Targets << " target tracks, " << m_nPass2Chains << " candidate chains ("
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
    std::cout << Name() << ":   pass 3: " << m_nPass3Chains << " selected chains with missing layers, " << nAdded << " layers added (L0-6:";
    for (const auto n : m_pass3Added)
    {
      std::cout << " " << n;
    }
    std::cout << ")" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
