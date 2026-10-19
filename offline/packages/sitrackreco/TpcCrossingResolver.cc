#include "TpcCrossingResolver.h"

#include "SiHitSeedData.h"

#include <tpctrackreco/TpcCrossingDecisionContainerv1.h>
#include <tpctrackreco/TpcCrossingDecisionv1.h>
#include <tpctrackreco/Tpc_AssembledTrack.h>
#include <tpctrackreco/Tpc_AssembledTrackContainer.h>
#include <tpctrackreco/Tpc_CrossingTrajectory.h>
#include <tpctrackreco/Tpc_CrossingTrajectoryContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <algorithm>
#include <format>
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

namespace
{
  // one crossing hypothesis of a track: its crossing trajectory and the best Si chain of it
  struct Hypothesis
  {
    unsigned int assembled = 0;
    short crossing = 0;
    int chain = -1;      // index in SiHitSeedEvent::chains
    double dz = 0;       // Si - TPC z offset at this crossing [cm]
    double tpcZ = 0;     // TPC z0 at this crossing [cm]
    double z = 0;        // collision z = Si z at the beam [cm]
    double match = 0;    // chain <-> trajectory
    double quality = 0;  // chain by itself
    double credit = 0;   // match * quality
    double vote = 0;     // credit * credit share of the track
    bool inMain = false; // voted for the main vertex of its crossing (pass 1)
    bool timed = false;  // the chain has INTT hits, i.e. a time (crossing) of its own
    double event = 0;    // leave-one-out support of the main vertex of its crossing
    double score = 0;
  };

  // vertex = 1-D cluster of weighted z in one crossing
  struct Vertex
  {
    double sumW = 0, sumWZ = 0;
    unsigned int n = 0;
    double mean() const { return sumW > 0 ? sumWZ / sumW : 0.0; }
  };

  // clusters (z, w) of one crossing; returns the vertices and the vertex index of every entry
  std::vector<Vertex> clusterZ(const std::vector<std::pair<double, double>>& zw, double window, std::vector<int>& index)
  {
    std::vector<std::size_t> order(zw.size());
    for (std::size_t i = 0; i < order.size(); ++i)
    {
      order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return zw[a].first < zw[b].first; });
    std::vector<Vertex> v;
    index.assign(zw.size(), -1);
    for (const std::size_t i : order)
    {
      if (v.empty() || std::fabs(zw[i].first - v.back().mean()) > window)
      {
        v.emplace_back();
      }
      v.back().sumW += zw[i].second;
      v.back().sumWZ += zw[i].second * zw[i].first;
      ++v.back().n;
      index[i] = static_cast<int>(v.size()) - 1;
    }
    return v;
  }

  int mainVertex(const std::vector<Vertex>& v)
  {
    int best = -1;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
      if (best < 0 || v[i].sumW > v[best].sumW)
      {
        best = static_cast<int>(i);
      }
    }
    return best;
  }

  constexpr unsigned char NoTier = std::numeric_limits<unsigned char>::max();

  double median(std::vector<double> v)
  {
    std::sort(v.begin(), v.end());
    return v.empty() ? 0.0 : (v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]));
  }
}  // namespace

TpcCrossingResolver::TpcCrossingResolver(const std::string& name)
  : SubsysReco(name)
{
}

int TpcCrossingResolver::InitRun(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  m_decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_outputNodeName);
  if (!m_decisions)
  {
    m_decisions = new TpcCrossingDecisionContainerv1();
    dstNode->addNode(new PHIODataNode<PHObject>(m_decisions, m_outputNodeName, "PHObject"));
  }
  std::cout << Name() << "::InitRun - "
            << (m_triggeredMode ? "triggered mode, crossing " + std::to_string(m_triggeredCrossing)
                                : "streaming mode, one hypothesis per crossing trajectory (" + m_trajectoryNodeName + " + " +
                                      m_siSeedNodeName + "): best chain -> credit -> vertex vote -> one collision per crossing")
            << " -> " << m_outputNodeName << std::endl;
  if (!m_triggeredMode && Verbosity() > 0)
  {
    std::cout << Name() << ":   credit = exp(-dz^2/2 " << m_dzSigma << "^2) exp(-" << m_rphiWeight << " chi2_rphi/n /2)"
              << " x (n_layers/7)^" << m_layerPower << " exp(-" << m_ddphiWeight << " <ddphi pull^2>/2) x " << m_noInttFactor
              << " if no INTT; chains chi2/ndf <= " << m_maxChainChi2 << "; vertices: window " << m_clusterZWindow
              << " cm, single-track Si z sigma " << m_siZSigma << " cm (support: x sqrt(1 + 1/n_other)), support scale " << m_supportScale << ", floor " << m_eventFloor
              << "; lone track: credit >= " << m_minLoneCredit << " and >= " << m_ambiguityRatio << " x its other hypotheses" << std::endl;
    std::cout << Name() << ":   collisions need a chain with INTT hits (MVTX-only chains only join them): "
              << (m_requireInttForCollision ? "yes" : "no") << std::endl;
    std::cout << Name() << ":   split collisions: "
              << (m_mergeSplitCollisions ? std::format("before pass 2, merged if |crossing difference| <= {}, |z difference| <= "
                                                       "min({} cm, {} sigma_z sqrt(1/n1 + 1/n2)) and chi2(one crossing) - "
                                                       "chi2(own crossings) < {} (dz sigma {} cm); tier 0 with >= {} tracks",
                                                       m_mergeMaxCrossings, m_mergeZWindow, m_mergeZNSigma, m_mergeMaxDeltaChi2, m_dzSigma,
                                                       m_mergeTier0MinTracks)
                                         : std::string("not merged"))
              << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCrossingResolver::process_event(PHCompositeNode* topNode)
{
  ++m_nEvents;
  m_decisions->Reset();
  auto* assembled = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_assembledNodeName);

  // ---- triggered: one crossing for every assembled track (as TpcCrossingFinder)
  if (m_triggeredMode)
  {
    const unsigned int n = assembled ? assembled->size() : 0;
    for (unsigned int i = 0; i < n; ++i)
    {
      const Tpc_AssembledTrack* t = assembled->get_track(i);
      if (!t)
      {
        continue;
      }
      auto* decision = new TpcCrossingDecisionv1();
      decision->set_assembled_track_id(t->get_track_id());
      decision->set_selected_crossing(m_triggeredCrossing);
      decision->set_selected_tier(0U);
      decision->set_number_of_available_crossings(1U);
      decision->set_number_of_tpc_valid_crossings(1U);
      decision->set_status(TpcCrossingStatus::SelectedByContainment);
      m_decisions->add_decision(decision);
      ++m_nTracks;
      ++m_nTier[0];
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto* trajectories = findNode::getClass<Tpc_CrossingTrajectoryContainer>(topNode, m_trajectoryNodeName);
  auto* seeds = findNode::getClass<SiHitSeedEvent>(topNode, m_siSeedNodeName);
  if (!trajectories)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_trajectoryNodeName << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // ---- tracks: TPC z0 of each of their crossing trajectories
  std::map<unsigned int, std::map<short, double>> tpcZ;  // assembled id -> crossing -> TPC z0
  std::map<unsigned int, std::map<short, int>> tpcIndex;  // assembled id -> crossing -> trajectory index
  for (unsigned int k = 0; k < trajectories->size(); ++k)
  {
    const Tpc_CrossingTrajectory* t = trajectories->get_trajectory(k);
    if (!t || t->get_fit_status() != 1 || !std::isfinite(t->get_z()))
    {
      continue;
    }
    tpcZ[t->get_source_assembled_track_id()][t->get_crossing()] = t->get_z();
    tpcIndex[t->get_source_assembled_track_id()][t->get_crossing()] = static_cast<int>(k);
  }

  // ---- 1. + 2. every trajectory takes its best chain; credit = match * quality
  std::vector<char> wasSelected;
  std::map<int, Hypothesis> bestOfTrajectory;  // trajectory index -> hypothesis
  const double sig2 = m_dzSigma * m_dzSigma;
  if (seeds)
  {
    wasSelected.resize(seeds->chains.size());
    for (std::size_t c = 0; c < seeds->chains.size(); ++c)
    {
      const SiHitChain& ch = seeds->chains[c];
      wasSelected[c] = ch.selected;
      if (ch.tpc_track_index < 0 || !std::isfinite(ch.tpc_chi2) || ch.tpc_chi2 > m_maxChainChi2 || !std::isfinite(ch.tpc_dz_offset))
      {
        continue;
      }
      const Tpc_CrossingTrajectory* t = trajectories->get_trajectory(static_cast<unsigned int>(ch.tpc_track_index));
      if (!t || t->get_fit_status() != 1 || !std::isfinite(t->get_z()) || std::fabs(t->get_z()) > m_maxVertexZ)
      {
        continue;
      }
      Hypothesis h;
      h.assembled = t->get_source_assembled_track_id();
      h.crossing = t->get_crossing();
      h.chain = static_cast<int>(c);
      h.dz = ch.tpc_dz_offset;
      h.tpcZ = t->get_z();
      h.z = h.tpcZ + h.dz;
      h.timed = ch.n_intt > 0;
      // match: z at this crossing and r*phi
      const double rphi = (ch.tpc_nres > 0 && std::isfinite(ch.tpc_chi2_rphi)) ? ch.tpc_chi2_rphi / ch.tpc_nres : 0.0;
      h.match = std::exp(-0.5 * h.dz * h.dz / sig2 - 0.5 * m_rphiWeight * rphi);
      // quality: layers, hardware delta-delta-phi, INTT time information
      double pull2 = 0;
      int npull = 0;
      for (const double p : ch.hw_ddphi_pull)
      {
        if (std::isfinite(p))
        {
          pull2 += p * p;
          ++npull;
        }
      }
      h.quality = std::pow(std::min(1.0, ch.n_layers() / 7.0), m_layerPower) *
                  std::exp(-0.5 * m_ddphiWeight * (npull > 0 ? pull2 / npull : 0.0)) * (ch.n_intt > 0 ? 1.0 : m_noInttFactor);
      h.credit = h.match * h.quality;
      auto [it, inserted] = bestOfTrajectory.try_emplace(ch.tpc_track_index, h);
      if (!inserted && h.credit > it->second.credit)
      {
        it->second = h;
      }
    }
  }

  // ---- 3. votes: credit * credit share of the track (at most one vote per track)
  std::map<unsigned int, std::vector<Hypothesis>> byTrack;
  for (const auto& [k, h] : bestOfTrajectory)
  {
    byTrack[h.assembled].push_back(h);
  }
  for (auto& [id, hyps] : byTrack)
  {
    double sum = 0;
    for (const auto& h : hyps)
    {
      sum += h.credit;
    }
    for (auto& h : hyps)
    {
      h.vote = sum > 0 ? h.credit * h.credit / sum : 0.0;
    }
    m_nHypotheses += hyps.size();
  }

  // ---- 4. pass 1: vertices of all hypotheses per crossing; the main vertex of a crossing is the
  // one with the largest vote sum (one collision per crossing), the others are the weak votes of
  // wrong hypotheses
  std::map<short, std::vector<Hypothesis*>> voters;
  for (auto& [id, hyps] : byTrack)
  {
    for (auto& h : hyps)
    {
      voters[h.crossing].push_back(&h);
    }
  }
  std::map<short, Vertex> main1;  // crossing -> main vertex (pass 1)
  unsigned int nCand1 = 0;
  for (auto& [crossing, list] : voters)
  {
    std::vector<std::pair<double, double>> zw;
    for (const Hypothesis* h : list)
    {
      zw.emplace_back(h->z, h->vote);
    }
    std::vector<int> index;
    const auto vtx = clusterZ(zw, m_clusterZWindow, index);
    nCand1 += vtx.size();
    const int m = mainVertex(vtx);
    main1[crossing] = vtx[m];
    for (std::size_t i = 0; i < list.size(); ++i)
    {
      list[i]->inMain = index[i] == m;
    }
  }

  // ---- 5. score: credit x leave-one-out support of the main vertex of its crossing
  const double sigH2 = m_siZSigma * m_siZSigma;
  for (auto& [id, hyps] : byTrack)
  {
    for (auto& h : hyps)
    {
      Vertex v = main1[h.crossing];
      unsigned int nOther = v.n;
      if (h.inMain)
      {
        v.sumW -= h.vote;
        v.sumWZ -= h.vote * h.z;
        nOther = nOther > 0 ? nOther - 1 : 0;
      }
      h.event = 0;
      if (v.sumW > 1e-6 && nOther > 0)
      {
        // track vs vertex of the other tracks: sigma^2 = sigma_track^2 (1 + 1 / n_other)
        const double d = h.z - v.mean();
        const double s2 = sigH2 * (1.0 + 1.0 / nOther);
        h.event = v.sumW / (v.sumW + m_supportScale) * std::exp(-0.5 * d * d / s2);
      }
      h.score = h.credit * (m_eventFloor + h.event);
    }
    std::sort(hyps.begin(), hyps.end(), [](const Hypothesis& a, const Hypothesis& b) { return a.score > b.score; });
  }

  // ---- ownership: best tracks first, each takes its best hypothesis with a free chain
  std::vector<unsigned int> order;
  for (const auto& [id, hyps] : byTrack)
  {
    if (!hyps.empty())
    {
      order.push_back(id);
    }
  }
  std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) { return byTrack[a].front().score > byTrack[b].front().score; });
  std::unordered_map<int, unsigned int> owner;  // blob -> assembled id
  auto chainFree = [&](const Hypothesis& h, unsigned int id)
  {
    unsigned int shared = 0;
    for (const int b : seeds->chains[h.chain].blob_ids)
    {
      const auto it = owner.find(b);
      shared += it != owner.end() && it->second != id;
    }
    return shared < m_minSharedBlobs;
  };
  auto take = [&](const Hypothesis& h, unsigned int id)
  {
    for (const int b : seeds->chains[h.chain].blob_ids)
    {
      owner.emplace(b, id);
    }
  };
  auto release = [&](const Hypothesis& h, unsigned int id)
  {
    for (const int b : seeds->chains[h.chain].blob_ids)
    {
      const auto it = owner.find(b);
      if (it != owner.end() && it->second == id)
      {
        owner.erase(it);
      }
    }
  };
  std::map<unsigned int, const Hypothesis*> winner;  // assembled id -> chosen hypothesis
  std::map<unsigned int, std::string> noChainReason;  // assembled id -> why it has no Si chain
  for (const unsigned int id : order)
  {
    for (const auto& h : byTrack[id])
    {
      if (chainFree(h, id))
      {
        winner[id] = &h;
        take(h, id);
        break;
      }
    }
    if (winner.count(id) == 0)
    {
      ++m_nLostByOwnership;
      noChainReason[id] = std::format("its {} Si chain(s) all taken by better tracks (ownership)", byTrack[id].size());
    }
  }

  // ---- 6. pass 2, one collision per crossing: vertices from the winners only (weight credit).
  //   supported: the winner is in the main vertex of its crossing together with another winner
  //   displaced : in a weaker vertex of its crossing -> moved to another hypothesis of the track
  //               that joins a supported main vertex, else no chain
  //   lone      : the only winner of its crossing -> moved like a displaced one if possible, else
  //               kept only if robust (credit >= setMinLoneCredit and >= setAmbiguityRatio x its
  //               other hypotheses, i.e. its other crossings are far away in z), else no chain
  // anchored: the main vertex of the crossing contains a winner with INTT hits. MVTX-only chains
  // have no time information (the MVTX integrates over many crossings, and a TPC-guided window
  // finds some MVTX-only chain near the TPC z of ANY crossing hypothesis), so they cannot make a
  // collision of their own, only join one that an INTT-timed chain establishes
  std::map<short, bool> anchored;
  struct SecondVertex
  {
    short crossing;
    double zMain, zSecond;
    unsigned int nTimed;
  };
  std::vector<SecondVertex> secondTimed;  // crossings with another vertex of >= 2 INTT-timed winners
  auto buildMain = [&](std::map<short, Vertex>& mainV, std::map<unsigned int, bool>& inMainV)
  {
    std::map<short, std::vector<unsigned int>> perCrossing;
    for (const auto& [id, h] : winner)
    {
      perCrossing[h->crossing].push_back(id);
    }
    mainV.clear();
    inMainV.clear();
    anchored.clear();
    secondTimed.clear();
    for (const auto& [crossing, ids] : perCrossing)
    {
      std::vector<std::pair<double, double>> zw;
      for (const unsigned int id : ids)
      {
        zw.emplace_back(winner[id]->z, winner[id]->credit);
      }
      std::vector<int> index;
      const auto vtx = clusterZ(zw, m_clusterZWindow, index);
      std::vector<unsigned int> timed(vtx.size(), 0);  // INTT-timed winners per vertex
      for (std::size_t i = 0; i < ids.size(); ++i)
      {
        timed[index[i]] += winner[ids[i]]->timed ? 1 : 0;
      }
      // main vertex: the largest one with an INTT-timed winner (the largest one if none or if
      // the INTT requirement is off)
      int m = -1;
      for (std::size_t j = 0; m_requireInttForCollision && j < vtx.size(); ++j)
      {
        if (timed[j] && (m < 0 || vtx[j].sumW > vtx[m].sumW))
        {
          m = static_cast<int>(j);
        }
      }
      if (m < 0)
      {
        m = mainVertex(vtx);
      }
      mainV[crossing] = vtx[m];
      anchored[crossing] = !m_requireInttForCollision || timed[m] > 0;
      for (std::size_t j = 0; j < vtx.size(); ++j)
      {
        if (static_cast<int>(j) != m && timed[j] >= 2)
        {
          secondTimed.push_back({crossing, vtx[m].mean(), vtx[j].mean(), timed[j]});
        }
      }
      for (std::size_t i = 0; i < ids.size(); ++i)
      {
        inMainV[ids[i]] = index[i] == m;
      }
    }
  };
  std::map<short, Vertex> main2;
  std::map<unsigned int, bool> inMain2;
  std::deque<Hypothesis> relabeled;                  // winners moved to the common crossing c* (stable addresses)
  std::map<unsigned int, unsigned int> mergedGroup;  // assembled id -> tracks of its merged collision
  buildMain(main2, inMain2);

  // ---- 6a. split collisions, one collision = one crossing also across crossings: a collision whose tracks were
  // split over neighbouring crossings (each track picks its crossing from its own dz, sigma ~0.3 cm
  // against a crossing step of ~0.8 - 1 cm) shows up as collisions at the same Si z in crossings
  // c, c +- 1, c +- 2. The Si z of a track does not depend on the crossing; only its TPC z (hence
  // dz) does, by one crossing step per crossing (taken from the trajectories, no step assumed).
  //   done BEFORE pass 2 drops displaced / lone tracks, so the tracks of a split collision are
  //   still there to be recovered
  //   candidates : main vertices of the winners, heaviest first; another one joins if
  //                |crossing difference| <= setMergeMaxCrossings and the Si vertices agree:
  //                |z difference| <= min(setMergeZWindow, setMergeZNSigma * sigma_z * sqrt(1/n1 + 1/n2)),
  //                sigma_z = setSiZSigma (single-track Si z at the beam)
  //   test       : chi2(c) = sum over the tracks of min((z_Si - z_TPC(c))^2 / sigma_dz^2, 9)
  //                (9 also for a track without a trajectory at c); joined if
  //                min_c chi2(c) - chi2(every track at its own crossing) < setMergeMaxDeltaChi2
  //                (two real collisions in neighbouring crossings cost ~ (step / sigma_dz)^2 per track)
  //   crossing   : c* = argmin chi2(c), i.e. the crossing of the combined dz (error sigma_dz / sqrt(n));
  //                the tracks of the group at other crossings move to c* with the same chain
  //                (dz recomputed at c*) unless their |dz| at c* is beyond 3 sigma_dz
  if (m_mergeSplitCollisions && seeds && !winner.empty())
  {
    buildMain(main2, inMain2);
    struct Collision
    {
      short crossing = 0;
      double z = 0, w = 0;
      std::vector<unsigned int> ids;
    };
    std::vector<Collision> cols;
    {
      std::map<short, Collision> byCrossing;
      for (const auto& [id, h] : winner)
      {
        if (inMain2[id])
        {
          byCrossing[h->crossing].ids.push_back(id);
        }
      }
      for (auto& [crossing, c] : byCrossing)
      {
        c.crossing = crossing;
        c.z = main2[crossing].mean();
        c.w = main2[crossing].sumW;
        cols.push_back(c);
      }
    }
    std::sort(cols.begin(), cols.end(), [](const Collision& a, const Collision& b) { return a.w > b.w; });

    constexpr double maxPull2 = 9.0;
    auto pull2 = [&](unsigned int id, short crossing)
    {
      const auto& zt = tpcZ.at(id);
      const auto f = zt.find(crossing);
      if (f == zt.end())
      {
        return maxPull2;
      }
      const double d = winner[id]->z - f->second;
      return std::min(d * d / sig2, maxPull2);
    };
    // best common crossing of a set of tracks: (crossing, chi2)
    auto bestCrossing = [&](const std::vector<unsigned int>& ids, short cmin, short cmax)
    {
      std::pair<short, double> best{cmin, std::numeric_limits<double>::max()};
      for (int c = cmin; c <= cmax; ++c)
      {
        double chi2 = 0;
        for (const unsigned int id : ids)
        {
          chi2 += pull2(id, static_cast<short>(c));
        }
        if (chi2 < best.second)
        {
          best = {static_cast<short>(c), chi2};
        }
      }
      return best;
    };
    auto ownChi2 = [&](const std::vector<unsigned int>& ids)
    {
      double chi2 = 0;
      for (const unsigned int id : ids)
      {
        chi2 += pull2(id, winner[id]->crossing);
      }
      return chi2;
    };

    std::vector<char> used(cols.size(), 0);
    for (std::size_t a = 0; a < cols.size(); ++a)
    {
      if (used[a])
      {
        continue;
      }
      used[a] = 1;
      std::vector<unsigned int> ids = cols[a].ids;
      short cmin = cols[a].crossing, cmax = cols[a].crossing;
      unsigned int nJoined = 0;
      for (std::size_t b = a + 1; b < cols.size(); ++b)
      {
        if (used[b] || std::abs(cols[b].crossing - cols[a].crossing) > m_mergeMaxCrossings ||
            std::fabs(cols[b].z - cols[a].z) >
                std::min(m_mergeZWindow, m_mergeZNSigma * m_siZSigma * std::sqrt(1.0 / cols[a].ids.size() + 1.0 / cols[b].ids.size())))
        {
          continue;
        }
        std::vector<unsigned int> trial = ids;
        trial.insert(trial.end(), cols[b].ids.begin(), cols[b].ids.end());
        const short lo = std::min(cmin, cols[b].crossing), hi = std::max(cmax, cols[b].crossing);
        const auto [c, chi2] = bestCrossing(trial, lo, hi);
        const double own = ownChi2(trial);
        if (Verbosity() > 1)
        {
          std::cout << Name() << ": split collision test: crossing " << cols[a].crossing << " (z " << cols[a].z << ", "
                    << cols[a].ids.size() << " tracks) + crossing " << cols[b].crossing << " (z " << cols[b].z << ", "
                    << cols[b].ids.size() << " tracks): best common crossing " << c << " chi2 " << chi2 << " vs own " << own
                    << (chi2 - own < m_mergeMaxDeltaChi2 ? " -> merged" : " -> kept apart") << std::endl;
        }
        if (chi2 - own < m_mergeMaxDeltaChi2)
        {
          ids = std::move(trial);
          cmin = lo;
          cmax = hi;
          used[b] = 1;
          ++nJoined;
        }
      }
      if (nJoined == 0)
      {
        continue;
      }
      m_nMergedCollisions += nJoined;
      const short cStar = bestCrossing(ids, cmin, cmax).first;
      for (const unsigned int id : ids)
      {
        const Hypothesis* w = winner[id];
        if (w->crossing != cStar)
        {
          if (pull2(id, cStar) >= maxPull2)
          {
            continue;  // no trajectory at c*, or |dz| beyond 3 sigma: stays where it is
          }
          Hypothesis r = *w;  // same chain, same Si z
          r.crossing = cStar;
          r.tpcZ = tpcZ.at(id).at(cStar);
          r.dz = r.z - r.tpcZ;
          const double f = std::exp(-0.5 * (r.dz * r.dz - w->dz * w->dz) / sig2);
          r.match *= f;
          r.credit *= f;
          relabeled.push_back(r);
          winner[id] = &relabeled.back();
          ++m_nMergeMoved;
        }
        mergedGroup[id] = static_cast<unsigned int>(ids.size());
      }
    }
  }
  buildMain(main2, inMain2);  // collisions after the merge, input of the classification
  m_nSecondTimed += secondTimed.size();
  if (Verbosity() > 1)
  {
    for (const auto& sv : secondTimed)
    {
      std::cout << Name() << ": crossing " << sv.crossing << " has a second vertex with " << sv.nTimed
                << " INTT-timed tracks at z " << sv.zSecond << " (main vertex z " << sv.zMain << "): suppressed by one collision per crossing" << std::endl;
    }
  }

  enum Kind
  {
    Supported = 0,
    Moved = 1,
    Lone = 2
  };
  std::map<unsigned int, int> kind;
  std::vector<unsigned int> byScore;
  for (const auto& [id, h] : winner)
  {
    byScore.push_back(id);
  }
  std::sort(byScore.begin(), byScore.end(), [&](unsigned int a, unsigned int b) { return winner[a]->score > winner[b]->score; });
  for (const unsigned int id : byScore)
  {
    const Hypothesis* w = winner[id];
    const Vertex& v = main2[w->crossing];
    if (inMain2[id] && v.n >= 2 && anchored[w->crossing])
    {
      kind[id] = Supported;
      continue;
    }
    // another hypothesis of the track joining a supported main vertex (without this track)
    const Hypothesis* alt = nullptr;
    for (const auto& h : byTrack[id])
    {
      if (&h == w)
      {
        continue;
      }
      const auto im = main2.find(h.crossing);
      if (im == main2.end() || im->second.n < 2 || !anchored[h.crossing] || std::fabs(h.z - im->second.mean()) > m_clusterZWindow ||
          !chainFree(h, id))
      {
        continue;
      }
      alt = &h;  // hypotheses are sorted by score: the first one found is the best
      break;
    }
    if (alt)
    {
      release(*w, id);
      winner[id] = alt;
      take(*alt, id);
      kind[id] = Moved;
      ++m_nMoved;
      continue;
    }
    if (inMain2[id] && v.n == 1)
    {
      double otherCredit = 0;
      for (const auto& h : byTrack[id])
      {
        if (&h != w)
        {
          otherCredit = std::max(otherCredit, h.credit);
        }
      }
      if ((w->timed || !m_requireInttForCollision) && w->credit >= m_minLoneCredit && w->credit >= m_ambiguityRatio * otherCredit)
      {
        kind[id] = Lone;
        ++m_nLoneKept;
        continue;
      }
    }
    // displaced, or a lone track that is not robust: no chain
    if (inMain2[id] && v.n == 1)
    {
      double otherCredit = 0;
      for (const auto& h : byTrack[id])
      {
        otherCredit = &h != w ? std::max(otherCredit, h.credit) : otherCredit;
      }
      if (m_requireInttForCollision && !w->timed)
      {
        noChainReason[id] = std::format("no INTT time: alone in crossing {} (z {:.2f}) with an MVTX-only chain", w->crossing, w->z);
        ++m_nUntimed;
      }
      else
      {
        noChainReason[id] = std::format("alone in crossing {} (z {:.2f}), not robust: credit {:.3f} (min {}), best other hypothesis {:.3f} (ratio {})",
                                        w->crossing, w->z, w->credit, m_minLoneCredit, otherCredit, m_ambiguityRatio);
      }
    }
    else if (inMain2[id] && !anchored[w->crossing])
    {
      noChainReason[id] = std::format("no INTT time: the collision of crossing {} (z {:.2f}, {} tracks) has only MVTX-only chains",
                                      w->crossing, v.mean(), v.n);
      ++m_nUntimed;
    }
    else
    {
      noChainReason[id] = std::format("displaced: crossing {} has a collision at z {:.2f} with {} tracks, its chain at z {:.2f}",
                                      w->crossing, v.mean(), v.n, w->z);
    }
    release(*w, id);
    winner.erase(id);
    ++m_nDropped;
  }

  buildMain(main2, inMain2);  // final collisions (statistics, TPC-only match)

  // ---- SI_HIT_SEED_EVENT: exactly the winning chains are selected, at their crossing
  if (seeds)
  {
    for (auto& ch : seeds->chains)
    {
      ch.selected = false;
    }
  }
  std::map<short, unsigned int> winnersPerCrossing;  // crossing -> Si-matched tracks
  for (const auto& [id, h] : winner)
  {
    SiHitChain& ch = seeds->chains[h->chain];
    ch.selected = true;
    ch.tpc_crossing = h->crossing;
    ch.tpc_dz_offset = h->dz;
    // trajectory of the decided crossing (differs from the chain's own one after a merge)
    const auto& ti = tpcIndex[h->assembled];
    if (const auto f = ti.find(h->crossing); f != ti.end())
    {
      ch.tpc_track_index = f->second;
    }
    m_nOtherChain += !wasSelected[h->chain];
    ++winnersPerCrossing[h->crossing];
  }
  auto allPerCrossing = winnersPerCrossing;  // + TPC-only matched tracks (below)

  // ---- decisions
  for (const auto& [id, zByCrossing] : tpcZ)
  {
    const auto& hyps = byTrack[id];
    auto* decision = new TpcCrossingDecisionv1();
    decision->set_assembled_track_id(id);
    decision->set_number_of_available_crossings(static_cast<unsigned short>(zByCrossing.size()));
    decision->set_number_of_tpc_valid_crossings(static_cast<unsigned short>(hyps.size()));

    const auto iw = winner.find(id);
    if (iw != winner.end())
    {
      const Hypothesis& best = *iw->second;
      const int k = kind[id];
      double second = 0, secondDz = std::numeric_limits<double>::quiet_NaN();
      // merged collision: the same chain at another crossing is not a competitor (the crossing was
      // decided by the combined dz of the collision)
      const auto ig = mergedGroup.find(id);
      const bool merged = ig != mergedGroup.end();
      for (const auto& h : hyps)
      {
        if (h.crossing != best.crossing && !(merged && h.chain == best.chain) && h.score > second)
        {
          second = h.score;
          secondDz = h.dz;
        }
      }
      const bool clear = second <= 0 || best.score >= m_ambiguityRatio * second;
      const unsigned char tier = (k == Supported && clear && (!merged || ig->second >= m_mergeTier0MinTracks)) ? 0 : 1;
      const auto bestCredit = std::max_element(hyps.begin(), hyps.end(), [](const Hypothesis& a, const Hypothesis& b) { return a.credit < b.credit; });
      m_nChangedByEvent += bestCredit->crossing != best.crossing;

      decision->set_selected_crossing(best.crossing);
      decision->set_selected_tier(tier);
      decision->set_selected_score(static_cast<float>(best.score));
      decision->set_tpc_z0(static_cast<float>(best.tpcZ));
      decision->set_silicon_vertex_z(static_cast<float>(k == Lone ? best.z : main2[best.crossing].mean()));
      decision->set_delta_z(static_cast<float>(best.dz));
      decision->set_best_abs_delta_z(static_cast<float>(std::fabs(best.dz)));
      decision->set_second_best_abs_delta_z(static_cast<float>(std::fabs(secondDz)));
      decision->set_number_of_vertex_compatible_crossings(
          static_cast<unsigned short>(std::count_if(hyps.begin(), hyps.end(), [](const Hypothesis& h) { return h.event > 0.5; })));
      decision->set_status(k == Lone ? TpcCrossingStatus::SelectedByVertexLoose
                                     : (tier == 0 ? TpcCrossingStatus::SelectedByVertex : TpcCrossingStatus::SelectedByVertexAmbiguous));
      ++m_nTier[tier];
      if (tier == 0)
      {
        m_dzTier0.push_back(best.dz);
      }

      if (Verbosity() > 1)
      {
        static const char* kindName[3] = {"supported", "moved", "lone"};
        std::cout << Name() << ": assembled " << id << " -> crossing " << best.crossing << " (" << kindName[k] << ", "
                  << main2[best.crossing].n << " tracks at z " << main2[best.crossing].mean() << ") chain "
                  << seeds->chains[best.chain].id << (wasSelected[best.chain] ? "" : " [not rank 0]") << " tier "
                  << static_cast<int>(tier) << " score " << best.score << " (other crossing " << second << "), hypotheses:";
        for (const auto& h : hyps)
        {
          std::cout << " [c " << h.crossing << " chain " << seeds->chains[h.chain].id << " L " << seeds->chains[h.chain].n_layers()
                    << " dz " << h.dz << " z " << h.z << " match " << h.match << " quality " << h.quality << " credit " << h.credit
                    << " vote " << h.vote << (h.inMain ? " main" : " weak") << " vertex z " << main1[h.crossing].mean() << " W "
                    << main1[h.crossing].sumW << " ev " << h.event << " score " << h.score << "]";
        }
        std::cout << std::endl;
      }
    }
    else
    {
      // no Si chain: final collisions of the track's crossings, tested with the TPC z0
      int nMatch = 0;
      short bestC = 0;
      double bestChi2 = std::numeric_limits<double>::max(), bestZ = 0;
      for (const auto& [crossing, zT] : zByCrossing)
      {
        const auto im = main2.find(crossing);
        if (im == main2.end() || im->second.n < 2 || im->second.sumW < m_minSupportTpcOnly)
        {
          continue;
        }
        const double d = zT - im->second.mean();
        const double chi2 = d * d / (m_tpcZSigma * m_tpcZSigma + sigH2 / im->second.n);  // vertex error sigma / sqrt(n)
        nMatch += chi2 < 4.0;
        if (chi2 < bestChi2)
        {
          bestChi2 = chi2;
          bestC = crossing;
          bestZ = im->second.mean();
        }
      }
      const bool matched = nMatch == 1;
      if (!matched)
      {
        // not resolved (tier 2): crossing of the trajectory closest to z = 0
        bestC = std::min_element(zByCrossing.begin(), zByCrossing.end(), [](const auto& a, const auto& b)
                                 { return std::fabs(a.second) < std::fabs(b.second); })
                    ->first;
      }
      const unsigned char tier = matched ? 1 : 2;
      decision->set_selected_crossing(bestC);
      decision->set_selected_tier(tier);
      decision->set_tpc_z0(static_cast<float>(zByCrossing.at(bestC)));
      if (matched)
      {
        decision->set_silicon_vertex_z(static_cast<float>(bestZ));
        decision->set_delta_z(static_cast<float>(bestZ - zByCrossing.at(bestC)));
        ++allPerCrossing[bestC];
      }
      decision->set_status(matched ? TpcCrossingStatus::SelectedByVertexLoose : TpcCrossingStatus::AmbiguousWithoutVertex);
      ++m_nTier[tier];
      m_nTpcOnly += matched;

      // why this track has no Si chain (End summary, Verbosity > 1 per track)
      std::string reason, category;
      const auto ir = noChainReason.find(id);
      if (ir != noChainReason.end())
      {
        reason = ir->second;
        category = reason.rfind("its ", 0) == 0     ? "lost by ownership"
                   : reason.rfind("no INTT", 0) == 0 ? "MVTX-only chain without an INTT-timed collision"
                   : reason.rfind("alone", 0) == 0   ? "alone in its crossing, not robust"
                                                     : "displaced by another collision in its crossing";
      }
      else
      {
        // no hypothesis: the SiHitSeedReco status of its crossing trajectories
        std::map<int, unsigned int> st;
        unsigned int nChainBad = 0;
        if (seeds)
        {
          for (const auto& info : seeds->tpc_seeds)
          {
            if (info.tpc_assembled_track_id == id)
            {
              ++st[info.status];
            }
          }
          for (const auto& ch : seeds->chains)
          {
            nChainBad += ch.tpc_assembled_track_id == id;
          }
        }
        static const char* stName[7] = {"chained", "too few TPC points", "TPC fit failed", "no prediction on the seed layers",
                                        "no hit in the anchor window", "no chain passing the layer requirements", "best chain fails the agreement cut"};
        static const int rankOf[7] = {6, 0, 1, 2, 3, 4, 5};  // most advanced status wins
        int best = -1;
        reason = std::format("no hypothesis: {} trajectories (", zByCrossing.size());
        for (const auto& [s, n] : st)
        {
          reason += std::format(" {} {},", stName[std::clamp(s, 0, 6)], n);
          best = (best < 0 || rankOf[std::clamp(s, 0, 6)] > rankOf[std::clamp(best, 0, 6)]) ? s : best;
        }
        reason += std::format(" ), {} chains, none with chi2/ndf <= {} and |TPC z0| <= {}", nChainBad, m_maxChainChi2, m_maxVertexZ);
        category = best < 0 ? "no hypothesis: no Si seeding input" : std::string("no hypothesis, best: ") + stName[std::clamp(best, 0, 6)];
      }
      ++m_noChainCategory[category];
      ++m_noChainTier[{category, tier}];
      if (Verbosity() > 1)
      {
        std::cout << Name() << ": assembled " << id << " -> no Si chain (" << reason << "); "
                  << (matched ? std::format("TPC-only match to the collision of crossing {} at z {:.2f} (tier 1)", bestC, bestZ)
                              : std::string("no unique TPC-only match (tier 2)"))
                  << std::endl;
      }
    }
    m_decisions->add_decision(decision);
    ++m_nTracks;
  }

  for (const auto& [crossing, n] : winnersPerCrossing)
  {
    ++m_multiplicity[n];
  }
  for (const auto& [crossing, n] : allPerCrossing)
  {
    ++m_multiplicityAll[n];
  }

  // assembled tracks without any crossing trajectory: no usable crossing
  const unsigned int n = assembled ? assembled->size() : 0;
  for (unsigned int i = 0; i < n; ++i)
  {
    const Tpc_AssembledTrack* t = assembled->get_track(i);
    if (!t || tpcZ.count(t->get_track_id()) != 0)
    {
      continue;
    }
    auto* decision = new TpcCrossingDecisionv1();
    decision->set_assembled_track_id(t->get_track_id());
    decision->set_selected_tier(NoTier);
    decision->set_status(TpcCrossingStatus::NoValidCrossing);
    m_decisions->add_decision(decision);
    ++m_nTier[3];
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << m_nEvents << " tracks " << tpcZ.size() << " (with a hypothesis " << order.size()
              << ", with a chain at the end " << winner.size() << "), hypotheses " << bestOfTrajectory.size()
              << " in " << voters.size() << " crossings (pass-1 vertices " << nCand1 << "), collisions " << winnersPerCrossing.size()
              << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCrossingResolver::End(PHCompositeNode* /*topNode*/)
{
  std::cout << Name() << "::End - " << m_nEvents << " events, decisions " << m_nTracks << ": tier 0 " << m_nTier[0]
            << ", tier 1 " << m_nTier[1] << " (lone kept " << m_nLoneKept << ", TPC-only match " << m_nTpcOnly << "), tier 2 "
            << m_nTier[2] << "; no crossing trajectory " << m_nTier[3];
  if (!m_triggeredMode)
  {
    std::cout << std::endl
              << Name() << "::End - hypotheses (crossing trajectories with a chain) " << m_nHypotheses
              << ", winner not the highest-credit hypothesis " << m_nChangedByEvent << ", chain other than SiHitSeedReco's rank 0 "
              << m_nOtherChain << ", left without a chain by the ownership " << m_nLostByOwnership << std::endl
              << Name() << "::End - one collision per crossing: moved to a supported vertex " << m_nMoved
              << ", lone kept (robust) " << m_nLoneKept << ", dropped (displaced / lone not robust) " << m_nDropped << std::endl
              << Name() << "::End - dropped because no INTT-timed chain in their collision (MVTX-only) " << m_nUntimed << std::endl
              << Name() << "::End - split collisions merged into the collision of a neighbouring crossing " << m_nMergedCollisions
              << ", tracks moved to the common crossing " << m_nMergeMoved << std::endl
              << Name() << "::End - crossings with a second vertex of >= 2 INTT-timed tracks (suppressed by one collision per crossing) "
              << m_nSecondTimed << std::endl
              << Name() << "::End - tracks without a Si chain, why (tier 1 = TPC-only match, drawn grey in SiTpc_EventDisplay / tier 2 = not used):" << std::endl;
    for (const auto& [cat, n] : m_noChainCategory)
    {
      const auto t1 = m_noChainTier.find({cat, 1});
      std::cout << Name() << "::End -   " << cat << ": " << n << " (tier 1 " << (t1 != m_noChainTier.end() ? t1->second : 0UL) << ")" << std::endl;
    }
    std::cout
              << Name() << "::End - Si-matched tracks per collision (crossing):";
    for (const auto& [mult, count] : m_multiplicity)
    {
      std::cout << " " << mult << ":" << count;
    }
    std::cout << std::endl
              << Name() << "::End - + TPC-only matched tracks per collision:";
    for (const auto& [mult, count] : m_multiplicityAll)
    {
      std::cout << " " << mult << ":" << count;
    }
    if (!m_dzTier0.empty())
    {
      std::vector<double> d = m_dzTier0;
      const double med = median(d);
      for (auto& x : d)
      {
        x = std::fabs(x - med);
      }
      std::cout << std::endl
                << Name() << "::End - Si - TPC dz of the tier-0 winners: median " << med << " cm, robust sigma "
                << 1.4826 * median(d) << " cm (setDzSigma in use " << m_dzSigma << ")";
    }
  }
  std::cout << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}