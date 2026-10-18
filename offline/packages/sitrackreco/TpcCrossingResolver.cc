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
#include <cmath>
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
              << " cm, sigma z " << m_siZSigma << " cm, support scale " << m_supportScale << ", floor " << m_eventFloor
              << "; lone track: credit >= " << m_minLoneCredit << " and >= " << m_ambiguityRatio << " x its other hypotheses" << std::endl;
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
  for (unsigned int k = 0; k < trajectories->size(); ++k)
  {
    const Tpc_CrossingTrajectory* t = trajectories->get_trajectory(k);
    if (!t || t->get_fit_status() != 1 || !std::isfinite(t->get_z()))
    {
      continue;
    }
    tpcZ[t->get_source_assembled_track_id()][t->get_crossing()] = t->get_z();
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
      if (h.inMain)
      {
        v.sumW -= h.vote;
        v.sumWZ -= h.vote * h.z;
      }
      h.event = 0;
      if (v.sumW > 1e-6)
      {
        const double d = h.z - v.mean();
        h.event = v.sumW / (v.sumW + m_supportScale) * std::exp(-0.5 * d * d / sigH2);
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
    m_nLostByOwnership += winner.count(id) == 0;
  }

  // ---- 6. pass 2, one collision per crossing: vertices from the winners only (weight credit).
  //   supported: the winner is in the main vertex of its crossing together with another winner
  //   displaced : in a weaker vertex of its crossing -> moved to another hypothesis of the track
  //               that joins a supported main vertex, else no chain
  //   lone      : the only winner of its crossing -> moved like a displaced one if possible, else
  //               kept only if robust (credit >= setMinLoneCredit and >= setAmbiguityRatio x its
  //               other hypotheses, i.e. its other crossings are far away in z), else no chain
  auto buildMain = [&](std::map<short, Vertex>& mainV, std::map<unsigned int, bool>& inMainV)
  {
    std::map<short, std::vector<unsigned int>> perCrossing;
    for (const auto& [id, h] : winner)
    {
      perCrossing[h->crossing].push_back(id);
    }
    mainV.clear();
    inMainV.clear();
    for (const auto& [crossing, ids] : perCrossing)
    {
      std::vector<std::pair<double, double>> zw;
      for (const unsigned int id : ids)
      {
        zw.emplace_back(winner[id]->z, winner[id]->credit);
      }
      std::vector<int> index;
      const auto vtx = clusterZ(zw, m_clusterZWindow, index);
      const int m = mainVertex(vtx);
      mainV[crossing] = vtx[m];
      for (std::size_t i = 0; i < ids.size(); ++i)
      {
        inMainV[ids[i]] = index[i] == m;
      }
    }
  };
  std::map<short, Vertex> main2;
  std::map<unsigned int, bool> inMain2;
  buildMain(main2, inMain2);

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
    if (inMain2[id] && v.n >= 2)
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
      if (im == main2.end() || im->second.n < 2 || std::fabs(h.z - im->second.mean()) > m_clusterZWindow || !chainFree(h, id))
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
      if (w->credit >= m_minLoneCredit && w->credit >= m_ambiguityRatio * otherCredit)
      {
        kind[id] = Lone;
        ++m_nLoneKept;
        continue;
      }
    }
    // displaced, or a lone track that is not robust: no chain
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
      for (const auto& h : hyps)
      {
        if (h.crossing != best.crossing && h.score > second)
        {
          second = h.score;
          secondDz = h.dz;
        }
      }
      const bool clear = second <= 0 || best.score >= m_ambiguityRatio * second;
      const unsigned char tier = (k == Supported && clear) ? 0 : 1;
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
        const double chi2 = d * d / (m_tpcZSigma * m_tpcZSigma + sigH2);
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
