#include "SiTpcOnlyCrossing.h"

#include "SiHitSeedData.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/TpcCrossingDecisionv1.h>
#include <tpctrackreco/Tpc_AssembledTrack.h>
#include <tpctrackreco/Tpc_AssembledTrackContainer.h>

#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertexMap_v1.h>
#include <globalvertex/SvtxVertex_v3.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace
{
  bool isSelectedStatus(unsigned char s)
  {
    using S = TpcCrossingStatus;
    return s == static_cast<unsigned char>(S::SelectedByVertex) || s == static_cast<unsigned char>(S::SelectedByVertexAmbiguous) ||
           s == static_cast<unsigned char>(S::SelectedByVertexLoose) || s == static_cast<unsigned char>(S::SelectedByContainment) ||
           s == static_cast<unsigned char>(S::SelectedByContainmentAmbiguous);
  }

  // assembled ids of the tracks with a selected Si chain
  std::set<unsigned int> siMatchedTracks(const SiHitSeedEvent* seeds)
  {
    std::set<unsigned int> out;
    if (seeds)
    {
      for (const auto& ch : seeds->chains)
      {
        if (ch.selected)
        {
          out.insert(ch.tpc_assembled_track_id);
        }
      }
    }
    return out;
  }
}  // namespace

//____________________________________________________________________________________________________
SiTpcOnlyCrossing::SiTpcOnlyCrossing(Mode mode, const std::string& name)
  : SubsysReco(name)
  , m_mode(mode)
{
}

//____________________________________________________________________________________________________
int SiTpcOnlyCrossing::InitRun(PHCompositeNode* topNode)
{
  if (m_mode == Mode::MakeVertices)
  {
    // output vertex map, created here so that TpcCrossingFinder finds it in its InitRun
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
    if (!findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName))
    {
      svtxNode->addNode(new PHIODataNode<PHObject>(new SvtxVertexMap_v1, m_vertexMapName, "PHObject"));
    }
    std::cout << Name() << "::InitRun - collisions of the Si-matched tracks (" << m_decisionName << " tier <= " << m_maxCollisionTier
              << ", >= " << m_minCollisionTracks << " tracks within " << m_clusterZWindow << " cm) -> " << m_vertexMapName << ", vertex z "
              << (m_vertexZ == VertexZ::TpcZ0Mean ? "mean TPC z0 of its tracks" : "mean Si z of its tracks") << std::endl;
  }
  else
  {
    std::cout << Name() << "::InitRun - tracks without a Si chain take the " << m_finderName << " decision in " << m_decisionName
              << ": vertex-matched -> tier " << m_vertexTier << (m_rejectAmbiguous ? " (ambiguous -> tier 2)" : "")
              << ", otherwise tier 2" << (m_addMissing ? "; tracks without a resolver decision are added" : "") << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcOnlyCrossing::process_event(PHCompositeNode* topNode)
{
  ++m_nEvents;
  return m_mode == Mode::MakeVertices ? makeVertices(topNode) : mergeDecisions(topNode);
}

//____________________________________________________________________________________________________
int SiTpcOnlyCrossing::makeVertices(PHCompositeNode* topNode)
{
  auto* vertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);
  if (!vertexMap)
  {
    std::cout << Name() << ": missing " << m_vertexMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  vertexMap->clear();

  const auto* decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionName);
  const auto* seeds = findNode::getClass<SiHitSeedEvent>(topNode, m_siSeedName);
  if (!decisions || !seeds)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_decisionName << " or " << m_siSeedName << " - no collisions" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Si-matched tracks per crossing: (Si z, TPC z0, assembled id)
  struct Member
  {
    double zSi = 0, zTpc = 0;
    unsigned int id = 0;
  };
  std::map<short, std::vector<Member>> byCrossing;
  for (const unsigned int id : siMatchedTracks(seeds))
  {
    const TpcCrossingDecision* d = decisions->get_decision(id);
    if (!d || d->get_selected_crossing() == SHRT_MAX || static_cast<unsigned int>(d->get_selected_tier()) > m_maxCollisionTier)
    {
      continue;
    }
    const double z0 = d->get_tpc_z0();
    const double dz = d->get_delta_z();
    if (!std::isfinite(z0) || !std::isfinite(dz) || !std::isfinite(d->get_silicon_vertex_z()))  // untimed lone: no collision
    {
      continue;
    }
    byCrossing[d->get_selected_crossing()].push_back({z0 + dz, z0, id});
  }

  // per crossing: 1-D clustering in Si z (sorted, a track joins if within the window of the running mean)
  for (auto& [crossing, members] : byCrossing)
  {
    std::sort(members.begin(), members.end(), [](const Member& a, const Member& b) { return a.zSi < b.zSi; });
    std::vector<std::vector<const Member*>> clusters;
    double sum = 0;
    for (const auto& m : members)
    {
      if (clusters.empty() || std::abs(m.zSi - sum / clusters.back().size()) > m_clusterZWindow)
      {
        clusters.emplace_back();
        sum = 0;
      }
      clusters.back().push_back(&m);
      sum += m.zSi;
    }
    for (const auto& c : clusters)
    {
      if (c.size() < m_minCollisionTracks)
      {
        continue;
      }
      double zSi = 0, zTpc = 0;
      for (const auto* m : c)
      {
        zSi += m->zSi;
        zTpc += m->zTpc;
      }
      zSi /= c.size();
      zTpc /= c.size();
      auto v = std::make_unique<SvtxVertex_v3>();
      v->set_x(0);
      v->set_y(0);
      v->set_z(static_cast<float>(m_vertexZ == VertexZ::TpcZ0Mean ? zTpc : zSi));
      v->set_t(0);
      v->set_chisq(0);
      v->set_ndof(0);
      v->set_beam_crossing(crossing);
      for (unsigned int a = 0; a < 3; ++a)
      {
        for (unsigned int b = 0; b < 3; ++b)
        {
          v->set_error(a, b, 0);
        }
      }
      v->set_error(2, 2, static_cast<float>(m_trackZSigma * m_trackZSigma / c.size()));
      for (const auto* m : c)
      {
        v->insert_track(m->id);
      }
      vertexMap->insert(v.release());
      ++m_nVertices;
      m_nVertexTracks += c.size();
      m_sumSiMinusTpc += (zSi - zTpc) * c.size();
      m_nSiMinusTpc += c.size();
      if (Verbosity() > 1)
      {
        std::cout << Name() << ": crossing " << crossing << " collision " << c.size() << " tracks, Si z " << zSi << " TPC z0 " << zTpc << std::endl;
      }
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": " << vertexMap->size() << " collisions -> " << m_vertexMapName << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcOnlyCrossing::mergeDecisions(PHCompositeNode* topNode)
{
  auto* decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionName);
  const auto* finder = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_finderName);
  const auto* seeds = findNode::getClass<SiHitSeedEvent>(topNode, m_siSeedName);
  const auto* assembled = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_assembledName);
  if (!decisions || !finder || !assembled)
  {
    std::cout << Name() << ": missing " << m_decisionName << ", " << m_finderName << " or " << m_assembledName << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  const std::set<unsigned int> withSi = siMatchedTracks(seeds);

  for (unsigned int i = 0; i < assembled->size(); ++i)
  {
    const Tpc_AssembledTrack* t = assembled->get_track(i);
    if (!t)
    {
      continue;
    }
    const unsigned int id = t->get_track_id();
    if (withSi.count(id) != 0)
    {
      continue;  // Si-matched: the resolver decision stays
    }
    TpcCrossingDecision* d = decisions->get_decision(id);
    if (!d && !m_addMissing)
    {
      continue;
    }
    ++m_nTpcOnly;
    const TpcCrossingDecision* f = finder->get_decision(id);
    if (!f)
    {
      ++m_nNoFinderDecision;
      continue;
    }
    const int oldTier = d ? static_cast<int>(d->get_selected_tier()) : -1;
    const short oldCrossing = d ? d->get_selected_crossing() : static_cast<short>(SHRT_MAX);
    const bool valid = isSelectedStatus(f->get_status()) && f->get_selected_crossing() != SHRT_MAX;
    const bool ambiguous = f->get_status() == static_cast<unsigned char>(TpcCrossingStatus::SelectedByVertexAmbiguous);
    unsigned char tier = 2;
    if (valid && f->get_selected_tier() == 0 && !(m_rejectAmbiguous && ambiguous))
    {
      tier = static_cast<unsigned char>(m_vertexTier);
    }

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
      ++m_nAdded;
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
    }
    d->set_selected_tier(tier);
    m_nChangedCrossing += valid && oldCrossing != SHRT_MAX && oldCrossing != f->get_selected_crossing();
    ++m_tierChange[{oldTier, static_cast<int>(tier)}];

    if (Verbosity() > 1)
    {
      std::cout << Name() << ": assembled " << id << " (no Si chain): resolver crossing " << oldCrossing << " tier " << oldTier << " -> "
                << (valid ? "crossing " + std::to_string(f->get_selected_crossing()) : std::string("no valid crossing")) << " tier "
                << static_cast<int>(tier) << " (finder tier " << static_cast<int>(f->get_selected_tier()) << ", status "
                << static_cast<int>(f->get_status()) << ", TPC z0 " << f->get_tpc_z0() << ", vertex z " << f->get_silicon_vertex_z()
                << ", dz " << f->get_delta_z() << ", 2nd best |dz| " << f->get_second_best_abs_delta_z() << ")" << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcOnlyCrossing::End(PHCompositeNode* /*topNode*/)
{
  if (m_mode == Mode::MakeVertices)
  {
    std::cout << Name() << "::End - " << m_nEvents << " events, " << m_nVertices << " collisions with " << m_nVertexTracks
              << " Si-matched tracks; mean Si z - TPC z0 of their tracks "
              << (m_nSiMinusTpc > 0 ? m_sumSiMinusTpc / m_nSiMinusTpc : std::numeric_limits<double>::quiet_NaN()) << " cm" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  std::cout << Name() << "::End - " << m_nEvents << " events, " << m_nTpcOnly << " tracks without a Si chain, no finder decision "
            << m_nNoFinderDecision << ", added " << m_nAdded << ", crossing changed " << m_nChangedCrossing << std::endl
            << Name() << "::End - tier resolver -> new (-1 = no resolver decision):";
  for (const auto& [k, n] : m_tierChange)
  {
    std::cout << " " << k.first << "->" << k.second << ": " << n;
  }
  std::cout << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
