#include "TpcSiliconCrossingRefiner.h"

#include "Full_PolyTrack.h"
#include "Full_PolyTrackContainer.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include "TpcCrossingDecisionContainerv1.h"
#include "TpcCrossingDecisionv1.h"
#include "Tpc_FittingTools.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

TpcSiliconCrossingRefiner::TpcSiliconCrossingRefiner(const std::string& name)
  : SubsysReco(name)
{
}

int TpcSiliconCrossingRefiner::InitRun(PHCompositeNode* topNode)
{
  return getNodes(topNode) == Fun4AllReturnCodes::EVENT_OK &&
                 createNodes(topNode) == Fun4AllReturnCodes::EVENT_OK
             ? Fun4AllReturnCodes::EVENT_OK
             : Fun4AllReturnCodes::ABORTRUN;
}

int TpcSiliconCrossingRefiner::getNodes(PHCompositeNode* topNode)
{
  m_fullTracks = findNode::getClass<Full_PolyTrackContainer>(topNode, m_fullTrackNodeName);
  m_tpcTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName);
  m_inputCrossings = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_inputCrossingNodeName);
  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_fullTracks || !m_tpcTracks || !m_inputCrossings || !m_geometry)
  {
    std::cerr << Name() << "::getNodes - missing "
              << (!m_fullTracks ? m_fullTrackNodeName
                                : (!m_tpcTracks ? m_tpcTrackNodeName
                                                : (!m_inputCrossings ? m_inputCrossingNodeName : "ActsGeometry")))
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcSiliconCrossingRefiner::createNodes(PHCompositeNode* topNode)
{
  if (m_outputCrossingNodeName == m_inputCrossingNodeName)
  {
    std::cerr << Name() << "::createNodes - input and output crossing nodes must differ" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  m_outputCrossings = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_outputCrossingNodeName);
  if (!m_outputCrossings)
  {
    m_outputCrossings = new TpcCrossingDecisionContainerv1();
    dstNode->addNode(new PHIODataNode<PHObject>(m_outputCrossings, m_outputCrossingNodeName, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

const Tpc_PolyTrack* TpcSiliconCrossingRefiner::findTpcTrack(const unsigned int track_id) const
{
  for (unsigned int i = 0; m_tpcTracks && i < m_tpcTracks->size(); ++i)
  {
    const Tpc_PolyTrack* track = m_tpcTracks->get_track(i);
    if (track && track->get_track_id() == track_id)
    {
      return track;
    }
  }
  return nullptr;
}

bool TpcSiliconCrossingRefiner::fitSiliconZ0(const Full_PolyTrack& track, double& z0) const
{
  std::vector<Tpc_FittingTools::FitPoint> points;
  points.reserve(track.size_silicon_states());
  for (unsigned int i = 0; i < track.size_silicon_states(); ++i)
  {
    const double x = track.get_state_x(i);
    const double y = track.get_state_y(i);
    const double z = track.get_state_z(i);
    const double r = std::hypot(x, y);
    if (std::isfinite(r) && std::isfinite(z))
    {
      points.emplace_back(r, z);
    }
  }
  if (points.size() < 2U)
  {
    return false;
  }

  const Tpc_FittingTools::LineFit fit = Tpc_FittingTools::fitLine(points);
  if (!fit.ok || !std::isfinite(fit.intercept))
  {
    return false;
  }
  z0 = fit.intercept;
  return true;
}

int TpcSiliconCrossingRefiner::process_event(PHCompositeNode* topNode)
{
  if ((!m_fullTracks || !m_tpcTracks || !m_inputCrossings || !m_outputCrossings || !m_geometry) &&
      (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
       createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK))
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  m_outputCrossings->Reset();

  const double crossing_z = m_geometry->get_drift_velocity() * m_crossingPeriodNs;
  if (!std::isfinite(crossing_z) || crossing_z <= 0.0)
  {
    std::cerr << Name() << "::process_event - invalid crossing z displacement" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  unsigned int refined = 0U;
  for (unsigned int i = 0; i < m_fullTracks->size(); ++i)
  {
    const Full_PolyTrack* full = m_fullTracks->get_track(i);
    if (!full)
    {
      continue;
    }
    const Tpc_PolyTrack* tpc = findTpcTrack(full->get_tpc_poly_track_id());
    const TpcCrossingDecision* input = m_inputCrossings->get_decision(full->get_source_assembled_track_id());
    if (!tpc || !input || tpc->size_cluster_keys() == 0U)
    {
      continue;
    }

    double silicon_z0 = 0.0;
    if (!fitSiliconZ0(*full, silicon_z0))
    {
      continue;
    }
    const double tpc_z0 = tpc->get_seed_z0();
    const short initial_crossing = input->get_selected_crossing();
    const unsigned int side = TpcDefs::getSide(tpc->get_cluster_key(0));
    if (!std::isfinite(tpc_z0) || side > 1U ||
        initial_crossing == std::numeric_limits<short>::max())
    {
      continue;
    }

    const double direction = side == 0U ? -1.0 : 1.0;
    auto* output = new TpcCrossingDecisionv1();
    output->set_assembled_track_id(full->get_source_assembled_track_id());

    struct RankedCandidate
    {
      TpcCrossingCandidate candidate;
      double delta_z{0.0};
      double abs_delta_z{0.0};
      double score{0.0};
    };
    std::vector<RankedCandidate> plausible_candidates;
    unsigned short ncompatible = 0U;
    for (unsigned int icandidate = 0; icandidate < input->get_number_of_candidates(); ++icandidate)
    {
      const TpcCrossingCandidate* finder_candidate = input->get_candidate(icandidate);
      if (!finder_candidate)
      {
        continue;
      }

      TpcCrossingCandidate candidate = *finder_candidate;
      candidate.is_selected = false;
      candidate.candidate_qa_bits &= ~IsSelected;
      const int candidate_delta = static_cast<int>(candidate.crossing) - static_cast<int>(initial_crossing);
      const double corrected_tpc_z0 = tpc_z0 + direction * static_cast<double>(candidate_delta) * crossing_z;
      const double delta_z = silicon_z0 - corrected_tpc_z0;
      const double abs_delta_z = std::fabs(delta_z);

      candidate.has_silicon_vertex = true;
      candidate.vertex_compatible = abs_delta_z <= m_maximumAbsDeltaZ;
      candidate.tpc_z_at_r0 = static_cast<float>(corrected_tpc_z0);
      candidate.closest_vertex_z = static_cast<float>(silicon_z0);
      candidate.closest_vertex_delta_z = static_cast<float>(delta_z);
      candidate.closest_vertex_abs_delta_z = static_cast<float>(abs_delta_z);
      candidate.candidate_qa_bits |= HasSiliconVertex;
      if (candidate.vertex_compatible)
      {
        candidate.candidate_qa_bits |= PassesVertexDz;
      }

      if (candidate.was_tested && candidate.tpc_valid &&
          candidate.confidence_tier != std::numeric_limits<unsigned char>::max() &&
          std::isfinite(candidate.confidence_score))
      {
        ncompatible += candidate.vertex_compatible ? 1U : 0U;
        plausible_candidates.push_back(
            {candidate, delta_z, abs_delta_z,
             abs_delta_z + m_crossingFinderScoreWeight * candidate.confidence_score});
      }
    }

    if (plausible_candidates.empty())
    {
      delete output;
      continue;
    }

    std::sort(plausible_candidates.begin(), plausible_candidates.end(),
              [](const RankedCandidate& lhs, const RankedCandidate& rhs) {
                if (lhs.candidate.vertex_compatible != rhs.candidate.vertex_compatible)
                {
                  return lhs.candidate.vertex_compatible;
                }
                if (lhs.candidate.confidence_tier != rhs.candidate.confidence_tier)
                {
                  return lhs.candidate.confidence_tier < rhs.candidate.confidence_tier;
                }
                if (lhs.score != rhs.score)
                {
                  return lhs.score < rhs.score;
                }
                return lhs.candidate.crossing < rhs.candidate.crossing;
              });

    RankedCandidate& best = plausible_candidates.front();
    best.candidate.is_selected = true;
    best.candidate.candidate_qa_bits |= IsSelected;
    best.candidate.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::Selected);
    const double second_abs_dz = plausible_candidates.size() > 1U &&
                                        plausible_candidates[1].candidate.vertex_compatible == best.candidate.vertex_compatible &&
                                        plausible_candidates[1].candidate.confidence_tier == best.candidate.confidence_tier
                                    ? plausible_candidates[1].abs_delta_z
                                    : std::numeric_limits<double>::infinity();
    for (const RankedCandidate& candidate : plausible_candidates)
    {
      output->add_candidate(candidate.candidate);
    }

    output->set_selected_crossing(best.candidate.crossing);
    output->set_tpc_z0(static_cast<float>(tpc_z0));
    output->set_silicon_vertex_z(static_cast<float>(silicon_z0));
    output->set_delta_z(static_cast<float>(best.delta_z));
    output->set_best_abs_delta_z(static_cast<float>(best.abs_delta_z));
    output->set_second_best_abs_delta_z(static_cast<float>(second_abs_dz));
    output->set_selected_score(static_cast<float>(best.score));
    // Tpc_PolyClusterizer accepts tiers 0 and 1 by default. Use tier 2 for an
    // incompatible result so the corrector does not silently materialize it.
    output->set_selected_tier(best.candidate.vertex_compatible
                                  ? std::min<unsigned char>(best.candidate.confidence_tier, 1U)
                                  : 2U);
    output->set_number_of_tested_crossings(input->get_number_of_tested_crossings());
    output->set_number_of_available_crossings(input->get_number_of_available_crossings());
    output->set_number_of_allowed_crossings(input->get_number_of_allowed_crossings());
    output->set_number_of_tpc_valid_crossings(
        static_cast<unsigned short>(std::min<std::size_t>(
            plausible_candidates.size(), std::numeric_limits<unsigned short>::max())));
    output->set_number_of_vertex_compatible_crossings(ncompatible);
    output->set_status(best.candidate.vertex_compatible
                           ? TpcCrossingStatus::SelectedBySiliconTrack
                           : TpcCrossingStatus::SiliconTrackIncompatible);
    m_outputCrossings->add_decision(output);
    ++refined;
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - input=" << m_fullTracks->size()
              << " refined=" << refined << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
