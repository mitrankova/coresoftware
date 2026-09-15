#include "TpcCrossingTrajectoryBuilder.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
#include "TpcDriftPolylineLookup.h"
#include "TpcTrackKalmanFitter.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyClusterContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phfield/PHFieldUtility.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

TpcCrossingTrajectoryBuilder::TpcCrossingTrajectoryBuilder(const std::string& name) : SubsysReco(name) {}
int TpcCrossingTrajectoryBuilder::getNodes(PHCompositeNode* topNode)
{
  m_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_trackNodeName);
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_clusterNodeName);
  m_decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_lookup = TpcDriftPolylineLookup::get(topNode);
  m_field = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
  if (!m_tracks || !m_clusters || !m_decisions || !m_hits || !m_lookup || !m_field)
  {
    std::cerr << Name() << "::getNodes - missing track, cluster, decision, hit, drift lookup, or field input" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrajectoryBuilder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dst = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst) { dst = new PHCompositeNode("DST"); topNode->addNode(dst); }
  m_trajectories = findNode::getClass<TpcCrossingTrajectoryContainer>(topNode, m_outputNodeName);
  if (!m_trajectories) { m_trajectories = new TpcCrossingTrajectoryContainer; dst->addNode(new PHIODataNode<PHObject>(m_trajectories, m_outputNodeName, "PHObject")); }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrajectoryBuilder::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  m_fitter = std::make_unique<FastFieldTrackFitter>(m_field);
  return createNodes(topNode);
}
bool TpcCrossingTrajectoryBuilder::addSiliconStates(TpcCrossingTrajectory& trajectory, const FastFieldTrackFitter::Result& fit) const
{
  const double theta = trajectory.get_state(TpcCrossingTrajectory::Theta);
  const double sinTheta = std::sin(theta);
  if (std::abs(sinTheta) < 1.e-9) return false;
  std::array<double, TpcTrackKalmanFitter::StateDim> state{{
      trajectory.get_state(TpcCrossingTrajectory::X), trajectory.get_state(TpcCrossingTrajectory::Y),
      trajectory.get_state(TpcCrossingTrajectory::Z), trajectory.get_state(TpcCrossingTrajectory::Phi),
      trajectory.get_state(TpcCrossingTrajectory::QOverP) / sinTheta, 1.0 / std::tan(theta)}};
  for (int layer = static_cast<int>(m_siliconRadii.size()) - 1; layer >= 0; --layer)
  {
    const double target = m_siliconRadii[layer];
    bool valid = false;
    for (unsigned int step = 0; step < 1600; ++step)
    {
      const double radius = std::hypot(state[TpcTrackKalmanFitter::X], state[TpcTrackKalmanFitter::Y]);
      if (radius <= target + 0.08) { valid = true; break; }
      state = TpcTrackKalmanFitter::propagate_state(state, -0.25, fit.propagationConfig);
      if (!std::isfinite(state[TpcTrackKalmanFitter::X])) break;
    }
    TpcCrossingTrajectory::LayerState output;
    output.layer = static_cast<unsigned int>(layer);
    output.x = state[TpcTrackKalmanFitter::X]; output.y = state[TpcTrackKalmanFitter::Y]; output.z = state[TpcTrackKalmanFitter::Z];
    output.phi = std::atan2(output.y, output.x); output.valid = valid;
    trajectory.add_layer_state(output);
  }
  return true;
}
int TpcCrossingTrajectoryBuilder::process_event(PHCompositeNode*)
{
  const auto begin = std::chrono::steady_clock::now();
  m_trajectories->Reset();
  std::map<TrkrDefs::cluskey, const Tpc_PolyCluster*> byKey;
  for (unsigned int i = 0; i < m_clusters->size(); ++i) if (const auto* c = m_clusters->get_cluster(i)) byKey[c->get_trkr_cluster_key()] = c;
  unsigned int referenceFits = 0, deltaBuilds = 0;
  double referenceFitSeconds = 0.0, responseSeconds = 0.0, crossingSeconds = 0.0;
  for (unsigned int i = 0; i < m_tracks->size(); ++i)
  {
    const auto* track = m_tracks->get_track(i);
    if (!track || !track->get_fit_status()) continue;
    const auto* decision = m_decisions->get_decision(track->get_source_assembled_track_id());
    if (!decision) continue;
    std::vector<const Tpc_PolyCluster*> clusters;
    for (const auto key : track->get_cluster_keys()) { const auto found = byKey.find(key); if (found != byKey.end()) clusters.push_back(found->second); }
    FastFieldTrackFitter::Result referenceFit;
    if (!m_fitter->fit(*track, clusters, referenceFit)) continue;
    ++referenceFits;
    referenceFitSeconds += referenceFit.fitSeconds;
    responseSeconds += referenceFit.responseSeconds;
    for (unsigned int ic = 0; ic < decision->get_number_of_candidates(); ++ic)
    {
      const auto* crossingCandidate = decision->get_candidate(ic);
      if (!crossingCandidate || !crossingCandidate->tpc_valid) continue;
      const auto crossingBegin = std::chrono::steady_clock::now();
      std::map<TrkrDefs::cluskey, std::array<double, 3>> displaced;
      for (const auto* cluster : clusters)
      {
        double weightSum = 0.0, x = 0.0, y = 0.0, z = 0.0;
        for (const auto& hitIndex : cluster->get_hit_indices())
        {
          auto* hitset = m_hits->findHitSet(hitIndex.first); auto* hit = hitset ? hitset->getHit(hitIndex.second) : nullptr;
          TpcDriftPolylineLookup::Point point;
          if (!hit || !m_lookup->getPosition(hitIndex.first, hitIndex.second, crossingCandidate->crossing, point)) continue;
          const double weight = hit->getAdc(); weightSum += weight; x += weight * point.x; y += weight * point.y; z += weight * point.z;
        }
        if (weightSum <= 0.0) continue;
        displaced[cluster->get_trkr_cluster_key()] = {x / weightSum, y / weightSum, z / weightSum};
      }
      if (displaced.size() != clusters.size()) continue;
      const auto update = m_fitter->linearUpdate(referenceFit, displaced);
      if (!update.valid) continue;
      auto* trajectory = new TpcCrossingTrajectory;
      trajectory->set_parent_track_id(track->get_track_id());
      trajectory->set_source_assembled_track_id(track->get_source_assembled_track_id());
      trajectory->set_crossing(crossingCandidate->crossing);
      trajectory->set_reference_crossing(decision->get_reference_crossing());
      for (unsigned int k = 0; k < TpcCrossingTrajectory::StateSize; ++k) { trajectory->set_delta(k, update.delta[k]); trajectory->set_state(k, referenceFit.state[k] + update.delta[k]); }
      for (unsigned int row = 0; row < TpcCrossingTrajectory::StateSize; ++row) for (unsigned int col = 0; col < TpcCrossingTrajectory::StateSize; ++col) trajectory->set_covariance(row, col, referenceFit.covariance[row * TpcCrossingTrajectory::StateSize + col]);
      trajectory->set_linear_chi2(update.chi2);
      addSiliconStates(*trajectory, referenceFit);
      const unsigned int validationHash = 2654435761U * (track->get_track_id() + 1U) + static_cast<unsigned int>(crossingCandidate->crossing);
      const bool validate = m_validationFraction > 0.0 &&
                            static_cast<double>(validationHash % 1000000U) / 1000000.0 < std::min(1.0, m_validationFraction);
      if (validate)
      {
        std::vector<TpcTrackPoint> candidatePoints;
        candidatePoints.reserve(clusters.size());
        for (const auto* cluster : clusters)
        {
          const auto& position = displaced.at(cluster->get_trkr_cluster_key());
          TpcTrackPoint point;
          point.track_id = static_cast<int>(track->get_track_id());
          point.layer = cluster->size_hits() ? static_cast<int>(TrkrDefs::getLayer(cluster->get_hit_index(0).first)) : 0;
          point.position = {position[0], position[1], position[2]};
          point.momentum = {track->get_px(), track->get_py(), track->get_pz()};
          point.detector = TpcTrackPoint::Detector::Tpc;
          point.cluster_key = cluster->get_trkr_cluster_key();
          candidatePoints.push_back(point);
        }
        FastFieldTrackFitter::Result candidateFit;
        if (m_fitter->fitMeasurements(*track, candidatePoints, candidateFit))
        {
          std::cout << Name() << " validation crossing=" << crossingCandidate->crossing;
          for (unsigned int k = 0; k < TpcCrossingTrajectory::StateSize; ++k)
            std::cout << " delta" << k << "_linear=" << update.delta[k]
                      << " delta" << k << "_full=" << candidateFit.state[k] - referenceFit.state[k];
          std::cout << std::endl;
          auto linearNative = referenceFit.nativeState;
          const double linearTheta = trajectory->get_state(TpcCrossingTrajectory::Theta);
          linearNative[TpcTrackKalmanFitter::X] = trajectory->get_state(TpcCrossingTrajectory::X);
          linearNative[TpcTrackKalmanFitter::Y] = trajectory->get_state(TpcCrossingTrajectory::Y);
          linearNative[TpcTrackKalmanFitter::Z] = trajectory->get_state(TpcCrossingTrajectory::Z);
          linearNative[TpcTrackKalmanFitter::Phi] = trajectory->get_state(TpcCrossingTrajectory::Phi);
          linearNative[TpcTrackKalmanFitter::QOverPt] = trajectory->get_state(TpcCrossingTrajectory::QOverP) / std::sin(linearTheta);
          linearNative[TpcTrackKalmanFitter::TanLambda] = 1.0 / std::tan(linearTheta);
          auto linearLayer = linearNative;
          auto fullLayer = candidateFit.nativeState;
          for (int layer = static_cast<int>(m_siliconRadii.size()) - 1; layer >= 0; --layer)
          {
            const double target = m_siliconRadii[layer];
            for (unsigned int step = 0; step < 1600 && std::hypot(linearLayer[0], linearLayer[1]) > target + 0.08; ++step)
              linearLayer = TpcTrackKalmanFitter::propagate_state(linearLayer, -0.25, referenceFit.propagationConfig);
            for (unsigned int step = 0; step < 1600 && std::hypot(fullLayer[0], fullLayer[1]) > target + 0.08; ++step)
              fullLayer = TpcTrackKalmanFitter::propagate_state(fullLayer, -0.25, candidateFit.propagationConfig);
            const double linearPhi = std::atan2(linearLayer[1], linearLayer[0]);
            const double fullPhi = std::atan2(fullLayer[1], fullLayer[0]);
            std::cout << Name() << " validation crossing=" << crossingCandidate->crossing << " si_layer=" << layer
                      << " delta_rphi_linear_minus_full=" << target * std::remainder(linearPhi - fullPhi, 2. * M_PI)
                      << " delta_z_linear_minus_full=" << linearLayer[2] - fullLayer[2] << std::endl;
          }
        }
      }
      m_trajectories->add(trajectory); ++deltaBuilds;
      crossingSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - crossingBegin).count();
      if (Verbosity() > 1)
      {
        std::cout << Name() << " crossing=" << crossingCandidate->crossing
                  << " delta_d0=" << update.delta[0] << " delta_z0=" << update.delta[2]
                  << " delta_phi=" << update.delta[3] << " delta_theta=" << update.delta[4]
                  << " delta_q_over_p=" << update.delta[5] << " linear_chi2=" << update.chi2 << std::endl;
      }
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << Name() << " reference_field_fits=" << referenceFits << " crossing_updates=" << deltaBuilds
              << " reference_fit_ms=" << (referenceFits ? 1.e3 * referenceFitSeconds / referenceFits : 0.)
              << " response_build_ms=" << (referenceFits ? 1.e3 * responseSeconds / referenceFits : 0.)
              << " crossing_update_us=" << (deltaBuilds ? 1.e6 * crossingSeconds / deltaBuilds : 0.)
              << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
