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
    for (unsigned int ic = 0; ic < decision->get_number_of_candidates(); ++ic)
    {
      const auto* crossingCandidate = decision->get_candidate(ic);
      if (!crossingCandidate || !crossingCandidate->tpc_valid) continue;
      std::vector<std::array<double, 3>> original, displaced;
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
        original.push_back({cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z()});
        displaced.push_back({x / weightSum, y / weightSum, z / weightSum});
      }
      if (original.size() != clusters.size()) continue;
      const auto delta = m_fitter->linearUpdate(referenceFit, original, displaced);
      auto* trajectory = new TpcCrossingTrajectory;
      trajectory->set_parent_track_id(track->get_track_id());
      trajectory->set_source_assembled_track_id(track->get_source_assembled_track_id());
      trajectory->set_crossing(crossingCandidate->crossing);
      trajectory->set_reference_crossing(decision->get_reference_crossing());
      for (unsigned int k = 0; k < TpcCrossingTrajectory::StateSize; ++k) { trajectory->set_delta(k, delta[k]); trajectory->set_state(k, referenceFit.state[k] + delta[k]); }
      for (unsigned int row = 0; row < TpcCrossingTrajectory::StateSize; ++row) for (unsigned int col = 0; col < TpcCrossingTrajectory::StateSize; ++col) trajectory->set_covariance(row, col, referenceFit.covariance[row * TpcCrossingTrajectory::StateSize + col]);
      trajectory->set_linear_chi2(referenceFit.ndf > 0 ? referenceFit.chi2 / referenceFit.ndf : referenceFit.chi2);
      addSiliconStates(*trajectory, referenceFit);
      m_trajectories->add(trajectory); ++deltaBuilds;
    }
  }
  if (Verbosity() > 0) std::cout << Name() << " reference_field_fits=" << referenceFits << " crossing_updates=" << deltaBuilds << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
