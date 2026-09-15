#include "TpcCrossingTrackFinalizer.h"
#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrackv1.h"
#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
#include "TpcCrossingClusterPosition.h"
#include "TpcSiliconMatchCandidate.h"
#include "TpcSiliconMatchCandidateContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyClusterv1.h"
#include "Tpc_PolyClusterContainer.h"
#include "Tpc_PolyClusterContainerv1.h"
#include "TpcDriftPolylineLookup.h"
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
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/ActsGeometry.h>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Definitions/Units.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

TpcCrossingTrackFinalizer::TpcCrossingTrackFinalizer(const std::string& name) : SubsysReco(name) {}
int TpcCrossingTrackFinalizer::getNodes(PHCompositeNode* topNode)
{
  m_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_trackNodeName);
  m_trajectories = findNode::getClass<TpcCrossingTrajectoryContainer>(topNode, m_trajectoryNodeName);
  m_candidates = findNode::getClass<TpcSiliconMatchCandidateContainer>(topNode, m_candidateNodeName);
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_clusterNodeName);
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_trkrClusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  m_lookup = TpcDriftPolylineLookup::get(topNode);
  m_field = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
  if (!m_tracks || !m_trajectories || !m_candidates || !m_clusters || !m_hits || !m_trkrClusters || !m_geometry || !m_lookup || !m_field)
  {
    std::cerr << Name() << "::getNodes - missing full-track input node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrackFinalizer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dst = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst) { dst = new PHCompositeNode("DST"); topNode->addNode(dst); }
  m_output = findNode::getClass<Full_PolyTrackContainer>(topNode, m_outputNodeName);
  if (!m_output)
  {
    m_output = new Full_PolyTrackContainerv1;
    dst->addNode(new PHIODataNode<PHObject>(m_output, m_outputNodeName, "PHObject"));
  }
  m_correctedClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_correctedClusterNodeName);
  if (!m_correctedClusters)
  {
    m_correctedClusters = new Tpc_PolyClusterContainerv1;
    dst->addNode(new PHIODataNode<PHObject>(m_correctedClusters, m_correctedClusterNodeName, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrackFinalizer::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  m_event = 0;
  m_fitter = std::make_unique<FastFieldTrackFitter>(m_field);
  return createNodes(topNode);
}
int TpcCrossingTrackFinalizer::process_event(PHCompositeNode*)
{
  const auto begin = std::chrono::steady_clock::now();
  ++m_event;
  m_output->Reset();
  m_correctedClusters->Reset();
  std::map<TrkrDefs::cluskey, const Tpc_PolyCluster*> clustersByKey;
  for (unsigned int i = 0; i < m_clusters->size(); ++i) if (const auto* cluster = m_clusters->get_cluster(i)) clustersByKey[cluster->get_trkr_cluster_key()] = cluster;
  unsigned int outputId = 0;
  unsigned int finalFits = 0;
  double finalFitSeconds = 0.0;
  for (unsigned int i = 0; i < m_candidates->size(); ++i)
  {
    const auto* candidate = m_candidates->get(i);
    if (!candidate || !candidate->get_selected()) continue;
    const Tpc_PolyTrack* parent = nullptr;
    for (unsigned int j = 0; j < m_tracks->size(); ++j)
    {
      const auto* value = m_tracks->get_track(j);
      if (value && value->get_track_id() == candidate->get_parent_track_id()) { parent = value; break; }
    }
    if (!parent) continue;
    const TpcCrossingTrajectory* trajectory = nullptr;
    for (unsigned int j = 0; j < m_trajectories->size(); ++j)
    {
      const auto* value = m_trajectories->get(j);
      if (value && value->get_parent_track_id() == candidate->get_parent_track_id() && value->get_crossing() == candidate->get_crossing()) { trajectory = value; break; }
    }
    if (!trajectory) continue;

    std::vector<const Tpc_PolyCluster*> fitClusters;
    for (const auto key : parent->get_cluster_keys())
    {
      const auto found = clustersByKey.find(key);
      if (found == clustersByKey.end()) continue;
      const auto* reference = found->second;
      auto* corrected = new Tpc_PolyClusterv1;
      corrected->set_event(m_event); corrected->set_cluster_id(m_correctedClusters->size());
      corrected->set_source_assembled_track_id(reference->get_source_assembled_track_id());
      corrected->set_trkr_cluster_key(reference->get_trkr_cluster_key()); corrected->set_side(reference->get_side());
      corrected->set_adc(reference->get_adc()); corrected->set_phi_width(reference->get_phi_width()); corrected->set_time_width(reference->get_time_width()); corrected->set_phase(reference->get_phase());
      std::array<double, 3> position{};
      std::vector<TpcCrossingClusterPosition::HitPosition> hitPositions;
      if (!TpcCrossingClusterPosition::get(*reference, *m_lookup, candidate->get_crossing(),
                                           trajectory->get_reference_crossing(), position, &hitPositions) ||
          hitPositions.size() != reference->size_hits())
      {
        delete corrected;
        continue;
      }
      double sx2 = 0.0, sy2 = 0.0, sz2 = 0.0;
      for (const auto& hitPosition : hitPositions)
      {
        corrected->add_hit(hitPosition.hitsetkey, hitPosition.hitkey,
                           hitPosition.candidate[0], hitPosition.candidate[1], hitPosition.candidate[2]);
        sx2 += (hitPosition.candidate[0] - position[0]) * (hitPosition.candidate[0] - position[0]);
        sy2 += (hitPosition.candidate[1] - position[1]) * (hitPosition.candidate[1] - position[1]);
        sz2 += (hitPosition.candidate[2] - position[2]) * (hitPosition.candidate[2] - position[2]);
      }
      const double inverseCount = 1.0 / static_cast<double>(hitPositions.size());
      corrected->set_centroid_x(position[0]); corrected->set_centroid_y(position[1]); corrected->set_centroid_z(position[2]);
      corrected->set_rms_x(std::sqrt(sx2 * inverseCount)); corrected->set_rms_y(std::sqrt(sy2 * inverseCount)); corrected->set_rms_z(std::sqrt(sz2 * inverseCount));
      m_correctedClusters->add_cluster(corrected); fitClusters.push_back(corrected);
    }
    std::vector<TpcTrackPoint> measurements;
    measurements.reserve(fitClusters.size() + candidate->get_silicon_cluster_keys().size());
    for (const auto* cluster : fitClusters)
    {
      TpcTrackPoint point;
      point.track_id = static_cast<int>(parent->get_track_id());
      point.layer = cluster->size_hits() ? static_cast<int>(TrkrDefs::getLayer(cluster->get_hit_index(0).first)) : 0;
      point.position = {cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z()};
      point.momentum = {parent->get_px(), parent->get_py(), parent->get_pz()};
      point.detector = TpcTrackPoint::Detector::Tpc;
      point.cluster_key = cluster->get_trkr_cluster_key();
      measurements.push_back(point);
    }
    for (const auto key : candidate->get_silicon_cluster_keys())
    {
      auto* cluster = m_trkrClusters->findCluster(key);
      const auto surface = cluster ? m_geometry->maps().getSurface(key, cluster) : nullptr;
      if (!cluster || !surface) continue;
      const auto global = m_geometry->getGlobalPosition(key, cluster);
      const Acts::Vector2 local(cluster->getLocalX() * Acts::UnitConstants::cm,
                                cluster->getLocalY() * Acts::UnitConstants::cm);
      const auto& context = m_geometry->geometry().geoContext;
      const Acts::Vector3 direction(1., 1., 1.);
      const auto origin = surface->localToGlobal(context, local, direction);
      const auto along0 = surface->localToGlobal(context, local + Acts::Vector2(Acts::UnitConstants::cm, 0.), direction) - origin;
      const auto along1 = surface->localToGlobal(context, local + Acts::Vector2(0., Acts::UnitConstants::cm), direction) - origin;
      const auto axis0 = along0.normalized();
      const auto axis1 = along1.normalized();
      const double sigma0 = std::max(1.e-6, static_cast<double>(cluster->getRPhiError()));
      const double sigma1 = std::max(1.e-6, static_cast<double>(cluster->getZError()));
      TpcTrackPoint point;
      point.track_id = static_cast<int>(parent->get_track_id());
      point.layer = static_cast<int>(TrkrDefs::getLayer(key));
      point.position = {global.x(), global.y(), global.z()};
      point.momentum = {parent->get_px(), parent->get_py(), parent->get_pz()};
      point.detector = TrkrDefs::getTrkrId(key) == TrkrDefs::mvtxId ? TpcTrackPoint::Detector::Mvtx : TpcTrackPoint::Detector::Intt;
      point.cluster_key = key;
      point.measurement_dimension = 2;
      point.has_measurement_model = true;
      point.measurement_projection = {axis0.x(), axis0.y(), axis0.z(), axis1.x(), axis1.y(), axis1.z(), 0., 0., 0.};
      point.measurement_covariance = {sigma0 * sigma0, 0., 0., 0., sigma1 * sigma1, 0., 0., 0., 1.};
      measurements.push_back(point);
    }
    FastFieldTrackFitter::Result finalFit;
    const bool allMeasurements = fitClusters.size() == parent->size_cluster_keys() &&
                                 measurements.size() == fitClusters.size() + candidate->get_silicon_cluster_keys().size();
    const bool fitOk = allMeasurements && m_fitter->fitMeasurements(*parent, measurements, finalFit);
    if (fitOk) { ++finalFits; finalFitSeconds += finalFit.fitSeconds; }

    auto* full = new Full_PolyTrackv1;
    full->set_event(m_event);
    full->set_track_id(outputId++);
    full->set_parent_tpc_track_id(parent->get_track_id());
    full->set_source_assembled_track_id(parent->get_source_assembled_track_id());
    full->set_crossing(candidate->get_crossing());
    full->set_status(1);
    full->set_n_mvtx(candidate->get_n_mvtx());
    full->set_n_intt(candidate->get_n_intt());
    full->set_score(candidate->get_score());
    full->set_max_abs_dz(candidate->get_max_abs_dz());
    full->set_max_abs_ddphi(candidate->get_max_abs_ddphi());
    full->set_fit_status(fitOk ? 1 : 0);
    full->set_chi2(fitOk ? finalFit.chi2 : parent->get_chi2()); full->set_ndf(fitOk ? finalFit.ndf : parent->get_ndf());
    const auto fallbackNative = std::array<double, 6>{{trajectory->get_state(0), trajectory->get_state(1), trajectory->get_state(2), trajectory->get_state(3), trajectory->get_state(4), trajectory->get_state(5)}};
    const auto fallbackState = FastFieldTrackFitter::externalState(fallbackNative);
    const auto& state = fitOk ? finalFit.state : fallbackState;
    full->set_x(state[0]); full->set_y(state[1]); full->set_z(state[2]);
    const double momentum = std::abs(state[5]) > 1.e-12 ? std::abs(1. / state[5]) : std::hypot(std::hypot(parent->get_px(), parent->get_py()), parent->get_pz());
    const double pt = momentum * std::sin(state[4]); full->set_px(pt * std::cos(state[3])); full->set_py(pt * std::sin(state[3])); full->set_pz(momentum * std::cos(state[4])); full->set_charge(state[5] < 0. ? -1. : 1.);
    for (unsigned int row = 0; row < 6; ++row) for (unsigned int col = 0; col < 6; ++col) full->set_cov(row, col, fitOk ? finalFit.covariance[row * 6 + col] : parent->get_cov(row, col));
    for (const auto key : parent->get_cluster_keys()) full->add_tpc_cluster_key(key);
    for (const auto key : candidate->get_silicon_cluster_keys()) full->add_silicon_cluster_key(key);
    m_output->add_track(full);
  }
  if (Verbosity() > 0) std::cout << Name() << " final_tracks=" << m_output->size() << " corrected_clusters=" << m_correctedClusters->size()
                                 << " combined_fits=" << finalFits << " final_fit_ms=" << (finalFits ? 1.e3 * finalFitSeconds / finalFits : 0.)
                                 << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
