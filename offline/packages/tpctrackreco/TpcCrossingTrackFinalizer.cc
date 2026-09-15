#include "TpcCrossingTrackFinalizer.h"
#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrackv1.h"
#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
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
  m_lookup = TpcDriftPolylineLookup::get(topNode);
  m_field = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
  if (!m_tracks || !m_trajectories || !m_candidates || !m_clusters || !m_hits || !m_lookup || !m_field)
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
      double sw = 0., sx = 0., sy = 0., sz = 0., sx2 = 0., sy2 = 0., sz2 = 0.;
      for (const auto& hitIndex : reference->get_hit_indices())
      {
        auto* hitset = m_hits->findHitSet(hitIndex.first); auto* hit = hitset ? hitset->getHit(hitIndex.second) : nullptr;
        TpcDriftPolylineLookup::Point point;
        if (!hit || !m_lookup->getPosition(hitIndex.first, hitIndex.second, candidate->get_crossing(), point)) continue;
        const double weight = hit->getAdc(); sw += weight; sx += weight * point.x; sy += weight * point.y; sz += weight * point.z;
        sx2 += weight * point.x * point.x; sy2 += weight * point.y * point.y; sz2 += weight * point.z * point.z;
        corrected->add_hit(hitIndex.first, hitIndex.second, point.x, point.y, point.z);
      }
      if (sw <= 0. || corrected->size_hits() != reference->size_hits()) { delete corrected; continue; }
      const double x = sx / sw, y = sy / sw, z = sz / sw;
      corrected->set_centroid_x(x); corrected->set_centroid_y(y); corrected->set_centroid_z(z);
      corrected->set_rms_x(std::sqrt(std::max(0., sx2 / sw - x * x))); corrected->set_rms_y(std::sqrt(std::max(0., sy2 / sw - y * y))); corrected->set_rms_z(std::sqrt(std::max(0., sz2 / sw - z * z)));
      m_correctedClusters->add_cluster(corrected); fitClusters.push_back(corrected);
    }
    FastFieldTrackFitter::Result finalFit;
    const bool fitOk = fitClusters.size() == parent->size_cluster_keys() && m_fitter->fit(*parent, fitClusters, finalFit);

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
    const auto& state = fitOk ? finalFit.state : std::array<double, 6>{{trajectory->get_state(0), trajectory->get_state(1), trajectory->get_state(2), trajectory->get_state(3), trajectory->get_state(4), trajectory->get_state(5)}};
    full->set_x(state[0]); full->set_y(state[1]); full->set_z(state[2]);
    const double momentum = std::abs(state[5]) > 1.e-12 ? std::abs(1. / state[5]) : std::hypot(std::hypot(parent->get_px(), parent->get_py()), parent->get_pz());
    const double pt = momentum * std::sin(state[4]); full->set_px(pt * std::cos(state[3])); full->set_py(pt * std::sin(state[3])); full->set_pz(momentum * std::cos(state[4])); full->set_charge(state[5] < 0. ? -1. : 1.);
    for (unsigned int row = 0; row < 6; ++row) for (unsigned int col = 0; col < 6; ++col) full->set_cov(row, col, fitOk ? finalFit.covariance[row * 6 + col] : parent->get_cov(row, col));
    for (const auto key : parent->get_cluster_keys()) full->add_tpc_cluster_key(key);
    for (const auto key : candidate->get_silicon_cluster_keys()) full->add_silicon_cluster_key(key);
    m_output->add_track(full);
  }
  if (Verbosity() > 0) std::cout << Name() << " final_tracks=" << m_output->size() << " corrected_clusters=" << m_correctedClusters->size() << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
