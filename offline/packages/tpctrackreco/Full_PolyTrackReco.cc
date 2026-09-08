#include "Full_PolyTrackReco.h"
#include "Full_PolyTrack.h"
#include "Full_PolyTrackContainer.h"
#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrackv1.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

Full_PolyTrackReco::Full_PolyTrackReco(const std::string& name) : SubsysReco(name) {}

int Full_PolyTrackReco::InitRun(PHCompositeNode* topNode)
{
  return (getNodes(topNode) == Fun4AllReturnCodes::EVENT_OK &&
          createNodes(topNode) == Fun4AllReturnCodes::EVENT_OK)
             ? Fun4AllReturnCodes::EVENT_OK : Fun4AllReturnCodes::ABORTRUN;
}

int Full_PolyTrackReco::getNodes(PHCompositeNode* topNode)
{
  m_inputTracks = findNode::getClass<Full_PolyTrackContainer>(topNode, m_inputNodeName);
  m_clusters = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterNodeName);
  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_inputTracks || !m_clusters || !m_geometry)
  {
    std::cerr << Name() << "::getNodes - missing "
              << (!m_inputTracks ? m_inputNodeName : (!m_clusters ? m_clusterNodeName : "ActsGeometry"))
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackReco::createNodes(PHCompositeNode* topNode)
{
  if (m_outputNodeName == m_inputNodeName)
  {
    std::cerr << Name() << "::createNodes - input and output node names must differ" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  m_outputTracks = findNode::getClass<Full_PolyTrackContainer>(topNode, m_outputNodeName);
  if (!m_outputTracks)
  {
    m_outputTracks = new Full_PolyTrackContainerv1();
    dstNode->addNode(new PHIODataNode<PHObject>(m_outputTracks, m_outputNodeName, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

bool Full_PolyTrackReco::fitSiliconClusters(const Full_PolyTrack& track,
                                             Tpc_FittingTools::FitResult& fit) const
{
  std::vector<Tpc_FittingTools::Point> points;
  points.reserve(track.size_cluster_keys());
  for (const auto key : track.get_cluster_keys())
  {
    TrkrCluster* cluster = m_clusters->findCluster(key);
    if (!cluster) continue;
    const auto global = m_geometry->getGlobalPosition(key, cluster);
    Tpc_FittingTools::Point point{global.x(), global.y(), global.z()};
    if (std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z)) points.push_back(point);
  }
  // The circle fitter unwraps its arc from the first point, so order outward
  // from the beam line irrespective of the matcher's layer traversal order.
  std::sort(points.begin(), points.end(), [](const auto& lhs, const auto& rhs) {
    return std::hypot(lhs.x, lhs.y) < std::hypot(rhs.x, rhs.y);
  });
  return Tpc_FittingTools::fit(points, fit);
}

void Full_PolyTrackReco::fillTrack(const Full_PolyTrack& in,
                                   const Tpc_FittingTools::FitResult& fit, const bool ok)
{
  auto* out = new Full_PolyTrackv1();
  out->set_event(in.get_event()); out->set_track_id(in.get_track_id());
  out->set_tpc_poly_track_id(in.get_tpc_poly_track_id());
  out->set_source_assembled_track_id(in.get_source_assembled_track_id());
  out->set_n_tpc_clusters(in.get_n_tpc_clusters());
  out->set_n_missing_layers(in.get_n_missing_layers());
  out->set_missing_layer_mask(in.get_missing_layer_mask());
  out->set_score(in.get_score()); out->set_charge(in.get_charge());
  out->set_seed_x(in.get_seed_x()); out->set_seed_y(in.get_seed_y()); out->set_seed_z(in.get_seed_z());
  out->set_seed_px(in.get_seed_px()); out->set_seed_py(in.get_seed_py()); out->set_seed_pz(in.get_seed_pz());

  if (ok)
  {
    double phi = fit.phi0;
    const double seedPhi = std::atan2(in.get_seed_py(), in.get_seed_px());
    if (std::cos(phi - seedPhi) < 0.0) phi += M_PI;
    // pt is the curvature-bearing quantity and is intentionally not refitted.
    const double pt = std::hypot(in.get_seed_px(), in.get_seed_py());
    const double tanTheta = std::tan(fit.theta);
    out->set_x(-fit.d0 * std::sin(phi)); out->set_y(fit.d0 * std::cos(phi)); out->set_z(fit.z0);
    out->set_px(pt * std::cos(phi)); out->set_py(pt * std::sin(phi));
    out->set_pz(std::fabs(tanTheta) > 1.e-12 ? pt / tanTheta : in.get_seed_pz());
    out->set_chi2(fit.chi2_xy + fit.chi2_z); out->set_ndf(fit.ndof_xy + fit.ndof_z);
    out->set_fit_status(1);
  }
  else
  {
    out->set_x(in.get_x()); out->set_y(in.get_y()); out->set_z(in.get_z());
    out->set_px(in.get_px()); out->set_py(in.get_py()); out->set_pz(in.get_pz());
    out->set_chi2(in.get_chi2()); out->set_ndf(in.get_ndf()); out->set_fit_status(0);
  }
  for (const auto key : in.get_cluster_keys()) out->add_cluster_key(key);
  for (unsigned int i = 0; i < in.size_silicon_states(); ++i)
  {
    out->add_silicon_state(in.get_state_layer(i), in.get_state_cluster_key(i),
      in.get_state_x(i), in.get_state_y(i), in.get_state_z(i),
      in.get_state_pred_x(i), in.get_state_pred_y(i), in.get_state_pred_z(i),
      in.get_state_rdphi(i), in.get_state_dz(i), in.get_state_chi2(i));
  }
  m_outputTracks->add_track(out);
}

int Full_PolyTrackReco::process_event(PHCompositeNode* topNode)
{
  if ((!m_inputTracks || !m_outputTracks || !m_clusters || !m_geometry) &&
      (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK || createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK))
    return Fun4AllReturnCodes::ABORTEVENT;
  m_outputTracks->Reset();
  unsigned int fitted = 0;
  for (unsigned int i = 0; i < m_inputTracks->size(); ++i)
  {
    const auto* in = m_inputTracks->get_track(i);
    if (!in) continue;
    Tpc_FittingTools::FitResult fit;
    const bool ok = fitSiliconClusters(*in, fit);
    fitted += ok ? 1U : 0U;
    fillTrack(*in, fit, ok);
  }
  if (Verbosity() > 0) std::cout << Name() << "::process_event - input=" << m_inputTracks->size()
    << " fitted=" << fitted << " output=" << m_outputTracks->size() << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
