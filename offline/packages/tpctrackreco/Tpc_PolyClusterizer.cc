#include "Tpc_PolyClusterizer.h"

#include "IdealPadMap.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include "Tpc_AssembledTrack.h"
#include "Tpc_AssembledTrackContainer.h"
#include "Tpc_PolyClusterContainerv1.h"
#include "Tpc_PolyClusterv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    while (phi > M_PI)
    {
      phi -= 2.0 * M_PI;
    }
    while (phi <= -M_PI)
    {
      phi += 2.0 * M_PI;
    }
    return phi;
  }

  double square(const double value)
  {
    return value * value;
  }
}  // namespace

Tpc_PolyClusterizer::Tpc_PolyClusterizer(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("TPC_ASSEMBLEDTRACKS")
  , m_outputNodeName("TPC_POLYCLUSTERS")
{
}

int Tpc_PolyClusterizer::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_driftLookup = TpcDriftPolylineLookup::getOrCreate(topNode, Verbosity());
  if (!m_driftLookup)
  {
    std::cerr << Name() << "::InitRun - failed to obtain " << TpcDriftPolylineLookup::NodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (!m_driftLookup->registerConfiguration(m_driftConfig, Name(), Verbosity()))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (!m_driftLookup->initialize(topNode, Verbosity()))
  {
    std::cerr << Name() << "::InitRun - failed to initialize shared TPC drift lookup" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_idealPadMap = m_driftLookup->idealPadMap();
  if (!m_idealPadMap)
  {
    std::cerr << Name() << "::InitRun - shared TPC drift lookup has no IdealPadMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyClusterizer::getNodes(PHCompositeNode* topNode)
{
  m_assembledTracks = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_inputNodeName);
  if (!m_assembledTracks)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  if (!m_crossingDecisions)
  {
    std::cerr << Name() << "::getNodes - missing " << m_crossingDecisionNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomContainerTpc = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!m_geomContainerTpc)
  {
    std::cerr << Name() << "::getNodes - missing TPCGEOMCONTAINER" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyClusterizer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_outputNodeName);
  if (!m_clusters)
  {
    m_clusters = new Tpc_PolyClusterContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_clusters, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool Tpc_PolyClusterizer::make_xyz_point(TrkrDefs::hitsetkey hsk,
                                         TrkrDefs::hitkey hk,
                                         const short crossing,
                                         Point& p) const
{
  if (!m_hits || !m_driftLookup)
  {
    return false;
  }

  TrkrHitSet* hitset = m_hits->findHitSet(hsk);
  if (!hitset)
  {
    return false;
  }
  TrkrHit* hit = hitset->getHit(hk);
  if (!hit)
  {
    return false;
  }

  const unsigned int layer = TrkrDefs::getLayer(hsk);
  const unsigned int hit_side = TpcDefs::getSide(hsk);
  const unsigned int pad = TpcDefs::getPad(hk);
  const unsigned int tbin = TpcDefs::getTBin(hk);

  TpcDriftPolylineLookup::Point drift_point;
  if (!m_driftLookup->getPosition(layer, hit_side, pad, tbin, crossing, drift_point))
  {
    return false;
  }

  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.layer = layer;
  p.side = hit_side;
  p.pad = pad;
  p.tbin = tbin;
  p.adc = hit->getAdc();
  p.x = drift_point.x;
  p.y = drift_point.y;
  p.z = drift_point.z;
  return true;
}

Tpc_PolyClusterizer::ClusterParameters
Tpc_PolyClusterizer::make_cluster_parameters(const std::vector<Point>& points,
                                             const Centroid& centroid,
                                             const int side) const
{
  ClusterParameters params;
  if (points.empty() || !centroid.ok || !m_idealPadMap)
  {
    return params;
  }

  std::set<unsigned int> pads;
  std::set<unsigned int> tbins;
  std::map<unsigned int, double> adc_by_pad;
  for (const Point& p : points)
  {
    params.adc += p.adc;
    pads.insert(p.pad);
    tbins.insert(p.tbin);
    adc_by_pad[p.pad] += p.adc;
  }

  params.phi_width = static_cast<unsigned int>(pads.size());
  params.time_width = static_cast<unsigned int>(tbins.size());

  unsigned int max_adc_pad = 0;
  double max_adc = -std::numeric_limits<double>::max();
  for (const auto& pad_adc : adc_by_pad)
  {
    if (pad_adc.second > max_adc)
    {
      max_adc = pad_adc.second;
      max_adc_pad = pad_adc.first;
    }
  }

  const unsigned int total_phibins = m_idealPadMap->get_total_phibins(centroid.layer);
  const double pad_phi_width = total_phibins > 0U ? 2.0 * M_PI / static_cast<double>(total_phibins) : 0.0;
  const double cluster_phi = std::atan2(centroid.y, centroid.x);
  const double max_adc_phi = m_idealPadMap->get_phi(static_cast<unsigned int>(side), centroid.layer, max_adc_pad);
  if (pad_phi_width > 0.0 && std::isfinite(cluster_phi) && std::isfinite(max_adc_phi))
  {
    params.phase = wrap_phi(cluster_phi - max_adc_phi) / pad_phi_width;
  }

  return params;
}

Tpc_PolyClusterizer::Centroid
Tpc_PolyClusterizer::make_centroid(const std::vector<Point>& points)
{
  Centroid c;
  if (points.empty())
  {
    return c;
  }

  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
  for (const Point& p : points)
  {
    sx += p.x;
    sy += p.y;
    sz += p.z;
  }

  const double n = static_cast<double>(points.size());
  c.x = sx / n;
  c.y = sy / n;
  c.z = sz / n;

  double sxx = 0.0;
  double syy = 0.0;
  double szz = 0.0;
  for (const Point& p : points)
  {
    const double dx = p.x - c.x;
    const double dy = p.y - c.y;
    const double dz = p.z - c.z;
    sxx += dx * dx;
    syy += dy * dy;
    szz += dz * dz;
  }

  c.rms_x = std::sqrt(sxx / n);
  c.rms_y = std::sqrt(syy / n);
  c.rms_z = std::sqrt(szz / n);
  c.layer = points.front().layer;
  c.ok = std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z);
  return c;
}

int Tpc_PolyClusterizer::process_event(PHCompositeNode* topNode)
{
  if (!m_assembledTracks || !m_clusters || !m_crossingDecisions || !m_driftLookup || !m_driftLookup->isInitialized())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  ActsGeometry* tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cerr << Name() << "::process_event - missing ActsGeometry, using RMS fallback for errors" << std::endl;
  }
  m_clusters->Reset();

  const unsigned int nassembled = m_assembledTracks->size();
  unsigned int nclusters = 0;
  std::map<TrkrDefs::hitsetkey, unsigned int> next_cluster_index_by_hitset;
  unsigned int nmissing_decision = 0;
  unsigned int nskipped_tier = 0;
  unsigned int nempty_points = 0;

  for (int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      for (unsigned int iassembled = 0; iassembled < nassembled; ++iassembled)
      {
        const Tpc_AssembledTrack* assembled = m_assembledTracks->get_track(iassembled);
        if (!assembled)
        {
          continue;
        }
        if (assembled->get_side() != side)
        {
          continue;
        }
        if (assembled->get_first_sector() % 12U != sector)
        {
          continue;
        }
        const TpcCrossingDecision* crossing_decision = m_crossingDecisions->get_decision(assembled->get_track_id());
        if (!crossing_decision)
        {
          ++nmissing_decision;
          continue;
        }
        const unsigned char selected_tier = crossing_decision->get_selected_tier();
        if (selected_tier > m_maxAcceptedTier)
        {
          ++nskipped_tier;
          continue;
        }
        const short selected_crossing = crossing_decision->get_selected_crossing();

        std::map<TrkrDefs::hitsetkey, std::vector<Point>> points_by_hitset;
        for (unsigned int ih = 0; ih < assembled->size_hit_indices(); ++ih)
        {
          const Tpc_AssembledTrack::HitIndex hi = assembled->get_hit_index(ih);
          if (TpcDefs::getSide(hi.first) != static_cast<unsigned int>(side))
          {
            continue;
          }

          Point p;
          if (make_xyz_point(hi.first, hi.second, selected_crossing, p))
          {
            points_by_hitset[p.hitsetkey].push_back(p);
          }
        }
        if (points_by_hitset.empty())
        {
          ++nempty_points;
          continue;
        }

        for (const auto& hitset_points : points_by_hitset)
        {
          const TrkrDefs::hitsetkey cluster_hitsetkey = hitset_points.first;
          const std::vector<Point>& points = hitset_points.second;
          const Centroid centroid = make_centroid(points);
          if (!centroid.ok)
          {
            continue;
          }

          const unsigned int cluster_index = next_cluster_index_by_hitset[cluster_hitsetkey]++;
          const TrkrDefs::cluskey trkr_cluster_key = TrkrDefs::genClusKey(cluster_hitsetkey, cluster_index);

          Tpc_PolyClusterv1* out = new Tpc_PolyClusterv1();
          out->set_event(m_event);
          out->set_cluster_id(m_clusters->size());
          out->set_source_assembled_track_id(assembled->get_track_id());
          out->set_trkr_cluster_key(trkr_cluster_key);
          out->set_side(side);
          out->set_centroid_x(centroid.x);
          out->set_centroid_y(centroid.y);
          out->set_centroid_z(centroid.z);

          double phi_error = std::hypot(centroid.rms_x, centroid.rms_y);
          double z_error = std::fabs(centroid.rms_z);
          PHG4TpcGeom* layergeom = m_geomContainerTpc ? m_geomContainerTpc->GetLayerCellGeom(centroid.layer) : nullptr;
          if (layergeom && tGeometry)
          {
            double adc_sum = 0.0;
            double iphi_sum = 0.0;
            double iphi2_sum = 0.0;
            double t_sum = 0.0;
            double t2_sum = 0.0;
            int phibinhi = -1;
            int phibinlo = std::numeric_limits<int>::max();
            int tbinhi = -1;
            int tbinlo = std::numeric_limits<int>::max();

            for (const Point& p : points)
            {
              if (p.adc <= 0.0)
              {
                continue;
              }

              const int iphi = static_cast<int>(p.pad);
              const int it = static_cast<int>(p.tbin);
              if(it >= layergeom->get_zbins())
              {
                continue;
              }
              const double adc = p.adc;
              phibinhi = std::max(iphi, phibinhi);
              phibinlo = std::min(iphi, phibinlo);
              tbinhi = std::max(it, tbinhi);
              tbinlo = std::min(it, tbinlo);

              iphi_sum += static_cast<double>(iphi) * adc;
              iphi2_sum += square(static_cast<double>(iphi)) * adc;

              const double t = layergeom->get_zcenter(it);
              t_sum += t * adc;
              t2_sum += square(t) * adc;
              adc_sum += adc;
            }

            const double drift_velocity = tGeometry->get_drift_velocity();
            if (adc_sum > 0.0 && std::isfinite(drift_velocity))
            {
              const double radius = layergeom->get_radius();
              const double clusiphi = iphi_sum / adc_sum;
              const double clust = t_sum / adc_sum;
              const double phi_cov = std::max(0.0, (iphi2_sum / adc_sum - square(clusiphi)) * square(layergeom->get_phistep()));
              const double t_cov = std::max(0.0, t2_sum / adc_sum - square(clust));
              const double phi_err_square = (phibinhi == phibinlo) ? 9.0 * (square(radius * layergeom->get_phistep()) / 12.0) : square(radius) * phi_cov / (adc_sum * 0.14);
              const double t_err_square = (tbinhi == tbinlo) ? 9.0 * (square(layergeom->get_zstep()) / 12.0) : t_cov / (adc_sum * 0.14);

              if (phi_err_square >= 0.0 && std::isfinite(phi_err_square))
              {
                phi_error = std::sqrt(phi_err_square);
              }
              const double z_err_square = t_err_square * square(drift_velocity);
              if (z_err_square >= 0.0 && std::isfinite(z_err_square))
              {
                z_error = std::sqrt(z_err_square);
              }
            }
          }

          out->set_rms_x(phi_error);
          out->set_rms_y(0.0);
          out->set_rms_z(z_error);

          const ClusterParameters params = make_cluster_parameters(points, centroid, static_cast<int>(points.front().side));
          out->set_adc(params.adc);
          out->set_phi_width(params.phi_width);
          out->set_time_width(params.time_width);
          out->set_phase(params.phase);
          for (const Point& p : points)
          {
            out->add_hit(p.hitsetkey, p.hitkey, p.x, p.y, p.z);
          }
          if (out->size_hits() == 0)
          {
            delete out;
            continue;
          }
          m_clusters->add_cluster(out);
          ++nclusters;
        }
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " assembled_tracks=" << nassembled
              << " poly_clusters=" << m_clusters->size()
              << " layer_clusters=" << nclusters
              << " missing_decisions=" << nmissing_decision
              << " skipped_tier=" << nskipped_tier
              << " empty_points=" << nempty_points << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
