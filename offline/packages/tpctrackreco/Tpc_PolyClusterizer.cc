#include "Tpc_PolyClusterizer.h"

#include "Tpc_AssembledTrack.h"
#include "Tpc_AssembledTrackContainer.h"
#include "IdealPadMap.h"
#include "Tpc_PolyClusterContainerv1.h"
#include "Tpc_PolyClusterv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

//#include <phgarfield/PHGarfield.h>
#include </sphenix/user/mitrankova/F4A/PHGarfield/install/include/phgarfield/PHGarfield.h>
#include <TPolyLine3D.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <thread>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    while (phi > M_PI) phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    return phi;
  }
}

Tpc_PolyClusterizer::Tpc_PolyClusterizer(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("TPC_ASSEMBLEDTRACKS")
  , m_outputNodeName("TPC_POLYCLUSTERS")
{
}

Tpc_PolyClusterizer::~Tpc_PolyClusterizer()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
  delete m_garfield;
  m_garfield = nullptr;
  for (PHGarfield* garfield : m_workerGarfields) delete garfield;
  m_workerGarfields.clear();
}

int Tpc_PolyClusterizer::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - failed to load IdealPadMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  delete m_garfield;
  //m_garfield = new PHGarfield(Name() + "_PHGarfield");

  const std::string electricFieldMap = "/sphenix/user/mitrankov/garf/include/sphenix_rossegger_garfield_field.root";

  m_garfield = new PHGarfield(Name() + "_PHGarfield", electricFieldMap, m_kEffSide0, m_kEffSide1);
  configure_garfield(m_garfield);
  if (m_garfield->InitRun(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << Name() << "::InitRun - PHGarfield InitRun failed" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (PHGarfield* garfield : m_workerGarfields) delete garfield;
  m_workerGarfields.clear();
  m_workerGarfields.reserve(24);
  for (unsigned int i = 0; i < 24; ++i)
  {
    PHGarfield* worker = new PHGarfield(Name() + "_PHGarfieldWorker" + std::to_string(i), electricFieldMap, m_kEffSide0, m_kEffSide1);
    configure_garfield(worker);
    if (worker->InitRun(topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      std::cerr << Name() << "::InitRun - PHGarfield worker InitRun failed" << std::endl;
      delete worker;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    m_workerGarfields.push_back(worker);
  }

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

void Tpc_PolyClusterizer::configure_garfield(PHGarfield* garfield) const
{
  if (!garfield) return;

  garfield->MoveTpc(m_tpcMove[0], m_tpcMove[1], m_tpcMove[2]);
  for (const auto& rotation : m_tpcRotations)
  {
    garfield->RotateTpc(rotation[0], rotation[1], rotation[2]);
  }
  garfield->SetCMVoltageDefault(m_cmVoltageDefault);
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
                                        PHGarfield* garfield,
                                        Point& p) const
{
  if (!m_hits || !m_idealPadMap || !garfield) return false;

  TrkrHitSet* hitset = m_hits->findHitSet(hsk);
  if (!hitset) return false;
  TrkrHit* hit = hitset->getHit(hk);
  if (!hit) return false;

  const unsigned int layer = TrkrDefs::getLayer(hsk);
  const unsigned int hit_side = TpcDefs::getSide(hsk);
  const unsigned int pad = TpcDefs::getPad(hk);
  const unsigned int tbin = TpcDefs::getTBin(hk);
  const double adc = hit->getAdc();
  if (layer < 7 || layer > 54) return false;
  if (hit_side >= 2U) return false;

  const double radius = m_idealPadMap->get_radius(layer);
  const double phi = m_idealPadMap->get_phi(hit_side, layer, pad);
  if (!std::isfinite(radius) || !std::isfinite(phi)) return false;

  const double corrected_tbin = static_cast<double>(tbin) - m_t0;
  const double target_time_ns = corrected_tbin * m_tpcAdcClock;
  if (target_time_ns <= 0.0 || !std::isfinite(target_time_ns)) return false;
  if (m_reverseDriftStepNs <= 0.0 || !std::isfinite(m_reverseDriftStepNs)) return false;

  const double x0 = radius * std::cos(phi);
  const double y0 = radius * std::sin(phi);
  const double z0 = (hit_side == 0U) ? m_startZSouth : m_startZNorth;

  TPolyLine3D* drift = garfield->ReverseDrift(x0, y0, z0, m_reverseDriftStepNs);
  if (!drift || drift->GetN() <= 0)
  {
    delete drift;
    return false;
  }

  const int npoints = drift->GetN();
  const Float_t* xyz = drift->GetP();
  if (!xyz || npoints <= 0)
  {
    delete drift;
    return false;
  }

  const double max_time_ns = static_cast<double>(npoints - 1) * m_reverseDriftStepNs;
  if (target_time_ns > max_time_ns)
  {
    delete drift;
    return false;
  }

  const double fbin = target_time_ns / m_reverseDriftStepNs;
  const int i0 = std::min(static_cast<int>(std::floor(fbin)), npoints - 1);
  const int i1 = std::min(i0 + 1, npoints - 1);
  const double frac = fbin - static_cast<double>(i0);

  const int idx0 = 3 * i0;
  const int idx1 = 3 * i1;
  const double x = xyz[idx0] + frac * (xyz[idx1] - xyz[idx0]);
  const double y = xyz[idx0 + 1] + frac * (xyz[idx1 + 1] - xyz[idx0 + 1]);
  const double z = xyz[idx0 + 2] + frac * (xyz[idx1 + 2] - xyz[idx0 + 2]);
  delete drift;

  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return false;
  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.layer = layer;
  p.side = hit_side;
  p.pad = pad;
  p.tbin = tbin;
  p.adc = adc;
  p.x = x;
  p.y = y;
  p.z = z;
  return true;
}

Tpc_PolyClusterizer::ClusterParameters
Tpc_PolyClusterizer::make_cluster_parameters(const std::vector<Point>& points,
                                            const Centroid& centroid,
                                            const int side) const
{
  ClusterParameters params;
  if (points.empty() || !centroid.ok || !m_idealPadMap) return params;

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
  if (points.empty()) return c;

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

void Tpc_PolyClusterizer::cluster_sector_side(unsigned int sector,
                                              int side,
                                              PHGarfield* garfield,
                                              std::vector<ClusterizedTrack>& output,
                                              unsigned int& nclusters) const
{
  output.clear();
  nclusters = 0;
  if (!m_assembledTracks || !garfield) return;

  const unsigned int nassembled = m_assembledTracks->size();
  for (unsigned int iassembled = 0; iassembled < nassembled; ++iassembled)
  {
    const Tpc_AssembledTrack* assembled = m_assembledTracks->get_track(iassembled);
    if (!assembled) continue;
    if (assembled->get_side() != side) continue;
    if (assembled->get_first_sector() % 12U != sector) continue;

    std::map<unsigned int, std::vector<Point>> points_by_layer;
    for (unsigned int ih = 0; ih < assembled->size_hit_indices(); ++ih)
    {
      const Tpc_AssembledTrack::HitIndex hi = assembled->get_hit_index(ih);
      if (TpcDefs::getSide(hi.first) != static_cast<unsigned int>(side)) continue;

      Point p;
      if (make_xyz_point(hi.first, hi.second, garfield, p)) points_by_layer[p.layer].push_back(p);
    }
    if (points_by_layer.empty()) continue;

    ClusterizedTrack out;
    out.source_assembled_track_id = assembled->get_track_id();
    out.side = side;

    for (const auto& layer_points : points_by_layer)
    {
      const unsigned int layer = layer_points.first;
      const std::vector<Point>& points = layer_points.second;
      const Centroid centroid = make_centroid(points);
      if (!centroid.ok) continue;

      const ClusterParameters params = make_cluster_parameters(points, centroid, static_cast<int>(points.front().side));
      out.layers.push_back(layer);
      out.centroids.push_back(centroid);
      out.parameters.push_back(params);
      out.cluster_points.push_back(points);
      ++nclusters;
    }

    if (!out.layers.empty()) output.push_back(out);
  }
}

int Tpc_PolyClusterizer::process_event(PHCompositeNode*)
{
  if (!m_assembledTracks || !m_clusters) return Fun4AllReturnCodes::EVENT_OK;
  m_clusters->Reset();

  const unsigned int nassembled = m_assembledTracks->size();
  std::vector<std::vector<ClusterizedTrack>> sector_outputs(24);
  std::vector<unsigned int> sector_nclusters(24, 0);
  std::vector<std::thread> workers;
  workers.reserve(24);

  for (int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      const unsigned int index = static_cast<unsigned int>(side) * 12U + sector;
      PHGarfield* garfield = index < m_workerGarfields.size() ? m_workerGarfields[index] : m_garfield;
      workers.push_back(std::thread(&Tpc_PolyClusterizer::cluster_sector_side, this,
                                    sector, side, garfield,
                                    std::ref(sector_outputs[index]),
                                    std::ref(sector_nclusters[index])));
    }
  }
  for (std::thread& worker : workers) worker.join();

  unsigned int nclusters = 0;
  for (unsigned int index = 0; index < sector_outputs.size(); ++index)
  {
    nclusters += sector_nclusters[index];
    for (const ClusterizedTrack& clustered : sector_outputs[index])
    {
      for (unsigned int ic = 0; ic < clustered.layers.size(); ++ic)
      {
        const Centroid& centroid = clustered.centroids[ic];
        Tpc_PolyClusterv1* out = new Tpc_PolyClusterv1();
        out->set_event(m_event);
        out->set_cluster_id(m_clusters->size());
        out->set_source_assembled_track_id(clustered.source_assembled_track_id);
        out->set_side(clustered.side);
        out->set_centroid_x(centroid.x);
        out->set_centroid_y(centroid.y);
        out->set_centroid_z(centroid.z);
        out->set_rms_x(centroid.rms_x);
        out->set_rms_y(centroid.rms_y);
        out->set_rms_z(centroid.rms_z);
        const ClusterParameters& params = clustered.parameters[ic];
        out->set_adc(params.adc);
        out->set_phi_width(params.phi_width);
        out->set_time_width(params.time_width);
        out->set_phase(params.phase);
        for (const Point& p : clustered.cluster_points[ic]) out->add_hit(p.hitsetkey, p.hitkey, p.x, p.y, p.z);
        if (out->size_hits() == 0)
        {
          delete out;
          continue;
        }
        m_clusters->add_cluster(out);
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " assembled_tracks=" << nassembled
              << " poly_clusters=" << m_clusters->size()
              << " layer_clusters=" << nclusters << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
