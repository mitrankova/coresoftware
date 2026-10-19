#include "Tpc_PolyTripletReco.h"

#include "IdealPadMap.h"
#include "Tpc_AssembledTrack.h"
#include "Tpc_AssembledTrackContainer.h"
#include "Tpc_CrossingTrajectoryContainerv1.h"
#include "Tpc_CrossingTrajectoryv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <phgarfield/PHGarfield.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TPolyLine3D.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>

namespace
{
  double unwrapNear(double phi, double ref)
  {
    while (phi - ref > M_PI)
    {
      phi -= 2.0 * M_PI;
    }
    while (phi - ref <= -M_PI)
    {
      phi += 2.0 * M_PI;
    }
    return phi;
  }
}  // namespace

Tpc_PolyTripletReco::Tpc_PolyTripletReco(const std::string& name)
  : SubsysReco(name)
{
}

Tpc_PolyTripletReco::~Tpc_PolyTripletReco()
{
  delete m_padMap;
}

int Tpc_PolyTripletReco::InitRun(PHCompositeNode* topNode)
{
  m_assembledTracks = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_inputNodeName);
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_geomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  m_garfield = findNode::getClass<PHGarfield>(topNode, "PHGarfield");
  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_assembledTracks || !m_hits || !m_geomContainer || !m_garfield || !m_geometry)
  {
    std::cerr << Name() << "::InitRun - missing input (" << m_inputNodeName
              << ", TRKR_HITSET, TPCGEOMCONTAINER, PHGarfield or ActsGeometry)" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  delete m_padMap;
  m_padMap = new IdealPadMap();
  if (m_padMap->load_from_cdb(Verbosity()) != 0 || !m_padMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - failed to load IdealPadMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // drift start z and timing, as in Tpc_PolyClusterizer
  PHG4TpcGeom* layergeom = m_geomContainer->GetLayerCellGeom(20);
  if (!layergeom)
  {
    std::cerr << Name() << "::InitRun - missing TPC geometry for layer 20" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  const double zStart = layergeom->get_max_driftlength() + layergeom->get_CM_halfwidth();
  m_startZ = {{-zStart, zStart}};
  if (m_garfield->GetIsNonDistortedMC())
  {
    if (!m_t0Override)
    {
      m_t0 = m_garfield->GetNonDistortedMCT0Tbins();
    }
    if (!m_tpcAdcClockOverride && layergeom->get_adc_clock() > 0.0)
    {
      m_tpcAdcClock = layergeom->get_adc_clock();
    }
  }
  m_t0Side = {{m_t0 + m_garfield->GetDeltaT0Tbins(0), m_t0 + m_garfield->GetDeltaT0Tbins(1)}};

  // first and last layer of every module (region)
  for (unsigned int i = 0; i < NLayers; ++i)
  {
    const unsigned int layer = FirstLayer + i;
    const int region = m_padMap->get_region(layer);
    m_moduleEdgeLayer[i] = (i == 0 || i + 1 == NLayers || m_padMap->get_region(layer - 1) != region ||
                            m_padMap->get_region(layer + 1) != region);
  }

  // drift polylines of every layer, built once here (not in the first event)
  m_lookup.assign(NLayers * 2 * 12 * m_nSamples, Polyline());
  m_layerState.fill(0);
  const auto tLookup = std::chrono::steady_clock::now();
  for (unsigned int i = 0; i < NLayers; ++i)
  {
    m_layerState[i] = buildLayer(FirstLayer + i) ? 1 : -1;
  }
  m_lookupMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tLookup).count();
  std::cout << Name() << "::InitRun - drift lookup: " << NLayers * 2 * 12 * m_nSamples << " polylines ("
            << m_nSamples << " per sector, substeps " << (m_driftSubsteps > 0 ? std::to_string(m_driftSubsteps) : std::string("default"))
            << ") built in " << m_lookupMs / 1000.0 << " s" << std::endl;

  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  m_trajectories = findNode::getClass<Tpc_CrossingTrajectoryContainer>(topNode, m_outputNodeName);
  if (!m_trajectories)
  {
    m_trajectories = new Tpc_CrossingTrajectoryContainerv1();
    dstNode->addNode(new PHIODataNode<PHObject>(m_trajectories, m_outputNodeName, "PHObject"));
  }

  std::cout << Name() << "::InitRun - " << (m_triggeredMode ? "triggered mode, crossing " + std::to_string(m_triggeredCrossing)
                                                             : "streaming mode, all INTT crossings with |z0| < " + std::to_string(m_maxVertexZ) + " cm")
            << "; edge pad margin " << m_edgePadMargin << ", min layer span " << m_minLayerSpan
            << "; t0 side0=" << m_t0Side[0] << " side1=" << m_t0Side[1] << " tbins; output "
            << m_outputNodeName << std::endl;
  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

// Reverse-drift polylines from m_nSamples pad positions per sector, one layer, both sides
bool Tpc_PolyTripletReco::buildLayer(unsigned int layer)
{
  const double radius = m_padMap->get_radius(layer);
  const unsigned int padsPerSector = m_padMap->get_pads_per_sector_for_layer(layer);
  if (!std::isfinite(radius) || padsPerSector == 0)
  {
    std::cerr << Name() << "::buildLayer - invalid geometry for layer " << layer << std::endl;
    return false;
  }
  for (unsigned int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      for (unsigned int s = 0; s < m_nSamples; ++s)
      {
        const auto localPad = static_cast<unsigned int>(std::llround(
            static_cast<double>(s) / (m_nSamples - 1) * static_cast<double>(padsPerSector - 1)));
        const double phi = m_padMap->get_phi(side, sector, layer, localPad);
        TPolyLine3D* line = m_driftSubsteps > 0
                                ? m_garfield->ReverseDrift(radius * std::cos(phi), radius * std::sin(phi), m_startZ[side],
                                                           m_driftStepNs, nullptr, m_driftSubsteps)
                                : m_garfield->ReverseDrift(radius * std::cos(phi), radius * std::sin(phi), m_startZ[side],
                                                           m_driftStepNs);
        if (!line || line->GetN() <= 0 || !line->GetP())
        {
          delete line;
          std::cerr << Name() << "::buildLayer - ReverseDrift failed, layer " << layer << " side " << side
                    << " sector " << sector << std::endl;
          return false;
        }
        Polyline& pl = polyline(layer, side, sector, s);
        pl.phi = phi;
        pl.points.resize(line->GetN());
        m_maxLookupTimeNs = std::max(m_maxLookupTimeNs, (line->GetN() - 1) * m_driftStepNs);
        const Float_t* xyz = line->GetP();
        for (int i = 0; i < line->GetN(); ++i)
        {
          const double x = xyz[3 * i], y = xyz[3 * i + 1];
          pl.points[i].dr = static_cast<float>(std::hypot(x, y) - radius);
          pl.points[i].dphi = static_cast<float>(unwrapNear(std::atan2(y, x), phi) - phi);
          pl.points[i].z = xyz[3 * i + 2];
        }
        delete line;
      }
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << Name() << "::buildLayer - built drift polylines for layer " << layer << std::endl;
  }
  return true;
}

// Position of a readout point (layer, continuous phi) after timeNs of reverse drift:
// linear in time along each polyline, linear in phi between the two neighbouring polylines.
bool Tpc_PolyTripletReco::drift(unsigned int layer, unsigned int side, unsigned int sector, double phi, double timeNs,
                                double& x, double& y, double& z) const
{
  if (timeNs <= 0.0 || !std::isfinite(timeNs) || !std::isfinite(phi))
  {
    return false;
  }

  // bracketing polylines in phi (clamped to the sector edges)
  const unsigned int ns = m_nSamples;
  std::array<double, MaxSamples> sphi{};
  for (unsigned int s = 0; s < ns; ++s)
  {
    sphi[s] = (s == 0) ? polyline(layer, side, sector, 0).phi : unwrapNear(polyline(layer, side, sector, s).phi, sphi[s - 1]);
  }
  phi = unwrapNear(phi, sphi[ns / 2]);
  const double dir = (sphi[ns - 1] >= sphi[0]) ? 1.0 : -1.0;
  unsigned int s0 = 0;
  while (s0 + 2 < ns && dir * (phi - sphi[s0 + 1]) > 0.0)
  {
    ++s0;
  }
  const unsigned int s1 = s0 + 1;
  const double denom = sphi[s1] - sphi[s0];
  const double f = denom != 0.0 ? std::clamp((phi - sphi[s0]) / denom, 0.0, 1.0) : 0.0;

  auto atTime = [&](const Polyline& pl, double& dr, double& dphi, double& pz)
  {
    if (pl.points.empty())
    {
      return false;
    }
    const double fbin = timeNs / m_driftStepNs;
    const auto i0 = static_cast<std::size_t>(fbin);
    if (i0 >= pl.points.size())
    {
      return false;
    }
    const std::size_t i1 = std::min(i0 + 1, pl.points.size() - 1);
    const double t = fbin - static_cast<double>(i0);
    const DriftPoint& a = pl.points[i0];
    const DriftPoint& b = pl.points[i1];
    dr = a.dr + t * (b.dr - a.dr);
    dphi = a.dphi + t * (b.dphi - a.dphi);
    pz = a.z + t * (b.z - a.z);
    return true;
  };

  double dr0 = 0, dphi0 = 0, z0 = 0, dr1 = 0, dphi1 = 0, z1 = 0;
  const bool ok0 = atTime(polyline(layer, side, sector, s0), dr0, dphi0, z0);
  const bool ok1 = atTime(polyline(layer, side, sector, s1), dr1, dphi1, z1);
  if (!ok0 && !ok1)
  {
    return false;
  }
  const double w = !ok1 ? 0.0 : (!ok0 ? 1.0 : f);
  const double r = m_padMap->get_radius(layer) + (1.0 - w) * dr0 + w * dr1;
  const double p = phi + (1.0 - w) * dphi0 + w * dphi1;
  z = (1.0 - w) * z0 + w * z1;
  x = r * std::cos(p);
  y = r * std::sin(p);
  return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

// crossing-independent part of a reference point: ADC-weighted pad phi and time bin on the layer
Tpc_PolyTripletReco::LayerMean Tpc_PolyTripletReco::layerMean(unsigned int layer, unsigned int side, const std::vector<Hit>& hits) const
{
  LayerMean m;
  m.layer = layer;
  if (hits.empty() || m_layerState[layer - FirstLayer] != 1)
  {
    return m;
  }
  const Hit& peak = *std::max_element(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.adc < b.adc; });
  const double phiPeak = m_padMap->get_phi(side, layer, peak.pad);
  double sw = 0, sphi = 0, stbin = 0;
  for (const Hit& h : hits)
  {
    const double w = h.adc > 0.0 ? h.adc : 1.0e-3;
    sw += w;
    sphi += w * unwrapNear(m_padMap->get_phi(side, layer, h.pad), phiPeak);
    stbin += w * h.tbin;
  }
  m.phi = sphi / sw;
  m.tbin = stbin / sw;
  m.sector = peak.pad / m_padMap->get_pads_per_sector_for_layer(layer);
  m.ok = m.sector < 12;
  return m;
}

// one reference point: the layer mean drifted back once with the crossing's time offset
bool Tpc_PolyTripletReco::makePoint(const LayerMean& m, unsigned int side, short crossing, RefPoint& p) const
{
  if (!m.ok)
  {
    return false;
  }
  const double timeNs = (m.tbin - m_t0Side[side]) * m_tpcAdcClock - crossing * m_crossingPeriodNs;
  p.layer = m.layer;
  return drift(m.layer, side, m.sector, m.phi, timeNs, p.x, p.y, p.z);
}

bool Tpc_PolyTripletReco::insideTpc(const RefPoint& p, unsigned int side) const
{
  const double r = std::hypot(p.x, p.y);
  const bool inVolume = r >= m_padMap->get_radius(FirstLayer) - m_radialTolerance &&
                        r <= m_padMap->get_radius(FirstLayer + NLayers - 1) + m_radialTolerance &&
                        std::fabs(p.z) <= m_tpcHalfLength + m_zTolerance;
  const bool rightSide = side == 0 ? p.z <= m_centralMembraneTolerance : p.z >= -m_centralMembraneTolerance;
  return inVolume && rightSide;
}

// INTT crossings of the event (time buckets of the INTT cluster sets)
std::set<short> Tpc_PolyTripletReco::inttCrossings() const
{
  std::set<short> crossings;
  if (!m_clusterMap)
  {
    return crossings;
  }
  for (const auto& hitsetkey : m_clusterMap->getHitSetKeys())
  {
    if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(hitsetkey)) == TrkrDefs::TrkrId::inttId)
    {
      crossings.insert(static_cast<short>(InttDefs::getTimeBucketId(hitsetkey)));
    }
  }
  return crossings;
}

int Tpc_PolyTripletReco::process_event(PHCompositeNode* topNode)
{
  const auto tStart = std::chrono::steady_clock::now();
  m_trajectories->Reset();

  std::set<short> crossings;
  if (m_triggeredMode)
  {
    crossings.insert(m_triggeredCrossing);
  }
  else
  {
    m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    crossings = inttCrossings();
  }

  for (unsigned int it = 0; it < m_assembledTracks->size(); ++it)
  {
    const Tpc_AssembledTrack* assembled = m_assembledTracks->get_track(it);
    if (!assembled)
    {
      continue;
    }
    const auto side = static_cast<unsigned int>(assembled->get_side());

    // hits of this track, by layer, and its time extrema
    std::map<unsigned int, std::vector<Hit>> byLayer;
    double minTbin = 1e9, maxTbin = -1e9;
    for (unsigned int ih = 0; ih < assembled->size_hit_indices(); ++ih)
    {
      const Tpc_AssembledTrack::HitIndex hi = assembled->get_hit_index(ih);
      const unsigned int layer = TrkrDefs::getLayer(hi.first);
      if (TpcDefs::getSide(hi.first) != side || layer < FirstLayer || layer >= FirstLayer + NLayers)
      {
        continue;
      }
      TrkrHitSet* hitset = m_hits->findHitSet(hi.first);
      TrkrHit* hit = hitset ? hitset->getHit(hi.second) : nullptr;
      if (!hit)
      {
        continue;
      }
      const auto tbin = static_cast<double>(TpcDefs::getTBin(hi.second));
      minTbin = std::min(minTbin, tbin);
      maxTbin = std::max(maxTbin, tbin);
      byLayer[layer].push_back({hi.first, hi.second, TpcDefs::getPad(hi.second), tbin, static_cast<double>(hit->getAdc())});
    }
    if (byLayer.empty())
    {
      continue;
    }
    ++m_nTracks;

    // usable layers (independent of the crossing): not the first / last layer of a module,
    // no hit on a sector-edge pad
    std::vector<unsigned int> usable;
    for (const auto& [layer, hits] : byLayer)
    {
      if (m_moduleEdgeLayer[layer - FirstLayer])
      {
        continue;
      }
      const unsigned int pps = m_padMap->get_pads_per_sector_for_layer(layer);
      const bool edgePad = std::any_of(hits.begin(), hits.end(), [&](const Hit& h)
                                       { const unsigned int local = h.pad % pps;
                                         return local < m_edgePadMargin || local + m_edgePadMargin >= pps; });
      if (!edgePad)
      {
        usable.push_back(layer);
      }
    }
    if (usable.size() < 3 || usable.back() - usable.front() < m_minLayerSpan)
    {
      ++m_nShortSpan;
      continue;
    }

    // crossing-independent: layer means and the radius of every usable layer
    std::vector<LayerMean> means;
    means.reserve(usable.size());
    for (const unsigned int layer : usable)
    {
      means.push_back(layerMean(layer, side, byLayer[layer]));
    }
    // middle-layer order for a (beginning, end) pair, cached (almost always the same pair)
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> middleOrder;
    const double tol = m_edgeToleranceTbins * m_tpcAdcClock;

    bool anyTrajectory = false;
    for (const short crossing : crossings)
    {
      ++m_nHypotheses;

      // streaming mode: drift-time window first (crossing-independent hit times, no drift needed)
      if (!m_triggeredMode)
      {
        const double tMin = (minTbin - m_t0Side[side]) * m_tpcAdcClock - crossing * m_crossingPeriodNs;
        const double tMax = (maxTbin - m_t0Side[side]) * m_tpcAdcClock - crossing * m_crossingPeriodNs;
        if (!(tMin > -tol && tMax <= m_maxLookupTimeNs + tol))
        {
          ++m_nTimeWindow;
          continue;
        }
      }

      // beginning: innermost usable layer that drifts; end: outermost one; middle: closest to
      // the radial midpoint
      std::array<RefPoint, 3> pts;
      bool okIn = false, okOut = false, okMid = false;
      std::size_t iIn = 0, iOut = means.size();
      for (; iIn < means.size() && !okIn; ++iIn)
      {
        ++m_nDrifts;
        okIn = makePoint(means[iIn], side, crossing, pts[0]);
      }
      for (; okIn && iOut-- > iIn && !okOut;)
      {
        ++m_nDrifts;
        okOut = makePoint(means[iOut], side, crossing, pts[2]);
      }
      if (okIn && okOut && pts[2].layer - pts[0].layer >= m_minLayerSpan)
      {
        auto [orderIt, inserted] = middleOrder.try_emplace({iIn - 1, iOut});
        std::vector<std::size_t>& order = orderIt->second;
        if (inserted)
        {
          const double rMid = 0.5 * (m_padMap->get_radius(pts[0].layer) + m_padMap->get_radius(pts[2].layer));
          for (std::size_t i = iIn; i < iOut; ++i)
          {
            order.push_back(i);
          }
          std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
                    { return std::fabs(m_padMap->get_radius(means[a].layer) - rMid) < std::fabs(m_padMap->get_radius(means[b].layer) - rMid); });
        }
        for (std::size_t k = 0; k < order.size() && !okMid; ++k)
        {
          ++m_nDrifts;
          okMid = makePoint(means[order[k]], side, crossing, pts[1]);
        }
      }
      if (!okMid)
      {
        ++m_nNoPoints;
        continue;
      }

      auto* trajectory = new Tpc_CrossingTrajectoryv1();
      fitTrajectory(pts, *trajectory);

      // streaming mode: physical crossing only
      if (!m_triggeredMode)
      {
        const bool volumeOk = insideTpc(pts[0], side) && insideTpc(pts[1], side) && insideTpc(pts[2], side);
        const bool vertexOk = trajectory->get_fit_status() == 1 &&
                              std::fabs(trajectory->get_z() - m_collisionZ) <= m_maxVertexZ;
        if (!volumeOk || !vertexOk)
        {
          ++(!volumeOk ? m_nOutsideTpc : m_nVertexZ);
          delete trajectory;
          continue;
        }
      }

      trajectory->set_event(m_event);
      trajectory->set_trajectory_id(m_trajectories->size());
      trajectory->set_source_assembled_track_id(assembled->get_track_id());
      trajectory->set_crossing(crossing);
      trajectory->set_side(static_cast<int>(side));
      for (unsigned int k = 0; k < 3; ++k)
      {
        trajectory->set_point(k, pts[k].layer, pts[k].x, pts[k].y, pts[k].z);
      }
      m_trajectories->add_trajectory(trajectory);
      ++m_nTrajectories;
      anyTrajectory = true;
    }
    m_nTracksWithTrajectory += anyTrajectory ? 1 : 0;
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event << " crossings=" << crossings.size()
              << " assembled=" << m_assembledTracks->size() << " trajectories=" << m_trajectories->size() << std::endl;
  }
  ++m_event;
  m_processMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tStart).count();
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyTripletReco::End(PHCompositeNode* /*topNode*/)
{
  std::cout << Name() << "::End - tracks with hits " << m_nTracks << ", with >= 1 trajectory " << m_nTracksWithTrajectory
            << ", usable span < " << m_minLayerSpan << " layers " << m_nShortSpan << std::endl
            << Name() << "::End - crossing hypotheses " << m_nHypotheses << ": trajectories " << m_nTrajectories
            << ", no 3 points " << m_nNoPoints;
  if (!m_triggeredMode)
  {
    std::cout << ", outside the drift time window " << m_nTimeWindow << ", outside the TPC " << m_nOutsideTpc
              << ", |z0| > " << m_maxVertexZ << " cm " << m_nVertexZ;
  }
  std::cout << std::endl
            << Name() << "::End - drift lookup " << m_lookupMs / 1000.0 << " s (InitRun); process_event "
            << (m_event > 0 ? m_processMs / m_event : 0.0) << " ms per event, " << m_nDrifts << " drift interpolations" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

// circle (sagitta) in xy + line in rz through the three reference points, world frame
void Tpc_PolyTripletReco::fitTrajectory(const std::array<RefPoint, 3>& pts, Tpc_CrossingTrajectoryv1& out) const
{
  double x[3], y[3], z[3];
  for (int k = 0; k < 3; ++k)
  {
    const Acts::Vector3 w = m_geometry->transformTpcEnvelopeToWorld(Acts::Vector3(pts[k].x, pts[k].y, pts[k].z));
    x[k] = w.x();
    y[k] = w.y();
    z[k] = w.z();
  }

  // ---- xy: kappa = 2 h / (|a| |b|), h = distance of P2 from the chord P1P3
  const double ax = x[1] - x[0], ay = y[1] - y[0];  // P1 -> P2
  const double bx = x[2] - x[1], by = y[2] - y[1];  // P2 -> P3
  const double cx = x[2] - x[0], cy = y[2] - y[0];  // chord P1 -> P3
  const double la = std::hypot(ax, ay), lb = std::hypot(bx, by), lc = std::hypot(cx, cy);
  if (!(la > 0.0 && lb > 0.0 && lc > 0.0))
  {
    return;
  }
  const double kappa = 2.0 * (ax * by - ay * bx) / (la * lb * lc);  // > 0: counter-clockwise

  // tangent at P1 = chord turned by half the subtended angle
  const double alpha = std::asin(std::clamp(0.5 * kappa * lc, -1.0, 1.0));
  const double ux = cx / lc, uy = cy / lc;
  const double tx = ux * std::cos(alpha) + uy * std::sin(alpha);
  const double ty = -ux * std::sin(alpha) + uy * std::cos(alpha);
  const double nx = -ty, ny = tx;  // left normal

  // closest approach to the beam (0,0), arc from P1 (negative = inward)
  const double wa = -(x[0] * nx + y[0] * ny);
  const double wb = -(x[0] * tx + y[0] * ty);
  const bool line = std::fabs(kappa) < m_minCurvature;
  const double phis = line ? 0.0 : std::atan2(wb * kappa, 1.0 - wa * kappa);
  const double s1 = line ? wb : std::sin(phis) / kappa;
  const double s2 = line ? 0.0 : 2.0 * std::pow(std::sin(0.5 * phis), 2) / kappa;
  const double pcaX = x[0] + tx * s1 + nx * s2;
  const double pcaY = y[0] + ty * s1 + ny * s2;
  const double dirX = tx * std::cos(phis) + nx * std::sin(phis);
  const double dirY = ty * std::cos(phis) + ny * std::sin(phis);

  // ---- rz: least-squares line z = z0 + slope * r
  double r[3], sr = 0, sz = 0, srr = 0, srz = 0;
  for (int k = 0; k < 3; ++k)
  {
    r[k] = std::hypot(x[k], y[k]);
    sr += r[k];
    sz += z[k];
    srr += r[k] * r[k];
    srz += r[k] * z[k];
  }
  const double det = 3.0 * srr - sr * sr;
  if (det <= 0.0)
  {
    return;
  }
  const double slope = (3.0 * srz - sr * sz) / det;
  const double z0 = (sz - slope * sr) / 3.0;
  double chi2 = 0.0;
  for (int k = 0; k < 3; ++k)
  {
    chi2 += std::pow(z[k] - z0 - slope * r[k], 2);  // [cm^2], no errors
  }
  const double pcaZ = z0 + slope * std::hypot(pcaX, pcaY);

  out.set_fit_status(1);
  out.set_x(pcaX);
  out.set_y(pcaY);
  out.set_z(pcaZ);
  out.set_phi(std::atan2(dirY, dirX));
  out.set_slope(slope);
  out.set_chi2(chi2);
  if (line)
  {
    const double norm = std::sqrt(1.0 + slope * slope);
    out.set_px(dirX / norm);
    out.set_py(dirY / norm);
    out.set_pz(slope / norm);
    out.set_charge(0.0);
    out.set_q_over_r(0.0);
  }
  else
  {
    const double pt = 0.003 * m_magneticFieldTesla / std::fabs(kappa);
    const double charge = kappa >= 0.0 ? -1.0 : 1.0;  // same convention as Tpc_PolyTrackReco
    out.set_px(pt * dirX);
    out.set_py(pt * dirY);
    out.set_pz(pt * slope);
    out.set_charge(charge);
    out.set_q_over_r(charge * std::fabs(kappa));
  }
}