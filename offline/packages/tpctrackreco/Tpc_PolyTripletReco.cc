#include "Tpc_PolyTripletReco.h"

#include "IdealPadMap.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include "Tpc_AssembledTrack.h"
#include "Tpc_AssembledTrackContainer.h"
#include "Tpc_PolyClusterContainerv1.h"
#include "Tpc_PolyClusterv1.h"
#include "Tpc_PolyTrackContainerv1.h"
#include "Tpc_PolyTrackv1.h"

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
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TPolyLine3D.h>

#include <algorithm>
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
  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_geomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  m_garfield = findNode::getClass<PHGarfield>(topNode, "PHGarfield");
  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_assembledTracks || !m_crossingDecisions || !m_hits || !m_geomContainer || !m_garfield || !m_geometry)
  {
    std::cerr << Name() << "::InitRun - missing input (" << m_inputNodeName << ", " << m_crossingDecisionNodeName
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

  // drift polylines are built per layer on first use
  m_lookup.assign(NLayers * 2 * 12 * NSamples, Polyline());
  m_layerState.fill(0);

  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  m_points = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_pointNodeName);
  if (!m_points)
  {
    m_points = new Tpc_PolyClusterContainerv1();
    dstNode->addNode(new PHIODataNode<PHObject>(m_points, m_pointNodeName, "PHObject"));
  }
  m_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_trackNodeName);
  if (!m_tracks)
  {
    m_tracks = new Tpc_PolyTrackContainerv1();
    dstNode->addNode(new PHIODataNode<PHObject>(m_tracks, m_trackNodeName, "PHObject"));
  }

  std::cout << Name() << "::InitRun - edge pad margin " << m_edgePadMargin << ", min layer span " << m_minLayerSpan
            << "; t0 side0=" << m_t0Side[0] << " side1=" << m_t0Side[1] << " tbins; output "
            << m_pointNodeName << ", " << m_trackNodeName << std::endl;
  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

// Reverse-drift polylines from NSamples pad positions per sector, one layer, both sides
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
      for (unsigned int s = 0; s < NSamples; ++s)
      {
        const auto localPad = static_cast<unsigned int>(std::llround(
            static_cast<double>(s) / (NSamples - 1) * static_cast<double>(padsPerSector - 1)));
        const double phi = m_padMap->get_phi(side, sector, layer, localPad);
        TPolyLine3D* line = m_garfield->ReverseDrift(radius * std::cos(phi), radius * std::sin(phi), m_startZ[side],
                                                     m_driftStepNs, nullptr, m_driftSubsteps);
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
  std::array<double, NSamples> sphi{};
  for (unsigned int s = 0; s < NSamples; ++s)
  {
    sphi[s] = (s == 0) ? polyline(layer, side, sector, 0).phi : unwrapNear(polyline(layer, side, sector, s).phi, sphi[s - 1]);
  }
  phi = unwrapNear(phi, sphi[NSamples / 2]);
  const double dir = (sphi[NSamples - 1] >= sphi[0]) ? 1.0 : -1.0;
  unsigned int s0 = 0;
  while (s0 + 2 < NSamples && dir * (phi - sphi[s0 + 1]) > 0.0)
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

// one reference point: ADC-weighted pad phi and time of the hits on the layer, drifted back once
bool Tpc_PolyTripletReco::makePoint(unsigned int layer, unsigned int side, short crossing,
                                    const std::vector<Hit>& hits, RefPoint& p)
{
  int& state = m_layerState[layer - FirstLayer];
  if (state == 0)
  {
    state = buildLayer(layer) ? 1 : -1;
  }
  if (state < 0 || hits.empty())
  {
    return false;
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
  const double timeNs = (stbin / sw - m_t0Side[side]) * m_tpcAdcClock - crossing * m_crossingPeriodNs;
  const unsigned int sector = peak.pad / m_padMap->get_pads_per_sector_for_layer(layer);
  p.layer = layer;
  return sector < 12 && drift(layer, side, sector, sphi / sw, timeNs, p.x, p.y, p.z);
}

int Tpc_PolyTripletReco::storePoint(unsigned int assembledId, unsigned int side, const RefPoint& p,
                                    const std::vector<Hit>& hits)
{
  const Hit& peak = *std::max_element(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.adc < b.adc; });
  auto* point = new Tpc_PolyClusterv1();
  point->set_event(m_event);
  point->set_cluster_id(m_points->size());
  point->set_source_assembled_track_id(assembledId);
  point->set_trkr_cluster_key(TrkrDefs::genClusKey(peak.hitsetkey, m_nextIndex[peak.hitsetkey]++));
  point->set_side(static_cast<int>(side));
  point->set_centroid_x(p.x);
  point->set_centroid_y(p.y);
  point->set_centroid_z(p.z);
  double adc = 0;
  for (const Hit& h : hits)
  {
    adc += std::max(h.adc, 0.0);
    point->add_hit(h.hitsetkey, h.hitkey, p.x, p.y, p.z);  // hit list only; xyz = reference point
  }
  point->set_adc(adc);
  const int index = static_cast<int>(m_points->size());
  m_points->add_cluster(point);
  return index;
}

int Tpc_PolyTripletReco::process_event(PHCompositeNode* /*topNode*/)
{
  m_points->Reset();
  m_tracks->Reset();
  m_nextIndex.clear();

  for (unsigned int it = 0; it < m_assembledTracks->size(); ++it)
  {
    const Tpc_AssembledTrack* assembled = m_assembledTracks->get_track(it);
    if (!assembled)
    {
      continue;
    }
    const TpcCrossingDecision* decision = m_crossingDecisions->get_decision(assembled->get_track_id());
    if (!decision || decision->get_selected_tier() > m_maxAcceptedTier)
    {
      continue;
    }
    const short crossing = decision->get_selected_crossing();
    const auto side = static_cast<unsigned int>(assembled->get_side());

    // hits of this track, by layer
    std::map<unsigned int, std::vector<Hit>> byLayer;
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
      byLayer[layer].push_back({hi.first, hi.second, TpcDefs::getPad(hi.second),
                                static_cast<double>(TpcDefs::getTBin(hi.second)), static_cast<double>(hit->getAdc())});
    }
    if (byLayer.empty())
    {
      continue;
    }
    ++m_nTracks;

    // usable layers: not the first / last layer of a module, no hit on a sector-edge pad
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

    // beginning: innermost usable layer that drifts; end: outermost one
    std::array<RefPoint, 3> pts;
    bool ok = false;
    std::size_t iIn = 0;
    for (; iIn < usable.size() && !ok; ++iIn)
    {
      ok = makePoint(usable[iIn], side, crossing, byLayer[usable[iIn]], pts[0]);
    }
    bool okOut = false;
    for (std::size_t i = usable.size(); ok && i-- > iIn && !okOut;)
    {
      okOut = makePoint(usable[i], side, crossing, byLayer[usable[i]], pts[2]);
    }
    if (!ok || !okOut)
    {
      ++m_nDriftFailed;
      continue;
    }
    if (pts[2].layer - pts[0].layer < m_minLayerSpan)
    {
      ++m_nShortSpan;
      continue;
    }

    // middle: usable layer closest to the radial midpoint of the two
    const double rMid = 0.5 * (m_padMap->get_radius(pts[0].layer) + m_padMap->get_radius(pts[2].layer));
    std::vector<unsigned int> middle;
    for (const unsigned int layer : usable)
    {
      if (layer > pts[0].layer && layer < pts[2].layer)
      {
        middle.push_back(layer);
      }
    }
    std::sort(middle.begin(), middle.end(), [&](unsigned int l1, unsigned int l2)
              { return std::fabs(m_padMap->get_radius(l1) - rMid) < std::fabs(m_padMap->get_radius(l2) - rMid); });
    bool okMid = false;
    for (std::size_t i = 0; i < middle.size() && !okMid; ++i)
    {
      okMid = makePoint(middle[i], side, crossing, byLayer[middle[i]], pts[1]);
    }
    if (!okMid)
    {
      ++m_nNoMiddle;
      continue;
    }

    ++m_nTriplets;
    std::array<int, 3> pointIndex{};
    for (int k = 0; k < 3; ++k)
    {
      pointIndex[k] = storePoint(assembled->get_track_id(), side, pts[k], byLayer[pts[k].layer]);
    }
    fillTrack(assembled->get_track_id(), pointIndex);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event << " assembled=" << m_assembledTracks->size()
              << " points=" << m_points->size() << " tracks=" << m_tracks->size() << std::endl;
  }
  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyTripletReco::End(PHCompositeNode* /*topNode*/)
{
  std::cout << Name() << "::End - tracks with hits " << m_nTracks << ": triplets " << m_nTriplets
            << ", usable span < " << m_minLayerSpan << " layers "
            << m_nShortSpan << ", no inner/outer point " << m_nDriftFailed << ", no middle point " << m_nNoMiddle
            << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

// circle (sagitta) in xy + line in rz through the three reference points, world frame
void Tpc_PolyTripletReco::fillTrack(unsigned int assembledId, const std::array<int, 3>& pointIndex)
{
  auto* out = new Tpc_PolyTrackv1();
  out->set_event(m_event);
  out->set_track_id(m_tracks->size());
  out->set_source_assembled_track_id(assembledId);

  double x[3], y[3], z[3];
  bool complete = true;
  for (int k = 0; k < 3; ++k)
  {
    if (pointIndex[k] < 0)
    {
      complete = false;
      continue;
    }
    const Tpc_PolyCluster* c = m_points->get_cluster(pointIndex[k]);
    out->add_cluster_key(c->get_trkr_cluster_key());
    const Acts::Vector3 w = m_geometry->transformTpcEnvelopeToWorld(
        Acts::Vector3(c->get_centroid_x(), c->get_centroid_y(), c->get_centroid_z()));
    x[k] = w.x();
    y[k] = w.y();
    z[k] = w.z();
  }

  if (complete)
  {
    // ---- xy: kappa = 2 h / (|a| |b|), h = distance of P2 from the chord P1P3
    const double ax = x[1] - x[0], ay = y[1] - y[0];  // P1 -> P2
    const double bx = x[2] - x[1], by = y[2] - y[1];  // P2 -> P3
    const double cx = x[2] - x[0], cy = y[2] - y[0];  // chord P1 -> P3
    const double la = std::hypot(ax, ay), lb = std::hypot(bx, by), lc = std::hypot(cx, cy);
    if (la > 0.0 && lb > 0.0 && lc > 0.0)
    {
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
      if (det > 0.0)
      {
        const double slope = (3.0 * srz - sr * sz) / det;
        const double z0 = (sz - slope * sr) / 3.0;
        double chi2 = 0.0;
        for (int k = 0; k < 3; ++k)
        {
          chi2 += std::pow(z[k] - z0 - slope * r[k], 2);  // [cm^2], no errors
        }
        const double pcaZ = z0 + slope * std::hypot(pcaX, pcaY);

        out->set_fit_status(1);
        out->set_x(pcaX);
        out->set_y(pcaY);
        out->set_z(pcaZ);
        out->set_seed_x0(pcaX);
        out->set_seed_y0(pcaY);
        out->set_seed_z0(pcaZ);
        out->set_seed_phi(std::atan2(dirY, dirX));
        out->set_seed_slope(slope);
        out->set_chi2(chi2);
        out->set_ndf(1.0);  // circle through 3 points is exact; rz line has 1 dof
        if (line)
        {
          const double norm = std::sqrt(1.0 + slope * slope);
          out->set_px(dirX / norm);
          out->set_py(dirY / norm);
          out->set_pz(slope / norm);
          out->set_charge(0.0);
          out->set_seed_q_over_r(0.0);
        }
        else
        {
          const double pt = 0.003 * m_magneticFieldTesla / std::fabs(kappa);
          const double charge = kappa >= 0.0 ? -1.0 : 1.0;  // same convention as Tpc_PolyTrackReco
          out->set_px(pt * dirX);
          out->set_py(pt * dirY);
          out->set_pz(pt * slope);
          out->set_charge(charge);
          out->set_seed_q_over_r(charge * std::fabs(kappa));
          out->set_helix_x0(x[0] + nx / kappa);
          out->set_helix_y0(y[0] + ny / kappa);
        }
      }
    }
  }

  m_tracks->add_track(out);
}