#include "SiHitSeedReco.h"

#include <tpctrackreco/Tpc_CrossingTrajectory.h>
#include <tpctrackreco/Tpc_CrossingTrajectoryContainer.h>
#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <intt/CylinderGeomIntt.h>
#include <mvtx/CylinderGeom_Mvtx.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <Eigen/Eigenvalues>

#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double twopi = 2.0 * pi;
constexpr std::array<int, 3> nstaves{{12, 16, 20}};
constexpr std::array<double, 3> mvtx_r{{2.523, 3.336, 4.148}};
constexpr std::array<double, 3> mvtx_phi0{{0.2382, 0.1937, 0.1481}};
constexpr std::array<double, 3> mvtx_tilt{{13.3 * pi / 180., 16.9 * pi / 180., 17.0 * pi / 180.}};
constexpr std::array<int, 4> intt_nladder{{12, 12, 16, 16}};
constexpr std::array<double, 4> intt_phi0{{-15 * pi / 180., 0., -11.25 * pi / 180., 0.}};
constexpr std::array<double, 4> intt_r{{7.188 - 0.0036, 7.732 - 0.0036, 9.680 - 0.0036, 10.262 - 0.0036}};
constexpr std::array<double, 9> mvtx_chip_z{{-12.06, -9.045, -6.03, -3.015, 0., 3.015, 6.03, 9.045, 12.06}};
constexpr double intt_inner_z = (8 * 1.6 + 0.2) / 2.;
constexpr double intt_outer_z = (8 * 1.6 + 0.2) + (5 * 2.0 + 0.2) / 2.;
constexpr std::array<double, 4> intt_sensor_z{{-intt_inner_z, -intt_outer_z, intt_inner_z, intt_outer_z}};
constexpr int nlayers = 7;

const char* statusName(int s)
{
  switch (s)
  {
  case SiTpcSeedInfo::Chained:
    return "chained";
  case SiTpcSeedInfo::TooFewClusters:
    return "too few TPC clusters";
  case SiTpcSeedInfo::FitFailed:
    return "TPC helix fit failed";
  case SiTpcSeedInfo::NoPrediction:
    return "no prediction on the seed layers";
  case SiTpcSeedInfo::NoAnchor:
    return "no hit in the anchor window";
  case SiTpcSeedInfo::NoChain:
    return "no chain passing the layer requirements";
  case SiTpcSeedInfo::Rejected:
    return "best chain fails the agreement cut";
  default:
    return "?";
  }
}
}  // namespace

SiHitSeedReco::SiHitSeedReco(const std::string& name)
  : SubsysReco(name)
{
  m_map.frame = &m_frame;
  m_map.alignment = &m_beamAlign;
}

int SiHitSeedReco::InitRun(PHCompositeNode* topNode)
{
  m_map.frame = &m_frame;
  m_map.alignment = &m_beamAlign;
  if (Verbosity() > 0)
  {
    const auto& c = m_frame.detectorCenterCm();
    std::cout << Name() << ": TPC-driven chaining, TPC from "
              << (m_tpcTrajectoryNodeName.empty() ? m_tpcTrackNodeName + " + " + m_tpcClusterNodeName
                                                  : m_tpcTrajectoryNodeName + " (crossing trajectories)")
              << " (>= " << m_minTpcClusters << " clusters, z offset " << m_tpcZOffset << " cm)" << std::endl
              << Name() << ":   frames: Uphi rotation " << m_frame.uphiRotation() << " rad, detector centre ("
              << 10 * c[0] << ", " << 10 * c[1] << ", " << 10 * c[2] << ") mm, INTT centre ("
              << 10 * m_frame.inttCenterCm()[0] << ", " << 10 * m_frame.inttCenterCm()[1] << ", " << 10 * m_frame.inttCenterCm()[2]
              << ") mm" << (m_frame.hasInttCenter() ? "" : " (= MVTX)") << ", Si halves on TPC beam axis "
              << (m_beamAlign.enabled() ? "on" : "off") << ", Bz " << m_bz << " T" << std::endl
              << Name() << ":   silicon hit positions on the " << (m_flatSensors ? "FLAT sensors (tilted MVTX staves)" : "nominal CYLINDERS (approximation)")
              << (m_applyAlignment ? " + per-sensor alignment" : ", no per-sensor alignment") << std::endl
              << Name() << ":   anchor window |r dphi| < " << m_anchorRphi << " cm, |dz| < " << m_anchorZ
              << " cm (max " << m_maxAnchors << " anchors); seed pairs 2->1" << (m_useInnerPair ? ", 1->0" : "") << std::endl
              << Name() << ":   seed window: z centre tanl_TPC * ds + " << m_seedZOffsetCm << " cm, +- " << m_seedZHalfCm << " cm; phi centre "
              << (m_seedPhiFromTpc ? "TPC curvature" : "0") << ", +- " << m_seedPhiHalf << " bins" << std::endl
              << Name() << ":   propagation (hardware coordinates): phi +- " << m_propPhiHalf << " bins, z +- "
              << m_propZHalf << " bins; chains >= " << m_minLayers << " layers, >= "
              << m_minMvtxLayers << " MVTX, >= " << m_minInttLayers << " INTT"
              << ((m_inttCrossingWindow >= 0 && !m_tpcTrajectoryNodeName.empty())
                      ? "; INTT hits only within +-" + std::to_string(m_inttCrossingWindow) + " time buckets of the trajectory crossing"
                      : std::string(""))
              << std::endl
              << Name() << ":   selection: " << (m_chainRanking == ChainRanking::AgreementFirst ? "agreeing chains first (chi2/ndf <= " + std::to_string(m_agreementMaxChi2) + "), then most layers" : std::string("most layers, then chi2/ndf"))
              << (m_exclusiveChains ? ", exclusive across TPC tracks (>= " + std::to_string(m_exclusiveMinShared) + " shared blobs)" : std::string(""))
              << " (sigma r*phi L0-6 " << m_sigRphiL[0] << " "
              << m_sigRphiL[1] << " " << m_sigRphiL[2] << " " << m_sigRphiL[3] << " " << m_sigRphiL[4] << " " << m_sigRphiL[5] << " "
              << m_sigRphiL[6] << ", z L0-6 " << m_sigZL[0] << " " << m_sigZL[1] << " " << m_sigZL[2] << " " << m_sigZL[3] << " "
              << m_sigZL[4] << " " << m_sigZL[5] << " " << m_sigZL[6] << " cm; INTT window centre "
              << (m_inttCentre == InttCentre::MvtxLine ? "straight MVTX line"
                  : (m_inttCentre == InttCentre::TpcOffset ? "TPC prediction + mean MVTX offset" : "TPC prediction + MVTX offset and slope"))
              << ", z model "
              << (m_zModel == TpcZModel::Fixed ? std::string("fixed")
                  : m_zModel == TpcZModel::Offset ? std::string("free offset")
                  : m_zModel == TpcZModel::OffsetSlope ? std::string("free offset + slope")
                  : "offset with prior sigma " + std::to_string(m_zOffsetSigma) + " cm")
              << (m_mvtxPropZFromTpc ? ", MVTX propagation z from the TPC inclination" : "")
              << ", INTT z half window " << (m_inttPropZHalf >= 0 ? m_inttPropZHalf : m_propZHalf) << " bins"
              << ", z weight " << m_zWeight << "; anchor z " << (m_anchorZ > 0 ? "window" : "not used") << ", z ordering weight " << m_anchorZWeight
              << ", max chi2/ndf " << (m_maxTpcChi2 > 0 ? std::to_string(m_maxTpcChi2) : std::string("off")) << ")"
              << std::endl;
  }
  if (m_useActsGeometry || m_compareActs)
  {
    m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    m_geomMvtx = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
    m_geomIntt = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
    if (!m_tGeometry || !m_geomMvtx || !m_geomIntt)
    {
      std::cerr << Name() << ": setUseActsGeometry / setCompareWithActsGeometry need ActsGeometry, CYLINDERGEOM_MVTX and CYLINDERGEOM_INTT"
                << " (TrackingInit before this module)" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    std::cout << Name() << ": silicon physics points from "
              << (m_useActsGeometry ? "the OFFICIAL tracking geometry (ActsGeometry), z shifted by " + std::to_string(m_actsZShift) + " cm"
                                    : std::string("own geometry"))
              << (m_compareActs ? ", compared with ActsGeometry at End" : "") << std::endl;
    std::cout << Name() << ": silicon HARDWARE coordinates (blobs, chaining) from "
              << (m_useActsGeometry && m_hwFromActs ? std::string("the official positions (Uphi about the MVTX / INTT centre)")
                                                    : std::string("the analytic stave / ladder model"))
              << std::endl;
  }
  if (m_inttPhysOffset[0] != 0.0 || m_inttPhysOffset[1] != 0.0)
  {
    std::cout << Name() << ": INTT physics points shifted by (" << 10 * m_inttPhysOffset[0] << ", " << 10 * m_inttPhysOffset[1] << ") mm" << std::endl;
  }
  if (std::any_of(m_layerZOffset.begin(), m_layerZOffset.end(), [](double v) { return v != 0.0; }))
  {
    std::cout << Name() << ": per-layer z offsets of the physics points [mm]:";
    for (int l = 0; l < nlayers; ++l)
    {
      std::cout << " L" << l << " " << 10 * m_layerZOffset[l];
    }
    std::cout << std::endl;
  }
  if (m_applyAlignment)
  {
    const int ret = loadAlignment();
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return ret;
    }
  }
  return createNodes(topNode);
}

int SiHitSeedReco::loadAlignment()
{
  std::unordered_map<TrkrDefs::hitsetkey, SiAlignmentTransform::Params> params;
  const std::string& source = m_alignmentFile;
  const int nread = source.empty() ? -1 : m_alignment.readFile(source, params, Verbosity());
  if (nread < 0)
  {
    std::cerr << Name() << ": cannot read silicon alignment parameters from '" << source
              << "' (set it with setSiAlignmentParamsFile)" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Ideal sensor centres/axes in the DETECTOR frame, from the same geometry that fillHits uses
  // (sensorPoint / sensorAxes: local (0, 0) of the sensor, local x along the sensor width, local
  // y along +z).  The axes matter only for the new rotation order.
  m_alignment.clear();
  for (const auto& [key, p] : params)
  {
    const int layer = TrkrDefs::getLayer(key);
    if (layer < 0 || layer >= nlayers)
    {
      continue;
    }
    const bool mvtx = layer <= 2;
    const int stave = mvtx ? MvtxDefs::getStaveId(key) : -1;
    const int chip = mvtx ? MvtxDefs::getChipId(key) : -1;
    const int ladderphi = mvtx ? -1 : InttDefs::getLadderPhiId(key);
    const int ladderz = mvtx ? -1 : InttDefs::getLadderZId(key);
    if ((mvtx && (stave >= nstaves[layer] || chip > 8)) ||
        (!mvtx && (ladderphi >= intt_nladder[layer - 3] || ladderz > 3)))
    {
      continue;
    }

    const Eigen::Vector3d centre = sensorPoint(layer, stave, chip, ladderphi, ladderz, 0.0, 0.0);
    const Eigen::Matrix3d ideal = sensorAxes(layer, stave, ladderphi);
    m_alignment.addSensor(SiAlignmentTransform::sensorId(key), p, centre, ideal);
  }

  std::cout << Name() << ": silicon alignment ON, " << m_alignment.size() << " sensors from " << source
            << (m_alignment.useNewSiliconRotationOrder() ? " (new" : " (old") << " rotation order)" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

Eigen::Vector3d SiHitSeedReco::alignPoint(TrkrDefs::hitsetkey hskey, const Eigen::Vector3d& p) const
{
  if (!m_applyAlignment)
  {
    return p;
  }
  const auto a = m_alignment.apply(SiAlignmentTransform::sensorId(hskey), {p.x(), p.y(), p.z()});
  return {a.x, a.y, a.z};
}

void SiHitSeedReco::alignHardware(TrkrDefs::hitsetkey hskey, int layer, double& z, double& phi) const
{
  if (!m_applyAlignment)
  {
    return;
  }
  // as before: nominal-cylinder point -> aligned -> back to (Uphi, Vz); only for the chaining
  const double rot = mvtx_phi0[2];
  const double p = phi + rot;
  const double r = layer <= 2 ? mvtx_r[layer] : intt_r[layer - 3];
  const Eigen::Vector3d a = alignPoint(hskey, {r * std::cos(p), r * std::sin(p), z});
  z = a.z();
  phi = std::fmod(std::atan2(a.y(), a.x()) - rot, twopi);
  if (phi < 0)
  {
    phi += twopi;
  }
}

bool SiHitSeedReco::actsGlobal(TrkrDefs::hitsetkey hskey, int layer, unsigned int row, unsigned int col, int ladderz, Eigen::Vector3d& g) const
{
  if (!m_tGeometry)
  {
    return false;
  }
  // local coordinates exactly as MvtxClusterizer / InttClusterizer (cm)
  double lx = 0;
  double lz = 0;
  if (layer <= 2)
  {
    auto* lg = dynamic_cast<CylinderGeom_Mvtx*>(m_geomMvtx->GetLayerGeom(layer));
    if (!lg)
    {
      return false;
    }
    const TVector3 loc = lg->get_local_coords_from_pixel(static_cast<int>(row), static_cast<int>(col));
    lx = loc.X();
    lz = loc.Z();
  }
  else
  {
    auto* lg = dynamic_cast<CylinderGeomIntt*>(m_geomIntt->GetLayerGeom(layer));
    if (!lg)
    {
      return false;
    }
    double loc[3] = {0.0, 0.0, 0.0};
    lg->find_strip_center_localcoords(ladderz, static_cast<int>(row), static_cast<int>(col), loc);
    lx = loc[1];
    lz = loc[2];
  }
  const auto surf = m_tGeometry->maps().getSiliconSurface(hskey);
  if (!surf)
  {
    return false;
  }
  const Acts::Vector2 local(lx * Acts::UnitConstants::cm, lz * Acts::UnitConstants::cm);
  const Acts::Vector3 glob = surf->localToGlobal(m_tGeometry->geometry().getGeoContext(), local, Acts::Vector3(1, 1, 1)) / Acts::UnitConstants::cm;
  g = Eigen::Vector3d(glob.x(), glob.y(), glob.z());
  return true;
}

// Sensor geometry in the DETECTOR frame (ideal sPHENIX frame before the detector centre).
//  MVTX (PHG4MvtxDetector, CylinderGeom_Mvtx): stave origin at R_layer (layer_nominal_radius)
//    and azimuth phi_s = phi0 + stave * 2 pi / N, stave rotated by phi_s + pi/2 + turbo tilt;
//    sensor centre at stave-local (+0.058128, -0.003) cm (loc_sensor_in_chip + chip/half-stave
//    offsets), hit local x (lx, sensor frame) along the stave-local x axis.  The sensor surface
//    spans r = 24.61 ... 27.93 mm on layer 0 (PHG4MvtxDefs rMin / rMax), not the cylinder 25.23 mm.
//  INTT: flat ladder tangent to R_layer at its azimuth, lx along +phi.
Eigen::Matrix3d SiHitSeedReco::sensorAxes(int layer, int stave, int ladderphi) const
{
  double a = 0;
  if (layer <= 2)
  {
    const double phis = mvtx_phi0[layer] + stave * twopi / nstaves[layer];
    a = phis + 0.5 * pi + (m_flatSensors ? mvtx_tilt[layer] : 0.0);
  }
  else
  {
    const int i = layer - 3;
    a = intt_phi0[i] + twopi * ladderphi / intt_nladder[i] + 0.5 * pi;
  }
  Eigen::Matrix3d m;
  m.col(0) = Eigen::Vector3d(std::cos(a), std::sin(a), 0.0);  // local x
  m.col(1) = Eigen::Vector3d::UnitZ();                       // local y (along z)
  m.col(2) = m.col(0).cross(m.col(1));                       // normal (outward)
  return m;
}

Eigen::Vector3d SiHitSeedReco::sensorPoint(int layer, int stave, int chip, int ladderphi, int ladderz, double lx, double ly) const
{
  const double z = longitudinalZ(layer, chip, ladderz, ly);
  if (!m_flatSensors)
  {
    // previous approximation: nominal cylinder
    const double phi = intrinsicPhi(layer, stave, ladderphi, lx) + mvtx_phi0[2];
    const double r = layer <= 2 ? mvtx_r[layer] : intt_r[layer - 3];
    return {r * std::cos(phi), r * std::sin(phi), z};
  }
  double phic = 0;  // azimuth of the stave / ladder origin
  double r = 0;
  double xs = lx;     // stave-local x
  double ys = 0.0;  // stave-local y (towards the beam for the MVTX after the rotation)
  if (layer <= 2)
  {
    phic = mvtx_phi0[layer] + stave * twopi / nstaves[layer];
    r = mvtx_r[layer];
    xs = lx + 0.058128;                // sensor centre in the chip (CylinderGeom_Mvtx)
    ys = -0.0005 - 0.02075 + 0.01825;  // sensor in chip + chip in module + half-stave in stave
  }
  else
  {
    const int i = layer - 3;
    phic = intt_phi0[i] + twopi * ladderphi / intt_nladder[i];
    r = intt_r[i];
  }
  const Eigen::Matrix3d ax = sensorAxes(layer, stave, ladderphi);
  const Eigen::Vector3d ex = ax.col(0);
  const Eigen::Vector3d ey(-ex.y(), ex.x(), 0.0);  // stave-local +y (Rz(a) * y)
  const Eigen::Vector3d p = r * Eigen::Vector3d(std::cos(phic), std::sin(phic), 0.0) + xs * ex + ys * ey;
  return {p.x(), p.y(), z};
}

int SiHitSeedReco::createNodes(PHCompositeNode* topNode)
{
  m_eventData = findNode::getClass<SiHitSeedEvent>(topNode, m_outputNodeName);
  if (m_eventData)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << ": DST node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_eventData = new SiHitSeedEvent();
  auto* node = new PHDataNode<PHObject>(m_eventData, m_outputNodeName, "PHObject");
  dstNode->addNode(node);

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": created transient node " << m_outputNodeName << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiHitSeedReco::process_event(PHCompositeNode* topNode)
{
  if (!m_eventData)
  {
    m_eventData = findNode::getClass<SiHitSeedEvent>(topNode, m_outputNodeName);
    if (!m_eventData)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  m_eventData->Reset();
  m_eventData->event = m_event++;
  if (!fillHits(topNode))
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (!m_eventData->hits.empty())
  {
    const auto [zlo, zhi] = std::minmax_element(m_eventData->hits.begin(), m_eventData->hits.end(),
                                                [](const auto& a, const auto& b) { return a.z < b.z; });
    const double span = std::max(0.01, zhi->z - zlo->z);
    m_zmin = zlo->z - std::max(0.02 * span, 0.01);
    m_zmax = zhi->z + std::max(0.02 * span, 0.01);
  }

  const double dz = (m_zmax - m_zmin) / m_nz;
  const double dphi = twopi / m_nphi;
  for (auto& h : m_eventData->hits)
  {
    h.zbin = (h.z - m_zmin) / dz;
    h.phibin = h.phi / dphi;
  }

  buildBlobs();
  buildHitGrid();

  // kept for the worker interface (rechain / recoverHits of SiSeedingOrchestrator)
  m_tpcInputs = collectTpcTracks(topNode);
  m_tpcInputPos.clear();
  for (std::size_t i = 0; i < m_tpcInputs.size(); ++i)
  {
    m_tpcInputPos[m_tpcInputs[i].index] = i;
  }
  for (const auto& t : m_tpcInputs)
  {
    processTpcTrack(t);
  }
  if (m_exclusiveChains)
  {
    resolveSharedChains();
  }
  if (m_measureWidths)
  {
    accumulateResidualWidths();
  }

  if (Verbosity() > 0)
  {
    std::array<unsigned int, 7> st{};
    for (const auto& s : m_eventData->tpc_seeds)
    {
      ++st[std::clamp(s.status, 0, 6)];
    }
    unsigned int nSel = 0;
    for (const auto& c : m_eventData->chains)
    {
      nSel += c.selected;
    }
    std::cout << Name() << ": event " << m_eventData->event
              << " hits=" << m_eventData->hits.size()
              << " blobs=" << m_eventData->blobs.size()
              << " tpc_tracks=" << m_eventData->tpc_seeds.size()
              << " chains=" << m_eventData->chains.size()
              << " selected=" << nSel
              << " | chained " << st[SiTpcSeedInfo::Chained]
              << ", <clusters " << st[SiTpcSeedInfo::TooFewClusters]
              << ", fit failed " << st[SiTpcSeedInfo::FitFailed]
              << ", no prediction " << st[SiTpcSeedInfo::NoPrediction]
              << ", no anchor " << st[SiTpcSeedInfo::NoAnchor]
              << ", no chain " << st[SiTpcSeedInfo::NoChain]
              << ", rejected " << st[SiTpcSeedInfo::Rejected] << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiHitSeedReco::End(PHCompositeNode*)
{
  if (m_hwActs > 0)
  {
    std::cout << Name() << ": hardware coordinates from the official positions for " << m_hwActs << " hits, analytic for " << m_hwAnalytic
              << "; analytic model - official per layer (dUphi [mrad], dVz [mm]: mean / rms; hits with |dUphi| > 50 mrad):" << std::endl;
    for (int l = 0; l < nlayers; ++l)
    {
      const auto& c = m_hwCmp[l];
      if (c[0] <= 0)
      {
        continue;
      }
      const double mu = c[1] / c[0], mz = c[3] / c[0];
      std::printf("%s:   layer %d  n %9.0f  dUphi %+8.2f / %7.2f   dVz %+7.3f / %6.3f   |dUphi|>50 mrad: %.0f (%.2f %%)\n", Name().c_str(), l, c[0], 1e3 * mu,
                  1e3 * std::sqrt(std::max(0.0, c[2] / c[0] - mu * mu)), 10 * mz, 10 * std::sqrt(std::max(0.0, c[4] / c[0] - mz * mz)), c[5],
                  100.0 * c[5] / c[0]);
    }
  }
  if (m_compareActs)
  {
    std::cout << Name() << ": own silicon hit positions (" << (m_flatSensors ? "flat sensors" : "cylinders")
              << (m_applyAlignment ? " + " + m_alignmentFile : std::string(", no alignment file"))
              << ", MVTX centre (" << 10 * m_frame.detectorCenterCm()[0] << ", " << 10 * m_frame.detectorCenterCm()[1] << ", " << 10 * m_frame.detectorCenterCm()[2]
              << ") mm, INTT centre (" << 10 * m_frame.inttCenterCm()[0] << ", " << 10 * m_frame.inttCenterCm()[1] << ", " << 10 * m_frame.inttCenterCm()[2]
              << ") mm) - OFFICIAL (ActsGeometry), global frame, mean / rms in mm:" << std::endl;
    std::printf("%s:   layer        n      dx  rms      dy  rms      dz  rms      dr  rms   r*dphi  rms\n", Name().c_str());
    for (int l = 0; l < nlayers; ++l)
    {
      const auto& c = m_actsCmp[l];
      if (c[0] <= 0)
      {
        continue;
      }
      std::printf("%s:   %5d %8.0f", Name().c_str(), l, c[0]);
      for (int k = 0; k < 5; ++k)
      {
        const double m = c[1 + k] / c[0];
        std::printf("  %+6.3f %5.3f", 10 * m, 10 * std::sqrt(std::max(0.0, c[6 + k] / c[0] - m * m)));
      }
      std::printf("\n");
    }
  }
  if (m_rechainCalls > 0)
  {
    std::cout << Name() << ": rechain (pass 2 on free blobs): " << m_rechainCalls << " calls, " << m_rechainInputs
              << " TPC inputs chained again, " << m_rechainChains << " new candidate chains (not in the pass-1 counts below)" << std::endl;
  }
  if (m_recoveryTried > 0)
  {
    unsigned long nAdded = 0;
    for (const auto n : m_recoveredPerLayer)
    {
      nAdded += n;
    }
    std::cout << Name() << ": hit recovery (pass 3): " << m_recoveryTried << " missing layers searched, " << nAdded << " added (L0-6:";
    for (const auto n : m_recoveredPerLayer)
    {
      std::cout << " " << n;
    }
    std::cout << "), no free blob in the window " << m_recoveryNoBlob << ", blob found but chain no longer agrees with the TPC "
              << m_recoveryBadChi2 << "; predictions: vertex-assisted " << m_recoveryVertexPred << ", TPC + offset fallback "
              << m_recoveryFallbackPred << ", none " << m_recoveryNoPred << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": " << m_event << " events, " << m_nChains << " chains stored; TPC tracks:";
    for (int s = 0; s < 7; ++s)
    {
      std::cout << " " << statusName(s) << " " << m_statusCount[s] << (s < 6 ? "," : "");
    }
    std::cout << std::endl;
    if (m_inttCentre != InttCentre::MvtxLine)
    {
      std::cout << Name() << ": INTT window from the TPC prediction ("
                << (m_inttCentre == InttCentre::TpcOffsetSlope ? "MVTX offset + slope" : "mean MVTX offset") << "): " << m_inttTpcUsed << " INTT searches, "
                << m_inttTpcFallback << " fell back to the straight MVTX line (no TPC prediction / MVTX offset)" << std::endl;
    }
    if (m_mvtxPropZFromTpc)
    {
      std::cout << Name() << ": MVTX propagations with the z centre from the TPC inclination: " << m_mvtxTpcZUsed << std::endl;
    }
    if (m_inttCrossingWindow >= 0 && !m_tpcTrajectoryNodeName.empty())
    {
      std::cout << Name() << ": INTT crossing window +-" << m_inttCrossingWindow << ": INTT hits in propagation windows "
                << m_inttTimeChecked << ", rejected by time bucket " << m_inttTimeRejected << std::endl;
    }
    if (m_measureWidths)
    {
      printResidualWidths();
    }
    if (m_exclusiveChains)
    {
      std::cout << Name() << ": exclusive chains (conflict = >= " << m_exclusiveMinShared << " shared blobs): kept rank 0 "
                << m_exclKept << ", moved to another chain " << m_exclReassigned << ", left without chain " << m_exclLost
                << " TPC tracks" << std::endl;
    }
    if (m_chainRanking == ChainRanking::AgreementFirst)
    {
      std::cout << Name() << ": agreement-first ranking (tpc_chi2 <= " << m_agreementMaxChi2 << "): rank 0 changed vs layers-first for "
                << m_rankChanged << " TPC tracks, no agreeing chain for " << m_noAgreeingChain << std::endl;
    }
    if (m_hwDdphiProp0 || m_hwDdphiCheck[0] || m_hwDdphiCheck[1])
    {
      std::cout << Name() << ": hw delta-delta-phi (" << m_hwDdphiNSigma << " sigma): layer-0 propagation windows "
                << m_hwDdphiProp0Used;
      for (int l = 0; l < 2; ++l)
      {
        if (m_hwDdphiCheck[l])
        {
          std::cout << ", layer " << l << " checked " << m_hwDdphiChecked[l]
                    << (m_hwDdphiReject ? " rejected " : " outside window ") << m_hwDdphiRejected[l];
        }
      }
      std::cout << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void SiHitSeedReco::accumulateResidualWidths()
{
  constexpr size_t kMaxPerLayer = 200000;
  for (const auto& c : m_eventData->chains)
  {
    if (!c.selected || !std::isfinite(c.tpc_chi2) || c.tpc_nres < static_cast<int>(c.layers.size()))
    {
      continue;
    }
    for (size_t i = 0; i < c.layers.size(); ++i)
    {
      const int l = c.layers[i];
      if (l < 0 || l > 6)
      {
        continue;
      }
      if (i < c.res_rdphi.size() && std::isfinite(c.res_rdphi[i]) && m_widthRphi[l].size() < kMaxPerLayer)
      {
        m_widthRphi[l].push_back(static_cast<float>(c.res_rdphi[i]));
      }
      if (i < c.res_dz_fit.size() && std::isfinite(c.res_dz_fit[i]) && m_widthZ[l].size() < kMaxPerLayer)
      {
        m_widthZ[l].push_back(static_cast<float>(c.res_dz_fit[i]));
      }
    }
  }
}

void SiHitSeedReco::printResidualWidths() const
{
  // robust width: 1.4826 * median |x - median|
  auto robust = [](std::vector<float> v, double& med, double& sig)
  {
    med = sig = std::numeric_limits<double>::quiet_NaN();
    if (v.size() < 5)
    {
      return;
    }
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    med = v[v.size() / 2];
    for (auto& x : v)
    {
      x = static_cast<float>(std::abs(x - med));
    }
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    sig = 1.4826 * v[v.size() / 2];
  };
  std::cout << Name() << ": residual widths of the selected chains to the TPC helix (robust sigma = 1.4826 MAD; z after the z model):" << std::endl;
  std::array<double, 7> sr{}, sz{};
  bool all = true;
  for (int l = 0; l < 7; ++l)
  {
    double mr = 0, mz = 0;
    robust(m_widthRphi[l], mr, sr[l]);
    robust(m_widthZ[l], mz, sz[l]);
    std::printf("%s:   layer %d  n %6zu  r*dphi median %+7.4f sigma %7.4f cm (in use %.4f) | z n %6zu median %+7.4f sigma %7.4f cm (in use %.4f)\n",
                Name().c_str(), l, m_widthRphi[l].size(), mr, sr[l], m_sigRphiL[l], m_widthZ[l].size(), mz, sz[l], m_sigZL[l]);
    all = all && std::isfinite(sr[l]) && std::isfinite(sz[l]);
  }
  std::cout << Name() << ":   to use them (the selection depends on the sigmas in use: iterate once or twice):" << std::endl;
  for (int l = 0; l < 7; ++l)
  {
    if (std::isfinite(sr[l]) && std::isfinite(sz[l]))
    {
      std::printf("      si_hit_seed->setTpcAgreementSigmaLayer(%d, %.3f, %.3f);\n", l, sr[l], sz[l]);
    }
    else
    {
      std::printf("      // layer %d: too few residuals\n", l);
    }
  }
  (void) all;
}

void SiHitSeedReco::resolveSharedChains()
{
  auto& chains = m_eventData->chains;
  auto& seeds = m_eventData->tpc_seeds;
  const bool dbgEvent = m_debugEvent >= 0 && static_cast<long>(m_eventData->event) == m_debugEvent;

  // same cuts as for the rank-0 selection in processTpcTrack
  auto passes = [&](const SiHitChain& c)
  {
    if (m_chainRanking == ChainRanking::AgreementFirst && !(std::isfinite(c.tpc_chi2) && c.tpc_chi2 <= m_agreementMaxChi2))
    {
      return false;
    }
    return m_maxTpcChi2 <= 0 || (std::isfinite(c.tpc_chi2) && c.tpc_chi2 <= m_maxTpcChi2);
  };
  auto chi2Key = [](const SiHitChain& c)
  {
    return std::isfinite(c.tpc_chi2) ? c.tpc_chi2 : std::numeric_limits<double>::max();
  };

  // chains of every TPC track (seed entry), in rank order
  std::vector<std::vector<int>> byTrack(seeds.size());
  std::unordered_map<int, size_t> seedOfTrack;
  for (size_t s = 0; s < seeds.size(); ++s)
  {
    seedOfTrack[seeds[s].tpc_track_index] = s;
  }
  for (size_t c = 0; c < chains.size(); ++c)
  {
    const auto it = seedOfTrack.find(chains[c].tpc_track_index);
    if (it != seedOfTrack.end())
    {
      byTrack[it->second].push_back(static_cast<int>(c));
    }
  }
  std::vector<size_t> order;
  for (size_t s = 0; s < seeds.size(); ++s)
  {
    auto& v = byTrack[s];
    std::sort(v.begin(), v.end(), [&](int a, int b) { return chains[a].rank < chains[b].rank; });
    if (seeds[s].selected_chain >= 0)
    {
      order.push_back(s);
    }
  }
  // best agreement first
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b)
  {
    return chi2Key(chains[seeds[a].selected_chain]) < chi2Key(chains[seeds[b].selected_chain]);
  });

  std::unordered_map<int, size_t> owner;  // blob -> seed entry that took it
  for (size_t s : order)
  {
    auto& info = seeds[s];
    const int oldSel = info.selected_chain;
    chains[oldSel].selected = false;
    int newSel = -1;
    size_t blocker = seeds.size();
    int blockShared = 0;
    for (int c : byTrack[s])
    {
      if (!passes(chains[c]) && c != oldSel)
      {
        continue;
      }
      std::unordered_map<size_t, int> sharedWith;
      for (int b : chains[c].blob_ids)
      {
        const auto it = owner.find(b);
        // crossing hypotheses of the same assembled track are alternatives, not competitors:
        // they may share silicon (TpcCrossingResolver picks one of them later)
        if (it != owner.end() && seeds[it->second].tpc_assembled_track_id != info.tpc_assembled_track_id)
        {
          ++sharedWith[it->second];
        }
      }
      size_t worst = seeds.size();
      int worstN = 0;
      for (const auto& [o, n] : sharedWith)
      {
        if (n > worstN)
        {
          worstN = n;
          worst = o;
        }
      }
      if (worstN >= m_exclusiveMinShared)
      {
        if (c == oldSel)
        {
          blocker = worst;
          blockShared = worstN;
        }
        continue;
      }
      newSel = c;
      break;
    }

    const bool dbg = dbgEvent && (m_debugTpcIndex < 0 || info.tpc_track_index == m_debugTpcIndex ||
                                  (blocker < seeds.size() && seeds[blocker].tpc_track_index == m_debugTpcIndex));
    if (newSel == oldSel)
    {
      ++m_exclKept;
    }
    else if (newSel >= 0)
    {
      ++m_exclReassigned;
    }
    else
    {
      ++m_exclLost;
    }
    if (dbg && newSel != oldSel)
    {
      std::printf("[SiSeedDebug] exclusive: TPC track %d (assembled %u) chain %d (rank 0, chi2/ndf %.3g, %d layers) shares %d blobs with the chain taken by TPC track %d (assembled %u, chi2/ndf %.3g)",
                  info.tpc_track_index, info.tpc_assembled_track_id, oldSel, chains[oldSel].tpc_chi2, chains[oldSel].n_layers(),
                  blockShared, blocker < seeds.size() ? seeds[blocker].tpc_track_index : -1,
                  blocker < seeds.size() ? seeds[blocker].tpc_assembled_track_id : 0u,
                  blocker < seeds.size() && seeds[blocker].selected_chain >= 0 ? chains[seeds[blocker].selected_chain].tpc_chi2 : NAN);
      if (newSel >= 0)
      {
        std::printf(" -> takes its chain %d (rank %d, chi2/ndf %.3g, %d layers)\n", newSel, chains[newSel].rank,
                    chains[newSel].tpc_chi2, chains[newSel].n_layers());
      }
      else
      {
        std::printf(" -> no other chain without conflict: Rejected\n");
      }
    }

    if (newSel >= 0)
    {
      chains[newSel].selected = true;
      info.selected_chain = newSel;
      for (int b : chains[newSel].blob_ids)
      {
        owner.emplace(b, s);
      }
    }
    else
    {
      info.selected_chain = -1;
      if (info.status == SiTpcSeedInfo::Chained)
      {
        --m_statusCount[SiTpcSeedInfo::Chained];
        ++m_statusCount[SiTpcSeedInfo::Rejected];
        info.status = SiTpcSeedInfo::Rejected;
      }
    }
  }

  if (!m_storeSelectedOnly)
  {
    return;
  }
  // prune the non-selected chains and re-index
  std::vector<int> newId(chains.size(), -1);
  std::vector<SiHitChain> kept;
  for (size_t c = 0; c < chains.size(); ++c)
  {
    if (chains[c].selected)
    {
      newId[c] = static_cast<int>(kept.size());
      kept.push_back(std::move(chains[c]));
      kept.back().id = newId[c];
    }
  }
  m_nChains -= chains.size() - kept.size();
  chains = std::move(kept);
  for (auto& info : seeds)
  {
    info.selected_chain = info.selected_chain >= 0 ? newId[info.selected_chain] : -1;
  }
}

bool SiHitSeedReco::debugThisTrack(const TpcTrackInput& in) const
{
  if (m_debugEvent < 0 || !m_eventData || static_cast<long>(m_eventData->event) != m_debugEvent)
  {
    return false;
  }
  if (m_debugTpcIndex >= 0 && in.index != m_debugTpcIndex)
  {
    return false;
  }
  if (m_debugAssembledId >= 0 && static_cast<long>(in.assembledId) != m_debugAssembledId)
  {
    return false;
  }
  return true;
}

void SiHitSeedReco::debugTrackHeader(const TpcTrackInput& in, const SiTpcReference::TpcReference* ref) const
{
  const double polyPt = std::hypot(in.polyPx, in.polyPy);
  std::printf("[SiSeedDebug] ==================================================\n");
  std::printf("[SiSeedDebug] event %llu  TPC_POLYTRACKS index %d  track id %u  assembled id %u  clusters %u\n",
              m_eventData->event, in.index, in.id, in.assembledId, in.nClusters);
  std::printf("[SiSeedDebug] Tpc_PolyTrack (its own fit, world frame): status %d  pt %.4f  phi %.5f  charge %+.0f  ref point (%.3f, %.3f, %.3f)\n",
              in.polyFitStatus, polyPt, std::atan2(in.polyPy, in.polyPx), in.polyCharge, in.polyX, in.polyY, in.polyZ);
  const auto& f = in.fit;
  std::printf("[SiSeedDebug] SiTpcHelixFit refit used here (raw cluster centroids, TPC frame): status %d  pt %.4f  phi %.5f  charge %+d  helicity %+d\n",
              f.status, f.pt, f.phi, f.charge, f.helicity);
  std::printf("[SiSeedDebug]   circle (%.3f, %.3f) R %.2f  pca (%.4f, %.4f)  dca %.4f  z0 %.3f  tanl %.4f  circle rms %.4f  z rms %.4f\n",
              f.cx, f.cy, f.R, f.pcaX, f.pcaY, f.dca, f.z0, f.tanl, f.circleRms, f.zRms);
  if (!ref)
  {
    std::printf("[SiSeedDebug]   no TPC reference (fit not usable for prediction)\n");
    return;
  }
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  std::printf("[SiSeedDebug] predictions (detector frame; Uphi bin %.4f rad, z bin %.4f cm this event):\n", dp, dz);
  for (int l = 0; l < nlayers; ++l)
  {
    for (int alt = 0; alt < 2; ++alt)
    {
      const auto& p = alt ? ref->alt[l] : ref->layers[l];
      if (!p.valid)
      {
        if (!alt)
        {
          std::printf("[SiSeedDebug]   layer %d: no prediction\n", l);
        }
        continue;
      }
      std::printf("[SiSeedDebug]   layer %d%s: half %d  Uphi %8.5f  z %8.3f  s %7.3f  x %7.3f y %7.3f  cosA %.3f\n",
                  l, alt ? " (other half)" : "", p.half, p.u, p.z, p.s, p.x, p.y, p.cosA);
    }
  }
}

void SiHitSeedReco::debugAnchors(const SiTpcReference::TpcReference& ref, int outer, const std::vector<int>& anchors) const
{
  const double R = m_radius[outer];
  std::printf("[SiSeedDebug]   blobs on layer %d within 3x anchor window (|r dphi| < %.2f, |dz| < %.2f cm); %zu anchors used (max %u):\n",
              outer, 3 * m_anchorRphi, 3 * m_anchorZ, anchors.size(), m_maxAnchors);
  int nPrinted = 0;
  for (const auto& b : m_eventData->blobs)
  {
    if (b.layer != outer)
    {
      continue;
    }
    double bestRphi = NAN, bestDz = NAN;
    for (const auto* p : {&ref.layers[outer], &ref.alt[outer]})
    {
      if (!p->valid)
      {
        continue;
      }
      const double rphi = R * wrapDelta(b.phi - p->u, twopi);
      const double dzc = b.z - p->z;
      const double zs = m_anchorZ > 0 ? m_anchorZ : 10.0;
      if (!std::isfinite(bestRphi) || std::hypot(rphi / m_anchorRphi, dzc / zs) < std::hypot(bestRphi / m_anchorRphi, bestDz / zs))
      {
        bestRphi = rphi;
        bestDz = dzc;
      }
    }
    const bool zCut = m_anchorZ > 0;
    if (!std::isfinite(bestRphi) || std::abs(bestRphi) > 3 * m_anchorRphi || (zCut && std::abs(bestDz) > 3 * m_anchorZ))
    {
      continue;
    }
    const auto it = std::find(anchors.begin(), anchors.end(), b.id);
    const bool inWin = std::abs(bestRphi) <= m_anchorRphi && (!zCut || std::abs(bestDz) <= m_anchorZ);
    std::printf("[SiSeedDebug]     blob %5d  Uphi %8.5f Vz %8.3f  nhits %2zu  r*dphi %+7.3f cm  dz %+7.3f cm  %s%s\n",
                b.id, b.phi, b.z, b.hit_ids.size(), bestRphi, bestDz, inWin ? "IN window" : "outside",
                it != anchors.end() ? (" -> anchor #" + std::to_string(it - anchors.begin())).c_str() : (inWin ? " (dropped: max anchors)" : ""));
    ++nPrinted;
  }
  if (nPrinted == 0)
  {
    std::printf("[SiSeedDebug]     (none)\n");
  }
}

void SiHitSeedReco::debugSeed(const SiBlob& seed, int inner, const SiTpcReference::TpcReference& ref,
                              const std::vector<Candidate>& cands) const
{
  double zc = 0, zh = 0, pc = 0;
  seedWindow(ref, seed.layer, inner, seed.z, zc, zh, pc);
  std::printf("[SiSeedDebug]   anchor blob %d (layer %d, Uphi %.5f, Vz %.3f): seed window on layer %d  dUphi centre %+.3f bins half %.2f, dz centre %+.3f bins half %.2f; %zu candidates (max %u used)\n",
              seed.id, seed.layer, seed.phi, seed.z, inner, pc, m_seedPhiHalf, zc, zh, cands.size(), m_maxSeedCand);
  for (size_t k = 0; k < cands.size(); ++k)
  {
    const auto& c = cands[k];
    std::printf("[SiSeedDebug]     cand #%zu blob %5d  res dUphi %+6.3f bins  res dz %+6.3f bins  score %.3f%s\n",
                k, c.blob, c.res_phi, c.res_z, c.score, k < m_maxSeedCand ? "" : "  (not used: max candidates)");
  }
  // inner-layer blobs near the window (blob centroid, relative to the anchor centroid)
  for (const auto& b : m_eventData->blobs)
  {
    if (b.layer != inner)
    {
      continue;
    }
    const double dph = wrapDelta(b.phibin - (seed.phibin + pc), m_nphi);
    const double dzb = b.zbin - (seed.zbin + zc);
    if (std::abs(dph) > 3 * m_seedPhiHalf || std::abs(dzb) > 3 * zh)
    {
      continue;
    }
    const bool isCand = std::any_of(cands.begin(), cands.end(), [&](const Candidate& c) { return c.blob == b.id; });
    std::printf("[SiSeedDebug]       near: blob %5d  Uphi %8.5f Vz %8.3f  dUphi-centre %+6.3f bins  dz-centre %+6.3f bins  %s\n",
                b.id, b.phi, b.z, dph, dzb, isCand ? "candidate" : "NOT a candidate (outside window, hit level)");
  }
}

void SiHitSeedReco::debugChain(const SiHitChain& ch, const char* tag) const
{
  std::cout << tag << ": layers " << ch.layers.size() << " (MVTX " << ch.n_mvtx << ", INTT " << ch.n_intt << "), blobs";
  for (size_t i = 0; i < ch.layers.size(); ++i)
  {
    std::cout << " L" << ch.layers[i] << ":" << ch.blob_ids[i];
  }
  std::cout << std::endl;
  for (const auto& st : ch.steps)
  {
    std::printf("[SiSeedDebug]       step kind %d  %d->%d  blob %5d  res Uphi %+6.3f bins (half %.2f)  res z %+6.3f bins (half %.2f)\n",
                st.kind, st.from_layer, st.to_layer, st.blob_id, st.res_phi_bins, st.phi_half_window_bins,
                st.res_z_bins, st.z_half_window_bins);
  }
}

bool SiHitSeedReco::hwDdphiTerms(int layer, const std::map<int, double>& uphiByLayer, double& dphi, double& dphiPred) const
{
  const int p = layer + 1;
  const auto iL = uphiByLayer.find(layer);
  const auto iP = uphiByLayer.find(p);
  if (iL == uphiByLayer.end() || iP == uphiByLayer.end())
  {
    return false;
  }
  const auto iPP = uphiByLayer.upper_bound(p);  // next outer layer of the chain
  if (iPP == uphiByLayer.end())
  {
    return false;
  }
  const double drPrev = m_radius[p] - m_radius[iPP->first];
  if (std::abs(drPrev) < 1e-3)
  {
    return false;
  }
  dphi = wrapDelta(iL->second - iP->second, twopi);
  dphiPred = wrapDelta(iP->second - iPP->second, twopi) * (m_radius[layer] - m_radius[p]) / drPrev;
  return true;
}

double SiHitSeedReco::wrapDelta(double x, double period)
{
  x = std::fmod(x + period / 2., period);
  if (x < 0)
  {
    x += period;
  }
  return x - period / 2.;
}

double SiHitSeedReco::circularMean(const std::vector<double>& values)
{
  double sx = 0, sy = 0;
  for (double p : values)
  {
    sx += std::cos(p);
    sy += std::sin(p);
  }
  double p = std::atan2(sy, sx);
  if (p < 0)
  {
    p += twopi;
  }
  return p;
}

std::pair<double, double> SiHitSeedReco::localCoordinates(int layer, unsigned int row, unsigned int col, int ladderz) const
{
  if (layer <= 2)
  {
    constexpr double row_pitch = 26.88e-4, col_pitch = 29.24e-4, sensor_x = 0.058128;
    constexpr double first_row = 0.5 * (512 * row_pitch - 37.44e-4 + 0.12 - row_pitch);
    return {first_row - row * row_pitch - sensor_x, (col + 0.5 - 512.) * col_pitch};
  }
  const int type = ladderz % 2;
  const int ncols = type == 0 ? 8 : 5;
  const double pitchz = type == 0 ? 1.6 : 2.0;
  return {(row + 0.5 - 128.) * 0.0078, (col + 0.5 - ncols / 2.) * pitchz};
}

double SiHitSeedReco::intrinsicPhi(int layer, int stave, int ladderphi, double lx) const
{
  double p = 0;
  if (layer <= 2)
  {
    p = mvtx_phi0[layer] - mvtx_phi0[2] + stave * twopi / nstaves[layer] + lx * std::cos(mvtx_tilt[layer]) / mvtx_r[layer];
  }
  else
  {
    const int i = layer - 3;
    p = intt_phi0[i] + twopi * ladderphi / intt_nladder[i] + lx / intt_r[i] - mvtx_phi0[2];
  }
  p = std::fmod(p, twopi);
  if (p < 0)
  {
    p += twopi;
  }
  return p;
}

double SiHitSeedReco::longitudinalZ(int layer, int chip, int ladderz, double ly) const
{
  return (layer <= 2 ? mvtx_chip_z.at(chip) : intt_sensor_z.at(ladderz)) + ly;
}

bool SiHitSeedReco::fillHits(PHCompositeNode* topNode)
{
  auto* cont = findNode::getClass<TrkrHitSetContainer>(topNode, m_hitNodeName);
  if (!cont)
  {
    std::cerr << Name() << ": missing " << m_hitNodeName << std::endl;
    return false;
  }

  for (auto det : {TrkrDefs::mvtxId, TrkrDefs::inttId})
  {
    auto range = cont->getHitSets(det);
    for (auto hsit = range.first; hsit != range.second; ++hsit)
    {
      const auto hskey = hsit->first;
      auto* hs = hsit->second;
      const int layer = TrkrDefs::getLayer(hskey);
      if (layer < 0 || layer > 6)
      {
        continue;
      }
      const int stave = layer <= 2 ? MvtxDefs::getStaveId(hskey) : -1;
      const int chip = layer <= 2 ? MvtxDefs::getChipId(hskey) : -1;
      const int ladderz = layer >= 3 ? InttDefs::getLadderZId(hskey) : -1;
      const int ladderphi = layer >= 3 ? InttDefs::getLadderPhiId(hskey) : -1;
      auto hits = hs->getHits();
      for (auto it = hits.first; it != hits.second; ++it)
      {
        const auto hk = it->first;
        const unsigned int row = layer <= 2 ? MvtxDefs::getRow(hk) : InttDefs::getRow(hk);
        const unsigned int col = layer <= 2 ? MvtxDefs::getCol(hk) : InttDefs::getCol(hk);
        const auto [lx, ly] = localCoordinates(layer, row, col, ladderz);
        SiHitPoint h;
        h.layer = layer;
        h.hitsetkey = hskey;
        h.hitkey = hk;
        h.row = row;
        h.col = col;
        h.stave = stave;
        h.chip = chip;
        h.ladderphi = ladderphi;
        h.ladderz = ladderz;
        h.adc = it->second ? it->second->getAdc() : 0;
        h.lx = lx;
        h.ly = ly;
        // hardware coordinates for blobs and chaining: analytic model, replaced below by the official
        // position of the hit (setHardwareFromActs, with setUseActsGeometry)
        h.z = longitudinalZ(layer, chip, ladderz, ly);
        h.phi = intrinsicPhi(layer, stave, ladderphi, lx);
        alignHardware(hskey, layer, h.z, h.phi);
        bool hwFromActs = false;
        // physics position on the sensor for the track fits
        Eigen::Vector3d pos = alignPoint(hskey, sensorPoint(layer, stave, chip, ladderphi, ladderz, lx, ly));
        Eigen::Vector3d g;
        if ((m_useActsGeometry || m_compareActs) && actsGlobal(hskey, layer, row, col, ladderz, g))
        {
          if (m_compareActs)
          {
            const auto own = m_frame.toGlobal(layer, {pos.x(), pos.y(), pos.z()});
            const double ro = std::hypot(own.x, own.y), ra = std::hypot(g.x(), g.y());
            double dphi = std::atan2(own.y, own.x) - std::atan2(g.y(), g.x());
            dphi = std::remainder(dphi, twopi);
            const std::array<double, 5> d{own.x - g.x(), own.y - g.y(), own.z - g.z(), ro - ra, ra * dphi};
            auto& c = m_actsCmp[layer];
            c[0] += 1;
            for (int k = 0; k < 5; ++k)
            {
              c[1 + k] += d[k];
              c[6 + k] += d[k] * d[k];
            }
          }
          if (m_useActsGeometry)
          {
            const auto dg = m_frame.toDetectorFromGlobal(layer, {g.x(), g.y(), g.z() + m_actsZShift});
            pos = Eigen::Vector3d(dg.x, dg.y, dg.z);
            if (m_hwFromActs && std::hypot(dg.x, dg.y) > 0)
            {
              // hardware coordinates of the same hit from its official position (same detector frame
              // as SiDetectorFrame::toDetector: angle about the layer's centre minus the Uphi rotation)
              double u = std::fmod(std::atan2(dg.y, dg.x) - m_frame.uphiRotation(), twopi);
              if (u < 0)
              {
                u += twopi;
              }
              const double du = std::remainder(h.phi - u, twopi);
              const double dz = h.z - dg.z;
              auto& c = m_hwCmp[layer];
              c[0] += 1;
              c[1] += du;
              c[2] += du * du;
              c[3] += dz;
              c[4] += dz * dz;
              c[5] += std::abs(du) > 0.05 ? 1 : 0;
              h.phi = u;
              h.z = dg.z;
              hwFromActs = true;
            }
          }
        }
        ++(hwFromActs ? m_hwActs : m_hwAnalytic);
        if (layer >= 3)
        {
          pos.x() += m_inttPhysOffset[0];
          pos.y() += m_inttPhysOffset[1];
        }
        h.px = pos.x();
        h.py = pos.y();
        h.pz = pos.z() + m_layerZOffset[layer];
        m_eventData->hits.push_back(h);
      }
    }
  }
  return true;
}

void SiHitSeedReco::buildBlobs()
{
  if (m_blobsFromPixels)
  {
    buildBlobsFromPixels();
    return;
  }

  using Cell = std::pair<int, int>;
  // per layer and INTT time bucket (MVTX: one group): hits of different bunch crossings are
  // different particles and never end up in the same blob
  std::array<std::map<int, std::map<Cell, std::vector<int>>>, nlayers> layerCells;
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  for (size_t i = 0; i < m_eventData->hits.size(); ++i)
  {
    auto& h = m_eventData->hits[i];
    const int iz = (int) std::floor((h.z - m_zmin) / dz);
    const int ip = (int) std::floor(h.phi / dp) % m_nphi;
    if (iz < 0 || iz >= m_nz)
    {
      continue;
    }
    const int bucket = (h.layer >= 3 && m_inttBlobsPerBucket) ? static_cast<int>(InttDefs::getTimeBucketId(h.hitsetkey)) : 0;
    layerCells[h.layer][bucket][{iz, ip}].push_back(i);
  }

  for (int layer = 6; layer >= 0; --layer)
  {
   for (auto& [bucket, remaining] : layerCells[layer])
   {
    while (!remaining.empty())
    {
      const Cell start = remaining.begin()->first;
      std::queue<Cell> q;
      q.push(start);
      std::vector<Cell> comp;
      std::vector<int> ids;
      while (!q.empty())
      {
        auto c = q.front();
        q.pop();
        auto it = remaining.find(c);
        if (it == remaining.end())
        {
          continue;
        }
        for (int id : it->second)
        {
          ids.push_back(id);
        }
        comp.push_back(c);
        remaining.erase(it);
        for (int a = -1; a <= 1; ++a)
        {
          for (int b = -1; b <= 1; ++b)
          {
            int ip = (c.second + b) % m_nphi;
            if (ip < 0)
            {
              ip += m_nphi;
            }
            const Cell n{c.first + a, ip};
            if (remaining.count(n))
            {
              q.push(n);
            }
          }
        }
      }
      std::sort(ids.begin(), ids.end());
      SiBlob blob;
      blob.id = m_eventData->blobs.size();
      blob.layer = layer;
      blob.hit_ids = ids;
      blob.cells = comp;
      std::vector<double> phis;
      double zsum = 0;
      for (int id : ids)
      {
        zsum += m_eventData->hits[id].z;
        phis.push_back(m_eventData->hits[id].phi);
        blob.adc_sum += m_eventData->hits[id].adc;
      }
      blob.z = zsum / ids.size();
      blob.phi = circularMean(phis);
      blob.zbin = (blob.z - m_zmin) / dz;
      blob.phibin = blob.phi / dp;
      m_eventData->blobs.push_back(std::move(blob));
    }
   }
  }
}

// Blobs = clusters of adjacent readout channels, as the hardware clusterizers do:
//   MVTX: same chip and strobe (hitsetkey), pixels touching in (row, col) incl. diagonals
//   INTT: same ladder, sensor and time bucket (hitsetkey), same z column, neighbouring strips
// so two particles only share a blob if their pixels / strips touch.
void SiHitSeedReco::buildBlobsFromPixels()
{
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  std::map<TrkrDefs::hitsetkey, std::vector<int>> bySet;
  for (size_t i = 0; i < m_eventData->hits.size(); ++i)
  {
    const auto& h = m_eventData->hits[i];
    const int iz = (int) std::floor((h.z - m_zmin) / dz);
    if (iz < 0 || iz >= m_nz)
    {
      continue;  // outside the z binning: not used by the finder
    }
    bySet[h.hitsetkey].push_back(static_cast<int>(i));
  }
  // outer layers first, as the grid blobs
  std::vector<std::pair<int, TrkrDefs::hitsetkey>> sets;
  for (const auto& [key, ids] : bySet)
  {
    sets.emplace_back(m_eventData->hits[ids.front()].layer, key);
  }
  std::stable_sort(sets.begin(), sets.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

  for (const auto& [layer, key] : sets)
  {
    const bool mvtx = layer <= 2;
    std::map<std::pair<int, int>, std::vector<int>> remaining;  // (row, col) -> hits
    for (const int id : bySet[key])
    {
      const auto& h = m_eventData->hits[id];
      remaining[{static_cast<int>(h.row), static_cast<int>(h.col)}].push_back(id);
    }
    while (!remaining.empty())
    {
      std::queue<std::pair<int, int>> q;
      q.push(remaining.begin()->first);
      std::vector<int> ids;
      while (!q.empty())
      {
        const auto px = q.front();
        q.pop();
        auto it = remaining.find(px);
        if (it == remaining.end())
        {
          continue;
        }
        ids.insert(ids.end(), it->second.begin(), it->second.end());
        remaining.erase(it);
        for (int a = -1; a <= 1; ++a)
        {
          for (int b = (mvtx ? -1 : 0); b <= (mvtx ? 1 : 0); ++b)  // INTT: same z column
          {
            const std::pair<int, int> n{px.first + a, px.second + b};
            if (remaining.count(n))
            {
              q.push(n);
            }
          }
        }
      }
      std::sort(ids.begin(), ids.end());
      SiBlob blob;
      blob.id = m_eventData->blobs.size();
      blob.layer = layer;
      blob.hit_ids = ids;
      std::vector<double> phis;
      double zsum = 0;
      std::set<std::pair<int, int>> cells;
      for (const int id : ids)
      {
        const auto& h = m_eventData->hits[id];
        zsum += h.z;
        phis.push_back(h.phi);
        blob.adc_sum += h.adc;
        cells.insert({(int) std::floor((h.z - m_zmin) / dz), ((int) std::floor(h.phi / dp) % m_nphi + m_nphi) % m_nphi});
      }
      blob.cells.assign(cells.begin(), cells.end());
      blob.z = zsum / ids.size();
      blob.phi = circularMean(phis);
      blob.zbin = (blob.z - m_zmin) / dz;
      blob.phibin = blob.phi / dp;
      m_eventData->blobs.push_back(std::move(blob));
    }
  }
}

// Hit -> blob map and a (layer, iz, iphi) cell index, so that every windowed search
// only visits the few cells inside its window instead of all hits of the layer.
void SiHitSeedReco::buildHitGrid(const std::vector<char>* excludedBlobs)
{
  const auto& hits = m_eventData->hits;
  m_hitBlob.assign(hits.size(), -1);
  for (const auto& b : m_eventData->blobs)
  {
    for (int id : b.hit_ids)
    {
      m_hitBlob[id] = b.id;
    }
  }

  const size_t ncell = static_cast<size_t>(nlayers) * m_nz * m_nphi;
  m_cellStart.assign(ncell + 1, 0);
  std::vector<int> cellOf(hits.size(), -1);
  for (size_t i = 0; i < hits.size(); ++i)
  {
    if (m_hitBlob[i] < 0)
    {
      continue;  // outside the z binning: not used by the finder
    }
    if (excludedBlobs && static_cast<std::size_t>(m_hitBlob[i]) < excludedBlobs->size() && (*excludedBlobs)[m_hitBlob[i]])
    {
      continue;  // blob not available in this search (worker interface: used by a kept chain)
    }
    const int iz = (int) std::floor(hits[i].zbin);
    const int ip = ((int) std::floor(hits[i].phibin) % m_nphi + m_nphi) % m_nphi;
    if (iz < 0 || iz >= m_nz)
    {
      continue;
    }
    cellOf[i] = (hits[i].layer * m_nz + iz) * m_nphi + ip;
    ++m_cellStart[cellOf[i] + 1];
  }
  for (size_t c = 0; c < ncell; ++c)
  {
    m_cellStart[c + 1] += m_cellStart[c];
  }
  m_cellHits.assign(m_cellStart.back(), -1);
  std::vector<int> fill(m_cellStart.begin(), m_cellStart.end() - 1);
  for (size_t i = 0; i < hits.size(); ++i)
  {
    if (cellOf[i] >= 0)
    {
      m_cellHits[fill[cellOf[i]]++] = static_cast<int>(i);
    }
  }
}

template <class F>
void SiHitSeedReco::forEachHitNear(int layer, double zbinCenter, double zHalf,
                                   double phibinCenter, double phiHalf, F&& f) const
{
  const int izlo = std::max(0, (int) std::floor(zbinCenter - zHalf));
  const int izhi = std::min(m_nz - 1, (int) std::floor(zbinCenter + zHalf));
  if (izlo > izhi)
  {
    return;
  }
  int iplo = (int) std::floor(phibinCenter - phiHalf);
  int iphi = (int) std::floor(phibinCenter + phiHalf);
  if (iphi - iplo + 1 >= m_nphi)
  {
    iplo = 0;
    iphi = m_nphi - 1;  // window covers the whole circle: visit each cell once
  }
  for (int iz = izlo; iz <= izhi; ++iz)
  {
    for (int ipr = iplo; ipr <= iphi; ++ipr)
    {
      const int ip = (ipr % m_nphi + m_nphi) % m_nphi;
      const int c = (layer * m_nz + iz) * m_nphi + ip;
      for (int k = m_cellStart[c]; k < m_cellStart[c + 1]; ++k)
      {
        f(m_cellHits[k]);
      }
    }
  }
}

// ============================================================================================
// TPC tracks
// ============================================================================================
std::vector<SiHitSeedReco::TpcTrackInput> SiHitSeedReco::collectTpcTracks(PHCompositeNode* topNode) const
{
  std::vector<TpcTrackInput> out;
  if (!m_tpcTrajectoryNodeName.empty())
  {
    // crossing trajectories: three points each, one trajectory per assembled track and crossing
    auto* trajectories = findNode::getClass<Tpc_CrossingTrajectoryContainer>(topNode, m_tpcTrajectoryNodeName);
    if (!trajectories)
    {
      if (Verbosity() > 0)
      {
        std::cout << Name() << ": missing " << m_tpcTrajectoryNodeName << " - no TPC-driven chaining in this event" << std::endl;
      }
      return out;
    }
    SiTpcHelixFit::Config cfg;
    cfg.beamX = m_beamAlign.referenceBeamX();
    cfg.beamY = m_beamAlign.referenceBeamY();
    cfg.bz = m_bz;
    for (unsigned int i = 0; i < trajectories->size(); ++i)
    {
      const Tpc_CrossingTrajectory* t = trajectories->get_trajectory(i);
      if (!t)
      {
        continue;
      }
      TpcTrackInput in;
      in.index = static_cast<int>(i);
      in.id = t->get_trajectory_id();
      in.assembledId = t->get_source_assembled_track_id();
      in.crossing = t->get_crossing();
      in.polyX = t->get_x();
      in.polyY = t->get_y();
      in.polyZ = t->get_z();
      in.polyPx = t->get_px();
      in.polyPy = t->get_py();
      in.polyPz = t->get_pz();
      in.polyCharge = t->get_charge();
      in.polyFitStatus = t->get_fit_status();
      std::vector<SiTpcHelixFit::Point> pts;
      for (unsigned int k = 0; k < Tpc_CrossingTrajectory::NPoints; ++k)
      {
        const auto p = m_beamAlign.tpcToTpcFrame(t->get_point_x(k), t->get_point_y(k), t->get_point_z(k));
        if (std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z))
        {
          pts.push_back({p.x, p.y, p.z + m_tpcZOffset, 1.0, 1.0, 1});
        }
      }
      in.nClusters = static_cast<unsigned int>(pts.size());
      if (pts.size() >= m_minTpcClusters)
      {
        std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b)
                  { return a.x * a.x + a.y * a.y < b.x * b.x + b.y * b.y; });
        in.fit = SiTpcHelixFit::fit(pts, cfg);
      }
      out.push_back(std::move(in));
    }
    return out;
  }

  auto* clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  auto* tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName);
  if (!clusters || !tracks)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_tpcClusterNodeName << " or " << m_tpcTrackNodeName
                << " - no TPC-driven chaining in this event" << std::endl;
    }
    return out;
  }
  // TPC = reference frame: clusters as measured (SiTpcBeamAlignment::tpcToTpcFrame is the identity)
  std::map<unsigned int, std::vector<SiTpcHelixFit::Point>> byTrack;
  for (unsigned int i = 0; i < clusters->size(); ++i)
  {
    const Tpc_PolyCluster* c = clusters->get_cluster(i);
    if (!c || !c->isValid())
    {
      continue;
    }
    const auto p = m_beamAlign.tpcToTpcFrame(c->get_centroid_x(), c->get_centroid_y(), c->get_centroid_z());
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
    {
      continue;
    }
    byTrack[c->get_source_assembled_track_id()].push_back({p.x, p.y, p.z + m_tpcZOffset, 1.0, 1.0, 1});
  }
  SiTpcHelixFit::Config cfg;
  cfg.beamX = m_beamAlign.referenceBeamX();
  cfg.beamY = m_beamAlign.referenceBeamY();
  cfg.bz = m_bz;
  for (unsigned int i = 0; i < tracks->size(); ++i)
  {
    const Tpc_PolyTrack* t = tracks->get_track(i);
    if (!t)
    {
      continue;
    }
    TpcTrackInput in;
    in.index = static_cast<int>(i);
    in.id = t->get_track_id();
    in.assembledId = t->get_source_assembled_track_id();
    in.polyX = t->get_x();
    in.polyY = t->get_y();
    in.polyZ = t->get_z();
    in.polyPx = t->get_px();
    in.polyPy = t->get_py();
    in.polyPz = t->get_pz();
    in.polyCharge = t->get_charge();
    in.polyFitStatus = t->get_fit_status();
    auto it = byTrack.find(in.assembledId);
    if (it != byTrack.end())
    {
      auto& pts = it->second;
      in.nClusters = static_cast<unsigned int>(pts.size());
      if (pts.size() >= m_minTpcClusters)
      {
        std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b)
                  { return a.x * a.x + a.y * a.y < b.x * b.x + b.y * b.y; });
        in.fit = SiTpcHelixFit::fit(pts, cfg);
      }
    }
    out.push_back(std::move(in));
  }
  return out;
}

// Blobs of `layer` with at least one hit in the anchor window around the TPC prediction
// (both clamshell translations near a boundary), closest first.
std::vector<int> SiHitSeedReco::anchorBlobs(const SiTpcReference::TpcReference& ref, int layer) const
{
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  std::unordered_map<int, double> metric;
  for (const auto* p : {&ref.layers[layer], &ref.alt[layer]})
  {
    if (!p->valid)
    {
      continue;
    }
    const double R = m_radius[layer];
    const bool zCut = m_anchorZ > 0;
    const double zCenterBins = zCut ? (p->z - m_zmin) / dz : 0.5 * m_nz;
    const double zHalfBins = zCut ? m_anchorZ / dz : m_nz;  // no z cut: whole z range
    const double zScale = zCut ? m_anchorZ : 10.0;          // ordering scale if no z window [cm]
    forEachHitNear(layer, zCenterBins, zHalfBins, p->u / dp, (m_anchorRphi / R) / dp, [&](int ih)
    {
      const int b = m_hitBlob[ih];
      if (b < 0)
      {
        return;
      }
      const auto& h = m_eventData->hits[ih];
      const double rdphi = R * SiTpcReference::wrapPi(h.phi - p->u);
      const double dzc = h.z - p->z;
      if (std::abs(rdphi) > m_anchorRphi || (zCut && std::abs(dzc) > m_anchorZ))
      {
        return;
      }
      const double m = (rdphi / m_anchorRphi) * (rdphi / m_anchorRphi) + m_anchorZWeight * (dzc / zScale) * (dzc / zScale);
      auto it = metric.find(b);
      if (it == metric.end() || m < it->second)
      {
        metric[b] = m;
      }
    });
  }
  std::vector<std::pair<double, int>> order;
  order.reserve(metric.size());
  for (const auto& [b, m] : metric)
  {
    order.emplace_back(m, b);
  }
  std::sort(order.begin(), order.end());
  std::vector<int> out;
  for (const auto& [m, b] : order)
  {
    if (out.size() >= m_maxAnchors)
    {
      break;
    }
    out.push_back(b);
  }
  return out;
}

// Seed window outer hit (layer `outer`, z = zOuter) -> layer `inner`, in bins of this event:
//   z  : centre zCenterBins relative to the outer hit, half width zHalfBins
//   phi: centre phiCenterBins relative to the outer hit (TPC curvature)
void SiHitSeedReco::seedWindow(const SiTpcReference::TpcReference& ref, int outer, int inner, double /*zOuter*/,
                               double& zCenterBins, double& zHalfBins, double& phiCenterBins) const
{
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  const auto& pOut = ref.layers[outer];
  const auto& pIn = ref.layers[inner];
  // TPC slope only: dz = tanl * ds, ds = transverse path between the two layers along the TPC
  // circle (< 0 going inward).  The TPC pca z is NOT used, so a TPC z offset / t0 does not bias it.
  const double ds = (pOut.valid && pIn.valid) ? (pIn.s - pOut.s) : (m_radius[inner] - m_radius[outer]) / std::max(1e-3, pOut.cosA);
  const double dzcm = ref.tanl * ds;
  zCenterBins = (dzcm + m_seedZOffsetCm) / dz;
  zHalfBins = m_seedZHalfCm / dz;
  phiCenterBins = (m_seedPhiFromTpc && pOut.valid && pIn.valid) ? SiTpcReference::wrapPi(pIn.u - pOut.u) / dp : 0.0;
}

std::vector<SiHitSeedReco::Candidate> SiHitSeedReco::seedCandidates(const SiBlob& seed, int inner,
                                                                   const SiTpcReference::TpcReference& ref) const
{
  std::unordered_map<int, Candidate> best;
  auto better = [](const Candidate& a, const Candidate& b)
  {
    return std::make_tuple(a.score, a.blob) < std::make_tuple(b.score, b.blob);
  };
  for (int oh : seed.hit_ids)
  {
    const auto& o = m_eventData->hits[oh];
    double zc = 0, zh = 0, pc = 0;
    seedWindow(ref, seed.layer, inner, o.z, zc, zh, pc);
    forEachHitNear(inner, o.zbin + zc, zh, o.phibin + pc, m_seedPhiHalf, [&](int ih)
    {
      const int b = m_hitBlob[ih];
      if (b < 0)
      {
        return;
      }
      const auto& h = m_eventData->hits[ih];
      const double dphi = wrapDelta(h.phibin - (o.phibin + pc), m_nphi);
      const double dzres = h.zbin - (o.zbin + zc);
      if (std::abs(dphi) > m_seedPhiHalf || std::abs(dzres) > zh)
      {
        return;
      }
      Candidate c;
      c.blob = b;
      c.hit = ih;
      c.layer = inner;
      c.outer_hit = oh;
      c.z = h.z;
      c.phi = h.phi;
      c.outer_z = o.z;
      c.outer_phi = o.phi;
      c.res_z = dzres;
      c.res_phi = dphi;
      c.z_center = zc;
      c.z_half = zh;
      c.phi_center = pc;
      c.score = (dphi / m_seedPhiHalf) * (dphi / m_seedPhiHalf) + (dzres / zh) * (dzres / zh);
      auto jt = best.find(b);
      if (jt == best.end() || better(c, jt->second))
      {
        best[b] = c;
      }
    });
  }
  std::vector<Candidate> out;
  out.reserve(best.size());
  for (auto& [_, c] : best)
  {
    out.push_back(c);
  }
  std::sort(out.begin(), out.end(), better);
  return out;
}

Eigen::Vector3d SiHitSeedReco::xyz(int layer, double z, double phi) const
{
  return {m_radius[layer] * std::cos(phi), m_radius[layer] * std::sin(phi), z};
}

SiHitSeedReco::Line3D SiHitSeedReco::fitLine(const std::vector<Eigen::Vector3d>& p) const
{
  Line3D l;
  if (p.size() == 2)
  {
    // two points: the line through them (same as the principal axis, without the eigen solver)
    l.c = 0.5 * (p[0] + p[1]);
    const Eigen::Vector3d d = p[1] - p[0];
    const double n = d.norm();
    l.d = n > 0 ? Eigen::Vector3d(d / n) : Eigen::Vector3d(0, 0, 1);
    return l;
  }
  l.c = Eigen::Vector3d::Zero();
  for (const auto& x : p)
  {
    l.c += x;
  }
  l.c /= p.size();
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto& x : p)
  {
    auto d = x - l.c;
    cov += d * d.transpose();
  }
  // closed-form 3x3 eigen decomposition (much faster than the iterative solver)
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(cov);
  l.d = es.eigenvectors().col(2).normalized();
  return l;
}

bool SiHitSeedReco::propagate(const Line3D& l, int layer, const Eigen::Vector3d& near, double& z, double& phi) const
{
  const double r = m_radius[layer], a = l.d.x() * l.d.x() + l.d.y() * l.d.y();
  if (a < 1e-12)
  {
    return false;
  }
  const double b = 2 * (l.c.x() * l.d.x() + l.c.y() * l.d.y());
  const double q = l.c.x() * l.c.x() + l.c.y() * l.c.y() - r * r;
  const double disc = b * b - 4 * a * q;
  if (disc < 0)
  {
    return false;
  }
  const double s = std::sqrt(disc);
  const Eigen::Vector3d p1 = l.c + ((-b - s) / (2 * a)) * l.d;
  const Eigen::Vector3d p2 = l.c + ((-b + s) / (2 * a)) * l.d;
  const auto& p = ((p1 - near).squaredNorm() < (p2 - near).squaredNorm()) ? p1 : p2;
  z = p.z();
  phi = std::atan2(p.y(), p.x());
  if (phi < 0)
  {
    phi += twopi;
  }
  return true;
}

// Seed pair -> chain, as before: straight line through the MVTX points, propagated to the missing
// MVTX layer and to the INTT layers.  Everything here is in hardware coordinates: a hit is the
// point (R_layer cos Uphi, R_layer sin Uphi, Vz) of the nominal cylinder in the detector frame
// (rotated by the common Uphi rotation, which a straight line does not care about), the same
// frame for all hits, so no conversion is needed.  The chain goes to physics space only
// afterwards (toPhysics).  No global hit mask.
bool SiHitSeedReco::buildChain(const SiBlob& seed, const Candidate& cand, SiHitChain& out,
                               const SiTpcReference::TpcReference* tpcRef) const
{
  // blobs taken by this chain (at most one per layer)
  std::array<int, nlayers> chainBlobs{};
  int nChainBlobs = 0;
  auto isUsed = [&](int b)
  {
    return std::find(chainBlobs.begin(), chainBlobs.begin() + nChainBlobs, b) != chainBlobs.begin() + nChainBlobs;
  };
  chainBlobs[nChainBlobs++] = seed.id;
  chainBlobs[nChainBlobs++] = cand.blob;

  struct A
  {
    int blob, hit;
    double z, phi;
  };
  std::map<int, A> acc;
  acc[seed.layer] = {seed.id, cand.outer_hit, cand.outer_z, cand.outer_phi};
  acc[cand.layer] = {cand.blob, cand.hit, cand.z, cand.phi};
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;

  SiChainStep s;
  s.kind = SiChainStep::Seed;
  s.from_layer = seed.layer;
  s.to_layer = cand.layer;
  s.blob_id = cand.blob;
  s.hit_id = cand.hit;
  s.ref_z = cand.outer_z;
  s.ref_phi = cand.outer_phi;
  s.pred_z = cand.outer_z + cand.z_center * dz;
  s.pred_phi = cand.outer_phi + cand.phi_center * dp;
  s.z = cand.z;
  s.phi = cand.phi;
  s.delta_z_bins = (cand.z - cand.outer_z) / dz;
  s.delta_phi_bins = wrapDelta(cand.phi - cand.outer_phi, twopi) / dp;
  s.res_z_bins = cand.res_z;
  s.res_phi_bins = cand.res_phi;
  s.z_half_window_bins = cand.z_half;
  s.phi_half_window_bins = m_seedPhiHalf;
  s.score = cand.score;
  out.steps.push_back(s);

  auto propagateTo = [&](int kind, int layer, const std::vector<int>& fitlayers) -> bool
  {
    std::vector<Eigen::Vector3d> pts;
    for (int l : fitlayers)
    {
      pts.push_back(xyz(l, acc.at(l).z, acc.at(l).phi));
    }
    if (pts.size() < 2)
    {
      return false;
    }
    const auto line = fitLine(pts);
    int ref = acc.begin()->first;
    double rd = 1e9;
    for (auto& [l, a] : acc)
    {
      const double d = std::abs(m_radius[l] - m_radius[layer]);
      if (d < rd)
      {
        rd = d;
        ref = l;
      }
    }
    double pz = 0, pp = 0;
    const bool lineOk = propagate(line, layer, xyz(ref, acc[ref].z, acc[ref].phi), pz, pp);
    double phiHalf = m_propPhiHalf;
    bool tpcGuided = false;
    if (kind == SiChainStep::InttPropagation && m_inttCentre != InttCentre::MvtxLine && tpcRef)
    {
      // TPC prediction on this layer + the chain's MVTX offsets to the TPC prediction
      // (mean, or straight-line fit in r / s extrapolated to this layer)
      double sRphi = 0, sZ = 0;
      int nOff = 0, nAlt = 0;
      std::vector<std::array<double, 4>> offs;  // r, r*phi offset, s, z offset
      for (const auto& [l, a] : acc)
      {
        if (l > 2)
        {
          continue;
        }
        const auto& p0 = tpcRef->layers[l];
        const auto& p1 = tpcRef->alt[l];
        const SiTpcReference::LayerPrediction* p = nullptr;
        if (p0.valid && p1.valid)
        {
          const bool useAlt = std::abs(wrapDelta(a.phi - p1.u, twopi)) < std::abs(wrapDelta(a.phi - p0.u, twopi));
          p = useAlt ? &p1 : &p0;
          nAlt += useAlt;
        }
        else if (p0.valid)
        {
          p = &p0;
        }
        if (!p)
        {
          continue;
        }
        const double oR = m_radius[l] * wrapDelta(a.phi - p->u, twopi);
        const double oZ = a.z - p->z;
        sRphi += oR;
        sZ += oZ;
        offs.push_back({m_radius[l], oR, p->s, oZ});
        ++nOff;
      }
      const auto& t0 = tpcRef->layers[layer];
      const auto& t1 = tpcRef->alt[layer];
      const SiTpcReference::LayerPrediction* tp = (t1.valid && 2 * nAlt > nOff) ? &t1 : (t0.valid ? &t0 : nullptr);
      if (nOff > 0 && tp)
      {
        double offR = sRphi / nOff, offZ = sZ / nOff;
        if (m_inttCentre == InttCentre::TpcOffsetSlope && nOff >= 2)
        {
          // least-squares lines offset(r) and offset(s), evaluated at this layer
          auto lineAt = [&](int ix, int iy, double x0, double fallback)
          {
            double mx = 0, my = 0;
            for (const auto& o : offs)
            {
              mx += o[ix];
              my += o[iy];
            }
            mx /= offs.size();
            my /= offs.size();
            double sxx = 0, sxy = 0;
            for (const auto& o : offs)
            {
              sxx += (o[ix] - mx) * (o[ix] - mx);
              sxy += (o[ix] - mx) * (o[iy] - my);
            }
            return sxx > 0.09 ? my + sxy / sxx * (x0 - mx) : fallback;  // lever arm > ~0.3 cm
          };
          offR = lineAt(0, 1, m_radius[layer], offR);
          offZ = lineAt(2, 3, tp->s, offZ);
        }
        pp = tp->u + offR / m_radius[layer];
        if (pp < 0)
        {
          pp += twopi;
        }
        else if (pp >= twopi)
        {
          pp -= twopi;
        }
        pz = tp->z + offZ;
        phiHalf = m_inttTpcPhiHalf;
        tpcGuided = true;
        ++m_inttTpcUsed;
      }
      else
      {
        ++m_inttTpcFallback;
      }
    }
    if (!lineOk && !tpcGuided)
    {
      return false;
    }
    if (layer == 0 && m_hwDdphiProp0 && acc.count(1) && acc.count(2))
    {
      // phi window from the step-to-step correlation instead of the straight line (z unchanged)
      // layer-0 Uphi is not known yet: only dUphi_pred(0) is used, the dUphi(0) entry is a dummy
      const std::map<int, double> uphi{{0, acc.at(1).phi}, {1, acc.at(1).phi}, {2, acc.at(2).phi}};
      double dummy = 0, pred = 0;
      if (hwDdphiTerms(0, uphi, dummy, pred))
      {
        pp = acc.at(1).phi + m_hwDdphiOffset[0] + m_hwDdphiSlope[0] * pred;
        phiHalf = m_hwDdphiNSigma * m_hwDdphiSigma[0] / dp;
        ++m_hwDdphiProp0Used;
      }
    }
    if (kind == SiChainStep::MvtxPropagation && m_mvtxPropZFromTpc && tpcRef && tpcRef->layers[layer].valid &&
        tpcRef->layers[ref].valid)
    {
      // inclination from the TPC: TPC z difference between this layer and the nearest chain layer
      pz = acc[ref].z + (tpcRef->layers[layer].z - tpcRef->layers[ref].z);
      ++m_mvtxTpcZUsed;
    }
    const double zHalf = (kind == SiChainStep::InttPropagation && m_inttPropZHalf >= 0) ? m_inttPropZHalf : m_propZHalf;
    int bestBlob = -1, bestHit = -1;
    double bestScore = 1e99, bz = 0, bp = 0, rz = 0, rp = 0;
    forEachHitNear(layer, (pz - m_zmin) / dz, zHalf, pp / dp, phiHalf, [&](int hid)
    {
      const int b = m_hitBlob[hid];
      if (b < 0 || isUsed(b))
      {
        return;
      }
      const auto& h = m_eventData->hits[hid];
      if (layer >= 3 && m_inttCrossingWindow >= 0 && !m_tpcTrajectoryNodeName.empty())
      {
        ++m_inttTimeChecked;
        const int bucket = static_cast<int>(InttDefs::getTimeBucketId(h.hitsetkey));
        if (std::abs(bucket - m_curCrossing) > m_inttCrossingWindow)
        {
          ++m_inttTimeRejected;
          return;
        }
      }
      const double rrphi = wrapDelta(h.phi - pp, twopi) / dp;
      const double rrz = (h.z - pz) / dz;
      if (std::abs(rrphi) > phiHalf || std::abs(rrz) > zHalf)
      {
        return;
      }
      const double sc = (rrphi / phiHalf) * (rrphi / phiHalf) + (rrz / zHalf) * (rrz / zHalf);
      if (sc < bestScore)
      {
        bestScore = sc;
        bestBlob = b;
        bestHit = hid;
        bz = h.z;
        bp = h.phi;
        rz = rrz;
        rp = rrphi;
      }
    });
    if (bestBlob < 0)
    {
      return false;
    }
    chainBlobs[nChainBlobs++] = bestBlob;
    acc[layer] = {bestBlob, bestHit, bz, bp};
    SiChainStep st;
    st.kind = kind;
    st.from_layer = ref;
    st.to_layer = layer;
    st.blob_id = bestBlob;
    st.hit_id = bestHit;
    st.ref_z = acc[ref].z;
    st.ref_phi = acc[ref].phi;
    st.pred_z = pz;
    st.pred_phi = pp;
    st.z = bz;
    st.phi = bp;
    st.delta_z_bins = (bz - acc[ref].z) / dz;
    st.delta_phi_bins = wrapDelta(bp - acc[ref].phi, twopi) / dp;
    st.res_z_bins = rz;
    st.res_phi_bins = rp;
    st.z_half_window_bins = m_propZHalf;
    st.phi_half_window_bins = phiHalf;
    st.score = bestScore;
    out.steps.push_back(st);
    return true;
  };

  for (int l = 2; l >= 0; --l)
  {
    if (!acc.count(l))
    {
      std::vector<int> fit;
      for (auto& [x, _] : acc)
      {
        if (x <= 2)
        {
          fit.push_back(x);
        }
      }
      propagateTo(SiChainStep::MvtxPropagation, l, fit);
    }
  }
  std::vector<int> mvtx;
  for (auto& [l, _] : acc)
  {
    if (l <= 2)
    {
      mvtx.push_back(l);
    }
  }
  for (int l = 3; l <= 6; ++l)
  {
    std::vector<int> fit = m_refitWithIntt ? std::vector<int>{} : mvtx;
    if (m_refitWithIntt)
    {
      for (auto& [x, _] : acc)
      {
        fit.push_back(x);
      }
    }
    propagateTo(SiChainStep::InttPropagation, l, fit);
  }

  // step-to-step delta-delta-phi check on the complete chain (layers 0, 1)
  double ddphiScore = 0;
  if (m_hwDdphiCheck[0] || m_hwDdphiCheck[1])
  {
    std::map<int, double> uphi;
    for (const auto& [l, a] : acc)
    {
      uphi[l] = a.phi;
    }
    for (int l = 0; l < 2; ++l)
    {
      double d = 0, pred = 0;
      if (!m_hwDdphiCheck[l] || !hwDdphiTerms(l, uphi, d, pred) || m_hwDdphiSigma[l] <= 0)
      {
        continue;
      }
      const double pull = (d - m_hwDdphiOffset[l] - m_hwDdphiSlope[l] * pred) / m_hwDdphiSigma[l];
      out.hw_ddphi_pull[l] = pull;
      ++m_hwDdphiChecked[l];
      if (std::abs(pull) > m_hwDdphiNSigma)
      {
        ++m_hwDdphiRejected[l];
        if (m_hwDdphiReject)
        {
          return false;
        }
      }
      ddphiScore += pull * pull;
    }
  }

  out.n_mvtx = 0;
  out.n_intt = 0;
  for (auto it = acc.rbegin(); it != acc.rend(); ++it)
  {
    out.layers.push_back(it->first);
    out.blob_ids.push_back(it->second.blob);
    out.point_hit_ids.push_back(it->second.hit);
    out.z.push_back(it->second.z);
    out.phi.push_back(it->second.phi);
    if (it->first <= 2)
    {
      ++out.n_mvtx;
    }
    else
    {
      ++out.n_intt;
    }
    for (int h : m_eventData->blobs[it->second.blob].hit_ids)
    {
      out.hit_ids.push_back(h);
    }
  }
  out.score = ddphiScore;
  for (auto& st : out.steps)
  {
    out.score += st.score;
  }
  return (int) out.layers.size() >= m_minLayers && out.n_mvtx >= m_minMvtxLayers && out.n_intt >= m_minInttLayers;
}

// Chain -> physics space (TPC frame) exactly as SiTrajectoryFitter::layerPoints, then residuals
// to the TPC helix at the radius of every point and the agreement chi2/ndf.
void SiHitSeedReco::toPhysics(const SiTpcHelixFit::Result& tpcFit, SiHitChain& chain) const
{
  chain.phys_x.clear();
  chain.phys_y.clear();
  chain.phys_z.clear();
  chain.res_rdphi.clear();
  chain.res_dz.clear();
  chain.res_dz_fit.clear();
  std::vector<double> pathS;  // TPC helix path length at every chain layer
  for (size_t i = 0; i < chain.layers.size(); ++i)
  {
    const int layer = chain.layers[i];
    const auto& blob = m_eventData->blobs[chain.blob_ids[i]];
    double sx = 0, sy = 0, sz = 0;
    for (int hid : blob.hit_ids)
    {
      // 3D hit position on the sensor (detector frame), not the nominal cylinder
      const auto& h = m_eventData->hits[hid];
      sx += h.px;
      sy += h.py;
      sz += h.pz;
    }
    const double n = std::max<size_t>(1, blob.hit_ids.size());
    const double phi = std::atan2(sy, sx);
    const SiDetectorFrame::Point det{sx / n, sy / n, sz / n};
    const auto glob = m_frame.toGlobal(layer, det);
    const auto tpc = m_beamAlign.siHalfToTpcAxis(m_beamAlign.siHalf(phi - m_frame.uphiRotation()), glob);
    chain.phys_x.push_back(tpc.x);
    chain.phys_y.push_back(tpc.y);
    chain.phys_z.push_back(tpc.z);

    double rdphi = std::numeric_limits<double>::quiet_NaN();
    double dzr = std::numeric_limits<double>::quiet_NaN();
    const double rad = std::hypot(tpc.x, tpc.y);
    SiTpcHelixFit::Pos hp;
    double s = 0;
    if (SiTpcHelixFit::atRadius(tpcFit, rad, tpc.z, hp, s))
    {
      rdphi = rad * SiTpcReference::wrapPi(std::atan2(tpc.y, tpc.x) - std::atan2(hp.y, hp.x));
      dzr = tpc.z - hp.z;
    }
    chain.res_rdphi.push_back(rdphi);
    chain.res_dz.push_back(dzr);
    pathS.push_back(s);
  }

  // ---- r*phi part
  double chi2r = 0;
  int nphi = 0;
  for (size_t i = 0; i < chain.layers.size(); ++i)
  {
    const double r = chain.res_rdphi[i];
    if (std::isfinite(r))
    {
      const double sig = std::max(1e-6, m_sigRphiL[std::clamp(chain.layers[i], 0, 6)]);
      chi2r += (r / sig) * (r / sig);
      ++nphi;
    }
  }

  // ---- z part: weighted fit of the model dz = a [+ b s]
  double S = 0, Ss = 0, Sz = 0, Sss = 0, Ssz = 0;
  int nz = 0;
  for (size_t i = 0; i < chain.layers.size(); ++i)
  {
    if (!std::isfinite(chain.res_dz[i]))
    {
      continue;
    }
    const double sig = std::max(1e-6, m_sigZL[std::clamp(chain.layers[i], 0, 6)]);
    const double w = 1.0 / (sig * sig);
    S += w;
    Ss += w * pathS[i];
    Sz += w * chain.res_dz[i];
    Sss += w * pathS[i] * pathS[i];
    Ssz += w * pathS[i] * chain.res_dz[i];
    ++nz;
  }
  double a = 0, b = 0;
  int nzpar = 0;
  if (m_zModel == TpcZModel::Offset && nz > 0)
  {
    a = Sz / S;
    nzpar = 1;
  }
  else if (m_zModel == TpcZModel::OffsetSlope && nz > 0)
  {
    const double D = S * Sss - Ss * Ss;
    if (nz >= 2 && D > 1e-9 * S * S)
    {
      b = (S * Ssz - Ss * Sz) / D;
      a = (Sz - b * Ss) / S;
      nzpar = 2;
    }
    else
    {
      a = Sz / S;
      nzpar = 1;
    }
  }
  const double priorW = 1.0 / (m_zOffsetSigma * m_zOffsetSigma);
  if (m_zModel == TpcZModel::OffsetPrior && nz > 0)
  {
    a = Sz / (S + priorW);  // offset pulled towards 0 by the prior
    nzpar = 0;              // one parameter, one prior constraint
  }
  // Fixed model: nothing is subtracted in the chi2, but the measured Si - TPC z offset is still
  // stored (TpcCrossingResolver needs it to move the chain to the neighbouring crossings)
  // (OffsetPrior: the measured offset too, not the one pulled by the prior)
  chain.tpc_dz_offset = ((m_zModel == TpcZModel::Fixed || m_zModel == TpcZModel::OffsetPrior) && nz > 0) ? Sz / S : a;
  chain.tpc_dz_slope = b;
  double chi2z = 0;
  for (size_t i = 0; i < chain.layers.size(); ++i)
  {
    if (std::isfinite(chain.res_dz[i]))
    {
      const double sig = std::max(1e-6, m_sigZL[std::clamp(chain.layers[i], 0, 6)]);
      const double d = (chain.res_dz[i] - a - b * pathS[i]) / sig;
      chi2z += d * d;
    }
  }
  if (m_zModel == TpcZModel::OffsetPrior && nz > 0)
  {
    chi2z += a * a * priorW;
  }
  chain.res_dz_fit.assign(chain.layers.size(), std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < chain.layers.size(); ++i)
  {
    if (std::isfinite(chain.res_dz[i]))
    {
      chain.res_dz_fit[i] = chain.res_dz[i] - a - b * pathS[i];
    }
  }
  chain.tpc_chi2_rphi = chi2r;
  chain.tpc_chi2_z = chi2z;
  const int ndfz = m_zWeight > 0 ? std::max(0, nz - nzpar) : 0;
  chain.tpc_ndf = nphi + ndfz;
  chain.tpc_nres = nphi;
  chain.tpc_chi2 = chain.tpc_ndf > 0 ? (chi2r + m_zWeight * chi2z) / chain.tpc_ndf : std::numeric_limits<double>::quiet_NaN();
  if (m_requireAllResiduals && nphi < static_cast<int>(chain.layers.size()))
  {
    // the TPC helix does not reach every chain layer: not this TPC track's chain
    chain.tpc_chi2 = std::numeric_limits<double>::quiet_NaN();
  }
}

void SiHitSeedReco::processTpcTrack(const TpcTrackInput& in)
{
  m_curCrossing = in.crossing;
  SiTpcSeedInfo info;
  info.tpc_track_index = in.index;
  info.tpc_track_id = in.id;
  info.tpc_assembled_track_id = in.assembledId;
  info.tpc_crossing = in.crossing;
  info.n_clusters = in.nClusters;
  const bool dbg = debugThisTrack(in);
  const bool dbgEvent = m_debugEvent >= 0 && static_cast<long>(m_eventData->event) == m_debugEvent;
  auto finish = [&](int status)
  {
    if (dbgEvent && !dbg)
    {
      std::printf("[SiSeedDebug] event %llu track index %4d assembled %6u clusters %3u  poly pt %.3f phi %+.4f z %+7.2f q %+.0f | refit pt %.3f phi %+.4f z0 %+7.2f | %s, anchors %d chains %d\n",
                  m_eventData->event, in.index, in.assembledId, in.nClusters, std::hypot(in.polyPx, in.polyPy),
                  std::atan2(in.polyPy, in.polyPx), in.polyZ, in.polyCharge, info.pt, info.phi, info.z0, statusName(status),
                  info.n_anchors, info.n_chains);
    }
    if (dbg)
    {
      std::cout << "[SiSeedDebug] => status " << statusName(status) << ", anchors " << info.n_anchors << ", chains "
                << info.n_chains << ", selected chain " << info.selected_chain << "\n"
                << "[SiSeedDebug] ==================================================" << std::endl;
    }
    info.status = status;
    ++m_statusCount[std::clamp(status, 0, 6)];
    m_eventData->tpc_seeds.push_back(info);
    if (Verbosity() > 1)
    {
      std::cout << Name() << ":   TPC track " << in.index << " (assembled " << in.assembledId << ", "
                << in.nClusters << " clusters, pt " << info.pt << ", z0 " << info.z0 << "): " << statusName(status)
                << ", anchors " << info.n_anchors << ", chains " << info.n_chains << ", selected chain "
                << info.selected_chain << std::endl;
    }
  };

  if (in.nClusters < m_minTpcClusters)
  {
    if (dbg)
    {
      debugTrackHeader(in, nullptr);
    }
    finish(SiTpcSeedInfo::TooFewClusters);
    return;
  }
  if (in.fit.fitted())
  {
    info.pt = in.fit.pt;
    info.phi = in.fit.phi;
    info.tanl = in.fit.tanl;
    info.z0 = in.fit.z0;
    info.charge = in.fit.charge;
  }
  SiTpcReference::TpcReference ref;
  const bool refOk = SiTpcReference::buildTpcReference(in.fit, m_map, ref, m_boundaryMargin);
  if (dbg)
  {
    debugTrackHeader(in, refOk ? &ref : nullptr);
  }
  if (!refOk)
  {
    finish(SiTpcSeedInfo::FitFailed);
    return;
  }

  // ---- chains of this TPC track (hits may be shared, identical chains once)
  std::vector<SiHitChain> chains;
  std::set<std::vector<int>> signatures;
  bool anyPrediction = false;
  std::vector<std::pair<int, int>> pairs{{2, 1}};
  if (m_useInnerPair)
  {
    pairs.emplace_back(1, 0);
  }
  for (const auto& [outer, inner] : pairs)
  {
    if (chains.size() >= m_maxChainsPerTrack)
    {
      break;
    }
    if (!ref.layers[outer].valid)
    {
      continue;
    }
    anyPrediction = true;
    const auto anchors = anchorBlobs(ref, outer);
    info.n_anchors += static_cast<int>(anchors.size());
    if (dbg)
    {
      std::cout << "[SiSeedDebug] --- seed pair " << outer << "->" << inner << std::endl;
      debugAnchors(ref, outer, anchors);
    }
    for (int a : anchors)
    {
      if (chains.size() >= m_maxChainsPerTrack)
      {
        if (dbg)
        {
          std::cout << "[SiSeedDebug]   chain cap " << m_maxChainsPerTrack << " reached, remaining anchors skipped" << std::endl;
        }
        break;
      }
      const auto& seed = m_eventData->blobs[a];
      const auto cands = seedCandidates(seed, inner, ref);
      if (dbg)
      {
        debugSeed(seed, inner, ref, cands);
      }
      const size_t ncand = std::min<size_t>(cands.size(), m_maxSeedCand);
      for (size_t k = 0; k < ncand && chains.size() < m_maxChainsPerTrack; ++k)
      {
        SiHitChain ch;
        const bool built = buildChain(seed, cands[k], ch, &ref);
        if (dbg)
        {
          std::cout << "[SiSeedDebug]     build anchor " << a << " + cand blob " << cands[k].blob << ": ";
          debugChain(ch, built ? "OK" : "FAILED requirements");
        }
        if (!built)
        {
          continue;
        }
        std::vector<int> sig = ch.blob_ids;
        std::sort(sig.begin(), sig.end());
        if (!signatures.insert(sig).second)
        {
          if (dbg)
          {
            std::cout << "[SiSeedDebug]       (duplicate of an earlier chain, dropped)" << std::endl;
          }
          continue;
        }
        ch.seed_outer_layer = outer;
        ch.anchor_blob = a;
        chains.push_back(std::move(ch));
      }
    }
  }
  info.n_chains = static_cast<int>(chains.size());
  if (!anyPrediction)
  {
    finish(SiTpcSeedInfo::NoPrediction);
    return;
  }
  if (info.n_anchors == 0)
  {
    finish(SiTpcSeedInfo::NoAnchor);
    return;
  }
  if (chains.empty())
  {
    finish(SiTpcSeedInfo::NoChain);
    return;
  }

  // ---- physics space, agreement with the TPC helix, ranking
  for (auto& ch : chains)
  {
    toPhysics(in.fit, ch);
  }
  // previous ordering: most layers, then tpc_chi2 (finite first), then score
  auto layersFirst = [](const SiHitChain& a, const SiHitChain& b)
  {
    if (a.n_layers() != b.n_layers())
    {
      return a.n_layers() > b.n_layers();
    }
    const bool fa = std::isfinite(a.tpc_chi2), fb = std::isfinite(b.tpc_chi2);
    if (fa != fb)
    {
      return fa;
    }
    if (fa && a.tpc_chi2 != b.tpc_chi2)
    {
      return a.tpc_chi2 < b.tpc_chi2;
    }
    return a.score < b.score;
  };
  const double agreeCut = m_agreementMaxChi2;
  auto agrees = [agreeCut](const SiHitChain& c)
  {
    return std::isfinite(c.tpc_chi2) && c.tpc_chi2 <= agreeCut;
  };
  // agreeing chains first (most layers among them), then the rest by tpc_chi2
  auto agreementFirst = [&](const SiHitChain& a, const SiHitChain& b)
  {
    const bool aa = agrees(a), ab = agrees(b);
    if (aa != ab)
    {
      return aa;
    }
    if (aa)
    {
      return layersFirst(a, b);
    }
    const bool fa = std::isfinite(a.tpc_chi2), fb = std::isfinite(b.tpc_chi2);
    if (fa != fb)
    {
      return fa;
    }
    if (fa && a.tpc_chi2 != b.tpc_chi2)
    {
      return a.tpc_chi2 < b.tpc_chi2;
    }
    if (a.n_layers() != b.n_layers())
    {
      return a.n_layers() > b.n_layers();
    }
    return a.score < b.score;
  };

  bool rankZeroAgrees = true;
  if (m_chainRanking == ChainRanking::AgreementFirst)
  {
    const auto oldBest = std::min_element(chains.begin(), chains.end(), layersFirst);
    const std::vector<int> oldBlobs = oldBest->blob_ids;
    std::sort(chains.begin(), chains.end(), agreementFirst);
    rankZeroAgrees = agrees(chains[0]);
    if (!rankZeroAgrees)
    {
      ++m_noAgreeingChain;
    }
    if (chains[0].blob_ids != oldBlobs)
    {
      ++m_rankChanged;
    }
  }
  else
  {
    std::sort(chains.begin(), chains.end(), layersFirst);
  }
  for (size_t r = 0; r < chains.size(); ++r)
  {
    chains[r].rank = static_cast<int>(r);
  }
  const bool pass = rankZeroAgrees &&
                    (m_maxTpcChi2 <= 0 || (std::isfinite(chains[0].tpc_chi2) && chains[0].tpc_chi2 <= m_maxTpcChi2));
  chains[0].selected = pass;

  // ---- store
  for (auto& ch : chains)
  {
    if (m_storeSelectedOnly && !m_exclusiveChains && !ch.selected)  // exclusive: pruned after arbitration
    {
      continue;
    }
    ch.id = static_cast<int>(m_eventData->chains.size());
    ch.tpc_track_index = in.index;
    ch.tpc_track_id = in.id;
    ch.tpc_assembled_track_id = in.assembledId;
    ch.tpc_crossing = in.crossing;
    if (ch.selected)
    {
      info.selected_chain = ch.id;
    }
    if (dbg)
    {
      std::cout << "[SiSeedDebug] stored chain " << ch.id << " rank " << ch.rank << (ch.selected ? " SELECTED" : "")
                << " anchor blob " << ch.anchor_blob << " seed " << ch.seed_outer_layer << "->" << ch.seed_outer_layer - 1
                << " chi2/ndf " << ch.tpc_chi2 << " (ndf " << ch.tpc_ndf << ", TPC residual on " << ch.tpc_nres << "/"
                << ch.layers.size() << " layers; chi2 r*phi " << ch.tpc_chi2_rphi << ", z " << ch.tpc_chi2_z
                << ") dz offset " << ch.tpc_dz_offset << " cm, slope " << ch.tpc_dz_slope << ", hw ddphi pull L0 " << ch.hw_ddphi_pull[0] << " L1 " << ch.hw_ddphi_pull[1] << std::endl;
      for (size_t i = 0; i < ch.layers.size(); ++i)
      {
        std::printf("[SiSeedDebug]     layer %d blob %5d  Uphi %8.5f Vz %8.3f | TPC frame x %8.3f y %8.3f z %8.3f | res r*dphi %+8.4f cm  dz %+8.4f cm\n",
                    ch.layers[i], ch.blob_ids[i], ch.phi[i], ch.z[i],
                    i < ch.phys_x.size() ? ch.phys_x[i] : NAN, i < ch.phys_y.size() ? ch.phys_y[i] : NAN,
                    i < ch.phys_z.size() ? ch.phys_z[i] : NAN,
                    i < ch.res_rdphi.size() ? ch.res_rdphi[i] : NAN, i < ch.res_dz.size() ? ch.res_dz[i] : NAN);
      }
    }
    if (Verbosity() > 2)
    {
      std::cout << Name() << ":     chain " << ch.id << " rank " << ch.rank << (ch.selected ? " SELECTED" : "")
                << " layers " << ch.n_layers() << " (MVTX " << ch.n_mvtx << ", INTT " << ch.n_intt << ") chi2/ndf "
                << ch.tpc_chi2 << " (ndf " << ch.tpc_ndf << ") dz offset " << ch.tpc_dz_offset << " cm, seed pair "
                << ch.seed_outer_layer << "->" << ch.seed_outer_layer - 1 << std::endl;
    }
    m_eventData->chains.push_back(std::move(ch));
    ++m_nChains;
  }
  finish(pass ? SiTpcSeedInfo::Chained : SiTpcSeedInfo::Rejected);
}
// ============================================================================================
// Worker interface (SiSeedingOrchestrator): pass 2 (rechain on free blobs), pass 3 (hit recovery)
// ============================================================================================
double SiHitSeedReco::tpcInputZ0(int tpcInputIndex) const
{
  const auto it = m_tpcInputPos.find(tpcInputIndex);
  if (it == m_tpcInputPos.end())
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const auto& in = m_tpcInputs[it->second];
  return in.fit.fitted() ? in.fit.z0 : std::numeric_limits<double>::quiet_NaN();
}

void SiHitSeedReco::updateHwDdphiPulls(SiHitChain& chain) const
{
  std::map<int, double> uphi;
  for (std::size_t i = 0; i < chain.layers.size() && i < chain.phi.size(); ++i)
  {
    uphi[chain.layers[i]] = chain.phi[i];
  }
  for (int l = 0; l < 2; ++l)
  {
    chain.hw_ddphi_pull[l] = std::numeric_limits<double>::quiet_NaN();
    double d = 0, pred = 0;
    if (!m_hwDdphiCheck[l] || m_hwDdphiSigma[l] <= 0 || !hwDdphiTerms(l, uphi, d, pred))
    {
      continue;
    }
    chain.hw_ddphi_pull[l] = (d - m_hwDdphiOffset[l] - m_hwDdphiSlope[l] * pred) / m_hwDdphiSigma[l];
  }
}


std::vector<int> SiHitSeedReco::rechain(const std::set<unsigned int>& assembledIds, const std::vector<char>& excludedBlobs,
                                        const RechainConfig& cfg, int passTag)
{
  std::set<int> inputs;
  for (const auto& in : m_tpcInputs)
  {
    if (assembledIds.count(in.assembledId) != 0)
    {
      inputs.insert(in.index);
    }
  }
  return rechainInputs(inputs, excludedBlobs, cfg, passTag);
}

std::vector<int> SiHitSeedReco::rechainInputs(const std::set<int>& tpcInputIndices, const std::vector<char>& excludedBlobs,
                                              const RechainConfig& cfg, int passTag)
{
  std::vector<int> out;
  if (!m_eventData || tpcInputIndices.empty())
  {
    return out;
  }
  ++m_rechainCalls;

  // pass-1 settings and run counters, restored at the end (pass-1 statistics stay pass-1 only)
  const int minLayers = m_minLayers, minMvtx = m_minMvtxLayers, minIntt = m_minInttLayers;
  const double anchorRphi = m_anchorRphi, anchorZ = m_anchorZ, seedZHalf = m_seedZHalfCm, seedPhiHalf = m_seedPhiHalf;
  const double propPhiHalf = m_propPhiHalf, propZHalf = m_propZHalf, inttPropZHalf = m_inttPropZHalf;
  const double inttTpcPhiHalf = m_inttTpcPhiHalf, agreeMax = m_agreementMaxChi2;
  const unsigned int maxChains = m_maxChainsPerTrack, maxAnchors = m_maxAnchors;
  const bool storeSelectedOnly = m_storeSelectedOnly;
  const auto statusCount = m_statusCount;
  const unsigned long nChainsRun = m_nChains, rankChanged = m_rankChanged, noAgreeing = m_noAgreeingChain;

  m_minLayers = cfg.minLayers.value_or(m_minLayers);
  m_minMvtxLayers = cfg.minMvtxLayers.value_or(m_minMvtxLayers);
  m_minInttLayers = cfg.minInttLayers.value_or(m_minInttLayers);
  m_anchorRphi = cfg.anchorRphiCm.value_or(m_anchorRphi);
  m_anchorZ = cfg.anchorZCm.value_or(m_anchorZ);
  m_seedZHalfCm = cfg.seedZHalfCm.value_or(m_seedZHalfCm);
  m_seedPhiHalf = cfg.seedPhiHalfBins.value_or(m_seedPhiHalf);
  m_propPhiHalf = cfg.propPhiHalfBins.value_or(m_propPhiHalf);
  m_propZHalf = cfg.propZHalfBins.value_or(m_propZHalf);
  m_inttPropZHalf = cfg.inttPropZHalfBins.value_or(m_inttPropZHalf);
  m_inttTpcPhiHalf = cfg.inttTpcPhiHalfBins.value_or(m_inttTpcPhiHalf);
  m_agreementMaxChi2 = cfg.agreementMaxChi2.value_or(m_agreementMaxChi2);
  m_maxChainsPerTrack = cfg.maxChainsPerTrack.value_or(m_maxChainsPerTrack);
  m_maxAnchors = cfg.maxAnchorsPerTrack.value_or(m_maxAnchors);
  m_storeSelectedOnly = false;  // every candidate is stored: the caller selects

  buildHitGrid(&excludedBlobs);
  const std::size_t nChains0 = m_eventData->chains.size();
  const std::size_t nSeeds0 = m_eventData->tpc_seeds.size();
  unsigned int nInputs = 0;
  for (const auto& in : m_tpcInputs)
  {
    if (tpcInputIndices.count(in.index) != 0)
    {
      processTpcTrack(in);
      ++nInputs;
    }
  }
  m_rechainInputs += nInputs;
  // tpc_seeds keeps one entry per TPC input (pass 1); the new chains carry their TPC input themselves
  m_eventData->tpc_seeds.erase(m_eventData->tpc_seeds.begin() + static_cast<std::ptrdiff_t>(nSeeds0), m_eventData->tpc_seeds.end());
  for (std::size_t c = nChains0; c < m_eventData->chains.size(); ++c)
  {
    auto& ch = m_eventData->chains[c];
    ch.pass = passTag;
    ch.selected = false;
    out.push_back(ch.id);
  }
  m_rechainChains += out.size();

  m_minLayers = minLayers;
  m_minMvtxLayers = minMvtx;
  m_minInttLayers = minIntt;
  m_anchorRphi = anchorRphi;
  m_anchorZ = anchorZ;
  m_seedZHalfCm = seedZHalf;
  m_seedPhiHalf = seedPhiHalf;
  m_propPhiHalf = propPhiHalf;
  m_propZHalf = propZHalf;
  m_inttPropZHalf = inttPropZHalf;
  m_inttTpcPhiHalf = inttTpcPhiHalf;
  m_agreementMaxChi2 = agreeMax;
  m_maxChainsPerTrack = maxChains;
  m_maxAnchors = maxAnchors;
  m_storeSelectedOnly = storeSelectedOnly;
  m_statusCount = statusCount;
  m_nChains = nChainsRun;
  m_rankChanged = rankChanged;
  m_noAgreeingChain = noAgreeing;

  buildHitGrid();  // full grid again
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": rechain pass " << passTag << ": " << nInputs << " TPC inputs, " << out.size()
              << " new candidate chains" << std::endl;
  }
  return out;
}

// Vertex-assisted prediction on `layer`: helix with the TPC curvature (synthetic points on the TPC helix
// between 25 and 75 cm, circle fit), centre refitted with the radius fixed to the chain's physics points
// and the beam (SiTpcHelixFit hybrid), z(s) from the MVTX points (+ INTT, small weight) and the vertex z.
// The crossing with the layer is searched at the radius where the DETECTOR-frame radius of the point
// equals the mean detector-frame radius of this event's hits on the layer (INTT: survey radii), and
// converted back to hardware (Uphi about the layer centre - rotation, Vz), as fillHits does for the hits.
bool SiHitSeedReco::predictWithVertex(const SiHitChain& chain, const SiTpcHelixFit::Result& tpcFit, int layer, double vertexZ,
                                      const RecoveryConfig& cfg, const std::array<double, 7>& layerRadius, double& z, double& phi) const
{
  if (!(m_useActsGeometry && m_hwFromActs) || !tpcFit.fitted() || layerRadius[layer] <= 0 || chain.phys_x.size() != chain.layers.size())
  {
    return false;
  }
  std::vector<SiTpcHelixFit::Point> pts;
  const double bx = m_beamAlign.referenceBeamX(), by = m_beamAlign.referenceBeamY();
  {
    SiTpcHelixFit::Point b;
    b.x = bx;
    b.y = by;
    b.z = std::isfinite(vertexZ) ? vertexZ : 0.0;
    b.wxy = 0;
    b.wz = std::isfinite(vertexZ) ? cfg.vertexZWeight : 0.0;
    b.wpos = cfg.beamPosWeight;
    pts.push_back(b);
  }
  std::vector<SiTpcHelixFit::Point> si;
  for (std::size_t i = 0; i < chain.layers.size(); ++i)
  {
    const bool intt = chain.layers[i] >= 3;
    SiTpcHelixFit::Point p;
    p.x = chain.phys_x[i];
    p.y = chain.phys_y[i];
    p.z = chain.phys_z[i];
    p.wxy = 0;
    p.wz = intt ? 1.0e-4 : 1.0;
    p.wpos = intt ? cfg.inttPosWeight : 1.0;
    si.push_back(p);
  }
  std::sort(si.begin(), si.end(), [](const auto& a, const auto& b) { return a.x * a.x + a.y * a.y < b.x * b.x + b.y * b.y; });
  pts.insert(pts.end(), si.begin(), si.end());
  unsigned int nTpc = 0;
  for (double rad = 25.0; rad <= 75.0; rad += 10.0)
  {
    SiTpcHelixFit::Pos hp;
    double s = 0;
    if (!SiTpcHelixFit::atRadius(tpcFit, rad, tpcFit.z0 + tpcFit.tanl * rad, hp, s))
    {
      continue;
    }
    SiTpcHelixFit::Point p;
    p.x = hp.x;
    p.y = hp.y;
    p.z = hp.z;
    p.wxy = 1;
    p.wz = 0;
    p.wpos = 0;
    pts.push_back(p);
    ++nTpc;
  }
  if (nTpc < 3)
  {
    return false;
  }
  SiTpcHelixFit::Config fc;
  fc.beamX = bx;
  fc.beamY = by;
  fc.bz = m_bz;
  fc.hybridMinPoints = 2;
  const auto r = SiTpcHelixFit::fit(pts, fc);
  if (!r.fitted() || !r.hybrid)
  {
    return false;
  }
  // nearest chain point: z reference for the crossing (the helix crosses the radius twice per turn)
  double zRef = chain.phys_z.front();
  double best = 1e9;
  for (std::size_t i = 0; i < chain.layers.size(); ++i)
  {
    const double d = std::abs(m_radius[chain.layers[i]] - m_radius[layer]);
    if (d < best)
    {
      best = d;
      zRef = chain.phys_z[i];
    }
  }
  const double target = layerRadius[layer];
  double rho = target;
  SiDetectorFrame::Point det;
  bool ok = false;
  for (int it = 0; it < 5; ++it)
  {
    SiTpcHelixFit::Pos hp;
    double s = 0;
    if (!SiTpcHelixFit::atRadius(r, rho, zRef, hp, s))
    {
      return false;
    }
    const auto g = m_beamAlign.siTpcFrameToGlobal(m_frame, {hp.x, hp.y, hp.z});
    det = m_frame.toDetectorFromGlobal(layer, g);
    const double rd = std::hypot(det.x, det.y);
    ok = true;
    if (std::abs(rd - target) < 1e-3)
    {
      break;
    }
    rho += target - rd;
  }
  if (!ok)
  {
    return false;
  }
  // physics detector frame -> hardware (inverse of fillHits: INTT physics offset, per-layer z offset)
  double hx = det.x, hy = det.y;
  if (layer >= 3)
  {
    hx -= m_inttPhysOffset[0];
    hy -= m_inttPhysOffset[1];
  }
  double u = std::fmod(std::atan2(hy, hx) - m_frame.uphiRotation(), twopi);
  if (u < 0)
  {
    u += twopi;
  }
  phi = u;
  z = det.z - m_layerZOffset[layer];
  return true;
}

std::array<unsigned int, 7> SiHitSeedReco::recoverHits(const std::vector<int>& chainIds, const std::vector<char>& excludedBlobs,
                                                       const RecoveryConfig& cfg, const std::vector<double>* vertexZ)
{
  std::array<unsigned int, 7> added{};
  if (!m_eventData || chainIds.empty())
  {
    return added;
  }
  auto& ev = *m_eventData;
  std::vector<char> taken(excludedBlobs);  // + blobs added by this call
  taken.resize(ev.blobs.size(), 0);
  buildHitGrid(&excludedBlobs);
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  const int savedCrossing = m_curCrossing;

  // mean detector-frame radius of this event's hits per layer (physics positions: INTT survey radii)
  std::array<double, 7> layerRadius{};
  {
    std::array<double, 7> n{};
    for (const auto& h : ev.hits)
    {
      if (h.layer >= 0 && h.layer < nlayers)
      {
        layerRadius[h.layer] += std::hypot(h.px, h.py);
        n[h.layer] += 1;
      }
    }
    for (int l = 0; l < nlayers; ++l)
    {
      layerRadius[l] = n[l] > 0 ? layerRadius[l] / n[l] : 0.0;
    }
  }

  for (const int cid : chainIds)
  {
    if (cid < 0 || cid >= static_cast<int>(ev.chains.size()))
    {
      continue;
    }
    SiHitChain& ch = ev.chains[cid];
    if (ch.layers.empty())
    {
      continue;
    }
    const auto ip = m_tpcInputPos.find(ch.tpc_track_index);
    if (ip == m_tpcInputPos.end())
    {
      continue;
    }
    const TpcTrackInput& in = m_tpcInputs[ip->second];
    SiTpcReference::TpcReference ref;
    if (!in.fit.fitted() || !SiTpcReference::buildTpcReference(in.fit, m_map, ref, m_boundaryMargin))
    {
      continue;
    }
    const double vz = (vertexZ && cid < static_cast<int>(vertexZ->size())) ? (*vertexZ)[cid] : std::numeric_limits<double>::quiet_NaN();
    m_curCrossing = ch.tpc_crossing;  // INTT time window of the (resolved) crossing
    for (const int b : ch.blob_ids)
    {
      if (b >= 0 && b < static_cast<int>(taken.size()))
      {
        taken[b] = 1;
      }
    }
    // agreement before, at the TPC input the chain belongs to now (the resolver may have moved it)
    SiHitChain before = ch;
    toPhysics(in.fit, before);
    double chi2Before = before.tpc_chi2;

    // layers to search
    auto has = [&](int l) { return std::find(ch.layers.begin(), ch.layers.end(), l) != ch.layers.end(); };
    std::vector<int> toSearch;
    if (cfg.mvtx)
    {
      for (int l = 0; l <= 2; ++l)
      {
        if (!has(l))
        {
          toSearch.push_back(l);
        }
      }
    }
    if (cfg.intt)
    {
      for (const auto& [a, b] : {std::pair<int, int>{3, 4}, std::pair<int, int>{5, 6}})
      {
        const bool ha = has(a), hb = has(b);
        if (cfg.inttEmptyPairsOnly && (ha || hb))
        {
          continue;
        }
        if (!ha)
        {
          toSearch.push_back(a);
        }
        if (!hb)
        {
          toSearch.push_back(b);
        }
      }
    }

    for (const int layer : toSearch)
    {
      if (has(layer))
      {
        continue;  // added by an earlier search of this chain
      }
      const bool intt = layer >= 3;
      ++m_recoveryTried;

      // nearest chain layer (reference of the step, and the 'near' point of the line propagation)
      std::size_t iNear = 0;
      for (std::size_t i = 1; i < ch.layers.size(); ++i)
      {
        if (std::abs(m_radius[ch.layers[i]] - m_radius[layer]) < std::abs(m_radius[ch.layers[iNear]] - m_radius[layer]))
        {
          iNear = i;
        }
      }

      double pz = 0, pp = 0;
      bool havePred = cfg.useVertex && predictWithVertex(ch, in.fit, layer, vz, cfg, layerRadius, pz, pp);
      if (havePred)
      {
        ++m_recoveryVertexPred;
      }
      if (!havePred && !intt)
      {
        // missing MVTX layer: straight line through the MVTX chain points, as the pass-1 propagation
        std::vector<Eigen::Vector3d> pts;
        for (std::size_t i = 0; i < ch.layers.size(); ++i)
        {
          if (ch.layers[i] <= 2)
          {
            pts.push_back(xyz(ch.layers[i], ch.z[i], ch.phi[i]));
          }
        }
        if (pts.size() >= 2)
        {
          havePred = propagate(fitLine(pts), layer, xyz(ch.layers[iNear], ch.z[iNear], ch.phi[iNear]), pz, pp);
        }
        m_recoveryFallbackPred += havePred;
      }
      if (!havePred)
      {
        // TPC prediction on this layer + mean offset of the chain points to the TPC prediction
        // (each point with the clamshell translation closer to it, as the pass-1 INTT window)
        double sR = 0, sZ = 0;
        int n = 0, nAlt = 0;
        for (std::size_t i = 0; i < ch.layers.size(); ++i)
        {
          const int l = ch.layers[i];
          const auto& p0 = ref.layers[l];
          const auto& p1 = ref.alt[l];
          const SiTpcReference::LayerPrediction* p = nullptr;
          if (p0.valid && p1.valid)
          {
            const bool useAlt = std::abs(wrapDelta(ch.phi[i] - p1.u, twopi)) < std::abs(wrapDelta(ch.phi[i] - p0.u, twopi));
            p = useAlt ? &p1 : &p0;
            nAlt += useAlt;
          }
          else if (p0.valid)
          {
            p = &p0;
          }
          if (!p)
          {
            continue;
          }
          sR += m_radius[l] * wrapDelta(ch.phi[i] - p->u, twopi);
          sZ += ch.z[i] - p->z;
          ++n;
        }
        const auto& t0 = ref.layers[layer];
        const auto& t1 = ref.alt[layer];
        const SiTpcReference::LayerPrediction* tp = (t1.valid && 2 * nAlt > n) ? &t1 : (t0.valid ? &t0 : nullptr);
        if (n > 0 && tp)
        {
          pp = tp->u + (sR / n) / m_radius[layer];
          pp = std::fmod(pp, twopi);
          if (pp < 0)
          {
            pp += twopi;
          }
          pz = tp->z + sZ / n;
          havePred = true;
          ++m_recoveryFallbackPred;
        }
      }
      if (!havePred)
      {
        ++m_recoveryNoPred;
        continue;
      }

      const double phiHalf = intt ? cfg.inttPhiHalfBins : cfg.mvtxPhiHalfBins;
      const double zHalf = intt ? cfg.inttZHalfBins : cfg.mvtxZHalfBins;
      int bestBlob = -1, bestHit = -1;
      double bestScore = 1e99, bz = 0, bp = 0, rz = 0, rp = 0;
      forEachHitNear(layer, (pz - m_zmin) / dz, zHalf, pp / dp, phiHalf, [&](int hid)
      {
        const int b = m_hitBlob[hid];
        if (b < 0 || taken[b])
        {
          return;
        }
        const auto& h = ev.hits[hid];
        if (intt && m_inttCrossingWindow >= 0 && !m_tpcTrajectoryNodeName.empty())
        {
          const int bucket = static_cast<int>(InttDefs::getTimeBucketId(h.hitsetkey));
          if (std::abs(bucket - m_curCrossing) > m_inttCrossingWindow)
          {
            return;
          }
        }
        const double rrphi = wrapDelta(h.phi - pp, twopi) / dp;
        const double rrz = (h.z - pz) / dz;
        if (std::abs(rrphi) > phiHalf || std::abs(rrz) > zHalf)
        {
          return;
        }
        const double sc = (rrphi / phiHalf) * (rrphi / phiHalf) + (rrz / zHalf) * (rrz / zHalf);
        if (sc < bestScore)
        {
          bestScore = sc;
          bestBlob = b;
          bestHit = hid;
          bz = h.z;
          bp = h.phi;
          rz = rrz;
          rp = rrphi;
        }
      });
      if (bestBlob < 0)
      {
        ++m_recoveryNoBlob;
        continue;
      }

      // trial chain with the new layer (per-layer vectors stay ordered outermost first)
      SiHitChain trial = ch;
      std::size_t pos = 0;
      while (pos < trial.layers.size() && trial.layers[pos] > layer)
      {
        ++pos;
      }
      const auto at = [pos](auto& v) { return v.begin() + static_cast<std::ptrdiff_t>(pos); };
      trial.layers.insert(at(trial.layers), layer);
      trial.blob_ids.insert(at(trial.blob_ids), bestBlob);
      trial.point_hit_ids.insert(at(trial.point_hit_ids), bestHit);
      trial.z.insert(at(trial.z), bz);
      trial.phi.insert(at(trial.phi), bp);
      for (const int h : ev.blobs[bestBlob].hit_ids)
      {
        trial.hit_ids.push_back(h);
      }
      ++(intt ? trial.n_intt : trial.n_mvtx);
      ++trial.n_recovered_layers;
      SiChainStep st;
      st.kind = intt ? SiChainStep::InttPropagation : SiChainStep::MvtxPropagation;
      st.from_layer = ch.layers[iNear];
      st.to_layer = layer;
      st.blob_id = bestBlob;
      st.hit_id = bestHit;
      st.ref_z = ch.z[iNear];
      st.ref_phi = ch.phi[iNear];
      st.pred_z = pz;
      st.pred_phi = pp;
      st.z = bz;
      st.phi = bp;
      st.delta_z_bins = (bz - ch.z[iNear]) / dz;
      st.delta_phi_bins = wrapDelta(bp - ch.phi[iNear], twopi) / dp;
      st.res_z_bins = rz;
      st.res_phi_bins = rp;
      st.z_half_window_bins = zHalf;
      st.phi_half_window_bins = phiHalf;
      st.score = bestScore;
      trial.steps.push_back(st);
      trial.score += bestScore;
      toPhysics(in.fit, trial);
      const bool agrees = std::isfinite(trial.tpc_chi2) && trial.tpc_chi2 <= cfg.maxChi2 &&
                          (!std::isfinite(chi2Before) || trial.tpc_chi2 <= chi2Before + cfg.maxChi2Increase);
      if (!agrees)
      {
        ++m_recoveryBadChi2;
        continue;
      }
      updateHwDdphiPulls(trial);
      chi2Before = trial.tpc_chi2;
      ch = std::move(trial);
      taken[bestBlob] = 1;
      ++added[layer];
      ++m_recoveredPerLayer[layer];
    }
  }

  m_curCrossing = savedCrossing;
  buildHitGrid();  // full grid again
  return added;
}
