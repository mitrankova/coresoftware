#include "Full_PolyTrackMatcher.h"

#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrackv1.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyClusterContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include "Tpc_FittingTools.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <TFile.h>
#include <TH3F.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>

namespace
{
  float finite_float(double value)
  {
    return std::isfinite(value) ? static_cast<float>(value) : std::numeric_limits<float>::quiet_NaN();
  }

  double square(double value)
  {
    return value * value;
  }
}

namespace
{
  double sagittaModelDerivative(double xrot, double x0, double invR)
  {
    const double dx = xrot - x0;
    const double dx2 = dx * dx;
    const double invR2 = invR * invR;
    const double invR3 = invR2 * invR;
    const double invR5 = invR3 * invR2;
    return -invR * dx - 0.5 * invR3 * dx2 * dx - 0.375 * invR5 * dx2 * dx2 * dx;
  }
}

Full_PolyTrackMatcher::Full_PolyTrackMatcher(const std::string& name)
  : SubsysReco(name)
{
}

Full_PolyTrackMatcher::~Full_PolyTrackMatcher()
{
  if (m_qaFile)
  {
    m_qaFile->Close();
    delete m_qaFile;
    m_qaFile = nullptr;
  }
}

int Full_PolyTrackMatcher::Init(PHCompositeNode*)
{
  createQaObjects();
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackMatcher::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackMatcher::getNodes(PHCompositeNode* topNode)
{
  m_tpcTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName);
  if (!m_tpcTracks)
  {
    std::cerr << Name() << "::getNodes - missing " << m_tpcTrackNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tpcClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  if (!m_tpcClusters)
  {
    std::cerr << Name() << "::getNodes - missing " << m_tpcClusterNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  if (!m_crossingDecisions && Verbosity() > 0)
  {
    std::cout << Name() << "::getNodes - optional " << m_crossingDecisionNodeName << " not found" << std::endl;
  }

  m_trkrClusters = findNode::getClass<TrkrClusterContainer>(topNode, m_trkrClusterNodeName);
  if (!m_trkrClusters)
  {
    std::cerr << Name() << "::getNodes - missing " << m_trkrClusterNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_actsGeometry)
  {
    std::cerr << Name() << "::getNodes - missing ActsGeometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackMatcher::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_fullTracks = findNode::getClass<Full_PolyTrackContainer>(topNode, m_outputNodeName);
  if (!m_fullTracks)
  {
    m_fullTracks = new Full_PolyTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_fullTracks, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void Full_PolyTrackMatcher::createQaObjects()
{
  if (!m_writeQA || m_qaFile)
  {
    return;
  }

  m_qaFile = new TFile(m_qaFileName.c_str(), "RECREATE");
  if (!m_qaFile || m_qaFile->IsZombie())
  {
    std::cerr << Name() << "::createQaObjects - cannot open QA file " << m_qaFileName << std::endl;
  }
}

double Full_PolyTrackMatcher::wrapPhi(double phi) const
{
  while (phi > M_PI) { phi -= 2.0 * M_PI; }
  while (phi <= -M_PI) { phi += 2.0 * M_PI; }
  return phi;
}

double Full_PolyTrackMatcher::unwrapPhiNear(double phi, double reference) const
{
  while (phi - reference > M_PI) { phi -= 2.0 * M_PI; }
  while (phi - reference <= -M_PI) { phi += 2.0 * M_PI; }
  return phi;
}

unsigned int Full_PolyTrackMatcher::layerBit(unsigned int layer) const
{
  return layer < 32U ? (1U << layer) : 0U;
}

bool Full_PolyTrackMatcher::getGlobalClusterPosition(TrkrDefs::cluskey key, TrkrCluster* cluster, SpacePoint& point) const
{
  if (!cluster || !m_actsGeometry)
  {
    return false;
  }

  const auto global = m_actsGeometry->getGlobalPosition(key, cluster);
  point.key = key;
  point.layer = TrkrDefs::getLayer(key);
  point.x = global.x();
  point.y = global.y();
  point.z = global.z();
  point.r = std::hypot(point.x, point.y);
  point.phi = std::atan2(point.y, point.x);

  const auto trkrid = TrkrDefs::getTrkrId(key);
  if (trkrid == TrkrDefs::inttId)
  {
    point.ladder = InttDefs::getLadderPhiId(key);
    point.sensor = InttDefs::getLadderZId(key);
  }
  else if (trkrid == TrkrDefs::mvtxId)
  {
    point.ladder = MvtxDefs::getStaveId(key);
    point.sensor = MvtxDefs::getChipId(key);
  }
  else
  {
    point.ladder = TrkrDefs::getPhiElement(key);
    point.sensor = TrkrDefs::getZElement(key);
  }

  return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) && point.r > 0.0;
}

std::vector<Full_PolyTrackMatcher::SpacePoint> Full_PolyTrackMatcher::collectSiliconClusters() const
{
  std::vector<SpacePoint> points;
  if (!m_trkrClusters)
  {
    return points;
  }

  auto physical_layer = [](TrkrDefs::TrkrId trkrid, unsigned int raw_layer) {
    if (trkrid == TrkrDefs::inttId && raw_layer < 3U)
    {
      return raw_layer + 3U;
    }
    return raw_layer;
  };

  auto collect_detector = [&](TrkrDefs::TrkrId trkrid) {
    const auto hitset_keys = m_trkrClusters->getHitSetKeys(trkrid);
    for (const TrkrDefs::hitsetkey hitset_key : hitset_keys)
    {
      const auto range = m_trkrClusters->getClusters(hitset_key);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        SpacePoint point;
        if (!getGlobalClusterPosition(iter->first, iter->second, point))
        {
          continue;
        }
        if (TrkrDefs::getTrkrId(point.key) != trkrid)
        {
          continue;
        }

        const unsigned int raw_layer = point.layer;
        point.layer = physical_layer(trkrid, raw_layer);
        if (trkrid == TrkrDefs::mvtxId && point.layer > 2U)
        {
          continue;
        }
        if (trkrid == TrkrDefs::inttId && (point.layer < 3U || point.layer > 6U))
        {
          continue;
        }
        points.push_back(point);
      }
    }
  };

  collect_detector(TrkrDefs::mvtxId);
  collect_detector(TrkrDefs::inttId);

  std::sort(points.begin(), points.end(), [](const SpacePoint& a, const SpacePoint& b) {
    if (a.r != b.r) { return a.r > b.r; }
    return a.layer > b.layer;
  });
  return points;
}

std::map<unsigned int, std::vector<const Tpc_PolyCluster*>> Full_PolyTrackMatcher::collectTpcClustersByTrack() const
{
  std::map<unsigned int, std::vector<const Tpc_PolyCluster*>> clusters_by_track;
  if (!m_tpcClusters)
  {
    return clusters_by_track;
  }

  for (unsigned int i = 0; i < m_tpcClusters->size(); ++i)
  {
    const Tpc_PolyCluster* cluster = m_tpcClusters->get_cluster(i);
    if (cluster)
    {
      clusters_by_track[cluster->get_source_assembled_track_id()].push_back(cluster);
    }
  }
  return clusters_by_track;
}

const TpcCrossingDecision* Full_PolyTrackMatcher::findCrossingDecision(unsigned int source_assembled_track_id) const
{
  if (!m_crossingDecisions)
  {
    return nullptr;
  }
  for (unsigned int i = 0; i < m_crossingDecisions->size(); ++i)
  {
    const TpcCrossingDecision* decision = m_crossingDecisions->get_decision(i);
    if (decision && decision->get_assembled_track_id() == source_assembled_track_id)
    {
      return decision;
    }
  }
  return nullptr;
}

Full_PolyTrackMatcher::TrajectoryState Full_PolyTrackMatcher::fitTrajectory(const std::vector<SpacePoint>& points) const
{
  TrajectoryState state;
  if (points.size() < 2)
  {
    return state;
  }

  std::vector<Tpc_FittingTools::FitPoint> phi_points;
  std::vector<Tpc_FittingTools::FitPoint> z_points;
  std::vector<Tpc_FittingTools::Point> xyz_points;
  std::vector<double> unwrapped_phi_values;
  phi_points.reserve(points.size());
  z_points.reserve(points.size());
  xyz_points.reserve(points.size());
  unwrapped_phi_values.reserve(points.size());

  double previous_phi = points.front().phi;
  for (const SpacePoint& point : points)
  {
    const double unwrapped_phi = unwrapPhiNear(point.phi, previous_phi);
    phi_points.emplace_back(point.r, unwrapped_phi);
    unwrapped_phi_values.push_back(unwrapped_phi);
    z_points.emplace_back(point.r, point.z);
    xyz_points.push_back({point.x, point.y, point.z});
    previous_phi = unwrapped_phi;
  }

  const Tpc_FittingTools::LineFit phi_line = Tpc_FittingTools::fitLine(phi_points);
  const Tpc_FittingTools::LineFit z_line = Tpc_FittingTools::fitLine(z_points);
  if (!phi_line.ok || !z_line.ok)
  {
    return state;
  }

  state.phi_slope = phi_line.slope;
  state.phi_intercept = phi_line.intercept;
  state.phi_theta = std::atan(phi_line.slope);
  state.phi_bline = phi_line.intercept;
  if (points.size() >= 3)
  {
    Tpc_FittingTools::FitResult circle_fit;
    if (Tpc_FittingTools::fit(xyz_points, circle_fit) && circle_fit.ok && !circle_fit.is_line)
    {
      const double circle_radius = std::hypot(points.back().x - circle_fit.cx, points.back().y - circle_fit.cy);
      if (std::isfinite(circle_fit.cx) && std::isfinite(circle_fit.cy) &&
          std::isfinite(circle_radius) && circle_radius > 0.0)
      {
        state.phi_circle_cx = circle_fit.cx;
        state.phi_circle_cy = circle_fit.cy;
        state.phi_circle_radius = circle_radius;
        state.phi_circle_ok = true;
      }
    }
  }

  if (m_useSagittaPhiFit && points.size() >= 3)
  {
    const Tpc_FittingTools::SagittaFit phi_sagitta = Tpc_FittingTools::fitSagitta(phi_points);
    if (phi_sagitta.ok &&
        std::isfinite(phi_sagitta.S) &&
        std::isfinite(phi_sagitta.x0) &&
        std::isfinite(phi_sagitta.invR) &&
        std::isfinite(phi_sagitta.theta) &&
        std::isfinite(phi_sagitta.b))
    {
      state.phi_S = phi_sagitta.S;
      state.phi_x0 = phi_sagitta.x0;
      state.phi_invR = phi_sagitta.invR;
      state.phi_theta = phi_sagitta.theta;
      state.phi_bline = phi_sagitta.b;
      state.phi_sagitta_ok = true;
    }
  }

  state.reference_r = points.back().r;
  state.reference_phi = previous_phi;
  state.reference_z = points.back().z;
  state.z_slope = z_line.slope;
  state.z_intercept = z_line.intercept;
  state.valid = std::isfinite(state.phi_slope) && std::isfinite(state.phi_intercept) &&
                std::isfinite(state.z_slope) && std::isfinite(state.z_intercept);
  return state;
}

bool Full_PolyTrackMatcher::tpcSeedXyAtZ(const Tpc_PolyTrack& track,
                                             double z,
                                             double arc_direction,
                                             bool use_straight_line,
                                             double& x,
                                             double& y) const
{
  if (track.get_fit_status() == 0 || !std::isfinite(z))
  {
    return false;
  }

  const double x0 = track.get_x();
  const double y0 = track.get_y();
  const double z0 = track.get_z();
  const double px = track.get_px();
  const double py = track.get_py();
  const double pz = track.get_pz();
  const double charge = track.get_charge();
  if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
      !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
      !std::isfinite(charge))
  {
    return false;
  }

  const double dz = z - z0;
  if (use_straight_line || std::fabs(charge * m_magneticFieldT) < 1.0e-12)
  {
    if (std::fabs(pz) < 1.0e-12)
    {
      return false;
    }
    x = x0 + arc_direction * px / pz * dz;
    y = y0 + arc_direction * py / pz * dz;
    return std::isfinite(x) && std::isfinite(y);
  }

  const double pt = std::hypot(px, py);
  if (pt <= 0.0 || std::fabs(pz) < 1.0e-12)
  {
    return false;
  }

  const double signed_radius = pt / (0.003 * charge * m_magneticFieldT);
  const double radius = std::fabs(signed_radius);
  if (!std::isfinite(radius) || radius <= 0.0)
  {
    return false;
  }

  const double tx = px / pt;
  const double ty = py / pt;
  const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
  const double xc = x0 + sign * radius * ty;
  const double yc = y0 - sign * radius * tx;
  const double phi0 = std::atan2(y0 - yc, x0 - xc);
  const double dzds = pz / pt;
  if (std::fabs(dzds) < 1.0e-12)
  {
    return false;
  }

  const double arc = arc_direction * dz / dzds;
  const double phi = phi0 - sign * arc / radius;
  x = xc + radius * std::cos(phi);
  y = yc + radius * std::sin(phi);
  return std::isfinite(x) && std::isfinite(y);
}

double Full_PolyTrackMatcher::tpcSeedResidual2(const Tpc_PolyTrack& track,
                                               const std::vector<const Tpc_PolyCluster*>& tpc_clusters,
                                               double arc_direction,
                                               bool use_straight_line) const
{
  if (tpc_clusters.empty())
  {
    return std::numeric_limits<double>::max();
  }

  double sum = 0.0;
  unsigned int n = 0;
  for (const Tpc_PolyCluster* cluster : tpc_clusters)
  {
    if (!cluster)
    {
      continue;
    }
    const double cx = cluster->get_centroid_x();
    const double cy = cluster->get_centroid_y();
    const double cz = cluster->get_centroid_z();
    if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(cz))
    {
      continue;
    }

    double x = 0.0;
    double y = 0.0;
    if (!tpcSeedXyAtZ(track, cz, arc_direction, use_straight_line, x, y))
    {
      continue;
    }

    const double dx = x - cx;
    const double dy = y - cy;
    sum += dx * dx + dy * dy;
    ++n;
  }

  return n > 0 ? sum / static_cast<double>(n) : std::numeric_limits<double>::max();
}

bool Full_PolyTrackMatcher::initializeTpcTrackState(const Tpc_PolyTrack& track,
                                                    const std::vector<const Tpc_PolyCluster*>& tpc_clusters,
                                                    TrajectoryState& state) const
{
  state = TrajectoryState();

  const double x0 = track.get_x();
  const double y0 = track.get_y();
  const double z0 = track.get_z();
  const double px = track.get_px();
  const double py = track.get_py();
  const double pz = track.get_pz();
  const double charge = track.get_charge();
  const double pt = std::hypot(px, py);
  if (track.get_fit_status() == 0 ||
      !std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
      !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
      !std::isfinite(charge) || !std::isfinite(pt) || pt <= 0.0 ||
      std::fabs(pz) < 1.0e-12)
  {
    return false;
  }

  const bool use_straight_line = std::fabs(charge * m_magneticFieldT) < 1.0e-12;
  const double forward_residual2 = tpcSeedResidual2(track, tpc_clusters, 1.0, use_straight_line);
  const double reverse_residual2 = tpcSeedResidual2(track, tpc_clusters, -1.0, use_straight_line);
  const double arc_direction = forward_residual2 <= reverse_residual2 ? 1.0 : -1.0;

  double x1 = 0.0;
  double y1 = 0.0;
  if (!tpcSeedXyAtZ(track, z0 + arc_direction, arc_direction, use_straight_line, x1, y1))
  {
    return false;
  }

  const double r0 = std::hypot(x0, y0);
  const double r1 = std::hypot(x1, y1);
  const double phi0 = std::atan2(y0, x0);
  const double phi1 = unwrapPhiNear(std::atan2(y1, x1), phi0);
  const double dr = r1 - r0;
  if (!std::isfinite(r0) || !std::isfinite(r1) || r0 <= 0.0 || r1 <= 0.0 ||
      !std::isfinite(phi0) || !std::isfinite(phi1) || std::fabs(dr) <= 1.0e-9)
  {
    return false;
  }

  state.tpc_seed_model = true;
  state.tpc_seed_straight_line = use_straight_line;
  state.tpc_seed_arc_direction = arc_direction;
  state.tpc_seed_x = x0;
  state.tpc_seed_y = y0;
  state.tpc_seed_z = z0;
  state.tpc_seed_px = px;
  state.tpc_seed_py = py;
  state.tpc_seed_pz = pz;
  state.tpc_seed_charge = charge;
  state.reference_r = r0;
  state.reference_phi = phi0;
  state.reference_z = z0;
  state.z_slope = arc_direction * pz / pt;
  state.z_intercept = z0 - state.z_slope * r0;
  state.phi_slope = (phi1 - phi0) / dr;
  state.phi_intercept = phi0 - state.phi_slope * r0;
  state.phi_theta = std::atan(state.phi_slope);
  state.phi_bline = state.phi_intercept;

  if (!use_straight_line)
  {
    const double signed_radius = pt / (0.003 * charge * m_magneticFieldT);
    const double radius = std::fabs(signed_radius);
    const double tx = px / pt;
    const double ty = py / pt;
    const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
    state.phi_circle_cx = x0 + sign * radius * ty;
    state.phi_circle_cy = y0 - sign * radius * tx;
    state.phi_circle_radius = radius;
    state.phi_circle_ok = std::isfinite(radius) && radius > 0.0 &&
                          std::isfinite(state.phi_circle_cx) && std::isfinite(state.phi_circle_cy);
  }

  state.valid = std::isfinite(state.reference_r) && state.reference_r > 0.0 &&
                std::isfinite(state.reference_phi) && std::isfinite(state.reference_z) &&
                std::isfinite(state.z_slope) && std::isfinite(state.z_intercept) &&
                (use_straight_line || state.phi_circle_ok);
  return state.valid;
}

bool Full_PolyTrackMatcher::initializeFixedRadiusState(const Tpc_PolyTrack& track,
                                                       const TrajectoryState& tpc_state,
                                                       const std::vector<SpacePoint>& silicon_points,
                                                       TrajectoryState& state) const
{
  if (!tpc_state.valid || silicon_points.empty())
  {
    return false;
  }

  const unsigned int outer_layer = m_siliconSearchLayers.empty() ? 6U : m_siliconSearchLayers.front();
  double outer_r = 0.0;
  for (const SpacePoint& point : silicon_points)
  {
    if (point.layer == outer_layer && point.r > outer_r)
    {
      outer_r = point.r;
    }
  }
  if (outer_r <= 0.0)
  {
    for (const SpacePoint& point : silicon_points)
    {
      outer_r = std::max(outer_r, point.r);
    }
  }
  if (outer_r <= 0.0)
  {
    return false;
  }

  SpacePoint crossing;
  crossing.layer = outer_layer;
  crossing.r = outer_r;
  if (!predictAtRadius(tpc_state, outer_r, std::hypot(track.get_px(), track.get_py()),
                       crossing.phi, crossing.z, crossing.x, crossing.y))
  {
    return false;
  }

  if (!initializeCircleThroughVertex(track, crossing, state))
  {
    return false;
  }
  state.reference_z = tpc_state.z_intercept + tpc_state.z_slope * state.reference_r;
  state.z_intercept = tpc_state.z_intercept;
  state.z_slope = tpc_state.z_slope;
  return true;
}

bool Full_PolyTrackMatcher::initializeCircleThroughVertex(const Tpc_PolyTrack& track,
                                                          const SpacePoint& crossing,
                                                          TrajectoryState& state) const
{
  const double pt = std::hypot(track.get_px(), track.get_py());
  const double q = track.get_charge();
  double radius = std::numeric_limits<double>::quiet_NaN();
  const double q_over_r = track.get_seed_q_over_r();
  if (std::isfinite(q_over_r) && std::fabs(q_over_r) > 1.0e-12)
  {
    radius = std::fabs(1.0 / q_over_r);
  }
  else if (std::isfinite(pt) && pt > 0.0 && std::isfinite(q) && std::fabs(q) > 0.0 &&
           std::isfinite(m_magneticFieldT) && std::fabs(m_magneticFieldT) > 0.0)
  {
    radius = 100.0 * pt / (0.3 * std::fabs(m_magneticFieldT) * std::fabs(q));
  }
  if (!std::isfinite(radius) || radius <= 0.0 || !std::isfinite(crossing.x) || !std::isfinite(crossing.y))
  {
    return false;
  }

  const double vx = 0.0;
  const double vy = 0.0;
  const double vz = 0.0;
  const double dx = crossing.x - vx;
  const double dy = crossing.y - vy;
  const double d = std::hypot(dx, dy);
  if (d <= 1.0e-9 || d > 2.0 * radius)
  {
    return false;
  }

  const double mx = 0.5 * (vx + crossing.x);
  const double my = 0.5 * (vy + crossing.y);
  const double h = std::sqrt(std::max(0.0, radius * radius - 0.25 * d * d));
  const double px = -dy / d;
  const double py = dx / d;
  const double cx_a = mx + h * px;
  const double cy_a = my + h * py;
  const double cx_b = mx - h * px;
  const double cy_b = my - h * py;

  auto bend_sign = [vx, vy, dx, dy](double cx, double cy)
  {
    const double rvx = vx - cx;
    const double rvy = vy - cy;
    const double tx1 = -rvy;
    const double ty1 = rvx;
    const double dot1 = tx1 * dx + ty1 * dy;
    const double tx = dot1 >= 0.0 ? tx1 : -tx1;
    const double ty = dot1 >= 0.0 ? ty1 : -ty1;
    const double cross = tx * (cy - vy) - ty * (cx - vx);
    return cross >= 0.0 ? 1.0 : -1.0;
  };

  const double wanted_bend = q >= 0.0 ? -1.0 : 1.0;
  const double bend_a = bend_sign(cx_a, cy_a);
  const bool use_a = bend_a == wanted_bend;
  state = TrajectoryState();
  state.phi_circle_cx = use_a ? cx_a : cx_b;
  state.phi_circle_cy = use_a ? cy_a : cy_b;
  state.phi_circle_radius = radius;
  state.phi_circle_ok = true;
  state.fixed_radius_circle = true;
  state.vertex_x = vx;
  state.vertex_y = vy;
  state.vertex_z = vz;
  state.reference_r = crossing.r;
  state.reference_phi = crossing.phi;
  state.reference_z = crossing.z;
  state.z_intercept = vz;
  state.z_slope = crossing.r > 1.0e-9 ? (crossing.z - vz) / crossing.r : 0.0;
  state.valid = true;
  return true;
}

void Full_PolyTrackMatcher::rotateFixedRadiusStateToHit(TrajectoryState& state,
                                                        const SpacePoint& point,
                                                        double pred_phi) const
{
  if (!state.fixed_radius_circle || !state.phi_circle_ok)
  {
    return;
  }

  const double dphi = wrapPhi(pred_phi - point.phi);
  const double c = std::cos(dphi);
  const double s = std::sin(dphi);
  const double cx = state.phi_circle_cx;
  const double cy = state.phi_circle_cy;
  state.phi_circle_cx = c * cx - s * cy;
  state.phi_circle_cy = s * cx + c * cy;
  state.reference_r = point.r;
  state.reference_phi = point.phi;
  state.reference_z = state.z_intercept + state.z_slope * point.r;
}

double Full_PolyTrackMatcher::predictSagittaPhi(double r, const TrajectoryState& state) const
{
  const double c = std::cos(state.phi_theta);
  const double s = std::sin(state.phi_theta);
  double yy = std::tan(state.phi_theta) * r;

  for (unsigned int iter = 0; iter < 25; ++iter)
  {
    const double xrot = c * r + s * yy;
    const double yrot = -s * r + c * yy;
    const double f = Tpc_FittingTools::sagittaModel(xrot, state.phi_S, state.phi_x0, state.phi_invR);
    const double g = yrot - f;
    const double df = sagittaModelDerivative(xrot, state.phi_x0, state.phi_invR);
    const double dg = c - df * s;
    if (std::fabs(dg) < 1.0e-12)
    {
      break;
    }
    const double step = g / dg;
    yy -= step;
    if (std::fabs(step) < 1.0e-10)
    {
      break;
    }
  }

  return state.phi_bline + yy;
}

double Full_PolyTrackMatcher::predictPhiAtRadiusNear(double r, const TrajectoryState& state, double reference_phi) const
{
  if (state.phi_circle_ok)
  {
    const double cx = state.phi_circle_cx;
    const double cy = state.phi_circle_cy;
    const double track_radius = state.phi_circle_radius;
    const double center_radius = std::hypot(cx, cy);
    if (center_radius > 1.0e-4 * std::max(r, track_radius) && track_radius > 0.0)
    {
      const double a = (r * r - track_radius * track_radius + center_radius * center_radius) / (2.0 * center_radius);
      const double h2 = r * r - a * a;
      if (h2 >= -1.0e-9)
      {
        const double h = std::sqrt(std::max(0.0, h2));
        const double ux = cx / center_radius;
        const double uy = cy / center_radius;
        const double bx = a * ux;
        const double by = a * uy;
        const double px = -uy;
        const double py = ux;
        const double phi_a = unwrapPhiNear(std::atan2(by + h * py, bx + h * px), reference_phi);
        const double phi_b = unwrapPhiNear(std::atan2(by - h * py, bx - h * px), reference_phi);
        return std::fabs(phi_a - reference_phi) <= std::fabs(phi_b - reference_phi) ? phi_a : phi_b;
      }
    }
  }

  if (!state.phi_sagitta_ok)
  {
    return unwrapPhiNear(state.phi_intercept + state.phi_slope * r, reference_phi);
  }

  const double c = std::cos(state.phi_theta);
  const double s = std::sin(state.phi_theta);
  double yy = unwrapPhiNear(reference_phi, state.phi_bline) - state.phi_bline;

  for (unsigned int iter = 0; iter < 25; ++iter)
  {
    const double xrot = c * r + s * yy;
    const double yrot = -s * r + c * yy;
    const double f = Tpc_FittingTools::sagittaModel(xrot, state.phi_S, state.phi_x0, state.phi_invR);
    const double g = yrot - f;
    const double df = sagittaModelDerivative(xrot, state.phi_x0, state.phi_invR);
    const double dg = c - df * s;
    if (std::fabs(dg) < 1.0e-12)
    {
      break;
    }
    const double step = g / dg;
    yy -= step;
    if (std::fabs(step) < 1.0e-10)
    {
      break;
    }
  }

  return unwrapPhiNear(state.phi_bline + yy, reference_phi);
}

double Full_PolyTrackMatcher::projectionStepSize(double pt) const
{
  const double min_step = std::max(m_projectionMinStepCm, 1.0e-3);
  const double max_step = std::max(m_projectionMaxStepCm, min_step);
  const double safe_pt = std::isfinite(pt) ? std::max(pt, 0.25) : 0.25;
  return std::min(max_step, std::max(min_step, max_step * safe_pt));
}

double Full_PolyTrackMatcher::pointTheta(const SpacePoint& point) const
{
  return std::atan2(point.r, point.z);
}


bool Full_PolyTrackMatcher::chainHasRequiredSiliconLayers(const Chain& chain) const
{
  for (const unsigned int required_layer : m_requiredSiliconLayers)
  {
    bool found = false;
    for (const ChainHit& hit : chain.hits)
    {
      if (hit.point.layer == required_layer)
      {
        found = true;
        break;
      }
    }
    if (!found)
    {
      return false;
    }
  }
  return true;
}

Full_PolyTrackMatcher::ResidualMatch Full_PolyTrackMatcher::myEventResidualMatch(unsigned int layer, const Chain& chain,
                                                                                double dphi, double dtheta, double pred_theta) const
{
  ResidualMatch out;

  if (std::fabs(dphi) > m_myEventAngularPrewindow ||
      std::fabs(dtheta) > m_myEventAngularPrewindow)
  {
    return out;
  }

  const double previous_dphi = chain.has_previous_residual ? chain.previous_dphi : 0.0;
  const double previous_dtheta = chain.has_previous_residual ? chain.previous_dtheta : 0.0;
  const double mean_phi = dynamicMeanPhi(layer, previous_dphi, chain.has_previous_residual);
  const double sigma_phi = dynamicSigmaPhi(chain.pt);
  const double mean_theta = dynamicMeanTheta(layer, previous_dtheta, chain.has_previous_residual);
  const double sigma_theta = dynamicSigmaTheta(chain.pt, pred_theta);
  if (sigma_phi <= 0.0 || sigma_theta <= 0.0)
  {
    return out;
  }

  out.sdphi = (dphi - mean_phi) / sigma_phi;
  out.sdtheta = (dtheta - mean_theta) / sigma_theta;

  double sigma_veto = m_phiWindowSigma;
  double sigma_inner = m_phiWindowSigma;
  if (layer == 0U)
  {
    sigma_veto = m_phiWindowSigmaInnerVeto;
    sigma_inner = m_phiWindowSigmaInner;
  }
  else if (layer == 1U)
  {
    sigma_veto = m_phiWindowSigmaSecond;
    sigma_inner = m_phiWindowSigmaSecond;
  }

  const double q = chain.charge < 0.0 ? -1.0 : 1.0;
  out.pass = out.sdphi * q > -sigma_veto &&
             out.sdphi * q < sigma_inner &&
             std::fabs(out.sdtheta) < m_thetaWindowSigma;
  out.chi2 = square(out.sdphi) + square(out.sdtheta);
  return out;
}

double Full_PolyTrackMatcher::slopeEta(double dzdr) const
{
  return std::asinh(dzdr);
}

double Full_PolyTrackMatcher::dynamicEtaWindow(double pt, double charge) const
{
  const double safe_pt = std::max(pt, 0.25);
  if (charge < 0.0)
  {
    return 0.045 + 0.0031 * std::exp(1.0 / safe_pt);
  }
  return 0.050 + 0.0064 * std::exp(1.1 / safe_pt);
}

bool Full_PolyTrackMatcher::siliconSlopeMatchesTpc(const Chain& chain, const SpacePoint& point) const
{
  if (chain.hits.empty())
  {
    return true;
  }

  const SpacePoint& previous = chain.hits.back().point;
  const double dr = point.r - previous.r;
  if (std::fabs(dr) < 1.0e-9)
  {
    return true;
  }

  const double silicon_eta = slopeEta((point.z - previous.z) / dr);
  const double tpc_eta = slopeEta(chain.state.z_slope);
  return std::fabs(silicon_eta - tpc_eta) < m_siliconSearchWindowFactor * dynamicEtaWindow(chain.pt, chain.charge);
}

double Full_PolyTrackMatcher::dynamicMeanPhi(unsigned int layer, double previous_dphi, bool has_previous) const
{
  if (!m_useDynamicResiduals || !has_previous || layer >= m_dynamicPhiMeanOffset.size())
  {
    return 0.0;
  }
  return m_dynamicPhiMeanOffset[layer] + m_dynamicPhiMeanSlope[layer] * previous_dphi;
}

double Full_PolyTrackMatcher::dynamicSigmaPhi(double pt) const
{
  const double safe_pt = std::max(pt, 0.25);
  const double dphi_window = std::max(0.05 + 0.0064 * std::exp(1.1 / safe_pt), 0.05);
  return std::max(m_siliconSearchWindowFactor * dphi_window / std::max(m_phiWindowSigma, 1.0e-9), 1.0e-9);
}

double Full_PolyTrackMatcher::dynamicMeanTheta(unsigned int layer, double previous_dtheta, bool has_previous) const
{
  if (!m_useDynamicResiduals || !has_previous || layer >= m_dynamicThetaMeanOffset.size())
  {
    return 0.0;
  }
  return m_dynamicThetaMeanOffset[layer] + m_dynamicThetaMeanSlope[layer] * previous_dtheta;
}

double Full_PolyTrackMatcher::dynamicSigmaTheta(double pt, double pred_theta) const
{
  const double safe_pt = std::max(pt, 0.25);
  const double deta_window = std::max(-0.014 + 0.0331 * std::exp(0.48 / safe_pt), 0.03);
  const double theta_window = deta_window * std::max(std::sin(pred_theta), 1.0e-9);
  return std::max(m_siliconSearchWindowFactor * theta_window / std::max(m_thetaWindowSigma, 1.0e-9), 1.0e-9);
}

double Full_PolyTrackMatcher::dynamicDzWindow(double pt) const
{
  const double safe_pt = std::max(pt, 0.25);
  return m_siliconSearchWindowFactor * (1.138 + 0.3919 * std::exp(0.84 / safe_pt));
}

bool Full_PolyTrackMatcher::predictTpcSeedAtRadius(const TrajectoryState& state,
                                                      double r,
                                                      double& pred_phi,
                                                      double& pred_z,
                                                      double& pred_x,
                                                      double& pred_y) const
{
  if (!state.tpc_seed_model || !state.valid || !std::isfinite(r) || r <= 0.0)
  {
    return false;
  }

  if (state.tpc_seed_straight_line || std::fabs(state.tpc_seed_charge * m_magneticFieldT) < 1.0e-12)
  {
    if (std::fabs(state.tpc_seed_pz) < 1.0e-12)
    {
      return false;
    }

    const double ax = state.tpc_seed_arc_direction * state.tpc_seed_px / state.tpc_seed_pz;
    const double ay = state.tpc_seed_arc_direction * state.tpc_seed_py / state.tpc_seed_pz;
    const double a = ax * ax + ay * ay;
    const double b = 2.0 * (state.tpc_seed_x * ax + state.tpc_seed_y * ay);
    const double c = state.tpc_seed_x * state.tpc_seed_x + state.tpc_seed_y * state.tpc_seed_y - r * r;
    if (a <= 1.0e-20)
    {
      return false;
    }

    const double disc = b * b - 4.0 * a * c;
    if (disc < -1.0e-9)
    {
      return false;
    }
    const double root = std::sqrt(std::max(0.0, disc));
    const double dz_a = (-b + root) / (2.0 * a);
    const double dz_b = (-b - root) / (2.0 * a);
    const double z_a = state.tpc_seed_z + dz_a;
    const double z_b = state.tpc_seed_z + dz_b;
    pred_z = std::fabs(z_a - state.reference_z) <= std::fabs(z_b - state.reference_z) ? z_a : z_b;
    pred_x = state.tpc_seed_x + ax * (pred_z - state.tpc_seed_z);
    pred_y = state.tpc_seed_y + ay * (pred_z - state.tpc_seed_z);
    pred_phi = wrapPhi(std::atan2(pred_y, pred_x));
    return std::isfinite(pred_phi) && std::isfinite(pred_z) && std::isfinite(pred_x) && std::isfinite(pred_y);
  }

  const double pt = std::hypot(state.tpc_seed_px, state.tpc_seed_py);
  if (pt <= 0.0 || std::fabs(state.tpc_seed_pz) < 1.0e-12 ||
      !state.phi_circle_ok || state.phi_circle_radius <= 0.0)
  {
    return false;
  }

  const double cx = state.phi_circle_cx;
  const double cy = state.phi_circle_cy;
  const double track_radius = state.phi_circle_radius;
  const double center_radius = std::hypot(cx, cy);
  if (center_radius <= 1.0e-9)
  {
    return false;
  }

  const double a = (r * r - track_radius * track_radius + center_radius * center_radius) / (2.0 * center_radius);
  const double h2 = r * r - a * a;
  if (h2 < -1.0e-9)
  {
    return false;
  }

  const double h = std::sqrt(std::max(0.0, h2));
  const double ux = cx / center_radius;
  const double uy = cy / center_radius;
  const double bx = a * ux;
  const double by = a * uy;
  const double px = -uy;
  const double py = ux;

  const double signed_radius = pt / (0.003 * state.tpc_seed_charge * m_magneticFieldT);
  const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
  const double phi0 = std::atan2(state.tpc_seed_y - cy, state.tpc_seed_x - cx);
  const double dzds = state.tpc_seed_pz / pt;

  auto candidate = [&](double x, double y, double& z, double& distance) {
    const double angle = std::atan2(y - cy, x - cx);
    double best_z = std::numeric_limits<double>::quiet_NaN();
    double best_distance = std::numeric_limits<double>::max();
    for (const double turn : {-1.0, 0.0, 1.0})
    {
      const double delta = unwrapPhiNear(angle, phi0) - phi0 + turn * 2.0 * M_PI;
      const double dz = -sign * state.tpc_seed_arc_direction * delta * dzds * track_radius;
      const double z_try = state.tpc_seed_z + dz;
      const double distance_try = std::fabs(z_try - state.reference_z);
      if (std::isfinite(z_try) && distance_try < best_distance)
      {
        best_z = z_try;
        best_distance = distance_try;
      }
    }
    z = best_z;
    distance = best_distance;
  };

  const double x_a = bx + h * px;
  const double y_a = by + h * py;
  const double x_b = bx - h * px;
  const double y_b = by - h * py;
  double z_a = 0.0;
  double z_b = 0.0;
  double distance_a = 0.0;
  double distance_b = 0.0;
  candidate(x_a, y_a, z_a, distance_a);
  candidate(x_b, y_b, z_b, distance_b);

  if (distance_a <= distance_b)
  {
    pred_x = x_a;
    pred_y = y_a;
    pred_z = z_a;
  }
  else
  {
    pred_x = x_b;
    pred_y = y_b;
    pred_z = z_b;
  }
  pred_phi = wrapPhi(std::atan2(pred_y, pred_x));
  return std::isfinite(pred_phi) && std::isfinite(pred_z) && std::isfinite(pred_x) && std::isfinite(pred_y);
}

bool Full_PolyTrackMatcher::predictAtRadius(const TrajectoryState& state, double r, double pt,
                                            double& pred_phi, double& pred_z,
                                            double& pred_x, double& pred_y) const
{
  if (state.tpc_seed_model)
  {
    return predictTpcSeedAtRadius(state, r, pred_phi, pred_z, pred_x, pred_y);
  }

  if (!state.valid || !std::isfinite(r) || r <= 0.0 || !std::isfinite(state.reference_r) || state.reference_r <= 0.0)
  {
    return false;
  }

  const double step_size = projectionStepSize(pt);
  const double direction = r >= state.reference_r ? 1.0 : -1.0;
  const unsigned int n_steps = std::max(1U, static_cast<unsigned int>(std::ceil(std::fabs(r - state.reference_r) / step_size)));
  double current_phi = state.reference_phi;

  for (unsigned int step = 1; step <= n_steps; ++step)
  {
    const double fraction = static_cast<double>(step) / static_cast<double>(n_steps);
    const double step_r = state.reference_r + direction * std::fabs(r - state.reference_r) * fraction;
    const double next_phi = predictPhiAtRadiusNear(step_r, state, current_phi);
    if (!std::isfinite(next_phi))
    {
      return false;
    }
    if (Verbosity() > 3)
    {
      std::cout << Name() << "::predictAtRadius - step=" << step
                << "/" << n_steps
                << " r=" << step_r
                << " phi=" << wrapPhi(next_phi)
                << " z=" << state.z_intercept + state.z_slope * step_r
                << " circle=" << state.phi_circle_ok
                << " sagitta=" << state.phi_sagitta_ok << std::endl;
    }
    current_phi = next_phi;
  }

  pred_phi = wrapPhi(current_phi);
  pred_z = state.z_intercept + state.z_slope * r;
  pred_x = r * std::cos(pred_phi);
  pred_y = r * std::sin(pred_phi);
  return std::isfinite(pred_phi) && std::isfinite(pred_z) && std::isfinite(pred_x) && std::isfinite(pred_y);
}

Full_PolyTrackMatcher::Chain Full_PolyTrackMatcher::extendWithHit(const Chain& chain, const SpacePoint& point,
                                                                  double pred_phi, double pred_z, double pred_theta,
                                                                  double pred_x, double pred_y,
                                                                  double dphi, double dtheta, double rdphi, double dz, double chi2) const
{
  Chain out = chain;
  ChainHit hit;
  hit.point = point;
  hit.pred_phi = pred_phi;
  hit.pred_z = pred_z;
  hit.pred_theta = pred_theta;
  hit.pred_x = pred_x;
  hit.pred_y = pred_y;
  hit.dphi = dphi;
  hit.dtheta = dtheta;
  hit.rdphi = rdphi;
  hit.dz = dz;
  hit.chi2 = chi2;
  out.hits.push_back(hit);
  out.chi2 += chi2;
  out.ndf += 2.0;
  out.previous_dphi = dphi;
  out.previous_dtheta = dtheta;
  out.has_previous_residual = true;
  out.score = out.chi2 + m_missingLayerPenalty * out.n_missing;
  return out;
}

Full_PolyTrackMatcher::Chain Full_PolyTrackMatcher::extendMissing(const Chain& chain, unsigned int layer) const
{
  Chain out = chain;
  out.missing_mask |= layerBit(layer);
  out.n_missing += 1U;
  out.score = out.chi2 + m_missingLayerPenalty * out.n_missing;
  return out;
}

void Full_PolyTrackMatcher::fillQaDphiCorrelation(const Tpc_PolyTrack& track, const Chain& chain,
                                                        const SpacePoint& point, double current_dphi)
{
  if (!m_writeQA || !m_qaFile || chain.hits.empty())
  {
    return;
  }
  const ChainHit& previous_hit = chain.hits.back();
  if (previous_hit.point.r <= 0.0 || !std::isfinite(current_dphi))
  {
    return;
  }

  const double previous_dphi = previous_hit.rdphi / previous_hit.point.r;
  const double pt = std::hypot(track.get_px(), track.get_py());
  if (!std::isfinite(previous_dphi) || !std::isfinite(pt))
  {
    return;
  }

  const unsigned int side = point.z >= 0.0 ? 1U : 0U;
  const int charge_bin = track.get_charge() < 0.0 ? -1 : 1;
  double phi = wrapPhi(point.phi);
  if (phi < 0.0)
  {
    phi += 2.0 * M_PI;
  }
  unsigned int phi_bin = static_cast<unsigned int>(m_qaPhiBins * phi / (2.0 * M_PI));
  if (phi_bin >= m_qaPhiBins)
  {
    phi_bin = m_qaPhiBins - 1U;
  }

  const double phi_bin_min = 2.0 * M_PI * static_cast<double>(phi_bin) / static_cast<double>(m_qaPhiBins);
  const double phi_bin_max = 2.0 * M_PI * static_cast<double>(phi_bin + 1U) / static_cast<double>(m_qaPhiBins);

  std::ostringstream key;
  key << "h_dphi_curr_prev_pt_l" << point.layer
      << "_side" << side
      << "_q" << (charge_bin < 0 ? "neg" : "pos")
      << "_phibin" << phi_bin;

  TH3F* hist = nullptr;
  auto found = m_hDphiCurrentVsPreviousVsPt.find(key.str());
  if (found != m_hDphiCurrentVsPreviousVsPt.end())
  {
    hist = found->second;
  }
  else
  {
    std::ostringstream title;
    title << "current dphi vs previous dphi vs pT layer " << point.layer
          << " side " << side
          << " q " << (charge_bin < 0 ? "-" : "+")
          << " phi bin " << phi_bin << " [" << phi_bin_min << ", " << phi_bin_max << ") rad"
          << ";#Delta#phi current layer [rad];#Delta#phi previous accepted layer [rad];p_{T} [GeV]";
    hist = new TH3F(key.str().c_str(), title.str().c_str(),
                    static_cast<int>(m_qaDphiBins), m_qaDphiMin, m_qaDphiMax,
                    static_cast<int>(m_qaDphiBins), m_qaDphiMin, m_qaDphiMax,
                    static_cast<int>(m_qaPtBins), m_qaPtMin, m_qaPtMax);
    hist->SetDirectory(nullptr);
    m_hDphiCurrentVsPreviousVsPt[key.str()] = hist;
  }

  hist->Fill(current_dphi, previous_dphi, pt);
}

std::vector<Full_PolyTrackMatcher::Chain> Full_PolyTrackMatcher::buildChains(
    const Tpc_PolyTrack& track,
    const std::vector<const Tpc_PolyCluster*>& tpc_clusters,
    const std::vector<SpacePoint>& silicon_points)
{
  Chain seed;
  seed.pt = std::hypot(track.get_px(), track.get_py());
  seed.charge = track.get_charge();
  if (tpc_clusters.size() < m_minTpcClustersForAssociation || seed.pt < m_minPtForAssociation)
  {
    if (Verbosity() > 1)
    {
      std::cout << Name() << "::buildChains - skip silicon association for source_track="
                << track.get_source_assembled_track_id()
                << " n_tpc_clusters=" << tpc_clusters.size()
                << " pt=" << seed.pt << std::endl;
    }
    return {};
  }
  const double tpc_pca_r = std::hypot(track.get_x(), track.get_y());
  if (std::isfinite(tpc_pca_r) && tpc_pca_r > m_maxTpcPcaRadiusForAssociation)
  {
    if (Verbosity() > 1)
    {
      std::cout << Name() << "::buildChains - skip silicon association for source_track="
                << track.get_source_assembled_track_id()
                << " tpc_pca_r=" << tpc_pca_r << std::endl;
    }
    return {};
  }
  if (!initializeTpcTrackState(track, tpc_clusters, seed.state))
  {
    if (Verbosity() > 1)
    {
      std::cout << Name() << "::buildChains - skip silicon association for source_track="
                << track.get_source_assembled_track_id()
                << " invalid TPC poly track parameters" << std::endl;
    }
    return {};
  }
  if (Verbosity() > 1)
  {
    std::cout << Name() << "::buildChains - source_track=" << track.get_source_assembled_track_id()
              << " pt=" << seed.pt
              << " seed_source=tpc_poly_track"
              << " circle=" << seed.state.phi_circle_ok
              << " sagitta=" << seed.state.phi_sagitta_ok
              << " fixed_tpc_seed=1"
              << " circle_center_r=" << std::hypot(seed.state.phi_circle_cx, seed.state.phi_circle_cy)
              << " circle_track_r=" << seed.state.phi_circle_radius
              << " reference_r=" << seed.state.reference_r
              << " reference_phi=" << wrapPhi(seed.state.reference_phi)
              << " reference_z=" << seed.state.reference_z << std::endl;
  }


  const std::vector<unsigned int> ordered_layers{6, 5, 4, 3, 2, 1, 0};
  if (Verbosity() > 1)
  {
    std::cout << "Silicon cluster ordering for track " << track.get_track_id() << std::endl;
    for (unsigned int layer : ordered_layers)
    {
      std::cout << " layer " << layer << ":";
      for (const SpacePoint& point : silicon_points)
      {
        if (point.layer != layer)
        {
          continue;
        }
        std::cout << " r=" << point.r << "("
                  << (TrkrDefs::getTrkrId(point.key) == TrkrDefs::mvtxId ? "MVTX" : "INTT")
                  << ")";
      }
      std::cout << std::endl;
    }
  }

  std::vector<unsigned int> active_layers;
  bool found_seed_window_start = false;
  for (unsigned int layer : ordered_layers)
  {
    bool layer_has_seed_window_candidate = false;
    for (const SpacePoint& point : silicon_points)
    {
      if (point.layer != layer)
      {
        continue;
      }

      double pred_phi = 0.0;
      double pred_z = 0.0;
      double pred_x = 0.0;
      double pred_y = 0.0;
      if (!predictAtRadius(seed.state, point.r, seed.pt, pred_phi, pred_z, pred_x, pred_y))
      {
        continue;
      }

      const double pred_theta = std::atan2(point.r, pred_z);
      const double dphi = wrapPhi(pred_phi - point.phi);
      const double dtheta = pred_theta - pointTheta(point);
      const ResidualMatch residual = myEventResidualMatch(layer, seed, dphi, dtheta, pred_theta);
      if (residual.pass)
      {
        layer_has_seed_window_candidate = true;
        break;
      }
    }

    if (!found_seed_window_start && !layer_has_seed_window_candidate)
    {
      continue;
    }
    found_seed_window_start = true;
    active_layers.push_back(layer);
  }

  if (active_layers.empty())
  {
    return {seed};
  }

  std::vector<Chain> chains{seed};
  for (unsigned int layer : active_layers)
  {
    std::vector<Chain> next_chains;
    for (const Chain& chain : chains)
    {
      std::vector<Chain> layer_extensions;
      for (const SpacePoint& point : silicon_points)
      {
        if (point.layer != layer)
        {
          continue;
        }

        double pred_phi = 0.0;
        double pred_z = 0.0;
        double pred_x = 0.0;
        double pred_y = 0.0;
        if (!predictAtRadius(chain.state, point.r, chain.pt, pred_phi, pred_z, pred_x, pred_y))
        {
          continue;
        }

        const double pred_theta = std::atan2(point.r, pred_z);
        const double dphi = wrapPhi(pred_phi - point.phi);
        const double dtheta = pred_theta - pointTheta(point);
        const double dz = point.z - pred_z;
        fillQaDphiCorrelation(track, chain, point, dphi);
        const ResidualMatch residual = myEventResidualMatch(layer, chain, dphi, dtheta, pred_theta);
        if (!residual.pass)
        {
          continue;
        }

        const double rdphi = point.r * dphi;
        const double chi2 = residual.chi2;
        layer_extensions.push_back(extendWithHit(chain, point, pred_phi, pred_z, pred_theta,
                                                 pred_x, pred_y, dphi, dtheta, rdphi, dz, chi2));
      }

      std::sort(layer_extensions.begin(), layer_extensions.end(), [](const Chain& a, const Chain& b) {
        if (a.hits.size() != b.hits.size()) { return a.hits.size() > b.hits.size(); }
        return a.score < b.score;
      });
      if (layer_extensions.size() > m_maxBranchesPerLayer)
      {
        layer_extensions.resize(m_maxBranchesPerLayer);
      }
      next_chains.insert(next_chains.end(), layer_extensions.begin(), layer_extensions.end());
      next_chains.push_back(extendMissing(chain, layer));
    }

    std::sort(next_chains.begin(), next_chains.end(), [](const Chain& a, const Chain& b) {
      if (a.hits.size() != b.hits.size()) { return a.hits.size() > b.hits.size(); }
      return a.score < b.score;
    });
    if (next_chains.size() > m_maxChains)
    {
      next_chains.resize(m_maxChains);
    }
    chains.swap(next_chains);
  }

  return chains;
}

const Full_PolyTrackMatcher::Chain* Full_PolyTrackMatcher::selectBestChain(const std::vector<Chain>& chains) const
{
  const Chain* best = nullptr;
  for (const Chain& chain : chains)
  {
    if (chain.hits.size() < m_minSiliconClusters)
    {
      continue;
    }
    if (!chainHasRequiredSiliconLayers(chain))
    {
      continue;
    }
    if (!best ||
        chain.hits.size() > best->hits.size() ||
        (chain.hits.size() == best->hits.size() && chain.score < best->score))
    {
      best = &chain;
    }
  }
  return best;
}

Full_PolyTrackMatcher::Chain Full_PolyTrackMatcher::attachClosestInttClusters(
    const Chain& mvtx_chain,
    const std::vector<SpacePoint>& silicon_points)
{
  Chain out = mvtx_chain;
  if (out.hits.empty() || !out.state.valid)
  {
    return out;
  }

  for (const unsigned int layer : m_inttMatchLayers)
  {
    const SpacePoint* best_point = nullptr;
    double best_pred_phi = 0.0;
    double best_pred_z = 0.0;
    double best_pred_theta = 0.0;
    double best_pred_x = 0.0;
    double best_pred_y = 0.0;
    double best_dphi = 0.0;
    double best_dtheta = 0.0;
    double best_rdphi = 0.0;
    double best_dz = 0.0;
    double best_chi2 = std::numeric_limits<double>::max();

    for (const SpacePoint& point : silicon_points)
    {
      if (point.layer != layer || static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(point.key)) != TrkrDefs::inttId)
      {
        continue;
      }

      double pred_phi = 0.0;
      double pred_z = 0.0;
      double pred_x = 0.0;
      double pred_y = 0.0;
      if (!predictAtRadius(out.state, point.r, out.pt, pred_phi, pred_z, pred_x, pred_y))
      {
        continue;
      }

      const double pred_theta = std::atan2(point.r, pred_z);
      const double dphi = wrapPhi(pred_phi - point.phi);
      const double dtheta = pred_theta - pointTheta(point);
      const double dz = point.z - pred_z;
      const ResidualMatch residual = myEventResidualMatch(layer, out, dphi, dtheta, pred_theta);
      if (!residual.pass)
      {
        continue;
      }

      const double rdphi = point.r * dphi;
      const double chi2 = residual.chi2;
      if (chi2 < best_chi2)
      {
        best_point = &point;
        best_pred_phi = pred_phi;
        best_pred_z = pred_z;
        best_pred_theta = pred_theta;
        best_pred_x = pred_x;
        best_pred_y = pred_y;
        best_dphi = dphi;
        best_dtheta = dtheta;
        best_rdphi = rdphi;
        best_dz = dz;
        best_chi2 = chi2;
      }
    }

    if (best_point)
    {
      out = extendWithHit(out, *best_point, best_pred_phi, best_pred_z, best_pred_theta,
                          best_pred_x, best_pred_y, best_dphi, best_dtheta, best_rdphi, best_dz, best_chi2);
    }
    else
    {
      out = extendMissing(out, layer);
    }
  }

  return out;
}

void Full_PolyTrackMatcher::fillTrack(const Tpc_PolyTrack& tpc_track, const Chain& chain)
{
  Full_PolyTrackv1* out = new Full_PolyTrackv1();
  out->set_event(m_event);
  out->set_track_id(m_fullTracks->size());
  out->set_tpc_poly_track_id(tpc_track.get_track_id());
  out->set_source_assembled_track_id(tpc_track.get_source_assembled_track_id());
  out->set_n_tpc_clusters(tpc_track.get_nclusters());
  out->set_n_missing_layers(chain.n_missing);
  out->set_missing_layer_mask(chain.missing_mask);
  out->set_fit_status((chain.state.valid || tpc_track.get_fit_status() != 0) ? 1 : 0);
  out->set_chi2(chain.chi2);
  out->set_ndf(chain.ndf);
  out->set_score(chain.score);
  out->set_charge(tpc_track.get_charge());
  out->set_seed_x(tpc_track.get_x());
  out->set_seed_y(tpc_track.get_y());
  out->set_seed_z(tpc_track.get_z());
  out->set_seed_px(tpc_track.get_px());
  out->set_seed_py(tpc_track.get_py());
  out->set_seed_pz(tpc_track.get_pz());

  if (chain.state.valid && !chain.hits.empty())
  {
    const double r0 = chain.hits.back().point.r;
    double phi = 0.0;
    double z = 0.0;
    double x = 0.0;
    double y = 0.0;
    predictAtRadius(chain.state, r0, chain.pt, phi, z, x, y);
    out->set_x(x);
    out->set_y(y);
    out->set_z(z);
    const double next_r = r0 + 1.0;
    double phi2 = 0.0;
    double z2 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    predictAtRadius(chain.state, next_r, chain.pt, phi2, z2, x2, y2);
    out->set_px(x2 - x);
    out->set_py(y2 - y);
    out->set_pz(z2 - z);
  }
  else
  {
    out->set_x(tpc_track.get_x());
    out->set_y(tpc_track.get_y());
    out->set_z(tpc_track.get_z());
    out->set_px(tpc_track.get_px());
    out->set_py(tpc_track.get_py());
    out->set_pz(tpc_track.get_pz());
  }

  std::vector<const ChainHit*> ordered_hits;
  ordered_hits.reserve(chain.hits.size());
  for (const ChainHit& hit : chain.hits)
  {
    ordered_hits.push_back(&hit);
  }
  std::sort(ordered_hits.begin(), ordered_hits.end(), [](const ChainHit* a, const ChainHit* b) {
    if (a->point.r != b->point.r) { return a->point.r < b->point.r; }
    return a->point.layer < b->point.layer;
  });

  if (Verbosity() > 1)
  {
    std::cout << "ASSOCIATED:";
    for (const ChainHit& hit : chain.hits)
    {
      std::cout << " [L" << hit.point.layer << " r=" << hit.point.r << "]";
    }
    std::cout << std::endl;

    std::cout << "STORED:";
    for (const ChainHit* hit : ordered_hits)
    {
      std::cout << " [L" << hit->point.layer << " r=" << hit->point.r << "]";
    }
    std::cout << std::endl;
  }

  for (const ChainHit* hit : ordered_hits)
  {
    out->add_cluster_key(hit->point.key);
    out->add_silicon_state(hit->point.layer, hit->point.key,
                           finite_float(hit->point.x), finite_float(hit->point.y), finite_float(hit->point.z),
                           finite_float(hit->pred_x), finite_float(hit->pred_y), finite_float(hit->pred_z),
                           finite_float(hit->rdphi), finite_float(hit->dz), finite_float(hit->chi2));
  }

  m_fullTracks->add_track(out);
}

int Full_PolyTrackMatcher::process_event(PHCompositeNode* topNode)
{
  if (!m_tpcTracks || !m_tpcClusters || !m_trkrClusters || !m_actsGeometry || !m_fullTracks)
  {
    if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  m_fullTracks->Reset();
  const std::vector<SpacePoint> silicon_points = collectSiliconClusters();
  const auto tpc_clusters_by_track = collectTpcClustersByTrack();

  for (unsigned int i = 0; i < m_tpcTracks->size(); ++i)
  {
    const Tpc_PolyTrack* tpc_track = m_tpcTracks->get_track(i);
    if (!tpc_track)
    {
      continue;
    }
    auto cluster_iter = tpc_clusters_by_track.find(tpc_track->get_source_assembled_track_id());
    const std::vector<const Tpc_PolyCluster*> empty_clusters;
    const std::vector<const Tpc_PolyCluster*>& tpc_clusters =
        cluster_iter != tpc_clusters_by_track.end() ? cluster_iter->second : empty_clusters;
    std::vector<Chain> chains = buildChains(*tpc_track, tpc_clusters, silicon_points);
    const Chain* best = selectBestChain(chains);
    Chain unmatched_chain;
    const Chain* selected_chain = best;
    if (!selected_chain)
    {
      selected_chain = chains.empty() ? &unmatched_chain : &chains.front();
    }

    fillTrack(*tpc_track, *selected_chain);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " tpc_tracks=" << m_tpcTracks->size()
              << " silicon_clusters=" << silicon_points.size()
              << " full_tracks=" << m_fullTracks->size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackMatcher::End(PHCompositeNode*)
{
  if (m_qaFile)
  {
    m_qaFile->cd();
    for (auto& item : m_hDphiCurrentVsPreviousVsPt)
    {
      if (item.second)
      {
        item.second->Write();
      }
    }
    m_qaFile->Close();
    delete m_qaFile;
    m_qaFile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
