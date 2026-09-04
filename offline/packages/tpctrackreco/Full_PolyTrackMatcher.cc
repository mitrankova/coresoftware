#include "Full_PolyTrackMatcher.h"

#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrackv1.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
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
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase/TrkrDefs.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TAxis.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>

namespace
{
  constexpr double SiliconBeamXIntercept = -0.0407;
  constexpr double SiliconBeamXSlope = -0.0015;
  constexpr double SiliconBeamYIntercept = 0.1645;
  constexpr double SiliconBeamYSlope = -0.0001;
  constexpr double TpcBeamXIntercept = -0.0103;
  constexpr double TpcBeamXSlope = -0.0013;
  constexpr double TpcBeamYIntercept = 0.1814;
  constexpr double TpcBeamYSlope = -0.0002;

  float finite_float(double value)
  {
    return std::isfinite(value) ? static_cast<float>(value) : std::numeric_limits<float>::quiet_NaN();
  }

  double square(double value)
  {
    return value * value;
  }
  constexpr double SeedPhiWindow = 0.15;
  constexpr double SeedThetaPreWindow = 4;

  double seedEtaWindow(double pt)
  {
    return std::max(0.03, -0.014 + 0.0331 * std::exp(0.48 / pt));
  }

  double seedThetaWindow(double pt, double theta)
  {
    return seedEtaWindow(pt) * std::max(std::sin(theta), 1.0e-9);
  }

  double thetaToEta(double theta)
  {
    if (!std::isfinite(theta))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double half_tan = std::tan(0.5 * theta);
    if (!std::isfinite(half_tan) || half_tan <= 0.0)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return -std::log(half_tan);
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

  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  if (!m_crossingDecisions && Verbosity() > 0)
  {
    std::cout << Name() << "::getNodes - optional " << m_crossingDecisionNodeName << " not found" << std::endl;
  }

  m_clusterCrossingAssoc = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, m_clusterCrossingAssocNodeName);
  if (!m_clusterCrossingAssoc && Verbosity() > 0)
  {
    std::cout << Name() << "::getNodes - optional " << m_clusterCrossingAssocNodeName << " not found" << std::endl;
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

  for (const TrkrDefs::TrkrId trkrid : {TrkrDefs::mvtxId, TrkrDefs::inttId})
  {
    for (const unsigned int layer : m_siliconSearchLayers)
    {
      const auto hitset_keys = m_trkrClusters->getHitSetKeys(trkrid, static_cast<uint8_t>(layer));
      for (const TrkrDefs::hitsetkey hitset_key : hitset_keys)
      {
        const auto range = m_trkrClusters->getClusters(hitset_key);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          SpacePoint point;
          if (getGlobalClusterPosition(iter->first, iter->second, point))
          {
            points.push_back(point);
          }
        }
      }
    }
  }
  std::sort(points.begin(), points.end(), [](const SpacePoint& a, const SpacePoint& b) {
    if (a.layer != b.layer) { return a.layer > b.layer; }
    return a.r > b.r;
  });
  return points;
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

Full_PolyTrackMatcher::TrajectoryState Full_PolyTrackMatcher::makeTpcReferenceTrajectory(const Tpc_PolyTrack& track) const
{
  TrajectoryState state = makeTpcSeedTrajectory(track);
  if (!state.valid || !m_trkrClusters)
  {
    return state;
  }

  std::vector<SpacePoint> points;
  points.reserve(track.size_cluster_keys() + 1U);
  for (unsigned int i = 0; i < track.size_cluster_keys(); ++i)
  {
    const TrkrDefs::cluskey key = track.get_cluster_key(i);
    TrkrCluster* cluster = m_trkrClusters->findCluster(key);
    SpacePoint point;
    if (cluster && getGlobalClusterPosition(key, cluster, point))
    {
      points.push_back(point);
    }
  }

  const double beam_z = std::isfinite(track.get_seed_z0()) ? track.get_seed_z0() : track.get_z();
  SpacePoint beam_point;
  beam_point.z = beam_z;
  beam_point.x = TpcBeamXIntercept + TpcBeamXSlope * beam_z;
  beam_point.y = TpcBeamYIntercept + TpcBeamYSlope * beam_z;
  beam_point.r = std::hypot(beam_point.x, beam_point.y);
  beam_point.phi = std::atan2(beam_point.y, beam_point.x);
  if (std::isfinite(beam_point.r) && beam_point.r > 0.0)
  {
    points.push_back(beam_point);
  }

  TrajectoryState fitted = fitTrajectory(points);
  if (!fitted.valid)
  {
    return state;
  }

  fitted.seed_x0 = state.seed_x0;
  fitted.seed_y0 = state.seed_y0;
  fitted.seed_z0 = state.seed_z0;
  fitted.seed_cx = state.seed_cx;
  fitted.seed_cy = state.seed_cy;
  fitted.seed_phi0 = state.seed_phi0;
  fitted.seed_slope = state.seed_slope;
  fitted.seed_q_over_r = fitted.phi_sagitta_ok ? fitted.phi_invR : state.seed_q_over_r;
  fitted.z_slope = state.seed_slope;
  fitted.z_intercept = state.seed_z0 - state.seed_slope * std::hypot(state.seed_x0, state.seed_y0);
  fitted.use_tpc_seed = false;
  fitted.valid = true;
  return fitted;
}

Full_PolyTrackMatcher::TrajectoryState Full_PolyTrackMatcher::makeSiliconSeedTrajectory(
    const TrajectoryState& tpc_state,
    const Tpc_PolyTrack& track,
    const SpacePoint& first_mvtx) const
{
  TrajectoryState state;
  const double q_over_r = std::fabs(tpc_state.seed_q_over_r) > 1.0e-12 ?
                            tpc_state.seed_q_over_r : track.get_seed_q_over_r();
  const double beam_z = first_mvtx.z;
  const double beam_x = SiliconBeamXIntercept + SiliconBeamXSlope * beam_z;
  const double beam_y = SiliconBeamYIntercept + SiliconBeamYSlope * beam_z;
  const double dx = first_mvtx.x - beam_x;
  const double dy = first_mvtx.y - beam_y;
  const double chord = std::hypot(dx, dy);
  if (!std::isfinite(chord) || chord <= 1.0e-9 || !std::isfinite(q_over_r))
  {
    return state;
  }

  state.seed_x0 = beam_x;
  state.seed_y0 = beam_y;
  state.seed_q_over_r = q_over_r;
  state.seed_slope = track.get_seed_slope();
  state.use_tpc_seed = true;
  state.use_silicon_seed = true;

  if (std::fabs(q_over_r) <= 1.0e-12)
  {
    state.seed_phi0 = std::atan2(dy, dx);
    const double arc = chord;
    state.seed_z0 = first_mvtx.z - state.seed_slope * arc;
    state.valid = std::isfinite(state.seed_phi0) && std::isfinite(state.seed_z0) &&
                  std::isfinite(state.seed_slope);
    return state;
  }

  const double radius = 1.0 / std::fabs(q_over_r);
  if (!std::isfinite(radius) || chord > 2.0 * radius)
  {
    return state;
  }

  const double mx = 0.5 * (beam_x + first_mvtx.x);
  const double my = 0.5 * (beam_y + first_mvtx.y);
  const double ux = dx / chord;
  const double uy = dy / chord;
  const double h = std::sqrt(std::max(radius * radius - 0.25 * chord * chord, 0.0));
  const double fit_sign = q_over_r > 0.0 ? -1.0 : 1.0;

  double best_residual = std::numeric_limits<double>::max();
  TrajectoryState best_state;
  for (const double side : {-1.0, 1.0})
  {
    TrajectoryState candidate = state;
    candidate.seed_cx = mx + side * h * (-uy);
    candidate.seed_cy = my + side * h * ux;
    const double start_angle = std::atan2(beam_y - candidate.seed_cy, beam_x - candidate.seed_cx);
    const double hit_angle = std::atan2(first_mvtx.y - candidate.seed_cy, first_mvtx.x - candidate.seed_cx);
    double dangle = unwrapPhiNear(hit_angle, start_angle) - start_angle;
    if (fit_sign * dangle < 0.0)
    {
      dangle += fit_sign * 2.0 * M_PI;
    }
    const double arc = std::fabs(radius * dangle);
    candidate.seed_phi0 = start_angle + fit_sign * 0.5 * M_PI;
    candidate.seed_z0 = first_mvtx.z - candidate.seed_slope * arc;
    candidate.valid = std::isfinite(candidate.seed_cx) && std::isfinite(candidate.seed_cy) &&
                      std::isfinite(candidate.seed_phi0) && std::isfinite(candidate.seed_z0);

    double pred_phi = 0.0;
    double pred_z = 0.0;
    double pred_x = 0.0;
    double pred_y = 0.0;
    if (candidate.valid && predictAtRadius(candidate, first_mvtx.r, pred_phi, pred_z, pred_x, pred_y))
    {
      const double residual = square(first_mvtx.r * wrapPhi(first_mvtx.phi - pred_phi)) + square(first_mvtx.z - pred_z);
      if (residual < best_residual)
      {
        best_residual = residual;
        best_state = candidate;
      }
    }
  }

  return best_state;
}

Full_PolyTrackMatcher::TrajectoryState Full_PolyTrackMatcher::correctSiliconSeedWithTwoHits(
    const TrajectoryState& seed_state,
    const ChainHit& outer_hit,
    const ChainHit& inner_hit) const
{
  TrajectoryState state = seed_state;
  if (!state.valid || std::fabs(state.seed_q_over_r) <= 1.0e-12)
  {
    return state;  // straight-line case: this correction only applies to curved seeds
  }

  const double radius = 1.0 / std::fabs(state.seed_q_over_r);   // curvature preserved, never refit
  const double x1 = outer_hit.point.x, y1 = outer_hit.point.y, z1 = outer_hit.point.z;
  const double x2 = inner_hit.point.x, y2 = inner_hit.point.y, z2 = inner_hit.point.z;
  const double dx = x2 - x1, dy = y2 - y1;
  const double chord = std::hypot(dx, dy);
  if (!std::isfinite(chord) || chord <= 1.0e-9 || chord > 2.0 * radius)
  {
    return state;  // hits inconsistent with this radius -- leave seed uncorrected
  }

  const double mx = 0.5 * (x1 + x2), my = 0.5 * (y1 + y2);
  const double ux = dx / chord, uy = dy / chord;
  const double h = std::sqrt(std::max(radius * radius - 0.25 * chord * chord, 0.0));

  double best_d2 = std::numeric_limits<double>::max();
  double cx = 0.0, cy = 0.0;
  for (const double side : {-1.0, 1.0})
  {
    const double ccx = mx + side * h * (-uy);
    const double ccy = my + side * h * ux;
    const double d2 = square(ccx - state.seed_cx) + square(ccy - state.seed_cy);
    if (d2 < best_d2) { best_d2 = d2; cx = ccx; cy = ccy; }
  }

  const double ox = state.seed_x0, oy = state.seed_y0;
  const double ocx = cx - ox, ocy = cy - oy;
  const double L = std::hypot(ocx, ocy);
  if (!std::isfinite(L) || L <= 1.0e-9)
  {
    return state;
  }
  const double dca = L - radius;
  const double fit_sign = state.seed_q_over_r > 0.0 ? -1.0 : 1.0;

  state.seed_x0 = ox + (ocx / L) * dca;
  state.seed_y0 = oy + (ocy / L) * dca;
  state.seed_cx = cx;
  state.seed_cy = cy;
  state.seed_phi0 = std::atan2(state.seed_y0 - cy, state.seed_x0 - cx) + fit_sign * 0.5 * M_PI;
  // Independent Si-side z-slope (theta) and z-anchor, derived from the two real hits
  // themselves -- NOT copied from track.get_seed_slope() (the TPC track's own slope), so
  // that delta_theta0/delta_z0 against the TPC trajectory are a real, non-circular check.
  const double start_angle = std::atan2(state.seed_y0 - cy, state.seed_x0 - cx);
  auto arcTo = [&](double px, double py) -> double
  {
    double dangle = unwrapPhiNear(std::atan2(py - cy, px - cx), start_angle) - start_angle;
    if (fit_sign * dangle < 0.0)
    {
      dangle += fit_sign * 2.0 * M_PI;
    }
    return std::fabs(radius * dangle);  // same arc convention as predictAtRadius()
  };
  const double s_outer = arcTo(x1, y1);
  const double s_inner = arcTo(x2, y2);
  if (std::isfinite(s_outer) && std::isfinite(s_inner) && std::fabs(s_outer - s_inner) > 1.0e-6)
  {
    state.seed_slope = (z1 - z2) / (s_outer - s_inner);
    state.seed_z0 = z2 - state.seed_slope * s_inner;
  }
  state.valid = std::isfinite(state.seed_x0) && std::isfinite(state.seed_y0) &&
                std::isfinite(state.seed_phi0) && std::isfinite(state.seed_slope) &&
                std::isfinite(state.seed_z0);
  return state;
}

Full_PolyTrackMatcher::TrajectoryState Full_PolyTrackMatcher::makeTpcSeedTrajectory(const Tpc_PolyTrack& track) const
{
  TrajectoryState state;
  if (track.get_fit_status() == 0)
  {
    return state;
  }

  state.seed_x0 = track.get_seed_x0();
  state.seed_y0 = track.get_seed_y0();
  state.seed_z0 = track.get_seed_z0();
  state.seed_cx = track.get_helix_x0();
  state.seed_cy = track.get_helix_y0();
  state.seed_phi0 = track.get_seed_phi();
  state.seed_slope = track.get_seed_slope();
  state.seed_q_over_r = track.get_seed_q_over_r();
  state.use_tpc_seed = true;
  state.valid = std::isfinite(state.seed_x0) && std::isfinite(state.seed_y0) &&
                std::isfinite(state.seed_z0) && std::isfinite(state.seed_phi0) &&
                std::isfinite(state.seed_slope) && std::isfinite(state.seed_q_over_r);
  if (std::fabs(state.seed_q_over_r) > 1.0e-12)
  {
    state.valid = state.valid && std::isfinite(state.seed_cx) && std::isfinite(state.seed_cy);
  }
  return state;
}



void Full_PolyTrackMatcher::updateSiliconTrajectoryHalfResidual(TrajectoryState& state, const ChainHit& hit) const
{
  if (!state.valid || !state.use_silicon_seed)
  {
    return;
  }

  const double shift = 0.5 * hit.dphi;
  const double bx = state.seed_x0;
  const double by = state.seed_y0;
  const double dcx = state.seed_cx - bx;
  const double dcy = state.seed_cy - by;
  const double c = std::cos(shift);
  const double s = std::sin(shift);
  if (std::isfinite(shift) && std::isfinite(dcx) && std::isfinite(dcy))
  {
    state.seed_cx = bx + c * dcx - s * dcy;
    state.seed_cy = by + s * dcx + c * dcy;
    state.seed_phi0 = wrapPhi(state.seed_phi0 + shift);
  }
  if (std::isfinite(hit.dz))
  {
    state.seed_z0 += 0.5 * hit.dz;
  }
}

void Full_PolyTrackMatcher::refitSiliconTrajectoryFromMvtx(Chain& chain) const
{
  if (!chain.state.valid || !chain.state.use_silicon_seed || chain.hits.empty())
  {
    return;
  }

  double sum_dphi = 0.0;
  double sum_dz = 0.0;
  unsigned int count = 0;
  for (const ChainHit& hit : chain.hits)
  {
    if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(hit.point.key)) != TrkrDefs::mvtxId)
    {
      continue;
    }
    double pred_phi = 0.0;
    double pred_z = 0.0;
    double pred_x = 0.0;
    double pred_y = 0.0;
    if (predictAtRadius(chain.state, hit.point.r, pred_phi, pred_z, pred_x, pred_y))
    {
      sum_dphi += wrapPhi(hit.point.phi - pred_phi);
      sum_dz += hit.point.z - pred_z;
      ++count;
    }
  }
  if (count == 0U)
  {
    return;
  }

  ChainHit correction;
  correction.dphi = sum_dphi / static_cast<double>(count);
  correction.dz = sum_dz / static_cast<double>(count);
  updateSiliconTrajectoryHalfResidual(chain.state, correction);
}

void Full_PolyTrackMatcher::refitSiliconTrajectoryFromHits(Chain& chain) const
{
  if (chain.hits.size() < 2U)
  {
    return;
  }

  std::vector<SpacePoint> fit_points;
  fit_points.reserve(chain.hits.size());
  for (const ChainHit& hit : chain.hits)
  {
    fit_points.push_back(hit.point);
  }

  TrajectoryState refit = fitTrajectory(fit_points);
  if (!refit.valid)
  {
    return;
  }

  refit.use_tpc_seed = false;
  refit.use_silicon_seed = true;
  refit.seed_z0 = chain.si_reference.seed_z0;
  chain.state = refit;
  chain.fit_points = fit_points;
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
  phi_points.reserve(points.size());
  z_points.reserve(points.size());

  const double ref_phi = points.front().phi;
  for (const SpacePoint& point : points)
  {
    phi_points.emplace_back(point.r, unwrapPhiNear(point.phi, ref_phi));
    z_points.emplace_back(point.r, point.z);
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

  state.z_slope = z_line.slope;
  state.z_intercept = z_line.intercept;
  state.valid = std::isfinite(state.phi_slope) && std::isfinite(state.phi_intercept) &&
                std::isfinite(state.z_slope) && std::isfinite(state.z_intercept);
  return state;
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

double Full_PolyTrackMatcher::pointTheta(const SpacePoint& point) const
{
  return std::atan2(point.r, point.z);
}

double Full_PolyTrackMatcher::trajectoryPhi0NearBeam(const TrajectoryState& state) const
{
  if (!state.valid)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (state.use_silicon_seed && state.use_tpc_seed)
  {
    return state.seed_phi0;
  }

  const double beam_z = std::isfinite(state.seed_z0) ? state.seed_z0 : 0.0;
  const double beam_x = TpcBeamXIntercept + TpcBeamXSlope * beam_z;
  const double beam_y = TpcBeamYIntercept + TpcBeamYSlope * beam_z;
  const double r_beam = std::hypot(beam_x, beam_y);
  const double r1 = std::max(r_beam + 0.01, 0.01);
  const double r2 = r1 + 0.01;

  double phi1 = 0.0;
  double z1 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  double phi2 = 0.0;
  double z2 = 0.0;
  double x2 = 0.0;
  double y2 = 0.0;
  if (predictAtRadius(state, r1, phi1, z1, x1, y1) &&
      predictAtRadius(state, r2, phi2, z2, x2, y2))
  {
    const double phi0 = std::atan2(y2 - y1, x2 - x1);
    if (std::isfinite(phi0))
    {
      return phi0;
    }
  }

  return state.use_tpc_seed ? state.seed_phi0 : state.phi_bline;
}

double Full_PolyTrackMatcher::trajectoryZ0NearBeam(const TrajectoryState& state) const
{
  if (!state.valid)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (state.use_tpc_seed)
  {
    return state.seed_z0;
  }
  return state.z_intercept;
}

double Full_PolyTrackMatcher::trajectoryTheta0NearBeam(const TrajectoryState& state) const
{
  if (!state.valid)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (state.use_tpc_seed)
  {
    return std::atan2(1.0, state.seed_slope);
  }
  return std::atan2(1.0, state.z_slope);
}

bool Full_PolyTrackMatcher::computeChainDcaMetrics(Chain& chain, const TrajectoryState& tpc_state) const
{
  const double tpc_phi0 = trajectoryPhi0NearBeam(tpc_state);
  const double silicon_phi0 = trajectoryPhi0NearBeam(chain.state);
  const double dphi0 = std::fabs(wrapPhi(silicon_phi0 - tpc_phi0));
  const double dphi0_flip = std::fabs(wrapPhi(silicon_phi0 + M_PI - tpc_phi0));
  chain.delta_phi0 = std::min(dphi0, dphi0_flip);
  chain.signed_delta_phi0 = wrapPhi(silicon_phi0 - tpc_phi0);

  const double tpc_z0 = trajectoryZ0NearBeam(tpc_state);
  const double silicon_z0 = trajectoryZ0NearBeam(chain.state);
  chain.delta_z0 = std::fabs(silicon_z0 - tpc_z0);
  chain.signed_delta_z0 = silicon_z0 - tpc_z0;

  const double tpc_theta0 = trajectoryTheta0NearBeam(tpc_state);
  const double silicon_theta0 = trajectoryTheta0NearBeam(chain.state);
  chain.delta_theta0 = std::fabs(silicon_theta0 - tpc_theta0);
  chain.signed_delta_theta0 = silicon_theta0 - tpc_theta0;

  const double tpc_eta0 = thetaToEta(tpc_theta0);
  const double silicon_eta0 = thetaToEta(silicon_theta0);
  chain.delta_eta0 = std::fabs(silicon_eta0 - tpc_eta0);

  const double phi_scale = SeedPhiWindow;
  const double theta_scale = seedThetaWindow(chain.pt, tpc_theta0);
  const double z_scale = dynamicDzWindow(chain.pt);
  const bool valid = std::isfinite(chain.delta_phi0) &&
                     std::isfinite(chain.delta_z0) &&
                     std::isfinite(chain.delta_theta0) &&
                     std::isfinite(chain.delta_eta0) &&
                     std::isfinite(phi_scale) && phi_scale > 0.0 &&
                     std::isfinite(theta_scale) && theta_scale > 0.0 &&
                     std::isfinite(z_scale) && z_scale > 0.0;
  if (!valid)
  {
    return false;
  }

  chain.dca_score = square(chain.delta_phi0 / phi_scale) +
                    square(chain.delta_theta0 / theta_scale) +
                    square(chain.delta_z0 / z_scale);
  return std::isfinite(chain.dca_score);
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
  (void) pt;
  return std::max(m_sigmaPhi, 1.0e-9);
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
  const double theta_window = seedThetaWindow(pt, pred_theta);
  return std::max(theta_window / std::max(m_thetaWindowSigma, 1.0e-9), 1.0e-9);
}

double Full_PolyTrackMatcher::dynamicDzWindow(double pt) const
{
  const double safe_pt = std::max(pt, 0.25);
  return 1.138 + 0.3919 * std::exp(0.84 / safe_pt);
}

bool Full_PolyTrackMatcher::predictAtRadius(const TrajectoryState& state, double r,
                                            double& pred_phi, double& pred_z,
                                            double& pred_x, double& pred_y) const
{
  if (!state.valid || !std::isfinite(r) || r <= 0.0)
  {
    return false;
  }

  if (state.use_tpc_seed)
  {
    double arc = std::numeric_limits<double>::quiet_NaN();
    if (std::fabs(state.seed_q_over_r) <= 1.0e-12)
    {
      const double tx = std::cos(state.seed_phi0);
      const double ty = std::sin(state.seed_phi0);
      const double b = state.seed_x0 * tx + state.seed_y0 * ty;
      const double c = state.seed_x0 * state.seed_x0 + state.seed_y0 * state.seed_y0 - r * r;
      const double discriminant = b * b - c;
      if (discriminant < 0.0)
      {
        return false;
      }
      const double root = std::sqrt(discriminant);
      const double s1 = -b + root;
      const double s2 = -b - root;
      if (s1 >= 0.0 && s2 >= 0.0)
      {
        arc = std::min(s1, s2);
      }
      else if (s1 >= 0.0)
      {
        arc = s1;
      }
      else if (s2 >= 0.0)
      {
        arc = s2;
      }
      else
      {
        arc = std::fabs(s1) < std::fabs(s2) ? s1 : s2;
      }
      pred_x = state.seed_x0 + arc * tx;
      pred_y = state.seed_y0 + arc * ty;
    }
    else
    {
      const double radius = 1.0 / std::fabs(state.seed_q_over_r);
      const double dc = std::hypot(state.seed_cx, state.seed_cy);
      if (!std::isfinite(radius) || radius <= 0.0 || dc <= 1.0e-12)
      {
        return false;
      }

      const double a = (r * r - radius * radius + dc * dc) / (2.0 * dc);
      const double h2 = r * r - a * a;
      if (h2 < -1.0e-9)
      {
        return false;
      }
      const double h = std::sqrt(std::max(h2, 0.0));
      const double ux = state.seed_cx / dc;
      const double uy = state.seed_cy / dc;
      const double base_x = a * ux;
      const double base_y = a * uy;
      const double cand_x[2] = {base_x - h * uy, base_x + h * uy};
      const double cand_y[2] = {base_y + h * ux, base_y - h * ux};
      const double start_angle = std::atan2(state.seed_y0 - state.seed_cy, state.seed_x0 - state.seed_cx);
      const double fit_sign = state.seed_q_over_r > 0.0 ? -1.0 : 1.0;

      double best_arc = std::numeric_limits<double>::max();
      unsigned int best = 0;
      for (unsigned int i = 0; i < 2; ++i)
      {
        double dangle = unwrapPhiNear(std::atan2(cand_y[i] - state.seed_cy, cand_x[i] - state.seed_cx), start_angle) - start_angle;
        if (fit_sign * dangle < 0.0)
        {
          dangle += fit_sign * 2.0 * M_PI;
        }
        const double candidate_arc = std::fabs(radius * dangle);
        if (candidate_arc < best_arc)
        {
          best_arc = candidate_arc;
          best = i;
        }
      }

      arc = best_arc;
      pred_x = cand_x[best];
      pred_y = cand_y[best];
    }

    pred_phi = std::atan2(pred_y, pred_x);
    pred_z = state.seed_z0 + state.seed_slope * arc;
    return std::isfinite(pred_phi) && std::isfinite(pred_z) && std::isfinite(pred_x) && std::isfinite(pred_y);
  }

  pred_phi = wrapPhi(state.phi_sagitta_ok ? predictSagittaPhi(r, state) : state.phi_intercept + state.phi_slope * r);
  pred_z = state.z_intercept + state.z_slope * r;
  pred_x = r * std::cos(pred_phi);
  pred_y = r * std::sin(pred_phi);
  return std::isfinite(pred_phi) && std::isfinite(pred_z) && std::isfinite(pred_x) && std::isfinite(pred_y);
}

Full_PolyTrackMatcher::Chain Full_PolyTrackMatcher::extendWithHit(const Chain& chain, const ChainHit& hit) const
{
  Chain out = chain;
  out.hits.push_back(hit);
  out.chi2 += hit.chi2;
  out.ndf += 2.0;
  out.previous_dphi = hit.dphi;
  out.previous_dtheta = hit.dtheta;
  out.has_previous_residual = true;
  if (!out.has_z_offset)
  {
    out.z_offset = hit.dz;
    out.has_z_offset = true;
  }
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

void Full_PolyTrackMatcher::fillQaResiduals(const Tpc_PolyTrack& track, const Chain& chain, const ChainHit& hit)
{
  if (!m_writeQA || !m_qaFile || chain.hits.empty())
  {
    return;
  }

  const ChainHit& previous_hit = chain.hits.back();
  const double pt = std::hypot(track.get_px(), track.get_py());
  if (!std::isfinite(pt) || !std::isfinite(hit.dphi) || !std::isfinite(hit.dtheta) ||
      !std::isfinite(previous_hit.dphi) || !std::isfinite(previous_hit.dtheta))
  {
    return;
  }

  const unsigned int side = hit.point.z >= 0.0 ? 1U : 0U;
  const int charge_bin = track.get_charge() < 0.0 ? -1 : 1;
  const auto angular_bin = [](double value, unsigned int bins) {
    if (bins <= 1U)
    {
      return 0U;
    }
    double wrapped = value;
    while (wrapped > M_PI) { wrapped -= 2.0 * M_PI; }
    while (wrapped <= -M_PI) { wrapped += 2.0 * M_PI; }
    unsigned int bin = static_cast<unsigned int>(bins * (wrapped + M_PI) / (2.0 * M_PI));
    if (bin >= bins)
    {
      bin = bins - 1U;
    }
    return bin;
  };
  const auto angular_bin_min = [](unsigned int bin, unsigned int bins) {
    return -M_PI + 2.0 * M_PI * static_cast<double>(bin) / static_cast<double>(bins);
  };
  const auto angular_bin_max = [](unsigned int bin, unsigned int bins) {
    return -M_PI + 2.0 * M_PI * static_cast<double>(bin + 1U) / static_cast<double>(bins);
  };

  const std::string charge_label = charge_bin < 0 ? "neg" : "pos";
  const double hit_phi = wrapPhi(hit.point.phi);
  const double hit_theta = pointTheta(hit.point);
  const double tpc_track_phi = wrapPhi(std::atan2(track.get_py(), track.get_px()));
  constexpr double tpc_residual_min = -0.2;
  constexpr double tpc_residual_max = 0.2;
  constexpr double si_residual_min = -0.2;
  constexpr double si_residual_max = 0.2;
  constexpr unsigned int current_previous_angular_bins = 1U;
  constexpr unsigned int standard_residual_angular_bins = 4U;
  const unsigned int current_previous_phi_bin = angular_bin(hit_phi, current_previous_angular_bins);
  const unsigned int current_previous_theta_bin = angular_bin(hit_theta, current_previous_angular_bins);
  const unsigned int standard_residual_phi_bin = angular_bin(hit_phi, standard_residual_angular_bins);
  const unsigned int standard_residual_theta_bin = angular_bin(hit_theta, standard_residual_angular_bins);

  auto fill_3d = [&](const std::string& prefix, const std::string& title_prefix,
                     double current, double previous, double min, double max) {
    std::ostringstream key;
    key << prefix << "_l" << hit.point.layer
        << "_side" << side
        << "_q" << charge_label
        << "_phibin" << current_previous_phi_bin
        << "_thetabin" << current_previous_theta_bin;

    TH3F* hist = nullptr;
    auto found = m_hResidualCurrentVsPreviousVsPt.find(key.str());
    if (found != m_hResidualCurrentVsPreviousVsPt.end())
    {
      hist = found->second;
    }
    else
    {
      std::ostringstream title;
      title << title_prefix << " layer " << hit.point.layer
            << " side " << side
            << " q " << (charge_bin < 0 ? "-" : "+")
            << " phi bin " << current_previous_phi_bin << " ["
            << angular_bin_min(current_previous_phi_bin, current_previous_angular_bins) << ", "
            << angular_bin_max(current_previous_phi_bin, current_previous_angular_bins) << ") rad"
            << " theta bin " << current_previous_theta_bin << " ["
            << angular_bin_min(current_previous_theta_bin, current_previous_angular_bins) << ", "
            << angular_bin_max(current_previous_theta_bin, current_previous_angular_bins) << ") rad"
            << ";current residual [rad];previous accepted layer residual [rad];p_{T} [GeV]";
      hist = new TH3F(key.str().c_str(), title.str().c_str(),
                      static_cast<int>(m_qaDphiBins), min, max,
                      static_cast<int>(m_qaDphiBins), min, max,
                      static_cast<int>(m_qaPtBins), m_qaPtMin, m_qaPtMax);
      hist->SetDirectory(nullptr);
      m_hResidualCurrentVsPreviousVsPt[key.str()] = hist;
    }
    hist->Fill(current, previous, pt);
  };

  fill_3d("h_si_dphi_curr_vs_prev_vs_pt", "fixed Si #Delta#phi current vs previous", hit.dphi, previous_hit.dphi, si_residual_min, si_residual_max);
  fill_3d("h_si_dtheta_curr_vs_prev_vs_pt", "fixed Si #Delta#theta current vs previous", hit.dtheta, previous_hit.dtheta, si_residual_min, si_residual_max);
  fill_3d("h_tpc_dphi_curr_vs_prev_vs_pt", "TPC #Delta#phi current vs previous", hit.tpc_dphi, previous_hit.tpc_dphi, tpc_residual_min, tpc_residual_max);
  fill_3d("h_tpc_dtheta_curr_vs_prev_vs_pt", "TPC #Delta#theta current vs previous", hit.tpc_dtheta, previous_hit.tpc_dtheta, tpc_residual_min, tpc_residual_max);

  auto fill_si_vs_tpc = [&](const std::string& key_base, const std::string& title, double tpc_residual, double si_residual,
                          double tpc_min, double tpc_max, double si_min, double si_max) {
    std::ostringstream key;
    key << key_base << "_l" << hit.point.layer;
    TH3F* hist = nullptr;
    auto found = m_hResidualSiVsTpc.find(key.str());
    if (found != m_hResidualSiVsTpc.end())
    {
      hist = found->second;
    }
    else
    {
      std::ostringstream hist_title;
      hist_title << title << " layer " << hit.point.layer
                 << ";TPC residual [rad];fixed Si residual [rad];TPC track #phi [rad]";
      hist = new TH3F(key.str().c_str(), hist_title.str().c_str(),
                      static_cast<int>(m_qaDphiBins), tpc_min, tpc_max,
                      static_cast<int>(m_qaDphiBins), si_min, si_max,
                      static_cast<int>(m_qaDphiBins), -M_PI, M_PI);
      hist->SetDirectory(nullptr);
      m_hResidualSiVsTpc[key.str()] = hist;
    }
    hist->Fill(tpc_residual, si_residual, tpc_track_phi);
  };

  fill_si_vs_tpc("h_dphi_si_vs_tpc", "fixed Si vs TPC #Delta#phi", hit.tpc_dphi, hit.dphi, tpc_residual_min, tpc_residual_max, si_residual_min, si_residual_max);
  fill_si_vs_tpc("h_dtheta_si_vs_tpc", "fixed Si vs TPC #Delta#theta", hit.tpc_dtheta, hit.dtheta, tpc_residual_min, tpc_residual_max, si_residual_min, si_residual_max);

  auto fill_1d = [&](const std::string& prefix, const std::string& title, double value) {
    std::ostringstream key;
    key << prefix << "_l" << hit.point.layer
        << "_side" << side
        << "_q" << charge_label
        << "_phibin" << standard_residual_phi_bin
        << "_thetabin" << standard_residual_theta_bin;
    TH1F* hist = nullptr;
    auto found = m_hStandardResidual.find(key.str());
    if (found != m_hStandardResidual.end())
    {
      hist = found->second;
    }
    else
    {
      std::ostringstream hist_title;
      hist_title << title << " layer " << hit.point.layer
                 << " side " << side << " q " << (charge_bin < 0 ? "-" : "+")
                 << " phi bin " << standard_residual_phi_bin << " ["
                 << angular_bin_min(standard_residual_phi_bin, standard_residual_angular_bins) << ", "
                 << angular_bin_max(standard_residual_phi_bin, standard_residual_angular_bins) << ") rad"
                 << " theta bin " << standard_residual_theta_bin << " ["
                 << angular_bin_min(standard_residual_theta_bin, standard_residual_angular_bins) << ", "
                 << angular_bin_max(standard_residual_theta_bin, standard_residual_angular_bins) << ") rad"
                 << ";standard residual;entries";
      hist = new TH1F(key.str().c_str(), hist_title.str().c_str(), 120, -12.0, 12.0);
      hist->SetDirectory(nullptr);
      m_hStandardResidual[key.str()] = hist;
    }
    hist->Fill(value);
  };

  fill_1d("h_si_sdphi", "fixed Si s_{#phi}", hit.sdphi);
  fill_1d("h_si_sdtheta", "fixed Si s_{#theta}", hit.sdtheta);
}

bool Full_PolyTrackMatcher::findInttCrossing(const Chain& chain, short& crossing) const
{
  if (!m_clusterCrossingAssoc)
  {
    return false;
  }

  std::map<short, unsigned int> counts;
  for (const ChainHit& hit : chain.hits)
  {
    if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(hit.point.key)) != TrkrDefs::inttId)
    {
      continue;
    }
    const auto range = m_clusterCrossingAssoc->getCrossings(hit.point.key);
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      ++counts[iter->second];
    }
  }

  if (counts.empty())
  {
    return false;
  }

  crossing = counts.begin()->first;
  unsigned int best_count = counts.begin()->second;
  for (const auto& item : counts)
  {
    if (item.second > best_count || (item.second == best_count && std::abs(item.first) < std::abs(crossing)))
    {
      crossing = item.first;
      best_count = item.second;
    }
  }
  return true;
}

void Full_PolyTrackMatcher::fillQaCrossing(const Tpc_PolyTrack& track, const Chain& chain)
{
  if (!m_writeQA || !m_qaFile)
  {
    return;
  }

  if (!m_hInttCrossingVsTpcCrossing)
  {
    m_hInttCrossingVsTpcCrossing = new TH2F("h_intt_crossing_vs_tpc_crossing",
                                            "INTT crossing vs TPC crossing;TPC crossing;INTT crossing",
                                            201, -100.5, 100.5, 201, -100.5, 100.5);
    m_hInttCrossingVsTpcCrossing->SetDirectory(nullptr);
    m_hDeltaCrossing = new TH1F("h_delta_crossing", "INTT - TPC crossing;#Delta crossing;tracks", 201, -100.5, 100.5);
    m_hDeltaCrossing->SetDirectory(nullptr);
    m_hCrossingStatus = new TH1F("h_crossing_status", "crossing availability/status;status;tracks", 5, 0.5, 5.5);
    m_hCrossingStatus->GetXaxis()->SetBinLabel(1, "both valid");
    m_hCrossingStatus->GetXaxis()->SetBinLabel(2, "TPC only");
    m_hCrossingStatus->GetXaxis()->SetBinLabel(3, "INTT only");
    m_hCrossingStatus->GetXaxis()->SetBinLabel(4, "neither valid");
    m_hCrossingStatus->GetXaxis()->SetBinLabel(5, "disagreement");
    m_hCrossingStatus->SetDirectory(nullptr);
  }

  const TpcCrossingDecision* decision = findCrossingDecision(track.get_source_assembled_track_id());
  const bool has_tpc = decision && decision->get_status() != static_cast<unsigned char>(TpcCrossingStatus::NoValidCrossing) &&
                       decision->get_status() != static_cast<unsigned char>(TpcCrossingStatus::Unknown);
  const short tpc_crossing = has_tpc ? decision->get_selected_crossing() : 0;

  short intt_crossing = 0;
  const bool has_intt = findInttCrossing(chain, intt_crossing);

  if (has_tpc && has_intt)
  {
    m_hCrossingStatus->Fill(1.0);
    m_hInttCrossingVsTpcCrossing->Fill(tpc_crossing, intt_crossing);
    m_hDeltaCrossing->Fill(static_cast<double>(intt_crossing) - static_cast<double>(tpc_crossing));
    if (intt_crossing != tpc_crossing)
    {
      m_hCrossingStatus->Fill(5.0);
    }
  }
  else if (has_tpc)
  {
    m_hCrossingStatus->Fill(2.0);
  }
  else if (has_intt)
  {
    m_hCrossingStatus->Fill(3.0);
  }
  else
  {
    m_hCrossingStatus->Fill(4.0);
  }
}

const Full_PolyTrackMatcher::SpacePoint* Full_PolyTrackMatcher::findBestMvtxCandidate(
    const Chain& chain,
    const std::vector<SpacePoint>& silicon_points,
    const std::set<TrkrDefs::cluskey>& used_keys,
    unsigned int layer,
    ChainHit& best_hit,
    double& best_chi2) const
{
  const SpacePoint* best_point = nullptr;
  best_chi2 = std::numeric_limits<double>::max();

  for (const SpacePoint& point : silicon_points)
  {
    if (point.layer != layer ||
        static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(point.key)) != TrkrDefs::mvtxId ||
        used_keys.find(point.key) != used_keys.end())
    {
      continue;
    }

    double pred_phi = 0.0;
    double pred_z = 0.0;
    double pred_x = 0.0;
    double pred_y = 0.0;
    if (!predictAtRadius(chain.si_reference, point.r, pred_phi, pred_z, pred_x, pred_y))
    {
      continue;
    }

    const double dphi = wrapPhi(point.phi - pred_phi);
    const double rdphi = point.r * dphi;
    const double dz = point.z - pred_z;
    if (std::fabs(rdphi) > m_looseRdphiWindow || std::fabs(dz) > m_looseDzWindow)
    {
      continue;
    }

    double tpc_pred_phi = 0.0;
    double tpc_pred_z = 0.0;
    double tpc_pred_x = 0.0;
    double tpc_pred_y = 0.0;
    if (!predictAtRadius(chain.tpc_reference, point.r, tpc_pred_phi, tpc_pred_z, tpc_pred_x, tpc_pred_y))
    {
      continue;
    }

    double si_ref_pred_phi = 0.0;
    double si_ref_pred_z = 0.0;
    double si_ref_pred_x = 0.0;
    double si_ref_pred_y = 0.0;
    double si_ref_dphi = 0.0;
    double si_ref_dtheta = 0.0;
    if (predictAtRadius(chain.si_reference, point.r, si_ref_pred_phi, si_ref_pred_z, si_ref_pred_x, si_ref_pred_y))
    {
      si_ref_dphi = wrapPhi(point.phi - si_ref_pred_phi);
      si_ref_dtheta = pointTheta(point) - std::atan2(point.r, si_ref_pred_z);
    }

    const double pred_theta = std::atan2(point.r, pred_z);
    const double dtheta = pointTheta(point) - pred_theta;
    const double tpc_pred_theta = std::atan2(point.r, tpc_pred_z);
    const double tpc_dphi = wrapPhi(point.phi - tpc_pred_phi);
    const double tpc_dtheta = pointTheta(point) - tpc_pred_theta;
    const double mean_phi = dynamicMeanPhi(layer, chain.previous_dphi, chain.has_previous_residual);
    const double mean_theta = dynamicMeanTheta(layer, chain.previous_dtheta, chain.has_previous_residual);
    const double sigma_phi = dynamicSigmaPhi(chain.pt);
    const double sigma_theta = dynamicSigmaTheta(chain.pt, pred_theta);
    const double sdphi = (dphi - mean_phi) / sigma_phi;
    const double sdtheta = (dtheta - mean_theta) / sigma_theta;

    if (!m_associationCalibrationMode &&
        (std::fabs(sdphi) > m_phiWindowSigma || std::fabs(sdtheta) > m_thetaWindowSigma))
    {
      continue;
    }

    const double chi2 = square(sdphi) + square(sdtheta);
    if (chi2 < best_chi2)
    {
      best_chi2 = chi2;
      best_point = &point;
      best_hit.point = point;
      best_hit.pred_phi = pred_phi;
      best_hit.pred_z = pred_z;
      best_hit.pred_theta = pred_theta;
      best_hit.pred_x = pred_x;
      best_hit.pred_y = pred_y;
      best_hit.dphi = dphi;
      best_hit.dtheta = dtheta;
      best_hit.rdphi = rdphi;
      best_hit.dz = dz;
      best_hit.tpc_pred_phi = tpc_pred_phi;
      best_hit.tpc_pred_z = tpc_pred_z;
      best_hit.tpc_pred_theta = tpc_pred_theta;
      best_hit.tpc_dphi = tpc_dphi;
      best_hit.tpc_dtheta = tpc_dtheta;
      best_hit.si_ref_pred_phi = si_ref_pred_phi;
      best_hit.si_ref_pred_z = si_ref_pred_z;
      best_hit.si_ref_dphi = si_ref_dphi;
      best_hit.si_ref_dtheta = si_ref_dtheta;
      best_hit.dynamic_mean_phi = mean_phi;
      best_hit.dynamic_mean_theta = mean_theta;
      best_hit.sigma_phi = sigma_phi;
      best_hit.sigma_theta = sigma_theta;
      best_hit.sdphi = sdphi;
      best_hit.sdtheta = sdtheta;
      best_hit.chi2 = chi2;
    }
  }

  return best_point;
}

std::vector<Full_PolyTrackMatcher::Chain> Full_PolyTrackMatcher::buildChains(
    const Tpc_PolyTrack& track,
    const std::vector<SpacePoint>& silicon_points)
{
  Chain tpc_seed;
  tpc_seed.pt = std::hypot(track.get_px(), track.get_py());
  tpc_seed.charge = track.get_charge();
  tpc_seed.state = makeTpcReferenceTrajectory(track);
  tpc_seed.tpc_reference = tpc_seed.state;
  if (!tpc_seed.state.valid)
  {
    return {};
  }

  std::vector<Chain> chains;
  for (const unsigned int seed_layer : m_matchLayers)
  {
    for (const SpacePoint& seed_point : silicon_points)
    {
      if (seed_point.layer != seed_layer ||
          static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(seed_point.key)) != TrkrDefs::mvtxId)
      {
        continue;
      }

      double tpc_pred_phi = 0.0;
      double tpc_pred_z = 0.0;
      double tpc_pred_x = 0.0;
      double tpc_pred_y = 0.0;
      if (!predictAtRadius(tpc_seed.state, seed_point.r, tpc_pred_phi, tpc_pred_z, tpc_pred_x, tpc_pred_y))
      {
        continue;
      }

      const double initial_dphi = wrapPhi(seed_point.phi - tpc_pred_phi);
      const double seed_theta = pointTheta(seed_point);
      const double tpc_seed_theta = std::atan2(seed_point.r, tpc_pred_z);
      const double phi_window = SeedPhiWindow+0.15;
      const double theta_window = SeedThetaPreWindow;
      //const double seed_eta = thetaToEta(seed_theta);
      if (std::fabs(initial_dphi) > phi_window  ||
          std::fabs(seed_theta - tpc_seed_theta) > theta_window)
      {
        continue;
      }
      /*const double tpc_seed_eta = thetaToEta(tpc_seed_theta);
      const double eta_window = seedEtaWindow(tpc_seed.pt);

      if (std::fabs(initial_dphi) > phi_window ||
          !std::isfinite(seed_eta) || !std::isfinite(tpc_seed_eta) ||
          std::fabs(seed_eta - tpc_seed_eta) > eta_window)
      {
        continue;
      }*/
     
      Chain chain;
      chain.pt = tpc_seed.pt;
      chain.charge = tpc_seed.charge;
      chain.tpc_reference = tpc_seed.state;
      chain.state = makeSiliconSeedTrajectory(tpc_seed.state, track, seed_point);
      if (!chain.state.valid)
      {
        continue;
      }
      chain.si_reference = chain.state;

      ChainHit seed_hit;
      double seed_chi2 = 0.0;
      std::set<TrkrDefs::cluskey> chain_used;
      const SpacePoint* accepted_seed = findBestMvtxCandidate(chain, silicon_points, chain_used,
                                                              seed_layer, seed_hit, seed_chi2);
      if (!accepted_seed || accepted_seed->key != seed_point.key)
      {
        double pred_phi = 0.0;
        double pred_z = 0.0;
        double pred_x = 0.0;
        double pred_y = 0.0;
        if (!predictAtRadius(chain.si_reference, seed_point.r, pred_phi, pred_z, pred_x, pred_y))
        {
          continue;
        }
        seed_hit.point = seed_point;
        seed_hit.pred_phi = pred_phi;
        seed_hit.pred_z = pred_z;
        seed_hit.pred_theta = std::atan2(seed_point.r, pred_z);
        seed_hit.pred_x = pred_x;
        seed_hit.pred_y = pred_y;
        seed_hit.dphi = wrapPhi(seed_point.phi - pred_phi);
        seed_hit.dtheta = pointTheta(seed_point) - seed_hit.pred_theta;
        seed_hit.rdphi = seed_point.r * seed_hit.dphi;
        seed_hit.dz = seed_point.z - pred_z;
        seed_hit.tpc_pred_phi = tpc_pred_phi;
        seed_hit.tpc_pred_z = tpc_pred_z;
        seed_hit.tpc_pred_theta = tpc_seed_theta;
        seed_hit.tpc_dphi = initial_dphi;
        seed_hit.tpc_dtheta = seed_theta - tpc_seed_theta;
        double si_ref_pred_phi = 0.0;
        double si_ref_pred_z = 0.0;
        double si_ref_pred_x = 0.0;
        double si_ref_pred_y = 0.0;
        if (predictAtRadius(chain.si_reference, seed_point.r, si_ref_pred_phi, si_ref_pred_z, si_ref_pred_x, si_ref_pred_y))
        {
          seed_hit.si_ref_pred_phi = si_ref_pred_phi;
          seed_hit.si_ref_pred_z = si_ref_pred_z;
          seed_hit.si_ref_dphi = wrapPhi(seed_point.phi - si_ref_pred_phi);
          seed_hit.si_ref_dtheta = pointTheta(seed_point) - std::atan2(seed_point.r, si_ref_pred_z);
        }
        seed_hit.dynamic_mean_phi = dynamicMeanPhi(seed_layer, chain.previous_dphi, chain.has_previous_residual);
        seed_hit.dynamic_mean_theta = dynamicMeanTheta(seed_layer, chain.previous_dtheta, chain.has_previous_residual);
        seed_hit.sigma_phi = dynamicSigmaPhi(chain.pt);
        seed_hit.sigma_theta = dynamicSigmaTheta(chain.pt, seed_hit.pred_theta);
        seed_hit.sdphi = (seed_hit.dphi - seed_hit.dynamic_mean_phi) / seed_hit.sigma_phi;
        seed_hit.sdtheta = (seed_hit.dtheta - seed_hit.dynamic_mean_theta) / seed_hit.sigma_theta;
        seed_hit.chi2 = square(seed_hit.sdphi) + square(seed_hit.sdtheta);
      }

      chain = extendWithHit(chain, seed_hit);
      chain_used.insert(seed_hit.point.key);

      for (const unsigned int layer : m_matchLayers)
      {
        if (layer >= seed_layer)
        {
          continue;
        }

        ChainHit best_hit;
        double best_chi2 = 0.0;
        const SpacePoint* best_point = findBestMvtxCandidate(chain, silicon_points, chain_used,
                                                             layer, best_hit, best_chi2);
        if (best_point)
        {
          fillQaResiduals(track, chain, best_hit);
          chain = extendWithHit(chain, best_hit);
          chain_used.insert(best_hit.point.key);

          if (chain.hits.size() >= 2)
          {
            const ChainHit& outer_hit = chain.hits[chain.hits.size() - 2];
            const ChainHit& inner_hit = chain.hits.back();
            TrajectoryState corrected = correctSiliconSeedWithTwoHits(chain.si_reference, outer_hit, inner_hit);
            if (corrected.valid)
            {
              chain.si_reference = corrected;
              chain.state = corrected;
            }
          }
        }
        else
        {
          chain = extendMissing(chain, layer);
        }
      }

      if (!computeChainDcaMetrics(chain, tpc_seed.state))
      {
        continue;
      }
      chains.push_back(chain);
      if (chains.size() >= m_maxChains)
      {
        return chains;
      }
    }
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
    if (!std::isfinite(chain.dca_score))
    {
      continue;
    }
    if (!best ||
        chain.hits.size() > best->hits.size() ||
        (chain.hits.size() == best->hits.size() && chain.dca_score < best->dca_score))
    {
      best = &chain;
    }
  }
  return best;
}

Full_PolyTrackMatcher::Chain Full_PolyTrackMatcher::attachClosestInttClusters(
    const Tpc_PolyTrack& track,
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
    ChainHit best_hit;
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
      if (!predictAtRadius(out.si_reference, point.r, pred_phi, pred_z, pred_x, pred_y))
      {
        continue;
      }

      const double dphi = wrapPhi(point.phi - pred_phi);
      const double rdphi = point.r * dphi;
      const double dz = point.z - pred_z;
      if (std::fabs(rdphi) > m_looseRdphiWindow || std::fabs(dz) > m_inttDzWindow)
      {
        continue;
      }

      const double pred_theta = std::atan2(point.r, pred_z);
      const double dtheta = pointTheta(point) - pred_theta;
      const double mean_phi = dynamicMeanPhi(layer, out.previous_dphi, out.has_previous_residual);
      const double mean_theta = dynamicMeanTheta(layer, out.previous_dtheta, out.has_previous_residual);
      const double sigma_phi = dynamicSigmaPhi(out.pt);
      const double sigma_theta = dynamicSigmaTheta(out.pt, pred_theta);
      const double sdphi = (dphi - mean_phi) / sigma_phi;
      const double sdtheta = (dtheta - mean_theta) / sigma_theta;
      if (!m_associationCalibrationMode &&
          (std::fabs(sdphi) > m_phiWindowSigma || std::fabs(sdtheta) > m_thetaWindowSigma))
      {
        continue;
      }

      const double chi2 = square(sdphi) + square(sdtheta);
      if (chi2 < best_chi2)
      {
        best_chi2 = chi2;
        best_point = &point;
        best_hit.point = point;
        best_hit.pred_phi = pred_phi;
        best_hit.pred_z = pred_z;
        best_hit.pred_theta = pred_theta;
        best_hit.pred_x = pred_x;
        best_hit.pred_y = pred_y;
        best_hit.dphi = dphi;
        best_hit.dtheta = dtheta;
        best_hit.rdphi = rdphi;
        best_hit.dz = dz;

        double tpc_pred_phi = 0.0;
        double tpc_pred_z = 0.0;
        double tpc_pred_x = 0.0;
        double tpc_pred_y = 0.0;
        if (predictAtRadius(out.tpc_reference, point.r, tpc_pred_phi, tpc_pred_z, tpc_pred_x, tpc_pred_y))
        {
          best_hit.tpc_pred_phi = tpc_pred_phi;
          best_hit.tpc_pred_z = tpc_pred_z;
          best_hit.tpc_pred_theta = std::atan2(point.r, tpc_pred_z);
          best_hit.tpc_dphi = wrapPhi(point.phi - tpc_pred_phi);
          best_hit.tpc_dtheta = pointTheta(point) - best_hit.tpc_pred_theta;
        }

        double si_ref_pred_phi = 0.0;
        double si_ref_pred_z = 0.0;
        double si_ref_pred_x = 0.0;
        double si_ref_pred_y = 0.0;
        if (predictAtRadius(out.si_reference, point.r, si_ref_pred_phi, si_ref_pred_z, si_ref_pred_x, si_ref_pred_y))
        {
          best_hit.si_ref_pred_phi = si_ref_pred_phi;
          best_hit.si_ref_pred_z = si_ref_pred_z;
          best_hit.si_ref_dphi = wrapPhi(point.phi - si_ref_pred_phi);
          best_hit.si_ref_dtheta = pointTheta(point) - std::atan2(point.r, si_ref_pred_z);
        }
        best_hit.dynamic_mean_phi = mean_phi;
        best_hit.dynamic_mean_theta = mean_theta;
        best_hit.sigma_phi = sigma_phi;
        best_hit.sigma_theta = sigma_theta;
        best_hit.sdphi = sdphi;
        best_hit.sdtheta = sdtheta;
        best_hit.chi2 = chi2;
      }
    }

    if (best_point)
    {
      fillQaResiduals(track, out, best_hit);
      out = extendWithHit(out, best_hit);
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
    predictAtRadius(chain.state, r0, phi, z, x, y);
    out->set_x(x);
    out->set_y(y);
    out->set_z(z);
    const double next_r = r0 + 1.0;
    double phi2 = 0.0;
    double z2 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    predictAtRadius(chain.state, next_r, phi2, z2, x2, y2);
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

  if (m_writeQA && m_qaFile && chain.hits.size() > 0U &&
      std::isfinite(chain.dca_score) && chain.dca_score < std::numeric_limits<double>::max())
  {
    if (!m_hChainDcaScore)
    {
      m_hChainDcaScore = new TH1F("h_chain_dca_score", "selected chain DCA score;DCA score;tracks", 200, 0.0, 200.0);
      m_hChainDcaScore->SetDirectory(nullptr);
      m_hChainDeltaEta0 = new TH1F("h_chain_delta_eta0", "selected chain #Delta#eta_{0};#Delta#eta_{0};tracks", 120, 0.0, 0.3);
      m_hChainDeltaEta0->SetDirectory(nullptr);
    }
    m_hChainDcaScore->Fill(chain.dca_score);
    m_hChainDeltaEta0->Fill(chain.delta_eta0);

    if (std::isfinite(chain.signed_delta_phi0) && std::isfinite(chain.signed_delta_z0) &&
        std::isfinite(chain.signed_delta_theta0) && std::isfinite(chain.pt))
    {
      const std::string charge_label = chain.charge < 0.0 ? "neg" : "pos";
      const auto fill_signed_profile = [&](const std::string& name_base, const std::string& title,
                                           double value, double y_min, double y_max) {
        const std::string key = name_base + "_" + charge_label;
        TProfile* prof = nullptr;
        auto found = m_hChainSignedResidualVsPt.find(key);
        if (found != m_hChainSignedResidualVsPt.end())
        {
          prof = found->second;
        }
        else
        {
          const std::string hist_name = "h_" + key + "_vs_pt";
          const std::string hist_title = title + " (q " + (chain.charge < 0.0 ? "-" : "+") +
                                         ");p_{T} [GeV];" + title;
          prof = new TProfile(hist_name.c_str(), hist_title.c_str(),
                              static_cast<int>(m_qaPtBins), m_qaPtMin, m_qaPtMax, y_min, y_max);
          prof->SetDirectory(nullptr);
          m_hChainSignedResidualVsPt[key] = prof;
        }
        prof->Fill(chain.pt, value);
      };

      fill_signed_profile("chain_signed_dphi0", "signed #Delta#phi_{0} (Si-TPC) [rad]",
                          chain.signed_delta_phi0, -0.3, 0.3);
      fill_signed_profile("chain_signed_dz0", "signed #Deltaz_{0} (Si-TPC) [cm]",
                          chain.signed_delta_z0, -2.0, 2.0);
      fill_signed_profile("chain_signed_dtheta0", "signed #Delta#theta_{0} (Si-TPC) [rad]",
                          chain.signed_delta_theta0, -0.3, 0.3);
    }
  }

  for (const ChainHit& hit : chain.hits)
  {
    out->add_cluster_key(hit.point.key);
    out->add_silicon_state(hit.point.layer, hit.point.key,
                           finite_float(hit.point.x), finite_float(hit.point.y), finite_float(hit.point.z),
                           finite_float(hit.pred_x), finite_float(hit.pred_y), finite_float(hit.pred_z),
                           finite_float(hit.rdphi), finite_float(hit.dz), finite_float(hit.chi2));
  }

  m_fullTracks->add_track(out);
}

bool Full_PolyTrackMatcher::passesTpcAssociationCuts(const Tpc_PolyTrack& tpc_track) const
{
  const double pt = std::hypot(tpc_track.get_px(), tpc_track.get_py());
  return tpc_track.get_nclusters() > m_tpcAssociationClusterCut &&
         std::isfinite(pt) && pt > m_tpcAssociationPtCut;
}

int Full_PolyTrackMatcher::process_event(PHCompositeNode* topNode)
{
  if (!m_tpcTracks || !m_trkrClusters || !m_actsGeometry || !m_fullTracks)
  {
    if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  m_fullTracks->Reset();
  const std::vector<SpacePoint> silicon_points = collectSiliconClusters();
  for (unsigned int i = 0; i < m_tpcTracks->size(); ++i)
  {
    const Tpc_PolyTrack* tpc_track = m_tpcTracks->get_track(i);
    if (!tpc_track)
    {
      continue;
    }
    if (!passesTpcAssociationCuts(*tpc_track))
    {
      Chain unmatched_chain;
      fillTrack(*tpc_track, unmatched_chain);
      continue;
    }

    std::vector<Chain> chains = buildChains(*tpc_track, silicon_points);
    const Chain* best = selectBestChain(chains);
    Chain unmatched_chain;
    const Chain* selected_chain = best;
    if (!selected_chain)
    {
      selected_chain = chains.empty() ? &unmatched_chain : &chains.front();
    }

    Chain output_chain = *selected_chain;
    if (best)
    {
      output_chain = attachClosestInttClusters(*tpc_track, *best, silicon_points);
      fillQaCrossing(*tpc_track, output_chain);
      refitSiliconTrajectoryFromHits(output_chain);
    }
    fillTrack(*tpc_track, output_chain);
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
    for (auto& item : m_hResidualCurrentVsPreviousVsPt)
    {
      if (item.second)
      {
        item.second->Write();
      }
    }
    for (auto& item : m_hResidualSiVsTpc)
    {
      if (item.second)
      {
        item.second->Write();
      }
    }
    for (auto& item : m_hStandardResidual)
    {
      if (item.second)
      {
        item.second->Write();
      }
    }
    for (auto& item : m_hChainSignedResidualVsPt)
    {
      if (item.second)
      {
        item.second->Write();
      }
    }
    if (m_hInttCrossingVsTpcCrossing)
    {
      m_hInttCrossingVsTpcCrossing->Write();
    }
    if (m_hDeltaCrossing)
    {
      m_hDeltaCrossing->Write();
    }
    if (m_hCrossingStatus)
    {
      m_hCrossingStatus->Write();
    }
    if (m_hChainDcaScore)
    {
      m_hChainDcaScore->Write();
    }
    if (m_hChainDeltaEta0)
    {
      m_hChainDeltaEta0->Write();
    }

    std::cout << Name() << "::End - association calibration summary" << std::endl;
    for (unsigned int layer = 0; layer < m_dynamicPhiMeanOffset.size(); ++layer)
    {
      double n_phi = 0.0;
      double n_phi_2 = 0.0;
      double n_phi_25 = 0.0;
      double n_phi_3 = 0.0;
      double n_theta = 0.0;
      double n_theta_2 = 0.0;
      double n_theta_25 = 0.0;
      double n_theta_3 = 0.0;

      for (const auto& item : m_hStandardResidual)
      {
        if (!item.second)
        {
          continue;
        }
        const std::string layer_token = "_l" + std::to_string(layer) + "_";
        if (item.first.find(layer_token) == std::string::npos)
        {
          continue;
        }
        TH1F* hist = item.second;
        for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
        {
          const double entries = hist->GetBinContent(bin);
          const double center = hist->GetBinCenter(bin);
          if (item.first.find("h_si_sdphi") == 0)
          {
            n_phi += entries;
            if (std::fabs(center) < 2.0) { n_phi_2 += entries; }
            if (std::fabs(center) < 2.5) { n_phi_25 += entries; }
            if (std::fabs(center) < 3.0) { n_phi_3 += entries; }
          }
          else if (item.first.find("h_si_sdtheta") == 0)
          {
            n_theta += entries;
            if (std::fabs(center) < 2.0) { n_theta_2 += entries; }
            if (std::fabs(center) < 2.5) { n_theta_25 += entries; }
            if (std::fabs(center) < 3.0) { n_theta_3 += entries; }
          }
        }
      }

      const double eff2 = n_phi > 0.0 ? n_phi_2 / n_phi : 0.0;
      const double eff25 = n_phi > 0.0 ? n_phi_25 / n_phi : 0.0;
      const double eff3 = n_phi > 0.0 ? n_phi_3 / n_phi : 0.0;
      const double theta_eff2 = n_theta > 0.0 ? n_theta_2 / n_theta : 0.0;
      const double theta_eff25 = n_theta > 0.0 ? n_theta_25 / n_theta : 0.0;
      const double theta_eff3 = n_theta > 0.0 ? n_theta_3 / n_theta : 0.0;

      std::cout << "Layer " << layer << ":\n"
                << "  phi mean offset " << m_dynamicPhiMeanOffset[layer] << "\n"
                << "  phi previous-layer slope " << m_dynamicPhiMeanSlope[layer] << "\n"
                << "  phi sigma " << m_sigmaPhi << "\n"
                << "  theta mean offset " << m_dynamicThetaMeanOffset[layer] << "\n"
                << "  theta previous-layer slope " << m_dynamicThetaMeanSlope[layer] << "\n"
                << "  theta sigma " << m_sigmaTheta << "\n"
                << "  efficiency within 2 sigma " << eff2 << " phi, " << theta_eff2 << " theta\n"
                << "  efficiency within 2.5 sigma " << eff25 << " phi, " << theta_eff25 << " theta\n"
                << "  efficiency within 3 sigma " << eff3 << " phi, " << theta_eff3 << " theta" << std::endl;
    }
    m_qaFile->Close();
    delete m_qaFile;
    m_qaFile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}