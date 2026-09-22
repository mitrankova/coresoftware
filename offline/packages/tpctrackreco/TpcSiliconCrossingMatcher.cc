#include "TpcSiliconCrossingMatcher.h"

#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
#include "TpcSiliconMatchCandidate.h"
#include "TpcSiliconMatchCandidateContainer.h"
#include "Tpc_FittingTools.h"
#include "TpcTrackKalmanFitter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phfield/PHFieldUtility.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>

namespace
{
  double square(const double value) { return value * value; }
  constexpr double SeedPhiWindow = 0.15;
  constexpr double SeedThetaPreWindow = 4.0;

  double seedEtaWindow(const double pt)
  {
    return std::max(0.03, -0.014 + 0.0331 * std::exp(0.48 / pt));
  }

  double seedThetaWindow(const double pt, const double theta)
  {
    return seedEtaWindow(pt) * std::max(std::sin(theta), 1.e-9);
  }

  double thetaToEta(const double theta)
  {
    const double halfTan = std::tan(0.5 * theta);
    return std::isfinite(halfTan) && halfTan > 0. ? -std::log(halfTan)
                                                 : std::numeric_limits<double>::quiet_NaN();
  }

  double sagittaModelDerivative(const double xrot, const double x0, const double invR)
  {
    const double dx = xrot - x0;
    const double dx2 = dx * dx;
    const double invR2 = invR * invR;
    const double invR3 = invR2 * invR;
    const double invR5 = invR3 * invR2;
    return -invR * dx - 0.5 * invR3 * dx2 * dx - 0.375 * invR5 * dx2 * dx2 * dx;
  }
}

TpcSiliconCrossingMatcher::TpcSiliconCrossingMatcher(const std::string& name) : SubsysReco(name) {}

int TpcSiliconCrossingMatcher::getNodes(PHCompositeNode* topNode)
{
  m_trajectories = findNode::getClass<TpcCrossingTrajectoryContainer>(topNode, m_trajectoryNodeName);
  m_clusters = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterNodeName);
  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_trajectories || !m_clusters || !m_geometry)
  {
    std::cerr << Name() << "::getNodes - missing trajectory, cluster, or ActsGeometry input" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcSiliconCrossingMatcher::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dst = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst) { dst = new PHCompositeNode("DST"); topNode->addNode(dst); }
  m_candidates = findNode::getClass<TpcSiliconMatchCandidateContainer>(topNode, m_outputNodeName);
  if (!m_candidates)
  {
    m_candidates = new TpcSiliconMatchCandidateContainer;
    dst->addNode(new PHIODataNode<PHObject>(m_candidates, m_outputNodeName, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcSiliconCrossingMatcher::InitRun(PHCompositeNode* topNode)
{
  if (!m_beamFrame.validate() || !(m_zSearchTimeBins > 0.) || !(m_tpcAdcClockNs > 0.))
  {
    std::cerr << Name() << "::InitRun - invalid beam frame or TPC time-bin search configuration" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  m_propagationConfig.magnetic_field = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
  m_propagationConfig.analytic_uniform_propagation = false;
  if (!m_propagationConfig.magnetic_field)
  {
    std::cerr << Name() << "::InitRun - missing magnetic field" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return createNodes(topNode);
}

double TpcSiliconCrossingMatcher::wrapPhi(double phi) const
{
  while (phi > M_PI) phi -= 2. * M_PI;
  while (phi <= -M_PI) phi += 2. * M_PI;
  return phi;
}

double TpcSiliconCrossingMatcher::unwrapPhiNear(double phi, const double reference) const
{
  while (phi - reference > M_PI) phi -= 2. * M_PI;
  while (phi - reference <= -M_PI) phi += 2. * M_PI;
  return phi;
}

double TpcSiliconCrossingMatcher::zSearchWindowCm() const
{
  return m_zSearchTimeBins * m_tpcAdcClockNs * m_geometry->get_drift_velocity();
}

std::vector<TpcSiliconCrossingMatcher::SpacePoint> TpcSiliconCrossingMatcher::collectSiliconClusters() const
{
  std::vector<SpacePoint> points;
  for (const auto detector : {TrkrDefs::TrkrId::mvtxId, TrkrDefs::TrkrId::inttId})
  {
    for (const auto hitsetkey : m_clusters->getHitSetKeys(detector))
    {
      const auto range = m_clusters->getClusters(hitsetkey);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        if (!iter->second) continue;
        const unsigned int layer = TrkrDefs::getLayer(iter->first);
        if (layer >= 7U) continue;
        const auto global = m_geometry->getGlobalPosition(iter->first, iter->second);
        const auto& beam = detector == TrkrDefs::mvtxId ? m_beamFrame.mvtxBeamLine()
                                                       : m_beamFrame.inttBeamLine();
        const auto centered = m_beamFrame.toBeamFrame(beam, global.x(), global.y(), global.z());
        SpacePoint point;
        point.key = iter->first;
        point.layer = layer;
        point.global_x = global.x();
        point.global_y = global.y();
        point.global_z = global.z();
        point.x = centered.x;
        point.y = centered.y;
        point.z = centered.z;
        point.r = std::hypot(point.x, point.y);
        point.phi = std::atan2(point.y, point.x);
        if (std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) && point.r > 0.)
          points.push_back(point);
      }
    }
  }
  return points;
}

TpcSiliconCrossingMatcher::TrajectoryState TpcSiliconCrossingMatcher::fitTrajectory(
    const std::vector<SpacePoint>& points) const
{
  TrajectoryState state;
  if (points.size() < 2U) return state;
  std::vector<Tpc_FittingTools::FitPoint> phiPoints;
  std::vector<Tpc_FittingTools::FitPoint> zPoints;
  phiPoints.reserve(points.size());
  zPoints.reserve(points.size());
  const double referencePhi = points.front().phi;
  for (const auto& point : points)
  {
    phiPoints.emplace_back(point.r, unwrapPhiNear(point.phi, referencePhi));
    zPoints.emplace_back(point.r, point.z);
  }
  const auto phiLine = Tpc_FittingTools::fitLine(phiPoints);
  const auto zLine = Tpc_FittingTools::fitLine(zPoints);
  if (!phiLine.ok || !zLine.ok) return state;
  state.phi_slope = phiLine.slope;
  state.phi_intercept = phiLine.intercept;
  state.phi_theta = std::atan(phiLine.slope);
  state.phi_bline = phiLine.intercept;
  if (points.size() >= 3U)
  {
    const auto sagitta = Tpc_FittingTools::fitSagitta(phiPoints);
    if (sagitta.ok && std::isfinite(sagitta.S) && std::isfinite(sagitta.x0) &&
        std::isfinite(sagitta.invR) && std::isfinite(sagitta.theta) && std::isfinite(sagitta.b))
    {
      state.phi_S = sagitta.S;
      state.phi_x0 = sagitta.x0;
      state.phi_invR = sagitta.invR;
      state.phi_theta = sagitta.theta;
      state.phi_bline = sagitta.b;
      state.phi_sagitta_ok = true;
    }
  }
  state.z_slope = zLine.slope;
  state.z_intercept = zLine.intercept;
  state.seed_slope = zLine.slope;
  state.seed_q_over_r = state.phi_sagitta_ok ? state.phi_invR : 0.;
  state.valid = std::isfinite(state.phi_slope) && std::isfinite(state.phi_intercept) &&
                std::isfinite(state.z_slope) && std::isfinite(state.z_intercept);
  return state;
}

TpcSiliconCrossingMatcher::TrajectoryState TpcSiliconCrossingMatcher::makeTpcReferenceTrajectory(
    const TpcCrossingTrajectory& trajectory) const
{
  std::vector<SpacePoint> points;
  for (unsigned int i = 0; i < trajectory.size_layer_states(); ++i)
  {
    const auto* layer = trajectory.get_layer_state(i);
    if (!layer || !layer->valid) continue;
    const auto centered = m_beamFrame.toBeamFrame(m_beamFrame.tpcBeamLine(), layer->x, layer->y, layer->z);
    SpacePoint point;
    point.layer = layer->layer;
    point.x = centered.x;
    point.y = centered.y;
    point.z = centered.z;
    point.r = std::hypot(point.x, point.y);
    point.phi = std::atan2(point.y, point.x);
    if (std::isfinite(point.z) && point.r > 0.) points.push_back(point);
  }
  std::sort(points.begin(), points.end(), [](const auto& lhs, const auto& rhs) { return lhs.r < rhs.r; });
  return fitTrajectory(points);
}

double TpcSiliconCrossingMatcher::predictSagittaPhi(const double r, const TrajectoryState& state) const
{
  const double c = std::cos(state.phi_theta);
  const double s = std::sin(state.phi_theta);
  double yy = std::tan(state.phi_theta) * r;
  for (unsigned int iter = 0; iter < 25U; ++iter)
  {
    const double xrot = c * r + s * yy;
    const double yrot = -s * r + c * yy;
    const double f = Tpc_FittingTools::sagittaModel(xrot, state.phi_S, state.phi_x0, state.phi_invR);
    const double derivative = c - sagittaModelDerivative(xrot, state.phi_x0, state.phi_invR) * s;
    if (std::abs(derivative) < 1.e-12) break;
    const double step = (yrot - f) / derivative;
    yy -= step;
    if (std::abs(step) < 1.e-10) break;
  }
  return state.phi_bline + yy;
}

bool TpcSiliconCrossingMatcher::predictAtRadius(const TrajectoryState& state, const double r,
                                                 double& predPhi, double& predZ,
                                                 double& predX, double& predY) const
{
  if (!state.valid || !std::isfinite(r) || r <= 0.) return false;
  if (!state.use_tpc_seed)
  {
    predPhi = wrapPhi(state.phi_sagitta_ok ? predictSagittaPhi(r, state)
                                          : state.phi_intercept + state.phi_slope * r);
    predZ = state.z_intercept + state.z_slope * r;
    predX = r * std::cos(predPhi);
    predY = r * std::sin(predPhi);
    return std::isfinite(predPhi) && std::isfinite(predZ);
  }

  double arc = std::numeric_limits<double>::quiet_NaN();
  if (std::abs(state.seed_q_over_r) <= 1.e-12)
  {
    const double tx = std::cos(state.seed_phi0);
    const double ty = std::sin(state.seed_phi0);
    const double b = state.seed_x0 * tx + state.seed_y0 * ty;
    const double c = square(state.seed_x0) + square(state.seed_y0) - square(r);
    const double discriminant = b * b - c;
    if (discriminant < 0.) return false;
    const double root = std::sqrt(discriminant);
    const double s1 = -b + root;
    const double s2 = -b - root;
    if (s1 >= 0. && s2 >= 0.) arc = std::min(s1, s2);
    else if (s1 >= 0.) arc = s1;
    else if (s2 >= 0.) arc = s2;
    else arc = std::abs(s1) < std::abs(s2) ? s1 : s2;
    predX = state.seed_x0 + arc * tx;
    predY = state.seed_y0 + arc * ty;
  }
  else
  {
    const double radius = 1. / std::abs(state.seed_q_over_r);
    const double dc = std::hypot(state.seed_cx, state.seed_cy);
    if (!std::isfinite(radius) || radius <= 0. || dc <= 1.e-12) return false;
    const double a = (r * r - radius * radius + dc * dc) / (2. * dc);
    const double h2 = r * r - a * a;
    if (h2 < -1.e-9) return false;
    const double h = std::sqrt(std::max(h2, 0.));
    const double ux = state.seed_cx / dc;
    const double uy = state.seed_cy / dc;
    const double baseX = a * ux;
    const double baseY = a * uy;
    const double candidateX[2] = {baseX - h * uy, baseX + h * uy};
    const double candidateY[2] = {baseY + h * ux, baseY - h * ux};
    const double startAngle = std::atan2(state.seed_y0 - state.seed_cy, state.seed_x0 - state.seed_cx);
    const double fitSign = state.seed_q_over_r > 0. ? -1. : 1.;
    double bestArc = std::numeric_limits<double>::max();
    unsigned int best = 0;
    for (unsigned int i = 0; i < 2U; ++i)
    {
      double delta = unwrapPhiNear(std::atan2(candidateY[i] - state.seed_cy,
                                             candidateX[i] - state.seed_cx), startAngle) - startAngle;
      if (fitSign * delta < 0.) delta += fitSign * 2. * M_PI;
      const double candidateArc = std::abs(radius * delta);
      if (candidateArc < bestArc) { bestArc = candidateArc; best = i; }
    }
    arc = bestArc;
    predX = candidateX[best];
    predY = candidateY[best];
  }
  predPhi = std::atan2(predY, predX);
  predZ = state.seed_z0 + state.seed_slope * arc;
  return std::isfinite(predPhi) && std::isfinite(predZ) && std::isfinite(predX) && std::isfinite(predY);
}

TpcSiliconCrossingMatcher::TrajectoryState TpcSiliconCrossingMatcher::makeSiliconSeedTrajectory(
    const TrajectoryState& tpcState, const SpacePoint& firstMvtx) const
{
  TrajectoryState state;
  const double qOverR = tpcState.seed_q_over_r;
  const double chord = std::hypot(firstMvtx.x, firstMvtx.y);
  if (!std::isfinite(chord) || chord <= 1.e-9 || !std::isfinite(qOverR)) return state;
  state.seed_q_over_r = qOverR;
  state.seed_slope = tpcState.seed_slope;
  state.use_tpc_seed = true;
  state.use_silicon_seed = true;
  if (std::abs(qOverR) <= 1.e-12)
  {
    state.seed_phi0 = std::atan2(firstMvtx.y, firstMvtx.x);
    state.seed_z0 = firstMvtx.z - state.seed_slope * chord;
    state.valid = std::isfinite(state.seed_phi0) && std::isfinite(state.seed_z0) && std::isfinite(state.seed_slope);
    return state;
  }
  const double radius = 1. / std::abs(qOverR);
  if (!std::isfinite(radius) || chord > 2. * radius) return state;
  const double mx = 0.5 * firstMvtx.x;
  const double my = 0.5 * firstMvtx.y;
  const double ux = firstMvtx.x / chord;
  const double uy = firstMvtx.y / chord;
  const double h = std::sqrt(std::max(radius * radius - 0.25 * chord * chord, 0.));
  const double fitSign = qOverR > 0. ? -1. : 1.;
  double bestArc = std::numeric_limits<double>::max();
  TrajectoryState bestState;
  for (const double side : {-1., 1.})
  {
    auto candidate = state;
    candidate.seed_cx = mx + side * h * (-uy);
    candidate.seed_cy = my + side * h * ux;
    const double start = std::atan2(-candidate.seed_cy, -candidate.seed_cx);
    double delta = unwrapPhiNear(std::atan2(firstMvtx.y - candidate.seed_cy,
                                           firstMvtx.x - candidate.seed_cx), start) - start;
    if (fitSign * delta < 0.) delta += fitSign * 2. * M_PI;
    const double arc = std::abs(radius * delta);
    candidate.seed_phi0 = start + fitSign * 0.5 * M_PI;
    candidate.seed_z0 = firstMvtx.z - candidate.seed_slope * arc;
    candidate.valid = std::isfinite(candidate.seed_cx) && std::isfinite(candidate.seed_cy) &&
                      std::isfinite(candidate.seed_phi0) && std::isfinite(candidate.seed_z0);
    if (candidate.valid && arc < bestArc) { bestArc = arc; bestState = candidate; }
  }
  return bestState;
}

TpcSiliconCrossingMatcher::TrajectoryState TpcSiliconCrossingMatcher::correctSiliconSeedWithTwoHits(
    const TrajectoryState& seed, const ChainHit& outer, const ChainHit& inner) const
{
  TrajectoryState state = seed;
  if (!state.valid || std::abs(state.seed_q_over_r) <= 1.e-12) return state;
  const double radius = 1. / std::abs(state.seed_q_over_r);
  const double dx = inner.point.x - outer.point.x;
  const double dy = inner.point.y - outer.point.y;
  const double chord = std::hypot(dx, dy);
  if (!std::isfinite(chord) || chord <= 1.e-9 || chord > 2. * radius) return state;
  const double mx = 0.5 * (outer.point.x + inner.point.x);
  const double my = 0.5 * (outer.point.y + inner.point.y);
  const double h = std::sqrt(std::max(radius * radius - 0.25 * chord * chord, 0.));
  double bestD2 = std::numeric_limits<double>::max();
  double cx = 0., cy = 0.;
  for (const double side : {-1., 1.})
  {
    const double candidateCx = mx + side * h * (-dy / chord);
    const double candidateCy = my + side * h * (dx / chord);
    const double d2 = square(candidateCx - state.seed_cx) + square(candidateCy - state.seed_cy);
    if (d2 < bestD2) { bestD2 = d2; cx = candidateCx; cy = candidateCy; }
  }
  const double length = std::hypot(cx - state.seed_x0, cy - state.seed_y0);
  if (!std::isfinite(length) || length <= 1.e-9) return state;
  const double dca = length - radius;
  state.seed_x0 += (cx - state.seed_x0) / length * dca;
  state.seed_y0 += (cy - state.seed_y0) / length * dca;
  state.seed_cx = cx;
  state.seed_cy = cy;
  const double fitSign = state.seed_q_over_r > 0. ? -1. : 1.;
  state.seed_phi0 = std::atan2(state.seed_y0 - cy, state.seed_x0 - cx) + fitSign * 0.5 * M_PI;
  const double start = std::atan2(state.seed_y0 - cy, state.seed_x0 - cx);
  const auto arcTo = [&](const SpacePoint& point)
  {
    double delta = unwrapPhiNear(std::atan2(point.y - cy, point.x - cx), start) - start;
    if (fitSign * delta < 0.) delta += fitSign * 2. * M_PI;
    return std::abs(radius * delta);
  };
  const double sOuter = arcTo(outer.point);
  const double sInner = arcTo(inner.point);
  if (std::isfinite(sOuter) && std::isfinite(sInner) && std::abs(sOuter - sInner) > 1.e-6)
  {
    state.seed_slope = (outer.point.z - inner.point.z) / (sOuter - sInner);
    state.seed_z0 = inner.point.z - state.seed_slope * sInner;
  }
  state.valid = std::isfinite(state.seed_phi0) && std::isfinite(state.seed_slope) && std::isfinite(state.seed_z0);
  return state;
}

TpcSiliconCrossingMatcher::TrajectoryState TpcSiliconCrossingMatcher::correctSiliconSeedWithAllHits(
    const TrajectoryState& seed, const std::vector<ChainHit>& hits) const
{
  TrajectoryState state = seed;
  if (!state.valid || std::abs(state.seed_q_over_r) <= 1.e-12 || hits.size() < 2U) return state;
  const double radius = 1. / std::abs(state.seed_q_over_r);
  double sxx = 0., sxy = 0., syy = 0., sxb = 0., syb = 0.;
  for (const auto& hit : hits)
  {
    const double dx = hit.point.x - state.seed_cx;
    const double dy = hit.point.y - state.seed_cy;
    const double distance = std::hypot(dx, dy);
    if (!std::isfinite(distance) || distance <= 1.e-9) continue;
    const double ux = dx / distance;
    const double uy = dy / distance;
    const double residual = distance - radius;
    sxx += ux * ux; sxy += ux * uy; syy += uy * uy;
    sxb += ux * residual; syb += uy * residual;
  }
  const double determinant = sxx * syy - sxy * sxy;
  if (!std::isfinite(determinant) || std::abs(determinant) <= 1.e-12) return state;
  const double cx = state.seed_cx + (sxb * syy - syb * sxy) / determinant;
  const double cy = state.seed_cy + (sxx * syb - sxy * sxb) / determinant;
  const double length = std::hypot(cx - state.seed_x0, cy - state.seed_y0);
  if (!std::isfinite(length) || length <= 1.e-9) return state;
  const double dca = length - radius;
  state.seed_x0 += (cx - state.seed_x0) / length * dca;
  state.seed_y0 += (cy - state.seed_y0) / length * dca;
  state.seed_cx = cx;
  state.seed_cy = cy;
  const double fitSign = state.seed_q_over_r > 0. ? -1. : 1.;
  state.seed_phi0 = std::atan2(state.seed_y0 - cy, state.seed_x0 - cx) + fitSign * 0.5 * M_PI;
  const double start = std::atan2(state.seed_y0 - cy, state.seed_x0 - cx);
  std::vector<double> arcs;
  arcs.reserve(hits.size());
  double arcMean = 0., zMean = 0.;
  for (const auto& hit : hits)
  {
    double delta = unwrapPhiNear(std::atan2(hit.point.y - cy, hit.point.x - cx), start) - start;
    if (fitSign * delta < 0.) delta += fitSign * 2. * M_PI;
    arcs.push_back(std::abs(radius * delta));
    arcMean += arcs.back();
    zMean += hit.point.z;
  }
  arcMean /= static_cast<double>(hits.size());
  zMean /= static_cast<double>(hits.size());
  double sss = 0., ssz = 0.;
  for (std::size_t i = 0; i < hits.size(); ++i)
  {
    sss += square(arcs[i] - arcMean);
    ssz += (arcs[i] - arcMean) * (hits[i].point.z - zMean);
  }
  if (sss > 1.e-6)
  {
    state.seed_slope = ssz / sss;
    state.seed_z0 = zMean - state.seed_slope * arcMean;
  }
  state.valid = std::isfinite(state.seed_phi0) && std::isfinite(state.seed_slope) && std::isfinite(state.seed_z0);
  return state;
}

double TpcSiliconCrossingMatcher::pointTheta(const SpacePoint& point) const
{
  return std::atan2(point.r, point.z);
}

double TpcSiliconCrossingMatcher::dynamicMeanPhi(const unsigned int layer, const double previous,
                                                  const bool hasPrevious) const
{
  if (!m_useDynamicResiduals || !hasPrevious || layer >= 7U) return 0.;
  return m_dynamicPhiMeanOffset[layer] + m_dynamicPhiMeanSlope[layer] * previous;
}

double TpcSiliconCrossingMatcher::dynamicSigmaPhi(const double) const
{
  return std::max(m_sigmaPhi, 1.e-9);
}

double TpcSiliconCrossingMatcher::dynamicMeanTheta(const unsigned int layer, const double previous,
                                                    const bool hasPrevious) const
{
  if (!m_useDynamicResiduals || !hasPrevious || layer >= 7U) return 0.;
  return m_dynamicThetaMeanOffset[layer] + m_dynamicThetaMeanSlope[layer] * previous;
}

double TpcSiliconCrossingMatcher::dynamicSigmaTheta(const double pt, const double predictedTheta) const
{
  return std::max(seedThetaWindow(pt, predictedTheta) / std::max(m_thetaWindowSigma, 1.e-9), 1.e-9);
}

double TpcSiliconCrossingMatcher::dynamicDzWindow(const double pt) const
{
  return 1.138 + 0.3919 * std::exp(0.84 / std::max(pt, 0.25));
}

bool TpcSiliconCrossingMatcher::matchToSurface(
    const TpcCrossingTrajectory& trajectory, const SpacePoint& point, SurfaceMatch& match) const
{
  auto* cluster = m_clusters->findCluster(point.key);
  const auto surface = cluster ? m_geometry->maps().getSurface(point.key, cluster) : nullptr;
  if (!cluster || !surface) return false;

  std::array<double, TpcTrackKalmanFitter::StateDim> startState{};
  for (unsigned int index = 0; index < TpcTrackKalmanFitter::StateDim; ++index)
    startState[index] = trajectory.get_state(index);
  std::array<double, TpcTrackKalmanFitter::StateDim> intersectionState{};
  const auto& context = m_geometry->geometry().geoContext;
  if (!TpcTrackKalmanFitter::propagate_to_surface(
          startState, m_propagationConfig, *surface, context, intersectionState, &match.path_length_cm))
    return false;

  match.intersection = {intersectionState[TpcTrackKalmanFitter::X],
                        intersectionState[TpcTrackKalmanFitter::Y],
                        intersectionState[TpcTrackKalmanFitter::Z]};
  const auto center = surface->center(context) / Acts::UnitConstants::cm;
  match.surface_center = {center.x(), center.y(), center.z()};

  const Acts::Vector2 local(cluster->getLocalX() * Acts::UnitConstants::cm,
                            cluster->getLocalY() * Acts::UnitConstants::cm);
  const auto tangent = TpcTrackKalmanFitter::state_tangent(intersectionState);
  const Acts::Vector3 direction(tangent.x, tangent.y, tangent.z);
  const auto origin = surface->localToGlobal(context, local, direction);
  const auto along0 = surface->localToGlobal(
      context, local + Acts::Vector2(Acts::UnitConstants::cm, 0.), direction) - origin;
  const auto along1 = surface->localToGlobal(
      context, local + Acts::Vector2(0., Acts::UnitConstants::cm), direction) - origin;
  if (!(along0.norm() > 0.) || !(along1.norm() > 0.)) return false;
  const auto axis0 = along0.normalized();
  const auto axis1 = along1.normalized();
  const Acts::Vector3 residual((point.global_x - match.intersection[0]) * Acts::UnitConstants::cm,
                               (point.global_y - match.intersection[1]) * Acts::UnitConstants::cm,
                               (point.global_z - match.intersection[2]) * Acts::UnitConstants::cm);
  match.local_residual_0 = residual.dot(axis0) / Acts::UnitConstants::cm;
  match.local_residual_1 = residual.dot(axis1) / Acts::UnitConstants::cm;
  return std::isfinite(match.local_residual_0) && std::isfinite(match.local_residual_1);
}

const TpcSiliconCrossingMatcher::SpacePoint* TpcSiliconCrossingMatcher::findBestMvtxCandidate(
    const TpcCrossingTrajectory& trajectory, const Chain& chain, const std::vector<SpacePoint>& points,
    const std::set<TrkrDefs::cluskey>& used, const unsigned int layer,
    ChainHit& bestHit, double& bestChi2, Counters& counters) const
{
  const SpacePoint* best = nullptr;
  bestChi2 = std::numeric_limits<double>::max();
  for (const auto& point : points)
  {
    if (point.layer != layer || TrkrDefs::getTrkrId(point.key) != TrkrDefs::mvtxId || used.count(point.key)) continue;
    ++counters.search[layer];
    double predPhi = 0., predZ = 0., predX = 0., predY = 0.;
    if (!predictAtRadius(chain.si_reference, point.r, predPhi, predZ, predX, predY)) continue;
    const double dz = point.z - predZ;
    if (std::abs(dz) > zSearchWindowCm()) continue;
    ++counters.zPass[layer];
    const double dphi = wrapPhi(point.phi - predPhi);
    const double rdphi = point.r * dphi;
    if (std::abs(rdphi) > m_looseRdphiWindow || std::abs(dz) > m_looseDzWindow) continue;
    const double predTheta = std::atan2(point.r, predZ);
    const double dtheta = pointTheta(point) - predTheta;
    const double meanPhi = dynamicMeanPhi(layer, chain.previous_dphi, chain.has_previous_residual);
    const double meanTheta = dynamicMeanTheta(layer, chain.previous_dtheta, chain.has_previous_residual);
    const double sigmaPhi = dynamicSigmaPhi(chain.pt);
    const double sigmaTheta = dynamicSigmaTheta(chain.pt, predTheta);
    double sdphi = (dphi - meanPhi) / sigmaPhi;
    double sdtheta = (dtheta - meanTheta) / sigmaTheta;
    if (!m_associationCalibrationMode && layer < 2U)
    {
      sdphi = (sdphi - m_vertexPhiMean[layer]) / std::max(m_vertexPhiSigma[layer], 1.e-9);
      sdtheta = (sdtheta - m_vertexThetaMean[layer]) / std::max(m_vertexThetaSigma[layer], 1.e-9);
    }
    if (!m_associationCalibrationMode &&
        (std::abs(sdphi) > m_phiWindowSigma || std::abs(sdtheta) > m_thetaWindowSigma)) continue;

    SurfaceMatch surfaceMatch;
    if (!matchToSurface(trajectory, point, surfaceMatch) ||
        std::abs(surfaceMatch.local_residual_0) > m_mvtxLocal0Window ||
        std::abs(surfaceMatch.local_residual_1) > m_mvtxLocal1Window) continue;
    ++counters.phiPass[layer];
    const double chi2 = square(surfaceMatch.local_residual_0 / m_mvtxLocal0Window) +
                        square(surfaceMatch.local_residual_1 / m_mvtxLocal1Window);
    if (chi2 < bestChi2)
    {
      bestChi2 = chi2;
      best = &point;
      bestHit.point = point;
      bestHit.pred_phi = predPhi;
      bestHit.pred_z = predZ;
      bestHit.dphi = dphi;
      bestHit.dtheta = dtheta;
      bestHit.rdphi = rdphi;
      bestHit.dz = dz;
      bestHit.ddphi = dphi - meanPhi;
      bestHit.chi2 = chi2;
      bestHit.surface_residual_0 = surfaceMatch.local_residual_0;
      bestHit.surface_residual_1 = surfaceMatch.local_residual_1;
    }
  }
  return best;
}

double TpcSiliconCrossingMatcher::trajectoryPhi0NearBeam(const TrajectoryState& state) const
{
  if (!state.valid) return std::numeric_limits<double>::quiet_NaN();
  if (state.use_silicon_seed && state.use_tpc_seed) return state.seed_phi0;
  double phi1 = 0., z1 = 0., x1 = 0., y1 = 0.;
  double phi2 = 0., z2 = 0., x2 = 0., y2 = 0.;
  if (predictAtRadius(state, 0.01, phi1, z1, x1, y1) &&
      predictAtRadius(state, 0.02, phi2, z2, x2, y2)) return std::atan2(y2 - y1, x2 - x1);
  return state.phi_bline;
}

double TpcSiliconCrossingMatcher::trajectoryZ0NearBeam(const TrajectoryState& state) const
{
  return state.use_tpc_seed ? state.seed_z0 : state.z_intercept;
}

double TpcSiliconCrossingMatcher::trajectoryTheta0NearBeam(const TrajectoryState& state) const
{
  return std::atan2(1., state.use_tpc_seed ? state.seed_slope : state.z_slope);
}

bool TpcSiliconCrossingMatcher::computeChainDcaMetrics(Chain& chain, const TrajectoryState& tpc) const
{
  const double tpcPhi = trajectoryPhi0NearBeam(tpc);
  const double siliconPhi = trajectoryPhi0NearBeam(chain.state);
  const double deltaPhi = std::min(std::abs(wrapPhi(siliconPhi - tpcPhi)),
                                   std::abs(wrapPhi(siliconPhi + M_PI - tpcPhi)));
  const double deltaZ = std::abs(trajectoryZ0NearBeam(chain.state) - trajectoryZ0NearBeam(tpc));
  const double tpcTheta = trajectoryTheta0NearBeam(tpc);
  const double siliconTheta = trajectoryTheta0NearBeam(chain.state);
  chain.delta_eta0 = std::abs(thetaToEta(siliconTheta) - thetaToEta(tpcTheta));
  const double thetaScale = seedThetaWindow(chain.pt, tpcTheta);
  const double zScale = dynamicDzWindow(chain.pt);
  if (!std::isfinite(deltaPhi) || !std::isfinite(deltaZ) || !std::isfinite(chain.delta_eta0) ||
      !(thetaScale > 0.) || !(zScale > 0.)) return false;
  chain.dca_score = square(deltaPhi / SeedPhiWindow) +
                    square((siliconTheta - tpcTheta) / thetaScale) + square(deltaZ / zScale);
  return std::isfinite(chain.dca_score);
}

std::vector<TpcSiliconCrossingMatcher::Chain> TpcSiliconCrossingMatcher::buildChains(
    const TpcCrossingTrajectory& trajectory, const std::vector<SpacePoint>& points, Counters& counters) const
{
  Chain tpcSeed;
  tpcSeed.pt = std::abs(trajectory.get_state(TpcCrossingTrajectory::QOverPt)) > 1.e-12
                   ? 1. / std::abs(trajectory.get_state(TpcCrossingTrajectory::QOverPt))
                   : 1.e12;
  tpcSeed.state = makeTpcReferenceTrajectory(trajectory);
  tpcSeed.tpc_reference = tpcSeed.state;
  if (!tpcSeed.state.valid) return {};
  std::vector<Chain> chains;
  for (const unsigned int seedLayer : m_matchLayers)
  {
    for (const auto& seedPoint : points)
    {
      if (seedPoint.layer != seedLayer || TrkrDefs::getTrkrId(seedPoint.key) != TrkrDefs::mvtxId) continue;
      ++counters.search[seedLayer];
      double tpcPhi = 0., tpcZ = 0., tpcX = 0., tpcY = 0.;
      if (!predictAtRadius(tpcSeed.state, seedPoint.r, tpcPhi, tpcZ, tpcX, tpcY)) continue;
      if (std::abs(seedPoint.z - tpcZ) > zSearchWindowCm()) continue;
      ++counters.zPass[seedLayer];
      const double initialDphi = wrapPhi(seedPoint.phi - tpcPhi);
      const double seedTheta = pointTheta(seedPoint);
      const double tpcSeedTheta = std::atan2(seedPoint.r, tpcZ);
      if (std::abs(initialDphi) > SeedPhiWindow + 0.15 ||
          std::abs(seedTheta - tpcSeedTheta) > SeedThetaPreWindow) continue;
      SurfaceMatch seedSurfaceMatch;
      if (!matchToSurface(trajectory, seedPoint, seedSurfaceMatch) ||
          std::abs(seedSurfaceMatch.local_residual_0) > m_mvtxLocal0Window ||
          std::abs(seedSurfaceMatch.local_residual_1) > m_mvtxLocal1Window) continue;
      ++counters.phiPass[seedLayer];
      Chain chain;
      chain.pt = tpcSeed.pt;
      chain.tpc_reference = tpcSeed.state;
      chain.state = makeSiliconSeedTrajectory(tpcSeed.state, seedPoint);
      if (!chain.state.valid) continue;
      chain.si_reference = chain.state;
      ChainHit seedHit;
      seedHit.point = seedPoint;
      seedHit.pred_phi = tpcPhi;
      seedHit.pred_z = tpcZ;
      seedHit.dphi = initialDphi;
      seedHit.dtheta = seedTheta - tpcSeedTheta;
      seedHit.rdphi = seedPoint.r * initialDphi;
      seedHit.dz = seedPoint.z - tpcZ;
      seedHit.surface_residual_0 = seedSurfaceMatch.local_residual_0;
      seedHit.surface_residual_1 = seedSurfaceMatch.local_residual_1;
      seedHit.chi2 = square(seedSurfaceMatch.local_residual_0 / m_mvtxLocal0Window) +
                     square(seedSurfaceMatch.local_residual_1 / m_mvtxLocal1Window);
      chain.chi2 = seedHit.chi2;
      chain.hits.push_back(seedHit);
      chain.has_previous_residual = true;
      std::set<TrkrDefs::cluskey> used{seedPoint.key};
      for (const unsigned int layer : m_matchLayers)
      {
        if (layer >= seedLayer) continue;
        ChainHit hit;
        double chi2 = 0.;
        if (findBestMvtxCandidate(trajectory, chain, points, used, layer, hit, chi2, counters))
        {
          chain.hits.push_back(hit);
          chain.chi2 += hit.chi2;
          chain.previous_dphi = hit.dphi;
          chain.previous_dtheta = hit.dtheta;
          chain.has_previous_residual = true;
          used.insert(hit.point.key);
          if (chain.hits.size() == 2U)
          {
            const auto corrected = correctSiliconSeedWithTwoHits(chain.si_reference, chain.hits[0], chain.hits[1]);
            if (corrected.valid) { chain.si_reference = corrected; chain.state = corrected; }
          }
          else if (chain.hits.size() >= 3U)
          {
            const auto refined = correctSiliconSeedWithAllHits(chain.si_reference, chain.hits);
            if (refined.valid) { chain.si_reference = refined; chain.state = refined; }
          }
        }
        else ++chain.n_missing;
      }
      chain.score = chain.chi2 + m_missingLayerPenalty * chain.n_missing;
      if (computeChainDcaMetrics(chain, tpcSeed.state)) chains.push_back(chain);
      if (chains.size() >= m_maxChains) return chains;
    }
  }
  return chains;
}

const TpcSiliconCrossingMatcher::Chain* TpcSiliconCrossingMatcher::selectBestChain(
    const std::vector<Chain>& chains) const
{
  const Chain* best = nullptr;
  for (const auto& chain : chains)
  {
    if (chain.hits.size() < m_minSiliconClusters || !std::isfinite(chain.dca_score) ||
        chain.dca_score >= m_maxChainDcaScore || !std::isfinite(chain.delta_eta0) ||
        chain.delta_eta0 >= m_maxChainDeltaEta) continue;
    if (!best || chain.hits.size() > best->hits.size() ||
        (chain.hits.size() == best->hits.size() && chain.dca_score < best->dca_score)) best = &chain;
  }
  return best;
}

TpcSiliconCrossingMatcher::Chain TpcSiliconCrossingMatcher::attachClosestInttClusters(
    const TpcCrossingTrajectory& trajectory, const Chain& mvtxChain, const std::vector<SpacePoint>& points, Counters& counters) const
{
  Chain output = mvtxChain;
  const TrajectoryState reference = mvtxChain.state;
  for (const unsigned int layer : m_inttMatchLayers)
  {
    const SpacePoint* best = nullptr;
    ChainHit bestHit;
    double bestScore = std::numeric_limits<double>::max();
    for (const auto& point : points)
    {
      if (point.layer != layer || TrkrDefs::getTrkrId(point.key) != TrkrDefs::inttId) continue;
      ++counters.search[layer];
      double predPhi = 0., predZ = 0., predX = 0., predY = 0.;
      if (!predictAtRadius(reference, point.r, predPhi, predZ, predX, predY)) continue;
      const double dz = point.z - predZ;
      if (std::abs(dz) > zSearchWindowCm()) continue;
      ++counters.zPass[layer];
      const double dphi = wrapPhi(point.phi - predPhi);
      const double rdphi = point.r * dphi;
      if (std::abs(rdphi) > m_inttRdphiWindow || std::abs(dz) > m_inttDzWindow) continue;
      SurfaceMatch surfaceMatch;
      if (!matchToSurface(trajectory, point, surfaceMatch) ||
          std::abs(surfaceMatch.local_residual_0) > m_inttLocal0Window) continue;
      ++counters.phiPass[layer];
      if (std::abs(surfaceMatch.local_residual_0) < bestScore)
      {
        bestScore = std::abs(surfaceMatch.local_residual_0);
        best = &point;
        bestHit.point = point;
        bestHit.pred_phi = predPhi;
        bestHit.pred_z = predZ;
        bestHit.dphi = dphi;
        bestHit.rdphi = rdphi;
        bestHit.dz = dz;
        bestHit.ddphi = dphi - dynamicMeanPhi(layer, mvtxChain.previous_dphi, mvtxChain.has_previous_residual);
        const double predTheta = std::atan2(point.r, predZ);
        bestHit.dtheta = pointTheta(point) - predTheta;
        bestHit.surface_residual_0 = surfaceMatch.local_residual_0;
        bestHit.surface_residual_1 = surfaceMatch.local_residual_1;
        // INTT is used only in its precise local measurement direction.
        bestHit.chi2 = square(surfaceMatch.local_residual_0 / m_inttLocal0Window);
      }
    }
    if (best) output.hits.push_back(bestHit);
  }
  return output;
}

int TpcSiliconCrossingMatcher::process_event(PHCompositeNode*)
{
  const auto begin = std::chrono::steady_clock::now();
  m_candidates->Reset();
  const auto siliconPoints = collectSiliconClusters();
  Counters counters;
  unsigned int rejectedPrinted = 0;
  unsigned int surfaceQaPrinted = 0;
  for (unsigned int i = 0; i < m_trajectories->size(); ++i)
  {
    const auto* trajectory = m_trajectories->get(i);
    if (!trajectory || !trajectory->isValid()) continue;
    ++counters.trajectoriesSeen;
    const auto chains = buildChains(*trajectory, siliconPoints, counters);
    const Chain* best = selectBestChain(chains);
    if (!best)
    {
      if (Verbosity() >= 10 && rejectedPrinted++ < m_maxRejectedTrajectoryPrints)
      {
        const auto tpcReference = makeTpcReferenceTrajectory(*trajectory);
        const SpacePoint* closest = nullptr;
        double closestMetric = std::numeric_limits<double>::max();
        double closestPredPhi = 0., closestPredZ = 0., closestDphi = 0., closestDz = 0.;
        for (const auto& point : siliconPoints)
        {
          if (TrkrDefs::getTrkrId(point.key) != TrkrDefs::mvtxId) continue;
          double predPhi = 0., predZ = 0., predX = 0., predY = 0.;
          if (!predictAtRadius(tpcReference, point.r, predPhi, predZ, predX, predY)) continue;
          const double dphi = wrapPhi(point.phi - predPhi);
          const double dz = point.z - predZ;
          const double metric = std::abs(point.r * dphi) / m_looseRdphiWindow +
                                std::abs(dz) / zSearchWindowCm();
          if (metric < closestMetric)
          {
            closestMetric = metric;
            closest = &point;
            closestPredPhi = predPhi;
            closestPredZ = predZ;
            closestDphi = dphi;
            closestDz = dz;
          }
        }
        std::cout << Name() << " rejected parent_track_id=" << trajectory->get_parent_track_id()
                  << " crossing=" << trajectory->get_crossing()
                  << " layer=" << (closest ? static_cast<int>(closest->layer) : -1)
                  << " predicted_phi=" << closestPredPhi
                  << " predicted_z=" << closestPredZ
                  << " best_delta_phi=" << closestDphi
                  << " best_delta_delta_phi=" << closestDphi
                  << " best_delta_z=" << closestDz
                  << " reason=no_valid_mvtx_chain" << std::endl;
      }
      continue;
    }
    Chain output = attachClosestInttClusters(*trajectory, *best, siliconPoints, counters);
    unsigned int nMvtx = 0, nIntt = 0;
    double maxDz = 0., maxDdphi = 0.;
    for (const auto& hit : output.hits)
    {
      if (TrkrDefs::getTrkrId(hit.point.key) == TrkrDefs::mvtxId) ++nMvtx;
      else if (TrkrDefs::getTrkrId(hit.point.key) == TrkrDefs::inttId) ++nIntt;
      maxDz = std::max(maxDz, std::abs(hit.dz));
      maxDdphi = std::max(maxDdphi, std::abs(hit.ddphi));
    }
    if (nMvtx == 1U) ++counters.candidates1Mvtx;
    if (nMvtx == 2U) ++counters.candidates2Mvtx;
    if (nMvtx >= 3U) ++counters.candidates3Mvtx;
    if (nIntt > 0U) ++counters.candidatesWithIntt;
    if (Verbosity() >= 10 && surfaceQaPrinted < m_maxSurfaceQaPrints)
    {
      const auto coarseReference = makeTpcReferenceTrajectory(*trajectory);
      for (const auto& hit : output.hits)
      {
        if (surfaceQaPrinted >= m_maxSurfaceQaPrints) break;
        double oldPhi = 0., oldZ = 0., oldX = 0., oldY = 0.;
        SurfaceMatch surfaceMatch;
        if (!predictAtRadius(coarseReference, hit.point.r, oldPhi, oldZ, oldX, oldY) ||
            !matchToSurface(*trajectory, hit.point, surfaceMatch)) continue;
        const auto& beam = TrkrDefs::getTrkrId(hit.point.key) == TrkrDefs::mvtxId
                               ? m_beamFrame.mvtxBeamLine() : m_beamFrame.inttBeamLine();
        const double oldGlobalX = oldX + beam.x0 + beam.dxdz * oldZ;
        const double oldGlobalY = oldY + beam.y0 + beam.dydz * oldZ;
        std::cout << Name() << " surface_match parent_track_id=" << trajectory->get_parent_track_id()
                  << " crossing=" << trajectory->get_crossing()
                  << " detector=" << (TrkrDefs::getTrkrId(hit.point.key) == TrkrDefs::mvtxId ? "MVTX" : "INTT")
                  << " layer=" << hit.point.layer << " cluster_key=" << hit.point.key
                  << " cluster_global_xyz=" << hit.point.global_x << "," << hit.point.global_y << "," << hit.point.global_z
                  << " cluster_beam_r=" << hit.point.r
                  << " old_radius_prediction_xyz=" << oldGlobalX << "," << oldGlobalY << "," << oldZ
                  << " surface_intersection_xyz=" << surfaceMatch.intersection[0] << ","
                  << surfaceMatch.intersection[1] << "," << surfaceMatch.intersection[2]
                  << " old_rdphi=" << hit.rdphi << " old_dz=" << hit.dz
                  << " surface_local_residual_0=" << surfaceMatch.local_residual_0
                  << " surface_local_residual_1=" << surfaceMatch.local_residual_1
                  << " surface_center_xyz=" << surfaceMatch.surface_center[0] << ","
                  << surfaceMatch.surface_center[1] << "," << surfaceMatch.surface_center[2]
                  << std::endl;
        ++surfaceQaPrinted;
      }
    }
    auto* candidate = new TpcSiliconMatchCandidate;
    candidate->set_parent_track_id(trajectory->get_parent_track_id());
    candidate->set_source_assembled_track_id(trajectory->get_source_assembled_track_id());
    candidate->set_crossing(trajectory->get_crossing());
    candidate->set_n_mvtx(nMvtx);
    candidate->set_n_intt(nIntt);
    candidate->set_score(output.score);
    candidate->set_max_abs_dz(maxDz);
    candidate->set_max_abs_ddphi(maxDdphi);
    for (const auto& hit : output.hits) candidate->add_silicon_cluster_key(hit.point.key);
    m_candidates->add(candidate);
  }
  if (Verbosity() > 0)
  {
    std::cout << Name() << " trajectories_seen=" << counters.trajectoriesSeen
              << " silicon_clusters=" << siliconPoints.size()
              << " z_search_time_bins=" << m_zSearchTimeBins
              << " z_search_window_cm=" << zSearchWindowCm();
    for (const unsigned int layer : m_matchLayers)
      std::cout << " mvtx_l" << layer << "_search=" << counters.search[layer]
                << " mvtx_l" << layer << "_z_pass=" << counters.zPass[layer]
                << " mvtx_l" << layer << "_ddphi_pass=" << counters.phiPass[layer];
    unsigned long long inttSearch = 0, inttZPass = 0, inttPhiPass = 0;
    for (const unsigned int layer : m_inttMatchLayers)
    {
      inttSearch += counters.search[layer];
      inttZPass += counters.zPass[layer];
      inttPhiPass += counters.phiPass[layer];
    }
    std::cout << " intt_search=" << inttSearch
              << " intt_z_pass=" << inttZPass
              << " intt_ddphi_pass=" << inttPhiPass
              << " candidates_1mvtx=" << counters.candidates1Mvtx
              << " candidates_2mvtx=" << counters.candidates2Mvtx
              << " candidates_3mvtx=" << counters.candidates3Mvtx
              << " candidates_with_intt=" << counters.candidatesWithIntt
              << " output_candidates=" << m_candidates->size()
              << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count()
              << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
