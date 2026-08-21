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
  const double dphi_window = 0.15;
  return std::max(dphi_window / std::max(m_phiWindowSigma, 1.0e-9), 1.0e-9);
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
  return std::max(theta_window / std::max(m_thetaWindowSigma, 1.0e-9), 1.0e-9);
}

double Full_PolyTrackMatcher::dynamicDzWindow(double pt) const
{
  const double safe_pt = std::max(pt, 0.25);
  return 1.138 + 0.3919 * std::exp(0.84 / safe_pt);
}

bool Full_PolyTrackMatcher::predictAtRadius(const TrajectoryState& state, double r, double pt,
                                            double& pred_phi, double& pred_z,
                                            double& pred_x, double& pred_y) const
{
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
  if (m_refitWithSiliconHits)
  {
    out.fit_points.push_back(point);
    std::sort(out.fit_points.begin(), out.fit_points.end(), [](const SpacePoint& a, const SpacePoint& b) {
      return a.r > b.r;
    });
    const TrajectoryState refit_state = fitTrajectory(out.fit_points);
    if (refit_state.valid)
    {
      out.state = refit_state;
    }
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
  seed.fit_points.reserve(tpc_clusters.size());
  for (const Tpc_PolyCluster* cluster : tpc_clusters)
  {
    if (!cluster)
    {
      continue;
    }
    SpacePoint point;
    point.key = cluster->get_trkr_cluster_key();
    point.layer = point.key != TrkrDefs::CLUSKEYMAX ? TrkrDefs::getLayer(point.key) : 0;
    point.x = cluster->get_centroid_x();
    point.y = cluster->get_centroid_y();
    point.z = cluster->get_centroid_z();
    point.r = std::hypot(point.x, point.y);
    point.phi = std::atan2(point.y, point.x);
    if (std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) && point.r > 0.0)
    {
      seed.fit_points.push_back(point);
    }
  }

  if (seed.fit_points.size() < 2)
  {
    SpacePoint p0;
    p0.x = track.get_x();
    p0.y = track.get_y();
    p0.z = track.get_z();
    p0.r = std::hypot(p0.x, p0.y);
    p0.phi = std::atan2(p0.y, p0.x);
    SpacePoint p1 = p0;
    const double pt = std::hypot(track.get_px(), track.get_py());
    if (pt > 0.0)
    {
      p1.x += 20.0 * track.get_px() / pt;
      p1.y += 20.0 * track.get_py() / pt;
      p1.z += 20.0 * track.get_pz() / std::max(pt, 1.0e-9);
      p1.r = std::hypot(p1.x, p1.y);
      p1.phi = std::atan2(p1.y, p1.x);
    }
    if (p0.r > 0.0 && p1.r > 0.0)
    {
      seed.fit_points.push_back(p0);
      seed.fit_points.push_back(p1);
    }
  }
  std::sort(seed.fit_points.begin(), seed.fit_points.end(), [](const SpacePoint& a, const SpacePoint& b) {
    return a.r > b.r;
  });

  seed.state = fitTrajectory(seed.fit_points);
  if (!seed.state.valid)
  {
    return {};
  }
  if (Verbosity() > 1)
  {
    std::cout << Name() << "::buildChains - source_track=" << track.get_source_assembled_track_id()
              << " pt=" << seed.pt
              << " fit_points=" << seed.fit_points.size()
              << " circle=" << seed.state.phi_circle_ok
              << " sagitta=" << seed.state.phi_sagitta_ok
              << " refit=" << m_refitWithSiliconHits
              << " circle_center_r=" << std::hypot(seed.state.phi_circle_cx, seed.state.phi_circle_cy)
              << " circle_track_r=" << seed.state.phi_circle_radius
              << " reference_r=" << seed.state.reference_r
              << " reference_phi=" << wrapPhi(seed.state.reference_phi)
              << " reference_z=" << seed.state.reference_z << std::endl;
  }


  std::vector<Chain> chains{seed};
  for (unsigned int layer : m_matchLayers)
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
        const double dphi = wrapPhi(point.phi - pred_phi);
        const double dtheta = pointTheta(point) - pred_theta;
        const double dz = point.z - pred_z;
        const double mean_phi = dynamicMeanPhi(layer, chain.previous_dphi, chain.has_previous_residual);
        const double sigma_phi = dynamicSigmaPhi(chain.pt);
        const double mean_theta = dynamicMeanTheta(layer, chain.previous_dtheta, chain.has_previous_residual);
        const double sigma_theta = dynamicSigmaTheta(chain.pt, pred_theta);
        const double sdphi = (dphi - mean_phi) / sigma_phi;
        const double sdtheta = (dtheta - mean_theta) / sigma_theta;
        if (std::fabs(sdphi) > m_phiWindowSigma || std::fabs(sdtheta) > m_thetaWindowSigma)
        {
          continue;
        }
        if (std::fabs(dz) > dynamicDzWindow(chain.pt))
        {
          continue;
        }

        fillQaDphiCorrelation(track, chain, point, dphi);
        if (chain.has_previous_residual && std::fabs(wrapPhi(dphi - chain.previous_dphi)) > m_deltaDeltaPhiWindow)
        {
          continue;
        }
        const double rdphi = point.r * dphi;
        const double chi2 = square(sdphi) + square(sdtheta);
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
    const std::vector<SpacePoint>& silicon_points) const
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
      const double dphi = wrapPhi(point.phi - pred_phi);
      const double dtheta = pointTheta(point) - pred_theta;
      const double dz = point.z - pred_z;
      const double mean_phi = dynamicMeanPhi(layer, out.previous_dphi, out.has_previous_residual);
      const double sigma_phi = dynamicSigmaPhi(out.pt);
      const double mean_theta = dynamicMeanTheta(layer, out.previous_dtheta, out.has_previous_residual);
      const double sigma_theta = dynamicSigmaTheta(out.pt, pred_theta);
      const double sdphi = (dphi - mean_phi) / sigma_phi;
      const double sdtheta = (dtheta - mean_theta) / sigma_theta;
      if (std::fabs(sdphi) > m_phiWindowSigma || std::fabs(sdtheta) > m_thetaWindowSigma)
      {
        continue;
      }
      if (std::fabs(dz) > dynamicDzWindow(out.pt))
      {
        continue;
      }

      const double rdphi = point.r * dphi;
      const double chi2 = square(sdphi) + square(sdtheta);
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

    Chain output_chain = *selected_chain;
    if (best)
    {
      output_chain = attachClosestInttClusters(*best, silicon_points);
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
