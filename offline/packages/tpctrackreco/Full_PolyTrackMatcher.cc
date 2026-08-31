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
            const double silicon_beam_x = SiliconBeamXIntercept + SiliconBeamXSlope * point.z;
            const double silicon_beam_y = SiliconBeamYIntercept + SiliconBeamYSlope * point.z;
            const double tpc_beam_x = TpcBeamXIntercept + TpcBeamXSlope * point.z;
            const double tpc_beam_y = TpcBeamYIntercept + TpcBeamYSlope * point.z;
            point.x += tpc_beam_x - silicon_beam_x;
            point.y += tpc_beam_y - silicon_beam_y;
            point.r = std::hypot(point.x, point.y);
            point.phi = std::atan2(point.y, point.x);
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
  if (!out.has_z_offset)
  {
    out.z_offset = dz;
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
    const std::vector<SpacePoint>& silicon_points)
{
  Chain seed;
  seed.pt = std::hypot(track.get_px(), track.get_py());
  seed.charge = track.get_charge();
  seed.state = makeTpcSeedTrajectory(track);
  if (!seed.state.valid)
  {
    return {};
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
        if (!predictAtRadius(chain.state, point.r, pred_phi, pred_z, pred_x, pred_y))
        {
          continue;
        }

        const double dz = point.z - pred_z;
        const double z_offset = chain.has_z_offset ? chain.z_offset : dz;
        const double z_residual = dz - z_offset;
        const double shifted_pred_z = pred_z + z_offset;
        const double pred_theta = std::atan2(point.r, shifted_pred_z);
        const double dphi = wrapPhi(point.phi - pred_phi);
        const double dtheta = pointTheta(point) - pred_theta;
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
        if (chain.has_z_offset && std::fabs(z_residual) > dynamicDzWindow(chain.pt))
        {
          continue;
        }

        const double rdphi = point.r * dphi;
        const double chi2 = square(sdphi) + square(sdtheta) + square(z_residual / std::max(m_sigmaDz, 1.0e-9));
        layer_extensions.push_back(extendWithHit(chain, point, pred_phi, pred_z, pred_theta,
                                                 pred_x, pred_y, dphi, dtheta, rdphi, dz, chi2));
        fillQaDphiCorrelation(track, chain, point, dphi);
      }

      std::sort(layer_extensions.begin(), layer_extensions.end(), [](const Chain& a, const Chain& b) {
        if (a.score != b.score) { return a.score < b.score; }
        return a.hits.size() > b.hits.size();
      });
      if (layer_extensions.size() > m_maxBranchesPerLayer)
      {
        layer_extensions.resize(m_maxBranchesPerLayer);
      }
      next_chains.insert(next_chains.end(), layer_extensions.begin(), layer_extensions.end());
      next_chains.push_back(extendMissing(chain, layer));
    }

    std::sort(next_chains.begin(), next_chains.end(), [](const Chain& a, const Chain& b) {
      if (a.score != b.score) { return a.score < b.score; }
      return a.hits.size() > b.hits.size();
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
        chain.score < best->score ||
        (chain.score == best->score && chain.hits.size() > best->hits.size()))
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

  const TrajectoryState tpc_state = out.state;
  const ChainHit* outer_mvtx_hit = nullptr;
  for (const ChainHit& hit : out.hits)
  {
    if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(hit.point.key)) == TrkrDefs::mvtxId &&
        (!outer_mvtx_hit || hit.point.r > outer_mvtx_hit->point.r))
    {
      outer_mvtx_hit = &hit;
    }
  }

  out.has_previous_residual = false;
  out.previous_dphi = 0.0;
  out.previous_dtheta = 0.0;


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
      if (!predictAtRadius(out.state, point.r, pred_phi, pred_z, pred_x, pred_y))
      {
        continue;
      }

      const double dz = point.z - pred_z;
      const double z_offset = out.has_z_offset ? out.z_offset : dz;
      const double z_residual = dz - z_offset;
      const double shifted_pred_z = pred_z + z_offset;
      const double pred_theta = std::atan2(point.r, shifted_pred_z);
      const double dphi = wrapPhi(point.phi - pred_phi);
      const double dtheta = pointTheta(point) - pred_theta;
      const double mean_phi = dynamicMeanPhi(layer, out.previous_dphi, out.has_previous_residual);
      const double sigma_phi = dynamicSigmaPhi(out.pt);
      const double mean_theta = dynamicMeanTheta(layer, out.previous_dtheta, out.has_previous_residual);
      const double sigma_theta = dynamicSigmaTheta(out.pt, pred_theta);
      if (outer_mvtx_hit)
      {
        double tpc_pred_phi = 0.0;
        double tpc_pred_z = 0.0;
        double tpc_pred_x = 0.0;
        double tpc_pred_y = 0.0;
        if (predictAtRadius(tpc_state, point.r, tpc_pred_phi, tpc_pred_z, tpc_pred_x, tpc_pred_y))
        {
          const double mvtx_phi = outer_mvtx_hit->point.phi;
          const double mvtx_theta = pointTheta(outer_mvtx_hit->point);
          const double tpc_pred_theta = std::atan2(point.r, tpc_pred_z + z_offset);
          const double dphi_to_tpc = wrapPhi(tpc_pred_phi - mvtx_phi);
          const double dphi_to_point = wrapPhi(point.phi - mvtx_phi);
          const double dtheta_to_tpc = tpc_pred_theta - mvtx_theta;
          const double dtheta_to_point = pointTheta(point) - mvtx_theta;
          const double phi_tolerance = m_phiWindowSigma * sigma_phi;
          const double theta_tolerance = m_thetaWindowSigma * sigma_theta;
          if ((std::fabs(dphi_to_tpc) > phi_tolerance && dphi_to_tpc * dphi_to_point < 0.0) ||
              (std::fabs(dtheta_to_tpc) > theta_tolerance && dtheta_to_tpc * dtheta_to_point < 0.0) ||
              std::fabs(dphi_to_point) > std::fabs(dphi_to_tpc) + phi_tolerance ||
              std::fabs(dtheta_to_point) > std::fabs(dtheta_to_tpc) + theta_tolerance)
          {
            continue;
          }
        }
      }
      const double sdphi = (dphi - mean_phi) / sigma_phi;
      const double sdtheta = (dtheta - mean_theta) / sigma_theta;
      if (std::fabs(sdphi) > m_phiWindowSigma || std::fabs(sdtheta) > m_thetaWindowSigma)
      {
        continue;
      }
      if (out.has_z_offset && std::fabs(z_residual) > dynamicDzWindow(out.pt))
      {
        continue;
      }

      const double rdphi = point.r * dphi;
      const double chi2 = square(rdphi / std::max(m_sigmaRdphi, 1.0e-9)) + square(z_residual / std::max(m_sigmaDz, 1.0e-9));
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