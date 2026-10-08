#include "SiHitSeedQA.h"
#include "SiHitSeedData.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <intt/CylinderGeomIntt.h>
#include <mvtx/CylinderGeom_Mvtx.h>
#include <trackbase/ActsGeometry.h>

#include <Acts/Definitions/Units.hpp>

#include <TVector3.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2D.h>
#include <TList.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

double wrapDelta(double x)
{
  x = std::fmod(x + kPi, kTwoPi);
  if (x < 0)
  {
    x += kTwoPi;
  }
  return x - kPi;
}

double median(std::vector<double> v)
{
  if (v.empty())
  {
    return kNaN;
  }
  const size_t m = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + m, v.end());
  double med = v[m];
  if (v.size() % 2 == 0)
  {
    med = 0.5 * (med + *std::max_element(v.begin(), v.begin() + m));
  }
  return med;
}
}  // namespace

SiHitSeedQA::SiHitSeedQA(const std::string& out, const std::string& name)
  : SubsysReco(name)
  , m_out(out)
{
}

int SiHitSeedQA::Init(PHCompositeNode*)
{
  m_file = new TFile(m_out.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cerr << Name() << ": cannot create " << m_out << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_profDir = m_file->mkdir("vertex_profiles");
  m_file->cd();

#define BRANCH(tree, name, variable) tree->Branch(name, &variable)

  m_evt = new TTree("events", "event summary");
  BRANCH(m_evt, "event", event);
  BRANCH(m_evt, "nhits", nhits);
  BRANCH(m_evt, "nclusters", nclusters);
  BRANCH(m_evt, "nchains", nchains);
  BRANCH(m_evt, "vertex_z", vertex_z);                  // SiHitSeedReco [cm]: tracklet z0 peak (window independent)
  BRANCH(m_evt, "vertex_z_linefit", vertex_z_linefit);  // SiHitSeedReco [cm]: dz=0 crossing of accepted seed links
  BRANCH(m_evt, "vertex_tracklet_npairs", vtx_tr_npairs);
  BRANCH(m_evt, "vertex_tracklet_npeak", vtx_tr_npeak);
  BRANCH(m_evt, "dz_vs_z_intercept_bins", fit_a);
  BRANCH(m_evt, "dz_vs_z_slope_bins_per_cm", fit_b);
  BRANCH(m_evt, "n_vertex_links", nvertex);
  BRANCH(m_evt, "n_vertex_chains", n_vertex_chains);
  BRANCH(m_evt, "n_vertex_tracklets", n_vertex_tracklets);
  BRANCH(m_evt, "vertex_phi_sectors", vertex_phi_sectors);
  BRANCH(m_evt, "vertex_phi_coverage_ok", vertex_phi_coverage_ok);
  BRANCH(m_evt, "vertex_half_agreement_ok", vertex_half_agreement_ok);
  BRANCH(m_evt, "vertex_selected", vertex_selected);
  BRANCH(m_evt, "vertex_half_dx", vertex_half_dx);
  BRANCH(m_evt, "vertex_half_dy", vertex_half_dy);
  BRANCH(m_evt, "final_vtx_x", final_vtx_x);
  BRANCH(m_evt, "final_vtx_y", final_vtx_y);
  BRANCH(m_evt, "final_vtx_z", final_vtx_z);
  BRANCH(m_evt, "final_vtx_x_err", final_vtx_x_err);
  BRANCH(m_evt, "final_vtx_y_err", final_vtx_y_err);
  BRANCH(m_evt, "final_vtx_z_err", final_vtx_z_err);
  BRANCH(m_evt, "final_vtx_ok", final_vtx_ok);
  // QA vertex fits [mm, detector-centred frame]: all links, clamshell half A, half B
  branchVertex(m_evt, "vtx", m_vtxAll);
  branchVertex(m_evt, "vtxA", m_vtxA);
  branchVertex(m_evt, "vtxB", m_vtxB);
  // Stage-7 independent estimator: one 3-MVTX straight tracklet per accepted chain.
  branchVertex(m_evt, "trk_vtx", m_trkVtxAll);
  branchVertex(m_evt, "trk_vtxA", m_trkVtxA);
  branchVertex(m_evt, "trk_vtxB", m_trkVtxB);
  // standard vertices: untouched [cm, global frame] and in the detector frame [mm]
  BRANCH(m_evt, "svx_x", svx_x);
  BRANCH(m_evt, "svx_y", svx_y);
  BRANCH(m_evt, "svx_z", svx_z);
  BRANCH(m_evt, "svx_ex", svx_ex);
  BRANCH(m_evt, "svx_ey", svx_ey);
  BRANCH(m_evt, "svx_ez", svx_ez);
  BRANCH(m_evt, "svx_chi2", svx_chi2);
  BRANCH(m_evt, "svx_ndof", svx_ndof);
  BRANCH(m_evt, "svx_ntracks", svx_ntracks);
  BRANCH(m_evt, "svx_det_x", svx_det_x);
  BRANCH(m_evt, "svx_det_y", svx_det_y);
  BRANCH(m_evt, "svx_det_z", svx_det_z);
  // detector frame in the global frame (accumulated fit) and per-event fit quality
  BRANCH(m_evt, "det_T_x", det_T_x);
  BRANCH(m_evt, "det_T_y", det_T_y);
  BRANCH(m_evt, "det_T_z", det_T_z);
  BRANCH(m_evt, "det_rot_x", det_rot_x);
  BRANCH(m_evt, "det_rot_y", det_rot_y);
  BRANCH(m_evt, "det_rot_z", det_rot_z);
  BRANCH(m_evt, "det_fit_rms", det_fit_rms);
  BRANCH(m_evt, "det_fit_npoints", det_fit_npoints);
  BRANCH(m_evt, "evt_fit_rms", evt_fit_rms);
  BRANCH(m_evt, "evt_fit_npoints", evt_fit_npoints);

  m_hits = new TTree("hits", "all silicon hits in notebook coordinates");
  BRANCH(m_hits, "event", event);
  BRANCH(m_hits, "hitsetkey", hitsetkey);
  BRANCH(m_hits, "hitkey", hitkey);
  BRANCH(m_hits, "layer", layer);
  BRANCH(m_hits, "row", row);
  BRANCH(m_hits, "col", col);
  BRANCH(m_hits, "stave", stave);
  BRANCH(m_hits, "chip", chip);
  BRANCH(m_hits, "ladderphi", ladderphi);
  BRANCH(m_hits, "ladderz", ladderz);
  BRANCH(m_hits, "adc", adc);
  BRANCH(m_hits, "lx", lx);
  BRANCH(m_hits, "ly", ly);
  BRANCH(m_hits, "z", z);
  BRANCH(m_hits, "phi", phi);

  m_clus = new TTree("clusters", "all silicon clusters in notebook coordinates");
  BRANCH(m_clus, "event", event);
  BRANCH(m_clus, "cluskey", cluskey);
  BRANCH(m_clus, "layer", layer);
  BRANCH(m_clus, "stave", stave);
  BRANCH(m_clus, "chip", chip);
  BRANCH(m_clus, "ladderphi", ladderphi);
  BRANCH(m_clus, "ladderz", ladderz);
  BRANCH(m_clus, "adc", adc);
  BRANCH(m_clus, "lx", lx);
  BRANCH(m_clus, "ly", ly);
  BRANCH(m_clus, "z", z);
  BRANCH(m_clus, "phi", phi);

  m_steps = new TTree("steps", "accepted chain links");
  BRANCH(m_steps, "event", event);
  BRANCH(m_steps, "chain_id", chain_id);
  BRANCH(m_steps, "kind", kind);
  BRANCH(m_steps, "from_layer", from_layer);
  BRANCH(m_steps, "to_layer", to_layer);
  BRANCH(m_steps, "ref_z", ref_z);
  BRANCH(m_steps, "ref_phi", ref_phi);
  BRANCH(m_steps, "pred_z", pred_z);
  BRANCH(m_steps, "pred_phi", pred_phi);
  BRANCH(m_steps, "z", z);
  BRANCH(m_steps, "phi", phi);
  BRANCH(m_steps, "delta_z_bins", delta_z);
  BRANCH(m_steps, "delta_phi_bins", delta_phi);
  BRANCH(m_steps, "res_z_bins", res_z);
  BRANCH(m_steps, "res_phi_bins", res_phi);
  BRANCH(m_steps, "score", score);

  m_chains = new TTree("chains", "reconstructed silicon hit chains");
  BRANCH(m_chains, "event", event);
  BRANCH(m_chains, "chain_id", chain_id);
  BRANCH(m_chains, "layers", vlayer);
  BRANCH(m_chains, "z", vz);
  BRANCH(m_chains, "phi", vphi);
  BRANCH(m_chains, "n_mvtx", chain_nmvtx);
  BRANCH(m_chains, "n_intt", chain_nintt);
  BRANCH(m_chains, "eta", chain_eta);
  BRANCH(m_chains, "chain_score", chain_score);
  BRANCH(m_chains, "max_prop_res_z_bins", chain_max_res_z);
  BRANCH(m_chains, "max_prop_res_phi_bins", chain_max_res_phi);
  BRANCH(m_chains, "passes_vertex_cuts", chain_pass_vertex_cuts);

#undef BRANCH

  m_dz = new TH2D(
      "h_seed_delta_z_vs_z",
      "seed #Delta z vs outer z;outer nominal V_{z} [cm];#Delta z [bins]",
      120, -20, 20, 160, -20, 20);

  m_dzcm = new TH2D(
      "h_seed_delta_z_cm_vs_z",
      "seed #Delta z vs outer z (cm, independent of the per-event bin width);outer nominal V_{z} [cm];#Delta z = z_{in} - z_{out} [cm]",
      120, -20, 20, 160, -8, 8);

  m_dphi = new TH2D(
      "h_seed_delta_phi_vs_phi",
      "seed #Delta U_{#phi} vs outer U_{#phi};outer U_{#phi} [rad];#Delta U_{#phi} [bins]",
      180, 0, 2.0 * M_PI, 120, -10, 10);

  return Fun4AllReturnCodes::EVENT_OK;
}

void SiHitSeedQA::branchVertex(TTree* t, const std::string& p, Vertex& v)
{
  auto B = [&](const std::string& n, auto& var) { t->Branch((p + "_" + n).c_str(), &var); };
  B("nlinks", v.nlinks);
  B("ok", v.ok);
  B("z", v.z);
  B("z_err", v.z_err);
  B("z_median", v.z_median);
  B("dz_intercept", v.dz_intercept);
  B("dz_slope", v.dz_slope);
  B("dz_rms", v.dz_rms);
  B("z_nlinks", v.z_nlinks);
  B("z_nbins", v.z_nbins);
  B("x", v.x);
  B("y", v.y);
  B("x_err", v.x_err);
  B("y_err", v.y_err);
  B("xy_corr", v.xy_corr);
  B("xy_const", v.xy_const);
  B("xy_rms", v.xy_rms);
  B("xy_nlinks", v.xy_nlinks);
}


double SiHitSeedQA::chainEta(const SiHitChain& c) const
{
  int iInner = -1;
  int iOuter = -1;
  double rInner = std::numeric_limits<double>::max();
  double rOuter = -1.0;

  const size_t n = std::min(c.layers.size(), c.z.size());
  for (size_t i = 0; i < n; ++i)
  {
    const int l = c.layers[i];
    if (l < 0 || l > 2)
    {
      continue;
    }

    const double r = m_radius[l];
    if (r < rInner)
    {
      rInner = r;
      iInner = static_cast<int>(i);
    }
    if (r > rOuter)
    {
      rOuter = r;
      iOuter = static_cast<int>(i);
    }
  }

  if (iInner < 0 || iOuter < 0 || iInner == iOuter)
  {
    return kNaN;
  }

  const double dr = rOuter - rInner;
  if (std::abs(dr) < 1e-9)
  {
    return kNaN;
  }

  return std::asinh((c.z[iOuter] - c.z[iInner]) / dr);
}

double SiHitSeedQA::chainMaxResidualZ(const SiHitChain& c) const
{
  double out = 0.0;
  for (const auto& s : c.steps)
  {
    if (s.kind == SiChainStep::Seed)
    {
      continue;
    }
    out = std::max(out, std::abs(s.res_z_bins));
  }
  return out;
}

double SiHitSeedQA::chainMaxResidualPhi(const SiHitChain& c) const
{
  double out = 0.0;
  for (const auto& s : c.steps)
  {
    if (s.kind == SiChainStep::Seed)
    {
      continue;
    }
    out = std::max(out, std::abs(s.res_phi_bins));
  }
  return out;
}

bool SiHitSeedQA::chainPassesVertexCuts(const SiHitChain& c) const
{
  if (c.n_mvtx < m_vtxMinMvtxLayers)
  {
    return false;
  }

  if (c.n_intt < m_vtxMinInttLayers)
  {
    return false;
  }

  if (m_vtxMaxAbsEta > 0)
  {
    const double eta = chainEta(c);
    if (!std::isfinite(eta) || std::abs(eta) > m_vtxMaxAbsEta)
    {
      return false;
    }
  }

  if (m_vtxMaxChainScore > 0 && c.score > m_vtxMaxChainScore)
  {
    return false;
  }

  if (m_vtxMaxResidualZBins > 0 &&
      chainMaxResidualZ(c) > m_vtxMaxResidualZBins)
  {
    return false;
  }

  if (m_vtxMaxResidualPhiBins > 0 &&
      chainMaxResidualPhi(c) > m_vtxMaxResidualPhiBins)
  {
    return false;
  }

  return true;
}

bool SiHitSeedQA::makeMvtxTracklet(const SiHitChain& c, Tracklet& t) const
{
  std::vector<Eigen::Vector3d> pts;

  const size_t n = std::min({c.layers.size(), c.z.size(), c.phi.size()});
  for (size_t i = 0; i < n; ++i)
  {
    const int l = c.layers[i];
    if (l < 0 || l > 2)
    {
      continue;
    }

    const double r = m_radius[l];
    const double ph = c.phi[i];
    pts.emplace_back(r * std::cos(ph), r * std::sin(ph), c.z[i]);
  }

  if (pts.size() != 3)
  {
    return false;
  }

  Eigen::Vector3d center = Eigen::Vector3d::Zero();
  for (const auto& p : pts)
  {
    center += p;
  }
  center /= static_cast<double>(pts.size());

  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto& p : pts)
  {
    const Eigen::Vector3d d = p - center;
    cov += d * d.transpose();
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
  if (es.info() != Eigen::Success)
  {
    return false;
  }

  Eigen::Vector3d dir = es.eigenvectors().col(2).normalized();

  // Give the direction a deterministic outward sign.
  if (center.x() * dir.x() + center.y() * dir.y() < 0)
  {
    dir = -dir;
  }

  double s2 = 0.0;
  for (const auto& p : pts)
  {
    const Eigen::Vector3d q = p - center;
    const Eigen::Vector3d perp = q - q.dot(dir) * dir;
    s2 += perp.squaredNorm();
  }

  const double dxy = std::hypot(dir.x(), dir.y());
  if (dxy < 1e-12)
  {
    return false;
  }

  t.point = {center.x(), center.y(), center.z()};
  t.dir = {dir.x(), dir.y(), dir.z()};
  t.eta = std::asinh(dir.z() / dxy);
  t.phi = std::atan2(center.y(), center.x());
  if (t.phi < 0)
  {
    t.phi += kTwoPi;
  }
  t.rms = std::sqrt(s2 / pts.size());
  return true;
}

void SiHitSeedQA::fitTrackletVertex(const std::vector<Tracklet>& tracklets, Vertex& v) const
{
  resetVertex(v);
  v.nlinks = static_cast<int>(tracklets.size());

  if (tracklets.size() < 3)
  {
    return;
  }

  std::vector<char> keep(tracklets.size(), 1);
  Eigen::Vector3d vertex = Eigen::Vector3d::Zero();
  Eigen::Matrix3d vertexCov = Eigen::Matrix3d::Zero();
  double rms = 0.0;
  int nkeep = 0;

  for (int iter = 0; iter <= m_vtxIter; ++iter)
  {
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
    Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
    nkeep = 0;

    for (size_t i = 0; i < tracklets.size(); ++i)
    {
      if (!keep[i])
      {
        continue;
      }

      const Eigen::Vector3d p(tracklets[i].point[0],
                              tracklets[i].point[1],
                              tracklets[i].point[2]);
      const Eigen::Vector3d d(tracklets[i].dir[0],
                              tracklets[i].dir[1],
                              tracklets[i].dir[2]);

      const Eigen::Matrix3d P =
          Eigen::Matrix3d::Identity() - d * d.transpose();

      A += P;
      rhs += P * p;
      ++nkeep;
    }

    if (nkeep < 3 || std::abs(A.determinant()) < 1e-12)
    {
      return;
    }

    vertex = A.ldlt().solve(rhs);

    std::vector<double> residuals;
    residuals.reserve(nkeep);
    double s2 = 0.0;

    for (size_t i = 0; i < tracklets.size(); ++i)
    {
      if (!keep[i])
      {
        continue;
      }

      const Eigen::Vector3d p(tracklets[i].point[0],
                              tracklets[i].point[1],
                              tracklets[i].point[2]);
      const Eigen::Vector3d d(tracklets[i].dir[0],
                              tracklets[i].dir[1],
                              tracklets[i].dir[2]);
      const Eigen::Vector3d q = vertex - p;
      const double r = (q - q.dot(d) * d).norm();
      residuals.push_back(r);
      s2 += r * r;
    }

    rms = residuals.empty() ? 0.0 : std::sqrt(s2 / residuals.size());

    // Use a robust MAD scale for rejection so the first iteration is not
    // inflated by a few fake tracklets.
    const double med = median(residuals);
    std::vector<double> ad;
    ad.reserve(residuals.size());
    for (double r : residuals)
    {
      ad.push_back(std::abs(r - med));
    }
    double robustSigma = 1.4826 * median(ad);
    if (!std::isfinite(robustSigma) || robustSigma <= 1e-9)
    {
      robustSigma = rms;
    }

    vertexCov = rms * rms * A.inverse();

    if (iter == m_vtxIter || robustSigma <= 0)
    {
      break;
    }

    for (size_t i = 0; i < tracklets.size(); ++i)
    {
      const Eigen::Vector3d p(tracklets[i].point[0],
                              tracklets[i].point[1],
                              tracklets[i].point[2]);
      const Eigen::Vector3d d(tracklets[i].dir[0],
                              tracklets[i].dir[1],
                              tracklets[i].dir[2]);
      const Eigen::Vector3d q = vertex - p;
      const double r = (q - q.dot(d) * d).norm();

      keep[i] = std::abs(r - med) <= m_vtxClip * robustSigma;
    }
  }

  const double cr = std::cos(m_rotation);
  const double sr = std::sin(m_rotation);

  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  R(0, 0) = cr;  R(0, 1) = -sr;
  R(1, 0) = sr;  R(1, 1) = cr;

  const Eigen::Vector3d g = R * vertex;
  const Eigen::Matrix3d cg = R * vertexCov * R.transpose();

  v.x = 10.0 * g.x() + m_center[0];
  v.y = 10.0 * g.y() + m_center[1];
  v.z = 10.0 * g.z() + m_center[2];

  v.x_err = 10.0 * std::sqrt(std::max(0.0, cg(0, 0)));
  v.y_err = 10.0 * std::sqrt(std::max(0.0, cg(1, 1)));
  v.z_err = 10.0 * std::sqrt(std::max(0.0, cg(2, 2)));

  v.z_median = v.z;
  v.xy_corr = (cg(0, 0) > 0 && cg(1, 1) > 0)
                  ? cg(0, 1) / std::sqrt(cg(0, 0) * cg(1, 1))
                  : kNaN;
  v.xy_rms = 10.0 * rms;
  v.dz_rms = 10.0 * rms;
  v.xy_nlinks = nkeep;
  v.z_nlinks = nkeep;
  v.z_nbins = 0;
  v.ok = 3;
}

int SiHitSeedQA::countPhiSectors(const std::vector<double>& phis) const
{
  if (m_vtxPhiSectors <= 0)
  {
    return 0;
  }

  std::vector<char> occupied(m_vtxPhiSectors, 0);
  for (double p : phis)
  {
    p = std::fmod(p, kTwoPi);
    if (p < 0)
    {
      p += kTwoPi;
    }

    int ib = static_cast<int>(std::floor(p / kTwoPi * m_vtxPhiSectors));
    ib = std::clamp(ib, 0, m_vtxPhiSectors - 1);
    occupied[ib] = 1;
  }

  int n = 0;
  for (char x : occupied)
  {
    n += x ? 1 : 0;
  }
  return n;
}

void SiHitSeedQA::resetVertex(Vertex& v) const
{
  v.z = v.z_err = v.z_median = v.dz_intercept = v.dz_slope = v.dz_rms = kNaN;
  v.x = v.y = v.x_err = v.y_err = v.xy_corr = v.xy_const = v.xy_rms = kNaN;
  v.a_cm = v.b_cm = v.c_cm = v.A_cm = v.B_cm = kNaN;
  v.z_nlinks = v.z_nbins = v.xy_nlinks = v.ok = v.nlinks = 0;
}

void SiHitSeedQA::fitVertex(const std::vector<Link>& links, Vertex& v) const
{
  resetVertex(v);
  v.nlinks = static_cast<int>(links.size());
  std::vector<double> zAtR0;
  zAtR0.reserve(links.size());
  for (const auto& l : links)
  {
    // straight line through the two hits, extrapolated to R = 0
    const double rf = m_radius[l.from_layer], rt = m_radius[l.to_layer];
    zAtR0.push_back(l.z_from - l.dz * rf / (rt - rf));
  }
  const double zMedCm = median(zAtR0);
  v.z_median = 10.0 * zMedCm + m_center[2];
  fitVertexZ(links, zMedCm, v);
  fitVertexXY(links, v);
}

// ---------------------------------------------------------------------------
// Vertex z: binned profile of dz vs z_outer, weighted straight-line fit of the bin means
// (weight = entries), iterated with rejection of links far from the line.
// ---------------------------------------------------------------------------
void SiHitSeedQA::fitVertexZ(const std::vector<Link>& links, double zSeed, Vertex& v) const
{
  // Start from links whose own straight-line extrapolation to R = 0 lies within
  // m_vtxSeedWindow of the robust (median) estimate: random combinations far from the
  // line would otherwise pull the bin means towards dz = 0 and bias the crossing.
  std::vector<char> keep(links.size(), 1);
  if (std::isfinite(zSeed) && m_vtxSeedWindow > 0)
  {
    for (size_t i = 0; i < links.size(); ++i)
    {
      const double rf = m_radius[links[i].from_layer], rt = m_radius[links[i].to_layer];
      keep[i] = std::abs(links[i].z_from - links[i].dz * rf / (rt - rf) - zSeed) <= m_vtxSeedWindow;
    }
  }
  const double bw = (m_vtxZhi - m_vtxZlo) / m_vtxNz;
  double a = 0, b = 0, rms = 0;
  Eigen::Matrix2d cov = Eigen::Matrix2d::Zero();
  bool ok = false;

  for (int it = 0; it <= m_vtxIter; ++it)
  {
    std::vector<double> n(m_vtxNz, 0.0), sy(m_vtxNz, 0.0);
    for (size_t i = 0; i < links.size(); ++i)
    {
      const int ib = (int) std::floor((links[i].z_from - m_vtxZlo) / bw);
      if (!keep[i] || ib < 0 || ib >= m_vtxNz)
      {
        continue;
      }
      n[ib] += 1;
      sy[ib] += links[i].dz;
    }
    // weighted LSQ on bin means: dz_mean = a + b * z_center
    Eigen::Matrix2d A = Eigen::Matrix2d::Zero();
    Eigen::Vector2d rhs = Eigen::Vector2d::Zero();
    int nb = 0;
    for (int ib = 0; ib < m_vtxNz; ++ib)
    {
      if (n[ib] < m_vtxMinPerBin)
      {
        continue;
      }
      const double x = m_vtxZlo + (ib + 0.5) * bw, y = sy[ib] / n[ib], w = n[ib];
      A(0, 0) += w;
      A(0, 1) += w * x;
      A(1, 1) += w * x * x;
      rhs(0) += w * y;
      rhs(1) += w * x * y;
      ++nb;
    }
    A(1, 0) = A(0, 1);
    if (nb < 2 || std::abs(A.determinant()) < 1e-12)
    {
      ok = false;
      break;
    }
    const Eigen::Vector2d p = A.ldlt().solve(rhs);
    a = p(0);
    b = p(1);
    // rms of the individual links around the line, then reject outliers
    double s2 = 0;
    int nk = 0;
    for (size_t i = 0; i < links.size(); ++i)
    {
      if (keep[i])
      {
        const double r = links[i].dz - (a + b * links[i].z_from);
        s2 += r * r;
        ++nk;
      }
    }
    rms = nk > 2 ? std::sqrt(s2 / (nk - 2)) : 0.0;
    cov = rms * rms * A.inverse();  // var(bin mean) = rms^2 / n
    v.z_nlinks = nk;
    v.z_nbins = nb;
    ok = true;
    if (it == m_vtxIter || rms <= 0)
    {
      break;
    }
    for (size_t i = 0; i < links.size(); ++i)
    {
      keep[i] = std::abs(links[i].dz - (a + b * links[i].z_from)) <= m_vtxClip * rms;
    }
  }

  if (!ok || std::abs(b) < 1e-6)
  {
    return;
  }
  const double zv = -a / b;  // cm
  const double da = -1.0 / b, db = a / (b * b);
  const double zerr = std::sqrt(std::max(0.0, da * da * cov(0, 0) + db * db * cov(1, 1) + 2 * da * db * cov(0, 1)));
  v.a_cm = a;
  v.b_cm = b;
  v.z = 10.0 * zv + m_center[2];
  v.z_err = 10.0 * zerr;
  v.dz_intercept = 10.0 * a;  // dz[mm] = intercept[mm] + slope * z_outer[mm] (nominal z, no offset)
  v.dz_slope = b;
  v.dz_rms = 10.0 * rms;
  v.ok |= 1;
}

// ---------------------------------------------------------------------------
// Vertex xy: unbinned fit of dphi/k = c + A cos(phi_out) + B sin(phi_out); x0 = -B, y0 = A
// in the Uphi frame, then rotated to the output frame.
// ---------------------------------------------------------------------------
void SiHitSeedQA::fitVertexXY(const std::vector<Link>& links, Vertex& v) const
{
  std::vector<double> yv(links.size());
  for (size_t i = 0; i < links.size(); ++i)
  {
    const double k = 1.0 / m_radius[links[i].to_layer] - 1.0 / m_radius[links[i].from_layer];
    yv[i] = links[i].dphi / k;
  }
  const int np = m_vtxXYConst ? 3 : 2;
  auto basis = [&](double ph)
  {
    Eigen::VectorXd f(np);
    if (m_vtxXYConst)
    {
      f << 1.0, std::cos(ph), std::sin(ph);
    }
    else
    {
      f << std::cos(ph), std::sin(ph);
    }
    return f;
  };
  std::vector<char> keep(links.size(), 1);
  Eigen::VectorXd p = Eigen::VectorXd::Zero(np);
  Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(np, np);
  double rms = 0;
  bool ok = false;

  for (int it = 0; it <= m_vtxIter; ++it)
  {
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(np, np);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(np);
    int nk = 0;
    for (size_t i = 0; i < links.size(); ++i)
    {
      if (!keep[i])
      {
        continue;
      }
      const Eigen::VectorXd f = basis(links[i].phi_from);
      A += f * f.transpose();
      rhs += f * yv[i];
      ++nk;
    }
    if (nk < np + 3 || std::abs(A.determinant()) < 1e-9)
    {
      ok = false;
      break;
    }
    p = A.ldlt().solve(rhs);
    double s2 = 0;
    for (size_t i = 0; i < links.size(); ++i)
    {
      if (keep[i])
      {
        const double r = yv[i] - basis(links[i].phi_from).dot(p);
        s2 += r * r;
      }
    }
    rms = std::sqrt(s2 / (nk - np));
    cov = rms * rms * A.inverse();
    v.xy_nlinks = nk;
    ok = true;
    if (it == m_vtxIter || rms <= 0)
    {
      break;
    }
    for (size_t i = 0; i < links.size(); ++i)
    {
      keep[i] = std::abs(yv[i] - basis(links[i].phi_from).dot(p)) <= m_vtxClip * rms;
    }
  }
  if (!ok)
  {
    return;
  }
  const int iA = m_vtxXYConst ? 1 : 0, iB = iA + 1;
  v.c_cm = m_vtxXYConst ? p(0) : 0.0;
  v.A_cm = p(iA);
  v.B_cm = p(iB);
  // Uphi frame [cm]: xu = -B, yu = A ; covariance of (xu, yu)
  const double xu = -p(iB), yu = p(iA);
  Eigen::Matrix2d cu;
  cu << cov(iB, iB), -cov(iB, iA), -cov(iA, iB), cov(iA, iA);
  // rotate to the output frame and convert to mm
  const double cr = std::cos(m_rotation), sr = std::sin(m_rotation);
  Eigen::Matrix2d R;
  R << cr, -sr, sr, cr;
  const Eigen::Vector2d g = R * Eigen::Vector2d(xu, yu);
  const Eigen::Matrix2d cg = R * cu * R.transpose();
  v.x = 10.0 * g(0) + m_center[0];
  v.y = 10.0 * g(1) + m_center[1];
  v.x_err = 10.0 * std::sqrt(cg(0, 0));
  v.y_err = 10.0 * std::sqrt(cg(1, 1));
  v.xy_corr = cg(0, 1) / std::sqrt(cg(0, 0) * cg(1, 1));
  v.xy_const = 10.0 * v.c_cm;
  v.xy_rms = 10.0 * rms;
  v.ok |= 2;
}

void SiHitSeedQA::saveEventProfiles(unsigned long long evt, const char* label, const std::vector<Link>& links, const Vertex& v)
{
  m_profDir->cd();
  const TString tag = Form("evt%llu_%s", evt, label);
  auto* pz = new TProfile("p_" + tag + "_dz_vs_z",
                          Form("event %llu %s: #Delta z vs outer z, v_{z} = %.2f mm;outer nominal V_{z} [cm];#Delta z = z_{in} - z_{out} [cm]",
                               evt, label, v.z),
                          m_vtxNz, m_vtxZlo, m_vtxZhi);
  auto* pxy = new TProfile("p_" + tag + "_dphik_vs_phi",
                           Form("event %llu %s: #Delta#phi/k vs outer U_{#phi}, v_{x} = %.2f, v_{y} = %.2f mm;outer U_{#phi} [rad];#Delta#phi / (1/R_{in} - 1/R_{out}) [cm]",
                                evt, label, v.x, v.y),
                           36, 0.0, kTwoPi);
  // the same links as 2D histograms (every link, before any outlier rejection)
  auto* hz = new TH2D("h_" + tag + "_dz_vs_z",
                      Form("event %llu %s: #Delta z vs outer z, v_{z} = %.2f mm;outer nominal V_{z} [cm];#Delta z = z_{in} - z_{out} [cm]",
                           evt, label, v.z),
                      m_vtxNz, m_vtxZlo, m_vtxZhi, m_h2dzN, m_h2dzLo, m_h2dzHi);
  auto* hxy = new TH2D("h_" + tag + "_dphik_vs_phi",
                       Form("event %llu %s: #Delta#phi/k vs outer U_{#phi}, v_{x} = %.2f, v_{y} = %.2f mm;outer U_{#phi} [rad];#Delta#phi / (1/R_{in} - 1/R_{out}) [cm]",
                            evt, label, v.x, v.y),
                       36, 0.0, kTwoPi, m_h2xyN, m_h2xyLo, m_h2xyHi);
  hz->SetDirectory(nullptr);
  hxy->SetDirectory(nullptr);
  hz->SetOption("COLZ");
  hxy->SetOption("COLZ");
  hz->SetStats(false);
  hxy->SetStats(false);
  for (const auto& l : links)
  {
    const double dphik = l.dphi / (1.0 / m_radius[l.to_layer] - 1.0 / m_radius[l.from_layer]);
    pz->Fill(l.z_from, l.dz);
    pxy->Fill(l.phi_from, dphik);
    hz->Fill(l.z_from, l.dz);
    hxy->Fill(l.phi_from, dphik);
  }
  // fitted curves on the profiles and on the 2D histograms (black, drawn over the COLZ map)
  if (v.ok & 1)
  {
    auto* f = new TF1("f_" + tag + "_dz", "[0]+[1]*x", m_vtxZlo, m_vtxZhi);
    f->SetParameters(v.a_cm, v.b_cm);
    pz->GetListOfFunctions()->Add(f);
    auto* f2 = static_cast<TF1*>(f->Clone("f_" + tag + "_dz_2d"));
    f2->SetLineColor(kBlack);
    f2->SetLineWidth(2);
    hz->GetListOfFunctions()->Add(f2);
  }
  if (v.ok & 2)
  {
    auto* f = new TF1("f_" + tag + "_xy", "[0]+[1]*cos(x)+[2]*sin(x)", 0.0, kTwoPi);
    f->SetParameters(v.c_cm, v.A_cm, v.B_cm);
    pxy->GetListOfFunctions()->Add(f);
    auto* f2 = static_cast<TF1*>(f->Clone("f_" + tag + "_xy_2d"));
    f2->SetLineColor(kBlack);
    f2->SetLineWidth(2);
    hxy->GetListOfFunctions()->Add(f2);
  }
  pz->Write();
  pxy->Write();
  hz->Write();
  hxy->Write();
  delete pz;  // also deletes the attached TF1
  delete pxy;
  delete hz;
  delete hxy;
  m_file->cd();
}

bool SiHitSeedQA::hitToGlobal(const SiHitPoint& h, double& gx, double& gy, double& gz) const
{
  if (!m_geom)
  {
    return false;
  }
  auto surf = m_geom->maps().getSiliconSurface(h.hitsetkey);
  if (!surf)
  {
    return false;
  }
  double locx = 0, locz = 0;
  if (h.layer <= 2)
  {
    auto* lg = m_mvtxGeom ? dynamic_cast<CylinderGeom_Mvtx*>(m_mvtxGeom->GetLayerGeom(h.layer)) : nullptr;
    if (!lg)
    {
      return false;
    }
    const TVector3 local = lg->get_local_coords_from_pixel(static_cast<int>(h.row), static_cast<int>(h.col));
    locx = local.X();
    locz = local.Z();
  }
  else
  {
    auto* lg = m_inttGeom ? dynamic_cast<CylinderGeomIntt*>(m_inttGeom->GetLayerGeom(h.layer)) : nullptr;
    if (!lg)
    {
      return false;
    }
    double lc[3] = {0, 0, 0};
    lg->find_strip_center_localcoords(h.ladderz, static_cast<int>(h.row), static_cast<int>(h.col), lc);
    locx = lc[1];
    locz = lc[2];
  }
  const Acts::Vector2 local2D(locx * Acts::UnitConstants::cm, locz * Acts::UnitConstants::cm);
  const Acts::Vector3 g = surf->localToGlobal(m_geom->geometry().getGeoContext(), local2D, Acts::Vector3(1, 1, 1)) / Acts::UnitConstants::cm;
  gx = g.x();
  gy = g.y();
  gz = g.z();
  return true;
}

// Rigid transform global = R * det + T (Kabsch), accumulated over events.
// det point of a hit: (R_layer cos(Uphi + m_rotation), R_layer sin(...), nominal Vz) [cm]
void SiHitSeedQA::fitFrame(const SiHitSeedEvent& d)
{
  evt_fit_npoints = 0;
  evt_fit_rms = kNaN;
  std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> pts;  // (det, global)
  for (const auto& c : d.chains)
  {
    for (int hid : c.point_hit_ids)
    {
      if ((int) pts.size() >= m_frameMaxPoints || hid < 0 || hid >= (int) d.hits.size())
      {
        continue;
      }
      const auto& h = d.hits[hid];
      double gx, gy, gz;
      if (!hitToGlobal(h, gx, gy, gz))
      {
        continue;
      }
      const double r = m_radius[h.layer], ph = h.phi + m_rotation;
      pts.emplace_back(Eigen::Vector3d(r * std::cos(ph), r * std::sin(ph), h.z), Eigen::Vector3d(gx, gy, gz));
    }
  }
  for (const auto& [dp, gp] : pts)
  {
    m_fN += 1;
    for (int i = 0; i < 3; ++i)
    {
      m_fSd[i] += dp(i);
      m_fSg[i] += gp(i);
      for (int j = 0; j < 3; ++j)
      {
        m_fSdg[3 * i + j] += dp(i) * gp(j);
      }
    }
  }
  evt_fit_npoints = static_cast<int>(pts.size());
  if (m_fN < 10)
  {
    return;
  }
  const Eigen::Vector3d cd(m_fSd[0] / m_fN, m_fSd[1] / m_fN, m_fSd[2] / m_fN);
  const Eigen::Vector3d cg(m_fSg[0] / m_fN, m_fSg[1] / m_fN, m_fSg[2] / m_fN);
  Eigen::Matrix3d H;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      H(i, j) = m_fSdg[3 * i + j] / m_fN - cd(i) * cg(j);  // covariance of det (rows) and global (cols)
    }
  }
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d R = svd.matrixV() * svd.matrixU().transpose();
  if (R.determinant() < 0)
  {
    Eigen::Matrix3d V = svd.matrixV();
    V.col(2) *= -1;
    R = V * svd.matrixU().transpose();
  }
  const Eigen::Vector3d T = cg - R * cd;
  for (int i = 0; i < 3; ++i)
  {
    m_T[i] = T(i);
    for (int j = 0; j < 3; ++j)
    {
      m_R[3 * i + j] = R(i, j);
    }
  }
  m_frameOk = true;
  // per-event residual rms of this event's points with the accumulated transform
  double s2 = 0;
  for (const auto& [dp, gp] : pts)
  {
    s2 += (gp - (R * dp + T)).squaredNorm();
  }
  if (!pts.empty())
  {
    evt_fit_rms = 10.0 * std::sqrt(s2 / pts.size());
  }
  det_fit_rms = evt_fit_rms;
  det_fit_npoints = static_cast<int>(m_fN);
  det_T_x = 10.0 * T.x() - m_center[0];  // global position of the detector-frame origin [mm]
  det_T_y = 10.0 * T.y() - m_center[1];
  det_T_z = 10.0 * T.z() - m_center[2];
  // small-angle rotation angles of R [mrad]
  det_rot_x = 1000.0 * 0.5 * (R(2, 1) - R(1, 2));
  det_rot_y = 1000.0 * 0.5 * (R(0, 2) - R(2, 0));
  det_rot_z = 1000.0 * 0.5 * (R(1, 0) - R(0, 1));
}

int SiHitSeedQA::process_event(PHCompositeNode* topNode)
{
  const auto* d = findNode::getClass<SiHitSeedEvent>(topNode, m_inputNodeName);
  if (!d)
  {
    std::cerr << Name() << ": missing input node " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  event = d->event;
  nhits = static_cast<int>(d->hits.size());
  nclusters = static_cast<int>(d->clusters.size());
  nchains = static_cast<int>(d->chains.size());
  vertex_z = d->vertex_z;
  vertex_z_linefit = d->vertex_z_linefit;
  vtx_tr_npairs = d->vertex_tracklet_npairs;
  vtx_tr_npeak = d->vertex_tracklet_npeak;
  fit_a = d->dz_vs_z_intercept_bins;
  fit_b = d->dz_vs_z_slope_bins_per_cm;
  nvertex = d->n_vertex_links;

  // ---- detector <-> global frame, and the standard vertices (untouched + detector frame)
  m_geom = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  m_mvtxGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  m_inttGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  fitFrame(*d);
  if (!m_frameOk)
  {
    det_T_x = det_T_y = det_T_z = det_rot_x = det_rot_y = det_rot_z = det_fit_rms = kNaN;
  }
  for (auto* v : {&svx_x, &svx_y, &svx_z, &svx_ex, &svx_ey, &svx_ez, &svx_chi2, &svx_ndof, &svx_det_x, &svx_det_y, &svx_det_z})
  {
    v->clear();
  }
  svx_ntracks.clear();
  if (auto* vmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName))
  {
    for (auto it = vmap->begin(); it != vmap->end(); ++it)
    {
      const SvtxVertex* sv = it->second;
      if (!sv)
      {
        continue;
      }
      svx_x.push_back(sv->get_x());
      svx_y.push_back(sv->get_y());
      svx_z.push_back(sv->get_z());
      svx_ex.push_back(std::sqrt(std::max(0.0, static_cast<double>(sv->get_error(0, 0)))));
      svx_ey.push_back(std::sqrt(std::max(0.0, static_cast<double>(sv->get_error(1, 1)))));
      svx_ez.push_back(std::sqrt(std::max(0.0, static_cast<double>(sv->get_error(2, 2)))));
      svx_chi2.push_back(sv->get_chisq());
      svx_ndof.push_back(sv->get_ndof());
      svx_ntracks.push_back(static_cast<int>(sv->size_tracks()));
      if (m_frameOk)
      {
        // det = R^T (global - T), then cm -> mm and the optional detector-centre offset
        const double gx = sv->get_x() - m_T[0], gy = sv->get_y() - m_T[1], gz = sv->get_z() - m_T[2];
        svx_det_x.push_back(10.0 * (m_R[0] * gx + m_R[3] * gy + m_R[6] * gz) + m_center[0]);
        svx_det_y.push_back(10.0 * (m_R[1] * gx + m_R[4] * gy + m_R[7] * gz) + m_center[1]);
        svx_det_z.push_back(10.0 * (m_R[2] * gx + m_R[5] * gy + m_R[8] * gz) + m_center[2]);
      }
      else
      {
        svx_det_x.push_back(kNaN);
        svx_det_y.push_back(kNaN);
        svx_det_z.push_back(kNaN);
      }
    }
  }

  // ---- cleaned vertex sample: first select chains, then build links and 3-MVTX tracklets
  double bnd = std::fmod(m_boundary, kTwoPi);
  if (bnd < 0)
  {
    bnd += kTwoPi;
  }

  std::vector<Link> links, linksA, linksB;
  std::vector<Tracklet> tracklets, trackletsA, trackletsB;
  std::vector<double> acceptedPhi;

  n_vertex_chains = 0;

  for (const auto& c : d->chains)
  {
    if (!chainPassesVertexCuts(c))
    {
      continue;
    }

    ++n_vertex_chains;

    Tracklet tr;
    if (makeMvtxTracklet(c, tr))
    {
      tracklets.push_back(tr);
      acceptedPhi.push_back(tr.phi);

      double u = std::fmod(tr.phi - bnd, kTwoPi);
      if (u < 0)
      {
        u += kTwoPi;
      }
      (u < kPi ? trackletsA : trackletsB).push_back(tr);
    }

    for (const auto& s : c.steps)
    {
      const int bit =
          s.kind == SiChainStep::Seed ? 1 :
          (s.kind == SiChainStep::MvtxPropagation ? 2 : 4);

      if (!(m_vtxStepMask & bit) ||
          s.from_layer < 0 ||
          s.to_layer < 0 ||
          s.from_layer == s.to_layer)
      {
        continue;
      }

      if (m_vtxSeedFromLayer >= 0 &&
          m_vtxSeedToLayer >= 0 &&
          (s.from_layer != m_vtxSeedFromLayer ||
           s.to_layer != m_vtxSeedToLayer))
      {
        continue;
      }

      const Link l{
          s.from_layer,
          s.to_layer,
          s.ref_z,
          s.z - s.ref_z,
          s.ref_phi,
          wrapDelta(s.phi - s.ref_phi)};

      links.push_back(l);

      // If no 3-MVTX tracklet was available, still let the link contribute
      // to the azimuthal-coverage diagnostic.
      if (c.n_mvtx < 3)
      {
        acceptedPhi.push_back(l.phi_from);
      }

      double u = std::fmod(l.phi_from - bnd, kTwoPi);
      if (u < 0)
      {
        u += kTwoPi;
      }
      (u < kPi ? linksA : linksB).push_back(l);
    }
  }

  // Original link estimator, but now after chain-level cleanup.
  fitVertex(links, m_vtxAll);
  fitVertex(linksA, m_vtxA);
  fitVertex(linksB, m_vtxB);

  // Stage 7: one measurement per accepted 3-MVTX chain.
  fitTrackletVertex(tracklets, m_trkVtxAll);
  fitTrackletVertex(trackletsA, m_trkVtxA);
  fitTrackletVertex(trackletsB, m_trkVtxB);

  n_vertex_tracklets = static_cast<int>(tracklets.size());

  // Stage 5: phi coverage.
  vertex_phi_sectors = countPhiSectors(acceptedPhi);
  vertex_phi_coverage_ok =
      (m_vtxMinPhiSectors <= 0 ||
       vertex_phi_sectors >= m_vtxMinPhiSectors) ? 1 : 0;

  // Choose which estimator defines the final event vertex / quality.
  const Vertex& qAll = m_useTrackletVertex ? m_trkVtxAll : m_vtxAll;
  const Vertex& qA   = m_useTrackletVertex ? m_trkVtxA   : m_vtxA;
  const Vertex& qB   = m_useTrackletVertex ? m_trkVtxB   : m_vtxB;

  vertex_half_dx = kNaN;
  vertex_half_dy = kNaN;
  vertex_half_agreement_ok = 1;

  // Stage 6: require both halves and compare their independent vertices.
  if (m_vtxMaxHalfDx > 0 || m_vtxMaxHalfDy > 0)
  {
    if ((qA.ok & 2) && (qB.ok & 2))
    {
      vertex_half_dx = qA.x - qB.x;
      vertex_half_dy = qA.y - qB.y;

      if (m_vtxMaxHalfDx > 0 &&
          std::abs(vertex_half_dx) > m_vtxMaxHalfDx)
      {
        vertex_half_agreement_ok = 0;
      }

      if (m_vtxMaxHalfDy > 0 &&
          std::abs(vertex_half_dy) > m_vtxMaxHalfDy)
      {
        vertex_half_agreement_ok = 0;
      }
    }
    else
    {
      vertex_half_agreement_ok = 0;
    }
  }

  final_vtx_x = qAll.x;
  final_vtx_y = qAll.y;
  final_vtx_z = qAll.z;
  final_vtx_x_err = qAll.x_err;
  final_vtx_y_err = qAll.y_err;
  final_vtx_z_err = qAll.z_err;
  final_vtx_ok = qAll.ok;

  vertex_selected =
      ((qAll.ok & 3) == 3 &&
       vertex_phi_coverage_ok &&
       vertex_half_agreement_ok) ? 1 : 0;
  if (m_nSavedProfiles < m_nSaveProfiles && !links.empty())
  {
    saveEventProfiles(event, "all", links, m_vtxAll);
    saveEventProfiles(event, "halfA", linksA, m_vtxA);
    saveEventProfiles(event, "halfB", linksB, m_vtxB);
    ++m_nSavedProfiles;
  }
  m_evt->Fill();

  if (Verbosity() > 0)
  {
    auto print = [&](const char* name, const Vertex& v)
    {
      std::cout << "  " << name << " (" << v.nlinks << " links): z = " << v.z << " +- " << v.z_err
                << " mm (median " << v.z_median << "), x = " << v.x << " +- " << v.x_err
                << ", y = " << v.y << " +- " << v.y_err << " mm" << std::endl;
    };
    std::cout << Name() << ": event " << event << " vertex [mm, detector frame]:" << std::endl;
    print("all  ", m_vtxAll);
    print("halfA", m_vtxA);
    print("halfB", m_vtxB);
    print("trk  ", m_trkVtxAll);
    print("trk A", m_trkVtxA);
    print("trk B", m_trkVtxB);
    std::cout << "  vertex cleanup: chains=" << n_vertex_chains
              << ", tracklets=" << n_vertex_tracklets
              << ", phi sectors=" << vertex_phi_sectors
              << ", phiOK=" << vertex_phi_coverage_ok
              << ", halfOK=" << vertex_half_agreement_ok
              << ", selected=" << vertex_selected
              << ", dHalf=(" << vertex_half_dx << ", " << vertex_half_dy << ") mm"
              << std::endl;
    std::cout << "  detector frame in global [mm]: T = (" << det_T_x << ", " << det_T_y << ", " << det_T_z
              << "), rot [mrad] = (" << det_rot_x << ", " << det_rot_y << ", " << det_rot_z << "), rms "
              << det_fit_rms << " mm from " << det_fit_npoints << " points" << std::endl;
    for (size_t i = 0; i < svx_z.size(); ++i)
    {
      std::cout << "  standard vertex " << i << ": global [cm] (" << svx_x[i] << ", " << svx_y[i] << ", " << svx_z[i]
                << ") -> detector [mm] (" << svx_det_x[i] << ", " << svx_det_y[i] << ", " << svx_det_z[i] << "), "
                << svx_ntracks[i] << " tracks" << std::endl;
    }
  }

  for (const auto& h : d->hits)
  {
    hitsetkey = h.hitsetkey;
    hitkey = h.hitkey;
    layer = h.layer;
    row = static_cast<int>(h.row);
    col = static_cast<int>(h.col);
    stave = h.stave;
    chip = h.chip;
    ladderphi = h.ladderphi;
    ladderz = h.ladderz;
    adc = h.adc;
    lx = h.lx;
    ly = h.ly;
    z = h.z;
    phi = h.phi;
    m_hits->Fill();
  }

  for (const auto& c : d->clusters)
  {
    cluskey = c.key;
    layer = c.layer;
    stave = c.stave;
    chip = c.chip;
    ladderphi = c.ladderphi;
    ladderz = c.ladderz;
    adc = c.adc;
    lx = c.lx;
    ly = c.ly;
    z = c.z;
    phi = c.phi;
    m_clus->Fill();
  }

  for (const auto& c : d->chains)
  {
    chain_id = c.id;
    vlayer = c.layers;
    vz = c.z;
    vphi = c.phi;
    chain_nmvtx = c.n_mvtx;
    chain_nintt = c.n_intt;
    chain_eta = chainEta(c);
    chain_score = c.score;
    chain_max_res_z = chainMaxResidualZ(c);
    chain_max_res_phi = chainMaxResidualPhi(c);
    chain_pass_vertex_cuts = chainPassesVertexCuts(c) ? 1 : 0;
    m_chains->Fill();

    for (const auto& s : c.steps)
    {
      kind = s.kind;
      from_layer = s.from_layer;
      to_layer = s.to_layer;
      ref_z = s.ref_z;
      ref_phi = s.ref_phi;
      pred_z = s.pred_z;
      pred_phi = s.pred_phi;
      z = s.z;
      phi = s.phi;
      delta_z = s.delta_z_bins;
      delta_phi = s.delta_phi_bins;
      res_z = s.res_z_bins;
      res_phi = s.res_phi_bins;
      score = s.score;
      m_steps->Fill();

      if (kind == SiChainStep::Seed)
      {
        m_dz->Fill(ref_z, delta_z);
        m_dzcm->Fill(ref_z, z - ref_z);
        m_dphi->Fill(ref_phi, delta_phi);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SiHitSeedQA::End(PHCompositeNode*)
{
  if (m_file)
  {
    m_file->cd();
    m_file->Write();
    m_file->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}