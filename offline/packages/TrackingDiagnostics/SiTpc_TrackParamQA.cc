#include "SiTpc_TrackParamQA.h"

#include <sitrackreco/SiTpcHelixFit.h>
#include <sitrackreco/SiTpc_Track.h>
#include <sitrackreco/SiTpc_TrackContainer.h>
#include <sitrackreco/Si_Trajectory.h>
#include <sitrackreco/Si_TrajectoryContainer.h>

#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>

namespace
{
  constexpr double kPi = 3.14159265358979323846;

  double wrapPi(double a)
  {
    a = std::fmod(a + kPi, 2.0 * kPi);
    if (a < 0)
    {
      a += 2.0 * kPi;
    }
    return a - kPi;
  }

  bool isCircle(const SiTpcHelixFit::Result& r)
  {
    return r.status == SiTpcHelixFit::Ok && r.R < SiTpcHelixFit::kStraightRadius;
  }

  // unit direction of motion at a point (x, y) on the helix
  void directionAt(const SiTpcHelixFit::Result& r, double x, double y, double& tx, double& ty)
  {
    if (!isCircle(r))
    {
      tx = r.tx;
      ty = r.ty;
      return;
    }
    const double a = std::atan2(y - r.cy, x - r.cx);
    tx = -r.helicity * std::sin(a);
    ty = r.helicity * std::cos(a);
  }

  // helix with its pca moved to the closest approach to (bx, by): pca, phi, dca, z0 updated
  SiTpcHelixFit::Result moveToBeam(const SiTpcHelixFit::Result& in, double bx, double by)
  {
    SiTpcHelixFit::Result r = in;
    double ds = 0;
    if (isCircle(in))
    {
      const double dx = bx - in.cx, dy = by - in.cy;
      const double d = std::hypot(dx, dy);
      if (d <= 0)
      {
        return r;
      }
      r.pcaX = in.cx + in.R * dx / d;
      r.pcaY = in.cy + in.R * dy / d;
      const double a0 = std::atan2(in.pcaY - in.cy, in.pcaX - in.cx);
      const double a1 = std::atan2(r.pcaY - in.cy, r.pcaX - in.cx);
      ds = in.helicity * in.R * wrapPi(a1 - a0);
      directionAt(in, r.pcaX, r.pcaY, r.tx, r.ty);
    }
    else
    {
      ds = (bx - in.pcaX) * in.tx + (by - in.pcaY) * in.ty;
      r.pcaX = in.pcaX + ds * in.tx;
      r.pcaY = in.pcaY + ds * in.ty;
    }
    r.z0 = in.z0 + in.tanl * ds;
    r.phi = std::atan2(r.ty, r.tx);
    // positive if the beam is to the right of the direction of motion
    r.dca = (bx - r.pcaX) * r.ty - (by - r.pcaY) * r.tx;
    return r;
  }

  struct TpcPoint
  {
    double x = 0, y = 0, z = 0, r = 0;
  };
}  // namespace

SiTpc_TrackParamQA::SiTpc_TrackParamQA(const std::string& name, const std::string& outfile)
  : SubsysReco(name)
  , m_outfilename(outfile)
{
}

SiTpc_TrackParamQA::~SiTpc_TrackParamQA() = default;

TH1F* SiTpc_TrackParamQA::h1(const std::string& n, const std::string& t, int nb, double lo, double hi)
{
  auto* h = new TH1F(n.c_str(), t.c_str(), nb, lo, hi);
  m_all.push_back(h);
  return h;
}

TH2F* SiTpc_TrackParamQA::h2(const std::string& n, const std::string& t, int nx, double xlo, double xhi, int ny,
                             double ylo, double yhi)
{
  auto* h = new TH2F(n.c_str(), t.c_str(), nx, xlo, xhi, ny, ylo, yhi);
  m_all.push_back(h);
  return h;
}

void SiTpc_TrackParamQA::book()
{
  const char* cat[kNCat] = {"all", "pos", "neg", "halfA", "halfB"};
  const char* catTitle[kNCat] = {"all", "TPC q>0", "TPC q<0", "Si half A", "Si half B"};
  for (int c = 0; c < kNCat; ++c)
  {
    const std::string p = std::string("h_") + cat[c] + "_";
    const std::string t = std::string(catTitle[c]) + ": ";
    m_hDphi[c] = h1(p + "dphi", t + "#phi_{Si} - #phi_{TPC} at the pca;#Delta#phi [rad]", 600, -0.3, 0.3);
    m_hDeta[c] = h1(p + "deta", t + "#eta_{Si} - #eta_{TPC};#Delta#eta", 600, -0.3, 0.3);
    m_hDtanl[c] = h1(p + "dtanl", t + "tan#lambda_{Si} - tan#lambda_{TPC};#Deltatan#lambda", 600, -0.3, 0.3);
    m_hDz0[c] = h1(p + "dz0", t + "z0_{Si} - z0_{TPC};#Deltaz0 [cm]", 1000, -25, 25);
    m_hDdca[c] = h1(p + "ddca", t + "dca_{Si} - dca_{TPC};#Deltadca [cm]", 800, -2, 2);
    m_hDqopt[c] = h1(p + "dqopt", t + "(q/p_{T})_{Si} - (q/p_{T})_{TPC} (Si circle);#Delta q/p_{T} [1/GeV]", 800, -4, 4);
    m_hDptrel[c] = h1(p + "dptrel", t + "(p_{T,Si} - p_{T,TPC}) / p_{T,TPC} (Si circle);#Deltap_{T}/p_{T}", 800, -2, 4);
  }

  const char* var[5] = {"dphi", "deta", "dz0", "ddca", "dqopt"};
  const char* varT[5] = {"#Delta#phi [rad]", "#Delta#eta", "#Deltaz0 [cm]", "#Deltadca [cm]", "#Delta q/p_{T} [1/GeV]"};
  const int vnb[5] = {300, 300, 500, 400, 400};
  const double vlo[5] = {-0.3, -0.3, -25, -2, -4};
  const double vhi[5] = {0.3, 0.3, 25, 2, 4};
  const char* xs[5] = {"tpcpt", "tpcphi", "tpceta", "tpcz0", "tpcqopt"};
  const char* xT[5] = {"p_{T,TPC} [GeV]", "#phi_{TPC} [rad]", "#eta_{TPC}", "z0_{TPC} [cm]", "(q/p_{T})_{TPC} [1/GeV]"};
  const int xnb[5] = {50, 72, 48, 120, 80};
  const double xlo[5] = {0, -kPi, -1.2, -120, -4};
  const double xhi[5] = {10, kPi, 1.2, 60, 4};
  for (int v = 0; v < 5; ++v)
  {
    for (int x = 0; x < 5; ++x)
    {
      m_h2Var[v][x] = h2(std::string("h2_") + var[v] + "_vs_" + xs[x], std::string(";") + xT[x] + ";" + varT[v], xnb[x],
                         xlo[x], xhi[x], vnb[v], vlo[v], vhi[v]);
    }
    for (int c = 1; c < kNCat; ++c)
    {
      m_h2VarPhiCat[v][c - 1] =
          h2(std::string("h2_") + var[v] + "_vs_tpcphi_" + cat[c], std::string(catTitle[c]) + ";#phi_{TPC} [rad];" + varT[v],
             72, -kPi, kPi, vnb[v], vlo[v], vhi[v]);
    }
  }

  m_hCorr[0] = h2("h2_corr_pt", ";p_{T,TPC} [GeV];p_{T,Si} [GeV]", 100, 0, 10, 100, 0, 10);
  m_hCorr[1] = h2("h2_corr_qopt", ";(q/p_{T})_{TPC};(q/p_{T})_{Si}", 160, -4, 4, 160, -4, 4);
  m_hCorr[2] = h2("h2_corr_phi", ";#phi_{TPC};#phi_{Si}", 144, -kPi, kPi, 144, -kPi, kPi);
  m_hCorr[3] = h2("h2_corr_eta", ";#eta_{TPC};#eta_{Si}", 120, -1.5, 1.5, 120, -1.5, 1.5);
  m_hCorr[4] = h2("h2_corr_z0", ";z0_{TPC} [cm];z0_{Si} [cm]", 180, -120, 60, 120, -30, 30);
  m_hCorr[5] = h2("h2_corr_dca", ";dca_{TPC} [cm];dca_{Si} [cm]", 100, -2, 2, 100, -2, 2);
  m_hCharge = h2("h2_charge", ";TPC charge;Si charge (0 = straight line)", 3, -1.5, 1.5, 3, -1.5, 1.5);

  const int nR = std::min<int>(kMaxRadii, static_cast<int>(m_radii.size()));
  for (int k = 0; k < nR; ++k)
  {
    const std::string p = "h_R" + std::to_string(k) + "_";
    const std::string p2 = "h2_R" + std::to_string(k) + "_";
    const std::string t = "r = " + std::to_string(m_radii[k]).substr(0, 6) + " cm: ";
    m_hRRdphi[k] = h1(p + "rdphi", t + "Si - TPC helix position;r#Delta#phi [cm]", 800, -2, 2);
    m_hRDz[k] = h1(p + "dz", t + "Si - TPC helix position;#Deltaz [cm]", 1000, -25, 25);
    m_hRDdir[k] = h1(p + "ddir", t + "Si - TPC direction;#Delta#phi_{dir} [rad]", 600, -0.3, 0.3);
    m_hRRdphiPhi[k] = h2(p2 + "rdphi_vs_phi", t + ";#phi [rad];r#Delta#phi [cm]", 72, -kPi, kPi, 400, -2, 2);
    m_hRDzZ[k] = h2(p2 + "dz_vs_z", t + ";z_{TPC} [cm];#Deltaz [cm]", 60, -30, 30, 500, -25, 25);
    m_hRRdphiQopt[k] = h2(p2 + "rdphi_vs_qopt", t + ";(q/p_{T})_{TPC} [1/GeV];r#Delta#phi [cm]", 80, -4, 4, 400, -2, 2);
  }

  for (int L = 0; L < 7; ++L)
  {
    const std::string p = "h_L" + std::to_string(L) + "_";
    const std::string p2 = "h2_L" + std::to_string(L) + "_";
    const std::string t = "layer " + std::to_string(L) + ": Si point - TPC-only helix";
    const double zr = L <= 2 ? 15 : 25;
    m_hLRdphi[L] = h1(p + "rdphi", t + ";r#Delta#phi [cm]", 800, -2, 2);
    m_hLDz[L] = h1(p + "dz", t + ";#Deltaz [cm]", 1000, -25, 25);
    m_hLRdphiPhi[L] = h2(p2 + "rdphi_vs_phi", t + ";#phi_{point} [rad];r#Delta#phi [cm]", 72, -kPi, kPi, 400, -2, 2);
    m_hLDzZ[L] = h2(p2 + "dz_vs_z", t + ";z_{point} [cm];#Deltaz [cm]", 60, -zr, zr, 500, -25, 25);
    m_hLRdphiQopt[L] = h2(p2 + "rdphi_vs_qopt", t + ";(q/p_{T})_{TPC} [1/GeV];r#Delta#phi [cm]", 80, -4, 4, 400, -2, 2);
    m_hLDzPhi[L] = h2(p2 + "dz_vs_phi", t + ";#phi_{point} [rad];#Deltaz [cm]", 72, -kPi, kPi, 500, -25, 25);
  }
}

int SiTpc_TrackParamQA::Init(PHCompositeNode* /*topNode*/)
{
  m_file = TFile::Open(m_outfilename.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cerr << Name() << ": cannot open " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_file->cd();
  book();

  if (m_writeTree)
  {
    m_tree = new TTree("params", "Si-only vs TPC-only parameters of matched pairs (Delta = Si - TPC)");
    m_tree->Branch("event", &t_event);
    m_tree->Branch("track_id", &t_trackId);
    m_tree->Branch("si_index", &t_siIndex);
    m_tree->Branch("tpc_index", &t_tpcIndex);
    m_tree->Branch("si_chain", &t_siChain);
    m_tree->Branch("si_status", &t_siStatus);
    m_tree->Branch("half", &t_half);
    m_tree->Branch("n_mvtx", &t_nMvtx);
    m_tree->Branch("n_intt", &t_nIntt);
    m_tree->Branch("n_tpc_clusters", &t_nTpcClus);
    m_tree->Branch("match_chi2", &t_matchChi2);
    m_tree->Branch("match_dphi", &t_matchDphi);
    m_tree->Branch("match_deta", &t_matchDeta);
    m_tree->Branch("match_dz0", &t_matchDz0);
    m_tree->Branch("fitted_tpc_z_offset", &t_fittedZOffset);
    m_tree->Branch("si_phi", &t_siPhi);
    m_tree->Branch("si_eta", &t_siEta);
    m_tree->Branch("si_tanl", &t_siTanl);
    m_tree->Branch("si_z0", &t_siZ0);
    m_tree->Branch("si_dca", &t_siDca);
    m_tree->Branch("si_pt", &t_siPt);
    m_tree->Branch("si_qopt", &t_siQopt);
    m_tree->Branch("si_R", &t_siR);
    m_tree->Branch("si_charge", &t_siCharge);
    m_tree->Branch("tpc_phi", &t_tpcPhi);
    m_tree->Branch("tpc_eta", &t_tpcEta);
    m_tree->Branch("tpc_tanl", &t_tpcTanl);
    m_tree->Branch("tpc_z0", &t_tpcZ0);
    m_tree->Branch("tpc_dca", &t_tpcDca);
    m_tree->Branch("tpc_pt", &t_tpcPt);
    m_tree->Branch("tpc_qopt", &t_tpcQopt);
    m_tree->Branch("tpc_R", &t_tpcR);
    m_tree->Branch("tpc_charge", &t_tpcCharge);
    m_tree->Branch("dphi", &t_dphi);
    m_tree->Branch("deta", &t_deta);
    m_tree->Branch("dtanl", &t_dtanl);
    m_tree->Branch("dz0", &t_dz0);
    m_tree->Branch("ddca", &t_ddca);
    m_tree->Branch("dqopt", &t_dqopt);
    m_tree->Branch("dptrel", &t_dptrel);
    m_tree->Branch("nR", &t_nR);
    m_tree->Branch("R", t_R, "R[nR]/F");
    m_tree->Branch("R_phi", t_Rphi, "R_phi[nR]/F");
    m_tree->Branch("R_z", t_Rz, "R_z[nR]/F");
    m_tree->Branch("R_rdphi", t_Rrdphi, "R_rdphi[nR]/F");
    m_tree->Branch("R_dz", t_Rdz, "R_dz[nR]/F");
    m_tree->Branch("R_ddir", t_Rddir, "R_ddir[nR]/F");
    m_tree->Branch("p_layer", &t_pLayer);
    m_tree->Branch("p_r", &t_pR);
    m_tree->Branch("p_phi", &t_pPhi);
    m_tree->Branch("p_z", &t_pZ);
    m_tree->Branch("p_rdphi", &t_pRdphi);
    m_tree->Branch("p_dz", &t_pDz);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": beam (" << beamX() << ", " << beamY() << ") cm, TPC z offset " << m_tpcZOffset
              << " cm" << (m_useFittedZOffset ? " + fitted per track" : "") << ", max match chi2 " << m_maxMatchChi2
              << ", output " << m_outfilename << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void SiTpc_TrackParamQA::resetTree()
{
  t_nR = 0;
  t_pLayer.clear();
  t_pR.clear();
  t_pPhi.clear();
  t_pZ.clear();
  t_pRdphi.clear();
  t_pDz.clear();
}

int SiTpc_TrackParamQA::process_event(PHCompositeNode* topNode)
{
  ++m_nEvents;
  auto* tracks = findNode::getClass<SiTpc_TrackContainer>(topNode, m_trackNodeName);
  auto* siTraj = findNode::getClass<Si_TrajectoryContainer>(topNode, m_siTrajNodeName);
  auto* clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  auto* tpcTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName);
  if (!tracks || !siTraj || !clusters || !tpcTracks)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing input (" << m_trackNodeName << (tracks ? " ok" : " missing") << ", "
                << m_siTrajNodeName << (siTraj ? " ok" : " missing") << ", " << m_tpcClusterNodeName
                << (clusters ? " ok" : " missing") << ", " << m_tpcTrackNodeName << (tpcTracks ? " ok" : " missing")
                << ")" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // TPC clusters per assembled track (as in SiTpcTrackMatcher)
  std::map<unsigned int, std::vector<TpcPoint>> clustersById;
  for (unsigned int i = 0; i < clusters->size(); ++i)
  {
    const Tpc_PolyCluster* c = clusters->get_cluster(i);
    if (!c || !c->isValid())
    {
      continue;
    }
    const auto b = m_alignment.tpcToTpcFrame(c->get_centroid_x(), c->get_centroid_y(), c->get_centroid_z());
    if (!std::isfinite(b.x) || !std::isfinite(b.y) || !std::isfinite(b.z))
    {
      continue;
    }
    clustersById[c->get_source_assembled_track_id()].push_back({b.x, b.y, b.z, std::hypot(b.x, b.y)});
  }

  SiTpcHelixFit::Config tpcCfg;
  tpcCfg.beamX = beamX();
  tpcCfg.beamY = beamY();
  tpcCfg.bz = m_bz;

  const double bx = beamX(), by = beamY();
  const int nR = std::min<int>(kMaxRadii, static_cast<int>(m_radii.size()));
  unsigned long nUsedEvent = 0;

  for (unsigned int it = 0; it < tracks->size(); ++it)
  {
    const SiTpc_Track* m = tracks->get_track(it);
    if (!m)
    {
      continue;
    }
    ++m_nTracks;
    const Si_Trajectory* s = siTraj->get_trajectory(m->get_si_trajectory_index());
    if (!s || !s->is_fitted())
    {
      ++m_nNoSi;
      continue;
    }
    const Tpc_PolyTrack* tt = tpcTracks->get_track(m->get_tpc_track_index());
    const auto cit = tt ? clustersById.find(tt->get_source_assembled_track_id()) : clustersById.end();
    if (cit == clustersById.end() || cit->second.size() < m_minTpcClusters)
    {
      ++m_nNoTpc;
      continue;
    }

    // ---- selection
    if (!(m->get_match_chi2() <= m_maxMatchChi2) || s->get_n_mvtx() < m_minMvtx || s->get_n_intt() < m_minIntt ||
        s->size_points() < m_minSiPoints || (m_requireSiCircle && s->get_fit_status() != Si_Trajectory::Ok))
    {
      ++m_nSel;
      continue;
    }

    // ---- TPC-only refit (identical to the matcher), TPC z corrected
    std::vector<TpcPoint> tp = cit->second;
    std::sort(tp.begin(), tp.end(), [](const auto& a, const auto& b) { return a.r < b.r; });
    const double zCorr = m_tpcZOffset + (m_useFittedZOffset ? m->get_tpc_z_offset() : 0.0);
    std::vector<SiTpcHelixFit::Point> pts;
    pts.reserve(tp.size());
    for (const auto& p : tp)
    {
      pts.push_back({p.x, p.y, p.z - zCorr, 1.0, 1.0, 1});
    }
    const SiTpcHelixFit::Result tpcRaw = SiTpcHelixFit::fit(pts, tpcCfg);
    if (!tpcRaw.fitted())
    {
      ++m_nNoTpc;
      continue;
    }
    if (tpcRaw.pt < m_minTpcPt || tpcRaw.pt > m_maxTpcPt)
    {
      ++m_nSel;
      continue;
    }

    // ---- both helices at the same reference point
    const SiTpcHelixFit::Result tpc = moveToBeam(tpcRaw, bx, by);
    const SiTpcHelixFit::Result si = moveToBeam(SiTpcHelixFit::fromTrack(*s), bx, by);
    const bool siCircle = isCircle(si);

    const double tpcQopt = tpc.charge / tpc.pt;
    const double siQopt = siCircle ? si.charge / si.pt : 0.0;
    const double dphi = wrapPi(si.phi - tpc.phi);
    const double deta = std::asinh(si.tanl) - std::asinh(tpc.tanl);
    const double dtanl = si.tanl - tpc.tanl;
    const double dz0 = si.z0 - tpc.z0;
    const double ddca = si.dca - tpc.dca;
    const double dqopt = siCircle ? siQopt - tpcQopt : NAN;
    const double dptrel = siCircle ? (si.pt - tpc.pt) / tpc.pt : NAN;
    const double tpcEta = std::asinh(tpc.tanl);

    // clamshell half of the innermost Si point
    int half = 0;
    {
      unsigned int kin = s->size_points();
      double rmin = 1e9;
      for (unsigned int k = 0; k < s->size_points(); ++k)
      {
        const double r = std::hypot(s->get_point_x(k), s->get_point_y(k));
        if (r < rmin)
        {
          rmin = r;
          kin = k;
        }
      }
      if (kin < s->size_points())
      {
        half = m_alignment.siHalfOfGlobal(m_frame, {s->get_point_x(kin), s->get_point_y(kin), s->get_point_z(kin)});
      }
    }

    bool inCat[kNCat] = {true, tpc.charge > 0, tpc.charge < 0, half == SiTpcBeamAlignment::SiHalfA,
                         half == SiTpcBeamAlignment::SiHalfB};
    for (int c = 0; c < kNCat; ++c)
    {
      if (!inCat[c])
      {
        continue;
      }
      m_hDphi[c]->Fill(dphi);
      m_hDeta[c]->Fill(deta);
      m_hDtanl[c]->Fill(dtanl);
      m_hDz0[c]->Fill(dz0);
      m_hDdca[c]->Fill(ddca);
      if (siCircle)
      {
        m_hDqopt[c]->Fill(dqopt);
        m_hDptrel[c]->Fill(dptrel);
      }
    }
    const double vals[5] = {dphi, deta, dz0, ddca, dqopt};
    const double xv[5] = {tpc.pt, tpc.phi, tpcEta, tpc.z0, tpcQopt};
    for (int v = 0; v < 5; ++v)
    {
      if (!std::isfinite(vals[v]))
      {
        continue;
      }
      for (int x = 0; x < 5; ++x)
      {
        m_h2Var[v][x]->Fill(xv[x], vals[v]);
      }
      for (int c = 1; c < kNCat; ++c)
      {
        if (inCat[c])
        {
          m_h2VarPhiCat[v][c - 1]->Fill(tpc.phi, vals[v]);
        }
      }
    }
    if (siCircle)
    {
      m_hCorr[0]->Fill(tpc.pt, si.pt);
      m_hCorr[1]->Fill(tpcQopt, siQopt);
    }
    m_hCorr[2]->Fill(tpc.phi, si.phi);
    m_hCorr[3]->Fill(tpcEta, std::asinh(si.tanl));
    m_hCorr[4]->Fill(tpc.z0, si.z0);
    m_hCorr[5]->Fill(tpc.dca, si.dca);
    m_hCharge->Fill(tpc.charge, siCircle ? si.charge : 0);

    // ---- tree: scalars
    resetTree();
    t_event = m_nEvents - 1;
    t_trackId = m->get_id();
    t_siIndex = m->get_si_trajectory_index();
    t_tpcIndex = m->get_tpc_track_index();
    t_siChain = s->get_chain_id();
    t_siStatus = s->get_fit_status();
    t_half = half;
    t_nMvtx = s->get_n_mvtx();
    t_nIntt = s->get_n_intt();
    t_nTpcClus = tp.size();
    t_matchChi2 = m->get_match_chi2();
    t_matchDphi = m->get_match_dphi();
    t_matchDeta = m->get_match_deta();
    t_matchDz0 = m->get_match_dz0();
    t_fittedZOffset = m->get_tpc_z_offset();
    t_siPhi = si.phi;
    t_siEta = std::asinh(si.tanl);
    t_siTanl = si.tanl;
    t_siZ0 = si.z0;
    t_siDca = si.dca;
    t_siPt = siCircle ? si.pt : NAN;
    t_siQopt = siQopt;
    t_siR = si.R;
    t_siCharge = siCircle ? si.charge : 0;
    t_tpcPhi = tpc.phi;
    t_tpcEta = tpcEta;
    t_tpcTanl = tpc.tanl;
    t_tpcZ0 = tpc.z0;
    t_tpcDca = tpc.dca;
    t_tpcPt = tpc.pt;
    t_tpcQopt = tpcQopt;
    t_tpcR = tpc.R;
    t_tpcCharge = tpc.charge;
    t_dphi = dphi;
    t_deta = deta;
    t_dtanl = dtanl;
    t_dz0 = dz0;
    t_ddca = ddca;
    t_dqopt = dqopt;
    t_dptrel = dptrel;

    // ---- at the reference radii
    for (int k = 0; k < nR; ++k)
    {
      const double rad = m_radii[k];
      SiTpcHelixFit::Pos ps, pt;
      double ss = 0, st = 0;
      if (!SiTpcHelixFit::atRadius(tpc, rad, tpc.z0 + tpc.tanl * rad, pt, st) ||
          !SiTpcHelixFit::atRadius(si, rad, pt.z, ps, ss))
      {
        continue;
      }
      const double phT = std::atan2(pt.y, pt.x);
      const double rdphi = rad * wrapPi(std::atan2(ps.y, ps.x) - phT);
      const double dz = ps.z - pt.z;
      double txs = 0, tys = 0, txt = 0, tyt = 0;
      directionAt(si, ps.x, ps.y, txs, tys);
      directionAt(tpc, pt.x, pt.y, txt, tyt);
      const double ddir = wrapPi(std::atan2(tys, txs) - std::atan2(tyt, txt));
      m_hRRdphi[k]->Fill(rdphi);
      m_hRDz[k]->Fill(dz);
      m_hRDdir[k]->Fill(ddir);
      m_hRRdphiPhi[k]->Fill(phT, rdphi);
      m_hRDzZ[k]->Fill(pt.z, dz);
      m_hRRdphiQopt[k]->Fill(tpcQopt, rdphi);
      t_R[t_nR] = rad;
      t_Rphi[t_nR] = phT;
      t_Rz[t_nR] = pt.z;
      t_Rrdphi[t_nR] = rdphi;
      t_Rdz[t_nR] = dz;
      t_Rddir[t_nR] = ddir;
      ++t_nR;
    }

    // ---- every Si point against the TPC-only helix
    for (unsigned int k = 0; k < s->size_points(); ++k)
    {
      const int L = s->get_point_layer(k);
      const double x = s->get_point_x(k), y = s->get_point_y(k), z = s->get_point_z(k);
      const double r = std::hypot(x, y);
      SiTpcHelixFit::Pos pt;
      double st = 0;
      if (L < 0 || L > 6 || !SiTpcHelixFit::atRadius(tpc, r, z, pt, st))
      {
        continue;
      }
      const double phi = std::atan2(y, x);
      const double rdphi = r * wrapPi(phi - std::atan2(pt.y, pt.x));
      const double dz = z - pt.z;
      m_hLRdphi[L]->Fill(rdphi);
      m_hLDz[L]->Fill(dz);
      m_hLRdphiPhi[L]->Fill(phi, rdphi);
      m_hLDzZ[L]->Fill(z, dz);
      m_hLRdphiQopt[L]->Fill(tpcQopt, rdphi);
      m_hLDzPhi[L]->Fill(phi, dz);
      t_pLayer.push_back(L);
      t_pR.push_back(r);
      t_pPhi.push_back(phi);
      t_pZ.push_back(z);
      t_pRdphi.push_back(rdphi);
      t_pDz.push_back(dz);
    }

    if (m_tree)
    {
      m_tree->Fill();
    }
    ++m_nUsed;
    ++nUsedEvent;

    if (Verbosity() > 1)
    {
      std::cout << Name() << ": ev " << t_event << " match " << t_trackId << " (si " << t_siIndex << ", tpc "
                << t_tpcIndex << ") half " << half << " nMVTX/nINTT " << t_nMvtx << "/" << t_nIntt
                << " | TPC pt " << tpc.pt << " q " << tpc.charge << " phi " << tpc.phi << " eta " << tpcEta
                << " z0 " << tpc.z0 << " dca " << tpc.dca << " | Si-TPC dphi " << dphi << " deta " << deta
                << " dz0 " << dz0 << " ddca " << ddca << " dq/pt " << dqopt << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << m_nEvents - 1 << " matched pairs " << tracks->size() << ", used "
              << nUsedEvent << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpc_TrackParamQA::End(PHCompositeNode* /*topNode*/)
{
  std::cout << Name() << ": events " << m_nEvents << ", matched pairs " << m_nTracks << ", used " << m_nUsed
            << " (no Si fit " << m_nNoSi << ", no TPC fit " << m_nNoTpc << ", failed selection " << m_nSel << ")"
            << std::endl;
  if (m_nUsed > 0)
  {
    auto core = [](TH1* h)
    {
      // mean / rms of the central part (+-2 rms around the peak, twice): robust to the pedestal
      double c = h->GetBinCenter(h->GetMaximumBin()), w = h->GetRMS();
      for (int i = 0; i < 2; ++i)
      {
        h->GetXaxis()->SetRangeUser(c - 2 * w, c + 2 * w);
        c = h->GetMean();
        w = std::max(h->GetRMS(), h->GetXaxis()->GetBinWidth(1));
      }
      h->GetXaxis()->SetRange(0, 0);
      return std::make_pair(c, w);
    };
    std::cout << Name() << ": core mean / rms (Si - TPC), all pairs:" << std::endl;
    const std::pair<const char*, TH1*> list[] = {
        {"dphi [mrad]", m_hDphi[0]}, {"deta [1e-3]", m_hDeta[0]}, {"dz0 [mm]", m_hDz0[0]}, {"ddca [mm]", m_hDdca[0]},
        {"dq/pt [1/GeV]", m_hDqopt[0]}};
    const double scale[] = {1e3, 1e3, 10, 10, 1};
    for (int i = 0; i < 5; ++i)
    {
      const auto [c, w] = core(list[i].second);
      std::cout << "    " << list[i].first << ": " << scale[i] * c << " / " << scale[i] * w << std::endl;
    }
    for (int L = 0; L < 7; ++L)
    {
      if (m_hLRdphi[L]->GetEntries() < 1)
      {
        continue;
      }
      const auto [c1, w1] = core(m_hLRdphi[L]);
      const auto [c2, w2] = core(m_hLDz[L]);
      std::cout << "    layer " << L << " (" << m_hLRdphi[L]->GetEntries() << " points): r*dphi " << 1e4 * c1 << " / "
                << 1e4 * w1 << " um, dz " << 10 * c2 << " / " << 10 * w2 << " mm" << std::endl;
    }
  }
  if (m_file)
  {
    m_file->cd();
    for (auto* h : m_all)
    {
      h->Write();
    }
    if (m_tree)
    {
      m_tree->Write();
    }
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
