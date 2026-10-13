// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKINGDIAGNOSTICS_SITPC_TRACKPARAMQA_H
#define TRACKINGDIAGNOSTICS_SITPC_TRACKPARAMQA_H

#include <sitrackreco/SiDetectorFrame.h>
#include <sitrackreco/SiTpcBeamAlignment.h>

#include <fun4all/SubsysReco.h>

#include <array>
#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TH1;
class TH1F;
class TH2F;
class TTree;

// SiTpc_TrackParamQA
//
// Old seeding (SiHitSeedReco -> SiTrajectoryFitter -> SiTpcTrackMatcher): for every MATCHED
// pair in SITPC_TRACKS, compares the Si-only trajectory (SI_TRAJECTORY) with the TPC-only
// helix of the same pair.
//
// Inputs: SITPC_TRACKS, SI_TRAJECTORY, TPC_POLYCLUSTERS, TPC_POLYTRACKS.
//
// 1. TPC-only helix: refitted from the TPC clusters exactly as SiTpcTrackMatcher does it
//    (SiTpcHelixFit, TPC frame, same beam position -> use the same setBeamPositionCm as in the
//    matcher).  Si-only helix: the stored Si_Trajectory fit.
// 2. Both helices are re-expressed at the SAME reference point (the beam position of this
//    module): pca, phi (direction at the pca), signed dca, z0 (z at the pca), tanl, eta,
//    q/pt.  So Delta dca and Delta z0 do not depend on which beam each fitter used.
// 3. Both helices are evaluated at fixed reference radii (setReferenceRadii, default MVTX L0,
//    MVTX L2, INTT L3, INTT L6, TPC inner 30 cm): r*Delta phi of the position, Delta z,
//    Delta of the direction.  The inner radii are where the TPC is extrapolated: this is the
//    number the TPC-driven matcher needs.
// 4. Every Si point of the trajectory against the TPC-only helix at the point radius
//    (per layer): this is the old-chain equivalent of the "probe" of SiTpcHitChainMatcher.
//
// Every Delta = Si - TPC.  TPC z: as measured minus setTpcZOffset (and minus the per-track
// fitted offset of the combined fit if setUseFittedTpcZOffset(true)).
//
// Categories for the 1D histograms: all, q+ (TPC charge), q-, half A, half B (clamshell half of
// the innermost Si point; SiTpcBeamAlignment boundary, SiDetectorFrame rotation / centre).
//
// Output file:
//   1D  h_<cat>_dphi, _deta, _dtanl, _dz0, _ddca, _dqopt, _dptrel           cat = all,pos,neg,halfA,halfB
//   2D  h2_<var>_vs_<x>   var = dphi,deta,dz0,ddca,dqopt ; x = tpcpt, tpcphi, tpceta, tpcz0, tpcqopt
//       h2_<var>_vs_tpcphi_<pos|neg|halfA|halfB>
//       h2_corr_<par>     Si vs TPC: pt, qopt, phi, eta, z0, dca
//       h2_charge         Si charge vs TPC charge (0 = straight-line Si fit)
//   at radius k:  h_R<k>_rdphi, _dz, _ddir, h2_R<k>_rdphi_vs_phi, h2_R<k>_dz_vs_z, h2_R<k>_rdphi_vs_qopt
//   per layer L:  h_L<L>_rdphi, _dz, h2_L<L>_rdphi_vs_phi, h2_L<L>_dz_vs_z, h2_L<L>_rdphi_vs_qopt,
//                 h2_L<L>_dz_vs_phi
//   tree "params": one entry per matched pair, all parameters, deltas, at-radius and per-point
//                  arrays - for any cut not foreseen here.
//
// Reading it:
//   - dphi / deta are the matching variables: their histograms are cut by the match windows
//     (in the macro 0.2 / 0.1).  Look at their CORE (width, mean) not the tails.
//     dz0, ddca, dqopt, dtanl, at-radius and per-point residuals are NOT cut on.
//   - accidental matches give a flat pedestal in dz0 / ddca / dqopt: fit core + flat.
//   - geometry: mean r*dphi at a radius vs phi = const (rotation) or a sine (shift);
//     same for both charges.  TPC distortion: opposite sign for q+ / q-, grows with 1/pt
//     (h2_..._rdphi_vs_qopt).  t0 / crossing: Delta z0 constant, h2_dz0_vs_tpcz0 shows the
//     out-of-time TPC tracks.  Drift velocity: Delta z proportional to z (h2_L<L>_dz_vs_z).
//   - half A vs half B: different mean dphi / ddca -> clamshell alignment.
class SiTpc_TrackParamQA : public SubsysReco
{
 public:
  explicit SiTpc_TrackParamQA(const std::string& name = "SiTpc_TrackParamQA",
                              const std::string& outfile = "sitpc_trackparam_qa.root");
  ~SiTpc_TrackParamQA() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  // ---- nodes
  void setOutputFileName(const std::string& s) { m_outfilename = s; }
  void setTrackNodeName(const std::string& s) { m_trackNodeName = s; }
  void setSiTrajectoryNodeName(const std::string& s) { m_siTrajNodeName = s; }
  void setTpcClusterNodeName(const std::string& s) { m_tpcClusterNodeName = s; }
  void setTpcTrackNodeName(const std::string& s) { m_tpcTrackNodeName = s; }

  // ---- frame (same as SiTpcTrackMatcher / SiTrajectoryFitter)
  void setBeamAlignment(const SiTpcBeamAlignment& a) { m_alignment = a; }
  void setDetectorFrame(const SiDetectorFrame& f) { m_frame = f; }
  void setUphiRotation(double rad) { m_frame.setUphiRotation(rad); }
  void setDetectorCenterMm(double x, double y) { m_frame.setDetectorCenterMm(x, y, 0.0); }
  // reference point for the perigee parameters and the TPC refit (as SiTpcTrackMatcher)
  void setBeamPositionCm(double x, double y)
  {
    m_beamX = x;
    m_beamY = y;
    m_beamFromAlignment = false;
  }
  double beamX() const { return m_beamFromAlignment ? m_alignment.referenceBeamX() : m_beamX; }
  double beamY() const { return m_beamFromAlignment ? m_alignment.referenceBeamY() : m_beamY; }
  void setBz(double tesla) { m_bz = tesla; }
  void setMinTpcClusters(unsigned int n) { m_minTpcClusters = n; }  // same as the matcher

  // ---- TPC z
  void setTpcZOffset(double cm) { m_tpcZOffset = cm; }               // z_tpc -> z_tpc - offset
  void setUseFittedTpcZOffset(bool v) { m_useFittedZOffset = v; }    // also - SiTpc_Track::get_tpc_z_offset

  // ---- selection of the matched pairs
  void setMaxMatchChi2(double v) { m_maxMatchChi2 = v; }  // default: all accepted pairs
  void setMinTpcPt(double v) { m_minTpcPt = v; }
  void setMaxTpcPt(double v) { m_maxTpcPt = v; }
  void setMinSiMvtx(unsigned int n) { m_minMvtx = n; }
  void setMinSiIntt(unsigned int n) { m_minIntt = n; }
  void setMinSiPoints(unsigned int n) { m_minSiPoints = n; }
  void setRequireSiCircle(bool v) { m_requireSiCircle = v; }  // only Si fits with curvature

  // ---- reference radii [cm] (max 8)
  void setReferenceRadii(const std::vector<double>& r) { m_radii = r; }

  void setWriteTree(bool v) { m_writeTree = v; }

 private:
  struct Perigee
  {
    bool ok = false;
    bool circle = false;
    double phi = 0, dca = 0, z0 = 0, tanl = 0, eta = 0;
    double pt = 0, qopt = 0;
    int charge = 0;
    double R = 0, pcaX = 0, pcaY = 0;
  };
  static constexpr int kMaxRadii = 8;
  static constexpr int kNCat = 5;  // all, pos, neg, halfA, halfB

  void book();
  TH1F* h1(const std::string& n, const std::string& t, int nb, double lo, double hi);
  TH2F* h2(const std::string& n, const std::string& t, int nx, double xlo, double xhi, int ny, double ylo, double yhi);
  void resetTree();

  std::string m_outfilename;
  std::string m_trackNodeName{"SITPC_TRACKS"};
  std::string m_siTrajNodeName{"SI_TRAJECTORY"};
  std::string m_tpcClusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_tpcTrackNodeName{"TPC_POLYTRACKS"};

  SiTpcBeamAlignment m_alignment;
  SiDetectorFrame m_frame;
  double m_beamX{0}, m_beamY{0};
  bool m_beamFromAlignment{true};
  double m_bz{1.4};
  unsigned int m_minTpcClusters{10};

  double m_tpcZOffset{0};
  bool m_useFittedZOffset{false};

  double m_maxMatchChi2{1e30};
  double m_minTpcPt{0}, m_maxTpcPt{1e30};
  unsigned int m_minMvtx{0}, m_minIntt{0}, m_minSiPoints{0};
  bool m_requireSiCircle{false};

  std::vector<double> m_radii{2.523, 4.148, 7.1844, 10.2584, 30.0};
  bool m_writeTree{true};

  TFile* m_file{nullptr};
  TTree* m_tree{nullptr};
  std::vector<TH1*> m_all;

  // 1D per category
  std::array<TH1F*, kNCat> m_hDphi{}, m_hDeta{}, m_hDtanl{}, m_hDz0{}, m_hDdca{}, m_hDqopt{}, m_hDptrel{};
  // 2D vs x (0 tpcpt, 1 tpcphi, 2 tpceta, 3 tpcz0, 4 tpcqopt), var (0 dphi,1 deta,2 dz0,3 ddca,4 dqopt)
  std::array<std::array<TH2F*, 5>, 5> m_h2Var{};
  // var vs tpc phi per category (pos, neg, halfA, halfB)
  std::array<std::array<TH2F*, 4>, 5> m_h2VarPhiCat{};
  std::array<TH2F*, 6> m_hCorr{};
  TH2F* m_hCharge{nullptr};
  // at radius
  std::array<TH1F*, kMaxRadii> m_hRRdphi{}, m_hRDz{}, m_hRDdir{};
  std::array<TH2F*, kMaxRadii> m_hRRdphiPhi{}, m_hRDzZ{}, m_hRRdphiQopt{};
  // per layer
  std::array<TH1F*, 7> m_hLRdphi{}, m_hLDz{};
  std::array<TH2F*, 7> m_hLRdphiPhi{}, m_hLDzZ{}, m_hLRdphiQopt{}, m_hLDzPhi{};

  // counters
  unsigned long m_nEvents{0}, m_nTracks{0}, m_nUsed{0}, m_nNoSi{0}, m_nNoTpc{0}, m_nSel{0};

  // ---- tree
  unsigned int t_event{0}, t_trackId{0}, t_siIndex{0}, t_tpcIndex{0};
  int t_siChain{-1}, t_siStatus{0}, t_half{0};
  unsigned int t_nMvtx{0}, t_nIntt{0}, t_nTpcClus{0};
  float t_matchChi2{0}, t_matchDphi{0}, t_matchDeta{0}, t_matchDz0{0}, t_fittedZOffset{0};
  float t_siPhi{0}, t_siEta{0}, t_siTanl{0}, t_siZ0{0}, t_siDca{0}, t_siPt{0}, t_siQopt{0}, t_siR{0};
  int t_siCharge{0};
  float t_tpcPhi{0}, t_tpcEta{0}, t_tpcTanl{0}, t_tpcZ0{0}, t_tpcDca{0}, t_tpcPt{0}, t_tpcQopt{0}, t_tpcR{0};
  int t_tpcCharge{0};
  float t_dphi{0}, t_deta{0}, t_dtanl{0}, t_dz0{0}, t_ddca{0}, t_dqopt{0}, t_dptrel{0};
  int t_nR{0};
  float t_R[kMaxRadii]{}, t_Rphi[kMaxRadii]{}, t_Rz[kMaxRadii]{}, t_Rrdphi[kMaxRadii]{}, t_Rdz[kMaxRadii]{},
      t_Rddir[kMaxRadii]{};
  std::vector<int> t_pLayer;
  std::vector<float> t_pR, t_pPhi, t_pZ, t_pRdphi, t_pDz;
};

#endif
