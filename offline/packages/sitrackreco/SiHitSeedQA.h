#ifndef SIHITSEEDQA_H
#define SIHITSEEDQA_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;
class TH2D;
class TDirectory;
class SiHitSeedEvent;
class ActsGeometry;
class PHG4CylinderGeomContainer;
struct SiHitPoint;
struct SiHitChain;

// QA of SiHitSeedReco + per-event vertex estimates from the accepted chain links.
//
// Vertex z  : profile of dz = z_inner - z_outer vs outer z, binned; straight-line fit of the
//             bin means; vertex z = z where the fitted dz = 0.  For a straight track from
//             (0,0,zv): dz = (z_outer - zv)(R_in/R_out - 1), so every layer pair crosses dz = 0
//             at z_outer = zv.  The fit starts from links whose own extrapolation to R = 0 is
//             within setVertexSeedWindow of the median (vtx*_z_median), then rejects outliers.
// Vertex xy : straight track from (x0,y0): dphi = phi_in - phi_out = k (y0 cos phi - x0 sin phi),
//             k = 1/R_in - 1/R_out.  Fit  dphi/k = c + A cos(phi_out) + B sin(phi_out)
//             -> y0 = A, x0 = -B   (in the intrinsic Uphi frame, then rotated, see below).
//
// OUTPUT FRAME: millimetres, nominal sPHENIX frame centred on the detector:
//   z : nominal longitudinal coordinate (MVTX chip / INTT sensor centres, midplane at 0)
//   x,y: Uphi frame rotated by setUphiToGlobalRotation (default +phi0 of MVTX layer-2 stave 0,
//        0.1481 rad, i.e. global phi = Uphi + 0.1481 if the stave phases are global azimuths)
//   plus an optional offset setDetectorCenterMm (survey/alignment), default (0,0,0).
// No alignment or survey is applied; nominal radii and the approximate localX signs enter.
//
// Standard vertices (SvtxVertexMap "SiliconSvtxVertexMap", PHSimpleVertexFinder):
//   svx_x/y/z, errors, chi2/ndof, ntracks : untouched values in the GLOBAL (aligned) frame [cm]
//   svx_det_x/y/z                          : the same vertices in the DETECTOR frame [mm]
// Detector frame <-> global frame: rigid transform global = R * det + T fitted (Kabsch) between
// the detector-frame position of chain hits (x = R_layer cos(Uphi + rotation), y = R_layer
// sin(...), z = nominal Vz; same frame as vtx_*) and the global position of the same hits from
// the aligned Acts geometry.  The fit is accumulated over all events so far (detector placement
// is constant) and stored as det_T_x/y/z [mm] = position of the nominal detector centre in the
// global frame, det_rot_x/y/z [mrad] = small rotation angles, det_fit_rms [mm].
// The beam position in the detector frame is then the vertex x/y in the detector frame
// (vtx_* from our links, or svx_det_* from the standard vertex), averaged over events.
//
// Three vertices per event: all links ("vtx_*"), and each clamshell half separately
// ("vtxA_*", "vtxB_*"), split by the outer hit's Uphi:
//   half A: Uphi in [B, B+pi),  half B: Uphi in [B+pi, B+2pi),  B = setClamshellBoundary.
class SiHitSeedQA : public SubsysReco
{
 public:
  explicit SiHitSeedQA(
      const std::string& out = "si_hit_seed_qa.root",
      const std::string& name = "SiHitSeedQA");

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setInputNodeName(const std::string& s) { m_inputNodeName = s; }

  // Which accepted chain steps enter the vertex fits (bit mask):
  // 1 = seed links (default), 2 = MVTX propagation steps, 4 = INTT propagation steps.
  void setVertexStepMask(int mask) { m_vtxStepMask = mask; }
  void setVertexZProfile(int nbins, double zloCm, double zhiCm) { m_vtxNz = nbins; m_vtxZlo = zloCm; m_vtxZhi = zhiCm; }
  void setVertexMinEntriesPerBin(int n) { m_vtxMinPerBin = n; }
  void setVertexClipSigma(double s) { m_vtxClip = s; }
  void setVertexIterations(int n) { m_vtxIter = n; }
  // Links enter the dz-vs-z fit only if their own extrapolation to R=0 is within this
  // distance [cm] of the median of all links (0 = use all links).
  void setVertexSeedWindow(double cm) { m_vtxSeedWindow = cm; }
  // Include the constant term c in the dphi/k fit (default true).
  void setVertexXYFitConstant(bool v) { m_vtxXYConst = v; }
  // Clamshell boundary in Uphi [rad]; the second boundary is at +pi (same as SiHitSeedDisplay).
  void setClamshellBoundary(double uphi) { m_boundary = uphi; }
  // Output frame: global phi = Uphi + rotation; then + detector-centre offset [mm].
  void setUphiToGlobalRotation(double rad) { m_rotation = rad; }
  void setDetectorCenterMm(double xMm, double yMm, double zMm) { m_center = {xMm, yMm, zMm}; }
  // Save the per-event dz-vs-z and dphi/k-vs-phi profiles (with fitted curves) for the first N events.
  void setSaveEventProfiles(unsigned int n) { m_nSaveProfiles = n; }
  // Per-event 2D histograms saved next to the profiles (same events, same x binning):
  //   h_evtN_<label>_dz_vs_z      : dz = z_in - z_out [cm] vs outer z [cm]
  //   h_evtN_<label>_dphik_vs_phi : dphi/k [cm] vs outer Uphi [rad]
  // y ranges [cm] and number of y bins:
  void setEventHist2DdzRange(int nbins, double loCm, double hiCm) { m_h2dzN = nbins; m_h2dzLo = loCm; m_h2dzHi = hiCm; }
  void setEventHist2DdphikRange(int nbins, double loCm, double hiCm) { m_h2xyN = nbins; m_h2xyLo = loCm; m_h2xyHi = hiCm; }
  void setStandardVertexMapName(const std::string& s) { m_vertexMapName = s; }
  // Max chain hits per event used in the detector <-> global transform fit.
  void setFrameFitMaxPoints(int n) { m_frameMaxPoints = n; }

  // ------------------------------------------------------------------------
  // Vertex-sample cleanup.  These cuts act only in SiHitSeedQA; the pattern
  // recognition itself stays loose, so the cut ladder can be scanned without
  // re-running SiHitSeedReco with different reconstruction thresholds.
  // ------------------------------------------------------------------------
  void setVertexMinMvtxLayers(int n) { m_vtxMinMvtxLayers = n; }
  void setVertexMinInttLayers(int n) { m_vtxMinInttLayers = n; }

  // Straight-line eta from the innermost/outermost MVTX points.
  // <= 0 disables the cut.
  void setVertexMaxAbsEta(double v) { m_vtxMaxAbsEta = v; }

  // Chain score / propagation residual quality. <= 0 disables the cut.
  void setVertexMaxChainScore(double v) { m_vtxMaxChainScore = v; }
  void setVertexMaxResidualZBins(double v) { m_vtxMaxResidualZBins = v; }
  void setVertexMaxResidualPhiBins(double v) { m_vtxMaxResidualPhiBins = v; }

  // Restrict the link estimator to a particular ordered layer pair.
  // (-1,-1) disables this requirement.
  void setVertexSeedLayerPair(int fromLayer, int toLayer)
  {
    m_vtxSeedFromLayer = fromLayer;
    m_vtxSeedToLayer = toLayer;
  }

  // Event-level azimuthal coverage after chain cuts.
  // Example: (12,6) means at least 6 occupied 30-degree sectors.
  void setVertexPhiCoverage(int nsectors, int minOccupied)
  {
    m_vtxPhiSectors = nsectors;
    m_vtxMinPhiSectors = minOccupied;
  }

  // Event-level half-A / half-B consistency [mm]. <=0 disables each cut.
  void setVertexHalfAgreement(double dxMm, double dyMm)
  {
    m_vtxMaxHalfDx = dxMm;
    m_vtxMaxHalfDy = dyMm;
  }

  // Stage 7: build one straight 3-MVTX tracklet per accepted chain and fit
  // their common 3D point.  Both link and tracklet vertices are always saved;
  // this switch chooses which one drives final_vtx_* and vertex_selected.
  void setUseTrackletVertex(bool v) { m_useTrackletVertex = v; }

 private:

  struct Tracklet
  {
    std::array<double, 3> point{{0, 0, 0}};  // cm, intrinsic Uphi frame
    std::array<double, 3> dir{{0, 0, 1}};    // unit vector
    double eta = 0;
    double phi = 0;
    double rms = 0;                          // cm, 3-hit line residual
  };

  struct Link
  {
    int from_layer, to_layer;
    double z_from, dz, phi_from, dphi;  // cm, cm, rad (Uphi), rad
  };
  // One vertex estimate.  Lengths in mm, output frame (except the *_cm fit internals).
  struct Vertex
  {
    double z, z_err, z_median;
    double dz_intercept, dz_slope, dz_rms;  // dz[mm] = intercept[mm] + slope * z_outer[mm]
    int z_nlinks, z_nbins;
    double x, y, x_err, y_err, xy_corr;
    double xy_const, xy_rms;  // [mm]
    int xy_nlinks;
    int ok;  // bit 1: z fit ok, bit 2: xy fit ok
    int nlinks;
    // fit internals in the Uphi frame [cm], used for the saved profiles
    double a_cm, b_cm, c_cm, A_cm, B_cm;
  };

  double chainEta(const SiHitChain& c) const;
  double chainMaxResidualZ(const SiHitChain& c) const;
  double chainMaxResidualPhi(const SiHitChain& c) const;
  bool chainPassesVertexCuts(const SiHitChain& c) const;

  bool makeMvtxTracklet(const SiHitChain& c, Tracklet& t) const;
  void fitTrackletVertex(const std::vector<Tracklet>& tracklets, Vertex& v) const;

  int countPhiSectors(const std::vector<double>& phis) const;
  void resetVertex(Vertex& v) const;
  void fitVertex(const std::vector<Link>& links, Vertex& v) const;
  void fitVertexZ(const std::vector<Link>& links, double zSeedCm, Vertex& v) const;
  void fitVertexXY(const std::vector<Link>& links, Vertex& v) const;
  void branchVertex(TTree* t, const std::string& prefix, Vertex& v);
  void saveEventProfiles(unsigned long long evt, const char* label, const std::vector<Link>& links, const Vertex& v);
  bool hitToGlobal(const SiHitPoint& h, double& gx, double& gy, double& gz) const;
  void fitFrame(const SiHitSeedEvent& d);  // updates m_R, m_T and the det_* branches

  std::string m_vertexMapName = "SiliconSvtxVertexMap";
  int m_frameMaxPoints = 5000;
  ActsGeometry* m_geom = nullptr;
  PHG4CylinderGeomContainer* m_mvtxGeom = nullptr;
  PHG4CylinderGeomContainer* m_inttGeom = nullptr;
  // accumulated sums for the Kabsch fit [cm]
  double m_fN = 0;
  std::array<double, 3> m_fSd{{0, 0, 0}}, m_fSg{{0, 0, 0}};
  std::array<double, 9> m_fSdg{{0, 0, 0, 0, 0, 0, 0, 0, 0}};  // sum d_i g_j
  double m_fSdd = 0, m_fSgg = 0;
  std::array<double, 9> m_R{{1, 0, 0, 0, 1, 0, 0, 0, 1}};  // global = R det + T  [cm]
  std::array<double, 3> m_T{{0, 0, 0}};
  bool m_frameOk = false;

  // standard vertex branches
  std::vector<double> svx_x, svx_y, svx_z, svx_ex, svx_ey, svx_ez, svx_chi2, svx_ndof;
  std::vector<int> svx_ntracks;
  std::vector<double> svx_det_x, svx_det_y, svx_det_z;
  // frame branches
  double det_T_x = 0, det_T_y = 0, det_T_z = 0, det_rot_x = 0, det_rot_y = 0, det_rot_z = 0;
  double det_fit_rms = 0, evt_fit_rms = 0;
  int det_fit_npoints = 0, evt_fit_npoints = 0;

  std::string m_inputNodeName = "SI_HIT_SEED_EVENT";
  std::string m_out;

  TFile* m_file = nullptr;
  TDirectory* m_profDir = nullptr;
  TTree* m_evt = nullptr;
  TTree* m_hits = nullptr;
  TTree* m_clus = nullptr;
  TTree* m_steps = nullptr;
  TTree* m_chains = nullptr;
  TH2D* m_dz = nullptr;
  TH2D* m_dzcm = nullptr;
  TH2D* m_dphi = nullptr;

  // vertex fit settings
  int m_vtxStepMask = 1;
  int m_vtxNz = 40;
  double m_vtxZlo = -25.0;
  double m_vtxZhi = 25.0;
  int m_vtxMinPerBin = 2;
  double m_vtxClip = 2.5;
  int m_vtxIter = 3;
  double m_vtxSeedWindow = 3.0;
  bool m_vtxXYConst = true;
  double m_boundary = 0.0;
  double m_rotation = 0.1481;
  std::array<double, 3> m_center{{0.0, 0.0, 0.0}};
  unsigned int m_nSaveProfiles = 5;
  unsigned int m_nSavedProfiles = 0;
  int m_h2dzN = 160;
  double m_h2dzLo = -8.0, m_h2dzHi = 8.0;
  int m_h2xyN = 100;
  double m_h2xyLo = -2.0, m_h2xyHi = 2.0;

  // vertex sample cleanup
  int m_vtxMinMvtxLayers = 0;
  int m_vtxMinInttLayers = 0;
  double m_vtxMaxAbsEta = -1.0;
  double m_vtxMaxChainScore = -1.0;
  double m_vtxMaxResidualZBins = -1.0;
  double m_vtxMaxResidualPhiBins = -1.0;
  int m_vtxSeedFromLayer = -1;
  int m_vtxSeedToLayer = -1;

  int m_vtxPhiSectors = 12;
  int m_vtxMinPhiSectors = 0;
  double m_vtxMaxHalfDx = -1.0;
  double m_vtxMaxHalfDy = -1.0;
  bool m_useTrackletVertex = false;

  // nominal layer radii [cm], same as SiHitSeedReco
  const std::array<double, 7> m_radius{{2.523, 3.336, 4.148, 7.188 - 0.0036, 7.732 - 0.0036, 9.680 - 0.0036, 10.262 - 0.0036}};

  unsigned long long event = 0;
  unsigned long long hitsetkey = 0;
  unsigned long long hitkey = 0;
  unsigned long long cluskey = 0;

  int layer = 0;
  int row = 0;
  int col = 0;
  int stave = 0;
  int chip = 0;
  int ladderphi = 0;
  int ladderz = 0;
  int chain_id = 0;
  int kind = 0;
  int from_layer = 0;
  int to_layer = 0;

  float adc = 0;

  double lx = 0;
  double ly = 0;
  double z = 0;
  double phi = 0;
  double ref_z = 0;
  double ref_phi = 0;
  double pred_z = 0;
  double pred_phi = 0;
  double delta_z = 0;
  double delta_phi = 0;
  double res_z = 0;
  double res_phi = 0;
  double score = 0;
  double vertex_z = 0;
  double vertex_z_linefit = 0;
  int vtx_tr_npairs = 0;
  int vtx_tr_npeak = 0;
  double fit_a = 0;
  double fit_b = 0;

  int nhits = 0;
  int nclusters = 0;
  int nchains = 0;
  int nvertex = 0;

  // event-level vertices: all links, clamshell half A, clamshell half B
  Vertex m_vtxAll{}, m_vtxA{}, m_vtxB{};

  // cleanup / final-estimator QA
  int n_vertex_chains = 0;
  int n_vertex_tracklets = 0;
  int vertex_phi_sectors = 0;
  int vertex_phi_coverage_ok = 0;
  int vertex_half_agreement_ok = 0;
  int vertex_selected = 0;
  double vertex_half_dx = 0;
  double vertex_half_dy = 0;

  double final_vtx_x = 0;
  double final_vtx_y = 0;
  double final_vtx_z = 0;
  double final_vtx_x_err = 0;
  double final_vtx_y_err = 0;
  double final_vtx_z_err = 0;
  int final_vtx_ok = 0;

  // chain-tree QA
  int chain_nmvtx = 0;
  int chain_nintt = 0;
  int chain_pass_vertex_cuts = 0;
  double chain_eta = 0;
  double chain_score = 0;
  double chain_max_res_z = 0;
  double chain_max_res_phi = 0;

  // tracklet vertices: all / half A / half B
  Vertex m_trkVtxAll{}, m_trkVtxA{}, m_trkVtxB{};

  std::vector<int> vlayer;
  std::vector<double> vz;
  std::vector<double> vphi;
};

#endif