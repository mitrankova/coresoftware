#ifndef SITPC_TRACKRESIDUALS_H
#define SITPC_TRACKRESIDUALS_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class SiHitSeedEvent;
struct SiHitChain;
class SiTpc_TrackContainer;
class SvtxVertexMap;
class TpcCrossingDecisionContainer;
class Tpc_PolyClusterContainer;
class TFile;
class TH3F;
class TTree;

// Residual QA of the matched Si + TPC tracks (SiTpcTrackMatcher), in the style of
// Tpc_PolyClusterResiduals.
//
// Tree "residuals": one entry per Si+TPC track
//   track / association / match variables (dphi, deta, dz0, chi2, Si-only and TPC-only phi,
//   eta, z0, TPC-only pt), the combined fit (pt, eta, phi, charge, dca, z0, tanl, R,
//   tpc_z_offset, rms), and per fit point (Si layer points and TPC clusters):
//     point_type (0 = Si, 1 = TPC), layer, source (Si point index / TPC cluster index),
//     x, y, z, r, phi of the point; state_* = fitted helix at the point's radius;
//     delta_phi = phi_point - phi_state, residual_rphi = r * delta_phi,
//     residual_z = z_point - [tpc_z_offset for TPC points] - z_state;
//     cluster_adc, cluster_pad_size (TPC points, from TPC_POLYCLUSTERS if present).
//   Frame: TPC frame of the matcher (Si halves on the TPC beam axis) (TPC z as measured, offset in residual_z).
//   Bunch crossing of the track (both from the same track):
//     tpc_crossing          crossing used for the TPC (TPC_CROSSING_DECISIONS, by assembled id;
//                           -9999 without decision)
//     tpc_crossing_tier     tier of that decision (0, 1, 2; -1 without decision)
//     intt_crossing         most frequent INTT time bucket of the INTT hits of the track's Si chain
//                           (SI_HIT_SEED_EVENT; -9999 if the chain has no INTT hit)
//     intt_n_hits           number of INTT hits of the chain
//     intt_crossing_spread  max - min INTT time bucket of those hits
//     crossing_diff         intt_crossing - tpc_crossing (-9999 if either is missing)
//     n_tracks_crossing     Si+TPC tracks of this event (all valid ones, before the tree cuts)
//                           with the same tpc_crossing, this track included (0 without crossing)
//     n_si_chains_crossing  selected Si chains of SI_HIT_SEED_EVENT with the same crossing
//                           (TpcCrossingResolver: Si-matched TPC tracks of that collision)
//   Collision vertex of the track (TpcCrossingResolver decision, Si frame, NaN without one):
//     track_vertex_z        BEFORE clustering: this track's own z at the beam from its Si chain,
//                           = decision tpc_z0 + delta_z (TPC z0 at the crossing + measured Si - TPC
//                           offset); for TPC-only matches (no own chain) it is the vertex z itself
//     vertex_z              AFTER clustering: decided collision z (silicon_vertex_z) = mean of the
//                           chain z of the tracks of that crossing's vertex (lone track: its own z)
//     vertex_dz             delta_z of the decision (Si - TPC z at the chosen crossing)
//     vertex_score          selected_score of the decision (NaN for TPC-only matches)
//     vertex_status         TpcCrossingStatus of the decision (10 = AmbiguousWithVertex: lone track
//                           with an MVTX-only chain, kept without a collision: vertex_* NaN,
//                           track_vertex_z = its own Si z)
//     vertex_own_z          1 if track_vertex_z comes from the track's own Si chain (tier 0, lone),
//                           0 for TPC-only matches (z taken from the vertex)
//     vertex_n_si           tracks of this event at this crossing with their own Si z
//     vertex_z_mean         unweighted mean of their track_vertex_z
//     vertex_z_rms          sample standard deviation of those z (NaN if fewer than 2)
//     vertex_z_err          setVertexTrackZSigma / sqrt(vertex_n_si) (expected uncertainty)
//     si_z0 - vertex_z      (plot it) Si-only fit z0 versus the decided vertex
//   Pair-finder vertices of the same crossing (SiTpc_VertexFinder; NaN / 0 if absent or the node is
//   missing), for the z-alignment / DCA_z study (plot_z_alignment.C):
//     pv_sitpc_z, pv_sitpc_ez, pv_sitpc_n   SiTpcPolyVertexMap vertex of this crossing nearest in z
//                                           to the track z0: z, sigma_z, number of tracks
//     pv_sitpc_in                           1 if this track is one of its tracks
//     pv_sitpc_z_excl                       vertex z WITHOUT this track (unweighted leave-one-out
//                                           (n z - z0) / (n - 1); = z if the track is not in it)
//     pv_si_z, pv_si_n                      SiPolyVertexMap vertex of this crossing nearest to si_z0
//   -> prompt DCA_z = z0 - pv_sitpc_z_excl (Si+TPC track), si_z0 - pv_si_z (Si-only track)
// Tree "vertices": one entry per (event, decided crossing with a vertex), i.e. per collision
//   event, crossing, vertex_z, n_tracks (all Si+TPC tracks), n_si (own Si z), n_tier0,
//   z_mean, z_rms, z_err (as above), and per track of the collision (vectors):
//   track_z (track_vertex_z), track_own_z, track_tier, track_eta, track_pt, track_phi, track_charge
//   -> multiplicity per collision: one entry per collision (not weighted by n as n_tracks_crossing);
//      dN/deta per collision = sum of track_eta / number of entries.
//   Per track also:
//     n_mvtx, n_intt, n_tpc          number of fit points per detector
//     x0, y0, z0                     point of closest approach to the beam (TPC frame)
//     res_rphi_layer[55], res_z_layer[55]
//                                    residual of the fit at each layer, index = layer
//                                    (0-2 MVTX, 3-6 INTT, 7-54 TPC); NaN if the track has no
//                                    point on that layer, mean if it has several
//     dphi_layer[55]                 delta_phi (rad) of the fit at each layer, same convention
//     dphi_prev_layer[55]            delta_phi of the previous (next outer) layer with a point
//     prev_layer[55]                 that previous layer index (-1 if none)
//     ddphi_layer[55]                delta-delta-phi = dphi_layer - dphi_prev_layer (NaN if no previous)
//     rms_rphi_{mvtx,intt,tpc}, rms_z_{mvtx,intt,tpc}
//                                    rms of the residuals per detector (INTT z: strip length)
// Tree "matching": every candidate Si-TPC pair stored by the matcher
//   (si_index, tpc_index, dphi, deta, dz0, chi2, accepted) - to tune the match windows.
//
// Delta-delta-phi 3D histograms (same chain logic as TpcSiliconCrossingMatcher: the chain is
// followed from the outside in, and each layer's delta_phi is compared with the delta_phi of the
// previous, i.e. next outer, layer that has a point). Separate histograms per charge
// (index 0 = negative, 1 = positive) since the phi residual sign follows the bending direction:
//   h3_dphi_vs_prev_l{L}_{neg,pos}   x = dphi(prev layer), y = dphi(L), z = pT
//                                    one per Si layer L = 0..6; profiling y vs x in pT slices
//                                    gives the offset/slope of TpcSiliconCrossingMatcher::setDynamicPhiMean
//   h3_ddphi_vs_layer_{neg,pos}      x = ddphi, y = layer (0-54), z = pT
//   h3_ddphi_vs_dphi_prev_{neg,pos}  x = dphi(prev layer), y = ddphi, z = layer (Si layers 0-6)
// The pT axis uses the combined fit pT; straight (non-curved) fits are not filled.
// NOTE: delta_phi above is the residual to the combined Si+TPC fit, which contains the point.
//
// Delta-delta-phi in silicon HARDWARE coordinates (SiHitSeedReco chain of the track, read from
// SI_HIT_SEED_EVENT via si_chain_id; skipped if that node is absent). Uses the aligned intrinsic
// Uphi / Vz of the representative hit per layer, chain order outermost layer first, nominal
// layer radii of SiHitSeedReco. No fit is involved, so these are the quantities to tune the
// chaining windows with. Axis pT = TPC-only pT, charge = TPC charge (what the chaining sees).
//   step into layer L from the previous (outer) chain layer P:
//     hw_dphi_step[L]  = Uphi(L) - Uphi(P)                                [rad]
//     hw_dphi_pred[L]  = hw_dphi_step[P] * (r_L - r_P) / (r_P - r_PP)      [rad]  (constant dphi/dr)
//     hw_ddphi[L]      = hw_dphi_step[L] - hw_dphi_pred[L]                 [rad]
//     hw_dz_step[L]    = Vz(L) - Vz(P)                                     [cm]
//   chain step residuals (exactly what SiHitSeedReco cuts on):
//     hw_res_phi[L] = Uphi(L) - pred_phi, hw_res_z[L] = Vz(L) - pred_z     [rad], [cm]
//     hw_step_kind[L] = SiChainStep::Kind (0 seed, 1 MVTX propagation, 2 INTT propagation)
//   histograms (per charge, _neg/_pos):
//     h3_hw_dphi_step_vs_layer   x = hw_dphi_step, y = L (0-6), z = pT  -> pT-dependent step window centre/width
//     h3_hw_dphi_vs_pred_l{L}    x = hw_dphi_pred, y = hw_dphi_step, z = pT, L = 0..6
//     h3_hw_ddphi_vs_layer       x = hw_ddphi, y = L, z = pT
//     h3_hw_res_phi_vs_layer     x = hw_res_phi, y = L, z = pT
//     h3_hw_res_z_vs_layer       x = hw_res_z, y = L, z = pT
class SiTpc_TrackResiduals : public SubsysReco
{
 public:
  explicit SiTpc_TrackResiduals(const std::string& name = "SiTpc_TrackResiduals",
                                const std::string& outfilename = "sitpc_track_residuals.root");
  ~SiTpc_TrackResiduals() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setTrackNodeName(const std::string& n) { m_trackNodeName = n; }
  void setTpcClusterNodeName(const std::string& n) { m_clusterNodeName = n; }
  void setSiSeedNodeName(const std::string& n) { m_siSeedNodeName = n; }
  void setCrossingDecisionNodeName(const std::string& n) { m_crossingNodeName = n; }
  void setOutputFileName(const std::string& n) { m_outfilename = n; }
  void setMinPt(double v) { m_minPt = v; }
  void setMaxPt(double v) { m_maxPt = v; }
  void setMinTpcClusters(unsigned int v) { m_minTpcClusters = v; }
  void setMinSiPoints(unsigned int v) { m_minSiPoints = v; }
  void setWriteMatchingTree(bool v) { m_writeMatching = v; }
  // pair-finder vertex maps (SiTpc_VertexFinder) for the pv_* branches, "" = not used
  void setSiTpcVertexMapName(const std::string& n) { m_sitpcVertexMapName = n; }
  void setSiVertexMapName(const std::string& n) { m_siVertexMapName = n; }
  // collision vertices: "vertices" tree, and the single-track Si z resolution at the beam [cm]
  // used for vertex_z_err (measured ~0.22 cm from the spread of multi-track vertices)
  void setWriteVertexTree(bool v) { m_writeVertices = v; }
  void setVertexTrackZSigma(double v) { m_vertexTrackZSigma = v; }

  // delta-delta-phi 3D histograms
  void setWriteDdphiHists(bool v) { m_writeDdphiHists = v; }
  void setDdphiDphiAxis(int nbins, double maxAbs)
  {
    m_dphiBins = nbins;
    m_dphiMax = maxAbs;
  }
  void setDdphiAxis(int nbins, double maxAbs)
  {
    m_ddphiBins = nbins;
    m_ddphiMax = maxAbs;
  }
  void setDdphiPtAxis(int nbins, double ptMin, double ptMax)
  {
    m_ptBins = nbins;
    m_ptMin = ptMin;
    m_ptMax = ptMax;
  }

  // delta-delta-phi 3D histograms in silicon hardware coordinates
  void setWriteHardwareDdphiHists(bool v) { m_writeHwHists = v; }
  void setHardwareDphiStepAxis(int nbins, double maxAbs)
  {
    m_hwStepBins = nbins;
    m_hwStepMax = maxAbs;
  }
  void setHardwareDdphiAxis(int nbins, double maxAbs)
  {
    m_hwDdphiBins = nbins;
    m_hwDdphiMax = maxAbs;
  }
  void setHardwareResidualAxes(int nbins, double phiMaxAbs, double zMaxAbsCm)
  {
    m_hwResBins = nbins;
    m_hwResPhiMax = phiMaxAbs;
    m_hwResZMax = zMaxAbsCm;
  }

 private:
  void reset_tree_values();
  void book_ddphi_hists();
  void fill_ddphi_hists();
  void book_hw_hists();
  void fill_hw(const SiHitChain& chain, double tpcPt, int tpcCharge);

  std::string m_outfilename;
  std::string m_trackNodeName{"SITPC_TRACKS"};
  std::string m_clusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_siSeedNodeName{"SI_HIT_SEED_EVENT"};
  std::string m_crossingNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_sitpcVertexMapName{"SiTpcPolyVertexMap"};
  std::string m_siVertexMapName{"SiPolyVertexMap"};
  SvtxVertexMap* m_sitpcVertices{nullptr};
  SvtxVertexMap* m_siVertices{nullptr};
  // pair-finder vertices of the track's crossing
  double m_pvSiTpcZ{0}, m_pvSiTpcEz{0}, m_pvSiTpcZExcl{0}, m_pvSiZ{0};
  double m_pvSiTpcX{0}, m_pvSiTpcY{0}, m_pvSiTpcEx{0}, m_pvSiTpcEy{0}, m_pvSiTpcDcaXY{0};
  int m_pvSiTpcN{0}, m_pvSiTpcIn{0}, m_pvSiN{0};

  double m_minPt{0.0};
  double m_maxPt{1.0e30};
  unsigned int m_minTpcClusters{0};
  unsigned int m_minSiPoints{0};
  bool m_writeMatching{true};
  bool m_writeVertices{true};
  double m_vertexTrackZSigma{0.22};

  unsigned int m_evt{0};
  TFile* m_outfile{nullptr};
  TTree* m_tree{nullptr};
  TTree* m_matchTree{nullptr};
  TTree* m_vertexTree{nullptr};
  SiTpc_TrackContainer* m_tracks{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};
  SiHitSeedEvent* m_siSeeds{nullptr};
  TpcCrossingDecisionContainer* m_decisions{nullptr};

  // ---- residuals tree
  unsigned int m_event{0};
  unsigned int m_trackId{0};
  int m_siChainId{-1};
  unsigned int m_siTrajIndex{0};
  unsigned int m_tpcTrackId{0};
  unsigned int m_tpcAssembledId{0};
  int m_fitStatus{0};
  unsigned int m_nSi{0};
  unsigned int m_nTpc{0};
  double m_matchDphi{0}, m_matchDeta{0}, m_matchDz0{0}, m_matchChi2{0};
  double m_siPhi{0}, m_siEta{0}, m_siZ0{0};
  double m_tpcPhi{0}, m_tpcEta{0}, m_tpcZ0{0}, m_tpcPt{0};
  int m_tpcCharge{0};
  double m_pt{0}, m_eta{0}, m_theta{0}, m_phi{0}, m_dca{0}, m_z0{0}, m_tanl{0};
  int m_charge{0};
  double m_R{0}, m_tpcZOffset{0}, m_circleRms{0}, m_zRms{0}, m_pcaX{0}, m_pcaY{0};
  unsigned int m_nMvtx{0}, m_nIntt{0}, m_nTpcLayers{0};
  static constexpr int kNoCrossing = -9999;
  int m_tpcCrossing{kNoCrossing};
  int m_tpcCrossingTier{-1};
  int m_inttCrossing{kNoCrossing};
  int m_inttNHits{0};
  int m_inttCrossingSpread{0};
  int m_crossingDiff{kNoCrossing};
  int m_nTracksCrossing{0};
  int m_nSiChainsCrossing{0};
  std::map<int, int> m_tracksPerCrossing;    // this event: tpc_crossing -> Si+TPC tracks
  std::map<int, int> m_siChainsPerCrossing;  // this event: tpc_crossing -> selected Si chains
  // collision vertex of the track
  double m_trackVertexZ{0}, m_vertexZ{0}, m_vertexDz{0}, m_vertexScore{0};
  int m_vertexStatus{0}, m_vertexOwnZ{0}, m_vertexNSi{0};
  double m_vertexZMean{0}, m_vertexZRms{0}, m_vertexZErr{0};
  // this event: crossing -> collision vertex summary (from all valid Si+TPC tracks)
  struct VertexSummary
  {
    double z{0};  // decided vertex z
    int nTracks{0}, nSi{0}, nTier0{0};
    double zMean{0}, zRms{0}, zErr{0};
    std::vector<float> trackZ, eta, pt, phi;
    std::vector<int> ownZ, tier, charge;
  };
  std::map<int, VertexSummary> m_vertices;
  // ---- vertices tree
  int m_vCrossing{0}, m_vNTracks{0}, m_vNSi{0}, m_vNTier0{0};
  double m_vZ{0}, m_vZMean{0}, m_vZRms{0}, m_vZErr{0};
  std::vector<float> m_vTrackZ, m_vEta, m_vPt, m_vPhi;
  std::vector<int> m_vOwnZ, m_vTier, m_vCharge;
  static constexpr int kNLayers = 55;    // 0-2 MVTX, 3-6 INTT, 7-54 TPC
  static constexpr int kNSiLayers = 7;   // 0-2 MVTX, 3-6 INTT
  double m_resRPhiLayer[kNLayers]{};
  double m_resZLayer[kNLayers]{};
  double m_dphiLayer[kNLayers]{};
  double m_dphiPrevLayer[kNLayers]{};
  int m_prevLayer[kNLayers]{};
  double m_ddphiLayer[kNLayers]{};
  // hardware chain, index = Si layer
  int m_hwHasChain{0};
  int m_hwPrevLayer[kNSiLayers]{};
  int m_hwStepKind[kNSiLayers]{};
  double m_hwUphi[kNSiLayers]{};
  double m_hwVz[kNSiLayers]{};
  double m_hwDphiStep[kNSiLayers]{};
  double m_hwDphiPred[kNSiLayers]{};
  double m_hwDdphi[kNSiLayers]{};
  double m_hwDzStep[kNSiLayers]{};
  double m_hwResPhi[kNSiLayers]{};
  double m_hwResZ[kNSiLayers]{};
  double m_rmsRPhi[3]{};  // MVTX, INTT, TPC
  double m_rmsZ[3]{};
  std::vector<int> m_pointType;
  std::vector<int> m_layer;
  std::vector<unsigned int> m_source;
  std::vector<double> m_x, m_y, m_z, m_r, m_pointPhi;
  std::vector<double> m_stateX, m_stateY, m_stateZ, m_stateR, m_statePhi;
  std::vector<double> m_deltaPhi, m_residualRPhi, m_residualZ;
  std::vector<double> m_clusterAdc;
  std::vector<unsigned int> m_clusterPadSize;

  // ---- matching tree
  unsigned int m_cSi{0}, m_cTpc{0};
  float m_cDphi{0}, m_cDeta{0}, m_cDz0{0}, m_cChi2{0};
  int m_cAccepted{0};

  // ---- delta-delta-phi 3D histograms, [charge] with 0 = negative, 1 = positive
  bool m_writeDdphiHists{true};
  int m_dphiBins{100};
  double m_dphiMax{0.02};   // rad
  int m_ddphiBins{100};
  double m_ddphiMax{0.01};  // rad
  int m_ptBins{20};
  double m_ptMin{0.0};
  double m_ptMax{5.0};      // GeV
  TH3F* m_h3DphiVsPrev[kNSiLayers][2]{};
  TH3F* m_h3DdphiVsLayer[2]{};
  TH3F* m_h3DdphiVsDphiPrev[2]{};

  // ---- hardware-coordinate delta-delta-phi 3D histograms, [charge]
  bool m_writeHwHists{true};
  int m_hwStepBins{100};
  double m_hwStepMax{0.1};      // rad
  int m_hwDdphiBins{100};
  double m_hwDdphiMax{0.02};    // rad
  int m_hwResBins{100};
  double m_hwResPhiMax{0.02};   // rad
  double m_hwResZMax{0.5};      // cm
  TH3F* m_h3HwDphiStepVsLayer[2]{};
  TH3F* m_h3HwDphiVsPred[kNSiLayers][2]{};
  TH3F* m_h3HwDdphiVsLayer[2]{};
  TH3F* m_h3HwResPhiVsLayer[2]{};
  TH3F* m_h3HwResZVsLayer[2]{};
};

#endif