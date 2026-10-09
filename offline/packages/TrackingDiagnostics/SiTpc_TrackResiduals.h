#ifndef SITPC_TRACKRESIDUALS_H
#define SITPC_TRACKRESIDUALS_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class SiTpc_TrackContainer;
class Tpc_PolyClusterContainer;
class TFile;
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
//   Frame: beam-axis frame of the matcher (TPC z as measured, offset in residual_z).
//   Per track also:
//     n_mvtx, n_intt, n_tpc          number of fit points per detector
//     x0, y0, z0                     point of closest approach to the beam (beam-axis frame)
//     res_rphi_layer[55], res_z_layer[55]
//                                    residual of the fit at each layer, index = layer
//                                    (0-2 MVTX, 3-6 INTT, 7-54 TPC); NaN if the track has no
//                                    point on that layer, mean if it has several
//     rms_rphi_{mvtx,intt,tpc}, rms_z_{mvtx,intt,tpc}
//                                    rms of the residuals per detector (INTT z: strip length)
// Tree "matching": every candidate Si-TPC pair stored by the matcher
//   (si_index, tpc_index, dphi, deta, dz0, chi2, accepted) - to tune the match windows.
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
  void setOutputFileName(const std::string& n) { m_outfilename = n; }
  void setMinPt(double v) { m_minPt = v; }
  void setMaxPt(double v) { m_maxPt = v; }
  void setMinTpcClusters(unsigned int v) { m_minTpcClusters = v; }
  void setMinSiPoints(unsigned int v) { m_minSiPoints = v; }
  void setWriteMatchingTree(bool v) { m_writeMatching = v; }

 private:
  void reset_tree_values();

  std::string m_outfilename;
  std::string m_trackNodeName{"SITPC_TRACKS"};
  std::string m_clusterNodeName{"TPC_POLYCLUSTERS"};

  double m_minPt{0.0};
  double m_maxPt{1.0e30};
  unsigned int m_minTpcClusters{0};
  unsigned int m_minSiPoints{0};
  bool m_writeMatching{true};

  unsigned int m_evt{0};
  TFile* m_outfile{nullptr};
  TTree* m_tree{nullptr};
  TTree* m_matchTree{nullptr};
  SiTpc_TrackContainer* m_tracks{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};

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
  static constexpr int kNLayers = 55;  // 0-2 MVTX, 3-6 INTT, 7-54 TPC
  double m_resRPhiLayer[kNLayers]{};
  double m_resZLayer[kNLayers]{};
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
};

#endif