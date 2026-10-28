#pragma once

// SiSeedComparisonQA
//
// 1) Cluster vs hit (tree "si_points"): for every layer point of the own Si chains (SI_TRAJECTORY, built from
//    hits) the official TRKR_CLUSTER cluster sharing most of its hits (TRKR_CLUSTERHITASSOC):
//      - hit bookkeeping: chain hits in the layer, sensors, clusters touched, hits shared with the best cluster,
//        hits of the best cluster not in the chain (merge), chain hits outside it (split)
//      - position: cluster - chain point, in the own frame (r*phi, z, r) [cm], with the official cluster
//        errors (ClusterErrorPara) for scale
//    Official -> own frame: x, y - SiBeam (+ INTT physics offset for the INTT), z + setSiActsZShiftCm.
//
// 2) Seed bias (tree "seed_bias"): every own Si trajectory refitted with SiTpcHelixFit twice, with the own points
//    and with the matched official cluster positions (same fit, same config as SiTrajectoryFitter): the change
//    of d0, phi, z0, tanl, pT is what using clusters instead of hits does to the Si seed.  For the trajectories
//    in SITPC_TRACKS the Si+TPC (hybrid) refit of SiTpcTrackMatcher is repeated the same way.
//    ref_* = refit with the own points - stored value: ~0 if the configs here match the macro.
//
// 3) Pattern recognition vs the standard silicon seeding (trees "own_chains", "std_seeds"), when
//    SiliconTrackSeedContainer exists (Tracking_Reco_SiliconSeed_run2pp + TrackSeedTrackMapConverter):
//      own chain  -> official seed with most shared clusters (n shared, n own, n official, parameters of both)
//      std seed   -> own chain with most shared clusters, and the closest own Si+TPC track (dphi, deta)
//    Parameters: pT, phi, eta, d0 to the beam line of each frame, z0 (official frame), crossing
//    (own: INTT time bucket of the chain hits / resolver decision; official: seed crossing).
//
// 4) Si-only vertices (tree "si_vertices"): own SiPolyVertexMap vs official SiliconSvtxVertexMap, matched in z
//    within the same crossing, positions compared in the official frame.

#include <fun4all/SubsysReco.h>

#include <array>
#include <memory>
#include <string>

class PHCompositeNode;
class TFile;
class TTree;

class SiSeedComparisonQA : public SubsysReco
{
 public:
  explicit SiSeedComparisonQA(const std::string& name = "SiSeedComparisonQA",
                              const std::string& fileName = "si_seed_comparison_qa.root");
  ~SiSeedComparisonQA() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // ---- own nodes
  void setSiTrajectoryNodeName(const std::string& s) { m_trajNode = s; }
  void setSiTpcTrackNodeName(const std::string& s) { m_siTpcNode = s; }
  void setOwnSiVertexNodeName(const std::string& s) { m_ownVertexNode = s; }
  void setCrossingDecisionNodeName(const std::string& s) { m_crossingNode = s; }
  // ---- standard nodes
  void setStdSeedNodeName(const std::string& s) { m_stdSeedNode = s; }
  void setStdTrackMapName(const std::string& s) { m_stdTrackNode = s; }
  void setStdVertexMapName(const std::string& s) { m_stdVertexNode = s; }

  // ---- frames (same numbers as SiHitSeedReco / SiTpcActsFitter) [cm]
  // shift of the own Si frame: own = official - SiBeam (SiTpcBeamAlignment SiFull, 0 without the correction)
  void setSiBeamSpotCm(double x, double y) { m_siBeam = {x, y}; }
  // beam line in the official frame (d0 of the standard seeds; own frame: this - SiBeam)
  void setOfficialBeamLineCm(double x0, double y0, double dxdz = 0, double dydz = 0) { m_beam = {x0, y0, dxdz, dydz}; }
  void setSiActsZShiftCm(double dz) { m_zShift = dz; }
  void setInttPhysicsOffsetCm(double dx, double dy) { m_inttOffset = {dx, dy}; }

  // ---- refits (SiTpcHelixFit), as in the macro
  void setSiFitConfig(double minLeverArm, double zWeightMvtx, double zWeightIntt, double bz = 1.4)
  {
    m_siLever = minLeverArm;
    m_siZw = {zWeightMvtx, zWeightIntt};
    m_bz = bz;
  }
  void setSiTpcFitConfig(std::array<double, 3> wxy, std::array<double, 3> wz, bool hybrid, double hybridWMvtx,
                         double hybridWIntt, unsigned int hybridMinPoints, bool tpcZOffset)
  {
    m_stWxy = wxy;
    m_stWz = wz;
    m_hybrid = hybrid;
    m_hybridW = {hybridWMvtx, hybridWIntt};
    m_hybridMin = hybridMinPoints;
    m_tpcZOffset = tpcZOffset;
  }

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;

  std::string m_fileName;
  std::string m_trajNode{"SI_TRAJECTORY"};
  std::string m_siTpcNode{"SITPC_TRACKS"};
  std::string m_ownVertexNode{"SiPolyVertexMap"};
  std::string m_crossingNode{"TPC_CROSSING_DECISIONS"};
  std::string m_stdSeedNode{"SiliconTrackSeedContainer"};
  std::string m_stdTrackNode{"SiliconSvtxTrackMap"};
  std::string m_stdVertexNode{"SiliconSvtxVertexMap"};

  std::array<double, 2> m_siBeam{0, 0};
  std::array<double, 4> m_beam{0, 0, 0, 0};
  double m_zShift{0};
  std::array<double, 2> m_inttOffset{0, 0};

  double m_siLever{3.0};
  std::array<double, 2> m_siZw{1.0, 1.0e-4};
  double m_bz{1.4};
  std::array<double, 3> m_stWxy{0, 0, 1};
  std::array<double, 3> m_stWz{1.0, 1.0e-4, 1.0e-4};
  bool m_hybrid{true};
  std::array<double, 2> m_hybridW{1.0, 0.04};
  unsigned int m_hybridMin{2};
  bool m_tpcZOffset{true};
};
