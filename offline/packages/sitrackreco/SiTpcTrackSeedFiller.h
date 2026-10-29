#pragma once

// SiTpcTrackSeedFiller
//
// Puts the own reconstruction into the standard sPHENIX seed containers, in place of the standard silicon
// seeding (Tracking_Reco_SiliconSeed_run2pp) and of PHSiliconTpcTrackMatching, so that the rest of the official
// chain (TPOT matching, PHTpcDeltaZCorrection, PHActsTrkFitter, PHSimpleVertexFinder, TrackResiduals) runs on it:
//
//   SiliconTrackSeedContainer : one TrackSeed_v2 per own Si trajectory (SI_TRAJECTORY)
//       cluster keys: per layer the official TRKR_CLUSTER cluster that contains most chain hits
//       (TRKR_CLUSTERHITASSOC, with the recovery of MVTX clusters removed by MvtxClusterPruner as duplicates of
//       the neighbouring strobe), else the closest official cluster of the sensor to the chain point;
//       optionally (setUseProjectedChainPoints) a new cluster at the projected chain point.
//       parameters: the own Si trajectory helix, moved to the official frame (+ Si beam spot, z - Si acts z shift):
//       circle centre (X0, Y0), q/R (charge from the bending and the sign of B_z), z0 and phi at the point of the
//       circle closest to the origin (TrackSeedHelper::get_xyz), slope = tan(lambda).
//   SvtxTrackSeedContainer    : one SvtxTrackSeed_v3 per own Si+TPC track (SITPC_TRACKS): Si seed index of its
//       trajectory, TPC seed index of its TPC poly track in TpcTrackSeedContainer (TpcPolyTrackSeedConverter,
//       must run before), crossing = crossing estimate = TPC_CROSSING_DECISIONS of its assembled TPC track.
//       Optionally (setAddUnmatchedTpcSeeds) TPC-only entries for the TPC seeds without Si+TPC track, as
//       PHSiliconTpcTrackMatching does.
//
// The TPC clusters and seeds come from the official converters (TpcPolyClusterTrkrClusterConverter,
// TpcPolyTrackSeedConverter): the TPC is in the official frame, the own Si-TPC beam-line matching is not used.

#include <fun4all/SubsysReco.h>

#include <array>
#include <memory>
#include <string>

class PHCompositeNode;

class SiTpcTrackSeedFiller : public SubsysReco
{
 public:
  explicit SiTpcTrackSeedFiller(const std::string& name = "SiTpcTrackSeedFiller");
  ~SiTpcTrackSeedFiller() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // ---- input nodes
  void setSiTpcTrackNodeName(const std::string& s) { m_siTpcTrackNodeName = s; }
  void setSiTrajectoryNodeName(const std::string& s) { m_siTrajectoryNodeName = s; }
  void setTpcPolyTrackNodeName(const std::string& s) { m_tpcPolyTrackNodeName = s; }
  void setCrossingDecisionNodeName(const std::string& s) { m_crossingNodeName = s; }
  void setClusterNodeName(const std::string& s) { m_clusterNodeName = s; }
  void setTpcSeedContainerName(const std::string& s) { m_tpcSeedName = s; }

  // ---- output nodes
  void setSiliconSeedContainerName(const std::string& s) { m_siSeedName = s; }
  void setSvtxSeedContainerName(const std::string& s) { m_svtxSeedName = s; }

  // ---- own frame -> official frame [cm]: Si_own = Si_official (z + siActsZShift) - SiBeam
  void setSiBeamSpotCm(double x, double y) { m_siBeam = {x, y}; }
  void setSiActsZShiftCm(double dz) { m_siActsZShift = dz; }

  // ---- options
  // write a new cluster (into the cluster container) at the projected chain point of layers without official cluster
  void setUseProjectedChainPoints(bool b) { m_useProjected = b; }
  // TPC-only SvtxTrackSeeds for the TPC seeds without own Si+TPC track
  void setAddUnmatchedTpcSeeds(bool b) { m_addUnmatchedTpc = b; }
  // Si seeds: minimum number of MVTX / Si clusters
  void setMinMvtx(unsigned int n) { m_minMvtx = n; }
  void setMinSi(unsigned int n) { m_minSi = n; }

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;

  std::string m_siTpcTrackNodeName{"SITPC_TRACKS"};
  std::string m_siTrajectoryNodeName{"SI_TRAJECTORY"};
  std::string m_tpcPolyTrackNodeName{"TPC_POLYTRACKS"};
  std::string m_crossingNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_clusterNodeName{"TRKR_CLUSTER"};
  std::string m_tpcSeedName{"TpcTrackSeedContainer"};
  std::string m_siSeedName{"SiliconTrackSeedContainer"};
  std::string m_svtxSeedName{"SvtxTrackSeedContainer"};

  std::array<double, 2> m_siBeam{0, 0};
  double m_siActsZShift{0};

  bool m_useProjected{false};
  bool m_addUnmatchedTpc{false};
  unsigned int m_minMvtx{2};
  unsigned int m_minSi{2};
};
