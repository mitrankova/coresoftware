#pragma once

// SiTpcTrackSeedFiller
//
// One module from the own Si chains to the standard sPHENIX seed containers.  It replaces
// SiTrajectoryFitter + SiTpcTrackMatcher (and the standard Si seeding + PHSiliconTpcTrackMatching), so
// the official chain (TPOT matching, PHTpcDeltaZCorrection, PHActsTrkFitter, PHSimpleVertexFinder,
// TrackResiduals) runs directly on the own pattern recognition.
//
// Input
//   SI_HIT_SEED_EVENT          SiHitSeedReco: TPC-driven Si chains.  A selected chain belongs to one TPC track
//                              (SiHitChain::tpc_assembled_track_id), so the Si-TPC association is already done.
//   TpcTrackSeedContainer      TpcPolyTrackSeedConverter (register it before this module)
//   TPC_POLYCLUSTERS           cluster key -> assembled TPC track id of the TPC seeds
//   TPC_CROSSING_DECISIONS     crossing of the assembled TPC track
//   TRKR_CLUSTER (+ TRKR_CLUSTERHITASSOC)   official Si clusters, TPC clusters of TpcPolyClusterTrkrClusterConverter
//
// Per selected chain whose TPC track has a TPC seed and a crossing:
//   Si clusters : per chain layer the official cluster that contains most hits of the chain blob
//                 (TRKR_CLUSTERHITASSOC; MVTX clusters removed by MvtxClusterPruner as duplicates are taken
//                 from the neighbouring strobe), else the closest official cluster of the sensor to the
//                 chain point (SiHitChain::phys_*, own frame -> official with setSiBeamSpotCm /
//                 setSiActsZShiftCm), optionally a new cluster at the projected chain point.
//   Si seed     : hybrid helix fit (SiTpcHelixFit, as SiTpcTrackMatcher::setHybridXY) in the OFFICIAL frame on
//                 the official cluster positions: curvature from the TPC clusters of the TPC seed (with the
//                 crossing), circle centre refitted to the Si clusters with that radius; z(s) with the MVTX /
//                 INTT / TPC weights and a free TPC z offset.  Stored as TrackSeed_v2: circle centre, q/R
//                 (charge from the bending and the sign of B_z), z0 and phi at the point closest to the origin
//                 (TrackSeedHelper::get_xyz = perigee of PHActsTrkFitter), tan(lambda), crossing.
//   SvtxTrackSeed_v3: Si seed index, TPC seed index, crossing (= crossing estimate).
//   Optionally TPC-only SvtxTrackSeeds for TPC seeds without selected chain (setAddUnmatchedTpcSeeds).
//     PHActsTrkFitter fits a seed WITHOUT silicon seed at crossing 0 (official behaviour), which in streaming
//     data puts the TPC clusters off in z by the drift of the real crossing.  setTpcOnlyCrossing:
//       FromDecision (default): the TPC-only seed gets a CROSSING CARRIER silicon seed: a TrackSeed_v2 without
//         cluster keys, with the TPC seed helix (X0, Y0, q/R, phi, tan(lambda)), z0 corrected for the decided
//         crossing and that crossing.  PHActsTrkFitter then fits the TPC (+TPOT) clusters at the decided crossing
//         (its start position = the carrier's perigee).  Only tracks whose decision tier is <= setTpcOnlyMaxTier;
//         the others as setTpcOnlyFallbackZero says.  The carrier seeds have no clusters: count Si tracks by
//         their MVTX / INTT clusters, not by the presence of a silicon seed.
//       Zero: plain TPC-only SvtxTrackSeed (no silicon seed), fitted at crossing 0, as PHSiliconTpcTrackMatching.

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
  void setSiSeedEventNodeName(const std::string& s) { m_siSeedEventName = s; }
  void setTpcClusterNodeName(const std::string& s) { m_tpcClusterNodeName = s; }
  void setCrossingDecisionNodeName(const std::string& s) { m_crossingNodeName = s; }
  void setClusterNodeName(const std::string& s) { m_clusterNodeName = s; }
  void setTpcSeedContainerName(const std::string& s) { m_tpcSeedName = s; }

  // ---- output nodes
  void setSiliconSeedContainerName(const std::string& s) { m_siSeedName = s; }
  void setSvtxSeedContainerName(const std::string& s) { m_svtxSeedName = s; }

  // ---- own frame of the chain points -> official frame [cm]: Si_own = Si_official (z + siActsZShift) - SiBeam.
  // Only used for the position match of chain layers without hit-associated cluster; the Si beam spot is also the
  // reference point of the seed fit.
  void setSiBeamSpotCm(double x, double y) { m_siBeam = {x, y}; }
  void setSiActsZShiftCm(double dz) { m_siActsZShift = dz; }

  // ---- seed fit weights (defaults = the SiTpcTrackMatcher settings of the macro)
  void setXYWeights(double mvtx, double intt, double tpc) { m_wxy = {mvtx, intt, tpc}; }
  void setZWeights(double mvtx, double intt, double tpc) { m_wz = {mvtx, intt, tpc}; }
  void setHybridXY(bool on, double wMvtx = 1.0, double wIntt = 0.04, unsigned int minPoints = 2)
  {
    m_hybrid = on;
    m_hybridW = {wMvtx, wIntt};
    m_hybridMinPoints = minPoints;
  }
  void setFitTpcZOffset(bool b) { m_fitTpcZOffset = b; }

  // ---- options
  // all chains instead of the selected one per TPC track (one SvtxTrackSeed per TPC seed in any case: best rank)
  void setUseSelectedChainsOnly(bool b) { m_selectedOnly = b; }
  // new cluster (into the cluster container) at the projected chain point of layers without official cluster
  void setUseProjectedChainPoints(bool b) { m_useProjected = b; }
  // TPC-only SvtxTrackSeeds for TPC seeds without chain
  void setAddUnmatchedTpcSeeds(bool b) { m_addUnmatchedTpc = b; }
  enum class TpcOnlyCrossing
  {
    Zero,
    FromDecision
  };
  void setTpcOnlyCrossing(TpcOnlyCrossing m) { m_tpcOnlyCrossing = m; }
  void setTpcOnlyMaxTier(unsigned int t) { m_tpcOnlyMaxTier = t; }
  // FromDecision, tracks without a usable decision (tier above setTpcOnlyMaxTier / no crossing): true = plain
  // crossing-0 TPC-only seed, false (default) = not added
  void setTpcOnlyFallbackZero(bool b) { m_tpcOnlyFallbackZero = b; }
  void setMinMvtx(unsigned int n) { m_minMvtx = n; }
  void setMinSi(unsigned int n) { m_minSi = n; }

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;

  std::string m_siSeedEventName{"SI_HIT_SEED_EVENT"};
  std::string m_tpcClusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_crossingNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_clusterNodeName{"TRKR_CLUSTER"};
  std::string m_tpcSeedName{"TpcTrackSeedContainer"};
  std::string m_siSeedName{"SiliconTrackSeedContainer"};
  std::string m_svtxSeedName{"SvtxTrackSeedContainer"};

  std::array<double, 2> m_siBeam{0, 0};
  double m_siActsZShift{0};

  std::array<double, 3> m_wxy{0.0, 0.0, 1.0};      // MVTX, INTT, TPC
  std::array<double, 3> m_wz{1.0, 1.0e-4, 1.0e-4};  // MVTX, INTT, TPC
  bool m_hybrid{true};
  std::array<double, 2> m_hybridW{1.0, 0.04};  // MVTX, INTT
  unsigned int m_hybridMinPoints{2};
  bool m_fitTpcZOffset{true};

  bool m_selectedOnly{true};
  bool m_useProjected{false};
  bool m_addUnmatchedTpc{false};
  TpcOnlyCrossing m_tpcOnlyCrossing{TpcOnlyCrossing::FromDecision};
  unsigned int m_tpcOnlyMaxTier{1};
  bool m_tpcOnlyFallbackZero{false};
  unsigned int m_minMvtx{2};
  unsigned int m_minSi{2};
};