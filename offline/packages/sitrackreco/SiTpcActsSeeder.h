#pragma once

// SiTpcActsSeeder
//
// Prepares the own Si+TPC tracks (SITPC_TRACKS from SiTpcTrackMatcher) for the OFFICIAL ACTS fitter
// (PHActsTrkFitter), i.e. writes what PHSiliconTpcTrackMatching + the standard clustering would have written:
//
//   - TrkrClusterContainer (default "TRKR_CLUSTER_SITPC", give it to PHActsTrkFitter::setTrkrClusterContainerName)
//       Si : copies of the official TRKR_CLUSTER clusters of the Si chain (same keys).  Per layer the cluster that
//            contains most chain hits (TRKR_CLUSTERHITASSOC, with the MVTX strobe-duplicate recovery), else the
//            closest official cluster of the sensor, else the chain point projected on the official surface
//            (new cluster, key index above the official ones of that sensor).
//            INTT physics offset (setInttPhysicsOffsetCm): added to the INTT local x (in-plane component).
//       TPC: one TrkrClusterv5 per poly cluster of the track (key = the poly cluster trkr key), written such that
//            the official chain (ActsGeometry::getGlobalPosition with the crossing of the seed, no TPC distortion
//            correction containers) gives back the own position of the cluster in the official frame:
//              position = TPC_raw - TpcBeamLine + SiBeam (x, y),  TPC_raw - siActsZShift - per-track TPC z offset (z)
//                         + setTpcExtraShiftCm
//            The official TPC surfaces sit at the pad-row radius, the own (PHGarfield corrected) clusters do not:
//            the cluster is moved to the surface along the own Si+TPC helix (TpcPlacement::AlongTrack, what the
//            official TpcClusterMover does) or along the surface normal (TpcPlacement::Normal).  The local z is
//            converted to the cluster time with the ActsGeometry drift parameters and the seed crossing, so the
//            crossing correction of the official fitter gives back exactly this z.
//            Errors: the poly-cluster base errors (rms_x: r*phi, rms_z: z) and sizes, the official fitter applies
//            ClusterErrorPara::get_clusterv5_modified_error on top (as for the standard clusters).
//
//   - SiliconTrackSeedContainer / TpcTrackSeedContainer (fixed names required by PHActsTrkFitter; appended to if
//     the standard seeding made them) with TrackSeed_v2 (circle centre, q/R, z0, tan(lambda), phi in the official
//     frame from the own Si+TPC fit, charge from the bending and the sign of B_z) and the cluster keys above.
//
//   - SvtxTrackSeedContainer (default "SiTpcSvtxTrackSeedContainer", PHActsTrkFitter::set_svtx_seed_map_name) with
//     SvtxTrackSeed_v3: Si seed index, TPC seed index, crossing (= crossing estimate) of the own track.
//
// The TPC distortion correction containers must NOT be loaded (G4TPC::ENABLE_*_CORRECTIONS = false): the own
// positions are already corrected (PHGarfield).  InitRun / the first event print a warning if they are there.
//
// Optional QA tree (setQaFileName): per TPC cluster the move to the surface and the closure
// |official global position - intended position|, per track the seed parameters.

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <array>
#include <memory>
#include <string>

class PHCompositeNode;

class SiTpcActsSeeder : public SubsysReco
{
 public:
  enum class TpcPlacement
  {
    AlongTrack = 0,  // move the cluster to the official surface along the own helix (as TpcClusterMover)
    Normal = 1       // project the cluster on the official surface along its normal
  };

  explicit SiTpcActsSeeder(const std::string& name = "SiTpcActsSeeder");
  ~SiTpcActsSeeder() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // ---- input nodes
  void setSiTpcTrackNodeName(const std::string& s) { m_siTpcTrackNodeName = s; }
  void setSiTrajectoryNodeName(const std::string& s) { m_siTrajectoryNodeName = s; }
  void setTpcClusterNodeName(const std::string& s) { m_tpcClusterNodeName = s; }
  void setCrossingDecisionNodeName(const std::string& s) { m_crossingNodeName = s; }
  void setSiClusterNodeName(const std::string& s) { m_siClusterNodeName = s; }

  // ---- output nodes
  void setClusterContainerName(const std::string& s) { m_clusterOutName = s; }
  void setSvtxSeedContainerName(const std::string& s) { m_svtxSeedName = s; }

  // ---- frames [cm], as SiTpcActsFitter
  void setSiBeamSpotCm(double x, double y) { m_siBeam = {x, y}; }
  void setTpcBeamLineCm(double x, double y) { m_tpcLine = {x, y}; }
  void setSiActsZShiftCm(double dz) { m_siActsZShift = dz; }
  void setTpcExtraShiftCm(double dx, double dy, double dz) { m_tpcExtraShift = {dx, dy, dz}; }
  void setApplyTrackTpcZOffset(bool b) { m_applyTrackTpcZOffset = b; }
  void setInttPhysicsOffsetCm(double dx, double dy) { m_inttOffset = {dx, dy}; }

  // ---- measurements
  void setTpcPlacement(TpcPlacement p) { m_tpcPlacement = p; }
  // drop TPC clusters crossed at an angle (track direction vs surface normal) above this [rad]
  void setMaxTpcCrossingAngle(double a) { m_maxTpcCrossingAngle = a; }
  // drop TPC clusters moved by more than this to reach the official surface [cm]
  void setMaxTpcMoveCm(double d) { m_maxTpcMove = d; }
  // use the projected chain point when no official Si cluster is found (else the layer is skipped)
  void setUseProjectedChainPoints(bool b) { m_useProjected = b; }

  // ---- track selection
  void setMinMvtx(unsigned int n) { m_minMvtx = n; }
  void setMinSi(unsigned int n) { m_minSi = n; }
  void setMinTpc(unsigned int n) { m_minTpc = n; }
  void setMinPt(double pt) { m_minPt = pt; }

  // ---- QA
  void setQaFileName(const std::string& s) { m_qaFileName = s; }

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;

  int getNodes(PHCompositeNode* topNode);
  int createNodes(PHCompositeNode* topNode);
  void checkDistortionContainers(PHCompositeNode* topNode);

  std::string m_siTpcTrackNodeName{"SITPC_TRACKS"};
  std::string m_siTrajectoryNodeName{"SI_TRAJECTORY"};
  std::string m_tpcClusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_crossingNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_siClusterNodeName{"TRKR_CLUSTER"};

  std::string m_clusterOutName{"TRKR_CLUSTER_SITPC"};
  std::string m_svtxSeedName{"SiTpcSvtxTrackSeedContainer"};

  std::array<double, 2> m_siBeam{0, 0};
  std::array<double, 2> m_tpcLine{0, 0};
  double m_siActsZShift{0};
  std::array<double, 3> m_tpcExtraShift{0, 0, 0};
  bool m_applyTrackTpcZOffset{true};
  std::array<double, 2> m_inttOffset{0, 0};

  TpcPlacement m_tpcPlacement{TpcPlacement::AlongTrack};
  double m_maxTpcCrossingAngle{1.2};
  double m_maxTpcMove{2.0};
  bool m_useProjected{true};

  unsigned int m_minMvtx{2};
  unsigned int m_minSi{2};
  unsigned int m_minTpc{10};
  double m_minPt{0.1};

  std::string m_qaFileName;
};
