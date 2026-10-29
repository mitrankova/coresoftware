#pragma once

// SiTpcActsFitter
//
// ACTS Kalman fit (direct navigation) of the Si+TPC tracks found by the own pattern recognition
// (SITPC_TRACKS from SiTpcTrackMatcher), with
//   - Si  : the Si chain measured on the OFFICIAL ACTS silicon surfaces.  Per layer the chain hits are
//           mapped to the official TRKR_CLUSTER cluster that contains most of them (TRKR_CLUSTERHITASSOC);
//           without such a cluster the chain point itself is projected on the official surface of its
//           sensor (SiMeasurement::ChainPoints does that for every layer).
//   - TPC : one ARTIFICIAL plane surface per TPC poly cluster, placed through the (PHGarfield corrected)
//           cluster position, normal = radial direction (pad-row geometry), local axes (r*phi, z).  The
//           measurement is the cluster itself (local (0, 0)).  The material of the official TPC surface at
//           that position is copied onto the plane, the non-sensitive official material surfaces (beam
//           pipe, field cages, ...) are put into the navigation sequence as in PHActsTrkFitter.
//
// Frames (all in cm, see the macro for the numbers):
//   The fit runs in the OFFICIAL global frame (official ACTS geometry + alignment), i.e. the silicon is
//   used as it is and the TPC is brought into the silicon frame:
//     own pattern-recognition frame:  Si_own  = Si_official (+ z: setSiActsZShiftCm) - SiBeam
//                                     TPC_own = TPC_raw - TpcBeamLine               (TPC_raw = poly centroid)
//     both own frames have the beam at (0, 0), so the TPC in the official frame is
//       TpcFrame::BeamAxis      : x, y = TPC_raw - TpcBeamLine + SiBeam
//                                 z    = TPC_raw - siActsZShift - (per track TPC z offset, optional)
//       TpcFrame::OfficialWorld : official TPC envelope -> world transform (no beam-axis matching)
//     + setTpcExtraShiftCm (fine tuning, e.g. from the QA fits).
//   The beam is then at SiBeam for both the Si and the TPC.  The fitted perigee is a line along z through
//   (SiBeam x, y), so d0 / z0 of the output are relative to the beam axis.
//   INTT: setInttPhysicsOffsetCm moves the INTT surfaces (transient alignment transforms, reset after every
//   event) by the same offset as SiHitSeedReco::setInttPhysicsOffsetCm.
//
// Output: SvtxTrackMap (default "SvtxTrackMap") with SvtxTrack_v4 (states on all measurements, Si/TPC seeds
// attached, so PHSimpleVertexFinder / KFParticle can run on it) + two TrackSeedContainers.
// Optional QA file: per track tree (full fit + optional TPC-only and Si-only fits at the same beam perigee:
// d0_tpc - d0_si vs phi measures the remaining Si - TPC offset) and per measurement residual tree.

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <array>
#include <memory>
#include <string>

class PHCompositeNode;

class SiTpcActsFitter : public SubsysReco
{
 public:
  enum class TpcFrame
  {
    BeamAxis = 0,
    OfficialWorld = 1
  };

  enum class TpcErrors
  {
    Standard = 0,  // ClusterErrorPara::get_clusterv5_modified_error (as PHActsTrkFitter) on the poly-cluster base errors
    Constant = 1   // setTpcErrorCm per module
  };

  enum class SiMeasurement
  {
    OfficialClusters = 0,  // official TRKR_CLUSTER clusters, fallback: projected chain point
    ChainPoints = 1        // the chain points projected on the official surfaces
  };

  explicit SiTpcActsFitter(const std::string& name = "SiTpcActsFitter");
  ~SiTpcActsFitter() override;

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
  void setTrackMapName(const std::string& s) { m_trackMapName = s; }
  void setSiliconSeedContainerName(const std::string& s) { m_siSeedName = s; }
  void setTpcSeedContainerName(const std::string& s) { m_tpcSeedName = s; }

  // ---- frames [cm]
  // Si beam spot used by the own reconstruction (SiTpcBeamAlignment SiFull x0, y0)
  void setSiBeamSpotCm(double x, double y) { m_siBeam = {x, y}; }
  // TPC beam line used by the own reconstruction (SiTpcBeamAlignment Tpc x0, y0)
  void setTpcBeamLineCm(double x, double y) { m_tpcLine = {x, y}; }
  // the z added to the official Si positions by SiHitSeedReco::setActsZShiftCm (+0.666 cm = -D_z)
  void setSiActsZShiftCm(double dz) { m_siActsZShift = dz; }
  void setTpcFrame(TpcFrame f) { m_tpcFrame = f; }
  // additional TPC translation in the official frame
  void setTpcExtraShiftCm(double dx, double dy, double dz) { m_tpcExtraShift = {dx, dy, dz}; }
  // per-track TPC z offset fitted by SiTpcTrackMatcher (setFitTpcZOffset): subtract it from the TPC z
  void setApplyTrackTpcZOffset(bool b) { m_applyTrackTpcZOffset = b; }
  // INTT relative to the MVTX, as SiHitSeedReco::setInttPhysicsOffsetCm
  void setInttPhysicsOffsetCm(double dx, double dy) { m_inttOffset = {dx, dy}; }

  // ---- measurements
  void setSiMeasurement(SiMeasurement m) { m_siMeasurement = m; }
  // errors of projected chain points [cm] (official clusters use ClusterErrorPara)
  void setMvtxChainPointErrorCm(double rphi, double z) { m_mvtxChainErr = {rphi, z}; }
  void setInttChainPointErrorCm(double rphi, double z) { m_inttChainErr = {rphi, z}; }
  // TPC cluster errors [cm]: sigma_rphi, sigma_z (per TPC module: inner, mid, outer)
  void setTpcErrorCm(double rphi, double z) { m_tpcErr = {{{rphi, z}, {rphi, z}, {rphi, z}}}; }
  void setTpcErrorCm(int module, double rphi, double z)
  {
    if (module >= 0 && module < 3)
    {
      m_tpcErr[module] = {rphi, z};
    }
  }
  void setTpcErrorMode(TpcErrors m) { m_tpcErrorMode = m; }
  // scale of the standard TPC errors (r*phi, z)
  void setTpcErrorScale(double rphi, double z) { m_tpcErrorScale = {rphi, z}; }
  // drop TPC clusters crossed at an angle (track direction vs radial normal) above this [rad]
  void setMaxTpcCrossingAngle(double a) { m_maxTpcCrossingAngle = a; }
  // copy the official TPC surface material onto the artificial planes
  void setCopyTpcMaterial(bool b) { m_copyTpcMaterial = b; }
  // include the non-sensitive official material surfaces in the navigation sequence
  void setUseMaterialSurfaces(bool b) { m_useMaterialSurfaces = b; }

  // ---- track selection
  void setMinMvtx(unsigned int n) { m_minMvtx = n; }
  void setMinSi(unsigned int n) { m_minSi = n; }
  void setMinTpc(unsigned int n) { m_minTpc = n; }
  void setMinPt(double pt) { m_minPt = pt; }

  // ---- fit
  void setMultipleScattering(bool b) { m_multipleScattering = b; }
  void setEnergyLoss(bool b) { m_energyLoss = b; }
  // initial covariance: sigma(d0, z0) [cm], sigma(phi, theta) [rad], relative sigma(q/p)
  void setSeedSigmas(double pos, double angle, double relQop)
  {
    m_seedSigmaPos = pos;
    m_seedSigmaAngle = angle;
    m_seedSigmaRelQop = relQop;
  }
  void setFillSvtxTrackStates(bool b) { m_fillStates = b; }
  // KF outlier rejection: predicted chi2 (2 dof) above the cut -> measurement not used; <= 0: off.
  // Only meaningful once the measurement errors are right (the pulls of the QA).
  void setTpcOutlierChi2(double c) { m_tpcOutlierChi2 = c; }
  void setSiOutlierChi2(double c) { m_siOutlierChi2 = c; }
  // a failed fit is repeated without the material surfaces (counted in End, flag in the QA tree)
  void setRetryWithoutMaterial(bool b) { m_retryWithoutMaterial = b; }
  // print measurements, navigation sequence and the Acts VERBOSE log of the first n failed fits
  void setDebugFailedFits(int n) { m_debugFailedFits = n; }
  // propagation step limit per fit (Acts default 1000)
  void setMaxSteps(unsigned int n) { m_maxSteps = n; }

  // ---- QA
  void setQaFileName(const std::string& s) { m_qaFileName = s; }
  // also fit TPC only and Si only (same perigee) for every track (QA tree only)
  void setQaSubFits(bool b) { m_qaSubFits = b; }

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;

  int getNodes(PHCompositeNode* topNode);
  int createNodes(PHCompositeNode* topNode);

  std::string m_siTpcTrackNodeName{"SITPC_TRACKS"};
  std::string m_siTrajectoryNodeName{"SI_TRAJECTORY"};
  std::string m_tpcClusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_crossingNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_siClusterNodeName{"TRKR_CLUSTER"};

  std::string m_trackMapName{"SvtxTrackMap"};
  std::string m_siSeedName{"SiTpcActsSiliconSeedContainer"};
  std::string m_tpcSeedName{"SiTpcActsTpcSeedContainer"};

  std::array<double, 2> m_siBeam{0, 0};
  std::array<double, 2> m_tpcLine{0, 0};
  double m_siActsZShift{0};
  TpcFrame m_tpcFrame{TpcFrame::BeamAxis};
  std::array<double, 3> m_tpcExtraShift{0, 0, 0};
  bool m_applyTrackTpcZOffset{true};
  std::array<double, 2> m_inttOffset{0, 0};

  SiMeasurement m_siMeasurement{SiMeasurement::OfficialClusters};
  std::array<double, 2> m_mvtxChainErr{0.0010, 0.0010};
  std::array<double, 2> m_inttChainErr{0.0025, 0.46};
  std::array<std::array<double, 2>, 3> m_tpcErr{{{0.030, 0.060}, {0.030, 0.060}, {0.030, 0.060}}};
  double m_maxTpcCrossingAngle{1.2};
  bool m_copyTpcMaterial{true};
  bool m_useMaterialSurfaces{true};

  unsigned int m_minMvtx{2};
  unsigned int m_minSi{2};
  unsigned int m_minTpc{10};
  double m_minPt{0.1};

  bool m_multipleScattering{true};
  bool m_energyLoss{true};
  double m_seedSigmaPos{0.1};
  double m_seedSigmaAngle{0.02};
  double m_seedSigmaRelQop{0.2};
  bool m_fillStates{true};
  unsigned int m_maxSteps{10000};
  double m_tpcOutlierChi2{0};
  double m_siOutlierChi2{0};
  bool m_retryWithoutMaterial{true};
  int m_debugFailedFits{0};
  TpcErrors m_tpcErrorMode{TpcErrors::Standard};
  std::array<double, 2> m_tpcErrorScale{1.0, 1.0};

  std::string m_qaFileName;
  bool m_qaSubFits{false};
};
