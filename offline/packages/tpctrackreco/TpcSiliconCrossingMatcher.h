#ifndef TPCTRACKRECO_TPCSILICONCROSSINGMATCHER_H
#define TPCTRACKRECO_TPCSILICONCROSSINGMATCHER_H

#include "BeamFrameTransform.h"
#include "TpcTrackFit.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <limits>
#include <set>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class TpcCrossingTrajectory;
class TpcCrossingTrajectoryContainer;
class TpcSiliconMatchCandidateContainer;
class TrkrClusterContainer;

class TpcSiliconCrossingMatcher : public SubsysReco
{
 public:
  explicit TpcSiliconCrossingMatcher(const std::string& name = "TpcSiliconCrossingMatcher");
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setTrajectoryNodeName(const std::string& value) { m_trajectoryNodeName = value; }
  void setClusterNodeName(const std::string& value) { m_clusterNodeName = value; }
  void setOutputNodeName(const std::string& value) { m_outputNodeName = value; }
  void setTpcBeamLine(double x0, double dxdz, double y0, double dydz) { m_beamFrame.setTpcBeamLine({x0, dxdz, y0, dydz}); }
  void setMvtxBeamLine(double x0, double dxdz, double y0, double dydz) { m_beamFrame.setMvtxBeamLine({x0, dxdz, y0, dydz}); }
  void setInttBeamLine(double x0, double dxdz, double y0, double dydz) { m_beamFrame.setInttBeamLine({x0, dxdz, y0, dydz}); }
  void setZSearchTimeBins(double value) { m_zSearchTimeBins = value; }
  void setTpcAdcClockNs(double value) { m_tpcAdcClockNs = value; }
  void setLooseWindow(double rdphi, double dz) { m_looseRdphiWindow = rdphi; m_looseDzWindow = dz; }
  void setInttWindow(double rdphi, double dz) { m_inttRdphiWindow = rdphi; m_inttDzWindow = dz; }
  void setInttRdphiWindow(double value) { m_inttRdphiWindow = value; }
  void setInttDzWindow(double value) { m_inttDzWindow = value; }
  void setPhiThetaWindowSigma(double phi, double theta) { m_phiWindowSigma = phi; m_thetaWindowSigma = theta; }
  void setAngularResidualSigma(double phi, double theta) { m_sigmaPhi = phi; m_sigmaTheta = theta; }
  void setUseDynamicResiduals(bool value) { m_useDynamicResiduals = value; }
  void setAssociationCalibrationMode(bool value) { m_associationCalibrationMode = value; }
  void setDynamicPhiMean(unsigned int layer, double offset, double slope)
  { if (layer < 7U) { m_dynamicPhiMeanOffset[layer] = offset; m_dynamicPhiMeanSlope[layer] = slope; } }
  void setDynamicThetaMean(unsigned int layer, double offset, double slope)
  { if (layer < 7U) { m_dynamicThetaMeanOffset[layer] = offset; m_dynamicThetaMeanSlope[layer] = slope; } }
  void setMinSiliconClusters(unsigned int value) { m_minSiliconClusters = value; }
  void setApplyChainDcaCut(bool value) { m_applyChainDcaCut = value; }
  void setMaxChainDcaScore(double value) { m_maxChainDcaScore = value; }
  void setMaxChainDeltaEta(double value) { m_maxChainDeltaEta = value; }
  void setMaxChains(unsigned int value) { m_maxChains = value; }
  void setMaxBranchesPerLayer(unsigned int value) { m_maxBranchesPerLayer = value; }
  void setMidpointCompatibilityLimits(double rdphi, double dz, double phi, double tanLambda)
  {
    m_maxMidpointAbsRdphi = rdphi;
    m_maxMidpointAbsDz = dz;
    m_maxMidpointAbsPhi = phi;
    m_maxMidpointAbsTanLambda = tanLambda;
  }
  void setMaxMidpointScore(double value) { m_maxMidpointScore = value; }
  void setMaxRejectedTrajectoryPrints(unsigned int value) { m_maxRejectedTrajectoryPrints = value; }
  void setSurfaceResidualWindows(double mvtxLocal0, double mvtxLocal1, double inttLocal0)
  { m_mvtxLocal0Window = mvtxLocal0; m_mvtxLocal1Window = mvtxLocal1; m_inttLocal0Window = inttLocal0; }
  void setMaxSurfaceQaPrints(unsigned int value) { m_maxSurfaceQaPrints = value; }

 private:
  struct SpacePoint
  {
    TrkrDefs::cluskey key{TrkrDefs::CLUSKEYMAX};
    unsigned int layer{0};
    double x{0.};
    double y{0.};
    double z{0.};
    double r{0.};
    double phi{0.};
    double global_x{0.};
    double global_y{0.};
    double global_z{0.};
  };

  struct SurfaceMatch
  {
    std::array<double, 3> intersection{};
    std::array<double, 3> surface_center{};
    double local_residual_0{0.};
    double local_residual_1{0.};
    double path_length_cm{0.};
  };

  struct TrajectoryState
  {
    double phi_intercept{0.};
    double phi_slope{0.};
    double phi_S{0.};
    double phi_x0{0.};
    double phi_invR{0.};
    double phi_theta{0.};
    double phi_bline{0.};
    double z_intercept{0.};
    double z_slope{0.};
    double seed_x0{0.};
    double seed_y0{0.};
    double seed_z0{0.};
    double seed_cx{0.};
    double seed_cy{0.};
    double seed_phi0{0.};
    double seed_slope{0.};
    double seed_q_over_r{0.};
    bool use_tpc_seed{false};
    bool use_silicon_seed{false};
    bool phi_sagitta_ok{false};
    bool valid{false};
  };

  struct ChainHit
  {
    SpacePoint point;
    double pred_phi{0.};
    double pred_z{0.};
    double dphi{0.};
    double dtheta{0.};
    double rdphi{0.};
    double dz{0.};
    double ddphi{0.};
    double chi2{0.};
    double surface_residual_0{0.};
    double surface_residual_1{0.};
  };

  struct Chain
  {
    std::vector<ChainHit> hits;
    TrajectoryState tpc_reference;
    TrajectoryState si_reference;
    TrajectoryState state;
    double chi2{0.};
    double score{0.};
    double midpoint_score{std::numeric_limits<double>::max()};
    double dca_score{std::numeric_limits<double>::max()};
    double delta_eta0{std::numeric_limits<double>::max()};
    double previous_dphi{0.};
    double previous_dtheta{0.};
    bool has_previous_residual{false};
    unsigned int n_missing{0};
    double pt{0.};
    double r_si_outer{0.};
    double r_tpc_inner{0.};
    double r_match{0.};
    double midpoint_delta_rphi{0.};
    double midpoint_delta_z{0.};
    double midpoint_delta_phi{0.};
    double midpoint_delta_tan_lambda{0.};
    bool midpoint_compatible{false};
    std::array<double, 4> midpoint_tpc{};
    std::array<double, 4> midpoint_si{};
  };

  struct Counters
  {
    unsigned int trajectoriesSeen{0};
    std::array<unsigned long long, 7> search{};
    std::array<unsigned long long, 7> zPass{};
    std::array<unsigned long long, 7> phiPass{};
    unsigned int l2Seeds{0};
    unsigned int candidates1Mvtx{0};
    unsigned int candidates2Mvtx{0};
    unsigned int candidates3Mvtx{0};
    unsigned int candidatesWithIntt{0};
  };

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  std::vector<SpacePoint> collectSiliconClusters() const;
  TrajectoryState makeTpcReferenceTrajectory(const TpcCrossingTrajectory&) const;
  TrajectoryState fitTrajectory(const std::vector<SpacePoint>&) const;
  TrajectoryState makeSiliconSeedTrajectory(const TrajectoryState&, const SpacePoint&) const;
  TrajectoryState correctSiliconSeedWithTwoHits(const TrajectoryState&, const ChainHit&, const ChainHit&) const;
  TrajectoryState correctSiliconSeedWithAllHits(const TrajectoryState&, const std::vector<ChainHit>&) const;
  bool predictAtRadius(const TrajectoryState&, double, double&, double&, double&, double&) const;
  bool matchToSurface(const TpcCrossingTrajectory&, const SpacePoint&, SurfaceMatch&) const;
  std::vector<ChainHit> findMvtxCandidates(const Chain&, const std::vector<SpacePoint>&,
                                           const std::set<TrkrDefs::cluskey>&, unsigned int,
                                           Counters&) const;
  std::vector<Chain> buildChains(const TpcCrossingTrajectory&, const std::vector<SpacePoint>&, Counters&) const;
  const Chain* selectBestChain(const std::vector<Chain>&) const;
  Chain attachClosestInttClusters(const Chain&, const std::vector<SpacePoint>&, Counters&) const;
  bool computeChainDcaMetrics(Chain&, const TrajectoryState&) const;
  bool computeMidpointMatch(Chain&, const TpcCrossingTrajectory&) const;
  bool propagateTpcToRadius(const TpcCrossingTrajectory&, double, std::array<double, 6>&) const;
  double wrapPhi(double) const;
  double unwrapPhiNear(double, double) const;
  double predictSagittaPhi(double, const TrajectoryState&) const;
  double pointTheta(const SpacePoint&) const;
  double trajectoryPhi0NearBeam(const TrajectoryState&) const;
  double trajectoryZ0NearBeam(const TrajectoryState&) const;
  double trajectoryTheta0NearBeam(const TrajectoryState&) const;
  double dynamicMeanPhi(unsigned int, double, bool) const;
  double dynamicSigmaPhi(double) const;
  double dynamicMeanTheta(unsigned int, double, bool) const;
  double dynamicSigmaTheta(double, double) const;
  double dynamicDzWindow(double) const;
  double zSearchWindowCm() const;

  std::string m_trajectoryNodeName{"TPC_CROSSING_TRAJECTORIES"};
  std::string m_clusterNodeName{"TRKR_CLUSTER"};
  std::string m_outputNodeName{"TPC_SILICON_MATCH_CANDIDATES"};
  BeamFrameTransform m_beamFrame;
  double m_zSearchTimeBins{2.0};
  double m_tpcAdcClockNs{56.881262};
  double m_looseRdphiWindow{0.15};
  double m_looseDzWindow{0.5};
  double m_inttRdphiWindow{0.25};
  double m_inttDzWindow{1.0};
  double m_sigmaPhi{0.5};
  double m_sigmaTheta{0.2};
  double m_phiWindowSigma{0.7};
  double m_thetaWindowSigma{0.2};
  double m_missingLayerPenalty{1.0};
  double m_maxChainDcaScore{5.0};
  double m_maxChainDeltaEta{0.2};
  double m_mvtxLocal0Window{2};
  double m_mvtxLocal1Window{3};
  double m_inttLocal0Window{2};
  double m_inttLocal1Window{3};
  bool m_applyChainDcaCut{true};
  bool m_useDynamicResiduals{true};
  bool m_associationCalibrationMode{true};
  std::array<double, 7> m_dynamicPhiMeanOffset{};
  std::array<double, 7> m_dynamicPhiMeanSlope{{1., 1., 1., 1., 1., 1., 1.}};
  std::array<double, 7> m_dynamicThetaMeanOffset{};
  std::array<double, 7> m_dynamicThetaMeanSlope{{1., 1., 1., 1., 1., 1., 1.}};
  std::array<double, 3> m_vertexPhiMean{{0.00720456, 0.00108784, 0.}};
  std::array<double, 3> m_vertexPhiSigma{{0.56235, 1.17511, 1.}};
  std::array<double, 3> m_vertexThetaMean{{0.000410664, 0.00323673, 0.}};
  std::array<double, 3> m_vertexThetaSigma{{0.209156, 0.193329, 1.}};
  unsigned int m_minSiliconClusters{0};
  unsigned int m_maxChains{256};
  unsigned int m_maxBranchesPerLayer{8};
  double m_maxMidpointAbsRdphi{2.0};
  double m_maxMidpointAbsDz{3.0};
  double m_maxMidpointAbsPhi{0.2};
  double m_maxMidpointAbsTanLambda{0.25};
  double m_maxMidpointScore{16.0};
  unsigned int m_maxRejectedTrajectoryPrints{10};
  unsigned int m_maxSurfaceQaPrints{20};
  TpcKalmanConfig m_propagationConfig;
  std::vector<unsigned int> m_matchLayers{2, 1, 0};
  std::vector<unsigned int> m_inttMatchLayers{3, 4, 5, 6};
  TpcCrossingTrajectoryContainer* m_trajectories{nullptr};
  TrkrClusterContainer* m_clusters{nullptr};
  ActsGeometry* m_geometry{nullptr};
  TpcSiliconMatchCandidateContainer* m_candidates{nullptr};
};

#endif
