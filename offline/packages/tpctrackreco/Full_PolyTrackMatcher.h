// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_FULLPOLYTRACKMATCHER_H
#define TPCTRACKRECO_FULLPOLYTRACKMATCHER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

class ActsGeometry;
class Full_PolyTrackContainer;
class PHCompositeNode;
class TFile;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TpcCrossingDecision;
class TpcCrossingDecisionContainer;
class Tpc_PolyTrack;
class Tpc_PolyTrackContainer;
class TrkrClusterCrossingAssoc;
class TrkrCluster;
class TrkrClusterContainer;

class Full_PolyTrackMatcher : public SubsysReco
{
 public:
  explicit Full_PolyTrackMatcher(const std::string& name = "Full_PolyTrackMatcher");
  ~Full_PolyTrackMatcher() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setTpcTrackNodeName(const std::string& n) { m_tpcTrackNodeName = n; }
  void setTpcClusterNodeName(const std::string& n) { m_tpcClusterNodeName = n; }
  void setCrossingDecisionNodeName(const std::string& n) { m_crossingDecisionNodeName = n; }
  void setClusterCrossingAssocNodeName(const std::string& n) { m_clusterCrossingAssocNodeName = n; }
  void setTrkrClusterNodeName(const std::string& n) { m_trkrClusterNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setQAFileName(const std::string& n) { m_qaFileName = n; }
  void setWriteQA(bool v) { m_writeQA = v; }
  void setLooseWindow(double rdphi, double dz)
  {
    m_looseRdphiWindow = rdphi;
    m_looseDzWindow = dz;
  }
  void setInttDzWindow(double dz) { m_inttDzWindow = dz; }
  void setSeedThetaWindow(double theta) { m_seedThetaWindow = theta; }
  void setResidualSigma(double rdphi, double dz)
  {
    m_sigmaRdphi = rdphi;
    m_sigmaDz = dz;
  }
  void setMissingLayerPenalty(double v) { m_missingLayerPenalty = v; }
  void setMaxBranchesPerLayer(unsigned int v) { m_maxBranchesPerLayer = v; }
  void setMaxChains(unsigned int v) { m_maxChains = v; }
  void setMinSiliconClusters(unsigned int v) { m_minSiliconClusters = v; }
  void setTpcAssociationClusterCut(unsigned int v) { m_tpcAssociationClusterCut = v; }
  void setTpcAssociationPtCut(double v) { m_tpcAssociationPtCut = v; }
  void setTpcAssociationCuts(unsigned int nclusters, double pt)
  {
    m_tpcAssociationClusterCut = nclusters;
    m_tpcAssociationPtCut = pt;
  }
  void setUseSagittaPhiFit(bool v) { m_useSagittaPhiFit = v; }
  void setPhiThetaWindowSigma(double phi, double theta)
  {
    m_phiWindowSigma = phi;
    m_thetaWindowSigma = theta;
  }
  void setUseDynamicAssociation(bool v) { m_useDynamicResiduals = v; }
  void setAssociationCalibrationMode(bool v) { m_associationCalibrationMode = v; }
  void setAngularResidualSigma(double phi, double theta)
  {
    m_sigmaPhi = phi;
    m_sigmaTheta = theta;
  }
  void setUseDynamicResiduals(bool v) { m_useDynamicResiduals = v; }
  void setDynamicPhiMean(unsigned int layer, double offset, double slope)
  {
    if (layer >= m_dynamicPhiMeanOffset.size())
    {
      return;
    }
    m_dynamicPhiMeanOffset[layer] = offset;
    m_dynamicPhiMeanSlope[layer] = slope;
  }
  void setDynamicThetaMean(unsigned int layer, double offset, double slope)
  {
    if (layer >= m_dynamicThetaMeanOffset.size())
    {
      return;
    }
    m_dynamicThetaMeanOffset[layer] = offset;
    m_dynamicThetaMeanSlope[layer] = slope;
  }
  void setDynamicResidualMeanSlope(double phi, double theta)
  {
    m_dynamicPhiMeanSlope.fill(phi);
    m_dynamicThetaMeanSlope.fill(theta);
  }

 private:
  struct SpacePoint
  {
    TrkrDefs::cluskey key{TrkrDefs::CLUSKEYMAX};
    unsigned int layer{0};
    unsigned int ladder{0};
    unsigned int sensor{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double r{0.0};
    double phi{0.0};
  };

  struct TrajectoryState
  {
    double phi_intercept{0.0};
    double phi_slope{0.0};
    double phi_S{0.0};
    double phi_x0{0.0};
    double phi_invR{0.0};
    double phi_theta{0.0};
    double phi_bline{0.0};
    double z_intercept{0.0};
    double z_slope{0.0};
    double seed_x0{0.0};
    double seed_y0{0.0};
    double seed_z0{0.0};
    double seed_cx{0.0};
    double seed_cy{0.0};
    double seed_phi0{0.0};
    double seed_slope{0.0};
    double seed_q_over_r{0.0};
    bool use_tpc_seed{false};
    bool use_silicon_seed{false};
    bool phi_sagitta_ok{false};
    bool valid{false};
  };

  struct ChainHit
  {
    SpacePoint point;
    double pred_phi{0.0};
    double pred_z{0.0};
    double pred_theta{0.0};
    double pred_x{0.0};
    double pred_y{0.0};
    double dphi{0.0};
    double dtheta{0.0};
    double rdphi{0.0};
    double dz{0.0};
    double tpc_pred_phi{0.0};
    double tpc_pred_z{0.0};
    double tpc_pred_theta{0.0};
    double tpc_dphi{0.0};
    double tpc_dtheta{0.0};
    double si_ref_pred_phi{0.0};
    double si_ref_pred_z{0.0};
    double si_ref_dphi{0.0};
    double si_ref_dtheta{0.0};
    double dynamic_mean_phi{0.0};
    double dynamic_mean_theta{0.0};
    double sigma_phi{0.0};
    double sigma_theta{0.0};
    double sdphi{0.0};
    double sdtheta{0.0};
    double chi2{0.0};
  };

  struct Chain
  {
    std::vector<SpacePoint> fit_points;
    std::vector<ChainHit> hits;
    TrajectoryState tpc_reference;
    TrajectoryState si_reference;
    TrajectoryState state;
    double chi2{0.0};
    double ndf{0.0};
    double score{0.0};
    double delta_phi0{std::numeric_limits<double>::max()};
    double delta_z0{std::numeric_limits<double>::max()};
    double delta_theta0{std::numeric_limits<double>::max()};
    double signed_delta_phi0{std::numeric_limits<double>::quiet_NaN()};
    double signed_delta_z0{std::numeric_limits<double>::quiet_NaN()};
    double signed_delta_theta0{std::numeric_limits<double>::quiet_NaN()};
    double delta_eta0{std::numeric_limits<double>::max()};
    double dca_score{std::numeric_limits<double>::max()};
    unsigned int missing_mask{0};
    unsigned int n_missing{0};
    double previous_dphi{0.0};
    double previous_dtheta{0.0};
    double z_offset{0.0};
    bool has_z_offset{false};
    bool has_previous_residual{false};
    double pt{0.0};
    double charge{0.0};
  };


  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  void createQaObjects();
  void fillQaResiduals(const Tpc_PolyTrack& track, const Chain& chain, const ChainHit& hit);
  void fillQaCrossing(const Tpc_PolyTrack& track, const Chain& chain);
  std::vector<SpacePoint> collectSiliconClusters() const;
  const TpcCrossingDecision* findCrossingDecision(unsigned int source_assembled_track_id) const;
  bool findInttCrossing(const Chain& chain, short& crossing) const;
  bool getGlobalClusterPosition(TrkrDefs::cluskey key, TrkrCluster* cluster, SpacePoint& point) const;
  TrajectoryState fitTrajectory(const std::vector<SpacePoint>& points) const;
  TrajectoryState makeTpcReferenceTrajectory(const Tpc_PolyTrack& track) const;
  TrajectoryState makeTpcSeedTrajectory(const Tpc_PolyTrack& track) const;
  TrajectoryState makeSiliconSeedTrajectory(const TrajectoryState& tpc_state,
                                            const Tpc_PolyTrack& track,
                                            const SpacePoint& first_mvtx) const;
  TrajectoryState correctSiliconSeedWithTwoHits(const TrajectoryState& seed_state,
                                                const ChainHit& outer_hit,
                                                const ChainHit& inner_hit) const;
  void updateSiliconTrajectoryHalfResidual(TrajectoryState& state, const ChainHit& hit) const;
  void refitSiliconTrajectoryFromMvtx(Chain& chain) const;
  void refitSiliconTrajectoryFromHits(Chain& chain) const;
  const SpacePoint* findBestMvtxCandidate(const Chain& chain,
                                          const std::vector<SpacePoint>& silicon_points,
                                          const std::set<TrkrDefs::cluskey>& used_keys,
                                          unsigned int layer,
                                          ChainHit& best_hit,
                                          double& best_chi2) const;
  bool predictAtRadius(const TrajectoryState& state, double r,
                       double& pred_phi, double& pred_z,
                       double& pred_x, double& pred_y) const;
  std::vector<Chain> buildChains(const Tpc_PolyTrack& track,
                                 const std::vector<SpacePoint>& silicon_points);
  Chain extendWithHit(const Chain& chain, const ChainHit& hit) const;
  Chain extendMissing(const Chain& chain, unsigned int layer) const;
  const Chain* selectBestChain(const std::vector<Chain>& chains) const;
  Chain attachClosestInttClusters(const Tpc_PolyTrack& track, const Chain& mvtx_chain,
                                 const std::vector<SpacePoint>& silicon_points);
  void fillTrack(const Tpc_PolyTrack& tpc_track, const Chain& chain);
  bool passesTpcAssociationCuts(const Tpc_PolyTrack& tpc_track) const;
  double wrapPhi(double phi) const;
  double unwrapPhiNear(double phi, double reference) const;
  double predictSagittaPhi(double r, const TrajectoryState& state) const;
  double pointTheta(const SpacePoint& point) const;
  double trajectoryPhi0NearBeam(const TrajectoryState& state) const;
  double trajectoryZ0NearBeam(const TrajectoryState& state) const;
  double trajectoryTheta0NearBeam(const TrajectoryState& state) const;
  bool computeChainDcaMetrics(Chain& chain, const TrajectoryState& tpc_state) const;
  double dynamicMeanPhi(unsigned int layer, double previous_dphi, bool has_previous) const;
  double dynamicSigmaPhi(double pt) const;
  double dynamicMeanTheta(unsigned int layer, double previous_dtheta, bool has_previous) const;
  double dynamicSigmaTheta(double pt, double pred_theta) const;
  double dynamicDzWindow(double pt) const;
  unsigned int layerBit(unsigned int layer) const;

  std::string m_tpcTrackNodeName{"TPC_POLYTRACKS"};
  std::string m_tpcClusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_crossingDecisionNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_clusterCrossingAssocNodeName{"TRKR_CLUSTERCROSSINGASSOC"};
  std::string m_trkrClusterNodeName{"TRKR_CLUSTER"};
  std::string m_outputNodeName{"FULL_POLYTRACKS"};
  std::string m_qaFileName{"full_polytrack_matcher_qa.root"};

  Tpc_PolyTrackContainer* m_tpcTracks{nullptr};
  TpcCrossingDecisionContainer* m_crossingDecisions{nullptr};
  TrkrClusterCrossingAssoc* m_clusterCrossingAssoc{nullptr};
  TrkrClusterContainer* m_trkrClusters{nullptr};
  ActsGeometry* m_actsGeometry{nullptr};
  Full_PolyTrackContainer* m_fullTracks{nullptr};

  bool m_writeQA{true};
  double m_looseRdphiWindow{0.2};
  double m_looseDzWindow{0.5};
  double m_inttDzWindow{1.0};
  double m_seedThetaWindow{0.02};
  double m_sigmaRdphi{0.15};
  double m_sigmaDz{1.0};
  double m_sigmaPhi{0.015};
  double m_sigmaTheta{0.02};
  double m_phiWindowSigma{3.0};
  double m_thetaWindowSigma{3.0};
  double m_missingLayerPenalty{6.0};
  bool m_useSagittaPhiFit{true};
  bool m_useDynamicResiduals{true};
  bool m_associationCalibrationMode{true};
  std::array<double, 7> m_dynamicPhiMeanOffset{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  std::array<double, 7> m_dynamicPhiMeanSlope{{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
  std::array<double, 7> m_dynamicThetaMeanOffset{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
  std::array<double, 7> m_dynamicThetaMeanSlope{{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}};
  unsigned int m_maxBranchesPerLayer{8};
  unsigned int m_maxChains{256};
  unsigned int m_minSiliconClusters{0};
  unsigned int m_tpcAssociationClusterCut{20};
  double m_tpcAssociationPtCut{0.1};
  std::vector<unsigned int> m_matchLayers{2, 1, 0};
  std::vector<unsigned int> m_inttMatchLayers{3, 4, 5, 6};
  std::vector<unsigned int> m_siliconSearchLayers{6, 5, 4, 3, 2, 1, 0};
  unsigned int m_event{0};

  TFile* m_qaFile{nullptr};
  unsigned int m_qaDphiBins{32};
  unsigned int m_qaPtBins{24};
  double m_qaDphiMin{-0.08};
  double m_qaDphiMax{0.08};
  double m_qaDthetaMin{-0.08};
  double m_qaDthetaMax{0.08};
  double m_qaPtMin{0.0};
  double m_qaPtMax{5.0};
  std::map<std::string, TH3F*> m_hResidualCurrentVsPreviousVsPt;
  std::map<std::string, TH3F*> m_hResidualSiVsTpc;
  std::map<std::string, TH1F*> m_hStandardResidual;
  std::map<std::string, TProfile*> m_hChainSignedResidualVsPt;
  TH2F* m_hInttCrossingVsTpcCrossing{nullptr};
  TH1F* m_hDeltaCrossing{nullptr};
  TH1F* m_hCrossingStatus{nullptr};
  TH1F* m_hChainDcaScore{nullptr};
  TH1F* m_hChainDeltaEta0{nullptr};
};

#endif