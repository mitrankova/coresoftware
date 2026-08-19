// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_FULLPOLYTRACKMATCHER_H
#define TPCTRACKRECO_FULLPOLYTRACKMATCHER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <map>
#include <string>
#include <vector>

class ActsGeometry;
class Full_PolyTrackContainer;
class PHCompositeNode;
class TFile;
class TH2F;
class TTree;
class TpcCrossingDecision;
class TpcCrossingDecisionContainer;
class Tpc_PolyCluster;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrack;
class Tpc_PolyTrackContainer;
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
  void setTrkrClusterNodeName(const std::string& n) { m_trkrClusterNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setQAFileName(const std::string& n) { m_qaFileName = n; }
  void setWriteQA(bool v) { m_writeQA = v; }
  void setLooseWindow(double rdphi, double dz)
  {
    m_looseRdphiWindow = rdphi;
    m_looseDzWindow = dz;
  }
  void setResidualSigma(double rdphi, double dz)
  {
    m_sigmaRdphi = rdphi;
    m_sigmaDz = dz;
  }
  void setMissingLayerPenalty(double v) { m_missingLayerPenalty = v; }
  void setMaxBranchesPerLayer(unsigned int v) { m_maxBranchesPerLayer = v; }
  void setMaxChains(unsigned int v) { m_maxChains = v; }
  void setMinSiliconClusters(unsigned int v) { m_minSiliconClusters = v; }

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
    double z_intercept{0.0};
    double z_slope{0.0};
    bool valid{false};
  };

  struct ChainHit
  {
    SpacePoint point;
    double pred_phi{0.0};
    double pred_z{0.0};
    double pred_x{0.0};
    double pred_y{0.0};
    double rdphi{0.0};
    double dz{0.0};
    double chi2{0.0};
  };

  struct Chain
  {
    std::vector<SpacePoint> fit_points;
    std::vector<ChainHit> hits;
    TrajectoryState state;
    double chi2{0.0};
    double ndf{0.0};
    double score{0.0};
    unsigned int missing_mask{0};
    unsigned int n_missing{0};
  };

  struct QaRow
  {
    unsigned int event{0};
    unsigned int tpc_track_id{0};
    unsigned int source_assembled_track_id{0};
    int crossing{0};
    unsigned int layer{0};
    unsigned int ladder{0};
    unsigned int sensor{0};
    unsigned long long cluster_key{0};
    float x{0.0F};
    float y{0.0F};
    float z{0.0F};
    float r{0.0F};
    float phi{0.0F};
    float pred_x{0.0F};
    float pred_y{0.0F};
    float pred_z{0.0F};
    float pred_phi{0.0F};
    float rdphi{0.0F};
    float dz{0.0F};
    float chi2{0.0F};
    float chain_chi2{0.0F};
    float chain_ndf{0.0F};
    float chain_score{0.0F};
    unsigned int chain_silicon_clusters{0};
    unsigned int chain_missing_layers{0};
    unsigned int chain_missing_mask{0};
    int accepted{0};
    int selected{0};
  };

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  void createQaObjects();
  void fillQaRow(const QaRow&);
  void fillQaForChain(unsigned int tpc_track_id, unsigned int source_id, int crossing,
                      const Chain& chain, bool selected);
  void fillQaForUnmatchedSeed(const Tpc_PolyTrack& tpc_track, int crossing, const Chain& chain);
  std::vector<SpacePoint> collectSiliconClusters() const;
  std::map<unsigned int, std::vector<const Tpc_PolyCluster*>> collectTpcClustersByTrack() const;
  const TpcCrossingDecision* findCrossingDecision(unsigned int source_assembled_track_id) const;
  bool getGlobalClusterPosition(TrkrDefs::cluskey key, TrkrCluster* cluster, SpacePoint& point) const;
  TrajectoryState fitTrajectory(const std::vector<SpacePoint>& points) const;
  bool predictAtRadius(const TrajectoryState& state, double r,
                       double& pred_phi, double& pred_z,
                       double& pred_x, double& pred_y) const;
  std::vector<Chain> buildChains(const Tpc_PolyTrack& track,
                                 const std::vector<const Tpc_PolyCluster*>& tpc_clusters,
                                 const std::vector<SpacePoint>& silicon_points);
  Chain extendWithHit(const Chain& chain, const SpacePoint& point,
                      double pred_phi, double pred_z,
                      double pred_x, double pred_y,
                      double rdphi, double dz, double chi2) const;
  Chain extendMissing(const Chain& chain, unsigned int layer) const;
  const Chain* selectBestChain(const std::vector<Chain>& chains) const;
  void fillTrack(const Tpc_PolyTrack& tpc_track, const Chain& chain);
  double wrapPhi(double phi) const;
  double unwrapPhiNear(double phi, double reference) const;
  unsigned int layerBit(unsigned int layer) const;

  std::string m_tpcTrackNodeName{"TPC_POLYTRACKS"};
  std::string m_tpcClusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_crossingDecisionNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_trkrClusterNodeName{"TRKR_CLUSTER"};
  std::string m_outputNodeName{"FULL_POLYTRACKS"};
  std::string m_qaFileName{"full_polytrack_matcher_qa.root"};

  Tpc_PolyTrackContainer* m_tpcTracks{nullptr};
  Tpc_PolyClusterContainer* m_tpcClusters{nullptr};
  TpcCrossingDecisionContainer* m_crossingDecisions{nullptr};
  TrkrClusterContainer* m_trkrClusters{nullptr};
  ActsGeometry* m_actsGeometry{nullptr};
  Full_PolyTrackContainer* m_fullTracks{nullptr};

  bool m_writeQA{true};
  double m_looseRdphiWindow{1.0};
  double m_looseDzWindow{5.0};
  double m_sigmaRdphi{0.15};
  double m_sigmaDz{1.0};
  double m_missingLayerPenalty{6.0};
  unsigned int m_maxBranchesPerLayer{8};
  unsigned int m_maxChains{256};
  unsigned int m_minSiliconClusters{1};
  std::vector<unsigned int> m_matchLayers{6, 5, 4, 3, 2, 1, 0};
  unsigned int m_event{0};

  TFile* m_qaFile{nullptr};
  TTree* m_candidateTree{nullptr};
  QaRow m_qaRow;
  std::map<unsigned int, TH2F*> m_hRdphiVsZ;
  std::map<unsigned int, TH2F*> m_hDzVsZ;
  std::map<unsigned int, TH2F*> m_hRdphiVsPhi;
  std::map<unsigned int, TH2F*> m_hDzVsPhi;
};

#endif
