#ifndef FULL_POLYCLUSTERRESIDUALS_H
#define FULL_POLYCLUSTERRESIDUALS_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class Full_PolyTrackContainer;
class Tpc_PolyTrackVertexContainer;
class PHCompositeNode;
class TFile;
class TTree;
class Tpc_PolyClusterContainer;

class Full_PolyClusterResiduals : public SubsysReco
{
 public:
  explicit Full_PolyClusterResiduals(const std::string& name = "Full_PolyClusterResiduals",
                                     const std::string& outfilename = "full_polycluster_residuals.root");
  ~Full_PolyClusterResiduals() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setClusterNodeName(const std::string& n) { m_clusterNodeName = n; }
  void setFull_PolyTrackNodeName(const std::string& n) { m_fullTrackNodeName = n; }
  void setTpc_PolyTrackVertexNodeName(const std::string& n) { m_trackVertexNodeName = n; }
  void setOutputFileName(const std::string& n) { m_outfilename = n; }
  void setMagneticFieldTesla(double b) { m_magneticFieldTesla = b; }
  void setMinPt(double v) { m_minPt = v; }
  void setMaxPt(double v) { m_maxPt = v; }
  void setMinTpcClusters(unsigned int v) { m_minTpcClusters = v; }
  void setMaxTpcClusters(unsigned int v) { m_maxTpcClusters = v; }
  void setUseStraightLineTracks(bool v) { m_useStraightLineTracks = v; }

 private:
  bool get_nodes(PHCompositeNode* topNode);
  void reset_tree_values();

  std::string m_outfilename;
  std::string m_clusterNodeName;
  std::string m_fullTrackNodeName;
  std::string m_trackVertexNodeName;

  double m_magneticFieldTesla{1.4};
  double m_minPt{0.0};
  double m_maxPt{1.0e30};
  unsigned int m_minTpcClusters{0};
  unsigned int m_maxTpcClusters{0xffffffffu};
  bool m_useStraightLineTracks{false};

  unsigned int m_evt{0};
  TFile* m_outfile{nullptr};
  TTree* m_tree{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};
  Full_PolyTrackContainer* m_fullTracks{nullptr};
  Tpc_PolyTrackVertexContainer* m_trackVertices{nullptr};

  unsigned int m_event{0};
  unsigned int m_fullTrackId{0};
  unsigned int m_tpcPolyTrackId{0};
  unsigned int m_sourceClusterId{0};
  unsigned int m_sourceAssembledTrackId{0};
  int m_side{0};
  unsigned int m_ntpcClusters{0};
  unsigned int m_nsiliconClusters{0};
  int m_fitStatus{0};
  double m_pt{0.0};
  double m_px{0.0};
  double m_py{0.0};
  double m_pz{0.0};
  double m_eta{0.0};
  double m_theta{0.0};
  double m_charge{0.0};
  double m_chi2{0.0};
  double m_ndf{0.0};
  double m_quality{0.0};
  double m_score{0.0};
  double m_vertexX{0.0};
  double m_vertexY{0.0};
  double m_vertexZ{0.0};
  double m_vertexR{0.0};
  double m_pcaX{0.0};
  double m_pcaY{0.0};
  double m_pcaZ{0.0};
  double m_zDCA{0.0};
  double m_rDCA{0.0};
  double m_rDCAZero{0.0};
  double m_R{0.0};
  double m_rzSlope{0.0};
  std::vector<unsigned int> m_clusterIndex;
  std::vector<unsigned int> m_sector;
  std::vector<unsigned int> m_layer;
  std::vector<double> m_clusterX;
  std::vector<double> m_clusterY;
  std::vector<double> m_clusterZ;
  std::vector<double> m_clusterR;
  std::vector<double> m_clusterPhi;
  std::vector<double> m_clusterAdc;
  std::vector<unsigned int> m_clusterPadSize;
  std::vector<double> m_stateX;
  std::vector<double> m_stateY;
  std::vector<double> m_stateZ;
  std::vector<double> m_stateZDca;
  std::vector<double> m_stateR;
  std::vector<double> m_statePhi;
  std::vector<double> m_deltaPhi;
  std::vector<double> m_residualRPhi;
  std::vector<double> m_residualZ;
};

#endif
