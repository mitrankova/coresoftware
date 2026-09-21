#ifndef TRACKINGDIAGNOSTICS_FULLPOLYTRACKRESIDUALS_H
#define TRACKINGDIAGNOSTICS_FULLPOLYTRACKRESIDUALS_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>
#include <vector>

class ActsGeometry;
class Full_PolyTrackContainer;
class PHCompositeNode;
class PHField;
class TFile;
class TTree;
class Tpc_PolyClusterContainer;
class TrkrClusterContainer;

class Full_PolyTrackResiduals : public SubsysReco
{
 public:
  explicit Full_PolyTrackResiduals(const std::string& name = "Full_PolyTrackResiduals",
                                   const std::string& outfilename = "full_polytrack_residuals.root");
  ~Full_PolyTrackResiduals() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setFullTrackNodeName(const std::string& value) { m_fullTrackNodeName = value; }
  void setTpcClusterNodeName(const std::string& value) { m_tpcClusterNodeName = value; }
  void setTrkrClusterNodeName(const std::string& value) { m_trkrClusterNodeName = value; }
  void setOutputFileName(const std::string& value) { m_outfilename = value; }
  void setMinPt(double value) { m_minPt = value; }
  void setMaxPt(double value) { m_maxPt = value; }
  void setPropagationStep(double value) { m_propagationStep = value; }

 private:
  bool getNodes(PHCompositeNode*);
  void resetTreeValues();

  std::string m_outfilename;
  std::string m_fullTrackNodeName{"FULL_POLYTRACKS"};
  std::string m_tpcClusterNodeName{"TPC_POLYCLUSTERS_CROSSING_CORRECTED"};
  std::string m_trkrClusterNodeName{"TRKR_CLUSTER"};
  double m_minPt{0.0};
  double m_maxPt{1.0e30};
  double m_propagationStep{0.25};
  unsigned int m_evt{0};

  TFile* m_outfile{nullptr};
  TTree* m_tree{nullptr};
  Full_PolyTrackContainer* m_fullTracks{nullptr};
  Tpc_PolyClusterContainer* m_tpcClusters{nullptr};
  TrkrClusterContainer* m_trkrClusters{nullptr};
  ActsGeometry* m_actsGeometry{nullptr};
  const PHField* m_field{nullptr};

  unsigned int m_event{0};
  unsigned int m_fullTrackId{0};
  unsigned int m_parentTpcTrackId{0};
  unsigned int m_sourceAssembledTrackId{0};
  short m_crossing{0};
  unsigned int m_ntpcClusters{0};
  unsigned int m_ninttClusters{0};
  unsigned int m_nmvtxClusters{0};
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
  std::vector<std::uint64_t> m_clusterKey;
  std::vector<unsigned int> m_detector;
  std::vector<unsigned int> m_layer;
  std::vector<double> m_clusterX;
  std::vector<double> m_clusterY;
  std::vector<double> m_clusterZ;
  std::vector<double> m_clusterR;
  std::vector<double> m_clusterPhi;
  std::vector<double> m_stateX;
  std::vector<double> m_stateY;
  std::vector<double> m_stateZ;
  std::vector<double> m_stateR;
  std::vector<double> m_statePhi;
  std::vector<double> m_deltaPhi;
  std::vector<double> m_residualRPhi;
  std::vector<double> m_residualZ;
};

#endif
