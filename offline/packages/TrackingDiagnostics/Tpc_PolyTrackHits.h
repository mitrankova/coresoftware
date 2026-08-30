#ifndef TPC_POLYTRACKHITS_H
#define TPC_POLYTRACKHITS_H

#include <fun4all/SubsysReco.h>

#include <string>

class IdealPadMap;
class PHCompositeNode;
class TFile;
class TTree;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrackContainer;
class TrkrHitSetContainer;

class Tpc_PolyTrackHits : public SubsysReco
{
 public:
  explicit Tpc_PolyTrackHits(const std::string& name = "Tpc_PolyTrackHits",
                             const std::string& outfilename = "tpc_polytrack_hits.root");
  ~Tpc_PolyTrackHits() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setClusterNodeName(const std::string& n) { m_clusterNodeName = n; }
  void setTpc_PolyTrackNodeName(const std::string& n) { m_finalTrackNodeName = n; }
  void setHitSetNodeName(const std::string& n) { m_hitSetNodeName = n; }
  void setOutputFileName(const std::string& n) { m_outfilename = n; }
  void setMinPt(double v) { m_minPt = v; }
  void setMinTpcClusters(unsigned int v) { m_minTpcClusters = v; }

 private:
  bool get_nodes(PHCompositeNode* topNode);
  void reset_tree_values();

  std::string m_outfilename;
  std::string m_clusterNodeName;
  std::string m_finalTrackNodeName;
  std::string m_hitSetNodeName;

  double m_minPt{0.1};
  unsigned int m_minTpcClusters{10};

  unsigned int m_evt{0};
  TFile* m_outfile{nullptr};
  TTree* m_tree{nullptr};
  IdealPadMap* m_idealPadMap{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};
  Tpc_PolyTrackContainer* m_finalTracks{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};

  unsigned int m_event{0};
  unsigned int m_polyTrackId{0};
  unsigned int m_sourceAssembledTrackId{0};
  unsigned int m_polyClusterId{0};
  unsigned int m_ntpcClusters{0};
  double m_pt{0.0};
  unsigned int m_pad{0};
  unsigned int m_timebin{0};
  unsigned int m_layer{0};
  unsigned int m_sector{0};
  unsigned int m_side{0};
  double m_phi{0.0};
  double m_radius{0.0};
  double m_x{0.0};
  double m_y{0.0};
  double m_z{0.0};
  double m_adc{0.0};
};

#endif
