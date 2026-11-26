#ifndef TPC_POLYTRACKHITS_H
#define TPC_POLYTRACKHITS_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <limits>
#include <string>

class IdealPadMap;
class PHCompositeNode;
class TFile;
class TTree;
class TpcCrossingDecisionContainer;
class Tpc_PolyCluster;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrackContainer;
class TrkrHit;
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
  void setCrossingDecisionNodeName(const std::string& n) { m_crossingDecisionNodeName = n; }
  void setTpc_PolyTrackNodeName(const std::string& n) { m_finalTrackNodeName = n; }
  void setHitSetNodeName(const std::string& n) { m_hitSetNodeName = n; }
  void setOutputFileName(const std::string& n) { m_outfilename = n; }

  // Optional track selection. Defaults keep ALL found tracks.
  // Hits of tracks rejected by these cuts are not written at all
  // (they are NOT relabelled as "no track" hits).
  void setMinPt(double v) { m_minPt = v; }
  void setMinTpcClusters(unsigned int v) { m_minTpcClusters = v; }

  // Write TPC hits that are not on any track (track fields = short max). Default: on.
  void setSaveUnassignedHits(bool v) { m_saveUnassignedHits = v; }

  // Sentinel used for every quantity that is undefined for a given hit
  static constexpr short kUndefShort = std::numeric_limits<short>::max();
  static constexpr unsigned int kUndefUInt = static_cast<unsigned int>(std::numeric_limits<short>::max());
  static constexpr double kUndefDouble = static_cast<double>(std::numeric_limits<short>::max());

 private:
  bool get_nodes(PHCompositeNode* topNode);
  void reset_tree_values();
  // Fills hit-level fields (geometry + adc). Returns false if the hit has no valid pad geometry.
  bool fill_hit_fields(uint32_t hitsetkey, uint32_t hitkey, const TrkrHit* hit);
  // Fills cluster-level fields and the hit position stored in the cluster
  void fill_cluster_fields(const Tpc_PolyCluster* cluster, unsigned int ihit);

  static uint64_t make_hit_id(uint32_t hitsetkey, uint32_t hitkey)
  {
    return (static_cast<uint64_t>(hitsetkey) << 32U) | static_cast<uint64_t>(hitkey);
  }

  std::string m_outfilename;
  std::string m_clusterNodeName;
  std::string m_crossingDecisionNodeName;
  std::string m_finalTrackNodeName;
  std::string m_hitSetNodeName;

  double m_minPt{-1.0};
  unsigned int m_minTpcClusters{0};
  bool m_saveUnassignedHits{true};

  unsigned int m_evt{0};
  TFile* m_outfile{nullptr};
  TTree* m_tree{nullptr};
  IdealPadMap* m_idealPadMap{nullptr};
  TpcCrossingDecisionContainer* m_crossingDecisions{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};
  Tpc_PolyTrackContainer* m_finalTracks{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};

  unsigned int m_event{0};
  unsigned int m_polyTrackId{0};
  unsigned int m_sourceAssembledTrackId{0};
  short m_crossing{0};
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
  double m_polyClusterX{0.0};
  double m_polyClusterY{0.0};
  double m_polyClusterZ{0.0};
  double m_polyClusterRmsX{0.0};
  double m_polyClusterRmsY{0.0};
  double m_polyClusterRmsZ{0.0};
  double m_polyClusterRPhiError{0.0};
  double m_polyClusterZError{0.0};
  double m_x{0.0};
  double m_y{0.0};
  double m_z{0.0};
  double m_adc{0.0};
};

#endif