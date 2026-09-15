#ifndef TPCTRACKRECO_TPCCROSSINGTRACKFINALIZER_H
#define TPCTRACKRECO_TPCCROSSINGTRACKFINALIZER_H
#include "FastFieldTrackFitter.h"
#include <fun4all/SubsysReco.h>
#include <memory>
#include <string>
class Full_PolyTrackContainer;
class PHCompositeNode;
class TpcCrossingTrajectoryContainer;
class TpcSiliconMatchCandidateContainer;
class Tpc_PolyTrackContainer;
class Tpc_PolyClusterContainer;
class TpcDriftPolylineLookup;
class TrkrHitSetContainer;
class PHField;
class TpcCrossingTrackFinalizer : public SubsysReco
{
 public:
  explicit TpcCrossingTrackFinalizer(const std::string& name = "TpcCrossingTrackFinalizer");
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  void setTrackNodeName(const std::string& value) { m_trackNodeName = value; }
  void setTrajectoryNodeName(const std::string& value) { m_trajectoryNodeName = value; }
  void setCandidateNodeName(const std::string& value) { m_candidateNodeName = value; }
  void setOutputNodeName(const std::string& value) { m_outputNodeName = value; }
  void setClusterNodeName(const std::string& value) { m_clusterNodeName = value; }
  void setCorrectedClusterNodeName(const std::string& value) { m_correctedClusterNodeName = value; }
 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  std::string m_trackNodeName{"TPC_POLYTRACKS"};
  std::string m_trajectoryNodeName{"TPC_CROSSING_TRAJECTORIES"};
  std::string m_candidateNodeName{"TPC_SILICON_MATCH_CANDIDATES"};
  std::string m_outputNodeName{"FULL_POLYTRACKS"};
  std::string m_clusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_correctedClusterNodeName{"TPC_POLYCLUSTERS_CROSSING_CORRECTED"};
  Tpc_PolyTrackContainer* m_tracks{nullptr};
  TpcCrossingTrajectoryContainer* m_trajectories{nullptr};
  TpcSiliconMatchCandidateContainer* m_candidates{nullptr};
  Full_PolyTrackContainer* m_output{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};
  Tpc_PolyClusterContainer* m_correctedClusters{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};
  TpcDriftPolylineLookup* m_lookup{nullptr};
  const PHField* m_field{nullptr};
  std::unique_ptr<FastFieldTrackFitter> m_fitter;
  unsigned int m_event{0};
};
#endif
