#ifndef TPCTRACKRECO_TPCCROSSINGTRAJECTORYBUILDER_H
#define TPCTRACKRECO_TPCCROSSINGTRAJECTORYBUILDER_H
#include "FastFieldTrackFitter.h"
#include <fun4all/SubsysReco.h>
#include <array>
#include <memory>
#include <string>
class PHCompositeNode;
class PHField;
class TpcCrossingDecisionContainer;
class TpcCrossingTrajectory;
class TpcCrossingTrajectoryContainer;
class TpcDriftPolylineLookup;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrackContainer;
class TrkrHitSetContainer;
class TpcCrossingTrajectoryBuilder : public SubsysReco
{
 public:
  explicit TpcCrossingTrajectoryBuilder(const std::string& name = "TpcCrossingTrajectoryBuilder");
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  void setTrackNodeName(const std::string& value) { m_trackNodeName = value; }
  void setClusterNodeName(const std::string& value) { m_clusterNodeName = value; }
  void setDecisionNodeName(const std::string& value) { m_decisionNodeName = value; }
  void setOutputNodeName(const std::string& value) { m_outputNodeName = value; }
  void setValidationFraction(double value) { m_validationFraction = value; }
 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool addSiliconStates(TpcCrossingTrajectory&, const FastFieldTrackFitter::Result&) const;
  std::string m_trackNodeName{"TPC_POLYTRACKS"};
  std::string m_clusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_decisionNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_outputNodeName{"TPC_CROSSING_TRAJECTORIES"};
  Tpc_PolyTrackContainer* m_tracks{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};
  TpcCrossingDecisionContainer* m_decisions{nullptr};
  TpcCrossingTrajectoryContainer* m_trajectories{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};
  TpcDriftPolylineLookup* m_lookup{nullptr};
  const PHField* m_field{nullptr};
  std::unique_ptr<FastFieldTrackFitter> m_fitter;
  std::array<float, 7> m_siliconRadii{{2.5F, 3.5F, 4.5F, 7.2F, 8.0F, 9.0F, 10.0F}};
  double m_validationFraction{0.0};
};
#endif
