// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCCROSSINGTRAJECTORYBUILDER_H
#define TPCTRACKRECO_TPCCROSSINGTRAJECTORYBUILDER_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class PHCompositeNode;
class TpcCrossingDecisionContainer;
class TpcCrossingTrajectoryContainer;
class Tpc_PolyTrackContainer;

class TpcCrossingTrajectoryBuilder : public SubsysReco
{
 public:
  explicit TpcCrossingTrajectoryBuilder(const std::string& name = "TpcCrossingTrajectoryBuilder");
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setTrackNodeName(const std::string& value) { m_trackNodeName = value; }
  void setDecisionNodeName(const std::string& value) { m_decisionNodeName = value; }
  void setOutputNodeName(const std::string& value) { m_outputNodeName = value; }

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  std::string m_trackNodeName{"TPC_POLYTRACKS"};
  std::string m_decisionNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_outputNodeName{"TPC_CROSSING_TRAJECTORIES"};
  Tpc_PolyTrackContainer* m_tracks{nullptr};
  TpcCrossingDecisionContainer* m_decisions{nullptr};
  TpcCrossingTrajectoryContainer* m_trajectories{nullptr};
  std::array<float, 7> m_siliconRadii{{2.5F, 3.5F, 4.5F, 7.2F, 8.0F, 9.0F, 10.0F}};
};

#endif
