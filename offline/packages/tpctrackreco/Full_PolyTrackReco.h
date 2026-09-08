// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_FULLPOLYTRACKRECO_H
#define TPCTRACKRECO_FULLPOLYTRACKRECO_H
#include "Tpc_FittingTools.h"
#include <fun4all/SubsysReco.h>
#include <string>
class ActsGeometry;
class Full_PolyTrack;
class Full_PolyTrackContainer;
class PHCompositeNode;
class TrkrClusterContainer;
/** Refit matched silicon clusters while retaining TPC-seed curvature. */
class Full_PolyTrackReco : public SubsysReco
{
 public:
  explicit Full_PolyTrackReco(const std::string& name = "Full_PolyTrackReco");
  ~Full_PolyTrackReco() override = default;
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  void setInputNodeName(const std::string& name) { m_inputNodeName = name; }
  void setOutputNodeName(const std::string& name) { m_outputNodeName = name; }
 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool fitSiliconClusters(const Full_PolyTrack&, Tpc_FittingTools::FitResult&) const;
  void fillTrack(const Full_PolyTrack&, const Tpc_FittingTools::FitResult&, bool);
  std::string m_inputNodeName{"FULL_POLYTRACKS"};
  std::string m_outputNodeName{"REFIT_FULL_POLYTRACKS"};
  std::string m_clusterNodeName{"TRKR_CLUSTER"};
  Full_PolyTrackContainer* m_inputTracks{nullptr};
  Full_PolyTrackContainer* m_outputTracks{nullptr};
  TrkrClusterContainer* m_clusters{nullptr};
  ActsGeometry* m_geometry{nullptr};
};
#endif
