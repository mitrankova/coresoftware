// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCSILICONCROSSINGREFINER_H
#define TPCTRACKRECO_TPCSILICONCROSSINGREFINER_H

#include <fun4all/SubsysReco.h>

#include <string>

class ActsGeometry;
class Full_PolyTrack;
class Full_PolyTrackContainer;
class PHCompositeNode;
class TpcCrossingDecisionContainer;
class Tpc_PolyTrack;
class Tpc_PolyTrackContainer;

/** Refine the TPC crossing from spatial silicon--TPC z compatibility.
 *
 * INTT positions already attached to Full_PolyTrack may contribute to the
 * silicon spatial fit. No INTT timing or crossing association is read.
 */
class TpcSiliconCrossingRefiner : public SubsysReco
{
 public:
  explicit TpcSiliconCrossingRefiner(const std::string& name = "TpcSiliconCrossingRefiner");
  ~TpcSiliconCrossingRefiner() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setFullTrackNodeName(const std::string& n) { m_fullTrackNodeName = n; }
  void setTpcTrackNodeName(const std::string& n) { m_tpcTrackNodeName = n; }
  void setInputCrossingNodeName(const std::string& n) { m_inputCrossingNodeName = n; }
  void setOutputCrossingNodeName(const std::string& n) { m_outputCrossingNodeName = n; }
  void setCrossingPeriodNs(double v) { m_crossingPeriodNs = v; }
  void setMaximumAbsDeltaZ(double v) { m_maximumAbsDeltaZ = v; }
  void setCrossingFinderScoreWeight(double v) { m_crossingFinderScoreWeight = v; }

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  const Tpc_PolyTrack* findTpcTrack(unsigned int track_id) const;
  bool fitSiliconZ0(const Full_PolyTrack& track, double& z0) const;

  std::string m_fullTrackNodeName{"FULL_POLYTRACKS"};
  std::string m_tpcTrackNodeName{"TPC_POLYTRACKS"};
  std::string m_inputCrossingNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_outputCrossingNodeName{"TPC_SILICON_CROSSING_DECISIONS"};
  Full_PolyTrackContainer* m_fullTracks{nullptr};
  Tpc_PolyTrackContainer* m_tpcTracks{nullptr};
  TpcCrossingDecisionContainer* m_inputCrossings{nullptr};
  TpcCrossingDecisionContainer* m_outputCrossings{nullptr};
  ActsGeometry* m_geometry{nullptr};
  double m_crossingPeriodNs{106.56};
  double m_maximumAbsDeltaZ{5.0};
  double m_crossingFinderScoreWeight{1.0};
};

#endif
