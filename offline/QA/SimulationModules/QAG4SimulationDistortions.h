// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONDISTORTIONS_H
#define QAG4SIMULATIONDISTORTIONS_H

#include <fun4all/SubsysReco.h>
#include <tpc/TpcGlobalPositionWrapper.h>
#include <trackbase/TrkrDefs.h>

#include <math.h>
#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class TrkrClusterContainer;
class SvtxTrack;
class ActsGeometry;

class QAG4SimulationDistortions : public SubsysReco
{
 public:
  QAG4SimulationDistortions(const std::string& name = "QAG4SimulationDistortions");

  ~QAG4SimulationDistortions() override;

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode *topNode) override;

  //! track map name
  void set_trackmap_name( const std::string& value )
  { m_trackmapname = value; }

 private:

  //! track map name
  std::string m_trackmapname = "SvtxSiliconMMTrackMap";

  std::string get_histo_prefix()
  {
    return std::string("h_") + Name() + std::string("_");
  }

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);
  bool checkTrack(SvtxTrack* track);
  SvtxTrackMap* m_trackMap = nullptr;
  TrkrClusterContainer* m_clusterContainer = nullptr;
  ActsGeometry* m_tGeometry = nullptr;

  //! tpc global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  int m_event = 0;
  float m_tanAlpha = NAN;
  float m_tanBeta = NAN;
  float m_drphi = NAN;
  float m_dz = NAN;
  float m_clusR = NAN;
  float m_clusPhi = NAN;
  float m_clusZ = NAN;
  float m_statePhi = NAN;
  float m_stateZ = NAN;
  float m_stateR = NAN;
  float m_stateRPhiErr = NAN;
  float m_stateZErr = NAN;
  float m_clusRPhiErr = NAN;
  float m_clusZErr = NAN;
  TrkrDefs::cluskey m_cluskey = TrkrDefs::CLUSKEYMAX;

  ///@name counters
  //@{
  int m_total_tracks = 0;
  int m_accepted_tracks = 0;

  int m_total_states = 0;
  int m_accepted_states = 0;
  //@}


};

#endif  // QAG4SIMULATIONDISTORTIONS_H
