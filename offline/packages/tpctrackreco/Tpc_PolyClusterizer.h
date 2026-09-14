// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCPOLYCLUSTERIZER_H
#define TPCTRACKRECO_TPCPOLYCLUSTERIZER_H

#include "TpcDriftPolylineLookup.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>

class IdealPadMap;
class PHCompositeNode;
class PHG4TpcGeomContainer;
class Tpc_AssembledTrackContainer;
class Tpc_PolyClusterContainer;
class TpcCrossingDecisionContainer;
class TrkrHitSetContainer;

class Tpc_PolyClusterizer : public SubsysReco
{
 public:
  explicit Tpc_PolyClusterizer(const std::string& name = "Tpc_PolyClusterizer");
  ~Tpc_PolyClusterizer() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setCrossingDecisionNodeName(const std::string& n) { m_crossingDecisionNodeName = n; }
  void setMaxAcceptedTier(unsigned char v) { m_maxAcceptedTier = v; }

  void setT0(double v) { m_driftConfig.t0 = v; m_driftConfig.t0Override = true; }
  void setTpcAdcClock(double v) { m_driftConfig.tpcAdcClock = v; m_driftConfig.tpcAdcClockOverride = true; }
  void setCrossingPeriodNs(double v) { m_driftConfig.crossingPeriodNs = v; m_driftConfig.crossingPeriodNsOverride = true; }
  void setReverseDriftStepNs(double v) { m_driftConfig.reverseDriftStepNs = v; m_driftConfig.reverseDriftStepNsOverride = true; }
  void setKEffSide0(double v) { m_driftConfig.kEffSide0 = v; m_driftConfig.kEffSide0Override = true; }
  void setKEffSide1(double v) { m_driftConfig.kEffSide1 = v; m_driftConfig.kEffSide1Override = true; }
  void setField3DCoefficientFile(const std::string& n) { m_driftConfig.field3DCoefficientFile = n; m_driftConfig.field3DCoefficientFileOverride = true; }
  void setElectricFieldMap(const std::string& n) { m_driftConfig.electricFieldMap = n; m_driftConfig.electricFieldMapOverride = true; }
  void setElectricFieldMap3DSide0(const std::string& n) { m_driftConfig.field3DSide0 = n; m_driftConfig.field3DSide0Override = true; }
  void setElectricFieldMap3DSide1(const std::string& n) { m_driftConfig.field3DSide1 = n; m_driftConfig.field3DSide1Override = true; }
  void setFrameElectricFieldMap3DSide0(const std::string& n) { m_driftConfig.framesSide0 = n; m_driftConfig.framesSide0Override = true; }
  void setFrameElectricFieldMap3DSide1(const std::string& n) { m_driftConfig.framesSide1 = n; m_driftConfig.framesSide1Override = true; }
  void setCMVoltageDefault(double v) { m_driftConfig.cmVoltageDefault = v; m_driftConfig.cmVoltageDefaultOverride = true; }
  void setUseSurveyGeometry(bool v) { m_driftConfig.useSurveyGeometry = v; m_driftConfig.useSurveyGeometryOverride = true; }
  void setMoveTpc(double x, double y, double z) { m_driftConfig.tpcMove = {{x, y, z}}; m_driftConfig.tpcMoveOverride = true; }
  void setRotateTpc(unsigned int index, double x, double y, double z)
  {
    if (index < m_driftConfig.tpcRotations.size())
    {
      m_driftConfig.tpcRotations[index] = {{x, y, z}};
      m_driftConfig.tpcRotationOverride[index] = true;
    }
  }
  void setStartZ(double south_z, double north_z)
  {
    m_driftConfig.startZSouth = south_z;
    m_driftConfig.startZNorth = north_z;
    m_driftConfig.startZOverride = true;
  }
  void setFrameChargeScale(double v) { m_driftConfig.frameChargeScale = v; m_driftConfig.frameChargeScaleOverride = true; }
  void setFieldCageVoltageOffsets(double ifcSouth, double ifcNorth, double ofcSouth, double ofcNorth)
  {
    m_driftConfig.fieldCageVoltageOffsets = {{ifcSouth, ifcNorth, ofcSouth, ofcNorth}};
    m_driftConfig.fieldCageVoltageOffsetsOverride = true;
  }
  void setUse2DElectricFieldMap(bool v) { m_driftConfig.use2DElectricFieldMap = v; m_driftConfig.use2DElectricFieldMapOverride = true; }

 private:
  struct Point
  {
    TrkrDefs::hitsetkey hitsetkey{0};
    TrkrDefs::hitkey hitkey{0};
    unsigned int layer{0};
    unsigned int side{0};
    unsigned int pad{0};
    unsigned int tbin{0};
    double adc{0.0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  struct Centroid
  {
    bool ok{false};
    unsigned int layer{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double rms_x{0.0};
    double rms_y{0.0};
    double rms_z{0.0};
  };

  struct ClusterParameters
  {
    double adc{0.0};
    unsigned int phi_width{0};
    unsigned int time_width{0};
    double phase{0.0};
  };

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool make_xyz_point(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, short crossing, Point& p) const;
  ClusterParameters make_cluster_parameters(const std::vector<Point>& points, const Centroid& centroid, int side) const;
  static Centroid make_centroid(const std::vector<Point>& points);

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  std::string m_crossingDecisionNodeName{"TPC_CROSSING_DECISIONS"};
  unsigned char m_maxAcceptedTier{1};

  Tpc_AssembledTrackContainer* m_assembledTracks{nullptr};
  Tpc_PolyClusterContainer* m_clusters{nullptr};
  TpcCrossingDecisionContainer* m_crossingDecisions{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};
  TpcDriftPolylineLookup* m_driftLookup{nullptr};
  const IdealPadMap* m_idealPadMap{nullptr};
  PHG4TpcGeomContainer* m_geomContainerTpc{nullptr};

  TpcDriftPolylineLookup::Config m_driftConfig;
  unsigned int m_event{0};
  double m_t0{8};
  double m_tpcAdcClock{56.881262};
  double m_crossingPeriodNs {106.56};
  double m_reverseDriftStepNs{56.881262};


  //! starting z position for primary electron backward drift
  /**
   * quoted values must be kept consistent with _max_driftlength + _CM_halfwidth
   * as defined in offline/packages/trackbase/ActsGeometry.h
   */
  double m_startZSouth{-102.605};
  double m_startZNorth{102.605};
  double m_kEffSide0{0.0};
  double m_kEffSide1{0.0};
  double m_cmVoltageDefault{375.0};
  bool use_survey_geometry = true;
  bool m_kEffSide0Override{false};
  bool m_kEffSide1Override{false};
  bool m_field3DCoefficientFileOverride{false};
  bool m_electricFieldMapOverride{false};
  bool m_field3DSide0Override{false};
  bool m_field3DSide1Override{false};
  bool m_framesSide0Override{false};
  bool m_framesSide1Override{false};
  bool m_use2DElectricFieldMap{false};
  std::array<double, 3> m_tpcMove{{0.0, 0.0, 0.0}};                                             //{{-0.16775, -0.0337, -0.71365}};
  std::array<std::array<double, 3>, 2> m_tpcRotations{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};  //{{{{0.0, 0.01485 / 10.0, 0.0}}, {{0.0298 / 8.0, 0.0, 0.0}}}};
  double m_frameChargeScale{-180.0};
  std::array<double, 4> m_fieldCageVoltageOffsets{{211.0, 0.0, 0.0, 0.0}};
};
#endif
