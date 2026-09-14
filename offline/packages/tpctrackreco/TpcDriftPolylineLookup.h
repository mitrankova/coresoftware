#ifndef TPCTRACKRECO_TPCDRIFTPOLYLINELOOKUP_H
#define TPCTRACKRECO_TPCDRIFTPOLYLINELOOKUP_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class IdealPadMap;
class PHCompositeNode;
class PHGarfield;
class TpcDriftPolylineLookupInit;

/**
 * Run-scoped, non-persistent TPC reverse-drift lookup.
 *
 * The lookup is stored once on RUN/TPC_DRIFT_LOOKUP as a transient PHDataNode.
 * Its PHGarfield setup and sampling reproduce the pre-refactor
 * Tpc_PolyClusterizer implementation. Configuration is supplied exactly once by
 * TpcDriftPolylineLookupInit during InitRun; consumers are read-only afterward.
 */
class TpcDriftPolylineLookup
{
 public:
  static constexpr const char* NodeName = "TPC_DRIFT_LOOKUP";
  static constexpr unsigned int NPhiSamples = 24;

  struct Point
  {
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  struct Config
  {
    double t0{8.0};
    double tpcAdcClock{56.881262};
    double crossingPeriodNs{106.56};
    double reverseDriftStepNs{56.881262};
    double startZSouth{-102.605};
    double startZNorth{102.605};
    double kEffSide0{0.0};
    double kEffSide1{0.0};
    double cmVoltageDefault{375.0};
    double frameChargeScale{-180.0};
    bool useSurveyGeometry{false};
    bool use2DElectricFieldMap{false};

    std::array<double, 3> tpcMove{{0.0, 0.0, 0.0}};
    std::array<std::array<double, 3>, 2> tpcRotations{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
    std::array<double, 4> fieldCageVoltageOffsets{{211.0, 0.0, 0.0, 0.0}};

    std::string electricFieldMap;
    std::string field3DCoefficientFile;
    std::string field3DSide0;
    std::string field3DSide1;
    std::string framesSide0;
    std::string framesSide1;

    bool t0Override{false};
    bool tpcAdcClockOverride{false};
    bool crossingPeriodNsOverride{false};
    bool reverseDriftStepNsOverride{false};
    bool startZOverride{false};
    bool kEffSide0Override{false};
    bool kEffSide1Override{false};
    bool cmVoltageDefaultOverride{false};
    bool frameChargeScaleOverride{false};
    bool useSurveyGeometryOverride{false};
    bool use2DElectricFieldMapOverride{false};
    bool tpcMoveOverride{false};
    std::array<bool, 2> tpcRotationOverride{{false, false}};
    bool fieldCageVoltageOffsetsOverride{false};
    bool electricFieldMapOverride{false};
    bool field3DCoefficientFileOverride{false};
    bool field3DSide0Override{false};
    bool field3DSide1Override{false};
    bool framesSide0Override{false};
    bool framesSide1Override{false};
  };

  TpcDriftPolylineLookup();
  ~TpcDriftPolylineLookup();

  TpcDriftPolylineLookup(const TpcDriftPolylineLookup&) = delete;
  TpcDriftPolylineLookup& operator=(const TpcDriftPolylineLookup&) = delete;

  static TpcDriftPolylineLookup* get(PHCompositeNode* topNode);
  bool isInitialized() const { return m_initialized; }

  bool getPosition(unsigned int layer,
                   unsigned int side,
                   unsigned int pad,
                   unsigned int tbin,
                   short crossing,
                   Point& point) const;

  bool getPosition(TrkrDefs::hitsetkey hitsetkey,
                   TrkrDefs::hitkey hitkey,
                   short crossing,
                   Point& point) const;

  bool getDelta(unsigned int layer,
                unsigned int side,
                unsigned int pad,
                unsigned int tbin,
                short crossing,
                short referenceCrossing,
                Point& delta) const;

  bool getDelta(TrkrDefs::hitsetkey hitsetkey,
                TrkrDefs::hitkey hitkey,
                short crossing,
                short referenceCrossing,
                Point& delta) const;

  bool hasEntry(unsigned int layer,
                unsigned int side,
                unsigned int pad,
                unsigned int tbin,
                short crossing) const;

  double crossingToDriftTimeNs(unsigned int tbin, short crossing) const;
  double maxLookupTimeNs() const { return m_maxLookupTimeNs; }
  double t0() const { return m_config.t0; }
  double tpcAdcClock() const { return m_config.tpcAdcClock; }
  double crossingPeriodNs() const { return m_config.crossingPeriodNs; }
  double reverseDriftStepNs() const { return m_config.reverseDriftStepNs; }

  const IdealPadMap* idealPadMap() const { return m_idealPadMap.get(); }

 private:
  friend class TpcDriftPolylineLookupInit;

  static TpcDriftPolylineLookup* getOrCreate(PHCompositeNode* topNode, int verbosity = 0);
  bool configure(const Config& config, int verbosity = 0);
  bool initialize(PHCompositeNode* topNode, int verbosity = 0);

  static constexpr unsigned int FirstLayer = 7;
  static constexpr unsigned int LastLayer = 54;
  static constexpr unsigned int NLayers = LastLayer - FirstLayer + 1;
  static constexpr unsigned int NSides = 2;
  static constexpr unsigned int NSectors = 12;
  static constexpr unsigned int NPolylines = NLayers * NSides * NSectors * NPhiSamples;

  struct DriftPoint
  {
    float delta_r{0.0F};
    float delta_phi{0.0F};
    float z{0.0F};
  };

  struct DriftPolyline
  {
    double phi{0.0};
    std::uint32_t offset{0};
    std::uint32_t count{0};
  };

  bool load_cdb_inputs(int verbosity);
  void configure_garfield(PHGarfield* garfield) const;
  bool build_drift_lookup(int verbosity);
  bool sample_drift_lookup(unsigned int layer,
                           unsigned int side,
                           unsigned int pad,
                           unsigned int tbin,
                           short crossing,
                           Point& point) const;

  const DriftPoint& driftPoint(const DriftPolyline& polyline, std::size_t index) const;
  static unsigned int drift_lookup_index(unsigned int layer_index,
                                         unsigned int side,
                                         unsigned int sector,
                                         unsigned int sample);

  Config m_config;
  std::unique_ptr<IdealPadMap> m_idealPadMap;
  std::unique_ptr<PHGarfield> m_garfield;
  std::array<DriftPolyline, NPolylines> m_driftLookup{};
  std::vector<DriftPoint> m_driftPoints;

  std::string m_electricFieldMap;
  std::string m_field3DSide0;
  std::string m_field3DSide1;
  std::string m_framesSide0;
  std::string m_framesSide1;
  double m_kEffSide0{0.0};
  double m_kEffSide1{0.0};
  double m_startZSouth{-102.605};
  double m_startZNorth{102.605};
  std::array<double, 3> m_tpcMove{{0.0, 0.0, 0.0}};
  std::array<std::array<double, 3>, 2> m_tpcRotations{{{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}}};
  double m_maxLookupTimeNs{0.0};
  bool m_configured{false};
  bool m_initialized{false};
};

/**
 * Fun4All InitRun-only configurator for the shared TPC drift lookup.
 *
 * Register this subsystem before TpcCrossingFinder and Tpc_PolyClusterizer.
 * All PHGarfield/lookup settings are made here exactly once.
 */
class TpcDriftPolylineLookupInit : public SubsysReco
{
 public:
  explicit TpcDriftPolylineLookupInit(const std::string& name = "TpcDriftPolylineLookupInit");
  ~TpcDriftPolylineLookupInit() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setT0(double v) { m_config.t0 = v; m_config.t0Override = true; }
  void setTpcAdcClock(double v) { m_config.tpcAdcClock = v; m_config.tpcAdcClockOverride = true; }
  void setCrossingPeriodNs(double v) { m_config.crossingPeriodNs = v; m_config.crossingPeriodNsOverride = true; }
  void setReverseDriftStepNs(double v) { m_config.reverseDriftStepNs = v; m_config.reverseDriftStepNsOverride = true; }
  void setKEffSide0(double v) { m_config.kEffSide0 = v; m_config.kEffSide0Override = true; }
  void setKEffSide1(double v) { m_config.kEffSide1 = v; m_config.kEffSide1Override = true; }
  void setField3DCoefficientFile(const std::string& n) { m_config.field3DCoefficientFile = n; m_config.field3DCoefficientFileOverride = true; }
  void setElectricFieldMap(const std::string& n) { m_config.electricFieldMap = n; m_config.electricFieldMapOverride = true; }
  void setElectricFieldMap3DSide0(const std::string& n) { m_config.field3DSide0 = n; m_config.field3DSide0Override = true; }
  void setElectricFieldMap3DSide1(const std::string& n) { m_config.field3DSide1 = n; m_config.field3DSide1Override = true; }
  void setFrameElectricFieldMap3DSide0(const std::string& n) { m_config.framesSide0 = n; m_config.framesSide0Override = true; }
  void setFrameElectricFieldMap3DSide1(const std::string& n) { m_config.framesSide1 = n; m_config.framesSide1Override = true; }
  void setCMVoltageDefault(double v) { m_config.cmVoltageDefault = v; m_config.cmVoltageDefaultOverride = true; }
  void setUseSurveyGeometry(bool v) { m_config.useSurveyGeometry = v; m_config.useSurveyGeometryOverride = true; }
  void setMoveTpc(double x, double y, double z) { m_config.tpcMove = {{x, y, z}}; m_config.tpcMoveOverride = true; }
  void setRotateTpc(unsigned int index, double x, double y, double z)
  {
    if (index < m_config.tpcRotations.size())
    {
      m_config.tpcRotations[index] = {{x, y, z}};
      m_config.tpcRotationOverride[index] = true;
    }
  }
  void setStartZ(double south_z, double north_z)
  {
    m_config.startZSouth = south_z;
    m_config.startZNorth = north_z;
    m_config.startZOverride = true;
  }
  void setFrameChargeScale(double v) { m_config.frameChargeScale = v; m_config.frameChargeScaleOverride = true; }
  void setFieldCageVoltageOffsets(double ifcSouth, double ifcNorth, double ofcSouth, double ofcNorth)
  {
    m_config.fieldCageVoltageOffsets = {{ifcSouth, ifcNorth, ofcSouth, ofcNorth}};
    m_config.fieldCageVoltageOffsetsOverride = true;
  }
  void setUse2DElectricFieldMap(bool v) { m_config.use2DElectricFieldMap = v; m_config.use2DElectricFieldMapOverride = true; }

 private:
  TpcDriftPolylineLookup::Config m_config;
};

#endif  // TPCTRACKRECO_TPCDRIFTPOLYLINELOOKUP_H
