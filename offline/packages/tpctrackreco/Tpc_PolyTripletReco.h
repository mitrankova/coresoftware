// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCPOLYTRIPLETRECO_H
#define TPCTRACKRECO_TPCPOLYTRIPLETRECO_H

// Fast TPC seeds for the silicon seeding, run before the full polyline clusterizer.
//
// For every assembled track (with an accepted crossing decision) three points along the track:
//   beginning : innermost usable layer
//   end       : outermost usable layer
//   middle    : usable layer closest to the radial midpoint of the two
// usable layer = not the first or last layer of a module, and no hit of the track on that layer
// within setEdgePadMargin pads of a sector edge (truncated clusters would bias the circle).
// Tracks whose usable layers span less than setMinLayerSpan layers get no points.
//   reference point = ADC-weighted mean pad phi and time bin of the track's hits on that layer,
//                     drifted back once through the PHGarfield polylines
//                     (drift lookup built per layer on first use)
// Track model:
//   xy: circle through the three points, curvature from the sagitta (kappa = 2 h / (a b)),
//       so a straight track is simply kappa = 0 (no infinite radius anywhere)
//   rz: straight line z = z0 + slope * r
//
// Output (same classes as the full chain, own node names, so the full chain can run later):
//   TPC_TRIPLETPOINTS : Tpc_PolyClusterContainer, one reference point per track and layer
//   TPC_TRIPLETTRACKS : Tpc_PolyTrackContainer
// SiHitSeedReco reads them with setTpcClusterNodeName / setTpcTrackNodeName / setMinTpcClusters(3).

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <map>
#include <string>
#include <vector>

class ActsGeometry;
class IdealPadMap;
class PHCompositeNode;
class PHGarfield;
class PHG4TpcGeomContainer;
class Tpc_AssembledTrackContainer;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrackContainer;
class TpcCrossingDecisionContainer;
class TrkrHitSetContainer;

class Tpc_PolyTripletReco : public SubsysReco
{
 public:
  explicit Tpc_PolyTripletReco(const std::string& name = "Tpc_PolyTripletReco");
  ~Tpc_PolyTripletReco() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  // pads from a sector edge whose hits make a layer unusable (1 = only the edge pads)
  void setEdgePadMargin(unsigned int v) { m_edgePadMargin = v; }
  // minimum (end - beginning) layer distance for a track to get points
  void setMinLayerSpan(unsigned int v) { m_minLayerSpan = v; }
  void setPointNodeName(const std::string& n) { m_pointNodeName = n; }
  void setTrackNodeName(const std::string& n) { m_trackNodeName = n; }
  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setCrossingDecisionNodeName(const std::string& n) { m_crossingDecisionNodeName = n; }
  void setMaxAcceptedTier(unsigned char v) { m_maxAcceptedTier = v; }
  void setMagneticFieldTesla(double v) { m_magneticFieldTesla = v; }
  // |kappa| below this [1/cm] is stored as a straight line (charge 0)
  void setMinCurvature(double v) { m_minCurvature = v; }
  // same meaning and defaults as in Tpc_PolyClusterizer
  void setT0(double v) { m_t0 = v; m_t0Override = true; }
  void setTpcAdcClock(double v) { m_tpcAdcClock = v; m_tpcAdcClockOverride = true; }
  void setCrossingPeriodNs(double v) { m_crossingPeriodNs = v; }

 private:
  static constexpr unsigned int FirstLayer = 7;
  static constexpr unsigned int NLayers = 48;
  static constexpr unsigned int NSamples = 24;  // drift polylines per sector

  struct DriftPoint
  {
    float dr{0.0F};
    float dphi{0.0F};
    float z{0.0F};
  };
  struct Polyline
  {
    double phi{0.0};
    std::vector<DriftPoint> points;
  };
  struct Hit
  {
    TrkrDefs::hitsetkey hitsetkey{0};
    TrkrDefs::hitkey hitkey{0};
    unsigned int pad{0};
    double tbin{0};
    double adc{0};
  };
  struct RefPoint
  {
    unsigned int layer{0};
    double x{0}, y{0}, z{0};
  };

  bool buildLayer(unsigned int layer);  // drift polylines of one layer, both sides
  bool drift(unsigned int layer, unsigned int side, unsigned int sector, double phi, double timeNs,
             double& x, double& y, double& z) const;
  bool makePoint(unsigned int layer, unsigned int side, short crossing, const std::vector<Hit>& hits, RefPoint& p);
  int storePoint(unsigned int assembledId, unsigned int side, const RefPoint& p, const std::vector<Hit>& hits);
  void fillTrack(unsigned int assembledId, const std::array<int, 3>& pointIndex);
  Polyline& polyline(unsigned int layer, unsigned int side, unsigned int sector, unsigned int sample)
  {
    return m_lookup[(((layer - FirstLayer) * 2 + side) * 12 + sector) * NSamples + sample];
  }
  const Polyline& polyline(unsigned int layer, unsigned int side, unsigned int sector, unsigned int sample) const
  {
    return m_lookup[(((layer - FirstLayer) * 2 + side) * 12 + sector) * NSamples + sample];
  }

  unsigned int m_edgePadMargin{1};
  std::array<bool, 48> m_moduleEdgeLayer{};
  unsigned int m_minLayerSpan{16};
  std::string m_inputNodeName{"TPC_ASSEMBLEDTRACKS"};
  std::string m_crossingDecisionNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_pointNodeName{"TPC_TRIPLETPOINTS"};
  std::string m_trackNodeName{"TPC_TRIPLETTRACKS"};
  unsigned char m_maxAcceptedTier{1};
  double m_magneticFieldTesla{1.4};
  double m_minCurvature{1.0e-7};

  double m_t0{8};
  bool m_t0Override{false};
  double m_tpcAdcClock{56.881262};
  bool m_tpcAdcClockOverride{false};
  double m_crossingPeriodNs{106.56};
  double m_driftStepNs{56.881262};
  unsigned int m_driftSubsteps{10};
  std::array<double, 2> m_t0Side{{0.0, 0.0}};
  std::array<double, 2> m_startZ{{-102.605, 102.605}};

  std::vector<Polyline> m_lookup;
  std::array<int, NLayers> m_layerState{};  // 0 not built, 1 built, -1 failed

  // counters for the End summary
  unsigned long m_nTracks{0};
  unsigned long m_nTriplets{0};
  unsigned long m_nShortSpan{0};
  unsigned long m_nNoMiddle{0};
  unsigned long m_nDriftFailed{0};

  std::map<TrkrDefs::hitsetkey, unsigned int> m_nextIndex;  // cluster key index per hitset (per event)

  Tpc_AssembledTrackContainer* m_assembledTracks{nullptr};
  TpcCrossingDecisionContainer* m_crossingDecisions{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};
  PHG4TpcGeomContainer* m_geomContainer{nullptr};
  PHGarfield* m_garfield{nullptr};
  ActsGeometry* m_geometry{nullptr};
  IdealPadMap* m_padMap{nullptr};
  Tpc_PolyClusterContainer* m_points{nullptr};
  Tpc_PolyTrackContainer* m_tracks{nullptr};
  unsigned int m_event{0};
};
#endif