// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCPOLYTRIPLETRECO_H
#define TPCTRACKRECO_TPCPOLYTRIPLETRECO_H

// Fast TPC seeds for the silicon seeding, run before the full polyline clusterizer.
//
// Crossing hypotheses per assembled track:
//   triggered mode (default): one crossing, setTriggeredCrossing (default 0), no further checks
//   streaming mode          : every INTT crossing of the event (INTT time buckets in TRKR_CLUSTER),
//                             no silicon vertex needed; a crossing is kept if
//                               - the track's earliest / latest hit times fit in the drift range,
//                               - the three reference points are inside the TPC on the track's side,
//                               - the trajectory z at the beam is within setMaxVertexZ of setCollisionZ
// One Tpc_CrossingTrajectory per accepted (assembled track, crossing); it keeps its assembled
// track id (get_source_assembled_track_id), its crossing (get_crossing) and its three points.
//
// Three points along the track, for each crossing:
//   beginning : innermost usable layer
//   end       : outermost usable layer
//   middle    : usable layer closest to the radial midpoint of the two
// usable layer = not the first or last layer of a module, and no hit of the track on that layer
// within setEdgePadMargin pads of a sector edge (truncated clusters would bias the circle).
// Tracks whose usable layers span less than setMinLayerSpan layers get no points.
//   reference point = ADC-weighted mean pad phi and time bin of the track's hits on that layer,
//                     drifted back once through the PHGarfield polylines with the crossing's
//                     time offset (drift lookup built per layer on first use)
// Track model:
//   xy: circle through the three points, curvature from the sagitta (kappa = 2 h / (a b)),
//       so a straight track is simply kappa = 0 (no infinite radius anywhere)
//   rz: straight line z = z0 + slope * r
//
// Output: TPC_CROSSINGTRAJECTORIES (Tpc_CrossingTrajectoryContainer). Tpc_PolyCluster / Tpc_PolyTrack
// stay for the full clusterizer and fit, which can run later.
// SiHitSeedReco reads it with setTpcTrajectoryNodeName("TPC_CROSSINGTRAJECTORIES") and
// setMinTpcClusters(3).

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>

class ActsGeometry;
class IdealPadMap;
class PHCompositeNode;
class PHGarfield;
class PHG4TpcGeomContainer;
class Tpc_AssembledTrackContainer;
class Tpc_CrossingTrajectoryContainer;
class Tpc_CrossingTrajectoryv1;
class TrkrClusterContainer;
class TrkrHitSetContainer;

class Tpc_PolyTripletReco : public SubsysReco
{
 public:
  explicit Tpc_PolyTripletReco(const std::string& name = "Tpc_PolyTripletReco");
  ~Tpc_PolyTripletReco() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  // crossing hypotheses (see above)
  void setTriggeredMode(bool v) { m_triggeredMode = v; }
  void setTriggeredCrossing(short v) { m_triggeredCrossing = v; }
  // streaming mode: |z at the beam - collision z| window [cm]
  void setMaxVertexZ(double v) { m_maxVertexZ = v; }
  void setCollisionZ(double v) { m_collisionZ = v; }

  // pads from a sector edge whose hits make a layer unusable (1 = only the edge pads)
  void setEdgePadMargin(unsigned int v) { m_edgePadMargin = v; }
  // minimum (end - beginning) layer distance for a track to get points
  void setMinLayerSpan(unsigned int v) { m_minLayerSpan = v; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
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
  bool insideTpc(const RefPoint& p, unsigned int side) const;
  void fitTrajectory(const std::array<RefPoint, 3>& pts, Tpc_CrossingTrajectoryv1& out) const;
  std::set<short> inttCrossings() const;
  Polyline& polyline(unsigned int layer, unsigned int side, unsigned int sector, unsigned int sample)
  {
    return m_lookup[(((layer - FirstLayer) * 2 + side) * 12 + sector) * NSamples + sample];
  }
  const Polyline& polyline(unsigned int layer, unsigned int side, unsigned int sector, unsigned int sample) const
  {
    return m_lookup[(((layer - FirstLayer) * 2 + side) * 12 + sector) * NSamples + sample];
  }

  bool m_triggeredMode{true};
  short m_triggeredCrossing{0};
  double m_maxVertexZ{20.0};
  double m_collisionZ{0.0};
  // TPC volume for the streaming-mode checks (as TpcCrossingFinder) [cm]
  double m_tpcHalfLength{105.5};
  double m_radialTolerance{2.0};
  double m_zTolerance{1.0};
  double m_centralMembraneTolerance{1.0};
  double m_edgeToleranceTbins{1.0};

  unsigned int m_edgePadMargin{1};
  std::array<bool, NLayers> m_moduleEdgeLayer{};
  unsigned int m_minLayerSpan{16};
  std::string m_inputNodeName{"TPC_ASSEMBLEDTRACKS"};
  std::string m_outputNodeName{"TPC_CROSSINGTRAJECTORIES"};
  double m_magneticFieldTesla{1.4};
  double m_minCurvature{1.0e-7};

  double m_t0{8};
  bool m_t0Override{false};
  double m_tpcAdcClock{56.881262};
  bool m_tpcAdcClockOverride{false};
  double m_crossingPeriodNs{106.56};
  double m_driftStepNs{56.881262};
  unsigned int m_driftSubsteps{10};
  double m_maxLookupTimeNs{0.0};  // longest drift polyline built so far
  std::array<double, 2> m_t0Side{{0.0, 0.0}};
  std::array<double, 2> m_startZ{{-102.605, 102.605}};

  std::vector<Polyline> m_lookup;
  std::array<int, NLayers> m_layerState{};  // 0 not built, 1 built, -1 failed

  // counters for the End summary
  unsigned long m_nTracks{0};
  unsigned long m_nShortSpan{0};
  unsigned long m_nTracksWithTrajectory{0};
  unsigned long m_nHypotheses{0};
  unsigned long m_nTrajectories{0};
  unsigned long m_nNoPoints{0};
  unsigned long m_nTimeWindow{0};
  unsigned long m_nOutsideTpc{0};
  unsigned long m_nVertexZ{0};

  Tpc_AssembledTrackContainer* m_assembledTracks{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};
  TrkrClusterContainer* m_clusterMap{nullptr};
  PHG4TpcGeomContainer* m_geomContainer{nullptr};
  PHGarfield* m_garfield{nullptr};
  ActsGeometry* m_geometry{nullptr};
  IdealPadMap* m_padMap{nullptr};
  Tpc_CrossingTrajectoryContainer* m_trajectories{nullptr};
  unsigned int m_event{0};
};
#endif
