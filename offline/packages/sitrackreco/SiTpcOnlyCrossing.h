// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCONLYCROSSING_H
#define SITRACKRECO_SITPCONLYCROSSING_H

// SiTpcOnlyCrossing: crossing (collision) of the TPC tracks WITHOUT a Si chain, decided with the logic of
// TpcCrossingFinder against the collisions found by the Si seeding.
//
// TpcCrossingResolver decides the crossing of a track without Si chain from the z0 of its 3-point crossing
// trajectories only (tier 1 if exactly one collision matches, else tier 2 at the crossing closest to z = 0).
// Here the standard TpcCrossingFinder (tpctrackreco, unchanged) does it instead, once the collisions are known.
// Two instances of this module around a TpcCrossingFinder (after SiSeedingOrchestrator, before
// Tpc_PolyClusterizer):
//
//   SiSeedingOrchestrator             TPC_CROSSING_DECISIONS, SI_HIT_SEED_EVENT (selected chains)
//   SiTpcOnlyCrossing(MakeVertices)   collisions of the Si-matched tracks -> SvtxVertexMap (setVertexMapName):
//                                     per crossing, the Si-matched tracks (selected chain, decision tier <=
//                                     setMaxCollisionTier, finite Si z) are clustered in z (setClusterZWindow);
//                                     clusters with >= setMinCollisionTracks tracks are vertices.  Vertex z in the
//                                     TPC frame of TpcCrossingFinder: the mean TPC z0 of its tracks (setVertexZ:
//                                     TpcZ0Mean, default) or their mean Si z (SiZ).
//   TpcCrossingFinder                 every assembled track, every INTT crossing: drift-time window, containment,
//                                     TPC z0 from its own Garfield points, closest vertex of the crossing ->
//                                     setOutputNodeName (e.g. TPC_CROSSING_DECISIONS_TPCONLY)
//   SiTpcOnlyCrossing(MergeDecisions) for every track WITHOUT a selected Si chain the decision of TPC_CROSSING_DECISIONS
//                                     is replaced by the TpcCrossingFinder one: crossing, z0, vertex z, dz, status.
//                                     Tier: TpcCrossingFinder tier 0 (|dz| <= its setMaxVertexDz) -> setVertexTier
//                                     (default 1: fitted as TPC-only by SiTpcTrackSeedFiller with tpcOnlyMaxTier 1);
//                                     with setRejectAmbiguous a close rival (SelectedByVertexAmbiguous) -> tier 2;
//                                     tier 1 / 2 of the finder (vertex too far / no vertex) -> tier 2.
//                                     No valid crossing at all: tier 2, crossing unchanged.
//                                     Tracks with a Si chain are never touched.
//                                     setAddMissing: also tracks without any resolver decision (no crossing
//                                     trajectory, e.g. span < Tpc_PolyTripletReco::setMinLayerSpan).
#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class PHCompositeNode;

class SiTpcOnlyCrossing : public SubsysReco
{
 public:
  enum class Mode
  {
    MakeVertices,
    MergeDecisions
  };
  enum class VertexZ
  {
    TpcZ0Mean,  // mean decision tpc_z0 of the collision's tracks (TPC frame, same drift as TpcCrossingFinder)
    SiZ         // mean Si z at the beam (tpc_z0 + delta_z)
  };

  explicit SiTpcOnlyCrossing(Mode mode, const std::string& name = "SiTpcOnlyCrossing");
  ~SiTpcOnlyCrossing() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // ---- nodes
  void setDecisionNodeName(const std::string& s) { m_decisionName = s; }     // resolver decisions (modified by Merge)
  void setFinderDecisionNodeName(const std::string& s) { m_finderName = s; }  // TpcCrossingFinder output
  void setSiSeedNodeName(const std::string& s) { m_siSeedName = s; }
  void setAssembledTrackNodeName(const std::string& s) { m_assembledName = s; }
  void setVertexMapName(const std::string& s) { m_vertexMapName = s; }        // MakeVertices output = finder input

  // ---- MakeVertices
  void setMaxCollisionTier(unsigned int t) { m_maxCollisionTier = t; }
  void setMinCollisionTracks(unsigned int n) { m_minCollisionTracks = n; }
  void setClusterZWindow(double cm) { m_clusterZWindow = cm; }
  void setVertexZ(VertexZ v) { m_vertexZ = v; }
  void setTrackZSigma(double cm) { m_trackZSigma = cm; }  // per-track z resolution for the vertex error

  // ---- MergeDecisions
  void setVertexTier(unsigned int t) { m_vertexTier = t; }
  void setRejectAmbiguous(bool b) { m_rejectAmbiguous = b; }
  void setAddMissing(bool b) { m_addMissing = b; }

 private:
  int makeVertices(PHCompositeNode* topNode);
  int mergeDecisions(PHCompositeNode* topNode);

  Mode m_mode;
  std::string m_decisionName{"TPC_CROSSING_DECISIONS"};
  std::string m_finderName{"TPC_CROSSING_DECISIONS_TPCONLY"};
  std::string m_siSeedName{"SI_HIT_SEED_EVENT"};
  std::string m_assembledName{"TPC_ASSEMBLEDTRACKS"};
  std::string m_vertexMapName{"SiSeedCollisionVertexMap"};

  unsigned int m_maxCollisionTier{1};
  unsigned int m_minCollisionTracks{2};
  double m_clusterZWindow{0.7};
  VertexZ m_vertexZ{VertexZ::TpcZ0Mean};
  double m_trackZSigma{0.5};

  unsigned int m_vertexTier{1};
  bool m_rejectAmbiguous{false};
  bool m_addMissing{false};

  // counters
  unsigned long m_nEvents{0};
  unsigned long m_nVertices{0};
  unsigned long m_nVertexTracks{0};
  double m_sumSiMinusTpc{0};      // mean Si z - TPC z0 of the vertex tracks (End: information)
  unsigned long m_nSiMinusTpc{0};
  unsigned long m_nTpcOnly{0};          // tracks without Si chain looked at
  unsigned long m_nNoFinderDecision{0};
  unsigned long m_nAdded{0};
  unsigned long m_nChangedCrossing{0};  // crossing differs from the resolver one
  std::map<std::pair<int, int>, unsigned long> m_tierChange;  // (resolver tier, new tier) -> count
};

#endif
