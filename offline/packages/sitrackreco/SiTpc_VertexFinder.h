// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCVERTEXFINDER_H
#define SITRACKRECO_SITPCVERTEXFINDER_H

// Collision vertices from one track source, per bunch crossing, with PairVertexFinder (the
// algorithm of Tpc_PolyTrackVertexFinder: track pairs as PHSimpleVertexFinder, outlier-pair
// removal, weighted vertex fit, absorption of unused tracks).  Register it once per source:
//
//   Source::Si    SI_TRAJECTORY (Si_Trajectory, SiTrajectoryFitter)   -> SiPolyVertexMap
//   Source::Tpc   TPC_POLYTRACKS (Tpc_PolyTrack, Tpc_PolyTrackReco)   -> TpcPolyVertexMap
//   Source::SiTpc SITPC_TRACKS (SiTpc_Track, SiTpcTrackMatcher)       -> SiTpcPolyVertexMap
//
// Output: SvtxVertexMap_v1 of SvtxVertex_v3 (x, y, z, 3x3 error, chi2, ndof, beam crossing);
// the track ids in a vertex are the ids of the source container:
//   Si: Si_Trajectory::get_id   Tpc: Tpc_PolyTrack::get_track_id   SiTpc: SiTpc_Track::get_id
//
// Crossing of a track (setGroupByCrossing, default on: vertices are found per crossing):
//   TPC_CROSSING_DECISIONS (TpcCrossingResolver) of its assembled TPC track
//     Si   : SI_HIT_SEED_EVENT chain (get_chain_id) -> tpc_assembled_track_id; without a
//            decision the chain's tpc_crossing
//     Tpc  : Tpc_PolyTrack::get_source_assembled_track_id
//     SiTpc: SiTpc_Track::get_tpc_assembled_track_id
//   only decisions with tier <= setMaxDecisionTier (default 1); tracks without a usable
//   crossing are skipped with setRequireCrossing (default on), otherwise grouped together with
//   beam crossing = SHRT_MAX.
// Frames: Si and SiTpc in the TPC frame of SiTpcBeamAlignment (beam at (0, 0) with the
// alignment); Tpc in the frame of the poly tracks (set the TPC beam spot with setBeamSpot).
// z of the Tpc vertices is the TPC z at the decided crossing (Tpc_PolyClusterizer drift).

#include "PairVertexFinder.h"

#include <fun4all/SubsysReco.h>

#include <cstddef>
#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class SvtxVertexMap;
class TpcCrossingDecisionContainer;

class SiTpc_VertexFinder : public SubsysReco
{
 public:
  enum class Source
  {
    Si,
    Tpc,
    SiTpc
  };

  //! defaults (node names, resolutions, cuts) are set for the source
  explicit SiTpc_VertexFinder(const std::string& name = "SiTpc_VertexFinder", Source source = Source::SiTpc);
  ~SiTpc_VertexFinder() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  // ---- nodes
  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setVertexMapName(const std::string& n) { m_vertexMapName = n; }
  void setCrossingDecisionNodeName(const std::string& n) { m_decisionNodeName = n; }
  void setSiSeedNodeName(const std::string& n) { m_siSeedNodeName = n; }

  // ---- crossings
  void setGroupByCrossing(bool v) { m_groupByCrossing = v; }
  void setRequireCrossing(bool v) { m_requireCrossing = v; }
  void setMaxDecisionTier(int t) { m_maxDecisionTier = t; }

  // ---- track selection
  void setTrackPtCut(double v) { m_minPt = v; }               // tracks with a measured pt
  void setMinTpcClusters(unsigned int v) { m_minTpcClusters = v; }  // Tpc, SiTpc
  void setMinSiPoints(unsigned int v) { m_minSiPoints = v; }        // Si, SiTpc
  void setMinMvtxPoints(unsigned int v) { m_minMvtx = v; }          // Si, SiTpc (MVTX points of the track)
  //    SiTpc: |phi_Si - phi_TPC| of the Si-TPC match below this (wrong chains / secondaries have large
  //    values; <= 0: no cut)
  void setMaxMatchDphi(double rad) { m_maxMatchDphi = rad; }
  void setTrackQualityCut(double chi2ndf) { m_maxChi2Ndf = chi2ndf; }  // Tpc fit chi2/ndf
  void setUseStraightLines(bool v) { m_useLines = v; }             // fits without curvature
  void setMagneticFieldTesla(double v) { m_bz = v; }                // Tpc helix radius from pt
  //    resolution of straight-line tracks (Si fits without curvature: the line misses the beam by
  //    ~ r1 r2 / 2R, ~0.1 cm at pt 0.2 GeV for MVTX-only chains) [cm]
  void setLineTrackResolution(double sigmaRPhi, double sigmaZ)
  {
    m_lineSigmaRPhi = sigmaRPhi;
    m_lineSigmaZ = sigmaZ;
  }
  //    SiTpc: track position from its silicon points, curvature (radius, helicity) and tanl from
  //    the Si+TPC fit (default on: the Si+TPC circle may be fitted with TPC points only, see
  //    SiTpcTrackMatcher::setXYWeights, which gives TPC precision at the beam)
  void setSiTpcPositionFromSi(bool v) { m_siTpcPositionFromSi = v; }
  //    shift added to every output vertex [cm] (common frame for the comparison of sources, e.g.
  //    TPC: minus the TPC beam line in x, y and the Si - TPC z offset in z)
  void setVertexShiftCm(double dx, double dy, double dz)
  {
    m_shift[0] = dx;
    m_shift[1] = dy;
    m_shift[2] = dz;
  }

  // ---- beam
  void setBeamSpot(double x, double y)
  {
    m_finder.config().beamX = x;
    m_finder.config().beamY = y;
  }
  void setBeamSpotSigma(double sx, double sy)
  {
    m_finder.config().beamSigmaX = sx;
    m_finder.config().beamSigmaY = sy;
  }
  void setUseBeamSpotConstraint(bool v) { m_finder.config().useBeamSpotConstraint = v; }
  void setMaxDcaXY(double v) { m_finder.config().maxDcaXY = v; }
  void setMaxAbsZ(double v) { m_finder.config().maxAbsZ = v; }

  // ---- pairs (PHSimpleVertexFinder-like)
  void setDcaCut(double v) { m_finder.config().pairDcaCut = v; }
  void setOutlierPairCut(double v) { m_finder.config().outlierPairCut = v; }
  void setRetryLooseDcaCut(bool v) { m_finder.config().retryLooseDcaCut = v; }
  void setMinTracksPerVertex(unsigned int v) { m_finder.config().minTracksPerVertex = v; }

  // ---- fit
  void setTrackResolution(double sigmaRPhi, double sigmaZ)
  {
    m_finder.config().sigmaRPhi = sigmaRPhi;
    m_finder.config().sigmaZ = sigmaZ;
  }
  void setMultipleScatteringTerm(double k) { m_finder.config().msTerm = k; }
  void setUseTrackCovariance(bool v) { m_finder.config().useTrackCovariance = v; }  // Tpc only
  void setMaxTrackChi2(double v) { m_finder.config().maxTrackChi2 = v; }
  void setMaxIterations(unsigned int v) { m_finder.config().maxIterations = v; }
  void setAbsorbUnusedTracks(bool v) { m_finder.config().absorbUnusedTracks = v; }

  PairVertexFinder::Config& finderConfig() { return m_finder.config(); }

 private:
  static constexpr short kNoCrossing = 32767;

  bool createNodes(PHCompositeNode* topNode);
  void collectSi(PHCompositeNode* topNode);
  void collectTpc(PHCompositeNode* topNode);
  void collectSiTpc(PHCompositeNode* topNode);
  // crossing of an assembled TPC track from the decisions; false if none usable
  bool crossingOf(unsigned int assembledId, short& crossing) const;
  void add(PairVertexFinder::TrackModel m, bool hasCrossing, short crossing);

  Source m_source;
  std::string m_inputNodeName;
  std::string m_vertexMapName;
  std::string m_decisionNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_siSeedNodeName{"SI_HIT_SEED_EVENT"};

  bool m_groupByCrossing{true};
  bool m_requireCrossing{true};
  int m_maxDecisionTier{1};

  double m_minPt{0.0};
  unsigned int m_minTpcClusters{10};
  unsigned int m_minSiPoints{3};
  unsigned int m_minMvtx{2};
  double m_maxMatchDphi{-1.0};
  unsigned long m_nRejMvtx{0};
  unsigned long m_nRejMatch{0};
  double m_maxChi2Ndf{1.0e9};
  bool m_useLines{true};
  double m_bz{1.4};
  double m_lineSigmaRPhi{0.1};
  double m_lineSigmaZ{0.1};
  bool m_siTpcPositionFromSi{true};
  double m_shift[3]{0.0, 0.0, 0.0};

  PairVertexFinder m_finder;
  SvtxVertexMap* m_vertexMap{nullptr};
  TpcCrossingDecisionContainer* m_decisions{nullptr};

  // this event: crossing -> track models
  std::map<short, std::vector<PairVertexFinder::TrackModel>> m_groups;

  // statistics
  unsigned long m_nEvents{0};
  unsigned long m_nInput{0};
  unsigned long m_nSelected{0};
  unsigned long m_nNoCrossing{0};
  unsigned long m_nSiPosition{0};
  unsigned long m_nGroups{0};
  unsigned long m_nVertices{0};
  unsigned long m_nTracksInVertices{0};
  std::map<std::size_t, unsigned long> m_multiplicity;  // tracks per vertex -> vertices
  std::vector<double> m_sigmaZ;                         // vertex z errors (End: median)
};

#endif
