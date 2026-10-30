// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_TPCCROSSINGRESOLVER_H
#define SITRACKRECO_TPCCROSSINGRESOLVER_H

// Chooses one bunch crossing AND one Si chain per assembled TPC track and writes the crossing as
// TPC_CROSSING_DECISIONS (same container as TpcCrossingFinder), so Tpc_PolyClusterizer and
// everything after it run unchanged with the resolved crossing.
//
// Flow: Tpc_PolyTripletReco -> SiHitSeedReco -> TpcCrossingResolver -> SiTrajectoryFitter
//       (fits the chains selected here) -> Tpc_PolyClusterizer ... -> SiTpcTrackMatcher.
// SiHitSeedReco must keep all its candidate chains (setStoreSelectedOnly(false), the default).
//
// Triggered mode: every assembled track gets setTriggeredCrossing (default 0), tier 0
//                 (identical to TpcCrossingFinder in triggered mode); the chains selected by
//                 SiHitSeedReco are kept.
//
// Streaming mode (a crossing hypothesis = one crossing trajectory of Tpc_PolyTripletReco, i.e. a
// crossing with INTT hits; no other crossings are made here):
//   1. best chain  : every crossing trajectory takes its best Si chain (highest credit) among the
//                    chains SiHitSeedReco found for it with tpc_chi2 <= setMaxChainChi2
//   2. credit      : credit = match * quality, both in (0, 1]
//                    match   : how well the chain fits the trajectory at its crossing
//                              exp(-dz^2 / 2 sigma_dz^2) * exp(-w_rphi * chi2_rphi/n / 2)
//                              dz = Si - TPC z offset at this crossing (tpc_dz_offset)
//                    quality : how good the chain is by itself
//                              (n_layers / 7)^p * exp(-w_ddphi * <ddphi pull^2> / 2)
//                              * (setNoInttFactor if the chain has no INTT hit)
//   3. vote        : each hypothesis votes for the collision z (Si z at the beam = TPC z0 + dz)
//                    in its crossing with  vote = credit * credit / sum of the track's credits,
//                    so a track has at most one vote and an ambiguous one splits it
//   4. vertices    : per crossing, 1-D clustering of the vote-weighted z; one crossing = one
//                    collision: its main vertex is the cluster with the largest vote sum, the
//                    other clusters are the weak votes of wrong hypotheses
//   5. decision    : score = credit * (floor + W(-h) / (W(-h) + supportScale)
//                    * exp(-(z - z_main)^2 / 2 sigma^2)), sigma^2 = sigma_z^2 (1 + 1/n_other),
//                    sigma_z = setSiZSigma (single-track Si z), W(-h) = vote sum of the main vertex
//                    of the crossing without the hypothesis' own vote (a hypothesis in a weak
//                    cluster sees the main vertex far away: credit * floor);
//                    tracks in order of their best score take their best hypothesis whose chain
//                    shares < setMinSharedBlobs blobs with the chains already taken
//   6. one collision per crossing (pass 2, vertices from the winners only, weight = credit):
//                    supported: in the main vertex of its crossing with >= 1 other winner
//                    displaced: in a weaker vertex of its crossing -> moved to another hypothesis
//                               of the track that joins a supported main vertex, else no chain
//                    lone     : the only winner of its crossing -> moved as above if possible,
//                               else kept only if robust: credit >= setMinLoneCredit and
//                               >= setAmbiguityRatio x its other hypotheses (its other crossings
//                               are far away in z; hypotheses with the same Si hits are not
//                               counted with setLoneIgnoreSameSi, default on), else no chain
//   INTT time: with setRequireInttForCollision (default on) a collision exists only if a winner
//                    with INTT hits is in it (the main vertex of a crossing is the largest one
//                    with such a winner); MVTX-only chains carry no time (MVTX integration, and the
//                    TPC-guided search finds an MVTX-only chain near the TPC z of any crossing), so
//                    they only join such collisions: a lone MVTX-only winner, or a collision made
//                    of MVTX-only winners only, is moved (as displaced) or gets no chain; a lone
//                    MVTX-only winner that passes the lone robustness test keeps its chain without a
//                    collision (setKeepUntimedLone, default on): tier 1, status AmbiguousWithVertex,
//                    silicon_vertex_z NaN, not in the collision counts
//   re-selection: a track that would end without a chain may take another of its candidate chains
//                    (any chain of its crossing trajectories with chi2/ndf <= cut, not only the
//                    best) that joins a supported INTT-timed collision (setVertexChainReselect)
//   split collisions (done at the start of pass 2, before any track is dropped): collisions in
//                    crossings within setMergeMaxCrossings whose Si vertices agree
//                    (|dz_vtx| <= min(setMergeZWindow, setMergeZNSigma sigma_z sqrt(1/n1+1/n2)))
//                    are one collision if a common crossing fits the dz of all their tracks
//                    (chi2 increase < setMergeMaxDeltaChi2); its tracks move to that crossing
//                    (same chain, dz recomputed), so one collision is in exactly one crossing
//   monitoring: crossings with a second vertex of >= 2 INTT-timed winners are counted (End) and
//                    printed (Verbosity > 1): two collisions in one crossing, the second suppressed
//   tier 0: supported, best score >= setAmbiguityRatio x the best of another crossing (for a
//           merged collision: the same chain at another crossing does not count, and the
//           collision needs >= setMergeTier0MinTracks tracks)
//   tier 1: supported but ambiguous, moved, merged (small), lone kept, or no chain and the TPC z0 matches exactly
//           one collision (>= 2 tracks, credit sum >= setMinSupportTpcOnly) of its crossings
//   tier 2: the rest (skipped by Tpc_PolyClusterizer with its default max tier 1)
//   SI_HIT_SEED_EVENT: exactly the winning chain of every track is selected, with tpc_crossing /
//   tpc_dz_offset of the winning hypothesis (SiTrajectoryFitter, matcher, display).
//   Needs the measured Si - TPC offset in tpc_dz_offset (SiHitSeedReco stores it with every z
//   model) and, in streaming mode, SiHitSeedReco::setInttCrossingWindow so that the chain of a
//   hypothesis has INTT hits of its own crossing.

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class TpcCrossingDecisionContainer;

class TpcCrossingResolver : public SubsysReco
{
 public:
  explicit TpcCrossingResolver(const std::string& name = "TpcCrossingResolver");
  ~TpcCrossingResolver() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setTriggeredMode(bool v) { m_triggeredMode = v; }
  void setTriggeredCrossing(short v) { m_triggeredCrossing = v; }

  void setAssembledTrackNodeName(const std::string& n) { m_assembledNodeName = n; }
  void setTrajectoryNodeName(const std::string& n) { m_trajectoryNodeName = n; }
  void setSiSeedNodeName(const std::string& n) { m_siSeedNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }

  // 1. chains considered for a trajectory: agreement chi2/ndf cut
  void setMaxChainChi2(double v) { m_maxChainChi2 = v; }
  // 2. match: Si - TPC z resolution at the silicon [cm] and weight of the r*phi chi2/n
  void setDzSigma(double cm) { m_dzSigma = cm; }
  void setRphiChi2Weight(double v) { m_rphiWeight = v; }
  //    quality: power of n_layers / 7, weight of the mean hw delta-delta-phi pull^2,
  //    factor for chains without INTT hit (no time information)
  void setLayerPower(double v) { m_layerPower = v; }
  void setDdphiWeight(double v) { m_ddphiWeight = v; }
  void setNoInttFactor(double v) { m_noInttFactor = v; }
  // 4. vertices: a vote joins a cluster if within this distance of its weighted mean [cm]
  void setClusterZWindow(double cm) { m_clusterZWindow = cm; }
  // 5. single-track Si z resolution at the beam [cm] (measured ~0.22); the support of a track by
  //    the vertex of n_other other tracks uses sigma^2 = sigma_z^2 (1 + 1/n_other), the TPC-only
  //    match sigma_tpc^2 + sigma_z^2 / n, the split-collision z test sigma_z sqrt(1/n1 + 1/n2);
  //    setTpcZSigma: TPC-only z0 resolution [cm]
  void setSiZSigma(double cm) { m_siZSigma = cm; }
  void setTpcZSigma(double cm) { m_tpcZSigma = cm; }
  //    event factor W / (W + scale); floor added to it (a lone track keeps credit * floor)
  void setSupportScale(double v) { m_supportScale = v; }
  void setEventFloor(double v) { m_eventFloor = v; }
  //    ownership: chains sharing at least this many blobs conflict
  void setMinSharedBlobs(unsigned int n) { m_minSharedBlobs = n; }
  //    tier 0 if best >= ratio * best hypothesis in another crossing; also the lone-track test
  void setAmbiguityRatio(double v) { m_ambiguityRatio = v; }
  // 6. a lone track (only winner of its crossing) is kept only with at least this credit
  void setMinLoneCredit(double v) { m_minLoneCredit = v; }
  //    the ambiguity ratio of a lone track ignores its hypotheses with the same Si hits (they differ
  //    only in the crossing, chosen by dz; the lone track is tier 1 anyway)
  void setLoneIgnoreSameSi(bool v) { m_loneIgnoreSameSi = v; }
  //    TPC-only tracks: minimum collision support (credit sum of its tracks) to be matched
  void setMinSupportTpcOnly(double v) { m_minSupportTpcOnly = v; }
  //    false: no TPC-only match here, tracks without a Si chain get tier 2 at the crossing closest to z = 0
  //    (SiSeedingOrchestrator decides them later with TpcCrossingFinder against the final collisions)
  void setTpcOnlyMatch(bool v) { m_tpcOnlyMatch = v; }
  // collisions need a winner with INTT hits; MVTX-only chains only join them
  void setRequireInttForCollision(bool v) { m_requireInttForCollision = v; }
  //    a lone winner with an MVTX-only chain that passes the lone robustness test keeps its Si
  //    association without making a collision: tier 1, status AmbiguousWithVertex, crossing from
  //    its dz, silicon_vertex_z NaN (its own z = tpc_z0 + delta_z), not counted as a collision;
  //    credit threshold: setMinUntimedLoneCredit (negative = setMinLoneCredit)
  void setKeepUntimedLone(bool v) { m_keepUntimedLone = v; }
  //    a track that would lose its chain (displaced, lone not robust) takes another of its
  //    candidate chains (not only the best one per crossing trajectory) that joins a supported
  //    INTT-timed collision within the cluster window (tier 1, moved)
  void setVertexChainReselect(bool v) { m_vertexChainReselect = v; }
  void setMinUntimedLoneCredit(double v) { m_minUntimedLoneCredit = v; }
  // 7. split collisions: merge on/off, max crossing distance, max Si z distance [cm],
  //    max chi2 increase of the common crossing (sigma = setDzSigma)
  void setMergeSplitCollisions(bool v) { m_mergeSplitCollisions = v; }
  void setMergeMaxCrossings(int n) { m_mergeMaxCrossings = n; }
  void setMergeZWindow(double cm) { m_mergeZWindow = cm; }  // upper limit of the z test [cm]
  void setMergeZNSigma(double n) { m_mergeZNSigma = n; }
  //    a merged track is tier 0 (if otherwise clear) when its collision has at least this many tracks
  void setMergeTier0MinTracks(unsigned int n) { m_mergeTier0MinTracks = n; }
  void setMergeMaxDeltaChi2(double v) { m_mergeMaxDeltaChi2 = v; }
  // hypotheses with |TPC z0| beyond this are not used [cm]
  void setMaxVertexZ(double cm) { m_maxVertexZ = cm; }

 private:
  bool m_triggeredMode{true};
  short m_triggeredCrossing{0};

  std::string m_assembledNodeName{"TPC_ASSEMBLEDTRACKS"};
  std::string m_trajectoryNodeName{"TPC_CROSSINGTRAJECTORIES"};
  std::string m_siSeedNodeName{"SI_HIT_SEED_EVENT"};
  std::string m_outputNodeName{"TPC_CROSSING_DECISIONS"};

  double m_maxChainChi2{3.0};
  double m_dzSigma{0.6};
  double m_rphiWeight{1.0};
  double m_layerPower{1.0};
  double m_ddphiWeight{1.0};
  double m_noInttFactor{0.5};
  double m_clusterZWindow{0.5};
  double m_siZSigma{0.22};
  double m_tpcZSigma{0.5};
  double m_supportScale{1.0};
  double m_eventFloor{0.05};
  unsigned int m_minSharedBlobs{2};
  double m_ambiguityRatio{2.0};
  double m_minLoneCredit{0.3};
  bool m_loneIgnoreSameSi{true};
  double m_minSupportTpcOnly{1.0};
  bool m_tpcOnlyMatch{true};
  double m_maxVertexZ{25.0};
  bool m_requireInttForCollision{true};
  bool m_keepUntimedLone{true};
  bool m_vertexChainReselect{true};
  double m_minUntimedLoneCredit{-1.0};
  bool m_mergeSplitCollisions{true};
  int m_mergeMaxCrossings{2};
  double m_mergeZWindow{1.0};
  double m_mergeZNSigma{3.0};
  unsigned int m_mergeTier0MinTracks{3};
  double m_mergeMaxDeltaChi2{3.0};

  TpcCrossingDecisionContainer* m_decisions{nullptr};
  unsigned long m_nEvents{0};
  unsigned long m_nTracks{0};
  unsigned long m_nTier[4]{0, 0, 0, 0};  // tier 0, 1, 2, no crossing trajectory
  unsigned long m_nHypotheses{0};        // crossing trajectories with a chain
  unsigned long m_nChangedByEvent{0};    // winner is not the hypothesis with the highest credit
  unsigned long m_nOtherChain{0};        // winner chain is not the one SiHitSeedReco selected
  unsigned long m_nLostByOwnership{0};   // tracks with hypotheses left without a chain
  unsigned long m_nTpcOnly{0};           // tier 1 from the TPC-only match
  unsigned long m_nMoved{0};             // pass 2: moved to another hypothesis at a supported vertex
  unsigned long m_nLoneKept{0};          // pass 2: lone and robust
  unsigned long m_nLoneSameSi{0};        // of those: kept only because same-Si hypotheses are ignored
  unsigned long m_nDropped{0};           // pass 2: displaced or lone and not robust -> no chain
  unsigned long m_nReselected{0};        // pass 2: another candidate chain joining a timed collision
  unsigned long m_nUntimedKept{0};       // pass 2: lone MVTX-only, Si association kept without a collision
  unsigned long m_nUntimed{0};           // pass 2: MVTX-only winners without an INTT-timed collision -> no chain
  unsigned long m_nMergedCollisions{0};  // 7: collisions merged into a neighbouring crossing's one
  unsigned long m_nMergeMoved{0};        // 7: tracks moved to the common crossing
  unsigned long m_nSecondTimed{0};       // crossings with a second vertex of >= 2 INTT-timed winners
  std::map<unsigned int, unsigned long> m_multiplicity;     // Si-matched tracks per collision (crossing) -> count
  std::map<unsigned int, unsigned long> m_multiplicityAll;  // + TPC-only matched tracks
  std::vector<double> m_dzTier0;                           // dz of the tier-0 winners (End: width)
  std::map<std::string, unsigned long> m_noChainCategory;                     // why tracks have no Si chain
  std::map<std::pair<std::string, int>, unsigned long> m_noChainTier;         // (why, tier) -> count
};
#endif