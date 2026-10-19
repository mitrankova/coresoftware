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
//                    * exp(-(z - z_main)^2 / 2 sigma_z^2)), W(-h) = vote sum of the main vertex
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
//                               are far away in z), else no chain
//   7. split collisions: collisions at the same Si z in crossings within setMergeMaxCrossings
//                    are one collision if a common crossing fits the dz of all their tracks
//                    (chi2 increase < setMergeMaxDeltaChi2); its tracks move to that crossing
//                    (same chain, dz recomputed), so one collision is in exactly one crossing
//   tier 0: supported, best score >= setAmbiguityRatio x the best of another crossing
//   tier 1: supported but ambiguous, moved, lone kept, or no chain and the TPC z0 matches exactly
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
  // 5. collision z resolution of a track with a chain / TPC-only [cm]
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
  //    TPC-only tracks: minimum collision support (credit sum of its tracks) to be matched
  void setMinSupportTpcOnly(double v) { m_minSupportTpcOnly = v; }
  // 7. split collisions: merge on/off, max crossing distance, max Si z distance [cm],
  //    max chi2 increase of the common crossing (sigma = setDzSigma)
  void setMergeSplitCollisions(bool v) { m_mergeSplitCollisions = v; }
  void setMergeMaxCrossings(int n) { m_mergeMaxCrossings = n; }
  void setMergeZWindow(double cm) { m_mergeZWindow = cm; }
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
  double m_siZSigma{0.15};
  double m_tpcZSigma{0.5};
  double m_supportScale{1.0};
  double m_eventFloor{0.05};
  unsigned int m_minSharedBlobs{2};
  double m_ambiguityRatio{2.0};
  double m_minLoneCredit{0.3};
  double m_minSupportTpcOnly{1.0};
  double m_maxVertexZ{25.0};
  bool m_mergeSplitCollisions{true};
  int m_mergeMaxCrossings{2};
  double m_mergeZWindow{0.5};
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
  unsigned long m_nDropped{0};           // pass 2: displaced or lone and not robust -> no chain
  unsigned long m_nMergedCollisions{0};  // 7: collisions merged into a neighbouring crossing's one
  unsigned long m_nMergeMoved{0};        // 7: tracks moved to the common crossing
  std::map<unsigned int, unsigned long> m_multiplicity;     // Si-matched tracks per collision (crossing) -> count
  std::map<unsigned int, unsigned long> m_multiplicityAll;  // + TPC-only matched tracks
  std::vector<double> m_dzTier0;                           // dz of the tier-0 winners (End: width)
  std::map<std::string, unsigned long> m_noChainCategory;                     // why tracks have no Si chain
  std::map<std::pair<std::string, int>, unsigned long> m_noChainTier;         // (why, tier) -> count
};
#endif