// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SISEEDINGORCHESTRATOR_H
#define SITRACKRECO_SISEEDINGORCHESTRATOR_H

#include "SiHitSeedReco.h"
#include "SiTrajectoryFitter.h"
#include "TpcCrossingResolver.h"

#include <fun4all/SubsysReco.h>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;
class TpcCrossingFinder;

// SiSeedingOrchestrator: ONE Fun4All module for the TPC-driven silicon seeding AND every crossing / vertex
// decision of the event (nothing downstream changes a crossing).
//
// It owns the workers and runs them itself (they are NOT registered with Fun4All):
//   finder()        SiHitSeedReco        TPC-driven Si chaining (pass 1, beam-constrained fallback), rechain (pass 2),
//                                        hit recovery (pass 3)
//   resolver()      TpcCrossingResolver  one crossing + one chain per assembled track (pass 1)
//   tpcOnlyFinder() TpcCrossingFinder    crossing of the tracks WITHOUT a Si chain against the final collisions
//   fitter()        SiTrajectoryFitter   Si-only trajectories of the final selected chains (SI_TRAJECTORY, display only)
// Configure them through these accessors with their usual setters; the shared frame setters below
// set the finder and the fitter at once.  Output nodes: SI_HIT_SEED_EVENT, TPC_CROSSING_DECISIONS (final),
// the collision vertex map (setCollisionVertexMapName), SI_TRAJECTORY (with setRunTrajectoryFit).
//
// process_event:
//   PASS 1  finder.process_event + resolver.process_event (exactly the previous chain).
//           Each assembled track ends with a GOOD chain, a WEAK chain or no chain:
//             good = selected, >= setGoodMinLayers layers, >= setGoodMinMvtxLayers MVTX,
//                    INTT hit if setGoodRequireIntt, chi2/ndf <= setGoodMaxChi2,
//                    decision tier <= setGoodMaxTier
//             weak = selected, not good
//   VERTEX  bootstrap vertices from the GOOD chains only: per crossing, 1-D clusters of the Si z at
//           the beam (TPC z0 of the chain's trajectory + Si - TPC dz, as the resolver), window
//           setVertexClusterWindow, >= setVertexMinTracks tracks.
//   MASK    blobs of the good chains (and of the weak ones unless setPass2RedoWeak) are removed.
//   PASS 2  (setPass2) finder.rechainInputs on the free blobs for the tracks without a chain (+ the
//           weak ones with setPass2RedoWeak), with the looser cuts of pass2Config() (e.g. 2 MVTX layers,
//           wider z windows for large |eta|).  Only their TPC inputs (crossing trajectories) that can
//           change: not those without a usable TPC prediction in pass 1 (setPass2SkipUnpredictable) and,
//           in streaming mode, only at crossings with a bootstrap vertex (setPass2VertexCrossingsOnly).  Selection is VERTEX-CONSTRAINED: a candidate counts
//           only if its Si z at the beam is within setPass2VertexNSigma of a bootstrap vertex of ITS
//           crossing (sigma^2 = setTrackZSigma^2 (1 + 1/n_vertex)).  Per track the best candidate
//           (score = (layers/7) exp(-chi2/2) exp(-pull^2/2), x setNoInttFactor without INTT) wins if
//           it is >= setPass2AmbiguityRatio x the best of any other crossing; tracks in score order,
//           a candidate sharing >= setPass2MinSharedBlobs blobs with an earlier pass-2 winner is
//           skipped.  A weak track takes the pass-2 chain only if it has more layers (or as many
//           and a lower chi2/ndf).  The winner is selected in SI_HIT_SEED_EVENT and its decision
//           rewritten (crossing, tier setPass2Tier, Si vertex z, dz, SelectedByVertexLoose);
//           triggered mode: the decision is not touched (crossing fixed).
//   PASS 3  (setPass3) finder.recoverHits for every selected chain with < 7 layers, on the blobs not
//           used by any selected chain, windows / cuts of pass3Config().  Default: only the INTT pairs
//           (3, 4) / (5, 6) in which the chain has no hit, predicted with the vertex-assisted helix (TPC
//           curvature, centre from the chain + beam, z from the MVTX + the chain's bootstrap vertex).
//   COLLISIONS  (streaming) the final collisions: per crossing, 1-D clusters (setVertexClusterWindow, the same
//           clustering as the bootstrap vertices) of the Si z at the beam of every Si-matched track (selected chain,
//           decision tier <= setCollisionMaxTier, timed collision: finite silicon_vertex_z), >= setCollisionMinTracks.
//           Written as SvtxVertexMap (setCollisionVertexMapName): beam crossing, z = mean TPC z0 of its tracks (the
//           TPC z0 frame of the decisions and of TpcCrossingFinder; mean Si z - TPC z0 is printed at End), error
//           setTpcZSigma^2 / n, the assembled track ids as tracks.
//   TPC-ONLY (streaming, setTpcOnlyCrossing, default on) tracks without a selected Si chain: the resolver's own
//           3-point TPC-only match is switched off (TpcCrossingResolver::setTpcOnlyMatch(false)); tpcOnlyFinder()
//           (TpcCrossingFinder, unchanged: every INTT crossing, drift-time window, containment, TPC z0 from its
//           Garfield points, closest collision of the crossing) decides them against the collision map.  Its
//           decision replaces theirs in TPC_CROSSING_DECISIONS: finder tier 0 (|dz| <= its setMaxVertexDz) ->
//           setTpcOnlyVertexTier (default 1: fitted as TPC-only by SiTpcTrackSeedFiller), a close rival with
//           setTpcOnlyRejectAmbiguous -> 2, otherwise tier 2; no valid crossing: tier 2, crossing unchanged.
//           setTpcOnlyAddMissing: also tracks without a crossing trajectory (no resolver decision).
//   DIAGNOSTICS  Verbosity 2: INTT time buckets of every bootstrap vertex; Verbosity 1: vertices in
//           neighbouring crossings with the time buckets of the smaller one (End: summary).
//   FIT     fitter.process_event (setRunTrajectoryFit) on the final selection (only the event display uses it).
class SiSeedingOrchestrator : public SubsysReco
{
 public:
  explicit SiSeedingOrchestrator(const std::string& name = "SiSeedingOrchestrator");
  ~SiSeedingOrchestrator() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  // ---- workers (configure with their own setters)
  SiHitSeedReco& finder() { return *m_finder; }
  TpcCrossingResolver& resolver() { return *m_resolver; }
  SiTrajectoryFitter& fitter() { return *m_fitter; }
  TpcCrossingFinder& tpcOnlyFinder() { return *m_tpcOnlyFinder; }

  // ---- shared frame: finder + fitter
  void setUphiRotation(double rad)
  {
    m_finder->setUphiRotation(rad);
    m_fitter->setUphiRotation(rad);
  }
  void setDetectorCenterMm(double x, double y, double z = 0.0)
  {
    m_finder->setDetectorCenterMm(x, y, z);
    m_fitter->setDetectorCenterMm(x, y, z);
  }
  void setInttCenterMm(double x, double y, double z = 0.0)
  {
    m_finder->setInttCenterMm(x, y, z);
    m_fitter->setInttCenterMm(x, y, z);
  }
  void setBeamAlignment(const SiTpcBeamAlignment& a)
  {
    m_finder->setBeamAlignment(a);
    m_fitter->setBeamAlignment(a);
  }
  // triggered (crossing 0) or streaming: finder needs nothing, resolver and the pass-2 decisions do
  void setTriggeredMode(bool v)
  {
    m_triggeredMode = v;
    m_resolver->setTriggeredMode(v);
  }
  // decision container (resolver output, rewritten by pass 2)
  void setDecisionNodeName(const std::string& n)
  {
    m_decisionNodeName = n;
    m_resolver->setOutputNodeName(n);
  }
  // false: no resolver (e.g. TpcCrossingFinder made the decisions before; the finder's own
  // selection is pass 1).  Default true.
  void setRunResolver(bool v) { m_runResolver = v; }
  void setRunTrajectoryFit(bool v) { m_runFit = v; }

  // ---- pass-1 classification
  void setGoodMinLayers(int n) { m_goodMinLayers = n; }
  void setGoodMinMvtxLayers(int n) { m_goodMinMvtx = n; }
  void setGoodRequireIntt(bool v) { m_goodRequireIntt = v; }
  void setGoodMaxChi2(double v) { m_goodMaxChi2 = v; }
  void setGoodMaxTier(int t) { m_goodMaxTier = t; }

  // ---- bootstrap vertices (good chains)
  void setVertexClusterWindow(double cm) { m_vertexWindow = cm; }
  void setVertexMinTracks(unsigned int n) { m_vertexMinTracks = n; }
  // single-track Si z resolution at the beam [cm] (resolver setSiZSigma)
  void setTrackZSigma(double cm) { m_trackZSigma = cm; }

  // ---- pass 2
  void setPass2(bool v) { m_pass2 = v; }
  void setPass2RedoWeak(bool v) { m_pass2RedoWeak = v; }
  SiHitSeedReco::RechainConfig& pass2Config() { return m_pass2Cfg; }
  // shortcuts of pass2Config()
  void setPass2MinLayers(int n) { m_pass2Cfg.minLayers = n; }
  void setPass2MinMvtxLayers(int n) { m_pass2Cfg.minMvtxLayers = n; }
  // large |eta|: wider z windows (anchor z [cm], seed z [cm], MVTX propagation z [bins])
  void setPass2HighEtaWindows(double anchorZCm, double seedZHalfCm, double propZHalfBins)
  {
    m_pass2Cfg.anchorZCm = anchorZCm;
    m_pass2Cfg.seedZHalfCm = seedZHalfCm;
    m_pass2Cfg.propZHalfBins = propZHalfBins;
  }
  void setPass2MaxChi2(double v) { m_pass2MaxChi2 = v; }
  void setPass2VertexNSigma(double v) { m_pass2VertexNSigma = v; }
  void setPass2AmbiguityRatio(double v) { m_pass2AmbiguityRatio = v; }
  void setPass2MinSharedBlobs(unsigned int n) { m_pass2MinShared = n; }
  void setPass2Tier(unsigned int t) { m_pass2Tier = static_cast<unsigned char>(t); }
  // skip TPC inputs whose pass-1 status was too few TPC points / TPC fit failed / no prediction (default on)
  void setPass2SkipUnpredictable(bool v) { m_pass2SkipUnpredictable = v; }
  // streaming: chain again only TPC inputs at crossings with a bootstrap vertex (default on)
  void setPass2VertexCrossingsOnly(bool v) { m_pass2VertexCrossingsOnly = v; }
  void setNoInttFactor(double v) { m_noInttFactor = v; }

  // ---- pass 3
  void setPass3(bool v) { m_pass3 = v; }
  SiHitSeedReco::RecoveryConfig& pass3Config() { return m_pass3Cfg; }

  // ---- final collisions (streaming)
  void setCollisionVertexMapName(const std::string& n) { m_collisionMapName = n; }
  void setCollisionMaxTier(unsigned int t) { m_collisionMaxTier = t; }
  void setCollisionMinTracks(unsigned int n) { m_collisionMinTracks = n; }
  // TPC z0 resolution [cm], only for the error of the collision z
  void setTpcZSigma(double cm) { m_tpcZSigma = cm; }

  // ---- crossing of the tracks without Si chain (streaming)
  void setTpcOnlyCrossing(bool v) { m_tpcOnly = v; }
  void setTpcOnlyVertexTier(unsigned int t) { m_tpcOnlyVertexTier = t; }
  void setTpcOnlyRejectAmbiguous(bool v) { m_tpcOnlyRejectAmbiguous = v; }
  void setTpcOnlyAddMissing(bool v) { m_tpcOnlyAddMissing = v; }
  void setTpcOnlyDecisionNodeName(const std::string& n) { m_tpcOnlyDecisionName = n; }
  void setAssembledTrackNodeName(const std::string& n) { m_assembledNodeName = n; }

 private:
  struct BootVertex
  {
    short crossing = 0;
    double z = 0;
    unsigned int n = 0;
    std::vector<int> members;  // chain indices (good chains)
  };
  // one entry of the per-crossing z clustering (bootstrap vertices, final collisions)
  struct ZEntry
  {
    short crossing = 0;
    double z = 0;   // clustered z (Si z at the beam)
    int index = 0;  // caller's index (chain or track)
  };
  // per crossing: sorted 1-D clusters, a value joins while within setVertexClusterWindow of the running mean;
  // clusters with >= minTracks entries, as lists of entries
  std::vector<std::vector<ZEntry>> clusterZ(std::vector<ZEntry> entries, unsigned int minTracks) const;
  // COLLISIONS + TPC-ONLY stages
  void makeCollisions(PHCompositeNode* topNode, const std::map<unsigned int, int>& selectedOf);
  void decideTpcOnly(PHCompositeNode* topNode, const std::map<unsigned int, int>& selectedOf);

  std::unique_ptr<SiHitSeedReco> m_finder;
  std::unique_ptr<TpcCrossingResolver> m_resolver;
  std::unique_ptr<SiTrajectoryFitter> m_fitter;
  std::unique_ptr<TpcCrossingFinder> m_tpcOnlyFinder;

  bool m_triggeredMode = true;  // as TpcCrossingResolver
  bool m_runResolver = true;
  bool m_runFit = true;
  std::string m_decisionNodeName = "TPC_CROSSING_DECISIONS";

  int m_goodMinLayers = 3;
  int m_goodMinMvtx = 2;
  bool m_goodRequireIntt = false;
  double m_goodMaxChi2 = 3.0;
  int m_goodMaxTier = 1;

  double m_vertexWindow = 0.7;
  unsigned int m_vertexMinTracks = 2;
  double m_trackZSigma = 0.22;

  bool m_pass2 = true;
  bool m_pass2RedoWeak = true;
  SiHitSeedReco::RechainConfig m_pass2Cfg;
  double m_pass2MaxChi2 = 3.0;
  double m_pass2VertexNSigma = 3.0;
  double m_pass2AmbiguityRatio = 1.5;
  unsigned int m_pass2MinShared = 2;
  unsigned char m_pass2Tier = 1;
  bool m_pass2SkipUnpredictable = true;
  bool m_pass2VertexCrossingsOnly = true;
  double m_noInttFactor = 0.5;

  bool m_pass3 = true;
  SiHitSeedReco::RecoveryConfig m_pass3Cfg;

  std::string m_collisionMapName = "SiSeedCollisionVertexMap";
  unsigned int m_collisionMaxTier = 1;
  unsigned int m_collisionMinTracks = 2;
  double m_tpcZSigma = 0.5;
  bool m_tpcOnly = true;
  unsigned int m_tpcOnlyVertexTier = 1;
  bool m_tpcOnlyRejectAmbiguous = false;
  bool m_tpcOnlyAddMissing = false;
  std::string m_tpcOnlyDecisionName = "TPC_CROSSING_DECISIONS_TPCONLY";
  std::string m_assembledNodeName = "TPC_ASSEMBLEDTRACKS";
  bool m_tpcOnlyActive = false;  // InitRun: streaming, resolver and setTpcOnlyCrossing

  // run counters
  unsigned long m_nEvents = 0;
  unsigned long m_nGood = 0;
  unsigned long m_nWeak = 0;
  unsigned long m_nNone = 0;
  unsigned long m_nVertices = 0;
  unsigned long m_nVertexTracks = 0;
  unsigned long m_nPass2Targets = 0;
  unsigned long m_nPass2Chains = 0;
  unsigned long m_nPass2Inputs = 0;          // TPC inputs chained again
  unsigned long m_nPass2SkippedNoPred = 0;   // ... skipped: no usable TPC prediction in pass 1
  unsigned long m_nPass2SkippedNoVertex = 0; // ... skipped: crossing without a bootstrap vertex
  unsigned long m_nPass2NoVertex = 0;     // candidate chains not compatible with a bootstrap vertex
  unsigned long m_nPass2Ambiguous = 0;    // tracks with compatible candidates in two crossings
  unsigned long m_nPass2Shared = 0;       // tracks whose candidates were all taken by better pass-2 tracks
  unsigned long m_nPass2Rescued = 0;      // tracks without a chain that got one
  unsigned long m_nPass2Improved = 0;     // weak tracks that took a better pass-2 chain
  unsigned long m_nPass2WeakKept = 0;     // weak tracks that kept their pass-1 chain
  unsigned long m_nPass3Chains = 0;
  unsigned long m_nPass3WithVertex = 0;
  // bootstrap vertices in neighbouring crossings
  unsigned long m_nAdjacentPairs = 0;
  unsigned long m_adjSmallOwn = 0;        // INTT hits of the smaller vertex at its own crossing
  unsigned long m_adjSmallNeighbour = 0;  // ... at the other vertex's crossing
  unsigned long m_adjSmallOther = 0;      // ... elsewhere
  std::vector<double> m_adjDz;            // z(small) - z(big)
  std::array<unsigned long, 7> m_pass3Added{};
  // final collisions, TPC-only crossings
  unsigned long m_nCollisions = 0;
  unsigned long m_nCollisionTracks = 0;
  double m_sumSiMinusTpc = 0;
  unsigned long m_nTpcOnlyTracks = 0;     // tracks without a Si chain decided by tpcOnlyFinder()
  unsigned long m_nTpcOnlyNoDecision = 0; // ... without a finder decision
  unsigned long m_nTpcOnlyAdded = 0;      // ... without a resolver decision, added
  unsigned long m_nTpcOnlyMoved = 0;      // ... crossing changed
  std::map<std::pair<int, int>, unsigned long> m_tpcOnlyTiers;  // (old tier, new tier), -1 = no decision before
};

#endif
