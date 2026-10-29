// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SISEEDINGORCHESTRATOR_H
#define SITRACKRECO_SISEEDINGORCHESTRATOR_H

#include "SiHitSeedReco.h"
#include "SiTrajectoryFitter.h"
#include "TpcCrossingResolver.h"

#include <fun4all/SubsysReco.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;

// SiSeedingOrchestrator: ONE Fun4All module for the whole TPC-driven silicon seeding, iterative.
//
// It owns the three workers and runs them itself (they are NOT registered with Fun4All):
//   finder()   SiHitSeedReco        TPC-driven Si chaining (pass 1), rechain (pass 2), hit recovery (pass 3)
//   resolver() TpcCrossingResolver  one crossing + one chain per assembled track (pass 1)
//   fitter()   SiTrajectoryFitter   Si-only trajectories of the final selected chains (SI_TRAJECTORY)
// Configure them through these accessors with their usual setters; the shared frame setters below
// set the finder and the fitter at once.  Output nodes as before: SI_HIT_SEED_EVENT,
// TPC_CROSSING_DECISIONS, SI_TRAJECTORY, so everything downstream (Tpc_PolyClusterizer,
// SiTpcTrackSeedFiller, SiTpc_EventDisplay) is unchanged.
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
//   PASS 2  (setPass2) finder.rechain on the free blobs for the tracks without a chain (+ the weak
//           ones with setPass2RedoWeak), with the looser cuts of pass2Config() (e.g. 2 MVTX layers,
//           wider z windows for large |eta|).  Selection is VERTEX-CONSTRAINED: a candidate counts
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
//           used by any selected chain, windows / cuts of pass3Config().
//   FIT     fitter.process_event (setRunTrajectoryFit) on the final selection.
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
  void setNoInttFactor(double v) { m_noInttFactor = v; }

  // ---- pass 3
  void setPass3(bool v) { m_pass3 = v; }
  SiHitSeedReco::RecoveryConfig& pass3Config() { return m_pass3Cfg; }

 private:
  struct BootVertex
  {
    short crossing = 0;
    double z = 0;
    unsigned int n = 0;
  };

  std::unique_ptr<SiHitSeedReco> m_finder;
  std::unique_ptr<TpcCrossingResolver> m_resolver;
  std::unique_ptr<SiTrajectoryFitter> m_fitter;

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
  double m_noInttFactor = 0.5;

  bool m_pass3 = true;
  SiHitSeedReco::RecoveryConfig m_pass3Cfg;

  // run counters
  unsigned long m_nEvents = 0;
  unsigned long m_nGood = 0;
  unsigned long m_nWeak = 0;
  unsigned long m_nNone = 0;
  unsigned long m_nVertices = 0;
  unsigned long m_nVertexTracks = 0;
  unsigned long m_nPass2Targets = 0;
  unsigned long m_nPass2Chains = 0;
  unsigned long m_nPass2NoVertex = 0;     // candidate chains not compatible with a bootstrap vertex
  unsigned long m_nPass2Ambiguous = 0;    // tracks with compatible candidates in two crossings
  unsigned long m_nPass2Shared = 0;       // tracks whose candidates were all taken by better pass-2 tracks
  unsigned long m_nPass2Rescued = 0;      // tracks without a chain that got one
  unsigned long m_nPass2Improved = 0;     // weak tracks that took a better pass-2 chain
  unsigned long m_nPass2WeakKept = 0;     // weak tracks that kept their pass-1 chain
  unsigned long m_nPass3Chains = 0;
  std::array<unsigned long, 7> m_pass3Added{};
};

#endif
