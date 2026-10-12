// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCHITCHAINMATCHER_H
#define SITRACKRECO_SITPCHITCHAINMATCHER_H

#include "SiDetectorFrame.h"
#include "SiHitGeometryTable.h"
#include "SiHitIndex.h"
#include "SiTpcBeamAlignment.h"
#include "SiTpcChainMath.h"
#include "SiTpcHelixFit.h"
#include "SiTpcHitMatch.h"

#include <fun4all/SubsysReco.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;
class SiHitSeedEvent;
class SiTpcHitMatchContainer;

// SiTpcHitChainMatcher
//
// TPC-driven silicon chaining on HITS (not clusters), for AuAu where the TPC trajectory is known.
// Runs next to the independent silicon chaining (SiHitSeedReco -> SiTrajectoryFitter ->
// SiTpcTrackMatcher), which is left untouched and serves as the cross-check.
//
// Input : SI_HIT_SEED_EVENT   decoded, aligned silicon hits (SiHitSeedReco; only the hit list is used)
//         TPC_POLYTRACKS + TPC_POLYCLUSTERS   TPC tracks; a TPC-only SiTpcHelixFit gives the reference
// Output: SI_TPC_HITMATCHES   SiTpcHitMatchContainerv1 under DST/SVTX (best chain per TPC track,
//                             or all chains with setStoreAllChains)
//
// Per event
//   1. hit index: every hit -> (u, z, r) from SiHitGeometryTable (table lookups, no trig),
//      per-layer azimuth buckets with hits sorted in z (SiHitIndex).
//   2. per TPC track: TPC-only helix, translated into the silicon detector frame and intersected
//      with the 7 layer cylinders once (SiTpcChainMath::buildTpcReference).  Within
//      setClamshellBoundaryMargin of a clamshell boundary the other half's translation is
//      searched too (the halves are displaced by ~1 mm).  The silicon chain itself lives in the
//      TPC frame, so it stays a circle when it crosses the boundary.
//   3. layer walk (default order 2,1,0,3,4,5,6 as TpcSiliconCrossingMatcher):
//      - entry layer (chain still empty): window around the TPC prediction
//            |dphi_tpc| < entryPhi, |dtheta_tpc| < entryTheta
//      - later layers: the silicon chain is rebuilt from the attached groups with the TPC radius,
//        helicity and tan(lambda) (SiTpcChainMath::fitSiChain); its prediction gives
//            dphi, dtheta  and  ddphi = dphi - (a_L + b_L dphi_prev),  ddtheta likewise
//        window |ddphi| < ddPhi, |ddtheta| < ddTheta
//      - polar angles are measured from the TPC vertex: theta = atan2(s, z - z0_TPC); the theta
//        window is turned into a z range once per branch and layer, the hit test itself is
//        two subtractions.
//      - every hit in the window seeds a GROUP: all geometric neighbours are attached
//            MVTX: same chip and strobe, |drow| <= 1, |dcol| <= 1 (8-connectivity), plus the
//                  abutting chip of the same stave across the chip gap (boundary columns)
//            INTT: same sensor and time bucket, |drow| <= 1, same column (optional |dcol| <= 1)
//        growth is bounded by the footprint expected from the TPC incidence angles plus a slack:
//            n_rphi <= ceil(1 + t tan(alpha_inc)/pitch_rphi) + slack
//            n_z    <= ceil(1 + t |tanl|/cos(alpha_inc)/pitch_z) + slack   (MVTX)
//        the group centroid is the measurement.
//      - ONE index query per track and layer (union of all branch windows); the groups found are
//        then tested per branch with centroid arithmetic only.
//      - windows ACCEPT; candidates are RANKED by physical resolutions (setScoreSigma,
//        setEntryScoreSigma; INTT z = strip / sqrt(12)), so a loose window does not decide which
//        group wins.
//      - every accepted group opens a branch (setMaxBranchesPerLayer: MVTX 8, INTT 1 = closest
//        group per INTT layer as in TpcSiliconCrossingMatcher; at most setMaxChains branches);
//        a layer without any accepted group is a missing layer.  After every layer a branch is
//        dropped if one with at least as many layers has a chi2 lower by setBranchScoreMargin
//        (combinatorics control, not a residual cut).
//      - optional (default off): a second group on the adjacent (overlapping) stave / ladder,
//        compatible with the main group, recorded as an Overlap step (not used in the chain).
//   4. the best branch (n layers, then n MVTX layers, then chi2) is stored with all per-step
//      residuals, so the delta-delta cuts can be tuned from the output.
//
// All cuts are LOOSE by default (see the defaults below): tune ddphi / ddtheta from SiTpcHitChainQA.
class SiTpcHitChainMatcher : public SubsysReco
{
 public:
  explicit SiTpcHitChainMatcher(const std::string& name = "SiTpcHitChainMatcher");
  ~SiTpcHitChainMatcher() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  // ---- nodes
  void setHitEventNodeName(const std::string& s) { m_hitEventNodeName = s; }
  void setTpcClusterNodeName(const std::string& s) { m_tpcClusterNodeName = s; }
  void setTpcTrackNodeName(const std::string& s) { m_tpcTrackNodeName = s; }
  void setOutputNodeName(const std::string& s) { m_outputNodeName = s; }

  // ---- frames (use the same objects as SiTrajectoryFitter / SiTpcTrackMatcher)
  void setUphiRotation(double rad) { m_frame.setUphiRotation(rad); }
  void setDetectorCenterMm(double x, double y, double z = 0.0) { m_frame.setDetectorCenterMm(x, y, z); }
  void setBeamAlignment(const SiTpcBeamAlignment& a) { m_alignment = a; }
  void setApplyBeamAlignment(bool v) { m_alignment.setEnabled(v); }
  void setBz(double tesla) { m_bz = tesla; }
  // within this r*phi distance [cm] of a clamshell boundary the window is searched with the
  // translations of both halves (the halves are displaced by ~1 mm)
  void setClamshellBoundaryMargin(double cm) { m_boundaryMargin = cm; }
  void setMinTpcClusters(unsigned int n) { m_minTpcClusters = n; }

  // ---- hit geometry (SiHitGeometryTable)
  void setHitModelPlanar(bool v) { m_geometry.setModel(v ? SiHitGeometryTable::Planar : SiHitGeometryTable::Linear); }
  void setMvtxTiltSign(int s) { m_geometry.setMvtxTiltSign(s); }
  void setMvtxThicknessCm(double t) { m_geometry.setMvtxThicknessCm(t); }
  void setInttThicknessCm(double t) { m_geometry.setInttThicknessCm(t); }
  void setBucketRphiCm(double w) { m_index.setBucketRphiCm(w); }

  // ---- layer order of the walk (default 2,1,0,3,4,5,6)
  void setLayerOrder(const std::vector<int>& order) { m_layerOrder = order; }

  // ---- windows [rad].  det = 0: MVTX, 1: INTT
  void setEntryWindow(int det, double dphi, double dtheta)
  {
    if (det == 0 || det == 1)
    {
      m_entryPhi[det] = dphi;
      m_entryTheta[det] = dtheta;
    }
  }
  void setDeltaDeltaWindow(int det, double ddphi, double ddtheta)
  {
    if (det == 0 || det == 1)
    {
      m_ddPhi[det] = ddphi;
      m_ddTheta[det] = ddtheta;
    }
  }
  // expected dphi at `layer` given the previous step's dphi: mean = offset + slope * dphi_prev
  // (default 0, 0: ddphi = dphi; the old TpcSiliconCrossingMatcher default slope was 1)
  void setDeltaDeltaPhiMean(int layer, double offset, double slope)
  {
    if (layer >= 0 && layer < 7)
    {
      m_ddPhiOffset[layer] = offset;
      m_ddPhiSlope[layer] = slope;
    }
  }
  void setDeltaDeltaThetaMean(int layer, double offset, double slope)
  {
    if (layer >= 0 && layer < 7)
    {
      m_ddThetaOffset[layer] = offset;
      m_ddThetaSlope[layer] = slope;
    }
  }

  // ---- ranking (NOT acceptance): resolutions used for the candidate score [cm].
  // score = (ddphi R / sigRphi)^2 + (ddtheta / sigTheta)^2, sigTheta = sigZ sin^2(theta) / s.
  // INTT: sigZ is at least strip length / sqrt(12).  Entry layer: TPC extrapolation resolutions.
  void setScoreSigma(int det, double sigRphi, double sigZ)
  {
    if (det == 0 || det == 1)
    {
      m_sigRphi[det] = sigRphi;
      m_sigZ[det] = sigZ;
    }
  }
  void setEntryScoreSigma(int det, double sigRphi, double sigZ)
  {
    if (det == 0 || det == 1)
    {
      m_entrySigRphi[det] = sigRphi;
      m_entrySigZ[det] = sigZ;
    }
  }

  // ---- groups
  void setGroupExtentSlack(double pitches) { m_groupSlack = pitches; }
  void setMaxGroupHits(unsigned int n) { m_maxGroupHits = n; }
  void setInttColumnNeighbours(bool v) { m_inttColumnNeighbours = v; }
  void setAttachOverlapGroups(bool v) { m_attachOverlap = v; }
  // overlap group = adjacent stave / ladder, within these distances of the main group [cm]
  // (z: MVTX, INTT; INTT z is a strip).  Recorded as an Overlap step, not used in the chain.
  void setOverlapTolerance(double rphi, double zMvtx, double zIntt)
  {
    m_overlapTolRphi = rphi;
    m_overlapTolZ[0] = zMvtx;
    m_overlapTolZ[1] = zIntt;
  }

  // ---- combinatorics and selection
  // branches opened per branch and layer; det = 0 MVTX (default 8), 1 INTT (default 1: closest
  // compatible group per INTT layer, as attachClosestInttClusters in TpcSiliconCrossingMatcher)
  void setMaxBranchesPerLayer(unsigned int n) { m_maxBranchesPerLayer[0] = n; }
  void setMaxBranchesPerLayer(int det, unsigned int n)
  {
    if (det == 0 || det == 1)
    {
      m_maxBranchesPerLayer[det] = n;
    }
  }
  void setMaxChains(unsigned int n) { m_maxChains = n; }
  void setMissingLayerPenalty(double v) { m_missingPenalty = v; }
  // after every layer, branches with chi2 > best chi2 of this track + margin are dropped
  // (combinatorics control in score units, not a residual cut; large value = off)
  void setBranchScoreMargin(double v) { m_branchScoreMargin = v; }
  void setMinLayers(unsigned int n) { m_minLayers = n; }
  void setMinMvtxLayers(unsigned int n) { m_minMvtxLayers = n; }
  void setStoreAllChains(bool v) { m_storeAll = v; }
  // weights of the group z in the chain z0 (INTT z = strip centre)
  void setZWeights(double mvtx, double intt)
  {
    m_zWeightMvtx = mvtx;
    m_zWeightIntt = intt;
  }

  const SiHitGeometryTable& geometry() const { return m_geometry; }

 private:
  struct TrackInput
  {
    unsigned int index = 0;
    unsigned int id = 0;
    unsigned int sourceId = 0;
    SiTpcHelixFit::Result fit;
  };

  struct Group
  {
    int layer = -1;
    std::vector<int> members;  // indices into m_index.hits(layer)
    double u = 0, z = 0, r = 0, x = 0, y = 0;
    int flags = 0;
    float extentRphi = 0, extentZ = 0;
    TrkrDefs::hitsetkey hitsetkey = TrkrDefs::HITSETKEYMAX;
    int phiId = -1;
    int zId = -1;
  };

  struct Limits
  {
    double expRphi = 1, expZ = 1;
    double maxRphi = 1, maxZ = 1;
  };

  struct Step
  {
    int group = -1;
    SiTpcHitMatchStep rec;
  };

  struct Branch
  {
    std::vector<Step> steps;
    std::vector<SiTpcChainMath::SiChainPoint> points;
    SiTpcChainMath::SiChainFit si;
    double chi2 = 0;
    unsigned int nMissing = 0;
    unsigned int nMvtx = 0;
    unsigned int nIntt = 0;
    bool hasPrev = false;
    double prevDphi = 0;
    double prevDtheta = 0;
    unsigned int nLayers() const { return nMvtx + nIntt; }
  };

  struct Window
  {
    bool valid = false;
    bool entry = true;
    const SiTpcChainMath::LayerPrediction* tp = nullptr;
    double uRef = 0, sRef = 0, zRef = 0, thetaRef = 0;  // reference curve on the layer
    double meanPhi = 0, meanTheta = 0;                 // delta-delta model
    double wPhi = 0, wTheta = 0;                       // half widths [rad]
    double uc = 0, zlo = 0, zhi = 0;                   // search box
  };

  struct Candidate
  {
    int group = -1;
    SiTpcHitMatchStep rec;
    double score = 0;
  };

  struct Counters
  {
    std::array<unsigned long long, 7> visited{};
    std::array<unsigned long long, 7> inWindow{};
    std::array<unsigned long long, 7> groupsBuilt{};
    std::array<unsigned long long, 7> groupsAccepted{};
    std::array<unsigned long long, 7> attached{};
    unsigned long long branches = 0;
    unsigned int tpcTracks = 0;
    unsigned int references = 0;
    unsigned int matched = 0;
  };

  int createNodes(PHCompositeNode*);
  std::vector<TrackInput> collectTpcTracks(PHCompositeNode*) const;
  void processTrack(const TrackInput&, Counters&);
  Window makeWindow(const SiTpcChainMath::TpcReference&, const SiTpcChainMath::LayerPrediction&, const Branch&) const;
  void collectGroups(int layer, const std::vector<Window>&, const Limits&, std::vector<int>& groups, Counters&);
  bool evaluateGroup(const SiTpcChainMath::TpcReference&, const Window&, int group, const Limits&, Candidate&) const;
  Limits groupLimits(const SiTpcChainMath::TpcReference&, const SiTpcChainMath::LayerPrediction&) const;
  int groupOf(int layer, int k, const Limits&, Counters&);
  bool neighbours(int layer, const SiIndexedHit& a, const SiIndexedHit& b, bool& crossChip) const;
  bool extendBranch(const SiTpcChainMath::TpcReference&, const SiTpcChainMath::LayerPrediction&, const Branch&,
                    const Candidate& main, const std::vector<Candidate>& overlaps, Branch& out) const;
  void storeBranch(const TrackInput&, const SiTpcChainMath::TpcReference&, const Branch&, unsigned int rank);
  int detIndex(int layer) const { return SiHitGeometryTable::isMvtx(layer) ? 0 : 1; }

  // nodes
  std::string m_hitEventNodeName = "SI_HIT_SEED_EVENT";
  std::string m_tpcClusterNodeName = "TPC_POLYCLUSTERS";
  std::string m_tpcTrackNodeName = "TPC_POLYTRACKS";
  std::string m_outputNodeName = "SI_TPC_HITMATCHES";
  SiTpcHitMatchContainer* m_container = nullptr;

  // frames, geometry
  SiDetectorFrame m_frame;
  SiTpcBeamAlignment m_alignment;
  SiTpcChainMath::FrameMap m_map;
  double m_bz = 1.4;
  double m_boundaryMargin = 0.3;
  unsigned int m_minTpcClusters = 10;
  SiHitGeometryTable m_geometry;
  SiHitIndex m_index;

  // walk
  std::vector<int> m_layerOrder{2, 1, 0, 3, 4, 5, 6};
  // LOOSE defaults [rad]; index 0 = MVTX, 1 = INTT (INTT theta: 1.6 / 2.0 cm strips)
  std::array<double, 2> m_entryPhi{{0.10, 0.10}};
  std::array<double, 2> m_entryTheta{{0.10, 0.30}};
  std::array<double, 2> m_ddPhi{{0.05, 0.05}};
  std::array<double, 2> m_ddTheta{{0.05, 0.30}};
  std::array<double, 7> m_ddPhiOffset{};
  std::array<double, 7> m_ddPhiSlope{};
  std::array<double, 7> m_ddThetaOffset{};
  std::array<double, 7> m_ddThetaSlope{};
  double m_zAcceptanceMargin = 0.5;  // [cm] beyond the active half length

  // groups
  double m_groupSlack = 3.0;  // [pitches]
  unsigned int m_maxGroupHits = 64;
  bool m_inttColumnNeighbours = false;
  bool m_attachOverlap = false;  // off until validated on data with real stave overlaps
  double m_overlapTolRphi = 0.03;
  std::array<double, 2> m_overlapTolZ{{0.3, 2.5}};
  // ranking resolutions [cm]; INTT z: strip length / sqrt(12) (sigZ <= 0 -> strip only)
  std::array<double, 2> m_sigRphi{{0.005, 0.005}};
  std::array<double, 2> m_sigZ{{0.005, 0.0}};
  std::array<double, 2> m_entrySigRphi{{0.05, 0.05}};
  std::array<double, 2> m_entrySigZ{{0.10, 0.0}};

  // combinatorics / selection
  std::array<unsigned int, 2> m_maxBranchesPerLayer{{8, 1}};
  unsigned int m_maxChains = 64;
  double m_missingPenalty = 1.0;
  double m_branchScoreMargin = 50.0;
  unsigned int m_minLayers = 1;
  unsigned int m_minMvtxLayers = 1;
  bool m_storeAll = false;
  double m_zWeightMvtx = 1.0;
  double m_zWeightIntt = 0.01;

  // per event / per track work areas
  std::array<std::vector<uint32_t>, 7> m_stamp;
  std::array<std::vector<int>, 7> m_groupOfHit;
  uint32_t m_curStamp = 0;
  std::vector<Group> m_groups;
  unsigned int m_nextMatchId = 0;

  // run summary
  unsigned long m_nEvents = 0;
  double m_timeIndex = 0, m_timeTpc = 0, m_timeMatch = 0;
};

#endif
