#ifndef SIHITSEEDRECO_H
#define SIHITSEEDRECO_H

#include "SiAlignmentTransform.h"
#include "SiDetectorFrame.h"
#include "SiHitSeedData.h"
#include "SiTpcBeamAlignment.h"
#include "SiTpcHelixFit.h"
#include "SiTpcReference.h"

#include <fun4all/SubsysReco.h>
#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <map>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHG4CylinderGeomContainer;

// SiHitSeedReco: TPC-driven silicon chaining in HARDWARE coordinates.
//
// Input : TRKR_HITSET (MVTX + INTT hits), TPC_POLYTRACKS + TPC_POLYCLUSTERS
// Output: SI_HIT_SEED_EVENT (SiHitSeedEvent, transient)
//
// Per event
//   1. hits -> (Uphi, Vz) on the nominal layer cylinders (optional MVTX/INTT alignment),
//      (z, Uphi) binning, blobs (8-connected cells), cell index.               [unchanged]
//   2. per TPC track: TPC-only helix (SiTpcHelixFit on its poly clusters) -> silicon DETECTOR
//      frame (SiDetectorFrame + SiTpcBeamAlignment, same objects as SiTrajectoryFitter) ->
//      prediction (Uphi, Vz) on every silicon layer (SiTpcReference).
//   3. for the seed pairs (2 -> 1) and (1 -> 0):
//        ANCHOR window on the outer seed layer around the TPC prediction, large enough for the
//        Si-TPC misalignment and the TPC extrapolation (setAnchorWindow, r*phi and z in cm);
//        every blob with a hit inside is an anchor.
//        (the TPC pca z enters only here, as the reference position of the anchor window)
//        SEED window outer hit -> inner layer, relative to the outer hit, centre from the TPC:
//          dUphi = Uphi_TPC(inner) - Uphi_TPC(outer)                 (TPC curvature)
//          dz    = tanl_TPC * (s_TPC(inner) - s_TPC(outer))          (TPC slope only, no pca z)
//        half widths setSeedPhiHalfWindowBins [Uphi bins], setSeedZHalfWindowCm [cm].
//        Chain building exactly as before, in hardware coordinates: straight line through the
//        MVTX points (R_layer, Uphi, Vz), propagation to the missing MVTX layer and the INTT,
//        windows in Uphi / z bins (setPropagationPhiHalfWindowBins / setPropagationZHalfWindowBins).
//        INTT window centre (default): TPC prediction + straight-line fit of the MVTX offsets to
//        it (setInttCentre).  Candidate chains of one TPC track may share hits.
//   4. every chain -> PHYSICS space (TPC frame, as SiTrajectoryFitter: blob centroid per layer ->
//      detector -> global -> its clamshell half onto the TPC beam axis), residuals to the TPC
//      helix at the radius of each point, agreement chi2/ndf with per-layer sigmas, free common
//      z offset, z weight 0.25.
//      Ranking per TPC track (default AgreementFirst): chains with chi2/ndf <= 3 first, most
//      layers among them; layer-0/1 hardware delta-delta-phi pull^2 as tie-breaker.
//   5. exclusive arbitration across TPC tracks (default on): a chain sharing >= 2 blobs with a
//      chain already taken by a better-agreeing TPC track is not used; the track falls back to
//      its next chain or is Rejected.
//
// Defaults (run-25 Au+Au tuning, Oct 2026): anchor window 1 cm r*phi x 5 cm z (z not used for
// ordering); seed window +-2 Uphi bins, +-1 cm; propagation +-2 Uphi bins, +-5 z bins; chains
// >= 3 layers, >= 2 MVTX; sigmas per layer from the End() residual-width printout; agreement
// cut 3; INTT centre TpcOffsetSlope; ddphi check layers 0 and 1 at 2 sigma, score only.
// Frames, beam alignment, Bz and debug settings are NOT defaulted here: set them in the macro.
class SiHitSeedReco : public SubsysReco
{
 public:
  explicit SiHitSeedReco(const std::string& name = "SiHitSeedReco");
  ~SiHitSeedReco() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  // ---- nodes
  void setHitNodeName(const std::string& s) { m_hitNodeName = s; }
  void setTpcClusterNodeName(const std::string& s) { m_tpcClusterNodeName = s; }
  void setTpcTrackNodeName(const std::string& s) { m_tpcTrackNodeName = s; }
  // non-empty: TPC input = Tpc_CrossingTrajectoryContainer of this name (Tpc_PolyTripletReco,
  // one trajectory per assembled track and crossing) instead of the poly clusters + tracks
  void setTpcTrajectoryNodeName(const std::string& s) { m_tpcTrajectoryNodeName = s; }
  void setOutputNodeName(const std::string& s) { m_outputNodeName = s; }

  // ---- chain requirements
  void setMinLayers(int v) { m_minLayers = v; }
  void setMinMvtxLayers(int v) { m_minMvtxLayers = v; }
  void setMinInttLayers(int v) { m_minInttLayers = v; }
  // crossing trajectories (streaming): an INTT hit is used for a trajectory only if its time
  // bucket is within +- n of the trajectory crossing, so the chain of a crossing hypothesis gets
  // only INTT hits of that crossing.  < 0: off (default; keep off in triggered mode, where the
  // trajectory crossing is not an INTT time bucket)
  void setInttCrossingWindow(int n) { m_inttCrossingWindow = n; }
  // true (default): INTT blobs are built per time bucket (hits of different crossings are never
  // merged into one blob); false: previous behaviour, all buckets of a strip region together
  void setInttBlobsPerTimeBucket(bool v) { m_inttBlobsPerBucket = v; }
  // true (default): blobs = clusters of touching pixels (MVTX, same chip / strobe) and
  // neighbouring strips (INTT, same sensor / z column / time bucket), as the hardware
  // clusterizers; false: previous behaviour, connected cells of the (z, Uphi) search grid
  // (cells of ~0.35 cm x 0.035 rad: hits of different particles within that are merged)
  void setBlobsFromPixelAdjacency(bool v) { m_blobsFromPixels = v; }
  void setRefitWithIntt(bool v) { m_refitWithIntt = v; }

  // ---------------------------------------------------------------------------
  // Silicon alignment (no ACTS): same parameters and transform as
  // trackbase/AlignmentTransformation, applied to every MVTX/INTT hit before (Uphi, Vz) are used.
  // The aligned 3D point is projected back onto the nominal layer cylinder.
  // ---------------------------------------------------------------------------
  void setApplySiAlignment(bool v) { m_applyAlignment = v; }
  // Physics points only (chaining always uses the hardware Uphi, Vz):
  // true (default): silicon hits are placed on the real FLAT sensors (MVTX staves tilted by the
  // turbo angle, sensor offset 0.058 cm in the chip, INTT flat ladders), as in PHG4MvtxDetector /
  // CylinderGeom_Mvtx.  false: previous approximation, hits on the nominal cylinder R_layer at
  // Uphi = stave phi + lx cos(tilt) / R (MVTX: up to ~2 mm wrong in r, ~0.7 mm in r*phi).
  void setUseFlatSensorGeometry(bool v) { m_flatSensors = v; }
  // Physics points from the OFFICIAL tracking geometry (ActsGeometry node: MVTX displacement,
  // INTT survey geometry, TRACKINGALIGNMENT from the CDB), with the local pixel / strip
  // coordinates of CylinderGeom_Mvtx / CylinderGeomIntt: the same global positions as the
  // TRKR_CLUSTER clusters.  The frame centres are subtracted (detector frame = global - centre),
  // so SiDetectorFrame::toGlobal gives back the official global position.  setUseFlatSensorGeometry
  // and the alignment file then affect only the chaining.
  void setUseActsGeometry(bool v) { m_useActsGeometry = v; }
  // z added to the official positions (setUseActsGeometry): the official MVTX z contains the global
  // displacement D_z = -0.666 cm, which the own geometry (and the TPC z tuning so far) does not.
  // +0.666 keeps the Si physics z in the same z frame as before (MVTX = nominal z) while the INTT
  // keeps its correct z relative to the MVTX.  0: true global z.
  void setActsZShiftCm(double dz) { m_actsZShift = dz; }
  // per-layer table at End: own physics position (flat sensors / cylinders) - official position
  void setCompareWithActsGeometry(bool v) { m_compareActs = v; }
  void setSiAlignmentParamsFile(const std::string& f) { m_alignmentFile = f; }
  void setUseNewSiliconRotationOrder(bool v) { m_alignment.setUseNewSiliconRotationOrder(v); }
  const SiAlignmentTransform& alignment() const { return m_alignment; }

  // ---------------------------------------------------------------------------
  // TPC -> silicon frames: keep identical to SiTrajectoryFitter / SiTpcTrackMatcher / display
  // ---------------------------------------------------------------------------
  void setUphiRotation(double rad) { m_frame.setUphiRotation(rad); }
  void setDetectorCenterMm(double x, double y, double z = 0.0) { m_frame.setDetectorCenterMm(x, y, z); }
  // INTT centre (layers 3-6); not set: the INTT uses the detector (MVTX) centre.  With the
  // official TRACKINGALIGNMENT file: setDetectorCenterMm(D) + setInttCenterMm(0, 0, 0)
  void setInttCenterMm(double x, double y, double z = 0.0) { m_frame.setInttCenterMm(x, y, z); }
  void setBeamAlignment(const SiTpcBeamAlignment& a) { m_beamAlign = a; }
  void setApplyBeamAlignment(bool v) { m_beamAlign.setEnabled(v); }
  void setBz(double tesla) { m_bz = tesla; }
  void setMinTpcClusters(unsigned int n) { m_minTpcClusters = n; }
  // added to every TPC cluster z before the TPC helix fit [cm] (TPC t0 / drift offset)
  void setTpcZOffset(double cm) { m_tpcZOffset = cm; }
  // within this r*phi distance [cm] of a clamshell boundary the anchor window is opened with
  // the translations of both halves
  void setClamshellBoundaryMargin(double cm) { m_boundaryMargin = cm; }

  // ---------------------------------------------------------------------------
  // Anchor window on the outer seed layer, around the TPC prediction [cm].
  // Large: it only limits the search space (Si-TPC misalignment, TPC extrapolation, t0).
  // ---------------------------------------------------------------------------
  // zCm <= 0: no z cut at all (the TPC z is not used for the anchor; whole z range searched)
  void setAnchorWindow(double rphiCm, double zCm)
  {
    m_anchorRphi = rphiCm;
    m_anchorZ = zCm;
  }
  // weight of the z term in the anchor ordering metric (rdphi/rphiWin)^2 + w (dz/zWin)^2, which
  // decides which anchors are tried first and which survive setMaxAnchorsPerTrack; 0 = r*phi only
  void setAnchorZWeight(double w) { m_anchorZWeight = std::max(0.0, w); }
  void setMaxAnchorsPerTrack(unsigned int n) { m_maxAnchors = n; }
  // seed pairs: (2 -> 1) always; (1 -> 0) on/off (default on)
  void setUseInnerSeedPair(bool v) { m_useInnerPair = v; }

  // ---------------------------------------------------------------------------
  // Seed window outer hit -> inner layer: dz = tanl_TPC * ds_TPC (TPC slope only; independent
  // of the TPC pca z / t0), dUphi from the TPC curvature (see class comment)
  // ---------------------------------------------------------------------------
  void setSeedZHalfWindowCm(double cm) { m_seedZHalfCm = cm; }
  void setSeedZOffsetCm(double cm) { m_seedZOffsetCm = cm; }
  void setSeedPhiHalfWindowBins(double v) { m_seedPhiHalf = v; }
  // false: phi window centred on dUphi = 0 (straight track) instead of the TPC curvature
  void setSeedPhiFromTpcCurvature(bool v) { m_seedPhiFromTpc = v; }
  void setMaxSeedCandidatesPerAnchor(unsigned int n) { m_maxSeedCand = n; }

  // ---------------------------------------------------------------------------
  // Propagation windows (missing MVTX layer, INTT), as before: Uphi bins, z bins of the event
  // ---------------------------------------------------------------------------
  void setPhiHalfWindowBins(double v) { m_seedPhiHalf = v; m_propPhiHalf = v; }  // both at once
  void setPropagationPhiHalfWindowBins(double v) { m_propPhiHalf = v; }
  void setPropagationZHalfWindowBins(double v) { m_propZHalf = v; }
  // INTT propagation z half window [z bins]; < 0 (default): same as setPropagationZHalfWindowBins.
  // Should cover half an INTT strip length + the prediction error.
  void setInttPropagationZHalfWindowBins(double v) { m_inttPropZHalf = v; }
  // true: z centre of the MVTX propagation = z of the nearest chain hit + the TPC z difference
  // between the two layers (inclination from the TPC, offset from the chain) instead of the
  // straight line through the chain hits (r*phi centre unchanged). Default false.
  void setMvtxPropagationZFromTpc(bool v) { m_mvtxPropZFromTpc = v; }

  // ---------------------------------------------------------------------------
  // Step-to-step delta-delta-phi constraint in hardware coordinates (layers 0 and 1), from the
  // fits of h3_hw_dphi_vs_pred_l{L} (SiTpc_TrackResiduals + FitHwDphiCorrelation.C):
  //   dUphi(L)      = Uphi(L) - Uphi(L+1)
  //   dUphi_pred(L) = dUphi(L+1) * (r_L - r_{L+1}) / (r_{L+1} - r_PP)   PP = next outer chain layer
  //   expected dUphi(L) = offset_L + slope_L * dUphi_pred(L),  width sigma_L  [rad]
  // Used in two places (defaults: check on for layers 0 and 1 at 2 sigma, score only;
  // layer-0 propagation window off):
  //   setHwDdphiLayer0Propagation(true): when layer 0 is reached by MVTX propagation from layers
  //     2 and 1, its Uphi window is centred on the expected dUphi(0) with half width
  //     nSigma * sigma_0 (instead of the straight line and setPropagationPhiHalfWindowBins);
  //     z window unchanged.
  //   setHwDdphiCheck(l0, l1): after the chain is complete (INTT attached), the pull of layer 0
  //     (needs layers 1, 2) and/or layer 1 (needs layer 2 and an outer layer: INTT, or none) is
  //     computed and stored in SiHitChain::hw_ddphi_pull; with setHwDdphiReject(true) a chain
  //     with |pull| > nSigma is dropped, otherwise pull^2 is added to the chain score.
  // ---------------------------------------------------------------------------
  void setHwDdphiParams(int layer, double offset, double slope, double sigma)
  {
    if (layer == 0 || layer == 1)
    {
      m_hwDdphiOffset[layer] = offset;
      m_hwDdphiSlope[layer] = slope;
      m_hwDdphiSigma[layer] = sigma;
    }
  }
  void setHwDdphiNSigma(double v) { m_hwDdphiNSigma = v; }
  void setHwDdphiLayer0Propagation(bool v) { m_hwDdphiProp0 = v; }
  void setHwDdphiCheck(bool layer0, bool layer1)
  {
    m_hwDdphiCheck[0] = layer0;
    m_hwDdphiCheck[1] = layer1;
  }
  void setHwDdphiReject(bool v) { m_hwDdphiReject = v; }

  // ---------------------------------------------------------------------------
  // Selection (physics space): agreement with the TPC helix
  //   chi2/ndf = [sum (rdphi/sigRphi)^2 + sum ((dz - offset)/sigZ)^2] / ndf
  //   offset = weighted mean dz if setTpcAgreementFreeZOffset(true) (TPC z depends on t0)
  // ---------------------------------------------------------------------------
  // sigma per Si layer (0-2 MVTX, 3-6 INTT); the residual to the TPC helix is dominated by the
  // TPC extrapolation, so it depends on the layer radius - measure it with the End() printout
  // (see setMeasureResidualWidths) and set it per layer.
  void setTpcAgreementSigma(double rphiCm, double zMvtxCm, double zInttCm)
  {
    setTpcAgreementSigma(rphiCm, rphiCm, zMvtxCm, zInttCm);
  }
  void setTpcAgreementSigma(double rphiMvtxCm, double rphiInttCm, double zMvtxCm, double zInttCm)
  {
    for (int l = 0; l < 7; ++l)
    {
      m_sigRphiL[l] = l <= 2 ? rphiMvtxCm : rphiInttCm;
      m_sigZL[l] = l <= 2 ? zMvtxCm : zInttCm;
    }
  }
  void setTpcAgreementSigmaLayer(int layer, double rphiCm, double zCm)
  {
    if (layer >= 0 && layer < 7)
    {
      m_sigRphiL[layer] = rphiCm;
      m_sigZL[layer] = zCm;
    }
  }
  // End(): robust width (1.4826 * MAD) of the r*phi and z residuals per layer of the SELECTED
  // chains (all layers with a TPC residual, finite tpc_chi2), printed with the setter lines to
  // paste (default on; needs Verbosity >= 1).  The selection depends on the sigmas in use, so
  // iterate once or twice.
  void setMeasureResidualWidths(bool v) { m_measureWidths = v; }

  // ---------------------------------------------------------------------------
  // INTT search window centre (propagation to layers 3-6):
  //   MvtxLine: straight line through the MVTX hits (behaviour before Oct 2026): right
  //     position and direction, no field curvature.
  //   TpcOffset: TPC prediction on that INTT layer + the chain's MEAN MVTX offset to the TPC
  //     prediction (r*phi in cm, z in cm).  Wrong whenever the MVTX segment's direction differs
  //     from the extrapolated TPC direction.
  //   TpcOffsetSlope (default): TPC prediction + a straight-line fit of the MVTX offsets to
  //     the TPC prediction, r*phi offset vs r and z offset vs path length s, extrapolated to the
  //     INTT layer.  TPC = curvature, MVTX hits = position and direction.  Offset only if the
  //     MVTX lever arm is too short (one layer).
  //   All TPC modes fall back to the straight line without a TPC prediction on that layer.
  // half width in Uphi bins for the TPC-guided window (z half width: setPropagationZHalfWindowBins)
  // ---------------------------------------------------------------------------
  enum class InttCentre
  {
    MvtxLine,
    TpcOffset,
    TpcOffsetSlope
  };
  void setInttCentre(InttCentre c) { m_inttCentre = c; }
  // kept for macros of the previous version: true = TpcOffset, false = MvtxLine
  void setInttFromTpcPrediction(bool v) { m_inttCentre = v ? InttCentre::TpcOffset : InttCentre::MvtxLine; }
  void setInttTpcPhiHalfWindowBins(double v) { m_inttTpcPhiHalf = v; }
  // z model of the Si - TPC agreement (the TPC z can be shifted by the crossing / t0 and its
  // tanl can be off):
  //   Fixed:       dz used as is
  //   Offset:      a free z offset is fitted and subtracted (default; same as FreeZOffset(true))
  //   OffsetSlope: a free offset and a free slope in path length (dz = a + b*s) are fitted, so a
  //                TPC tanl error does not count either; needs >= 3 z points to constrain anything
  //   OffsetPrior: the TPC z error at the silicon is mostly a common offset (crossing, t0, 3-point
  //                extrapolation) while the TPC inclination (tanl) is precise.  A z offset is
  //                fitted with a Gaussian prior of width setTpcAgreementZOffsetSigma:
  //                chi2_z = sum ((dz_i - a) / sigma_z,i)^2 + (a / sigma_a)^2.
  //                sigma_z,i (setTpcAgreementSigmaLayer) are then the widths AROUND the offset
  //                (End printout with this model), so a chain with the wrong inclination fails
  //                even if its mean z agrees, and the offset itself is limited by sigma_a.
  enum class TpcZModel
  {
    Fixed,
    Offset,
    OffsetSlope,
    OffsetPrior
  };
  void setTpcAgreementZModel(TpcZModel m) { m_zModel = m; }
  // OffsetPrior: width of the prior on the Si - TPC z offset [cm] (TpcCrossingResolver End: robust
  // sigma of the tier-0 dz)
  void setTpcAgreementZOffsetSigma(double cm) { m_zOffsetSigma = std::max(1e-3, cm); }
  void setTpcAgreementFreeZOffset(bool v) { m_zModel = v ? TpcZModel::Offset : TpcZModel::Fixed; }
  // weight of the z part relative to r*phi:
  //   tpc_chi2 = (chi2_rphi + w * chi2_z) / (n_rphi + [w > 0] * (n_z - n_zpar))
  // w = 1: z as r*phi; 0 < w < 1: z counts less (default 0.25); w = 0: r*phi only (z only stored)
  void setTpcAgreementZWeight(double w) { m_zWeight = std::max(0.0, w); }
  // true (default): if the TPC helix does not reach the radius of some chain layer (e.g. a
  // large-dca / low-pt TPC track with hits on the inner MVTX layers), the chain cannot belong to
  // that TPC track: tpc_chi2 = NaN.  false: previous behaviour, those layers are skipped and the
  // chi2/ndf is computed from the rest (can be ndf 1 and a meaningless tiny chi2).
  void setRequireTpcResidualAllLayers(bool v) { m_requireAllResiduals = v; }
  void setMaxTpcChi2(double v) { m_maxTpcChi2 = v; }  // <= 0: no cut

  // ---------------------------------------------------------------------------
  // Ranking of the candidate chains of one TPC track (rank 0 = the one passed on):
  //   LayersFirst (behaviour before Oct 2026): most Si layers, then tpc_chi2, then score.
  //     A longer chain made of a neighbouring track's hits wins over the correct shorter one.
  //   AgreementFirst (default): chains with a finite tpc_chi2 <= setAgreementMaxChi2 (default 3) come
  //     first; among them most layers, then tpc_chi2, then score.  The rest follow, ordered by
  //     tpc_chi2.  If no chain agrees, rank 0 is the best-chi2 chain and it is NOT selected.
  // setMaxTpcChi2 still applies on top (rank 0 selected only if it also passes it).
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Exclusive chains across TPC tracks (default ON; off = behaviour before Oct 2026: every TPC track picks
  // its chain on its own, so the same Si hits can be the selected chain of several TPC tracks,
  // and a TPC track can take the chain of a nearby track instead of its own).
  // On: after all TPC tracks of the event are chained, the selections are arbitrated.  TPC
  // tracks are processed in order of the tpc_chi2 of their selected chain (best agreement
  // first); each takes the first chain in its own rank order that shares fewer than
  // setExclusiveMinSharedBlobs blobs with chains already taken by other TPC tracks and passes
  // the same cuts as the rank-0 selection (setMaxTpcChi2, agreement cut with AgreementFirst).
  // A TPC track left without such a chain gets status Rejected.  Chain ranks are unchanged;
  // only the "selected" flags move.  Works with setStoreSelectedOnly (pruned after arbitration).
  // ---------------------------------------------------------------------------
  void setExclusiveChains(bool v) { m_exclusiveChains = v; }
  void setExclusiveMinSharedBlobs(int n) { m_exclusiveMinShared = std::max(1, n); }

  enum class ChainRanking
  {
    LayersFirst,
    AgreementFirst
  };
  void setChainRanking(ChainRanking r) { m_chainRanking = r; }
  void setAgreementMaxChi2(double v) { m_agreementMaxChi2 = v; }
  void setMaxChainsPerTrack(unsigned int n) { m_maxChainsPerTrack = n; }
  // keep only the selected chain of every TPC track in SI_HIT_SEED_EVENT (default: keep all)
  void setStoreSelectedOnly(bool v) { m_storeSelectedOnly = v; }

  // ---------------------------------------------------------------------------
  // Debug dump of the chaining of selected TPC tracks (printed regardless of Verbosity):
  //   setDebugEvent(n): event n of THIS module (0 = first event SiHitSeedReco processes; events
  //     aborted earlier in the chain, e.g. by SkipFirstN, are not counted).  -1 = off.
  //   setDebugTpcTrack(i): only TPC_POLYTRACKS index i in that event; -1 = every TPC track.
  //   setDebugTpcAssembledId(id): only the track with this source assembled id (alternative).
  //   Every other TPC track of the debug event gets a one-line summary (index, pt, phi, z, status)
  //   so the track of interest can be identified first.
  // Per track: Tpc_PolyTrack parameters next to the SiTpcHelixFit refit used here, the predicted
  // (Uphi, z) on every Si layer, every blob on the outer seed layer within 3x the anchor window
  // (in / out, order), per anchor the seed candidates and the inner-layer blobs within 3x the
  // seed window, the result of every chain build, and per stored chain the residuals to the
  // TPC helix per layer, chi2/ndf, rank and selection.
  // ---------------------------------------------------------------------------
  void setDebugEvent(long n) { m_debugEvent = n; }
  void setDebugTpcTrack(int index) { m_debugTpcIndex = index; }
  void setDebugTpcAssembledId(long id) { m_debugAssembledId = id; }

  // ---- binning: number of bins.  The z range is recomputed every event from the hits.
  void setZBins(int n) { m_nz = n; }
  void setPhiBins(int n) { m_nphi = n; }

 private:
  struct Line3D
  {
    Eigen::Vector3d c;
    Eigen::Vector3d d;
  };
  struct Candidate
  {
    int blob = -1;
    int hit = -1;
    int layer = -1;
    int outer_hit = -1;
    double z = 0, phi = 0, outer_z = 0, outer_phi = 0;
    double res_z = 0, res_phi = 0, z_center = 0, z_half = 0, phi_center = 0, score = 0;
  };
  struct TpcTrackInput
  {
    int index = -1;
    unsigned int id = 0;
    unsigned int assembledId = 0;
    short crossing = 0;
    unsigned int nClusters = 0;
    SiTpcHelixFit::Result fit;
    // Tpc_PolyTrack's own parameters (for the debug dump)
    double polyX = 0, polyY = 0, polyZ = 0, polyPx = 0, polyPy = 0, polyPz = 0, polyCharge = 0;
    int polyFitStatus = 0;
  };

  int createNodes(PHCompositeNode*);
  int loadAlignment();
  // ideal detector-frame position of local point (lx, ly) of a sensor, and the sensor axes
  Eigen::Vector3d sensorPoint(int layer, int stave, int chip, int ladderphi, int ladderz, double lx, double ly) const;
  Eigen::Matrix3d sensorAxes(int layer, int stave, int ladderphi) const;
  // official global position of a pixel / strip [cm]; false if not available
  bool actsGlobal(TrkrDefs::hitsetkey hskey, int layer, unsigned int row, unsigned int col, int ladderz, Eigen::Vector3d& g) const;
  // ideal point -> aligned point (unchanged without alignment)
  Eigen::Vector3d alignPoint(TrkrDefs::hitsetkey hskey, const Eigen::Vector3d& p) const;
  // hardware coordinates (Uphi, Vz) with the alignment applied on the nominal cylinder (chaining)
  void alignHardware(TrkrDefs::hitsetkey hskey, int layer, double& z, double& phi) const;
  bool fillHits(PHCompositeNode*);
  void buildBlobs();
  void buildBlobsFromPixels();
  void buildHitGrid();
  std::vector<TpcTrackInput> collectTpcTracks(PHCompositeNode*) const;
  void processTpcTrack(const TpcTrackInput&);
  void resolveSharedChains();
  bool debugThisTrack(const TpcTrackInput&) const;
  void debugTrackHeader(const TpcTrackInput&, const SiTpcReference::TpcReference*) const;
  void debugAnchors(const SiTpcReference::TpcReference&, int outer, const std::vector<int>& anchors) const;
  void debugSeed(const SiBlob& seed, int inner, const SiTpcReference::TpcReference&, const std::vector<Candidate>& cands) const;
  void debugChain(const SiHitChain& ch, const char* tag) const;

  static double wrapDelta(double x, double period);
  static double circularMean(const std::vector<double>& values);
  std::pair<double, double> localCoordinates(int layer, unsigned int row, unsigned int col, int ladderz) const;
  double intrinsicPhi(int layer, int stave, int ladderphi, double lx) const;
  double longitudinalZ(int layer, int chip, int ladderz, double ly) const;
  Eigen::Vector3d xyz(int layer, double z, double phi) const;
  Line3D fitLine(const std::vector<Eigen::Vector3d>& points) const;
  bool propagate(const Line3D&, int layer, const Eigen::Vector3d& near, double& z, double& phi) const;

  template <class F>
  void forEachHitNear(int layer, double zbinCenter, double zHalf,
                      double phibinCenter, double phiHalf, F&& f) const;

  std::vector<int> anchorBlobs(const SiTpcReference::TpcReference&, int layer) const;
  void seedWindow(const SiTpcReference::TpcReference&, int outer, int inner, double zOuter,
                  double& zCenterBins, double& zHalfBins, double& phiCenterBins) const;
  std::vector<Candidate> seedCandidates(const SiBlob& seed, int innerLayer, const SiTpcReference::TpcReference&) const;
  bool buildChain(const SiBlob& seed, const Candidate& cand, SiHitChain& out,
                  const SiTpcReference::TpcReference* tpcRef = nullptr) const;
  void accumulateResidualWidths();
  void printResidualWidths() const;
  // dUphi(L) and dUphi_pred(L) for L = 0, 1 from the hardware hits per layer (layer -> Uphi);
  // false if L, L+1 or the next outer layer is missing
  bool hwDdphiTerms(int layer, const std::map<int, double>& uphiByLayer, double& dphi, double& dphiPred) const;
  void toPhysics(const SiTpcHelixFit::Result& tpcFit, SiHitChain& chain) const;

  // nodes
  std::string m_hitNodeName = "TRKR_HITSET";
  std::string m_tpcClusterNodeName = "TPC_POLYCLUSTERS";
  std::string m_tpcTrackNodeName = "TPC_POLYTRACKS";
  std::string m_tpcTrajectoryNodeName;
  std::string m_outputNodeName = "SI_HIT_SEED_EVENT";
  unsigned long long m_event = 0;
  SiHitSeedEvent* m_eventData = nullptr;

  int m_minLayers = 3;
  int m_minMvtxLayers = 2;
  int m_minInttLayers = 0;
  int m_inttCrossingWindow = -1;                // < 0: off
  bool m_inttBlobsPerBucket = true;
  bool m_blobsFromPixels = true;
  int m_curCrossing = 0;                        // crossing of the TPC input being chained
  mutable unsigned long m_inttTimeChecked = 0;  // INTT hits in a propagation window (time cut on)
  mutable unsigned long m_inttTimeRejected = 0; // ... of them outside the crossing window
  bool m_refitWithIntt = false;

  // silicon (hardware) alignment
  bool m_applyAlignment = false;
  bool m_flatSensors = true;
  bool m_useActsGeometry = false;
  bool m_compareActs = false;
  double m_actsZShift = 0.0;
  ActsGeometry* m_tGeometry = nullptr;
  PHG4CylinderGeomContainer* m_geomMvtx = nullptr;
  PHG4CylinderGeomContainer* m_geomIntt = nullptr;
  // own - official, per layer: n, sums of dx, dy, dz, dr, r*dphi and their squares
  std::array<std::array<double, 11>, 7> m_actsCmp{};
  std::string m_alignmentFile = "mvtx_align_params.txt";  // used only with setApplySiAlignment(true)
  SiAlignmentTransform m_alignment;

  // TPC -> silicon frames
  SiDetectorFrame m_frame;
  SiTpcBeamAlignment m_beamAlign;
  SiTpcReference::FrameMap m_map;
  double m_bz = 1.4;
  unsigned int m_minTpcClusters = 10;
  double m_tpcZOffset = 0.0;
  double m_boundaryMargin = 0.3;

  // anchor
  double m_anchorRphi = 1.0;  // [cm]
  double m_anchorZ = 5.0;     // [cm], <= 0: no z cut
  double m_anchorZWeight = 0.0;
  unsigned int m_maxAnchors = 16;
  bool m_useInnerPair = true;

  // seed window
  double m_seedZHalfCm = 1.0;
  double m_seedZOffsetCm = 0.0;
  double m_seedPhiHalf = 2.0;  // [Uphi bins]
  bool m_seedPhiFromTpc = true;
  unsigned int m_maxSeedCand = 4;

  // propagation windows
  double m_propPhiHalf = 2.0;  // [Uphi bins]
  double m_propZHalf = 5.0;    // [z bins]

  // step-to-step delta-delta-phi, layers 0 and 1 [rad]; defaults from FitHwDphiCorrelation.C on
  // run 76905 Au+Au without the TPC track-count cut (both charges)
  std::array<double, 2> m_hwDdphiOffset{{-0.00149, 0.00082}};
  std::array<double, 2> m_hwDdphiSlope{{1.4445, 1.7818}};
  std::array<double, 2> m_hwDdphiSigma{{0.0208, 0.0214}};
  double m_hwDdphiNSigma = 2.0;
  bool m_hwDdphiProp0 = false;
  std::array<bool, 2> m_hwDdphiCheck{{true, true}};
  bool m_hwDdphiReject = false;
  mutable std::array<unsigned long, 2> m_hwDdphiChecked{};
  mutable std::array<unsigned long, 2> m_hwDdphiRejected{};
  mutable unsigned long m_hwDdphiProp0Used = 0;

  // selection
  // [cm] per layer, robust residual widths of selected chains (End() printout, z model Offset)
  std::array<double, 7> m_sigRphiL{{0.310, 0.340, 0.378, 0.376, 0.389, 0.403, 0.439}};
  std::array<double, 7> m_sigZL{{0.277, 0.143, 0.287, 1.025, 0.990, 1.057, 0.874}};
  bool m_measureWidths = true;
  std::array<std::vector<float>, 7> m_widthRphi;
  std::array<std::vector<float>, 7> m_widthZ;

  // INTT window from the TPC prediction
  InttCentre m_inttCentre = InttCentre::TpcOffsetSlope;
  double m_inttTpcPhiHalf = 2.0;  // [Uphi bins]
  mutable unsigned long m_inttTpcUsed = 0;
  mutable unsigned long m_inttTpcFallback = 0;
  TpcZModel m_zModel = TpcZModel::Offset;
  double m_zOffsetSigma = 0.35;
  double m_inttPropZHalf = -1.0;
  bool m_mvtxPropZFromTpc = false;
  mutable unsigned long m_mvtxTpcZUsed = 0;
  double m_zWeight = 0.25;
  bool m_requireAllResiduals = true;
  double m_maxTpcChi2 = 0.0;
  ChainRanking m_chainRanking = ChainRanking::AgreementFirst;
  double m_agreementMaxChi2 = 3.0;
  // run counters for the ranking
  unsigned long m_rankChanged = 0;      // AgreementFirst picked a different rank-0 than LayersFirst would
  unsigned long m_noAgreeingChain = 0;  // AgreementFirst: TPC tracks with chains but none agreeing
  unsigned int m_maxChainsPerTrack = 32;
  bool m_storeSelectedOnly = false;

  // exclusive chains across TPC tracks
  bool m_exclusiveChains = true;
  int m_exclusiveMinShared = 2;
  unsigned long m_exclReassigned = 0;  // TPC tracks moved to another of their chains
  unsigned long m_exclLost = 0;        // TPC tracks left without a chain
  unsigned long m_exclKept = 0;        // TPC tracks that kept their rank-0 chain

  // debug dump
  long m_debugEvent = -1;
  int m_debugTpcIndex = -1;
  long m_debugAssembledId = -1;

  // binning
  int m_nz = 120;
  int m_nphi = 180;
  double m_zmin = -20.0;
  double m_zmax = 20.0;

  // per-event hit lookup: hits of each (layer, iz, iphi) cell, CSR layout
  std::vector<int> m_hitBlob;    // blob id of each hit, -1 if not in a blob
  std::vector<int> m_cellStart;  // size 7*nz*nphi + 1
  std::vector<int> m_cellHits;

  // run counters
  std::array<unsigned long, 7> m_statusCount{};
  unsigned long m_nChains = 0;

  const std::array<double, 7> m_radius{{2.523, 3.336, 4.148, 7.188 - 0.0036, 7.732 - 0.0036, 9.680 - 0.0036, 10.262 - 0.0036}};
};

#endif