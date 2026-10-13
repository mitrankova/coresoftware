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

#include <array>
#include <string>
#include <vector>

class PHCompositeNode;

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
//        Hits are NOT masked: chains of one TPC track may share hits, and so may chains of
//        different TPC tracks.  Identical chains (same blobs) are stored once per TPC track.
//   4. every chain -> PHYSICS space (TPC frame, as SiTrajectoryFitter: blob centroid per layer ->
//      detector -> global -> its clamshell half onto the TPC beam axis), residuals to the TPC
//      helix at the radius of each point, agreement chi2/ndf (optional free common z offset).
//      Ranking per TPC track: most silicon layers first, then the smallest chi2/ndf.
//      Rank 0 is SELECTED if its chi2/ndf <= setMaxTpcChi2 (<= 0: no cut).
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
  void setOutputNodeName(const std::string& s) { m_outputNodeName = s; }

  // ---- chain requirements
  void setMinLayers(int v) { m_minLayers = v; }
  void setMinMvtxLayers(int v) { m_minMvtxLayers = v; }
  void setMinInttLayers(int v) { m_minInttLayers = v; }
  void setRefitWithIntt(bool v) { m_refitWithIntt = v; }

  // ---------------------------------------------------------------------------
  // Silicon alignment (no ACTS): same parameters and transform as
  // trackbase/AlignmentTransformation, applied to every MVTX/INTT hit before (Uphi, Vz) are used.
  // The aligned 3D point is projected back onto the nominal layer cylinder.
  // ---------------------------------------------------------------------------
  void setApplySiAlignment(bool v) { m_applyAlignment = v; }
  void setSiAlignmentParamsFile(const std::string& f) { m_alignmentFile = f; }
  void setUseNewSiliconRotationOrder(bool v) { m_alignment.setUseNewSiliconRotationOrder(v); }
  const SiAlignmentTransform& alignment() const { return m_alignment; }

  // ---------------------------------------------------------------------------
  // TPC -> silicon frames: keep identical to SiTrajectoryFitter / SiTpcTrackMatcher / display
  // ---------------------------------------------------------------------------
  void setUphiRotation(double rad) { m_frame.setUphiRotation(rad); }
  void setDetectorCenterMm(double x, double y, double z = 0.0) { m_frame.setDetectorCenterMm(x, y, z); }
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
  void setAnchorWindow(double rphiCm, double zCm)
  {
    m_anchorRphi = rphiCm;
    m_anchorZ = zCm;
  }
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

  // ---------------------------------------------------------------------------
  // Selection (physics space): agreement with the TPC helix
  //   chi2/ndf = [sum (rdphi/sigRphi)^2 + sum ((dz - offset)/sigZ)^2] / ndf
  //   offset = weighted mean dz if setTpcAgreementFreeZOffset(true) (TPC z depends on t0)
  // ---------------------------------------------------------------------------
  void setTpcAgreementSigma(double rphiCm, double zMvtxCm, double zInttCm)
  {
    m_sigRphi = rphiCm;
    m_sigZ[0] = zMvtxCm;
    m_sigZ[1] = zInttCm;
  }
  void setTpcAgreementFreeZOffset(bool v) { m_freeZOffset = v; }
  void setMaxTpcChi2(double v) { m_maxTpcChi2 = v; }  // <= 0: no cut
  void setMaxChainsPerTrack(unsigned int n) { m_maxChainsPerTrack = n; }
  // keep only the selected chain of every TPC track in SI_HIT_SEED_EVENT (default: keep all)
  void setStoreSelectedOnly(bool v) { m_storeSelectedOnly = v; }

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
    unsigned int nClusters = 0;
    SiTpcHelixFit::Result fit;
  };

  int createNodes(PHCompositeNode*);
  int loadAlignment();
  void alignPoint(TrkrDefs::hitsetkey hskey, int layer, double& z, double& phi) const;
  bool fillHits(PHCompositeNode*);
  void buildBlobs();
  void buildHitGrid();
  std::vector<TpcTrackInput> collectTpcTracks(PHCompositeNode*) const;
  void processTpcTrack(const TpcTrackInput&);

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
  bool buildChain(const SiBlob& seed, const Candidate& cand, SiHitChain& out) const;
  void toPhysics(const SiTpcHelixFit::Result& tpcFit, SiHitChain& chain) const;

  // nodes
  std::string m_hitNodeName = "TRKR_HITSET";
  std::string m_tpcClusterNodeName = "TPC_POLYCLUSTERS";
  std::string m_tpcTrackNodeName = "TPC_POLYTRACKS";
  std::string m_outputNodeName = "SI_HIT_SEED_EVENT";
  unsigned long long m_event = 0;
  SiHitSeedEvent* m_eventData = nullptr;

  int m_minLayers = 3;
  int m_minMvtxLayers = 2;
  int m_minInttLayers = 0;
  bool m_refitWithIntt = false;

  // silicon (hardware) alignment
  bool m_applyAlignment = false;
  std::string m_alignmentFile;
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
  double m_anchorZ = 3.0;     // [cm]
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

  // selection
  double m_sigRphi = 0.1;                       // [cm]
  std::array<double, 2> m_sigZ{{0.2, 0.6}};      // [cm] MVTX, INTT
  bool m_freeZOffset = true;
  double m_maxTpcChi2 = 0.0;
  unsigned int m_maxChainsPerTrack = 32;
  bool m_storeSelectedOnly = false;

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
