// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCTRACKMATCHER_H
#define SITRACKRECO_SITPCTRACKMATCHER_H

#include "SiTpcBeamAlignment.h"

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class Si_TrajectoryContainer;
class SiTpc_TrackContainer;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrackContainer;

// SiTpcTrackMatcher
//
// Pairs silicon trajectories (SI_TRAJECTORY, SiTrajectoryFitter) with TPC poly tracks
// (TPC_POLYTRACKS + their TPC_POLYCLUSTERS) and refits every pair with all its points (full Si +
// TPC trajectory).  Output: SiTpc_TrackContainerv1 "SITPC_TRACKS" under DST/SVTX.
//
// Frame: TPC frame (SiTpcBeamAlignment).  The TPC clusters are used as measured; the Si
// trajectories are already in it (SiTrajectoryFitter moves both clamshell halves onto the TPC
// beam axis - use the same alignment object there).  Beam = TPC beam axis.
//
// 1. TPC-only helix: SiTpcHelixFit on the track's clusters (same model as the Si fit).
// 2. Association
//    FromSiChains (default): the chaining is TPC-driven (SiHitSeedReco), so every Si chain already
//      belongs to a TPC track (SiHitChain::tpc_track_index, read from SI_HIT_SEED_EVENT, which
//      must be in the same job).  A pair is accepted if its chain was SELECTED by SiHitSeedReco
//      (and, with setApplyMatchWindowToChains(true), also inside the phi/eta windows below).
//      One Si+TPC track per TPC track; Si hits may be shared by different TPC tracks.
//    PhiEtaWindow: independent matching at the pca as before,
//      dphi = phi_Si - phi_TPC, deta = eta_Si - eta_TPC, dz0 = z0_Si - z0_TPC (optional),
//      candidate if |dphi| < dphiWindow and |deta| < detaWindow [and |dz0| < dz0Window, same
//      charge if requested]; chi2 = sum of the squared normalised differences; best chi2 first,
//      every Si and TPC track used at most once.
//    In both modes the matching variables are stored for every pair (candidate list, "accepted"
//    flag) so SiTpc_TrackResiduals can write its matching tree.
// 3. Refit: SiTpcHelixFit on Si layer points + TPC clusters with per-detector weights (xy and z
//    separately; INTT z is a strip length, so its z weight is small by default), optional free
//    TPC z offset, optional beam constraint.
class SiTpcTrackMatcher : public SubsysReco
{
 public:
  explicit SiTpcTrackMatcher(const std::string& name = "SiTpcTrackMatcher");
  ~SiTpcTrackMatcher() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  enum class Association
  {
    FromSiChains,
    PhiEtaWindow
  };
  void setAssociation(Association a) { m_association = a; }
  // FromSiChains: also require the pair to be inside the phi/eta(/z0) windows (default false)
  void setApplyMatchWindowToChains(bool v) { m_windowOnChains = v; }

  // ---- nodes
  void setSiTrajectoryNodeName(const std::string& s) { m_siTrajNodeName = s; }
  void setSiSeedNodeName(const std::string& s) { m_siSeedNodeName = s; }
  void setTpcClusterNodeName(const std::string& s) { m_tpcClusterNodeName = s; }
  void setTpcTrackNodeName(const std::string& s) { m_tpcTrackNodeName = s; }
  void setOutputNodeName(const std::string& s) { m_outputNodeName = s; }

  // ---- frame (same object as in SiTrajectoryFitter / SiTpc_EventDisplay)
  void setBeamAlignment(const SiTpcBeamAlignment& a) { m_alignment = a; }
  void setApplyBeamAlignment(bool v) { m_alignment.setEnabled(v); }
  // Beam position for the perigee / beam constraint.  Default: TPC beam axis at z = 0 from the
  // alignment object.  Calling this overrides it.
  void setBeamPositionCm(double x, double y)
  {
    m_beamX = x;
    m_beamY = y;
    m_beamFromAlignment = false;
  }
  double beamX() const { return m_beamFromAlignment ? m_alignment.referenceBeamX() : m_beamX; }
  double beamY() const { return m_beamFromAlignment ? m_alignment.referenceBeamY() : m_beamY; }
  void setBz(double tesla) { m_bz = tesla; }

  // ---- matching windows (PhiEtaWindow; FromSiChains only with setApplyMatchWindowToChains)
  void setMatchWindow(double dphi, double deta)
  {
    m_dphiWindow = dphi;
    m_detaWindow = deta;
  }
  void setMatchWindowDz0(double cm) { m_dz0Window = cm; }  // <= 0: z0 not used (default)
  void setRequireSameCharge(bool v) { m_sameCharge = v; }   // only if both are circle fits
  void setCandidateWindowScale(double v) { m_candScale = v; }
  void setMinTpcClusters(unsigned int n) { m_minTpcClusters = n; }

  // ---- combined refit
  void setXYWeights(double mvtx, double intt, double tpc)
  {
    m_wxy[0] = mvtx;
    m_wxy[1] = intt;
    m_wxy[2] = tpc;
  }
  void setZWeights(double mvtx, double intt, double tpc)
  {
    m_wz[0] = mvtx;
    m_wz[1] = intt;
    m_wz[2] = tpc;
  }
  void setFitTpcZOffset(bool v) { m_fitTpcZOffset = v; }
  // Hybrid xy: curvature from the circle fit with the xy weights above (with setXYWeights(0, 0, 1)
  // the TPC-only circle: TPC pt unchanged), circle CENTRE refitted with that radius to the silicon
  // points (weights wMvtx, wIntt ~ 1/sigma_rphi^2 relative to the MVTX: MVTX ~5 um, INTT ~25 um ->
  // 0.04).  The Si+TPC dca, pca, phi and z0 then have silicon precision; the TPC stays free in
  // position (its distortions / offsets do not move the track at the beam).  Needs at least
  // minPoints silicon points with weight > 0; otherwise the plain fit is kept.
  void setHybridXY(bool on, double wMvtx = 1.0, double wIntt = 0.04, unsigned int minPoints = 2)
  {
    m_hybrid = on;
    m_hybridW[0] = wMvtx;
    m_hybridW[1] = wIntt;
    m_hybridMinPoints = minPoints;
  }
  void setBeamConstraint(bool on, double weight = 1.0)
  {
    m_beamConstraint = on;
    m_beamWeight = weight;
  }

 private:
  int createNodes(PHCompositeNode*);

  Association m_association = Association::FromSiChains;
  bool m_windowOnChains = false;

  std::string m_siTrajNodeName = "SI_TRAJECTORY";
  std::string m_siSeedNodeName = "SI_HIT_SEED_EVENT";
  std::string m_tpcClusterNodeName = "TPC_POLYCLUSTERS";
  std::string m_tpcTrackNodeName = "TPC_POLYTRACKS";
  std::string m_outputNodeName = "SITPC_TRACKS";
  SiTpc_TrackContainer* m_container = nullptr;

  SiTpcBeamAlignment m_alignment;
  double m_beamX = 0.0;
  double m_beamY = 0.0;
  bool m_beamFromAlignment = true;
  double m_bz = 1.4;

  double m_dphiWindow = 0.05;
  double m_detaWindow = 0.05;
  double m_dz0Window = 0.0;
  bool m_sameCharge = false;
  double m_candScale = 3.0;
  unsigned int m_minTpcClusters = 10;

  double m_wxy[3] = {10.0, 10.0, 1.0};  // MVTX, INTT, TPC
  double m_wz[3] = {10.0, 0.01, 1.0};
  bool m_fitTpcZOffset = true;
  bool m_hybrid = false;
  double m_hybridW[2] = {1.0, 0.04};  // MVTX, INTT
  unsigned int m_hybridMinPoints = 2;
  unsigned long m_nFitted = 0;
  unsigned long m_nHybrid = 0;
  bool m_beamConstraint = false;
  double m_beamWeight = 1.0;

  unsigned long m_nEvents = 0;
  unsigned long m_nMatched = 0;
};

#endif
