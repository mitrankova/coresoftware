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
// Matches silicon trajectories (SI_TRAJECTORY, SiTrajectoryFitter) to TPC poly tracks
// (TPC_POLYTRACKS + their TPC_POLYCLUSTERS) in phi and eta, and refits every matched pair
// with all its points.  Output: SiTpc_TrackContainerv1 "SITPC_TRACKS" under DST/SVTX.
//
// Frame: beam-axis frame (SiTpcBeamAlignment).  The Si trajectories are already in it (use the
// same alignment object in SiTrajectoryFitter); the TPC clusters are moved into it here.
//
// 1. TPC-only helix: SiTpcHelixFit on the track's clusters (same model as the Si fit), so
//    both sides are described the same way at their point of closest approach to the beam.
// 2. Matching variables, all at the pca:
//      dphi = phi_Si - phi_TPC (direction of motion, wrapped),  deta = eta_Si - eta_TPC,
//      dz0  = z0_Si - z0_TPC (optional: TPC z depends on the bunch crossing / t0).
//    A pair is a candidate if |dphi| < dphiWindow and |deta| < detaWindow [and |dz0| <
//    dz0Window, and same charge if requested and both are circle fits].
//    chi2 = (dphi/dphiWindow)^2 + (deta/detaWindow)^2 [+ (dz0/dz0Window)^2].
//    Assignment: best chi2 first, every Si and TPC track used at most once.
// 3. Refit: SiTpcHelixFit on Si layer points + TPC clusters with per-detector weights (xy and z
//    separately; INTT z is a strip length, so its z weight is small by default), optional free
//    TPC z offset, optional beam constraint.
//
// All pairs within setCandidateWindowScale x the windows are stored in the container (with an
// "accepted" flag) to tune the windows, see SiTpc_TrackResiduals.
class SiTpcTrackMatcher : public SubsysReco
{
 public:
  explicit SiTpcTrackMatcher(const std::string& name = "SiTpcTrackMatcher");
  ~SiTpcTrackMatcher() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  // ---- nodes
  void setSiTrajectoryNodeName(const std::string& s) { m_siTrajNodeName = s; }
  void setTpcClusterNodeName(const std::string& s) { m_tpcClusterNodeName = s; }
  void setTpcTrackNodeName(const std::string& s) { m_tpcTrackNodeName = s; }
  void setOutputNodeName(const std::string& s) { m_outputNodeName = s; }

  // ---- frame (same object as in SiTrajectoryFitter / SiTpc_EventDisplay)
  void setBeamAlignment(const SiTpcBeamAlignment& a) { m_alignment = a; }
  void setApplyBeamAlignment(bool v) { m_alignment.setEnabled(v); }
  void setBeamPositionCm(double x, double y)
  {
    m_beamX = x;
    m_beamY = y;
  }
  void setBz(double tesla) { m_bz = tesla; }

  // ---- matching
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
  void setBeamConstraint(bool on, double weight = 1.0)
  {
    m_beamConstraint = on;
    m_beamWeight = weight;
  }

 private:
  int createNodes(PHCompositeNode*);

  std::string m_siTrajNodeName = "SI_TRAJECTORY";
  std::string m_tpcClusterNodeName = "TPC_POLYCLUSTERS";
  std::string m_tpcTrackNodeName = "TPC_POLYTRACKS";
  std::string m_outputNodeName = "SITPC_TRACKS";
  SiTpc_TrackContainer* m_container = nullptr;

  SiTpcBeamAlignment m_alignment;
  double m_beamX = 0.0;
  double m_beamY = 0.0;
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
  bool m_beamConstraint = false;
  double m_beamWeight = 1.0;

  unsigned long m_nEvents = 0;
  unsigned long m_nMatched = 0;
};

#endif
