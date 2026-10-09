#ifndef SITPC_EVENTDISPLAY_H
#define SITPC_EVENTDISPLAY_H

#include <sitrackreco/SiDetectorFrame.h>
#include <sitrackreco/SiTpcBeamAlignment.h>

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class SiHitSeedEvent;
class Si_TrajectoryContainer;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrackContainer;
class Tpc_PolyTrackVertexContainer;

// Event display of the silicon seeds / trajectories together with the TPC poly clusters and
// poly tracks.  One TDirectory per event (events/event_NNNNNN), only canvases are written.
//
// Inputs (each optional, the display draws whatever is present):
//   SI_HIT_SEED_EVENT      SiHitSeedEvent (transient: SiHitSeedReco must run in the same job)
//   SI_TRAJECTORY          Si_TrajectoryContainer (SiTrajectoryFitter)
//   TPC_POLYCLUSTERS       Tpc_PolyClusterContainer
//   TPC_POLYTRACKS         Tpc_PolyTrackContainer
//   TPC_POLYTRACKVERTICES  Tpc_PolyTrackVertexContainer (pca + collision vertices, z0 selection)
//
// ---- DETECTOR coordinates (as SiHitSeedDisplay): axes nominal Vz, intrinsic Uphi, layer.
//   c3_evtN_det_hits_chains_traj : one histogram (BOX2Z) with all Si hits over the full Uphi
//                                  turn [B, B+2pi) (B = setClamshellBoundary), every chain as a
//                                  solid line, every Si trajectory as a dashed line (crossing of
//                                  the fitted helix with each nominal layer), and the vertex
//                                  plane at layer -1 (R = 0): tracklet vertex_z (blue) and each
//                                  trajectory's z0 (dotted, chain colour).
//   Lines crossing the Uphi seam at B are split there.
//   setDrawUnassociatedHits(false) keeps only the chain hits in the histogram.
//
// ---- PHYSICS space (as Tpc_PolyClusterDisplay): axes z, x, y in the global frame [cm].
//   c3_evtN_phys_z_x_y      TPC poly clusters + poly tracks + pca / collision vertices,
//                           Si hits (chain colour; others grey), Si trajectories.
//                           Each trajectory is solid in its own detector; setTrajectoryExtension
//                           selects which ones are also drawn (dashed) through the other one.
//   c3_evtN_phys_si_z_x_y   zoom on the silicon
//   c_evtN_phys_xy          x-y projection, everything
//   c_evtN_phys_xy_si       x-y zoom on the silicon with nominal layers around the detector
//                           centre, beam and detector-centre markers
//   c_evtN_phys_z_r         z vs radius, everything
//
// Coordinates of the physics views: BEAM-AXIS frame (SiTpcBeamAlignment), beam at (0, 0).
//   Si : hit (Uphi, Vz) -> detector frame -> detector shift (SiDetectorFrame) -> beam axis of
//        its clamshell half.  Keep the same settings as SiTrajectoryFitter.
//   TPC: clusters, poly tracks, pca and collision vertices -> beam axis with the TPC beam line.
// setApplyBeamAlignment(false) shows the old frames (Si global frame, TPC as reconstructed).
// Colours: Si seed and its trajectory share a colour (chain id); a TPC poly cluster and its
// poly track share a colour (assembled track id).  Si hits are squares, TPC clusters circles.
class SiTpc_EventDisplay : public SubsysReco
{
 public:
  SiTpc_EventDisplay(const std::string& name = "SiTpc_EventDisplay",
                     const std::string& outfilename = "sitpc_event_display.root",
                     unsigned int maxEventDisplays = 5);
  ~SiTpc_EventDisplay() override;

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // ---- node names
  void setSiSeedNodeName(const std::string& s) { m_siSeedNodeName = s; }
  void setSiTrajectoryNodeName(const std::string& s) { m_siTrajNodeName = s; }
  void setTpcClusterNodeName(const std::string& s) { m_tpcClusterNodeName = s; }
  void setTpcTrackNodeName(const std::string& s) { m_tpcTrackNodeName = s; }
  void setTpcVertexNodeName(const std::string& s) { m_tpcVertexNodeName = s; }

  // ---- what to draw
  void setDrawDetectorView(bool v) { m_drawDetector = v; }
  void setDrawPhysicsView(bool v) { m_drawPhysics = v; }
  void setDrawProjections(bool v) { m_drawProjections = v; }  // 2D x-y and z-r canvases
  void setDrawUnassociatedHits(bool v) { m_drawUnassociated = v; }
  void setDrawTpc(bool v) { m_drawTpc = v; }
  void setMaxEventDisplays(unsigned int n) { m_maxEventDisplays = n; }
  // Only save events with at least one Si trajectory / one TPC poly track.
  void setRequireSiTrajectory(bool v) { m_requireSiTraj = v; }
  void setRequireTpcTrack(bool v) { m_requireTpcTrack = v; }

  // ---- Si frame (must match SiTrajectoryFitter)
  void setUphiRotation(double rad) { m_frame.setUphiRotation(rad); }
  void setDetectorCenterMm(double x, double y, double z = 0.0) { m_frame.setDetectorCenterMm(x, y, z); }
  void setBeamPositionCm(double x, double y) { m_beamX = x; m_beamY = y; }
  // Beam-axis alignment of TPC and the two Si clamshell halves (SiTpcBeamAlignment, default:
  // values of the vertex QA PDFs).  Must be the same object as in SiTrajectoryFitter: the Si
  // trajectories are fitted in the beam-axis frame.  With it on, the beam is at (0, 0).
  void setBeamAlignment(const SiTpcBeamAlignment& a) { m_alignment = a; }
  void setApplyBeamAlignment(bool v) { m_alignment.setEnabled(v); }

  // ---- How far trajectories are drawn in the physics views (TPC + Si together).
  // The own-detector part is always solid; the extension into the other detector is dashed.
  //   ExtendNone      : Si only inside the silicon, TPC only over its clusters
  //   ExtendSiIntoTpc : Si trajectories extrapolated out through the TPC (default)
  //   ExtendTpcIntoSi : TPC poly tracks extrapolated in through the silicon to their pca
  //   ExtendBoth      : both
  enum TrajectoryExtension
  {
    ExtendNone = 0,
    ExtendSiIntoTpc = 1,
    ExtendTpcIntoSi = 2,
    ExtendBoth = 3
  };
  void setTrajectoryExtension(int mode) { m_extension = mode; }

  // ---- Si trajectory selection
  void setSiOnlyGoodFits(bool v) { m_siOnlyGood = v; }
  // Circle fits below this pt are drawn dotted (or hidden with setDrawLowPtSi(false)).
  // Straight-line fits (no pt) are never cut.  Chains without a usable fit are marked with
  // open black diamonds on their hits; Verbosity 1 prints the reason for each of them.
  void setSiMinPt(double gev) { m_siMinPt = gev; }
  void setDrawLowPtSi(bool v) { m_drawLowPtSi = v; }
  void setSiMinPoints(unsigned int n) { m_siMinPoints = n; }
  // Radius [cm] up to which Si trajectories are extrapolated (ExtendSiIntoTpc / ExtendBoth).
  void setSiExtrapolationRadius(double r) { m_siExtrapR = r; }
  // Radius [cm] down to which TPC tracks are extrapolated (ExtendTpcIntoSi / ExtendBoth);
  // the extension also stops at the track's closest approach to the beam axis.
  void setTpcExtrapolationRadius(double r) { m_tpcExtrapR = r; }

  // ---- detector view
  void setZBins(int n) { m_nz = n; }
  void setPhiBins(int n) { m_nphi = n; }  // full circle
  // Uphi [rad] where the Uphi axis starts (the seam of the display).
  void setClamshellBoundary(double uphi) { m_boundary = uphi; }

  // ---- physics view (as Tpc_PolyClusterDisplay)
  void setZRange(double zmin, double zmax)
  {
    m_zmin = zmin;
    m_zmax = zmax;
  }
  void setXYRange(double xymax) { m_xymax = xymax; }
  // TPC poly tracks (and their clusters) are drawn only if their vertex z0 is in this range
  // (applied when TPC_POLYTRACKVERTICES exists).
  void setTrackVertexZRange(double zmin, double zmax)
  {
    m_trackVertexZMin = zmin;
    m_trackVertexZMax = zmax;
  }
  void setMagneticFieldTesla(double b) { m_bz = b; }
  void setUseStraightLineTracks(bool v) { m_useStraightLineTracks = v; }

 private:
  void get_nodes(PHCompositeNode* topNode);

  std::string m_outfilename;
  std::string m_siSeedNodeName = "SI_HIT_SEED_EVENT";
  std::string m_siTrajNodeName = "SI_TRAJECTORY";
  std::string m_tpcClusterNodeName = "TPC_POLYCLUSTERS";
  std::string m_tpcTrackNodeName = "TPC_POLYTRACKS";
  std::string m_tpcVertexNodeName = "TPC_POLYTRACKVERTICES";
  unsigned int m_maxEventDisplays;
  unsigned int m_evt = 0;
  unsigned int m_eventsSaved = 0;

  TFile* m_outfile = nullptr;
  SiHitSeedEvent* m_siSeeds = nullptr;
  Si_TrajectoryContainer* m_siTraj = nullptr;
  Tpc_PolyClusterContainer* m_tpcClusters = nullptr;
  Tpc_PolyTrackContainer* m_tpcTracks = nullptr;
  Tpc_PolyTrackVertexContainer* m_tpcVertices = nullptr;

  bool m_drawDetector = true;
  bool m_drawPhysics = true;
  bool m_drawProjections = true;
  bool m_drawUnassociated = true;
  bool m_drawTpc = true;
  bool m_requireSiTraj = false;
  bool m_requireTpcTrack = false;

  SiDetectorFrame m_frame;
  SiTpcBeamAlignment m_alignment;
  double m_beamX = 0.0;
  double m_beamY = 0.0;

  bool m_siOnlyGood = true;
  double m_siMinPt = 0.0;
  bool m_drawLowPtSi = true;
  unsigned int m_siMinPoints = 3;
  double m_siExtrapR = 80.0;
  double m_tpcExtrapR = 0.0;
  int m_extension = ExtendSiIntoTpc;

  int m_nz = 120;
  int m_nphi = 180;
  double m_boundary = 0.0;

  double m_zmin = -102.0;
  double m_zmax = 102.0;
  double m_xymax = 85.0;
  double m_trackVertexZMin = -20.0;
  double m_trackVertexZMax = 20.0;
  double m_bz = 1.4;
  bool m_useStraightLineTracks = false;
};

#endif
