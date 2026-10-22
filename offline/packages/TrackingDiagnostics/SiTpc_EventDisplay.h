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
class SiTpc_TrackContainer;
class TpcCrossingDecisionContainer;

// Event display of the silicon seeds / trajectories together with the TPC poly clusters and
// poly tracks.  One TDirectory per event (events/event_NNNNNN), only canvases are written.
//
// Inputs (each optional, the display draws whatever is present):
//   SI_HIT_SEED_EVENT      SiHitSeedEvent (transient: SiHitSeedReco must run in the same job)
//   SI_TRAJECTORY          Si_TrajectoryContainer (SiTrajectoryFitter)
//   TPC_POLYCLUSTERS       Tpc_PolyClusterContainer
//   TPC_POLYTRACKS         Tpc_PolyTrackContainer
//   TPC_POLYTRACKVERTICES  Tpc_PolyTrackVertexContainer (pca + collision vertices, z0 selection)
//   TPC_CROSSING_DECISIONS TpcCrossingDecisionContainer (crossing of every assembled track;
//                          without it every TPC track is in crossing 0)
//
// ---- DETECTOR coordinates (as SiHitSeedDisplay): axes nominal Vz, intrinsic Uphi, layer.
//   c3_evtN_det_hits_chains_traj : one histogram (BOX2Z) with all Si hits over the full Uphi
//                                  turn [B, B+2pi) (B = setClamshellBoundary), every TPC-driven
//                                  chain selected by SiHitSeedReco as a solid line (candidates
//                                  dotted with setDrawCandidateChains), every Si trajectory as a
//                                  dashed line (crossing of the fitted helix with each nominal
//                                  layer), and the vertex plane at layer -1 (R = 0): the pca z of
//                                  each chain's TPC seed (triangle) and each trajectory's z0
//                                  (dotted, chain colour).
//   Lines crossing the Uphi seam at B are split there.
//   setDrawUnassociatedHits(false) keeps only the chain hits in the histogram.
//
// ---- PHYSICAL (phi, eta, r) coordinates of the hits: same content as det_hits_chains_traj.
//   c3_evtN_phys_phi_eta_r_hits_chains_traj : one histogram (BOX2Z) with all Si hits at their
//                                  position in the TPC frame (same transform as the physics views),
//                                  axes phi, eta, r.  Every TPC-driven chain as a solid line through
//                                  its hits (candidates dotted with setDrawCandidateChains), every Si
//                                  trajectory as a dashed line (helix sampled from
//                                  setPhysViewMinRadius to the edge of the silicon), chains without a
//                                  usable fit marked with open black diamonds.
//     r   = distance from the beam axis (beam marker position, see setBeamPositionCm)
//     phi = azimuth around the beam axis, axis [P, P+2pi) (P = setPhysPhiBoundary, default -pi)
//     eta = asinh((z - z_ref) / r), z_ref = z of the first TPC collision vertex
//           (setEtaReferenceFromVertex, default) or setEtaReferenceZ (also the fallback).
//           With z_ref at the collision, a primary track lies at constant eta.
//   setDrawPhysHitsView(false) switches it off; setDrawUnassociatedHits applies as well.
//
// ---- PHYSICS space (as Tpc_PolyClusterDisplay): axes z, x, y in the global frame [cm].
//   c3_evtN_phys_z_x_y_crossingC   one canvas per bunch crossing C (negative: crossing_mC):
//                           the TPC poly clusters + poly tracks + pca of the tracks resolved to
//                           crossing C (TPC_CROSSING_DECISIONS), the Si chains / trajectories /
//                           Si+TPC tracks of those tracks (Si hits not in any chain are not drawn).
//                           Collision vertices only if the event has a single crossing.
//                           Each trajectory is solid in its own detector; setTrajectoryExtension
//                           selects which ones are also drawn (dashed) through the other one.
//   c3_evtN_phys_z_x_y_allcrossings   the same, all crossings in one canvas
//   c3_evtN_phys_si_z_x_y   zoom on the silicon, all crossings, with the event summary and all
//                           unassociated Si hits (grey)
//   Optional (setDrawSiTpcTracks): fitted Si+TPC tracks as solid lines in all physics views.
//
// Coordinates of the physics views: TPC frame (SiTpcBeamAlignment), beam on the TPC beam axis.
//   TPC: clusters, poly tracks, pca and collision vertices as reconstructed (reference).
//   Si : hit (Uphi, Vz) -> detector frame -> detector shift (SiDetectorFrame) -> its clamshell
//        half moved onto the TPC beam axis.  Keep the same settings as SiTrajectoryFitter.
// setApplyBeamAlignment(false) shows the Si in the global frame (only the detector shift).
// Colours: Si seed and its trajectory share a colour (chain id); a TPC poly cluster and its
// poly track share a colour (assembled track id).  Si hits are squares, TPC clusters circles.
//
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
  void setCrossingDecisionNodeName(const std::string& s) { m_crossingNodeName = s; }

  // ---- what to draw
  void setDrawDetectorView(bool v) { m_drawDetector = v; }
  void setDrawPhysicsView(bool v) { m_drawPhysics = v; }
  void setDrawUnassociatedHits(bool v) { m_drawUnassociated = v; }
  void setDrawTpc(bool v) { m_drawTpc = v; }
  void setMaxEventDisplays(unsigned int n) { m_maxEventDisplays = n; }
  // Only save events with at least one Si trajectory / one TPC poly track.
  void setRequireSiTrajectory(bool v) { m_requireSiTraj = v; }
  void setRequireTpcTrack(bool v) { m_requireTpcTrack = v; }
  // Fitted Si+TPC tracks (SiTpcTrackMatcher, node SITPC_TRACKS) as solid lines in the physics
  // views, from the pca to the edge of the view (default off; colour default black).
  // Note: with a free TPC z offset in the refit, the drawn TPC clusters are at their measured
  // z, i.e. displaced in z from the line by the fitted offset.
  void setDrawSiTpcTracks(bool v) { m_drawSiTpcTracks = v; }
  void setSiTpcTrackNodeName(const std::string& s) { m_siTpcTrackNodeName = s; }
  void setSiTpcTrackColor(int c) { m_siTpcTrackColor = c; }
  // With SITPC_TRACKS present: TPC clusters / tracks matched to a Si trajectory are drawn in the
  // colour of that Si chain, unmatched TPC in grey (default).  false: colour by assembled id.
  void setColorTpcByMatch(bool v) { m_colorTpcByMatch = v; }
  // Colour only Si chains with at least this many MVTX layers, and their Si trajectories and the
  // TPC tracks matched to them; the others are grey like unmatched TPC tracks (default 0: every
  // chain in colour). Si+TPC tracks keep setSiTpcTrackColor.
  void setColorMinMvtxLayers(int n) { m_colorMinMvtx = n; }
  // false: grey objects (unmatched TPC tracks, chains below setColorMinMvtxLayers and what is
  // matched to them, incl. their Si+TPC tracks) are not drawn in the physics views (default true)
  void setDrawGreyTracks(bool v) { m_drawGreyTracks = v; }

  // ---- TPC-driven Si chains (SiHitSeedReco): only the selected chain of every TPC track
  // (default) or also all candidate chains (dotted in the detector view)
  void setDrawCandidateChains(bool v) { m_drawCandidateChains = v; }

  // ---- Si frame (must match SiTrajectoryFitter)
  void setUphiRotation(double rad) { m_frame.setUphiRotation(rad); }
  void setDetectorCenterMm(double x, double y, double z = 0.0) { m_frame.setDetectorCenterMm(x, y, z); }
  // INTT centre (layers 3-6); not set: the INTT uses the detector (MVTX) centre
  void setInttCenterMm(double x, double y, double z = 0.0) { m_frame.setInttCenterMm(x, y, z); }
  // Beam marker position.  Default: TPC beam axis at z = 0 from the alignment object.
  void setBeamPositionCm(double x, double y)
  {
    m_beamX = x;
    m_beamY = y;
    m_beamFromAlignment = false;
  }
  // Si clamshell halves moved onto the TPC beam axis (SiTpcBeamAlignment, default: values of
  // the vertex QA PDFs).  Must be the same object as in SiTrajectoryFitter: the Si
  // trajectories are fitted in the TPC frame.
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
  // Radius [cm] at which Si trajectories end (ExtendSiIntoTpc / ExtendBoth): solid up to the
  // silicon (13 cm) or to this radius if smaller, dashed beyond.  Without Si extension they
  // end at 13 cm.
  void setSiExtrapolationRadius(double r) { m_siExtrapR = r; }
  // Radius [cm] down to which TPC tracks are extrapolated (ExtendTpcIntoSi / ExtendBoth);
  // the extension also stops at the track's closest approach to the beam axis.
  void setTpcExtrapolationRadius(double r) { m_tpcExtrapR = r; }

  // ---- detector view
  void setZBins(int n) { m_nz = n; }
  void setPhiBins(int n) { m_nphi = n; }  // full circle
  // Uphi [rad] where the Uphi axis starts (the seam of the display).
  void setClamshellBoundary(double uphi) { m_boundary = uphi; }

  // ---- physical (phi, eta, r) view of the Si hits / chains / trajectories
  void setDrawPhysHitsView(bool v) { m_drawPhysHits = v; }
  // phi [rad] where the phi axis starts (seam of the display), default -pi
  void setPhysPhiBoundary(double phi) { m_physPhiBoundary = phi; }
  // eta = asinh((z - z_ref) / r): z_ref from the first TPC collision vertex if present (default),
  // otherwise / with setEtaReferenceFromVertex(false) the fixed value setEtaReferenceZ [cm]
  void setEtaReferenceZ(double z) { m_etaRefZ = z; }
  void setEtaReferenceFromVertex(bool v) { m_etaRefFromVertex = v; }
  void setEtaBins(int n) { m_neta = n; }
  void setRadiusBins(int n) { m_nr = n; }
  // Si trajectories are drawn from this radius [cm] outward (eta is ill defined near r = 0)
  void setPhysViewMinRadius(double r) { m_physMinR = r; }

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
  std::string m_crossingNodeName = "TPC_CROSSING_DECISIONS";
  int m_colorMinMvtx = 0;
  bool m_drawGreyTracks = true;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt = 0;
  unsigned int m_eventsSaved = 0;

  TFile* m_outfile = nullptr;
  SiHitSeedEvent* m_siSeeds = nullptr;
  Si_TrajectoryContainer* m_siTraj = nullptr;
  Tpc_PolyClusterContainer* m_tpcClusters = nullptr;
  Tpc_PolyTrackContainer* m_tpcTracks = nullptr;
  Tpc_PolyTrackVertexContainer* m_tpcVertices = nullptr;
  SiTpc_TrackContainer* m_siTpcTracks = nullptr;
  TpcCrossingDecisionContainer* m_decisions = nullptr;
  std::string m_siTpcTrackNodeName = "SITPC_TRACKS";
  bool m_drawSiTpcTracks = false;
  bool m_colorTpcByMatch = true;
  int m_siTpcTrackColor = 1;  // kBlack
  bool m_drawCandidateChains = false;

  bool m_drawDetector = true;
  bool m_drawPhysics = true;
  bool m_drawUnassociated = true;
  bool m_drawTpc = true;
  bool m_requireSiTraj = false;
  bool m_requireTpcTrack = false;

  SiDetectorFrame m_frame;
  SiTpcBeamAlignment m_alignment;
  double m_beamX = 0.0;
  double m_beamY = 0.0;
  bool m_beamFromAlignment = true;
  double beamX() const { return m_beamFromAlignment ? m_alignment.referenceBeamX() : m_beamX; }
  double beamY() const { return m_beamFromAlignment ? m_alignment.referenceBeamY() : m_beamY; }

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

  bool m_drawPhysHits = true;
  double m_physPhiBoundary = -3.14159265358979323846;
  double m_etaRefZ = 0.0;
  bool m_etaRefFromVertex = true;
  int m_neta = 100;
  int m_nr = 52;
  double m_physMinR = 2.0;

  double m_zmin = -102.0;
  double m_zmax = 102.0;
  double m_xymax = 85.0;
  double m_trackVertexZMin = -20.0;
  double m_trackVertexZMax = 20.0;
  double m_bz = 1.4;
  bool m_useStraightLineTracks = false;
};

#endif