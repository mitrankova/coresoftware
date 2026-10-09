// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITRAJECTORYFITTER_H
#define SITRACKRECO_SITRAJECTORYFITTER_H

#include "SiDetectorFrame.h"
#include "SiTpcBeamAlignment.h"
#include "SiTpcHelixFit.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class SiHitSeedEvent;
class Si_Trajectory;
class Si_TrajectoryContainer;

// SiTrajectoryFitter
//
// Input : SiHitSeedEvent (transient, from SiHitSeedReco), node "SI_HIT_SEED_EVENT".
// Output: Si_TrajectoryContainerv1 of Si_Trajectoryv1, persistent node "SI_TRAJECTORY"
//         under DST/SVTX.
//
// For every chain:
//   1. every hit on the chain -> DETECTOR frame   (SiDetectorFrame::toDetector)
//   2. -> GLOBAL frame: shift by the detector centre (default +6, -1, 0 mm), own step
//   2b. -> BEAM-AXIS frame: per clamshell half, the measured beam line is moved onto the
//       z axis (SiTpcBeamAlignment; setBeamAlignment / setApplyBeamAlignment)
//   3. one fit point per layer = centroid of that layer's hits (projected back to the
//      layer radius, i.e. the azimuth of the centroid on the nominal cylinder)
//   4. xy: Taubin circle fit (or a straight line if the radial lever arm is too short to
//      measure curvature, see setMinCurvatureLeverArm), perigee w.r.t. the beam position;
//      then a straight-line fit z(s).
class SiTrajectoryFitter : public SubsysReco
{
 public:
  explicit SiTrajectoryFitter(const std::string& name = "SiTrajectoryFitter");
  ~SiTrajectoryFitter() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& s) { m_inputNodeName = s; }
  void setOutputNodeName(const std::string& s) { m_outputNodeName = s; }

  // step 1: global phi = Uphi + rotation [rad]
  void setUphiRotation(double rad) { m_frame.setUphiRotation(rad); }
  // step 2: position of the detector-frame origin in the global frame
  void setDetectorCenterMm(double x, double y, double z = 0.0) { m_frame.setDetectorCenterMm(x, y, z); }
  // Beam-axis alignment of the two clamshell halves (SiTpcBeamAlignment, default: values of
  // si_alignment_new3.pdf).  With it on, the beam is at (0, 0): keep setBeamPositionCm(0, 0).
  void setBeamAlignment(const SiTpcBeamAlignment& a) { m_alignment = a; }
  void setApplyBeamAlignment(bool v) { m_alignment.setEnabled(v); }
  const SiTpcBeamAlignment& beamAlignment() const { return m_alignment; }
  // beam position (cm, frame of the fit points) used for the perigee
  void setBeamPositionCm(double x, double y) { m_beamX = x; m_beamY = y; }
  // solenoid field [T], signed: only used for pt and charge
  void setBz(double tesla) { m_bz = tesla; }

  void setUseIntt(bool v) { m_useIntt = v; }
  void setMinPoints(unsigned int n) { m_minPoints = n; }
  // The curvature of a short track is not measurable: three MVTX points span only 1.6 cm, so
  // a 30 um position error already gives a 1 GeV track any pt between ~0.3 and ~2 GeV.
  // If the radial span of the fit points is below this [cm], a straight line is fitted
  // instead (fit status StraightLine, pt = NaN, charge = 0).  Default 3 cm: MVTX-only
  // chains become straight lines, chains reaching the INTT get a circle.  0 = always circle.
  void setMinCurvatureLeverArm(double cm) { m_minCurvatureLeverArm = cm; }
  // Add the beam position (setBeamPositionCm) to the xy fit as an extra point with this
  // weight relative to one layer point.  Only for tracks from the primary vertex.
  void setBeamConstraint(bool on, double weight = 1.0)
  {
    m_beamConstraint = on;
    m_beamWeight = weight;
  }

  // Weights of the MVTX and INTT points in the z(s) fit.  INTT z is the centre of a 1.6 / 2.0 cm
  // long strip (sigma ~0.5 cm), MVTX z is a pixel (~10 um): with equal weights the INTT
  // dominates the eta of the trajectory.  Default 1 / 0.01.
  void setZWeights(double mvtx, double intt)
  {
    m_zWeightMvtx = mvtx;
    m_zWeightIntt = intt;
  }

  // Verbosity 1: per-event count of fits per status.  Verbosity 2: also every chain that was
  // not fitted or got a circle fit with pt below this [GeV] (layers, hits, radii).
  void setReportPt(double gev) { m_reportPt = gev; }

  const SiDetectorFrame& frame() const { return m_frame; }

  // --- fit primitives (public so they can be unit-tested / reused)
  struct Circle
  {
    double x = 0, y = 0, r = 0;
    bool ok = false;
  };
  // Taubin algebraic circle fit (N. Chernov), >= 3 points.
  static Circle fitCircleTaubin(const std::vector<double>& x, const std::vector<double>& y);
  static Circle fitCircleTaubin(const std::vector<double>& x, const std::vector<double>& y,
                                const std::vector<double>& w);
  // radius used to store a straight-line fit in the circle parametrization [cm]
  static constexpr double kStraightRadius = SiTpcHelixFit::kStraightRadius;

 private:
  struct LayerPoint
  {
    int layer = -1;
    SiDetectorFrame::Point global;
    unsigned int nhits = 0;
  };

  int createNodes(PHCompositeNode*);
  std::vector<LayerPoint> layerPoints(const SiHitSeedEvent& event, const std::vector<int>& hitIds) const;
  void fit(const std::vector<LayerPoint>& points, Si_Trajectory& trajectory) const;

  std::string m_inputNodeName = "SI_HIT_SEED_EVENT";
  std::string m_outputNodeName = "SI_TRAJECTORY";
  Si_TrajectoryContainer* m_container = nullptr;

  SiDetectorFrame m_frame;
  SiTpcBeamAlignment m_alignment;
  double m_beamX = 0.0;
  double m_beamY = 0.0;
  double m_bz = 1.4;
  bool m_useIntt = true;
  unsigned int m_minPoints = 3;
  double m_minCurvatureLeverArm = 3.0;
  bool m_beamConstraint = false;
  double m_beamWeight = 1.0;
  double m_reportPt = 0.2;
  double m_zWeightMvtx = 1.0;
  double m_zWeightIntt = 0.01;
};

#endif
