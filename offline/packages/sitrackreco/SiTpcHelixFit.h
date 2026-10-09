// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCHELIXFIT_H
#define SITRACKRECO_SITPCHELIXFIT_H

#include <cmath>
#include <vector>

// Helix fit shared by SiTrajectoryFitter (silicon only), SiTpcTrackMatcher (TPC only, and
// the combined Si + TPC refit) and the residual QA, so every track uses the same model.
//
//   xy : weighted Taubin circle fit, or a straight line if the radial lever arm of the points
//        is shorter than Config::minLeverArm (curvature not measurable);
//        optional beam-spot constraint (beam position as an extra point with beamWeight)
//   pca: point of closest approach of the circle to the beam position
//   z  : weighted straight-line fit z(s) = z0 + tanl * s, s = transverse path length from
//        the pca (s > 0 outward); optionally a free z offset for the points of group 1
//        (TPC: its z depends on the assumed bunch crossing / t0)
//
// Conventions (same as Si_Trajectory): helicity +1 = counter-clockwise going outward;
// charge = -helicity * sign(Bz); phi = direction of motion at the pca; dca signed, positive
// if the beam is to the right of the direction of motion.  Units cm, rad, GeV, T.
namespace SiTpcHelixFit
{
  struct Point
  {
    double x = 0;
    double y = 0;
    double z = 0;
    double wxy = 1;  // weight in the circle fit (0 = not used in xy)
    double wz = 1;   // weight in the z(s) fit  (0 = not used in z)
    int group = 0;   // 0 = silicon, 1 = TPC (gets the optional z offset)
  };

  struct Config
  {
    double beamX = 0;
    double beamY = 0;
    double bz = 1.4;
    double minLeverArm = 0;    // [cm] radial span below which a straight line is fitted
    double beamWeight = 0;     // > 0: beam position added to the xy fit with this weight
    bool groupZOffset = false;  // free z offset for group 1
  };

  // fit status values, identical to Si_Trajectory::FitStatus
  enum Status
  {
    NotFitted = 0,
    Ok = 1,
    TooFewPoints = 2,
    Degenerate = 3,
    StraightLine = 4
  };

  // radius used to store a straight line in the circle parametrization [cm]
  constexpr double kStraightRadius = 1.0e6;

  struct Result
  {
    int status = NotFitted;
    double cx = 0, cy = 0, R = 0;  // circle (straight line: very large R)
    int helicity = 0;
    double circleRms = 0;          // rms of the xy distances to the circle (points with wxy > 0)
    double pcaX = 0, pcaY = 0;
    double tx = 0, ty = 0;         // unit direction of motion at the pca
    double phi = 0, dca = 0;
    double pt = 0;                 // NaN for a straight line
    int charge = 0;
    double z0 = 0, tanl = 0;
    double zOffset = 0;            // fitted offset of group 1 (0 if not fitted)
    double zRms = 0;               // rms of the z residuals (points with wz > 0, offset applied)
    std::vector<double> s;         // path length of every input point
    bool fitted() const { return status == Ok || status == StraightLine; }
  };

  struct Circle
  {
    double x = 0, y = 0, r = 0;
    bool ok = false;
  };
  Circle fitCircleTaubin(const std::vector<double>& x, const std::vector<double>& y,
                         const std::vector<double>& w);

  // pts ordered outward (the helicity is taken from the first and last point)
  Result fit(const std::vector<Point>& pts, const Config& cfg, unsigned int minPoints = 3);

  struct Pos
  {
    double x = 0, y = 0, z = 0;
  };
  // point on the fitted helix at path length s (group 1: add zOffset yourself)
  Pos at(const Result& r, double s);
  // crossing of the helix with the circle of radius rad around (0, 0); of the (up to two)
  // crossings within half a turn the one with z closest to zRef is returned
  bool atRadius(const Result& r, double rad, double zRef, Pos& out, double& sOut);

  // Result rebuilt from a stored track (Si_Trajectory, SiTpc_Track: same getters)
  template <class Track>
  Result fromTrack(const Track& t)
  {
    Result r;
    r.status = t.get_fit_status();
    r.cx = t.get_circle_x();
    r.cy = t.get_circle_y();
    r.R = t.get_radius();
    r.helicity = t.get_helicity();
    r.circleRms = t.get_circle_rms();
    r.pcaX = t.get_pca_x();
    r.pcaY = t.get_pca_y();
    r.phi = t.get_phi();
    r.tx = std::cos(r.phi);
    r.ty = std::sin(r.phi);
    r.dca = t.get_dca();
    r.pt = t.get_pt();
    r.charge = t.get_charge();
    r.z0 = t.get_z0();
    r.tanl = t.get_tanl();
    r.zRms = t.get_z_rms();
    return r;
  }
}  // namespace SiTpcHelixFit

#endif
