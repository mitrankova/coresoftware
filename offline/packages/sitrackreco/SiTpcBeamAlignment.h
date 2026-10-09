// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCBEAMALIGNMENT_H
#define SITRACKRECO_SITPCBEAMALIGNMENT_H

#include "SiDetectorFrame.h"

#include <array>
#include <cmath>

// Alignment of the TPC and of the two silicon clamshell halves with respect to the beam axis.
//
// Each part has its own measured beam line (from the vertex x vs z / y vs z fits):
//     x_beam(z) = x0 + dxdz * z ,   y_beam(z) = y0 + dydz * z
// and the correction moves that line onto the z axis:
//     x' = x - x_beam(z) ,  y' = y - y_beam(z) ,  z' = z
// (shear form; the slopes are ~1 mrad, so this equals the rotation to O(1e-6)).
// After the correction the beam is at (0, 0) for every part, so TPC and Si are in one frame.
//
// Silicon transformation chain, every step separate:
//   1. intrinsic (Uphi, Vz)  -> DETECTOR frame        SiDetectorFrame::toDetector
//   2. DETECTOR frame        -> GLOBAL frame          SiDetectorFrame::toGlobal: the
//                                                     independent +6, -1 mm detector shift
//   3. GLOBAL frame          -> BEAM-AXIS frame       this class, per clamshell half
// The Si beam lines are therefore given in the GLOBAL (translated) frame, i.e. measured
// after step 2 (si_alignment_new3.pdf pages 9-10, summary on the last page).  If the
// detector shift is changed, the Si beam lines have to be re-measured.
// TPC: TPC coordinates -> BEAM-AXIS frame with the TPC beam line.
//
// Clamshell halves as in SiHitSeedQA: half A for Uphi in [B, B+pi), half B for
// [B+pi, B+2pi), B = clamshell boundary in intrinsic Uphi (default 0).
//
// All lengths in cm.  Header only.
class SiTpcBeamAlignment
{
 public:
  enum Part
  {
    Tpc = 0,
    SiHalfA = 1,
    SiHalfB = 2,
    NParts = 3
  };

  struct BeamLine
  {
    double x0 = 0;    // [cm] at z = 0
    double y0 = 0;    // [cm] at z = 0
    double dxdz = 0;  // [rad]
    double dydz = 0;  // [rad]
    double x(double z) const { return x0 + dxdz * z; }
    double y(double z) const { return y0 + dydz * z; }
  };

  using Point = SiDetectorFrame::Point;

  // Values from the vertex QA (linear fits of the Crystal-Ball slice centres):
  //   TPC   : VertexDCAAlignmentQA_All.pdf                       x(z), y(z) [cm]
  //   Si A/B: si_alignment_new3.pdf, "SiHitSeed half A/B, global frame (translated)",
  //           pages 9-10 (summary on the last page)               x(z), y(z) [mm]
  // measured with the Si detector shift (+6, -1) mm and clamshell boundary B = 0.
  SiTpcBeamAlignment()
  {
    setBeamLineCm(Tpc, -0.106835, 0.267482, 8.478616e-04, -8.041254e-04);
    setBeamLineMm(SiHalfA, -0.533300, 1.772693, 7.094553e-04, -2.743455e-04);
    setBeamLineMm(SiHalfB, 0.595174, 1.608226, 2.343532e-03, 6.926381e-04);
  }

  // ---- configuration
  void setEnabled(bool v) { m_enabled = v; }  // false: identity (only the detector shift remains)
  bool enabled() const { return m_enabled; }
  void setBeamLineCm(int part, double x0, double y0, double dxdz, double dydz)
  {
    if (part >= 0 && part < NParts)
    {
      m_line[part] = {x0, y0, dxdz, dydz};
    }
  }
  void setBeamLineMm(int part, double x0, double y0, double dxdz, double dydz)
  {
    setBeamLineCm(part, 0.1 * x0, 0.1 * y0, dxdz, dydz);
  }
  const BeamLine& beamLine(int part) const { return m_line[part]; }
  // Clamshell boundary in intrinsic Uphi [rad], must match the one used for the Si fits.
  void setClamshellBoundary(double uphi) { m_boundary = uphi; }
  double clamshellBoundary() const { return m_boundary; }

  // ---- clamshell half of a Si point
  int siHalf(double uphi) const
  {
    constexpr double twopi = 2.0 * M_PI;
    double u = std::fmod(uphi - m_boundary, twopi);
    if (u < 0)
    {
      u += twopi;
    }
    return u < M_PI ? SiHalfA : SiHalfB;
  }
  // half of a point given in the GLOBAL or BEAM-AXIS frame (Uphi from its azimuth around the
  // detector centre; the alignment moves points by < 1 mm, far from the boundary for hits)
  int siHalfOfGlobal(const SiDetectorFrame& frame, const Point& g) const
  {
    const Point d = frame.toDetectorFromGlobal(g);
    return siHalf(std::atan2(d.y, d.x) - frame.uphiRotation());
  }

  // ---- step 3 and its inverse for one part
  Point toBeamAxis(int part, const Point& p) const
  {
    if (!m_enabled)
    {
      return p;
    }
    const BeamLine& l = m_line[part];
    return {p.x - l.x(p.z), p.y - l.y(p.z), p.z};
  }
  Point fromBeamAxis(int part, const Point& p) const
  {
    if (!m_enabled)
    {
      return p;
    }
    const BeamLine& l = m_line[part];
    return {p.x + l.x(p.z), p.y + l.y(p.z), p.z};
  }

  // ---- full chains
  // TPC point (TPC/global coordinates) -> beam-axis frame
  Point tpcToBeamAxis(const Point& p) const { return toBeamAxis(Tpc, p); }
  Point tpcToBeamAxis(double x, double y, double z) const { return toBeamAxis(Tpc, {x, y, z}); }

  // Si hit on a layer: intrinsic (Uphi, Vz) -> detector -> global (shift) -> beam axis
  Point siToBeamAxis(const SiDetectorFrame& frame, int layer, double uphi, double vz) const
  {
    return toBeamAxis(siHalf(uphi), frame.toGlobal(frame.toDetector(layer, uphi, vz)));
  }
  // Si point already in the global frame (after the shift), half from its own Uphi
  Point siGlobalToBeamAxis(const SiDetectorFrame& frame, const Point& g) const
  {
    return toBeamAxis(siHalfOfGlobal(frame, g), g);
  }
  // Inverse for display: beam-axis frame -> global frame (half from the point's Uphi)
  Point siBeamAxisToGlobal(const SiDetectorFrame& frame, const Point& b) const
  {
    return fromBeamAxis(siHalfOfGlobal(frame, b), b);
  }

 private:
  bool m_enabled = true;
  double m_boundary = 0.0;
  std::array<BeamLine, NParts> m_line{};
};

#endif
