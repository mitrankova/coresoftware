// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCBEAMALIGNMENT_H
#define SITRACKRECO_SITPCBEAMALIGNMENT_H

#include "SiDetectorFrame.h"

#include <array>
#include <cmath>

// Beam-axis alignment of the two silicon clamshell halves to the TPC.
//
// ---- 1. The beam axes (measured, one per part)
// Every part has its own beam axis, measured from its vertices (linear fits of the vertex
// x vs z and y vs z):
//     x_axis(z) = x0 + dxdz * z ,   y_axis(z) = y0 + dydz * z
//   TPC      : VertexDCAAlignmentQA_All.pdf                                  [cm]
//   Si half A: si_alignment_new3.pdf, half A, global frame (translated), p. 9 [mm]
//   Si half B: si_alignment_new3.pdf, half B, global frame (translated), p.10 [mm]
// The TPC beam axis is the REFERENCE.
//
// ---- 2. The transformation (Si half -> TPC beam axis)
// Each Si clamshell half is moved so that its beam axis coincides with the TPC beam axis:
//     dx(z) = x_half(z) - x_tpc(z) ,   dy(z) = y_half(z) - y_tpc(z)
//     x' = x - dx(z) ,  y' = y - dy(z) ,  z' = z
// (shear form; the slope differences are ~1 mrad, equal to the rotation to O(1e-6)).
// The TPC is not transformed (identity): the common frame is the TPC frame, the beam is on
// the TPC beam axis (x_tpc(z), y_tpc(z)), NOT at (0, 0).  Fits use the TPC beam axis as
// beam position (referenceBeamX/Y).
//
// Silicon chain, every step separate:
//   a. intrinsic (Uphi, Vz) -> DETECTOR frame   SiDetectorFrame::toDetector
//   b. DETECTOR frame       -> GLOBAL frame     SiDetectorFrame::toGlobal (+6, -1 mm shift)
//   c. GLOBAL frame         -> TPC frame        this class: half A or B onto the TPC beam axis
// The Si beam axes are measured in the GLOBAL frame (after step b).  If the detector shift
// is changed, they have to be re-measured.
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

  SiTpcBeamAlignment()
  {
    // ---- 1. beam axes (defaults = vertex QA PDFs, Si shift (+6, -1) mm, boundary B = 0)
    //                    x0          y0          dx/dz          dy/dz
    setTpcBeamAxisCm(-0.106835, 0.267482, 8.478616e-04, -8.041254e-04);
    setSiBeamAxisMm(SiHalfA, -0.533300, 1.772693, 7.094553e-04, -2.743455e-04);
    setSiBeamAxisMm(SiHalfB, 0.595174, 1.608226, 2.343532e-03, 6.926381e-04);
  }

  // =========================================================================================
  // 1. beam axes
  // =========================================================================================
  void setTpcBeamAxisCm(double x0, double y0, double dxdz, double dydz)
  {
    m_axis[Tpc] = {x0, y0, dxdz, dydz};
  }
  void setTpcBeamAxisMm(double x0, double y0, double dxdz, double dydz)
  {
    setTpcBeamAxisCm(0.1 * x0, 0.1 * y0, dxdz, dydz);
  }
  void setSiBeamAxisCm(int half, double x0, double y0, double dxdz, double dydz)
  {
    if (half == SiHalfA || half == SiHalfB)
    {
      m_axis[half] = {x0, y0, dxdz, dydz};
    }
  }
  void setSiBeamAxisMm(int half, double x0, double y0, double dxdz, double dydz)
  {
    setSiBeamAxisCm(half, 0.1 * x0, 0.1 * y0, dxdz, dydz);
  }
  // generic form (part = Tpc, SiHalfA, SiHalfB)
  void setBeamLineCm(int part, double x0, double y0, double dxdz, double dydz)
  {
    if (part >= 0 && part < NParts)
    {
      m_axis[part] = {x0, y0, dxdz, dydz};
    }
  }
  void setBeamLineMm(int part, double x0, double y0, double dxdz, double dydz)
  {
    setBeamLineCm(part, 0.1 * x0, 0.1 * y0, dxdz, dydz);
  }

  const BeamLine& tpcBeamAxis() const { return m_axis[Tpc]; }
  const BeamLine& siBeamAxis(int half) const { return m_axis[half]; }
  const BeamLine& beamLine(int part) const { return m_axis[part]; }

  // Beam position in the common (TPC) frame = TPC beam axis.  Used as beam / perigee
  // reference by the fits.
  double referenceBeamX(double z = 0.0) const { return m_axis[Tpc].x(z); }
  double referenceBeamY(double z = 0.0) const { return m_axis[Tpc].y(z); }

  // ---- other configuration
  // false: no Si correction (Si stays in the GLOBAL frame after the detector shift)
  void setEnabled(bool v) { m_enabled = v; }
  bool enabled() const { return m_enabled; }
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
  // half of a point given in the GLOBAL or TPC frame (Uphi from its azimuth around the
  // detector centre; the correction moves points by ~1 mm, far from the boundary for hits)
  int siHalfOfGlobal(const SiDetectorFrame& frame, const Point& g) const
  {
    const Point d = frame.toDetectorFromGlobal(g);
    return siHalf(std::atan2(d.y, d.x) - frame.uphiRotation());
  }

  // =========================================================================================
  // 2. transformation: Si half beam axis -> TPC beam axis
  // =========================================================================================
  // Offset of the Si half beam axis from the TPC beam axis at z: dx = x_half - x_tpc, dy = ...
  double siOffsetX(int half, double z) const { return m_axis[half].x(z) - m_axis[Tpc].x(z); }
  double siOffsetY(int half, double z) const { return m_axis[half].y(z) - m_axis[Tpc].y(z); }

  // Si point of one half: GLOBAL frame -> TPC frame
  Point siHalfToTpcAxis(int half, const Point& g) const
  {
    if (!m_enabled || (half != SiHalfA && half != SiHalfB))
    {
      return g;
    }
    return {g.x - siOffsetX(half, g.z), g.y - siOffsetY(half, g.z), g.z};
  }
  // inverse: TPC frame -> GLOBAL frame
  Point siHalfFromTpcAxis(int half, const Point& t) const
  {
    if (!m_enabled || (half != SiHalfA && half != SiHalfB))
    {
      return t;
    }
    return {t.x + siOffsetX(half, t.z), t.y + siOffsetY(half, t.z), t.z};
  }

  // TPC point: the TPC frame is the reference, identity (kept so the chain is explicit)
  Point tpcToTpcFrame(const Point& p) const { return p; }
  Point tpcToTpcFrame(double x, double y, double z) const { return {x, y, z}; }

  // ---- full Si chains
  // Si hit: intrinsic (Uphi, Vz) -> a. detector -> b. global (shift) -> c. TPC frame
  Point siToTpcFrame(const SiDetectorFrame& frame, int layer, double uphi, double vz) const
  {
    const Point det = frame.toDetector(layer, uphi, vz);  // a
    const Point glob = frame.toGlobal(det);               // b
    return siHalfToTpcAxis(siHalf(uphi), glob);           // c
  }
  // Si point already in the GLOBAL frame, half from its own Uphi
  Point siGlobalToTpcFrame(const SiDetectorFrame& frame, const Point& g) const
  {
    return siHalfToTpcAxis(siHalfOfGlobal(frame, g), g);
  }
  // inverse for display: TPC frame -> GLOBAL frame (half from the point's Uphi)
  Point siTpcFrameToGlobal(const SiDetectorFrame& frame, const Point& t) const
  {
    return siHalfFromTpcAxis(siHalfOfGlobal(frame, t), t);
  }

 private:
  bool m_enabled = true;
  double m_boundary = 0.0;
  std::array<BeamLine, NParts> m_axis{};
};

#endif
