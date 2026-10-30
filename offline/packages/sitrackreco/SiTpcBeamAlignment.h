// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCBEAMALIGNMENT_H
#define SITRACKRECO_SITPCBEAMALIGNMENT_H

#include "SiDetectorFrame.h"

#include <array>
#include <cmath>
#include <iostream>

// Beam-centred common transverse frame for the TPC and the silicon.
//
// ---- 1. The beam lines (measured, one per part)
//     x_axis(z) = x0 + dxdz * z ,   y_axis(z) = y0 + dydz * z
//   Tpc     : TPC polytracks, pre-ACTS, TPC frame (no TPC survey transform)        [cm]
//   SiHalfA : silicon half A, Uphi in [B, B+pi)                                     [mm/cm]
//   SiHalfB : silicon half B, Uphi in [B+pi, B+2pi)                                 [mm/cm]
//   SiFull  : whole silicon (both halves together), used when per-half mode is off  [mm/cm]
//   B = clamshell boundary in intrinsic Uphi (default 0).  With the default rotation
//   (global phi = Uphi + 0.1481) half A is the upper half (global phi in [8.5, 188.5) deg).
//
// ---- 2. Which frame the Si beam lines must be measured in
// The Si chain in SiHitSeedReco is
//   a. intrinsic (Uphi, Vz)                       hit
//   b. per-sensor alignment (optional)            SiAlignmentTransform, setApplySiAlignment(true)
//   c. DETECTOR -> GLOBAL                         SiDetectorFrame::toGlobal (+ detector centre)
//   d. GLOBAL -> beam-centred frame               this class
// The Si beam lines must be fitted from vertices reconstructed with EXACTLY the same a-c:
// same alignment file (on/off), same rotation order, same detector centre, same boundary.
// The detector centre then cancels in the output (the beam is at (0,0) whatever centre is
// used), but a beam line measured with a different centre or without the alignment file
// gives a wrong offset.  Use shiftBeamLineMm() only as a PROVISIONAL bridge between geometries.
//
// The TRACKINGALIGNMENT silicon constants are relative to the ACTS geometry with the MVTX
// survey displacement applied (MakeActsGeometry::makeMvtxMapPairs, "Global displacement"),
// i.e. the MVTX wrapper at D = (6.2955, -1.06874, -6.66289) mm.  For the per-sensor file to be
// applied consistently, the SiDetectorFrame centre x, y must be D_x, D_y (not 6, -1).  Keep the
// centre z at 0 until SiHitSeedReco seeds with the same z as its residuals (z-alignment step).
//
// ABI note: this class is passed BY VALUE between libsitrackreco and libTrackingDiagnostics
// (SiHitSeedReco, SiTrajectoryFitter, SiTpcTrackMatcher, SiTpc_EventDisplay).  After changing
// this header rebuild ALL libraries that include it.
//
// ---- 3. The transformation (BeamLine mode: x, y only, tilts ignored; World mode: see 5.)
// centerBeam = true  (default):  every part is moved so that its own beam is at (0, 0)
//     TPC : x' = x - x0_tpc            y' = y - y0_tpc
//     Si  : x' = x - x0_si             y' = y - y0_si       (x0_si = half A/B or SiFull)
// centerBeam = false            :  legacy mode, TPC unchanged, Si moved onto the TPC beam
//     Si  : x' = x - x0_si + x0_tpc    y' = y - y0_si + y0_tpc
// The fitted slopes are kept (for the later tilt step) but are NOT applied.
//
// ---- 4. Per-half vs common silicon offset
// perHalf = true  (default): each half gets its own offset (absorbs a residual half-to-half
//                             misalignment, makes the two halves agree BY CONSTRUCTION).
// perHalf = false           : both halves get the SiFull offset; use this to test whether the
//                             per-sensor alignment file alone makes the halves consistent.
//
// ---- 5. TPC frame: beam line (legacy) or world (setTpcFrame)
// TpcFrame::BeamLine (default): TPC output = raw poly position (TPC envelope coordinates) - TPC beam line.
//     The seeder TPC frame then differs from the ACTS TPC frame by the envelope -> world transform E
//     (TPC survey placement / tilt from Tracking_Geometry), and the TPC is aligned to the Si only at z = 0.
// TpcFrame::World: TPC output = E'(raw) - Si beam (SiFull x0, y0), z + tpcZShift, where E' is the
//     envelope -> world transform of ActsGeometry, i.e. EXACTLY what the ACTS fit sees (including the
//     TpcGlobalAlignment correction).  The TPC and the Si are then in one frame, the same one in the
//     seeder and in ACTS; the TPC beam line is not used.  E' is copied from ActsGeometry by
//     SiHitSeedReco::InitRun (setTpcToWorld); until then the BeamLine mode is used.
//     tpcZShift = the Si z shift of the seeder (SiHitSeedReco::setActsZShiftCm) so that TPC and Si z agree.
//
// All lengths in cm internally.  Header only.
class SiTpcBeamAlignment
{
 public:
  enum Part
  {
    Tpc = 0,
    SiHalfA = 1,
    SiHalfB = 2,
    SiFull = 3,
    NParts = 4
  };

  enum class TpcFrame
  {
    BeamLine = 0,
    World = 1
  };

  struct BeamLine
  {
    double x0 = 0;    // [cm] at z = 0
    double y0 = 0;    // [cm] at z = 0
    double dxdz = 0;  // [rad], metadata only
    double dydz = 0;  // [rad], metadata only
    double x(double z) const { return x0 + dxdz * z; }
    double y(double z) const { return y0 + dydz * z; }
  };

  using Point = SiDetectorFrame::Point;

  SiTpcBeamAlignment()
  {
    // TPC polytracks, pre-ACTS, run 79513                       [cm]
    setBeamLineCm(Tpc, -0.106835, 0.267482, 8.478616e-04, -8.041254e-04);
    // UNALIGNED silicon fits (no per-sensor alignment), global frame with centre (+6,-1,0) mm.
    // Re-measure with setApplySiAlignment(true) and centre D before using them for real.  [mm]
    setBeamLineMm(SiHalfA, -0.533300, 1.772693, 7.094553e-04, -2.743455e-04);
    setBeamLineMm(SiHalfB, 0.595174, 1.608226, 2.343532e-03, 6.926381e-04);
    setBeamLineMm(SiFull, 0.5 * (-0.533300 + 0.595174), 0.5 * (1.772693 + 1.608226), 0., 0.);
  }

  // =========================================================================================
  // 1. beam lines
  // =========================================================================================
  void setBeamLineCm(int part, double x0, double y0, double dxdz, double dydz)
  {
    if (validPart(part))
    {
      m_axis[part] = {x0, y0, dxdz, dydz};
    }
  }
  void setBeamLineMm(int part, double x0, double y0, double dxdz, double dydz)
  {
    setBeamLineCm(part, 0.1 * x0, 0.1 * y0, dxdz, dydz);
  }
  // Shift an already set beam line (e.g. old unaligned fit + alignment-file translation).
  void shiftBeamLineMm(int part, double dx, double dy)
  {
    if (validPart(part))
    {
      m_axis[part].x0 += 0.1 * dx;
      m_axis[part].y0 += 0.1 * dy;
    }
  }

  // backward-compatible setters
  void setTpcBeamAxisCm(double x0, double y0, double dxdz, double dydz) { setBeamLineCm(Tpc, x0, y0, dxdz, dydz); }
  void setTpcBeamAxisMm(double x0, double y0, double dxdz, double dydz) { setBeamLineMm(Tpc, x0, y0, dxdz, dydz); }
  void setSiBeamAxisCm(int half, double x0, double y0, double dxdz, double dydz)
  {
    if (validSi(half))
    {
      setBeamLineCm(half, x0, y0, dxdz, dydz);
    }
  }
  void setSiBeamAxisMm(int half, double x0, double y0, double dxdz, double dydz)
  {
    setSiBeamAxisCm(half, 0.1 * x0, 0.1 * y0, dxdz, dydz);
  }

  const BeamLine& beamLine(int part) const { return m_axis[validPart(part) ? part : Tpc]; }
  const BeamLine& tpcBeamAxis() const { return m_axis[Tpc]; }
  const BeamLine& siBeamAxis(int half) const { return m_axis[siPart(half)]; }

  // =========================================================================================
  // configuration
  // =========================================================================================
  void setEnabled(bool on) { m_enabled = on; }
  bool enabled() const { return m_enabled; }

  void setCenterBeam(bool on) { m_centerBeam = on; }
  bool centerBeam() const { return m_centerBeam; }

  void setTpcFrame(TpcFrame f) { m_tpcFrame = f; }
  TpcFrame tpcFrame() const { return m_tpcFrame; }

  // envelope -> world of the TPC: p' = R p + t (R row-major), plus z shift [cm]
  void setTpcToWorld(const std::array<double, 9>& r, const std::array<double, 3>& t, double zShift)
  {
    m_tpcR = r;
    m_tpcT = t;
    m_tpcZShift = zShift;
    m_tpcToWorldSet = true;
  }
  bool tpcToWorldSet() const { return m_tpcToWorldSet; }
  bool tpcWorldActive() const { return m_tpcFrame == TpcFrame::World && m_tpcToWorldSet; }

  void setSiPerHalf(bool on) { m_perHalf = on; }
  bool siPerHalf() const { return m_perHalf; }

  void setClamshellBoundary(double b) { m_boundary = b; }
  double clamshellBoundary() const { return m_boundary; }

  // Beam position in the OUTPUT frame (used as the helix-fit beam constraint).
  // Disabled or legacy mode: the TPC frame is the output, so the beam sits at the TPC beam.
  // World mode is always beam-centred.
  double referenceBeamX(double /*z*/ = 0.) const { return (m_enabled && centred()) ? 0. : m_axis[Tpc].x0; }
  double referenceBeamY(double /*z*/ = 0.) const { return (m_enabled && centred()) ? 0. : m_axis[Tpc].y0; }

  // =========================================================================================
  // halves
  // =========================================================================================
  int siHalf(double uphi) const
  {
    constexpr double pi = 3.14159265358979323846;
    constexpr double twopi = 2 * pi;
    double u = std::fmod(uphi - m_boundary, twopi);
    if (u < 0)
    {
      u += twopi;
    }
    return u < pi ? SiHalfA : SiHalfB;
  }
  int siHalfOfGlobal(const SiDetectorFrame& frame, const Point& g) const
  {
    const auto d = frame.toDetectorFromGlobal(g);
    return siHalf(std::atan2(d.y, d.x) - frame.uphiRotation());
  }

  // Si offset actually used for a half: its own line, or SiFull in common mode.
  double siX0(int half) const { return m_axis[siPart(half)].x0; }
  double siY0(int half) const { return m_axis[siPart(half)].y0; }

  // Offset of the Si beam from the TPC beam, in the input coordinates [cm].
  double siOffsetX(int half, double /*z*/ = 0.) const { return siX0(half) - m_axis[Tpc].x0; }
  double siOffsetY(int half, double /*z*/ = 0.) const { return siY0(half) - m_axis[Tpc].y0; }

  // =========================================================================================
  // transformations
  // =========================================================================================
  // Si point (GLOBAL frame, after per-sensor alignment and detector centre) -> output frame
  Point siHalfToTpcAxis(int half, const Point& g) const
  {
    if (!m_enabled || !validHalf(half))
    {
      return g;
    }
    return {g.x - siX0(half) + refX(), g.y - siY0(half) + refY(), g.z};
  }
  Point siHalfFromTpcAxis(int half, const Point& p) const
  {
    if (!m_enabled || !validHalf(half))
    {
      return p;
    }
    return {p.x + siX0(half) - refX(), p.y + siY0(half) - refY(), p.z};
  }

  // TPC point (raw poly position, TPC envelope coordinates) -> output frame
  Point tpcToTpcFrame(const Point& p) const
  {
    if (!m_enabled)
    {
      return p;
    }
    if (tpcWorldActive())
    {
      const Point w = tpcEnvelopeToWorld(p);
      return {w.x - m_axis[SiFull].x0, w.y - m_axis[SiFull].y0, w.z + m_tpcZShift};
    }
    if (!centred())
    {
      return p;
    }
    return {p.x - m_axis[Tpc].x0, p.y - m_axis[Tpc].y0, p.z};
  }
  Point tpcToTpcFrame(double x, double y, double z) const { return tpcToTpcFrame(Point{x, y, z}); }
  Point tpcFromTpcFrame(const Point& p) const
  {
    if (!m_enabled)
    {
      return p;
    }
    if (tpcWorldActive())
    {
      return tpcWorldToEnvelope({p.x + m_axis[SiFull].x0, p.y + m_axis[SiFull].y0, p.z - m_tpcZShift});
    }
    if (!centred())
    {
      return p;
    }
    return {p.x + m_axis[Tpc].x0, p.y + m_axis[Tpc].y0, p.z};
  }

  // the copied ActsGeometry transform (identity until setTpcToWorld)
  Point tpcEnvelopeToWorld(const Point& p) const
  {
    const auto& r = m_tpcR;
    return {r[0] * p.x + r[1] * p.y + r[2] * p.z + m_tpcT[0],
            r[3] * p.x + r[4] * p.y + r[5] * p.z + m_tpcT[1],
            r[6] * p.x + r[7] * p.y + r[8] * p.z + m_tpcT[2]};
  }
  Point tpcWorldToEnvelope(const Point& w) const
  {
    // R is a rotation: inverse = transpose
    const auto& r = m_tpcR;
    const double x = w.x - m_tpcT[0], y = w.y - m_tpcT[1], z = w.z - m_tpcT[2];
    return {r[0] * x + r[3] * y + r[6] * z, r[1] * x + r[4] * y + r[7] * z, r[2] * x + r[5] * y + r[8] * z};
  }

  // Full Si chain without per-sensor alignment (intrinsic -> detector -> global -> output).
  Point siToTpcFrame(const SiDetectorFrame& frame, int layer, double uphi, double vz) const
  {
    const Point det = frame.toDetector(layer, uphi, vz);
    return siHalfToTpcAxis(siHalf(uphi), frame.toGlobal(layer, det));
  }
  Point siGlobalToTpcFrame(const SiDetectorFrame& frame, const Point& g) const
  {
    return siHalfToTpcAxis(siHalfOfGlobal(frame, g), g);
  }

  // Inverse for display.  The half of an output point is ambiguous near the boundary when
  // the two halves have different offsets; iterate until the half is self-consistent.
  Point siTpcFrameToGlobal(const SiDetectorFrame& frame, const Point& p) const
  {
    int h = siHalfOfGlobal(frame, p);
    for (int i = 0; i < 2; ++i)
    {
      const Point g = siHalfFromTpcAxis(h, p);
      const int next = siHalfOfGlobal(frame, g);
      if (next == h)
      {
        return g;
      }
      h = next;
    }
    return siHalfFromTpcAxis(h, p);
  }

  // =========================================================================================
  // printout (call from the macro to check what was configured)
  // =========================================================================================
  void print(std::ostream& os = std::cout) const
  {
    static const char* names[NParts] = {"TPC", "Si half A", "Si half B", "Si full"};
    os << "SiTpcBeamAlignment: " << (m_enabled ? "ON" : "OFF")
       << ", " << (centred() ? "beam-centred frame" : "TPC frame (legacy)")
       << ", TPC " << (m_tpcFrame == TpcFrame::World ? (m_tpcToWorldSet ? "WORLD (ActsGeometry envelope -> world, minus Si beam)"
                                                                          : "WORLD requested, transform not loaded yet (beam-line mode until SiHitSeedReco::InitRun)")
                                                      : "beam line")
       << ", Si offsets " << (m_perHalf ? "per half" : "common (SiFull)")
       << ", clamshell boundary " << m_boundary << " rad" << std::endl;
    for (int p = 0; p < NParts; ++p)
    {
      os << "  " << names[p] << ": x0 = " << 10. * m_axis[p].x0 << " mm, y0 = " << 10. * m_axis[p].y0
         << " mm (dx/dz = " << 1e3 * m_axis[p].dxdz << ", dy/dz = " << 1e3 * m_axis[p].dydz << " mrad, not applied)"
         << std::endl;
    }
    for (int h : {SiHalfA, SiHalfB})
    {
      os << "  applied shift " << names[h] << ": (" << -10. * (siX0(h) - refX()) << ", " << -10. * (siY0(h) - refY())
         << ") mm" << std::endl;
    }
    if (tpcWorldActive())
    {
      const Point o = tpcToTpcFrame(Point{0., 0., 0.});
      const Point n = tpcToTpcFrame(Point{0., 0., 100.});
      os << "  TPC: envelope (0,0,0) -> (" << 10. * o.x << ", " << 10. * o.y << ", " << 10. * o.z << ") mm, (0,0,1000 mm) -> ("
         << 10. * n.x << ", " << 10. * n.y << ", " << 10. * n.z << ") mm, z shift " << 10. * m_tpcZShift << " mm" << std::endl;
    }
    else
    {
      os << "  applied shift TPC: (" << (centred() ? -10. * m_axis[Tpc].x0 : 0.) << ", "
         << (centred() ? -10. * m_axis[Tpc].y0 : 0.) << ") mm" << std::endl;
    }
  }

 private:
  static bool validPart(int p) { return p >= 0 && p < NParts; }
  static bool validHalf(int h) { return h == SiHalfA || h == SiHalfB; }
  static bool validSi(int h) { return validHalf(h) || h == SiFull; }
  int siPart(int half) const { return (m_perHalf && validHalf(half)) ? half : SiFull; }
  bool centred() const { return m_centerBeam || m_tpcFrame == TpcFrame::World; }
  double refX() const { return centred() ? 0. : m_axis[Tpc].x0; }
  double refY() const { return centred() ? 0. : m_axis[Tpc].y0; }

  bool m_enabled = true;
  bool m_centerBeam = true;
  bool m_perHalf = true;
  double m_boundary = 0.0;
  std::array<BeamLine, NParts> m_axis{};

  TpcFrame m_tpcFrame = TpcFrame::BeamLine;
  bool m_tpcToWorldSet = false;
  std::array<double, 9> m_tpcR{1., 0., 0., 0., 1., 0., 0., 0., 1.};
  std::array<double, 3> m_tpcT{0., 0., 0.};
  double m_tpcZShift = 0.;
};

#endif