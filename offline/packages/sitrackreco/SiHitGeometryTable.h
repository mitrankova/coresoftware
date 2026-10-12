// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SIHITGEOMETRYTABLE_H
#define SITRACKRECO_SIHITGEOMETRYTABLE_H

#include "SiDetectorFrame.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// SiHitGeometryTable
//
// Per-run lookup tables that turn a decoded silicon hit (SiHitPoint from SiHitSeedReco: layer,
// row, col, aligned intrinsic Uphi and Vz) into the physics-space position used by
// SiTpcHitChainMatcher, plus the detector constants that define which hits are neighbours
// (pitches, sensitive thickness, chip gap, layer half lengths).
//
// Frame: the silicon DETECTOR frame of SiDetectorFrame (centred on the barrel), intrinsic
// azimuth u = Uphi (global phi = u + rotation).  All lengths in cm.
//
// Position model
//   Linear (default)  exactly the SiHitSeedReco coordinates:
//                       u = SiHitPoint::phi, z = SiHitPoint::z, r = nominal layer radius.
//                     The TPC-driven and the independent chaining then see identical hit
//                     positions, so the hit-level cross-check is clean.
//   Planar            the sensor is a flat plane at the nominal radius, tilted by the stave tilt
//                     (MVTX) or not tilted (INTT).  SiHitSeedReco uses the first-order azimuth
//                     lx*cos(tilt)/R; the planar model adds the remainder per row:
//                       du(row) = atan2(lx cos(tilt), R + s lx sin(tilt)) - lx cos(tilt)/R
//                       r(row)  = sqrt(R^2 + 2 s R lx sin(tilt) + lx^2)
//                     s = setMvtxTiltSign (+1/-1), the orientation of the stave tilt.
//                     VALIDATE the sign against the ACTS surface positions before using the
//                     planar model for tuning.  The alignment applied by SiHitSeedReco is kept;
//                     the correction is added on top (exact to first order, alignment angles
//                     are mrad).
//
// No trig per hit: everything is a table lookup by (layer, row), filled once in build().
class SiHitGeometryTable
{
 public:
  enum Model
  {
    Linear = 0,
    Planar = 1
  };

  static constexpr int kNLayers = 7;
  static constexpr int kMvtxRows = 512;
  static constexpr int kMvtxCols = 1024;
  static constexpr int kInttRows = 256;

  struct Position
  {
    float u = 0;  // intrinsic azimuth [0, 2pi)
    float z = 0;
    float r = 0;
  };

  SiHitGeometryTable() { build(); }

  // ---- configuration (call build() afterwards; the matcher does it in InitRun)
  void setModel(Model m) { m_model = m; }
  void setMvtxTiltSign(int s) { m_tiltSign = s >= 0 ? 1 : -1; }
  void setMvtxThicknessCm(double t) { m_mvtxThickness = t; }
  void setInttThicknessCm(double t) { m_inttThickness = t; }
  Model model() const { return m_model; }

  void build()
  {
    for (int layer = 0; layer < kNLayers; ++layer)
    {
      const bool mvtx = isMvtx(layer);
      const int nrows = mvtx ? kMvtxRows : kInttRows;
      const double R = radius(layer);
      const double tilt = mvtx ? kMvtxTilt[layer] : 0.0;
      auto& du = m_du[layer];
      auto& rr = m_r[layer];
      du.assign(nrows, 0.f);
      rr.assign(nrows, static_cast<float>(R));
      m_maxDr[layer] = 0;
      if (m_model != Planar)
      {
        continue;
      }
      const double c = std::cos(tilt), s = m_tiltSign * std::sin(tilt);
      for (int row = 0; row < nrows; ++row)
      {
        const double lx = localX(layer, row);
        const double exact = std::atan2(lx * c, R + lx * s);
        du[row] = static_cast<float>(exact - lx * c / R);
        rr[row] = static_cast<float>(std::sqrt(R * R + 2 * R * lx * s + lx * lx));
        m_maxDr[layer] = std::max(m_maxDr[layer], std::abs(rr[row] - R));
      }
    }
  }

  // ---- hit position (uphi, vz: aligned values from SiHitPoint)
  Position position(int layer, unsigned int row, double uphi, double vz) const
  {
    Position p;
    const auto& du = m_du[layer];
    const auto& rr = m_r[layer];
    const bool ok = row < du.size();
    double u = uphi + (ok ? du[row] : 0.0);
    if (u >= kTwoPi)
    {
      u -= kTwoPi;
    }
    else if (u < 0)
    {
      u += kTwoPi;
    }
    p.u = static_cast<float>(u);
    p.z = static_cast<float>(vz);
    p.r = ok ? rr[row] : static_cast<float>(radius(layer));
    return p;
  }

  // ---- detector constants
  static bool isMvtx(int layer) { return layer >= 0 && layer <= 2; }
  static bool isIntt(int layer) { return layer >= 3 && layer <= 6; }
  static double radius(int layer) { return SiDetectorFrame::layerRadius(layer); }

  // pitch along r*phi (MVTX rows, INTT strips) [cm]
  static double pitchRphi(int layer) { return isMvtx(layer) ? kMvtxRowPitch : kInttStripPitch; }
  // pitch along z: MVTX column; INTT strip length (ladderz % 2 == 0: 1.6 cm, else 2.0 cm)
  static double pitchZ(int layer, int ladderz = 0)
  {
    if (isMvtx(layer))
    {
      return kMvtxColPitch;
    }
    return (ladderz % 2 == 0) ? 1.6 : 2.0;
  }
  // number of MVTX staves / INTT ladders in phi
  static int nPhiSegments(int layer)
  {
    constexpr std::array<int, kNLayers> n{{12, 16, 20, 12, 12, 16, 16}};
    return (layer >= 0 && layer < kNLayers) ? n[layer] : 0;
  }
  double thickness(int layer) const { return isMvtx(layer) ? m_mvtxThickness : m_inttThickness; }
  static double staveTilt(int layer) { return isMvtx(layer) ? kMvtxTilt[layer] : 0.0; }
  // gap between abutting MVTX chips along z (chip pitch - active length) [cm]
  static double mvtxChipGap() { return kMvtxChipPitchZ - kMvtxCols * kMvtxColPitch; }
  // half length of the active area along z [cm]
  static double halfLengthZ(int layer) { return isMvtx(layer) ? 13.57 : 23.2; }
  // largest |r_hit - R_layer| over the rows (0 in the Linear model)
  double maxDr(int layer) const { return (layer >= 0 && layer < kNLayers) ? m_maxDr[layer] : 0.0; }

  // local x of a row, identical to SiHitSeedReco::localCoordinates
  static double localX(int layer, int row)
  {
    if (isMvtx(layer))
    {
      constexpr double sensor_x = 0.058128;
      constexpr double first_row = 0.5 * (kMvtxRows * kMvtxRowPitch - 37.44e-4 + 0.12 - kMvtxRowPitch);
      return first_row - row * kMvtxRowPitch - sensor_x;
    }
    return (row + 0.5 - 128.) * kInttStripPitch;
  }

  static constexpr double kTwoPi = 2.0 * M_PI;

 private:
  // same constants as SiHitSeedReco
  static constexpr double kMvtxRowPitch = 26.88e-4;
  static constexpr double kMvtxColPitch = 29.24e-4;
  static constexpr double kMvtxChipPitchZ = 3.015;
  static constexpr double kInttStripPitch = 0.0078;
  static constexpr std::array<double, 3> kMvtxTilt{{13.3 * M_PI / 180., 16.9 * M_PI / 180., 17.0 * M_PI / 180.}};

  Model m_model = Linear;
  int m_tiltSign = 1;
  double m_mvtxThickness = 25e-4;   // sensitive epi layer [cm]
  double m_inttThickness = 320e-4;  // [cm]
  std::array<std::vector<float>, kNLayers> m_du;
  std::array<std::vector<float>, kNLayers> m_r;
  std::array<double, kNLayers> m_maxDr{};
};

#endif
