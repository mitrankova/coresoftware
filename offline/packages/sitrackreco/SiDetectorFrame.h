// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SIDETECTORFRAME_H
#define SITRACKRECO_SIDETECTORFRAME_H

#include <array>
#include <cmath>

// Two-step silicon coordinate chain, kept deliberately separate:
//
//   1. hit (Uphi, Vz) -> DETECTOR frame      (idealized, centred on the silicon barrel)
//        x_det = R_layer cos(Uphi + rotation)
//        y_det = R_layer sin(Uphi + rotation)
//        z_det = Vz
//      Uphi, Vz are the nominal intrinsic coordinates filled by SiHitSeedReco
//      (SiHitPoint::phi, SiHitPoint::z); rotation = phi0 of MVTX layer-2 stave 0.
//
//   2. DETECTOR frame -> GLOBAL frame        (global alignment of the whole detector)
//        x_glob = x_det + center
//      center = position of the detector-frame origin in the global frame
//      (default +6 mm, -1 mm, 0 mm, i.e. the same convention as
//      SiHitSeedQA::setDetectorCenterMm).
//
// Header-only, no ROOT dictionary.  All lengths in cm.
class SiDetectorFrame
{
 public:
  struct Point
  {
    double x = 0;
    double y = 0;
    double z = 0;
  };

  static constexpr int kNLayers = 7;

  // nominal layer radii [cm], same as SiHitSeedReco
  static constexpr std::array<double, kNLayers> kRadius{
      {2.523, 3.336, 4.148, 7.188 - 0.0036, 7.732 - 0.0036, 9.680 - 0.0036, 10.262 - 0.0036}};

  static double layerRadius(int layer)
  {
    return (layer >= 0 && layer < kNLayers) ? kRadius[layer] : NAN;
  }

  // ---- configuration
  void setUphiRotation(double rad) { m_rotation = rad; }
  void setDetectorCenterCm(double x, double y, double z) { m_center = {x, y, z}; }
  void setDetectorCenterMm(double x, double y, double z) { m_center = {0.1 * x, 0.1 * y, 0.1 * z}; }

  double uphiRotation() const { return m_rotation; }
  const std::array<double, 3>& detectorCenterCm() const { return m_center; }

  // ---- step 1: intrinsic (Uphi, Vz) on a layer -> detector frame
  Point toDetector(int layer, double uphi, double vz) const
  {
    const double r = layerRadius(layer);
    const double p = uphi + m_rotation;
    return {r * std::cos(p), r * std::sin(p), vz};
  }

  // ---- step 2: detector frame -> global frame
  Point toGlobal(const Point& det) const
  {
    return {det.x + m_center[0], det.y + m_center[1], det.z + m_center[2]};
  }

  Point toDetectorFromGlobal(const Point& glob) const
  {
    return {glob.x - m_center[0], glob.y - m_center[1], glob.z - m_center[2]};
  }

 private:
  double m_rotation = 0.1481;                    // [rad]
  std::array<double, 3> m_center{{0.6, -0.1, 0.0}};  // [cm] = (+6, -1, 0) mm
};

#endif
