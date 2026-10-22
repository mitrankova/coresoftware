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
//        x_glob = x_det + center(layer)
//      center = position of the detector-frame origin in the global frame
//      (default +6 mm, -1 mm, 0 mm, i.e. the same convention as
//      SiHitSeedQA::setDetectorCenterMm).
//      MVTX (layers 0-2) and INTT (layers 3-6) can have different centres: in the sPHENIX
//      geometry only the MVTX is displaced (MakeActsGeometry::makeMvtxMapPairs, MVTX_ALIGNMENT
//      global displacement D = (6.2955, -1.06874, -6.66289) mm; PHActsSiliconSeeding: MVTX
//      layers centred at (0.6, -0.1) cm, INTT at (0, 0)), and the TRACKINGALIGNMENT constants
//      are small corrections on top of that.  setDetectorCenterMm sets the MVTX centre and,
//      until setInttCenterMm is called, the INTT centre too (previous behaviour).
//      Use with the official alignment file:  setDetectorCenterMm(D)  +  setInttCenterMm(0, 0, 0).
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
  // INTT (layers 3-6) centre; not set: the INTT uses the MVTX centre
  void setInttCenterCm(double x, double y, double z)
  {
    m_inttCenter = {x, y, z};
    m_inttCenterSet = true;
  }
  void setInttCenterMm(double x, double y, double z) { setInttCenterCm(0.1 * x, 0.1 * y, 0.1 * z); }

  double uphiRotation() const { return m_rotation; }
  const std::array<double, 3>& detectorCenterCm() const { return m_center; }  // MVTX
  const std::array<double, 3>& inttCenterCm() const { return m_inttCenterSet ? m_inttCenter : m_center; }
  bool hasInttCenter() const { return m_inttCenterSet; }
  // centre used for a layer (< 0: MVTX centre)
  const std::array<double, 3>& centerCm(int layer) const { return layer >= 3 ? inttCenterCm() : m_center; }

  // ---- step 1: intrinsic (Uphi, Vz) on a layer -> detector frame
  Point toDetector(int layer, double uphi, double vz) const
  {
    const double r = layerRadius(layer);
    const double p = uphi + m_rotation;
    return {r * std::cos(p), r * std::sin(p), vz};
  }

  // ---- step 2: detector frame -> global frame (layer: MVTX / INTT centre; without: MVTX)
  Point toGlobal(int layer, const Point& det) const
  {
    const auto& c = centerCm(layer);
    return {det.x + c[0], det.y + c[1], det.z + c[2]};
  }
  Point toGlobal(const Point& det) const { return toGlobal(-1, det); }

  Point toDetectorFromGlobal(int layer, const Point& glob) const
  {
    const auto& c = centerCm(layer);
    return {glob.x - c[0], glob.y - c[1], glob.z - c[2]};
  }
  Point toDetectorFromGlobal(const Point& glob) const { return toDetectorFromGlobal(-1, glob); }

 private:
  double m_rotation = 0.1481;                    // [rad]
  std::array<double, 3> m_center{{0.6, -0.1, 0.0}};  // [cm] = (+6, -1, 0) mm, MVTX
  std::array<double, 3> m_inttCenter{{0.0, 0.0, 0.0}};  // [cm], used only if set
  bool m_inttCenterSet = false;
};

#endif
