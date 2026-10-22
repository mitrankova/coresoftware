// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCREFERENCE_H
#define SITRACKRECO_SITPCREFERENCE_H

#include "SiDetectorFrame.h"
#include "SiTpcBeamAlignment.h"
#include "SiTpcHelixFit.h"

#include <algorithm>
#include <array>
#include <cmath>

// TPC track -> silicon HARDWARE (detector-frame) coordinates.  Used by SiHitSeedReco to put the
// anchor window of the TPC-driven silicon chaining on the outer seed layer.  Free of Fun4All.
//
// Frames
//   TPC helix : TPC frame (SiTpcHelixFit on the TPC clusters, beam = TPC beam axis).
//   Si hits   : silicon DETECTOR frame (SiDetectorFrame): layers are centred cylinders, the hit
//               coordinates (Uphi, Vz) of SiHitSeedReco live there.
//   A TPC-frame point goes to the detector frame of clamshell half `half` as
//     det = SiDetectorFrame::toDetectorFromGlobal(SiTpcBeamAlignment::siHalfFromTpcAxis(half, tpc))
//   i.e. a pure translation that depends on the clamshell half and (weakly) on z.
//
// Per layer the helix is intersected with the nominal layer cylinder once (closed-form
// circle-circle intersection, iterated for the half / z dependence of the translation).
// Near a clamshell boundary the prediction with the other half's translation is kept as well
// (the halves are displaced by ~1 mm).
namespace SiTpcReference
{
  constexpr double kPi = M_PI;
  constexpr double kTwoPi = 2.0 * M_PI;

  inline double wrapPi(double a)
  {
    a = std::fmod(a + kPi, kTwoPi);
    if (a < 0)
    {
      a += kTwoPi;
    }
    return a - kPi;
  }
  inline double wrap2Pi(double a)
  {
    a = std::fmod(a, kTwoPi);
    return a < 0 ? a + kTwoPi : a;
  }

  // ---- TPC frame <-> silicon detector frame
  struct FrameMap
  {
    const SiDetectorFrame* frame = nullptr;
    const SiTpcBeamAlignment* alignment = nullptr;

    // layer: MVTX / INTT detector centre (SiDetectorFrame::centerCm); < 0: MVTX
    SiDetectorFrame::Point tpcToDet(double x, double y, double z, int half, int layer = -1) const
    {
      const auto g = alignment->siHalfFromTpcAxis(half, {x, y, z});
      return frame->toDetectorFromGlobal(layer, g);
    }
    // detector frame point of clamshell half `half` -> TPC frame
    SiDetectorFrame::Point detToTpc(double x, double y, double z, int half, int layer = -1) const
    {
      const auto g = frame->toGlobal(layer, {x, y, z});
      return alignment->siHalfToTpcAxis(half, g);
    }
    // det = tpc + (dx, dy) at fixed z, half and layer
    void offset(int half, double z, double& dx, double& dy, int layer = -1) const
    {
      const auto p = tpcToDet(0, 0, z, half, layer);
      dx = p.x;
      dy = p.y;
    }
    // intrinsic Uphi of a detector-frame point (same convention as SiHitPoint::phi)
    double uphi(double xd, double yd) const { return wrap2Pi(std::atan2(yd, xd) - frame->uphiRotation()); }
    int half(double u) const { return alignment->siHalf(u); }
    int otherHalf(int h) const { return h == SiTpcBeamAlignment::SiHalfA ? SiTpcBeamAlignment::SiHalfB : SiTpcBeamAlignment::SiHalfA; }
    // r*phi distance [cm] of azimuth u on radius rl to the nearest clamshell boundary (B or B + pi)
    double boundaryDistance(double u, double rl) const
    {
      const double b = alignment->clamshellBoundary();
      return rl * std::min(std::abs(wrapPi(u - b)), std::abs(wrapPi(u - b - kPi)));
    }
  };

  // ---- prediction of the TPC helix on one silicon layer, detector frame
  struct LayerPrediction
  {
    bool valid = false;
    int layer = -1;
    int half = 0;
    double R = 0;                // nominal layer radius
    double x = 0, y = 0, z = 0;  // detector frame
    double u = 0;                // intrinsic azimuth Uphi
    double s = 0;                // transverse path length from the TPC pca
    double cosA = 1;             // incidence: cos of the angle between direction and radial
  };

  struct TpcReference
  {
    bool valid = false;
    // TPC frame helix (SiTpcHelixFit::Result)
    double cx = 0, cy = 0, R = 0;
    int helicity = 0;
    double pcaX = 0, pcaY = 0;
    double z0 = 0, tanl = 0;
    double phi = 0, dca = 0, pt = 0;
    int charge = 0;
    std::array<LayerPrediction, SiDetectorFrame::kNLayers> layers{};
    // prediction with the other clamshell half's translation, valid only near a boundary
    std::array<LayerPrediction, SiDetectorFrame::kNLayers> alt{};
  };

  // Intersections of the circle (ccx, ccy, R) with the circle of radius rl around (0, 0).
  inline int circleIntersections(double ccx, double ccy, double R, double rl, double px[2], double py[2])
  {
    const double d = std::hypot(ccx, ccy);
    if (!(d > 0) || d > rl + R || d < std::abs(rl - R))
    {
      return 0;
    }
    const double a = (rl * rl - R * R + d * d) / (2 * d);
    const double h = std::sqrt(std::max(0.0, rl * rl - a * a));
    const double ux = ccx / d, uy = ccy / d;
    px[0] = a * ux - h * uy;
    py[0] = a * uy + h * ux;
    px[1] = a * ux + h * uy;
    py[1] = a * uy - h * ux;
    return 2;
  }

  // forceHalf >= 0: use that clamshell half's translation instead of the half of the point
  inline bool predictTpcLayer(const TpcReference& ref, const FrameMap& map, int layer, double rl, LayerPrediction& out,
                              int forceHalf = -1)
  {
    out = LayerPrediction();
    if (!ref.valid || !(rl > 0))
    {
      return false;
    }
    // first guess of the half: TPC-frame point at ~rl along the pca direction
    const double tx0 = std::cos(ref.phi), ty0 = std::sin(ref.phi);
    const auto p0 = map.tpcToDet(ref.pcaX + rl * tx0, ref.pcaY + rl * ty0, ref.z0 + ref.tanl * rl, SiTpcBeamAlignment::SiHalfA, layer);
    int half = forceHalf >= 0 ? forceHalf : map.half(map.uphi(p0.x, p0.y));
    double zg = ref.z0 + ref.tanl * rl;
    double px = 0, py = 0, s = 0, ccx = 0, ccy = 0;
    bool found = false;
    for (int iter = 0; iter < 3; ++iter)
    {
      double dx = 0, dy = 0;
      map.offset(half, zg, dx, dy, layer);
      ccx = ref.cx + dx;
      ccy = ref.cy + dy;
      const double qx = ref.pcaX + dx, qy = ref.pcaY + dy;
      double cxs[2], cys[2];
      if (circleIntersections(ccx, ccy, ref.R, rl, cxs, cys) == 0)
      {
        return false;
      }
      const double a0 = std::atan2(qy - ccy, qx - ccx);
      double best = 1e30;
      for (int k = 0; k < 2; ++k)
      {
        const double sk = ref.R * ref.helicity * wrapPi(std::atan2(cys[k] - ccy, cxs[k] - ccx) - a0);
        if (sk > 0 && sk < best)
        {
          best = sk;
          px = cxs[k];
          py = cys[k];
        }
      }
      if (best > 1e29)
      {
        return false;  // not on the outgoing branch
      }
      s = best;
      found = true;
      zg = ref.z0 + ref.tanl * s;
      const int newHalf = forceHalf >= 0 ? forceHalf : map.half(map.uphi(px, py));
      if (newHalf == half && iter > 0)
      {
        break;
      }
      half = newHalf;
    }
    if (!found)
    {
      return false;
    }
    out.layer = layer;
    out.half = half;
    out.R = rl;
    out.x = px;
    out.y = py;
    out.s = s;
    out.z = ref.z0 + ref.tanl * s;
    out.u = map.uphi(px, py);
    const double tx = ref.helicity * (-(py - ccy)) / ref.R;
    const double ty = ref.helicity * (px - ccx) / ref.R;
    out.cosA = tx * px / rl + ty * py / rl;
    if (out.cosA < 1e-3)
    {
      return false;  // grazing / looping: no useful prediction
    }
    out.valid = true;
    return true;
  }

  // Reference from a TPC-only SiTpcHelixFit result.  Straight-line fits are rejected (the TPC
  // measures the curvature).
  inline bool buildTpcReference(const SiTpcHelixFit::Result& fit, const FrameMap& map, TpcReference& ref,
                                double boundaryMargin = 0.3)
  {
    ref = TpcReference();
    if (fit.status != SiTpcHelixFit::Ok || !(fit.R > 0) || fit.R >= 0.5 * SiTpcHelixFit::kStraightRadius ||
        (fit.helicity != 1 && fit.helicity != -1))
    {
      return false;
    }
    ref.cx = fit.cx;
    ref.cy = fit.cy;
    ref.R = fit.R;
    ref.helicity = fit.helicity;
    ref.pcaX = fit.pcaX;
    ref.pcaY = fit.pcaY;
    ref.z0 = fit.z0;
    ref.tanl = fit.tanl;
    ref.phi = fit.phi;
    ref.dca = fit.dca;
    ref.pt = fit.pt;
    ref.charge = fit.charge;
    ref.valid = true;
    for (int layer = 0; layer < SiDetectorFrame::kNLayers; ++layer)
    {
      const double rl = SiDetectorFrame::layerRadius(layer);
      auto& p = ref.layers[layer];
      predictTpcLayer(ref, map, layer, rl, p);
      if (p.valid && map.boundaryDistance(p.u, rl) < boundaryMargin)
      {
        predictTpcLayer(ref, map, layer, rl, ref.alt[layer], map.otherHalf(p.half));
      }
    }
    return true;
  }
}  // namespace SiTpcReference

#endif
