// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCCHAINMATH_H
#define SITRACKRECO_SITPCCHAINMATH_H

#include "SiDetectorFrame.h"
#include "SiTpcBeamAlignment.h"
#include "SiTpcHelixFit.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// Geometry used by SiTpcHitChainMatcher, kept free of Fun4All so it can be unit tested.
//
// Frames
//   TPC helix  : TPC frame (SiTpcHelixFit on TPC clusters, beam = TPC beam axis).
//   hits       : silicon DETECTOR frame (SiDetectorFrame), where the layers are centred
//                cylinders and the hit tables are static; all window tests happen there.
//   Si chain   : TPC frame (the two clamshell halves are displaced by ~1 mm).  A TPC-frame point goes there as
//                  det = SiDetectorFrame::toDetectorFromGlobal(SiTpcBeamAlignment::siHalfFromTpcAxis(half, tpc))
//                i.e. a pure translation that depends on the clamshell half and (weakly) on z.
//
// TPC reference (AuAu: the TPC trajectory is known)
//   For every layer the helix is intersected with the layer cylinder once per track
//   (closed-form circle-circle intersection, two iterations for the half / z dependence of the
//   translation).  Stored per layer: point, u, path length s, z, the local incidence of the
//   track on the cylinder and the derivatives du/dr, dz/dr, ds/dr used to move the prediction
//   to the radius of a hit that is not exactly on the nominal cylinder.
//
// Silicon chain (TPC curvature and TPC theta)
//   The chain is NOT fitted freely: radius, helicity and tan(lambda) are those of the TPC
//   helix.  The silicon hits only anchor it:
//     xy : 1 point  -> the circle of TPC radius through the beam point and the point
//          2 points -> the circle of TPC radius through both points
//          >=3      -> least-squares centre with the radius fixed (Gauss-Newton)
//     z  : z(s) = z0 + tanl_TPC * s, z0 = weighted mean of z_i - tanl_TPC * s_i
//   Polar angles are measured from the TPC vertex (z0 of the TPC helix on the beam axis):
//     theta(point) = atan2(s, z - z0_TPC).
namespace SiTpcChainMath
{
  constexpr double kPi = M_PI;
  constexpr double kTwoPi = 2.0 * M_PI;

  inline double wrapPi(double a)
  {
    while (a > kPi)
    {
      a -= kTwoPi;
    }
    while (a <= -kPi)
    {
      a += kTwoPi;
    }
    return a;
  }
  inline double wrap2Pi(double a)
  {
    while (a >= kTwoPi)
    {
      a -= kTwoPi;
    }
    while (a < 0)
    {
      a += kTwoPi;
    }
    return a;
  }

  // ---- TPC frame <-> silicon detector frame
  struct FrameMap
  {
    const SiDetectorFrame* frame = nullptr;
    const SiTpcBeamAlignment* alignment = nullptr;

    SiDetectorFrame::Point tpcToDet(double x, double y, double z, int half) const
    {
      const auto g = alignment->siHalfFromTpcAxis(half, {x, y, z});
      return frame->toDetectorFromGlobal(g);
    }
    // inverse: detector frame point of clamshell half `half` -> TPC frame
    SiDetectorFrame::Point detToTpc(double x, double y, double z, int half) const
    {
      const auto g = frame->toGlobal({x, y, z});
      return alignment->siHalfToTpcAxis(half, g);
    }
    // det = tpc + (dx, dy) at fixed z and half
    void offset(int half, double z, double& dx, double& dy) const
    {
      const auto p = tpcToDet(0, 0, z, half);
      dx = p.x;
      dy = p.y;
    }
    double uphi(double xd, double yd) const { return wrap2Pi(std::atan2(yd, xd) - frame->uphiRotation()); }
    int half(double u) const { return alignment->siHalf(u); }
    int otherHalf(int h) const { return h == SiTpcBeamAlignment::SiHalfA ? SiTpcBeamAlignment::SiHalfB : SiTpcBeamAlignment::SiHalfA; }
    // r*phi distance [cm] of azimuth u on radius rl to the nearest clamshell boundary (B or B + pi)
    double boundaryDistance(double u, double rl) const
    {
      const double b = alignment->clamshellBoundary();
      return rl * std::min(std::abs(wrapPi(u - b)), std::abs(wrapPi(u - b - kPi)));
    }
    // TPC beam axis at z, in the detector frame
    void beamDet(double z, int half, double& bx, double& by) const
    {
      const auto& axis = alignment->tpcBeamAxis();
      const auto p = tpcToDet(axis.x(z), axis.y(z), z, half);
      bx = p.x;
      by = p.y;
    }
  };

  // ---- TPC reference
  struct LayerPrediction
  {
    bool valid = false;
    int layer = -1;
    int half = 0;
    double R = 0;               // nominal layer radius
    double x = 0, y = 0, z = 0;  // detector frame
    double u = 0;               // intrinsic azimuth
    double s = 0;               // transverse path from the TPC pca
    double tx = 0, ty = 0;      // unit direction of motion (xy)
    double cosA = 1, sinA = 0;  // incidence: direction vs radial, cos = t.n, sin = n x t
    double dudr = 0, dzdr = 0, dsdr = 1;
    double cx = 0, cy = 0;      // TPC circle centre translated into the detector frame
  };

  struct TpcReference
  {
    bool valid = false;
    // TPC frame helix (SiTpcHelixFit::Result)
    double cx = 0, cy = 0, R = 0;
    int helicity = 0;
    double pcaX = 0, pcaY = 0;
    double z0 = 0, tanl = 0, theta = 0;  // theta = atan2(1, tanl)
    double phi = 0, dca = 0, pt = 0;
    int charge = 0;
    std::array<LayerPrediction, SiDetectorFrame::kNLayers> layers{};
    // Near a clamshell boundary the half of the predicted point is ambiguous (the halves are
    // displaced by ~1 mm): alt = prediction with the other half's translation (valid only there).
    std::array<LayerPrediction, SiDetectorFrame::kNLayers> alt{};
  };

  // Intersection of the circle (ccx, ccy, R) with the circle of radius rl around (0, 0).
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
    auto p0 = map.tpcToDet(ref.pcaX + rl * tx0, ref.pcaY + rl * ty0, ref.z0 + ref.tanl * rl, SiTpcBeamAlignment::SiHalfA);
    int half = forceHalf >= 0 ? forceHalf : map.half(map.uphi(p0.x, p0.y));
    double zg = ref.z0 + ref.tanl * rl;
    double px = 0, py = 0, s = 0, ccx = 0, ccy = 0;
    bool found = false;
    for (int iter = 0; iter < 3; ++iter)
    {
      double dx = 0, dy = 0;
      map.offset(half, zg, dx, dy);
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
    out.cx = ccx;
    out.cy = ccy;
    out.tx = ref.helicity * (-(py - ccy)) / ref.R;
    out.ty = ref.helicity * (px - ccx) / ref.R;
    const double nx = px / rl, ny = py / rl;
    out.cosA = out.tx * nx + out.ty * ny;
    out.sinA = nx * out.ty - ny * out.tx;
    if (out.cosA < 1e-3)
    {
      return false;  // grazing / looping: no useful prediction
    }
    out.dsdr = 1.0 / out.cosA;
    out.dudr = out.sinA / (rl * out.cosA);
    out.dzdr = ref.tanl / out.cosA;
    out.valid = true;
    return true;
  }

  // Build the reference from a TPC-only SiTpcHelixFit result.  Straight-line fits are rejected
  // (the TPC measures curvature, and the chain needs it).
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
    ref.theta = std::atan2(1.0, fit.tanl);
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

  // ---- silicon chain with TPC radius / helicity / tanl
  // The chain lives in the TPC frame: the two clamshell halves are displaced w.r.t. each other
  // (~1 mm, SiTpcBeamAlignment), so a chain crossing the half boundary is only straight / circular
  // in the common frame.  Points are moved there with FrameMap::detToTpc (half of the point),
  // predictions are moved back to the detector frame of the predicted layer point.
  struct SiChainPoint
  {
    double x = 0, y = 0, z = 0;  // TPC frame
    double r = 0;
    double wz = 1;  // weight in the z0 mean (INTT strips: small)
  };

  struct SiChainFit
  {
    bool valid = false;
    double cx = 0, cy = 0;      // circle centre (radius = TPC radius), TPC frame
    double pcaX = 0, pcaY = 0;  // closest point of the circle to the beam point
    double bx = 0, by = 0;      // beam point used (TPC beam axis)
    double z0 = 0;              // z at the pca
    double dca = 0;             // |beam - centre| - R
  };

  // centre of the circle of radius R through a -> b (motion from a to b) with helicity hel
  inline bool circleThrough(double ax, double ay, double bx, double by, double R, int hel, double& cx, double& cy)
  {
    const double dx = bx - ax, dy = by - ay;
    const double chord = std::hypot(dx, dy);
    if (!(chord > 1e-9) || chord > 2 * R)
    {
      return false;
    }
    const double h = std::sqrt(std::max(0.0, R * R - 0.25 * chord * chord));
    cx = 0.5 * (ax + bx) + hel * h * (-dy / chord);
    cy = 0.5 * (ay + by) + hel * h * (dx / chord);
    return true;
  }

  inline double arcFromPca(const SiChainFit& f, double R, double x, double y)
  {
    const double c = std::hypot(x - f.pcaX, y - f.pcaY);
    return 2 * R * std::asin(std::min(1.0, c / (2 * R)));
  }

  // pts: any order.  previous: last fit of this chain (start point for >= 3 points), may be null.
  inline bool fitSiChain(std::vector<SiChainPoint> pts, double R, int hel, double tanl, double bx, double by,
                         const SiChainFit* previous, SiChainFit& out)
  {
    out = SiChainFit();
    out.bx = bx;
    out.by = by;
    if (pts.empty() || !(R > 0))
    {
      return false;
    }
    std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b) { return a.r < b.r; });
    double cx = 0, cy = 0;
    if (pts.size() == 1)
    {
      if (!circleThrough(bx, by, pts[0].x, pts[0].y, R, hel, cx, cy))
      {
        return false;
      }
    }
    else
    {
      if (!circleThrough(pts.front().x, pts.front().y, pts.back().x, pts.back().y, R, hel, cx, cy))
      {
        return false;
      }
      if (pts.size() >= 3)
      {
        if (previous && previous->valid)
        {
          cx = previous->cx;
          cy = previous->cy;
        }
        for (int iter = 0; iter < 3; ++iter)
        {
          double sxx = 0, sxy = 0, syy = 0, sxb = 0, syb = 0;
          for (const auto& p : pts)
          {
            const double dx = p.x - cx, dy = p.y - cy;
            const double d = std::hypot(dx, dy);
            if (!(d > 1e-9))
            {
              continue;
            }
            const double ux = dx / d, uy = dy / d, res = d - R;
            sxx += ux * ux;
            sxy += ux * uy;
            syy += uy * uy;
            sxb += ux * res;
            syb += uy * res;
          }
          const double det = sxx * syy - sxy * sxy;
          if (!(std::abs(det) > 1e-12))
          {
            break;
          }
          cx += (sxb * syy - syb * sxy) / det;
          cy += (sxx * syb - sxy * sxb) / det;
        }
      }
    }
    out.cx = cx;
    out.cy = cy;
    const double vx = bx - cx, vy = by - cy;
    const double vl = std::hypot(vx, vy);
    if (!(vl > 1e-9))
    {
      return false;
    }
    out.pcaX = cx + R * vx / vl;
    out.pcaY = cy + R * vy / vl;
    out.dca = vl - R;
    double sw = 0, swz = 0;
    for (const auto& p : pts)
    {
      const double s = arcFromPca(out, R, p.x, p.y);
      sw += p.wz;
      swz += p.wz * (p.z - tanl * s);
    }
    if (!(sw > 0))
    {
      return false;
    }
    out.z0 = swz / sw;
    out.valid = std::isfinite(out.cx) && std::isfinite(out.cy) && std::isfinite(out.z0);
    return out.valid;
  }

  struct SiPrediction
  {
    bool valid = false;
    double x = 0, y = 0;    // detector frame
    double xt = 0, yt = 0;  // TPC frame
    double z = 0, u = 0, s = 0;
  };

  // Chain point on the layer of tp (radius tp.R, translation of tp's clamshell half at tp.z);
  // of the two crossings the one nearest to the TPC prediction on that layer is taken.
  inline SiPrediction predictSiChain(const SiChainFit& f, double R, double tanl, const LayerPrediction& tp,
                                     const FrameMap& map)
  {
    SiPrediction p;
    if (!f.valid || !tp.valid)
    {
      return p;
    }
    double dx = 0, dy = 0;
    map.offset(tp.half, tp.z, dx, dy);
    double xs[2], ys[2];
    if (circleIntersections(f.cx + dx, f.cy + dy, R, tp.R, xs, ys) == 0)
    {
      return p;
    }
    const int k = (std::hypot(xs[0] - tp.x, ys[0] - tp.y) <= std::hypot(xs[1] - tp.x, ys[1] - tp.y)) ? 0 : 1;
    p.x = xs[k];
    p.y = ys[k];
    p.xt = p.x - dx;
    p.yt = p.y - dy;
    p.s = arcFromPca(f, R, p.xt, p.yt);
    p.z = f.z0 + tanl * p.s;
    p.u = map.uphi(p.x, p.y);
    p.valid = true;
    return p;
  }

  // direction of motion of a circle (centre cx, cy, helicity hel) at (x, y)
  inline void tangent(double cx, double cy, int hel, double x, double y, double& tx, double& ty)
  {
    const double dx = x - cx, dy = y - cy;
    const double d = std::hypot(dx, dy);
    tx = d > 0 ? hel * (-dy) / d : 0;
    ty = d > 0 ? hel * dx / d : 0;
  }

  inline double thetaFromVertex(double s, double z, double zVertex) { return std::atan2(s, z - zVertex); }
}  // namespace SiTpcChainMath

#endif
