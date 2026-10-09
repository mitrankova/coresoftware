#include "SiTpcHelixFit.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
  constexpr double kPi = 3.14159265358979323846;

  double wrapPi(double a)
  {
    a = std::fmod(a + kPi, 2.0 * kPi);
    if (a < 0)
    {
      a += 2.0 * kPi;
    }
    return a - kPi;
  }

  // solves the symmetric n x n system (n <= 3) A p = b by Gaussian elimination
  bool solve(int n, double A[3][3], double b[3], double p[3])
  {
    for (int c = 0; c < n; ++c)
    {
      int piv = c;
      for (int r = c + 1; r < n; ++r)
      {
        if (std::abs(A[r][c]) > std::abs(A[piv][c]))
        {
          piv = r;
        }
      }
      if (std::abs(A[piv][c]) < 1e-12)
      {
        return false;
      }
      std::swap(A[c], A[piv]);
      std::swap(b[c], b[piv]);
      for (int r = c + 1; r < n; ++r)
      {
        const double f = A[r][c] / A[c][c];
        for (int k = c; k < n; ++k)
        {
          A[r][k] -= f * A[c][k];
        }
        b[r] -= f * b[c];
      }
    }
    for (int r = n - 1; r >= 0; --r)
    {
      double v = b[r];
      for (int k = r + 1; k < n; ++k)
      {
        v -= A[r][k] * p[k];
      }
      p[r] = v / A[r][r];
    }
    return true;
  }
}  // namespace

namespace SiTpcHelixFit
{
  // Weighted Taubin algebraic circle fit, Newton iteration on the characteristic polynomial
  // (N. Chernov, "Circular and linear regression", CircleFitByTaubin).
  Circle fitCircleTaubin(const std::vector<double>& x, const std::vector<double>& y,
                         const std::vector<double>& w)
  {
    Circle out;
    const std::size_t n = std::min({x.size(), y.size(), w.size()});
    std::size_t nUsed = 0;
    double sw = 0, mx = 0, my = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (w[i] <= 0)
      {
        continue;
      }
      ++nUsed;
      sw += w[i];
      mx += w[i] * x[i];
      my += w[i] * y[i];
    }
    if (nUsed < 3 || sw <= 0)
    {
      return out;
    }
    mx /= sw;
    my /= sw;

    double Mxx = 0, Myy = 0, Mxy = 0, Mxz = 0, Myz = 0, Mzz = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (w[i] <= 0)
      {
        continue;
      }
      const double xi = x[i] - mx, yi = y[i] - my, zi = xi * xi + yi * yi;
      Mxy += w[i] * xi * yi;
      Mxx += w[i] * xi * xi;
      Myy += w[i] * yi * yi;
      Mxz += w[i] * xi * zi;
      Myz += w[i] * yi * zi;
      Mzz += w[i] * zi * zi;
    }
    Mxx /= sw;
    Myy /= sw;
    Mxy /= sw;
    Mxz /= sw;
    Myz /= sw;
    Mzz /= sw;

    const double Mz = Mxx + Myy;
    const double covXY = Mxx * Myy - Mxy * Mxy;
    const double varZ = Mzz - Mz * Mz;
    const double A3 = 4.0 * Mz;
    const double A2 = -3.0 * Mz * Mz - Mzz;
    const double A1 = varZ * Mz + 4.0 * covXY * Mz - Mxz * Mxz - Myz * Myz;
    const double A0 = Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy) - varZ * covXY;
    const double A22 = A2 + A2;
    const double A33 = A3 + A3 + A3;

    double xn = 0.0, yn = A0;
    for (int iter = 0; iter < 99; ++iter)
    {
      const double dy = A1 + xn * (A22 + A33 * xn);
      if (dy == 0)
      {
        break;
      }
      const double xnew = xn - yn / dy;
      if (xnew == xn || !std::isfinite(xnew))
      {
        break;
      }
      const double ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
      if (std::abs(ynew) >= std::abs(yn))
      {
        break;
      }
      xn = xnew;
      yn = ynew;
    }

    const double det = xn * xn - xn * Mz + covXY;
    if (std::abs(det) < 1e-30)
    {
      return out;  // collinear points
    }
    const double xc = (Mxz * (Myy - xn) - Myz * Mxy) / det / 2.0;
    const double yc = (Myz * (Mxx - xn) - Mxz * Mxy) / det / 2.0;

    out.x = xc + mx;
    out.y = yc + my;
    out.r = std::sqrt(xc * xc + yc * yc + Mz);
    out.ok = std::isfinite(out.x) && std::isfinite(out.y) && std::isfinite(out.r);
    return out;
  }

  Result fit(const std::vector<Point>& pts, const Config& cfg, unsigned int minPoints)
  {
    Result res;
    const std::size_t n = pts.size();
    std::size_t nxy = 0;
    for (const auto& p : pts)
    {
      nxy += p.wxy > 0;
    }
    if (nxy < std::max(3U, minPoints))
    {
      res.status = TooFewPoints;
      return res;
    }

    std::vector<double> xs, ys, ws;
    double rmin = std::numeric_limits<double>::max(), rmax = 0;
    int iFirst = -1, iLast = -1;  // innermost / outermost point used in xy (input order)
    for (std::size_t i = 0; i < n; ++i)
    {
      if (pts[i].wxy <= 0)
      {
        continue;
      }
      xs.push_back(pts[i].x);
      ys.push_back(pts[i].y);
      ws.push_back(pts[i].wxy);
      const double r = std::hypot(pts[i].x - cfg.beamX, pts[i].y - cfg.beamY);
      rmin = std::min(rmin, r);
      rmax = std::max(rmax, r);
      if (iFirst < 0)
      {
        iFirst = static_cast<int>(i);
      }
      iLast = static_cast<int>(i);
    }
    const double x1 = pts[iFirst].x, y1 = pts[iFirst].y, x2 = pts[iLast].x, y2 = pts[iLast].y;
    if (cfg.beamWeight > 0)
    {
      xs.push_back(cfg.beamX);
      ys.push_back(cfg.beamY);
      ws.push_back(cfg.beamWeight);
    }
    const bool straight = (rmax - rmin) < cfg.minLeverArm;

    res.s.assign(n, 0.0);
    if (!straight)
    {
      const Circle c = fitCircleTaubin(xs, ys, ws);
      if (!c.ok)
      {
        res.status = Degenerate;
        return res;
      }
      const double ax = x1 - c.x, ay = y1 - c.y, bx = x2 - c.x, by = y2 - c.y;
      const int h = (ax * by - ay * bx) >= 0 ? +1 : -1;
      const double dx = cfg.beamX - c.x, dy = cfg.beamY - c.y;
      const double dist = std::hypot(dx, dy);
      if (dist <= 0)
      {
        res.status = Degenerate;
        return res;
      }
      const double ux = dx / dist, uy = dy / dist;  // centre -> beam
      res.pcaX = c.x + c.r * ux;
      res.pcaY = c.y + c.r * uy;
      res.tx = -h * uy;
      res.ty = h * ux;
      res.cx = c.x;
      res.cy = c.y;
      res.R = c.r;
      res.helicity = h;
      res.dca = ((res.pcaX - cfg.beamX) * res.ty - (res.pcaY - cfg.beamY) * res.tx >= 0 ? -1.0 : 1.0) * std::abs(dist - c.r);
      res.pt = 0.003 * std::abs(cfg.bz) * c.r;
      res.charge = cfg.bz == 0 ? 0 : -h * (cfg.bz > 0 ? 1 : -1);
      const double a0 = std::atan2(res.pcaY - c.y, res.pcaX - c.x);
      double ss = 0, sw = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        res.s[i] = c.r * h * wrapPi(std::atan2(pts[i].y - c.y, pts[i].x - c.x) - a0);
        if (pts[i].wxy > 0)
        {
          const double d = std::hypot(pts[i].x - c.x, pts[i].y - c.y) - c.r;
          ss += d * d;
          sw += 1;
        }
      }
      res.circleRms = std::sqrt(ss / sw);
    }
    else
    {
      // weighted total least squares line (principal axis)
      double sw = 0, mx = 0, my = 0;
      for (std::size_t i = 0; i < xs.size(); ++i)
      {
        sw += ws[i];
        mx += ws[i] * xs[i];
        my += ws[i] * ys[i];
      }
      mx /= sw;
      my /= sw;
      double sxx = 0, syy = 0, sxy = 0;
      for (std::size_t i = 0; i < xs.size(); ++i)
      {
        const double dx = xs[i] - mx, dy = ys[i] - my;
        sxx += ws[i] * dx * dx;
        syy += ws[i] * dy * dy;
        sxy += ws[i] * dx * dy;
      }
      const double ang = 0.5 * std::atan2(2 * sxy, sxx - syy);
      double tx = std::cos(ang), ty = std::sin(ang);
      if ((x2 - x1) * tx + (y2 - y1) * ty < 0)  // outward
      {
        tx = -tx;
        ty = -ty;
      }
      const double proj = (cfg.beamX - mx) * tx + (cfg.beamY - my) * ty;
      res.pcaX = mx + proj * tx;
      res.pcaY = my + proj * ty;
      res.tx = tx;
      res.ty = ty;
      // stored as a circle of very large radius, helicity +1 (centre to the left)
      res.cx = res.pcaX - kStraightRadius * ty;
      res.cy = res.pcaY + kStraightRadius * tx;
      res.R = kStraightRadius;
      res.helicity = +1;
      res.dca = ((res.pcaX - cfg.beamX) * ty - (res.pcaY - cfg.beamY) * tx >= 0 ? -1.0 : 1.0) *
                std::hypot(res.pcaX - cfg.beamX, res.pcaY - cfg.beamY);
      res.pt = std::numeric_limits<double>::quiet_NaN();
      res.charge = 0;
      double ss = 0, sn = 0;
      for (std::size_t i = 0; i < n; ++i)
      {
        res.s[i] = (pts[i].x - res.pcaX) * tx + (pts[i].y - res.pcaY) * ty;
        if (pts[i].wxy > 0)
        {
          const double d = (pts[i].x - res.pcaX) * (-ty) + (pts[i].y - res.pcaY) * tx;
          ss += d * d;
          sn += 1;
        }
      }
      res.circleRms = std::sqrt(ss / sn);
    }
    res.phi = std::atan2(res.ty, res.tx);

    // ---- z(s): weighted linear least squares in (z0, tanl[, offset of group 1])
    std::size_t nz0 = 0, nz1 = 0;
    for (const auto& p : pts)
    {
      if (p.wz > 0)
      {
        (p.group == 1 ? nz1 : nz0)++;
      }
    }
    const bool useOffset = cfg.groupZOffset && nz0 >= 1 && nz1 >= 1 && nz0 + nz1 >= 3;
    const int np = useOffset ? 3 : 2;
    double A[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    double b[3] = {0, 0, 0};
    for (std::size_t i = 0; i < n; ++i)
    {
      if (pts[i].wz <= 0)
      {
        continue;
      }
      const double g[3] = {1.0, res.s[i], (useOffset && pts[i].group == 1) ? 1.0 : 0.0};
      for (int r = 0; r < np; ++r)
      {
        for (int k = 0; k < np; ++k)
        {
          A[r][k] += pts[i].wz * g[r] * g[k];
        }
        b[r] += pts[i].wz * g[r] * pts[i].z;
      }
    }
    double p[3] = {0, 0, 0};
    if (!solve(np, A, b, p))
    {
      res.status = Degenerate;
      return res;
    }
    res.z0 = p[0];
    res.tanl = p[1];
    res.zOffset = useOffset ? p[2] : 0.0;
    double zz = 0, zn = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
      if (pts[i].wz <= 0)
      {
        continue;
      }
      const double d = pts[i].z - (res.z0 + res.tanl * res.s[i] + (pts[i].group == 1 ? res.zOffset : 0.0));
      zz += d * d;
      zn += 1;
    }
    res.zRms = std::sqrt(zz / zn);
    res.status = straight ? StraightLine : Ok;
    return res;
  }

  Pos at(const Result& r, double s)
  {
    if (r.R >= kStraightRadius)
    {
      return {r.pcaX + s * r.tx, r.pcaY + s * r.ty, r.z0 + r.tanl * s};
    }
    const double a0 = std::atan2(r.pcaY - r.cy, r.pcaX - r.cx);
    const double a = a0 + r.helicity * s / r.R;
    return {r.cx + r.R * std::cos(a), r.cy + r.R * std::sin(a), r.z0 + r.tanl * s};
  }

  bool atRadius(const Result& r, double rad, double zRef, Pos& out, double& sOut)
  {
    if (!r.fitted() || rad <= 0)
    {
      return false;
    }
    double cand[2];
    int nc = 0;
    if (r.R >= kStraightRadius)
    {
      // |pca + s t| = rad
      const double bq = r.pcaX * r.tx + r.pcaY * r.ty;
      const double cq = r.pcaX * r.pcaX + r.pcaY * r.pcaY - rad * rad;
      const double disc = bq * bq - cq;
      if (disc < 0)
      {
        return false;
      }
      cand[nc++] = -bq + std::sqrt(disc);
      cand[nc++] = -bq - std::sqrt(disc);
    }
    else
    {
      // intersection of the track circle with the circle of radius rad around (0, 0)
      const double d = std::hypot(r.cx, r.cy);
      if (d <= 0 || d > rad + r.R || d < std::abs(rad - r.R))
      {
        return false;
      }
      const double a = (rad * rad - r.R * r.R + d * d) / (2 * d);
      const double h = std::sqrt(std::max(0.0, rad * rad - a * a));
      const double ux = r.cx / d, uy = r.cy / d;
      const double a0 = std::atan2(r.pcaY - r.cy, r.pcaX - r.cx);
      for (int sgn : {-1, 1})
      {
        const double px = a * ux - sgn * h * uy, py = a * uy + sgn * h * ux;
        cand[nc++] = r.R * r.helicity * wrapPi(std::atan2(py - r.cy, px - r.cx) - a0);
      }
    }
    // The circle of radius rad is crossed once on the outgoing branch (s > 0) and once behind
    // the pca (s < 0, opposite side of the beam).  Points of a track from the beam are on the
    // outgoing branch; at small radius both crossings can have a similar z, so the outgoing
    // one is preferred and z only decides between crossings on the same side.
    constexpr double kBackTolerance = 0.5;  // [cm] allows |dca|-size negative s
    bool anyForward = false;
    for (int i = 0; i < nc; ++i)
    {
      anyForward |= cand[i] >= -kBackTolerance;
    }
    double best = std::numeric_limits<double>::max();
    for (int i = 0; i < nc; ++i)
    {
      if (anyForward && cand[i] < -kBackTolerance)
      {
        continue;
      }
      const double dz = std::abs(r.z0 + r.tanl * cand[i] - zRef);
      if (dz < best)
      {
        best = dz;
        sOut = cand[i];
      }
    }
    out = at(r, sOut);
    return true;
  }
}  // namespace SiTpcHelixFit
