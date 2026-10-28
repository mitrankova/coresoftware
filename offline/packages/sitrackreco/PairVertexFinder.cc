#include "PairVertexFinder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstddef>
#include <map>
#include <numeric>
#include <set>
#include <utility>

namespace
{
  constexpr double kPi = M_PI;

  double wrap_phi(double phi)
  {
    while (phi > kPi)
    {
      phi -= 2.0 * kPi;
    }
    while (phi <= -kPi)
    {
      phi += 2.0 * kPi;
    }
    return phi;
  }

  double median(std::vector<double> v)
  {
    if (v.empty())
    {
      return 0.0;
    }
    const std::size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(n), v.end());
    double med = v[n];
    if (v.size() % 2 == 0)
    {
      med = 0.5 * (med + *std::max_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(n)));
    }
    return med;
  }

  //! simple union-find
  struct DisjointSet
  {
    std::vector<unsigned int> parent;
    explicit DisjointSet(unsigned int n)
      : parent(n)
    {
      std::iota(parent.begin(), parent.end(), 0U);
    }
    unsigned int find(unsigned int i)
    {
      while (parent[i] != i)
      {
        parent[i] = parent[parent[i]];
        i = parent[i];
      }
      return i;
    }
    void unite(unsigned int a, unsigned int b)
    {
      a = find(a);
      b = find(b);
      if (a != b)
      {
        parent[b] = a;
      }
    }
  };
}  // namespace

//____________________________________________________________________________
bool PairVertexFinder::fromPerigee(unsigned int id, double pcaX, double pcaY, double z0, double phi, double tanl,
                                   double cx, double cy, double radius, int helicity, double pt, TrackModel& m)
{
  if (!std::isfinite(pcaX) || !std::isfinite(pcaY) || !std::isfinite(z0) || !std::isfinite(phi) || !std::isfinite(tanl))
  {
    return false;
  }
  m.track_id = id;
  m.x0 = pcaX;
  m.y0 = pcaY;
  m.z0 = z0;
  m.phi0 = phi;
  m.tanl = tanl;
  m.p = (std::isfinite(pt) && pt > 0.0) ? pt * std::sqrt(1.0 + tanl * tanl) : 0.0;
  const bool helix = std::isfinite(cx) && std::isfinite(cy) && std::isfinite(radius) && radius > 0.0 && radius < 1.0e4 &&
                     (helicity == 1 || helicity == -1);
  if (!helix)
  {
    m.is_line = true;
    m.dir = Eigen::Vector3d(std::cos(phi), std::sin(phi), tanl).normalized();
    return true;
  }
  m.is_line = false;
  m.xc = cx;
  m.yc = cy;
  m.radius = radius;
  m.h = helicity;
  return true;
}

//____________________________________________________________________________
bool PairVertexFinder::fromPointsWithCurvature(unsigned int id, const std::vector<Eigen::Vector3d>& points, double radius, int helicity,
                                               double tanl, double pt, TrackModel& m)
{
  if (points.size() < 2 || !std::isfinite(tanl))
  {
    return false;
  }
  // innermost and outermost point
  const auto byR = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b) { return a.head<2>().norm() < b.head<2>().norm(); };
  const Eigen::Vector3d p1 = *std::min_element(points.begin(), points.end(), byR);
  const Eigen::Vector3d p2 = *std::max_element(points.begin(), points.end(), byR);
  const Eigen::Vector2d chord = p2.head<2>() - p1.head<2>();
  const double L = chord.norm();
  if (L < 0.1)
  {
    return false;  // points on (nearly) one layer
  }
  const Eigen::Vector2d u = chord / L;
  m.track_id = id;
  m.tanl = tanl;
  m.p = (std::isfinite(pt) && pt > 0.0) ? pt * std::sqrt(1.0 + tanl * tanl) : 0.0;
  m.x0 = p1.x();
  m.y0 = p1.y();
  const bool helix = std::isfinite(radius) && radius > 0.5 * L && radius < 1.0e4 && (helicity == 1 || helicity == -1);
  // arc length of every point from p1, then z_ref
  auto arc = [&](const Eigen::Vector3d& q) -> double
  {
    if (!helix)
    {
      return (q.head<2>() - p1.head<2>()).dot(u);
    }
    const Eigen::Vector2d c(m.xc, m.yc);
    const Eigen::Vector2d a = p1.head<2>() - c, b = q.head<2>() - c;
    const double dphi = std::atan2(a.x() * b.y() - a.y() * b.x(), a.dot(b));  // signed angle a -> b
    return m.h * radius * dphi;
  };
  if (helix)
  {
    const Eigen::Vector2d mid = 0.5 * (p1.head<2>() + p2.head<2>());
    const Eigen::Vector2d left(-u.y(), u.x());
    const Eigen::Vector2d c = mid + helicity * std::sqrt(radius * radius - 0.25 * L * L) * left;  // centre left of the motion for CCW
    m.is_line = false;
    m.xc = c.x();
    m.yc = c.y();
    m.radius = radius;
    m.h = helicity;
    const Eigen::Vector2d r1 = (p1.head<2>() - c) / radius;
    const Eigen::Vector2d t = helicity * Eigen::Vector2d(-r1.y(), r1.x());  // direction of motion at p1
    m.phi0 = std::atan2(t.y(), t.x());
  }
  else
  {
    m.is_line = true;
    m.phi0 = std::atan2(u.y(), u.x());
    m.dir = Eigen::Vector3d(u.x(), u.y(), tanl).normalized();
  }
  // z at p1 from every point, weighted by its z resolution: MVTX pixels (r < kMaxPixelR) ~0.003 cm,
  // INTT strips ~0.46 cm (weight kStripZWeight ~ (0.003/0.46)^2); with equal weights the INTT strips
  // moved z_ref by ~1-2 mm (Si+TPC vertex z much worse than the Si-only one)
  constexpr double kMaxPixelR = 6.0;       // [cm] MVTX layers r < ~5 cm, INTT > ~7 cm
  constexpr double kStripZWeight = 1.0e-4;
  double sumWZ = 0, sumW = 0;
  for (const auto& q : points)
  {
    const double w = q.head<2>().norm() < kMaxPixelR ? 1.0 : kStripZWeight;
    sumWZ += w * (q.z() - tanl * arc(q));
    sumW += w;
  }
  m.z0 = sumW > 0 ? sumWZ / sumW : std::numeric_limits<double>::quiet_NaN();
  return std::isfinite(m.z0) && std::isfinite(m.phi0);
}

//____________________________________________________________________________
bool PairVertexFinder::acceptAtBeam(const TrackModel& m) const
{
  const LinearTrack lin = linearize(m, Eigen::Vector3d(m_cfg.beamX, m_cfg.beamY, 0.0));
  if (!lin.ok)
  {
    return false;
  }
  const double dcaxy = std::hypot(lin.point.x() - m_cfg.beamX, lin.point.y() - m_cfg.beamY);
  return dcaxy <= m_cfg.maxDcaXY && std::fabs(lin.point.z()) <= m_cfg.maxAbsZ;
}

//____________________________________________________________________________
// Point of the track closest to `target` in the transverse plane, and the unit tangent there.
// Residuals are projected perpendicular to the tangent, so the remaining along-track offset
// does not matter.
PairVertexFinder::LinearTrack PairVertexFinder::linearize(const TrackModel& m, const Eigen::Vector3d& target) const
{
  LinearTrack lin;
  if (m.is_line)
  {
    const Eigen::Vector3d p0(m.x0, m.y0, m.z0);
    const Eigen::Vector3d diff = target - p0;
    const double dxy2 = m.dir.x() * m.dir.x() + m.dir.y() * m.dir.y();
    const double s = dxy2 > 1.0e-12 ? (diff.x() * m.dir.x() + diff.y() * m.dir.y()) / dxy2 : diff.dot(m.dir);
    lin.point = p0 + s * m.dir;
    lin.dir = m.dir;
    lin.ok = lin.point.allFinite();
    return lin;
  }
  // helix: angle of the target around the circle centre
  const double psi = std::atan2(target.y() - m.yc, target.x() - m.xc);
  const double phi_t = psi + m.h * 0.5 * kPi;  // tangent angle there
  const double dphi = wrap_phi(phi_t - m.phi0);
  const double s = m.h * m.radius * dphi;  // transverse arc length from the reference point
  lin.point = Eigen::Vector3d(m.xc + m.radius * std::cos(psi), m.yc + m.radius * std::sin(psi), m.z0 + s * m.tanl);
  lin.dir = Eigen::Vector3d(std::cos(phi_t), std::sin(phi_t), m.tanl).normalized();
  lin.ok = lin.point.allFinite() && lin.dir.allFinite();
  return lin;
}

//____________________________________________________________________________
// 3x3 weight matrix W = B C^-1 B^T, B = [u w] spans the plane normal to the track; u is
// transverse, w = dir x u. C is the 2x2 position covariance of the track at the PCA in that plane.
Eigen::Matrix3d PairVertexFinder::weightMatrix(const TrackModel& m, const LinearTrack& lin) const
{
  const Eigen::Vector3d& d = lin.dir;
  Eigen::Vector3d u(-d.y(), d.x(), 0.0);
  if (u.norm() < 1.0e-9)
  {
    u = Eigen::Vector3d(1.0, 0.0, 0.0);
  }
  u.normalize();
  const Eigen::Vector3d w = d.cross(u).normalized();

  Eigen::Matrix<double, 3, 2> B;
  B.col(0) = u;
  B.col(1) = w;

  Eigen::Matrix2d C;
  if (m_cfg.useTrackCovariance && m.has_cov)
  {
    C = B.transpose() * m.cov * B;
  }
  else
  {
    const double srphi = m.sigmaRPhi > 0.0 ? m.sigmaRPhi : m_cfg.sigmaRPhi;
    const double sz = m.sigmaZ > 0.0 ? m.sigmaZ : m_cfg.sigmaZ;
    C.setZero();
    C(0, 0) = srphi * srphi;
    C(1, 1) = sz * sz * w.z() * w.z();
  }
  if (m.p > 0.0 && m_cfg.msTerm > 0.0)
  {
    const double ms = m_cfg.msTerm / m.p;
    C(0, 0) += ms * ms;
    C(1, 1) += ms * ms;
  }
  C(0, 0) = std::max(C(0, 0), 1.0e-8);
  C(1, 1) = std::max(C(1, 1), 1.0e-8);
  return B * C.inverse() * B.transpose();
}

//____________________________________________________________________________
double PairVertexFinder::dcaTwoLines(const LinearTrack& t1, const LinearTrack& t2, Eigen::Vector3d& pca1, Eigen::Vector3d& pca2)
{
  const Eigen::Vector3d w0 = t1.point - t2.point;
  const double b = t1.dir.dot(t2.dir);
  const double d = t1.dir.dot(w0);
  const double e = t2.dir.dot(w0);
  const double den = 1.0 - b * b;  // directions are unit vectors
  double s = 0.0;
  double t = 0.0;
  if (den < 1.0e-12)
  {
    t = e;  // parallel: project point 1 onto line 2
  }
  else
  {
    s = (b * e - d) / den;
    t = (e - b * d) / den;
  }
  pca1 = t1.point + s * t1.dir;
  pca2 = t2.point + t * t2.dir;
  return (pca1 - pca2).norm();
}

//____________________________________________________________________________
std::vector<PairVertexFinder::TrackPair> PairVertexFinder::findPairs(const std::vector<TrackModel>& models,
                                                                      const std::vector<LinearTrack>& lin, double dcacut) const
{
  std::vector<TrackPair> pairs;
  const auto n = static_cast<unsigned int>(lin.size());
  const double loose = 3.0 * dcacut;
  for (unsigned int i = 0; i < n; ++i)
  {
    for (unsigned int j = i + 1; j < n; ++j)
    {
      Eigen::Vector3d p1;
      Eigen::Vector3d p2;
      double dca = dcaTwoLines(lin[i], lin[j], p1, p2);
      if (dca > loose)
      {
        continue;
      }
      // refine: re-linearise both tracks at the pair midpoint
      for (int k = 0; k < 2; ++k)
      {
        const Eigen::Vector3d mid = 0.5 * (p1 + p2);
        const LinearTrack li = linearize(models[i], mid);
        const LinearTrack lj = linearize(models[j], mid);
        if (!li.ok || !lj.ok)
        {
          break;
        }
        dca = dcaTwoLines(li, lj, p1, p2);
      }
      if (dca > dcacut)
      {
        continue;
      }
      const Eigen::Vector3d mid = 0.5 * (p1 + p2);
      if (std::hypot(mid.x() - m_cfg.beamX, mid.y() - m_cfg.beamY) > m_cfg.maxDcaXY || std::fabs(mid.z()) > m_cfg.maxAbsZ)
      {
        continue;
      }
      TrackPair pair;
      pair.a = i;
      pair.b = j;
      pair.dca = dca;
      pair.pca_a = p1;
      pair.pca_b = p2;
      pairs.push_back(pair);
    }
  }
  return pairs;
}

//____________________________________________________________________________
std::vector<std::vector<unsigned int>> PairVertexFinder::connectedComponents(unsigned int nnodes, const std::vector<TrackPair>& pairs)
{
  DisjointSet ds(nnodes);
  std::vector<bool> used(nnodes, false);
  for (const auto& p : pairs)
  {
    ds.unite(p.a, p.b);
    used[p.a] = true;
    used[p.b] = true;
  }
  std::vector<std::vector<unsigned int>> comps;
  std::vector<int> root_to_comp(nnodes, -1);
  for (unsigned int i = 0; i < nnodes; ++i)
  {
    if (!used[i])
    {
      continue;
    }
    const unsigned int r = ds.find(i);
    if (root_to_comp[r] < 0)
    {
      root_to_comp[r] = static_cast<int>(comps.size());
      comps.emplace_back();
    }
    comps[static_cast<std::size_t>(root_to_comp[r])].push_back(i);
  }
  return comps;
}

//____________________________________________________________________________
// For every candidate, keep the pairs whose midpoint lies within the outlier cut of the
// candidate median. Rejected pairs are re-examined (they may belong to a second vertex merged
// into the same candidate), but tracks already claimed by a kept pair are not reused.
std::vector<PairVertexFinder::TrackPair> PairVertexFinder::removeOutlierPairs(unsigned int nnodes, const std::vector<TrackPair>& pairs) const
{
  std::vector<TrackPair> result;
  std::vector<TrackPair> remaining = pairs;
  std::set<unsigned int> claimed;
  for (int pass = 0; pass < 5 && !remaining.empty(); ++pass)
  {
    if (pass > 0)
    {
      remaining.erase(std::remove_if(remaining.begin(), remaining.end(), [&claimed](const TrackPair& p)
                                     { return claimed.count(p.a) != 0 || claimed.count(p.b) != 0; }),
                      remaining.end());
    }
    DisjointSet ds(nnodes);
    for (const auto& p : remaining)
    {
      ds.unite(p.a, p.b);
    }
    std::map<unsigned int, std::vector<const TrackPair*>> by_root;
    for (const auto& p : remaining)
    {
      by_root[ds.find(p.a)].push_back(&p);
    }
    std::vector<TrackPair> rejected;
    std::vector<TrackPair> kept_this_pass;
    for (const auto& [root, plist] : by_root)
    {
      std::vector<double> xs;
      std::vector<double> ys;
      std::vector<double> zs;
      for (const auto* p : plist)
      {
        const Eigen::Vector3d mid = p->midpoint();
        xs.push_back(mid.x());
        ys.push_back(mid.y());
        zs.push_back(mid.z());
      }
      const Eigen::Vector3d med(median(xs), median(ys), median(zs));
      for (const auto* p : plist)
      {
        if ((p->midpoint() - med).norm() < m_cfg.outlierPairCut)
        {
          kept_this_pass.push_back(*p);
        }
        else
        {
          rejected.push_back(*p);
        }
      }
    }
    for (const auto& p : kept_this_pass)
    {
      claimed.insert(p.a);
      claimed.insert(p.b);
      result.push_back(p);
    }
    remaining.swap(rejected);
  }
  return result;
}

//____________________________________________________________________________
double PairVertexFinder::trackChi2(const TrackModel& m, const Eigen::Vector3d& vtx) const
{
  const LinearTrack lin = linearize(m, vtx);
  if (!lin.ok)
  {
    return 1.0e30;
  }
  const Eigen::Vector3d r = vtx - lin.point;
  return r.dot(weightMatrix(m, lin) * r);
}

//____________________________________________________________________________
PairVertexFinder::VertexFit PairVertexFinder::fitVertex(const std::vector<TrackModel>& models, std::vector<unsigned int> tracks,
                                                        const Eigen::Vector3d& seed) const
{
  VertexFit fit;
  Eigen::Vector3d vtx = seed;

  // one Gauss-Newton solve: tracks linearised at `at`
  auto solve = [&](const std::vector<unsigned int>& active, const Eigen::Vector3d& at, Eigen::Vector3d& out, Eigen::Matrix3d& A) -> bool
  {
    A.setZero();
    Eigen::Vector3d b = Eigen::Vector3d::Zero();
    unsigned int nused = 0;
    for (const unsigned int it : active)
    {
      const TrackModel& m = models[it];
      const LinearTrack lin = linearize(m, at);
      if (!lin.ok)
      {
        continue;
      }
      const Eigen::Matrix3d W = weightMatrix(m, lin);
      A += W;
      b += W * lin.point;
      ++nused;
    }
    if (m_cfg.useBeamSpotConstraint)
    {
      const double wx = 1.0 / (m_cfg.beamSigmaX * m_cfg.beamSigmaX);
      const double wy = 1.0 / (m_cfg.beamSigmaY * m_cfg.beamSigmaY);
      A(0, 0) += wx;
      A(1, 1) += wy;
      b(0) += wx * m_cfg.beamX;
      b(1) += wy * m_cfg.beamY;
    }
    if (nused < 2)
    {
      return false;
    }
    const Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
    if (!lu.isInvertible())
    {
      return false;
    }
    out = lu.solve(b);
    return out.allFinite();
  };

  Eigen::Matrix3d A;
  while (true)
  {
    if (tracks.size() < m_cfg.minTracksPerVertex)
    {
      return fit;
    }
    for (unsigned int iter = 0; iter < m_cfg.maxIterations; ++iter)
    {
      Eigen::Vector3d next;
      if (!solve(tracks, vtx, next, A))
      {
        return fit;
      }
      const double shift = (next - vtx).norm();
      vtx = next;
      if (shift < m_cfg.convergenceTol)
      {
        break;
      }
    }
    // drop the worst track while it is incompatible
    double worst_chi2 = -1.0;
    std::size_t worst = 0;
    for (std::size_t k = 0; k < tracks.size(); ++k)
    {
      const double c = trackChi2(models[tracks[k]], vtx);
      if (c > worst_chi2)
      {
        worst_chi2 = c;
        worst = k;
      }
    }
    if (worst_chi2 > m_cfg.maxTrackChi2 && tracks.size() > m_cfg.minTracksPerVertex)
    {
      tracks.erase(tracks.begin() + static_cast<std::ptrdiff_t>(worst));
      continue;
    }
    if (worst_chi2 > m_cfg.maxTrackChi2)
    {
      return fit;  // only the minimum number of tracks left and still incompatible
    }
    break;
  }

  Eigen::Vector3d dummy;
  if (!solve(tracks, vtx, dummy, A))
  {
    return fit;
  }
  fit.position = vtx;
  fit.covariance = A.inverse();
  fit.chi2 = 0.0;
  for (const unsigned int it : tracks)
  {
    fit.chi2 += trackChi2(models[it], vtx);
  }
  if (m_cfg.useBeamSpotConstraint)
  {
    const double dx = (vtx.x() - m_cfg.beamX) / m_cfg.beamSigmaX;
    const double dy = (vtx.y() - m_cfg.beamY) / m_cfg.beamSigmaY;
    fit.chi2 += dx * dx + dy * dy;
  }
  fit.ndf = 2 * static_cast<int>(tracks.size()) - 3 + (m_cfg.useBeamSpotConstraint ? 2 : 0);
  fit.tracks = std::move(tracks);
  fit.ok = fit.position.allFinite() && fit.covariance.allFinite();
  return fit;
}

//____________________________________________________________________________
PairVertexFinder::VertexFit PairVertexFinder::fit(const std::vector<TrackModel>& models, const Eigen::Vector3d& seed) const
{
  std::vector<unsigned int> tracks(models.size());
  for (unsigned int i = 0; i < tracks.size(); ++i)
  {
    tracks[i] = i;
  }
  return fitVertex(models, std::move(tracks), seed);
}

//____________________________________________________________________________
std::vector<PairVertexFinder::VertexFit> PairVertexFinder::find(const std::vector<TrackModel>& models, Stats* stats) const
{
  std::vector<VertexFit> vertices;
  if (models.size() < std::max(2U, m_cfg.minTracksPerVertex))
  {
    return vertices;
  }
  const auto n = static_cast<unsigned int>(models.size());

  // 2. linearise at the beam spot and find compatible pairs
  const Eigen::Vector3d beam(m_cfg.beamX, m_cfg.beamY, 0.0);
  std::vector<LinearTrack> lin;
  lin.reserve(n);
  for (const auto& m : models)
  {
    lin.push_back(linearize(m, beam));
  }
  std::vector<TrackPair> pairs = findPairs(models, lin, m_cfg.pairDcaCut);
  bool loose = false;
  if (pairs.empty() && m_cfg.retryLooseDcaCut)
  {
    pairs = findPairs(models, lin, 3.0 * m_cfg.pairDcaCut);
    loose = true;
  }

  // 3. candidates = connected components after outlier-pair removal
  const std::vector<TrackPair> clean = removeOutlierPairs(n, pairs);
  const auto comps = connectedComponents(n, clean);
  if (stats)
  {
    stats->pairs = static_cast<unsigned int>(pairs.size());
    stats->cleanPairs = static_cast<unsigned int>(clean.size());
    stats->candidates = static_cast<unsigned int>(comps.size());
    stats->looseRetry = loose;
  }

  // 4. fit candidates, seeded at the median pair midpoint
  for (const auto& comp : comps)
  {
    if (comp.size() < m_cfg.minTracksPerVertex)
    {
      continue;
    }
    const std::set<unsigned int> members(comp.begin(), comp.end());
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
    for (const auto& p : clean)
    {
      if (members.count(p.a) != 0)
      {
        const Eigen::Vector3d mid = p.midpoint();
        xs.push_back(mid.x());
        ys.push_back(mid.y());
        zs.push_back(mid.z());
      }
    }
    const Eigen::Vector3d seed(median(xs), median(ys), median(zs));
    VertexFit fit = fitVertex(models, comp, seed);
    if (fit.ok)
    {
      vertices.push_back(std::move(fit));
    }
  }

  // 5. unused tracks join their most compatible vertex, which is refitted
  if (m_cfg.absorbUnusedTracks && !vertices.empty())
  {
    std::vector<bool> used(n, false);
    for (const auto& v : vertices)
    {
      for (const unsigned int it : v.tracks)
      {
        used[it] = true;
      }
    }
    std::vector<bool> changed(vertices.size(), false);
    for (unsigned int it = 0; it < n; ++it)
    {
      if (used[it])
      {
        continue;
      }
      double best = m_cfg.maxTrackChi2;
      int ibest = -1;
      for (std::size_t iv = 0; iv < vertices.size(); ++iv)
      {
        const double c = trackChi2(models[it], vertices[iv].position);
        if (c < best)
        {
          best = c;
          ibest = static_cast<int>(iv);
        }
      }
      if (ibest >= 0)
      {
        vertices[static_cast<std::size_t>(ibest)].tracks.push_back(it);
        changed[static_cast<std::size_t>(ibest)] = true;
      }
    }
    for (std::size_t iv = 0; iv < vertices.size(); ++iv)
    {
      if (!changed[iv])
      {
        continue;
      }
      VertexFit refit = fitVertex(models, vertices[iv].tracks, vertices[iv].position);
      if (refit.ok)
      {
        vertices[iv] = std::move(refit);
      }
    }
  }

  std::sort(vertices.begin(), vertices.end(), [](const VertexFit& a, const VertexFit& b) { return a.tracks.size() > b.tracks.size(); });
  return vertices;
}
