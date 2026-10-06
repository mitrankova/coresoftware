#include "Tpc_PolyTrackVertexFinder.h"

#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include "Tpc_PolyTrackVertexContainer.h"

#include <globalvertex/SvtxVertexMap_v1.h>
#include <globalvertex/SvtxVertex_v3.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

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

  bool good(double x)
  {
    return std::isfinite(x) && std::fabs(x) < 1.0e30;
  }

  double median(std::vector<double> v)
  {
    if (v.empty())
    {
      return 0.0;
    }
    const std::size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + n, v.end());
    double med = v[n];
    if (v.size() % 2 == 0)
    {
      med = 0.5 * (med + *std::max_element(v.begin(), v.begin() + n));
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
Tpc_PolyTrackVertexFinder::Tpc_PolyTrackVertexFinder(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________
int Tpc_PolyTrackVertexFinder::InitRun(PHCompositeNode* topNode)
{
  if (!createNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
bool Tpc_PolyTrackVertexFinder::getNodes(PHCompositeNode* topNode)
{
  m_polyTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_inputNodeName);
  if (!m_polyTracks)
  {
    std::cerr << Name() << "::getNodes - missing Tpc_PolyTrackContainer node "
              << m_inputNodeName << std::endl;
    return false;
  }

  if (m_fillPolyVertexContainer)
  {
    m_polyVertices = findNode::getClass<Tpc_PolyTrackVertexContainer>(topNode, m_polyVertexNodeName);
    if (!m_polyVertices)
    {
      std::cerr << Name() << "::getNodes - missing " << m_polyVertexNodeName
                << " (run Tpc_PolyTrackVertexer first); collision fields will not be filled"
                << std::endl;
    }
  }
  return true;
}

//____________________________________________________________________________
bool Tpc_PolyTrackVertexFinder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  PHNodeIterator dstIter(dstNode);
  auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);
  if (!m_vertexMap)
  {
    m_vertexMap = new SvtxVertexMap_v1;
    auto* node = new PHIODataNode<PHObject>(m_vertexMap, m_vertexMapName, "PHObject");
    svtxNode->addNode(node);
    if (Verbosity() > 0)
    {
      std::cout << Name() << "::createNodes - created " << m_vertexMapName << std::endl;
    }
  }
  return true;
}

//____________________________________________________________________________
bool Tpc_PolyTrackVertexFinder::buildModel(const Tpc_PolyTrack* trk, unsigned int index,
                                          TrackModel& model) const
{
  if (!trk || trk->get_fit_status() <= 0)
  {
    return false;
  }
  if (trk->get_nclusters() < m_minClusters)
  {
    return false;
  }
  if (trk->get_ndf() > 0.0 && trk->get_chi2() / trk->get_ndf() > m_maxChi2Ndf)
  {
    return false;
  }

  const double x = trk->get_x();
  const double y = trk->get_y();
  const double z = trk->get_z();
  const double px = trk->get_px();
  const double py = trk->get_py();
  const double pz = trk->get_pz();
  const double q = trk->get_charge();
  if (!good(x) || !good(y) || !good(z) || !good(px) || !good(py) || !good(pz) || !good(q))
  {
    return false;
  }

  const double pt = std::hypot(px, py);
  if (pt <= 1.0e-9)
  {
    return false;
  }

  model.index = index;
  model.track_id = trk->get_track_id();
  model.nclusters = trk->get_nclusters();
  model.x0 = x;
  model.y0 = y;
  model.z0 = z;
  model.phi0 = std::atan2(py, px);
  model.tanl = pz / pt;
  model.is_line = m_zeroField || std::fabs(q) < 0.5 || std::fabs(m_magneticFieldTesla) < 1.0e-9;

  if (model.is_line)
  {
    // Line3D fits store a point and a direction in (x,y,z),(px,py,pz)
    model.dir = Eigen::Vector3d(px, py, pz).normalized();
    model.p = 0.0;
  }
  else
  {
    if (pt < m_minPt)
    {
      return false;
    }
    model.p = pt * std::sqrt(1.0 + model.tanl * model.tanl);

    const double r_mom = pt / (0.003 * std::fabs(m_magneticFieldTesla));
    const double sphi = std::sin(model.phi0);
    const double cphi = std::cos(model.phi0);

    // preferred: circle centre stored by Tpc_PolyTrackReco, if consistent
    bool use_stored = false;
    const double cx = trk->get_helix_x0();
    const double cy = trk->get_helix_y0();
    if (good(cx) && good(cy) && (std::fabs(cx) + std::fabs(cy)) > 0.0)
    {
      const double dx = cx - x;
      const double dy = cy - y;
      const double rg = std::hypot(dx, dy);
      const double along = dx * cphi + dy * sphi;  // should vanish: centre is normal to tangent
      if (rg > 0.0 && std::fabs(along) < 0.05 * rg && std::fabs(rg - r_mom) < 0.2 * r_mom)
      {
        use_stored = true;
        model.xc = cx;
        model.yc = cy;
        model.radius = rg;
        // centre - P0 = R h (-sin phi0, cos phi0)
        model.h = (-dx * sphi + dy * cphi) >= 0.0 ? 1.0 : -1.0;
      }
    }

    if (!use_stored)
    {
      // positive charge in +Bz bends clockwise
      model.radius = r_mom;
      model.h = (q * m_magneticFieldTesla > 0.0) ? -1.0 : 1.0;
      model.xc = x - model.radius * model.h * sphi;
      model.yc = y + model.radius * model.h * cphi;
    }

    if (!good(model.radius) || model.radius <= 0.0)
    {
      return false;
    }
  }

  // optional measured covariance (x,y,z block of the 6x6)
  if (trk->get_cov(0, 0) > 0.0 && trk->get_cov(1, 1) > 0.0 && trk->get_cov(2, 2) > 0.0)
  {
    model.has_cov = true;
    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        model.cov(i, j) = trk->get_cov(i, j);
      }
    }
  }

  // loose consistency with the beam spot
  const LinearTrack lin = linearize(model, Eigen::Vector3d(m_beamX, m_beamY, 0.0));
  if (!lin.ok)
  {
    return false;
  }
  const double dcaxy = std::hypot(lin.point.x() - m_beamX, lin.point.y() - m_beamY);
  if (dcaxy > m_maxDcaXY || std::fabs(lin.point.z()) > m_maxAbsZ0)
  {
    return false;
  }

  return true;
}

//____________________________________________________________________________
// Point of the track closest to `target` in the transverse plane, and the
// unit tangent there. Residuals are projected perpendicular to the tangent,
// so the remaining along-track offset does not matter.
Tpc_PolyTrackVertexFinder::LinearTrack
Tpc_PolyTrackVertexFinder::linearize(const TrackModel& m, const Eigen::Vector3d& target) const
{
  LinearTrack lin;

  if (m.is_line)
  {
    const Eigen::Vector3d p0(m.x0, m.y0, m.z0);
    const Eigen::Vector3d diff = target - p0;
    const double dxy2 = m.dir.x() * m.dir.x() + m.dir.y() * m.dir.y();
    double s = 0.0;
    if (dxy2 > 1.0e-12)
    {
      s = (diff.x() * m.dir.x() + diff.y() * m.dir.y()) / dxy2;
    }
    else
    {
      s = diff.dot(m.dir);
    }
    lin.point = p0 + s * m.dir;
    lin.dir = m.dir;
    lin.ok = lin.point.allFinite();
    return lin;
  }

  // helix: angle of the target around the circle centre
  const double psi = std::atan2(target.y() - m.yc, target.x() - m.xc);
  const double phi_t = psi + m.h * 0.5 * kPi;              // tangent angle there
  const double dphi = wrap_phi(phi_t - m.phi0);
  const double s = m.h * m.radius * dphi;                  // transverse arc length

  lin.point = Eigen::Vector3d(m.xc + m.radius * std::cos(psi),
                              m.yc + m.radius * std::sin(psi),
                              m.z0 + s * m.tanl);
  lin.dir = Eigen::Vector3d(std::cos(phi_t), std::sin(phi_t), m.tanl).normalized();
  lin.ok = lin.point.allFinite() && lin.dir.allFinite();
  return lin;
}

//____________________________________________________________________________
// 3x3 weight matrix W = B C^-1 B^T, B = [u w] spans the plane normal to the
// track; u is transverse, w = dir x u. C is the 2x2 position covariance of
// the track at the PCA in that plane.
Eigen::Matrix3d Tpc_PolyTrackVertexFinder::weightMatrix(const TrackModel& m,
                                                        const LinearTrack& lin) const
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
  if (m_useTrackCovariance && m.has_cov)
  {
    C = B.transpose() * m.cov * B;
  }
  else
  {
    // position error: sigma_rphi along u, sigma_z along z projected onto w
    C.setZero();
    C(0, 0) = m_sigmaRPhi * m_sigmaRPhi;
    C(1, 1) = m_sigmaZ * m_sigmaZ * w.z() * w.z();
  }

  if (m.p > 0.0 && m_msTerm > 0.0)
  {
    const double ms = m_msTerm / m.p;
    C(0, 0) += ms * ms;
    C(1, 1) += ms * ms;
  }
  C(0, 0) = std::max(C(0, 0), 1.0e-8);
  C(1, 1) = std::max(C(1, 1), 1.0e-8);

  return B * C.inverse() * B.transpose();
}

//____________________________________________________________________________
double Tpc_PolyTrackVertexFinder::dcaTwoLines(const LinearTrack& t1, const LinearTrack& t2,
                                             Eigen::Vector3d& pca1, Eigen::Vector3d& pca2)
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
std::vector<Tpc_PolyTrackVertexFinder::TrackPair>
Tpc_PolyTrackVertexFinder::findPairs(const std::vector<LinearTrack>& lin, double dcacut) const
{
  std::vector<TrackPair> pairs;
  const unsigned int n = lin.size();
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

      // refine: re-linearise both helices at the pair midpoint
      for (int k = 0; k < 2; ++k)
      {
        const Eigen::Vector3d mid = 0.5 * (p1 + p2);
        const LinearTrack li = linearize(m_models[i], mid);
        const LinearTrack lj = linearize(m_models[j], mid);
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
      if (std::hypot(mid.x() - m_beamX, mid.y() - m_beamY) > m_maxDcaXY ||
          std::fabs(mid.z()) > m_maxAbsZ0)
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
std::vector<std::vector<unsigned int>>
Tpc_PolyTrackVertexFinder::connectedComponents(unsigned int nnodes,
                                               const std::vector<TrackPair>& pairs)
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
    comps[root_to_comp[r]].push_back(i);
  }
  return comps;
}

//____________________________________________________________________________
// For every candidate, keep the pairs whose midpoint lies within the outlier
// cut of the candidate median. Rejected pairs are re-examined (they may
// belong to a second vertex merged into the same candidate), but tracks
// already claimed by a kept pair are not reused.
std::vector<Tpc_PolyTrackVertexFinder::TrackPair>
Tpc_PolyTrackVertexFinder::removeOutlierPairs(const std::vector<TrackPair>& pairs) const
{
  std::vector<TrackPair> result;
  std::vector<TrackPair> remaining = pairs;
  std::set<unsigned int> claimed;
  const unsigned int n = m_models.size();

  for (int pass = 0; pass < 5 && !remaining.empty(); ++pass)
  {
    // drop pairs touching tracks claimed in earlier passes
    if (pass > 0)
    {
      remaining.erase(std::remove_if(remaining.begin(), remaining.end(),
                                     [&claimed](const TrackPair& p)
                                     { return claimed.count(p.a) || claimed.count(p.b); }),
                      remaining.end());
    }

    DisjointSet ds(n);
    for (const auto& p : remaining)
    {
      ds.unite(p.a, p.b);
    }

    // median midpoint per component
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
        if ((p->midpoint() - med).norm() < m_outlierPairCut)
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
double Tpc_PolyTrackVertexFinder::trackChi2(const TrackModel& model, const Eigen::Vector3d& vtx) const
{
  const LinearTrack lin = linearize(model, vtx);
  if (!lin.ok)
  {
    return 1.0e30;
  }
  const Eigen::Vector3d r = vtx - lin.point;
  return r.dot(weightMatrix(model, lin) * r);
}

//____________________________________________________________________________
Tpc_PolyTrackVertexFinder::VertexFit
Tpc_PolyTrackVertexFinder::fitVertex(std::vector<unsigned int> tracks,
                                     const Eigen::Vector3d& seed) const
{
  VertexFit fit;
  Eigen::Vector3d vtx = seed;

  // one Gauss-Newton solve: tracks linearised at `at`
  auto solve = [&](const std::vector<unsigned int>& active, const Eigen::Vector3d& at,
                   Eigen::Vector3d& out, Eigen::Matrix3d& A) -> bool
  {
    A.setZero();
    Eigen::Vector3d b = Eigen::Vector3d::Zero();
    unsigned int nused = 0;
    for (const unsigned int it : active)
    {
      const TrackModel& m = m_models[it];
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
    if (m_useBeamSpotConstraint)
    {
      const double wx = 1.0 / (m_beamSigmaX * m_beamSigmaX);
      const double wy = 1.0 / (m_beamSigmaY * m_beamSigmaY);
      A(0, 0) += wx;
      A(1, 1) += wy;
      b(0) += wx * m_beamX;
      b(1) += wy * m_beamY;
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
    if (tracks.size() < m_minTracksPerVertex)
    {
      return fit;
    }

    // iterate to convergence with re-linearisation at the current vertex
    bool converged = false;
    for (unsigned int iter = 0; iter < m_maxIterations; ++iter)
    {
      Eigen::Vector3d next;
      if (!solve(tracks, vtx, next, A))
      {
        return fit;
      }
      const double shift = (next - vtx).norm();
      vtx = next;
      if (shift < m_convergenceTol)
      {
        converged = true;
        break;
      }
    }
    if (!converged && Verbosity() > 1)
    {
      std::cout << Name() << "::fitVertex - not converged after " << m_maxIterations
                << " iterations, ntracks " << tracks.size() << std::endl;
    }

    // drop the worst track while it is incompatible
    double worst_chi2 = -1.0;
    std::size_t worst = 0;
    for (std::size_t k = 0; k < tracks.size(); ++k)
    {
      const double c = trackChi2(m_models[tracks[k]], vtx);
      if (c > worst_chi2)
      {
        worst_chi2 = c;
        worst = k;
      }
    }
    if (worst_chi2 > m_maxTrackChi2 && tracks.size() > m_minTracksPerVertex)
    {
      tracks.erase(tracks.begin() + static_cast<std::ptrdiff_t>(worst));
      continue;
    }
    if (worst_chi2 > m_maxTrackChi2)
    {
      return fit;  // only minimum number of tracks left and still incompatible
    }
    break;
  }

  // final quantities at the converged vertex
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
    fit.chi2 += trackChi2(m_models[it], vtx);
  }
  if (m_useBeamSpotConstraint)
  {
    const double dx = (vtx.x() - m_beamX) / m_beamSigmaX;
    const double dy = (vtx.y() - m_beamY) / m_beamSigmaY;
    fit.chi2 += dx * dx + dy * dy;
  }
  fit.ndf = 2 * static_cast<int>(tracks.size()) - 3 + (m_useBeamSpotConstraint ? 2 : 0);
  fit.tracks = std::move(tracks);
  fit.ok = fit.position.allFinite() && fit.covariance.allFinite();
  return fit;
}

//____________________________________________________________________________
int Tpc_PolyTrackVertexFinder::process_event(PHCompositeNode* topNode)
{
  if (!m_polyTracks || !m_vertexMap)
  {
    if (!getNodes(topNode) || !createNodes(topNode))
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }
  ++m_nEvents;
  m_vertexMap->Reset();
  m_models.clear();

  // 1. track models
  const unsigned int ntracks = m_polyTracks->size();
  for (unsigned int i = 0; i < ntracks; ++i)
  {
    TrackModel model;
    if (buildModel(m_polyTracks->get_track(i), i, model))
    {
      m_models.push_back(model);
    }
  }

  std::vector<VertexFit> vertices;
  if (m_models.size() >= m_minTracksPerVertex)
  {
    // 2. linearise at the beam spot and find compatible pairs
    const Eigen::Vector3d beam(m_beamX, m_beamY, 0.0);
    std::vector<LinearTrack> lin;
    lin.reserve(m_models.size());
    for (const auto& m : m_models)
    {
      lin.push_back(linearize(m, beam));
    }

    std::vector<TrackPair> pairs = findPairs(lin, m_pairDcaCut);
    if (pairs.empty())
    {
      pairs = findPairs(lin, 3.0 * m_pairDcaCut);
    }

    // 3. candidates = connected components after outlier-pair removal
    const std::vector<TrackPair> clean = removeOutlierPairs(pairs);
    const auto comps = connectedComponents(m_models.size(), clean);

    if (Verbosity() > 1)
    {
      std::cout << Name() << " selected tracks " << m_models.size()
                << " pairs " << pairs.size() << " clean pairs " << clean.size()
                << " candidates " << comps.size() << std::endl;
    }

    // 4. fit candidates
    for (const auto& comp : comps)
    {
      if (comp.size() < m_minTracksPerVertex)
      {
        continue;
      }
      const std::set<unsigned int> members(comp.begin(), comp.end());
      std::vector<double> xs;
      std::vector<double> ys;
      std::vector<double> zs;
      for (const auto& p : clean)
      {
        if (members.count(p.a))
        {
          const Eigen::Vector3d mid = p.midpoint();
          xs.push_back(mid.x());
          ys.push_back(mid.y());
          zs.push_back(mid.z());
        }
      }
      const Eigen::Vector3d seed(median(xs), median(ys), median(zs));
      VertexFit fit = fitVertex(comp, seed);
      if (fit.ok)
      {
        vertices.push_back(std::move(fit));
      }
    }

    // 5. absorb unused tracks into their most compatible vertex, then refit
    if (m_absorbUnusedTracks && !vertices.empty())
    {
      std::vector<bool> used(m_models.size(), false);
      for (const auto& v : vertices)
      {
        for (const unsigned int it : v.tracks)
        {
          used[it] = true;
        }
      }
      std::vector<bool> changed(vertices.size(), false);
      for (unsigned int it = 0; it < m_models.size(); ++it)
      {
        if (used[it])
        {
          continue;
        }
        double best = m_maxTrackChi2;
        int ibest = -1;
        for (std::size_t iv = 0; iv < vertices.size(); ++iv)
        {
          const double c = trackChi2(m_models[it], vertices[iv].position);
          if (c < best)
          {
            best = c;
            ibest = static_cast<int>(iv);
          }
        }
        if (ibest >= 0)
        {
          vertices[ibest].tracks.push_back(it);
          changed[ibest] = true;
        }
      }
      for (std::size_t iv = 0; iv < vertices.size(); ++iv)
      {
        if (!changed[iv])
        {
          continue;
        }
        VertexFit refit = fitVertex(vertices[iv].tracks, vertices[iv].position);
        if (refit.ok)
        {
          vertices[iv] = std::move(refit);
        }
      }
    }

    // largest vertex first
    std::sort(vertices.begin(), vertices.end(),
              [](const VertexFit& a, const VertexFit& b)
              { return a.tracks.size() > b.tracks.size(); });
  }

  // 6. output
  writeVertices(vertices);

  if (!vertices.empty())
  {
    ++m_nEventsWithVertex;
    m_nVertices += vertices.size();
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - input tracks " << ntracks
              << " selected " << m_models.size()
              << " vertices " << vertices.size() << std::endl;
    if (Verbosity() > 1)
    {
      for (const auto& v : vertices)
      {
        std::cout << "   vertex (" << v.position.x() << ", " << v.position.y() << ", "
                  << v.position.z() << ") +- (" << std::sqrt(v.covariance(0, 0)) << ", "
                  << std::sqrt(v.covariance(1, 1)) << ", " << std::sqrt(v.covariance(2, 2))
                  << ") ntracks " << v.tracks.size() << " chi2/ndf " << v.chi2 << "/"
                  << v.ndf << std::endl;
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
void Tpc_PolyTrackVertexFinder::writeVertices(const std::vector<VertexFit>& vertices)
{
  for (const auto& v : vertices)
  {
    auto svtx = std::make_unique<SvtxVertex_v3>();
    svtx->set_x(v.position.x());
    svtx->set_y(v.position.y());
    svtx->set_z(v.position.z());
    svtx->set_t(0.0);
    svtx->set_chisq(v.chi2);
    svtx->set_ndof(static_cast<unsigned int>(std::max(0, v.ndf)));
    svtx->set_beam_crossing(m_beamCrossing);
    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        svtx->set_error(i, j, v.covariance(i, j));
      }
    }
    for (const unsigned int it : v.tracks)
    {
      svtx->insert_track(m_models[it].track_id);
    }
    m_vertexMap->insert(svtx.release());  // the map assigns the vertex id
  }

  if (m_fillPolyVertexContainer && m_polyVertices)
  {
    m_polyVertices->clear_collision_vertices();
    for (const auto& v : vertices)
    {
      m_polyVertices->add_collision_vertex(v.position.x(), v.position.y(), v.position.z(),
                                           std::sqrt(std::max(0.0, v.covariance(2, 2))),
                                           static_cast<unsigned int>(v.tracks.size()));
    }
    m_polyVertices->set_collision_min_clusters(m_minClusters);
    m_polyVertices->set_collision_vertex_valid(vertices.empty() ? 0 : 1);
  }
}

//____________________________________________________________________________
int Tpc_PolyTrackVertexFinder::End(PHCompositeNode* /*topNode*/)
{
  std::cout << Name() << "::End - events " << m_nEvents
            << " with vertex " << m_nEventsWithVertex
            << " total vertices " << m_nVertices << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
