// PHSimpleGeomProp.cc

#include "PHSimpleGeomProp.h"

#include "PHGhostRejection.h"

#include <ffamodules/CDBInterface.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <phfield/PHFieldUtility.h>
#include <phfield/PHFieldConfig.h>

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TSystem.h>

#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <syncstream>

#include <omp.h>

namespace
{
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }

  struct CircleFitResult
  {
    double xc = 0.;
    double yc = 0.;
    double R = 0.;
    bool ok = false;
  };

  // ------------------------------------------------------------
  // Simple 3x3 solver (Gaussian elimination)
  // ------------------------------------------------------------
  bool solve3x3(const double A[3][3], const double b[3], double x[3])
  {
    double M[3][3];
    double rhs[3];
    for (int i = 0; i < 3; ++i)
    {
      rhs[i] = b[i];
      for (int j = 0; j < 3; ++j) M[i][j] = A[i][j];
    }

    // Forward elimination
    for (int k = 0; k < 3; ++k)
    {
      // pivot
      int piv = k;
      double maxAbs = std::fabs(M[k][k]);
      for (int i = k + 1; i < 3; ++i)
      {
        double val = std::fabs(M[i][k]);
        if (val > maxAbs)
        {
          maxAbs = val;
          piv = i;
        }
      }

      if (maxAbs < 1e-20) return false;

      if (piv != k)
      {
        for (int j = 0; j < 3; ++j) std::swap(M[k][j], M[piv][j]);
        std::swap(rhs[k], rhs[piv]);
      }

      const double diag = M[k][k];
      for (int i = k + 1; i < 3; ++i)
      {
        const double f = M[i][k] / diag;
        rhs[i] -= f * rhs[k];
        for (int j = k; j < 3; ++j) M[i][j] -= f * M[k][j];
      }
    }

    // Back substitution
    for (int i = 2; i >= 0; --i)
    {
      double sum = rhs[i];
      for (int j = i + 1; j < 3; ++j) sum -= M[i][j] * x[j];
      if (std::fabs(M[i][i]) < 1e-20) return false;
      x[i] = sum / M[i][i];
    }
    return true;
  }

  // ------------------------------------------------------------
  // Sagitta model in rotated frame: f(u') = sagitta at u'
  // ------------------------------------------------------------
  void sagittaModelCpp(
      double S, double x0, double invR,
      const std::vector<double>& up,
      std::vector<double>& f_out)
  {
    const size_t n = up.size();
    f_out.resize(n);

    for (size_t i = 0; i < n; ++i)
    {
      const double dx  = up[i] - x0;
      const double dx2 = dx*dx;
      const double dx4 = dx2*dx2;
      const double dx6 = dx4*dx2;

      const double invR3 = invR*invR*invR;
      const double invR5 = invR3*invR*invR;

      f_out[i] =
        S
        - 0.5   * invR  * dx2
        - 0.125 * invR3 * dx4
        - 0.0625* invR5 * dx6;
    }
  }

  // ------------------------------------------------------------
  // Sagitta-based circle fit in XY
  // ------------------------------------------------------------
  bool fitCircleSagittaXY(
    const std::vector<std::pair<double,double>>& pts,
    double& xc_out, double& yc_out, double& R_out)
  {
    const size_t n = pts.size();
    if (n < 3) return false;

    // 1) straight line y = m x + b
    double Sx = 0., Sy = 0., Sxx = 0., Sxy = 0.;
    for (size_t i = 0; i < n; ++i)
    {
      const double x = pts[i].first;
      const double y = pts[i].second;
      Sx  += x;
      Sy  += y;
      Sxx += x*x;
      Sxy += x*y;
    }

    const double denom = n*Sxx - Sx*Sx;
    double m = 0., b = 0.;
    if (std::fabs(denom) > 1e-20)
    {
      m = (n*Sxy - Sx*Sy) / denom;
      b = (Sy - m*Sx) / n;
    }
    else
    {
      m = 0.;
      b = (n > 0) ? (Sy / n) : 0.;
    }

    const double theta = std::atan(m);
    const double cth   = std::cos(theta);
    const double sth   = std::sin(theta);

    // 2) rotate to (u', v') with chord horizontal
    std::vector<double> up(n), vp(n);
    for (size_t i = 0; i < n; ++i)
    {
      const double x = pts[i].first;
      const double y = pts[i].second;

      const double x_shift = x;
      const double y_shift = y - b; // line passes through origin

      // xr = cos θ * x + sin θ * (y-b)
      // yr = -sin θ * x + cos θ * (y-b)
      up[i] =  cth * x_shift + sth * y_shift;
      vp[i] = -sth * x_shift + cth * y_shift;
    }

    // 3) initial guesses: use max v' as sagitta apex
    double S0   = vp[0];
    double x0_0 = up[0];
    for (size_t i = 1; i < n; ++i)
    {
      if (vp[i] > S0)
      {
        S0   = vp[i];
        x0_0 = up[i];
      }
    }

    double S    = S0;
    double x0   = x0_0;
    double invR = 1.0 / 100.0; // R ~ 100 cm

    // 4) Gauss–Newton with numerical Jacobian
    const int    maxIter   = 50;
    const double tolDelta  = 1e-8;
    const double eps_param = 1e-6;

    std::vector<double> f(n), r(n);
    int nfev = 0;

    for (int iter = 0; iter < maxIter; ++iter)
    {
      // model + residuals
      sagittaModelCpp(S,   x0,   invR, up, f);
      nfev += 1;

      double cost = 0.0;
      for (size_t i = 0; i < n; ++i)
      {
        r[i] = f[i] - vp[i];
        cost += 0.5 * r[i] * r[i];
      }

      double JTJ[3][3] = { {0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.} };
      double JTr[3]    = { 0., 0., 0. };

      // perturb each parameter
      std::vector<double> fS_plus(n), fS_minus(n);
      std::vector<double> fx0_plus(n), fx0_minus(n);
      std::vector<double> fR_plus(n),  fR_minus(n);

      sagittaModelCpp(S + eps_param, x0, invR, up, fS_plus);
      sagittaModelCpp(S - eps_param, x0, invR, up, fS_minus);
      sagittaModelCpp(S, x0 + eps_param, invR, up, fx0_plus);
      sagittaModelCpp(S, x0 - eps_param, invR, up, fx0_minus);
      sagittaModelCpp(S, x0, invR + eps_param, up, fR_plus);
      sagittaModelCpp(S, x0, invR - eps_param, up, fR_minus);
      nfev += 6;

      for (size_t i = 0; i < n; ++i)
      {
        const double dfdS    = (fS_plus[i]  - fS_minus[i])  / (2.0 * eps_param);
        const double dfdx0   = (fx0_plus[i] - fx0_minus[i]) / (2.0 * eps_param);
        const double dfdinvR = (fR_plus[i]  - fR_minus[i])  / (2.0 * eps_param);

        const double J[3] = { dfdS, dfdx0, dfdinvR };
        const double ri   = r[i];

        for (int a = 0; a < 3; ++a)
        {
          JTr[a] += J[a] * ri;
          for (int j = 0; j < 3; ++j)
          {
            JTJ[a][j] += J[a] * J[j];
          }
        }
      }

      for (int a = 0; a < 3; ++a) JTr[a] = -JTr[a];

      double delta[3] = {0.,0.,0.};
      if (!solve3x3(JTJ, JTr, delta)) break;

      const double alpha = 1.0;
      S    += alpha * delta[0];
      x0   += alpha * delta[1];
      invR += alpha * delta[2];

      if (!std::isfinite(S) || !std::isfinite(x0) || !std::isfinite(invR))
        break;

      if (invR == 0.) invR = 1.0 / 100.0;

      const double maxDelta =
        std::max(std::fabs(delta[0]),
        std::max(std::fabs(delta[1]), std::fabs(delta[2])));

      if (maxDelta < tolDelta) break;
    }

    const double R = (invR != 0.0) ? (1.0 / invR) : 1e9;
    if (!std::isfinite(R) || R == 0.) return false;
    
    const double Xc_prime = x0;
    const double Yc_prime = S - R; // signed curvature

    const double xc_shift =  cth * Xc_prime - sth * Yc_prime;
    const double yc_shift =  sth * Xc_prime + cth * Yc_prime;

    const double xc = xc_shift;
    const double yc = yc_shift + b;  // undo shift in y

    if (!std::isfinite(xc) || !std::isfinite(yc)) return false;

    xc_out = xc;
    yc_out = yc;
    R_out =R;
    return true;
  }

  // ------------------------------------------------------------
  // Wrapper: sagitta-based circle fit returning CircleFitResult
  // ------------------------------------------------------------
  /*CircleFitResult sagitta_circle_fit(
      const std::vector<double>& xs,
      const std::vector<double>& ys)
  {
    CircleFitResult res;
    const size_t n = xs.size();
    if (n < 3 || ys.size() != n) return res;

    std::vector<std::pair<double,double>> pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
      pts.emplace_back(xs[i], ys[i]);
    }

    double xc = 0., yc = 0., R = 0.;
    if (!fitCircleSagittaXY(pts, xc, yc, R)) return res;

    res.xc = xc;
    res.yc = yc;
    res.R  = R;
    res.ok = true;
    return res;
  }*/

} // namespace



/*
  CircleFitResult weighted_circle_fit(
      const std::vector<double> &xs,
      const std::vector<double> &ys,
      const std::vector<double> &w)
  {
    CircleFitResult res;

    const size_t n = xs.size();
    if (n < 3 || ys.size() != n || w.size() != n)
    {
      return res;
    }

    double Sw = 0.;
    double Sx = 0., Sy = 0.;
    double Sxx = 0., Syy = 0., Sxy = 0.;
    double Sx3 = 0., Sy3 = 0.;
    double Sx1y2 = 0., Sx2y1 = 0.;

    for (size_t i = 0; i < n; ++i)
    {
      const double wi = w[i];
      if (wi <= 0.)
      {
        continue;
      }

      const double x = xs[i];
      const double y = ys[i];
      const double x2 = x * x;
      const double y2 = y * y;

      Sw += wi;
      Sx += wi * x;
      Sy += wi * y;
      Sxx += wi * x2;
      Syy += wi * y2;
      Sxy += wi * x * y;
      Sx3 += wi * x2 * x;
      Sy3 += wi * y2 * y;
      Sx1y2 += wi * x * y2;
      Sx2y1 += wi * x2 * y;
    }

    if (Sw <= 0.)
    {
      return res;
    }

    Eigen::Matrix3d A;
    Eigen::Vector3d b;

    A << Sxx, Sxy, Sx,
        Sxy, Syy, Sy,
        Sx, Sy, Sw;

    b << -(Sx3 + Sx1y2),
        -(Sx2y1 + Sy3),
        -(Sxx + Syy);

    const double det = A.determinant();
    if (std::abs(det) < 1e-12)
    {
      return res;
    }

    Eigen::Vector3d p = A.fullPivLu().solve(b);
    const double a = p(0);
    const double c = p(1);
    const double d = p(2);

    const double xc = -a / 2.;
    const double yc = -c / 2.;
    const double R2 = xc * xc + yc * yc - d;
    if (R2 <= 0.)
    {
      return res;
    }

    res.xc = xc;
    res.yc = yc;
    res.R = std::sqrt(R2);
    res.ok = std::isfinite(res.R);
    return res;
  }

  // ------------------------------------------------------------------
  // Robust circle fit with Tukey weight function (Eq. 3.8 in the figs)
  //
  //   w(Δ/RMS) = 1 - ( (Δ/RMS) / C_T )^n  for |Δ|/RMS < C_T
  //           = 0                         otherwise
  //
  // Δ is the radial residual: (sqrt((x-xc)^2 + (y-yc)^2) - R)
  // ------------------------------------------------------------------
  CircleFitResult robust_tukey_circle_fit(
      const std::vector<double> &xs,
      const std::vector<double> &ys,
      double C_T = 3.0,
      double n = 1.5,
      int max_iter = 10)
  {
    const size_t N = xs.size();
    CircleFitResult result;

    if (N < 3 || ys.size() != N)
    {
      return result;
    }

    // start with equal weights
    std::vector<double> w(N, 1.0);

    CircleFitResult current = weighted_circle_fit(xs, ys, w);
    if (!current.ok)
    {
      return current;
    }

    for (int it = 0; it < max_iter; ++it)
    {
      // compute residuals
      double sumw = 0.;
      double sumw_res2 = 0.;
      std::vector<double> res(N, 0.);

      for (size_t i = 0; i < N; ++i)
      {
        const double dx = xs[i] - current.xc;
        const double dy = ys[i] - current.yc;
        const double r = std::sqrt(dx * dx + dy * dy);
        const double delta = r - current.R;

        res[i] = delta;
        const double wi = w[i];
        sumw += wi;
        sumw_res2 += wi * delta * delta;
      }

      if (sumw <= 0.)
      {
        break;
      }

      const double rms = std::sqrt(sumw_res2 / std::max<double>(1., sumw - 3.));
      if (!std::isfinite(rms) || rms <= 0.)
      {
        break;
      }

      // update weights with Tukey function
      bool changed = false;
      for (size_t i = 0; i < N; ++i)
      {
        const double u = std::abs(res[i]) / rms;
        double wi_old = w[i];
        double wi_new = 0.;

        if (u < C_T)
        {
          const double ratio = u / C_T;
          wi_new = 1.0 - std::pow(ratio, n);
          if (wi_new < 0.)
          {
            wi_new = 0.;
          }
        }
        else
        {
          wi_new = 0.;
        }

        if (std::abs(wi_new - wi_old) > 1e-4)
        {
          changed = true;
        }
        w[i] = wi_new;
      }

      if (!changed)
      {
        break;
      }

      CircleFitResult next = weighted_circle_fit(xs, ys, w);
      if (!next.ok)
      {
        break;
      }
      current = next;
    }

    return current;
  }

}  // namespace*/

// ------------------------------------------------------------------
// ctor / dtor / End
// ------------------------------------------------------------------

PHSimpleGeomProp::PHSimpleGeomProp(const std::string &name)
    : SubsysReco(name)
{
}

int PHSimpleGeomProp::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

// ------------------------------------------------------------------
// InitRun: grab nodes, set threads
// ------------------------------------------------------------------

int PHSimpleGeomProp::InitRun(PHCompositeNode *topNode)
{
  int ret = get_nodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  std::cout << "PHSimpleGeomProp::InitRun - m_num_threads: " << m_num_threads << std::endl;
  if (m_num_threads >= 1)
  {
    omp_set_num_threads(m_num_threads);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

// ------------------------------------------------------------------
// get_nodes: ActsGeometry, TPC geometry, clusters, seeds
// ------------------------------------------------------------------

int PHSimpleGeomProp::get_nodes(PHCompositeNode *topNode)
{
  // Acts geometry
  m_tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tgeometry)
  {
    std::cout << PHWHERE << " No Acts tracking geometry, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc global position wrapper (for distortions)
  m_globalPositionWrapper.loadNodes(topNode);

  // clusters
  if (_use_truth_clusters)
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_TRUTH");
  }
  else
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  if (!_cluster_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // TPC seeds
  _track_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_track_map)
  {
    std::cerr << PHWHERE << " ERROR: Can't find TrackSeedContainer TpcTrackSeedContainer" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // TPC geometry
  auto geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!geom_container)
  {
    std::cerr << PHWHERE << "ERROR: Can't find node TPCGEOMCONTAINER" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  radii.clear();
  for (int i = 7; i <= 54; ++i)
  {
    radii.push_back(geom_container->GetLayerCellGeom(i)->get_radius());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

// ------------------------------------------------------------------
// main event function
// ------------------------------------------------------------------

int PHSimpleGeomProp::process_event(PHCompositeNode *topNode)
{
  if (_n_iteration != 0)
  {
    _iteration_map =
        findNode::getClass<TrkrClusterIterationMapv1>(topNode, "CLUSTER_ITERATION_MAP");
    if (!_iteration_map)
    {
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  PHTimer timer("GeomPropTimer");
  timer.restart();

  if (_max_seeds > 0 && _track_map->size() > _max_seeds)
  {
    std::cout << "number of TPC seeds: " << _track_map->size() << std::endl;
    std::cout << PHWHERE << "number of TPC seeds > " << _max_seeds << " aborting event."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const auto globalPositions = PrepareKDTrees();
  if (Verbosity())
  {
    std::cout << "PHSimpleGeomProp::process_event - PrepareKDTrees time: "
              << timer.elapsed() << " ms" << std::endl;
  }

  std::vector<std::vector<TrkrDefs::cluskey>> new_chains;
  std::vector<TrackSeed_v2> unused_tracks;

  timer.restart();
#pragma omp parallel
  {
    if (Verbosity())
    {
      std::osyncstream(std::cout)
          << "PHSimpleGeomProp -"
          << " num_threads: " << omp_get_num_threads()
          << " this thread: " << omp_get_thread_num()
          << std::endl;
    }

    std::vector<std::vector<TrkrDefs::cluskey>> local_chains;
    std::vector<TrackSeed_v2> local_unused;

#pragma omp for schedule(static)
    for (size_t track_it = 0; track_it < _track_map->size(); ++track_it)
    {
      auto track = _track_map->get(track_it);
      if (!track)
      {
        continue;
      }

      // keep non-TPC seeds unchanged
      const bool is_tpc = std::any_of(
          track->begin_cluster_keys(), track->end_cluster_keys(),
          [](const TrkrDefs::cluskey &key)
          { return TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId; });

      if (!is_tpc)
      {
        local_unused.emplace_back(*track);
        continue;
      }

      // copy seed cluster keys
      std::vector<TrkrDefs::cluskey> seed_keys;
      std::copy(track->begin_cluster_keys(), track->end_cluster_keys(),
                std::back_inserter(seed_keys));

      if (seed_keys.size() < 3)
      {
        continue;
      }

      // propagate both directions using robust circle geometry
      auto chain_in = PropagateTrackSimple(track, seed_keys,
                                           PropagationDirection::Inward,
                                           globalPositions);
      auto chain_out = PropagateTrackSimple(track, seed_keys,
                                            PropagationDirection::Outward,
                                            globalPositions);

      std::vector<TrkrDefs::cluskey> full_chain;
      if (chain_in.size() >= chain_out.size())
      {
        full_chain = std::move(chain_in);
      }
      else
      {
        full_chain = std::move(chain_out);
      }

      if (full_chain.size() < 3)
      {
        continue;
      }

      std::sort(full_chain.begin(), full_chain.end());
      full_chain.erase(std::unique(full_chain.begin(), full_chain.end()),
                       full_chain.end());

      local_chains.push_back(std::move(full_chain));
    }  // seed loop

    std::sort(local_chains.begin(), local_chains.end());
    local_chains.erase(std::unique(local_chains.begin(), local_chains.end()),
                       local_chains.end());

#pragma omp critical
    {
      new_chains.reserve(new_chains.size() + local_chains.size());
      new_chains.insert(new_chains.end(),
                        std::make_move_iterator(local_chains.begin()),
                        std::make_move_iterator(local_chains.end()));

      unused_tracks.reserve(unused_tracks.size() + local_unused.size());
      unused_tracks.insert(unused_tracks.end(),
                           std::make_move_iterator(local_unused.begin()),
                           std::make_move_iterator(local_unused.end()));
    }
  }  // parallel

  if (Verbosity())
  {
    std::cout << "PHSimpleGeomProp::process_event - first seed loop time: "
              << timer.elapsed() << " ms" << std::endl;
  }

  // global deduplication & simple cleaning
  timer.restart();
  std::sort(new_chains.begin(), new_chains.end());
  new_chains.erase(std::unique(new_chains.begin(), new_chains.end()),
                   new_chains.end());

  new_chains.erase(
      std::remove_if(new_chains.begin(), new_chains.end(),
                     [](const auto &c)
                     { return c.size() < 3; }),
      new_chains.end());

  if (Verbosity())
  {
    std::cout << "PHSimpleGeomProp::process_event - cleaned chains: "
              << new_chains.size() << std::endl;
  }

  // build seeds with circle+line fits (no Kalman)
  std::vector<TrackSeed_v2> new_seeds;
  new_seeds.reserve(new_chains.size());
  for (const auto &chain : new_chains)
  {
    new_seeds.push_back(makeSeedFromChain(chain, globalPositions));
  }

  // dummy chi2 for ghost rejection (you can compute from residuals if needed)
  std::vector<float> trackChi2(new_seeds.size(), 1.0);

  // reset track map and publish
  _track_map->Reset();

  if (m_ghostrejection)
  {
    rejectAndPublishSeeds(new_seeds, globalPositions, trackChi2);
  }
  else
  {
    publishSeeds(new_seeds);
  }

  // also publish unused non-TPC seeds
  publishSeeds(unused_tracks);

  if (Verbosity())
  {
    std::cout << "PHSimpleGeomProp::process_event - total time: "
              << timer.elapsed() << " ms" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

// ------------------------------------------------------------------
// global position helper
// ------------------------------------------------------------------

Acts::Vector3 PHSimpleGeomProp::getGlobalPosition(
    TrkrDefs::cluskey key, TrkrCluster *cluster) const
{
  return _pp_mode
             ? m_tgeometry->getGlobalPosition(key, cluster)
             : m_globalPositionWrapper.getGlobalPositionDistortionCorrected(key, cluster, 0);
}

// ------------------------------------------------------------------
// PrepareKDTrees: identical logic to KF version
// ------------------------------------------------------------------

PHSimpleGeomProp::PositionMap PHSimpleGeomProp::PrepareKDTrees()
{
  PositionMap globalPositions;

  std::vector<std::vector<std::vector<double>>> kdhits;
  kdhits.resize(58);

  if (!_cluster_map)
  {
    std::cout << "WARNING: cluster map not provided" << std::endl;
    return globalPositions;
  }

  for (const auto &hitsetkey : _cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto range = _cluster_map->getClusters(hitsetkey);
    for (TrkrClusterContainer::ConstIterator it = range.first; it != range.second; ++it)
    {
      const auto &[cluskey, cluster] = *it;
      if (!cluster)
      {
        continue;
      }

      // skip hits used in a previous iteration
      if (_n_iteration && _iteration_map &&
          _iteration_map->getIteration(cluskey) > 0)
      {
        continue;
      }

      const auto gpos = getGlobalPosition(cluskey, cluster);
      globalPositions.emplace(cluskey, gpos);

      const int layer = TrkrDefs::getLayer(cluskey);
      std::vector<double> kdhit{gpos.x(), gpos.y(), gpos.z(), 0.};
      const uint64_t key = cluskey;
      std::memcpy(&kdhit[3], &key, sizeof(key));

      kdhits[layer].push_back(std::move(kdhit));
    }
  }

  _ptclouds.resize(kdhits.size());
  _kdtrees.resize(kdhits.size());

  for (size_t l = 0; l < kdhits.size(); ++l)
  {
    _ptclouds[l] = std::make_shared<KDPointCloud<double>>();
    _ptclouds[l]->pts = std::move(kdhits[l]);

    _kdtrees[l] = std::make_shared<KDTree_t>(
        3, *(_ptclouds[l]),
        nanoflann::KDTreeSingleIndexAdaptorParams(10));
    _kdtrees[l]->buildIndex();
  }

  return globalPositions;
}

// ------------------------------------------------------------------
// Geometry-only propagation using robust circle fit + line in (r,z)
// ------------------------------------------------------------------

std::vector<TrkrDefs::cluskey> PHSimpleGeomProp::PropagateTrackSimple(
    TrackSeed *track,
    const std::vector<TrkrDefs::cluskey> &seed_keys,
    PropagationDirection direction,
    const PositionMap &globalPositions) const
{
  std::vector<TrkrDefs::cluskey> result;

  if (!track || seed_keys.size() < 3)
  {
    return result;
  }

  // build xy and rz arrays
  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<std::pair<double, double> > rz_pts;
  std::vector<std::pair<double, double> > allPts;

  xs.reserve(seed_keys.size());
  ys.reserve(seed_keys.size());
  rz_pts.reserve(seed_keys.size());

  for (auto key : seed_keys)
  {
    const auto &pos = globalPositions.at(key);
    const double x = pos(0);
    const double y = pos(1);
    const double z = pos(2);
    const double r = std::sqrt(square(x) + square(y));

    xs.push_back(x);
    ys.push_back(y);
    allPts.push_back(std::make_pair(x,y));
    rz_pts.emplace_back(r, z);
    std::cout<<"PropagateTrackSimple: seed key " << key
             << " pos = ( " << x << "   " << y << "   " << z << " ) r=" << r
             << std::endl;
  }

  // robust circle fit (Tukey)
  //const CircleFitResult circle =
 //     robust_tukey_circle_fit(xs, ys, m_tukey_C, m_tukey_n, m_circle_max_iter);

  // sagitta-based geometric circle fit
 /* const CircleFitResult circle = sagitta_circle_fit(xs, ys);
  if (!circle.ok)
  {
    
    return seed_keys;
  }
  const double R = circle.R;
  const double X0 = circle.xc;
  const double Y0 = circle.yc;
  */
double R=0, X0=0, Y0=0;
double xc_tmp=0, yc_tmp=0, R_tmp=0;
    if (fitCircleSagittaXY(allPts, xc_tmp, yc_tmp, R_tmp))
    {
      X0 = xc_tmp;
      Y0 = yc_tmp;
      R  = R_tmp;
    }


  // line fit in (r,z): z = A + B r
  const auto [B, A] = TrackFitUtils::line_fit(rz_pts);

  if (Verbosity() > 1)
  {
    std::cout << "PropagateTrackSimple: circle xc=  " << X0
              << " yc=  " << Y0 << " R=  " << R
              << ", z(r)=  " << A << "  +  " << B << "  r" << std::endl;
  }

  // map layer -> seed cluskey
  std::map<unsigned int, TrkrDefs::cluskey> layer_to_key;
  for (auto key : seed_keys)
  {
    layer_to_key[TrkrDefs::getLayer(key)] = key;
  }

  // starting layer and position: innermost or outermost
  unsigned int current_layer = 0;
  TrkrDefs::cluskey current_key = 0;

  if (direction == PropagationDirection::Outward)
  {
    current_key = seed_keys.back();
  }
  else
  {
    current_key = seed_keys.front();
  }

  current_layer = TrkrDefs::getLayer(current_key);
  auto pos0 = globalPositions.at(current_key);
  double curr_x = pos0(0);
  double curr_y = pos0(1);
  double curr_z = pos0(2);

  if (Verbosity() > 1)
  {
    std::cout << "PropagateTrackSimple: start at layer " << current_layer
              << " pos = (" << curr_x << ", " << curr_y << ", " << curr_z
              << ")" << std::endl;
  }

  result = seed_keys;

  int max_missing_layers = 3;
  int missing_layers = 0;

  for (int step = 0; step <= _max_propagation_steps; ++step)
  {
    int next_layer = (direction == PropagationDirection::Outward)
                         ? int(current_layer) + 1
                         : int(current_layer) - 1;

    if (next_layer < 7 || next_layer > 54)
    {
      if (Verbosity() > 1)
      {
        std::cout << "PropagateTrackSimple: reached layer " << next_layer
                  << " (outside TPC), stop" << std::endl;
      }
      break;
    }

    const double target_radius = radii[next_layer - 7];

    // if we already have a seed hit in this layer, just move there
    auto it_layer = layer_to_key.find(next_layer);
    if (it_layer != layer_to_key.end())
    {
      current_layer = next_layer;
      current_key = it_layer->second;

      auto pos = globalPositions.at(current_key);
      curr_x = pos(0);
      curr_y = pos(1);
      curr_z = pos(2);

      if (std::find(result.begin(), result.end(), current_key) == result.end())
      {
        result.push_back(current_key);
      }

      missing_layers = 0;
      continue;
    }

    // project circle onto cylinder of radius target_radius
    auto inter = TrackFitUtils::circle_circle_intersection(
        target_radius, std::fabs(R), X0, Y0);

    const double xplus = std::get<0>(inter);
    const double yplus = std::get<1>(inter);
    const double xminus = std::get<2>(inter);
    const double yminus = std::get<3>(inter);

    // choose intersection closer to previous point
    const double dxp = xplus - curr_x;
    const double dyp = yplus - curr_y;
    const double dxm = xminus - curr_x;
    const double dym = yminus - curr_y;
    const double d2p = dxp * dxp + dyp * dyp;
    const double d2m = dxm * dxm + dym * dym;

    const double pred_x = (d2p < d2m) ? xplus : xminus;
    const double pred_y = (d2p < d2m) ? yplus : yminus;
    const double pred_r = std::sqrt(square(pred_x) + square(pred_y));
    const double pred_z = A + B * pred_r;

    if (Verbosity() > 2)
    {
      std::cout << "PropagateTrackSimple: step " << step
                << " layer " << next_layer
                << " predicted (" << pred_x << ", " << pred_y << ", " << pred_z
                << ")" << std::endl;
    }

    // KD nearest cluster
    double query_pt[3] = {pred_x, pred_y, pred_z};
    std::vector<long unsigned int> index_out(1);
    std::vector<double> distance_out(1);

    int n_results =
        _kdtrees[next_layer]->knnSearch(query_pt, 1, &index_out[0], &distance_out[0]);

    if (!n_results)
    {
      ++missing_layers;
      current_layer = next_layer;
      if (missing_layers > max_missing_layers)
      {
        break;
      }
      continue;
    }

    const std::vector<double> &point =
        _ptclouds[next_layer]->pts[index_out[0]];
    TrkrDefs::cluskey closest_ckey = (*((int64_t *)&point[3]));
    TrkrCluster *clusterCandidate = _cluster_map->findCluster(closest_ckey);
    if (!clusterCandidate)
    {
      ++missing_layers;
      current_layer = next_layer;
      if (missing_layers > max_missing_layers)
      {
        break;
      }
      continue;
    }

    const auto cand_pos = globalPositions.at(closest_ckey);
    const double cand_x = cand_pos(0);
    const double cand_y = cand_pos(1);
    const double cand_z = cand_pos(2);

    const double dx = cand_x - pred_x;
    const double dy = cand_y - pred_y;
    const double dz = cand_z - pred_z;
    const double dr = std::sqrt(dx * dx + dy * dy + dz * dz);

    if (dr > _max_dist)
    {
      ++missing_layers;
      current_layer = next_layer;
      if (missing_layers > max_missing_layers)
      {
        break;
      }
      continue;
    }

    // accept
    if (std::find(result.begin(), result.end(), closest_ckey) == result.end())
    {
      result.push_back(closest_ckey);
    }

    current_layer = next_layer;
    current_key = closest_ckey;
    curr_x = cand_x;
    curr_y = cand_y;
    curr_z = cand_z;
    missing_layers = 0;
  }

  if (Verbosity() > 1)
  {
    std::cout << "PropagateTrackSimple: final chain size " << result.size()
              << " layers: ";
    for (auto c : result)
    {
      std::cout << (int)TrkrDefs::getLayer(c) << " ";
    }
    std::cout << std::endl;
  }

  return result;
}

// ------------------------------------------------------------------
// Build TrackSeed_v2 from chain using circle+line fits
// ------------------------------------------------------------------

TrackSeed_v2 PHSimpleGeomProp::makeSeedFromChain(
    const std::vector<TrkrDefs::cluskey> &chain,
    const PositionMap &globalPositions) const
{
  TrackSeed_v2 seed;

  if (chain.size() < 3)
  {
    return seed;
  }

  for (auto key : chain)
  {
    seed.insert_cluster_key(key);
  }

  // build local position map
  PositionMap local;
  for (auto key : chain)
  {
    local.emplace(key, globalPositions.at(key));
  }

  // let TrackSeedHelper fill standard track params
  TrackSeedHelper::circleFitByTaubin(&seed, local, 7, 55);
  TrackSeedHelper::lineFit(&seed, local, 7, 55);
  seed.set_phi(TrackSeedHelper::get_phi(&seed, local));

  // DO NOT call set_charge (TrackSeed_v2 doesn't have it)
  // If you really want to enforce sign consistency, you can do:
  // int q = seed.get_charge();
  // if (q != 0) seed.set_qOverR(std::abs(seed.get_qOverR()) * q);

  return seed;
}



// ------------------------------------------------------------------
// RemoveBadClusters: keep only chains with >=3 hits (simple for now)
// ------------------------------------------------------------------

std::vector<std::vector<TrkrDefs::cluskey>>
PHSimpleGeomProp::RemoveBadClusters(
    const std::vector<std::vector<TrkrDefs::cluskey>> &chains,
    const PositionMap & /*globalPositions*/) const
{
  std::vector<std::vector<TrkrDefs::cluskey>> clean;
  std::copy_if(chains.begin(), chains.end(), std::back_inserter(clean),
               [](const auto &c)
               { return c.size() >= 3; });
  return clean;
}

// ------------------------------------------------------------------
// Ghost rejection + publishing (same strategy as KF module)
// ------------------------------------------------------------------

void PHSimpleGeomProp::rejectAndPublishSeeds(
    std::vector<TrackSeed_v2> &seeds,
    const PositionMap &positions,
    std::vector<float> &trackChi2)
{
  PHTimer timer("GeomPropGhostTimer");

  timer.restart();
  PHGhostRejection rejector(Verbosity(), seeds);

  // you can expose these as setters if you like
  rejector.set_phi_cut(0.02);
  rejector.set_eta_cut(0.05);
  rejector.set_x_cut(0.3);
  rejector.set_y_cut(0.3);
  rejector.set_z_cut(0.5);

  for (unsigned int itrack = 0; itrack < seeds.size(); ++itrack)
  {
    rejector.cut_from_clusters(itrack);
  }

#pragma omp parallel for schedule(static)
  for (unsigned int itrack = 0; itrack < seeds.size(); ++itrack)
  {
    if (rejector.is_rejected(itrack))
    {
      continue;
    }

    auto &seed = seeds[itrack];

    PositionMap local;
    std::transform(seed.begin_cluster_keys(), seed.end_cluster_keys(),
                   std::inserter(local, local.end()),
                   [&positions](const auto &key)
                   { return std::make_pair(key, positions.at(key)); });

    TrackSeedHelper::circleFitByTaubin(&seed, local, 7, 55);
    TrackSeedHelper::lineFit(&seed, local, 7, 55);
    seed.set_phi(TrackSeedHelper::get_phi(&seed, local));

    const int q = seed.get_charge();
    seed.set_qOverR(std::abs(seed.get_qOverR()) * q);
  }

  if (Verbosity())
  {
    std::cout << "PHSimpleGeomProp::rejectAndPublishSeeds - circle fit time: "
              << timer.elapsed() << " ms" << std::endl;
  }

  timer.restart();
  rejector.find_ghosts(trackChi2);
  if (Verbosity())
  {
    std::cout << "PHSimpleGeomProp::rejectAndPublishSeeds - ghost rejection: "
              << timer.elapsed() << " ms" << std::endl;
  }

  timer.restart();
  for (unsigned int itrack = 0; itrack < seeds.size(); ++itrack)
  {
    if (rejector.is_rejected(itrack))
    {
      if (Verbosity() > 0)
      {
        std::cout << " Seed " << ((int)itrack)
                  << " rejected. Not getting published." << std::endl;
      }
      continue;
    }

    auto &seed = seeds[itrack];
    _track_map->insert(&seed);

    if (Verbosity() > 0)
    {
      std::cout << "Publishing seed " << ((int)itrack)
                << " q " << seed.get_charge()
                << " qOverR " << seed.get_qOverR()
                << " x " << TrackSeedHelper::get_x(&seed)
                << " y " << TrackSeedHelper::get_y(&seed)
                << " z " << TrackSeedHelper::get_z(&seed)
                << " pT " << seed.get_pt()
                << " eta " << seed.get_eta()
                << " phi " << seed.get_phi()
                << std::endl;
    }
  }

  if (Verbosity())
  {
    std::cout << "PHSimpleGeomProp::rejectAndPublishSeeds - publication: "
              << timer.elapsed() << " ms" << std::endl;
  }
}

// ------------------------------------------------------------------
// publishSeeds: simple insertion
// ------------------------------------------------------------------

void PHSimpleGeomProp::publishSeeds(const std::vector<TrackSeed_v2> &seeds)
{
  for (const auto &seed : seeds)
  {
    _track_map->insert(&seed);
  }
}
