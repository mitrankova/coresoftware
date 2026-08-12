#include "PHGarfieldRossegger.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNamed.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace
{
  constexpr double epsilon0 = 8.8541878128e-12;
  constexpr double cmToM = 1.0e-2;
  constexpr double mToCm = 100.0;
  constexpr double elementaryCharge = 1.602176634e-19;

  unsigned int src_index(unsigned int ir, unsigned int ip, unsigned int iz, unsigned int np, unsigned int nz)
  { return (ir * np + ip) * nz + iz; }

  unsigned int fld_index(unsigned int ir, unsigned int ip, unsigned int iz, unsigned int np, unsigned int nz)
  { return (ir * np + ip) * nz + iz; }

  unsigned int mode_index(unsigned int in, unsigned int il, unsigned int nl)
  { return in * nl + il; }

  std::vector<double> edges(double lo, double hi, unsigned int n)
  {
    std::vector<double> out(n + 1);
    const double step = (hi - lo) / static_cast<double>(n);
    for (unsigned int i = 0; i <= n; ++i) { out[i] = lo + step * static_cast<double>(i); }
    return out;
  }

  std::vector<double> centers(const std::vector<double>& e)
  {
    std::vector<double> out(e.size() - 1);
    for (unsigned int i = 0; i + 1 < e.size(); ++i) { out[i] = 0.5 * (e[i] + e[i + 1]); }
    return out;
  }

  std::vector<double> phi_edges(unsigned int n)
  {
    std::vector<double> out(n + 1);
    const double step = 2.0 * std::numbers::pi / static_cast<double>(n);
    for (unsigned int i = 0; i <= n; ++i) { out[i] = static_cast<double>(i) * step; }
    return out;
  }

  std::vector<double> phi_centers(unsigned int n)
  {
    std::vector<double> out(n);
    const double step = 2.0 * std::numbers::pi / static_cast<double>(n);
    for (unsigned int i = 0; i < n; ++i) { out[i] = (static_cast<double>(i) + 0.5) * step; }
    return out;
  }

  double jprime(unsigned int m, double x)
  {
    return m == 0 ? -std::cyl_bessel_j(1, x) : 0.5 * (std::cyl_bessel_j(m - 1, x) - std::cyl_bessel_j(m + 1, x));
  }

  double yprime(unsigned int m, double x)
  {
    return m == 0 ? -std::cyl_neumann(1, x) : 0.5 * (std::cyl_neumann(m - 1, x) - std::cyl_neumann(m + 1, x));
  }

  void scale_to_cm(std::vector<double>& values)
  {
    for (double& value : values) { value *= mToCm; }
  }
}  // namespace

PHGarfieldRossegger::PHGarfieldRossegger(const std::string& name)
  : SubsysReco(name)
{
}

void PHGarfieldRossegger::setGeometryCm(double inner_radius, double outer_radius, double half_length)
{
  m_aCm = inner_radius; m_bCm = outer_radius; m_lCm = half_length;
}

void PHGarfieldRossegger::setSourceRadiusCm(double min_radius, double max_radius)
{
  m_sourceRMinCm = min_radius; m_sourceRMaxCm = max_radius;
}

void PHGarfieldRossegger::setDensity(double reference_density_nC_per_m3, double k_eff, double radial_power_alpha)
{
  m_rhoReferenceNCPerM3 = reference_density_nC_per_m3; m_kEff = k_eff; m_radialPowerAlpha = radial_power_alpha;
}

void PHGarfieldRossegger::setPhiModulation(double m1_amplitude, double m1_phase, double m12_amplitude, double m12_phase)
{
  m_m1Amplitude = m1_amplitude; m_m1Phase = m1_phase; m_m12Amplitude = m12_amplitude; m_m12Phase = m12_phase;
}

void PHGarfieldRossegger::setSourceGrid(unsigned int nr, unsigned int nphi, unsigned int nz)
{
  m_nrSource = nr; m_nphiSource = nphi; m_nzSource = nz;
}

void PHGarfieldRossegger::setObservationGrid(unsigned int nr, unsigned int nphi, unsigned int nz)
{
  m_nrObs = nr; m_nphiObs = nphi; m_nzObs = nz;
}

void PHGarfieldRossegger::setModeTruncation(unsigned int m_phi_max, unsigned int n_radial_modes, unsigned int n_longitudinal_modes)
{
  m_mPhiMax = m_phi_max; m_nRadialModes = n_radial_modes; m_nLongitudinalModes = n_longitudinal_modes;
}

void PHGarfieldRossegger::setRadialJob(unsigned int job_index, unsigned int n_jobs)
{
  m_jobIndex = job_index; m_nJobs = std::max(1U, n_jobs);
}

int PHGarfieldRossegger::Init(PHCompositeNode*)
{
  return InitRun(nullptr);
}

int PHGarfieldRossegger::InitRun(PHCompositeNode*)
{
  if (m_done) { return Fun4AllReturnCodes::EVENT_OK; }
  const int status = calculate();
  m_done = status == Fun4AllReturnCodes::EVENT_OK;
  return status;
}

int PHGarfieldRossegger::process_event(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHGarfieldRossegger::validateConfig() const
{
  const double a = m_aCm * cmToM;
  const double b = m_bCm * cmToM;
  const double smin = m_sourceRMinCm * cmToM;
  const double smax = m_sourceRMaxCm * cmToM;
  if (!(0.0 < a && a < smin && smin < smax && smax <= b))
  {
    std::cout << PHWHERE << " invalid Rossegger geometry/source radii" << std::endl;
    return false;
  }
  if (m_lCm <= 0.0 || m_nrSource == 0 || m_nphiSource == 0 || m_nzSource == 0 ||
      m_nrObs == 0 || m_nphiObs == 0 || m_nzObs == 0 || m_nRadialModes == 0 || m_nLongitudinalModes == 0)
  {
    std::cout << PHWHERE << " invalid Rossegger grid or mode configuration" << std::endl;
    return false;
  }
  if (m_jobIndex >= m_nJobs)
  {
    std::cout << PHWHERE << " job index outside number of jobs" << std::endl;
    return false;
  }
  return true;
}

unsigned int PHGarfieldRossegger::effectivePhiMax() const
{
  const bool uniform = std::fabs(m_m1Amplitude) < 1.0e-15 && std::fabs(m_m12Amplitude) < 1.0e-15;
  return (m_autoAxisymmetric && uniform) ? 0 : m_mPhiMax;
}

std::pair<unsigned int, unsigned int> PHGarfieldRossegger::radialRange() const
{
  return {(m_nrObs * m_jobIndex) / m_nJobs, (m_nrObs * (m_jobIndex + 1)) / m_nJobs};
}

double PHGarfieldRossegger::radialBoundaryFunction(double kval, unsigned int mode_m) const
{
  const double a = m_aCm * cmToM;
  const double b = m_bCm * cmToM;
  return std::cyl_bessel_j(mode_m, kval * b) * std::cyl_neumann(mode_m, kval * a) -
         std::cyl_neumann(mode_m, kval * b) * std::cyl_bessel_j(mode_m, kval * a);
}

double PHGarfieldRossegger::radialBasis(double kval, unsigned int mode_m, double radius_m) const
{
  const double a = m_aCm * cmToM;
  return std::cyl_bessel_j(mode_m, kval * radius_m) * std::cyl_neumann(mode_m, kval * a) -
         std::cyl_neumann(mode_m, kval * radius_m) * std::cyl_bessel_j(mode_m, kval * a);
}

double PHGarfieldRossegger::radialBasisDerivative(double kval, unsigned int mode_m, double radius_m) const
{
  const double a = m_aCm * cmToM;
  return kval * (jprime(mode_m, kval * radius_m) * std::cyl_neumann(mode_m, kval * a) -
                 yprime(mode_m, kval * radius_m) * std::cyl_bessel_j(mode_m, kval * a));
}

std::vector<double> PHGarfieldRossegger::findRadialRoots(unsigned int mode_m, unsigned int n_roots) const
{
  const double width = (m_bCm - m_aCm) * cmToM;
  const double step = 0.20 * std::numbers::pi / width;
  double left_k = std::max(1.0e-6, 0.10 * std::numbers::pi / width);
  double left_val = radialBoundaryFunction(left_k, mode_m);
  const double max_k = static_cast<double>(n_roots + mode_m + 20) * std::numbers::pi / width;
  std::vector<double> roots;

  for (double right_k = left_k + step; right_k <= max_k && roots.size() < n_roots; right_k += step)
  {
    const double right_val = radialBoundaryFunction(right_k, mode_m);
    if (std::isfinite(left_val) && std::isfinite(right_val) && left_val * right_val < 0.0)
    {
      double lo = right_k - step;
      double hi = right_k;
      double flo = radialBoundaryFunction(lo, mode_m);
      for (unsigned int iter = 0; iter < 120; ++iter)
      {
        const double mid = 0.5 * (lo + hi);
        const double fmid = radialBoundaryFunction(mid, mode_m);
        if (std::fabs(fmid) < 1.0e-14 || std::fabs(hi - lo) < 1.0e-12) { lo = mid; hi = mid; break; }
        if (flo * fmid <= 0.0) { hi = mid; }
        else { lo = mid; flo = fmid; }
      }
      const double root = 0.5 * (lo + hi);
      if (roots.empty() || std::fabs(root - roots.back()) > 1.0e-8) { roots.push_back(root); }
    }
    left_k = right_k;
    left_val = right_val;
  }
  if (roots.size() != n_roots)
  {
    throw std::runtime_error(std::format("Found {} roots for m={}, requested {}", roots.size(), mode_m, n_roots));
  }
  return roots;
}

std::vector<double> PHGarfieldRossegger::legendreRootsAndWeights(unsigned int n_points, std::vector<double>& weights) const
{
  std::vector<double> roots(n_points);
  weights.assign(n_points, 0.0);
  const unsigned int half = (n_points + 1) / 2;
  for (unsigned int i = 0; i < half; ++i)
  {
    double z = std::cos(std::numbers::pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(n_points) + 0.5));
    double old_z = 0.0;
    double pp = 0.0;
    while (std::fabs(z - old_z) > 1.0e-15)
    {
      double p1 = 1.0;
      double p2 = 0.0;
      for (unsigned int j = 1; j <= n_points; ++j)
      {
        const double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * static_cast<double>(j) - 1.0) * z * p2 - (static_cast<double>(j) - 1.0) * p3) / static_cast<double>(j);
      }
      pp = static_cast<double>(n_points) * (z * p1 - p2) / (z * z - 1.0);
      old_z = z;
      z = old_z - p1 / pp;
    }
    roots[i] = -z;
    roots[n_points - 1 - i] = z;
    const double weight = 2.0 / ((1.0 - z * z) * pp * pp);
    weights[i] = weight;
    weights[n_points - 1 - i] = weight;
  }
  return roots;
}

int PHGarfieldRossegger::calculate()
{
  if (!validateConfig()) { return Fun4AllReturnCodes::ABORTRUN; }
  try
  {
    const double a = m_aCm * cmToM;
    const double b = m_bCm * cmToM;
    const double len = m_lCm * cmToM;
    const double smin = m_sourceRMinCm * cmToM;
    const double smax = m_sourceRMaxCm * cmToM;
    const unsigned int mmax = effectivePhiMax();
    const auto [r_begin, r_end] = radialRange();

    std::cout << Name() << " Rossegger field calculation" << std::endl;
    std::cout << "  radial job " << m_jobIndex << "/" << m_nJobs << " bins [" << r_begin << ", " << r_end << ")" << std::endl;

    const auto rse = edges(smin, smax, m_nrSource);
    const auto pse = phi_edges(m_nphiSource);
    const auto zse = edges(0.0, len, m_nzSource);
    const auto roe = edges(a, b, m_nrObs);
    const auto poe = phi_edges(m_nphiObs);
    const auto zoe = edges(0.0, len, m_nzObs);
    const auto rs = centers(rse);
    const auto ps = phi_centers(m_nphiSource);
    const auto zs = centers(zse);
    const auto ro = centers(roe);
    const auto po = phi_centers(m_nphiObs);
    const auto zo = centers(zoe);

    std::vector<double> volumes(m_nrSource * m_nphiSource * m_nzSource, 0.0);
    std::vector<double> raw(volumes.size(), 0.0);
    double shape_sum = 0.0;
    double volume_sum = 0.0;
    for (unsigned int irs = 0; irs < m_nrSource; ++irs)
    {
      const double rvol = 0.5 * (rse[irs + 1] * rse[irs + 1] - rse[irs] * rse[irs]);
      const double rshape = std::pow(smin / rs[irs], m_radialPowerAlpha);
      for (unsigned int ips = 0; ips < m_nphiSource; ++ips)
      {
        const double pshape = 1.0 + m_m1Amplitude * std::cos(ps[ips] - m_m1Phase) + m_m12Amplitude * std::cos(12.0 * (ps[ips] - m_m12Phase));
        if (pshape < 0.0) { throw std::runtime_error("Phi modulation makes rho negative"); }
        for (unsigned int izs = 0; izs < m_nzSource; ++izs)
        {
          const unsigned int idx = src_index(irs, ips, izs, m_nphiSource, m_nzSource);
          volumes[idx] = rvol * (pse[ips + 1] - pse[ips]) * (zse[izs + 1] - zse[izs]);
          raw[idx] = rshape * pshape;
          shape_sum += raw[idx] * volumes[idx];
          volume_sum += volumes[idx];
        }
      }
    }
    const double target_rho = m_kEff * m_rhoReferenceNCPerM3 * 1.0e-9;
    const double shape_mean = shape_sum / volume_sum;
    std::vector<double> rho(raw.size(), 0.0);
    double total_charge = 0.0;
    for (unsigned int idx = 0; idx < rho.size(); ++idx)
    {
      rho[idx] = target_rho * raw[idx] / shape_mean;
      total_charge += rho[idx] * volumes[idx];
    }
    std::cout << "  mean rho = " << target_rho << " C/m^3, charge = " << total_charge << " C, ions = " << total_charge / elementaryCharge << std::endl;

    std::vector<std::vector<double>> km(mmax + 1);
    std::vector<std::vector<double>> norms(mmax + 1);
    std::vector<double> wg;
    const auto xg = legendreRootsAndWeights(192, wg);
    std::vector<double> rg(xg.size(), 0.0);
    const double jac = 0.5 * (b - a);
    for (unsigned int ig = 0; ig < xg.size(); ++ig) { rg[ig] = jac * xg[ig] + 0.5 * (b + a); }
    for (unsigned int im = 0; im <= mmax; ++im)
    {
      km[im] = findRadialRoots(im, m_nRadialModes);
      norms[im].assign(m_nRadialModes, 0.0);
      for (unsigned int in = 0; in < m_nRadialModes; ++in)
      {
        double sum = 0.0;
        for (unsigned int ig = 0; ig < xg.size(); ++ig)
        {
          const double rb = radialBasis(km[im][in], im, rg[ig]);
          sum += wg[ig] * rg[ig] * rb * rb;
        }
        norms[im][in] = jac * sum;
      }
    }

    std::vector<double> q(m_nLongitudinalModes, 0.0);
    for (unsigned int il = 0; il < m_nLongitudinalModes; ++il) { q[il] = static_cast<double>(il + 1) * std::numbers::pi / len; }

    std::vector<std::vector<double>> mc(mmax + 1), ms(mmax + 1);
    for (unsigned int im = 0; im <= mmax; ++im)
    {
      mc[im].assign(m_nRadialModes * m_nLongitudinalModes, 0.0);
      ms[im].assign(m_nRadialModes * m_nLongitudinalModes, 0.0);
      for (unsigned int in = 0; in < m_nRadialModes; ++in)
      {
        for (unsigned int il = 0; il < m_nLongitudinalModes; ++il)
        {
          double cproj = 0.0;
          double sproj = 0.0;
          for (unsigned int irs = 0; irs < m_nrSource; ++irs)
          {
            const double rb = radialBasis(km[im][in], im, rs[irs]);
            for (unsigned int ips = 0; ips < m_nphiSource; ++ips)
            {
              const double cp = std::cos(static_cast<double>(im) * ps[ips]);
              const double sp = std::sin(static_cast<double>(im) * ps[ips]);
              for (unsigned int izs = 0; izs < m_nzSource; ++izs)
              {
                const unsigned int idx = src_index(irs, ips, izs, m_nphiSource, m_nzSource);
                const double charge = rho[idx] * volumes[idx];
                const double zfactor = std::sin(zs[izs] * q[il]);
                cproj += charge * rb * cp * zfactor;
                sproj += charge * rb * sp * zfactor;
              }
            }
          }
          mc[im][mode_index(in, il, m_nLongitudinalModes)] = cproj;
          ms[im][mode_index(in, il, m_nLongitudinalModes)] = im == 0 ? 0.0 : sproj;
        }
      }
    }

    const unsigned int fsize = m_nrObs * m_nphiObs * m_nzObs;
    std::vector<double> phi(fsize, 0.0), er(fsize, 0.0), ep(fsize, 0.0), ez(fsize, 0.0);
    for (unsigned int iro = r_begin; iro < r_end; ++iro)
    {
      for (unsigned int im = 0; im <= mmax; ++im)
      {
        const double anorm = im == 0 ? 2.0 * std::numbers::pi : std::numbers::pi;
        std::vector<double> rb(m_nRadialModes, 0.0), drb(m_nRadialModes, 0.0);
        for (unsigned int in = 0; in < m_nRadialModes; ++in)
        {
          rb[in] = radialBasis(km[im][in], im, ro[iro]);
          drb[in] = radialBasisDerivative(km[im][in], im, ro[iro]);
        }
        for (unsigned int ipo = 0; ipo < m_nphiObs; ++ipo)
        {
          const double cp = std::cos(static_cast<double>(im) * po[ipo]);
          const double sp = std::sin(static_cast<double>(im) * po[ipo]);
          for (unsigned int izo = 0; izo < m_nzObs; ++izo)
          {
            const unsigned int idx = fld_index(iro, ipo, izo, m_nphiObs, m_nzObs);
            for (unsigned int in = 0; in < m_nRadialModes; ++in)
            {
              for (unsigned int il = 0; il < m_nLongitudinalModes; ++il)
              {
                const double kval = km[im][in];
                const double qval = q[il];
                const unsigned int midx = mode_index(in, il, m_nLongitudinalModes);
                const double denom = epsilon0 * norms[im][in] * (len / 2.0) * anorm * (kval * kval + qval * qval);
                const double ac = mc[im][midx] / denom;
                const double as = ms[im][midx] / denom;
                const double amp = ac * cp + as * sp;
                const double sz = std::sin(zo[izo] * qval);
                const double qcz = qval * std::cos(zo[izo] * qval);
                phi[idx] += rb[in] * sz * amp;
                er[idx] -= drb[in] * sz * amp;
                ez[idx] -= rb[in] * qcz * amp;
                if (im > 0) { ep[idx] += (static_cast<double>(im) / ro[iro]) * rb[in] * sz * (ac * sp - as * cp); }
              }
            }
          }
        }
      }
    }

    writeGarfieldRootFile(roe, zoe, er, ez, r_begin, r_end);
    if (m_writeField3D) { writeField3DRootFile(rse, roe, pse, poe, zse, zoe, rho, phi, er, ep, ez, r_begin, r_end); }
    if (m_verifyOutput && !verifyOutput()) { return Fun4AllReturnCodes::ABORTRUN; }
  }
  catch (const std::exception& e)
  {
    std::cout << PHWHERE << " " << Name() << " failed: " << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGarfieldRossegger::writeGarfieldRootFile(const std::vector<double>& r_edges_m, const std::vector<double>& z_edges_m, const std::vector<double>& er, const std::vector<double>& ez, unsigned int r_begin, unsigned int r_end) const
{
  auto rcm = r_edges_m;
  auto zcm = z_edges_m;
  scale_to_cm(rcm);
  scale_to_cm(zcm);

  TFile file(m_garfieldOutputFile.c_str(), "RECREATE");
  if (!file.IsOpen() || file.IsZombie()) { throw std::runtime_error("Could not create " + m_garfieldOutputFile); }
  file.mkdir("QA");
  file.cd("QA");
  TH2D hEr("hErDefault", "Axisymmetric radial field for PHGarfield", m_nrObs, rcm.data(), m_nzObs, zcm.data());
  TH2D hEz("hEzDefault", "Axisymmetric local longitudinal field for PHGarfield", m_nrObs, rcm.data(), m_nzObs, zcm.data());
  hEr.GetXaxis()->SetTitle("r [cm]"); hEr.GetYaxis()->SetTitle("|z| [cm]"); hEr.GetZaxis()->SetTitle("E_{r} [V/m]");
  hEz.GetXaxis()->SetTitle("r [cm]"); hEz.GetYaxis()->SetTitle("|z| [cm]"); hEz.GetZaxis()->SetTitle("E_{z}^{local} [V/m]");

  for (unsigned int ir = r_begin; ir < r_end; ++ir)
  {
    for (unsigned int iz = 0; iz < m_nzObs; ++iz)
    {
      double ersum = 0.0;
      double ezsum = 0.0;
      for (unsigned int ip = 0; ip < m_nphiObs; ++ip)
      {
        const unsigned int idx = fld_index(ir, ip, iz, m_nphiObs, m_nzObs);
        ersum += er[idx];
        ezsum += ez[idx];
      }
      hEr.SetBinContent(ir + 1, iz + 1, ersum / static_cast<double>(m_nphiObs));
      hEz.SetBinContent(ir + 1, iz + 1, ezsum / static_cast<double>(m_nphiObs));
    }
  }
  hEr.Write();
  hEz.Write();
  file.cd();
  const std::string meta = std::format("{{\"format\":\"PHGarfield axisymmetric correction map\",\"histograms\":{{\"Er\":\"QA/hErDefault\",\"Ez\":\"QA/hEzDefault\"}},\"field_units\":\"V/m\",\"job_index\":{},\"n_jobs\":{},\"radial_bin_begin\":{},\"radial_bin_end\":{}}}", m_jobIndex, m_nJobs, r_begin, r_end);
  TNamed("metadata_json", meta.c_str()).Write();
  file.Close();
  std::cout << Name() << " wrote " << m_garfieldOutputFile << std::endl;
}

void PHGarfieldRossegger::writeField3DRootFile(const std::vector<double>& r_source_edges_m, const std::vector<double>& r_obs_edges_m, const std::vector<double>& phi_source_edges, const std::vector<double>& phi_obs_edges, const std::vector<double>& z_source_edges_m, const std::vector<double>& z_obs_edges_m, const std::vector<double>& rho, const std::vector<double>& potential, const std::vector<double>& er, const std::vector<double>& ephi, const std::vector<double>& ez, unsigned int r_begin, unsigned int r_end) const
{
  auto rscm = r_source_edges_m;
  auto rocm = r_obs_edges_m;
  auto zscm = z_source_edges_m;
  auto zocm = z_obs_edges_m;
  scale_to_cm(rscm); scale_to_cm(rocm); scale_to_cm(zscm); scale_to_cm(zocm);
  TFile file(m_field3DOutputFile.c_str(), "RECREATE");
  if (!file.IsOpen() || file.IsZombie()) { throw std::runtime_error("Could not create " + m_field3DOutputFile); }
  TH3D hRho("hRho", "Positive-ion density", m_nrSource, rscm.data(), m_nphiSource, phi_source_edges.data(), m_nzSource, zscm.data());
  TH3D hPhi("hPhi", "Space-charge potential", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEr("hEr", "Radial electric field", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEphi("hEphi", "Azimuthal electric field", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  TH3D hEz("hEz", "Longitudinal electric field", m_nrObs, rocm.data(), m_nphiObs, phi_obs_edges.data(), m_nzObs, zocm.data());
  for (unsigned int ir = 0; ir < m_nrSource; ++ir)
    for (unsigned int ip = 0; ip < m_nphiSource; ++ip)
      for (unsigned int iz = 0; iz < m_nzSource; ++iz)
        hRho.SetBinContent(ir + 1, ip + 1, iz + 1, rho[src_index(ir, ip, iz, m_nphiSource, m_nzSource)]);
  for (unsigned int ir = r_begin; ir < r_end; ++ir)
    for (unsigned int ip = 0; ip < m_nphiObs; ++ip)
      for (unsigned int iz = 0; iz < m_nzObs; ++iz)
      {
        const unsigned int idx = fld_index(ir, ip, iz, m_nphiObs, m_nzObs);
        hPhi.SetBinContent(ir + 1, ip + 1, iz + 1, potential[idx]);
        hEr.SetBinContent(ir + 1, ip + 1, iz + 1, er[idx]);
        hEphi.SetBinContent(ir + 1, ip + 1, iz + 1, ephi[idx]);
        hEz.SetBinContent(ir + 1, ip + 1, iz + 1, ez[idx]);
      }
  hRho.Write(); hPhi.Write(); hEr.Write(); hEphi.Write(); hEz.Write();
  TNamed("metadata_json", "Rossegger modal solution for sPHENIX TPC").Write();
  file.Close();
  std::cout << Name() << " wrote " << m_field3DOutputFile << std::endl;
}

bool PHGarfieldRossegger::verifyOutput() const
{
  TFile file(m_garfieldOutputFile.c_str(), "READ");
  if (!file.IsOpen() || file.IsZombie()) { return false; }
  return file.Get("QA/hErDefault") && file.Get("QA/hEzDefault") && file.Get("metadata_json");
}
