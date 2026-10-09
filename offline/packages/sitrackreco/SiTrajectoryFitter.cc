#include "SiTrajectoryFitter.h"

#include "SiHitSeedData.h"
#include "Si_Trajectoryv1.h"
#include "Si_TrajectoryContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>

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
}  // namespace

SiTrajectoryFitter::SiTrajectoryFitter(const std::string& name)
  : SubsysReco(name)
{
}

int SiTrajectoryFitter::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    const auto& c = m_frame.detectorCenterCm();
    std::cout << Name() << ": Uphi rotation " << m_frame.uphiRotation() << " rad"
              << ", detector centre (" << 10 * c[0] << ", " << 10 * c[1] << ", " << 10 * c[2] << ") mm"
              << ", beam (" << m_beamX << ", " << m_beamY << ") cm"
              << ", Bz " << m_bz << " T"
              << ", use INTT " << m_useIntt << ", min points " << m_minPoints
              << ", straight line if radial lever arm < " << m_minCurvatureLeverArm << " cm"
              << ", beam constraint " << (m_beamConstraint ? "on" : "off") << std::endl;
    if (m_alignment.enabled())
    {
      const char* names[] = {"TPC", "Si half A", "Si half B"};
      for (int part = SiTpcBeamAlignment::SiHalfA; part <= SiTpcBeamAlignment::SiHalfB; ++part)
      {
        const auto& l = m_alignment.beamLine(part);
        std::cout << Name() << ":   beam-axis alignment " << names[part] << ": x0 " << 10 * l.x0 << " mm, y0 "
                  << 10 * l.y0 << " mm, dx/dz " << 1e3 * l.dxdz << " mrad, dy/dz " << 1e3 * l.dydz << " mrad" << std::endl;
      }
    }
    else
    {
      std::cout << Name() << ":   beam-axis alignment off" << std::endl;
    }
  }
  return createNodes(topNode);
}

int SiTrajectoryFitter::createNodes(PHCompositeNode* topNode)
{
  m_container = findNode::getClass<Si_TrajectoryContainer>(topNode, m_outputNodeName);
  if (m_container)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << ": DST node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator dstIter(dstNode);
  auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_container = new Si_TrajectoryContainerv1();
  svtxNode->addNode(new PHIODataNode<PHObject>(m_container, m_outputNodeName, "PHObject"));

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": created DST/SVTX/" << m_outputNodeName << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTrajectoryFitter::process_event(PHCompositeNode* topNode)
{
  auto* event = findNode::getClass<SiHitSeedEvent>(topNode, m_inputNodeName);
  if (!event)
  {
    std::cerr << Name() << ": missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!m_container)
  {
    m_container = findNode::getClass<Si_TrajectoryContainer>(topNode, m_outputNodeName);
    if (!m_container)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  m_container->Reset();

  std::array<unsigned int, 5> nStatus{};  // NotFitted, Ok, TooFewPoints, Degenerate, StraightLine
  unsigned int nLowPt = 0;
  for (const auto& chain : event->chains)
  {
    auto* trajectory = new Si_Trajectoryv1();
    trajectory->set_id(m_container->size());
    trajectory->set_chain_id(chain.id);

    for (int hid : chain.hit_ids)
    {
      if (hid < 0 || hid >= static_cast<int>(event->hits.size()))
      {
        continue;
      }
      const auto& h = event->hits[hid];
      if (!m_useIntt && h.layer > 2)
      {
        continue;
      }
      trajectory->add_hit_index(h.hitsetkey, h.hitkey);
    }

    const auto points = layerPoints(*event, chain.hit_ids);
    unsigned int nMvtx = 0;
    unsigned int nIntt = 0;
    for (const auto& p : points)
    {
      trajectory->add_point(p.layer, p.global.x, p.global.y, p.global.z, p.nhits);
      (p.layer <= 2 ? nMvtx : nIntt)++;
    }
    trajectory->set_n_mvtx(nMvtx);
    trajectory->set_n_intt(nIntt);

    fit(points, *trajectory);
    const int status = trajectory->get_fit_status();
    ++nStatus[std::clamp(status, 0, 4)];
    const bool lowPt = status == Si_Trajectory::Ok && trajectory->get_pt() < m_reportPt;
    nLowPt += lowPt;
    // Verbosity 2: every chain that is not fitted or is fitted with a low pt
    if (Verbosity() > 2 || (Verbosity() > 1 && (!trajectory->is_fitted() || lowPt)))
    {
      std::cout << Name() << ":   chain " << chain.id << " status " << status
                << " pt " << trajectory->get_pt() << " circle_rms " << trajectory->get_circle_rms() << " cm, layers:";
      for (const auto& p : points)
      {
        std::cout << " " << p.layer << "(" << p.nhits << " hits, r=" << std::hypot(p.global.x, p.global.y) << ")";
      }
      std::cout << std::endl;
      if (Verbosity() > 2)
      {
        trajectory->identify();
      }
    }
    m_container->add_trajectory(trajectory);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << event->event << " chains=" << event->chains.size()
              << " circle=" << nStatus[Si_Trajectory::Ok] << " (pt<" << m_reportPt << " GeV: " << nLowPt << ")"
              << " straight=" << nStatus[Si_Trajectory::StraightLine]
              << " too_few_points=" << nStatus[Si_Trajectory::TooFewPoints]
              << " degenerate=" << nStatus[Si_Trajectory::Degenerate] << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

// hits -> detector frame -> global frame, then one centroid per layer (on the layer cylinder
// in the detector frame, so the shift is applied to a point that lies on the nominal layer).
std::vector<SiTrajectoryFitter::LayerPoint> SiTrajectoryFitter::layerPoints(
    const SiHitSeedEvent& event, const std::vector<int>& hitIds) const
{
  struct Sum
  {
    double x = 0, y = 0, z = 0;
    unsigned int n = 0;
  };
  std::map<int, Sum> sums;  // ordered by layer = ordered outward
  for (int hid : hitIds)
  {
    if (hid < 0 || hid >= static_cast<int>(event.hits.size()))
    {
      continue;
    }
    const auto& h = event.hits[hid];
    if (h.layer < 0 || h.layer >= SiDetectorFrame::kNLayers || (!m_useIntt && h.layer > 2))
    {
      continue;
    }
    // step 1: detector frame
    const auto det = m_frame.toDetector(h.layer, h.phi, h.z);
    auto& s = sums[h.layer];
    s.x += det.x;
    s.y += det.y;
    s.z += det.z;
    ++s.n;
  }

  std::vector<LayerPoint> out;
  out.reserve(sums.size());
  for (const auto& [layer, s] : sums)
  {
    const double r = SiDetectorFrame::layerRadius(layer);
    const double phi = std::atan2(s.y, s.x);
    const SiDetectorFrame::Point det{r * std::cos(phi), r * std::sin(phi), s.z / s.n};
    LayerPoint p;
    p.layer = layer;
    // step 2: detector shift (global frame); step 3: beam-axis alignment of its clamshell half
    p.global = m_alignment.toBeamAxis(m_alignment.siHalf(phi - m_frame.uphiRotation()), m_frame.toGlobal(det));
    p.nhits = s.n;
    out.push_back(p);
  }
  return out;
}

void SiTrajectoryFitter::fit(const std::vector<LayerPoint>& points, Si_Trajectory& t) const
{
  const unsigned int n = points.size();
  if (n < std::max(3U, m_minPoints))
  {
    t.set_fit_status(Si_Trajectory::TooFewPoints);
    return;
  }

  std::vector<double> xs, ys, ws;
  double rmin = std::numeric_limits<double>::max(), rmax = 0;
  for (const auto& p : points)
  {
    xs.push_back(p.global.x);
    ys.push_back(p.global.y);
    ws.push_back(1.0);
    const double r = std::hypot(p.global.x - m_beamX, p.global.y - m_beamY);
    rmin = std::min(rmin, r);
    rmax = std::max(rmax, r);
  }
  // optional beam-spot constraint: the beam position enters the xy fit as an extra point
  if (m_beamConstraint)
  {
    xs.push_back(m_beamX);
    ys.push_back(m_beamY);
    ws.push_back(m_beamWeight);
  }
  const bool straight = (rmax - rmin) < m_minCurvatureLeverArm;

  // ---- xy model.  Both models end with: pca, unit tangent (tx, ty) at the pca along the
  // direction of motion (outward), and the transverse path length s of every point.
  double pcaX = 0, pcaY = 0, tx = 0, ty = 0;
  std::vector<double> s(n);
  if (!straight)
  {
    const Circle c = fitCircleTaubin(xs, ys, ws);
    if (!c.ok)
    {
      t.set_fit_status(Si_Trajectory::Degenerate);
      return;
    }
    // helicity from the turning sense going outward (points are ordered by layer)
    const double ax = xs.front() - c.x, ay = ys.front() - c.y;
    const double bx = xs[n - 1] - c.x, by = ys[n - 1] - c.y;
    const int h = (ax * by - ay * bx) >= 0 ? +1 : -1;

    const double dx = m_beamX - c.x, dy = m_beamY - c.y;
    const double dist = std::hypot(dx, dy);
    if (dist <= 0)
    {
      t.set_fit_status(Si_Trajectory::Degenerate);
      return;
    }
    const double ux = dx / dist, uy = dy / dist;  // centre -> beam
    pcaX = c.x + c.r * ux;
    pcaY = c.y + c.r * uy;
    tx = -h * uy;
    ty = h * ux;

    double ss = 0;
    const double a0 = std::atan2(pcaY - c.y, pcaX - c.x);
    for (unsigned int i = 0; i < n; ++i)
    {
      const double d = std::hypot(xs[i] - c.x, ys[i] - c.y) - c.r;
      ss += d * d;
      s[i] = c.r * h * wrapPi(std::atan2(ys[i] - c.y, xs[i] - c.x) - a0);
    }
    t.set_circle_x(c.x);
    t.set_circle_y(c.y);
    t.set_radius(c.r);
    t.set_helicity(h);
    t.set_circle_rms(std::sqrt(ss / n));
    t.set_dca(((pcaX - m_beamX) * ty - (pcaY - m_beamY) * tx >= 0 ? -1.0 : 1.0) * std::abs(dist - c.r));
    t.set_pt(0.003 * std::abs(m_bz) * c.r);  // GeV, R in cm, B in T
    t.set_charge(m_bz == 0 ? 0 : -h * (m_bz > 0 ? 1 : -1));
  }
  else
  {
    // straight line: weighted total least squares (principal axis)
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
    tx = std::cos(ang);
    ty = std::sin(ang);
    if ((xs[n - 1] - xs.front()) * tx + (ys[n - 1] - ys.front()) * ty < 0)  // outward
    {
      tx = -tx;
      ty = -ty;
    }
    const double proj = (m_beamX - mx) * tx + (m_beamY - my) * ty;
    pcaX = mx + proj * tx;
    pcaY = my + proj * ty;

    double ss = 0;
    for (unsigned int i = 0; i < n; ++i)
    {
      const double d = (xs[i] - pcaX) * (-ty) + (ys[i] - pcaY) * tx;
      ss += d * d;
      s[i] = (xs[i] - pcaX) * tx + (ys[i] - pcaY) * ty;
    }
    // stored as a circle of very large radius (helicity +1: centre to the left), so every
    // consumer can keep using the circle parametrization
    t.set_circle_x(pcaX - kStraightRadius * ty);
    t.set_circle_y(pcaY + kStraightRadius * tx);
    t.set_radius(kStraightRadius);
    t.set_helicity(+1);
    t.set_circle_rms(std::sqrt(ss / n));
    t.set_dca(((pcaX - m_beamX) * ty - (pcaY - m_beamY) * tx >= 0 ? -1.0 : 1.0) *
              std::hypot(pcaX - m_beamX, pcaY - m_beamY));
    t.set_pt(NAN);  // curvature not measured
    t.set_charge(0);
  }
  t.set_pca_x(pcaX);
  t.set_pca_y(pcaY);
  t.set_phi(std::atan2(ty, tx));

  // ---- line z(s)
  double sS = 0, sZ = 0, sSS = 0, sSZ = 0;
  for (unsigned int i = 0; i < n; ++i)
  {
    const double z = points[i].global.z;
    sS += s[i];
    sZ += z;
    sSS += s[i] * s[i];
    sSZ += s[i] * z;
  }
  const double det = n * sSS - sS * sS;
  if (std::abs(det) < 1e-12)
  {
    t.set_fit_status(Si_Trajectory::Degenerate);
    return;
  }
  const double tanl = (n * sSZ - sS * sZ) / det;
  const double z0 = (sZ - tanl * sS) / n;
  t.set_tanl(tanl);
  t.set_z0(z0);

  double zz = 0;
  for (unsigned int i = 0; i < n; ++i)
  {
    const double d = points[i].global.z - (z0 + tanl * s[i]);
    zz += d * d;
  }
  t.set_z_rms(std::sqrt(zz / n));

  t.set_fit_status(straight ? Si_Trajectory::StraightLine : Si_Trajectory::Ok);
}

// Taubin algebraic circle fit, Newton iteration on the characteristic polynomial
// (N. Chernov, "Circular and linear regression", CircleFitByTaubin).
SiTrajectoryFitter::Circle SiTrajectoryFitter::fitCircleTaubin(const std::vector<double>& x,
                                                                const std::vector<double>& y)
{
  return fitCircleTaubin(x, y, std::vector<double>(std::min(x.size(), y.size()), 1.0));
}

// Weighted version: the moments are weighted averages.
SiTrajectoryFitter::Circle SiTrajectoryFitter::fitCircleTaubin(const std::vector<double>& x,
                                                                const std::vector<double>& y,
                                                                const std::vector<double>& w)
{
  Circle out;
  const std::size_t n = std::min({x.size(), y.size(), w.size()});
  if (n < 3)
  {
    return out;
  }

  double sw = 0, mx = 0, my = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    sw += w[i];
    mx += w[i] * x[i];
    my += w[i] * y[i];
  }
  if (sw <= 0)
  {
    return out;
  }
  mx /= sw;
  my /= sw;

  double Mxx = 0, Myy = 0, Mxy = 0, Mxz = 0, Myz = 0, Mzz = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
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
