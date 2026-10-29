#include "SiTpc_EventDisplay.h"

#include "tpctrackreco/TpcCrossingDecision.h"
#include "tpctrackreco/TpcCrossingDecisionContainer.h"
#include "tpctrackreco/Tpc_PolyCluster.h"
#include "tpctrackreco/Tpc_PolyClusterContainer.h"
#include "tpctrackreco/Tpc_PolyTrack.h"
#include "tpctrackreco/Tpc_PolyTrackContainer.h"
#include "tpctrackreco/Tpc_PolyTrackVertex.h"
#include "tpctrackreco/Tpc_PolyTrackVertexContainer.h"

#include <sitrackreco/SiHitSeedData.h>
#include <sitrackreco/Si_Trajectory.h>
#include <sitrackreco/Si_TrajectoryContainer.h>
#include <sitrackreco/SiTpcHelixFit.h>
#include <sitrackreco/SiTpc_Track.h>
#include <sitrackreco/SiTpc_TrackContainer.h>

#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackSeed.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace
{
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2.0 * kPi;
  constexpr double kVertexLayer = -1.0;  // detector view: "layer" of the vertex plane (R = 0)
  constexpr double kSiZoomR = 13.0;      // physics zoom on the silicon [cm]

  // =====================================================================
  // generic helpers
  // =====================================================================
  int palette_color(const int i)
  {
    static const int colors[] = {
        kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
        kOrange + 7, kViolet + 1, kAzure + 1, kPink + 7, kTeal + 3};
    const int n = sizeof(colors) / sizeof(colors[0]);
    return colors[((i % n) + n) % n];
  }

  // Objects drawn on a canvas get kCanDelete: deleting the canvas frees them.
  template <class T>
  T* owned(T* obj)
  {
    obj->SetBit(kCanDelete);
    return obj;
  }

  void write_and_delete(TCanvas* c)
  {
    c->Modified();
    c->Update();
    c->Write();
    delete c;
  }

  struct XYZ
  {
    double x = 0, y = 0, z = 0;
  };

  double radius_of(const XYZ& p) { return std::hypot(p.x, p.y); }

  // 3D polyline (z, x, y), as Tpc_PolyClusterDisplay
  void draw_line_zxy(const std::vector<XYZ>& pts, const int color, const int style, const int width)
  {
    if (pts.size() < 2)
    {
      return;
    }
    auto* line = owned(new TPolyLine3D(static_cast<int>(pts.size())));
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
      line->SetPoint(static_cast<int>(i), pts[i].z, pts[i].x, pts[i].y);
    }
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    line->Draw("same");
  }

  void draw_markers_zxy(const std::vector<XYZ>& pts, const int color, const int style, const double size)
  {
    if (pts.empty())
    {
      return;
    }
    auto* pm = owned(new TPolyMarker3D(static_cast<int>(pts.size())));
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
      pm->SetPoint(static_cast<int>(i), pts[i].z, pts[i].x, pts[i].y);
    }
    pm->SetMarkerColor(color);
    pm->SetMarkerStyle(style);
    pm->SetMarkerSize(size);
    pm->Draw("same");
  }

  // =====================================================================
  // detector view helpers (as SiHitSeedDisplay)
  // =====================================================================
  struct P3
  {
    double z, phi, layer;
  };

  double wrapDelta(double x)
  {
    x = std::fmod(x + kPi, kTwoPi);
    if (x < 0)
    {
      x += kTwoPi;
    }
    return x - kPi;
  }

  double wrap2pi(double x)
  {
    x = std::fmod(x, kTwoPi);
    return x < 0 ? x + kTwoPi : x;
  }

  // Display coordinate of a phi value: u in [b, b+2pi).
  double toU(double phi, double b)
  {
    return b + wrap2pi(phi - b);
  }

  // Splits a (z, phi, layer) path at the display seam u = b (= b + 2pi), so every piece lies
  // inside [b, b+2pi): the path is unwrapped continuously, cut where it crosses b + 2k*pi, and
  // each piece is shifted by a multiple of 2pi into the display range.
  std::vector<std::vector<P3>> splitAtSeam(const std::vector<P3>& in, double b)
  {
    std::vector<std::vector<P3>> pieces;
    if (in.size() < 2)
    {
      return pieces;
    }
    std::vector<P3> p = in;
    p[0].phi = toU(in[0].phi, b);
    for (size_t i = 1; i < p.size(); ++i)
    {
      p[i].phi = p[i - 1].phi + wrapDelta(in[i].phi - p[i - 1].phi);
    }
    for (size_t i = 1; i < p.size(); ++i)
    {
      const P3 a = p[i - 1];
      const P3 c = p[i];
      std::vector<double> ts;
      const double lo = std::min(a.phi, c.phi), hi = std::max(a.phi, c.phi);
      for (int k = (int) std::ceil((lo - b) / kTwoPi); b + k * kTwoPi < hi; ++k)
      {
        const double edge = b + k * kTwoPi;
        if (edge > lo && edge < hi)
        {
          ts.push_back((edge - a.phi) / (c.phi - a.phi));
        }
      }
      ts.push_back(1.0);
      double t0 = 0.0;
      for (double t1 : ts)
      {
        if (t1 - t0 < 1e-12)
        {
          continue;
        }
        auto at = [&](double t)
        { return P3{a.z + t * (c.z - a.z), a.phi + t * (c.phi - a.phi), a.layer + t * (c.layer - a.layer)}; };
        const P3 s0 = at(t0), s1 = at(t1);
        const double mid = 0.5 * (s0.phi + s1.phi);
        const double offset = toU(mid, b) - mid;
        const P3 u0{s0.z, s0.phi + offset, s0.layer}, u1{s1.z, s1.phi + offset, s1.layer};
        if (!pieces.empty() && std::abs(pieces.back().back().phi - u0.phi) < 1e-9 &&
            std::abs(pieces.back().back().z - u0.z) < 1e-9)
        {
          pieces.back().push_back(u1);
        }
        else
        {
          pieces.push_back({u0, u1});
        }
        t0 = t1;
      }
    }
    return pieces;
  }

  TPolyLine3D* make_det_line(const std::vector<P3>& piece, const int color, const int style)
  {
    auto* line = owned(new TPolyLine3D(static_cast<int>(piece.size())));
    for (size_t i = 0; i < piece.size(); ++i)
    {
      line->SetPoint(static_cast<int>(i), piece[i].z, piece[i].phi, piece[i].layer);
    }
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(style == 1 ? 3 : 2);
    return line;
  }

  // physical (phi, eta, r) view: P3 is reused as {z = eta, phi, layer = r}, so splitAtSeam cuts
  // the lines at the phi seam; drawn with x = phi, y = eta, z = r
  TPolyLine3D* make_phys_line(const std::vector<P3>& piece, const int color, const int style, const int width)
  {
    auto* line = owned(new TPolyLine3D(static_cast<int>(piece.size())));
    for (size_t i = 0; i < piece.size(); ++i)
    {
      line->SetPoint(static_cast<int>(i), piece[i].phi, piece[i].z, piece[i].layer);
    }
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    return line;
  }

  TH3F* make_det_hist(const std::string& name, const std::string& title, int nz, double zmin, double zmax,
                      int nphi, double phimin, double phimax)
  {
    // layer axis also covers layer -1: the vertex plane (R = 0)
    auto* h = owned(new TH3F(name.c_str(), title.c_str(), nz, zmin, zmax, nphi, phimin, phimax, 8, -1.5, 6.5));
    h->SetDirectory(nullptr);
    h->SetStats(false);
    return h;
  }

  TCanvas* make_det_canvas(const std::string& name, const std::string& title)
  {
    auto* c = new TCanvas(name.c_str(), title.c_str(), 1150, 850);
    c->SetRightMargin(0.20);
    c->SetBottomMargin(0.16);
    c->SetTheta(24);
    c->SetPhi(35);
    return c;
  }

  // =====================================================================
  // Si trajectory helix (global frame)
  // =====================================================================
  struct Helix
  {
    bool ok = false;
    double cx = 0, cy = 0, R = 0, a0 = 0, z0 = 0, tanl = 0;
    int h = 0;
    XYZ at(double s) const
    {
      const double a = a0 + h * s / R;
      return {cx + R * std::cos(a), cy + R * std::sin(a), z0 + tanl * s};
    }
  };

  Helix helix_of(const Si_Trajectory& t)
  {
    Helix hx;
    hx.cx = t.get_circle_x();
    hx.cy = t.get_circle_y();
    hx.R = t.get_radius();
    hx.h = t.get_helicity();
    hx.z0 = t.get_z0();
    hx.tanl = t.get_tanl();
    hx.ok = std::isfinite(hx.cx) && std::isfinite(hx.cy) && std::isfinite(hx.R) && hx.R > 0 && hx.h != 0 &&
            std::isfinite(hx.z0) && std::isfinite(hx.tanl);
    if (hx.ok)
    {
      hx.a0 = std::atan2(t.get_pca_y() - hx.cy, t.get_pca_x() - hx.cx);
    }
    return hx;
  }

  // Samples the helix from the pca outward with a fixed step, while its distance from (ox, oy)
  // <= rmax and |z| <= zmax, at most half a turn (loopers).  A straight-line fit is stored as
  // a very large circle, so the path length is also capped.
  constexpr double kHelixStep = 0.05;  // [cm]
  std::vector<XYZ> sample_helix(const Helix& hx, const double ox, const double oy,
                                const double rmax, const double zmax)
  {
    std::vector<XYZ> out;
    if (!hx.ok)
    {
      return out;
    }
    const double smax = std::min(kPi * hx.R, 4.0 * rmax + 100.0);
    for (double s = 0; s <= smax; s += kHelixStep)
    {
      const XYZ p = hx.at(s);
      if (std::hypot(p.x - ox, p.y - oy) > rmax || std::fabs(p.z) > zmax)
      {
        break;
      }
      out.push_back(p);
    }
    return out;
  }

  // First crossing of the helix (from the pca outward) with the circle of radius r around
  // (ox, oy); false if it does not reach it within half a turn.
  // First crossing of the helix (from the pca outward) with the circle of radius r, where the
  // radius of a helix point p is |toFrame(p)|_xy, e.g. its radius in the Si detector frame.
  template <class Map>
  bool helix_at_radius_mapped(const Helix& hx, Map toFrame, const double r, XYZ& out)
  {
    if (!hx.ok)
    {
      return false;
    }
    auto dist = [&](double s)
    {
      const XYZ p = toFrame(hx.at(s));
      return std::hypot(p.x, p.y) - r;
    };
    const double smax = std::min(kPi * hx.R, 4.0 * r + 100.0);
    const double ds = kHelixStep;
    double s0 = 0;
    if (dist(0) >= 0)
    {
      return false;
    }
    for (double s1 = ds; s1 <= smax; s1 += ds)
    {
      if (dist(s1) >= 0)
      {
        for (int it = 0; it < 40; ++it)
        {
          const double sm = 0.5 * (s0 + s1);
          (dist(sm) < 0 ? s0 : s1) = sm;
        }
        out = hx.at(0.5 * (s0 + s1));
        return true;
      }
      s0 = s1;
    }
    return false;
  }

  // =====================================================================
  // TPC poly tracks (as Tpc_PolyClusterDisplay)
  // =====================================================================
  bool track_vertex_z_selected(const Tpc_PolyTrackVertex* vtx, const double zmin, const double zmax)
  {
    if (!vtx)
    {
      return false;
    }
    const double z0 = vtx->get_z0();
    return std::isfinite(z0) && z0 >= zmin && z0 <= zmax;
  }

  bool poly_track_xy_at_z(const Tpc_PolyTrack* trk, const double z, const double bz,
                          const double arc_direction, const bool use_straight_line, double& x, double& y)
  {
    if (!trk || trk->get_fit_status() == 0 || !std::isfinite(z))
    {
      return false;
    }
    const double x0 = trk->get_x(), y0 = trk->get_y(), z0 = trk->get_z();
    const double px = trk->get_px(), py = trk->get_py(), pz = trk->get_pz();
    const double charge = trk->get_charge();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) || !std::isfinite(px) ||
        !std::isfinite(py) || !std::isfinite(pz) || !std::isfinite(charge))
    {
      return false;
    }
    const double dz = z - z0;
    if (use_straight_line || std::fabs(charge * bz) < 1.0e-12)
    {
      if (std::fabs(pz) < 1.0e-12)
      {
        return false;
      }
      x = x0 + arc_direction * px / pz * dz;
      y = y0 + arc_direction * py / pz * dz;
      return std::isfinite(x) && std::isfinite(y);
    }
    const double pt = std::hypot(px, py);
    if (pt <= 0.0 || std::fabs(pz) < 1.0e-12)
    {
      return false;
    }
    const double signed_radius = pt / (0.003 * charge * bz);
    const double radius = std::fabs(signed_radius);
    if (!std::isfinite(radius) || radius <= 0.0)
    {
      return false;
    }
    const double tx = px / pt, ty = py / pt;
    const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
    const double xc = x0 + sign * radius * ty;
    const double yc = y0 - sign * radius * tx;
    const double phi0 = std::atan2(y0 - yc, x0 - xc);
    const double dzds = pz / pt;
    if (std::fabs(dzds) < 1.0e-12)
    {
      return false;
    }
    const double arc = arc_direction * dz / dzds;
    const double phi = phi0 - sign * arc / radius;
    x = xc + radius * std::cos(phi);
    y = yc + radius * std::sin(phi);
    return std::isfinite(x) && std::isfinite(y);
  }

  double cluster_line_residual2(const Tpc_PolyTrack* trk, const std::vector<XYZ>& clusters, const double bz,
                                const double arc_direction, const bool use_straight_line)
  {
    double sum = 0.0;
    unsigned int n = 0;
    for (const auto& c : clusters)
    {
      double x = 0, y = 0;
      if (!poly_track_xy_at_z(trk, c.z, bz, arc_direction, use_straight_line, x, y))
      {
        continue;
      }
      sum += (x - c.x) * (x - c.x) + (y - c.y) * (y - c.y);
      ++n;
    }
    return n > 0 ? sum / n : std::numeric_limits<double>::max();
  }

  // Extension of a TPC poly track inward: from its innermost cluster z (zIn), moving in z
  // away from the outermost cluster (zOut), until the radius stops decreasing (closest approach
  // to the beam axis), drops below rmin, or |z| leaves zabs.  The first point is at zIn so the
  // extension joins the solid part.
  std::vector<XYZ> extend_poly_track_inward(const Tpc_PolyTrack* trk, const double zIn, const double zOut,
                                            const double rmin, const double zabs, const double bz,
                                            const double arc_direction, const bool use_straight_line)
  {
    std::vector<XYZ> out;
    const double dir = zIn <= zOut ? -1.0 : 1.0;
    const double dz = std::clamp(std::fabs(zOut - zIn) / 400.0, 0.002, 0.05);  // fine enough for the Si zoom
    double prevR = std::numeric_limits<double>::max();
    for (int i = 0; i < 5000; ++i)
    {
      const double z = zIn + dir * dz * i;
      if (std::fabs(z) > zabs)
      {
        break;
      }
      double x = 0, y = 0;
      if (!poly_track_xy_at_z(trk, z, bz, arc_direction, use_straight_line, x, y))
      {
        break;
      }
      const double r = std::hypot(x, y);
      if (r > prevR + 1e-6)
      {
        break;  // passed the closest approach
      }
      out.push_back({x, y, z});
      if (r < rmin)
      {
        break;
      }
      prevR = r;
    }
    return out;
  }

  // Keeps the parts of a curve inside r <= rmax and |z| <= zabs (3D lines are not clipped by
  // the TH3 frame), split into connected pieces.
  std::vector<std::vector<XYZ>> clip_curve(const std::vector<XYZ>& pts, const double rmax, const double zabs)
  {
    std::vector<std::vector<XYZ>> out;
    bool inside = false;
    for (const auto& p : pts)
    {
      const bool in = radius_of(p) <= rmax && std::fabs(p.z) <= zabs;
      if (in)
      {
        if (!inside)
        {
          out.emplace_back();
        }
        out.back().push_back(p);
      }
      inside = in;
    }
    return out;
  }

  std::vector<XYZ> sample_poly_track(const Tpc_PolyTrack* trk, const double zlo, const double zhi,
                                     const double xymax, const double bz, const double arc_direction,
                                     const bool use_straight_line)
  {
    std::vector<XYZ> out;
    const unsigned int nsteps = 80;
    for (unsigned int i = 0; i <= nsteps; ++i)
    {
      const double z = zlo + (zhi - zlo) * i / static_cast<double>(nsteps);
      double x = 0, y = 0;
      if (!poly_track_xy_at_z(trk, z, bz, arc_direction, use_straight_line, x, y))
      {
        continue;
      }
      if (std::fabs(x) > xymax || std::fabs(y) > xymax)
      {
        continue;
      }
      out.push_back({x, y, z});
    }
    return out;
  }
}  // namespace

// =====================================================================
// SiTpc_EventDisplay
// =====================================================================
SiTpc_EventDisplay::SiTpc_EventDisplay(const std::string& name,
                                       const std::string& outfilename,
                                       const unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_maxEventDisplays(maxEventDisplays)
{
}

SiTpc_EventDisplay::~SiTpc_EventDisplay()
{
  delete m_outfile;
}

int SiTpc_EventDisplay::Init(PHCompositeNode* /*unused*/)
{
  // cppcheck-suppress publicAllocationError
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open output file " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_outfile->mkdir("events");
  const auto& c = m_frame.detectorCenterCm();
  std::cout << Name() << "::Init - writing up to " << m_maxEventDisplays << " events to " << m_outfilename
            << " (Si detector centre " << 10 * c[0] << ", " << 10 * c[1] << " mm)" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void SiTpc_EventDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_siSeeds = findNode::getClass<SiHitSeedEvent>(topNode, m_siSeedNodeName);
  m_siTraj = findNode::getClass<Si_TrajectoryContainer>(topNode, m_siTrajNodeName);
  m_tpcClusters = m_drawTpc ? findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName) : nullptr;
  m_tpcTracks = m_drawTpc ? findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName) : nullptr;
  m_tpcVertices = m_drawTpc ? findNode::getClass<Tpc_PolyTrackVertexContainer>(topNode, m_tpcVertexNodeName) : nullptr;
  m_siTpcTracks = findNode::getClass<SiTpc_TrackContainer>(topNode, m_siTpcTrackNodeName);
  m_decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingNodeName);
  m_actsTracks = m_drawActsTracks ? findNode::getClass<SvtxTrackMap>(topNode, m_actsTrackMapName) : nullptr;

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": nodes Si seeds " << (m_siSeeds ? "yes" : "no")
              << ", Si trajectories " << (m_siTraj ? "yes" : "no")
              << ", TPC clusters " << (m_tpcClusters ? "yes" : "no")
              << ", TPC tracks " << (m_tpcTracks ? "yes" : "no")
              << ", TPC vertices " << (m_tpcVertices ? "yes" : "no")
              << ", Si+TPC tracks " << (m_siTpcTracks ? "yes" : "no")
              << ", ACTS tracks " << (m_actsTracks ? "yes" : "no")
              << ", crossing decisions " << (m_decisions ? "yes" : "no") << std::endl;
  }
}

int SiTpc_EventDisplay::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  if (!m_outfile || m_eventsSaved >= m_maxEventDisplays)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  get_nodes(topNode);
  if (!m_siSeeds && !m_siTraj && !m_tpcClusters)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const auto& center = m_frame.detectorCenterCm();

  // colour of a Si chain (and of everything matched to it): grey if it has fewer MVTX layers
  // than setColorMinMvtxLayers
  std::map<int, int> chainNMvtx;  // chain id -> MVTX layers
  if (m_siSeeds)
  {
    for (const auto& ch : m_siSeeds->chains)
    {
      chainNMvtx[ch.id] = ch.n_mvtx;
    }
  }
  const int grey = static_cast<int>(kGray + 1);
  auto chainColor = [&](int chainId) -> int
  {
    if (m_colorMinMvtx > 0)
    {
      const auto it = chainNMvtx.find(chainId);
      if (it == chainNMvtx.end() || it->second < m_colorMinMvtx)
      {
        return grey;
      }
    }
    return palette_color(chainId);
  };

  // ==================================================================
  // Si trajectories
  // ==================================================================
  struct SiTrajDraw
  {
    const Si_Trajectory* t = nullptr;
    Helix helix;
    int color = 0;
    bool lowPt = false;  // circle fit below setSiMinPt: drawn dotted (setDrawLowPtSi)
  };
  std::vector<SiTrajDraw> siTrajs;
  std::map<int, int> siChainStatus;  // chain id -> fit status of its trajectory
  unsigned int nSiLowPt = 0;
  if (m_siTraj)
  {
    for (unsigned int i = 0; i < m_siTraj->size(); ++i)
    {
      const Si_Trajectory* t = m_siTraj->get_trajectory(i);
      if (!t)
      {
        continue;
      }
      siChainStatus[t->get_chain_id()] = t->get_fit_status();
      if (t->size_points() < m_siMinPoints || (m_siOnlyGood && !t->is_fitted()))
      {
        continue;
      }
      // straight-line fits (MVTX only) have no pt: they are not affected by the pt cut
      const bool lowPt = t->get_fit_status() == Si_Trajectory::Ok && t->get_pt() < m_siMinPt;
      nSiLowPt += lowPt;
      if (lowPt && !m_drawLowPtSi)
      {
        continue;
      }
      siTrajs.push_back({t, helix_of(*t), chainColor(t->get_chain_id()), lowPt});
    }
  }

  // ==================================================================
  // TPC poly clusters / tracks / vertices (selection as Tpc_PolyClusterDisplay)
  // ==================================================================
  const bool extendSi = m_extension == ExtendSiIntoTpc || m_extension == ExtendBoth;
  const bool extendTpc = m_extension == ExtendTpcIntoSi || m_extension == ExtendBoth;
  const double zabs = std::max(std::fabs(m_zmin), std::fabs(m_zmax));

  // TPC colours: a TPC track matched to a Si trajectory (SITPC_TRACKS) takes the colour of the
  // Si chain, unmatched TPC tracks are grey; without SITPC_TRACKS (or setColorTpcByMatch(false))
  // the colour comes from the assembled track id as before.
  std::map<unsigned int, int> tpcMatchedChain;  // assembled track id -> Si chain id
  if (m_siTpcTracks)
  {
    for (unsigned int i = 0; i < m_siTpcTracks->size(); ++i)
    {
      if (const SiTpc_Track* t = m_siTpcTracks->get_track(i))
      {
        tpcMatchedChain[t->get_tpc_assembled_track_id()] = t->get_si_chain_id();
      }
    }
  }
  const bool colorByMatch = m_colorTpcByMatch && m_siTpcTracks;
  auto tpcColor = [&](unsigned int assembledId)
  {
    if (!colorByMatch)
    {
      return palette_color(static_cast<int>(assembledId));
    }
    const auto it = tpcMatchedChain.find(assembledId);
    return it != tpcMatchedChain.end() ? chainColor(it->second) : grey;
  };
  // crossing of an assembled track (TPC_CROSSING_DECISIONS; 0 without decisions)
  auto crossingOf = [&](unsigned int assembledId) -> short
  {
    const TpcCrossingDecision* d = m_decisions ? m_decisions->get_decision(assembledId) : nullptr;
    return d ? d->get_selected_crossing() : static_cast<short>(0);
  };

  struct TpcLine
  {
    unsigned int id = 0;  // assembled track id
    int color = 0;
    std::vector<XYZ> pts;  // over the clusters (solid)
    std::vector<XYZ> ext;  // inward extension through the silicon (dashed), if requested
  };
  std::map<unsigned int, std::vector<XYZ>> tpcClustersById;  // assembled track id -> clusters
  std::vector<TpcLine> tpcLines;
  struct PcaDraw
  {
    unsigned int id = 0;  // assembled track id
    int color = 0;
    XYZ p;
  };
  std::vector<PcaDraw> tpcPca;
  std::vector<XYZ> tpcCollision;
  if (m_tpcClusters)
  {
    std::map<unsigned int, const Tpc_PolyTrackVertex*> vtxById;
    const unsigned int nvtx = m_tpcVertices ? m_tpcVertices->size() : 0;
    for (unsigned int i = 0; i < nvtx; ++i)
    {
      if (const auto* v = m_tpcVertices->get_vertex(i))
      {
        vtxById[v->get_source_assembled_track_id()] = v;
      }
    }
    for (unsigned int i = 0; i < m_tpcClusters->size(); ++i)
    {
      const Tpc_PolyCluster* c = m_tpcClusters->get_cluster(i);
      if (!c || !c->isValid())
      {
        continue;
      }
      const unsigned int id = c->get_source_assembled_track_id();
      if (m_tpcVertices)
      {
        const auto it = vtxById.find(id);
        if (it == vtxById.end() || !track_vertex_z_selected(it->second, m_trackVertexZMin, m_trackVertexZMax))
        {
          continue;
        }
      }
      const XYZ p{c->get_centroid_x(), c->get_centroid_y(), c->get_centroid_z()};
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) || p.z < m_zmin || p.z > m_zmax)
      {
        continue;
      }
      tpcClustersById[id].push_back(p);
    }

    const unsigned int ntrk = m_tpcTracks ? m_tpcTracks->size() : 0;
    for (unsigned int i = 0; i < ntrk; ++i)
    {
      const Tpc_PolyTrack* trk = m_tpcTracks->get_track(i);
      if (!trk)
      {
        continue;
      }
      const auto it = tpcClustersById.find(trk->get_source_assembled_track_id());
      if (it == tpcClustersById.end() || it->second.empty())
      {
        continue;
      }
      double zlo = std::numeric_limits<double>::max(), zhi = -zlo;
      for (const auto& p : it->second)
      {
        zlo = std::min(zlo, p.z);
        zhi = std::max(zhi, p.z);
      }
      if (zlo == zhi)
      {
        zlo -= 0.1;
        zhi += 0.1;
      }
      const bool straight = m_useStraightLineTracks || std::fabs(trk->get_charge() * m_bz) < 1.0e-12;
      const double fwd = cluster_line_residual2(trk, it->second, m_bz, 1.0, straight);
      const double rev = cluster_line_residual2(trk, it->second, m_bz, -1.0, straight);
      const double arc = fwd <= rev ? 1.0 : -1.0;
      TpcLine line;
      line.id = trk->get_source_assembled_track_id();
      line.color = tpcColor(trk->get_source_assembled_track_id());
      line.pts = sample_poly_track(trk, zlo, zhi, m_xymax, m_bz, arc, straight);
      if (line.pts.size() < 2)
      {
        continue;
      }
      if (extendTpc)
      {
        // inner end of the track = the end of the solid line closer to the beam axis
        const bool innerIsFront = radius_of(line.pts.front()) <= radius_of(line.pts.back());
        const double zIn = innerIsFront ? line.pts.front().z : line.pts.back().z;
        const double zOut = innerIsFront ? line.pts.back().z : line.pts.front().z;
        line.ext = extend_poly_track_inward(trk, zIn, zOut, m_tpcExtrapR, zabs, m_bz, arc, straight);
      }
      tpcLines.push_back(std::move(line));
    }

    for (unsigned int i = 0; i < nvtx; ++i)
    {
      const auto* v = m_tpcVertices->get_vertex(i);
      if (!v || !v->get_pca_valid() || !track_vertex_z_selected(v, m_trackVertexZMin, m_trackVertexZMax))
      {
        continue;
      }
      const XYZ p{v->get_pca_x(), v->get_pca_y(), v->get_pca_z()};
      if (std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z))
      {
        tpcPca.push_back({v->get_source_assembled_track_id(), tpcColor(v->get_source_assembled_track_id()), p});
      }
    }
    const unsigned int ncoll = (m_tpcVertices && m_tpcVertices->get_collision_vertex_valid())
                                   ? m_tpcVertices->get_collision_vertex_count()
                                   : 0;
    for (unsigned int i = 0; i < ncoll; ++i)
    {
      const XYZ p{m_tpcVertices->get_collision_x(i), m_tpcVertices->get_collision_y(i), m_tpcVertices->get_collision_z(i)};
      if (std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z))
      {
        tpcCollision.push_back(p);
      }
    }
  }

  // TPC: reference frame of the beam-axis alignment, clusters / tracks / vertices are drawn as
  // reconstructed (SiTpcBeamAlignment::tpcToTpcFrame is the identity).  Only the Si is moved.

  if ((m_requireSiTraj && siTrajs.empty()) || (m_requireTpcTrack && tpcLines.empty()))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // ==================================================================
  // Si hits: global positions, split by chain
  // ==================================================================
  // chains drawn: the chain selected for every TPC track (setDrawCandidateChains(true): all
  // candidate chains; they may share hits, a shared hit is drawn in each chain's colour)
  std::vector<const SiHitChain*> drawnChains;
  std::vector<std::pair<const SiHitChain*, double>> chainSeedZ0;  // chain, z0 of its TPC seed
  if (m_siSeeds)
  {
    std::map<int, double> seedZ0;  // TPC track index -> pca z of the TPC seed
    for (const auto& ts : m_siSeeds->tpc_seeds)
    {
      seedZ0[ts.tpc_track_index] = ts.z0;
    }
    for (const auto& chain : m_siSeeds->chains)
    {
      if (!chain.selected && !m_drawCandidateChains)
      {
        continue;
      }
      drawnChains.push_back(&chain);
      const auto it = seedZ0.find(chain.tpc_track_index);
      if (it != seedZ0.end() && std::isfinite(it->second))
      {
        chainSeedZ0.emplace_back(&chain, it->second);
      }
    }
  }
  std::map<int, short> chainCrossing;  // chain id -> crossing of its TPC trajectory
  for (const auto* cp : drawnChains)
  {
    chainCrossing[cp->id] = cp->tpc_crossing;
  }
  std::vector<XYZ> siFreeHits;
  std::map<int, std::vector<XYZ>> siChainHits;  // chain id -> hits
  std::vector<char> onChain;
  if (m_siSeeds)
  {
    onChain.assign(m_siSeeds->hits.size(), 0);
    for (const auto* cp : drawnChains)
    {
      const auto& chain = *cp;
      for (int hid : chain.hit_ids)
      {
        if (hid >= 0 && hid < static_cast<int>(m_siSeeds->hits.size()))
        {
          onChain[hid] = 1;
        }
      }
    }
    for (std::size_t i = 0; i < m_siSeeds->hits.size(); ++i)
    {
      const auto& h = m_siSeeds->hits[i];
      if (h.layer < 0 || h.layer >= SiDetectorFrame::kNLayers)
      {
        continue;
      }
      const auto g = m_alignment.siToTpcFrame(m_frame, h.layer, h.phi, h.z);
      if (!onChain[i])
      {
        siFreeHits.push_back({g.x, g.y, g.z});
      }
    }
    for (const auto* cp : drawnChains)
    {
      const auto& chain = *cp;
      for (int hid : chain.hit_ids)
      {
        if (hid < 0 || hid >= static_cast<int>(m_siSeeds->hits.size()))
        {
          continue;
        }
        const auto& h = m_siSeeds->hits[hid];
        const auto g = m_alignment.siToTpcFrame(m_frame, h.layer, h.phi, h.z);
        siChainHits[chain.id].push_back({g.x, g.y, g.z});
      }
    }
  }
  // chains without a usable trajectory: no trajectory at all, or a failed fit
  std::map<int, std::vector<XYZ>> siUnfittedHits;  // chain id -> hits
  unsigned int nSiUnfitted = 0;
  for (const auto& [id, pts] : siChainHits)
  {
    const auto it = siChainStatus.find(id);
    const bool fitted = it != siChainStatus.end() &&
                        (it->second == Si_Trajectory::Ok || it->second == Si_Trajectory::StraightLine);
    if (fitted)
    {
      continue;
    }
    ++nSiUnfitted;
    siUnfittedHits[id] = pts;
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": chain " << id << " not fitted ("
                << (it == siChainStatus.end() ? std::string("no trajectory in ") + m_siTrajNodeName
                                              : "fit status " + std::to_string(it->second))
                << "), " << pts.size() << " hits" << std::endl;
    }
  }
  std::string siTpcSummary;
  if (!m_siTpcTracks)
  {
    siTpcSummary = "no " + m_siTpcTrackNodeName + " node (SiTpcTrackMatcher not run before the display?)";
  }
  else
  {
    siTpcSummary = std::format("{} Si+TPC tracks ({}), TPC coloured by {}", m_siTpcTracks->size(),
                               m_actsTracks ? std::format("ACTS fit drawn: {} tracks, failed fits dotted", m_actsTracks->size())
                                            : (m_drawSiTpcTracks ? "drawn" : "lines off: setDrawSiTpcTracks(true)"),
                               colorByMatch ? "matched Si chain, unmatched grey" : "assembled id");
  }
  const std::string siSummary = std::format("{} chains: {} drawn, {} below {:.2f} GeV{}, {} not fitted (#diamond)",
                                            siChainHits.size(), siTrajs.size() - (m_drawLowPtSi ? nSiLowPt : 0),
                                            nSiLowPt, m_siMinPt, m_drawLowPtSi ? " (dotted)" : " (hidden)", nSiUnfitted);

  // ==================================================================
  // output directory
  // ==================================================================
  TDirectory* eventsTop = m_outfile->GetDirectory("events");
  if (!eventsTop)
  {
    eventsTop = m_outfile->mkdir("events");
  }
  TDirectory* eventDir = eventsTop->mkdir(std::format("event_{:06}", m_evt).c_str());
  if (!eventDir)
  {
    std::cerr << Name() << " - failed to create event directory" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  eventDir->cd();
  const std::string tag = std::format("evt{:06}", m_evt);
  const std::string title = std::format("event {}", m_evt);
  std::vector<TObject*> legendStyles;  // legend entries do not own their style objects

  // ==================================================================
  // 1. DETECTOR coordinates: (Vz, Uphi, layer), one canvas
  // ==================================================================
  if (m_drawDetector && m_siSeeds && !m_siSeeds->hits.empty())
  {
    // Si trajectories in detector coordinates: helix crossing with each nominal layer
    // (circle of radius R_layer around the detector centre), converted global -> detector.
    struct DetTraj
    {
      std::vector<P3> layers;  // outermost first, like the chains
      P3 vertex{0, 0, kVertexLayer};
      bool hasVertex = false;
      int color = 0;
    };
    std::vector<DetTraj> detTrajs;
    for (const auto& d : siTrajs)
    {
      DetTraj dt;
      dt.color = d.color;
      for (int l = SiDetectorFrame::kNLayers - 1; l >= 0; --l)
      {
        // the trajectory is in the TPC frame: undo the move of the clamshell half onto the TPC
        // beam axis and the detector shift to get the crossing with the nominal layer
        auto toDet = [&](const XYZ& b)
        {
          const auto dp = m_frame.toDetectorFromGlobal(l, m_alignment.siTpcFrameToGlobal(m_frame, {b.x, b.y, b.z}));
          return XYZ{dp.x, dp.y, dp.z};
        };
        XYZ g;
        if (!helix_at_radius_mapped(d.helix, toDet, SiDetectorFrame::layerRadius(l), g))
        {
          continue;
        }
        const XYZ det = toDet(g);
        dt.layers.push_back({det.z, wrap2pi(std::atan2(det.y, det.x) - m_frame.uphiRotation()), static_cast<double>(l)});
      }
      if (!dt.layers.empty() && d.helix.ok)
      {
        dt.vertex = {d.helix.z0 - center[2], dt.layers.back().phi, kVertexLayer};
        dt.hasVertex = true;
      }
      detTrajs.push_back(std::move(dt));
    }

    // TPC seed pca z (detector frame z = TPC frame z - detector centre z) at layer -1
    const auto [zlo, zhi] = std::minmax_element(m_siSeeds->hits.begin(), m_siSeeds->hits.end(),
                                                [](const auto& a, const auto& b) { return a.z < b.z; });
    double zlow = zlo->z, zhigh = zhi->z;
    for (const auto& [cp, z0] : chainSeedZ0)
    {
      zlow = std::min(zlow, z0 - center[2]);
      zhigh = std::max(zhigh, z0 - center[2]);
    }
    for (const auto& dt : detTrajs)
    {
      if (dt.hasVertex)
      {
        zlow = std::min(zlow, dt.vertex.z);
        zhigh = std::max(zhigh, dt.vertex.z);
      }
    }
    const double span = std::max(0.01, zhigh - zlow);
    const double zmin = zlow - std::max(0.02 * span, 0.01);
    const double zmax = zhigh + std::max(0.02 * span, 0.01);

    // Uphi axis: one full turn [b, b+2pi), b snapped to a phi-bin edge.
    const double dphi = kTwoPi / m_nphi;
    const double b = std::round(wrap2pi(m_boundary) / dphi) * dphi;

    // One histogram with all hits (chain hits and, if enabled, the others).
    TH3F* hHits = make_det_hist("h3_" + tag + "_det_hits",
                                title + " Si hits, solid TPC-driven chains / dashed Si trajectories;"
                                        "nominal V_{z} [cm];intrinsic U_{#phi} [rad];layer",
                                m_nz, zmin, zmax, m_nphi, b, b + kTwoPi);
    int nAssoc = 0, nFree = 0;
    for (std::size_t i = 0; i < m_siSeeds->hits.size(); ++i)
    {
      const auto& hit = m_siSeeds->hits[i];
      if (!onChain[i] && !m_drawUnassociated)
      {
        continue;
      }
      hHits->Fill(hit.z, toU(hit.phi, b), hit.layer);
      ++(onChain[i] ? nAssoc : nFree);
    }

    auto* c = make_det_canvas("c3_" + tag + "_det_hits_chains_traj", hHits->GetTitle());
    hHits->Draw("BOX2Z");
    auto drawPath = [&](const std::vector<P3>& pts, int color, int style)
    {
      for (const auto& piece : splitAtSeam(pts, b))
      {
        make_det_line(piece, color, style)->Draw("same");
      }
    };
    for (const auto* cp : drawnChains)
    {
      const auto& chain = *cp;
      std::vector<P3> pts;
      for (size_t i = 0; i < chain.z.size(); ++i)
      {
        pts.push_back({chain.z[i], chain.phi[i], static_cast<double>(chain.layers[i])});
      }
      // candidate (not selected) chains thin dotted
      if (chain.selected)
      {
        drawPath(pts, chainColor(chain.id), 1);
      }
      else
      {
        for (const auto& piece : splitAtSeam(pts, b))
        {
          auto* l = make_det_line(piece, chainColor(chain.id), 3);
          l->SetLineWidth(1);
          l->Draw("same");
        }
      }
    }
    for (const auto& dt : detTrajs)
    {
      drawPath(dt.layers, dt.color, 2);
      if (dt.hasVertex)
      {
        drawPath({dt.layers.back(), dt.vertex}, dt.color, 3);
      }
    }
    // vertex plane: pca z of the TPC seed of every drawn chain (chain colour, at the chain's
    // innermost Uphi), joined to the innermost chain point
    for (const auto& [cp, z0] : chainSeedZ0)
    {
      if (cp->z.empty())
      {
        continue;
      }
      const P3 vtx{z0 - center[2], cp->phi.back(), kVertexLayer};
      drawPath({{cp->z.back(), cp->phi.back(), static_cast<double>(cp->layers.back())}, vtx}, chainColor(cp->id), 3);
      auto* m = owned(new TPolyMarker3D(1));
      m->SetPoint(0, vtx.z, toU(vtx.phi, b), kVertexLayer);
      m->SetMarkerStyle(23);
      m->SetMarkerColor(chainColor(cp->id));
      m->SetMarkerSize(1.6);
      m->Draw("same");
    }
    auto* legend = owned(new TLegend(0.12, 0.76, 0.68, 0.92));
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    auto* chainStyle = new TLine();
    chainStyle->SetLineColor(kRed + 1);
    chainStyle->SetLineWidth(3);
    auto* trajStyle = new TLine();
    trajStyle->SetLineColor(kRed + 1);
    trajStyle->SetLineStyle(2);
    trajStyle->SetLineWidth(2);
    auto* vtxStyle = new TMarker(0, 0, 23);
    vtxStyle->SetMarkerColor(kRed + 1);
    vtxStyle->SetMarkerSize(1.6);
    legendStyles.insert(legendStyles.end(), {chainStyle, trajStyle, vtxStyle});
    unsigned int nSelected = 0;
    for (const auto* cp : drawnChains)
    {
      nSelected += cp->selected;
    }
    legend->AddEntry(chainStyle, std::format("Hit chains: {} selected{} ({} chain hits, {} other hits)", nSelected,
                                             m_drawCandidateChains ? std::format(", {} candidates dotted", drawnChains.size() - nSelected) : std::string(),
                                             nAssoc, nFree).c_str(), "l");
    legend->AddEntry(trajStyle, std::format("Si trajectories: {} (dotted: to z0 at layer -1)", detTrajs.size()).c_str(), "l");
    legend->AddEntry(vtxStyle, "TPC seed pca z (layer -1, chain colour)", "p");
    legend->Draw();
    write_and_delete(c);
  }

  // ==================================================================
  // 1b. PHYSICAL coordinates of the hits: (phi, eta, r), one canvas
  //     same content as det_hits_chains_traj, hits at their TPC-frame position
  // ==================================================================
  if (m_drawPhysHits && m_siSeeds && !m_siSeeds->hits.empty())
  {
    const double ox = beamX(), oy = beamY();
    const bool zrefFromVertex = m_etaRefFromVertex && !tpcCollision.empty();
    const double zref = zrefFromVertex ? tpcCollision.front().z : m_etaRefZ;
    const double rmaxView = kSiZoomR;

    // TPC-frame point -> {eta, phi, r} (P3 layout {z, phi, layer}); eta NaN at r = 0
    auto toPhys = [&](const XYZ& g) -> P3
    {
      const double dx = g.x - ox, dy = g.y - oy;
      const double r = std::hypot(dx, dy);
      const double eta = r > 0 ? std::asinh((g.z - zref) / r) : std::numeric_limits<double>::quiet_NaN();
      return {eta, std::atan2(dy, dx), r};
    };
    auto hitPhys = [&](int hid) -> P3
    {
      const auto& h = m_siSeeds->hits[hid];
      const auto g = m_alignment.siToTpcFrame(m_frame, h.layer, h.phi, h.z);
      return toPhys({g.x, g.y, g.z});
    };
    auto validHit = [&](int hid)
    {
      if (hid < 0 || hid >= static_cast<int>(m_siSeeds->hits.size()))
      {
        return false;
      }
      const int l = m_siSeeds->hits[hid].layer;
      return l >= 0 && l < SiDetectorFrame::kNLayers;
    };

    // hits drawn and the eta range from them
    std::vector<P3> hitPts;
    int nAssoc = 0, nFree = 0;
    double etaLo = std::numeric_limits<double>::max(), etaHi = -etaLo;
    for (std::size_t i = 0; i < m_siSeeds->hits.size(); ++i)
    {
      if (!validHit(static_cast<int>(i)) || (!onChain[i] && !m_drawUnassociated))
      {
        continue;
      }
      const P3 p = hitPhys(static_cast<int>(i));
      if (!std::isfinite(p.z))
      {
        continue;
      }
      hitPts.push_back(p);
      etaLo = std::min(etaLo, p.z);
      etaHi = std::max(etaHi, p.z);
      ++(onChain[i] ? nAssoc : nFree);
    }
    if (hitPts.empty())
    {
      etaLo = -1.0;
      etaHi = 1.0;
    }
    const double etaSpan = std::max(0.05, etaHi - etaLo);
    const double etaMin = etaLo - std::max(0.03 * etaSpan, 0.05);
    const double etaMax = etaHi + std::max(0.03 * etaSpan, 0.05);

    // phi axis: one full turn [bp, bp+2pi), bp snapped to a phi-bin edge
    const double dphi = kTwoPi / m_nphi;
    const double bp = std::round(m_physPhiBoundary / dphi) * dphi;

    TH3F* hPhys = owned(new TH3F(("h3_" + tag + "_phys_phi_eta_r_hits").c_str(),
                                 std::format("{} Si hits, solid TPC-driven chains / dashed Si trajectories "
                                             "(#eta from z_{{ref}} = {:.2f} cm{});"
                                             "#phi [rad];#eta;r [cm]",
                                             title, zref, zrefFromVertex ? ", collision vertex" : "")
                                     .c_str(),
                                 m_nphi, bp, bp + kTwoPi, m_neta, etaMin, etaMax, m_nr, 0.0, rmaxView));
    hPhys->SetDirectory(nullptr);
    hPhys->SetStats(false);
    for (const auto& p : hitPts)
    {
      hPhys->Fill(toU(p.phi, bp), p.z, p.layer);
    }

    auto* c = make_det_canvas("c3_" + tag + "_phys_phi_eta_r_hits_chains_traj", hPhys->GetTitle());
    hPhys->Draw("BOX2Z");

    // lines: keep only points with finite eta inside the axis range (3D lines are not clipped
    // by the frame), split into connected pieces, then at the phi seam
    auto drawPhysPath = [&](const std::vector<P3>& pts, int color, int style, int width)
    {
      std::vector<std::vector<P3>> runs;
      bool inside = false;
      for (const auto& p : pts)
      {
        const bool in = std::isfinite(p.z) && p.z >= etaMin && p.z <= etaMax && p.layer <= rmaxView;
        if (in)
        {
          if (!inside)
          {
            runs.emplace_back();
          }
          runs.back().push_back(p);
        }
        inside = in;
      }
      for (const auto& run : runs)
      {
        for (const auto& piece : splitAtSeam(run, bp))
        {
          make_phys_line(piece, color, style, width)->Draw("same");
        }
      }
    };

    // chains through their hits (hit_ids, outermost first like the chain arrays)
    for (const auto* cp : drawnChains)
    {
      std::vector<P3> pts;
      for (int hid : cp->hit_ids)
      {
        if (validHit(hid))
        {
          pts.push_back(hitPhys(hid));
        }
      }
      if (cp->selected)
      {
        drawPhysPath(pts, chainColor(cp->id), 1, 3);
      }
      else
      {
        drawPhysPath(pts, chainColor(cp->id), 3, 1);
      }
    }

    // Si trajectories: helix (TPC frame) sampled from the pca to the edge of the silicon,
    // drawn from setPhysViewMinRadius outward; low-pt fits dotted
    unsigned int nTrajDrawn = 0;
    for (const auto& d : siTrajs)
    {
      std::vector<P3> pts;
      for (const auto& g : sample_helix(d.helix, ox, oy, rmaxView, zabs))
      {
        const P3 p = toPhys(g);
        if (p.layer >= m_physMinR)
        {
          pts.push_back(p);
        }
      }
      if (pts.size() < 2)
      {
        continue;
      }
      ++nTrajDrawn;
      drawPhysPath(pts, d.color, d.lowPt ? 3 : 2, 2);
    }

    // chains without a usable fit: open black diamonds on their hits
    for (const auto* cp : drawnChains)
    {
      if (siUnfittedHits.find(cp->id) == siUnfittedHits.end())
      {
        continue;
      }
      std::vector<P3> pts;
      for (int hid : cp->hit_ids)
      {
        if (validHit(hid))
        {
          const P3 p = hitPhys(hid);
          if (std::isfinite(p.z))
          {
            pts.push_back(p);
          }
        }
      }
      if (pts.empty())
      {
        continue;
      }
      auto* m = owned(new TPolyMarker3D(static_cast<int>(pts.size())));
      for (std::size_t i = 0; i < pts.size(); ++i)
      {
        m->SetPoint(static_cast<int>(i), toU(pts[i].phi, bp), pts[i].z, pts[i].layer);
      }
      m->SetMarkerStyle(27);
      m->SetMarkerColor(kBlack);
      m->SetMarkerSize(1.4);
      m->Draw("same");
    }

    auto* legend = owned(new TLegend(0.12, 0.78, 0.68, 0.92));
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    auto* chainStyle = new TLine();
    chainStyle->SetLineColor(kRed + 1);
    chainStyle->SetLineWidth(3);
    auto* trajStyle = new TLine();
    trajStyle->SetLineColor(kRed + 1);
    trajStyle->SetLineStyle(2);
    trajStyle->SetLineWidth(2);
    auto* unfitStyle = new TMarker(0, 0, 27);
    unfitStyle->SetMarkerColor(kBlack);
    unfitStyle->SetMarkerSize(1.4);
    legendStyles.insert(legendStyles.end(), {chainStyle, trajStyle, unfitStyle});
    unsigned int nSelected = 0;
    for (const auto* cp : drawnChains)
    {
      nSelected += cp->selected;
    }
    legend->AddEntry(chainStyle, std::format("Hit chains: {} selected{} ({} chain hits, {} other hits)", nSelected,
                                             m_drawCandidateChains ? std::format(", {} candidates dotted", drawnChains.size() - nSelected) : std::string(),
                                             nAssoc, nFree).c_str(), "l");
    legend->AddEntry(trajStyle, std::format("Si trajectories: {} (r #geq {:.1f} cm)", nTrajDrawn, m_physMinR).c_str(), "l");
    legend->AddEntry(unfitStyle, std::format("chains not fitted: {}", siUnfittedHits.size()).c_str(), "p");
    legend->Draw();
    write_and_delete(c);
  }

  // ==================================================================
  // 2. PHYSICS space
  // ==================================================================
  // Si trajectory curves: inside the silicon (solid) and, if requested, the extension out
  // through the TPC (dashed), starting where the solid part ends.
  std::vector<std::vector<XYZ>> siCurveIn, siCurveExt;
  for (const auto& d : siTrajs)
  {
    // one sampling, split where the trajectory leaves the silicon; the extension starts at
    // the last silicon point so the dashed part joins the solid one
    // solid inside the silicon (up to the silicon radius, or less if the extrapolation radius
    // is smaller), dashed beyond; extendSi: the trajectory ends exactly at m_siExtrapR
    const double rEnd = extendSi ? m_siExtrapR : kSiZoomR;
    const double rSolid = std::min(kSiZoomR, rEnd);
    const auto full = sample_helix(d.helix, 0.0, 0.0, rEnd, zabs);
    std::size_t nIn = 0;
    while (nIn < full.size() && radius_of(full[nIn]) <= rSolid)
    {
      ++nIn;
    }
    siCurveIn.emplace_back(full.begin(), full.begin() + static_cast<long>(nIn));
    std::vector<XYZ> ext;
    if (extendSi && nIn > 0 && nIn < full.size())
    {
      ext.assign(full.begin() + static_cast<long>(nIn - 1), full.end());
    }
    siCurveExt.push_back(std::move(ext));
  }

  // In the silicon zoom only the part of the TPC tracks inside the zoom box is drawn.
  const double zoomZ = 30.0;
  auto tpcCurves = [&](const TpcLine& l, const bool zoom, const bool extension)
  {
    const auto& pts = extension ? l.ext : l.pts;
    return zoom ? clip_curve(pts, kSiZoomR, zoomZ) : std::vector<std::vector<XYZ>>{pts};
  };
  // Si+TPC fitted tracks (optional): from the pca outward through the whole tracker.
  // With the ACTS track map: the ACTS tracks instead (solid), Si+TPC tracks without ACTS fit dotted.
  std::vector<std::vector<XYZ>> siTpcCurves;
  std::vector<unsigned int> siTpcIds;  // assembled track id of every curve
  std::vector<int> siTpcStyles;        // line style of every curve
  std::vector<int> siTpcColors;        // line colour of every curve
  std::set<unsigned int> actsAssembled;  // assembled ids with an ACTS track
  if (m_actsTracks)
  {
    // TPC cluster key -> assembled track id (the TPC seed of an ACTS track holds its poly cluster keys)
    std::map<TrkrDefs::cluskey, unsigned int> assembledOfKey;
    if (m_tpcClusters)
    {
      for (unsigned int i = 0; i < m_tpcClusters->size(); ++i)
      {
        const Tpc_PolyCluster* c = m_tpcClusters->get_cluster(i);
        if (c && c->get_trkr_cluster_key() != TrkrDefs::CLUSKEYMAX)
        {
          assembledOfKey[c->get_trkr_cluster_key()] = c->get_source_assembled_track_id();
        }
      }
    }
    // assembled TPC track id -> per-track TPC z offset of its Si+TPC track
    std::map<unsigned int, double> zOffsetOfAssembled;
    if (m_siTpcTracks)
    {
      for (unsigned int i = 0; i < m_siTpcTracks->size(); ++i)
      {
        const SiTpc_Track* t = m_siTpcTracks->get_track(i);
        if (t && t->isValid())
        {
          zOffsetOfAssembled[t->get_tpc_assembled_track_id()] = t->get_tpc_z_offset();
        }
      }
    }
    for (const auto& [key, trk] : *m_actsTracks)
    {
      if (!trk)
      {
        continue;
      }
      // assembled id: majority of the TPC clusters of the track
      std::map<unsigned int, int> votes;
      if (const TrackSeed* tpcSeed = trk->get_tpc_seed())
      {
        for (auto it = tpcSeed->begin_cluster_keys(); it != tpcSeed->end_cluster_keys(); ++it)
        {
          const auto f = assembledOfKey.find(*it);
          if (f != assembledOfKey.end())
          {
            ++votes[f->second];
          }
        }
      }
      if (votes.empty())
      {
        continue;
      }
      const unsigned int assembledId = std::max_element(votes.begin(), votes.end(), [](const auto& a, const auto& b)
                                                        { return a.second < b.second; })
                                           ->first;
      // line through the fitted states, in the display frame.  Si states: + frame shift; TPC states also
      // minus the TPC-only shifts of the fitter (per-track TPC z offset, extra TPC shift).
      double tpcZOffset = 0;
      if (m_actsTpcZOffset)
      {
        const auto f = zOffsetOfAssembled.find(assembledId);
        if (f != zOffsetOfAssembled.end())
        {
          tpcZOffset = f->second;
        }
      }
      struct State
      {
        double s, x, y, z, px, py, pz;
      };
      std::vector<State> states;
      states.push_back({0.0, trk->get_x() + m_actsShift[0], trk->get_y() + m_actsShift[1], trk->get_z() + m_actsShift[2],
                        trk->get_px(), trk->get_py(), trk->get_pz()});
      for (auto it = trk->begin_states(); it != trk->end_states(); ++it)
      {
        const SvtxTrackState* st = it->second;
        if (!st || st->get_pathlength() <= 0)
        {
          continue;  // the state at path length 0 holds no position (perigee = track position)
        }
        const bool tpc = st->get_cluskey() != std::numeric_limits<TrkrDefs::cluskey>::max() &&
                         TrkrDefs::getTrkrId(st->get_cluskey()) == TrkrDefs::tpcId;
        const double dx = m_actsShift[0] - (tpc ? m_actsTpcExtra[0] : 0.0);
        const double dy = m_actsShift[1] - (tpc ? m_actsTpcExtra[1] : 0.0);
        const double dz = m_actsShift[2] + (tpc ? tpcZOffset - m_actsTpcExtra[2] : 0.0);
        states.push_back({st->get_pathlength(), st->get_x() + dx, st->get_y() + dy, st->get_z() + dz, st->get_px(), st->get_py(), st->get_pz()});
      }
      std::sort(states.begin(), states.end(), [](const State& a, const State& b) { return a.s < b.s; });
      // between two states: the circle tangent to the fitted direction at the first state through the
      // second one (transverse plane), z linear in the transverse arc length
      std::vector<XYZ> pts;
      bool outside = false;
      auto add = [&](const XYZ& v)
      {
        if (std::hypot(v.x, v.y) > m_xymax || std::fabs(v.z) > zabs)
        {
          outside = true;
          return;
        }
        pts.push_back(v);
      };
      for (std::size_t k = 0; k + 1 < states.size() && !outside; ++k)
      {
        const State& a = states[k];
        const State& b = states[k + 1];
        const double pt = std::hypot(a.px, a.py);
        const double cx = b.x - a.x;
        const double cy = b.y - a.y;
        const double chord = std::hypot(cx, cy);
        if (!(pt > 0) || chord < 1e-6)
        {
          add({a.x, a.y, a.z});
          continue;
        }
        const double tx = a.px / pt;
        const double ty = a.py / pt;
        const double kappa = 2.0 * (tx * cy - ty * cx) / (chord * chord);  // signed curvature [1/cm]
        const double half = 0.5 * kappa * chord;
        double len = chord;  // transverse arc length
        bool circle = std::fabs(kappa) * chord > 1e-6;
        if (circle)
        {
          if (std::fabs(half) < 1.0 && (tx * cx + ty * cy) > 0)
          {
            len = 2.0 * std::asin(half) / kappa;
          }
          else
          {
            circle = false;  // next state behind the fitted direction: straight connection
          }
        }
        const double psi0 = std::atan2(ty, tx);
        for (double sv = 0; sv < len && !outside; sv += kHelixStep)
        {
          const double z = a.z + (b.z - a.z) * sv / len;
          if (circle)
          {
            add({a.x + (std::sin(psi0 + kappa * sv) - std::sin(psi0)) / kappa, a.y - (std::cos(psi0 + kappa * sv) - std::cos(psi0)) / kappa, z});
          }
          else
          {
            add({a.x + cx * sv / len, a.y + cy * sv / len, z});
          }
        }
      }
      if (!outside && !states.empty())
      {
        const State& b = states.back();
        add({b.x, b.y, b.z});
      }
      siTpcCurves.push_back(std::move(pts));
      siTpcIds.push_back(assembledId);
      siTpcStyles.push_back(1);
      siTpcColors.push_back(m_actsTrackColor);
      actsAssembled.insert(assembledId);
    }
  }
  if (m_siTpcTracks && (m_actsTracks ? m_drawSiTpcForFailedActs : m_drawSiTpcTracks))
  {
    for (unsigned int i = 0; i < m_siTpcTracks->size(); ++i)
    {
      const SiTpc_Track* t = m_siTpcTracks->get_track(i);
      if (!t || !t->isValid())
      {
        continue;
      }
      if (m_actsTracks && actsAssembled.contains(t->get_tpc_assembled_track_id()))
      {
        continue;  // drawn as ACTS track
      }
      const auto hx = SiTpcHelixFit::fromTrack(*t);
      const double smax = std::min(kPi * hx.R, 4.0 * m_xymax + 100.0);
      std::vector<XYZ> pts;
      for (double sv = 0; sv <= smax; sv += kHelixStep)
      {
        const auto p = SiTpcHelixFit::at(hx, sv);
        if (std::hypot(p.x, p.y) > m_xymax || std::fabs(p.z) > zabs)
        {
          break;
        }
        pts.push_back({p.x, p.y, p.z});
      }
      siTpcCurves.push_back(std::move(pts));
      siTpcIds.push_back(t->get_tpc_assembled_track_id());
      siTpcStyles.push_back(m_actsTracks ? 3 : 1);  // dotted: no ACTS fit for this Si+TPC track
      siTpcColors.push_back(m_siTpcTrackColor);
    }
  }
  auto siTpcSegments = [&](const std::vector<XYZ>& pts, const bool zoom)
  {
    return zoom ? clip_curve(pts, kSiZoomR, zoomZ) : std::vector<std::vector<XYZ>>{pts};
  };

  const int siHitStyle = 21;  // squares
  const int tpcClusterStyle = 20;  // circles

  auto chainXing = [&](int chainId)
  {
    const auto it = chainCrossing.find(chainId);
    return it != chainCrossing.end() ? it->second : static_cast<short>(0);
  };
  // xing: draw only what belongs to this crossing; nullopt: all crossings
  // freeSiHits: also the Si hits not in any chain (grey)
  auto drawEverything3D = [&](const bool zoom, const std::optional<short> xing, const bool collisions, const bool freeSiHits)
  {
    auto keep = [&](short crossing) { return !xing || crossing == *xing; };
    auto shown = [&](int color) { return m_drawGreyTracks || color != grey; };
    for (const auto& l : tpcLines)
    {
      if (!keep(crossingOf(l.id)) || !shown(l.color))
      {
        continue;
      }
      for (const auto& seg : tpcCurves(l, zoom, false))
      {
        draw_line_zxy(seg, l.color, 1, 3);
      }
      for (const auto& seg : tpcCurves(l, zoom, true))
      {
        draw_line_zxy(seg, l.color, 2, 2);
      }
    }
    for (const auto& [id, pts] : tpcClustersById)
    {
      if (keep(crossingOf(id)) && shown(tpcColor(id)))
      {
        draw_markers_zxy(pts, tpcColor(id), tpcClusterStyle, zoom ? 0.8 : 1.0);
      }
    }
    for (const auto& v : tpcPca)
    {
      if (keep(crossingOf(v.id)) && shown(v.color))
      {
        draw_markers_zxy({v.p}, v.color, 20, 1.2);
      }
    }
    if (collisions)
    {
      draw_markers_zxy(tpcCollision, kBlack, 29, 2.0);
    }
    if (freeSiHits)
    {
      draw_markers_zxy(siFreeHits, kGray + 1, siHitStyle, zoom ? 0.4 : 0.3);
    }
    for (const auto& [id, pts] : siChainHits)
    {
      if (keep(chainXing(id)) && shown(chainColor(id)))
      {
        draw_markers_zxy(pts, chainColor(id), siHitStyle, zoom ? 0.8 : 0.6);
      }
    }
    for (const auto& [id, pts] : siUnfittedHits)
    {
      if (keep(chainXing(id)))
      {
        draw_markers_zxy(pts, kBlack, 27, zoom ? 1.4 : 1.0);
      }
    }
    for (std::size_t i = 0; i < siTpcCurves.size(); ++i)
    {
      if (!keep(crossingOf(siTpcIds[i])) || !shown(tpcColor(siTpcIds[i])))
      {
        continue;
      }
      for (const auto& seg : siTpcSegments(siTpcCurves[i], zoom))
      {
        draw_line_zxy(seg, siTpcColors[i], siTpcStyles[i], 2);
      }
    }
    for (std::size_t i = 0; i < siTrajs.size(); ++i)
    {
      if (!keep(chainXing(siTrajs[i].t->get_chain_id())) || !shown(siTrajs[i].color))
      {
        continue;
      }
      const bool low = siTrajs[i].lowPt;
      if (!zoom)
      {
        draw_line_zxy(siCurveExt[i], siTrajs[i].color, low ? 3 : 2, low ? 1 : 2);
      }
      draw_line_zxy(siCurveIn[i], siTrajs[i].color, low ? 3 : 1, low ? 2 : 3);
    }
  };

  if (m_drawPhysics)
  {
    // crossings present in the event: of the drawn TPC tracks, Si trajectories, Si+TPC tracks
    std::set<short> crossings;
    for (const auto& l : tpcLines)
    {
      crossings.insert(crossingOf(l.id));
    }
    for (const auto& d : siTrajs)
    {
      crossings.insert(chainXing(d.t->get_chain_id()));
    }
    for (const unsigned int id : siTpcIds)
    {
      crossings.insert(crossingOf(id));
    }
    const bool singleCrossing = crossings.size() <= 1;

    // full detector, one canvas per crossing
    for (const short xing : crossings)
    {
      unsigned int nTpc = 0, nSi = 0, nSiTpc = 0;
      for (const auto& l : tpcLines)
      {
        nTpc += crossingOf(l.id) == xing;
      }
      for (const auto& d : siTrajs)
      {
        nSi += chainXing(d.t->get_chain_id()) == xing;
      }
      for (const unsigned int id : siTpcIds)
      {
        nSiTpc += crossingOf(id) == xing;
      }
      // collision z of this crossing: silicon_vertex_z of the tier 0 / 1 decisions of its tracks
      // (TpcCrossingResolver: one collision per crossing, so normally a single value)
      std::map<long, unsigned int> vertexVotes;  // z in units of 0.01 cm -> tracks
      {
        std::set<unsigned int> ids;
        for (const auto& l : tpcLines)
        {
          ids.insert(l.id);
        }
        ids.insert(siTpcIds.begin(), siTpcIds.end());
        for (const unsigned int id : ids)
        {
          const TpcCrossingDecision* d = m_decisions ? m_decisions->get_decision(id) : nullptr;
          if (d && d->get_selected_crossing() == xing && d->get_selected_tier() <= 1 && std::isfinite(d->get_silicon_vertex_z()))
          {
            ++vertexVotes[std::lround(100.0 * d->get_silicon_vertex_z())];
          }
        }
      }
      std::string vtxText = "no vertex";
      if (!vertexVotes.empty())
      {
        vtxText.clear();
        for (const auto& [zc, n] : vertexVotes)
        {
          vtxText += std::format("{}z_{{vtx}} = {:.2f} cm ({} trk)", vtxText.empty() ? "" : ", ", 0.01 * static_cast<double>(zc), n);
        }
      }
      const std::string xtag = xing < 0 ? std::format("m{}", -xing) : std::format("{}", xing);
      const std::string xtitle = std::format("{} crossing {}, {}: {} TPC tracks, {} Si trajectories, {} Si+TPC tracks",
                                             title, xing, vtxText, nTpc, nSi, nSiTpc);
      auto* c = new TCanvas(("c3_" + tag + "_phys_z_x_y_crossing" + xtag).c_str(), xtitle.c_str(), 1200, 900);
      auto* h3 = owned(new TH3D(("h3_" + tag + "_phys_z_x_y_crossing" + xtag).c_str(),
                                (xtitle + " - TPC (circles) + Si (squares);z [cm];x [cm];y [cm]").c_str(),
                                204, m_zmin, m_zmax, 170, -m_xymax, m_xymax, 170, -m_xymax, m_xymax));
      h3->SetDirectory(nullptr);
      h3->SetStats(false);
      h3->Draw();
      drawEverything3D(false, xing, singleCrossing, false);
      write_and_delete(c);
    }
    // full detector, all crossings in one canvas (same content as the per-crossing canvases)
    {
      const std::string xtitle = std::format("{} all {} crossings: {} TPC tracks, {} Si trajectories, {} Si+TPC tracks",
                                             title, crossings.size(), tpcLines.size(), siTrajs.size(), siTpcIds.size());
      auto* c = new TCanvas(("c3_" + tag + "_phys_z_x_y_allcrossings").c_str(), xtitle.c_str(), 1200, 900);
      auto* h3 = owned(new TH3D(("h3_" + tag + "_phys_z_x_y_allcrossings").c_str(),
                                (xtitle + " - TPC (circles) + Si (squares);z [cm];x [cm];y [cm]").c_str(),
                                204, m_zmin, m_zmax, 170, -m_xymax, m_xymax, 170, -m_xymax, m_xymax));
      h3->SetDirectory(nullptr);
      h3->SetStats(false);
      h3->Draw();
      drawEverything3D(false, std::nullopt, singleCrossing, false);
      write_and_delete(c);
    }
    // silicon zoom, all crossings
    {
      auto* c = new TCanvas(("c3_" + tag + "_phys_si_z_x_y").c_str(), (title + " silicon").c_str(), 1200, 900);
      auto* h3 = owned(new TH3D(("h3_" + tag + "_phys_si_z_x_y").c_str(),
                                (title + " silicon zoom;z [cm];x [cm];y [cm]").c_str(),
                                60, -30, 30, 52, -kSiZoomR, kSiZoomR, 52, -kSiZoomR, kSiZoomR));
      h3->SetDirectory(nullptr);
      h3->SetStats(false);
      h3->Draw();
      drawEverything3D(true, std::nullopt, singleCrossing, true);
      draw_markers_zxy({{beamX(), beamY(), 0.0}}, kBlack, 29, 1.5);
      auto* tx = owned(new TLatex(0.02, 0.97, std::format("#star beam (TPC beam axis)   shift ({:.1f}, {:.1f}) mm, Si halves on TPC beam axis {}   {} crossings   {}",
                                                          10 * center[0], 10 * center[1], m_alignment.enabled() ? "on" : "off",
                                                          crossings.size(), siSummary)
                                                  .c_str()));
      tx->SetNDC();
      tx->SetTextSize(0.022);
      tx->Draw();
      auto* tx2 = owned(new TLatex(0.02, 0.94, siTpcSummary.c_str()));
      tx2->SetNDC();
      tx2->SetTextSize(0.022);
      tx2->Draw();
      write_and_delete(c);
    }
  }

  for (auto* o : legendStyles)
  {
    delete o;
  }

  unsigned int nTpcClusters = 0;
  for (const auto& [id, pts] : tpcClustersById)
  {
    nTpcClusters += pts.size();
  }
  std::cout << Name() << " - saved event " << m_evt << ": Si " << siSummary << "; " << siTpcSummary << "; TPC clusters " << nTpcClusters
            << ", TPC poly tracks " << tpcLines.size() << ", TPC collision vertices " << tpcCollision.size() << std::endl;
  ++m_eventsSaved;
  // keys and header on disk after every saved event: the file stays readable if the job does
  // not reach End (crash, kill, time limit)
  m_outfile->Write();
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpc_EventDisplay::End(PHCompositeNode* /*unused*/)
{
  if (m_outfile)
  {
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
  }
  std::cout << Name() << "::End - events seen: " << m_evt << ", events written: " << m_eventsSaved
            << ", output file: " << m_outfilename << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}