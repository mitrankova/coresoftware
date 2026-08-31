#include "Full_PolyTrackDisplay.h"

#include "tpctrackreco/Full_PolyTrack.h"
#include "tpctrackreco/Full_PolyTrackContainer.h"
#include "tpctrackreco/TpcCrossingDecision.h"
#include "tpctrackreco/TpcCrossingDecisionContainer.h"
#include "tpctrackreco/Tpc_PolyCluster.h"
#include "tpctrackreco/Tpc_PolyClusterContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace
{
  struct DisplayPoint
  {
    double x{std::numeric_limits<double>::quiet_NaN()};
    double y{std::numeric_limits<double>::quiet_NaN()};
    double z{std::numeric_limits<double>::quiet_NaN()};
    unsigned int layer{0};
    TrkrDefs::cluskey key{TrkrDefs::CLUSKEYMAX};
  };

  int crossing_color(const int crossing)
  {
    static const int colors[] = {
        kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
        kCyan + 2, kOrange + 7, kViolet + 1, kAzure + 1, kPink + 7,
        kTeal + 3, kSpring + 5};
    const int ncolors = static_cast<int>(sizeof(colors) / sizeof(colors[0]));
    const unsigned int index = static_cast<unsigned int>(((crossing % ncolors) + ncolors) % ncolors);
    return colors[index % (sizeof(colors) / sizeof(colors[0]))];
  }

  int track_color(const unsigned int track_index)
  {
    static const int colors[] = {
        kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
        kOrange + 7, kViolet + 1, kAzure + 1, kPink + 7, kTeal + 3,
        kSpring + 5, kYellow + 2};
    return colors[track_index % (sizeof(colors) / sizeof(colors[0]))];
  }

  int layer_marker_style(const unsigned int layer)
  {
    if (layer <= 2U)
    {
      return 20;
    }
    if (layer <= 6U)
    {
      return 21;
    }
    return 24;
  }

  bool finite_point(const DisplayPoint& point)
  {
    return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z);
  }

  TPolyMarker3D* make_cluster_marker(const DisplayPoint& point, const int color)
  {
    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, point.z, point.x, point.y);
    marker->SetMarkerColor(color);
    marker->SetMarkerStyle(layer_marker_style(point.layer));
    marker->SetMarkerSize(point.layer <= 2U ? 1.5 : 1.7);
    return marker;
  }

  TPolyMarker3D* make_unused_silicon_seed_marker(const DisplayPoint& point, const double marker_size)
  {
    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, point.z, point.x, point.y);
    marker->SetMarkerColor(kGray + 1);
    marker->SetMarkerStyle(24);
    marker->SetMarkerSize(marker_size);
    return marker;
  }

  TPolyMarker3D* make_tpc_cluster_marker(const DisplayPoint& point, const int color)
  {
    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, point.z, point.x, point.y);
    marker->SetMarkerColor(color);
    marker->SetMarkerStyle(24);
    marker->SetMarkerSize(0.9);
    return marker;
  }

  TPolyMarker3D* make_track_seed_marker(const Full_PolyTrack* trk, const int color)
  {
    if (!trk)
    {
      return nullptr;
    }
    const double x = trk->get_seed_x();
    const double y = trk->get_seed_y();
    const double z = trk->get_seed_z();
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
    {
      return nullptr;
    }
    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, z, x, y);
    marker->SetMarkerColor(color);
    marker->SetMarkerStyle(29);
    marker->SetMarkerSize(1.4);
    return marker;
  }

  TPolyLine3D* make_cluster_chain_line(std::vector<DisplayPoint> points, const int color)
  {
    points.erase(std::remove_if(points.begin(), points.end(), [](const DisplayPoint& point) {
                   return !finite_point(point);
                 }),
                 points.end());
    std::sort(points.begin(), points.end(), [](const DisplayPoint& a, const DisplayPoint& b) {
      const double ra = std::hypot(a.x, a.y);
      const double rb = std::hypot(b.x, b.y);
      return ra > rb;
    });
    if (points.size() < 2)
    {
      return nullptr;
    }

    TPolyLine3D* line = new TPolyLine3D(static_cast<int>(points.size()));
    for (unsigned int i = 0; i < points.size(); ++i)
    {
      line->SetPoint(static_cast<int>(i), points[i].z, points[i].x, points[i].y);
    }
    line->SetLineColor(color);
    line->SetLineWidth(2);
    return line;
  }


  bool point_group_z_range(const std::vector<DisplayPoint>& points,
                           const double display_zmin,
                           const double display_zmax,
                           double& zmin,
                           double& zmax)
  {
    if (points.empty())
    {
      return false;
    }

    zmin = std::numeric_limits<double>::max();
    zmax = -std::numeric_limits<double>::max();
    for (const DisplayPoint& point : points)
    {
      if (!finite_point(point))
      {
        continue;
      }
      zmin = std::min(zmin, point.z);
      zmax = std::max(zmax, point.z);
    }

    if (zmin == std::numeric_limits<double>::max() || zmax == -std::numeric_limits<double>::max())
    {
      return false;
    }

    zmin = std::max(zmin, display_zmin);
    zmax = std::min(zmax, display_zmax);
    if (zmin > zmax)
    {
      return false;
    }
    if (zmin == zmax)
    {
      zmin = std::max(zmin - 0.1, display_zmin);
      zmax = std::min(zmax + 0.1, display_zmax);
    }
    return zmin < zmax;
  }

  void extend_z_range_to_tpc_seed_point(const Full_PolyTrack* trk,
                                          const double display_zmin,
                                          const double display_zmax,
                                          double& zmin,
                                          double& zmax)
  {
    if (!trk)
    {
      return;
    }
    const double z = trk->get_seed_z();
    if (!std::isfinite(z))
    {
      return;
    }
    const double clipped_z = std::max(display_zmin, std::min(display_zmax, z));
    zmin = std::min(zmin, clipped_z);
    zmax = std::max(zmax, clipped_z);
  }

  bool tpc_seed_xy_at_z(const Full_PolyTrack* trk,
                          const double z,
                          const double magnetic_field_tesla,
                          const double arc_direction,
                          const bool use_straight_line,
                          double& x,
                          double& y)
  {
    if (!trk || trk->get_fit_status() == 0 || !std::isfinite(z))
    {
      return false;
    }

    const double x0 = trk->get_seed_x();
    const double y0 = trk->get_seed_y();
    const double z0 = trk->get_seed_z();
    const double px = trk->get_seed_px();
    const double py = trk->get_seed_py();
    const double pz = trk->get_seed_pz();
    const double charge = trk->get_charge();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(charge))
    {
      return false;
    }

    const double dz = z - z0;
    if (use_straight_line || std::fabs(charge * magnetic_field_tesla) < 1.0e-12)
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

    const double signed_radius = pt / (0.003 * charge * magnetic_field_tesla);
    const double radius = std::fabs(signed_radius);
    if (!std::isfinite(radius) || radius <= 0.0)
    {
      return false;
    }

    const double tx = px / pt;
    const double ty = py / pt;
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

  double tpc_seed_line_residual2(const Full_PolyTrack* trk,
                              const std::vector<DisplayPoint>& points,
                              const double magnetic_field_tesla,
                              const double arc_direction,
                              const bool use_straight_line)
  {
    if (!trk || points.empty())
    {
      return std::numeric_limits<double>::max();
    }

    double sum = 0.0;
    unsigned int n = 0;
    for (const DisplayPoint& point : points)
    {
      if (!finite_point(point))
      {
        continue;
      }

      double x = 0.0;
      double y = 0.0;
      if (!tpc_seed_xy_at_z(trk, point.z, magnetic_field_tesla, arc_direction, use_straight_line, x, y))
      {
        continue;
      }

      const double dx = x - point.x;
      const double dy = y - point.y;
      sum += dx * dx + dy * dy;
      ++n;
    }

    return n > 0 ? sum / static_cast<double>(n) : std::numeric_limits<double>::max();
  }

  TPolyLine3D* make_tpc_seed_fit_line(const Full_PolyTrack* trk,
                                        const double zmin,
                                        const double zmax,
                                        const double xymax,
                                        const double magnetic_field_tesla,
                                        const double arc_direction,
                                        const bool use_straight_line,
                                        const int color)
  {
    if (!trk || trk->get_fit_status() == 0)
    {
      return nullptr;
    }

    std::vector<double> zs;
    std::vector<double> xs;
    std::vector<double> ys;
    const unsigned int nsteps = 80;
    zs.reserve(nsteps + 1);
    xs.reserve(nsteps + 1);
    ys.reserve(nsteps + 1);

    for (unsigned int istep = 0; istep <= nsteps; ++istep)
    {
      const double f = static_cast<double>(istep) / static_cast<double>(nsteps);
      const double z = zmin + f * (zmax - zmin);
      double x = 0.0;
      double y = 0.0;
      if (!tpc_seed_xy_at_z(trk, z, magnetic_field_tesla, arc_direction, use_straight_line, x, y))
      {
        continue;
      }
      if (std::fabs(x) > xymax || std::fabs(y) > xymax)
      {
        continue;
      }
      zs.push_back(z);
      xs.push_back(x);
      ys.push_back(y);
    }

    if (zs.size() < 2)
    {
      return nullptr;
    }
    TPolyLine3D* line = new TPolyLine3D(static_cast<int>(zs.size()));
    for (unsigned int i = 0; i < zs.size(); ++i)
    {
      line->SetPoint(static_cast<int>(i), zs[i], xs[i], ys[i]);
    }
    line->SetLineColor(color);
    line->SetLineStyle(2);
    line->SetLineWidth(3);
    return line;
  }

}  // namespace

Full_PolyTrackDisplay::Full_PolyTrackDisplay(const std::string& name,
                                             const std::string& outfilename,
                                             const std::string& fullTrackNodeName,
                                             const unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_fullTrackNodeName(fullTrackNodeName)
  , m_tpcClusterNodeName("TPC_POLYCLUSTERS")
  , m_trkrClusterNodeName("TRKR_CLUSTER")
  , m_crossingDecisionNodeName("TPC_CROSSING_DECISIONS")
  , m_maxEventDisplays(maxEventDisplays)
  , m_evt(0)
  , m_eventsSaved(0)
  , m_minMvtxHits(3)
  , m_minInttHits(0)
  , m_zmin(-102.0)
  , m_zmax(102.0)
  , m_xymax(85.0)
  , m_magneticFieldTesla(1.4)
  , m_unusedSiliconSeedMarkerSize(0.25)
  , m_minTrackPt(0.1)
  , m_drawTrackLines(true)
  , m_useStraightLineTracks(false)
  , m_drawTpcOnlyFullPolyTracks(true)
  , m_drawUnusedSiliconSeeds(false)
  , m_outfile(nullptr)
  , m_fullTracks(nullptr)
  , m_tpcClusters(nullptr)
  , m_trkrClusters(nullptr)
  , m_actsGeometry(nullptr)
  , m_crossingDecisions(nullptr)
{
}

Full_PolyTrackDisplay::~Full_PolyTrackDisplay()
{
  if (m_outfile)
  {
    delete m_outfile;
    m_outfile = nullptr;
  }
}

int Full_PolyTrackDisplay::Init(PHCompositeNode* /*unused*/)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open output file " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_outfile->mkdir("events");
  return Fun4AllReturnCodes::EVENT_OK;
}

bool Full_PolyTrackDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_fullTracks = findNode::getClass<Full_PolyTrackContainer>(topNode, m_fullTrackNodeName);
  if (!m_fullTracks)
  {
    std::cerr << Name() << " - missing " << m_fullTrackNodeName << std::endl;
    return false;
  }

  m_tpcClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  if (!m_tpcClusters)
  {
    std::cerr << Name() << " - missing " << m_tpcClusterNodeName
              << ", drawing without Tpc_PolyCluster points" << std::endl;
  }

  m_trkrClusters = findNode::getClass<TrkrClusterContainer>(topNode, m_trkrClusterNodeName);
  if (!m_trkrClusters)
  {
    std::cerr << Name() << " - missing " << m_trkrClusterNodeName
              << ", drawing only stored Full_PolyTrack states" << std::endl;
  }

  m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_actsGeometry && m_trkrClusters)
  {
    std::cerr << Name() << " - missing ActsGeometry"
              << ", cannot resolve cluster keys without stored states" << std::endl;
  }

  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  if (!m_crossingDecisions)
  {
    std::cerr << Name() << " - missing " << m_crossingDecisionNodeName
              << ", coloring by source track id fallback" << std::endl;
  }
  return true;
}

int Full_PolyTrackDisplay::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  if (!get_nodes(topNode))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  if (!m_outfile || m_eventsSaved >= m_maxEventDisplays)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  std::map<unsigned int, int> crossing_by_source_id;
  const unsigned int ndecisions = m_crossingDecisions ? m_crossingDecisions->size() : 0;
  for (unsigned int idecision = 0; idecision < ndecisions; ++idecision)
  {
    const TpcCrossingDecision* decision = m_crossingDecisions->get_decision(idecision);
    if (decision)
    {
      crossing_by_source_id[decision->get_assembled_track_id()] = decision->get_selected_crossing();
    }
  }

  TDirectory* eventsTop = m_outfile->GetDirectory("events");
  if (!eventsTop)
  {
    eventsTop = m_outfile->mkdir("events");
  }
  eventsTop->cd();

  TDirectory* eventDir = eventsTop->mkdir(std::format("event_{:06}", m_evt).c_str());
  if (!eventDir)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  eventDir->cd();

  TH3D* h3 = new TH3D(std::format("h3_evt{:06}_full_polytrack_clusters", m_evt).c_str(),
                      std::format("event {} Full_PolyTrack selected clusters;z [cm];x [cm];y [cm]", m_evt).c_str(),
                      204, m_zmin, m_zmax,
                      170, -m_xymax, m_xymax,
                      170, -m_xymax, m_xymax);
  h3->SetStats(false);
  h3->SetDirectory(nullptr);

  std::vector<TPolyMarker3D*> cluster_markers;
  std::vector<TPolyMarker3D*> tpc_cluster_markers;
  std::vector<TPolyMarker3D*> seed_markers;
  std::vector<TPolyMarker3D*> unused_silicon_seed_markers;
  std::vector<TPolyLine3D*> lines;
  std::vector<TPolyLine3D*> fit_lines;
  std::map<int, std::vector<TPolyMarker3D*> > cluster_track_markers_by_crossing;
  std::map<int, std::vector<TPolyMarker3D*> > tpc_cluster_track_markers_by_crossing;
  std::map<int, std::vector<TPolyMarker3D*> > seed_track_markers_by_crossing;
  std::map<int, std::vector<TPolyLine3D*> > track_lines_by_crossing;
  std::map<int, std::vector<TPolyLine3D*> > fit_lines_by_crossing;
  std::map<int, bool> crossings_to_draw;
  std::set<TrkrDefs::cluskey> used_silicon_keys;

  unsigned int ntracks_selected = 0;
  unsigned int nclusters_selected = 0;
  unsigned int nclusters_plotted = 0;
  unsigned int ntpc_clusters_selected = 0;
  unsigned int ntpc_clusters_plotted = 0;
  unsigned int nunused_silicon_seeds_selected = 0;
  unsigned int nunused_silicon_seeds_plotted = 0;
  std::map<unsigned int, std::vector<const Tpc_PolyCluster*> > tpc_clusters_by_source_id;
  const unsigned int ntpc_clusters_total = m_tpcClusters ? m_tpcClusters->size() : 0;
  for (unsigned int icluster = 0; icluster < ntpc_clusters_total; ++icluster)
  {
    const Tpc_PolyCluster* cluster = m_tpcClusters->get_cluster(icluster);
    if (!cluster || !cluster->isValid())
    {
      continue;
    }
    tpc_clusters_by_source_id[cluster->get_source_assembled_track_id()].push_back(cluster);
  }
  const unsigned int ntracks = m_fullTracks->size();
  for (unsigned int itrack = 0; itrack < ntracks; ++itrack)
  {
    const Full_PolyTrack* trk = m_fullTracks->get_track(itrack);
    if (!trk || !trk->isValid())
    {
      continue;
    }

    const double track_pt = std::hypot(trk->get_seed_px(), trk->get_seed_py());
    if (!std::isfinite(track_pt) || track_pt < m_minTrackPt)
    {
      continue;
    }

    const unsigned int source_id = trk->get_source_assembled_track_id();

    std::vector<DisplayPoint> points;
    std::vector<DisplayPoint> tpc_points;
    const auto tpc_cluster_iter = tpc_clusters_by_source_id.find(source_id);
    if (tpc_cluster_iter != tpc_clusters_by_source_id.end())
    {
      for (const Tpc_PolyCluster* cluster : tpc_cluster_iter->second)
      {
        if (!cluster)
        {
          continue;
        }
        DisplayPoint point;
        point.key = cluster->get_trkr_cluster_key();
        point.layer = point.key != TrkrDefs::CLUSKEYMAX ? TrkrDefs::getLayer(point.key) : 999U;
        point.x = cluster->get_centroid_x();
        point.y = cluster->get_centroid_y();
        point.z = cluster->get_centroid_z();
        tpc_points.push_back(point);
      }
    }
    const unsigned int nstates = trk->size_silicon_states();
    for (unsigned int istate = 0; istate < nstates; ++istate)
    {
      DisplayPoint point;
      point.key = trk->get_state_cluster_key(istate);
      point.layer = trk->get_state_layer(istate);
      point.x = trk->get_state_x(istate);
      point.y = trk->get_state_y(istate);
      point.z = trk->get_state_z(istate);
      points.push_back(point);
    }

    if (points.empty() && m_trkrClusters && m_actsGeometry)
    {
      const unsigned int nkeys = trk->size_cluster_keys();
      for (unsigned int ikey = 0; ikey < nkeys; ++ikey)
      {
        const TrkrDefs::cluskey key = trk->get_cluster_key(ikey);
        TrkrCluster* cluster = m_trkrClusters->findCluster(key);
        if (!cluster)
        {
          continue;
        }
        const auto global = m_actsGeometry->getGlobalPosition(key, cluster);
        DisplayPoint point;
        point.key = key;
        point.layer = TrkrDefs::getLayer(key);
        point.x = global.x();
        point.y = global.y();
        point.z = global.z();
        points.push_back(point);
      }
    }

    unsigned int nmvtx_hits = 0;
    unsigned int nintt_hits = 0;
    for (const DisplayPoint& point : points)
    {
      if (point.key != TrkrDefs::CLUSKEYMAX)
      {
        const TrkrDefs::TrkrId trkrid = static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(point.key));
        if (trkrid == TrkrDefs::mvtxId)
        {
          ++nmvtx_hits;
        }
        else if (trkrid == TrkrDefs::inttId)
        {
          ++nintt_hits;
        }
      }
      else if (point.layer <= 2U)
      {
        ++nmvtx_hits;
      }
      else if (point.layer >= 3U && point.layer <= 6U)
      {
        ++nintt_hits;
      }
    }
    if (nmvtx_hits < m_minMvtxHits || nintt_hits < m_minInttHits)
    {
      continue;
    }

    if (!m_drawTpcOnlyFullPolyTracks && points.empty())
    {
      continue;
    }

    for (const DisplayPoint& point : points)
    {
      if (point.key != TrkrDefs::CLUSKEYMAX)
      {
        used_silicon_keys.insert(point.key);
      }
    }

    ++ntracks_selected;

    const auto crossing_iter = crossing_by_source_id.find(source_id);
    const int crossing = crossing_iter != crossing_by_source_id.end()
                             ? crossing_iter->second
                             : static_cast<int>(source_id);
    const int color = crossing_color(crossing);
    const int per_crossing_color = track_color(source_id);
    crossings_to_draw[crossing] = true;

    TPolyMarker3D* seed_marker = make_track_seed_marker(trk, color);
    if (seed_marker)
    {
      seed_markers.push_back(seed_marker);
      TPolyMarker3D* seed_track_marker = make_track_seed_marker(trk, per_crossing_color);
      if (seed_track_marker)
      {
        seed_track_markers_by_crossing[crossing].push_back(seed_track_marker);
      }
    }

    for (const DisplayPoint& point : tpc_points)
    {
      ++ntpc_clusters_selected;
      if (!finite_point(point))
      {
        continue;
      }
      if (point.z < m_zmin || point.z > m_zmax ||
          std::fabs(point.x) > m_xymax || std::fabs(point.y) > m_xymax)
      {
        continue;
      }
      TPolyMarker3D* marker = make_tpc_cluster_marker(point, color);
      tpc_cluster_markers.push_back(marker);
      tpc_cluster_track_markers_by_crossing[crossing].push_back(make_tpc_cluster_marker(point, per_crossing_color));
      ++ntpc_clusters_plotted;
    }

    for (const DisplayPoint& point : points)
    {
      ++nclusters_selected;
      if (!finite_point(point))
      {
        continue;
      }
      if (point.z < m_zmin || point.z > m_zmax ||
          std::fabs(point.x) > m_xymax || std::fabs(point.y) > m_xymax)
      {
        continue;
      }
      TPolyMarker3D* marker = make_cluster_marker(point, color);
      cluster_markers.push_back(marker);
      cluster_track_markers_by_crossing[crossing].push_back(make_cluster_marker(point, per_crossing_color));
      ++nclusters_plotted;
    }

    if (m_drawTrackLines)
    {
      std::vector<DisplayPoint> line_points;
      line_points.reserve(tpc_points.size() + points.size());
      line_points.insert(line_points.end(), tpc_points.begin(), tpc_points.end());
      line_points.insert(line_points.end(), points.begin(), points.end());

      TPolyLine3D* line = make_cluster_chain_line(line_points, color);
      if (line)
      {
        lines.push_back(line);
      }

      TPolyLine3D* track_line = make_cluster_chain_line(line_points, per_crossing_color);
      if (track_line)
      {
        track_lines_by_crossing[crossing].push_back(track_line);
      }

      double fit_zmin = m_zmin;
      double fit_zmax = m_zmax;
      if (point_group_z_range(tpc_points, m_zmin, m_zmax, fit_zmin, fit_zmax))
      {
        extend_z_range_to_tpc_seed_point(trk, m_zmin, m_zmax, fit_zmin, fit_zmax);
        const bool use_straight_line = m_useStraightLineTracks || std::fabs(trk->get_charge() * m_magneticFieldTesla) < 1.0e-12;
        const double forward_residual2 = tpc_seed_line_residual2(trk, tpc_points, m_magneticFieldTesla, 1.0, use_straight_line);
        const double reverse_residual2 = tpc_seed_line_residual2(trk, tpc_points, m_magneticFieldTesla, -1.0, use_straight_line);
        const double arc_direction = forward_residual2 <= reverse_residual2 ? 1.0 : -1.0;

        TPolyLine3D* fit_line = make_tpc_seed_fit_line(trk, fit_zmin, fit_zmax, m_xymax,
                                                       m_magneticFieldTesla, arc_direction,
                                                       use_straight_line, color);
        if (fit_line)
        {
          fit_lines.push_back(fit_line);
        }

        TPolyLine3D* crossing_fit_line = make_tpc_seed_fit_line(trk, fit_zmin, fit_zmax, m_xymax,
                                                                m_magneticFieldTesla, arc_direction,
                                                                use_straight_line, per_crossing_color);
        if (crossing_fit_line)
        {
          fit_lines_by_crossing[crossing].push_back(crossing_fit_line);
        }
      }
    }
  }

  if (m_drawUnusedSiliconSeeds && m_trkrClusters && m_actsGeometry)
  {
    for (const TrkrDefs::TrkrId trkrid : {TrkrDefs::mvtxId, TrkrDefs::inttId})
    {
      for (unsigned int layer = 0; layer <= 6U; ++layer)
      {
        const auto hitset_keys = m_trkrClusters->getHitSetKeys(trkrid, static_cast<uint8_t>(layer));
        for (const TrkrDefs::hitsetkey hitset_key : hitset_keys)
        {
          const auto range = m_trkrClusters->getClusters(hitset_key);
          for (auto iter = range.first; iter != range.second; ++iter)
          {
            const TrkrDefs::cluskey key = iter->first;
            TrkrCluster* cluster = iter->second;
            if (!cluster || used_silicon_keys.find(key) != used_silicon_keys.end())
            {
              continue;
            }
            const auto global = m_actsGeometry->getGlobalPosition(key, cluster);
            DisplayPoint point;
            point.key = key;
            point.layer = TrkrDefs::getLayer(key);
            point.x = global.x();
            point.y = global.y();
            point.z = global.z();
            ++nunused_silicon_seeds_selected;
            if (!finite_point(point))
            {
              continue;
            }
            if (point.z < m_zmin || point.z > m_zmax ||
                std::fabs(point.x) > m_xymax || std::fabs(point.y) > m_xymax)
            {
              continue;
            }
            unused_silicon_seed_markers.push_back(make_unused_silicon_seed_marker(point, m_unusedSiliconSeedMarkerSize));
            ++nunused_silicon_seeds_plotted;
          }
        }
      }
    }
  }

  TCanvas* c3 = new TCanvas(std::format("c3_evt{:06}_full_polytrack_clusters", m_evt).c_str(),
                            std::format("event {} Full_PolyTrack clusters", m_evt).c_str(),
                            1200, 900);
  h3->Draw();
  for (TPolyMarker3D* marker : unused_silicon_seed_markers)
  {
    if (marker)
    {
      marker->Draw("same");
    }
  }
  for (TPolyLine3D* line : lines)
  {
    if (line)
    {
      line->Draw("same");
    }
  }
  for (TPolyMarker3D* marker : tpc_cluster_markers)
  {
    if (marker)
    {
      marker->Draw("same");
    }
  }
  for (TPolyMarker3D* marker : cluster_markers)
  {
    if (marker)
    {
      marker->Draw("same");
    }
  }
  for (TPolyMarker3D* marker : seed_markers)
  {
    if (marker)
    {
      marker->Draw("same");
    }
  }
  for (TPolyLine3D* line : fit_lines)
  {
    if (line)
    {
      line->Draw("same");
    }
  }
  c3->Write();
  for (const auto& crossing_item : crossings_to_draw)
  {
    const int crossing = crossing_item.first;
    TH3D* hcross = new TH3D(std::format("h3_evt{:06}_full_polytrack_clusters_crossing_{}", m_evt, crossing).c_str(),
                            std::format("event {} Full_PolyTrack clusters crossing {};z [cm];x [cm];y [cm]", m_evt, crossing).c_str(),
                            204, m_zmin, m_zmax,
                            170, -m_xymax, m_xymax,
                            170, -m_xymax, m_xymax);
    hcross->SetStats(false);
    hcross->SetDirectory(nullptr);

    TCanvas* ccross = new TCanvas(std::format("c3_evt{:06}_full_polytrack_clusters_crossing_{}", m_evt, crossing).c_str(),
                                  std::format("event {} Full_PolyTrack clusters crossing {}", m_evt, crossing).c_str(),
                                  1200, 900);
    hcross->Draw();
    for (TPolyMarker3D* marker : unused_silicon_seed_markers)
    {
      if (marker)
      {
        marker->Draw("same");
      }
    }
    for (TPolyLine3D* line : track_lines_by_crossing[crossing])
    {
      if (line)
      {
        line->Draw("same");
      }
    }
    for (TPolyMarker3D* marker : tpc_cluster_track_markers_by_crossing[crossing])
    {
      if (marker)
      {
        marker->Draw("same");
      }
    }
    for (TPolyMarker3D* marker : cluster_track_markers_by_crossing[crossing])
    {
      if (marker)
      {
        marker->Draw("same");
      }
    }
    for (TPolyMarker3D* marker : seed_track_markers_by_crossing[crossing])
    {
      if (marker)
      {
        marker->Draw("same");
      }
    }
    for (TPolyLine3D* line : fit_lines_by_crossing[crossing])
    {
      if (line)
      {
        line->Draw("same");
      }
    }
    ccross->Write();
  }


  std::cout << Name() << " - saved event " << m_evt
            << " full_tracks=" << ntracks
            << " selected_tracks=" << ntracks_selected
            << " si_clusters=" << nclusters_selected
            << " si_plotted=" << nclusters_plotted
            << " tpc_clusters=" << ntpc_clusters_selected
            << " tpc_plotted=" << ntpc_clusters_plotted
            << " unused_si_seeds=" << nunused_silicon_seeds_selected
            << " unused_si_plotted=" << nunused_silicon_seeds_plotted
            << " lines=" << lines.size()
            << " fit_lines=" << fit_lines.size()
            << " crossing_canvases=" << crossings_to_draw.size()
            << " crossing_decisions=" << ndecisions << std::endl;

  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackDisplay::End(PHCompositeNode* /*unused*/)
{
  if (m_outfile)
  {
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
