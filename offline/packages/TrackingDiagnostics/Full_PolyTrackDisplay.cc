#include "Full_PolyTrackDisplay.h"

#include <tpctrackreco/Full_PolyTrack.h>
#include <tpctrackreco/Full_PolyTrackContainer.h>
#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>

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
#include <format>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace
{
  struct Point
  {
    double x{0.0};
    double y{0.0};
    double z{0.0};
    unsigned int layer{0};
    bool silicon{false};
  };

  int crossingColor(const int crossing)
  {
    static const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2,
                                 kMagenta + 1, kCyan + 2, kOrange + 7,
                                 kViolet + 1, kAzure + 1, kPink + 7};
    constexpr int count = sizeof(colors) / sizeof(colors[0]);
    return colors[((crossing % count) + count) % count];
  }

  int trackColor(const unsigned int index)
  {
    static const int colors[] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                                 kCyan + 2, kOrange + 7, kViolet + 1, kAzure + 1};
    return colors[index % (sizeof(colors) / sizeof(colors[0]))];
  }

  bool visible(const Point& point, const double zmin, const double zmax, const double xymax)
  {
    return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) &&
           point.z >= zmin && point.z <= zmax &&
           std::abs(point.x) <= xymax && std::abs(point.y) <= xymax;
  }

  TPolyMarker3D* marker(const Point& point, const int color)
  {
    auto* output = new TPolyMarker3D(1);
    output->SetPoint(0, point.z, point.x, point.y);
    output->SetMarkerColor(color);
    output->SetMarkerStyle(point.silicon ? (point.layer <= 2U ? 20 : 21) : 24);
    output->SetMarkerSize(point.silicon ? 1.5 : 0.9);
    return output;
  }

  bool xyAtZ(const Full_PolyTrack& track, const double z, const double field,
             const bool straight, const double direction, double& x, double& y)
  {
    const double x0 = track.get_x();
    const double y0 = track.get_y();
    const double z0 = track.get_z();
    const double px = track.get_px();
    const double py = track.get_py();
    const double pz = track.get_pz();
    const double charge = track.get_charge();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        std::abs(pz) < 1.e-12)
    {
      return false;
    }
    const double dz = z - z0;
    if (straight || !std::isfinite(charge) || std::abs(charge * field) < 1.e-12)
    {
      x = x0 + direction * px * dz / pz;
      y = y0 + direction * py * dz / pz;
      return std::isfinite(x) && std::isfinite(y);
    }
    const double pt = std::hypot(px, py);
    if (pt <= 0.0) return false;
    const double signedRadius = pt / (0.003 * charge * field);
    const double radius = std::abs(signedRadius);
    const double sign = signedRadius > 0.0 ? 1.0 : -1.0;
    const double xc = x0 + sign * radius * py / pt;
    const double yc = y0 - sign * radius * px / pt;
    const double phi0 = std::atan2(y0 - yc, x0 - xc);
    const double arc = direction * dz * pt / pz;
    const double phi = phi0 - sign * arc / radius;
    x = xc + radius * std::cos(phi);
    y = yc + radius * std::sin(phi);
    return std::isfinite(x) && std::isfinite(y);
  }

  double residual2(const Full_PolyTrack& track, const std::vector<Point>& points,
                   const double field, const bool straight, const double direction)
  {
    double sum = 0.0;
    unsigned int count = 0;
    for (const auto& point : points)
    {
      double x = 0.0, y = 0.0;
      if (!xyAtZ(track, point.z, field, straight, direction, x, y)) continue;
      const double dx = x - point.x;
      const double dy = y - point.y;
      sum += dx * dx + dy * dy;
      ++count;
    }
    return count ? sum / count : std::numeric_limits<double>::max();
  }

  TPolyLine3D* fitLine(const Full_PolyTrack& track, const double zmin, const double zmax,
                       const double xymax, const double field, const bool straight,
                       const double direction, const int color)
  {
    std::vector<Point> points;
    for (unsigned int i = 0; i <= 100; ++i)
    {
      Point point;
      point.z = zmin + (zmax - zmin) * static_cast<double>(i) / 100.0;
      if (xyAtZ(track, point.z, field, straight, direction, point.x, point.y) &&
          std::abs(point.x) <= xymax && std::abs(point.y) <= xymax) points.push_back(point);
    }
    if (points.size() < 2) return nullptr;
    auto* line = new TPolyLine3D(points.size());
    for (unsigned int i = 0; i < points.size(); ++i) line->SetPoint(i, points[i].z, points[i].x, points[i].y);
    line->SetLineColor(color);
    line->SetLineWidth(3);
    return line;
  }
}

Full_PolyTrackDisplay::Full_PolyTrackDisplay(const std::string& name,
                                             const std::string& outfilename,
                                             const std::string& fullTrackNodeName,
                                             const unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_fullTrackNodeName(fullTrackNodeName)
  , m_maxEventDisplays(maxEventDisplays)
{
}

Full_PolyTrackDisplay::~Full_PolyTrackDisplay()
{
  delete m_outfile;
}

int Full_PolyTrackDisplay::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie()) return Fun4AllReturnCodes::ABORTRUN;
  m_outfile->mkdir("events");
  return Fun4AllReturnCodes::EVENT_OK;
}

bool Full_PolyTrackDisplay::getNodes(PHCompositeNode* topNode)
{
  m_fullTracks = findNode::getClass<Full_PolyTrackContainer>(topNode, m_fullTrackNodeName);
  m_tpcClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  m_trkrClusters = findNode::getClass<TrkrClusterContainer>(topNode, m_trkrClusterNodeName);
  m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_fullTracks) std::cerr << Name() << " - missing " << m_fullTrackNodeName << std::endl;
  if (!m_tpcClusters) std::cerr << Name() << " - missing " << m_tpcClusterNodeName << std::endl;
  if (!m_trkrClusters || !m_actsGeometry) std::cerr << Name() << " - silicon cluster position input missing" << std::endl;
  return m_fullTracks && m_tpcClusters && m_trkrClusters && m_actsGeometry;
}

int Full_PolyTrackDisplay::process_event(PHCompositeNode* topNode)
{
  ++m_event;
  if (m_eventsSaved >= m_maxEventDisplays || !getNodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;

  std::map<TrkrDefs::cluskey, const Tpc_PolyCluster*> tpcByKey;
  for (unsigned int i = 0; i < m_tpcClusters->size(); ++i)
  {
    const auto* cluster = m_tpcClusters->get_cluster(i);
    if (cluster) tpcByKey[cluster->get_trkr_cluster_key()] = cluster;
  }

  std::map<int, std::vector<TPolyMarker3D*>> markersByCrossing;
  std::map<int, std::vector<TPolyLine3D*>> linesByCrossing;
  std::vector<TPolyMarker3D*> allMarkers;
  std::vector<TPolyLine3D*> allLines;
  unsigned int selected = 0;
  for (unsigned int i = 0; i < m_fullTracks->size(); ++i)
  {
    const auto* track = m_fullTracks->get_track(i);
    if (!track || !track->isValid() || track->get_n_mvtx() < m_minMvtxHits ||
        track->get_n_intt() < m_minInttHits || std::hypot(track->get_px(), track->get_py()) < m_minTrackPt) continue;
    std::vector<Point> points;
    for (const auto key : track->get_tpc_cluster_keys())
    {
      const auto found = tpcByKey.find(key);
      if (found == tpcByKey.end()) continue;
      const auto* cluster = found->second;
      points.push_back({cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z(), TrkrDefs::getLayer(key), false});
    }
    for (const auto key : track->get_silicon_cluster_keys())
    {
      auto* cluster = m_trkrClusters->findCluster(key);
      if (!cluster) continue;
      const auto global = m_actsGeometry->getGlobalPosition(key, cluster);
      points.push_back({global.x(), global.y(), global.z(), TrkrDefs::getLayer(key), true});
    }
    const int crossing = track->get_crossing();
    const int color = crossingColor(crossing);
    const int individualColor = trackColor(i);
    for (const auto& point : points)
    {
      if (!visible(point, m_zmin, m_zmax, m_xymax)) continue;
      allMarkers.push_back(marker(point, color));
      markersByCrossing[crossing].push_back(marker(point, individualColor));
    }
    const bool straight = m_useStraightLineTracks || std::abs(track->get_charge() * m_magneticFieldTesla) < 1.e-12;
    const double direction = residual2(*track, points, m_magneticFieldTesla, straight, 1.0) <=
                                     residual2(*track, points, m_magneticFieldTesla, straight, -1.0) ? 1.0 : -1.0;
    if (auto* line = fitLine(*track, m_zmin, m_zmax, m_xymax, m_magneticFieldTesla, straight, direction, color)) allLines.push_back(line);
    if (auto* line = fitLine(*track, m_zmin, m_zmax, m_xymax, m_magneticFieldTesla, straight, direction, individualColor)) linesByCrossing[crossing].push_back(line);
    ++selected;
  }

  auto* events = m_outfile->GetDirectory("events");
  events->cd();
  auto* eventDir = events->mkdir(std::format("event_{:06}", m_event).c_str());
  eventDir->cd();
  auto draw = [&](const std::string& suffix, const std::string& title,
                  const std::vector<TPolyMarker3D*>& markers, const std::vector<TPolyLine3D*>& lines)
  {
    auto* hist = new TH3D(("h3_" + suffix).c_str(), (title + ";z [cm];x [cm];y [cm]").c_str(),
                          204, m_zmin, m_zmax, 170, -m_xymax, m_xymax, 170, -m_xymax, m_xymax);
    hist->SetStats(false); hist->SetDirectory(nullptr);
    auto* canvas = new TCanvas(("c3_" + suffix).c_str(), title.c_str(), 1200, 900);
    hist->Draw();
    for (auto* line : lines) if (line) line->Draw("same");
    for (auto* value : markers) if (value) value->Draw("same");
    canvas->Write();
  };
  draw(std::format("evt{:06}_full_polytracks", m_event), std::format("event {} FULL_POLYTRACKS", m_event), allMarkers, allLines);
  for (const auto& entry : markersByCrossing)
  {
    draw(std::format("evt{:06}_full_polytracks_crossing_{}", m_event, entry.first),
         std::format("event {} FULL_POLYTRACKS crossing {}", m_event, entry.first), entry.second, linesByCrossing[entry.first]);
  }
  std::cout << Name() << " - saved event " << m_event << " selected_tracks=" << selected << std::endl;
  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackDisplay::End(PHCompositeNode*)
{
  if (m_outfile) { m_outfile->Close(); delete m_outfile; m_outfile = nullptr; }
  return Fun4AllReturnCodes::EVENT_OK;
}
