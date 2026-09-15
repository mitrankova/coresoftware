#include "Full_PolyTrackDisplay.h"

#include <tpctrackreco/Full_PolyTrack.h>
#include <tpctrackreco/Full_PolyTrackContainer.h>
#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/TpcTrackKalmanFitter.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phfield/PHFieldUtility.h>
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
    TrkrDefs::cluskey key{TrkrDefs::CLUSKEYMAX};
    TrkrDefs::TrkrId detector{TrkrDefs::tpcId};
    bool silicon{false};
  };

  const char* detectorName(const TrkrDefs::TrkrId detector)
  {
    if (detector == TrkrDefs::mvtxId) return "MVTX";
    if (detector == TrkrDefs::inttId) return "INTT";
    if (detector == TrkrDefs::tpcId) return "TPC";
    return "UNKNOWN";
  }

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

  TPolyLine3D* fittedTrajectory(const Full_PolyTrack& track,
                                const std::vector<Point>& measurements,
                                const PHField* field,
                                const bool straight,
                                const double fallbackField,
                                const double zmin,
                                const double zmax,
                                const double xymax,
                                const int color,
                                unsigned int& attempted,
                                unsigned int& failed,
                                unsigned int& samples)
  {
    const double px = track.get_px();
    const double py = track.get_py();
    const double pz = track.get_pz();
    const double pt = std::hypot(px, py);
    if (!(pt > 0.) || measurements.empty()) return nullptr;
    std::array<double, TpcTrackKalmanFitter::StateDim> state{{
        track.get_x(), track.get_y(), track.get_z(), std::atan2(py, px),
        straight ? 0. : track.get_charge() / pt, pz / pt}};
    for (const double value : state) if (!std::isfinite(value)) return nullptr;
    TpcKalmanConfig config;
    config.magnetic_field = straight ? nullptr : field;
    config.analytic_uniform_propagation = straight;
    config.bfield_t = straight ? 0. : fallbackField;
    const auto propagate = [&](const std::array<double, TpcTrackKalmanFitter::StateDim>& input,
                               const double ds)
    {
      ++attempted;
      return TpcTrackKalmanFitter::propagate_state(input, ds, config);
    };
    const auto plus = propagate(state, 0.25);
    const auto minus = propagate(state, -0.25);
    const double radius0 = std::hypot(state[TpcTrackKalmanFitter::X], state[TpcTrackKalmanFitter::Y]);
    const double radiusPlus = std::hypot(plus[TpcTrackKalmanFitter::X], plus[TpcTrackKalmanFitter::Y]);
    const double radiusMinus = std::hypot(minus[TpcTrackKalmanFitter::X], minus[TpcTrackKalmanFitter::Y]);
    if (!std::isfinite(radiusPlus) || !std::isfinite(radiusMinus)) { ++failed; return nullptr; }
    const double step = radiusPlus >= radiusMinus ? 0.25 : -0.25;
    double maximumRadius = 0.;
    for (const auto& point : measurements) maximumRadius = std::max(maximumRadius, std::hypot(point.x, point.y));
    std::vector<Point> points;
    if (visible({state[0], state[1], state[2]}, zmin, zmax, xymax))
      points.push_back({state[0], state[1], state[2]});
    double previousRadius = radius0;
    for (unsigned int i = 0; i < 2000 && previousRadius <= maximumRadius + 0.25; ++i)
    {
      const auto next = propagate(state, step);
      const double radius = std::hypot(next[TpcTrackKalmanFitter::X], next[TpcTrackKalmanFitter::Y]);
      if (!std::isfinite(radius) || !std::isfinite(next[TpcTrackKalmanFitter::Z])) { ++failed; break; }
      if (radius + 1.e-4 < previousRadius) { ++failed; break; }
      state = next;
      previousRadius = radius;
      Point point{state[0], state[1], state[2]};
      if (visible(point, zmin, zmax, xymax)) points.push_back(point);
    }
    if (points.size() < 2) return nullptr;
    samples += points.size();
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
  m_field = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
  if (!m_fullTracks) std::cerr << Name() << " - missing " << m_fullTrackNodeName << std::endl;
  if (!m_tpcClusters) std::cerr << Name() << " - missing " << m_tpcClusterNodeName << std::endl;
  if (!m_trkrClusters || !m_actsGeometry) std::cerr << Name() << " - silicon cluster position input missing" << std::endl;
  if (!m_field && !m_useStraightLineTracks) std::cerr << Name() << " - magnetic field input missing" << std::endl;
  return m_fullTracks && m_tpcClusters && m_trkrClusters && m_actsGeometry &&
         (m_field || m_useStraightLineTracks);
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
  unsigned int propagationAttempted = 0;
  unsigned int propagationFailed = 0;
  unsigned int trajectorySamples = 0;
  bool dumpedTrack = false;
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
      points.push_back({cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z(),
                        TrkrDefs::getLayer(key), key, TrkrDefs::tpcId, false});
    }
    for (const auto key : track->get_silicon_cluster_keys())
    {
      auto* cluster = m_trkrClusters->findCluster(key);
      if (!cluster) continue;
      const auto global = m_actsGeometry->getGlobalPosition(key, cluster);
      points.push_back({global.x(), global.y(), global.z(), TrkrDefs::getLayer(key), key,
                        static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(key)), true});
    }
    if (Verbosity() >= 10 && !dumpedTrack)
    {
      std::cout << Name() << " display_track parent=" << track->get_parent_tpc_track_id()
                << " order=stored_before_sort" << std::endl;
      for (unsigned int index = 0; index < points.size(); ++index)
      {
        const auto& point = points[index];
        std::cout << "  idx=" << index << " det=" << detectorName(point.detector)
                  << " layer=" << point.layer << " cluster_key=" << point.key
                  << " x=" << point.x << " y=" << point.y << " z=" << point.z
                  << " r=" << std::hypot(point.x, point.y)
                  << " phi=" << std::atan2(point.y, point.x) << std::endl;
      }
    }
    std::sort(points.begin(), points.end(), [](const Point& lhs, const Point& rhs)
    {
      return std::hypot(lhs.x, lhs.y) > std::hypot(rhs.x, rhs.y);
    });
    if (Verbosity() >= 10 && !dumpedTrack)
    {
      std::cout << Name() << " display_track parent=" << track->get_parent_tpc_track_id()
                << " order=radius_outer_to_inner" << std::endl;
      for (unsigned int index = 0; index < points.size(); ++index)
      {
        const auto& point = points[index];
        std::cout << "  idx=" << index << " det=" << detectorName(point.detector)
                  << " layer=" << point.layer << " cluster_key=" << point.key
                  << " x=" << point.x << " y=" << point.y << " z=" << point.z
                  << " r=" << std::hypot(point.x, point.y)
                  << " phi=" << std::atan2(point.y, point.x) << std::endl;
      }
    }
    const int crossing = track->get_crossing();
    const int color = crossingColor(crossing);
    const int individualColor = trackColor(i);
    if (m_drawMeasurements) for (const auto& point : points)
    {
      if (!visible(point, m_zmin, m_zmax, m_xymax)) continue;
      allMarkers.push_back(marker(point, color));
      markersByCrossing[crossing].push_back(marker(point, individualColor));
    }
    const bool straight = m_useStraightLineTracks || std::abs(track->get_charge() * m_magneticFieldTesla) < 1.e-12;
    const unsigned int attemptsBefore = propagationAttempted;
    const unsigned int failuresBefore = propagationFailed;
    const unsigned int samplesBefore = trajectorySamples;
    if (m_drawFittedTrajectory)
    {
      if (auto* line = fittedTrajectory(*track, points, m_field, straight, m_magneticFieldTesla,
                                        m_zmin, m_zmax, m_xymax, color,
                                        propagationAttempted, propagationFailed, trajectorySamples))
      {
        allLines.push_back(line);
        auto* individualLine = dynamic_cast<TPolyLine3D*>(line->Clone());
        if (individualLine)
        {
          individualLine->SetLineColor(individualColor);
          linesByCrossing[crossing].push_back(individualLine);
        }
      }
    }
    if (Verbosity() >= 10 && !dumpedTrack)
    {
      unsigned int nTpc = 0, nIntt = 0, nMvtx = 0;
      double minRadius = std::numeric_limits<double>::max(), maxRadius = 0.;
      for (const auto& point : points)
      {
        if (point.detector == TrkrDefs::tpcId) ++nTpc;
        else if (point.detector == TrkrDefs::inttId) ++nIntt;
        else if (point.detector == TrkrDefs::mvtxId) ++nMvtx;
        const double radius = std::hypot(point.x, point.y);
        minRadius = std::min(minRadius, radius);
        maxRadius = std::max(maxRadius, radius);
      }
      std::cout << Name() << " display_track_summary parent=" << track->get_parent_tpc_track_id()
                << " tpc_markers=" << nTpc << " intt_markers=" << nIntt
                << " mvtx_markers=" << nMvtx
                << " fitted_trajectory_samples=" << (trajectorySamples - samplesBefore)
                << " min_radius=" << minRadius << " max_radius=" << maxRadius
                << " propagation_attempted=" << (propagationAttempted - attemptsBefore)
                << " propagation_failed=" << (propagationFailed - failuresBefore) << std::endl;
      dumpedTrack = true;
    }
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
  std::cout << Name() << " - saved event " << m_event << " selected_tracks=" << selected
            << " display_propagation_attempted=" << propagationAttempted
            << " display_propagation_failed=" << propagationFailed
            << " fitted_trajectory_samples=" << trajectorySamples << std::endl;
  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackDisplay::End(PHCompositeNode*)
{
  if (m_outfile) { m_outfile->Close(); delete m_outfile; m_outfile = nullptr; }
  return Fun4AllReturnCodes::EVENT_OK;
}
