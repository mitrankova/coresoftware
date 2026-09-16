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
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <utility>
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

  class NamedPolyLine3D : public TPolyLine3D
  {
   public:
    NamedPolyLine3D(const int count, std::string name)
      : TPolyLine3D(count)
      , m_name(std::move(name))
    {
    }

    const char* GetName() const override { return m_name.c_str(); }

   private:
    std::string m_name;
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

  using NativeState = std::array<double, TpcTrackKalmanFitter::StateDim>;

  bool finiteState(const NativeState& state)
  {
    return std::all_of(state.begin(), state.end(), [](const double value) { return std::isfinite(value); });
  }

  TPolyLine3D* fieldTrajectory(const Full_PolyTrack& track,
                               const NativeState& initialState,
                               const std::vector<Point>& measurements,
                               const PHField* field,
                               const char* stateName,
                               const double zmin,
                               const double zmax,
                               const double xymax,
                               const int color,
                               const int verbosity,
                               unsigned int& attempted,
                               unsigned int& succeeded,
                               unsigned int& failed,
                               unsigned int& reportedFailures,
                               unsigned int& samples)
  {
    if (!field || !finiteState(initialState) || measurements.empty()) return nullptr;

    double maximumMeasurementRadius = 0.;
    for (const auto& point : measurements)
    {
      maximumMeasurementRadius = std::max(maximumMeasurementRadius, std::hypot(point.x, point.y));
    }
    const double minimumRadius = 2.0;
    const double maximumRadius = std::min(xymax, std::max(80.0, maximumMeasurementRadius + 1.0));

    TpcKalmanConfig config;
    config.magnetic_field = field;
    config.analytic_uniform_propagation = false;

    const auto propagate = [&](const NativeState& input, const double ds)
    {
      ++attempted;
      const auto output = TpcTrackKalmanFitter::propagate_state(input, ds, config);
      if (finiteState(output))
      {
        ++succeeded;
      }
      else
      {
        ++failed;
      }
      return output;
    };

    const auto reportFailure = [&](const char* reason, const double targetRadius,
                                   const NativeState& state)
    {
      if (verbosity < 5 || reportedFailures >= 5) return;
      ++reportedFailures;
      std::cout << "Full_PolyTrackDisplay propagation_failure parent_track_id="
                << track.get_parent_tpc_track_id() << " crossing=" << track.get_crossing()
                << " state=" << stateName << " target_radius=" << targetRadius
                << " status=" << reason
                << " x=" << state[TpcTrackKalmanFitter::X]
                << " y=" << state[TpcTrackKalmanFitter::Y]
                << " z=" << state[TpcTrackKalmanFitter::Z]
                << " Phi=" << state[TpcTrackKalmanFitter::Phi]
                << " QOverPt=" << state[TpcTrackKalmanFitter::QOverPt]
                << " TanLambda=" << state[TpcTrackKalmanFitter::TanLambda] << std::endl;
    };

    const auto sampleBranch = [&](const int radialDirection)
    {
      std::vector<NativeState> branch;
      branch.push_back(initialState);
      const double radius0 = std::hypot(initialState[TpcTrackKalmanFitter::X],
                                        initialState[TpcTrackKalmanFitter::Y]);
      const auto plus = propagate(initialState, 0.5);
      const auto minus = propagate(initialState, -0.5);
      const double plusAdvance = finiteState(plus)
          ? radialDirection * (std::hypot(plus[TpcTrackKalmanFitter::X], plus[TpcTrackKalmanFitter::Y]) - radius0)
          : -std::numeric_limits<double>::infinity();
      const double minusAdvance = finiteState(minus)
          ? radialDirection * (std::hypot(minus[TpcTrackKalmanFitter::X], minus[TpcTrackKalmanFitter::Y]) - radius0)
          : -std::numeric_limits<double>::infinity();
      if (!(std::max(plusAdvance, minusAdvance) > 1.e-5))
      {
        reportFailure("no_radial_direction", radialDirection > 0 ? maximumRadius : minimumRadius,
                      initialState);
        return branch;
      }

      const double step = plusAdvance >= minusAdvance ? 0.5 : -0.5;
      NativeState state = initialState;
      double previousRadius = radius0;
      for (unsigned int index = 0; index < 2000; ++index)
      {
        const double targetRadius = radialDirection > 0 ? maximumRadius : minimumRadius;
        if ((radialDirection > 0 && previousRadius >= maximumRadius) ||
            (radialDirection < 0 && previousRadius <= minimumRadius)) break;
        const auto next = propagate(state, step);
        if (!finiteState(next))
        {
          reportFailure("non_finite_state", targetRadius, state);
          break;
        }
        const double radius = std::hypot(next[TpcTrackKalmanFitter::X],
                                         next[TpcTrackKalmanFitter::Y]);
        const double dx = next[TpcTrackKalmanFitter::X] - state[TpcTrackKalmanFitter::X];
        const double dy = next[TpcTrackKalmanFitter::Y] - state[TpcTrackKalmanFitter::Y];
        const double dz = next[TpcTrackKalmanFitter::Z] - state[TpcTrackKalmanFitter::Z];
        if (!std::isfinite(radius) || radialDirection * (radius - previousRadius) < -1.e-4 ||
            dx * dx + dy * dy + dz * dz < 1.e-12)
        {
          --succeeded;
          ++failed;
          reportFailure("trajectory_sanity_check", targetRadius, next);
          break;
        }
        branch.push_back(next);
        state = next;
        previousRadius = radius;
      }
      return branch;
    };

    auto outward = sampleBranch(+1);
    auto inward = sampleBranch(-1);
    std::vector<NativeState> ordered;
    ordered.reserve(outward.size() + inward.size());
    for (auto iterator = outward.rbegin(); iterator != outward.rend(); ++iterator)
    {
      ordered.push_back(*iterator);
    }
    for (auto iterator = std::next(inward.begin()); iterator != inward.end(); ++iterator)
    {
      ordered.push_back(*iterator);
    }

    std::vector<NativeState> visibleStates;
    for (const auto& state : ordered)
    {
      if (visible({state[TpcTrackKalmanFitter::X], state[TpcTrackKalmanFitter::Y],
                   state[TpcTrackKalmanFitter::Z]}, zmin, zmax, xymax))
      {
        visibleStates.push_back(state);
      }
    }
    if (visibleStates.size() < 2) return nullptr;

    if (verbosity >= 10)
    {
      std::cout << "Full_PolyTrackDisplay native_state parent=" << track.get_parent_tpc_track_id()
                << " crossing=" << track.get_crossing() << " source=" << stateName
                << " frame=global_detector"
                << " x=" << initialState[TpcTrackKalmanFitter::X]
                << " y=" << initialState[TpcTrackKalmanFitter::Y]
                << " z=" << initialState[TpcTrackKalmanFitter::Z]
                << " Phi=" << initialState[TpcTrackKalmanFitter::Phi]
                << " QOverPt=" << initialState[TpcTrackKalmanFitter::QOverPt]
                << " TanLambda=" << initialState[TpcTrackKalmanFitter::TanLambda] << std::endl;
      for (unsigned int index = 0; index < visibleStates.size(); ++index)
      {
        const auto& state = visibleStates[index];
        std::cout << "  sample=" << index
                  << " radius=" << std::hypot(state[TpcTrackKalmanFitter::X],
                                               state[TpcTrackKalmanFitter::Y])
                  << " x=" << state[TpcTrackKalmanFitter::X]
                  << " y=" << state[TpcTrackKalmanFitter::Y]
                  << " z=" << state[TpcTrackKalmanFitter::Z]
                  << " phi=" << state[TpcTrackKalmanFitter::Phi]
                  << " QOverPt=" << state[TpcTrackKalmanFitter::QOverPt]
                  << " TanLambda=" << state[TpcTrackKalmanFitter::TanLambda] << std::endl;
      }
    }

    samples += visibleStates.size();
    auto* line = new NamedPolyLine3D(
        visibleStates.size(),
        std::format("{}_parent_{}", stateName, track.get_parent_tpc_track_id()));
    for (unsigned int index = 0; index < visibleStates.size(); ++index)
    {
      const auto& state = visibleStates[index];
      line->SetPoint(index, state[TpcTrackKalmanFitter::Z],
                     state[TpcTrackKalmanFitter::X], state[TpcTrackKalmanFitter::Y]);
    }
    line->SetLineColor(color);
    line->SetLineWidth(stateName == std::string("final_field_trajectory") ? 3 : 2);
    line->SetLineStyle(stateName == std::string("final_field_trajectory") ? 1 : 2);
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
  if (!m_field) std::cerr << Name() << " - magnetic field input missing" << std::endl;
  return m_fullTracks && m_tpcClusters && m_trkrClusters && m_actsGeometry && m_field;
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
  unsigned int propagationSucceeded = 0;
  unsigned int propagationFailed = 0;
  unsigned int reportedFailures = 0;
  unsigned int trajectorySamples = 0;
  bool dumpedTrack = false;
  for (unsigned int i = 0; i < m_fullTracks->size(); ++i)
  {
    const auto* track = m_fullTracks->get_track(i);
    if (!track || (!track->has_final_native_state() && !track->has_fast_native_state()) ||
        track->get_n_mvtx() < m_minMvtxHits || track->get_n_intt() < m_minInttHits) continue;
    const bool useFinalForSelection = track->has_final_native_state();
    const double selectionQOverPt = useFinalForSelection
        ? track->get_final_native_state(TpcTrackKalmanFitter::QOverPt)
        : track->get_fast_native_state(TpcTrackKalmanFitter::QOverPt);
    if (!(std::abs(selectionQOverPt) > 1.e-12) ||
        1. / std::abs(selectionQOverPt) < m_minTrackPt) continue;
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
    const unsigned int attemptsBefore = propagationAttempted;
    const unsigned int successesBefore = propagationSucceeded;
    const unsigned int failuresBefore = propagationFailed;
    const unsigned int samplesBefore = trajectorySamples;
    const auto drawTrajectory = [&](const NativeState& state, const char* stateName)
    {
      if (auto* line = fieldTrajectory(*track, state, points, m_field, stateName,
                                       m_zmin, m_zmax, m_xymax, color, Verbosity(),
                                       propagationAttempted, propagationSucceeded,
                                       propagationFailed, reportedFailures, trajectorySamples))
      {
        allLines.push_back(line);
        auto* individualLine = new NamedPolyLine3D(
            *static_cast<NamedPolyLine3D*>(line));
        if (individualLine)
        {
          individualLine->SetLineColor(individualColor);
          linesByCrossing[crossing].push_back(individualLine);
        }
      }
    };
    if (m_drawFinalFieldTrajectory && track->has_final_native_state())
    {
      NativeState state{};
      for (unsigned int index = 0; index < state.size(); ++index)
        state[index] = track->get_final_native_state(index);
      drawTrajectory(state, "final_field_trajectory");
    }
    if (m_drawFastFieldTrajectory && track->has_fast_native_state())
    {
      NativeState state{};
      for (unsigned int index = 0; index < state.size(); ++index)
        state[index] = track->get_fast_native_state(index);
      drawTrajectory(state, "fast_field_trajectory");
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
                << " crossing=" << track->get_crossing()
                << " state_frame=global_detector"
                << " tpc_markers=" << nTpc << " intt_markers=" << nIntt
                << " mvtx_markers=" << nMvtx
                << " field_trajectory_samples=" << (trajectorySamples - samplesBefore)
                << " min_radius=" << minRadius << " max_radius=" << maxRadius
                << " propagation_attempts=" << (propagationAttempted - attemptsBefore)
                << " propagation_success=" << (propagationSucceeded - successesBefore)
                << " propagation_failures=" << (propagationFailed - failuresBefore) << std::endl;
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
    for (auto* line : lines) if (line)
    {
      line->Write();
      line->Draw("same");
    }
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
            << " display_propagation_attempts=" << propagationAttempted
            << " display_propagation_success=" << propagationSucceeded
            << " display_propagation_failures=" << propagationFailed
            << " field_trajectory_samples=" << trajectorySamples << std::endl;
  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackDisplay::End(PHCompositeNode*)
{
  if (m_outfile) { m_outfile->Close(); delete m_outfile; m_outfile = nullptr; }
  return Fun4AllReturnCodes::EVENT_OK;
}
