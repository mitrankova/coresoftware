#include "Full_PolyTrackResiduals.h"

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

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <vector>

namespace
{
  using NativeState = std::array<double, TpcTrackKalmanFitter::StateDim>;

  struct Measurement
  {
    TrkrDefs::cluskey key{TrkrDefs::CLUSKEYMAX};
    TrkrDefs::TrkrId detector{TrkrDefs::tpcId};
    unsigned int layer{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  bool finiteState(const NativeState& state)
  {
    return std::all_of(state.begin(), state.end(),
                       [](const double value) { return std::isfinite(value); });
  }

  double wrapPhi(double phi)
  {
    const double pi = std::acos(-1.0);
    while (phi > pi) phi -= 2.0 * pi;
    while (phi <= -pi) phi += 2.0 * pi;
    return phi;
  }

  std::vector<NativeState> sampleFinalFieldTrajectory(
      const NativeState& initialState, const PHField* field,
      const double minimumRadius, const double maximumRadius, const double stepSize)
  {
    std::vector<NativeState> trajectory;
    if (!field || !finiteState(initialState) || !(stepSize > 0.0)) return trajectory;

    TpcKalmanConfig config;
    config.magnetic_field = field;
    config.analytic_uniform_propagation = false;
    const auto sampleBranch = [&](const int radialDirection)
    {
      std::vector<NativeState> branch{initialState};
      const double radius0 = std::hypot(initialState[TpcTrackKalmanFitter::X],
                                        initialState[TpcTrackKalmanFitter::Y]);
      const auto plus = TpcTrackKalmanFitter::propagate_state(initialState, stepSize, config);
      const auto minus = TpcTrackKalmanFitter::propagate_state(initialState, -stepSize, config);
      const auto advance = [&](const NativeState& state)
      {
        return finiteState(state)
            ? radialDirection * (std::hypot(state[TpcTrackKalmanFitter::X],
                                             state[TpcTrackKalmanFitter::Y]) - radius0)
            : -std::numeric_limits<double>::infinity();
      };
      const double plusAdvance = advance(plus);
      const double minusAdvance = advance(minus);
      if (!(std::max(plusAdvance, minusAdvance) > 1.e-6)) return branch;

      const double step = plusAdvance >= minusAdvance ? stepSize : -stepSize;
      NativeState state = initialState;
      double previousRadius = radius0;
      for (unsigned int index = 0; index < 4000; ++index)
      {
        if ((radialDirection > 0 && previousRadius >= maximumRadius) ||
            (radialDirection < 0 && previousRadius <= minimumRadius)) break;
        const auto next = TpcTrackKalmanFitter::propagate_state(state, step, config);
        if (!finiteState(next)) break;
        const double radius = std::hypot(next[TpcTrackKalmanFitter::X],
                                         next[TpcTrackKalmanFitter::Y]);
        if (!std::isfinite(radius) || radialDirection * (radius - previousRadius) < -1.e-4) break;
        branch.push_back(next);
        state = next;
        previousRadius = radius;
      }
      return branch;
    };

    auto outward = sampleBranch(+1);
    auto inward = sampleBranch(-1);
    trajectory.reserve(outward.size() + inward.size());
    trajectory.insert(trajectory.end(), outward.begin(), outward.end());
    trajectory.insert(trajectory.end(), std::next(inward.begin()), inward.end());
    return trajectory;
  }

  bool stateAtRadius(const std::vector<NativeState>& trajectory,
                     const double targetRadius, const double referenceZ, NativeState& output)
  {
    double bestDz = std::numeric_limits<double>::max();
    bool found = false;
    for (std::size_t index = 1; index < trajectory.size(); ++index)
    {
      const auto& first = trajectory[index - 1];
      const auto& second = trajectory[index];
      const double r1 = std::hypot(first[TpcTrackKalmanFitter::X], first[TpcTrackKalmanFitter::Y]);
      const double r2 = std::hypot(second[TpcTrackKalmanFitter::X], second[TpcTrackKalmanFitter::Y]);
      if ((targetRadius - r1) * (targetRadius - r2) > 0.0) continue;
      const double dr = r2 - r1;
      const double fraction = std::fabs(dr) > 1.e-12
          ? std::clamp((targetRadius - r1) / dr, 0.0, 1.0) : 0.0;
      NativeState candidate{};
      for (unsigned int component = 0; component < candidate.size(); ++component)
        candidate[component] = first[component] + fraction * (second[component] - first[component]);
      const double dz = std::fabs(candidate[TpcTrackKalmanFitter::Z] - referenceZ);
      if (dz < bestDz)
      {
        bestDz = dz;
        output = candidate;
        found = true;
      }
    }
    return found;
  }
}  // namespace

Full_PolyTrackResiduals::Full_PolyTrackResiduals(const std::string& name,
                                                 const std::string& outfilename)
  : SubsysReco(name)
  , m_outfilename(outfilename)
{
}

Full_PolyTrackResiduals::~Full_PolyTrackResiduals()
{
  delete m_outfile;
}

int Full_PolyTrackResiduals::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open output file " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_tree = new TTree("residuals", "Full poly track cluster residuals");
  m_tree->Branch("event", &m_event, "event/i");
  m_tree->Branch("full_track_id", &m_fullTrackId, "full_track_id/i");
  m_tree->Branch("parent_tpc_track_id", &m_parentTpcTrackId, "parent_tpc_track_id/i");
  m_tree->Branch("source_assembled_track_id", &m_sourceAssembledTrackId, "source_assembled_track_id/i");
  m_tree->Branch("crossing", &m_crossing, "crossing/S");
  m_tree->Branch("ntpc_clusters", &m_ntpcClusters, "ntpc_clusters/i");
  m_tree->Branch("nintt_clusters", &m_ninttClusters, "nintt_clusters/i");
  m_tree->Branch("nmvtx_clusters", &m_nmvtxClusters, "nmvtx_clusters/i");
  m_tree->Branch("fit_status", &m_fitStatus, "fit_status/I");
  m_tree->Branch("pt", &m_pt, "pt/D");
  m_tree->Branch("px", &m_px, "px/D");
  m_tree->Branch("py", &m_py, "py/D");
  m_tree->Branch("pz", &m_pz, "pz/D");
  m_tree->Branch("eta", &m_eta, "eta/D");
  m_tree->Branch("theta", &m_theta, "theta/D");
  m_tree->Branch("charge", &m_charge, "charge/D");
  m_tree->Branch("chi2", &m_chi2, "chi2/D");
  m_tree->Branch("ndf", &m_ndf, "ndf/D");
  m_tree->Branch("quality", &m_quality, "quality/D");
  m_tree->Branch("cluster_key", &m_clusterKey);
  m_tree->Branch("detector", &m_detector);
  m_tree->Branch("layer", &m_layer);
  m_tree->Branch("cluster_x", &m_clusterX);
  m_tree->Branch("cluster_y", &m_clusterY);
  m_tree->Branch("cluster_z", &m_clusterZ);
  m_tree->Branch("cluster_r", &m_clusterR);
  m_tree->Branch("cluster_phi", &m_clusterPhi);
  m_tree->Branch("state_x", &m_stateX);
  m_tree->Branch("state_y", &m_stateY);
  m_tree->Branch("state_z", &m_stateZ);
  m_tree->Branch("state_r", &m_stateR);
  m_tree->Branch("state_phi", &m_statePhi);
  m_tree->Branch("delta_phi", &m_deltaPhi);
  m_tree->Branch("residual_rphi", &m_residualRPhi);
  m_tree->Branch("residual_z", &m_residualZ);
  return Fun4AllReturnCodes::EVENT_OK;
}

bool Full_PolyTrackResiduals::getNodes(PHCompositeNode* topNode)
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

void Full_PolyTrackResiduals::resetTreeValues()
{
  m_event = m_evt;
  m_fullTrackId = m_parentTpcTrackId = m_sourceAssembledTrackId = 0;
  m_crossing = 0;
  m_ntpcClusters = m_ninttClusters = m_nmvtxClusters = 0;
  m_fitStatus = 0;
  m_pt = m_px = m_py = m_pz = m_eta = m_theta = m_charge = m_chi2 = m_ndf =
      m_quality = std::numeric_limits<double>::quiet_NaN();
  m_clusterKey.clear(); m_detector.clear(); m_layer.clear();
  m_clusterX.clear(); m_clusterY.clear(); m_clusterZ.clear(); m_clusterR.clear(); m_clusterPhi.clear();
  m_stateX.clear(); m_stateY.clear(); m_stateZ.clear(); m_stateR.clear(); m_statePhi.clear();
  m_deltaPhi.clear(); m_residualRPhi.clear(); m_residualZ.clear();
}

int Full_PolyTrackResiduals::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  if (!getNodes(topNode) || !m_tree) return Fun4AllReturnCodes::EVENT_OK;

  std::map<TrkrDefs::cluskey, const Tpc_PolyCluster*> tpcByKey;
  for (unsigned int index = 0; index < m_tpcClusters->size(); ++index)
  {
    const auto* cluster = m_tpcClusters->get_cluster(index);
    if (cluster && cluster->isValid()) tpcByKey[cluster->get_trkr_cluster_key()] = cluster;
  }

  unsigned int tracksFilled = 0;
  unsigned int residualsFilled = 0;
  for (unsigned int index = 0; index < m_fullTracks->size(); ++index)
  {
    const auto* track = m_fullTracks->get_track(index);
    if (!track || !track->isValid() || !track->has_final_native_state()) continue;
    NativeState initialState{};
    for (unsigned int component = 0; component < initialState.size(); ++component)
      initialState[component] = track->get_final_native_state(component);
    if (!finiteState(initialState)) continue;
    const double qOverPt = initialState[TpcTrackKalmanFitter::QOverPt];
    if (!(std::fabs(qOverPt) > 1.e-12)) continue;
    const double pt = 1.0 / std::fabs(qOverPt);
    if (pt < m_minPt || pt > m_maxPt) continue;

    std::vector<Measurement> measurements;
    for (const auto key : track->get_tpc_cluster_keys())
    {
      const auto found = tpcByKey.find(key);
      if (found == tpcByKey.end()) continue;
      const auto* cluster = found->second;
      measurements.push_back({key, TrkrDefs::tpcId, TrkrDefs::getLayer(key),
                              cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z()});
    }
    for (const auto key : track->get_silicon_cluster_keys())
    {
      auto* cluster = m_trkrClusters->findCluster(key);
      if (!cluster) continue;
      const auto global = m_actsGeometry->getGlobalPosition(key, cluster);
      measurements.push_back({key, static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(key)),
                              TrkrDefs::getLayer(key), global.x(), global.y(), global.z()});
    }
    if (measurements.empty()) continue;

    double minimumRadius = std::numeric_limits<double>::max();
    double maximumRadius = 0.0;
    for (const auto& measurement : measurements)
    {
      const double radius = std::hypot(measurement.x, measurement.y);
      minimumRadius = std::min(minimumRadius, radius);
      maximumRadius = std::max(maximumRadius, radius);
    }
    minimumRadius = std::max(0.0, minimumRadius - m_propagationStep);
    maximumRadius += m_propagationStep;
    const auto trajectory = sampleFinalFieldTrajectory(initialState, m_field, minimumRadius,
                                                       maximumRadius, m_propagationStep);
    if (trajectory.size() < 2) continue;

    resetTreeValues();
    m_fullTrackId = track->get_track_id();
    m_parentTpcTrackId = track->get_parent_tpc_track_id();
    m_sourceAssembledTrackId = track->get_source_assembled_track_id();
    m_crossing = track->get_crossing();
    m_fitStatus = track->get_fit_status();
    m_pt = pt; m_px = track->get_px(); m_py = track->get_py(); m_pz = track->get_pz();
    m_eta = std::isfinite(m_pz) ? std::asinh(m_pz / pt) : std::numeric_limits<double>::quiet_NaN();
    m_theta = std::isfinite(m_pz) ? std::atan2(pt, m_pz) : std::numeric_limits<double>::quiet_NaN();
    m_charge = track->get_charge(); m_chi2 = track->get_chi2(); m_ndf = track->get_ndf();
    m_quality = std::isfinite(m_chi2) && std::isfinite(m_ndf) && m_ndf > 0.0
        ? m_chi2 / m_ndf : std::numeric_limits<double>::quiet_NaN();

    for (const auto& measurement : measurements)
    {
      if (measurement.detector == TrkrDefs::tpcId) ++m_ntpcClusters;
      else if (measurement.detector == TrkrDefs::inttId) ++m_ninttClusters;
      else if (measurement.detector == TrkrDefs::mvtxId) ++m_nmvtxClusters;
    }

    for (const auto& measurement : measurements)
    {
      NativeState state{};
      const double clusterRadius = std::hypot(measurement.x, measurement.y);
      if (!stateAtRadius(trajectory, clusterRadius, measurement.z, state)) continue;
      const double clusterPhi = std::atan2(measurement.y, measurement.x);
      const double stateX = state[TpcTrackKalmanFitter::X];
      const double stateY = state[TpcTrackKalmanFitter::Y];
      const double stateZ = state[TpcTrackKalmanFitter::Z];
      const double statePhi = std::atan2(stateY, stateX);
      const double deltaPhi = wrapPhi(clusterPhi - statePhi);
      m_clusterKey.push_back(static_cast<std::uint64_t>(measurement.key));
      m_detector.push_back(static_cast<unsigned int>(measurement.detector));
      m_layer.push_back(measurement.layer);
      m_clusterX.push_back(measurement.x); m_clusterY.push_back(measurement.y);
      m_clusterZ.push_back(measurement.z); m_clusterR.push_back(clusterRadius);
      m_clusterPhi.push_back(clusterPhi); m_stateX.push_back(stateX); m_stateY.push_back(stateY);
      m_stateZ.push_back(stateZ); m_stateR.push_back(std::hypot(stateX, stateY));
      m_statePhi.push_back(statePhi); m_deltaPhi.push_back(deltaPhi);
      m_residualRPhi.push_back(clusterRadius * deltaPhi); m_residualZ.push_back(measurement.z - stateZ);
      ++residualsFilled;
    }
    if (!m_clusterKey.empty())
    {
      m_tree->Fill();
      ++tracksFilled;
    }
  }
  if (Verbosity() > 0)
    std::cout << Name() << "::process_event - event " << m_evt
              << " full_tracks=" << m_fullTracks->size() << " tracks_filled=" << tracksFilled
              << " residuals=" << residualsFilled << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackResiduals::End(PHCompositeNode*)
{
  if (m_outfile)
  {
    m_outfile->cd();
    if (m_tree) m_tree->Write();
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
    m_tree = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
