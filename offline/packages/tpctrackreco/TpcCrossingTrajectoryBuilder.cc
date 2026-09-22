#include "TpcCrossingTrajectoryBuilder.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include "TpcCrossingClusterPosition.h"
#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
#include "TpcDriftPolylineLookup.h"
#include "TpcTrackKalmanFitter.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyClusterContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phfield/PHFieldUtility.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TpcDefs.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace
{
  enum class FullFitOracleStatus
  {
    Valid,
    FitFailed,
    NonFinite,
    BadCovariance,
    DiscontinuousSolution
  };

  const char* oracleStatusName(const FullFitOracleStatus status)
  {
    switch (status)
    {
    case FullFitOracleStatus::Valid: return "Valid";
    case FullFitOracleStatus::FitFailed: return "FitFailed";
    case FullFitOracleStatus::NonFinite: return "NonFinite";
    case FullFitOracleStatus::BadCovariance: return "BadCovariance";
    case FullFitOracleStatus::DiscontinuousSolution: return "DiscontinuousSolution";
    }
    return "Unknown";
  }

  FullFitOracleStatus classifyOracle(const FastFieldTrackFitter::Result& reference,
                                     const FastFieldTrackFitter::Result& candidate,
                                     double& continuityMaxPull)
  {
    continuityMaxPull = 0.;
    if (!candidate.fitSuccess || !candidate.valid) return FullFitOracleStatus::FitFailed;
    if (!std::isfinite(candidate.chi2)) return FullFitOracleStatus::NonFinite;
    for (const double value : candidate.nativeState)
      if (!std::isfinite(value)) return FullFitOracleStatus::NonFinite;
    for (const double value : candidate.covariance)
      if (!std::isfinite(value)) return FullFitOracleStatus::NonFinite;
    for (unsigned int index = 0; index < FastFieldTrackFitter::StateSize; ++index)
      if (!(candidate.covariance[7 * index] > 0.)) return FullFitOracleStatus::BadCovariance;

    // This is only a QA continuity diagnostic. A 50-sigma jump in any of the
    // angular/curvature parameters is reported, never clipped or rejected by
    // reconstruction. It combines both fit covariances rather than imposing a
    // lone absolute QOverPt cut.
    const std::array<unsigned int, 3> continuityIndices{{TpcTrackKalmanFitter::Phi,
                                                         TpcTrackKalmanFitter::QOverPt,
                                                         TpcTrackKalmanFitter::TanLambda}};
    for (const unsigned int index : continuityIndices)
    {
      double delta = candidate.nativeState[index] - reference.nativeState[index];
      if (index == TpcTrackKalmanFitter::Phi) delta = std::remainder(delta, 2. * M_PI);
      const double variance = reference.covariance[7 * index] + candidate.covariance[7 * index];
      if (!(variance > 0.) || !std::isfinite(variance)) return FullFitOracleStatus::BadCovariance;
      continuityMaxPull = std::max(continuityMaxPull, std::abs(delta) / std::sqrt(variance));
    }
    return continuityMaxPull > 50. ? FullFitOracleStatus::DiscontinuousSolution
                                   : FullFitOracleStatus::Valid;
  }

  double relativeDifference(const double lhs, const double rhs)
  {
    constexpr double floor = 1.e-12;
    return std::abs(lhs - rhs) / std::max({std::abs(lhs), std::abs(rhs), floor});
  }

  void validateLongitudinalJacobian(
      const unsigned int parentTrackId,
      const short int referenceCrossing,
      const FastFieldTrackFitter::Result& fit,
      const std::map<TrkrDefs::cluskey, const Tpc_PolyCluster*>& clusters,
      const double epsilonZ,
      const double epsilonTanLambda)
  {
    if (!(epsilonZ > 0.) || !(epsilonTanLambda > 0.) || fit.measurements.empty()) return;

    auto zPlus = fit.nativeState;
    auto zMinus = fit.nativeState;
    auto tanLPlus = fit.nativeState;
    auto tanLMinus = fit.nativeState;
    zPlus[TpcTrackKalmanFitter::Z] += epsilonZ;
    zMinus[TpcTrackKalmanFitter::Z] -= epsilonZ;
    tanLPlus[TpcTrackKalmanFitter::TanLambda] += epsilonTanLambda;
    tanLMinus[TpcTrackKalmanFitter::TanLambda] -= epsilonTanLambda;

    const std::size_t count = fit.measurements.size();
    const std::set<std::size_t> printed{{0, count / 4, count / 2, (3 * count) / 4, count - 1}};
    double maxAbsZ = 0., maxRelZ = 0., maxAbsTanL = 0., maxRelTanL = 0.;
    bool finite = true;
    for (std::size_t i = 0; i < count; ++i)
    {
      const double path = i < fit.pathS.size() ? fit.pathS[i] - fit.pathS.front() : 0.;
      const auto predictedZPlus = TpcTrackKalmanFitter::propagate_state(zPlus, path, fit.propagationConfig);
      const auto predictedZMinus = TpcTrackKalmanFitter::propagate_state(zMinus, path, fit.propagationConfig);
      const auto predictedTanLPlus = TpcTrackKalmanFitter::propagate_state(tanLPlus, path, fit.propagationConfig);
      const auto predictedTanLMinus = TpcTrackKalmanFitter::propagate_state(tanLMinus, path, fit.propagationConfig);
      std::array<double, 3> fdZ{};
      std::array<double, 3> fdTanL{};
      for (unsigned int row = 0; row < 3; ++row)
      {
        fdZ[row] = (predictedZPlus[row] - predictedZMinus[row]) / (2. * epsilonZ);
        fdTanL[row] = (predictedTanLPlus[row] - predictedTanLMinus[row]) / (2. * epsilonTanLambda);
        const double jacZ = fit.measurements[i].jacobian[6 * row + TpcTrackKalmanFitter::Z];
        const double jacTanL = fit.measurements[i].jacobian[6 * row + TpcTrackKalmanFitter::TanLambda];
        finite = finite && std::isfinite(fdZ[row]) && std::isfinite(fdTanL[row]) &&
                 std::isfinite(jacZ) && std::isfinite(jacTanL);
        maxAbsZ = std::max(maxAbsZ, std::abs(jacZ - fdZ[row]));
        maxRelZ = std::max(maxRelZ, relativeDifference(jacZ, fdZ[row]));
        maxAbsTanL = std::max(maxAbsTanL, std::abs(jacTanL - fdTanL[row]));
        maxRelTanL = std::max(maxRelTanL, relativeDifference(jacTanL, fdTanL[row]));
      }
      if (printed.count(i) != 0)
      {
        const auto cluster = clusters.find(fit.measurements[i].key);
        const int layer = cluster != clusters.end() && cluster->second->size_hits()
                              ? static_cast<int>(TrkrDefs::getLayer(cluster->second->get_hit_index(0).first))
                              : -1;
        const auto& position = fit.measurements[i].reference;
        std::cout << "LongitudinalJacobian parent_track_id=" << parentTrackId
                  << " reference_crossing=" << referenceCrossing
                  << " measurement_index=" << i << " layer=" << layer
                  << " radius=" << std::hypot(position[0], position[1])
                  << " xyz=" << position[0] << "," << position[1] << "," << position[2];
        static const std::array<const char*, 3> axes{{"x", "y", "z"}};
        for (unsigned int row = 0; row < 3; ++row)
        {
          const double jac = fit.measurements[i].jacobian[6 * row + TpcTrackKalmanFitter::Z];
          std::cout << " J_Z_" << axes[row] << "=" << jac
                    << " FD_Z_" << axes[row] << "=" << fdZ[row]
                    << " diff_Z_" << axes[row] << "=" << jac - fdZ[row];
        }
        for (unsigned int row = 0; row < 3; ++row)
        {
          const double jac = fit.measurements[i].jacobian[6 * row + TpcTrackKalmanFitter::TanLambda];
          std::cout << " J_TanL_" << axes[row] << "=" << jac
                    << " FD_TanL_" << axes[row] << "=" << fdTanL[row]
                    << " diff_TanL_" << axes[row] << "=" << jac - fdTanL[row];
        }
        std::cout << std::endl;
      }
    }
    std::cout << "LongitudinalJacobianSummary parent_track_id=" << parentTrackId
              << " reference_crossing=" << referenceCrossing
              << " epsilon_Z_cm=" << epsilonZ
              << " epsilon_TanLambda=" << epsilonTanLambda
              << " max_abs_diff_Z_column=" << maxAbsZ
              << " max_rel_diff_Z_column=" << maxRelZ
              << " max_abs_diff_TanLambda_column=" << maxAbsTanL
              << " max_rel_diff_TanLambda_column=" << maxRelTanL
              << " finite=" << finite << std::endl;
  }
}

TpcCrossingTrajectoryBuilder::TpcCrossingTrajectoryBuilder(const std::string& name) : SubsysReco(name) {}
int TpcCrossingTrajectoryBuilder::getNodes(PHCompositeNode* topNode)
{
  m_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_trackNodeName);
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_clusterNodeName);
  m_decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_lookup = TpcDriftPolylineLookup::get(topNode);
  m_field = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
  if (!m_tracks || !m_clusters || !m_decisions || !m_hits || !m_lookup || !m_field)
  {
    std::cerr << Name() << "::getNodes - missing track, cluster, decision, hit, drift lookup, or field input" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrajectoryBuilder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dst = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst) { dst = new PHCompositeNode("DST"); topNode->addNode(dst); }
  m_trajectories = findNode::getClass<TpcCrossingTrajectoryContainer>(topNode, m_outputNodeName);
  if (!m_trajectories) { m_trajectories = new TpcCrossingTrajectoryContainer; dst->addNode(new PHIODataNode<PHObject>(m_trajectories, m_outputNodeName, "PHObject")); }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrajectoryBuilder::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  m_fitter = std::make_unique<FastFieldTrackFitter>(m_field);
  return createNodes(topNode);
}
bool TpcCrossingTrajectoryBuilder::addSiliconStates(TpcCrossingTrajectory& trajectory, const FastFieldTrackFitter::Result& fit,
    const std::array<double, FastFieldTrackFitter::StateSize>& candidateState) const
{
  auto state = candidateState;
  for (int layer = static_cast<int>(m_trajectorySamplingRadii.size()) - 1; layer >= 0; --layer)
  {
    const double target = m_trajectorySamplingRadii[layer];
    bool valid = false;
    for (unsigned int step = 0; step < 1600; ++step)
    {
      const double radius = std::hypot(state[TpcTrackKalmanFitter::X], state[TpcTrackKalmanFitter::Y]);
      if (radius <= target + 0.08) { valid = true; break; }
      state = TpcTrackKalmanFitter::propagate_state(state, -0.25, fit.propagationConfig);
      if (!std::isfinite(state[TpcTrackKalmanFitter::X])) break;
    }
    TpcCrossingTrajectory::LayerState output;
    output.layer = static_cast<unsigned int>(layer);
    output.x = state[TpcTrackKalmanFitter::X]; output.y = state[TpcTrackKalmanFitter::Y]; output.z = state[TpcTrackKalmanFitter::Z];
    output.phi = std::atan2(output.y, output.x); output.valid = valid;
    trajectory.add_layer_state(output);
  }
  return true;
}
int TpcCrossingTrajectoryBuilder::process_event(PHCompositeNode*)
{
  const auto begin = std::chrono::steady_clock::now();
  m_trajectories->Reset();
  std::map<TrkrDefs::cluskey, const Tpc_PolyCluster*> byKey;
  for (unsigned int i = 0; i < m_clusters->size(); ++i) if (const auto* c = m_clusters->get_cluster(i)) byKey[c->get_trkr_cluster_key()] = c;
  unsigned int referenceFits = 0, deltaBuilds = 0;
  unsigned int longitudinalJacobianValidations = 0;
  unsigned int fullValidationAttempted = 0, fullValidationValid = 0;
  unsigned int fullValidationFailed = 0, fullValidationNonFinite = 0;
  unsigned int fullValidationBadCovariance = 0, fullValidationDiscontinuous = 0;
  unsigned int nearValidationAttempted = 0, nearValidationValid = 0;
  double referenceFitSeconds = 0.0, responseSeconds = 0.0, crossingSeconds = 0.0;
  unsigned int inputTracks = 0;
  unsigned int rejectedFitStatus = 0;
  unsigned int rejectedMinPt = 0;
  unsigned int rejectedMinTpcClusters = 0;
  unsigned int selectedTracks = 0;
  for (unsigned int i = 0; i < m_tracks->size(); ++i)
  {
    const auto* track = m_tracks->get_track(i);
    ++inputTracks;

    // Require a valid TPC PolyTrack fit.
    if (!track || !track->get_fit_status())
    {
      ++rejectedFitStatus;
      continue;
    }

    // Require at least 18 TPC clusters, i.e. ntpc_clusters > 17.
    const unsigned int nTpcClusters = track->size_cluster_keys();
    if (nTpcClusters < m_minTpcClusters)
    {
      ++rejectedMinTpcClusters;
      continue;
    }

    // Cheap TPC PolyTrack pT preselection.
    // This is intentionally applied before the expensive FastFieldTrackFitter.
    const double pt = std::hypot(track->get_px(), track->get_py());
    if (!std::isfinite(pt) || pt <= m_minPt)
    {
      ++rejectedMinPt;
      continue;
    }

    ++selectedTracks;

    const auto* decision =
        m_decisions->get_decision(track->get_source_assembled_track_id());
    if (!decision)
    {
      continue;
    }
    std::vector<const Tpc_PolyCluster*> clusters;
    for (const auto key : track->get_cluster_keys()) { const auto found = byKey.find(key); if (found != byKey.end()) clusters.push_back(found->second); }
    FastFieldTrackFitter::Result referenceFit;
    if (!m_fitter->fit(*track, clusters, referenceFit)) continue;
    ++referenceFits;
    referenceFitSeconds += referenceFit.fitSeconds;
    responseSeconds += referenceFit.responseSeconds;
    if (m_validateLongitudinalJacobian &&
        longitudinalJacobianValidations < m_longitudinalJacobianValidationTracks)
    {
      validateLongitudinalJacobian(track->get_track_id(), decision->get_reference_crossing(),
                                   referenceFit, byKey, m_longitudinalJacobianEpsilonZ,
                                   m_longitudinalJacobianEpsilonTanLambda);
      ++longitudinalJacobianValidations;
    }
    for (unsigned int ic = 0; ic < decision->get_number_of_candidates(); ++ic)
    {
      const auto* crossingCandidate = decision->get_candidate(ic);
      if (!crossingCandidate || !crossingCandidate->tpc_valid) continue;
      const auto crossingBegin = std::chrono::steady_clock::now();
      std::map<TrkrDefs::cluskey, std::array<double, 3>> displaced;
      bool positionsOk = true;
      for (const auto* cluster : clusters)
      {
        std::array<double, 3> position{};
        if (!TpcCrossingClusterPosition::get(*cluster, *m_lookup, crossingCandidate->crossing,
                                             decision->get_reference_crossing(), position))
        {
          positionsOk = false;
          break;
        }
        displaced[cluster->get_trkr_cluster_key()] = position;

        const double dx = position[0] - cluster->get_centroid_x();
        const double dy = position[1] - cluster->get_centroid_y();
        const double dz = position[2] - cluster->get_centroid_z();
        const double magnitude = std::sqrt(dx * dx + dy * dy + dz * dz);
        static bool dumpedIdentityFailure = false;
        if (!dumpedIdentityFailure && Verbosity() >= 10 &&
            crossingCandidate->crossing == decision->get_reference_crossing() && magnitude > 1.e-7)
        {
          dumpedIdentityFailure = true;
          const auto first = cluster->get_hit_index(0);
          std::vector<TpcCrossingClusterPosition::HitPosition> hitPositions;
          TpcCrossingClusterPosition::get(*cluster, *m_lookup, crossingCandidate->crossing,
                                          decision->get_reference_crossing(), position, &hitPositions);
          std::cout << Name() << " identity_failure parent_track_id=" << track->get_track_id()
                    << " cluster_id=" << cluster->get_cluster_id()
                    << " layer=" << TrkrDefs::getLayer(first.first) << " side=" << cluster->get_side()
                    << " raw_hits=" << cluster->size_hits()
                    << " reference_xyz=" << cluster->get_centroid_x() << "," << cluster->get_centroid_y() << "," << cluster->get_centroid_z()
                    << " candidate_xyz=" << position[0] << "," << position[1] << "," << position[2]
                    << " delta_xyz=" << dx << "," << dy << "," << dz << " delta_magnitude=" << magnitude << std::endl;
          for (const auto& hitPosition : hitPositions)
          {
            auto* hitset = m_hits->findHitSet(hitPosition.hitsetkey);
            auto* hit = hitset ? hitset->getHit(hitPosition.hitkey) : nullptr;
            std::cout << "  hitsetkey=" << hitPosition.hitsetkey << " hitkey=" << hitPosition.hitkey
                      << " pad=" << TpcDefs::getPad(hitPosition.hitkey) << " tbin=" << TpcDefs::getTBin(hitPosition.hitkey)
                      << " adc=" << (hit ? hit->getAdc() : 0.0)
                      << " p_ref=" << hitPosition.reference[0] << "," << hitPosition.reference[1] << "," << hitPosition.reference[2]
                      << " p_candidate=" << hitPosition.candidate[0] << "," << hitPosition.candidate[1] << "," << hitPosition.candidate[2]
                      << " delta_hit=" << hitPosition.delta[0] << "," << hitPosition.delta[1] << "," << hitPosition.delta[2] << std::endl;
          }
        }
      }
      if (!positionsOk || displaced.size() != clusters.size()) continue;
      const auto update = m_fitter->linearUpdate(referenceFit, displaced);
      if (!update.valid) continue;
      auto* trajectory = new TpcCrossingTrajectory;
      trajectory->set_parent_track_id(track->get_track_id());
      trajectory->set_source_assembled_track_id(track->get_source_assembled_track_id());
      trajectory->set_crossing(crossingCandidate->crossing);
      trajectory->set_reference_crossing(decision->get_reference_crossing());
      for (unsigned int k = 0; k < TpcCrossingTrajectory::StateSize; ++k) { trajectory->set_delta(k, update.delta[k]); trajectory->set_state(k, update.state[k]); }
      for (unsigned int row = 0; row < TpcCrossingTrajectory::StateSize; ++row) for (unsigned int col = 0; col < TpcCrossingTrajectory::StateSize; ++col) trajectory->set_covariance(row, col, referenceFit.covariance[row * TpcCrossingTrajectory::StateSize + col]);
      trajectory->set_linear_chi2(update.chi2);
      addSiliconStates(*trajectory, referenceFit, update.state);
      const unsigned int validationHash = 2654435761U * (track->get_track_id() + 1U) + static_cast<unsigned int>(crossingCandidate->crossing);
      const bool validate = m_validationFraction > 0.0 &&
                            static_cast<double>(validationHash % 1000000U) / 1000000.0 < std::min(1.0, m_validationFraction);
      if (validate)
      {
        std::vector<TpcTrackPoint> candidatePoints;
        candidatePoints.reserve(clusters.size());
        for (const auto* cluster : clusters)
        {
          const auto& position = displaced.at(cluster->get_trkr_cluster_key());
          TpcTrackPoint point;
          point.track_id = static_cast<int>(track->get_track_id());
          point.layer = cluster->size_hits() ? static_cast<int>(TrkrDefs::getLayer(cluster->get_hit_index(0).first)) : 0;
          point.position = {position[0], position[1], position[2]};
          point.momentum = {track->get_px(), track->get_py(), track->get_pz()};
          point.detector = TpcTrackPoint::Detector::Tpc;
          point.cluster_key = cluster->get_trkr_cluster_key();
          candidatePoints.push_back(point);
        }
        ++fullValidationAttempted;
        const int crossingOffset = static_cast<int>(crossingCandidate->crossing) -
                                   static_cast<int>(decision->get_reference_crossing());
        const bool nearReference = std::abs(crossingOffset) <= 3;
        if (nearReference) ++nearValidationAttempted;
        FastFieldTrackFitter::Result candidateFit;
        m_fitter->fitMeasurements(*track, candidatePoints, referenceFit.nativeState, candidateFit);
        double continuityMaxPull = 0.;
        const auto oracleStatus = classifyOracle(referenceFit, candidateFit, continuityMaxPull);
        switch (oracleStatus)
        {
        case FullFitOracleStatus::Valid:
          ++fullValidationValid;
          if (nearReference) ++nearValidationValid;
          break;
        case FullFitOracleStatus::FitFailed: ++fullValidationFailed; break;
        case FullFitOracleStatus::NonFinite: ++fullValidationNonFinite; break;
        case FullFitOracleStatus::BadCovariance: ++fullValidationBadCovariance; break;
        case FullFitOracleStatus::DiscontinuousSolution: ++fullValidationDiscontinuous; break;
        }
        std::cout << Name() << " validation parent_track_id=" << track->get_track_id()
                  << " reference_crossing=" << decision->get_reference_crossing()
                  << " candidate_crossing=" << crossingCandidate->crossing
                  << " crossing_offset=" << crossingOffset
                  << " near_reference=" << nearReference
                  << " oracle_status=" << oracleStatusName(oracleStatus)
                  << " full_fit_success=" << candidateFit.fitSuccess
                  << " full_chi2=" << candidateFit.chi2
                  << " full_ndf=" << candidateFit.ndf
                  << " full_chi2_ndf=" << (candidateFit.ndf > 0 ? candidateFit.chi2 / candidateFit.ndf : -1.)
                  << " n_measurements=" << candidateFit.nMeasurements
                  << " n_accepted=" << candidateFit.nAccepted
                  << " continuity_max_pull=" << continuityMaxPull;
        if (!candidateFit.fitSuccess)
        {
          std::cout << " fit_message=\"" << candidateFit.fitMessage << "\"" << std::endl;
        }
        else
        {
          static const std::array<const char*, 6> nativeNames{{"X", "Y", "Z", "Phi", "QOverPt", "TanLambda"}};
          for (unsigned int k = 0; k < TpcCrossingTrajectory::StateSize; ++k)
          {
            double fullDelta = candidateFit.nativeState[k] - referenceFit.nativeState[k];
            if (k == TpcTrackKalmanFitter::Phi) fullDelta = std::remainder(fullDelta, 2. * M_PI);
            std::cout << " delta" << nativeNames[k] << "_linear=" << update.delta[k]
                      << " delta" << nativeNames[k] << "_full=" << fullDelta;
          }
          const auto linearExternal = FastFieldTrackFitter::externalState(update.state);
          const auto roundTripNative = FastFieldTrackFitter::nativeState(linearExternal);
          double roundTripError = 0.0;
          for (unsigned int k = 0; k < TpcCrossingTrajectory::StateSize; ++k)
          {
            const double difference = k == TpcTrackKalmanFitter::Phi
                ? std::remainder(roundTripNative[k] - update.state[k], 2. * M_PI)
                : roundTripNative[k] - update.state[k];
            roundTripError = std::max(roundTripError, std::abs(difference));
          }
          assert(roundTripError < 1.e-10);
          std::cout << " deltaTheta_linear=" << linearExternal[4] - referenceFit.state[4]
                    << " deltaTheta_full=" << candidateFit.state[4] - referenceFit.state[4]
                    << " delta(q/p)_linear=" << linearExternal[5] - referenceFit.state[5]
                    << " delta(q/p)_full=" << candidateFit.state[5] - referenceFit.state[5]
                    << " A_condition=" << referenceFit.informationCondition
                    << " LDLT_ok=" << referenceFit.informationSolveOk
                    << " A_eigenvalues=";
          for (const auto value : referenceFit.informationEigenvalues) std::cout << value << ",";
          double qOverPtResponse2 = 0.0, tanLambdaResponse2 = 0.0;
          for (const auto& response : referenceFit.measurements)
            for (unsigned int column = 0; column < 3; ++column)
            {
              qOverPtResponse2 += response.response[3 * TpcTrackKalmanFitter::QOverPt + column] * response.response[3 * TpcTrackKalmanFitter::QOverPt + column];
              tanLambdaResponse2 += response.response[3 * TpcTrackKalmanFitter::TanLambda + column] * response.response[3 * TpcTrackKalmanFitter::TanLambda + column];
            }
          std::cout << " QOverPt_response_norm=" << std::sqrt(qOverPtResponse2)
                    << " TanLambda_response_norm=" << std::sqrt(tanLambdaResponse2)
                    << " reference_QOverPt=" << referenceFit.nativeState[TpcTrackKalmanFitter::QOverPt]
                    << " candidate_QOverPt=" << candidateFit.nativeState[TpcTrackKalmanFitter::QOverPt]
                    << " reference_QOverPt_variance=" << referenceFit.covariance[7 * TpcTrackKalmanFitter::QOverPt]
                    << " candidate_QOverPt_variance=" << candidateFit.covariance[7 * TpcTrackKalmanFitter::QOverPt]
                    << " reference_Phi=" << referenceFit.nativeState[TpcTrackKalmanFitter::Phi]
                    << " candidate_Phi=" << candidateFit.nativeState[TpcTrackKalmanFitter::Phi]
                    << " reference_TanLambda=" << referenceFit.nativeState[TpcTrackKalmanFitter::TanLambda]
                    << " candidate_TanLambda=" << candidateFit.nativeState[TpcTrackKalmanFitter::TanLambda]
                    << " reference_TanLambda_variance=" << referenceFit.covariance[7 * TpcTrackKalmanFitter::TanLambda]
                    << " candidate_TanLambda_variance=" << candidateFit.covariance[7 * TpcTrackKalmanFitter::TanLambda]
                    << " linear_chi2=" << update.chi2 << std::endl;
          if (crossingCandidate->crossing == decision->get_reference_crossing())
          {
            const double maxStateDelta = *std::max_element(update.delta.begin(), update.delta.end(),
                [](double lhs, double rhs) { return std::abs(lhs) < std::abs(rhs); });
            const bool identityOk = update.maxMeasurementDelta < 1.e-9 && update.rhsNorm < 1.e-7 && std::abs(maxStateDelta) < 1.e-9;
            std::cout << Name() << " reference_identity max_delta_m=" << update.maxMeasurementDelta
                      << " rhs_norm=" << update.rhsNorm << " max_delta_state=" << std::abs(maxStateDelta)
                      << " ok=" << identityOk << std::endl;
            assert(identityOk);
          }
          if (oracleStatus == FullFitOracleStatus::Valid)
          {
            auto linearLayer = update.state;
            auto fullLayer = candidateFit.nativeState;
            for (int layer = static_cast<int>(m_trajectorySamplingRadii.size()) - 1; layer >= 0; --layer)
            {
              const double target = m_trajectorySamplingRadii[layer];
              for (unsigned int step = 0; step < 1600 && std::hypot(linearLayer[0], linearLayer[1]) > target + 0.08; ++step)
                linearLayer = TpcTrackKalmanFitter::propagate_state(linearLayer, -0.25, referenceFit.propagationConfig);
              for (unsigned int step = 0; step < 1600 && std::hypot(fullLayer[0], fullLayer[1]) > target + 0.08; ++step)
                fullLayer = TpcTrackKalmanFitter::propagate_state(fullLayer, -0.25, candidateFit.propagationConfig);
              const double linearPhi = std::atan2(linearLayer[1], linearLayer[0]);
              const double fullPhi = std::atan2(fullLayer[1], fullLayer[0]);
              std::cout << Name() << " validation crossing=" << crossingCandidate->crossing << " si_layer=" << layer
                        << " delta_rphi_linear_minus_full=" << target * std::remainder(linearPhi - fullPhi, 2. * M_PI)
                        << " delta_z_linear_minus_full=" << linearLayer[2] - fullLayer[2] << std::endl;
            }
          }
        }
      }
      m_trajectories->add(trajectory); ++deltaBuilds;
      crossingSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - crossingBegin).count();
      if (Verbosity() > 1)
      {
        std::cout << Name() << " crossing=" << crossingCandidate->crossing
                  << " delta_d0=" << update.delta[0] << " delta_z0=" << update.delta[2]
                  << " delta_phi=" << update.delta[3] << " delta_q_over_pt=" << update.delta[4]
                  << " delta_tan_lambda=" << update.delta[5] << " linear_chi2=" << update.chi2 << std::endl;
      }
    }
  }
  if (m_validationFraction > 0.)
  {
    std::cout << Name() << " full_validation_attempted=" << fullValidationAttempted
              << " full_validation_valid=" << fullValidationValid
              << " full_validation_failed=" << fullValidationFailed
              << " full_validation_nonfinite=" << fullValidationNonFinite
              << " full_validation_bad_covariance=" << fullValidationBadCovariance
              << " full_validation_discontinuous=" << fullValidationDiscontinuous
              << " near_validation_attempted=" << nearValidationAttempted
              << " near_validation_valid=" << nearValidationValid << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << Name()
              << " input_tracks=" << inputTracks
              << " rejected_fit_status=" << rejectedFitStatus
              << " rejected_min_pt=" << rejectedMinPt
              << " rejected_min_tpc_clusters=" << rejectedMinTpcClusters
              << " selected_tracks=" << selectedTracks
              << " reference_field_fits=" << referenceFits
              << " crossing_updates=" << deltaBuilds
              << " reference_fit_ms="
              << (referenceFits ? 1.e3 * referenceFitSeconds / referenceFits : 0.)
              << " response_build_ms="
              << (referenceFits ? 1.e3 * responseSeconds / referenceFits : 0.)
              << " crossing_update_us="
              << (deltaBuilds ? 1.e6 * crossingSeconds / deltaBuilds : 0.)
              << " seconds="
              << std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - begin).count()
              << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
