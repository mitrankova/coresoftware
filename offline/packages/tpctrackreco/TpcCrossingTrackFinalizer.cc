#include "TpcCrossingTrackFinalizer.h"
#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrackv1.h"
#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
#include "TpcCrossingClusterPosition.h"
#include "TpcSiliconMatchCandidate.h"
#include "TpcSiliconMatchCandidateContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyClusterv1.h"
#include "Tpc_PolyClusterContainer.h"
#include "Tpc_PolyClusterContainerv1.h"
#include "TpcDriftPolylineLookup.h"
#include "TpcTrackHelixFitter.h"
#include "TpcTrackKalmanFitter.h"
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
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/ActsGeometry.h>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Definitions/Units.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace
{
  using NativeState = std::array<double, FastFieldTrackFitter::StateSize>;
  using NativeCovariance = std::array<double, FastFieldTrackFitter::StateSize *
                                               FastFieldTrackFitter::StateSize>;

  const char* finalFitStatusName(const FullPolyTrackFitStatus status)
  {
    switch (status)
    {
    case FullPolyTrackFitStatus::FinalFitAccepted: return "FinalFitAccepted";
    case FullPolyTrackFitStatus::FinalFitFailedFallback: return "FinalFitFailedFallback";
    case FullPolyTrackFitStatus::FinalFitNonFiniteFallback: return "FinalFitNonFiniteFallback";
    case FullPolyTrackFitStatus::FinalFitBadCovarianceFallback: return "FinalFitBadCovarianceFallback";
    case FullPolyTrackFitStatus::FinalFitDiscontinuousFallback: return "FinalFitDiscontinuousFallback";
    case FullPolyTrackFitStatus::Unknown: return "Unknown";
    }
    return "Unknown";
  }

  FullPolyTrackFitStatus validateFinalFit(
      const bool fitCalled,
      const NativeState& fastNative,
      const NativeCovariance& fastCovariance,
      const FastFieldTrackFitter::Result& finalFit,
      const double maximumContinuityPull,
      const std::array<bool, 3>& continuityEnabled,
      double& continuityMaxPull)
  {
    continuityMaxPull = std::numeric_limits<double>::quiet_NaN();
    if (!fitCalled || !finalFit.fitSuccess || !finalFit.valid || finalFit.nAccepted == 0)
      return FullPolyTrackFitStatus::FinalFitFailedFallback;
    if (!std::isfinite(finalFit.chi2))
      return FullPolyTrackFitStatus::FinalFitNonFiniteFallback;
    for (const double value : finalFit.nativeState)
      if (!std::isfinite(value)) return FullPolyTrackFitStatus::FinalFitNonFiniteFallback;
    for (const double value : finalFit.covariance)
      if (!std::isfinite(value)) return FullPolyTrackFitStatus::FinalFitNonFiniteFallback;

    for (unsigned int index = 0; index < FastFieldTrackFitter::StateSize; ++index)
      if (!(finalFit.covariance[7 * index] > 0.))
        return FullPolyTrackFitStatus::FinalFitBadCovarianceFallback;

    continuityMaxPull = 0.;
    const std::array<unsigned int, 3> continuityIndices{{
        TpcTrackKalmanFitter::Phi,
        TpcTrackKalmanFitter::QOverPt,
        TpcTrackKalmanFitter::TanLambda}};
    constexpr double twoPi = 6.28318530717958647692;
    for (std::size_t continuityIndex = 0; continuityIndex < continuityIndices.size(); ++continuityIndex)
    {
      if (!continuityEnabled[continuityIndex]) continue;
      const unsigned int index = continuityIndices[continuityIndex];
      double delta = finalFit.nativeState[index] - fastNative[index];
      if (index == TpcTrackKalmanFitter::Phi) delta = std::remainder(delta, twoPi);
      // The trajectory covariance describes the reference fast fit. It is used
      // as the QA reference covariance without attempting to transport it.
      const double variance = fastCovariance[7 * index] + finalFit.covariance[7 * index];
      if (!(variance > 0.) || !std::isfinite(variance))
        return FullPolyTrackFitStatus::FinalFitBadCovarianceFallback;
      continuityMaxPull = std::max(continuityMaxPull, std::abs(delta) / std::sqrt(variance));
    }
    return continuityMaxPull > maximumContinuityPull
        ? FullPolyTrackFitStatus::FinalFitDiscontinuousFallback
        : FullPolyTrackFitStatus::FinalFitAccepted;
  }

  bool finiteState(const NativeState& state)
  {
    return std::all_of(state.begin(), state.end(), [](const double value)
    { return std::isfinite(value); });
  }

  double positionDistance(const NativeState& state, const TpcTrackPoint& target)
  {
    return std::hypot(std::hypot(state[TpcTrackKalmanFitter::X] - target.position.x,
                                 state[TpcTrackKalmanFitter::Y] - target.position.y),
                      state[TpcTrackKalmanFitter::Z] - target.position.z);
  }

  bool transportToMeasurement(
      const NativeState& state, const TpcTrackPoint& target,
      const TpcKalmanConfig& config, const Acts::Surface* surface,
      const Acts::GeometryContext& geoContext, NativeState& output,
      double& pathCm, double& distanceCm)
  {
    pathCm = 0.;
    distanceCm = std::numeric_limits<double>::quiet_NaN();
    if (target.detector != TpcTrackPoint::Detector::Tpc && surface &&
        TpcTrackKalmanFitter::propagate_to_surface(
            state, config, *surface, geoContext, output, &pathCm) && finiteState(output))
    {
      distanceCm = positionDistance(output, target);
      return distanceCm <= 1.;
    }

    pathCm = 0.;

    constexpr double stepCm = -0.25;
    constexpr unsigned int maximumSteps = 1600;
    const double targetRadius = std::hypot(target.position.x, target.position.y);
    NativeState previous = state;
    double previousPath = 0.;
    NativeState current = state;
    bool bracketed = std::hypot(current[TpcTrackKalmanFitter::X],
                               current[TpcTrackKalmanFitter::Y]) <= targetRadius;
    for (unsigned int step = 0; step < maximumSteps && !bracketed; ++step)
    {
      previous = current;
      previousPath = pathCm;
      current = TpcTrackKalmanFitter::propagate_state(previous, stepCm, config);
      pathCm += stepCm;
      if (!finiteState(current)) return false;
      bracketed = std::hypot(current[TpcTrackKalmanFitter::X],
                             current[TpcTrackKalmanFitter::Y]) <= targetRadius;
    }
    if (!bracketed) return false;

    double lowPath = pathCm;
    double highPath = previousPath;
    for (unsigned int iteration = 0; iteration < 16; ++iteration)
    {
      const double width = highPath - lowPath;
      const double leftPath = lowPath + 0.25 * width;
      const double rightPath = highPath - 0.25 * width;
      const auto left = TpcTrackKalmanFitter::propagate_state(state, leftPath, config);
      const auto right = TpcTrackKalmanFitter::propagate_state(state, rightPath, config);
      if (!finiteState(left) || !finiteState(right)) return false;
      if (positionDistance(left, target) < positionDistance(right, target))
        highPath = 0.5 * (lowPath + highPath);
      else
        lowPath = 0.5 * (lowPath + highPath);
    }
    pathCm = 0.5 * (lowPath + highPath);
    output = TpcTrackKalmanFitter::propagate_state(state, pathCm, config);
    if (!finiteState(output)) return false;
    distanceCm = positionDistance(output, target);
    return distanceCm <= 1.;
  }

  double median(std::vector<double> values)
  {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const auto middle = values.begin() + values.size() / 2;
    std::nth_element(values.begin(), middle, values.end());
    if (values.size() % 2 != 0) return *middle;
    return 0.5 * (*middle + *std::max_element(values.begin(), middle));
  }
}

TpcCrossingTrackFinalizer::TpcCrossingTrackFinalizer(const std::string& name) : SubsysReco(name) {}
int TpcCrossingTrackFinalizer::getNodes(PHCompositeNode* topNode)
{
  m_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_trackNodeName);
  m_trajectories = findNode::getClass<TpcCrossingTrajectoryContainer>(topNode, m_trajectoryNodeName);
  m_candidates = findNode::getClass<TpcSiliconMatchCandidateContainer>(topNode, m_candidateNodeName);
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_clusterNodeName);
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_trkrClusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  m_lookup = TpcDriftPolylineLookup::get(topNode);
  m_field = PHFieldUtility::GetFieldMapNode(nullptr, topNode, Verbosity());
  if (!m_tracks || !m_trajectories || !m_candidates || !m_clusters || !m_hits || !m_trkrClusters || !m_geometry || !m_lookup || !m_field)
  {
    std::cerr << Name() << "::getNodes - missing full-track input node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrackFinalizer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dst = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst) { dst = new PHCompositeNode("DST"); topNode->addNode(dst); }
  m_output = findNode::getClass<Full_PolyTrackContainer>(topNode, m_outputNodeName);
  if (!m_output)
  {
    m_output = new Full_PolyTrackContainerv1;
    dst->addNode(new PHIODataNode<PHObject>(m_output, m_outputNodeName, "PHObject"));
  }
  m_correctedClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_correctedClusterNodeName);
  if (!m_correctedClusters)
  {
    m_correctedClusters = new Tpc_PolyClusterContainerv1;
    dst->addNode(new PHIODataNode<PHObject>(m_correctedClusters, m_correctedClusterNodeName, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrackFinalizer::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  m_event = 0;
  m_fitter = std::make_unique<FastFieldTrackFitter>(m_field);
  return createNodes(topNode);
}
int TpcCrossingTrackFinalizer::process_event(PHCompositeNode*)
{
  const auto begin = std::chrono::steady_clock::now();
  ++m_event;
  m_output->Reset();
  m_correctedClusters->Reset();
  std::map<TrkrDefs::cluskey, const Tpc_PolyCluster*> clustersByKey;
  for (unsigned int i = 0; i < m_clusters->size(); ++i) if (const auto* cluster = m_clusters->get_cluster(i)) clustersByKey[cluster->get_trkr_cluster_key()] = cluster;
  unsigned int outputId = 0;
  unsigned int selectedEqualReference = 0;
  unsigned int selectedNonreference = 0;
  unsigned int selectedOffsetAbs1 = 0;
  unsigned int selectedOffsetAbs2 = 0;
  unsigned int selectedOffsetAbs3 = 0;
  unsigned int selectedOffsetAbsGt3 = 0;
  unsigned int combinedFitAttempted = 0;
  unsigned int combinedFitCalls = 0;
  unsigned int combinedFitAccepted = 0;
  unsigned int combinedFitFailed = 0;
  unsigned int combinedFitNonfinite = 0;
  unsigned int combinedFitBadCovariance = 0;
  unsigned int combinedFitDiscontinuous = 0;
  unsigned int combinedFitFastFallback = 0;
  unsigned int combinedFitTransportFailed = 0;
  unsigned int finalFitQaPrinted = 0;
  std::vector<double> acceptedChi2Ndf;
  std::map<int, unsigned int> crossingOffsets;
  double finalFitSeconds = 0.0;
  for (unsigned int i = 0; i < m_candidates->size(); ++i)
  {
    const auto* candidate = m_candidates->get(i);
    if (!candidate || !candidate->get_selected()) continue;
    const Tpc_PolyTrack* parent = nullptr;
    for (unsigned int j = 0; j < m_tracks->size(); ++j)
    {
      const auto* value = m_tracks->get_track(j);
      if (value && value->get_track_id() == candidate->get_parent_track_id()) { parent = value; break; }
    }
    if (!parent) continue;
    const TpcCrossingTrajectory* trajectory = nullptr;
    for (unsigned int j = 0; j < m_trajectories->size(); ++j)
    {
      const auto* value = m_trajectories->get(j);
      if (value && value->get_parent_track_id() == candidate->get_parent_track_id() && value->get_crossing() == candidate->get_crossing()) { trajectory = value; break; }
    }
    if (!trajectory) continue;

    const short referenceCrossing = trajectory->get_reference_crossing();
    const short selectedCrossing = trajectory->get_crossing();
    const int crossingOffset = static_cast<int>(selectedCrossing) -
                               static_cast<int>(referenceCrossing);
    ++crossingOffsets[crossingOffset];
    if (crossingOffset == 0) ++selectedEqualReference;
    else ++selectedNonreference;
    switch (std::abs(crossingOffset))
    {
    case 0: break;
    case 1: ++selectedOffsetAbs1; break;
    case 2: ++selectedOffsetAbs2; break;
    case 3: ++selectedOffsetAbs3; break;
    default: ++selectedOffsetAbsGt3; break;
    }

    std::vector<const Tpc_PolyCluster*> fitClusters;
    for (const auto key : parent->get_cluster_keys())
    {
      const auto found = clustersByKey.find(key);
      if (found == clustersByKey.end()) continue;
      const auto* reference = found->second;
      auto* corrected = new Tpc_PolyClusterv1;
      corrected->set_event(m_event); corrected->set_cluster_id(m_correctedClusters->size());
      corrected->set_source_assembled_track_id(reference->get_source_assembled_track_id());
      corrected->set_trkr_cluster_key(reference->get_trkr_cluster_key()); corrected->set_side(reference->get_side());
      corrected->set_adc(reference->get_adc()); corrected->set_phi_width(reference->get_phi_width()); corrected->set_time_width(reference->get_time_width()); corrected->set_phase(reference->get_phase());
      std::array<double, 3> position{};
      std::vector<TpcCrossingClusterPosition::HitPosition> hitPositions;
      if (!TpcCrossingClusterPosition::get(*reference, *m_lookup, candidate->get_crossing(),
                                           trajectory->get_reference_crossing(), position, &hitPositions) ||
          hitPositions.size() != reference->size_hits())
      {
        delete corrected;
        continue;
      }
      double sx2 = 0.0, sy2 = 0.0, sz2 = 0.0;
      for (const auto& hitPosition : hitPositions)
      {
        corrected->add_hit(hitPosition.hitsetkey, hitPosition.hitkey,
                           hitPosition.candidate[0], hitPosition.candidate[1], hitPosition.candidate[2]);
        sx2 += (hitPosition.candidate[0] - position[0]) * (hitPosition.candidate[0] - position[0]);
        sy2 += (hitPosition.candidate[1] - position[1]) * (hitPosition.candidate[1] - position[1]);
        sz2 += (hitPosition.candidate[2] - position[2]) * (hitPosition.candidate[2] - position[2]);
      }
      const double inverseCount = 1.0 / static_cast<double>(hitPositions.size());
      corrected->set_centroid_x(position[0]); corrected->set_centroid_y(position[1]); corrected->set_centroid_z(position[2]);
      corrected->set_rms_x(std::sqrt(sx2 * inverseCount)); corrected->set_rms_y(std::sqrt(sy2 * inverseCount)); corrected->set_rms_z(std::sqrt(sz2 * inverseCount));
      m_correctedClusters->add_cluster(corrected); fitClusters.push_back(corrected);
    }
    std::vector<TpcTrackPoint> measurements;
    measurements.reserve(fitClusters.size() + candidate->get_silicon_cluster_keys().size());
    for (const auto* cluster : fitClusters)
    {
      TpcTrackPoint point;
      point.track_id = static_cast<int>(parent->get_track_id());
      point.layer = cluster->size_hits() ? static_cast<int>(TrkrDefs::getLayer(cluster->get_hit_index(0).first)) : 0;
      point.position = {cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z()};
      point.momentum = {parent->get_px(), parent->get_py(), parent->get_pz()};
      point.detector = TpcTrackPoint::Detector::Tpc;
      point.cluster_key = cluster->get_trkr_cluster_key();
      measurements.push_back(point);
    }
    for (const auto key : candidate->get_silicon_cluster_keys())
    {
      auto* cluster = m_trkrClusters->findCluster(key);
      const auto surface = cluster ? m_geometry->maps().getSurface(key, cluster) : nullptr;
      if (!cluster || !surface) continue;
      const auto global = m_geometry->getGlobalPosition(key, cluster);
      const Acts::Vector2 local(cluster->getLocalX() * Acts::UnitConstants::cm,
                                cluster->getLocalY() * Acts::UnitConstants::cm);
      const auto& context = m_geometry->geometry().geoContext;
      const Acts::Vector3 direction(1., 1., 1.);
      const auto origin = surface->localToGlobal(context, local, direction);
      const auto along0 = surface->localToGlobal(context, local + Acts::Vector2(Acts::UnitConstants::cm, 0.), direction) - origin;
      const auto along1 = surface->localToGlobal(context, local + Acts::Vector2(0., Acts::UnitConstants::cm), direction) - origin;
      const auto axis0 = along0.normalized();
      const auto axis1 = along1.normalized();
      const bool isMvtx = TrkrDefs::getTrkrId(key) == TrkrDefs::mvtxId;
      const auto& misalignmentSigma = isMvtx
          ? m_mvtxMisalignmentSigma : m_inttMisalignmentSigma;
      const double clusterSigma0 = std::max(1.e-6, static_cast<double>(cluster->getRPhiError()));
      const double clusterSigma1 = std::max(1.e-6, static_cast<double>(cluster->getZError()));
      const double variance0 = clusterSigma0 * clusterSigma0 +
                               misalignmentSigma[0] * misalignmentSigma[0];
      const double variance1 = clusterSigma1 * clusterSigma1 +
                               misalignmentSigma[1] * misalignmentSigma[1];
      TpcTrackPoint point;
      point.track_id = static_cast<int>(parent->get_track_id());
      point.layer = static_cast<int>(TrkrDefs::getLayer(key));
      point.position = {global.x(), global.y(), global.z()};
      point.momentum = {parent->get_px(), parent->get_py(), parent->get_pz()};
      point.detector = isMvtx ? TpcTrackPoint::Detector::Mvtx : TpcTrackPoint::Detector::Intt;
      point.cluster_key = key;
      point.measurement_dimension = 2;
      point.has_measurement_model = true;
      point.measurement_projection = {axis0.x(), axis0.y(), axis0.z(), axis1.x(), axis1.y(), axis1.z(), 0., 0., 0.};
      point.measurement_covariance = {variance0, 0., 0., 0., variance1, 0., 0., 0., 1.};
      measurements.push_back(point);
    }
    const NativeState fastNative{{
        trajectory->get_state(TpcCrossingTrajectory::X),
        trajectory->get_state(TpcCrossingTrajectory::Y),
        trajectory->get_state(TpcCrossingTrajectory::Z),
        trajectory->get_state(TpcCrossingTrajectory::Phi),
        trajectory->get_state(TpcCrossingTrajectory::QOverPt),
        trajectory->get_state(TpcCrossingTrajectory::TanLambda)}};
    NativeCovariance fastCovariance{};
    for (unsigned int row = 0; row < FastFieldTrackFitter::StateSize; ++row)
      for (unsigned int column = 0; column < FastFieldTrackFitter::StateSize; ++column)
        fastCovariance[FastFieldTrackFitter::StateSize * row + column] =
            trajectory->get_covariance(row, column);

    NativeState transportedFastNative = fastNative;
    double transportPathCm = std::numeric_limits<double>::quiet_NaN();
    double transportDistanceCm = std::numeric_limits<double>::quiet_NaN();
    bool transportSucceeded = false;
    if (!measurements.empty())
    {
      auto orderedMeasurements = measurements;
      // Keep this ordering identical to FastFieldTrackFitter::fitMeasurementsImpl.
      TpcTrackHelixFitter::order_points(orderedMeasurements, TpcTrackPointOrder::Radius);
      const auto& target = orderedMeasurements.front();
      TrkrCluster* targetCluster = target.detector == TpcTrackPoint::Detector::Tpc
          ? nullptr : m_trkrClusters->findCluster(target.cluster_key);
      const auto targetSurface = targetCluster
          ? m_geometry->maps().getSurface(target.cluster_key, targetCluster) : nullptr;
      transportSucceeded = transportToMeasurement(
          fastNative, target, m_fitter->makePropagationConfig(), targetSurface.get(),
          m_geometry->geometry().geoContext, transportedFastNative,
          transportPathCm, transportDistanceCm);
    }

    FastFieldTrackFitter::Result finalFit;
    const bool allMeasurements = fitClusters.size() == parent->size_cluster_keys() &&
                                 measurements.size() == fitClusters.size() + candidate->get_silicon_cluster_keys().size();
    ++combinedFitAttempted;
    bool fitReturned = false;
    if (allMeasurements && transportSucceeded)
    {
      ++combinedFitCalls;
      fitReturned = m_fitter->fitMeasurements(
          *parent, measurements, transportedFastNative, finalFit);
      finalFitSeconds += finalFit.fitSeconds;
    }
    else if (allMeasurements)
    {
      ++combinedFitTransportFailed;
    }
    double continuityMaxPull = std::numeric_limits<double>::quiet_NaN();
    const auto finalFitStatus = validateFinalFit(
        fitReturned, transportedFastNative, fastCovariance, finalFit,
        m_finalFitContinuityMaxPull, m_finalFitContinuityIndices,
        continuityMaxPull);
    const bool finalFitAccepted =
        finalFitStatus == FullPolyTrackFitStatus::FinalFitAccepted;
    switch (finalFitStatus)
    {
    case FullPolyTrackFitStatus::FinalFitAccepted: ++combinedFitAccepted; break;
    case FullPolyTrackFitStatus::FinalFitFailedFallback: ++combinedFitFailed; break;
    case FullPolyTrackFitStatus::FinalFitNonFiniteFallback: ++combinedFitNonfinite; break;
    case FullPolyTrackFitStatus::FinalFitBadCovarianceFallback: ++combinedFitBadCovariance; break;
    case FullPolyTrackFitStatus::FinalFitDiscontinuousFallback: ++combinedFitDiscontinuous; break;
    case FullPolyTrackFitStatus::Unknown: ++combinedFitFailed; break;
    }
    if (!finalFitAccepted) ++combinedFitFastFallback;
    if (finalFitAccepted && finalFit.ndf > 0)
      acceptedChi2Ndf.push_back(finalFit.chi2 / static_cast<double>(finalFit.ndf));

    std::array<double, FastFieldTrackFitter::StateSize> finalMinusFast{};
    finalMinusFast.fill(std::numeric_limits<double>::quiet_NaN());
    if (fitReturned)
    {
      for (unsigned int index = 0; index < finalMinusFast.size(); ++index)
      {
        if (!std::isfinite(finalFit.nativeState[index])) continue;
        finalMinusFast[index] = finalFit.nativeState[index] - transportedFastNative[index];
      }
      if (std::isfinite(finalMinusFast[TpcTrackKalmanFitter::Phi]))
        finalMinusFast[TpcTrackKalmanFitter::Phi] =
            std::remainder(finalMinusFast[TpcTrackKalmanFitter::Phi],
                           6.28318530717958647692);
    }

    const double fastFinalPositionDistanceCm = fitReturned && finiteState(finalFit.nativeState)
        ? std::hypot(std::hypot(
              finalFit.nativeState[TpcTrackKalmanFitter::X] - transportedFastNative[TpcTrackKalmanFitter::X],
              finalFit.nativeState[TpcTrackKalmanFitter::Y] - transportedFastNative[TpcTrackKalmanFitter::Y]),
              finalFit.nativeState[TpcTrackKalmanFitter::Z] - transportedFastNative[TpcTrackKalmanFitter::Z])
        : std::numeric_limits<double>::quiet_NaN();
    unsigned int nTpcFit = 0, nMvtxFit = 0, nInttFit = 0;
    double chi2Tpc = 0., chi2Mvtx = 0., chi2Intt = 0.;
    const auto numberFitMeasurements = std::min(
        finalFit.measurementChi2.size(), finalFit.measurementDetector.size());
    for (std::size_t index = 0; index < numberFitMeasurements; ++index)
    {
      const auto detector = static_cast<TpcTrackPoint::Detector>(finalFit.measurementDetector[index]);
      if (detector == TpcTrackPoint::Detector::Tpc)
      { ++nTpcFit; chi2Tpc += finalFit.measurementChi2[index]; }
      else if (detector == TpcTrackPoint::Detector::Mvtx)
      { ++nMvtxFit; chi2Mvtx += finalFit.measurementChi2[index]; }
      else if (detector == TpcTrackPoint::Detector::Intt)
      { ++nInttFit; chi2Intt += finalFit.measurementChi2[index]; }
    }

    const bool printFinalFitQa = Verbosity() >= 10 ||
        (Verbosity() >= 5 &&
         (finalFitQaPrinted < m_maxFinalFitQaTracks || !finalFitAccepted));
    if (printFinalFitQa)
    {
      ++finalFitQaPrinted;
      std::cout << Name() << " crossing_selection"
                << " parent=" << parent->get_track_id()
                << " source_assembled=" << parent->get_source_assembled_track_id()
                << " reference=" << referenceCrossing
                << " selected=" << selectedCrossing
                << " offset=" << crossingOffset << std::endl;
      std::cout << Name() << " final_fit_qa"
                << " parent_track_id=" << parent->get_track_id()
                << " reference_crossing=" << referenceCrossing
                << " selected_crossing=" << selectedCrossing
                << " crossing_offset=" << crossingOffset
                << " final_fit_status=" << finalFitStatusName(finalFitStatus)
                << " final_fit_success=" << finalFit.fitSuccess
                << " chi2=" << finalFit.chi2
                << " ndf=" << finalFit.ndf
                << " nAccepted=" << finalFit.nAccepted
                << " n_tpc=" << nTpcFit
                << " n_mvtx_fit=" << nMvtxFit
                << " n_intt_fit=" << nInttFit
                << " chi2_tpc=" << chi2Tpc
                << " chi2_mvtx=" << chi2Mvtx
                << " chi2_intt=" << chi2Intt
                << " transport_path_cm=" << transportPathCm
                << " transport_distance_cm=" << transportDistanceCm
                << " fast_final_position_distance_cm=" << fastFinalPositionDistanceCm
                << " deltaX_final_minus_fast=" << finalMinusFast[TpcTrackKalmanFitter::X]
                << " deltaY_final_minus_fast=" << finalMinusFast[TpcTrackKalmanFitter::Y]
                << " deltaZ_final_minus_fast=" << finalMinusFast[TpcTrackKalmanFitter::Z]
                << " deltaPhi_final_minus_fast=" << finalMinusFast[TpcTrackKalmanFitter::Phi]
                << " deltaQOverPt_final_minus_fast=" << finalMinusFast[TpcTrackKalmanFitter::QOverPt]
                << " deltaTanLambda_final_minus_fast=" << finalMinusFast[TpcTrackKalmanFitter::TanLambda]
                << " continuity_max_pull=" << continuityMaxPull << std::endl;
    }

    auto* full = new Full_PolyTrackv1;
    full->set_event(m_event);
    full->set_track_id(outputId++);
    full->set_parent_tpc_track_id(parent->get_track_id());
    full->set_source_assembled_track_id(parent->get_source_assembled_track_id());
    full->set_crossing(candidate->get_crossing());
    full->set_status(1);
    full->set_n_mvtx(candidate->get_n_mvtx());
    full->set_n_intt(candidate->get_n_intt());
    full->set_score(candidate->get_score());
    full->set_max_abs_dz(candidate->get_max_abs_dz());
    full->set_max_abs_ddphi(candidate->get_max_abs_ddphi());
    full->set_fit_status(static_cast<int>(finalFitStatus));
    full->set_chi2(finalFitAccepted ? finalFit.chi2 : parent->get_chi2());
    full->set_ndf(finalFitAccepted ? finalFit.ndf : parent->get_ndf());
    const auto& fallbackNative = transportSucceeded ? transportedFastNative : fastNative;
    const auto& finalNative = finalFitAccepted ? finalFit.nativeState : fallbackNative;
    for (unsigned int index = 0; index < finalNative.size(); ++index)
    {
      full->set_final_native_state(index, finalNative[index]);
      // The fast state is defined at the innermost combined-fit measurement.
      full->set_fast_native_state(index, fallbackNative[index]);
    }
    const auto fastExternal = FastFieldTrackFitter::externalState(fallbackNative);
    const auto& state = finalFitAccepted ? finalFit.state : fastExternal;
    full->set_x(state[0]); full->set_y(state[1]); full->set_z(state[2]);
    const double momentum = std::abs(state[5]) > 1.e-12
        ? std::abs(1. / state[5])
        : std::hypot(std::hypot(parent->get_px(), parent->get_py()), parent->get_pz());
    const double pt = momentum * std::sin(state[4]);
    full->set_px(pt * std::cos(state[3]));
    full->set_py(pt * std::sin(state[3]));
    full->set_pz(momentum * std::cos(state[4]));
    full->set_charge(state[5] < 0. ? -1. : 1.);
    for (unsigned int row = 0; row < FastFieldTrackFitter::StateSize; ++row)
      for (unsigned int column = 0; column < FastFieldTrackFitter::StateSize; ++column)
        full->set_cov(row, column, finalFitAccepted
            ? finalFit.covariance[FastFieldTrackFitter::StateSize * row + column]
            : fastCovariance[FastFieldTrackFitter::StateSize * row + column]);
    for (const auto key : parent->get_cluster_keys()) full->add_tpc_cluster_key(key);
    for (const auto key : candidate->get_silicon_cluster_keys()) full->add_silicon_cluster_key(key);
    m_output->add_track(full);
  }
  if (Verbosity() > 0)
  {
    std::cout << Name()
              << " final_tracks=" << m_output->size()
              << " selected_equal_reference=" << selectedEqualReference
              << " selected_nonreference=" << selectedNonreference
              << " selected_offset_abs_1=" << selectedOffsetAbs1
              << " selected_offset_abs_2=" << selectedOffsetAbs2
              << " selected_offset_abs_3=" << selectedOffsetAbs3
              << " selected_offset_abs_gt3=" << selectedOffsetAbsGt3
              << " combined_fit_attempted=" << combinedFitAttempted
              << " combined_fit_accepted=" << combinedFitAccepted
              << " combined_fit_failed=" << combinedFitFailed
              << " combined_fit_nonfinite=" << combinedFitNonfinite
              << " combined_fit_bad_covariance=" << combinedFitBadCovariance
              << " combined_fit_discontinuous=" << combinedFitDiscontinuous
              << " combined_fit_fast_fallback=" << combinedFitFastFallback
              << " combined_fit_transport_failed=" << combinedFitTransportFailed
              << " accepted_median_chi2_ndf=" << median(acceptedChi2Ndf)
              << " corrected_clusters=" << m_correctedClusters->size()
              << " final_fit_ms=" << (combinedFitCalls
                  ? 1.e3 * finalFitSeconds / combinedFitCalls : 0.)
              << " seconds="
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - begin).count()
              << std::endl;
    std::cout << Name() << " crossing_offsets";
    for (const auto& [offset, count] : crossingOffsets)
      std::cout << " " << offset << ":" << count;
    std::cout << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
