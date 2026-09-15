#include "FastFieldTrackFitter.h"
#include "TpcTrackHelixFitter.h"
#include "TpcTrackKalmanFitter.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyTrack.h"
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace
{
  using Matrix6 = Eigen::Matrix<double, 6, 6>;
  using Matrix36 = Eigen::Matrix<double, 3, 6, Eigen::RowMajor>;
  using Matrix63 = Eigen::Matrix<double, 6, 3, Eigen::RowMajor>;
  using Matrix3 = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
  using Vector6 = Eigen::Matrix<double, 6, 1>;
  using Vector3 = Eigen::Matrix<double, 3, 1>;
  double wrap_phi(const double value) { return std::remainder(value, 2.0 * M_PI); }
}

FastFieldTrackFitter::FastFieldTrackFitter(const PHField* field) : m_field(field) {}

std::array<double, FastFieldTrackFitter::StateSize> FastFieldTrackFitter::externalState(
    const std::array<double, StateSize>& native)
{
  const double qOverPt = native[TpcTrackKalmanFitter::QOverPt];
  const double tanLambda = native[TpcTrackKalmanFitter::TanLambda];
  return {native[0], native[1], native[2], native[3], std::atan2(1.0, tanLambda),
          qOverPt / std::sqrt(1.0 + tanLambda * tanLambda)};
}

std::array<double, FastFieldTrackFitter::StateSize> FastFieldTrackFitter::nativeState(
    const std::array<double, StateSize>& external)
{
  const double theta = external[4];
  return {external[0], external[1], external[2], wrap_phi(external[3]),
          external[5] / std::sin(theta), 1.0 / std::tan(theta)};
}

bool FastFieldTrackFitter::fitMeasurements(const Tpc_PolyTrack& track,
                                            const std::vector<TpcTrackPoint>& input,
                                            Result& output) const
{
  output = Result{};
  if (input.size() < 5 || !m_field) return false;
  auto points = input;
  TpcTrackHelixFitter::order_points(points, TpcTrackPointOrder::Radius);
  TpcKalmanConfig config;
  config.magnetic_field = m_field;
  config.analytic_uniform_propagation = false;
  config.point_order = TpcTrackPointOrder::Input;
  const auto fitBegin = std::chrono::steady_clock::now();
  TpcKalmanResult fit;
  const int charge = track.get_charge() < 0 ? -1 : 1;
  if (!TpcTrackKalmanFitter::fit(points, charge, config, fit) || fit.states_smoothed.empty()) return false;
  output.fitSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - fitBegin).count();
  output.nativeState = fit.states_smoothed.front();
  output.state = externalState(output.nativeState);
  output.covariance = fit.covs_smoothed.front();
  output.pathS = fit.path_s;
  output.propagationConfig = config;
  output.chi2 = fit.chi2;
  output.ndf = fit.ndof;
  output.valid = true;
  return true;
}

bool FastFieldTrackFitter::fit(const Tpc_PolyTrack& track,
                               const std::vector<const Tpc_PolyCluster*>& clusters,
                               Result& output) const
{
  std::vector<std::pair<const Tpc_PolyCluster*, TpcTrackPoint>> ordered;
  ordered.reserve(clusters.size());
  for (const auto* cluster : clusters)
  {
    if (!cluster) continue;
    TpcTrackPoint point;
    point.track_id = static_cast<int>(track.get_track_id());
    point.layer = cluster->size_hits() ? static_cast<int>(TrkrDefs::getLayer(cluster->get_hit_index(0).first)) : 0;
    point.position = {cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z()};
    point.momentum = {track.get_px(), track.get_py(), track.get_pz()};
    ordered.emplace_back(cluster, point);
  }
  std::sort(ordered.begin(), ordered.end(), [](const auto& lhs, const auto& rhs)
  {
    return std::hypot(lhs.second.position.x, lhs.second.position.y) < std::hypot(rhs.second.position.x, rhs.second.position.y);
  });
  std::vector<TpcTrackPoint> points;
  points.reserve(ordered.size());
  for (const auto& item : ordered) points.push_back(item.second);
  if (!fitMeasurements(track, points, output)) return false;

  const auto responseBegin = std::chrono::steady_clock::now();
  output.measurements.reserve(points.size());
  const double varR = output.propagationConfig.meas_sigma_r_cm * output.propagationConfig.meas_sigma_r_cm;
  const double varRPhi = output.propagationConfig.meas_sigma_rphi_cm * output.propagationConfig.meas_sigma_rphi_cm;
  const double varZ = output.propagationConfig.meas_sigma_z_cm * output.propagationConfig.meas_sigma_z_cm;
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    MeasurementResponse response;
    response.key = ordered[i].first->get_trkr_cluster_key();
    response.reference = {points[i].position.x, points[i].position.y, points[i].position.z};
    const double path = i < output.pathS.size() ? output.pathS[i] - output.pathS.front() : 0.0;
    const auto predicted = TpcTrackKalmanFitter::propagate_state(output.nativeState, path, output.propagationConfig);
    response.prediction = {predicted[0], predicted[1], predicted[2]};
    // The response must differentiate the same full-field propagation used for
    // the nominal prediction. The Kalman fitter may use its local-uniform-Bz
    // Jacobian approximation internally, but that is not the derivative of
    // propagate_state in a nonuniform field.
    auto responseConfig = output.propagationConfig;
    responseConfig.rkn_fast_field_jacobian = false;
    const auto fullJacobian = TpcTrackKalmanFitter::propagation_jacobian(output.nativeState, path, responseConfig);
    for (unsigned int row = 0; row < 3; ++row) for (unsigned int col = 0; col < 6; ++col) response.jacobian[row * 6 + col] = fullJacobian[row * 6 + col];
    const double radius = std::hypot(points[i].position.x, points[i].position.y);
    const double c = radius > 0. ? points[i].position.x / radius : 1.;
    const double s = radius > 0. ? points[i].position.y / radius : 0.;
    Matrix3 covariance = Matrix3::Zero();
    covariance(0, 0) = varR * c * c + varRPhi * s * s;
    covariance(1, 1) = varR * s * s + varRPhi * c * c;
    covariance(0, 1) = covariance(1, 0) = (varR - varRPhi) * s * c;
    covariance(2, 2) = varZ;
    const Matrix3 weight = covariance.inverse();
    for (unsigned int row = 0; row < 3; ++row) for (unsigned int col = 0; col < 3; ++col) response.weight[row * 3 + col] = weight(row, col);
    output.measurements.push_back(response);
  }
  Matrix6 normal = Matrix6::Zero();
  const std::array<double, 6> priorSigma{{output.propagationConfig.initial_sigma_pos_cm,
      output.propagationConfig.initial_sigma_pos_cm, output.propagationConfig.initial_sigma_pos_cm,
      output.propagationConfig.initial_sigma_phi, output.propagationConfig.initial_sigma_qop_t,
      output.propagationConfig.initial_sigma_tanl}};
  for (unsigned int i = 0; i < 6; ++i)
  {
    if (!(priorSigma[i] > 0.0) || !std::isfinite(priorSigma[i])) { output.valid = false; return false; }
    normal(i, i) = 1.0 / (priorSigma[i] * priorSigma[i]);
  }
  for (const auto& response : output.measurements)
  {
    const Eigen::Map<const Matrix36> jacobian(response.jacobian.data());
    const Eigen::Map<const Matrix3> weight(response.weight.data());
    normal.noalias() += jacobian.transpose() * weight * jacobian;
  }
  const Eigen::LDLT<Matrix6> decomposition(normal);
  if (decomposition.info() != Eigen::Success || decomposition.vectorD().cwiseAbs().minCoeff() < 1.e-14)
  {
    output.valid = false;
    return false;
  }
  output.informationSolveOk = true;
  const Eigen::SelfAdjointEigenSolver<Matrix6> eigenSolver(normal);
  if (eigenSolver.info() == Eigen::Success)
  {
    for (unsigned int i = 0; i < 6; ++i) output.informationEigenvalues[i] = eigenSolver.eigenvalues()(i);
    const double smallest = eigenSolver.eigenvalues().minCoeff();
    output.informationCondition = smallest > 0.0 ? eigenSolver.eigenvalues().maxCoeff() / smallest : std::numeric_limits<double>::infinity();
  }
  for (auto& response : output.measurements)
  {
    const Eigen::Map<const Matrix36> jacobian(response.jacobian.data());
    const Eigen::Map<const Matrix3> weight(response.weight.data());
    const Matrix63 gain = decomposition.solve(jacobian.transpose() * weight);
    if (!gain.allFinite()) { output.valid = false; return false; }
    for (unsigned int row = 0; row < 6; ++row) for (unsigned int col = 0; col < 3; ++col) response.response[row * 3 + col] = gain(row, col);
  }
  output.responseSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - responseBegin).count();
  return true;
}

FastFieldTrackFitter::Update FastFieldTrackFitter::linearUpdate(
    const Result& reference,
    const std::map<TrkrDefs::cluskey, std::array<double, 3>>& displaced) const
{
  Update output;
  if (!reference.valid || reference.measurements.empty()) return output;
  Vector6 nativeDelta = Vector6::Zero();
  Vector6 rhs = Vector6::Zero();
  for (const auto& response : reference.measurements)
  {
    const auto found = displaced.find(response.key);
    if (found == displaced.end()) return output;
    const Eigen::Map<const Matrix63> gain(response.response.data());
    const Vector3 deltaMeasurement(found->second[0] - response.reference[0],
                                   found->second[1] - response.reference[1],
                                   found->second[2] - response.reference[2]);
    nativeDelta.noalias() += gain * deltaMeasurement;
    const Eigen::Map<const Matrix36> jacobian(response.jacobian.data());
    const Eigen::Map<const Matrix3> weight(response.weight.data());
    rhs.noalias() += jacobian.transpose() * weight * deltaMeasurement;
    output.maxMeasurementDelta = std::max(output.maxMeasurementDelta, deltaMeasurement.cwiseAbs().maxCoeff());
  }
  if (!nativeDelta.allFinite()) return output;
  std::array<double, StateSize> updatedNative = reference.nativeState;
  for (unsigned int i = 0; i < StateSize; ++i) updatedNative[i] += nativeDelta(i);
  updatedNative[TpcTrackKalmanFitter::Phi] = wrap_phi(updatedNative[TpcTrackKalmanFitter::Phi]);
  for (unsigned int i = 0; i < StateSize; ++i) output.delta[i] = nativeDelta(i);
  output.delta[TpcTrackKalmanFitter::Phi] = wrap_phi(updatedNative[TpcTrackKalmanFitter::Phi] - reference.nativeState[TpcTrackKalmanFitter::Phi]);
  output.state = updatedNative;
  output.rhsNorm = rhs.norm();
  for (const auto& response : reference.measurements)
  {
    const auto& measurement = displaced.at(response.key);
    const Eigen::Map<const Matrix36> jacobian(response.jacobian.data());
    const Eigen::Map<const Matrix3> weight(response.weight.data());
    const Vector3 deltaMeasurement(measurement[0] - response.reference[0], measurement[1] - response.reference[1], measurement[2] - response.reference[2]);
    const Vector3 error = deltaMeasurement - jacobian * nativeDelta;
    output.chi2 += (error.transpose() * weight * error)(0, 0);
  }
  output.valid = true;
  return output;
}
