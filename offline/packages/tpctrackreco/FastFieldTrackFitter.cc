#include "FastFieldTrackFitter.h"
#include "TpcTrackKalmanFitter.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyTrack.h"
#include <algorithm>
#include <cmath>

FastFieldTrackFitter::FastFieldTrackFitter(const PHField* field) : m_field(field) {}

bool FastFieldTrackFitter::fit(const Tpc_PolyTrack& track,
                               const std::vector<const Tpc_PolyCluster*>& clusters,
                               Result& output) const
{
  std::vector<TpcTrackPoint> points;
  points.reserve(clusters.size());
  for (const auto* cluster : clusters)
  {
    if (!cluster) continue;
    TpcTrackPoint point;
    point.track_id = static_cast<int>(track.get_track_id());
    point.layer = cluster->size_hits() ? static_cast<int>(TrkrDefs::getLayer(cluster->get_hit_index(0).first)) : 0;
    point.position = {cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z()};
    point.momentum = {track.get_px(), track.get_py(), track.get_pz()};
    points.push_back(point);
  }
  if (points.size() < 5 || !m_field) return false;
  TpcKalmanConfig config;
  config.magnetic_field = m_field;
  config.analytic_uniform_propagation = false;
  config.point_order = TpcTrackPointOrder::Radius;
  TpcKalmanResult fit;
  const int charge = track.get_charge() < 0 ? -1 : 1;
  if (!TpcTrackKalmanFitter::fit(points, charge, config, fit) || fit.states_smoothed.empty()) return false;
  const auto& state = fit.states_smoothed.front();
  const auto& covariance = fit.covs_smoothed.front();
  const double qOverPt = state[TpcTrackKalmanFitter::QOverPt];
  const double tanLambda = state[TpcTrackKalmanFitter::TanLambda];
  output.state = {state[TpcTrackKalmanFitter::X], state[TpcTrackKalmanFitter::Y],
                  state[TpcTrackKalmanFitter::Z], state[TpcTrackKalmanFitter::Phi],
                  std::atan2(1.0, tanLambda),
                  qOverPt / std::sqrt(1.0 + tanLambda * tanLambda)};
  output.covariance = covariance;
  output.propagationConfig = config;
  output.chi2 = fit.chi2;
  output.ndf = fit.ndof;
  output.valid = true;
  return true;
}

std::array<double, FastFieldTrackFitter::StateSize> FastFieldTrackFitter::linearUpdate(
    const Result& reference, const std::vector<std::array<double, 3>>& original,
    const std::vector<std::array<double, 3>>& displaced) const
{
  std::array<double, StateSize> delta{};
  const std::size_t count = std::min(original.size(), displaced.size());
  if (!reference.valid || count == 0) return delta;
  double sumR2 = 0.0, sumRdPhi = 0.0, sumR2Dphi = 0.0, sumDz = 0.0, sumDr = 0.0;
  for (std::size_t i = 0; i < count; ++i)
  {
    const double x = original[i][0], y = original[i][1];
    const double dx = displaced[i][0] - x, dy = displaced[i][1] - y;
    const double r2 = x * x + y * y;
    const double dphi = r2 > 1.e-12 ? (x * dy - y * dx) / r2 : 0.0;
    delta[0] += dx; delta[1] += dy; sumDz += displaced[i][2] - original[i][2];
    sumRdPhi += std::sqrt(r2) * dphi; sumR2Dphi += r2 * dphi; sumR2 += r2;
    sumDr += r2 > 1.e-12 ? (x * dx + y * dy) / std::sqrt(r2) : 0.0;
  }
  const double inv = 1.0 / static_cast<double>(count);
  delta[0] *= inv; delta[1] *= inv; delta[2] = sumDz * inv;
  delta[3] = sumRdPhi * inv / std::max(std::sqrt(sumR2 * inv), 1.e-6);
  delta[4] = -sumDr * inv * std::sin(reference.state[4]) / std::max(std::sqrt(sumR2 * inv), 1.e-6);
  delta[5] = reference.state[5] * sumR2Dphi / std::max(sumR2, 1.e-9);
  return delta;
}
