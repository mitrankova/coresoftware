#ifndef TPCTRACKRECO_FASTFIELDTRACKFITTER_H
#define TPCTRACKRECO_FASTFIELDTRACKFITTER_H

#include "TpcTrackFit.h"
#include <array>
#include <vector>

class PHField;
class Tpc_PolyCluster;
class Tpc_PolyTrack;

class FastFieldTrackFitter
{
 public:
  static constexpr unsigned int StateSize = 6;
  struct Result
  {
    bool valid{false};
    std::array<double, StateSize> state{}; // x,y,z,phi,theta,q/p
    std::array<double, StateSize * StateSize> covariance{};
    TpcKalmanConfig propagationConfig;
    double chi2{0.0};
    int ndf{-1};
  };

  explicit FastFieldTrackFitter(const PHField* field);
  bool fit(const Tpc_PolyTrack&, const std::vector<const Tpc_PolyCluster*>&, Result&) const;
  std::array<double, StateSize> linearUpdate(
      const Result&, const std::vector<std::array<double, 3>>& reference,
      const std::vector<std::array<double, 3>>& displaced) const;

 private:
  const PHField* m_field{nullptr};
};
#endif
