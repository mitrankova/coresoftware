#ifndef TPCTRACKRECO_FASTFIELDTRACKFITTER_H
#define TPCTRACKRECO_FASTFIELDTRACKFITTER_H
#include "TpcTrackFit.h"
#include <trackbase/TrkrDefs.h>
#include <array>
#include <map>
#include <string>
#include <vector>
class PHField;
class Tpc_PolyCluster;
class Tpc_PolyTrack;
class FastFieldTrackFitter
{
 public:
  static constexpr unsigned int StateSize = 6;
  struct MeasurementResponse
  {
    TrkrDefs::cluskey key{TrkrDefs::CLUSKEYMAX};
    std::array<double, 3> reference{};
    std::array<double, 3> prediction{};
    std::array<double, 18> jacobian{}; // d(x,y,z)_i / dX_reference, row-major 3x6
    std::array<double, 9> weight{};    // measurement inverse covariance, 3x3
    std::array<double, 18> response{}; // A^-1 J_i^T W_i, row-major 6x3
  };
  struct Result
  {
    bool valid{false};
    std::array<double, StateSize> state{};       // x,y,z,phi,theta,q/p
    std::array<double, StateSize> nativeState{}; // x,y,z,phi,q/pt,tan(lambda)
    std::array<double, StateSize * StateSize> covariance{};
    std::vector<MeasurementResponse> measurements;
    std::vector<double> pathS;
    TpcKalmanConfig propagationConfig;
    double chi2{0.0};
    int ndf{-1};
    double fitSeconds{0.0};
    double responseSeconds{0.0};
    std::array<double, StateSize> informationEigenvalues{};
    double informationCondition{0.0};
    bool informationSolveOk{false};
    bool fitSuccess{false};
    std::string fitMessage;
    std::size_t nMeasurements{0};
    std::size_t nAccepted{0};
  };
  struct Update
  {
    bool valid{false};
    std::array<double, StateSize> delta{};  // native: x,y,z,phi,q/pt,tan(lambda)
    std::array<double, StateSize> state{};
    double rhsNorm{0.0};
    double maxMeasurementDelta{0.0};
    double chi2{0.0};
  };
  explicit FastFieldTrackFitter(const PHField* field);
  bool fit(const Tpc_PolyTrack&, const std::vector<const Tpc_PolyCluster*>&, Result&) const;
  bool fitMeasurements(const Tpc_PolyTrack&, const std::vector<TpcTrackPoint>&, Result&) const;
  bool fitMeasurements(const Tpc_PolyTrack&, const std::vector<TpcTrackPoint>&,
                       const std::array<double, StateSize>& initialNativeState, Result&) const;
  Update linearUpdate(const Result&, const std::map<TrkrDefs::cluskey, std::array<double, 3>>&) const;
  static std::array<double, StateSize> externalState(const std::array<double, StateSize>& native);
  static std::array<double, StateSize> nativeState(const std::array<double, StateSize>& external);
 private:
  bool fitMeasurementsImpl(const Tpc_PolyTrack&, const std::vector<TpcTrackPoint>&,
                           const std::array<double, StateSize>* initialNativeState, Result&) const;
  const PHField* m_field{nullptr};
};
#endif
