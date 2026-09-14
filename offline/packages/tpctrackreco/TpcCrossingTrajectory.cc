#include "TpcCrossingTrajectory.h"

#include <algorithm>

ClassImp(TpcCrossingTrajectory)

void TpcCrossingTrajectory::identify(std::ostream& os) const
{
  os << "TpcCrossingTrajectory parent=" << m_parentTrackId << " crossing=" << m_crossing
     << " linear chi2=" << m_linearChi2 << std::endl;
}

void TpcCrossingTrajectory::Reset()
{
  m_parentTrackId = 0;
  m_crossing = 0;
  m_referenceCrossing = 0;
  m_state.fill(0.F);
  m_delta.fill(0.F);
  m_linearChi2 = nan();
  m_layerStates.fill(LayerState{});
  m_nLayerStates = 0;
}

int TpcCrossingTrajectory::isValid() const
{
  return std::isfinite(m_state[X]) && std::isfinite(m_state[Y]) &&
         std::isfinite(m_state[Z]) && std::isfinite(m_state[Phi]) &&
         std::isfinite(m_state[Theta]) && std::isfinite(m_state[QOverP]);
}
