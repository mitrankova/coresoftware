#include "TpcSiliconMatchCandidate.h"

#include <cmath>

ClassImp(TpcSiliconMatchCandidate)

void TpcSiliconMatchCandidate::identify(std::ostream& os) const
{
  os << "TpcSiliconMatchCandidate parent=" << m_parentTrackId << " crossing=" << m_crossing
     << " silicon=" << m_siliconKeys.size() << " direction_score=" << m_tpcSiMatchScore
     << " si_internal_score=" << m_siInternalScore
     << " usable=" << m_tpcSiCompatible
     << " selected=" << m_selected << std::endl;
}

void TpcSiliconMatchCandidate::Reset()
{
  m_parentTrackId = 0;
  m_sourceAssembledTrackId = 0;
  m_crossing = 0;
  m_score = m_siInternalScore = m_tpcSiMatchScore = m_maxAbsDz = m_maxAbsDdphi =
      std::numeric_limits<float>::quiet_NaN();
  m_rSiOuter = m_rTpcInner = m_rMatch = m_midpointDeltaRdphi = m_midpointDeltaZ =
      m_rDirectionMatch = m_outerMvtxDeltaPhi = m_outerMvtxDeltaTanLambda =
          std::numeric_limits<float>::quiet_NaN();
  m_nMvtx = m_nIntt = 0;
  m_tpcSiCompatible = false;
  m_selected = false;
  m_siliconKeys.clear();
}

int TpcSiliconMatchCandidate::isValid() const
{
  return m_tpcSiCompatible && std::isfinite(m_tpcSiMatchScore) &&
         std::isfinite(m_siInternalScore) &&
         std::isfinite(m_rDirectionMatch) &&
         std::isfinite(m_outerMvtxDeltaPhi) &&
         std::isfinite(m_outerMvtxDeltaTanLambda) && !m_siliconKeys.empty();
}
