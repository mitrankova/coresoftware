#include "TpcSiliconMatchCandidate.h"

#include <cmath>

ClassImp(TpcSiliconMatchCandidate)

void TpcSiliconMatchCandidate::identify(std::ostream& os) const
{
  os << "TpcSiliconMatchCandidate parent=" << m_parentTrackId << " crossing=" << m_crossing
     << " silicon=" << m_siliconKeys.size() << " midpoint_score=" << m_tpcSiMidpointScore
     << " si_internal_score=" << m_siInternalScore
     << " selected=" << m_selected << std::endl;
}

void TpcSiliconMatchCandidate::Reset()
{
  m_parentTrackId = 0;
  m_sourceAssembledTrackId = 0;
  m_crossing = 0;
  m_score = m_siInternalScore = m_tpcSiMidpointScore = m_maxAbsDz = m_maxAbsDdphi =
      std::numeric_limits<float>::quiet_NaN();
  m_nMvtx = m_nIntt = 0;
  m_selected = false;
  m_siliconKeys.clear();
}

int TpcSiliconMatchCandidate::isValid() const
{
  return std::isfinite(m_tpcSiMidpointScore) && !m_siliconKeys.empty();
}
