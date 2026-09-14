#include "Full_PolyTrackv1.h"

#include <cmath>

void Full_PolyTrackv1::identify(std::ostream& os) const
{
  os << "Full_PolyTrackv1 id=" << m_trackId << " parent=" << m_parentTpcTrackId
     << " crossing=" << m_crossing << " TPC=" << m_tpcKeys.size()
     << " MVTX=" << m_nMvtx << " INTT=" << m_nIntt << " score=" << m_score << std::endl;
}

void Full_PolyTrackv1::Reset()
{
  m_event = m_trackId = m_parentTpcTrackId = m_sourceAssembledTrackId = 0;
  m_crossing = 0; m_status = 0; m_nMvtx = m_nIntt = 0;
  m_score = m_maxAbsDz = m_maxAbsDdphi = nan(); m_fitStatus = 0;
  m_chi2 = m_ndf = m_x = m_y = m_z = m_px = m_py = m_pz = m_charge = nan();
  m_cov.fill(nan()); m_tpcKeys.clear(); m_siliconKeys.clear();
}

int Full_PolyTrackv1::isValid() const
{
  return m_fitStatus != 0 && std::isfinite(m_x) && std::isfinite(m_y) &&
         std::isfinite(m_z) && std::isfinite(m_px) && std::isfinite(m_py) && std::isfinite(m_pz);
}

double Full_PolyTrackv1::get_cov(unsigned int i, unsigned int j) const
{ return i < 6 && j < 6 ? m_cov[6 * i + j] : nan(); }

void Full_PolyTrackv1::set_cov(unsigned int i, unsigned int j, double value)
{ if (i < 6 && j < 6) { m_cov[6 * i + j] = value; m_cov[6 * j + i] = value; } }
