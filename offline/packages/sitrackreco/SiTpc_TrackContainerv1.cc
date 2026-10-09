#include "SiTpc_TrackContainerv1.h"

#include "SiTpc_Track.h"

SiTpc_TrackContainerv1::~SiTpc_TrackContainerv1()
{
  Reset();
}

void SiTpc_TrackContainerv1::identify(std::ostream& os) const
{
  os << "SiTpc_TrackContainerv1 with " << m_tracks.size() << " Si+TPC tracks, "
     << m_cand_si.size() << " candidate pairs" << std::endl;
}

void SiTpc_TrackContainerv1::Reset()
{
  for (auto* track : m_tracks)
  {
    delete track;
  }
  m_tracks.clear();
  m_cand_si.clear();
  m_cand_tpc.clear();
  m_cand_dphi.clear();
  m_cand_deta.clear();
  m_cand_dz0.clear();
  m_cand_chi2.clear();
  m_cand_accepted.clear();
}

int SiTpc_TrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* SiTpc_TrackContainerv1::CloneMe() const
{
  auto* copy = new SiTpc_TrackContainerv1();
  for (const auto* track : m_tracks)
  {
    copy->m_tracks.push_back(static_cast<SiTpc_Track*>(track->CloneMe()));
  }
  copy->m_cand_si = m_cand_si;
  copy->m_cand_tpc = m_cand_tpc;
  copy->m_cand_dphi = m_cand_dphi;
  copy->m_cand_deta = m_cand_deta;
  copy->m_cand_dz0 = m_cand_dz0;
  copy->m_cand_chi2 = m_cand_chi2;
  copy->m_cand_accepted = m_cand_accepted;
  return copy;
}
