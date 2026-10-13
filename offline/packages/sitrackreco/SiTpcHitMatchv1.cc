#include "SiTpcHitMatchv1.h"

void SiTpcHitMatchv1::identify(std::ostream& os) const
{
  os << "SiTpcHitMatchv1: id=" << m_id << " tpc_track_index=" << m_tpc_track_index
     << " crossing=" << m_crossing << " rank=" << m_rank << " selected=" << m_selected
     << " n_mvtx_layers=" << m_n_mvtx_layers << " n_intt_layers=" << m_n_intt_layers
     << " n_missing=" << m_n_missing << " chi2=" << m_chi2
     << " steps=" << m_step_layer.size() << " hits=" << m_hits.size()
     << " tpc_pt=" << m_tpc_pt << " delta_phi_dir=" << m_delta_phi_dir
     << " delta_z0=" << m_delta_z0 << std::endl;
}

void SiTpcHitMatchv1::Reset()
{
  *this = SiTpcHitMatchv1();
}

int SiTpcHitMatchv1::isValid() const
{
  return m_hits.empty() ? 0 : 1;
}

void SiTpcHitMatchv1::add_step(const SiTpcHitMatchStep& s)
{
  m_step_layer.push_back(static_cast<int8_t>(s.layer));
  m_step_n_hits.push_back(static_cast<uint16_t>(s.n_hits));
  m_step_flags.push_back(static_cast<uint8_t>(s.flags));
  m_step_hitsetkey.push_back(s.hitsetkey);
  m_step_u.push_back(s.u);
  m_step_z.push_back(s.z);
  m_step_r.push_back(s.r);
  m_step_dphi_tpc.push_back(s.dphi_tpc);
  m_step_dtheta_tpc.push_back(s.dtheta_tpc);
  m_step_dz_tpc.push_back(s.dz_tpc);
  m_step_dphi.push_back(s.dphi);
  m_step_dtheta.push_back(s.dtheta);
  m_step_ddphi.push_back(s.ddphi);
  m_step_ddtheta.push_back(s.ddtheta);
  m_step_chi2.push_back(s.chi2);
  m_step_extent_rphi.push_back(s.extent_rphi);
  m_step_extent_z.push_back(s.extent_z);
  m_step_exp_extent_rphi.push_back(s.exp_extent_rphi);
  m_step_exp_extent_z.push_back(s.exp_extent_z);
}

SiTpcHitMatchStep SiTpcHitMatchv1::get_step(unsigned int i) const
{
  SiTpcHitMatchStep s;
  if (i >= m_step_layer.size())
  {
    return s;
  }
  s.layer = m_step_layer[i];
  s.n_hits = m_step_n_hits[i];
  s.flags = m_step_flags[i];
  s.hitsetkey = m_step_hitsetkey[i];
  s.u = m_step_u[i];
  s.z = m_step_z[i];
  s.r = m_step_r[i];
  s.dphi_tpc = m_step_dphi_tpc[i];
  s.dtheta_tpc = m_step_dtheta_tpc[i];
  s.dz_tpc = m_step_dz_tpc[i];
  s.dphi = m_step_dphi[i];
  s.dtheta = m_step_dtheta[i];
  s.ddphi = m_step_ddphi[i];
  s.ddtheta = m_step_ddtheta[i];
  s.chi2 = m_step_chi2[i];
  s.extent_rphi = m_step_extent_rphi[i];
  s.extent_z = m_step_extent_z[i];
  s.exp_extent_rphi = m_step_exp_extent_rphi[i];
  s.exp_extent_z = m_step_exp_extent_z[i];
  return s;
}