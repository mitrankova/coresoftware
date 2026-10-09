#include "SiTpc_Trackv1.h"

SiTpc_Trackv1::SiTpc_Trackv1()
{
  Reset();
}

void SiTpc_Trackv1::identify(std::ostream& os) const
{
  os << "SiTpc_Trackv1: id=" << m_id << " si_chain=" << m_si_chain_id << " (traj " << m_si_trajectory_index << ")"
     << " tpc_track=" << m_tpc_track_id << " (assembled " << m_tpc_assembled_track_id << ")"
     << " points si/tpc=" << m_n_si_points << "/" << m_n_tpc_points << " status=" << m_fit_status << std::endl;
  os << "  match: dphi=" << m_match_dphi << " deta=" << m_match_deta << " dz0=" << m_match_dz0
     << " chi2=" << m_match_chi2 << "  (si phi/eta/z0 " << m_si_phi << "/" << m_si_eta << "/" << m_si_z0
     << ", tpc phi/eta/z0/pt " << m_tpc_phi << "/" << m_tpc_eta << "/" << m_tpc_z0 << "/" << m_tpc_pt << ")" << std::endl;
  os << "  fit: c=(" << m_circle_x << ", " << m_circle_y << ") R=" << m_radius << " h=" << m_helicity
     << " pca=(" << m_pca_x << ", " << m_pca_y << ") phi=" << m_phi << " dca=" << m_dca << " pt=" << m_pt
     << " q=" << m_charge << " z0=" << m_z0 << " tanl=" << m_tanl << " tpc_z_offset=" << m_tpc_z_offset
     << " rms xy/z=" << m_circle_rms << "/" << m_z_rms << std::endl;
}

void SiTpc_Trackv1::Reset()
{
  m_id = 0;
  m_si_trajectory_index = 0;
  m_si_chain_id = -1;
  m_tpc_track_index = 0;
  m_tpc_track_id = 0;
  m_tpc_assembled_track_id = 0;
  m_n_si_points = 0;
  m_n_tpc_points = 0;
  m_fit_status = 0;
  m_si_phi = NAN;
  m_si_eta = NAN;
  m_si_z0 = NAN;
  m_tpc_phi = NAN;
  m_tpc_eta = NAN;
  m_tpc_z0 = NAN;
  m_tpc_pt = NAN;
  m_tpc_charge = 0;
  m_match_dphi = NAN;
  m_match_deta = NAN;
  m_match_dz0 = NAN;
  m_match_chi2 = NAN;
  m_circle_x = NAN;
  m_circle_y = NAN;
  m_radius = NAN;
  m_helicity = 0;
  m_circle_rms = NAN;
  m_pca_x = NAN;
  m_pca_y = NAN;
  m_phi = NAN;
  m_dca = NAN;
  m_pt = NAN;
  m_charge = 0;
  m_z0 = NAN;
  m_tanl = NAN;
  m_tpc_z_offset = 0;
  m_z_rms = NAN;
  m_point_type.clear();
  m_point_layer.clear();
  m_point_source.clear();
  m_point_x.clear();
  m_point_y.clear();
  m_point_z.clear();
}

int SiTpc_Trackv1::isValid() const
{
  return (m_fit_status == 1 || m_fit_status == 4) ? 1 : 0;  // SiTpcHelixFit::Ok / StraightLine
}
