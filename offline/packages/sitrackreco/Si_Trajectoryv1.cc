#include "Si_Trajectoryv1.h"

Si_Trajectoryv1::Si_Trajectoryv1()
{
  Reset();
}

void Si_Trajectoryv1::identify(std::ostream& os) const
{
  os << "Si_Trajectoryv1:"
     << " id=" << m_id
     << " chain_id=" << m_chain_id
     << " n_mvtx=" << m_n_mvtx
     << " n_intt=" << m_n_intt
     << " status=" << m_fit_status
     << " points=" << m_point_layer.size()
     << " hits=" << m_hit_indices.size()
     << std::endl;
  os << "  circle: c=(" << m_circle_x << ", " << m_circle_y << ") R=" << m_radius
     << " h=" << m_helicity << " rms=" << m_circle_rms << std::endl;
  os << "  perigee: pca=(" << m_pca_x << ", " << m_pca_y << ") phi=" << m_phi
     << " dca=" << m_dca << " z0=" << m_z0 << " tanl=" << m_tanl
     << " pt=" << m_pt << " q=" << m_charge << " z_rms=" << m_z_rms << std::endl;
  for (unsigned int i = 0; i < m_point_layer.size(); ++i)
  {
    os << "  layer " << m_point_layer[i] << ": (" << m_point_x[i] << ", " << m_point_y[i]
       << ", " << m_point_z[i] << ") nhits=" << m_point_nhits[i] << std::endl;
  }
}

void Si_Trajectoryv1::Reset()
{
  m_id = 0;
  m_chain_id = -1;
  m_n_mvtx = 0;
  m_n_intt = 0;
  m_fit_status = NotFitted;

  m_point_layer.clear();
  m_point_x.clear();
  m_point_y.clear();
  m_point_z.clear();
  m_point_nhits.clear();

  m_circle_x = NAN;
  m_circle_y = NAN;
  m_radius = NAN;
  m_helicity = 0;
  m_circle_rms = NAN;

  m_z0 = NAN;
  m_tanl = NAN;
  m_z_rms = NAN;

  m_pca_x = NAN;
  m_pca_y = NAN;
  m_phi = NAN;
  m_dca = NAN;
  m_pt = NAN;
  m_charge = 0;

  m_hit_indices.clear();
}

int Si_Trajectoryv1::isValid() const
{
  return (m_fit_status == Ok || m_fit_status == StraightLine) ? 1 : 0;
}
