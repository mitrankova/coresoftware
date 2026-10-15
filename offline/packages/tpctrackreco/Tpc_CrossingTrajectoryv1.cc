#include "Tpc_CrossingTrajectoryv1.h"

void Tpc_CrossingTrajectoryv1::identify(std::ostream& os) const
{
  os << "Tpc_CrossingTrajectoryv1:"
     << " event=" << m_event
     << " id=" << m_trajectory_id
     << " assembled=" << m_source_assembled_track_id
     << " crossing=" << m_crossing
     << " side=" << m_side
     << " fit_status=" << m_fit_status
     << " layers=(" << m_point_layer[0] << "," << m_point_layer[1] << "," << m_point_layer[2] << ")"
     << " pos=(" << m_x << "," << m_y << "," << m_z << ")"
     << " mom=(" << m_px << "," << m_py << "," << m_pz << ")"
     << " charge=" << m_charge
     << " qOverR=" << m_q_over_r
     << std::endl;
}
