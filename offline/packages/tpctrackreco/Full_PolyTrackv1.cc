#include "Full_PolyTrackv1.h"

#include <cmath>

Full_PolyTrackv1::Full_PolyTrackv1()
{
  Reset();
}

void Full_PolyTrackv1::identify(std::ostream& os) const
{
  os << "Full_PolyTrackv1:"
     << " event=" << m_event
     << " track_id=" << m_track_id
     << " tpc_poly_track_id=" << m_tpc_poly_track_id
     << " source_assembled_track_id=" << m_source_assembled_track_id
     << " n_tpc_clusters=" << m_n_tpc_clusters
     << " n_silicon_clusters=" << m_silicon_cluster_keys.size()
     << " n_missing_layers=" << m_n_missing_layers
     << " chi2=" << m_chi2
     << " ndf=" << m_ndf
     << " score=" << m_score
     << " pos=(" << m_x << "," << m_y << "," << m_z << ")"
     << " mom=(" << m_px << "," << m_py << "," << m_pz << ")"
     << std::endl;
}

void Full_PolyTrackv1::Reset()
{
  m_event = 0;
  m_track_id = 0;
  m_tpc_poly_track_id = 0;
  m_source_assembled_track_id = 0;
  m_n_tpc_clusters = 0;
  m_n_missing_layers = 0;
  m_missing_layer_mask = 0;
  m_fit_status = 0;
  m_chi2 = 0.0;
  m_ndf = 0.0;
  m_score = 0.0;
  m_x = 0.0;
  m_y = 0.0;
  m_z = 0.0;
  m_px = 0.0;
  m_py = 0.0;
  m_pz = 0.0;
  m_charge = 0.0;
  m_seed_x = 0.0;
  m_seed_y = 0.0;
  m_seed_z = 0.0;
  m_seed_px = 0.0;
  m_seed_py = 0.0;
  m_seed_pz = 0.0;
  clear_cluster_keys();
  clear_silicon_states();
}

int Full_PolyTrackv1::isValid() const
{
  return (m_fit_status != 0 &&
          std::isfinite(m_x) && std::isfinite(m_y) && std::isfinite(m_z))
             ? 1
             : 0;
}

void Full_PolyTrackv1::add_silicon_state(unsigned int layer, TrkrDefs::cluskey key,
                                         float x, float y, float z,
                                         float pred_x, float pred_y, float pred_z,
                                         float rdphi, float dz, float chi2)
{
  m_state_layer.push_back(layer);
  m_state_cluster_key.push_back(key);
  m_state_x.push_back(x);
  m_state_y.push_back(y);
  m_state_z.push_back(z);
  m_state_pred_x.push_back(pred_x);
  m_state_pred_y.push_back(pred_y);
  m_state_pred_z.push_back(pred_z);
  m_state_rdphi.push_back(rdphi);
  m_state_dz.push_back(dz);
  m_state_chi2.push_back(chi2);
}

void Full_PolyTrackv1::clear_silicon_states()
{
  m_state_layer.clear();
  m_state_cluster_key.clear();
  m_state_x.clear();
  m_state_y.clear();
  m_state_z.clear();
  m_state_pred_x.clear();
  m_state_pred_y.clear();
  m_state_pred_z.clear();
  m_state_rdphi.clear();
  m_state_dz.clear();
  m_state_chi2.clear();
}
