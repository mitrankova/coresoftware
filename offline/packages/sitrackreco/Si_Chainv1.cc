#include "Si_Chainv1.h"

Si_Chainv1::Si_Chainv1()
{
  Reset();
}

void Si_Chainv1::identify(std::ostream& os) const
{
  os << "Si_Chainv1:"
     << " chain_id=" << m_chain_id
     << " n_mvtx_hits=" << m_n_mvtx_hits
     << " n_intt_hits=" << m_n_intt_hits
     << " vertex=" << m_vertex
     << " slope=" << m_slope
     << " phi=" << m_phi
     << " hit_indices=" << m_hit_indices.size()
     << std::endl;
}

void Si_Chainv1::Reset()
{
  m_chain_id = 0;
  m_n_mvtx_hits = 0;
  m_n_intt_hits = 0;
  m_vertex = 0.0;
  m_slope = 0.0;
  m_phi = 0.0;
  m_hit_indices.clear();
}

int Si_Chainv1::isValid() const
{
  return m_hit_indices.empty() ? 0 : 1;
}
