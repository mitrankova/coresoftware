// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SICHAINV1_H
#define SITRACKRECO_SICHAINV1_H

#include "Si_Chain.h"

#include <iostream>
#include <vector>

class Si_Chainv1 : public Si_Chain
{
 public:
  Si_Chainv1();
  ~Si_Chainv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Si_Chainv1(*this); }

  unsigned int get_chain_id() const override { return m_chain_id; }
  unsigned int get_n_mvtx_hits() const override { return m_n_mvtx_hits; }
  unsigned int get_n_intt_hits() const override { return m_n_intt_hits; }
  double get_vertex() const override { return m_vertex; }
  double get_slope() const override { return m_slope; }
  double get_phi() const override { return m_phi; }

  void set_chain_id(unsigned int value) override { m_chain_id = value; }
  void set_n_mvtx_hits(unsigned int value) override { m_n_mvtx_hits = value; }
  void set_n_intt_hits(unsigned int value) override { m_n_intt_hits = value; }
  void set_vertex(double value) override { m_vertex = value; }
  void set_slope(double value) override { m_slope = value; }
  void set_phi(double value) override { m_phi = value; }

  void add_hit_index(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::hitkey hitkey) override
  {
    m_hit_indices.emplace_back(hitsetkey, hitkey);
  }

  unsigned int size_hit_indices() const override
  {
    return static_cast<unsigned int>(m_hit_indices.size());
  }

  HitIndex get_hit_index(unsigned int index) const override
  {
    if (index >= m_hit_indices.size()) return {0, 0};
    return m_hit_indices[index];
  }

  const std::vector<HitIndex>& get_hit_indices() const override
  {
    return m_hit_indices;
  }

 private:
  unsigned int m_chain_id{0};
  unsigned int m_n_mvtx_hits{0};
  unsigned int m_n_intt_hits{0};
  double m_vertex{0.0};
  double m_slope{0.0};
  double m_phi{0.0};
  std::vector<HitIndex> m_hit_indices;

  ClassDefOverride(Si_Chainv1, 1)
};

#endif
