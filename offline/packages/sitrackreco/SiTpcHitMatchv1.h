// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCHITMATCHV1_H
#define SITRACKRECO_SITPCHITMATCHV1_H

#include "SiTpcHitMatch.h"

#include <cstdint>
#include <iostream>
#include <vector>

class SiTpcHitMatchv1 : public SiTpcHitMatch
{
 public:
  SiTpcHitMatchv1() = default;
  ~SiTpcHitMatchv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new SiTpcHitMatchv1(*this); }

  // ---- bookkeeping
  unsigned int get_id() const override { return m_id; }
  unsigned int get_tpc_track_index() const override { return m_tpc_track_index; }
  unsigned int get_tpc_track_id() const override { return m_tpc_track_id; }
  unsigned int get_tpc_assembled_track_id() const override { return m_tpc_assembled_track_id; }
  short get_crossing() const override { return m_crossing; }
  unsigned int get_rank() const override { return m_rank; }
  bool get_selected() const override { return m_selected; }

  void set_id(unsigned int v) override { m_id = v; }
  void set_tpc_track_index(unsigned int v) override { m_tpc_track_index = v; }
  void set_tpc_track_id(unsigned int v) override { m_tpc_track_id = v; }
  void set_tpc_assembled_track_id(unsigned int v) override { m_tpc_assembled_track_id = v; }
  void set_crossing(short v) override { m_crossing = v; }
  void set_rank(unsigned int v) override { m_rank = v; }
  void set_selected(bool v) override { m_selected = v; }

  // ---- TPC reference
  float get_tpc_pt() const override { return m_tpc_pt; }
  int get_tpc_charge() const override { return m_tpc_charge; }
  float get_tpc_phi() const override { return m_tpc_phi; }
  float get_tpc_tanl() const override { return m_tpc_tanl; }
  float get_tpc_z0() const override { return m_tpc_z0; }
  float get_tpc_dca() const override { return m_tpc_dca; }

  void set_tpc_pt(float v) override { m_tpc_pt = v; }
  void set_tpc_charge(int v) override { m_tpc_charge = v; }
  void set_tpc_phi(float v) override { m_tpc_phi = v; }
  void set_tpc_tanl(float v) override { m_tpc_tanl = v; }
  void set_tpc_z0(float v) override { m_tpc_z0 = v; }
  void set_tpc_dca(float v) override { m_tpc_dca = v; }

  // ---- summary
  unsigned int get_n_mvtx_layers() const override { return m_n_mvtx_layers; }
  unsigned int get_n_intt_layers() const override { return m_n_intt_layers; }
  unsigned int get_n_missing() const override { return m_n_missing; }
  float get_chi2() const override { return m_chi2; }
  float get_si_z0() const override { return m_si_z0; }
  float get_si_dca() const override { return m_si_dca; }
  float get_delta_phi_dir() const override { return m_delta_phi_dir; }
  float get_delta_z0() const override { return m_delta_z0; }
  float get_delta_dca() const override { return m_delta_dca; }

  void set_n_mvtx_layers(unsigned int v) override { m_n_mvtx_layers = v; }
  void set_n_intt_layers(unsigned int v) override { m_n_intt_layers = v; }
  void set_n_missing(unsigned int v) override { m_n_missing = v; }
  void set_chi2(float v) override { m_chi2 = v; }
  void set_si_z0(float v) override { m_si_z0 = v; }
  void set_si_dca(float v) override { m_si_dca = v; }
  void set_delta_phi_dir(float v) override { m_delta_phi_dir = v; }
  void set_delta_z0(float v) override { m_delta_z0 = v; }
  void set_delta_dca(float v) override { m_delta_dca = v; }

  // ---- steps
  void add_step(const SiTpcHitMatchStep& s) override;
  unsigned int size_steps() const override { return static_cast<unsigned int>(m_step_layer.size()); }
  SiTpcHitMatchStep get_step(unsigned int i) const override;

  // ---- hits
  void add_hit_index(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, unsigned int step) override
  {
    m_hits.emplace_back(hsk, hk);
    m_hit_step.push_back(static_cast<uint8_t>(step));
  }
  unsigned int size_hit_indices() const override { return static_cast<unsigned int>(m_hits.size()); }
  HitIndex get_hit_index(unsigned int i) const override { return i < m_hits.size() ? m_hits[i] : HitIndex{0, 0}; }
  unsigned int get_hit_step(unsigned int i) const override { return i < m_hit_step.size() ? m_hit_step[i] : 0; }
  const std::vector<HitIndex>& get_hit_indices() const override { return m_hits; }

 private:
  unsigned int m_id{0};
  unsigned int m_tpc_track_index{0};
  unsigned int m_tpc_track_id{0};
  unsigned int m_tpc_assembled_track_id{0};
  short m_crossing{0};
  unsigned int m_rank{0};
  bool m_selected{false};

  float m_tpc_pt{NAN};
  int m_tpc_charge{0};
  float m_tpc_phi{NAN};
  float m_tpc_tanl{NAN};
  float m_tpc_z0{NAN};
  float m_tpc_dca{NAN};

  unsigned int m_n_mvtx_layers{0};
  unsigned int m_n_intt_layers{0};
  unsigned int m_n_missing{0};
  float m_chi2{NAN};
  float m_si_z0{NAN};
  float m_si_dca{NAN};
  float m_delta_phi_dir{NAN};
  float m_delta_z0{NAN};
  float m_delta_dca{NAN};

  // steps, column-wise
  std::vector<int8_t> m_step_layer;
  std::vector<uint16_t> m_step_n_hits;
  std::vector<uint8_t> m_step_flags;
  std::vector<unsigned int> m_step_hitsetkey;
  std::vector<float> m_step_u;
  std::vector<float> m_step_z;
  std::vector<float> m_step_r;
  std::vector<float> m_step_dphi_tpc;
  std::vector<float> m_step_dtheta_tpc;
  std::vector<float> m_step_dz_tpc;
  std::vector<float> m_step_dphi;
  std::vector<float> m_step_dtheta;
  std::vector<float> m_step_ddphi;
  std::vector<float> m_step_ddtheta;
  std::vector<float> m_step_chi2;
  std::vector<float> m_step_extent_rphi;
  std::vector<float> m_step_extent_z;
  std::vector<float> m_step_exp_extent_rphi;
  std::vector<float> m_step_exp_extent_z;

  std::vector<HitIndex> m_hits;
  std::vector<uint8_t> m_hit_step;

  ClassDefOverride(SiTpcHitMatchv1, 1)
};

#endif
