// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_FULLPOLYTRACKV1_H
#define TPCTRACKRECO_FULLPOLYTRACKV1_H

#include "Full_PolyTrack.h"

#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <vector>

class Full_PolyTrackv1 : public Full_PolyTrack
{
 public:
  Full_PolyTrackv1();
  ~Full_PolyTrackv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Full_PolyTrackv1(*this); }

  unsigned int get_event() const override { return m_event; }
  unsigned int get_track_id() const override { return m_track_id; }
  unsigned int get_tpc_poly_track_id() const override { return m_tpc_poly_track_id; }
  unsigned int get_source_assembled_track_id() const override { return m_source_assembled_track_id; }
  unsigned int get_n_tpc_clusters() const override { return m_n_tpc_clusters; }
  unsigned int get_n_silicon_clusters() const override { return static_cast<unsigned int>(m_silicon_cluster_keys.size()); }
  unsigned int get_n_missing_layers() const override { return m_n_missing_layers; }
  unsigned int get_missing_layer_mask() const override { return m_missing_layer_mask; }
  int get_fit_status() const override { return m_fit_status; }
  double get_chi2() const override { return m_chi2; }
  double get_ndf() const override { return m_ndf; }
  double get_score() const override { return m_score; }
  double get_x() const override { return m_x; }
  double get_y() const override { return m_y; }
  double get_z() const override { return m_z; }
  double get_px() const override { return m_px; }
  double get_py() const override { return m_py; }
  double get_pz() const override { return m_pz; }
  double get_charge() const override { return m_charge; }
  double get_seed_x() const override { return m_seed_x; }
  double get_seed_y() const override { return m_seed_y; }
  double get_seed_z() const override { return m_seed_z; }
  double get_seed_px() const override { return m_seed_px; }
  double get_seed_py() const override { return m_seed_py; }
  double get_seed_pz() const override { return m_seed_pz; }

  void set_event(unsigned int v) override { m_event = v; }
  void set_track_id(unsigned int v) override { m_track_id = v; }
  void set_tpc_poly_track_id(unsigned int v) override { m_tpc_poly_track_id = v; }
  void set_source_assembled_track_id(unsigned int v) override { m_source_assembled_track_id = v; }
  void set_n_tpc_clusters(unsigned int v) override { m_n_tpc_clusters = v; }
  void set_n_missing_layers(unsigned int v) override { m_n_missing_layers = v; }
  void set_missing_layer_mask(unsigned int v) override { m_missing_layer_mask = v; }
  void set_fit_status(int v) override { m_fit_status = v; }
  void set_chi2(double v) override { m_chi2 = v; }
  void set_ndf(double v) override { m_ndf = v; }
  void set_score(double v) override { m_score = v; }
  void set_x(double v) override { m_x = v; }
  void set_y(double v) override { m_y = v; }
  void set_z(double v) override { m_z = v; }
  void set_px(double v) override { m_px = v; }
  void set_py(double v) override { m_py = v; }
  void set_pz(double v) override { m_pz = v; }
  void set_charge(double v) override { m_charge = v; }
  void set_seed_x(double v) override { m_seed_x = v; }
  void set_seed_y(double v) override { m_seed_y = v; }
  void set_seed_z(double v) override { m_seed_z = v; }
  void set_seed_px(double v) override { m_seed_px = v; }
  void set_seed_py(double v) override { m_seed_py = v; }
  void set_seed_pz(double v) override { m_seed_pz = v; }

  unsigned int size_cluster_keys() const override { return static_cast<unsigned int>(m_silicon_cluster_keys.size()); }
  TrkrDefs::cluskey get_cluster_key(unsigned int i) const override
  {
    return i < m_silicon_cluster_keys.size() ? m_silicon_cluster_keys[i] : TrkrDefs::CLUSKEYMAX;
  }
  const std::vector<TrkrDefs::cluskey>& get_cluster_keys() const override { return m_silicon_cluster_keys; }
  void add_cluster_key(TrkrDefs::cluskey key) override { m_silicon_cluster_keys.push_back(key); }
  void clear_cluster_keys() override { m_silicon_cluster_keys.clear(); }

  unsigned int size_silicon_states() const override { return static_cast<unsigned int>(m_state_layer.size()); }
  unsigned int get_state_layer(unsigned int i) const override { return i < m_state_layer.size() ? m_state_layer[i] : 0; }
  TrkrDefs::cluskey get_state_cluster_key(unsigned int i) const override { return i < m_state_cluster_key.size() ? m_state_cluster_key[i] : TrkrDefs::CLUSKEYMAX; }
  float get_state_x(unsigned int i) const override { return i < m_state_x.size() ? m_state_x[i] : 0.0F; }
  float get_state_y(unsigned int i) const override { return i < m_state_y.size() ? m_state_y[i] : 0.0F; }
  float get_state_z(unsigned int i) const override { return i < m_state_z.size() ? m_state_z[i] : 0.0F; }
  float get_state_pred_x(unsigned int i) const override { return i < m_state_pred_x.size() ? m_state_pred_x[i] : 0.0F; }
  float get_state_pred_y(unsigned int i) const override { return i < m_state_pred_y.size() ? m_state_pred_y[i] : 0.0F; }
  float get_state_pred_z(unsigned int i) const override { return i < m_state_pred_z.size() ? m_state_pred_z[i] : 0.0F; }
  float get_state_rdphi(unsigned int i) const override { return i < m_state_rdphi.size() ? m_state_rdphi[i] : 0.0F; }
  float get_state_dz(unsigned int i) const override { return i < m_state_dz.size() ? m_state_dz[i] : 0.0F; }
  float get_state_chi2(unsigned int i) const override { return i < m_state_chi2.size() ? m_state_chi2[i] : 0.0F; }
  void add_silicon_state(unsigned int layer, TrkrDefs::cluskey key,
                         float x, float y, float z,
                         float pred_x, float pred_y, float pred_z,
                         float rdphi, float dz, float chi2) override;
  void clear_silicon_states() override;

 private:
  unsigned int m_event{0};
  unsigned int m_track_id{0};
  unsigned int m_tpc_poly_track_id{0};
  unsigned int m_source_assembled_track_id{0};
  unsigned int m_n_tpc_clusters{0};
  unsigned int m_n_missing_layers{0};
  unsigned int m_missing_layer_mask{0};
  int m_fit_status{0};
  double m_chi2{0.0};
  double m_ndf{0.0};
  double m_score{0.0};
  double m_x{0.0};
  double m_y{0.0};
  double m_z{0.0};
  double m_px{0.0};
  double m_py{0.0};
  double m_pz{0.0};
  double m_charge{0.0};
  double m_seed_x{0.0};
  double m_seed_y{0.0};
  double m_seed_z{0.0};
  double m_seed_px{0.0};
  double m_seed_py{0.0};
  double m_seed_pz{0.0};
  std::vector<TrkrDefs::cluskey> m_silicon_cluster_keys;
  std::vector<unsigned int> m_state_layer;
  std::vector<TrkrDefs::cluskey> m_state_cluster_key;
  std::vector<float> m_state_x;
  std::vector<float> m_state_y;
  std::vector<float> m_state_z;
  std::vector<float> m_state_pred_x;
  std::vector<float> m_state_pred_y;
  std::vector<float> m_state_pred_z;
  std::vector<float> m_state_rdphi;
  std::vector<float> m_state_dz;
  std::vector<float> m_state_chi2;

  ClassDefOverride(Full_PolyTrackv1, 1)
};
#endif
