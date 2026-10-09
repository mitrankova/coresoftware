// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCTRACKV1_H
#define SITRACKRECO_SITPCTRACKV1_H

#include "SiTpc_Track.h"

#include <iostream>
#include <vector>

class SiTpc_Trackv1 : public SiTpc_Track
{
 public:
  SiTpc_Trackv1();
  ~SiTpc_Trackv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new SiTpc_Trackv1(*this); }

  unsigned int get_id() const override { return m_id; }
  void set_id(unsigned int v) override { m_id = v; }
  unsigned int get_si_trajectory_index() const override { return m_si_trajectory_index; }
  void set_si_trajectory_index(unsigned int v) override { m_si_trajectory_index = v; }
  int get_si_chain_id() const override { return m_si_chain_id; }
  void set_si_chain_id(int v) override { m_si_chain_id = v; }
  unsigned int get_tpc_track_index() const override { return m_tpc_track_index; }
  void set_tpc_track_index(unsigned int v) override { m_tpc_track_index = v; }
  unsigned int get_tpc_track_id() const override { return m_tpc_track_id; }
  void set_tpc_track_id(unsigned int v) override { m_tpc_track_id = v; }
  unsigned int get_tpc_assembled_track_id() const override { return m_tpc_assembled_track_id; }
  void set_tpc_assembled_track_id(unsigned int v) override { m_tpc_assembled_track_id = v; }
  unsigned int get_n_si_points() const override { return m_n_si_points; }
  void set_n_si_points(unsigned int v) override { m_n_si_points = v; }
  unsigned int get_n_tpc_points() const override { return m_n_tpc_points; }
  void set_n_tpc_points(unsigned int v) override { m_n_tpc_points = v; }
  int get_fit_status() const override { return m_fit_status; }
  void set_fit_status(int v) override { m_fit_status = v; }
  double get_si_phi() const override { return m_si_phi; }
  void set_si_phi(double v) override { m_si_phi = v; }
  double get_si_eta() const override { return m_si_eta; }
  void set_si_eta(double v) override { m_si_eta = v; }
  double get_si_z0() const override { return m_si_z0; }
  void set_si_z0(double v) override { m_si_z0 = v; }
  double get_tpc_phi() const override { return m_tpc_phi; }
  void set_tpc_phi(double v) override { m_tpc_phi = v; }
  double get_tpc_eta() const override { return m_tpc_eta; }
  void set_tpc_eta(double v) override { m_tpc_eta = v; }
  double get_tpc_z0() const override { return m_tpc_z0; }
  void set_tpc_z0(double v) override { m_tpc_z0 = v; }
  double get_tpc_pt() const override { return m_tpc_pt; }
  void set_tpc_pt(double v) override { m_tpc_pt = v; }
  int get_tpc_charge() const override { return m_tpc_charge; }
  void set_tpc_charge(int v) override { m_tpc_charge = v; }
  double get_match_dphi() const override { return m_match_dphi; }
  void set_match_dphi(double v) override { m_match_dphi = v; }
  double get_match_deta() const override { return m_match_deta; }
  void set_match_deta(double v) override { m_match_deta = v; }
  double get_match_dz0() const override { return m_match_dz0; }
  void set_match_dz0(double v) override { m_match_dz0 = v; }
  double get_match_chi2() const override { return m_match_chi2; }
  void set_match_chi2(double v) override { m_match_chi2 = v; }
  double get_circle_x() const override { return m_circle_x; }
  void set_circle_x(double v) override { m_circle_x = v; }
  double get_circle_y() const override { return m_circle_y; }
  void set_circle_y(double v) override { m_circle_y = v; }
  double get_radius() const override { return m_radius; }
  void set_radius(double v) override { m_radius = v; }
  int get_helicity() const override { return m_helicity; }
  void set_helicity(int v) override { m_helicity = v; }
  double get_circle_rms() const override { return m_circle_rms; }
  void set_circle_rms(double v) override { m_circle_rms = v; }
  double get_pca_x() const override { return m_pca_x; }
  void set_pca_x(double v) override { m_pca_x = v; }
  double get_pca_y() const override { return m_pca_y; }
  void set_pca_y(double v) override { m_pca_y = v; }
  double get_phi() const override { return m_phi; }
  void set_phi(double v) override { m_phi = v; }
  double get_dca() const override { return m_dca; }
  void set_dca(double v) override { m_dca = v; }
  double get_pt() const override { return m_pt; }
  void set_pt(double v) override { m_pt = v; }
  int get_charge() const override { return m_charge; }
  void set_charge(int v) override { m_charge = v; }
  double get_z0() const override { return m_z0; }
  void set_z0(double v) override { m_z0 = v; }
  double get_tanl() const override { return m_tanl; }
  void set_tanl(double v) override { m_tanl = v; }
  double get_tpc_z_offset() const override { return m_tpc_z_offset; }
  void set_tpc_z_offset(double v) override { m_tpc_z_offset = v; }
  double get_z_rms() const override { return m_z_rms; }
  void set_z_rms(double v) override { m_z_rms = v; }

  void add_point(int type, int layer, unsigned int source, float x, float y, float z) override
  {
    m_point_type.push_back(type);
    m_point_layer.push_back(layer);
    m_point_source.push_back(source);
    m_point_x.push_back(x);
    m_point_y.push_back(y);
    m_point_z.push_back(z);
  }
  unsigned int size_points() const override { return static_cast<unsigned int>(m_point_type.size()); }
  int get_point_type(unsigned int i) const override { return i < m_point_type.size() ? m_point_type[i] : -1; }
  int get_point_layer(unsigned int i) const override { return i < m_point_layer.size() ? m_point_layer[i] : -1; }
  unsigned int get_point_source(unsigned int i) const override { return i < m_point_source.size() ? m_point_source[i] : 0; }
  float get_point_x(unsigned int i) const override { return i < m_point_x.size() ? m_point_x[i] : NAN; }
  float get_point_y(unsigned int i) const override { return i < m_point_y.size() ? m_point_y[i] : NAN; }
  float get_point_z(unsigned int i) const override { return i < m_point_z.size() ? m_point_z[i] : NAN; }

 private:
  unsigned int m_id{0};
  unsigned int m_si_trajectory_index{0};
  int m_si_chain_id{-1};
  unsigned int m_tpc_track_index{0};
  unsigned int m_tpc_track_id{0};
  unsigned int m_tpc_assembled_track_id{0};
  unsigned int m_n_si_points{0};
  unsigned int m_n_tpc_points{0};
  int m_fit_status{0};
  double m_si_phi{NAN};
  double m_si_eta{NAN};
  double m_si_z0{NAN};
  double m_tpc_phi{NAN};
  double m_tpc_eta{NAN};
  double m_tpc_z0{NAN};
  double m_tpc_pt{NAN};
  int m_tpc_charge{0};
  double m_match_dphi{NAN};
  double m_match_deta{NAN};
  double m_match_dz0{NAN};
  double m_match_chi2{NAN};
  double m_circle_x{NAN};
  double m_circle_y{NAN};
  double m_radius{NAN};
  int m_helicity{0};
  double m_circle_rms{NAN};
  double m_pca_x{NAN};
  double m_pca_y{NAN};
  double m_phi{NAN};
  double m_dca{NAN};
  double m_pt{NAN};
  int m_charge{0};
  double m_z0{NAN};
  double m_tanl{NAN};
  double m_tpc_z_offset{0};
  double m_z_rms{NAN};

  std::vector<int> m_point_type;
  std::vector<int> m_point_layer;
  std::vector<unsigned int> m_point_source;
  std::vector<float> m_point_x;
  std::vector<float> m_point_y;
  std::vector<float> m_point_z;

  ClassDefOverride(SiTpc_Trackv1, 1)
};

#endif
