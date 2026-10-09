// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITRAJECTORYV1_H
#define SITRACKRECO_SITRAJECTORYV1_H

#include "Si_Trajectory.h"

#include <iostream>
#include <vector>

class Si_Trajectoryv1 : public Si_Trajectory
{
 public:
  Si_Trajectoryv1();
  ~Si_Trajectoryv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Si_Trajectoryv1(*this); }

  // bookkeeping
  unsigned int get_id() const override { return m_id; }
  int get_chain_id() const override { return m_chain_id; }
  unsigned int get_n_mvtx() const override { return m_n_mvtx; }
  unsigned int get_n_intt() const override { return m_n_intt; }
  int get_fit_status() const override { return m_fit_status; }

  void set_id(unsigned int v) override { m_id = v; }
  void set_chain_id(int v) override { m_chain_id = v; }
  void set_n_mvtx(unsigned int v) override { m_n_mvtx = v; }
  void set_n_intt(unsigned int v) override { m_n_intt = v; }
  void set_fit_status(int v) override { m_fit_status = v; }

  // fit points
  void add_point(int layer, float x, float y, float z, unsigned int nhits) override
  {
    m_point_layer.push_back(layer);
    m_point_x.push_back(x);
    m_point_y.push_back(y);
    m_point_z.push_back(z);
    m_point_nhits.push_back(nhits);
  }
  unsigned int size_points() const override { return static_cast<unsigned int>(m_point_layer.size()); }
  int get_point_layer(unsigned int i) const override { return i < m_point_layer.size() ? m_point_layer[i] : -1; }
  float get_point_x(unsigned int i) const override { return i < m_point_x.size() ? m_point_x[i] : NAN; }
  float get_point_y(unsigned int i) const override { return i < m_point_y.size() ? m_point_y[i] : NAN; }
  float get_point_z(unsigned int i) const override { return i < m_point_z.size() ? m_point_z[i] : NAN; }
  unsigned int get_point_nhits(unsigned int i) const override { return i < m_point_nhits.size() ? m_point_nhits[i] : 0; }

  // circle
  double get_circle_x() const override { return m_circle_x; }
  double get_circle_y() const override { return m_circle_y; }
  double get_radius() const override { return m_radius; }
  int get_helicity() const override { return m_helicity; }
  double get_circle_rms() const override { return m_circle_rms; }

  void set_circle_x(double v) override { m_circle_x = v; }
  void set_circle_y(double v) override { m_circle_y = v; }
  void set_radius(double v) override { m_radius = v; }
  void set_helicity(int v) override { m_helicity = v; }
  void set_circle_rms(double v) override { m_circle_rms = v; }

  // line
  double get_z0() const override { return m_z0; }
  double get_tanl() const override { return m_tanl; }
  double get_z_rms() const override { return m_z_rms; }

  void set_z0(double v) override { m_z0 = v; }
  void set_tanl(double v) override { m_tanl = v; }
  void set_z_rms(double v) override { m_z_rms = v; }

  // perigee
  double get_pca_x() const override { return m_pca_x; }
  double get_pca_y() const override { return m_pca_y; }
  double get_phi() const override { return m_phi; }
  double get_dca() const override { return m_dca; }
  double get_pt() const override { return m_pt; }
  int get_charge() const override { return m_charge; }

  void set_pca_x(double v) override { m_pca_x = v; }
  void set_pca_y(double v) override { m_pca_y = v; }
  void set_phi(double v) override { m_phi = v; }
  void set_dca(double v) override { m_dca = v; }
  void set_pt(double v) override { m_pt = v; }
  void set_charge(int v) override { m_charge = v; }

  // hits
  void add_hit_index(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::hitkey hitkey) override
  {
    m_hit_indices.emplace_back(hitsetkey, hitkey);
  }
  unsigned int size_hit_indices() const override { return static_cast<unsigned int>(m_hit_indices.size()); }
  HitIndex get_hit_index(unsigned int i) const override
  {
    if (i >= m_hit_indices.size()) return {0, 0};
    return m_hit_indices[i];
  }
  const std::vector<HitIndex>& get_hit_indices() const override { return m_hit_indices; }

 private:
  unsigned int m_id{0};
  int m_chain_id{-1};
  unsigned int m_n_mvtx{0};
  unsigned int m_n_intt{0};
  int m_fit_status{NotFitted};

  std::vector<int> m_point_layer;
  std::vector<float> m_point_x;
  std::vector<float> m_point_y;
  std::vector<float> m_point_z;
  std::vector<unsigned int> m_point_nhits;

  double m_circle_x{NAN};
  double m_circle_y{NAN};
  double m_radius{NAN};
  int m_helicity{0};
  double m_circle_rms{NAN};

  double m_z0{NAN};
  double m_tanl{NAN};
  double m_z_rms{NAN};

  double m_pca_x{NAN};
  double m_pca_y{NAN};
  double m_phi{NAN};
  double m_dca{NAN};
  double m_pt{NAN};
  int m_charge{0};

  std::vector<HitIndex> m_hit_indices;

  ClassDefOverride(Si_Trajectoryv1, 1)
};

#endif
