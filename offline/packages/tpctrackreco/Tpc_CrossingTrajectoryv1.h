// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCCROSSINGTRAJECTORYV1_H
#define TPCTRACKRECO_TPCCROSSINGTRAJECTORYV1_H

#include "Tpc_CrossingTrajectory.h"

#include <iostream>

class Tpc_CrossingTrajectoryv1 : public Tpc_CrossingTrajectory
{
 public:
  Tpc_CrossingTrajectoryv1() = default;
  ~Tpc_CrossingTrajectoryv1() override = default;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = Tpc_CrossingTrajectoryv1(); }
  int isValid() const override { return m_fit_status != 0 ? 1 : 0; }
  PHObject* CloneMe() const override { return new Tpc_CrossingTrajectoryv1(*this); }

  unsigned int get_event() const override { return m_event; }
  unsigned int get_trajectory_id() const override { return m_trajectory_id; }
  unsigned int get_source_assembled_track_id() const override { return m_source_assembled_track_id; }
  short get_crossing() const override { return m_crossing; }
  int get_side() const override { return m_side; }
  int get_fit_status() const override { return m_fit_status; }

  unsigned int get_point_layer(unsigned int i) const override { return i < NPoints ? m_point_layer[i] : 0; }
  double get_point_x(unsigned int i) const override { return i < NPoints ? m_point_x[i] : 0.0; }
  double get_point_y(unsigned int i) const override { return i < NPoints ? m_point_y[i] : 0.0; }
  double get_point_z(unsigned int i) const override { return i < NPoints ? m_point_z[i] : 0.0; }

  double get_x() const override { return m_x; }
  double get_y() const override { return m_y; }
  double get_z() const override { return m_z; }
  double get_px() const override { return m_px; }
  double get_py() const override { return m_py; }
  double get_pz() const override { return m_pz; }
  double get_charge() const override { return m_charge; }
  double get_q_over_r() const override { return m_q_over_r; }
  double get_phi() const override { return m_phi; }
  double get_slope() const override { return m_slope; }
  double get_chi2() const override { return m_chi2; }

  void set_event(unsigned int v) override { m_event = v; }
  void set_trajectory_id(unsigned int v) override { m_trajectory_id = v; }
  void set_source_assembled_track_id(unsigned int v) override { m_source_assembled_track_id = v; }
  void set_crossing(short v) override { m_crossing = v; }
  void set_side(int v) override { m_side = v; }
  void set_fit_status(int v) override { m_fit_status = v; }
  void set_point(unsigned int i, unsigned int layer, double x, double y, double z) override
  {
    if (i >= NPoints) return;
    m_point_layer[i] = layer;
    m_point_x[i] = x;
    m_point_y[i] = y;
    m_point_z[i] = z;
  }
  void set_x(double v) override { m_x = v; }
  void set_y(double v) override { m_y = v; }
  void set_z(double v) override { m_z = v; }
  void set_px(double v) override { m_px = v; }
  void set_py(double v) override { m_py = v; }
  void set_pz(double v) override { m_pz = v; }
  void set_charge(double v) override { m_charge = v; }
  void set_q_over_r(double v) override { m_q_over_r = v; }
  void set_phi(double v) override { m_phi = v; }
  void set_slope(double v) override { m_slope = v; }
  void set_chi2(double v) override { m_chi2 = v; }

 private:
  unsigned int m_event{0};
  unsigned int m_trajectory_id{0};
  unsigned int m_source_assembled_track_id{0};
  short m_crossing{0};
  int m_side{0};
  int m_fit_status{0};
  unsigned int m_point_layer[NPoints]{0, 0, 0};
  double m_point_x[NPoints]{0.0, 0.0, 0.0};
  double m_point_y[NPoints]{0.0, 0.0, 0.0};
  double m_point_z[NPoints]{0.0, 0.0, 0.0};
  double m_x{0.0};
  double m_y{0.0};
  double m_z{0.0};
  double m_px{0.0};
  double m_py{0.0};
  double m_pz{0.0};
  double m_charge{0.0};
  double m_q_over_r{0.0};
  double m_phi{0.0};
  double m_slope{0.0};
  double m_chi2{0.0};

  ClassDefOverride(Tpc_CrossingTrajectoryv1, 1)
};
#endif
