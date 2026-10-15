// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCCROSSINGTRAJECTORY_H
#define TPCTRACKRECO_TPCCROSSINGTRAJECTORY_H

#include <phool/PHObject.h>

#include <iostream>

// TPC trajectory of one assembled track under one bunch-crossing hypothesis
// (Tpc_PolyTripletReco): three reference points + circle (xy) / line (rz).
//   reference points: TPC frame (as Tpc_PolyCluster centroids), inner / middle / outer
//   track parameters: world frame, at the point of closest approach to the beam (as Tpc_PolyTrack)
class Tpc_CrossingTrajectory : public PHObject
{
 public:
  static constexpr unsigned int NPoints = 3;

  Tpc_CrossingTrajectory() = default;
  ~Tpc_CrossingTrajectory() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_CrossingTrajectory base class" << std::endl;
  }
  int isValid() const override { return 0; }

  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_trajectory_id() const { return 0; }
  virtual unsigned int get_source_assembled_track_id() const { return 0; }
  virtual short get_crossing() const { return 0; }
  virtual int get_side() const { return 0; }
  virtual int get_fit_status() const { return 0; }

  // reference points, i = 0 (inner), 1 (middle), 2 (outer)
  virtual unsigned int get_point_layer(unsigned int) const { return 0; }
  virtual double get_point_x(unsigned int) const { return 0.0; }
  virtual double get_point_y(unsigned int) const { return 0.0; }
  virtual double get_point_z(unsigned int) const { return 0.0; }

  // track parameters
  virtual double get_x() const { return 0.0; }
  virtual double get_y() const { return 0.0; }
  virtual double get_z() const { return 0.0; }
  virtual double get_px() const { return 0.0; }
  virtual double get_py() const { return 0.0; }
  virtual double get_pz() const { return 0.0; }
  virtual double get_charge() const { return 0.0; }
  virtual double get_q_over_r() const { return 0.0; }  // charge * |curvature| [1/cm]
  virtual double get_phi() const { return 0.0; }       // direction at the pca
  virtual double get_slope() const { return 0.0; }     // dz/dr
  virtual double get_chi2() const { return 0.0; }      // rz line [cm^2]

  virtual void set_event(unsigned int) {}
  virtual void set_trajectory_id(unsigned int) {}
  virtual void set_source_assembled_track_id(unsigned int) {}
  virtual void set_crossing(short) {}
  virtual void set_side(int) {}
  virtual void set_fit_status(int) {}
  virtual void set_point(unsigned int /*i*/, unsigned int /*layer*/, double /*x*/, double /*y*/, double /*z*/) {}
  virtual void set_x(double) {}
  virtual void set_y(double) {}
  virtual void set_z(double) {}
  virtual void set_px(double) {}
  virtual void set_py(double) {}
  virtual void set_pz(double) {}
  virtual void set_charge(double) {}
  virtual void set_q_over_r(double) {}
  virtual void set_phi(double) {}
  virtual void set_slope(double) {}
  virtual void set_chi2(double) {}

 private:
  ClassDefOverride(Tpc_CrossingTrajectory, 0)
};
#endif
