// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCTRACK_H
#define SITRACKRECO_SITPCTRACK_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>

// Matched silicon + TPC track (SiTpcTrackMatcher): one Si_Trajectory matched in phi / eta to
// one TPC poly track, refitted with the Si layer points and the TPC clusters together.
// The payload is implemented by versioned classes (SiTpc_Trackv1, ...).
//
// Frame: beam-axis frame of SiTpcBeamAlignment (beam at (0, 0)).  Units cm, rad, GeV.
// Helix conventions as Si_Trajectory / SiTpcHelixFit (helicity, charge, signed dca, z(s)).
// TPC z: the refit can have a free z offset for the TPC points (bunch crossing / t0 not
// known): TPC point z = z0 + tanl * s + tpc_z_offset.
class SiTpc_Track : public PHObject
{
 public:
  // point types
  enum PointType
  {
    SiPoint = 0,
    TpcPoint = 1
  };

  SiTpc_Track() = default;
  ~SiTpc_Track() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SiTpc_Track base class" << std::endl;
  }
  int isValid() const override { return 0; }

  // ---- bookkeeping / association
  virtual unsigned int get_id() const { return 0; }
  virtual unsigned int get_si_trajectory_index() const { return 0; }  // index in SI_TRAJECTORY
  virtual int get_si_chain_id() const { return -1; }
  virtual unsigned int get_tpc_track_index() const { return 0; }       // index in TPC_POLYTRACKS
  virtual unsigned int get_tpc_track_id() const { return 0; }          // Tpc_PolyTrack::get_track_id
  virtual unsigned int get_tpc_assembled_track_id() const { return 0; }
  virtual unsigned int get_n_si_points() const { return 0; }
  virtual unsigned int get_n_tpc_points() const { return 0; }
  virtual int get_fit_status() const { return 0; }  // SiTpcHelixFit::Status

  virtual void set_id(unsigned int) {}
  virtual void set_si_trajectory_index(unsigned int) {}
  virtual void set_si_chain_id(int) {}
  virtual void set_tpc_track_index(unsigned int) {}
  virtual void set_tpc_track_id(unsigned int) {}
  virtual void set_tpc_assembled_track_id(unsigned int) {}
  virtual void set_n_si_points(unsigned int) {}
  virtual void set_n_tpc_points(unsigned int) {}
  virtual void set_fit_status(int) {}

  // ---- matching: Si and TPC-only parameters at their pca, and the differences
  virtual double get_si_phi() const { return NAN; }
  virtual double get_si_eta() const { return NAN; }
  virtual double get_si_z0() const { return NAN; }
  virtual double get_tpc_phi() const { return NAN; }
  virtual double get_tpc_eta() const { return NAN; }
  virtual double get_tpc_z0() const { return NAN; }
  virtual double get_tpc_pt() const { return NAN; }
  virtual int get_tpc_charge() const { return 0; }
  virtual double get_match_dphi() const { return NAN; }  // si - tpc, wrapped
  virtual double get_match_deta() const { return NAN; }  // si - tpc
  virtual double get_match_dz0() const { return NAN; }   // si - tpc
  virtual double get_match_chi2() const { return NAN; }

  virtual void set_si_phi(double) {}
  virtual void set_si_eta(double) {}
  virtual void set_si_z0(double) {}
  virtual void set_tpc_phi(double) {}
  virtual void set_tpc_eta(double) {}
  virtual void set_tpc_z0(double) {}
  virtual void set_tpc_pt(double) {}
  virtual void set_tpc_charge(int) {}
  virtual void set_match_dphi(double) {}
  virtual void set_match_deta(double) {}
  virtual void set_match_dz0(double) {}
  virtual void set_match_chi2(double) {}

  // ---- combined fit
  virtual double get_circle_x() const { return NAN; }
  virtual double get_circle_y() const { return NAN; }
  virtual double get_radius() const { return NAN; }
  virtual int get_helicity() const { return 0; }
  virtual double get_circle_rms() const { return NAN; }
  virtual double get_pca_x() const { return NAN; }
  virtual double get_pca_y() const { return NAN; }
  virtual double get_phi() const { return NAN; }
  virtual double get_dca() const { return NAN; }
  virtual double get_pt() const { return NAN; }
  virtual int get_charge() const { return 0; }
  virtual double get_z0() const { return NAN; }
  virtual double get_tanl() const { return NAN; }
  virtual double get_tpc_z_offset() const { return 0; }
  virtual double get_z_rms() const { return NAN; }

  virtual void set_circle_x(double) {}
  virtual void set_circle_y(double) {}
  virtual void set_radius(double) {}
  virtual void set_helicity(int) {}
  virtual void set_circle_rms(double) {}
  virtual void set_pca_x(double) {}
  virtual void set_pca_y(double) {}
  virtual void set_phi(double) {}
  virtual void set_dca(double) {}
  virtual void set_pt(double) {}
  virtual void set_charge(int) {}
  virtual void set_z0(double) {}
  virtual void set_tanl(double) {}
  virtual void set_tpc_z_offset(double) {}
  virtual void set_z_rms(double) {}

  double get_eta() const { return std::asinh(get_tanl()); }
  double get_theta() const { return std::atan2(1.0, get_tanl()); }
  double get_p() const { return get_pt() * std::sqrt(1.0 + get_tanl() * get_tanl()); }

  // ---- fit points (beam-axis frame, TPC z as measured, i.e. without the offset)
  //   type SiPoint : layer 0-6, source = Si_Trajectory point index
  //   type TpcPoint: layer = TPC layer, source = index in TPC_POLYCLUSTERS
  virtual void add_point(int /*type*/, int /*layer*/, unsigned int /*source*/, float /*x*/, float /*y*/, float /*z*/) {}
  virtual unsigned int size_points() const { return 0; }
  virtual int get_point_type(unsigned int) const { return -1; }
  virtual int get_point_layer(unsigned int) const { return -1; }
  virtual unsigned int get_point_source(unsigned int) const { return 0; }
  virtual float get_point_x(unsigned int) const { return NAN; }
  virtual float get_point_y(unsigned int) const { return NAN; }
  virtual float get_point_z(unsigned int) const { return NAN; }

 private:
  ClassDefOverride(SiTpc_Track, 0)
};

#endif
