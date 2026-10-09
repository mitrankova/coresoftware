// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITRAJECTORY_H
#define SITRACKRECO_SITRAJECTORY_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

// Persistent silicon trajectory: one SiHitSeedReco chain, its hits converted to the
// detector frame, shifted to the global (aligned) frame, and fitted with a circle in xy
// plus a straight line in (s, z).  The payload is implemented by versioned classes
// (Si_Trajectoryv1, ...) so the on-disk schema can evolve without changing this interface.
//
// Units: cm, rad, GeV.  All positions are in the GLOBAL frame (after the detector shift).
//
// Conventions
//   helicity  +1: the trajectory turns counter-clockwise (seen from +z) while moving outward,
//             -1: clockwise.
//   charge    = -helicity * sign(Bz)   (set by the fitter from its configured field).
//   pca       point of closest approach of the circle to the beam position (xy).
//   phi       direction of motion at the pca.
//   dca       signed distance beam -> pca, positive if the beam is to the right of the
//             direction of motion (sign of (pca - beam) x t, t = unit tangent).
//   z0, tanl  z(s) = z0 + tanl * s, s = transverse arc length from the pca (s > 0 outward).
class Si_Trajectory : public PHObject
{
 public:
  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  enum FitStatus
  {
    NotFitted = 0,
    Ok = 1,
    TooFewPoints = 2,
    Degenerate = 3,
    // xy fitted with a straight line (lever arm too short for curvature): the circle is a
    // very large radius, pt = NaN, charge = 0.  All other parameters are valid.
    StraightLine = 4
  };

  bool is_fitted() const { return get_fit_status() == Ok || get_fit_status() == StraightLine; }

  Si_Trajectory() = default;
  ~Si_Trajectory() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Si_Trajectory base class" << std::endl;
  }
  int isValid() const override { return 0; }

  // ---- bookkeeping
  virtual unsigned int get_id() const { return 0; }
  virtual int get_chain_id() const { return -1; }
  virtual unsigned int get_n_mvtx() const { return 0; }
  virtual unsigned int get_n_intt() const { return 0; }
  virtual int get_fit_status() const { return NotFitted; }

  virtual void set_id(unsigned int) {}
  virtual void set_chain_id(int) {}
  virtual void set_n_mvtx(unsigned int) {}
  virtual void set_n_intt(unsigned int) {}
  virtual void set_fit_status(int) {}

  // ---- fit points: one per layer, global frame [cm]
  virtual void add_point(int /*layer*/, float /*x*/, float /*y*/, float /*z*/, unsigned int /*nhits*/) {}
  virtual unsigned int size_points() const { return 0; }
  virtual int get_point_layer(unsigned int) const { return -1; }
  virtual float get_point_x(unsigned int) const { return NAN; }
  virtual float get_point_y(unsigned int) const { return NAN; }
  virtual float get_point_z(unsigned int) const { return NAN; }
  virtual unsigned int get_point_nhits(unsigned int) const { return 0; }

  // ---- circle in xy
  virtual double get_circle_x() const { return NAN; }
  virtual double get_circle_y() const { return NAN; }
  virtual double get_radius() const { return NAN; }
  virtual int get_helicity() const { return 0; }
  virtual double get_circle_rms() const { return NAN; }  // rms of (|p - c| - R) [cm]

  virtual void set_circle_x(double) {}
  virtual void set_circle_y(double) {}
  virtual void set_radius(double) {}
  virtual void set_helicity(int) {}
  virtual void set_circle_rms(double) {}

  // ---- line in (s, z)
  virtual double get_z0() const { return NAN; }
  virtual double get_tanl() const { return NAN; }
  virtual double get_z_rms() const { return NAN; }

  virtual void set_z0(double) {}
  virtual void set_tanl(double) {}
  virtual void set_z_rms(double) {}

  // ---- perigee w.r.t. the beam position
  virtual double get_pca_x() const { return NAN; }
  virtual double get_pca_y() const { return NAN; }
  virtual double get_phi() const { return NAN; }
  virtual double get_dca() const { return NAN; }
  virtual double get_pt() const { return NAN; }
  virtual int get_charge() const { return 0; }

  virtual void set_pca_x(double) {}
  virtual void set_pca_y(double) {}
  virtual void set_phi(double) {}
  virtual void set_dca(double) {}
  virtual void set_pt(double) {}
  virtual void set_charge(int) {}

  // derived, no storage needed
  double get_eta() const { return std::asinh(get_tanl()); }
  double get_theta() const { return std::atan2(1.0, get_tanl()); }
  double get_p() const { return get_pt() * std::sqrt(1.0 + get_tanl() * get_tanl()); }

  // ---- raw hits of the chain
  virtual void add_hit_index(TrkrDefs::hitsetkey, TrkrDefs::hitkey) {}
  virtual unsigned int size_hit_indices() const { return 0; }
  virtual HitIndex get_hit_index(unsigned int) const { return {0, 0}; }
  virtual const std::vector<HitIndex>& get_hit_indices() const
  {
    static const std::vector<HitIndex> empty_indices;
    return empty_indices;
  }

 private:
  ClassDefOverride(Si_Trajectory, 0)
};

#endif
