// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCHITMATCH_H
#define SITRACKRECO_SITPCHITMATCH_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

// One attached hit group ("step") of a TPC-driven silicon hit chain.  Plain struct used by the
// interface; SiTpcHitMatchv1 stores the fields column-wise.
//
// Residual conventions (all angles in rad, intrinsic azimuth u, polar angles from the TPC vertex)
//   dphi_tpc, dtheta_tpc : group centroid - TPC helix at the centroid radius
//   dz_tpc               : z of the centroid - TPC helix z [cm]
//   dphi, dtheta         : group centroid - silicon chain prediction (chain built from the
//                          previous groups with TPC curvature and TPC theta); for the entry group
//                          they equal dphi_tpc / dtheta_tpc
//   ddphi, ddtheta       : dphi - (offset + slope * dphi_previous), same for theta
//                          (SiTpcHitChainMatcher::setDeltaDeltaPhiMean / ThetaMean)
struct SiTpcHitMatchStep
{
  enum Flag
  {
    Entry = 1,      // first group of the chain: tested against the TPC helix
    Overlap = 2,    // second measurement in the same layer (MVTX stave / INTT ladder overlap)
    Truncated = 4,  // neighbour growth stopped by the expected-footprint limit
    CrossChip = 8   // group spans two abutting MVTX chips
  };

  int layer = -1;
  int n_hits = 0;
  int flags = 0;
  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::HITSETKEYMAX;  // sensor of the first hit
  float u = NAN, z = NAN, r = NAN;
  float dphi_tpc = NAN, dtheta_tpc = NAN, dz_tpc = NAN;
  float dphi = NAN, dtheta = NAN;
  float ddphi = NAN, ddtheta = NAN;
  float chi2 = NAN;
  float extent_rphi = NAN, extent_z = NAN;          // group size [pitches]
  float exp_extent_rphi = NAN, exp_extent_z = NAN;  // expected from the TPC incidence [pitches]
};

// TPC-driven silicon hit chain for one TPC track (SiTpcHitChainMatcher).
// Persistent; payload in versioned classes (SiTpcHitMatchv1, ...).
class SiTpcHitMatch : public PHObject
{
 public:
  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  SiTpcHitMatch() = default;
  ~SiTpcHitMatch() override = default;

  void identify(std::ostream& os = std::cout) const override { os << "SiTpcHitMatch base class" << std::endl; }
  int isValid() const override { return 0; }

  // ---- bookkeeping
  virtual unsigned int get_id() const { return 0; }
  virtual unsigned int get_tpc_track_index() const { return 0; }        // index in TPC_POLYTRACKS
  virtual unsigned int get_tpc_track_id() const { return 0; }           // Tpc_PolyTrack::get_track_id
  virtual unsigned int get_tpc_assembled_track_id() const { return 0; }
  virtual short get_crossing() const { return 0; }  // 0 in AuAu
  virtual unsigned int get_rank() const { return 0; }  // 0 = best chain of this TPC track
  virtual bool get_selected() const { return false; }  // set by SiTpcHitChainResolver

  virtual void set_id(unsigned int) {}
  virtual void set_tpc_track_index(unsigned int) {}
  virtual void set_tpc_track_id(unsigned int) {}
  virtual void set_tpc_assembled_track_id(unsigned int) {}
  virtual void set_crossing(short) {}
  virtual void set_rank(unsigned int) {}
  virtual void set_selected(bool) {}

  // ---- TPC reference (TPC-only helix, TPC frame)
  virtual float get_tpc_pt() const { return NAN; }
  virtual int get_tpc_charge() const { return 0; }
  virtual float get_tpc_phi() const { return NAN; }
  virtual float get_tpc_tanl() const { return NAN; }
  virtual float get_tpc_z0() const { return NAN; }
  virtual float get_tpc_dca() const { return NAN; }

  virtual void set_tpc_pt(float) {}
  virtual void set_tpc_charge(int) {}
  virtual void set_tpc_phi(float) {}
  virtual void set_tpc_tanl(float) {}
  virtual void set_tpc_z0(float) {}
  virtual void set_tpc_dca(float) {}

  // ---- chain summary
  virtual unsigned int get_n_mvtx_layers() const { return 0; }
  virtual unsigned int get_n_intt_layers() const { return 0; }
  virtual unsigned int get_n_missing() const { return 0; }
  virtual float get_chi2() const { return NAN; }
  // silicon chain (TPC radius and tanl) in the detector frame
  virtual float get_si_z0() const { return NAN; }   // z at the pca to the beam
  virtual float get_si_dca() const { return NAN; }  // |beam - centre| - R
  // silicon chain circle, TPC frame: radius and helicity are the TPC ones, the centre is anchored
  // on the silicon groups; pca = closest point to the beam.  With get_si_z0 and get_tpc_tanl it
  // gives the full chain helix (e.g. for event displays).
  virtual float get_si_circle_x() const { return NAN; }
  virtual float get_si_circle_y() const { return NAN; }
  virtual float get_si_radius() const { return NAN; }
  virtual int get_si_helicity() const { return 0; }
  virtual float get_si_pca_x() const { return NAN; }
  virtual float get_si_pca_y() const { return NAN; }
  // TPC-Si compatibility, for tuning
  virtual float get_delta_phi_dir() const { return NAN; }  // direction si - tpc at the outermost chain layer
  virtual float get_delta_z0() const { return NAN; }       // si_z0 - TPC z0 (both at the beam)
  virtual float get_delta_dca() const { return NAN; }      // si_dca - TPC dca (|beam - centre| - R convention)

  virtual void set_n_mvtx_layers(unsigned int) {}
  virtual void set_n_intt_layers(unsigned int) {}
  virtual void set_n_missing(unsigned int) {}
  virtual void set_chi2(float) {}
  virtual void set_si_z0(float) {}
  virtual void set_si_dca(float) {}
  virtual void set_si_circle_x(float) {}
  virtual void set_si_circle_y(float) {}
  virtual void set_si_radius(float) {}
  virtual void set_si_helicity(int) {}
  virtual void set_si_pca_x(float) {}
  virtual void set_si_pca_y(float) {}
  virtual void set_delta_phi_dir(float) {}
  virtual void set_delta_z0(float) {}
  virtual void set_delta_dca(float) {}

  // ---- attached groups
  virtual void add_step(const SiTpcHitMatchStep&) {}
  virtual unsigned int size_steps() const { return 0; }
  virtual SiTpcHitMatchStep get_step(unsigned int) const { return {}; }

  // ---- hits (every hit of every group); step index of each hit
  virtual void add_hit_index(TrkrDefs::hitsetkey, TrkrDefs::hitkey, unsigned int /*step*/) {}
  virtual unsigned int size_hit_indices() const { return 0; }
  virtual HitIndex get_hit_index(unsigned int) const { return {0, 0}; }
  virtual unsigned int get_hit_step(unsigned int) const { return 0; }
  virtual const std::vector<HitIndex>& get_hit_indices() const
  {
    static const std::vector<HitIndex> empty;
    return empty;
  }

 private:
  ClassDefOverride(SiTpcHitMatch, 0)
};

#endif