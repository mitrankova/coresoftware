// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCTRACKCONTAINER_H
#define SITRACKRECO_SITPCTRACKCONTAINER_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>

class SiTpc_Track;

// Matched Si + TPC tracks, plus the list of all Si-TPC candidate pairs inside the loose
// candidate window of SiTpcTrackMatcher (for tuning the match windows).
class SiTpc_TrackContainer : public PHObject
{
 public:
  SiTpc_TrackContainer() = default;
  ~SiTpc_TrackContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SiTpc_TrackContainer base class" << std::endl;
  }
  int isValid() const override { return 0; }

  // ---- tracks (the container takes ownership)
  virtual unsigned int size() const { return 0; }
  virtual void add_track(SiTpc_Track*) {}
  virtual const SiTpc_Track* get_track(unsigned int) const { return nullptr; }
  virtual SiTpc_Track* get_track(unsigned int) { return nullptr; }

  // ---- candidate pairs
  // si/tpc: indices in SI_TRAJECTORY / TPC_POLYTRACKS; accepted: 1 if this pair became a track
  virtual void add_candidate(unsigned int /*si*/, unsigned int /*tpc*/, float /*dphi*/, float /*deta*/,
                             float /*dz0*/, float /*chi2*/, int /*accepted*/) {}
  virtual unsigned int size_candidates() const { return 0; }
  virtual unsigned int get_candidate_si(unsigned int) const { return 0; }
  virtual unsigned int get_candidate_tpc(unsigned int) const { return 0; }
  virtual float get_candidate_dphi(unsigned int) const { return NAN; }
  virtual float get_candidate_deta(unsigned int) const { return NAN; }
  virtual float get_candidate_dz0(unsigned int) const { return NAN; }
  virtual float get_candidate_chi2(unsigned int) const { return NAN; }
  virtual int get_candidate_accepted(unsigned int) const { return 0; }
  virtual void set_candidate_accepted(unsigned int, int) {}

 private:
  ClassDefOverride(SiTpc_TrackContainer, 0)
};

#endif
