// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_FULLPOLYTRACKCONTAINERV1_H
#define TPCTRACKRECO_FULLPOLYTRACKCONTAINERV1_H

#include "Full_PolyTrackContainer.h"

#include <iostream>
#include <vector>

class Full_PolyTrack;

class Full_PolyTrackContainerv1 : public Full_PolyTrackContainer
{
 public:
  Full_PolyTrackContainerv1();
  ~Full_PolyTrackContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_tracks.size()); }
  void add_track(Full_PolyTrack* trk) override { m_tracks.push_back(trk); }
  const Full_PolyTrack* get_track(unsigned int i) const override
  {
    return i < m_tracks.size() ? m_tracks[i] : nullptr;
  }
  Full_PolyTrack* get_track(unsigned int i) override
  {
    return i < m_tracks.size() ? m_tracks[i] : nullptr;
  }

 private:
  std::vector<Full_PolyTrack*> m_tracks;

  ClassDefOverride(Full_PolyTrackContainerv1, 1)
};
#endif
