#ifndef TPCTRACKRECO_FULLPOLYTRACKCONTAINERV1_H
#define TPCTRACKRECO_FULLPOLYTRACKCONTAINERV1_H
#include "Full_PolyTrackContainer.h"
#include <vector>
class Full_PolyTrack;
class Full_PolyTrackContainerv1 : public Full_PolyTrackContainer
{
 public:
  ~Full_PolyTrackContainerv1() override { Reset(); }
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return !m_tracks.empty(); }
  PHObject* CloneMe() const override;
  unsigned int size() const override { return m_tracks.size(); }
  void add_track(Full_PolyTrack* track) override { if (track) m_tracks.push_back(track); }
  Full_PolyTrack* get_track(unsigned int i) override { return i < m_tracks.size() ? m_tracks[i] : nullptr; }
  const Full_PolyTrack* get_track(unsigned int i) const override { return i < m_tracks.size() ? m_tracks[i] : nullptr; }
 private:
  std::vector<Full_PolyTrack*> m_tracks;
  ClassDefOverride(Full_PolyTrackContainerv1, 1)
};
#endif
