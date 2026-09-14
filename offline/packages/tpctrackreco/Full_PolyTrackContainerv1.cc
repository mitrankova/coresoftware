#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrack.h"
void Full_PolyTrackContainerv1::identify(std::ostream& os) const { os << "Full_PolyTrackContainerv1 size=" << m_tracks.size() << std::endl; }
void Full_PolyTrackContainerv1::Reset() { for (auto* track : m_tracks) delete track; m_tracks.clear(); }
PHObject* Full_PolyTrackContainerv1::CloneMe() const
{
  auto* copy = new Full_PolyTrackContainerv1;
  for (const auto* track : m_tracks) if (track) copy->add_track(static_cast<Full_PolyTrack*>(track->CloneMe()));
  return copy;
}
