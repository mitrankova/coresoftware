#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrack.h"

Full_PolyTrackContainerv1::Full_PolyTrackContainerv1()
{
  Reset();
}

Full_PolyTrackContainerv1::~Full_PolyTrackContainerv1()
{
  Reset();
}

void Full_PolyTrackContainerv1::identify(std::ostream& os) const
{
  os << "Full_PolyTrackContainerv1 with " << m_tracks.size() << " tracks" << std::endl;
}

void Full_PolyTrackContainerv1::Reset()
{
  for (auto*& track : m_tracks)
  {
    delete track;
  }
  m_tracks.clear();
}

int Full_PolyTrackContainerv1::isValid() const
{
  return m_tracks.empty() ? 0 : 1;
}

PHObject* Full_PolyTrackContainerv1::CloneMe() const
{
  Full_PolyTrackContainerv1* copy = new Full_PolyTrackContainerv1();
  for (auto* track : m_tracks)
  {
    if (track)
    {
      copy->m_tracks.push_back(static_cast<Full_PolyTrack*>(track->CloneMe()));
    }
  }
  return copy;
}
