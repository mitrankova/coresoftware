#include "SiTpcHitMatchContainerv1.h"

#include "SiTpcHitMatch.h"

SiTpcHitMatchContainerv1::~SiTpcHitMatchContainerv1()
{
  Reset();
}

void SiTpcHitMatchContainerv1::identify(std::ostream& os) const
{
  os << "SiTpcHitMatchContainerv1 with " << m_matches.size() << " TPC-driven silicon hit chains" << std::endl;
}

void SiTpcHitMatchContainerv1::Reset()
{
  for (auto* m : m_matches)
  {
    delete m;
  }
  m_matches.clear();
}

int SiTpcHitMatchContainerv1::isValid() const
{
  return m_matches.empty() ? 0 : 1;
}

PHObject* SiTpcHitMatchContainerv1::CloneMe() const
{
  auto* copy = new SiTpcHitMatchContainerv1();
  for (const auto* m : m_matches)
  {
    copy->m_matches.push_back(static_cast<SiTpcHitMatch*>(m->CloneMe()));
  }
  return copy;
}
