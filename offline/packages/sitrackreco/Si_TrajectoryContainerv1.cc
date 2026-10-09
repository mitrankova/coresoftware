#include "Si_TrajectoryContainerv1.h"

#include "Si_Trajectory.h"

Si_TrajectoryContainerv1::~Si_TrajectoryContainerv1()
{
  Reset();
}

void Si_TrajectoryContainerv1::identify(std::ostream& os) const
{
  os << "Si_TrajectoryContainerv1 with "
     << m_trajectories.size()
     << " silicon trajectories"
     << std::endl;
}

void Si_TrajectoryContainerv1::Reset()
{
  for (auto* trajectory : m_trajectories)
  {
    delete trajectory;
  }
  m_trajectories.clear();
}

int Si_TrajectoryContainerv1::isValid() const
{
  return m_trajectories.empty() ? 0 : 1;
}

PHObject* Si_TrajectoryContainerv1::CloneMe() const
{
  auto* copy = new Si_TrajectoryContainerv1();
  for (const auto* trajectory : m_trajectories)
  {
    copy->m_trajectories.push_back(static_cast<Si_Trajectory*>(trajectory->CloneMe()));
  }
  return copy;
}
