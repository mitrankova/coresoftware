#include "Tpc_CrossingTrajectoryContainerv1.h"
#include "Tpc_CrossingTrajectory.h"

Tpc_CrossingTrajectoryContainerv1::~Tpc_CrossingTrajectoryContainerv1()
{
  Reset();
}

void Tpc_CrossingTrajectoryContainerv1::identify(std::ostream& os) const
{
  os << "Tpc_CrossingTrajectoryContainerv1 with " << m_trajectories.size() << " trajectories" << std::endl;
}

void Tpc_CrossingTrajectoryContainerv1::Reset()
{
  for (auto* t : m_trajectories)
  {
    delete t;
  }
  m_trajectories.clear();
}

PHObject* Tpc_CrossingTrajectoryContainerv1::CloneMe() const
{
  auto* copy = new Tpc_CrossingTrajectoryContainerv1();
  for (const auto* t : m_trajectories)
  {
    if (t)
    {
      copy->m_trajectories.push_back(static_cast<Tpc_CrossingTrajectory*>(t->CloneMe()));
    }
  }
  return copy;
}
