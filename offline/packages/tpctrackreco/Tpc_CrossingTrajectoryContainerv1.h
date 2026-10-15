// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCCROSSINGTRAJECTORYCONTAINERV1_H
#define TPCTRACKRECO_TPCCROSSINGTRAJECTORYCONTAINERV1_H

#include "Tpc_CrossingTrajectoryContainer.h"

#include <iostream>
#include <vector>

class Tpc_CrossingTrajectory;

class Tpc_CrossingTrajectoryContainerv1 : public Tpc_CrossingTrajectoryContainer
{
 public:
  Tpc_CrossingTrajectoryContainerv1() = default;
  ~Tpc_CrossingTrajectoryContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return m_trajectories.empty() ? 0 : 1; }
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_trajectories.size()); }
  void add_trajectory(Tpc_CrossingTrajectory* t) override { m_trajectories.push_back(t); }
  const Tpc_CrossingTrajectory* get_trajectory(unsigned int i) const override
  {
    return i < m_trajectories.size() ? m_trajectories[i] : nullptr;
  }
  Tpc_CrossingTrajectory* get_trajectory(unsigned int i) override
  {
    return i < m_trajectories.size() ? m_trajectories[i] : nullptr;
  }

 private:
  std::vector<Tpc_CrossingTrajectory*> m_trajectories;

  ClassDefOverride(Tpc_CrossingTrajectoryContainerv1, 1)
};
#endif
