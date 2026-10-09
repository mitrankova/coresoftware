// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITRAJECTORYCONTAINERV1_H
#define SITRACKRECO_SITRAJECTORYCONTAINERV1_H

#include "Si_TrajectoryContainer.h"

#include <iostream>
#include <vector>

class Si_Trajectory;

class Si_TrajectoryContainerv1 : public Si_TrajectoryContainer
{
 public:
  Si_TrajectoryContainerv1() = default;
  ~Si_TrajectoryContainerv1() override;

  // owns raw pointers: no implicit copies, use CloneMe()
  Si_TrajectoryContainerv1(const Si_TrajectoryContainerv1&) = delete;
  Si_TrajectoryContainerv1& operator=(const Si_TrajectoryContainerv1&) = delete;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override
  {
    return static_cast<unsigned int>(m_trajectories.size());
  }

  void add_trajectory(Si_Trajectory* trajectory) override
  {
    m_trajectories.push_back(trajectory);
  }

  const Si_Trajectory* get_trajectory(unsigned int index) const override
  {
    if (index >= m_trajectories.size()) return nullptr;
    return m_trajectories[index];
  }

  Si_Trajectory* get_trajectory(unsigned int index) override
  {
    if (index >= m_trajectories.size()) return nullptr;
    return m_trajectories[index];
  }

 private:
  std::vector<Si_Trajectory*> m_trajectories;

  ClassDefOverride(Si_TrajectoryContainerv1, 1)
};

#endif
