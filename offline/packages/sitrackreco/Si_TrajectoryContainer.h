// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITRAJECTORYCONTAINER_H
#define SITRACKRECO_SITRAJECTORYCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>

class Si_Trajectory;

class Si_TrajectoryContainer : public PHObject
{
 public:
  Si_TrajectoryContainer() = default;
  ~Si_TrajectoryContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Si_TrajectoryContainer base class" << std::endl;
  }
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }
  // The container takes ownership of the trajectory.
  virtual void add_trajectory(Si_Trajectory*) {}
  virtual const Si_Trajectory* get_trajectory(unsigned int) const { return nullptr; }
  virtual Si_Trajectory* get_trajectory(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(Si_TrajectoryContainer, 0)
};

#endif
