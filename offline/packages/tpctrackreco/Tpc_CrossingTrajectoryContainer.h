// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCCROSSINGTRAJECTORYCONTAINER_H
#define TPCTRACKRECO_TPCCROSSINGTRAJECTORYCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>

class Tpc_CrossingTrajectory;

class Tpc_CrossingTrajectoryContainer : public PHObject
{
 public:
  Tpc_CrossingTrajectoryContainer() = default;
  ~Tpc_CrossingTrajectoryContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_CrossingTrajectoryContainer base class" << std::endl;
  }
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }
  virtual void add_trajectory(Tpc_CrossingTrajectory*) {}
  virtual const Tpc_CrossingTrajectory* get_trajectory(unsigned int) const { return nullptr; }
  virtual Tpc_CrossingTrajectory* get_trajectory(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(Tpc_CrossingTrajectoryContainer, 0)
};
#endif
