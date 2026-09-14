#ifndef TPCTRACKRECO_TPCCROSSINGTRAJECTORYCONTAINER_H
#define TPCTRACKRECO_TPCCROSSINGTRAJECTORYCONTAINER_H

#include "TpcCrossingTrajectory.h"

#include <phool/PHObject.h>

#include <iostream>
#include <vector>

class TpcCrossingTrajectoryContainer : public PHObject
{
 public:
  ~TpcCrossingTrajectoryContainer() override { Reset(); }
  void identify(std::ostream& os = std::cout) const override { os << "TpcCrossingTrajectoryContainer size=" << m_values.size() << std::endl; }
  void Reset() override { for (auto* value : m_values) delete value; m_values.clear(); }
  int isValid() const override { return !m_values.empty(); }
  PHObject* CloneMe() const override;
  unsigned int size() const { return m_values.size(); }
  TpcCrossingTrajectory* get(unsigned int i) { return i < m_values.size() ? m_values[i] : nullptr; }
  const TpcCrossingTrajectory* get(unsigned int i) const { return i < m_values.size() ? m_values[i] : nullptr; }
  void add(TpcCrossingTrajectory* value) { if (value) m_values.push_back(value); }

 private:
  std::vector<TpcCrossingTrajectory*> m_values;
  ClassDefOverride(TpcCrossingTrajectoryContainer, 1)
};

#endif
