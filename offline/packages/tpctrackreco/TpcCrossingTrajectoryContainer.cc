#include "TpcCrossingTrajectoryContainer.h"

ClassImp(TpcCrossingTrajectoryContainer)

PHObject* TpcCrossingTrajectoryContainer::CloneMe() const
{
  auto* output = new TpcCrossingTrajectoryContainer;
  for (const auto* value : m_values)
  {
    if (value) output->add(static_cast<TpcCrossingTrajectory*>(value->CloneMe()));
  }
  return output;
}
