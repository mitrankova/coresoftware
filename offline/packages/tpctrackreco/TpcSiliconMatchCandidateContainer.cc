#include "TpcSiliconMatchCandidateContainer.h"
ClassImp(TpcSiliconMatchCandidateContainer)
PHObject* TpcSiliconMatchCandidateContainer::CloneMe() const
{
  auto* output = new TpcSiliconMatchCandidateContainer;
  for (const auto* value : m_values) if (value) output->add(static_cast<TpcSiliconMatchCandidate*>(value->CloneMe()));
  return output;
}
