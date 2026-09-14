#ifndef TPCTRACKRECO_TPCSILICONMATCHCANDIDATECONTAINER_H
#define TPCTRACKRECO_TPCSILICONMATCHCANDIDATECONTAINER_H

#include "TpcSiliconMatchCandidate.h"
#include <phool/PHObject.h>
#include <iostream>
#include <vector>

class TpcSiliconMatchCandidateContainer : public PHObject
{
 public:
  ~TpcSiliconMatchCandidateContainer() override { Reset(); }
  void identify(std::ostream& os = std::cout) const override { os << "TpcSiliconMatchCandidateContainer size=" << m_values.size() << std::endl; }
  void Reset() override { for (auto* value : m_values) delete value; m_values.clear(); }
  int isValid() const override { return !m_values.empty(); }
  PHObject* CloneMe() const override;
  unsigned int size() const { return m_values.size(); }
  TpcSiliconMatchCandidate* get(unsigned int i) { return i < m_values.size() ? m_values[i] : nullptr; }
  const TpcSiliconMatchCandidate* get(unsigned int i) const { return i < m_values.size() ? m_values[i] : nullptr; }
  void add(TpcSiliconMatchCandidate* value) { if (value) m_values.push_back(value); }
 private:
  std::vector<TpcSiliconMatchCandidate*> m_values;
  ClassDefOverride(TpcSiliconMatchCandidateContainer, 1)
};

#endif
