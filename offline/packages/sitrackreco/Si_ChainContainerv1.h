// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SICHAINCONTAINERV1_H
#define SITRACKRECO_SICHAINCONTAINERV1_H

#include "Si_ChainContainer.h"

#include <iostream>
#include <vector>

class Si_Chain;

class Si_ChainContainerv1 : public Si_ChainContainer
{
 public:
  Si_ChainContainerv1() = default;
  ~Si_ChainContainerv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override
  {
    return static_cast<unsigned int>(m_chains.size());
  }

  void add_chain(Si_Chain* chain) override
  {
    m_chains.push_back(chain);
  }

  const Si_Chain* get_chain(unsigned int index) const override
  {
    if (index >= m_chains.size()) return nullptr;
    return m_chains[index];
  }

  Si_Chain* get_chain(unsigned int index) override
  {
    if (index >= m_chains.size()) return nullptr;
    return m_chains[index];
  }

 private:
  std::vector<Si_Chain*> m_chains;

  ClassDefOverride(Si_ChainContainerv1, 1)
};

#endif
