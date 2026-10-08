#include "Si_ChainContainerv1.h"

#include "Si_Chain.h"

Si_ChainContainerv1::~Si_ChainContainerv1()
{
  Reset();
}

void Si_ChainContainerv1::identify(std::ostream& os) const
{
  os << "Si_ChainContainerv1 with "
     << m_chains.size()
     << " silicon chains"
     << std::endl;
}

void Si_ChainContainerv1::Reset()
{
  for (auto* chain : m_chains)
  {
    delete chain;
  }
  m_chains.clear();
}

int Si_ChainContainerv1::isValid() const
{
  return m_chains.empty() ? 0 : 1;
}

PHObject* Si_ChainContainerv1::CloneMe() const
{
  auto* copy = new Si_ChainContainerv1();
  for (const auto* chain : m_chains)
  {
    copy->m_chains.push_back(static_cast<Si_Chain*>(chain->CloneMe()));
  }
  return copy;
}
