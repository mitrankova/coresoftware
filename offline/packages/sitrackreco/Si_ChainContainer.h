// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SICHAINCONTAINER_H
#define SITRACKRECO_SICHAINCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>

class Si_Chain;

class Si_ChainContainer : public PHObject
{
 public:
  Si_ChainContainer() = default;
  ~Si_ChainContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Si_ChainContainer base class" << std::endl;
  }
  int isValid() const override { return 0; }

  virtual unsigned int size() const { return 0; }
  virtual void add_chain(Si_Chain*) {}
  virtual const Si_Chain* get_chain(unsigned int) const { return nullptr; }
  virtual Si_Chain* get_chain(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(Si_ChainContainer, 0)
};

#endif
