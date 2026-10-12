// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCHITMATCHCONTAINER_H
#define SITRACKRECO_SITPCHITMATCHCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>

class SiTpcHitMatch;

// TPC-driven silicon hit chains of one event (SiTpcHitChainMatcher), node SI_TPC_HITMATCHES.
class SiTpcHitMatchContainer : public PHObject
{
 public:
  SiTpcHitMatchContainer() = default;
  ~SiTpcHitMatchContainer() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SiTpcHitMatchContainer base class" << std::endl;
  }
  int isValid() const override { return 0; }

  // the container takes ownership
  virtual unsigned int size() const { return 0; }
  virtual void add_match(SiTpcHitMatch*) {}
  virtual const SiTpcHitMatch* get_match(unsigned int) const { return nullptr; }
  virtual SiTpcHitMatch* get_match(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(SiTpcHitMatchContainer, 0)
};

#endif
