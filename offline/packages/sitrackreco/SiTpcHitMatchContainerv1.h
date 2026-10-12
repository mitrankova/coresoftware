// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCHITMATCHCONTAINERV1_H
#define SITRACKRECO_SITPCHITMATCHCONTAINERV1_H

#include "SiTpcHitMatchContainer.h"

#include <iostream>
#include <vector>

class SiTpcHitMatch;

class SiTpcHitMatchContainerv1 : public SiTpcHitMatchContainer
{
 public:
  SiTpcHitMatchContainerv1() = default;
  ~SiTpcHitMatchContainerv1() override;

  // owns raw pointers: no implicit copies, use CloneMe()
  SiTpcHitMatchContainerv1(const SiTpcHitMatchContainerv1&) = delete;
  SiTpcHitMatchContainerv1& operator=(const SiTpcHitMatchContainerv1&) = delete;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_matches.size()); }
  void add_match(SiTpcHitMatch* m) override { m_matches.push_back(m); }
  const SiTpcHitMatch* get_match(unsigned int i) const override { return i < m_matches.size() ? m_matches[i] : nullptr; }
  SiTpcHitMatch* get_match(unsigned int i) override { return i < m_matches.size() ? m_matches[i] : nullptr; }

 private:
  std::vector<SiTpcHitMatch*> m_matches;

  ClassDefOverride(SiTpcHitMatchContainerv1, 1)
};

#endif
