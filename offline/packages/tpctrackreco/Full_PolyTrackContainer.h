#ifndef TPCTRACKRECO_FULLPOLYTRACKCONTAINER_H
#define TPCTRACKRECO_FULLPOLYTRACKCONTAINER_H
#include <phool/PHObject.h>
#include <iostream>
class Full_PolyTrack;
class Full_PolyTrackContainer : public PHObject
{
 public:
  ~Full_PolyTrackContainer() override = default;
  void identify(std::ostream& os = std::cout) const override { os << "Full_PolyTrackContainer base class" << std::endl; }
  int isValid() const override { return 0; }
  virtual unsigned int size() const { return 0; }
  virtual void add_track(Full_PolyTrack*) {}
  virtual Full_PolyTrack* get_track(unsigned int) { return nullptr; }
  virtual const Full_PolyTrack* get_track(unsigned int) const { return nullptr; }
 private:
  ClassDefOverride(Full_PolyTrackContainer, 1)
};
#endif
