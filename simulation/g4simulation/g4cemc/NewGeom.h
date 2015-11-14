#ifndef NEWGEOM_H_
#define NEWGEOM_H_

#include "RawTowerDefs.h"

#include <phool/phool.h>
#include <phool/PHObject.h>
#include <iostream>
#include <map>

class NewGeom : public PHObject {

 public:

  virtual ~NewGeom() {}

  virtual void identify(std::ostream& os=std::cout) const { PHOOL_VIRTUAL_WARN("identify()"); }

  virtual RawTowerDefs::keytype get_id() const { PHOOL_VIRTUAL_WARN("get_id()"); return 0; }

  virtual void set_center_x( double ) { PHOOL_VIRTUAL_WARN("set_center_x()"); return ; }
  virtual void set_center_y( double ) { PHOOL_VIRTUAL_WARN("set_center_y()"); return ; }
  virtual void set_center_z( double ) { PHOOL_VIRTUAL_WARN("set_center_z()"); return ; }

  virtual void set_size_x( double ) { PHOOL_VIRTUAL_WARN("set_size_x()"); return ; }
  virtual void set_size_y( double ) { PHOOL_VIRTUAL_WARN("set_size_y()"); return ; }
  virtual void set_size_z( double ) { PHOOL_VIRTUAL_WARN("set_size_z()"); return ; }

  virtual double get_center_x() const { PHOOL_VIRTUAL_WARN("get_center_x()"); return -1; }
  virtual double get_center_y() const { PHOOL_VIRTUAL_WARN("get_center_y()"); return -1; }
  virtual double get_center_z() const { PHOOL_VIRTUAL_WARN("get_center_z()"); return -1; }

  virtual double get_size_x() const { PHOOL_VIRTUAL_WARN("get_size_x()"); return -1; }
  virtual double get_size_y() const { PHOOL_VIRTUAL_WARN("get_size_y()"); return -1; }
  virtual double get_size_z() const { PHOOL_VIRTUAL_WARN("get_size_z()"); return -1; }
  virtual double get_volume() const { PHOOL_VIRTUAL_WARN("get_volume()"); return -1; }

  virtual double get_center_radius() const { PHOOL_VIRTUAL_WARN("get_center_radius()"); return -1; }
  virtual double get_eta() const { PHOOL_VIRTUAL_WARN("get_eta()"); return -1; }
  virtual double get_phi() const { PHOOL_VIRTUAL_WARN("get_phi()"); return -1; }

 protected:
  NewGeom() {}

  ClassDef(NewGeom,1)

};

#endif /* NEWGEOM_H_ */
