// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SICHAIN_H
#define SITRACKRECO_SICHAIN_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <utility>
#include <vector>

// Persistent silicon hit chain. Payload is implemented by versioned classes
// so the on-disk schema can evolve without changing this interface.
class Si_Chain : public PHObject
{
 public:
  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;

  Si_Chain() = default;
  ~Si_Chain() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Si_Chain base class" << std::endl;
  }
  int isValid() const override { return 0; }

  virtual unsigned int get_chain_id() const { return 0; }
  virtual unsigned int get_n_mvtx_hits() const { return 0; }
  virtual unsigned int get_n_intt_hits() const { return 0; }

  // Straight-line parameters in cylindrical coordinates: vertex is the z
  // intercept at r=0, slope is dz/dr, and phi is in radians.
  virtual double get_vertex() const { return 0.0; }
  virtual double get_slope() const { return 0.0; }
  virtual double get_phi() const { return 0.0; }

  virtual void set_chain_id(unsigned int) {}
  virtual void set_n_mvtx_hits(unsigned int) {}
  virtual void set_n_intt_hits(unsigned int) {}
  virtual void set_vertex(double) {}
  virtual void set_slope(double) {}
  virtual void set_phi(double) {}

  virtual void add_hit_index(TrkrDefs::hitsetkey, TrkrDefs::hitkey) {}
  virtual unsigned int size_hit_indices() const { return 0; }
  virtual HitIndex get_hit_index(unsigned int) const { return {0, 0}; }
  virtual const std::vector<HitIndex>& get_hit_indices() const
  {
    static const std::vector<HitIndex> empty_indices;
    return empty_indices;
  }

 private:
  ClassDefOverride(Si_Chain, 0)
};

#endif
