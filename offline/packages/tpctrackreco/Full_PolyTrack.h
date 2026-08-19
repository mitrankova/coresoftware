// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_FULLPOLYTRACK_H
#define TPCTRACKRECO_FULLPOLYTRACK_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <vector>

class Full_PolyTrack : public PHObject
{
 public:
  Full_PolyTrack() = default;
  ~Full_PolyTrack() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Full_PolyTrack base class" << std::endl;
  }
  int isValid() const override { return 0; }

  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }
  virtual unsigned int get_tpc_poly_track_id() const { return 0; }
  virtual unsigned int get_source_assembled_track_id() const { return 0; }
  virtual unsigned int get_n_tpc_clusters() const { return 0; }
  virtual unsigned int get_n_silicon_clusters() const { return 0; }
  virtual unsigned int get_n_missing_layers() const { return 0; }
  virtual unsigned int get_missing_layer_mask() const { return 0; }
  virtual int get_fit_status() const { return 0; }
  virtual double get_chi2() const { return 0.0; }
  virtual double get_ndf() const { return 0.0; }
  virtual double get_score() const { return 0.0; }
  virtual double get_x() const { return 0.0; }
  virtual double get_y() const { return 0.0; }
  virtual double get_z() const { return 0.0; }
  virtual double get_px() const { return 0.0; }
  virtual double get_py() const { return 0.0; }
  virtual double get_pz() const { return 0.0; }
  virtual double get_charge() const { return 0.0; }
  virtual double get_seed_x() const { return 0.0; }
  virtual double get_seed_y() const { return 0.0; }
  virtual double get_seed_z() const { return 0.0; }
  virtual double get_seed_px() const { return 0.0; }
  virtual double get_seed_py() const { return 0.0; }
  virtual double get_seed_pz() const { return 0.0; }

  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}
  virtual void set_tpc_poly_track_id(unsigned int) {}
  virtual void set_source_assembled_track_id(unsigned int) {}
  virtual void set_n_tpc_clusters(unsigned int) {}
  virtual void set_n_missing_layers(unsigned int) {}
  virtual void set_missing_layer_mask(unsigned int) {}
  virtual void set_fit_status(int) {}
  virtual void set_chi2(double) {}
  virtual void set_ndf(double) {}
  virtual void set_score(double) {}
  virtual void set_x(double) {}
  virtual void set_y(double) {}
  virtual void set_z(double) {}
  virtual void set_px(double) {}
  virtual void set_py(double) {}
  virtual void set_pz(double) {}
  virtual void set_charge(double) {}
  virtual void set_seed_x(double) {}
  virtual void set_seed_y(double) {}
  virtual void set_seed_z(double) {}
  virtual void set_seed_px(double) {}
  virtual void set_seed_py(double) {}
  virtual void set_seed_pz(double) {}

  virtual unsigned int size_cluster_keys() const { return 0; }
  virtual TrkrDefs::cluskey get_cluster_key(unsigned int) const { return TrkrDefs::CLUSKEYMAX; }
  virtual const std::vector<TrkrDefs::cluskey>& get_cluster_keys() const
  {
    static const std::vector<TrkrDefs::cluskey> empty_keys;
    return empty_keys;
  }
  virtual void add_cluster_key(TrkrDefs::cluskey) {}
  virtual void clear_cluster_keys() {}

  virtual unsigned int size_silicon_states() const { return 0; }
  virtual unsigned int get_state_layer(unsigned int) const { return 0; }
  virtual TrkrDefs::cluskey get_state_cluster_key(unsigned int) const { return TrkrDefs::CLUSKEYMAX; }
  virtual float get_state_x(unsigned int) const { return 0.0F; }
  virtual float get_state_y(unsigned int) const { return 0.0F; }
  virtual float get_state_z(unsigned int) const { return 0.0F; }
  virtual float get_state_pred_x(unsigned int) const { return 0.0F; }
  virtual float get_state_pred_y(unsigned int) const { return 0.0F; }
  virtual float get_state_pred_z(unsigned int) const { return 0.0F; }
  virtual float get_state_rdphi(unsigned int) const { return 0.0F; }
  virtual float get_state_dz(unsigned int) const { return 0.0F; }
  virtual float get_state_chi2(unsigned int) const { return 0.0F; }
  virtual void add_silicon_state(unsigned int, TrkrDefs::cluskey,
                                 float, float, float,
                                 float, float, float,
                                 float, float, float) {}
  virtual void clear_silicon_states() {}

 private:
  ClassDefOverride(Full_PolyTrack, 0)
};
#endif
