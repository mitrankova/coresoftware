#ifndef TPCTRACKRECO_FULLPOLYTRACK_H
#define TPCTRACKRECO_FULLPOLYTRACK_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <limits>
#include <vector>

class Full_PolyTrack : public PHObject
{
 public:
  ~Full_PolyTrack() override = default;
  void identify(std::ostream& os = std::cout) const override { os << "Full_PolyTrack base class" << std::endl; }
  int isValid() const override { return 0; }

  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }
  virtual unsigned int get_parent_tpc_track_id() const { return 0; }
  virtual unsigned int get_source_assembled_track_id() const { return 0; }
  virtual short get_crossing() const { return 0; }
  virtual unsigned char get_status() const { return 0; }
  virtual unsigned int get_n_mvtx() const { return 0; }
  virtual unsigned int get_n_intt() const { return 0; }
  virtual double get_score() const { return nan(); }
  virtual double get_max_abs_dz() const { return nan(); }
  virtual double get_max_abs_ddphi() const { return nan(); }
  virtual int get_fit_status() const { return 0; }
  virtual double get_chi2() const { return nan(); }
  virtual double get_ndf() const { return nan(); }
  virtual double get_x() const { return nan(); }
  virtual double get_y() const { return nan(); }
  virtual double get_z() const { return nan(); }
  virtual double get_px() const { return nan(); }
  virtual double get_py() const { return nan(); }
  virtual double get_pz() const { return nan(); }
  virtual double get_charge() const { return nan(); }
  virtual double get_cov(unsigned int, unsigned int) const { return nan(); }
  virtual bool has_final_native_state() const { return false; }
  virtual double get_final_native_state(unsigned int) const { return nan(); }
  virtual bool has_fast_native_state() const { return false; }
  virtual double get_fast_native_state(unsigned int) const { return nan(); }

  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}
  virtual void set_parent_tpc_track_id(unsigned int) {}
  virtual void set_source_assembled_track_id(unsigned int) {}
  virtual void set_crossing(short) {}
  virtual void set_status(unsigned char) {}
  virtual void set_n_mvtx(unsigned int) {}
  virtual void set_n_intt(unsigned int) {}
  virtual void set_score(double) {}
  virtual void set_max_abs_dz(double) {}
  virtual void set_max_abs_ddphi(double) {}
  virtual void set_fit_status(int) {}
  virtual void set_chi2(double) {}
  virtual void set_ndf(double) {}
  virtual void set_x(double) {}
  virtual void set_y(double) {}
  virtual void set_z(double) {}
  virtual void set_px(double) {}
  virtual void set_py(double) {}
  virtual void set_pz(double) {}
  virtual void set_charge(double) {}
  virtual void set_cov(unsigned int, unsigned int, double) {}
  virtual void set_final_native_state(unsigned int, double) {}
  virtual void set_fast_native_state(unsigned int, double) {}

  virtual unsigned int size_tpc_cluster_keys() const { return 0; }
  virtual TrkrDefs::cluskey get_tpc_cluster_key(unsigned int) const { return TrkrDefs::CLUSKEYMAX; }
  virtual const std::vector<TrkrDefs::cluskey>& get_tpc_cluster_keys() const { return empty_keys(); }
  virtual void add_tpc_cluster_key(TrkrDefs::cluskey) {}
  virtual void clear_tpc_cluster_keys() {}
  virtual unsigned int size_silicon_cluster_keys() const { return 0; }
  virtual TrkrDefs::cluskey get_silicon_cluster_key(unsigned int) const { return TrkrDefs::CLUSKEYMAX; }
  virtual const std::vector<TrkrDefs::cluskey>& get_silicon_cluster_keys() const { return empty_keys(); }
  virtual void add_silicon_cluster_key(TrkrDefs::cluskey) {}
  virtual void clear_silicon_cluster_keys() {}

 protected:
  static double nan() { return std::numeric_limits<double>::quiet_NaN(); }
  static const std::vector<TrkrDefs::cluskey>& empty_keys()
  {
    static const std::vector<TrkrDefs::cluskey> keys;
    return keys;
  }

 private:
  ClassDefOverride(Full_PolyTrack, 1)
};

#endif
