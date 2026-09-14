#ifndef TPCTRACKRECO_FULLPOLYTRACKV1_H
#define TPCTRACKRECO_FULLPOLYTRACKV1_H

#include "Full_PolyTrack.h"

#include <array>

class Full_PolyTrackv1 : public Full_PolyTrack
{
 public:
  Full_PolyTrackv1() { Reset(); }
  ~Full_PolyTrackv1() override = default;
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new Full_PolyTrackv1(*this); }

#define FULL_GETSET(type, name, member) \
  type get_##name() const override { return member; } \
  void set_##name(type value) override { member = value; }
  FULL_GETSET(unsigned int, event, m_event)
  FULL_GETSET(unsigned int, track_id, m_trackId)
  FULL_GETSET(unsigned int, parent_tpc_track_id, m_parentTpcTrackId)
  FULL_GETSET(unsigned int, source_assembled_track_id, m_sourceAssembledTrackId)
  FULL_GETSET(short, crossing, m_crossing)
  FULL_GETSET(unsigned char, status, m_status)
  FULL_GETSET(unsigned int, n_mvtx, m_nMvtx)
  FULL_GETSET(unsigned int, n_intt, m_nIntt)
  FULL_GETSET(double, score, m_score)
  FULL_GETSET(double, max_abs_dz, m_maxAbsDz)
  FULL_GETSET(double, max_abs_ddphi, m_maxAbsDdphi)
  FULL_GETSET(int, fit_status, m_fitStatus)
  FULL_GETSET(double, chi2, m_chi2)
  FULL_GETSET(double, ndf, m_ndf)
  FULL_GETSET(double, x, m_x)
  FULL_GETSET(double, y, m_y)
  FULL_GETSET(double, z, m_z)
  FULL_GETSET(double, px, m_px)
  FULL_GETSET(double, py, m_py)
  FULL_GETSET(double, pz, m_pz)
  FULL_GETSET(double, charge, m_charge)
#undef FULL_GETSET

  double get_cov(unsigned int i, unsigned int j) const override;
  void set_cov(unsigned int i, unsigned int j, double value) override;
  unsigned int size_tpc_cluster_keys() const override { return m_tpcKeys.size(); }
  TrkrDefs::cluskey get_tpc_cluster_key(unsigned int i) const override { return i < m_tpcKeys.size() ? m_tpcKeys[i] : TrkrDefs::CLUSKEYMAX; }
  const std::vector<TrkrDefs::cluskey>& get_tpc_cluster_keys() const override { return m_tpcKeys; }
  void add_tpc_cluster_key(TrkrDefs::cluskey key) override { m_tpcKeys.push_back(key); }
  void clear_tpc_cluster_keys() override { m_tpcKeys.clear(); }
  unsigned int size_silicon_cluster_keys() const override { return m_siliconKeys.size(); }
  TrkrDefs::cluskey get_silicon_cluster_key(unsigned int i) const override { return i < m_siliconKeys.size() ? m_siliconKeys[i] : TrkrDefs::CLUSKEYMAX; }
  const std::vector<TrkrDefs::cluskey>& get_silicon_cluster_keys() const override { return m_siliconKeys; }
  void add_silicon_cluster_key(TrkrDefs::cluskey key) override { m_siliconKeys.push_back(key); }
  void clear_silicon_cluster_keys() override { m_siliconKeys.clear(); }

 private:
  unsigned int m_event{0}, m_trackId{0}, m_parentTpcTrackId{0}, m_sourceAssembledTrackId{0};
  short m_crossing{0};
  unsigned char m_status{0};
  unsigned int m_nMvtx{0}, m_nIntt{0};
  double m_score{nan()}, m_maxAbsDz{nan()}, m_maxAbsDdphi{nan()};
  int m_fitStatus{0};
  double m_chi2{nan()}, m_ndf{nan()};
  double m_x{nan()}, m_y{nan()}, m_z{nan()}, m_px{nan()}, m_py{nan()}, m_pz{nan()}, m_charge{nan()};
  std::array<double, 36> m_cov{};
  std::vector<TrkrDefs::cluskey> m_tpcKeys;
  std::vector<TrkrDefs::cluskey> m_siliconKeys;
  ClassDefOverride(Full_PolyTrackv1, 1)
};

#endif
