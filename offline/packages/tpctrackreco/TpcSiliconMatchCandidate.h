#ifndef TPCTRACKRECO_TPCSILICONMATCHCANDIDATE_H
#define TPCTRACKRECO_TPCSILICONMATCHCANDIDATE_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <limits>
#include <vector>

class TpcSiliconMatchCandidate : public PHObject
{
 public:
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new TpcSiliconMatchCandidate(*this); }

  unsigned int get_parent_track_id() const { return m_parentTrackId; }
  void set_parent_track_id(unsigned int value) { m_parentTrackId = value; }
  unsigned int get_source_assembled_track_id() const { return m_sourceAssembledTrackId; }
  void set_source_assembled_track_id(unsigned int value) { m_sourceAssembledTrackId = value; }
  short get_crossing() const { return m_crossing; }
  void set_crossing(short value) { m_crossing = value; }
  float get_score() const { return m_score; }
  void set_score(float value) { m_score = value; }
  float get_max_abs_dz() const { return m_maxAbsDz; }
  void set_max_abs_dz(float value) { m_maxAbsDz = value; }
  float get_max_abs_ddphi() const { return m_maxAbsDdphi; }
  void set_max_abs_ddphi(float value) { m_maxAbsDdphi = value; }
  unsigned int get_n_mvtx() const { return m_nMvtx; }
  void set_n_mvtx(unsigned int value) { m_nMvtx = value; }
  unsigned int get_n_intt() const { return m_nIntt; }
  void set_n_intt(unsigned int value) { m_nIntt = value; }
  bool get_selected() const { return m_selected; }
  void set_selected(bool value) { m_selected = value; }
  const std::vector<TrkrDefs::cluskey>& get_silicon_cluster_keys() const { return m_siliconKeys; }
  void add_silicon_cluster_key(TrkrDefs::cluskey key) { m_siliconKeys.push_back(key); }

 private:
  unsigned int m_parentTrackId{0};
  unsigned int m_sourceAssembledTrackId{0};
  short m_crossing{0};
  float m_score{std::numeric_limits<float>::quiet_NaN()};
  float m_maxAbsDz{std::numeric_limits<float>::quiet_NaN()};
  float m_maxAbsDdphi{std::numeric_limits<float>::quiet_NaN()};
  unsigned int m_nMvtx{0};
  unsigned int m_nIntt{0};
  bool m_selected{false};
  std::vector<TrkrDefs::cluskey> m_siliconKeys;
  ClassDefOverride(TpcSiliconMatchCandidate, 2)
};

#endif
