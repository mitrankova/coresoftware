// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCTRACKCONTAINERV1_H
#define SITRACKRECO_SITPCTRACKCONTAINERV1_H

#include "SiTpc_TrackContainer.h"

#include <iostream>
#include <vector>

class SiTpc_Track;

class SiTpc_TrackContainerv1 : public SiTpc_TrackContainer
{
 public:
  SiTpc_TrackContainerv1() = default;
  ~SiTpc_TrackContainerv1() override;

  // owns raw pointers: no implicit copies, use CloneMe()
  SiTpc_TrackContainerv1(const SiTpc_TrackContainerv1&) = delete;
  SiTpc_TrackContainerv1& operator=(const SiTpc_TrackContainerv1&) = delete;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override;

  unsigned int size() const override { return static_cast<unsigned int>(m_tracks.size()); }
  void add_track(SiTpc_Track* track) override { m_tracks.push_back(track); }
  const SiTpc_Track* get_track(unsigned int i) const override { return i < m_tracks.size() ? m_tracks[i] : nullptr; }
  SiTpc_Track* get_track(unsigned int i) override { return i < m_tracks.size() ? m_tracks[i] : nullptr; }

  void add_candidate(unsigned int si, unsigned int tpc, float dphi, float deta, float dz0, float chi2, int accepted) override
  {
    m_cand_si.push_back(si);
    m_cand_tpc.push_back(tpc);
    m_cand_dphi.push_back(dphi);
    m_cand_deta.push_back(deta);
    m_cand_dz0.push_back(dz0);
    m_cand_chi2.push_back(chi2);
    m_cand_accepted.push_back(accepted);
  }
  unsigned int size_candidates() const override { return static_cast<unsigned int>(m_cand_si.size()); }
  unsigned int get_candidate_si(unsigned int i) const override { return i < m_cand_si.size() ? m_cand_si[i] : 0; }
  unsigned int get_candidate_tpc(unsigned int i) const override { return i < m_cand_tpc.size() ? m_cand_tpc[i] : 0; }
  float get_candidate_dphi(unsigned int i) const override { return i < m_cand_dphi.size() ? m_cand_dphi[i] : NAN; }
  float get_candidate_deta(unsigned int i) const override { return i < m_cand_deta.size() ? m_cand_deta[i] : NAN; }
  float get_candidate_dz0(unsigned int i) const override { return i < m_cand_dz0.size() ? m_cand_dz0[i] : NAN; }
  float get_candidate_chi2(unsigned int i) const override { return i < m_cand_chi2.size() ? m_cand_chi2[i] : NAN; }
  int get_candidate_accepted(unsigned int i) const override { return i < m_cand_accepted.size() ? m_cand_accepted[i] : 0; }
  void set_candidate_accepted(unsigned int i, int v) override
  {
    if (i < m_cand_accepted.size())
    {
      m_cand_accepted[i] = v;
    }
  }

 private:
  std::vector<SiTpc_Track*> m_tracks;
  std::vector<unsigned int> m_cand_si;
  std::vector<unsigned int> m_cand_tpc;
  std::vector<float> m_cand_dphi;
  std::vector<float> m_cand_deta;
  std::vector<float> m_cand_dz0;
  std::vector<float> m_cand_chi2;
  std::vector<int> m_cand_accepted;

  ClassDefOverride(SiTpc_TrackContainerv1, 1)
};

#endif
