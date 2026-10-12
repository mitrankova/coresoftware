// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCHITCHAINQA_H
#define SITRACKRECO_SITPCHITCHAINQA_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TTree;

// SiTpcHitChainQA
//
// Tuning and cross-check trees for SiTpcHitChainMatcher (+ SiTpcHitChainResolver).
//
//   step  : one row per attached group: layer, group size vs expected footprint, residuals to the
//           TPC helix (dphi_tpc, dtheta_tpc, dz_tpc) and to the silicon chain (dphi, dtheta,
//           ddphi, ddtheta), flags (entry / overlap / truncated / cross-chip), plus TPC pt, tanl.
//           -> tune setEntryWindow / setDeltaDeltaWindow / setDeltaDeltaPhiMean / ThetaMean
//              and the group slack.
//   match : one row per chain: summary, TPC-Si compatibility, and the HIT-LEVEL cross-check
//           against the independent silicon chaining (SI_HIT_SEED_EVENT chains):
//             xc_chain     SiHitSeedReco chain sharing most hits (-1: none)
//             xc_shared    number of shared hits, xc_chain_hits its size, xc_jaccard
//             xc_ntouched  number of different SiHitSeedReco chains touched
//           and against SiTpcTrackMatcher (SITPC_TRACKS, if present):
//             st_chain     SiHitSeedReco chain matched by SiTpcTrackMatcher to the same TPC track
//             st_shared    hits shared with that chain, st_agree = (st_chain == xc_chain)
//   event : counts per event.
class SiTpcHitChainQA : public SubsysReco
{
 public:
  explicit SiTpcHitChainQA(const std::string& name = "SiTpcHitChainQA");
  ~SiTpcHitChainQA() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setOutputFile(const std::string& f) { m_outputFile = f; }
  void setMatchNodeName(const std::string& s) { m_matchNodeName = s; }
  void setHitEventNodeName(const std::string& s) { m_hitEventNodeName = s; }
  void setSiTpcTrackNodeName(const std::string& s) { m_siTpcTrackNodeName = s; }

 private:
  std::string m_outputFile = "SiTpcHitChainQA.root";
  std::string m_matchNodeName = "SI_TPC_HITMATCHES";
  std::string m_hitEventNodeName = "SI_HIT_SEED_EVENT";
  std::string m_siTpcTrackNodeName = "SITPC_TRACKS";

  TFile* m_file = nullptr;
  TTree* m_tEvent = nullptr;
  TTree* m_tMatch = nullptr;
  TTree* m_tStep = nullptr;
  unsigned long m_event = 0;

  // ---- branch buffers
  // event
  int b_nMatches = 0, b_nSelected = 0, b_nSiHits = 0, b_nSeedChains = 0, b_nSiTpc = 0;
  // match
  int b_id = 0, b_tpcIndex = 0, b_rank = 0, b_selected = 0, b_tpcCharge = 0;
  float b_tpcPt = 0, b_tpcPhi = 0, b_tpcTanl = 0, b_tpcZ0 = 0, b_tpcDca = 0;
  int b_nMvtx = 0, b_nIntt = 0, b_nMissing = 0, b_nSteps = 0, b_nHits = 0;
  float b_chi2 = 0, b_siZ0 = 0, b_siDca = 0, b_dPhiDir = 0, b_dZ0 = 0, b_dDca = 0;
  int b_xcChain = -1, b_xcShared = 0, b_xcChainHits = 0, b_xcTouched = 0;
  float b_xcJaccard = 0;
  int b_stChain = -1, b_stShared = 0, b_stAgree = 0;
  // step
  int s_index = 0, s_layer = 0, s_nHits = 0, s_flags = 0;
  float s_u = 0, s_z = 0, s_r = 0;
  float s_dphiTpc = 0, s_dthetaTpc = 0, s_dzTpc = 0;
  float s_dphi = 0, s_dtheta = 0, s_ddphi = 0, s_ddtheta = 0, s_chi2 = 0;
  float s_extRphi = 0, s_extZ = 0, s_expRphi = 0, s_expZ = 0;
};

#endif
