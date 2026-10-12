#include "SiTpcHitChainQA.h"

#include "SiHitSeedData.h"
#include "SiTpcHitMatch.h"
#include "SiTpcHitMatchContainer.h"
#include "SiTpc_Track.h"
#include "SiTpc_TrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

namespace
{
  uint64_t packHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk)
  {
    return (static_cast<uint64_t>(hsk) << 32) | static_cast<uint64_t>(hk);
  }
}  // namespace

SiTpcHitChainQA::SiTpcHitChainQA(const std::string& name)
  : SubsysReco(name)
{
}

int SiTpcHitChainQA::InitRun(PHCompositeNode*)
{
  if (m_file)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_file = new TFile(m_outputFile.c_str(), "RECREATE");

  m_tEvent = new TTree("event", "SiTpcHitChainMatcher event summary");
  m_tEvent->Branch("event", &m_event, "event/l");
  m_tEvent->Branch("n_matches", &b_nMatches, "n_matches/I");
  m_tEvent->Branch("n_selected", &b_nSelected, "n_selected/I");
  m_tEvent->Branch("n_si_hits", &b_nSiHits, "n_si_hits/I");
  m_tEvent->Branch("n_seed_chains", &b_nSeedChains, "n_seed_chains/I");
  m_tEvent->Branch("n_sitpc_tracks", &b_nSiTpc, "n_sitpc_tracks/I");

  m_tMatch = new TTree("match", "TPC-driven silicon hit chains");
  m_tMatch->Branch("event", &m_event, "event/l");
  m_tMatch->Branch("id", &b_id, "id/I");
  m_tMatch->Branch("tpc_index", &b_tpcIndex, "tpc_index/I");
  m_tMatch->Branch("rank", &b_rank, "rank/I");
  m_tMatch->Branch("selected", &b_selected, "selected/I");
  m_tMatch->Branch("tpc_pt", &b_tpcPt, "tpc_pt/F");
  m_tMatch->Branch("tpc_charge", &b_tpcCharge, "tpc_charge/I");
  m_tMatch->Branch("tpc_phi", &b_tpcPhi, "tpc_phi/F");
  m_tMatch->Branch("tpc_tanl", &b_tpcTanl, "tpc_tanl/F");
  m_tMatch->Branch("tpc_z0", &b_tpcZ0, "tpc_z0/F");
  m_tMatch->Branch("tpc_dca", &b_tpcDca, "tpc_dca/F");
  m_tMatch->Branch("n_mvtx_layers", &b_nMvtx, "n_mvtx_layers/I");
  m_tMatch->Branch("n_intt_layers", &b_nIntt, "n_intt_layers/I");
  m_tMatch->Branch("n_missing", &b_nMissing, "n_missing/I");
  m_tMatch->Branch("n_steps", &b_nSteps, "n_steps/I");
  m_tMatch->Branch("n_hits", &b_nHits, "n_hits/I");
  m_tMatch->Branch("chi2", &b_chi2, "chi2/F");
  m_tMatch->Branch("si_z0", &b_siZ0, "si_z0/F");
  m_tMatch->Branch("si_dca", &b_siDca, "si_dca/F");
  m_tMatch->Branch("delta_phi_dir", &b_dPhiDir, "delta_phi_dir/F");
  m_tMatch->Branch("delta_z0", &b_dZ0, "delta_z0/F");
  m_tMatch->Branch("delta_dca", &b_dDca, "delta_dca/F");
  m_tMatch->Branch("xc_chain", &b_xcChain, "xc_chain/I");
  m_tMatch->Branch("xc_shared", &b_xcShared, "xc_shared/I");
  m_tMatch->Branch("xc_chain_hits", &b_xcChainHits, "xc_chain_hits/I");
  m_tMatch->Branch("xc_ntouched", &b_xcTouched, "xc_ntouched/I");
  m_tMatch->Branch("xc_jaccard", &b_xcJaccard, "xc_jaccard/F");
  m_tMatch->Branch("st_chain", &b_stChain, "st_chain/I");
  m_tMatch->Branch("st_shared", &b_stShared, "st_shared/I");
  m_tMatch->Branch("st_agree", &b_stAgree, "st_agree/I");

  m_tStep = new TTree("step", "attached hit groups");
  m_tStep->Branch("event", &m_event, "event/l");
  m_tStep->Branch("id", &b_id, "id/I");
  m_tStep->Branch("rank", &b_rank, "rank/I");
  m_tStep->Branch("selected", &b_selected, "selected/I");
  m_tStep->Branch("tpc_pt", &b_tpcPt, "tpc_pt/F");
  m_tStep->Branch("tpc_tanl", &b_tpcTanl, "tpc_tanl/F");
  m_tStep->Branch("tpc_charge", &b_tpcCharge, "tpc_charge/I");
  m_tStep->Branch("n_layers", &b_nSteps, "n_layers/I");
  m_tStep->Branch("index", &s_index, "index/I");
  m_tStep->Branch("layer", &s_layer, "layer/I");
  m_tStep->Branch("n_hits", &s_nHits, "n_hits/I");
  m_tStep->Branch("flags", &s_flags, "flags/I");
  m_tStep->Branch("u", &s_u, "u/F");
  m_tStep->Branch("z", &s_z, "z/F");
  m_tStep->Branch("r", &s_r, "r/F");
  m_tStep->Branch("dphi_tpc", &s_dphiTpc, "dphi_tpc/F");
  m_tStep->Branch("dtheta_tpc", &s_dthetaTpc, "dtheta_tpc/F");
  m_tStep->Branch("dz_tpc", &s_dzTpc, "dz_tpc/F");
  m_tStep->Branch("dphi", &s_dphi, "dphi/F");
  m_tStep->Branch("dtheta", &s_dtheta, "dtheta/F");
  m_tStep->Branch("ddphi", &s_ddphi, "ddphi/F");
  m_tStep->Branch("ddtheta", &s_ddtheta, "ddtheta/F");
  m_tStep->Branch("chi2", &s_chi2, "chi2/F");
  m_tStep->Branch("extent_rphi", &s_extRphi, "extent_rphi/F");
  m_tStep->Branch("extent_z", &s_extZ, "extent_z/F");
  m_tStep->Branch("exp_extent_rphi", &s_expRphi, "exp_extent_rphi/F");
  m_tStep->Branch("exp_extent_z", &s_expZ, "exp_extent_z/F");
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpcHitChainQA::process_event(PHCompositeNode* topNode)
{
  auto* matches = findNode::getClass<SiTpcHitMatchContainer>(topNode, m_matchNodeName);
  if (!matches)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  auto* event = findNode::getClass<SiHitSeedEvent>(topNode, m_hitEventNodeName);
  auto* sitpc = findNode::getClass<SiTpc_TrackContainer>(topNode, m_siTpcTrackNodeName);

  // ---- independent chaining: hit -> SiHitSeedReco chain
  std::unordered_map<uint64_t, int> chainOfHit;
  std::vector<int> chainSize;
  if (event)
  {
    std::vector<int> chainOfIndex(event->hits.size(), -1);
    chainSize.assign(event->chains.size(), 0);
    for (std::size_t c = 0; c < event->chains.size(); ++c)
    {
      for (int h : event->chains[c].hit_ids)
      {
        if (h >= 0 && h < static_cast<int>(chainOfIndex.size()))
        {
          chainOfIndex[h] = static_cast<int>(c);
          ++chainSize[c];
        }
      }
    }
    chainOfHit.reserve(event->hits.size());
    for (std::size_t i = 0; i < event->hits.size(); ++i)
    {
      if (chainOfIndex[i] >= 0)
      {
        chainOfHit[packHit(event->hits[i].hitsetkey, event->hits[i].hitkey)] = chainOfIndex[i];
      }
    }
  }
  // ---- SiTpcTrackMatcher: TPC track index -> SiHitSeedReco chain id
  std::map<unsigned int, int> stChainOfTpc;
  if (sitpc)
  {
    for (unsigned int i = 0; i < sitpc->size(); ++i)
    {
      const SiTpc_Track* t = sitpc->get_track(i);
      if (t)
      {
        stChainOfTpc[t->get_tpc_track_index()] = t->get_si_chain_id();
      }
    }
  }

  b_nMatches = static_cast<int>(matches->size());
  b_nSelected = 0;
  for (unsigned int i = 0; i < matches->size(); ++i)
  {
    const SiTpcHitMatch* m = matches->get_match(i);
    if (!m)
    {
      continue;
    }
    b_id = static_cast<int>(m->get_id());
    b_tpcIndex = static_cast<int>(m->get_tpc_track_index());
    b_rank = static_cast<int>(m->get_rank());
    b_selected = m->get_selected() ? 1 : 0;
    b_nSelected += b_selected;
    b_tpcPt = m->get_tpc_pt();
    b_tpcCharge = m->get_tpc_charge();
    b_tpcPhi = m->get_tpc_phi();
    b_tpcTanl = m->get_tpc_tanl();
    b_tpcZ0 = m->get_tpc_z0();
    b_tpcDca = m->get_tpc_dca();
    b_nMvtx = static_cast<int>(m->get_n_mvtx_layers());
    b_nIntt = static_cast<int>(m->get_n_intt_layers());
    b_nMissing = static_cast<int>(m->get_n_missing());
    b_nSteps = static_cast<int>(m->size_steps());
    b_nHits = static_cast<int>(m->size_hit_indices());
    b_chi2 = m->get_chi2();
    b_siZ0 = m->get_si_z0();
    b_siDca = m->get_si_dca();
    b_dPhiDir = m->get_delta_phi_dir();
    b_dZ0 = m->get_delta_z0();
    b_dDca = m->get_delta_dca();

    // hit-level cross-check
    std::map<int, int> shared;
    for (const auto& h : m->get_hit_indices())
    {
      auto it = chainOfHit.find(packHit(h.first, h.second));
      if (it != chainOfHit.end())
      {
        ++shared[it->second];
      }
    }
    b_xcChain = -1;
    b_xcShared = 0;
    b_xcChainHits = 0;
    b_xcJaccard = 0;
    b_xcTouched = static_cast<int>(shared.size());
    for (const auto& [c, n] : shared)
    {
      if (n > b_xcShared)
      {
        b_xcShared = n;
        b_xcChain = c;
      }
    }
    if (b_xcChain >= 0)
    {
      b_xcChainHits = chainSize[b_xcChain];
      b_xcJaccard = static_cast<float>(b_xcShared) / static_cast<float>(b_nHits + b_xcChainHits - b_xcShared);
    }
    b_stChain = -1;
    b_stShared = 0;
    b_stAgree = 0;
    auto st = stChainOfTpc.find(m->get_tpc_track_index());
    if (st != stChainOfTpc.end())
    {
      b_stChain = st->second;
      auto sh = shared.find(b_stChain);
      b_stShared = sh != shared.end() ? sh->second : 0;
      b_stAgree = (b_stChain >= 0 && b_stChain == b_xcChain) ? 1 : 0;
    }
    m_tMatch->Fill();

    for (unsigned int k = 0; k < m->size_steps(); ++k)
    {
      const auto s = m->get_step(k);
      s_index = static_cast<int>(k);
      s_layer = s.layer;
      s_nHits = s.n_hits;
      s_flags = s.flags;
      s_u = s.u;
      s_z = s.z;
      s_r = s.r;
      s_dphiTpc = s.dphi_tpc;
      s_dthetaTpc = s.dtheta_tpc;
      s_dzTpc = s.dz_tpc;
      s_dphi = s.dphi;
      s_dtheta = s.dtheta;
      s_ddphi = s.ddphi;
      s_ddtheta = s.ddtheta;
      s_chi2 = s.chi2;
      s_extRphi = s.extent_rphi;
      s_extZ = s.extent_z;
      s_expRphi = s.exp_extent_rphi;
      s_expZ = s.exp_extent_z;
      m_tStep->Fill();
    }
  }
  b_nSiHits = event ? static_cast<int>(event->hits.size()) : -1;
  b_nSeedChains = event ? static_cast<int>(event->chains.size()) : -1;
  b_nSiTpc = sitpc ? static_cast<int>(sitpc->size()) : -1;
  m_tEvent->Fill();
  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpcHitChainQA::End(PHCompositeNode*)
{
  if (m_file)
  {
    m_file->cd();
    m_tEvent->Write();
    m_tMatch->Write();
    m_tStep->Write();
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
