#include "TpcSiliconCrossingResolver.h"
#include "TpcSiliconMatchCandidate.h"
#include "TpcSiliconMatchCandidateContainer.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
#include <chrono>

TpcSiliconCrossingResolver::TpcSiliconCrossingResolver(const std::string& name) : SubsysReco(name) {}
int TpcSiliconCrossingResolver::InitRun(PHCompositeNode* topNode)
{
  m_candidates = findNode::getClass<TpcSiliconMatchCandidateContainer>(topNode, m_candidateNodeName);
  m_decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  if (!m_candidates || !m_decisions)
  {
    std::cerr << Name() << "::InitRun - missing " << m_candidateNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcSiliconCrossingResolver::process_event(PHCompositeNode*)
{
  const auto begin = std::chrono::steady_clock::now();
  std::vector<TpcSiliconMatchCandidate*> ordered;
  std::set<unsigned int> candidateParents;
  for (unsigned int i = 0; i < m_candidates->size(); ++i)
  {
    auto* candidate = m_candidates->get(i);
    if (candidate && candidate->isValid()) { candidate->set_selected(false); ordered.push_back(candidate); candidateParents.insert(candidate->get_source_assembled_track_id()); }
  }
  for (const auto id : candidateParents)
  {
    if (auto* decision = m_decisions->get_decision(id)) decision->set_status(TpcCrossingStatus::NoSiliconMatch);
  }
  std::sort(ordered.begin(), ordered.end(), [](const auto* lhs, const auto* rhs)
  {
    if (lhs->get_n_mvtx() != rhs->get_n_mvtx()) return lhs->get_n_mvtx() > rhs->get_n_mvtx();
    const auto lhsTotal = lhs->get_n_mvtx() + lhs->get_n_intt();
    const auto rhsTotal = rhs->get_n_mvtx() + rhs->get_n_intt();
    if (lhsTotal != rhsTotal) return lhsTotal > rhsTotal;
    if (lhs->get_n_intt() != rhs->get_n_intt()) return lhs->get_n_intt() > rhs->get_n_intt();
    return lhs->get_score() < rhs->get_score();
  });
  std::set<unsigned int> usedTracks;
  std::set<TrkrDefs::cluskey> usedClusters;
  for (auto* candidate : ordered)
  {
    if (usedTracks.count(candidate->get_parent_track_id())) continue;
    bool conflict = false;
    for (const auto key : candidate->get_silicon_cluster_keys()) if (usedClusters.count(key)) { conflict = true; break; }
    if (conflict) continue;
    candidate->set_selected(true);
    if (auto* decision = m_decisions->get_decision(candidate->get_source_assembled_track_id()))
    {
      decision->set_selected_crossing(candidate->get_crossing());
      decision->set_selected_score(candidate->get_score());
      decision->set_status(TpcCrossingStatus::ResolvedBySilicon);
    }
    usedTracks.insert(candidate->get_parent_track_id());
    usedClusters.insert(candidate->get_silicon_cluster_keys().begin(), candidate->get_silicon_cluster_keys().end());
  }
  if (Verbosity() > 0) std::cout << Name() << " resolved_tracks=" << usedTracks.size() << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
