#include "TpcSiliconCrossingResolver.h"
#include "TpcSiliconMatchCandidate.h"
#include "TpcSiliconMatchCandidateContainer.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

TpcSiliconCrossingResolver::TpcSiliconCrossingResolver(const std::string& name) : SubsysReco(name) {}
int TpcSiliconCrossingResolver::InitRun(PHCompositeNode* topNode)
{
  m_candidates = findNode::getClass<TpcSiliconMatchCandidateContainer>(topNode, m_candidateNodeName);
  if (!m_candidates)
  {
    std::cerr << Name() << "::InitRun - missing " << m_candidateNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcSiliconCrossingResolver::process_event(PHCompositeNode*)
{
  std::vector<TpcSiliconMatchCandidate*> ordered;
  for (unsigned int i = 0; i < m_candidates->size(); ++i)
  {
    auto* candidate = m_candidates->get(i);
    if (candidate && candidate->isValid()) { candidate->set_selected(false); ordered.push_back(candidate); }
  }
  std::sort(ordered.begin(), ordered.end(), [](const auto* lhs, const auto* rhs) { return lhs->get_score() < rhs->get_score(); });
  std::set<unsigned int> usedTracks;
  std::set<TrkrDefs::cluskey> usedClusters;
  for (auto* candidate : ordered)
  {
    if (usedTracks.count(candidate->get_parent_track_id())) continue;
    bool conflict = false;
    for (const auto key : candidate->get_silicon_cluster_keys()) if (usedClusters.count(key)) { conflict = true; break; }
    if (conflict) continue;
    candidate->set_selected(true);
    usedTracks.insert(candidate->get_parent_track_id());
    usedClusters.insert(candidate->get_silicon_cluster_keys().begin(), candidate->get_silicon_cluster_keys().end());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
