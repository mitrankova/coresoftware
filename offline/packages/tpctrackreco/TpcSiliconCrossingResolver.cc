#include "TpcSiliconCrossingResolver.h"
#include "TpcSiliconMatchCandidate.h"
#include "TpcSiliconMatchCandidateContainer.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <set>
#include <vector>

TpcSiliconCrossingResolver::TpcSiliconCrossingResolver(const std::string& name)
  : SubsysReco(name)
{
}

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
  std::map<unsigned int, unsigned int> parentToSource;
  std::map<unsigned int, unsigned int> invalidByParent;
  std::map<unsigned int, unsigned int> conflictsByParent;
  for (unsigned int i = 0; i < m_candidates->size(); ++i)
  {
    auto* candidate = m_candidates->get(i);
    if (!candidate) continue;
    const unsigned int parent = candidate->get_parent_track_id();
    candidate->set_selected(false);
    candidateParents.insert(parent);
    parentToSource[parent] = candidate->get_source_assembled_track_id();
    if (!candidate->isValid())
    {
      ++invalidByParent[parent];
      continue;
    }
    ordered.push_back(candidate);
  }
  for (const auto parent : candidateParents)
  {
    const auto source = parentToSource.find(parent);
    if (source != parentToSource.end())
      if (auto* decision = m_decisions->get_decision(source->second))
        decision->set_status(TpcCrossingStatus::NoSiliconMatch);
  }

  std::sort(ordered.begin(), ordered.end(), [](const auto* lhs, const auto* rhs)
  {
    const auto lhsTotal = lhs->get_n_mvtx() + lhs->get_n_intt();
    const auto rhsTotal = rhs->get_n_mvtx() + rhs->get_n_intt();
    if (lhsTotal != rhsTotal) return lhsTotal > rhsTotal;
    if (lhs->get_direction_score() != rhs->get_direction_score())
      return lhs->get_direction_score() < rhs->get_direction_score();
    if (lhs->get_n_mvtx() != rhs->get_n_mvtx()) return lhs->get_n_mvtx() > rhs->get_n_mvtx();
    if (lhs->get_si_internal_score() != rhs->get_si_internal_score())
      return lhs->get_si_internal_score() < rhs->get_si_internal_score();
    if (lhs->get_crossing() != rhs->get_crossing()) return lhs->get_crossing() < rhs->get_crossing();
    return lhs->get_parent_track_id() < rhs->get_parent_track_id();
  });

  std::set<unsigned int> usedTracks;
  std::set<TrkrDefs::cluskey> usedClusters;
  for (auto* candidate : ordered)
  {
    const unsigned int parent = candidate->get_parent_track_id();
    if (usedTracks.count(parent)) continue;
    TrkrDefs::cluskey conflictingKey = TrkrDefs::CLUSKEYMAX;
    for (const auto key : candidate->get_silicon_cluster_keys())
    {
      if (usedClusters.count(key))
      {
        conflictingKey = key;
        break;
      }
    }
    if (conflictingKey != TrkrDefs::CLUSKEYMAX)
    {
      ++conflictsByParent[parent];
      if (Verbosity() >= 2)
      {
        std::cout << Name() << " skip_shared_cluster parent_track_id=" << parent
                  << " crossing=" << candidate->get_crossing()
                  << " cluster_key=" << conflictingKey
                  << " direction_score=" << candidate->get_direction_score() << std::endl;
      }
      continue;
    }

    candidate->set_selected(true);
    if (auto* decision = m_decisions->get_decision(candidate->get_source_assembled_track_id()))
    {
      decision->set_selected_crossing(candidate->get_crossing());
      decision->set_selected_score(candidate->get_direction_score());
      decision->set_status(TpcCrossingStatus::ResolvedBySilicon);
    }
    usedTracks.insert(parent);
    usedClusters.insert(candidate->get_silicon_cluster_keys().begin(),
                        candidate->get_silicon_cluster_keys().end());
  }

  if (Verbosity() >= 2)
  {
    for (const auto parent : candidateParents)
    {
      if (usedTracks.count(parent)) continue;
      std::cout << Name() << " unmatched parent_track_id=" << parent
                << " invalid_candidates=" << invalidByParent[parent]
                << " conflicting_candidates=" << conflictsByParent[parent]
                << " reason="
                << (conflictsByParent[parent] ? "no_nonconflicting_candidate"
                                             : "no_valid_candidate")
                << std::endl;
    }
  }
  if (Verbosity() > 0)
    std::cout << Name() << " resolved_tracks=" << usedTracks.size()
              << " seconds="
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count()
              << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
