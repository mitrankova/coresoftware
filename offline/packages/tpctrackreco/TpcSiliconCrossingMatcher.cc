#include "TpcSiliconCrossingMatcher.h"
#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
#include "TpcSiliconMatchCandidate.h"
#include "TpcSiliconMatchCandidateContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
  struct SiliconPoint { TrkrDefs::cluskey key; float z; float phi; };
  float angleDifference(float a, float b) { return std::remainder(a - b, 2.F * static_cast<float>(M_PI)); }
}

TpcSiliconCrossingMatcher::TpcSiliconCrossingMatcher(const std::string& name) : SubsysReco(name) {}

int TpcSiliconCrossingMatcher::getNodes(PHCompositeNode* topNode)
{
  m_trajectories = findNode::getClass<TpcCrossingTrajectoryContainer>(topNode, m_trajectoryNodeName);
  m_clusters = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterNodeName);
  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_trajectories || !m_clusters || !m_geometry)
  {
    std::cerr << Name() << "::getNodes - missing trajectory, cluster, or ActsGeometry input" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcSiliconCrossingMatcher::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dst = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst) { dst = new PHCompositeNode("DST"); topNode->addNode(dst); }
  m_candidates = findNode::getClass<TpcSiliconMatchCandidateContainer>(topNode, m_outputNodeName);
  if (!m_candidates)
  {
    m_candidates = new TpcSiliconMatchCandidateContainer;
    dst->addNode(new PHIODataNode<PHObject>(m_candidates, m_outputNodeName, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcSiliconCrossingMatcher::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  return createNodes(topNode);
}

int TpcSiliconCrossingMatcher::process_event(PHCompositeNode*)
{
  m_candidates->Reset();
  std::array<std::vector<SiliconPoint>, 7> byLayer;
  for (const auto detector : {TrkrDefs::TrkrId::mvtxId, TrkrDefs::TrkrId::inttId})
  {
    for (const auto hitsetkey : m_clusters->getHitSetKeys(detector))
    {
      const auto range = m_clusters->getClusters(hitsetkey);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        const auto key = iter->first;
        const unsigned int layer = TrkrDefs::getLayer(key);
        if (layer >= byLayer.size() || !iter->second) continue;
        const auto position = m_geometry->getGlobalPosition(key, iter->second);
        byLayer[layer].push_back({key, static_cast<float>(position.z()), static_cast<float>(std::atan2(position.y(), position.x()))});
      }
    }
  }

  for (unsigned int i = 0; i < m_trajectories->size(); ++i)
  {
    const auto* trajectory = m_trajectories->get(i);
    if (!trajectory || !trajectory->isValid()) continue;
    auto* candidate = new TpcSiliconMatchCandidate;
    candidate->set_parent_track_id(trajectory->get_parent_track_id());
    candidate->set_crossing(trajectory->get_crossing());
    float score = 0.F, maxDz = 0.F, maxDPhi = 0.F;
    unsigned int nMvtx = 0, nIntt = 0;
    for (unsigned int j = 0; j < trajectory->size_layer_states(); ++j)
    {
      const auto* state = trajectory->get_layer_state(j);
      if (!state || !state->valid || state->layer >= byLayer.size()) continue;
      const SiliconPoint* best = nullptr;
      float bestScore = std::numeric_limits<float>::max(), bestDz = 0.F, bestDPhi = 0.F;
      for (const auto& point : byLayer[state->layer])
      {
        const float dz = std::abs(point.z - state->z);
        const float dphi = std::abs(angleDifference(point.phi, state->phi));
        if (dz > m_maxDz || dphi > m_maxDPhi) continue;
        const float value = (dz / m_maxDz) * (dz / m_maxDz) + (dphi / m_maxDPhi) * (dphi / m_maxDPhi);
        if (value < bestScore) { bestScore = value; best = &point; bestDz = dz; bestDPhi = dphi; }
      }
      if (!best) continue;
      candidate->add_silicon_cluster_key(best->key);
      score += bestScore;
      maxDz = std::max(maxDz, bestDz);
      maxDPhi = std::max(maxDPhi, bestDPhi);
      if (state->layer <= 2U) ++nMvtx; else ++nIntt;
    }
    if (nMvtx < m_minMvtx || nIntt < m_minIntt) { delete candidate; continue; }
    candidate->set_n_mvtx(nMvtx);
    candidate->set_n_intt(nIntt);
    candidate->set_score(score);
    candidate->set_max_abs_dz(maxDz);
    candidate->set_max_abs_ddphi(maxDPhi);
    m_candidates->add(candidate);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
