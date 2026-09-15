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
#include <map>
#include <chrono>
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
  if (!m_beamFrame.validate())
  {
    std::cerr << Name() << "::InitRun - invalid BeamFrameTransform" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  return createNodes(topNode);
}

int TpcSiliconCrossingMatcher::process_event(PHCompositeNode*)
{
  const auto begin = std::chrono::steady_clock::now();
  m_candidates->Reset();
  using Bin = std::pair<int, int>;
  std::array<std::map<Bin, std::vector<SiliconPoint>>, 7> byLayer;
  unsigned int indexedClusters = 0;
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
        const auto& beam = layer <= 2U ? m_beamFrame.mvtxBeamLine() : m_beamFrame.inttBeamLine();
        const auto centered = m_beamFrame.toBeamFrame(beam, position.x(), position.y(), position.z());
        const float phi = std::atan2(centered.y, centered.x);
        const int phiBin = static_cast<int>(std::floor((phi + static_cast<float>(M_PI)) / m_maxDPhi));
        const int zBin = static_cast<int>(std::floor(centered.z / m_maxDz));
        byLayer[layer][{phiBin, zBin}].push_back({key, static_cast<float>(centered.z), phi});
        ++indexedClusters;
      }
    }
  }

  for (unsigned int i = 0; i < m_trajectories->size(); ++i)
  {
    const auto* trajectory = m_trajectories->get(i);
    if (!trajectory || !trajectory->isValid()) continue;
    auto* candidate = new TpcSiliconMatchCandidate;
    candidate->set_parent_track_id(trajectory->get_parent_track_id());
    candidate->set_source_assembled_track_id(trajectory->get_source_assembled_track_id());
    candidate->set_crossing(trajectory->get_crossing());
    float score = 0.F, maxDz = 0.F, maxDPhi = 0.F;
    unsigned int nMvtx = 0, nIntt = 0;
    float previousDPhi = 0.F;
    bool hasPrevious = false;
    for (unsigned int j = 0; j < trajectory->size_layer_states(); ++j)
    {
      const auto* state = trajectory->get_layer_state(j);
      if (!state || !state->valid || state->layer >= byLayer.size()) continue;
      const auto predicted = m_beamFrame.toBeamFrame(m_beamFrame.tpcBeamLine(), state->x, state->y, state->z);
      const float predictedPhi = std::atan2(predicted.y, predicted.x);
      const int phiBin = static_cast<int>(std::floor((predictedPhi + static_cast<float>(M_PI)) / m_maxDPhi));
      const int zBin = static_cast<int>(std::floor(predicted.z / m_maxDz));
      const SiliconPoint* best = nullptr;
      float bestScore = std::numeric_limits<float>::max(), bestDz = 0.F, bestDPhi = 0.F;
      for (int dp = -1; dp <= 1; ++dp) for (int dzBin = -1; dzBin <= 1; ++dzBin)
      {
        const auto found = byLayer[state->layer].find({phiBin + dp, zBin + dzBin});
        if (found == byLayer[state->layer].end()) continue;
        for (const auto& point : found->second)
        {
          const float dz = point.z - predicted.z;
          const float dphi = angleDifference(point.phi, predictedPhi);
          const float expected = hasPrevious ? static_cast<float>(m_phiOffset[state->layer] + m_phiSlope[state->layer] * previousDPhi) : 0.F;
          const float ddphi = angleDifference(dphi, expected);
          if (std::abs(dz) > m_maxDz || std::abs(ddphi) > m_maxDPhi) continue;
          const float value = (dz / m_maxDz) * (dz / m_maxDz) + (ddphi / m_maxDPhi) * (ddphi / m_maxDPhi);
          if (value < bestScore) { bestScore = value; best = &point; bestDz = std::abs(dz); bestDPhi = std::abs(ddphi); }
        }
      }
      if (!best) continue;
      candidate->add_silicon_cluster_key(best->key);
      score += bestScore;
      maxDz = std::max(maxDz, bestDz);
      maxDPhi = std::max(maxDPhi, bestDPhi);
      previousDPhi = angleDifference(best->phi, predictedPhi);
      hasPrevious = true;
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
  if (Verbosity() > 0) std::cout << Name() << " indexed_silicon_clusters=" << indexedClusters << " candidates=" << m_candidates->size() << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
