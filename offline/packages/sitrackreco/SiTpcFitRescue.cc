#include "SiTpcFitRescue.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackSeed_v3.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/sphenix_constants.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <set>

SiTpcFitRescue::SiTpcFitRescue(Mode mode, const std::string& name)
  : SubsysReco(name)
  , m_mode(mode)
{
}

int SiTpcFitRescue::InitRun(PHCompositeNode* topNode)
{
  if (m_mode == Mode::MakeSeeds)
  {
    // the rescue seed container must exist before the 2nd PHActsTrkFitter's InitRun (it looks for it there)
    if (!findNode::getClass<TrackSeedContainer>(topNode, m_rescueSeedName))
    {
      PHNodeIterator iter(topNode);
      auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
      if (!dstNode)
      {
        std::cerr << Name() << "::InitRun - DST node missing" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
      PHNodeIterator dstIter(dstNode);
      auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
      if (!svtxNode)
      {
        svtxNode = new PHCompositeNode("SVTX");
        dstNode->addNode(svtxNode);
      }
      svtxNode->addNode(new PHIODataNode<PHObject>(new TrackSeedContainer_v1, m_rescueSeedName, "PHObject"));
    }
    std::cout << Name() << "::InitRun - make seeds: Si+TPC seeds of " << m_svtxSeedName << " without a track in " << m_trackMapName
              << " -> TPC-only seeds with a crossing carrier silicon seed in " << m_rescueSeedName
              << (m_retryTpcOnly ? " (failed TPC-only seeds retried too)" : "") << std::endl;
  }
  else
  {
    std::cout << Name() << "::InitRun - merge tracks: " << m_rescueTrackMapName << " appended to " << m_trackMapName << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpcFitRescue::process_event(PHCompositeNode* topNode)
{
  ++m_nEvents;
  return m_mode == Mode::MakeSeeds ? makeSeeds(topNode) : mergeTracks(topNode);
}

int SiTpcFitRescue::makeSeeds(PHCompositeNode* topNode)
{
  auto* svtxSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_svtxSeedName);
  auto* siSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_siSeedName);
  auto* tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_tpcSeedName);
  auto* trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto* rescue = findNode::getClass<TrackSeedContainer>(topNode, m_rescueSeedName);
  auto* geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!svtxSeeds || !siSeeds || !tpcSeeds || !trackMap || !rescue || !geometry)
  {
    std::cerr << Name() << ": missing node (" << m_svtxSeedName << ", " << m_siSeedName << ", " << m_tpcSeedName << ", "
              << m_trackMapName << ", " << m_rescueSeedName << " or ActsGeometry)" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  rescue->Reset();

  // TPC seeds that have a fitted track
  std::set<std::size_t> fittedTpc;
  for (const auto& [key, track] : *trackMap)
  {
    if (!track)
    {
      continue;
    }
    if (const TrackSeed* t = track->get_tpc_seed())
    {
      fittedTpc.insert(tpcSeeds->find(t));
    }
    else if (track->get_tpc_seed_id() != USHRT_MAX)
    {
      fittedTpc.insert(track->get_tpc_seed_id());  // slim tracks
    }
  }

  const double driftVelocity = geometry->get_drift_velocity();
  const double zBunchSeparation = sphenix_constants::time_between_crossings * driftVelocity;
  unsigned int nEvent = 0;
  for (std::size_t i = 0; i < svtxSeeds->size(); ++i)
  {
    const TrackSeed* svtx = svtxSeeds->get(i);
    if (!svtx)
    {
      continue;
    }
    ++m_nSeeds;
    const unsigned int tpcIndex = svtx->get_tpc_seed_index();
    const TrackSeed* tpcSeed = tpcSeeds->get(tpcIndex);
    if (!tpcSeed || fittedTpc.contains(tpcIndex))
    {
      continue;
    }
    ++m_nFailed;
    const TrackSeed* siSeed = siSeeds->get(svtx->get_silicon_seed_index());
    const bool tpcOnly = !siSeed || siSeed->size_cluster_keys() == 0;
    if (tpcOnly)
    {
      ++m_nFailedTpcOnly;
      if (!m_retryTpcOnly)
      {
        continue;
      }
    }
    // crossing of the seed (silicon seed, else the SvtxTrackSeed, else the TPC seed) and the TPC side
    short crossing = siSeed ? siSeed->get_crossing() : SHRT_MAX;
    if (crossing == SHRT_MAX)
    {
      crossing = svtx->get_crossing();
    }
    if (crossing == SHRT_MAX)
    {
      crossing = tpcSeed->get_crossing();
    }
    int side = -1;
    for (auto it = tpcSeed->begin_cluster_keys(); it != tpcSeed->end_cluster_keys(); ++it)
    {
      if (TrkrDefs::getTrkrId(*it) == TrkrDefs::tpcId)
      {
        side = static_cast<int>(TpcDefs::getSide(TrkrDefs::getHitSetKeyFromClusKey(*it)));
        break;
      }
    }
    if (crossing == SHRT_MAX || side < 0 || !(driftVelocity > 0))
    {
      ++m_nNoCrossing;
      continue;
    }
    // TPC seed z0 is stored for crossing 0 (TpcPolyTrackSeedConverter): undo that for the carrier
    const double z0 = (side == 0) ? tpcSeed->get_Z0() - static_cast<double>(crossing) * zBunchSeparation
                                  : tpcSeed->get_Z0() + static_cast<double>(crossing) * zBunchSeparation;
    TrackSeed_v2 carrier;  // no cluster keys: TPC (+ TPOT) only fit at this crossing
    carrier.set_X0(tpcSeed->get_X0());
    carrier.set_Y0(tpcSeed->get_Y0());
    carrier.set_qOverR(tpcSeed->get_qOverR());
    carrier.set_Z0(static_cast<float>(z0));
    carrier.set_slope(tpcSeed->get_slope());
    carrier.set_phi(tpcSeed->get_phi());
    carrier.set_crossing(crossing);
    siSeeds->insert(&carrier);

    SvtxTrackSeed_v3 seed;
    seed.set_silicon_seed_index(static_cast<unsigned int>(siSeeds->size() - 1));
    seed.set_tpc_seed_index(tpcIndex);
    seed.set_crossing(crossing);
    seed.set_crossing_estimate(crossing);
    rescue->insert(&seed);
    ++m_nRescueSeeds;
    ++nEvent;
  }
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << m_nEvents << " " << trackMap->size() << " fitted tracks, " << nEvent
              << " Si+TPC seeds without a track -> TPC-only rescue seeds" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpcFitRescue::mergeTracks(PHCompositeNode* topNode)
{
  auto* trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto* rescued = findNode::getClass<SvtxTrackMap>(topNode, m_rescueTrackMapName);
  if (!trackMap || !rescued)
  {
    std::cerr << Name() << ": missing " << m_trackMapName << " or " << m_rescueTrackMapName << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  unsigned int next = 0;
  for (const auto& [key, track] : *trackMap)
  {
    next = std::max(next, key + 1);
  }
  unsigned int nEvent = 0;
  for (const auto& [key, track] : *rescued)
  {
    if (!track)
    {
      continue;
    }
    if (SvtxTrack* copy = trackMap->insertWithKey(track, next))
    {
      copy->set_id(next);
    }
    ++next;
    ++nEvent;
  }
  m_nMerged += nEvent;
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << m_nEvents << " " << nEvent << " rescued TPC-only tracks appended, " << m_trackMapName
              << " size " << trackMap->size() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpcFitRescue::End(PHCompositeNode* /*topNode*/)
{
  if (m_mode == Mode::MakeSeeds)
  {
    std::cout << Name() << "::End - " << m_nEvents << " events, " << m_nSeeds << " SvtxTrackSeeds, " << m_nFailed
              << " without a fitted track (" << m_nFailedTpcOnly << " TPC-only" << (m_retryTpcOnly ? ", retried" : ", not retried")
              << "), rescue seeds " << m_nRescueSeeds << ", no crossing / side " << m_nNoCrossing << std::endl;
  }
  else
  {
    std::cout << Name() << "::End - " << m_nEvents << " events, rescued TPC-only tracks appended " << m_nMerged << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
