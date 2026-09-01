#include "PHTrackCleaner.h"

#include "PHTrackCleaner.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>            // for cluskey, getLayer, TrkrId
#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>     // for sqrt, fabs, atan2, cos
#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for map
#include <set>       // for _Rb_tree_const_iterator
#include <utility>   // for pair, make_pair

//____________________________________________________________________________..
PHTrackCleaner::PHTrackCleaner(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHTrackCleaner::~PHTrackCleaner() = default;

//____________________________________________________________________________..
int PHTrackCleaner::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return ret;
}

//____________________________________________________________________________..
int PHTrackCleaner::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;
  }

  std::set<unsigned int> track_keep_list;
  std::set<unsigned int> track_delete_list;

  unsigned int good_track = 0;  // for diagnostic output only
  unsigned int ok_track = 0;    // tracks to keep

  auto get_tpc_group = [this](const SvtxTrack* track, unsigned int track_id)
  {
    if (!track)
    {
      return std::make_pair(2U, track_id);
    }

    auto tpc_seed = track->get_tpc_seed();
    const unsigned int tpc_index = _tpc_seed_map->find(tpc_seed);
    std::pair<unsigned int, unsigned int> tpc_group = std::make_pair(2U, track_id);
    if (tpc_seed)
    {
      if (tpc_index < _tpc_seed_map->size())
      {
        tpc_group = std::make_pair(0U, tpc_index);
        if (tpc_seed->get_tpc_seed_index() != UINT_MAX)
        {
          tpc_group = std::make_pair(1U, tpc_seed->get_tpc_seed_index());
        }
      }
      else if (tpc_seed->get_tpc_seed_index() < _tpc_seed_map->size())
      {
        const TrackSeed* indexed_tpc_seed = _tpc_seed_map->get(tpc_seed->get_tpc_seed_index());
        tpc_group = std::make_pair(0U, tpc_seed->get_tpc_seed_index());
        if (indexed_tpc_seed && indexed_tpc_seed->get_tpc_seed_index() != UINT_MAX)
        {
          tpc_group = std::make_pair(1U, indexed_tpc_seed->get_tpc_seed_index());
        }
      }
    }
    return tpc_group;
  };

  std::map<std::pair<unsigned int, unsigned int>, unsigned int> best_track_by_tpc_group;
  std::map<std::pair<unsigned int, unsigned int>, double> best_quality_by_tpc_group;

  // loop over the fitted tracks and keep only the best track per TPC group.
  // Crossing duplicates share the source assembled-track group.
  for (auto &it : *_track_map)
  {
    const unsigned int track_id = it.first;
    SvtxTrack* track = it.second;
    if (!track || track->get_ndf() <= min_ndf || track->get_ndf() == UINT_MAX)
    {
      continue;
    }

    const double qual = track->get_chisq() / track->get_ndf();
    if (!std::isfinite(qual) || qual >= quality_cut * 2)
    {
      continue;
    }

    const auto tpc_group = get_tpc_group(track, track_id);

    if (Verbosity() > 1)
    {
      unsigned int si_index = UINT_MAX;
      auto si_seed = track->get_silicon_seed();
      if (si_seed)
      {
        si_index = _silicon_seed_map->find(si_seed);
      }
      else
      {
        std::cout << "      no silicon seed found " << std::endl;
      }

      std::cout << "        track ID " << track_id << " TPC group type " << tpc_group.first
                << " id " << tpc_group.second << " si index " << si_index
                << " crossing " << track->get_crossing()
                << " chisq " << track->get_chisq() << " ndf " << track->get_ndf()
                << " chisq/ndf " << qual << std::endl;
    }

    const auto best_iter = best_quality_by_tpc_group.find(tpc_group);
    if (best_iter == best_quality_by_tpc_group.end() || qual < best_iter->second)
    {
      best_track_by_tpc_group[tpc_group] = track_id;
      best_quality_by_tpc_group[tpc_group] = qual;
    }
  }

  for (const auto& item : best_track_by_tpc_group)
  {
    const double qual = best_quality_by_tpc_group[item.first];
    track_keep_list.insert(item.second);
    ok_track++;
    if (qual < quality_cut)
    {
      good_track++;
    }

    if (Verbosity() > 1)
    {
      std::cout << "        best track for TPC group type " << item.first.first
                << " id " << item.first.second << " has track_id " << item.second
                << " chisq/ndf " << qual << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << " Number of good tracks with qual < " << quality_cut << "  is " << good_track << " OK tracks " << ok_track << std::endl;
  }

  // make a list of tracks that did not make the keep list
  for (auto &track_it : *_track_map)
  {
    auto id = track_it.first;

    auto set_it = track_keep_list.find(id);
    if (set_it == track_keep_list.end())
    {
      if (Verbosity() > 1)
      {
        std::cout << "    add id " << id << " to track_delete_list " << std::endl;
      }
      track_delete_list.insert(id);
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << " track_delete_list size " << track_delete_list.size() << std::endl;
  }

  // delete failed tracks
  for (unsigned int it : track_delete_list)
  {
    if (Verbosity() > 1)
    {
      std::cout << " erasing track ID " << it << std::endl;
    }
    _track_map->erase(it);
  }

  if (Verbosity() > 0)
  {
    std::cout << "Track map size after choosing best silicon match: " << _track_map->size() << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackCleaner::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackCleaner::GetNodes(PHCompositeNode *topNode)
{
  _tpc_seed_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_tpc_seed_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find TpcTrackSeedContainer: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _silicon_seed_map = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!_silicon_seed_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
