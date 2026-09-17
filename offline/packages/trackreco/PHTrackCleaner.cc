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
#include <limits>
#include <map>       // for map
#include <set>       // for _Rb_tree_const_iterator
#include <tuple>
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

  using SourceKey = std::pair<unsigned int, unsigned int>;
  using HypothesisKey = std::tuple<unsigned int, unsigned int, short int, unsigned int>;

  struct Candidate
  {
    unsigned int track_id{std::numeric_limits<unsigned int>::max()};
    SourceKey source_key{2U, std::numeric_limits<unsigned int>::max()};
    unsigned int source_id{std::numeric_limits<unsigned int>::max()};
    unsigned int tpc_seed_index{std::numeric_limits<unsigned int>::max()};
    unsigned int silicon_seed_index{std::numeric_limits<unsigned int>::max()};
    short int hypothesis_crossing{std::numeric_limits<short int>::max()};
    short int best_crossing{std::numeric_limits<short int>::max()};
    double quality{std::numeric_limits<double>::infinity()};
    bool has_silicon_seed{false};
    bool has_duplicate_identity{false};
  };

  const auto has_valid_best_crossing = [](const Candidate& candidate)
  {
    return candidate.best_crossing != std::numeric_limits<short int>::max();
  };

  const auto crossing_distance = [&has_valid_best_crossing](const Candidate& candidate)
  {
    return has_valid_best_crossing(candidate)
               ? std::abs(static_cast<int>(candidate.hypothesis_crossing) -
                          static_cast<int>(candidate.best_crossing))
               : std::numeric_limits<int>::max();
  };

  const auto crossing_distance_label = [&has_valid_best_crossing, &crossing_distance](const Candidate& candidate)
  {
    return candidate.has_duplicate_identity && has_valid_best_crossing(candidate)
               ? std::to_string(crossing_distance(candidate))
               : std::string("invalid");
  };

  const auto better_within_hypothesis = [](const Candidate& candidate, const Candidate& current)
  {
    return (candidate.has_silicon_seed && !current.has_silicon_seed) ||
           (candidate.has_silicon_seed == current.has_silicon_seed &&
            candidate.quality < current.quality);
  };

  const auto better_across_hypotheses =
      [&has_valid_best_crossing, &crossing_distance, &better_within_hypothesis](
          const Candidate& candidate, const Candidate& current)
  {
    if (candidate.has_duplicate_identity && current.has_duplicate_identity)
    {
      const bool candidate_best_valid = has_valid_best_crossing(candidate);
      const bool current_best_valid = has_valid_best_crossing(current);
      if (candidate_best_valid != current_best_valid)
      {
        return candidate_best_valid;
      }
      if (candidate_best_valid)
      {
        const int candidate_distance = crossing_distance(candidate);
        const int current_distance = crossing_distance(current);
        if (candidate_distance != current_distance)
        {
          return candidate_distance < current_distance;
        }
      }
    }

    return better_within_hypothesis(candidate, current);
  };

  std::map<HypothesisKey, Candidate> best_by_hypothesis;

  for (auto& [track_id, track] : *_track_map)
  {
    if (!track || track->get_ndf() <= min_ndf ||
        track->get_ndf() == std::numeric_limits<unsigned int>::max())
    {
      continue;
    }

    const double quality = track->get_chisq() / track->get_ndf();
    if (!std::isfinite(quality) || quality >= quality_cut * 2)
    {
      continue;
    }

    Candidate candidate;
    candidate.track_id = track_id;
    candidate.quality = quality;
    candidate.best_crossing = track->get_crossing();

    TrackSeed* tpc_seed = track->get_tpc_seed();
    if (tpc_seed)
    {
      const std::size_t tpc_index = _tpc_seed_map->find(tpc_seed);
      if (tpc_index < _tpc_seed_map->size())
      {
        candidate.tpc_seed_index = static_cast<unsigned int>(tpc_index);
      }

      candidate.source_id = tpc_seed->get_tpc_seed_index();
      candidate.hypothesis_crossing = tpc_seed->get_crossing();
      candidate.has_duplicate_identity =
          candidate.source_id != std::numeric_limits<unsigned int>::max() &&
          candidate.hypothesis_crossing != std::numeric_limits<short int>::max();
    }

    // A valid source ID groups crossing hypotheses from the same assembled
    // track. Legacy seeds fall back to their unique TPC container index.
    if (candidate.source_id != std::numeric_limits<unsigned int>::max())
    {
      candidate.source_key = {0U, candidate.source_id};
    }
    else if (candidate.tpc_seed_index != std::numeric_limits<unsigned int>::max())
    {
      candidate.source_key = {1U, candidate.tpc_seed_index};
      candidate.source_id = candidate.tpc_seed_index;
    }
    else
    {
      candidate.source_key = {2U, track_id};
      candidate.source_id = track_id;
    }

    TrackSeed* silicon_seed = track->get_silicon_seed();
    candidate.has_silicon_seed = silicon_seed != nullptr;
    if (silicon_seed)
    {
      const std::size_t silicon_index = _silicon_seed_map->find(silicon_seed);
      if (silicon_index < _silicon_seed_map->size())
      {
        candidate.silicon_seed_index = static_cast<unsigned int>(silicon_index);
      }
    }

    // Invalid hypothesis crossings must not merge distinct legacy TPC seeds.
    const unsigned int invalid_crossing_instance =
        candidate.hypothesis_crossing == std::numeric_limits<short int>::max()
            ? (candidate.tpc_seed_index != std::numeric_limits<unsigned int>::max()
                   ? candidate.tpc_seed_index
                   : track_id)
            : std::numeric_limits<unsigned int>::max();
    const HypothesisKey hypothesis_key{
        candidate.source_key.first,
        candidate.source_key.second,
        candidate.hypothesis_crossing,
        invalid_crossing_instance};

    if (Verbosity() > 1)
    {
      std::cout << "source=" << candidate.source_id
                << " hyp=" << candidate.hypothesis_crossing
                << " best=" << candidate.best_crossing
                << " delta=" << crossing_distance_label(candidate)
                << " track=" << candidate.track_id
                << " si=" << candidate.silicon_seed_index
                << " chi2/ndf=" << candidate.quality
                << std::endl;
    }

    const auto hypothesis_iter = best_by_hypothesis.find(hypothesis_key);
    if (hypothesis_iter == best_by_hypothesis.end() ||
        better_within_hypothesis(candidate, hypothesis_iter->second))
    {
      best_by_hypothesis[hypothesis_key] = candidate;
    }
  }

  std::map<SourceKey, Candidate> best_by_source;
  for (const auto& hypothesis_candidate : best_by_hypothesis)
  {
    const Candidate& candidate = hypothesis_candidate.second;
    if (Verbosity() > 1)
    {
      std::cout << "hypothesis winner: source=" << candidate.source_id
                << " hyp=" << candidate.hypothesis_crossing
                << " track=" << candidate.track_id
                << " best=" << candidate.best_crossing
                << " delta=" << crossing_distance_label(candidate)
                << " q=" << candidate.quality
                << std::endl;
    }

    const auto source_iter = best_by_source.find(candidate.source_key);
    if (source_iter == best_by_source.end() ||
        better_across_hypotheses(candidate, source_iter->second))
    {
      best_by_source[candidate.source_key] = candidate;
    }
  }

  std::set<unsigned int> track_keep_list;
  unsigned int good_track = 0;
  for (const auto& source_candidate : best_by_source)
  {
    const Candidate& candidate = source_candidate.second;
    track_keep_list.insert(candidate.track_id);
    if (candidate.quality < quality_cut)
    {
      ++good_track;
    }

    if (Verbosity() > 1)
    {
      std::cout << "source winner: source=" << candidate.source_id
                << " hyp=" << candidate.hypothesis_crossing
                << " best=" << candidate.best_crossing
                << " delta=" << crossing_distance_label(candidate)
                << " track=" << candidate.track_id
                << " q=" << candidate.quality
                << std::endl;
    }
  }

  std::set<unsigned int> track_delete_list;
  for (const auto& track_item : *_track_map)
  {
    if (!track_keep_list.contains(track_item.first))
    {
      track_delete_list.insert(track_item.first);
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << " Number of good tracks with qual < " << quality_cut
              << " is " << good_track
              << " OK tracks " << track_keep_list.size() << std::endl;
    std::cout << " track_delete_list size " << track_delete_list.size() << std::endl;
  }

  for (const unsigned int track_id : track_delete_list)
  {
    if (Verbosity() > 1)
    {
      std::cout << " erasing track ID " << track_id << std::endl;
    }
    _track_map->erase(track_id);
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
