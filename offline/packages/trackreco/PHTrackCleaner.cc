#include "PHTrackCleaner.h"

#include "PHTrackCleaner.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>            // for cluskey, getLayer, TrkrId
#include <trackbase/TpcDefs.h>
#include <tpc/TpcClusterZCrossingCorrection.h>
#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedHelper.h>

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
    short int silicon_crossing{std::numeric_limits<short int>::max()};
    double chisq{std::numeric_limits<double>::quiet_NaN()};
    unsigned int ndf{std::numeric_limits<unsigned int>::max()};
    double quality{std::numeric_limits<double>::infinity()};
    double z_residual{std::numeric_limits<double>::infinity()};
    float pt{std::numeric_limits<float>::quiet_NaN()};
    float eta{std::numeric_limits<float>::quiet_NaN()};
    float phi{std::numeric_limits<float>::quiet_NaN()};
    unsigned int ntpc{0};
    unsigned int nintt{0};
    unsigned int nmvtx{0};
    bool has_silicon_seed{false};
    bool has_duplicate_identity{false};
    bool has_valid_z_residual{false};
    bool tight_z_match{false};
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

  const auto z_residual_label = [](const Candidate& candidate)
  {
    return candidate.has_valid_z_residual
               ? std::to_string(candidate.z_residual)
               : std::string("INVALID");
  };

  const auto better_within_hypothesis = [](const Candidate& candidate, const Candidate& current)
  {
    return (candidate.has_silicon_seed && !current.has_silicon_seed) ||
           (candidate.has_silicon_seed == current.has_silicon_seed &&
            candidate.quality < current.quality);
  };

  const auto better_within_hypothesis_with_z =
      [&better_within_hypothesis](const Candidate& candidate, const Candidate& current)
  {
    if (candidate.tight_z_match != current.tight_z_match)
    {
      return candidate.tight_z_match;
    }
    return better_within_hypothesis(candidate, current);
  };

  const auto better_across_hypotheses_without_z =
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

  const auto better_across_hypotheses =
      [&better_across_hypotheses_without_z](const Candidate& candidate, const Candidate& current)
  {
    if (candidate.tight_z_match != current.tight_z_match)
    {
      return candidate.tight_z_match;
    }
    return better_across_hypotheses_without_z(candidate, current);
  };

  std::map<HypothesisKey, Candidate> best_by_hypothesis;
  std::map<HypothesisKey, Candidate> legacy_best_by_hypothesis;
  unsigned int candidate_total = 0;
  unsigned int candidate_tight_z = 0;

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
    candidate.chisq = track->get_chisq();
    candidate.ndf = track->get_ndf();
    candidate.quality = quality;
    candidate.pt = track->get_pt();
    candidate.eta = track->get_eta();
    candidate.phi = track->get_phi();
    candidate.best_crossing = track->get_crossing();

    const auto count_clusters = [&candidate](const TrackSeed* seed)
    {
      if (!seed)
      {
        return;
      }
      for (auto cluster_iter = seed->begin_cluster_keys();
           cluster_iter != seed->end_cluster_keys(); ++cluster_iter)
      {
        switch (TrkrDefs::getTrkrId(*cluster_iter))
        {
        case TrkrDefs::tpcId:
          ++candidate.ntpc;
          break;
        case TrkrDefs::inttId:
          ++candidate.nintt;
          break;
        case TrkrDefs::mvtxId:
          ++candidate.nmvtx;
          break;
        default:
          break;
        }
      }
    };

    TrackSeed* tpc_seed = track->get_tpc_seed();
    count_clusters(tpc_seed);
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
      count_clusters(silicon_seed);
      candidate.silicon_crossing = silicon_seed->get_crossing();
      const std::size_t silicon_index = _silicon_seed_map->find(silicon_seed);
      if (silicon_index < _silicon_seed_map->size())
      {
        candidate.silicon_seed_index = static_cast<unsigned int>(silicon_index);
      }
    }
    // Matcher's seed-level z consistency: correct the TPC seed z with the
    // candidate's matcher-selected best crossing, then compare to silicon.
    if (tpc_seed && silicon_seed && has_valid_best_crossing(candidate))
    {
      bool has_tpc_side = false;
      unsigned int tpc_side = 0;
      for (auto cluster_iter = tpc_seed->begin_cluster_keys();
           cluster_iter != tpc_seed->end_cluster_keys(); ++cluster_iter)
      {
        if (TrkrDefs::getTrkrId(*cluster_iter) == TrkrDefs::tpcId)
        {
          tpc_side = TpcDefs::getSide(*cluster_iter);
          has_tpc_side = true;
          break;
        }
      }

      if (has_tpc_side)
      {
        const double tpc_z = TrackSeedHelper::get_z(tpc_seed);
        const double silicon_z = TrackSeedHelper::get_z(silicon_seed);
        if (std::isfinite(tpc_z) && std::isfinite(silicon_z))
        {
          const double corrected_tpc_z =
              TpcClusterZCrossingCorrection::correctZ(
                  tpc_z, tpc_side, candidate.best_crossing);
          if (std::isfinite(corrected_tpc_z))
          {
            candidate.z_residual = std::abs(corrected_tpc_z - silicon_z);
            candidate.has_valid_z_residual = true;
            candidate.tight_z_match = candidate.z_residual < m_tight_z_match_max;
          }
        }
      }
    }

    ++candidate_total;
    if (candidate.tight_z_match)
    {
      ++candidate_tight_z;
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
                << " zres=" << z_residual_label(candidate)
                << " tightZ=" << candidate.tight_z_match
                << " si=" << candidate.silicon_seed_index
                << " chi2/ndf=" << candidate.quality
                << std::endl;
    }

    const auto legacy_hypothesis_iter = legacy_best_by_hypothesis.find(hypothesis_key);
    if (legacy_hypothesis_iter == legacy_best_by_hypothesis.end() ||
        better_within_hypothesis(candidate, legacy_hypothesis_iter->second))
    {
      legacy_best_by_hypothesis[hypothesis_key] = candidate;
    }

    const auto hypothesis_iter = best_by_hypothesis.find(hypothesis_key);
    if (hypothesis_iter == best_by_hypothesis.end() ||
        better_within_hypothesis_with_z(candidate, hypothesis_iter->second))
    {
      best_by_hypothesis[hypothesis_key] = candidate;
    }
  }

std::map<SourceKey, Candidate> quality_by_source;
  std::map<SourceKey, Candidate> best_by_source;
  std::map<SourceKey, Candidate> previous_stage2_by_source;
  std::map<SourceKey, Candidate> old_style_by_source;
  std::map<SourceKey, std::vector<Candidate>> hypothesis_winners_by_source;
  for (const auto& hypothesis_candidate : best_by_hypothesis)
  {
    const Candidate& candidate = hypothesis_candidate.second;
    hypothesis_winners_by_source[candidate.source_key].push_back(candidate);
    if (Verbosity() > 1)
    {
      std::cout << "hypothesis winner: source=" << candidate.source_id
                << " hyp=" << candidate.hypothesis_crossing
                << " track=" << candidate.track_id
                << " best=" << candidate.best_crossing
                << " delta=" << crossing_distance_label(candidate)
                << " q=" << candidate.quality
                << " zres=" << z_residual_label(candidate)
                << " tightZ=" << candidate.tight_z_match
                << std::endl;
    }
    const auto quality_iter = quality_by_source.find(candidate.source_key);

    if (quality_iter == quality_by_source.end() ||
        better_within_hypothesis(candidate, quality_iter->second))
    {
      quality_by_source[candidate.source_key] = candidate;
    }

    const auto source_iter = best_by_source.find(candidate.source_key);
    const auto previous_source_iter = previous_stage2_by_source.find(candidate.source_key);
    if (previous_source_iter == previous_stage2_by_source.end() ||
        better_across_hypotheses_without_z(candidate, previous_source_iter->second))
    {
      previous_stage2_by_source[candidate.source_key] = candidate;
    }

    if (source_iter == best_by_source.end() ||
        better_across_hypotheses(candidate, source_iter->second))
    {
      best_by_source[candidate.source_key] = candidate;
    }

    const auto old_style_iter = old_style_by_source.find(candidate.source_key);
    if (old_style_iter == old_style_by_source.end() ||
        better_within_hypothesis(candidate, old_style_iter->second))
    {
      old_style_by_source[candidate.source_key] = candidate;
    }
  }
  unsigned int hypothesis_winners_tight_z = 0;
  for (const auto& [key, candidate] : best_by_hypothesis)
  {
    if (candidate.tight_z_match)
    {
      ++hypothesis_winners_tight_z;
    }
  }

  unsigned int source_winners_tight_z = 0;
  unsigned int stage1_changed_by_tight_z = 0;
  unsigned int stage2_changed_by_tight_z = 0;
  for (const auto& [key, candidate] : best_by_hypothesis)
  {
    const auto legacy_iter = legacy_best_by_hypothesis.find(key);
    if (legacy_iter != legacy_best_by_hypothesis.end() &&
        candidate.track_id != legacy_iter->second.track_id)
    {
      ++stage1_changed_by_tight_z;
    }
  }
  for (const auto& [key, candidate] : best_by_source)
  {
    if (candidate.tight_z_match)
    {
      ++source_winners_tight_z;
    }
    const auto previous_iter = previous_stage2_by_source.find(key);
    if (previous_iter != previous_stage2_by_source.end() &&
        candidate.track_id != previous_iter->second.track_id)
    {
      ++stage2_changed_by_tight_z;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "PHTrackCleaner tight-z summary:"
              << " candidate_total=" << candidate_total
              << " candidate_tight_z=" << candidate_tight_z
              << " hypothesis_winners=" << best_by_hypothesis.size()
              << " hypothesis_winners_tight_z=" << hypothesis_winners_tight_z
              << " source_winners=" << best_by_source.size()
              << " source_winners_tight_z=" << source_winners_tight_z
              << " stage1_changed_by_tight_z=" << stage1_changed_by_tight_z
              << " stage2_changed_by_tight_z=" << stage2_changed_by_tight_z
              << std::endl;
  }

  if (Verbosity() > 1)
  {
    const auto crossing_label = [](const short int crossing)
    {
      return crossing == std::numeric_limits<short int>::max()
                 ? std::string("INVALID")
                 : std::to_string(crossing);
    };
    const auto seed_index_label = [](const Candidate& candidate)
    {
      return candidate.has_silicon_seed &&
                     candidate.silicon_seed_index != std::numeric_limits<unsigned int>::max()
                 ? std::to_string(candidate.silicon_seed_index)
                 : std::string("NONE");
    };
    const auto tpc_seed_index_label = [](const Candidate& candidate)
    {
      return candidate.tpc_seed_index != std::numeric_limits<unsigned int>::max()
                 ? std::to_string(candidate.tpc_seed_index)
                 : std::string("INVALID");
    };
    const auto print_candidate = [&crossing_distance_label, &z_residual_label, &crossing_label, &seed_index_label, &tpc_seed_index_label](
                                     const char* prefix, const Candidate& candidate)
    {
      std::cout << prefix
                << "source=" << candidate.source_id
                << " track=" << candidate.track_id
                << " tpc=" << tpc_seed_index_label(candidate)
                << " si=" << seed_index_label(candidate)
                << " siCross=" << (candidate.has_silicon_seed
                                         ? crossing_label(candidate.silicon_crossing)
                                         : std::string("NONE"))
                << " hyp=" << crossing_label(candidate.hypothesis_crossing)
                << " best=" << crossing_label(candidate.best_crossing)
                << " delta=" << crossing_distance_label(candidate)
                << " zres=" << z_residual_label(candidate)
                << " tightZ=" << candidate.tight_z_match
                << " chisq=" << candidate.chisq
                << " ndf=" << candidate.ndf
                << " q=" << candidate.quality
                << " pt=" << candidate.pt
                << " eta=" << candidate.eta
                << " phi=" << candidate.phi
                << " nTPC=" << candidate.ntpc
                << " nINTT=" << candidate.nintt
                << " nMVTX=" << candidate.nmvtx
                << std::endl;
    };

    for (const auto& [source_key, old_candidate] : old_style_by_source)
    {
      const Candidate& new_candidate = best_by_source.at(source_key);
      if (old_candidate.track_id == new_candidate.track_id)
      {
        continue;
      }

      std::cout << "CLEANER_DIFF source=" << old_candidate.source_id << std::endl;
      std::cout << "  candidates:" << std::endl;
      for (const Candidate& candidate : hypothesis_winners_by_source.at(source_key))
      {
        print_candidate("    ", candidate);
      }
      print_candidate("  OLDSTYLE -> ", old_candidate);
      print_candidate("  NEW      -> ", new_candidate);
    }
    if (Verbosity() > 1)
{
  for (const auto& [source_key, selected] : best_by_source)
  {
    const auto qiter = quality_by_source.find(source_key);
    if (qiter == quality_by_source.end())
    {
      continue;
    }

    const Candidate& quality = qiter->second;

    if (quality.track_id == selected.track_id)
    {
      continue;
    }

    const SvtxTrack* quality_track = _track_map->get(quality.track_id);
    const SvtxTrack* selected_track = _track_map->get(selected.track_id);

    std::cout
        << "DIAG_CLEAN_CHOICE"
        << " source=" << selected.source_id

        << " qualityTrack=" << quality.track_id
        << " qualityHyp=" << quality.hypothesis_crossing
        << " qualityBest=" << quality.best_crossing
        << " qualityDelta=" << crossing_distance_label(quality)
        << " qualityQ=" << quality.quality
        << " qualityZres=" << z_residual_label(quality)
        << " qualityTightZ=" << quality.tight_z_match
        << " qualityPt="
        << (quality_track
                ? std::hypot(quality_track->get_px(), quality_track->get_py())
                : -1.0)

        << " selectedTrack=" << selected.track_id
        << " selectedHyp=" << selected.hypothesis_crossing
        << " selectedBest=" << selected.best_crossing
        << " selectedDelta=" << crossing_distance_label(selected)
        << " selectedQ=" << selected.quality
        << " selectedZres=" << z_residual_label(selected)
        << " selectedTightZ=" << selected.tight_z_match
        << " selectedPt="
        << (selected_track
                ? std::hypot(selected_track->get_px(), selected_track->get_py())
                : -1.0)

        << std::endl;
  }
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
                << " zres=" << z_residual_label(candidate)
                << " tightZ=" << candidate.tight_z_match
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
