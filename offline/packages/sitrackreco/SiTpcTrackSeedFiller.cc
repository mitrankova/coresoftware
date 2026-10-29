#include "SiTpcTrackSeedFiller.h"

#include "SiTpc_Track.h"
#include "SiTpc_TrackContainer.h"
#include "Si_Trajectory.h"
#include "Si_TrajectoryContainer.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrDefs.h>

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
#include <phool/phool.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  constexpr double cm = Acts::UnitConstants::cm;  // 10 (Acts length unit = mm)
  constexpr double twopi = 2. * M_PI;
  constexpr double kStraightRadius = 1.0e4;  // [cm] radius used for trajectories without curvature

  // window of the position match chain point - official cluster [cm] (local 0, local 1)
  constexpr std::array<double, 2> kMvtxMatchWindow{0.03, 0.03};
  constexpr std::array<double, 2> kInttMatchWindow{0.05, 1.5};
  // errors of projected chain points [cm] (r*phi, z)
  constexpr std::array<double, 2> kMvtxChainErr{0.0010, 0.0010};
  constexpr std::array<double, 2> kInttChainErr{0.0025, 0.46};

  template <class T>
  T square(T x)
  {
    return x * x;
  }
}  // namespace

//____________________________________________________________________________________________________
struct SiTpcTrackSeedFiller::Impl
{
  ActsGeometry* geometry = nullptr;
  SiTpc_TrackContainer* siTpcTracks = nullptr;
  Si_TrajectoryContainer* siTrajectories = nullptr;
  Tpc_PolyTrackContainer* tpcPolyTracks = nullptr;
  TpcCrossingDecisionContainer* crossings = nullptr;
  TrkrClusterContainer* clusters = nullptr;
  TrkrClusterHitAssoc* hitAssoc = nullptr;
  TrackSeedContainer* tpcSeeds = nullptr;
  TrackSeedContainer* siSeeds = nullptr;
  TrackSeedContainer* svtxSeeds = nullptr;
  double bzSign = 1;

  std::unordered_map<TrkrDefs::hitsetkey, std::unordered_map<TrkrDefs::hitkey, TrkrDefs::cluskey>> hitToCluster;
  std::map<TrkrDefs::hitsetkey, uint32_t> nextIndex;

  unsigned long nEvents = 0, nTrajectories = 0, nSiSeeds = 0, nSiTpc = 0, nSvtxSeeds = 0, nTpcOnly = 0;
  unsigned long nSiOfficialAssoc = 0, nSiOfficialPosition = 0, nSiProjected = 0, nSiMissing = 0, nSiOtherStrobe = 0;
  std::map<std::string, unsigned long> rejected;

  TrkrDefs::cluskey clusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk);
  TrkrDefs::cluskey officialClusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, bool& otherStrobe);
  uint32_t freeIndex(TrkrDefs::hitsetkey hsk);
};

//____________________________________________________________________________________________________
TrkrDefs::cluskey SiTpcTrackSeedFiller::Impl::clusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk)
{
  if (!hitAssoc)
  {
    return TrkrDefs::CLUSKEYMAX;
  }
  auto it = hitToCluster.find(hsk);
  if (it == hitToCluster.end())
  {
    std::unordered_map<TrkrDefs::hitkey, TrkrDefs::cluskey> m;
    if (auto* map = hitAssoc->getClusterMap(hsk))
    {
      for (const auto& [ckey, hitkey] : *map)
      {
        m.emplace(hitkey, ckey);
      }
    }
    it = hitToCluster.emplace(hsk, std::move(m)).first;
  }
  const auto c = it->second.find(hk);
  return c == it->second.end() ? TrkrDefs::CLUSKEYMAX : c->second;
}

//____________________________________________________________________________________________________
TrkrDefs::cluskey SiTpcTrackSeedFiller::Impl::officialClusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, bool& otherStrobe)
{
  // MvtxClusterPruner removes the copy of a cluster duplicated in the next strobe from TRKR_CLUSTER without
  // updating TRKR_CLUSTERHITASSOC: then the cluster with the same pixel in the neighbouring strobe is used
  otherStrobe = false;
  const auto ckey = clusterOfHit(hsk, hk);
  if (ckey != TrkrDefs::CLUSKEYMAX && clusters->findCluster(ckey))
  {
    return ckey;
  }
  if (TrkrDefs::getTrkrId(hsk) != TrkrDefs::mvtxId)
  {
    return TrkrDefs::CLUSKEYMAX;
  }
  const auto layer = TrkrDefs::getLayer(hsk);
  const auto stave = MvtxDefs::getStaveId(hsk);
  const auto chip = MvtxDefs::getChipId(hsk);
  const int strobe = MvtxDefs::getStrobeId(hsk);
  for (const int ds : {-1, 1})
  {
    const auto other = clusterOfHit(MvtxDefs::genHitSetKey(layer, stave, chip, strobe + ds), hk);
    if (other != TrkrDefs::CLUSKEYMAX && clusters->findCluster(other))
    {
      otherStrobe = true;
      return other;
    }
  }
  return TrkrDefs::CLUSKEYMAX;
}

//____________________________________________________________________________________________________
uint32_t SiTpcTrackSeedFiller::Impl::freeIndex(TrkrDefs::hitsetkey hsk)
{
  // TrkrClusterContainerv4 stores the clusters of a hitset in a vector indexed by the cluster index:
  // new keys continue after the largest index of the hitset
  auto it = nextIndex.find(hsk);
  if (it == nextIndex.end())
  {
    uint32_t next = 0;
    const auto range = clusters->getClusters(hsk);
    for (auto c = range.first; c != range.second; ++c)
    {
      next = std::max(next, TrkrDefs::getClusIndex(c->first) + 1);
    }
    it = nextIndex.emplace(hsk, next).first;
  }
  return it->second++;
}

//____________________________________________________________________________________________________
SiTpcTrackSeedFiller::SiTpcTrackSeedFiller(const std::string& name)
  : SubsysReco(name)
  , m_impl(std::make_unique<Impl>())
{
}

SiTpcTrackSeedFiller::~SiTpcTrackSeedFiller() = default;

//____________________________________________________________________________________________________
int SiTpcTrackSeedFiller::InitRun(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
  d.geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  d.siTpcTracks = findNode::getClass<SiTpc_TrackContainer>(topNode, m_siTpcTrackNodeName);
  d.siTrajectories = findNode::getClass<Si_TrajectoryContainer>(topNode, m_siTrajectoryNodeName);
  d.tpcPolyTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcPolyTrackNodeName);
  d.crossings = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingNodeName);
  d.clusters = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterNodeName);
  d.hitAssoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  d.tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_tpcSeedName);

  bool ok = true;
  auto require = [&](const void* p, const std::string& what)
  {
    if (!p)
    {
      std::cout << PHWHERE << " " << Name() << ": missing " << what << std::endl;
      ok = false;
    }
  };
  require(d.geometry, "ActsGeometry");
  require(d.siTpcTracks, m_siTpcTrackNodeName);
  require(d.siTrajectories, m_siTrajectoryNodeName);
  require(d.tpcPolyTracks, m_tpcPolyTrackNodeName);
  require(d.crossings, m_crossingNodeName);
  require(d.clusters, m_clusterNodeName);
  require(d.tpcSeeds, m_tpcSeedName + " (register TpcPolyTrackSeedConverter before " + Name() + ")");
  if (!d.hitAssoc)
  {
    std::cout << Name() << ": no TRKR_CLUSTERHITASSOC - Si clusters matched by position only" << std::endl;
  }
  if (findNode::getClass<TrackSeedContainer>(topNode, m_siSeedName))
  {
    std::cout << PHWHERE << " " << Name() << ": " << m_siSeedName << " already exists (standard silicon seeding registered?)."
              << " The own Si seeds replace the standard seeding: switch it off." << std::endl;
    ok = false;
  }
  if (!ok)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // output containers
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " " << Name() << ": DST node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator dstIter(dstNode);
  auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }
  d.siSeeds = new TrackSeedContainer_v1;
  svtxNode->addNode(new PHIODataNode<PHObject>(d.siSeeds, m_siSeedName, "PHObject"));
  d.svtxSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_svtxSeedName);
  if (!d.svtxSeeds)
  {
    d.svtxSeeds = new TrackSeedContainer_v1;
    svtxNode->addNode(new PHIODataNode<PHObject>(d.svtxSeeds, m_svtxSeedName, "PHObject"));
  }

  // charge: from the sense of rotation and the field (positive charge, B_z > 0: clockwise).
  // With B_z > 0 this is the convention of TrackSeedHelper::circleFitByTaubin.
  {
    const auto& field = *d.geometry->geometry().magField;
    auto cache = field.makeCache(d.geometry->geometry().magFieldContext);
    const auto b = field.getField(Acts::Vector3(0, 0, 0), cache);
    d.bzSign = (b.ok() && b.value().z() < 0) ? -1. : 1.;
    std::cout << Name() << "::InitRun - B_z at the centre " << (b.ok() ? b.value().z() / Acts::UnitConstants::T : NAN) << " T" << std::endl;
  }

  std::cout << Name() << "::InitRun - Si beam spot (" << m_siBeam[0] << ", " << m_siBeam[1] << ") cm, Si acts z shift "
            << m_siActsZShift << " cm, projected chain points " << (m_useProjected ? "on" : "off")
            << ", TPC-only seeds " << (m_addUnmatchedTpc ? "on" : "off") << std::endl
            << "  output " << m_siSeedName << ", " << m_svtxSeedName << " (TPC seeds from " << m_tpcSeedName << ")" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcTrackSeedFiller::process_event(PHCompositeNode* /*topNode*/)
{
  auto& d = *m_impl;
  ++d.nEvents;
  d.siSeeds->Reset();
  d.svtxSeeds->Reset();
  d.hitToCluster.clear();
  d.nextIndex.clear();

  const auto gctx = d.geometry->geometry().getGeoContext();
  auto reject = [&](const std::string& why)
  {
    ++d.rejected[why];
  };

  // ---- crossing of a trajectory: of its Si+TPC track (if any)
  auto crossingOfAssembled = [&](unsigned int assembledId) -> short
  {
    const TpcCrossingDecision* decision = d.crossings->get_decision(assembledId);
    return decision ? decision->get_selected_crossing() : static_cast<short>(SHRT_MAX);
  };
  std::map<unsigned int, short> crossingOfTrajectory;
  for (unsigned int i = 0; i < d.siTpcTracks->size(); ++i)
  {
    if (const SiTpc_Track* t = d.siTpcTracks->get_track(i))
    {
      crossingOfTrajectory.emplace(t->get_si_trajectory_index(), crossingOfAssembled(t->get_tpc_assembled_track_id()));
    }
  }

  // ---- Si seeds: one per own Si trajectory
  std::map<unsigned int, unsigned int> siSeedOfTrajectory;
  for (unsigned int itraj = 0; itraj < d.siTrajectories->size(); ++itraj)
  {
    const Si_Trajectory* traj = d.siTrajectories->get_trajectory(itraj);
    if (!traj)
    {
      continue;
    }
    ++d.nTrajectories;

    // helix in the official frame
    double cx = traj->get_circle_x() + m_siBeam[0];
    double cy = traj->get_circle_y() + m_siBeam[1];
    double R = traj->get_radius();
    const double phiPca = traj->get_phi();
    const double pcaX = traj->get_pca_x() + m_siBeam[0];
    const double pcaY = traj->get_pca_y() + m_siBeam[1];
    const double tanl = traj->get_tanl();
    const double z0 = traj->get_z0() - m_siActsZShift;  // at the pca to the own beam axis (= Si beam spot)
    if (!std::isfinite(phiPca) || !std::isfinite(pcaX) || !std::isfinite(pcaY) || !std::isfinite(tanl) || !std::isfinite(z0))
    {
      reject("Si trajectory without valid fit");
      continue;
    }
    double omega = 1;  // +1: counter-clockwise seen from +z
    if (std::isfinite(R) && R > 0 && std::isfinite(cx) && std::isfinite(cy))
    {
      omega = ((pcaX - cx) * std::sin(phiPca) - (pcaY - cy) * std::cos(phiPca)) >= 0 ? 1. : -1.;
    }
    else
    {
      // no curvature: large circle through the pca, tangent to phi
      R = kStraightRadius;
      cx = pcaX - std::sin(phiPca) * R;
      cy = pcaY + std::cos(phiPca) * R;
    }
    const double r0 = std::hypot(cx, cy);
    if (!(r0 > 0))
    {
      reject("Si trajectory with degenerate circle");
      continue;
    }
    auto angleOf = [&](double x, double y)
    { return std::atan2(y - cy, x - cx); };
    // point of the circle closest to the origin: z and direction there
    const double ox = cx - R * cx / r0, oy = cy - R * cy / r0;
    const double ds = R * std::remainder(omega * (angleOf(ox, oy) - angleOf(pcaX, pcaY)), twopi);
    const double zOrigin = z0 + tanl * ds;
    const double rx = (ox - cx) / R, ry = (oy - cy) / R;
    const float seedPhi = static_cast<float>(std::atan2(omega * rx, -omega * ry));
    const int charge = -omega * d.bzSign > 0 ? 1 : -1;

    // chain points and hits per layer
    std::map<int, Acts::Vector3> chainPoint;  // official frame [cm]
    for (unsigned int k = 0; k < traj->size_points(); ++k)
    {
      const int layer = traj->get_point_layer(k);
      if (layer >= 0 && layer < 7 && !chainPoint.contains(layer))
      {
        chainPoint[layer] = Acts::Vector3(traj->get_point_x(k) + m_siBeam[0], traj->get_point_y(k) + m_siBeam[1], traj->get_point_z(k) - m_siActsZShift);
      }
    }
    std::map<int, std::vector<Si_Trajectory::HitIndex>> hitsByLayer;
    for (const auto& hit : traj->get_hit_indices())
    {
      const int layer = TrkrDefs::getLayer(hit.first);
      if (layer < 7)
      {
        hitsByLayer[layer].push_back(hit);
      }
    }

    struct NewCluster
    {
      TrkrDefs::cluskey key = TrkrDefs::CLUSKEYMAX;
      TrkrClusterv5 cluster;
    };
    std::vector<TrkrDefs::cluskey> keys;
    std::vector<NewCluster> newClusters;
    unsigned int nmvtx = 0, nintt = 0;
    unsigned long nAssoc = 0, nPosition = 0, nOther = 0, nProj = 0, nMiss = 0;
    for (const auto& [layer, hits] : hitsByLayer)
    {
      const bool isIntt = layer >= 3;
      std::map<TrkrDefs::hitsetkey, int> nPerSensor;
      for (const auto& hit : hits)
      {
        ++nPerSensor[hit.first];
      }
      const auto sensor = std::max_element(nPerSensor.begin(), nPerSensor.end(), [](const auto& a, const auto& b)
                                           { return a.second < b.second; })
                              ->first;
      const auto cp = chainPoint.find(layer);

      // 1) official cluster containing most chain hits
      TrkrDefs::cluskey best = TrkrDefs::CLUSKEYMAX;
      bool bestOtherStrobe = false;
      {
        std::map<TrkrDefs::cluskey, std::pair<int, bool>> nPerCluster;  // count, found in the other strobe only
        for (const auto& hit : hits)
        {
          bool otherStrobe = false;
          const auto ckey = d.officialClusterOfHit(hit.first, hit.second, otherStrobe);
          if (ckey != TrkrDefs::CLUSKEYMAX)
          {
            auto& e = nPerCluster.try_emplace(ckey, 0, true).first->second;
            ++e.first;
            e.second = e.second && otherStrobe;
          }
        }
        int nbest = 0;
        for (const auto& [ckey, e] : nPerCluster)
        {
          if (e.first > nbest)
          {
            nbest = e.first;
            best = ckey;
            bestOtherStrobe = e.second;
          }
        }
      }
      // 2) closest official cluster of the majority sensor
      bool byPosition = false;
      if (best == TrkrDefs::CLUSKEYMAX && cp != chainPoint.end())
      {
        if (const auto surface = d.geometry->maps().getSiliconSurface(sensor))
        {
          const Acts::Vector3 lt = surface->localToGlobalTransform(gctx).inverse() * (cp->second * cm);
          const auto& win = isIntt ? kInttMatchWindow : kMvtxMatchWindow;
          double bestDist = std::numeric_limits<double>::max();
          const auto range = d.clusters->getClusters(sensor);
          for (auto it = range.first; it != range.second; ++it)
          {
            const Acts::Vector2 lc = d.geometry->getLocalCoords(it->first, it->second) * cm;
            const double d0 = std::abs(lc.x() - lt.x()) / cm, d1 = std::abs(lc.y() - lt.y()) / cm;
            if (d0 > win[0] || d1 > win[1])
            {
              continue;
            }
            const double dist = square(d0 / win[0]) + square(d1 / win[1]);
            if (dist < bestDist)
            {
              bestDist = dist;
              best = it->first;
            }
          }
          byPosition = best != TrkrDefs::CLUSKEYMAX;
        }
      }
      if (best != TrkrDefs::CLUSKEYMAX)
      {
        keys.push_back(best);
        ++(byPosition ? nPosition : nAssoc);
        nOther += bestOtherStrobe;
        ++(isIntt ? nintt : nmvtx);
        continue;
      }
      // 3) optional: new cluster at the chain point projected on the official surface of its sensor
      if (m_useProjected && cp != chainPoint.end())
      {
        if (const auto surface = d.geometry->maps().getSiliconSurface(sensor))
        {
          const Acts::Vector3 lt = surface->localToGlobalTransform(gctx).inverse() * (cp->second * cm);
          NewCluster nc;
          nc.key = TrkrDefs::genClusKey(sensor, d.freeIndex(sensor));
          nc.cluster.setLocalX(static_cast<float>(lt.x() / cm));
          nc.cluster.setLocalY(static_cast<float>(lt.y() / cm));
          const auto& err = isIntt ? kInttChainErr : kMvtxChainErr;
          nc.cluster.setPhiError(static_cast<float>(err[0]));
          nc.cluster.setZError(static_cast<float>(err[1]));
          nc.cluster.setPhiSize(2);
          nc.cluster.setZSize(2);
          nc.cluster.setAdc(1);
          keys.push_back(nc.key);
          newClusters.push_back(std::move(nc));
          ++nProj;
          ++(isIntt ? nintt : nmvtx);
          continue;
        }
      }
      ++nMiss;
    }
    d.nSiMissing += nMiss;
    if (nmvtx < m_minMvtx || nmvtx + nintt < m_minSi)
    {
      reject("Si seed with too few clusters");
      continue;
    }
    d.nSiOfficialAssoc += nAssoc;
    d.nSiOfficialPosition += nPosition;
    d.nSiOtherStrobe += nOther;
    d.nSiProjected += nProj;
    for (const auto& nc : newClusters)
    {
      d.clusters->addClusterSpecifyKey(nc.key, new TrkrClusterv5(nc.cluster));
    }

    TrackSeed_v2 seed;
    seed.set_X0(static_cast<float>(cx));
    seed.set_Y0(static_cast<float>(cy));
    seed.set_qOverR(static_cast<float>(charge / R));  // TrackSeed_v2::get_charge: sign of q/R
    seed.set_Z0(static_cast<float>(zOrigin));
    seed.set_slope(static_cast<float>(tanl));
    seed.set_phi(seedPhi);
    const auto xing = crossingOfTrajectory.find(itraj);
    seed.set_crossing(xing == crossingOfTrajectory.end() ? static_cast<short>(SHRT_MAX) : xing->second);
    for (const auto key : keys)
    {
      seed.insert_cluster_key(key);
    }
    d.siSeeds->insert(&seed);
    siSeedOfTrajectory[itraj] = static_cast<unsigned int>(d.siSeeds->size() - 1);
    ++d.nSiSeeds;
  }

  // ---- TPC seed of a TPC poly track: TpcPolyTrackSeedConverter copies the poly-track cluster keys
  std::unordered_map<TrkrDefs::cluskey, unsigned int> tpcSeedOfKey;
  for (unsigned int i = 0; i < d.tpcSeeds->size(); ++i)
  {
    if (const TrackSeed* s = d.tpcSeeds->get(i))
    {
      for (auto it = s->begin_cluster_keys(); it != s->end_cluster_keys(); ++it)
      {
        if (TrkrDefs::getTrkrId(*it) == TrkrDefs::tpcId)
        {
          tpcSeedOfKey.emplace(*it, i);
        }
      }
    }
  }

  // ---- SvtxTrackSeeds: one per own Si+TPC track
  std::set<unsigned int> matchedTpcSeeds;
  for (unsigned int i = 0; i < d.siTpcTracks->size(); ++i)
  {
    const SiTpc_Track* t = d.siTpcTracks->get_track(i);
    if (!t || !t->isValid())
    {
      continue;
    }
    ++d.nSiTpc;
    const auto si = siSeedOfTrajectory.find(t->get_si_trajectory_index());
    if (si == siSeedOfTrajectory.end())
    {
      reject("Si+TPC track without Si seed");
      continue;
    }
    const Tpc_PolyTrack* poly = t->get_tpc_track_index() < d.tpcPolyTracks->size() ? d.tpcPolyTracks->get_track(t->get_tpc_track_index()) : nullptr;
    if (!poly)
    {
      reject("Si+TPC track without TPC poly track");
      continue;
    }
    std::map<unsigned int, int> votes;
    for (unsigned int k = 0; k < poly->size_cluster_keys(); ++k)
    {
      const auto f = tpcSeedOfKey.find(poly->get_cluster_key(k));
      if (f != tpcSeedOfKey.end())
      {
        ++votes[f->second];
      }
    }
    if (votes.empty())
    {
      reject("Si+TPC track without TPC seed (TpcPolyTrackSeedConverter cuts: > 20 clusters, pT > 0.1)");
      continue;
    }
    const unsigned int tpcIndex = std::max_element(votes.begin(), votes.end(), [](const auto& a, const auto& b)
                                                   { return a.second < b.second; })
                                      ->first;
    const short crossing = crossingOfAssembled(t->get_tpc_assembled_track_id());
    if (crossing == SHRT_MAX)
    {
      reject("Si+TPC track without crossing");
      continue;
    }
    SvtxTrackSeed_v3 svtx;
    svtx.set_silicon_seed_index(si->second);
    svtx.set_tpc_seed_index(tpcIndex);
    svtx.set_crossing(crossing);
    svtx.set_crossing_estimate(crossing);
    d.svtxSeeds->insert(&svtx);
    matchedTpcSeeds.insert(tpcIndex);
    ++d.nSvtxSeeds;
  }

  // ---- optional TPC-only entries (as PHSiliconTpcTrackMatching)
  if (m_addUnmatchedTpc)
  {
    for (unsigned int i = 0; i < d.tpcSeeds->size(); ++i)
    {
      if (d.tpcSeeds->get(i) && !matchedTpcSeeds.contains(i))
      {
        SvtxTrackSeed_v3 svtx;
        svtx.set_tpc_seed_index(i);
        d.svtxSeeds->insert(&svtx);
        ++d.nTpcOnly;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcTrackSeedFiller::End(PHCompositeNode* /*topNode*/)
{
  auto& d = *m_impl;
  std::cout << Name() << "::End - events " << d.nEvents << std::endl
            << "  Si trajectories " << d.nTrajectories << " -> Si seeds " << d.nSiSeeds << std::endl
            << "  Si clusters: official by hit association " << d.nSiOfficialAssoc << " (pruned MVTX duplicate, other strobe "
            << d.nSiOtherStrobe << "), by position " << d.nSiOfficialPosition << ", projected chain points " << d.nSiProjected
            << ", layers without cluster " << d.nSiMissing << std::endl
            << "  Si+TPC tracks " << d.nSiTpc << " -> SvtxTrackSeeds " << d.nSvtxSeeds << ", TPC-only SvtxTrackSeeds " << d.nTpcOnly << std::endl;
  for (const auto& [why, n] : d.rejected)
  {
    std::cout << "  rejected (" << why << "): " << n << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
