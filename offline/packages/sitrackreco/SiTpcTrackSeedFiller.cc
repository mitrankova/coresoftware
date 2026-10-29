#include "SiTpcTrackSeedFiller.h"

#include "SiHitSeedData.h"
#include "SiTpcHelixFit.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>

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
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  constexpr double cm = Acts::UnitConstants::cm;  // 10 (Acts length unit = mm)
  constexpr double twopi = 2. * M_PI;

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

  // detector index of the fit weights: 0 = MVTX, 1 = INTT, 2 = TPC
  int detIndex(unsigned int layer)
  {
    return layer <= 2 ? 0 : (layer <= 6 ? 1 : 2);
  }

  // key of a poly cluster in TRKR_CLUSTER, as TpcPolyClusterTrkrClusterConverter::getClusterKey
  TrkrDefs::cluskey polyClusterKey(const Tpc_PolyCluster* cluster)
  {
    if (!cluster || cluster->size_hits() == 0)
    {
      return TrkrDefs::CLUSKEYMAX;
    }
    const TrkrDefs::hitsetkey hitsetkey = cluster->get_hit_index(0).first;
    const TrkrDefs::cluskey key = cluster->get_trkr_cluster_key();
    if (key != TrkrDefs::CLUSKEYMAX && TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId && TrkrDefs::getHitSetKeyFromClusKey(key) == hitsetkey)
    {
      return key;
    }
    return TrkrDefs::genClusKey(hitsetkey, cluster->get_cluster_id());
  }
}  // namespace

//____________________________________________________________________________________________________
struct SiTpcTrackSeedFiller::Impl
{
  ActsGeometry* geometry = nullptr;
  SiHitSeedEvent* siSeedEvent = nullptr;
  Tpc_PolyClusterContainer* tpcClusters = nullptr;
  TpcCrossingDecisionContainer* crossings = nullptr;
  TrkrClusterContainer* clusters = nullptr;
  TrkrClusterHitAssoc* hitAssoc = nullptr;
  TrackSeedContainer* tpcSeeds = nullptr;
  TrackSeedContainer* siSeeds = nullptr;
  TrackSeedContainer* svtxSeeds = nullptr;
  double bzSign = 1;

  std::unordered_map<TrkrDefs::hitsetkey, std::unordered_map<TrkrDefs::hitkey, TrkrDefs::cluskey>> hitToCluster;
  std::map<TrkrDefs::hitsetkey, uint32_t> nextIndex;

  unsigned long nEvents = 0, nChains = 0, nTpcSeeds = 0, nSeeds = 0, nHybrid = 0, nTpcOnly = 0;
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
  d.tpcClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
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
  require(d.tpcClusters, m_tpcClusterNodeName);
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

  std::cout << Name() << "::InitRun - chains " << (m_selectedOnly ? "selected" : "all (best rank per TPC track)")
            << ", Si beam spot (" << m_siBeam[0] << ", " << m_siBeam[1] << ") cm, Si acts z shift " << m_siActsZShift << " cm" << std::endl
            << "  seed fit: xy weights MVTX/INTT/TPC " << m_wxy[0] << "/" << m_wxy[1] << "/" << m_wxy[2] << ", z " << m_wz[0] << "/"
            << m_wz[1] << "/" << m_wz[2] << ", hybrid " << (m_hybrid ? "on" : "off") << " (" << m_hybridW[0] << "/" << m_hybridW[1]
            << ", min " << m_hybridMinPoints << "), TPC z offset " << (m_fitTpcZOffset ? "free" : "fixed") << std::endl
            << "  projected chain points " << (m_useProjected ? "on" : "off") << ", TPC-only seeds " << (m_addUnmatchedTpc ? "on" : "off")
            << ", output " << m_siSeedName << ", " << m_svtxSeedName << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcTrackSeedFiller::process_event(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
  ++d.nEvents;
  d.siSeeds->Reset();
  d.svtxSeeds->Reset();
  d.hitToCluster.clear();
  d.nextIndex.clear();

  d.siSeedEvent = findNode::getClass<SiHitSeedEvent>(topNode, m_siSeedEventName);
  if (!d.siSeedEvent)
  {
    std::cout << PHWHERE << " " << Name() << ": missing " << m_siSeedEventName << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  const SiHitSeedEvent& event = *d.siSeedEvent;
  const auto gctx = d.geometry->geometry().getGeoContext();
  auto reject = [&](const std::string& why)
  {
    ++d.rejected[why];
  };

  // ---- TPC seed of an assembled TPC track: the poly-cluster keys of the seed
  std::unordered_map<TrkrDefs::cluskey, unsigned int> assembledOfKey;
  for (unsigned int i = 0; i < d.tpcClusters->size(); ++i)
  {
    const Tpc_PolyCluster* c = d.tpcClusters->get_cluster(i);
    if (c && c->isValid())
    {
      assembledOfKey.emplace(polyClusterKey(c), c->get_source_assembled_track_id());
    }
  }
  std::map<unsigned int, unsigned int> tpcSeedOfAssembled;
  for (unsigned int i = 0; i < d.tpcSeeds->size(); ++i)
  {
    const TrackSeed* s = d.tpcSeeds->get(i);
    if (!s)
    {
      continue;
    }
    std::map<unsigned int, int> votes;
    for (auto it = s->begin_cluster_keys(); it != s->end_cluster_keys(); ++it)
    {
      if (const auto f = assembledOfKey.find(*it); f != assembledOfKey.end())
      {
        ++votes[f->second];
      }
    }
    if (!votes.empty())
    {
      tpcSeedOfAssembled[std::max_element(votes.begin(), votes.end(), [](const auto& a, const auto& b)
                                          { return a.second < b.second; })
                             ->first] = i;
      ++d.nTpcSeeds;
    }
  }

  // ---- chains -> TPC seeds (one chain per TPC seed: selected, or best rank / most layers)
  std::map<unsigned int, const SiHitChain*> chainOfTpcSeed;
  for (const auto& chain : event.chains)
  {
    if (m_selectedOnly && !chain.selected)
    {
      continue;
    }
    ++d.nChains;
    const auto t = tpcSeedOfAssembled.find(chain.tpc_assembled_track_id);
    if (chain.tpc_track_index < 0 || t == tpcSeedOfAssembled.end())
    {
      reject("chain without TPC seed (TpcPolyTrackSeedConverter cuts: > 20 clusters, pT > 0.1)");
      continue;
    }
    auto& best = chainOfTpcSeed[t->second];
    if (!best || chain.rank < best->rank || (chain.rank == best->rank && chain.n_layers() > best->n_layers()))
    {
      best = &chain;
    }
  }

  std::set<unsigned int> matchedTpcSeeds;
  for (const auto& [tpcIndex, chainPtr] : chainOfTpcSeed)
  {
    const SiHitChain& chain = *chainPtr;
    const TrackSeed* tpcSeed = d.tpcSeeds->get(tpcIndex);
    const TpcCrossingDecision* decision = d.crossings->get_decision(chain.tpc_assembled_track_id);
    const short crossing = decision ? decision->get_selected_crossing() : static_cast<short>(SHRT_MAX);
    if (crossing == SHRT_MAX)
    {
      reject("chain without crossing");
      continue;
    }

    // ---- Si clusters, one per chain layer
    struct NewCluster
    {
      TrkrDefs::cluskey key = TrkrDefs::CLUSKEYMAX;
      TrkrClusterv5 cluster;
    };
    std::vector<TrkrDefs::cluskey> siKeys;
    std::vector<NewCluster> newClusters;
    unsigned int nmvtx = 0, nintt = 0;
    unsigned long nAssoc = 0, nPosition = 0, nOther = 0, nProj = 0, nMiss = 0;
    for (std::size_t il = 0; il < chain.layers.size(); ++il)
    {
      const int layer = chain.layers[il];
      const bool isIntt = layer >= 3;
      std::vector<const SiHitPoint*> hits;
      if (il < chain.blob_ids.size() && chain.blob_ids[il] >= 0 && chain.blob_ids[il] < static_cast<int>(event.blobs.size()))
      {
        for (const int hid : event.blobs[chain.blob_ids[il]].hit_ids)
        {
          if (hid >= 0 && hid < static_cast<int>(event.hits.size()))
          {
            hits.push_back(&event.hits[hid]);
          }
        }
      }
      if (hits.empty())
      {
        ++nMiss;
        continue;
      }
      std::map<TrkrDefs::hitsetkey, int> nPerSensor;
      for (const auto* h : hits)
      {
        ++nPerSensor[h->hitsetkey];
      }
      const auto sensor = std::max_element(nPerSensor.begin(), nPerSensor.end(), [](const auto& a, const auto& b)
                                           { return a.second < b.second; })
                              ->first;
      // chain point in the official frame
      const bool hasPoint = il < chain.phys_x.size() && il < chain.phys_y.size() && il < chain.phys_z.size();
      const Acts::Vector3 target = hasPoint ? Acts::Vector3(chain.phys_x[il] + m_siBeam[0], chain.phys_y[il] + m_siBeam[1], chain.phys_z[il] - m_siActsZShift)
                                            : Acts::Vector3::Zero();

      // 1) official cluster containing most hits of the chain blob
      TrkrDefs::cluskey best = TrkrDefs::CLUSKEYMAX;
      bool bestOtherStrobe = false;
      {
        std::map<TrkrDefs::cluskey, std::pair<int, bool>> nPerCluster;  // count, found in the other strobe only
        for (const auto* h : hits)
        {
          bool otherStrobe = false;
          const auto ckey = d.officialClusterOfHit(h->hitsetkey, h->hitkey, otherStrobe);
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
      if (best == TrkrDefs::CLUSKEYMAX && hasPoint)
      {
        if (const auto surface = d.geometry->maps().getSiliconSurface(sensor))
        {
          const Acts::Vector3 lt = surface->localToGlobalTransform(gctx).inverse() * (target * cm);
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
        siKeys.push_back(best);
        ++(byPosition ? nPosition : nAssoc);
        nOther += bestOtherStrobe;
        ++(isIntt ? nintt : nmvtx);
        continue;
      }
      // 3) optional: new cluster at the chain point projected on the official surface of its sensor
      if (m_useProjected && hasPoint)
      {
        if (const auto surface = d.geometry->maps().getSiliconSurface(sensor))
        {
          const Acts::Vector3 lt = surface->localToGlobalTransform(gctx).inverse() * (target * cm);
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
          siKeys.push_back(nc.key);
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
      reject("too few Si clusters");
      continue;
    }

    // ---- seed fit in the official frame: Si clusters + TPC clusters of the TPC seed
    for (const auto& nc : newClusters)
    {
      d.clusters->addClusterSpecifyKey(nc.key, new TrkrClusterv5(nc.cluster));
    }
    std::vector<std::pair<double, SiTpcHelixFit::Point>> pts;  // radius, point
    for (const auto key : siKeys)
    {
      TrkrCluster* c = d.clusters->findCluster(key);
      if (!c)
      {
        continue;
      }
      const Acts::Vector3 g = d.geometry->getGlobalPosition(key, c);
      const int det = detIndex(TrkrDefs::getLayer(key));
      SiTpcHelixFit::Point p{g.x(), g.y(), g.z(), m_wxy[det], m_wz[det], 0, m_hybrid ? m_hybridW[det] : 0.0};
      pts.emplace_back(std::hypot(g.x(), g.y()), p);
    }
    for (auto it = tpcSeed->begin_cluster_keys(); it != tpcSeed->end_cluster_keys(); ++it)
    {
      if (TrkrDefs::getTrkrId(*it) != TrkrDefs::tpcId)
      {
        continue;
      }
      TrkrCluster* c = d.clusters->findCluster(*it);
      if (!c)
      {
        continue;
      }
      const Acts::Vector3 g = d.geometry->getGlobalPosition(*it, c, crossing);
      if (!g.allFinite())
      {
        continue;
      }
      SiTpcHelixFit::Point p{g.x(), g.y(), g.z(), m_wxy[2], m_wz[2], 1, 0.0};
      pts.emplace_back(std::hypot(g.x(), g.y()), p);
    }
    std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b)
              { return a.first < b.first; });
    std::vector<SiTpcHelixFit::Point> fitPoints;
    fitPoints.reserve(pts.size());
    for (const auto& [r, p] : pts)
    {
      fitPoints.push_back(p);
    }
    SiTpcHelixFit::Config cfg;
    cfg.beamX = m_siBeam[0];
    cfg.beamY = m_siBeam[1];
    cfg.bz = 1.4 * d.bzSign;
    cfg.groupZOffset = m_fitTpcZOffset;
    cfg.hybridMinPoints = m_hybrid ? m_hybridMinPoints : 0;
    const auto r = SiTpcHelixFit::fit(fitPoints, cfg);
    if (r.status != SiTpcHelixFit::Ok || !(r.R > 0) || r.charge == 0)
    {
      reject("seed fit failed");
      continue;
    }
    d.nHybrid += r.hybrid;

    // z and direction at the point of the circle closest to the origin (TrackSeedHelper::get_xyz)
    const double r0 = std::hypot(r.cx, r.cy);
    if (!(r0 > 0))
    {
      reject("seed fit failed");
      continue;
    }
    const double h = r.helicity >= 0 ? 1. : -1.;  // +1: counter-clockwise going outward
    auto angleOf = [&](double x, double y)
    { return std::atan2(y - r.cy, x - r.cx); };
    const double ox = r.cx - r.R * r.cx / r0, oy = r.cy - r.R * r.cy / r0;
    const double ds = r.R * std::remainder(h * (angleOf(ox, oy) - angleOf(r.pcaX, r.pcaY)), twopi);
    const double zOrigin = r.z0 + r.tanl * ds;
    const double rx = (ox - r.cx) / r.R, ry = (oy - r.cy) / r.R;
    const float seedPhi = static_cast<float>(std::atan2(h * rx, -h * ry));

    TrackSeed_v2 seed;
    seed.set_X0(static_cast<float>(r.cx));
    seed.set_Y0(static_cast<float>(r.cy));
    seed.set_qOverR(static_cast<float>(r.charge / r.R));  // TrackSeed_v2::get_charge: sign of q/R
    seed.set_Z0(static_cast<float>(zOrigin));
    seed.set_slope(static_cast<float>(r.tanl));
    seed.set_phi(seedPhi);
    seed.set_crossing(crossing);
    for (const auto key : siKeys)
    {
      seed.insert_cluster_key(key);
    }
    d.siSeeds->insert(&seed);

    SvtxTrackSeed_v3 svtx;
    svtx.set_silicon_seed_index(static_cast<unsigned int>(d.siSeeds->size() - 1));
    svtx.set_tpc_seed_index(tpcIndex);
    svtx.set_crossing(crossing);
    svtx.set_crossing_estimate(crossing);
    d.svtxSeeds->insert(&svtx);
    matchedTpcSeeds.insert(tpcIndex);
    ++d.nSeeds;
    d.nSiOfficialAssoc += nAssoc;
    d.nSiOfficialPosition += nPosition;
    d.nSiOtherStrobe += nOther;
    d.nSiProjected += nProj;
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
  std::cout << Name() << "::End - events " << d.nEvents << ", chains " << d.nChains << ", TPC seeds " << d.nTpcSeeds
            << " -> Si seeds = Si+TPC SvtxTrackSeeds " << d.nSeeds << " (hybrid xy " << d.nHybrid << "), TPC-only SvtxTrackSeeds "
            << d.nTpcOnly << std::endl
            << "  Si clusters: official by hit association " << d.nSiOfficialAssoc << " (pruned MVTX duplicate, other strobe "
            << d.nSiOtherStrobe << "), by position " << d.nSiOfficialPosition << ", projected chain points " << d.nSiProjected
            << ", chain layers without cluster " << d.nSiMissing << std::endl;
  for (const auto& [why, n] : d.rejected)
  {
    std::cout << "  rejected (" << why << "): " << n << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
