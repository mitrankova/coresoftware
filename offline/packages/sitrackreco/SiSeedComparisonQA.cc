#include "SiSeedComparisonQA.h"

#include "SiTpcHelixFit.h"
#include "SiTpc_Track.h"
#include "SiTpc_TrackContainer.h"
#include "Si_Trajectory.h"
#include "Si_TrajectoryContainer.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

namespace
{
  constexpr float fnan = std::numeric_limits<float>::quiet_NaN();

  int detOfLayer(int layer) { return layer <= 2 ? 0 : (layer <= 6 ? 1 : 2); }

  // own Si chain, per layer: official cluster with most shared hits
  struct LayerMatch
  {
    int layer = -1;
    int nChainHits = 0;
    int nSensors = 0;
    int nClustersTouched = 0;
    int nShared = 0;
    int nClusterHits = 0;
    TrkrDefs::cluskey key = TrkrDefs::CLUSKEYMAX;
    TrkrCluster* cluster = nullptr;
    bool hasCluster = false;
    bool otherStrobe = false;  // the cluster copy of the neighbouring MVTX strobe (pruned duplicate)
    double cx = 0, cy = 0, cz = 0;  // cluster, own frame [cm]
  };

  // parameters of a track for the comparison (official frame)
  struct Params
  {
    float pt = fnan, phi = fnan, eta = fnan, d0 = fnan, z0 = fnan;
    int charge = 0, crossing = SHRT_MAX, nmvtx = 0, nintt = 0;
  };

  struct OwnChain
  {
    const Si_Trajectory* traj = nullptr;
    unsigned int index = 0;
    std::map<int, LayerMatch> layers;
    std::set<TrkrDefs::cluskey> clusters;
    Params par;
    int siTpcIndex = -1;  // index in SITPC_TRACKS
    int resolverCrossing = SHRT_MAX;
  };

  struct StdSeed
  {
    const TrackSeed* seed = nullptr;
    std::set<TrkrDefs::cluskey> clusters;
    Params par;
  };

  float wrap(double a) { return static_cast<float>(std::remainder(a, 2 * M_PI)); }
}  // namespace

//____________________________________________________________________________________________________
struct SiSeedComparisonQA::Impl
{
  TFile* file = nullptr;
  TTree* points = nullptr;
  TTree* bias = nullptr;
  TTree* ownChains = nullptr;
  TTree* stdSeeds = nullptr;
  TTree* vertices = nullptr;
  int event = 0;

  struct PointRow
  {
    int event = 0, traj = 0, chain = 0, layer = 0, in_sitpc = 0;
    int nchain = 0, nsensors = 0, nclusters = 0, nshared = 0, ncluster_hits = 0, has_cluster = 0, other_strobe = 0;
    int phisize = 0, zsize = 0;
    float px = fnan, py = fnan, pz = fnan, cx = fnan, cy = fnan, cz = fnan;
    float d_rphi = fnan, d_z = fnan, d_r = fnan, err_rphi = fnan, err_z = fnan, pt = fnan;
  } p;
  struct BiasRow
  {
    int event = 0, traj = 0, in_sitpc = 0, nsi = 0, nreplaced = 0, si_ok = 0, st_ok = 0;
    float pt = fnan, eta = fnan;
    // Si-only refit: clusters - own points; ref = own-point refit - stored
    float si_d0 = fnan, si_phi = fnan, si_z0 = fnan, si_tanl = fnan, si_relpt = fnan;
    float si_ref_d0 = fnan, si_ref_phi = fnan, si_ref_z0 = fnan;
    // Si+TPC (hybrid) refit
    float st_d0 = fnan, st_phi = fnan, st_z0 = fnan, st_tanl = fnan, st_relpt = fnan;
    float st_ref_d0 = fnan, st_ref_phi = fnan, st_ref_z0 = fnan;
  } b;
  struct ChainRow
  {
    int event = 0, traj = 0, in_sitpc = 0, nown = 0, nstd = 0, nshared = 0, std_found = 0;
    int own_crossing = 0, resolver_crossing = 0, std_crossing = 0, own_nmvtx = 0, own_nintt = 0, std_nmvtx = 0, std_nintt = 0;
    float own_pt = fnan, own_phi = fnan, own_eta = fnan, own_d0 = fnan, own_z0 = fnan;
    float std_pt = fnan, std_phi = fnan, std_eta = fnan, std_d0 = fnan, std_z0 = fnan;
  } c;
  struct SeedRow
  {
    int event = 0, seed = 0, nstd = 0, nown = 0, nshared = 0, own_found = 0, own_in_sitpc = 0, crossing = 0, nmvtx = 0, nintt = 0;
    float pt = fnan, phi = fnan, eta = fnan, d0 = fnan, z0 = fnan;
    float near_sitpc_dphi = fnan, near_sitpc_deta = fnan, near_sitpc_pt = fnan;
  } s;
  struct VertexRow
  {
    int event = 0, type = 0, crossing = 0, own_ntrk = 0, std_ntrk = 0, matched = 0;
    float own_x = fnan, own_y = fnan, own_z = fnan, std_x = fnan, std_y = fnan, std_z = fnan;
  } v;

  // summary counters
  std::map<int, std::array<double, 6>> layerSums;  // n, sum drphi, sum drphi^2, sum dz, sum dz^2, n identical
  std::map<int, std::array<long, 6>> layerCounts;   // points, no cluster, identical, split, merged, other strobe
  long nChains = 0, nChainsAll = 0, nChainsPartial = 0, nChainsNone = 0;
  long nChainsSiTpc = 0, nChainsSiTpcAll = 0, nChainsSiTpcPartial = 0, nChainsSiTpcNone = 0;
  long nStd = 0, nStdMatched = 0, nStdNearSiTpc = 0;
  long nOwnVtx = 0, nOwnVtxMatched = 0, nStdVtx = 0;

  std::unordered_map<TrkrDefs::hitsetkey, std::unordered_map<TrkrDefs::hitkey, TrkrDefs::cluskey>> hitToCluster;
  TrkrDefs::cluskey clusterOfHit(TrkrClusterHitAssoc* assoc, TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk)
  {
    if (!assoc)
    {
      return TrkrDefs::CLUSKEYMAX;
    }
    auto it = hitToCluster.find(hsk);
    if (it == hitToCluster.end())
    {
      std::unordered_map<TrkrDefs::hitkey, TrkrDefs::cluskey> m;
      if (auto* map = assoc->getClusterMap(hsk))
      {
        for (const auto& [ckey, hitkey] : *map)
        {
          m.emplace(hitkey, ckey);
        }
      }
      it = hitToCluster.emplace(hsk, std::move(m)).first;
    }
    const auto f = it->second.find(hk);
    return f == it->second.end() ? TrkrDefs::CLUSKEYMAX : f->second;
  }

  // official cluster of a chain hit.  MVTX: MvtxClusterPruner removes the copy of a cluster duplicated in the
  // next strobe (identical pixels) from TRKR_CLUSTER without updating TRKR_CLUSTERHITASSOC, so a chain hit
  // can point to a removed cluster: then the cluster with the same pixel in the neighbouring strobe is used.
  TrkrDefs::cluskey officialClusterOfHit(TrkrClusterHitAssoc* assoc, TrkrClusterContainer* clusters, TrkrDefs::hitsetkey hsk,
                                         TrkrDefs::hitkey hk, bool& otherStrobe)
  {
    otherStrobe = false;
    if (!clusters)
    {
      return TrkrDefs::CLUSKEYMAX;
    }
    const auto ckey = clusterOfHit(assoc, hsk, hk);
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
      const auto other = clusterOfHit(assoc, MvtxDefs::genHitSetKey(layer, stave, chip, strobe + ds), hk);
      if (other != TrkrDefs::CLUSKEYMAX && clusters->findCluster(other))
      {
        otherStrobe = true;
        return other;
      }
    }
    return TrkrDefs::CLUSKEYMAX;
  }
};

//____________________________________________________________________________________________________
SiSeedComparisonQA::SiSeedComparisonQA(const std::string& name, const std::string& fileName)
  : SubsysReco(name)
  , m_impl(std::make_unique<Impl>())
  , m_fileName(fileName)
{
}

SiSeedComparisonQA::~SiSeedComparisonQA() = default;

//____________________________________________________________________________________________________
int SiSeedComparisonQA::InitRun(PHCompositeNode* /*topNode*/)
{
  auto& d = *m_impl;
  d.file = new TFile(m_fileName.c_str(), "RECREATE");

  d.points = new TTree("si_points", "own chain point vs official cluster (own frame, cm)");
  auto& p = d.p;
  d.points->Branch("event", &p.event);
  d.points->Branch("traj", &p.traj);
  d.points->Branch("chain", &p.chain);
  d.points->Branch("layer", &p.layer);
  d.points->Branch("in_sitpc", &p.in_sitpc);
  d.points->Branch("nchain", &p.nchain);
  d.points->Branch("nsensors", &p.nsensors);
  d.points->Branch("nclusters", &p.nclusters);
  d.points->Branch("nshared", &p.nshared);
  d.points->Branch("ncluster_hits", &p.ncluster_hits);
  d.points->Branch("has_cluster", &p.has_cluster);
  d.points->Branch("other_strobe", &p.other_strobe);
  d.points->Branch("phisize", &p.phisize);
  d.points->Branch("zsize", &p.zsize);
  d.points->Branch("px", &p.px);
  d.points->Branch("py", &p.py);
  d.points->Branch("pz", &p.pz);
  d.points->Branch("cx", &p.cx);
  d.points->Branch("cy", &p.cy);
  d.points->Branch("cz", &p.cz);
  d.points->Branch("d_rphi", &p.d_rphi);
  d.points->Branch("d_z", &p.d_z);
  d.points->Branch("d_r", &p.d_r);
  d.points->Branch("err_rphi", &p.err_rphi);
  d.points->Branch("err_z", &p.err_z);
  d.points->Branch("pt", &p.pt);

  d.bias = new TTree("seed_bias", "refit with official clusters - refit with own points (own frame, cm, rad)");
  auto& b = d.b;
  d.bias->Branch("event", &b.event);
  d.bias->Branch("traj", &b.traj);
  d.bias->Branch("in_sitpc", &b.in_sitpc);
  d.bias->Branch("nsi", &b.nsi);
  d.bias->Branch("nreplaced", &b.nreplaced);
  d.bias->Branch("si_ok", &b.si_ok);
  d.bias->Branch("st_ok", &b.st_ok);
  d.bias->Branch("pt", &b.pt);
  d.bias->Branch("eta", &b.eta);
  d.bias->Branch("si_d0", &b.si_d0);
  d.bias->Branch("si_phi", &b.si_phi);
  d.bias->Branch("si_z0", &b.si_z0);
  d.bias->Branch("si_tanl", &b.si_tanl);
  d.bias->Branch("si_relpt", &b.si_relpt);
  d.bias->Branch("si_ref_d0", &b.si_ref_d0);
  d.bias->Branch("si_ref_phi", &b.si_ref_phi);
  d.bias->Branch("si_ref_z0", &b.si_ref_z0);
  d.bias->Branch("st_d0", &b.st_d0);
  d.bias->Branch("st_phi", &b.st_phi);
  d.bias->Branch("st_z0", &b.st_z0);
  d.bias->Branch("st_tanl", &b.st_tanl);
  d.bias->Branch("st_relpt", &b.st_relpt);
  d.bias->Branch("st_ref_d0", &b.st_ref_d0);
  d.bias->Branch("st_ref_phi", &b.st_ref_phi);
  d.bias->Branch("st_ref_z0", &b.st_ref_z0);

  d.ownChains = new TTree("own_chains", "own Si chain -> standard seed with most shared clusters (official frame)");
  auto& c = d.c;
  d.ownChains->Branch("event", &c.event);
  d.ownChains->Branch("traj", &c.traj);
  d.ownChains->Branch("in_sitpc", &c.in_sitpc);
  d.ownChains->Branch("nown", &c.nown);
  d.ownChains->Branch("nstd", &c.nstd);
  d.ownChains->Branch("nshared", &c.nshared);
  d.ownChains->Branch("std_found", &c.std_found);
  d.ownChains->Branch("own_crossing", &c.own_crossing);
  d.ownChains->Branch("resolver_crossing", &c.resolver_crossing);
  d.ownChains->Branch("std_crossing", &c.std_crossing);
  d.ownChains->Branch("own_nmvtx", &c.own_nmvtx);
  d.ownChains->Branch("own_nintt", &c.own_nintt);
  d.ownChains->Branch("std_nmvtx", &c.std_nmvtx);
  d.ownChains->Branch("std_nintt", &c.std_nintt);
  d.ownChains->Branch("own_pt", &c.own_pt);
  d.ownChains->Branch("own_phi", &c.own_phi);
  d.ownChains->Branch("own_eta", &c.own_eta);
  d.ownChains->Branch("own_d0", &c.own_d0);
  d.ownChains->Branch("own_z0", &c.own_z0);
  d.ownChains->Branch("std_pt", &c.std_pt);
  d.ownChains->Branch("std_phi", &c.std_phi);
  d.ownChains->Branch("std_eta", &c.std_eta);
  d.ownChains->Branch("std_d0", &c.std_d0);
  d.ownChains->Branch("std_z0", &c.std_z0);

  d.stdSeeds = new TTree("std_seeds", "standard silicon seed -> own chain with most shared clusters (official frame)");
  auto& s = d.s;
  d.stdSeeds->Branch("event", &s.event);
  d.stdSeeds->Branch("seed", &s.seed);
  d.stdSeeds->Branch("nstd", &s.nstd);
  d.stdSeeds->Branch("nown", &s.nown);
  d.stdSeeds->Branch("nshared", &s.nshared);
  d.stdSeeds->Branch("own_found", &s.own_found);
  d.stdSeeds->Branch("own_in_sitpc", &s.own_in_sitpc);
  d.stdSeeds->Branch("crossing", &s.crossing);
  d.stdSeeds->Branch("nmvtx", &s.nmvtx);
  d.stdSeeds->Branch("nintt", &s.nintt);
  d.stdSeeds->Branch("pt", &s.pt);
  d.stdSeeds->Branch("phi", &s.phi);
  d.stdSeeds->Branch("eta", &s.eta);
  d.stdSeeds->Branch("d0", &s.d0);
  d.stdSeeds->Branch("z0", &s.z0);
  d.stdSeeds->Branch("near_sitpc_dphi", &s.near_sitpc_dphi);
  d.stdSeeds->Branch("near_sitpc_deta", &s.near_sitpc_deta);
  d.stdSeeds->Branch("near_sitpc_pt", &s.near_sitpc_pt);

  d.vertices = new TTree("si_vertices", "own Si vertices vs standard silicon vertices (official frame); type 0 own, 1 standard only");
  auto& v = d.v;
  d.vertices->Branch("event", &v.event);
  d.vertices->Branch("type", &v.type);
  d.vertices->Branch("crossing", &v.crossing);
  d.vertices->Branch("own_ntrk", &v.own_ntrk);
  d.vertices->Branch("std_ntrk", &v.std_ntrk);
  d.vertices->Branch("matched", &v.matched);
  d.vertices->Branch("own_x", &v.own_x);
  d.vertices->Branch("own_y", &v.own_y);
  d.vertices->Branch("own_z", &v.own_z);
  d.vertices->Branch("std_x", &v.std_x);
  d.vertices->Branch("std_y", &v.std_y);
  d.vertices->Branch("std_z", &v.std_z);

  std::cout << Name() << ": own frame = official - (" << m_siBeam[0] << ", " << m_siBeam[1] << ") cm, z + " << m_zShift
            << " cm, INTT offset (" << m_inttOffset[0] << ", " << m_inttOffset[1] << ") cm; official beam (" << m_beam[0] << ", "
            << m_beam[1] << ") slopes (" << m_beam[2] << ", " << m_beam[3] << "); output " << m_fileName << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiSeedComparisonQA::process_event(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
  d.hitToCluster.clear();

  auto* trajs = findNode::getClass<Si_TrajectoryContainer>(topNode, m_trajNode);
  if (!trajs)
  {
    ++d.event;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  auto* siTpc = findNode::getClass<SiTpc_TrackContainer>(topNode, m_siTpcNode);
  auto* crossings = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingNode);
  auto* geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  auto* clusters = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  auto* assoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // own frame <-> official frame
  auto toOwn = [&](const Acts::Vector3& g, int layer)
  {
    const bool intt = layer >= 3;
    return Acts::Vector3(g.x() - m_siBeam[0] + (intt ? m_inttOffset[0] : 0.), g.y() - m_siBeam[1] + (intt ? m_inttOffset[1] : 0.),
                         g.z() + m_zShift);
  };
  // d0 to the beam line (official frame), Acts convention: (pca - beam) . (-sin phi, cos phi)
  auto d0Official = [&](double x, double y, double z, double phi)
  {
    const double bx = m_beam[0] + m_beam[2] * z, by = m_beam[1] + m_beam[3] * z;
    return static_cast<float>(-(x - bx) * std::sin(phi) + (y - by) * std::cos(phi));
  };

  // ---- SITPC track of every trajectory
  std::map<unsigned int, int> siTpcOfTraj;
  if (siTpc)
  {
    for (unsigned int i = 0; i < siTpc->size(); ++i)
    {
      if (const SiTpc_Track* t = siTpc->get_track(i))
      {
        siTpcOfTraj[t->get_si_trajectory_index()] = static_cast<int>(i);
      }
    }
  }

  // ---- own chains: per layer the official cluster with most shared hits
  std::vector<OwnChain> own;
  for (unsigned int it = 0; it < trajs->size(); ++it)
  {
    const Si_Trajectory* traj = trajs->get_trajectory(it);
    if (!traj)
    {
      continue;
    }
    OwnChain oc;
    oc.traj = traj;
    oc.index = it;
    const auto f = siTpcOfTraj.find(it);
    oc.siTpcIndex = f == siTpcOfTraj.end() ? -1 : f->second;

    std::map<int, std::vector<Si_Trajectory::HitIndex>> hitsByLayer;
    std::map<int, int> inttBuckets;
    for (const auto& hit : traj->get_hit_indices())
    {
      const int layer = TrkrDefs::getLayer(hit.first);
      if (layer < 7)
      {
        hitsByLayer[layer].push_back(hit);
      }
      if (TrkrDefs::getTrkrId(hit.first) == TrkrDefs::inttId)
      {
        ++inttBuckets[InttDefs::getTimeBucketId(hit.first)];
      }
    }
    for (const auto& [layer, hits] : hitsByLayer)
    {
      LayerMatch lm;
      lm.layer = layer;
      lm.nChainHits = static_cast<int>(hits.size());
      std::set<TrkrDefs::hitsetkey> sensors;
      std::map<TrkrDefs::cluskey, int> shared;
      std::set<TrkrDefs::cluskey> fromOtherStrobe;
      for (const auto& hit : hits)
      {
        sensors.insert(hit.first);
        bool otherStrobe = false;
        const auto ckey = d.officialClusterOfHit(assoc, clusters, hit.first, hit.second, otherStrobe);
        if (ckey != TrkrDefs::CLUSKEYMAX)
        {
          ++shared[ckey];
          if (otherStrobe)
          {
            fromOtherStrobe.insert(ckey);
          }
        }
      }
      lm.nSensors = static_cast<int>(sensors.size());
      lm.nClustersTouched = static_cast<int>(shared.size());
      for (const auto& [ckey, n] : shared)
      {
        if (n > lm.nShared)
        {
          lm.nShared = n;
          lm.key = ckey;
        }
      }
      if (lm.key != TrkrDefs::CLUSKEYMAX && clusters && geometry)
      {
        lm.cluster = clusters->findCluster(lm.key);
        if (lm.cluster)
        {
          const auto range = assoc->getHits(lm.key);
          lm.nClusterHits = static_cast<int>(std::distance(range.first, range.second));
          const Acts::Vector3 g = geometry->getGlobalPosition(lm.key, lm.cluster);
          const Acts::Vector3 o = toOwn(g, layer);
          lm.cx = o.x();
          lm.cy = o.y();
          lm.cz = o.z();
          lm.hasCluster = true;
          lm.otherStrobe = fromOtherStrobe.contains(lm.key);
          oc.clusters.insert(lm.key);
        }
      }
      oc.layers[layer] = lm;
    }
    // parameters, official frame
    if (traj->get_fit_status() == 1 || traj->get_fit_status() == 4)
    {
      const double x = traj->get_pca_x() + m_siBeam[0], y = traj->get_pca_y() + m_siBeam[1], z = traj->get_z0() - m_zShift;
      oc.par.pt = traj->get_pt();
      oc.par.phi = traj->get_phi();
      oc.par.eta = std::asinh(traj->get_tanl());
      oc.par.d0 = d0Official(x, y, z, traj->get_phi());
      oc.par.z0 = z;
      oc.par.charge = traj->get_charge();
    }
    oc.par.nmvtx = static_cast<int>(traj->get_n_mvtx());
    oc.par.nintt = static_cast<int>(traj->get_n_intt());
    if (!inttBuckets.empty())
    {
      oc.par.crossing = std::max_element(inttBuckets.begin(), inttBuckets.end(), [](const auto& a, const auto& b)
                                         { return a.second < b.second; })
                            ->first;
    }
    if (crossings && oc.siTpcIndex >= 0)
    {
      const SiTpc_Track* t = siTpc->get_track(oc.siTpcIndex);
      if (const TpcCrossingDecision* dec = crossings->get_decision(t->get_tpc_assembled_track_id()))
      {
        oc.resolverCrossing = dec->get_selected_crossing();
      }
    }
    own.push_back(std::move(oc));
  }

  // ---- 1) cluster vs chain point and 2) seed bias
  for (const auto& oc : own)
  {
    const Si_Trajectory* traj = oc.traj;
    for (unsigned int k = 0; k < traj->size_points(); ++k)
    {
      const int layer = traj->get_point_layer(k);
      const auto lmIt = oc.layers.find(layer);
      auto& p = d.p;
      p = Impl::PointRow();
      p.event = d.event;
      p.traj = static_cast<int>(oc.index);
      p.chain = traj->get_chain_id();
      p.layer = layer;
      p.in_sitpc = oc.siTpcIndex >= 0;
      p.px = traj->get_point_x(k);
      p.py = traj->get_point_y(k);
      p.pz = traj->get_point_z(k);
      p.pt = traj->get_pt();
      auto& cnt = d.layerCounts[layer];
      ++cnt[0];
      if (lmIt != oc.layers.end())
      {
        const auto& lm = lmIt->second;
        p.nchain = lm.nChainHits;
        p.nsensors = lm.nSensors;
        p.nclusters = lm.nClustersTouched;
        p.nshared = lm.nShared;
        p.ncluster_hits = lm.nClusterHits;
        p.has_cluster = lm.hasCluster;
        p.other_strobe = lm.otherStrobe;
        if (lm.hasCluster)
        {
          p.cx = lm.cx;
          p.cy = lm.cy;
          p.cz = lm.cz;
          const double phi = std::atan2(p.py, p.px);
          const double dx = lm.cx - p.px, dy = lm.cy - p.py;
          p.d_rphi = -dx * std::sin(phi) + dy * std::cos(phi);
          p.d_r = dx * std::cos(phi) + dy * std::sin(phi);
          p.d_z = lm.cz - p.pz;
          p.phisize = lm.cluster->getPhiSize();
          p.zsize = lm.cluster->getZSize();
          const auto err = ClusterErrorPara::get_clusterv5_modified_error(lm.cluster, std::hypot(lm.cx, lm.cy), lm.key);
          p.err_rphi = std::sqrt(err.first);
          p.err_z = std::sqrt(err.second);
          const bool identical = lm.nShared == lm.nChainHits && lm.nClusterHits == lm.nShared;
          cnt[2] += identical;
          cnt[3] += lm.nClustersTouched > 1;
          cnt[4] += lm.nClusterHits > lm.nShared;
          cnt[5] += lm.otherStrobe;
          auto& sums = d.layerSums[layer];
          sums[0] += 1;
          sums[1] += p.d_rphi;
          sums[2] += p.d_rphi * p.d_rphi;
          sums[3] += p.d_z;
          sums[4] += p.d_z * p.d_z;
          sums[5] += identical;
        }
        else
        {
          ++cnt[1];
        }
      }
      else
      {
        ++cnt[1];
      }
      d.points->Fill();
    }

    // ---- seed bias: same fit with the own points and with the cluster positions
    auto& b = d.b;
    b = Impl::BiasRow();
    b.event = d.event;
    b.traj = static_cast<int>(oc.index);
    b.in_sitpc = oc.siTpcIndex >= 0;
    b.nsi = static_cast<int>(traj->size_points());
    b.pt = traj->get_pt();
    b.eta = std::asinh(traj->get_tanl());
    std::vector<SiTpcHelixFit::Point> ownPts, cluPts;
    for (unsigned int k = 0; k < traj->size_points(); ++k)
    {
      const int layer = traj->get_point_layer(k);
      const double wz = layer <= 2 ? m_siZw[0] : m_siZw[1];
      SiTpcHelixFit::Point q{traj->get_point_x(k), traj->get_point_y(k), traj->get_point_z(k), 1.0, wz, 0};
      ownPts.push_back(q);
      const auto lmIt = oc.layers.find(layer);
      if (lmIt != oc.layers.end() && lmIt->second.hasCluster)
      {
        q.x = lmIt->second.cx;
        q.y = lmIt->second.cy;
        q.z = lmIt->second.cz;
        ++b.nreplaced;
      }
      cluPts.push_back(q);
    }
    SiTpcHelixFit::Config siCfg;
    siCfg.bz = m_bz;
    siCfg.minLeverArm = m_siLever;
    const auto ro = SiTpcHelixFit::fit(ownPts, siCfg);
    const auto rc = SiTpcHelixFit::fit(cluPts, siCfg);
    if (ro.fitted() && rc.fitted())
    {
      b.si_ok = 1;
      b.si_d0 = rc.dca - ro.dca;
      b.si_phi = wrap(rc.phi - ro.phi);
      b.si_z0 = rc.z0 - ro.z0;
      b.si_tanl = rc.tanl - ro.tanl;
      b.si_relpt = ro.status == SiTpcHelixFit::Ok && rc.status == SiTpcHelixFit::Ok ? rc.pt / ro.pt - 1 : fnan;
      if (traj->get_fit_status() == 1 || traj->get_fit_status() == 4)
      {
        b.si_ref_d0 = ro.dca - traj->get_dca();
        b.si_ref_phi = wrap(ro.phi - traj->get_phi());
        b.si_ref_z0 = ro.z0 - traj->get_z0();
      }
    }
    if (oc.siTpcIndex >= 0)
    {
      const SiTpc_Track* t = siTpc->get_track(oc.siTpcIndex);
      std::vector<SiTpcHelixFit::Point> o2, c2;
      for (unsigned int k = 0; k < t->size_points(); ++k)
      {
        const int layer = t->get_point_layer(k);
        const bool isSi = t->get_point_type(k) == SiTpc_Track::SiPoint;
        const int det = isSi ? detOfLayer(layer) : 2;
        SiTpcHelixFit::Point q{t->get_point_x(k), t->get_point_y(k), t->get_point_z(k), m_stWxy[det], m_stWz[det], isSi ? 0 : 1,
                               (m_hybrid && isSi) ? m_hybridW[det] : 0.0};
        o2.push_back(q);
        const auto lmIt = isSi ? oc.layers.find(layer) : oc.layers.end();
        if (lmIt != oc.layers.end() && lmIt->second.hasCluster)
        {
          q.x = lmIt->second.cx;
          q.y = lmIt->second.cy;
          q.z = lmIt->second.cz;
        }
        c2.push_back(q);
      }
      SiTpcHelixFit::Config stCfg;
      stCfg.bz = m_bz;
      stCfg.groupZOffset = m_tpcZOffset;
      stCfg.hybridMinPoints = m_hybrid ? m_hybridMin : 0;
      const auto so = SiTpcHelixFit::fit(o2, stCfg);
      const auto sc = SiTpcHelixFit::fit(c2, stCfg);
      if (so.fitted() && sc.fitted())
      {
        b.st_ok = 1;
        b.st_d0 = sc.dca - so.dca;
        b.st_phi = wrap(sc.phi - so.phi);
        b.st_z0 = sc.z0 - so.z0;
        b.st_tanl = sc.tanl - so.tanl;
        b.st_relpt = so.status == SiTpcHelixFit::Ok && sc.status == SiTpcHelixFit::Ok ? sc.pt / so.pt - 1 : fnan;
        b.st_ref_d0 = so.dca - t->get_dca();
        b.st_ref_phi = wrap(so.phi - t->get_phi());
        b.st_ref_z0 = so.z0 - t->get_z0();
      }
    }
    d.bias->Fill();
  }

  // ---- 3) standard silicon seeds
  // the standard seeding ran only if its converted track map exists too (SiTpcActsFitter creates an empty
  // SiliconTrackSeedContainer placeholder when the standard seeding is off).  Only seeds with a converted track
  // are standard seeds: TrackSeedTrackMapConverter converts every seed, and SiTpcActsSeeder appends its own
  // seeds (for PHActsTrkFitter) to the same container after the conversion.
  auto* stdSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_stdSeedNode);
  if (stdSeeds && findNode::getClass<SvtxTrackMap>(topNode, m_stdTrackNode))
  {
    std::map<const TrackSeed*, const SvtxTrack*> trackOfSeed;
    if (auto* tmap = findNode::getClass<SvtxTrackMap>(topNode, m_stdTrackNode))
    {
      for (const auto& [key, trk] : *tmap)
      {
        if (trk && trk->get_silicon_seed())
        {
          trackOfSeed[trk->get_silicon_seed()] = trk;
        }
      }
    }
    std::vector<StdSeed> stdList;
    std::map<TrkrDefs::cluskey, std::vector<int>> seedsOfCluster;
    for (std::size_t i = 0; i < stdSeeds->size(); ++i)
    {
      const TrackSeed* seed = stdSeeds->get(i);
      if (!seed || !trackOfSeed.contains(seed))
      {
        continue;
      }
      StdSeed ss;
      ss.seed = seed;
      for (auto it = seed->begin_cluster_keys(); it != seed->end_cluster_keys(); ++it)
      {
        const auto id = TrkrDefs::getTrkrId(*it);
        if (id == TrkrDefs::mvtxId || id == TrkrDefs::inttId)
        {
          ss.clusters.insert(*it);
          ss.par.nmvtx += id == TrkrDefs::mvtxId;
          ss.par.nintt += id == TrkrDefs::inttId;
          seedsOfCluster[*it].push_back(static_cast<int>(stdList.size()));
        }
      }
      ss.par.crossing = seed->get_crossing();
      const auto t = trackOfSeed.find(seed);
      if (t != trackOfSeed.end())
      {
        const SvtxTrack* trk = t->second;
        const double phi = std::atan2(trk->get_py(), trk->get_px());
        ss.par.pt = std::hypot(trk->get_px(), trk->get_py());
        ss.par.phi = phi;
        ss.par.eta = trk->get_eta();
        ss.par.d0 = d0Official(trk->get_x(), trk->get_y(), trk->get_z(), phi);
        ss.par.z0 = trk->get_z();
        ss.par.charge = trk->get_charge();
      }
      stdList.push_back(ss);
    }

    // own chain -> standard seed
    std::vector<int> bestOwnOfStd(stdList.size(), -1), bestOwnShared(stdList.size(), 0);
    for (std::size_t io = 0; io < own.size(); ++io)
    {
      const auto& oc = own[io];
      std::map<int, int> sharedWith;
      for (const auto key : oc.clusters)
      {
        const auto f = seedsOfCluster.find(key);
        if (f != seedsOfCluster.end())
        {
          for (const int is : f->second)
          {
            ++sharedWith[is];
          }
        }
      }
      int best = -1, nbest = 0;
      for (const auto& [is, n] : sharedWith)
      {
        if (n > nbest || (n == nbest && best >= 0 && std::abs(wrap(stdList[is].par.phi - oc.par.phi)) < std::abs(wrap(stdList[best].par.phi - oc.par.phi))))
        {
          nbest = n;
          best = is;
        }
        if (n > bestOwnShared[is])
        {
          bestOwnShared[is] = n;
          bestOwnOfStd[is] = static_cast<int>(io);
        }
      }
      auto& c = d.c;
      c = Impl::ChainRow();
      c.event = d.event;
      c.traj = static_cast<int>(oc.index);
      c.in_sitpc = oc.siTpcIndex >= 0;
      c.nown = static_cast<int>(oc.clusters.size());
      c.nshared = nbest;
      c.std_found = best >= 0;
      c.own_crossing = oc.par.crossing;
      c.resolver_crossing = oc.resolverCrossing;
      c.own_nmvtx = oc.par.nmvtx;
      c.own_nintt = oc.par.nintt;
      c.own_pt = oc.par.pt;
      c.own_phi = oc.par.phi;
      c.own_eta = oc.par.eta;
      c.own_d0 = oc.par.d0;
      c.own_z0 = oc.par.z0;
      c.std_crossing = SHRT_MAX;
      if (best >= 0)
      {
        const auto& ss = stdList[best];
        c.nstd = static_cast<int>(ss.clusters.size());
        c.std_crossing = ss.par.crossing;
        c.std_nmvtx = ss.par.nmvtx;
        c.std_nintt = ss.par.nintt;
        c.std_pt = ss.par.pt;
        c.std_phi = ss.par.phi;
        c.std_eta = ss.par.eta;
        c.std_d0 = ss.par.d0;
        c.std_z0 = ss.par.z0;
      }
      d.ownChains->Fill();

      const bool all = c.nown > 0 && nbest == c.nown && c.nstd == c.nown;
      const bool none = nbest == 0;
      ++d.nChains;
      d.nChainsAll += all;
      d.nChainsNone += none;
      d.nChainsPartial += !all && !none;
      if (c.in_sitpc)
      {
        ++d.nChainsSiTpc;
        d.nChainsSiTpcAll += all;
        d.nChainsSiTpcNone += none;
        d.nChainsSiTpcPartial += !all && !none;
      }
    }

    // standard seed -> own chain, closest own Si+TPC track
    for (std::size_t is = 0; is < stdList.size(); ++is)
    {
      const auto& ss = stdList[is];
      auto& s = d.s;
      s = Impl::SeedRow();
      s.event = d.event;
      s.seed = static_cast<int>(is);
      s.nstd = static_cast<int>(ss.clusters.size());
      s.nshared = bestOwnShared[is];
      s.own_found = bestOwnOfStd[is] >= 0 && bestOwnShared[is] >= 2;
      if (bestOwnOfStd[is] >= 0)
      {
        s.nown = static_cast<int>(own[bestOwnOfStd[is]].clusters.size());
        s.own_in_sitpc = own[bestOwnOfStd[is]].siTpcIndex >= 0;
      }
      s.crossing = ss.par.crossing;
      s.nmvtx = ss.par.nmvtx;
      s.nintt = ss.par.nintt;
      s.pt = ss.par.pt;
      s.phi = ss.par.phi;
      s.eta = ss.par.eta;
      s.d0 = ss.par.d0;
      s.z0 = ss.par.z0;
      if (siTpc && std::isfinite(ss.par.phi))
      {
        double bestScore = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < siTpc->size(); ++i)
        {
          const SiTpc_Track* t = siTpc->get_track(i);
          if (!t || t->get_fit_status() != 1)
          {
            continue;
          }
          const double dphi = wrap(t->get_phi() - ss.par.phi), deta = std::asinh(t->get_tanl()) - ss.par.eta;
          const double score = dphi * dphi + deta * deta;
          if (score < bestScore)
          {
            bestScore = score;
            s.near_sitpc_dphi = static_cast<float>(dphi);
            s.near_sitpc_deta = static_cast<float>(deta);
            s.near_sitpc_pt = t->get_pt();
          }
        }
      }
      d.stdSeeds->Fill();
      ++d.nStd;
      d.nStdMatched += s.own_found;
      d.nStdNearSiTpc += !s.own_found && std::abs(s.near_sitpc_dphi) < 0.05 && std::abs(s.near_sitpc_deta) < 0.05;
    }
  }

  // ---- 4) Si-only vertices
  auto* ownV = findNode::getClass<SvtxVertexMap>(topNode, m_ownVertexNode);
  auto* stdV = findNode::getClass<SvtxVertexMap>(topNode, m_stdVertexNode);
  if (ownV && stdV)
  {
    std::set<unsigned int> usedStd;
    for (const auto& [ko, vo] : *ownV)
    {
      auto& v = d.v;
      v = Impl::VertexRow();
      v.event = d.event;
      v.type = 0;
      v.crossing = vo->get_beam_crossing();
      v.own_ntrk = static_cast<int>(vo->size_tracks());
      v.own_x = vo->get_x() + m_siBeam[0];
      v.own_y = vo->get_y() + m_siBeam[1];
      v.own_z = vo->get_z() - m_zShift;
      double bestDz = 1.0;
      unsigned int best = UINT_MAX;
      for (const auto& [ks, vs] : *stdV)
      {
        if (usedStd.contains(ks) || vs->get_beam_crossing() != vo->get_beam_crossing())
        {
          continue;
        }
        const double dz = std::abs(vs->get_z() - v.own_z);
        if (dz < bestDz)
        {
          bestDz = dz;
          best = ks;
        }
      }
      ++d.nOwnVtx;
      if (best != UINT_MAX)
      {
        usedStd.insert(best);
        const SvtxVertex* vs = stdV->get(best);
        v.matched = 1;
        v.std_ntrk = static_cast<int>(vs->size_tracks());
        v.std_x = vs->get_x();
        v.std_y = vs->get_y();
        v.std_z = vs->get_z();
        ++d.nOwnVtxMatched;
      }
      d.vertices->Fill();
    }
    for (const auto& [ks, vs] : *stdV)
    {
      ++d.nStdVtx;
      if (usedStd.contains(ks))
      {
        continue;
      }
      auto& v = d.v;
      v = Impl::VertexRow();
      v.event = d.event;
      v.type = 1;
      v.crossing = vs->get_beam_crossing();
      v.std_ntrk = static_cast<int>(vs->size_tracks());
      v.std_x = vs->get_x();
      v.std_y = vs->get_y();
      v.std_z = vs->get_z();
      d.vertices->Fill();
    }
  }

  ++d.event;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiSeedComparisonQA::End(PHCompositeNode* /*topNode*/)
{
  auto& d = *m_impl;
  std::cout << Name() << "::End - " << d.event << " events" << std::endl;
  std::cout << "  own chain point vs official cluster (own frame) [um]:" << std::endl
            << "  layer   points  no-cluster  other-strobe  identical  split  merged     <d_rphi>  rms(d_rphi)     <d_z>  rms(d_z)" << std::endl;
  for (const auto& [layer, cnt] : d.layerCounts)
  {
    const auto& s = d.layerSums[layer];
    const double n = std::max(1., s[0]);
    const double mr = s[1] / n, mz = s[3] / n;
    std::cout << std::fixed << std::setprecision(1) << "  " << std::setw(5) << layer << std::setw(9) << cnt[0] << std::setw(12) << cnt[1] << std::setw(14) << cnt[5]
              << std::setw(11) << cnt[2] << std::setw(7) << cnt[3] << std::setw(8) << cnt[4] << std::setw(13) << 1e4 * mr << std::setw(13)
              << 1e4 * std::sqrt(std::max(0., s[2] / n - mr * mr)) << std::setw(10) << 1e4 * mz << std::setw(10)
              << 1e4 * std::sqrt(std::max(0., s[4] / n - mz * mz)) << std::endl;
  }
  std::cout << std::defaultfloat;
  if (d.nChains > 0 && d.nStd > 0)
  {
    std::cout << "  own chains " << d.nChains << ": identical standard seed " << d.nChainsAll << ", partial " << d.nChainsPartial
              << ", none " << d.nChainsNone << std::endl
              << "    of them in SITPC_TRACKS " << d.nChainsSiTpc << ": identical " << d.nChainsSiTpcAll << ", partial "
              << d.nChainsSiTpcPartial << ", none " << d.nChainsSiTpcNone << std::endl
              << "  standard seeds " << d.nStd << ": with an own chain (>= 2 shared clusters) " << d.nStdMatched
              << ", without one but an own Si+TPC track within 0.05 in phi and eta " << d.nStdNearSiTpc << std::endl;
  }
  if (d.nOwnVtx > 0 || d.nStdVtx > 0)
  {
    std::cout << "  Si vertices: own " << d.nOwnVtx << " (matched " << d.nOwnVtxMatched << "), standard " << d.nStdVtx << std::endl;
  }
  if (d.file)
  {
    d.file->cd();
    d.points->Write();
    d.bias->Write();
    d.ownChains->Write();
    d.stdSeeds->Write();
    d.vertices->Write();
    d.file->Close();
    delete d.file;
    d.file = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
