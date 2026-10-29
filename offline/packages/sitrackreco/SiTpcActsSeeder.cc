#include "SiTpcActsSeeder.h"

#include "SiTpc_Track.h"
#include "SiTpc_TrackContainer.h"
#include "Si_Trajectory.h"
#include "Si_TrajectoryContainer.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>
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

#include <TFile.h>
#include <TTree.h>

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

  // circle + straight line in (s, z) in the official frame, from the own Si+TPC fit
  struct Helix
  {
    double cx = 0, cy = 0, R = 0;
    double omega = 1;  // +1: counter-clockwise seen from +z
    double x0 = 0, y0 = 0, z0 = 0, tanl = 0;  // start point (closest to the beam axis) and z there

    double angleOf(double x, double y) const { return std::atan2(y - cy, x - cx); }

    // signed path length in xy from the start point along the direction of motion [cm]
    double arcLength(double x, double y) const
    {
      return R * std::remainder(omega * (angleOf(x, y) - angleOf(x0, y0)), twopi);
    }

    // unit direction of motion at a point on (or near) the circle
    Acts::Vector3 direction(double x, double y) const
    {
      const double rx = (x - cx) / R, ry = (y - cy) / R;
      return Acts::Vector3(-omega * ry, omega * rx, tanl).normalized();
    }
  };
}  // namespace

//____________________________________________________________________________________________________
struct SiTpcActsSeeder::Impl
{
  // nodes
  ActsGeometry* geometry = nullptr;
  SiTpc_TrackContainer* siTpcTracks = nullptr;
  Si_TrajectoryContainer* siTrajectories = nullptr;
  Tpc_PolyClusterContainer* tpcClusters = nullptr;
  TpcCrossingDecisionContainer* crossings = nullptr;
  TrkrClusterContainer* siClusters = nullptr;
  TrkrClusterHitAssoc* hitAssoc = nullptr;
  TrkrClusterContainer* outClusters = nullptr;
  TrackSeedContainer* siSeeds = nullptr;
  TrackSeedContainer* tpcSeeds = nullptr;
  TrackSeedContainer* svtxSeeds = nullptr;
  bool ownSiSeeds = false, ownTpcSeeds = false;
  double bzSign = 1;
  bool distortionWarned = false;

  // per event
  std::unordered_map<TrkrDefs::hitsetkey, std::unordered_map<TrkrDefs::hitkey, TrkrDefs::cluskey>> hitToCluster;
  std::map<TrkrDefs::hitsetkey, uint32_t> nextIndex;  // free cluster index per hitset for new keys
  std::map<TrkrDefs::cluskey, Acts::Vector3> written;   // TPC key -> intended position (duplicate check)

  // counters
  unsigned long nEvents = 0, nInput = 0, nSeeds = 0;
  unsigned long nSiOfficialAssoc = 0, nSiOfficialPosition = 0, nSiProjected = 0, nSiMissing = 0, nSiOtherStrobe = 0;
  unsigned long nTpc = 0, nTpcAngleCut = 0, nTpcMoveCut = 0, nTpcNoSurface = 0, nTpcDuplicateKey = 0;
  double sumMove = 0, maxMove = 0, sumClosure = 0, maxClosure = 0;
  unsigned long nClosure = 0;
  std::map<std::string, unsigned long> rejected;

  // QA
  TFile* qaFile = nullptr;
  TTree* clusterTree = nullptr;
  TTree* seedTree = nullptr;
  struct ClusterRow
  {
    int event = 0, seed = 0, layer = 0, side = 0, crossing = 0;
    float x = NAN, y = NAN, z = NAN;     // intended (own) position, official frame [cm]
    float sx = NAN, sy = NAN, sz = NAN;  // position on the official surface [cm]
    float gx = NAN, gy = NAN, gz = NAN;  // official ActsGeometry global position of the written cluster [cm]
    float move = NAN, dr = NAN, angle = NAN, closure = NAN, time = NAN;
  } clusterRow;
  struct SeedRow
  {
    int event = 0, seed = 0, crossing = 0, charge = 0, nmvtx = 0, nintt = 0, ntpc = 0, nsiProjected = 0;
    float pt = NAN, phi = NAN, tanl = NAN, z0 = NAN, X0 = NAN, Y0 = NAN, R = NAN, tpcZOffset = NAN;
  } seedRow;

  TrkrDefs::cluskey clusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk);
  TrkrDefs::cluskey officialClusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, bool& otherStrobe);
  uint32_t freeIndex(TrkrDefs::hitsetkey hsk);
};

//____________________________________________________________________________________________________
TrkrDefs::cluskey SiTpcActsSeeder::Impl::clusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk)
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
TrkrDefs::cluskey SiTpcActsSeeder::Impl::officialClusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, bool& otherStrobe)
{
  // MvtxClusterPruner removes the copy of a cluster duplicated in the next strobe from TRKR_CLUSTER without
  // updating TRKR_CLUSTERHITASSOC: then the cluster with the same pixel in the neighbouring strobe is used
  otherStrobe = false;
  const auto ckey = clusterOfHit(hsk, hk);
  if (ckey != TrkrDefs::CLUSKEYMAX && siClusters->findCluster(ckey))
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
    if (other != TrkrDefs::CLUSKEYMAX && siClusters->findCluster(other))
    {
      otherStrobe = true;
      return other;
    }
  }
  return TrkrDefs::CLUSKEYMAX;
}

//____________________________________________________________________________________________________
uint32_t SiTpcActsSeeder::Impl::freeIndex(TrkrDefs::hitsetkey hsk)
{
  // TrkrClusterContainerv4 stores the clusters of a hitset in a vector indexed by the cluster index:
  // new keys continue after the largest index used in the input containers
  auto it = nextIndex.find(hsk);
  if (it == nextIndex.end())
  {
    uint32_t next = 0;
    if (TrkrDefs::getTrkrId(hsk) != TrkrDefs::tpcId)
    {
      const auto range = siClusters->getClusters(hsk);
      for (auto c = range.first; c != range.second; ++c)
      {
        next = std::max(next, TrkrDefs::getClusIndex(c->first) + 1);
      }
    }
    it = nextIndex.emplace(hsk, next).first;
  }
  return it->second++;
}

//____________________________________________________________________________________________________
SiTpcActsSeeder::SiTpcActsSeeder(const std::string& name)
  : SubsysReco(name)
  , m_impl(std::make_unique<Impl>())
{
}

SiTpcActsSeeder::~SiTpcActsSeeder() = default;

//____________________________________________________________________________________________________
int SiTpcActsSeeder::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  auto& d = *m_impl;

  // charge of the seed: from the sense of rotation and the field (positive charge, B_z > 0: clockwise).
  // With B_z > 0 this is the convention of TrackSeedHelper::circleFitByTaubin.
  {
    const auto& field = *d.geometry->geometry().magField;
    auto cache = field.makeCache(d.geometry->geometry().magFieldContext);
    const auto b = field.getField(Acts::Vector3(0, 0, 0), cache);
    d.bzSign = (b.ok() && b.value().z() < 0) ? -1. : 1.;
    std::cout << Name() << "::InitRun - B_z at the centre " << (b.ok() ? b.value().z() / Acts::UnitConstants::T : NAN) << " T" << std::endl;
  }

  checkDistortionContainers(topNode);

  if (!m_qaFileName.empty())
  {
    d.qaFile = new TFile(m_qaFileName.c_str(), "RECREATE");
    d.clusterTree = new TTree("tpc_clusters", "SiTpcActsSeeder TPC clusters: own position, on the official surface, closure");
    auto& c = d.clusterRow;
    d.clusterTree->Branch("event", &c.event);
    d.clusterTree->Branch("seed", &c.seed);
    d.clusterTree->Branch("layer", &c.layer);
    d.clusterTree->Branch("side", &c.side);
    d.clusterTree->Branch("crossing", &c.crossing);
    d.clusterTree->Branch("x", &c.x);
    d.clusterTree->Branch("y", &c.y);
    d.clusterTree->Branch("z", &c.z);
    d.clusterTree->Branch("sx", &c.sx);
    d.clusterTree->Branch("sy", &c.sy);
    d.clusterTree->Branch("sz", &c.sz);
    d.clusterTree->Branch("gx", &c.gx);
    d.clusterTree->Branch("gy", &c.gy);
    d.clusterTree->Branch("gz", &c.gz);
    d.clusterTree->Branch("move", &c.move);
    d.clusterTree->Branch("dr", &c.dr);
    d.clusterTree->Branch("angle", &c.angle);
    d.clusterTree->Branch("closure", &c.closure);
    d.clusterTree->Branch("time", &c.time);

    d.seedTree = new TTree("seeds", "SiTpcActsSeeder seeds (official frame)");
    auto& s = d.seedRow;
    d.seedTree->Branch("event", &s.event);
    d.seedTree->Branch("seed", &s.seed);
    d.seedTree->Branch("crossing", &s.crossing);
    d.seedTree->Branch("charge", &s.charge);
    d.seedTree->Branch("nmvtx", &s.nmvtx);
    d.seedTree->Branch("nintt", &s.nintt);
    d.seedTree->Branch("ntpc", &s.ntpc);
    d.seedTree->Branch("nsi_projected", &s.nsiProjected);
    d.seedTree->Branch("pt", &s.pt);
    d.seedTree->Branch("phi", &s.phi);
    d.seedTree->Branch("tanl", &s.tanl);
    d.seedTree->Branch("z0", &s.z0);
    d.seedTree->Branch("X0", &s.X0);
    d.seedTree->Branch("Y0", &s.Y0);
    d.seedTree->Branch("R", &s.R);
    d.seedTree->Branch("tpc_z_offset", &s.tpcZOffset);
  }

  std::cout << Name() << "::InitRun - configuration" << std::endl
            << "  Si beam spot (" << m_siBeam[0] << ", " << m_siBeam[1] << ") cm, TPC beam line (" << m_tpcLine[0] << ", "
            << m_tpcLine[1] << ") cm, Si acts z shift " << m_siActsZShift << " cm" << std::endl
            << "  TPC extra shift (" << m_tpcExtraShift[0] << ", " << m_tpcExtraShift[1] << ", " << m_tpcExtraShift[2] << ") cm"
            << ", per-track TPC z offset " << (m_applyTrackTpcZOffset ? "applied" : "not applied")
            << ", INTT offset (" << m_inttOffset[0] << ", " << m_inttOffset[1] << ") cm (in-plane)" << std::endl
            << "  TPC placement " << (m_tpcPlacement == TpcPlacement::AlongTrack ? "along the own helix" : "along the surface normal")
            << ", max crossing angle " << m_maxTpcCrossingAngle << " rad, max move " << m_maxTpcMove << " cm" << std::endl
            << "  output clusters " << m_clusterOutName << ", seeds SiliconTrackSeedContainer"
            << (d.ownSiSeeds ? "" : " (appended)") << " / TpcTrackSeedContainer" << (d.ownTpcSeeds ? "" : " (appended)")
            << " / " << m_svtxSeedName << (m_qaFileName.empty() ? "" : ", QA " + m_qaFileName) << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
void SiTpcActsSeeder::checkDistortionContainers(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
  if (d.distortionWarned)
  {
    return;
  }
  for (const std::string name : {"TpcDistortionCorrectionContainerModuleEdge", "TpcDistortionCorrectionContainerStatic",
                                 "TpcDistortionCorrectionContainerAverage", "TpcDistortionCorrectionContainerFluctuation"})
  {
    PHNodeIterator iter(topNode);
    if (iter.findFirst("PHDataNode", name) || iter.findFirst("PHIODataNode", name))
    {
      std::cout << Name() << " - WARNING " << name << " is on the node tree: PHActsTrkFitter will apply the official TPC"
                << " distortion correction ON TOP of the own (PHGarfield) positions.  Set G4TPC::ENABLE_*_CORRECTIONS = false."
                << std::endl;
      d.distortionWarned = true;
    }
  }
}

//____________________________________________________________________________________________________
int SiTpcActsSeeder::getNodes(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
  d.geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  d.siTpcTracks = findNode::getClass<SiTpc_TrackContainer>(topNode, m_siTpcTrackNodeName);
  d.siTrajectories = findNode::getClass<Si_TrajectoryContainer>(topNode, m_siTrajectoryNodeName);
  d.tpcClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  d.crossings = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingNodeName);
  d.siClusters = findNode::getClass<TrkrClusterContainer>(topNode, m_siClusterNodeName);
  d.hitAssoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

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
  require(d.tpcClusters, m_tpcClusterNodeName);
  require(d.siClusters, m_siClusterNodeName);
  if (!d.crossings)
  {
    std::cout << Name() << ": no " << m_crossingNodeName << " - crossing 0 for every track" << std::endl;
  }
  if (!d.hitAssoc)
  {
    std::cout << Name() << ": no TRKR_CLUSTERHITASSOC - Si clusters matched by position only" << std::endl;
  }
  if (m_clusterOutName == m_siClusterNodeName)
  {
    std::cout << PHWHERE << " " << Name() << ": the output cluster container must not be the input " << m_siClusterNodeName << std::endl;
    ok = false;
  }
  return ok ? Fun4AllReturnCodes::EVENT_OK : Fun4AllReturnCodes::ABORTRUN;
}

//____________________________________________________________________________________________________
int SiTpcActsSeeder::createNodes(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
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
  auto* trkrNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "TRKR"));
  if (!trkrNode)
  {
    trkrNode = new PHCompositeNode("TRKR");
    dstNode->addNode(trkrNode);
  }

  d.outClusters = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterOutName);
  if (!d.outClusters)
  {
    d.outClusters = new TrkrClusterContainerv4;
    trkrNode->addNode(new PHIODataNode<PHObject>(d.outClusters, m_clusterOutName, "PHObject"));
  }

  auto seedContainer = [&](const std::string& name, bool& own) -> TrackSeedContainer*
  {
    auto* c = findNode::getClass<TrackSeedContainer>(topNode, name);
    own = (c == nullptr);
    if (!c)
    {
      c = new TrackSeedContainer_v1;
      svtxNode->addNode(new PHIODataNode<PHObject>(c, name, "PHObject"));
    }
    return c;
  };
  d.siSeeds = seedContainer("SiliconTrackSeedContainer", d.ownSiSeeds);
  d.tpcSeeds = seedContainer("TpcTrackSeedContainer", d.ownTpcSeeds);
  bool ownSvtx = false;
  d.svtxSeeds = seedContainer(m_svtxSeedName, ownSvtx);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcActsSeeder::process_event(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
  ++d.nEvents;
  if (d.nEvents == 1)
  {
    checkDistortionContainers(topNode);
  }
  d.outClusters->Reset();
  d.svtxSeeds->Reset();
  if (d.ownSiSeeds)
  {
    d.siSeeds->Reset();
  }
  if (d.ownTpcSeeds)
  {
    d.tpcSeeds->Reset();
  }
  d.hitToCluster.clear();
  d.nextIndex.clear();
  d.written.clear();

  // TPC: the poly-cluster keys use the indices 0..n-1 per hitset; duplicates get indices after them
  for (unsigned int i = 0; i < d.tpcClusters->size(); ++i)
  {
    const Tpc_PolyCluster* c = d.tpcClusters->get_cluster(i);
    if (c && c->get_trkr_cluster_key() != TrkrDefs::CLUSKEYMAX)
    {
      const auto hsk = TrkrDefs::getHitSetKeyFromClusKey(c->get_trkr_cluster_key());
      auto& next = d.nextIndex[hsk];
      next = std::max(next, TrkrDefs::getClusIndex(c->get_trkr_cluster_key()) + 1);
    }
  }

  const auto gctx = d.geometry->geometry().getGeoContext();  // the context of ActsGeometry::getGlobalPosition
  const bool shiftIntt = m_inttOffset[0] != 0 || m_inttOffset[1] != 0;
  const Acts::Vector3 inttShiftMm(m_inttOffset[0] * cm, m_inttOffset[1] * cm, 0);

  auto reject = [&](const std::string& why)
  {
    ++d.rejected[why];
  };

  for (unsigned int itrack = 0; itrack < d.siTpcTracks->size(); ++itrack)
  {
    const SiTpc_Track* track = d.siTpcTracks->get_track(itrack);
    if (!track)
    {
      continue;
    }
    ++d.nInput;

    // ---- selection
    if (track->get_fit_status() != 1)  // SiTpcHelixFit::Ok (StraightLine has no curvature)
    {
      reject("own fit not ok");
      continue;
    }
    const double pt = track->get_pt();
    if (track->get_charge() == 0 || !std::isfinite(pt) || pt < m_minPt || !(track->get_radius() > 0))
    {
      reject("charge / pt");
      continue;
    }
    if (track->get_si_trajectory_index() >= d.siTrajectories->size())
    {
      reject("no Si trajectory");
      continue;
    }
    const Si_Trajectory* traj = d.siTrajectories->get_trajectory(track->get_si_trajectory_index());
    if (!traj)
    {
      reject("no Si trajectory");
      continue;
    }

    short crossing = 0;
    if (d.crossings)
    {
      const TpcCrossingDecision* decision = d.crossings->get_decision(track->get_tpc_assembled_track_id());
      crossing = decision ? decision->get_selected_crossing() : SHRT_MAX;
    }
    if (crossing == SHRT_MAX)
    {
      reject("no crossing");
      continue;
    }

    // ---- helix in the official frame: own frame = official frame - Si beam spot (z: + siActsZShift)
    Helix h;
    h.cx = track->get_circle_x() + m_siBeam[0];
    h.cy = track->get_circle_y() + m_siBeam[1];
    h.R = track->get_radius();
    {
      const double vx = m_siBeam[0] - h.cx, vy = m_siBeam[1] - h.cy;
      const double v = std::hypot(vx, vy);
      if (!(v > 0))
      {
        reject("degenerate circle");
        continue;
      }
      h.x0 = h.cx + h.R * vx / v;
      h.y0 = h.cy + h.R * vy / v;
      const double tx = std::cos(track->get_phi()), ty = std::sin(track->get_phi());
      h.omega = ((h.x0 - h.cx) * ty - (h.y0 - h.cy) * tx) >= 0 ? 1. : -1.;
      h.tanl = track->get_tanl();
      h.z0 = track->get_z0() - m_siActsZShift;
    }
    if (!std::isfinite(h.cx) || !std::isfinite(h.cy) || !std::isfinite(h.tanl) || !std::isfinite(h.z0) || !(std::hypot(h.cx, h.cy) > 0))
    {
      reject("invalid helix");
      continue;
    }

    // ---- Si: chain points and hits per layer
    std::map<int, Acts::Vector3> chainPoint;  // own frame [cm]
    for (unsigned int k = 0; k < traj->size_points(); ++k)
    {
      const int layer = traj->get_point_layer(k);
      if (layer >= 0 && layer < 7 && !chainPoint.contains(layer))
      {
        chainPoint[layer] = Acts::Vector3(traj->get_point_x(k), traj->get_point_y(k), traj->get_point_z(k));
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

    // clusters are written to the output container only once the track is accepted
    struct NewCluster
    {
      TrkrDefs::cluskey key = TrkrDefs::CLUSKEYMAX;
      TrkrClusterv5 cluster;
    };
    std::vector<TrkrDefs::cluskey> siKeys;
    std::vector<NewCluster> siNew;
    unsigned int nmvtx = 0, nintt = 0, nProjected = 0;
    for (const auto& [layer, hits] : hitsByLayer)
    {
      const bool isIntt = layer >= 3;
      std::optional<Acts::Vector3> target;  // chain point, official frame [cm]
      if (const auto cp = chainPoint.find(layer); cp != chainPoint.end())
      {
        target = Acts::Vector3(cp->second.x() + m_siBeam[0], cp->second.y() + m_siBeam[1], cp->second.z() - m_siActsZShift);
      }
      std::map<TrkrDefs::hitsetkey, int> nPerSensor;
      for (const auto& hit : hits)
      {
        ++nPerSensor[hit.first];
      }
      const auto sensor = std::max_element(nPerSensor.begin(), nPerSensor.end(), [](const auto& a, const auto& b)
                                           { return a.second < b.second; })
                              ->first;

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
      if (best == TrkrDefs::CLUSKEYMAX && target)
      {
        if (const auto surface = d.geometry->maps().getSiliconSurface(sensor))
        {
          const Acts::Vector3 lt = surface->localToGlobalTransform(gctx).inverse() * (*target * cm);
          const auto& win = isIntt ? kInttMatchWindow : kMvtxMatchWindow;
          double bestDist = std::numeric_limits<double>::max();
          const auto range = d.siClusters->getClusters(sensor);
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
        ++(byPosition ? d.nSiOfficialPosition : d.nSiOfficialAssoc);
        d.nSiOtherStrobe += bestOtherStrobe;
        ++(isIntt ? nintt : nmvtx);
        continue;
      }
      // 3) chain point projected on the official surface of its sensor (new cluster)
      if (m_useProjected && target)
      {
        if (const auto surface = d.geometry->maps().getSiliconSurface(sensor))
        {
          const Acts::Vector3 lt = surface->localToGlobalTransform(gctx).inverse() * (*target * cm);
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
          siNew.push_back(std::move(nc));
          ++nProjected;
          ++(isIntt ? nintt : nmvtx);
          continue;
        }
      }
      ++d.nSiMissing;
    }
    if (nmvtx < m_minMvtx || nmvtx + nintt < m_minSi)
    {
      reject("too few Si clusters");
      continue;
    }

    // ---- TPC: own cluster positions in the official frame, moved to the official surface
    const double zOffset = m_applyTrackTpcZOffset ? track->get_tpc_z_offset() : 0.;
    std::vector<NewCluster> tpcNew;
    std::vector<Impl::ClusterRow> tpcRows;
    for (unsigned int ip = 0; ip < track->size_points(); ++ip)
    {
      if (track->get_point_type(ip) != SiTpc_Track::TpcPoint)
      {
        continue;
      }
      const unsigned int icluster = track->get_point_source(ip);
      const Tpc_PolyCluster* poly = icluster < d.tpcClusters->size() ? d.tpcClusters->get_cluster(icluster) : nullptr;
      if (!poly || !poly->isValid() || poly->size_hits() == 0)
      {
        continue;
      }
      const TrkrDefs::hitsetkey hitsetkey = poly->get_hit_index(0).first;
      const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
      const unsigned int side = TpcDefs::getSide(hitsetkey);
      const Acts::Vector3 pos(poly->get_centroid_x() - m_tpcLine[0] + m_siBeam[0] + m_tpcExtraShift[0],
                              poly->get_centroid_y() - m_tpcLine[1] + m_siBeam[1] + m_tpcExtraShift[1],
                              poly->get_centroid_z() - m_siActsZShift - zOffset + m_tpcExtraShift[2]);
      if (!pos.allFinite())
      {
        continue;
      }
      ++d.nTpc;

      // official surface (phi segment) and the position on it; the segment is searched again after the move
      const Acts::Vector3 dir = h.direction(pos.x(), pos.y());
      Acts::Vector3 onSurface = pos;
      TrkrDefs::subsurfkey subsurfkey = 0;
      std::shared_ptr<const Acts::Surface> surface;
      double angle = 0;
      bool ok = true;
      for (int iter = 0; iter < 2 && ok; ++iter)
      {
        surface = d.geometry->get_tpc_surface_from_coords(hitsetkey, onSurface, subsurfkey);
        if (!surface)
        {
          ok = false;
          break;
        }
        const Acts::Transform3 transform = surface->localToGlobalTransform(gctx);
        const Acts::Vector3 normal = transform.linear().col(2);
        const Acts::Vector3 center = transform.translation() / cm;
        const double nd = normal.dot(dir);
        angle = std::acos(std::min(1., std::abs(nd)));
        if (m_tpcPlacement == TpcPlacement::AlongTrack)
        {
          if (std::abs(nd) < 1e-3)
          {
            ok = false;
            break;
          }
          onSurface = pos + dir * (normal.dot(center - pos) / nd);
        }
        else
        {
          onSurface = pos + normal * normal.dot(center - pos);
        }
      }
      if (!ok || !surface)
      {
        ++d.nTpcNoSurface;
        continue;
      }
      if (angle > m_maxTpcCrossingAngle)
      {
        ++d.nTpcAngleCut;
        continue;
      }
      const double move = (onSurface - pos).norm();
      if (move > m_maxTpcMove)
      {
        ++d.nTpcMoveCut;
        continue;
      }

      // key: the poly cluster key; a new one if the cluster was already written at another position
      TrkrDefs::cluskey key = poly->get_trkr_cluster_key();
      if (key == TrkrDefs::CLUSKEYMAX || TrkrDefs::getTrkrId(key) != TrkrDefs::tpcId)
      {
        key = TrkrDefs::genClusKey(hitsetkey, d.freeIndex(hitsetkey));
      }
      if (const auto w = d.written.find(key); w != d.written.end() && (w->second - onSurface).norm() > 1e-4)
      {
        key = TrkrDefs::genClusKey(TrkrDefs::getHitSetKeyFromClusKey(key), d.freeIndex(TrkrDefs::getHitSetKeyFromClusKey(key)));
        ++d.nTpcDuplicateKey;
      }

      // local coordinates on the surface; local z -> cluster time (inverse of ActsGeometry::getLocalCoords with
      // the crossing of the seed, linear in the time)
      const Acts::Vector3 local = surface->localToGlobalTransform(gctx).inverse() * (onSurface * cm);
      NewCluster nc;
      nc.key = key;
      nc.cluster.setSubSurfKey(subsurfkey);
      nc.cluster.setLocalX(static_cast<float>(local.x() / cm));
      nc.cluster.setLocalY(0.f);
      const double zl0 = d.geometry->getLocalCoords(key, &nc.cluster, crossing)(1);
      nc.cluster.setLocalY(100.f);
      const double zl100 = d.geometry->getLocalCoords(key, &nc.cluster, crossing)(1);
      if (!(std::abs(zl100 - zl0) > 0))
      {
        ++d.nTpcNoSurface;
        continue;
      }
      const double time = 100. * (local.y() / cm - zl0) / (zl100 - zl0);
      nc.cluster.setLocalY(static_cast<float>(time));
      nc.cluster.setPhiError(static_cast<float>(std::hypot(poly->get_rms_x(), poly->get_rms_y())));
      nc.cluster.setZError(static_cast<float>(std::abs(poly->get_rms_z())));
      nc.cluster.setPhiSize(static_cast<char>(std::min(poly->get_phi_width(), 127U)));
      nc.cluster.setZSize(static_cast<char>(std::min(poly->get_time_width(), 127U)));
      nc.cluster.setEdge(0);
      nc.cluster.setOverlap(0);
      nc.cluster.setAdc(static_cast<unsigned int>(std::max(0., poly->get_adc())));
      nc.cluster.setMaxAdc(static_cast<uint16_t>(std::min(65535., std::max(0., poly->get_adc()))));

      // closure: the official global position of the written cluster
      const Acts::Vector3 global = d.geometry->getGlobalPosition(key, &nc.cluster, crossing);
      const double closure = (global - onSurface).norm();
      d.sumMove += move;
      d.maxMove = std::max(d.maxMove, move);
      if (std::isfinite(closure))
      {
        d.sumClosure += closure;
        d.maxClosure = std::max(d.maxClosure, closure);
        ++d.nClosure;
      }
      d.written[key] = onSurface;

      Impl::ClusterRow r;
      r.layer = static_cast<int>(layer);
      r.side = static_cast<int>(side);
      r.crossing = crossing;
      r.x = pos.x();
      r.y = pos.y();
      r.z = pos.z();
      r.sx = onSurface.x();
      r.sy = onSurface.y();
      r.sz = onSurface.z();
      r.gx = global.x();
      r.gy = global.y();
      r.gz = global.z();
      r.move = move;
      r.dr = std::hypot(onSurface.x(), onSurface.y()) - std::hypot(pos.x(), pos.y());
      r.angle = angle;
      r.closure = closure;
      r.time = time;
      tpcRows.push_back(r);
      tpcNew.push_back(std::move(nc));
    }
    if (tpcNew.size() < m_minTpc)
    {
      reject("too few TPC clusters");
      continue;
    }

    // ---- accepted: write the clusters
    for (const auto key : siKeys)
    {
      if (d.outClusters->findCluster(key))
      {
        continue;
      }
      if (const auto nc = std::find_if(siNew.begin(), siNew.end(), [key](const auto& c)
                                       { return c.key == key; });
          nc != siNew.end())
      {
        d.outClusters->addClusterSpecifyKey(key, new TrkrClusterv5(nc->cluster));
        continue;
      }
      TrkrCluster* official = d.siClusters->findCluster(key);
      auto* copy = new TrkrClusterv5;
      copy->CopyFrom(official);
      if (shiftIntt && TrkrDefs::getTrkrId(key) == TrkrDefs::inttId)
      {
        // INTT offset relative to the MVTX: in-plane component on the official surface
        if (const auto surface = d.geometry->maps().getSurface(key, official))
        {
          const Acts::Vector3 dl = surface->localToGlobalTransform(gctx).linear().transpose() * inttShiftMm;
          copy->setLocalX(static_cast<float>(copy->getLocalX() + dl.x() / cm));
          copy->setLocalY(static_cast<float>(copy->getLocalY() + dl.y() / cm));
        }
      }
      d.outClusters->addClusterSpecifyKey(key, copy);
    }
    for (const auto& nc : tpcNew)
    {
      if (!d.outClusters->findCluster(nc.key))
      {
        d.outClusters->addClusterSpecifyKey(nc.key, new TrkrClusterv5(nc.cluster));
      }
    }
    d.nSiProjected += nProjected;

    // ---- seeds: circle in the official frame, z0 and phi at the point closest to the origin
    // (TrackSeedHelper::get_xyz uses that point as the perigee of the official fit)
    const int charge = -h.omega * d.bzSign > 0 ? 1 : -1;
    const double r0 = std::hypot(h.cx, h.cy);  // > 0, checked with the helix
    const double ox = h.cx - h.R * h.cx / r0, oy = h.cy - h.R * h.cy / r0;  // point of the circle closest to (0, 0)
    const double zOrigin = h.z0 + h.tanl * h.arcLength(ox, oy);
    const Acts::Vector3 dirOrigin = h.direction(ox, oy);
    const float seedPhi = static_cast<float>(std::atan2(dirOrigin.y(), dirOrigin.x()));
    auto insertSeed = [&](TrackSeedContainer* container, const std::vector<TrkrDefs::cluskey>& keys) -> unsigned int
    {
      TrackSeed_v2 s;
      s.set_X0(static_cast<float>(h.cx));
      s.set_Y0(static_cast<float>(h.cy));
      s.set_qOverR(static_cast<float>(charge / h.R));  // TrackSeed_v2::get_charge: sign of q/R; get_pt: 1.4 T * R
      s.set_Z0(static_cast<float>(zOrigin));
      s.set_slope(static_cast<float>(h.tanl));
      s.set_phi(seedPhi);
      s.set_crossing(crossing);
      for (const auto key : keys)
      {
        s.insert_cluster_key(key);
      }
      container->insert(&s);
      return static_cast<unsigned int>(container->size() - 1);
    };
    std::vector<TrkrDefs::cluskey> tpcKeys;
    tpcKeys.reserve(tpcNew.size());
    for (const auto& nc : tpcNew)
    {
      tpcKeys.push_back(nc.key);
    }
    const unsigned int siIndex = insertSeed(d.siSeeds, siKeys);
    const unsigned int tpcIndex = insertSeed(d.tpcSeeds, tpcKeys);
    SvtxTrackSeed_v3 svtx;
    svtx.set_silicon_seed_index(siIndex);
    svtx.set_tpc_seed_index(tpcIndex);
    svtx.set_crossing(crossing);
    svtx.set_crossing_estimate(crossing);
    d.svtxSeeds->insert(&svtx);
    const int seedIndex = static_cast<int>(d.svtxSeeds->size() - 1);
    ++d.nSeeds;

    if (d.qaFile)
    {
      for (auto& r : tpcRows)
      {
        r.event = static_cast<int>(d.nEvents - 1);
        r.seed = seedIndex;
        d.clusterRow = r;
        d.clusterTree->Fill();
      }
      auto& s = d.seedRow;
      s.event = static_cast<int>(d.nEvents - 1);
      s.seed = seedIndex;
      s.crossing = crossing;
      s.charge = charge;
      s.nmvtx = static_cast<int>(nmvtx);
      s.nintt = static_cast<int>(nintt);
      s.ntpc = static_cast<int>(tpcNew.size());
      s.nsiProjected = static_cast<int>(nProjected);
      s.pt = static_cast<float>(pt);
      s.phi = seedPhi;
      s.tanl = static_cast<float>(h.tanl);
      s.z0 = static_cast<float>(zOrigin);
      s.X0 = static_cast<float>(h.cx);
      s.Y0 = static_cast<float>(h.cy);
      s.R = static_cast<float>(h.R);
      s.tpcZOffset = static_cast<float>(zOffset);
      d.seedTree->Fill();
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcActsSeeder::End(PHCompositeNode* /*topNode*/)
{
  auto& d = *m_impl;
  const unsigned long nWritten = d.nTpc - d.nTpcAngleCut - d.nTpcMoveCut - d.nTpcNoSurface;
  std::cout << Name() << "::End - events " << d.nEvents << ", input tracks " << d.nInput << ", seeds " << d.nSeeds << std::endl
            << "  Si clusters: official by hit association " << d.nSiOfficialAssoc << " (pruned MVTX duplicate, other strobe "
            << d.nSiOtherStrobe << "), by position " << d.nSiOfficialPosition << ", projected chain points " << d.nSiProjected
            << ", layers without cluster " << d.nSiMissing << std::endl
            << "  TPC clusters " << d.nTpc << ": no surface " << d.nTpcNoSurface << ", crossing angle cut " << d.nTpcAngleCut
            << ", move cut " << d.nTpcMoveCut << ", duplicate keys " << d.nTpcDuplicateKey << std::endl
            << "  TPC move to the official surface: mean " << (nWritten ? 1e4 * d.sumMove / static_cast<double>(nWritten) : 0.)
            << " um, max " << 1e4 * d.maxMove << " um" << std::endl
            << "  closure |official global - intended| : mean " << (d.nClosure ? 1e4 * d.sumClosure / static_cast<double>(d.nClosure) : 0.)
            << " um, max " << 1e4 * d.maxClosure << " um (should be ~0)" << std::endl;
  for (const auto& [why, n] : d.rejected)
  {
    std::cout << "  rejected (" << why << "): " << n << std::endl;
  }
  if (d.qaFile)
  {
    d.qaFile->cd();
    d.clusterTree->Write();
    d.seedTree->Write();
    d.qaFile->Close();
    delete d.qaFile;
    d.qaFile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
