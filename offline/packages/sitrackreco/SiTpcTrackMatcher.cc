#include "SiTpcTrackMatcher.h"

#include "SiHitSeedData.h"
#include "SiTpcHelixFit.h"
#include "SiTpc_TrackContainerv1.h"
#include "SiTpc_Trackv1.h"
#include "Si_Trajectory.h"
#include "Si_TrajectoryContainer.h"

#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <trackbase/TrkrDefs.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace
{
  constexpr double kPi = 3.14159265358979323846;

  double wrapPi(double a)
  {
    a = std::fmod(a + kPi, 2.0 * kPi);
    if (a < 0)
    {
      a += 2.0 * kPi;
    }
    return a - kPi;
  }

  // detector index for the weights: 0 = MVTX, 1 = INTT, 2 = TPC
  int detIndex(int layer)
  {
    return layer <= 2 ? 0 : (layer <= 6 ? 1 : 2);
  }

  struct SiInfo
  {
    unsigned int index = 0;
    const Si_Trajectory* t = nullptr;
    bool fitted = false;
    double phi = 0, eta = 0, z0 = 0;
    int charge = 0;  // 0 for straight-line fits
  };

  struct TpcPoint
  {
    unsigned int clusterIndex = 0;
    int layer = -1;
    double x = 0, y = 0, z = 0, r = 0;
  };

  struct TpcInfo
  {
    unsigned int index = 0;
    const Tpc_PolyTrack* t = nullptr;
    std::vector<TpcPoint> points;  // sorted outward
    SiTpcHelixFit::Result fit;
    double eta = 0;
  };

  struct Candidate
  {
    std::size_t si = 0, tpc = 0;  // indices in the local vectors
    double dphi = 0, deta = 0, dz0 = 0, chi2 = 0;
    int rank = 0;  // chain rank (FromSiChains)
    bool eligible = false;
    std::size_t store = 0;  // index in the container's candidate list
  };
}  // namespace

SiTpcTrackMatcher::SiTpcTrackMatcher(const std::string& name)
  : SubsysReco(name)
{
}

int SiTpcTrackMatcher::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": association "
              << (m_association == Association::FromSiChains
                      ? std::string("from the TPC-driven Si chains (") + m_siSeedNodeName + ", selected chains" +
                            (m_windowOnChains ? " inside the windows)" : ")")
                      : std::string("phi/eta windows"))
              << "; windows |dphi| < " << m_dphiWindow << " rad, |deta| < " << m_detaWindow;
    if (m_dz0Window > 0)
    {
      std::cout << ", |dz0| < " << m_dz0Window << " cm";
    }
    std::cout << (m_sameCharge ? ", same charge" : "") << "; refit weights xy (MVTX/INTT/TPC) " << m_wxy[0] << "/"
              << m_wxy[1] << "/" << m_wxy[2] << ", z " << m_wz[0] << "/" << m_wz[1] << "/" << m_wz[2]
              << ", TPC z offset " << (m_fitTpcZOffset ? "free" : "fixed") << ", Si halves on TPC beam axis "
              << (m_alignment.enabled() ? "on" : "off") << ", beam (" << beamX() << ", " << beamY() << ") cm" << std::endl;
  }
  return createNodes(topNode);
}

int SiTpcTrackMatcher::createNodes(PHCompositeNode* topNode)
{
  m_container = findNode::getClass<SiTpc_TrackContainer>(topNode, m_outputNodeName);
  if (m_container)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << ": DST node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator dstIter(dstNode);
  auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }
  m_container = new SiTpc_TrackContainerv1();
  svtxNode->addNode(new PHIODataNode<PHObject>(m_container, m_outputNodeName, "PHObject"));
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpcTrackMatcher::process_event(PHCompositeNode* topNode)
{
  ++m_nEvents;
  if (!m_container)
  {
    m_container = findNode::getClass<SiTpc_TrackContainer>(topNode, m_outputNodeName);
    if (!m_container)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  m_container->Reset();

  auto* siTraj = findNode::getClass<Si_TrajectoryContainer>(topNode, m_siTrajNodeName);
  auto* clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  auto* tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName);
  const bool fromChains = m_association == Association::FromSiChains;
  auto* seeds = fromChains ? findNode::getClass<SiHitSeedEvent>(topNode, m_siSeedNodeName) : nullptr;
  if (!siTraj || !clusters || !tracks || (fromChains && !seeds))
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing input (" << m_siTrajNodeName << " " << (siTraj ? "ok" : "missing") << ", "
                << m_tpcClusterNodeName << " " << (clusters ? "ok" : "missing") << ", " << m_tpcTrackNodeName << " "
                << (tracks ? "ok" : "missing");
      if (fromChains)
      {
        std::cout << ", " << m_siSeedNodeName << " " << (seeds ? "ok" : "missing");
      }
      std::cout << ")" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // ---- silicon side (chain mode: unfitted trajectories are kept, their points still go into
  // the combined refit; the matching variables are then NaN)
  std::vector<SiInfo> si;
  for (unsigned int i = 0; i < siTraj->size(); ++i)
  {
    const Si_Trajectory* t = siTraj->get_trajectory(i);
    if (!t || (!fromChains && !t->is_fitted()) || t->size_points() == 0)
    {
      continue;
    }
    SiInfo s;
    s.index = i;
    s.t = t;
    s.fitted = t->is_fitted();
    s.phi = s.fitted ? t->get_phi() : std::numeric_limits<double>::quiet_NaN();
    s.eta = s.fitted ? std::asinh(t->get_tanl()) : std::numeric_limits<double>::quiet_NaN();
    s.z0 = s.fitted ? t->get_z0() : std::numeric_limits<double>::quiet_NaN();
    s.charge = t->get_fit_status() == Si_Trajectory::Ok ? t->get_charge() : 0;
    si.push_back(s);
  }

  // ---- TPC side: clusters per assembled track, TPC frame, TPC-only helix
  std::map<unsigned int, std::vector<TpcPoint>> clustersById;
  for (unsigned int i = 0; i < clusters->size(); ++i)
  {
    const Tpc_PolyCluster* c = clusters->get_cluster(i);
    if (!c || !c->isValid())
    {
      continue;
    }
    // TPC = reference frame: clusters as measured
    const auto b = m_alignment.tpcToTpcFrame(c->get_centroid_x(), c->get_centroid_y(), c->get_centroid_z());
    if (!std::isfinite(b.x) || !std::isfinite(b.y) || !std::isfinite(b.z))
    {
      continue;
    }
    TpcPoint p;
    p.clusterIndex = i;
    p.layer = c->size_hits() > 0 ? static_cast<int>(TrkrDefs::getLayer(c->get_hit_index(0).first)) : -1;
    p.x = b.x;
    p.y = b.y;
    p.z = b.z;
    p.r = std::hypot(b.x, b.y);
    clustersById[c->get_source_assembled_track_id()].push_back(p);
  }

  SiTpcHelixFit::Config tpcCfg;
  tpcCfg.beamX = beamX();
  tpcCfg.beamY = beamY();
  tpcCfg.bz = m_bz;
  std::vector<TpcInfo> tpc;
  // assembled track id -> position in tpc (the Si chains may come from a different TPC track
  // container, e.g. TPC_TRIPLETTRACKS, so they are linked by assembled id, not by index)
  std::map<unsigned int, std::size_t> tpcByAssembledId;
  for (unsigned int i = 0; i < tracks->size(); ++i)
  {
    const Tpc_PolyTrack* t = tracks->get_track(i);
    if (!t)
    {
      continue;
    }
    const auto it = clustersById.find(t->get_source_assembled_track_id());
    if (it == clustersById.end() || it->second.size() < m_minTpcClusters)
    {
      continue;
    }
    TpcInfo info;
    info.index = i;
    info.t = t;
    info.points = it->second;
    std::sort(info.points.begin(), info.points.end(), [](const auto& a, const auto& b) { return a.r < b.r; });
    std::vector<SiTpcHelixFit::Point> pts;
    for (const auto& p : info.points)
    {
      pts.push_back({p.x, p.y, p.z, 1.0, 1.0, 1});
    }
    info.fit = SiTpcHelixFit::fit(pts, tpcCfg);
    if (!info.fit.fitted())
    {
      continue;
    }
    info.eta = std::asinh(info.fit.tanl);
    tpcByAssembledId[t->get_source_assembled_track_id()] = tpc.size();
    tpc.push_back(std::move(info));
  }

  // ---- candidates
  const bool useDz0 = m_dz0Window > 0;
  std::vector<Candidate> cands;
  auto makeCandidate = [&](std::size_t a, std::size_t b, Candidate& c)
  {
    c.si = a;
    c.tpc = b;
    c.dphi = wrapPi(si[a].phi - tpc[b].fit.phi);
    c.deta = si[a].eta - tpc[b].eta;
    c.dz0 = si[a].z0 - tpc[b].fit.z0;
    const double ndphi = c.dphi / m_dphiWindow, ndeta = c.deta / m_detaWindow;
    const double ndz0 = useDz0 ? c.dz0 / m_dz0Window : 0.0;
    c.chi2 = ndphi * ndphi + ndeta * ndeta + ndz0 * ndz0;
    const bool chargeOk = !m_sameCharge || si[a].charge == 0 || si[a].charge == tpc[b].fit.charge;
    const bool inWindow = std::abs(ndphi) <= 1 && std::abs(ndeta) <= 1 && std::abs(ndz0) <= 1 && chargeOk;
    const bool inCandWindow = std::abs(ndphi) <= m_candScale && std::abs(ndeta) <= m_candScale && std::abs(ndz0) <= m_candScale;
    return std::make_pair(inWindow, inCandWindow);
  };

  if (fromChains)
  {
    // every trajectory belongs to the TPC track of its chain
    for (std::size_t a = 0; a < si.size(); ++a)
    {
      const int cid = si[a].t->get_chain_id();
      if (cid < 0 || cid >= static_cast<int>(seeds->chains.size()) || seeds->chains[cid].id != cid)
      {
        continue;
      }
      const SiHitChain& chain = seeds->chains[cid];
      const auto it = tpcByAssembledId.find(chain.tpc_assembled_track_id);
      if (chain.tpc_track_index < 0 || it == tpcByAssembledId.end())
      {
        continue;
      }
      Candidate c;
      const auto [inWindow, inCandWindow] = makeCandidate(a, it->second, c);
      (void) inCandWindow;  // chain pairs are always stored
      c.rank = chain.rank;
      c.eligible = chain.selected && (!m_windowOnChains || inWindow);
      c.store = m_container->size_candidates();
      m_container->add_candidate(si[a].index, tpc[it->second].index, c.dphi, c.deta, c.dz0, c.chi2, 0);
      cands.push_back(c);
    }
  }
  else
  {
    for (std::size_t a = 0; a < si.size(); ++a)
    {
      for (std::size_t b = 0; b < tpc.size(); ++b)
      {
        Candidate c;
        const auto [inWindow, inCandWindow] = makeCandidate(a, b, c);
        if (!inCandWindow)
        {
          continue;
        }
        c.eligible = inWindow;
        c.store = m_container->size_candidates();
        m_container->add_candidate(si[a].index, tpc[b].index, c.dphi, c.deta, c.dz0, c.chi2, 0);
        cands.push_back(c);
      }
    }
  }

  // ---- assignment
  std::vector<const Candidate*> order;
  for (const auto& c : cands)
  {
    if (c.eligible)
    {
      order.push_back(&c);
    }
  }
  std::sort(order.begin(), order.end(), [](const auto* a, const auto* b)
  {
    if (a->rank != b->rank)
    {
      return a->rank < b->rank;
    }
    const bool fa = std::isfinite(a->chi2), fb = std::isfinite(b->chi2);
    if (fa != fb)
    {
      return fa;
    }
    return a->chi2 < b->chi2;
  });
  std::vector<char> siUsed(si.size(), 0), tpcUsed(tpc.size(), 0);

  SiTpcHelixFit::Config cfg;
  cfg.beamX = beamX();
  cfg.beamY = beamY();
  cfg.bz = m_bz;
  cfg.groupZOffset = m_fitTpcZOffset;
  cfg.beamWeight = m_beamConstraint ? m_beamWeight : 0.0;

  for (const Candidate* c : order)
  {
    // one Si+TPC track per TPC track; a Si trajectory is used once (window mode: its hits are
    // exclusive; chain mode: one trajectory per chain, chains of different TPC tracks may share hits)
    if (siUsed[c->si] || tpcUsed[c->tpc])
    {
      continue;
    }
    siUsed[c->si] = tpcUsed[c->tpc] = 1;
    m_container->set_candidate_accepted(c->store, 1);
    const SiInfo& s = si[c->si];
    const TpcInfo& t = tpc[c->tpc];

    auto* track = new SiTpc_Trackv1();
    track->set_id(m_container->size());
    track->set_si_trajectory_index(s.index);
    track->set_si_chain_id(s.t->get_chain_id());
    track->set_tpc_track_index(t.index);
    track->set_tpc_track_id(t.t->get_track_id());
    track->set_tpc_assembled_track_id(t.t->get_source_assembled_track_id());
    track->set_si_phi(s.phi);
    track->set_si_eta(s.eta);
    track->set_si_z0(s.z0);
    track->set_tpc_phi(t.fit.phi);
    track->set_tpc_eta(t.eta);
    track->set_tpc_z0(t.fit.z0);
    track->set_tpc_pt(t.fit.pt);
    track->set_tpc_charge(t.fit.charge);
    track->set_match_dphi(c->dphi);
    track->set_match_deta(c->deta);
    track->set_match_dz0(c->dz0);
    track->set_match_chi2(c->chi2);

    // combined points: Si layer points (inner -> outer), then TPC clusters (sorted outward)
    std::vector<SiTpcHelixFit::Point> pts;
    for (unsigned int k = 0; k < s.t->size_points(); ++k)
    {
      const int layer = s.t->get_point_layer(k);
      const int d = detIndex(layer);
      pts.push_back({s.t->get_point_x(k), s.t->get_point_y(k), s.t->get_point_z(k), m_wxy[d], m_wz[d], 0});
      track->add_point(SiTpc_Track::SiPoint, layer, k, s.t->get_point_x(k), s.t->get_point_y(k), s.t->get_point_z(k));
    }
    for (const auto& p : t.points)
    {
      pts.push_back({p.x, p.y, p.z, m_wxy[2], m_wz[2], 1});
      track->add_point(SiTpc_Track::TpcPoint, p.layer, p.clusterIndex, p.x, p.y, p.z);
    }
    track->set_n_si_points(s.t->size_points());
    track->set_n_tpc_points(t.points.size());

    const auto r = SiTpcHelixFit::fit(pts, cfg);
    track->set_fit_status(r.status);
    if (r.fitted())
    {
      track->set_circle_x(r.cx);
      track->set_circle_y(r.cy);
      track->set_radius(r.R);
      track->set_helicity(r.helicity);
      track->set_circle_rms(r.circleRms);
      track->set_pca_x(r.pcaX);
      track->set_pca_y(r.pcaY);
      track->set_phi(r.phi);
      track->set_dca(r.dca);
      track->set_pt(r.pt);
      track->set_charge(r.charge);
      track->set_z0(r.z0);
      track->set_tanl(r.tanl);
      track->set_tpc_z_offset(r.zOffset);
      track->set_z_rms(r.zRms);
    }
    if (Verbosity() > 2)
    {
      track->identify();
    }
    m_container->add_track(track);
  }
  m_nMatched += m_container->size();

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << m_nEvents << " Si trajectories " << si.size() << ", TPC tracks " << tpc.size()
              << ", candidates " << m_container->size_candidates() << " (eligible " << order.size() << ")"
              << ", Si+TPC tracks " << m_container->size() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
