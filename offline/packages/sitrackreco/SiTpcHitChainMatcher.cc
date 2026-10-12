#include "SiTpcHitChainMatcher.h"

#include "SiHitSeedData.h"
#include "SiTpcHitMatchContainerv1.h"
#include "SiTpcHitMatchv1.h"

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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>

namespace
{
  using SiTpcChainMath::wrapPi;
  double square(double x) { return x * x; }
  double clampTheta(double t) { return std::min(M_PI - 1e-4, std::max(1e-4, t)); }
  double cotan(double t) { return std::cos(t) / std::sin(t); }
  double seconds(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b)
  {
    return std::chrono::duration<double>(b - a).count();
  }
}  // namespace

SiTpcHitChainMatcher::SiTpcHitChainMatcher(const std::string& name)
  : SubsysReco(name)
{
}

int SiTpcHitChainMatcher::InitRun(PHCompositeNode* topNode)
{
  m_geometry.build();
  m_map.frame = &m_frame;
  m_map.alignment = &m_alignment;
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": hit model " << (m_geometry.model() == SiHitGeometryTable::Planar ? "planar" : "linear")
              << ", layer order";
    for (int l : m_layerOrder)
    {
      std::cout << " " << l;
    }
    std::cout << "\n  entry  |dphi| < " << m_entryPhi[0] << " / " << m_entryPhi[1] << ", |dtheta| < " << m_entryTheta[0]
              << " / " << m_entryTheta[1] << " rad (MVTX / INTT)"
              << "\n  dd     |ddphi| < " << m_ddPhi[0] << " / " << m_ddPhi[1] << ", |ddtheta| < " << m_ddTheta[0] << " / "
              << m_ddTheta[1] << " rad"
              << "\n  branches/layer " << m_maxBranchesPerLayer[0] << " / " << m_maxBranchesPerLayer[1] << ", chains " << m_maxChains << ", group slack "
              << m_groupSlack << " pitches, overlap groups " << (m_attachOverlap ? "on" : "off") << std::endl;
  }
  return createNodes(topNode);
}

int SiTpcHitChainMatcher::createNodes(PHCompositeNode* topNode)
{
  m_container = findNode::getClass<SiTpcHitMatchContainer>(topNode, m_outputNodeName);
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
  m_container = new SiTpcHitMatchContainerv1();
  svtxNode->addNode(new PHIODataNode<PHObject>(m_container, m_outputNodeName, "PHObject"));
  return Fun4AllReturnCodes::EVENT_OK;
}

// TPC-only helix per TPC track, as in SiTpcTrackMatcher (TPC frame, beam = TPC beam axis).
std::vector<SiTpcHitChainMatcher::TrackInput> SiTpcHitChainMatcher::collectTpcTracks(PHCompositeNode* topNode) const
{
  std::vector<TrackInput> out;
  auto* clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  auto* tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName);
  if (!clusters || !tracks)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_tpcClusterNodeName << " or " << m_tpcTrackNodeName << std::endl;
    }
    return out;
  }
  std::map<unsigned int, std::vector<SiTpcHelixFit::Point>> byTrack;
  for (unsigned int i = 0; i < clusters->size(); ++i)
  {
    const Tpc_PolyCluster* c = clusters->get_cluster(i);
    if (!c || !c->isValid())
    {
      continue;
    }
    const auto p = m_alignment.tpcToTpcFrame(c->get_centroid_x(), c->get_centroid_y(), c->get_centroid_z());
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
    {
      continue;
    }
    byTrack[c->get_source_assembled_track_id()].push_back({p.x, p.y, p.z, 1.0, 1.0, 1});
  }
  SiTpcHelixFit::Config cfg;
  cfg.beamX = m_alignment.referenceBeamX();
  cfg.beamY = m_alignment.referenceBeamY();
  cfg.bz = m_bz;
  for (unsigned int i = 0; i < tracks->size(); ++i)
  {
    const Tpc_PolyTrack* t = tracks->get_track(i);
    if (!t)
    {
      continue;
    }
    auto it = byTrack.find(t->get_source_assembled_track_id());
    if (it == byTrack.end() || it->second.size() < m_minTpcClusters)
    {
      continue;
    }
    auto& pts = it->second;
    std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b)
              { return a.x * a.x + a.y * a.y < b.x * b.x + b.y * b.y; });
    TrackInput in;
    in.index = i;
    in.id = t->get_track_id();
    in.sourceId = t->get_source_assembled_track_id();
    in.fit = SiTpcHelixFit::fit(pts, cfg);
    if (in.fit.status == SiTpcHelixFit::Ok)
    {
      out.push_back(std::move(in));
    }
  }
  return out;
}

int SiTpcHitChainMatcher::process_event(PHCompositeNode* topNode)
{
  ++m_nEvents;
  if (!m_container)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  m_container->Reset();
  m_nextMatchId = 0;

  auto* event = findNode::getClass<SiHitSeedEvent>(topNode, m_hitEventNodeName);
  if (!event)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_hitEventNodeName << " (run SiHitSeedReco first)" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const auto t0 = std::chrono::steady_clock::now();
  m_index.build(event->hits, m_geometry);
  for (int layer = 0; layer < 7; ++layer)
  {
    m_stamp[layer].assign(m_index.size(layer), 0);
    m_groupOfHit[layer].assign(m_index.size(layer), -1);
  }
  m_curStamp = 0;
  const auto t1 = std::chrono::steady_clock::now();
  const auto tracks = collectTpcTracks(topNode);
  const auto t2 = std::chrono::steady_clock::now();

  Counters counters;
  counters.tpcTracks = static_cast<unsigned int>(tracks.size());
  for (const auto& track : tracks)
  {
    processTrack(track, counters);
  }
  const auto t3 = std::chrono::steady_clock::now();
  m_timeIndex += seconds(t0, t1);
  m_timeTpc += seconds(t1, t2);
  m_timeMatch += seconds(t2, t3);

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << m_nEvents << " si_hits=" << m_index.size() << " tpc_tracks=" << counters.tpcTracks
              << " references=" << counters.references << " matched=" << counters.matched
              << " output=" << m_container->size() << " branches=" << counters.branches;
    for (int layer = 0; layer < 7; ++layer)
    {
      std::cout << " L" << layer << "[visit=" << counters.visited[layer] << " win=" << counters.inWindow[layer]
                << " grp=" << counters.groupsBuilt[layer] << " acc=" << counters.groupsAccepted[layer]
                << " att=" << counters.attached[layer] << "]";
    }
    std::cout << " t_index=" << seconds(t0, t1) << " t_tpc=" << seconds(t1, t2) << " t_match=" << seconds(t2, t3)
              << " s" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpcHitChainMatcher::End(PHCompositeNode*)
{
  if (Verbosity() > 0 && m_nEvents > 0)
  {
    std::cout << Name() << ": " << m_nEvents << " events, mean time per event: index " << m_timeIndex / m_nEvents
              << " s, TPC helices " << m_timeTpc / m_nEvents << " s, matching " << m_timeMatch / m_nEvents << " s"
              << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

// Expected footprint of one particle crossing on a sensor of `layer`, from the TPC incidence.
SiTpcHitChainMatcher::Limits SiTpcHitChainMatcher::groupLimits(const SiTpcChainMath::TpcReference& ref,
                                                              const SiTpcChainMath::LayerPrediction& tp) const
{
  Limits lim;
  const int layer = tp.layer;
  const double t = m_geometry.thickness(layer);
  const double alpha = std::asin(std::min(1.0, std::abs(tp.sinA)));
  // MVTX staves are tilted: the incidence on the sensor is bounded by |alpha| + tilt
  const double inc = std::min(1.45, alpha + SiHitGeometryTable::staveTilt(layer));
  lim.expRphi = 1.0 + t * std::tan(inc) / SiHitGeometryTable::pitchRphi(layer);
  lim.maxRphi = std::ceil(lim.expRphi) + m_groupSlack;
  if (SiHitGeometryTable::isMvtx(layer))
  {
    lim.expZ = 1.0 + t * std::abs(ref.tanl) / std::cos(inc) / SiHitGeometryTable::pitchZ(layer);
    lim.maxZ = std::ceil(lim.expZ) + m_groupSlack;
  }
  else
  {
    lim.expZ = 1.0;
    lim.maxZ = m_inttColumnNeighbours ? 2.0 : 1.0;
  }
  return lim;
}

bool SiTpcHitChainMatcher::neighbours(int layer, const SiIndexedHit& a, const SiIndexedHit& b, bool& crossChip) const
{
  crossChip = false;
  if (a.time != b.time)
  {
    return false;
  }
  const int drow = std::abs(static_cast<int>(a.row) - static_cast<int>(b.row));
  const int dcol = std::abs(static_cast<int>(a.col) - static_cast<int>(b.col));
  if (drow > 1)
  {
    return false;
  }
  if (a.hitsetkey == b.hitsetkey)
  {
    if (SiHitGeometryTable::isMvtx(layer))
    {
      return dcol <= 1;
    }
    return dcol == 0 || (m_inttColumnNeighbours && dcol == 1);
  }
  // MVTX: abutting chips of the same stave, facing boundary columns (col grows with z)
  if (SiHitGeometryTable::isMvtx(layer) && a.phiId == b.phiId && std::abs(a.zId - b.zId) == 1)
  {
    const SiIndexedHit& lo = a.zId < b.zId ? a : b;
    const SiIndexedHit& hi = a.zId < b.zId ? b : a;
    if (lo.col == SiHitGeometryTable::kMvtxCols - 1 && hi.col == 0)
    {
      crossChip = true;
      return true;
    }
  }
  return false;
}

// Group (connected neighbour set) containing hit k of `layer`; built once per track and layer.
int SiTpcHitChainMatcher::groupOf(int layer, int k, const Limits& lim, Counters& counters)
{
  if (m_stamp[layer][k] == m_curStamp)
  {
    return m_groupOfHit[layer][k];
  }
  const auto& hits = m_index.hits(layer);
  const double R = SiHitGeometryTable::radius(layer);
  const bool mvtx = SiHitGeometryTable::isMvtx(layer);
  const double pr = SiHitGeometryTable::pitchRphi(layer);
  const double gap = SiHitGeometryTable::mvtxChipGap();

  const int gid = static_cast<int>(m_groups.size());
  m_groups.emplace_back();
  Group& g = m_groups.back();
  g.layer = layer;
  g.hitsetkey = hits[k].hitsetkey;
  g.phiId = hits[k].phiId;
  g.zId = hits[k].zId;
  g.members.push_back(k);
  m_stamp[layer][k] = m_curStamp;
  m_groupOfHit[layer][k] = gid;
  ++counters.groupsBuilt[layer];

  const double u0 = hits[k].u;
  double uMin = 0, uMax = 0, zMin = hits[k].z, zMax = hits[k].z;
  bool anyCrossChip = false;
  std::vector<int> stack{k};
  while (!stack.empty())
  {
    const int q = stack.back();
    stack.pop_back();
    const SiIndexedHit& hq = hits[q];
    const double duN = 1.5 * pr / R;
    double dzN = 0.01;
    if (mvtx)
    {
      const bool edge = hq.col == 0 || hq.col == SiHitGeometryTable::kMvtxCols - 1;
      dzN = 1.5 * SiHitGeometryTable::pitchZ(layer) + (edge ? gap : 0.0);
    }
    else if (m_inttColumnNeighbours)
    {
      dzN = 1.01 * SiHitGeometryTable::pitchZ(layer, hq.zId);
    }
    m_index.forEachInWindow(layer, hq.u, duN, hq.z - dzN, hq.z + dzN, [&](int n)
    {
      if (m_stamp[layer][n] == m_curStamp)
      {
        return;
      }
      const SiIndexedHit& hn = hits[n];
      if (std::abs(wrapPi(hn.u - hq.u)) > duN)
      {
        return;
      }
      bool crossChip = false;
      if (!neighbours(layer, hq, hn, crossChip))
      {
        return;
      }
      const double ur = wrapPi(hn.u - u0);
      const double nuMin = std::min(uMin, ur), nuMax = std::max(uMax, ur);
      const double nzMin = std::min(zMin, static_cast<double>(hn.z)), nzMax = std::max(zMax, static_cast<double>(hn.z));
      const bool cc = anyCrossChip || crossChip;
      const double extR = (nuMax - nuMin) * R / pr + 1.0;
      const double extZ = mvtx ? std::max(0.0, nzMax - nzMin - (cc ? gap : 0.0)) / SiHitGeometryTable::pitchZ(layer) + 1.0
                               : 1.0;
      if (extR > lim.maxRphi + 1e-6 || (mvtx && extZ > lim.maxZ + 1e-6) || g.members.size() >= m_maxGroupHits)
      {
        g.flags |= SiTpcHitMatchStep::Truncated;
        return;
      }
      uMin = nuMin;
      uMax = nuMax;
      zMin = nzMin;
      zMax = nzMax;
      anyCrossChip = cc;
      m_stamp[layer][n] = m_curStamp;
      m_groupOfHit[layer][n] = gid;
      g.members.push_back(n);
      stack.push_back(n);
    });
  }
  // centroid (binary MVTX readout: unweighted)
  double su = 0, sz = 0, sr = 0;
  for (int m : g.members)
  {
    su += wrapPi(hits[m].u - u0);
    sz += hits[m].z;
    sr += hits[m].r;
  }
  const double n = static_cast<double>(g.members.size());
  g.u = SiTpcChainMath::wrap2Pi(u0 + su / n);
  g.z = sz / n;
  g.r = sr / n;
  const double phi = g.u + m_frame.uphiRotation();
  g.x = g.r * std::cos(phi);
  g.y = g.r * std::sin(phi);
  if (anyCrossChip)
  {
    g.flags |= SiTpcHitMatchStep::CrossChip;
  }
  g.extentRphi = static_cast<float>((uMax - uMin) * R / pr + 1.0);
  g.extentZ = mvtx ? static_cast<float>(std::max(0.0, zMax - zMin - (anyCrossChip ? gap : 0.0)) /
                                            SiHitGeometryTable::pitchZ(layer) +
                                        1.0)
                   : 1.f;
  return gid;
}

// Window of branch `br` on the layer of tp: centre from the TPC helix (entry) or from the
// silicon chain built with TPC curvature and theta, half widths, theta window as z range.
SiTpcHitChainMatcher::Window SiTpcHitChainMatcher::makeWindow(const SiTpcChainMath::TpcReference& ref,
                                                              const SiTpcChainMath::LayerPrediction& tp,
                                                              const Branch& br) const
{
  Window w;
  const int layer = tp.layer;
  const int det = detIndex(layer);
  w.tp = &tp;
  w.entry = br.steps.empty();
  w.uRef = tp.u;
  w.sRef = tp.s;
  w.zRef = tp.z;
  w.wPhi = m_entryPhi[det];
  w.wTheta = m_entryTheta[det];
  if (!w.entry)
  {
    const auto sp = SiTpcChainMath::predictSiChain(br.si, ref.R, ref.tanl, tp, m_map);
    if (!sp.valid)
    {
      return w;
    }
    w.uRef = sp.u;
    w.sRef = sp.s;
    w.zRef = sp.z;
    if (br.hasPrev)
    {
      w.meanPhi = m_ddPhiOffset[layer] + m_ddPhiSlope[layer] * br.prevDphi;
      w.meanTheta = m_ddThetaOffset[layer] + m_ddThetaSlope[layer] * br.prevDtheta;
    }
    w.wPhi = m_ddPhi[det];
    w.wTheta = m_ddTheta[det];
  }
  w.thetaRef = SiTpcChainMath::thetaFromVertex(w.sRef, w.zRef, ref.z0);
  w.uc = w.uRef + w.meanPhi;
  const double thc = w.thetaRef + w.meanTheta;
  // theta window -> z range on the nominal cylinder, widened for hits off the nominal radius
  const double dr = m_geometry.maxDr(layer);
  w.zlo = ref.z0 + w.sRef * cotan(clampTheta(thc + w.wTheta)) - std::abs(tp.dzdr) * dr - 1e-3;
  w.zhi = ref.z0 + w.sRef * cotan(clampTheta(thc - w.wTheta)) + std::abs(tp.dzdr) * dr + 1e-3;
  w.valid = true;
  return w;
}

// ONE index query per track and layer: the union of all branch windows.  Every hit inside seeds
// (or joins) a group; the groups are then tested per branch with centroid arithmetic only.
void SiTpcHitChainMatcher::collectGroups(int layer, const std::vector<Window>& windows, const Limits& lim,
                                         std::vector<int>& groups, Counters& counters)
{
  groups.clear();
  const Window* first = nullptr;
  for (const auto& w : windows)
  {
    if (w.valid)
    {
      first = &w;
      break;
    }
  }
  if (!first)
  {
    return;
  }
  const double dr = m_geometry.maxDr(layer);
  const double u0 = first->uc;
  double lo = 1e9, hi = -1e9, zlo = 1e9, zhi = -1e9;
  for (const auto& w : windows)
  {
    if (!w.valid)
    {
      continue;
    }
    const double du = wrapPi(w.uc - u0);
    const double ext = w.wPhi + std::abs(w.tp->dudr) * dr;
    lo = std::min(lo, du - ext);
    hi = std::max(hi, du + ext);
    zlo = std::min(zlo, w.zlo);
    zhi = std::max(zhi, w.zhi);
  }
  const double uc = u0 + 0.5 * (lo + hi);
  const double half = 0.5 * (hi - lo);
  const auto& hits = m_index.hits(layer);
  m_index.forEachInWindow(layer, uc, half, zlo, zhi, [&](int k)
  {
    ++counters.visited[layer];
    if (std::abs(wrapPi(hits[k].u - uc)) > half)
    {
      return;
    }
    ++counters.inWindow[layer];
    const int gid = groupOf(layer, k, lim, counters);
    if (std::find(groups.begin(), groups.end(), gid) == groups.end())
    {
      groups.push_back(gid);
    }
  });
}

// Residuals of group gid w.r.t. the TPC helix and w.r.t. window w; false if outside the window.
// Ranking score from physical resolutions (setScoreSigma / setEntryScoreSigma), not from the
// window widths: the windows only accept.
bool SiTpcHitChainMatcher::evaluateGroup(const SiTpcChainMath::TpcReference& ref, const Window& w, int gid,
                                         const Limits& lim, Candidate& c) const
{
  if (!w.valid)
  {
    return false;
  }
  const auto& tp = *w.tp;
  const int layer = tp.layer;
  const int det = detIndex(layer);
  const Group& g = m_groups[gid];
  const double R = tp.R;
  const double z0T = ref.z0;
  const double drg = g.r - R;

  double dphi = 0, dtheta = 0;
  const double sT = tp.s + tp.dsdr * drg;
  const double zT = tp.z + tp.dzdr * drg;
  const double dphiTpc = wrapPi(g.u - (tp.u + tp.dudr * drg));
  const double dthetaTpc = SiTpcChainMath::thetaFromVertex(sT, g.z, z0T) - ref.theta;
  double sS = sT;
  if (w.entry)
  {
    dphi = dphiTpc;
    dtheta = dthetaTpc;
  }
  else
  {
    sS = w.sRef + tp.dsdr * drg;
    const double zS = w.zRef + ref.tanl * tp.dsdr * drg;
    dphi = wrapPi(g.u - (w.uRef + tp.dudr * drg));
    dtheta = SiTpcChainMath::thetaFromVertex(sS, g.z, z0T) - SiTpcChainMath::thetaFromVertex(sS, zS, z0T);
  }
  const double ddphi = dphi - w.meanPhi;
  const double ddtheta = dtheta - w.meanTheta;
  if (std::abs(ddphi) > w.wPhi || std::abs(ddtheta) > w.wTheta)
  {
    return false;
  }
  // score: r*phi and z resolutions turned into angles at this point
  double sigRphi = w.entry ? m_entrySigRphi[det] : m_sigRphi[det];
  double sigZ = w.entry ? m_entrySigZ[det] : m_sigZ[det];
  if (!SiHitGeometryTable::isMvtx(layer))
  {
    sigZ = std::max(sigZ, SiHitGeometryTable::pitchZ(layer, g.zId) / std::sqrt(12.0));  // INTT strip
  }
  const double sinT = std::sin(w.thetaRef);
  const double sigTheta = std::max(1e-6, sigZ * sinT * sinT / std::max(sS, 1e-3));
  const double sigPhi = std::max(1e-7, sigRphi / R);

  c = Candidate();
  c.group = gid;
  c.score = square(ddphi / sigPhi) + square(ddtheta / sigTheta);
  SiTpcHitMatchStep& s = c.rec;
  s.layer = layer;
  s.n_hits = static_cast<int>(g.members.size());
  s.flags = g.flags | (w.entry ? SiTpcHitMatchStep::Entry : 0);
  s.hitsetkey = g.hitsetkey;
  s.u = static_cast<float>(g.u);
  s.z = static_cast<float>(g.z);
  s.r = static_cast<float>(g.r);
  s.dphi_tpc = static_cast<float>(dphiTpc);
  s.dtheta_tpc = static_cast<float>(dthetaTpc);
  s.dz_tpc = static_cast<float>(g.z - zT);
  s.dphi = static_cast<float>(dphi);
  s.dtheta = static_cast<float>(dtheta);
  s.ddphi = static_cast<float>(ddphi);
  s.ddtheta = static_cast<float>(ddtheta);
  s.chi2 = static_cast<float>(c.score);
  s.extent_rphi = g.extentRphi;
  s.extent_z = g.extentZ;
  s.exp_extent_rphi = static_cast<float>(lim.expRphi);
  s.exp_extent_z = static_cast<float>(lim.expZ);
  return true;
}

bool SiTpcHitChainMatcher::extendBranch(const SiTpcChainMath::TpcReference& ref, const SiTpcChainMath::LayerPrediction& tp,
                                        const Branch& br, const Candidate& main, const std::vector<Candidate>& overlaps,
                                        Branch& out) const
{
  out = br;
  const bool mvtx = SiHitGeometryTable::isMvtx(tp.layer);
  const double wz = mvtx ? m_zWeightMvtx : m_zWeightIntt;
  auto addGroup = [&](const Candidate& c, bool overlap)
  {
    const Group& g = m_groups[c.group];
    Step st;
    st.group = c.group;
    st.rec = c.rec;
    if (overlap)
    {
      st.rec.flags |= SiTpcHitMatchStep::Overlap;
    }
    out.steps.push_back(st);
    if (overlap)
    {
      return;  // recorded, not used in the chain (the main group of the layer anchors it)
    }
    // chain points live in the TPC frame (clamshell half of the group)
    const auto t = m_map.detToTpc(g.x, g.y, g.z, m_map.half(g.u));
    SiTpcChainMath::SiChainPoint p;
    p.x = t.x;
    p.y = t.y;
    p.z = t.z;
    p.r = g.r;
    p.wz = wz;
    out.points.push_back(p);
  };
  addGroup(main, false);
  for (const auto& o : overlaps)
  {
    addGroup(o, true);
  }

  const double bx = m_alignment.tpcBeamAxis().x(ref.z0);
  const double by = m_alignment.tpcBeamAxis().y(ref.z0);
  SiTpcChainMath::SiChainFit fit;
  if (!SiTpcChainMath::fitSiChain(out.points, ref.R, ref.helicity, ref.tanl, bx, by, &br.si, fit))
  {
    return false;
  }
  out.si = fit;
  out.chi2 += main.score;
  if (mvtx)
  {
    ++out.nMvtx;
  }
  else
  {
    ++out.nIntt;
  }
  if (!br.steps.empty())
  {
    out.hasPrev = true;
    out.prevDphi = main.rec.dphi;
    out.prevDtheta = main.rec.dtheta;
  }
  return true;
}

void SiTpcHitChainMatcher::processTrack(const TrackInput& track, Counters& counters)
{
  SiTpcChainMath::TpcReference ref;
  if (!SiTpcChainMath::buildTpcReference(track.fit, m_map, ref, m_boundaryMargin))
  {
    return;
  }
  ++counters.references;
  m_groups.clear();

  std::vector<Branch> branches(1);
  std::vector<Branch> next;
  std::vector<Candidate> cands;
  std::vector<Window> windows;
  std::vector<int> groups;
  for (int layer : m_layerOrder)
  {
    if (layer < 0 || layer >= 7)
    {
      continue;
    }
    const auto& tp = ref.layers[layer];
    if (!tp.valid || std::abs(tp.z) > SiHitGeometryTable::halfLengthZ(layer) + m_zAcceptanceMargin ||
        m_index.size(layer) == 0)
    {
      continue;  // outside the layer: neither searched nor counted as missing
    }
    ++m_curStamp;  // new group cache for this track and layer
    const Limits lim = groupLimits(ref, tp);
    const auto& alt = ref.alt[layer];
    const unsigned int maxBranches = m_maxBranchesPerLayer[detIndex(layer)];

    // windows of all branches (and, at a clamshell boundary, with the other half's translation)
    windows.clear();
    for (const auto& br : branches)
    {
      windows.push_back(makeWindow(ref, tp, br));
      windows.push_back(alt.valid ? makeWindow(ref, alt, br) : Window());
    }
    collectGroups(layer, windows, lim, groups, counters);

    next.clear();
    for (std::size_t ib = 0; ib < branches.size(); ++ib)
    {
      const Branch& br = branches[ib];
      cands.clear();
      for (int gid : groups)
      {
        Candidate best, c;
        bool ok = false;
        for (int v = 0; v < 2; ++v)
        {
          if (evaluateGroup(ref, windows[2 * ib + v], gid, lim, c) && (!ok || c.score < best.score))
          {
            best = c;
            ok = true;
          }
        }
        if (ok)
        {
          ++counters.groupsAccepted[layer];
          cands.push_back(best);
        }
      }
      std::sort(cands.begin(), cands.end(), [](const Candidate& x, const Candidate& y) { return x.score < y.score; });

      bool extended = false;
      for (std::size_t i = 0; i < cands.size() && i < maxBranches; ++i)
      {
        // overlap groups: a second measurement of the same particle on the neighbouring stave /
        // ladder of this layer.  Compatibility is required with the MAIN group itself (not only
        // with the window): adjacent phi segment, |d(r*phi)| and |dz| below setOverlapTolerance.
        std::vector<Candidate> overlaps;
        if (m_attachOverlap)
        {
          const Group& gm = m_groups[cands[i].group];
          const int nseg = SiHitGeometryTable::nPhiSegments(layer);
          const int det = detIndex(layer);
          std::vector<int> usedPhi{gm.phiId};
          for (std::size_t j = 0; j < cands.size(); ++j)
          {
            const Group& go = m_groups[cands[j].group];
            const int dseg = ((go.phiId - gm.phiId) % nseg + nseg) % nseg;
            if (j == i || (dseg != 1 && dseg != nseg - 1) ||
                std::find(usedPhi.begin(), usedPhi.end(), go.phiId) != usedPhi.end() ||
                std::abs(wrapPi(go.u - gm.u)) * tp.R > m_overlapTolRphi || std::abs(go.z - gm.z) > m_overlapTolZ[det])
            {
              continue;
            }
            usedPhi.push_back(go.phiId);
            overlaps.push_back(cands[j]);
          }
        }
        Branch ext;
        if (extendBranch(ref, tp, br, cands[i], overlaps, ext))
        {
          next.push_back(std::move(ext));
          extended = true;
        }
      }
      if (!extended)
      {
        Branch miss = br;
        if (!miss.steps.empty())
        {
          ++miss.nMissing;
          miss.chi2 += m_missingPenalty;
        }
        next.push_back(std::move(miss));
      }
    }
    std::stable_sort(next.begin(), next.end(), [](const Branch& a, const Branch& b)
                     {
                       if (a.nLayers() != b.nLayers()) return a.nLayers() > b.nLayers();
                       return a.chi2 < b.chi2;
                     });
    // combinatorics control (not a residual cut): drop a branch if another branch with at least
    // as many layers has a chi2 lower by more than the margin (next is sorted by layers, chi2)
    double bestChi2 = 1e300;
    std::size_t kept = 0;
    for (std::size_t ib = 0; ib < next.size(); ++ib)
    {
      if (next[ib].chi2 > bestChi2 + m_branchScoreMargin)
      {
        continue;
      }
      bestChi2 = std::min(bestChi2, next[ib].chi2);
      if (kept != ib)
      {
        next[kept] = std::move(next[ib]);
      }
      ++kept;
    }
    next.resize(kept);
    if (next.size() > m_maxChains)
    {
      next.resize(m_maxChains);
    }
    counters.branches += next.size();
    branches.swap(next);
  }

  std::vector<const Branch*> good;
  for (const auto& br : branches)
  {
    if (!br.steps.empty() && br.nLayers() >= m_minLayers && br.nMvtx >= m_minMvtxLayers)
    {
      good.push_back(&br);
    }
  }
  if (good.empty())
  {
    return;
  }
  std::stable_sort(good.begin(), good.end(), [](const Branch* a, const Branch* b)
                   {
                     if (a->nLayers() != b->nLayers()) return a->nLayers() > b->nLayers();
                     if (a->nMvtx != b->nMvtx) return a->nMvtx > b->nMvtx;
                     return a->chi2 < b->chi2;
                   });
  ++counters.matched;
  for (const auto& st : good.front()->steps)
  {
    ++counters.attached[st.rec.layer];
  }
  const std::size_t nStore = m_storeAll ? good.size() : 1;
  for (std::size_t i = 0; i < nStore; ++i)
  {
    storeBranch(track, ref, *good[i], static_cast<unsigned int>(i));
  }
}

void SiTpcHitChainMatcher::storeBranch(const TrackInput& track, const SiTpcChainMath::TpcReference& ref, const Branch& br,
                                       unsigned int rank)
{
  auto* m = new SiTpcHitMatchv1();
  m->set_id(m_nextMatchId++);
  m->set_tpc_track_index(track.index);
  m->set_tpc_track_id(track.id);
  m->set_tpc_assembled_track_id(track.sourceId);
  m->set_crossing(0);
  m->set_rank(rank);
  m->set_tpc_pt(static_cast<float>(ref.pt));
  m->set_tpc_charge(ref.charge);
  m->set_tpc_phi(static_cast<float>(ref.phi));
  m->set_tpc_tanl(static_cast<float>(ref.tanl));
  m->set_tpc_z0(static_cast<float>(ref.z0));
  m->set_tpc_dca(static_cast<float>(ref.dca));
  m->set_n_mvtx_layers(br.nMvtx);
  m->set_n_intt_layers(br.nIntt);
  m->set_n_missing(br.nMissing);
  m->set_chi2(static_cast<float>(br.chi2));
  m->set_si_z0(static_cast<float>(br.si.z0));
  m->set_si_dca(static_cast<float>(br.si.dca));
  m->set_si_circle_x(static_cast<float>(br.si.cx));
  m->set_si_circle_y(static_cast<float>(br.si.cy));
  m->set_si_radius(static_cast<float>(ref.R));
  m->set_si_helicity(ref.helicity);
  m->set_si_pca_x(static_cast<float>(br.si.pcaX));
  m->set_si_pca_y(static_cast<float>(br.si.pcaY));

  // TPC-Si compatibility at the outermost attached (non-overlap) layer
  const Step* outer = nullptr;
  for (const auto& st : br.steps)
  {
    if (!(st.rec.flags & SiTpcHitMatchStep::Overlap) && (!outer || st.rec.r > outer->rec.r))
    {
      outer = &st;
    }
  }
  if (outer)
  {
    const auto& tp = ref.layers[outer->rec.layer];
    const auto sp = SiTpcChainMath::predictSiChain(br.si, ref.R, ref.tanl, tp, m_map);
    if (sp.valid)
    {
      double tx = 0, ty = 0;
      SiTpcChainMath::tangent(br.si.cx, br.si.cy, ref.helicity, sp.xt, sp.yt, tx, ty);
      m->set_delta_phi_dir(static_cast<float>(std::atan2(tp.tx * ty - tp.ty * tx, tp.tx * tx + tp.ty * ty)));
    }
    const double tpcDca = std::hypot(br.si.bx - ref.cx, br.si.by - ref.cy) - ref.R;  // same beam point, TPC frame
    m->set_delta_dca(static_cast<float>(br.si.dca - tpcDca));
  }
  m->set_delta_z0(static_cast<float>(br.si.z0 - ref.z0));

  for (std::size_t i = 0; i < br.steps.size(); ++i)
  {
    const auto& st = br.steps[i];
    m->add_step(st.rec);
    const auto& hits = m_index.hits(st.rec.layer);
    for (int k : m_groups[st.group].members)
    {
      m->add_hit_index(hits[k].hitsetkey, hits[k].hitkey, static_cast<unsigned int>(i));
    }
  }
  m_container->add_match(m);
}