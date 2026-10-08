#include "SiHitSeedReco.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>

#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <tuple>
#include <unordered_map>

namespace
{
constexpr double pi = 3.14159265358979323846;
constexpr double twopi = 2.0 * pi;
constexpr std::array<int, 3> nstaves{{12, 16, 20}};
constexpr std::array<double, 3> mvtx_r{{2.523, 3.336, 4.148}};
constexpr std::array<double, 3> mvtx_phi0{{0.2382, 0.1937, 0.1481}};
constexpr std::array<double, 3> mvtx_tilt{{13.3 * pi / 180., 16.9 * pi / 180., 17.0 * pi / 180.}};
constexpr std::array<int, 4> intt_nladder{{12, 12, 16, 16}};
constexpr std::array<double, 4> intt_phi0{{-15 * pi / 180., 0., -11.25 * pi / 180., 0.}};
constexpr std::array<double, 4> intt_r{{7.188 - 0.0036, 7.732 - 0.0036, 9.680 - 0.0036, 10.262 - 0.0036}};
constexpr std::array<double, 9> mvtx_chip_z{{-12.06, -9.045, -6.03, -3.015, 0., 3.015, 6.03, 9.045, 12.06}};
constexpr double intt_inner_z = (8 * 1.6 + 0.2) / 2.;
constexpr double intt_outer_z = (8 * 1.6 + 0.2) + (5 * 2.0 + 0.2) / 2.;
constexpr std::array<double, 4> intt_sensor_z{{-intt_inner_z, -intt_outer_z, intt_inner_z, intt_outer_z}};
constexpr int nlayers = 7;
}  // namespace

SiHitSeedReco::SiHitSeedReco(const std::string& name)
  : SubsysReco(name)
{
}

int SiHitSeedReco::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": seed z window ";
    if (m_zMode == ZWindowMode::Band)
    {
      std::cout << "band +-" << m_zBandHalf << " bins, center " << m_zBandOffset << " - "
                << m_zBandShiftAtEdge << "*z/Zmax bins";
    }
    else if (m_zMode == ZWindowMode::Centimetre)
    {
      std::cout << "dz_cm in (" << m_zInterceptCm << " + " << m_zSlopeCm << "*z) +- " << m_zHalfCm << " cm";
    }
    else if (m_zMode == ZWindowMode::Vertex)
    {
      std::cout << "line through the event vertex, offset " << m_zVtxOffsetCm << " cm, +- " << m_zVtxHalfCm << " cm";
    }
    else
    {
      std::cout << "(" << m_zIntercept << " + " << m_zSlope << "*z) +- " << m_zHalf << " bins";
    }
    std::cout << "; seed phi +-" << m_seedPhiHalf << " bins; propagation phi +-" << m_propPhiHalf
              << ", z +-" << m_propZHalf << " bins" << std::endl;
  }
  return createNodes(topNode);
}

int SiHitSeedReco::createNodes(PHCompositeNode* topNode)
{
  m_eventData = findNode::getClass<SiHitSeedEvent>(topNode, m_outputNodeName);
  if (m_eventData)
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

  m_eventData = new SiHitSeedEvent();
  auto* node = new PHDataNode<PHObject>(m_eventData, m_outputNodeName, "PHObject");
  dstNode->addNode(node);

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": created transient node " << m_outputNodeName << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SiHitSeedReco::process_event(PHCompositeNode* topNode)
{
  if (!m_eventData)
  {
    m_eventData = findNode::getClass<SiHitSeedEvent>(topNode, m_outputNodeName);
    if (!m_eventData)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  m_eventData->Reset();
  m_eventData->event = m_event++;
  if (!fillHits(topNode))
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  fillClusters(topNode);

  if (!m_eventData->hits.empty())
  {
    const auto [zlo, zhi] = std::minmax_element(m_eventData->hits.begin(), m_eventData->hits.end(),
                                                [](const auto& a, const auto& b) { return a.z < b.z; });
    const double span = std::max(0.01, zhi->z - zlo->z);
    m_zmin = zlo->z - std::max(0.02 * span, 0.01);
    m_zmax = zhi->z + std::max(0.02 * span, 0.01);
  }

  const double dz = (m_zmax - m_zmin) / m_nz;
  const double dphi = twopi / m_nphi;
  for (auto& h : m_eventData->hits)
  {
    h.zbin = (h.z - m_zmin) / dz;
    h.phibin = h.phi / dphi;
  }

  buildBlobs();
  buildHitGrid();
  findTrackletVertex();
  findChains();
  fitVertexFromSeedLinks();

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << m_eventData->event
              << " hits=" << m_eventData->hits.size()
              << " clusters=" << m_eventData->clusters.size()
              << " blobs=" << m_eventData->blobs.size()
              << " chains=" << m_eventData->chains.size()
              << " vertex_z=" << m_eventData->vertex_z << " cm (tracklet, " << m_eventData->vertex_tracklet_npeak
              << "/" << m_eventData->vertex_tracklet_npairs << " pairs in peak) vertex_z_linefit="
              << m_eventData->vertex_z_linefit << " cm" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiHitSeedReco::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

double SiHitSeedReco::wrapDelta(double x, double period)
{
  x = std::fmod(x + period / 2., period);
  if (x < 0)
  {
    x += period;
  }
  return x - period / 2.;
}

double SiHitSeedReco::circularMean(const std::vector<double>& values)
{
  double sx = 0, sy = 0;
  for (double p : values)
  {
    sx += std::cos(p);
    sy += std::sin(p);
  }
  double p = std::atan2(sy, sx);
  if (p < 0)
  {
    p += twopi;
  }
  return p;
}

std::pair<double, double> SiHitSeedReco::localCoordinates(int layer, unsigned int row, unsigned int col, int ladderz) const
{
  if (layer <= 2)
  {
    constexpr double row_pitch = 26.88e-4, col_pitch = 29.24e-4, sensor_x = 0.058128;
    constexpr double first_row = 0.5 * (512 * row_pitch - 37.44e-4 + 0.12 - row_pitch);
    return {first_row - row * row_pitch - sensor_x, (col + 0.5 - 512.) * col_pitch};
  }
  const int type = ladderz % 2;
  const int ncols = type == 0 ? 8 : 5;
  const double pitchz = type == 0 ? 1.6 : 2.0;
  return {(row + 0.5 - 128.) * 0.0078, (col + 0.5 - ncols / 2.) * pitchz};
}

double SiHitSeedReco::intrinsicPhi(int layer, int stave, int ladderphi, double lx) const
{
  double p = 0;
  if (layer <= 2)
  {
    p = mvtx_phi0[layer] - mvtx_phi0[2] + stave * twopi / nstaves[layer] + lx * std::cos(mvtx_tilt[layer]) / mvtx_r[layer];
  }
  else
  {
    const int i = layer - 3;
    p = intt_phi0[i] + twopi * ladderphi / intt_nladder[i] + lx / intt_r[i] - mvtx_phi0[2];
  }
  p = std::fmod(p, twopi);
  if (p < 0)
  {
    p += twopi;
  }
  return p;
}

double SiHitSeedReco::longitudinalZ(int layer, int chip, int ladderz, double ly) const
{
  return (layer <= 2 ? mvtx_chip_z.at(chip) : intt_sensor_z.at(ladderz)) + ly;
}

bool SiHitSeedReco::fillHits(PHCompositeNode* topNode)
{
  auto* cont = findNode::getClass<TrkrHitSetContainer>(topNode, m_hitNodeName);
  if (!cont)
  {
    std::cerr << Name() << ": missing " << m_hitNodeName << std::endl;
    return false;
  }

  for (auto det : {TrkrDefs::mvtxId, TrkrDefs::inttId})
  {
    auto range = cont->getHitSets(det);
    for (auto hsit = range.first; hsit != range.second; ++hsit)
    {
      const auto hskey = hsit->first;
      auto* hs = hsit->second;
      const int layer = TrkrDefs::getLayer(hskey);
      if (layer < 0 || layer > 6)
      {
        continue;
      }
      const int stave = layer <= 2 ? MvtxDefs::getStaveId(hskey) : -1;
      const int chip = layer <= 2 ? MvtxDefs::getChipId(hskey) : -1;
      const int ladderz = layer >= 3 ? InttDefs::getLadderZId(hskey) : -1;
      const int ladderphi = layer >= 3 ? InttDefs::getLadderPhiId(hskey) : -1;
      auto hits = hs->getHits();
      for (auto it = hits.first; it != hits.second; ++it)
      {
        const auto hk = it->first;
        const unsigned int row = layer <= 2 ? MvtxDefs::getRow(hk) : InttDefs::getRow(hk);
        const unsigned int col = layer <= 2 ? MvtxDefs::getCol(hk) : InttDefs::getCol(hk);
        const auto [lx, ly] = localCoordinates(layer, row, col, ladderz);
        SiHitPoint h;
        h.layer = layer;
        h.hitsetkey = hskey;
        h.hitkey = hk;
        h.row = row;
        h.col = col;
        h.stave = stave;
        h.chip = chip;
        h.ladderphi = ladderphi;
        h.ladderz = ladderz;
        h.adc = it->second ? it->second->getAdc() : 0;
        h.lx = lx;
        h.ly = ly;
        h.z = longitudinalZ(layer, chip, ladderz, ly);
        h.phi = intrinsicPhi(layer, stave, ladderphi, lx);
        m_eventData->hits.push_back(h);
      }
    }
  }
  return true;
}

void SiHitSeedReco::fillClusters(PHCompositeNode* topNode)
{
  auto* cont = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterNodeName);
  if (!cont)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": no " << m_clusterNodeName << " node, clusters=0" << std::endl;
    }
    return;
  }
  for (auto det : {TrkrDefs::mvtxId, TrkrDefs::inttId})
  {
    auto keys = cont->getHitSetKeys(det);
    for (auto hskey : keys)
    {
      const int layer = TrkrDefs::getLayer(hskey);
      auto range = cont->getClusters(hskey);
      for (auto it = range.first; it != range.second; ++it)
      {
        auto* c = it->second;
        if (!c)
        {
          continue;
        }
        SiClusterPoint p;
        p.key = it->first;
        p.layer = layer;
        p.lx = c->getLocalX();
        p.ly = c->getLocalY();
        p.adc = c->getAdc();
        if (layer <= 2)
        {
          p.stave = MvtxDefs::getStaveId(hskey);
          p.chip = MvtxDefs::getChipId(hskey);
        }
        else
        {
          p.ladderz = InttDefs::getLadderZId(hskey);
          p.ladderphi = InttDefs::getLadderPhiId(hskey);
        }
        p.z = longitudinalZ(layer, p.chip, p.ladderz, p.ly);
        p.phi = intrinsicPhi(layer, p.stave, p.ladderphi, p.lx);
        m_eventData->clusters.push_back(p);
      }
    }
  }
}

void SiHitSeedReco::buildBlobs()
{
  using Cell = std::pair<int, int>;
  std::array<std::map<Cell, std::vector<int>>, nlayers> layerCells;
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  for (size_t i = 0; i < m_eventData->hits.size(); ++i)
  {
    auto& h = m_eventData->hits[i];
    const int iz = (int) std::floor((h.z - m_zmin) / dz);
    const int ip = (int) std::floor(h.phi / dp) % m_nphi;
    if (iz < 0 || iz >= m_nz)
    {
      continue;
    }
    layerCells[h.layer][{iz, ip}].push_back(i);
  }

  for (int layer = 6; layer >= 0; --layer)
  {
    auto& remaining = layerCells[layer];
    while (!remaining.empty())
    {
      const Cell start = remaining.begin()->first;
      std::queue<Cell> q;
      q.push(start);
      std::vector<Cell> comp;
      std::vector<int> ids;
      while (!q.empty())
      {
        auto c = q.front();
        q.pop();
        auto it = remaining.find(c);
        if (it == remaining.end())
        {
          continue;
        }
        for (int id : it->second)
        {
          ids.push_back(id);
        }
        comp.push_back(c);
        remaining.erase(it);
        for (int a = -1; a <= 1; ++a)
        {
          for (int b = -1; b <= 1; ++b)
          {
            int ip = (c.second + b) % m_nphi;
            if (ip < 0)
            {
              ip += m_nphi;
            }
            const Cell n{c.first + a, ip};
            if (remaining.count(n))
            {
              q.push(n);
            }
          }
        }
      }
      std::sort(ids.begin(), ids.end());
      SiBlob blob;
      blob.id = m_eventData->blobs.size();
      blob.layer = layer;
      blob.hit_ids = ids;
      blob.cells = comp;
      std::vector<double> phis;
      double zsum = 0;
      for (int id : ids)
      {
        zsum += m_eventData->hits[id].z;
        phis.push_back(m_eventData->hits[id].phi);
        blob.adc_sum += m_eventData->hits[id].adc;
      }
      blob.z = zsum / ids.size();
      blob.phi = circularMean(phis);
      blob.zbin = (blob.z - m_zmin) / dz;
      blob.phibin = blob.phi / dp;
      m_eventData->blobs.push_back(std::move(blob));
    }
  }
}

// Hit -> blob map and a (layer, iz, iphi) cell index, so that every windowed search
// only visits the few cells inside its window instead of all hits of the layer.
void SiHitSeedReco::buildHitGrid()
{
  const auto& hits = m_eventData->hits;
  m_hitBlob.assign(hits.size(), -1);
  for (const auto& b : m_eventData->blobs)
  {
    for (int id : b.hit_ids)
    {
      m_hitBlob[id] = b.id;
    }
  }

  const size_t ncell = static_cast<size_t>(nlayers) * m_nz * m_nphi;
  m_cellStart.assign(ncell + 1, 0);
  std::vector<int> cellOf(hits.size(), -1);
  for (size_t i = 0; i < hits.size(); ++i)
  {
    if (m_hitBlob[i] < 0)
    {
      continue;  // outside the z binning: not used by the finder
    }
    const int iz = (int) std::floor(hits[i].zbin);
    const int ip = ((int) std::floor(hits[i].phibin) % m_nphi + m_nphi) % m_nphi;
    if (iz < 0 || iz >= m_nz)
    {
      continue;
    }
    cellOf[i] = (hits[i].layer * m_nz + iz) * m_nphi + ip;
    ++m_cellStart[cellOf[i] + 1];
  }
  for (size_t c = 0; c < ncell; ++c)
  {
    m_cellStart[c + 1] += m_cellStart[c];
  }
  m_cellHits.assign(m_cellStart.back(), -1);
  std::vector<int> fill(m_cellStart.begin(), m_cellStart.end() - 1);
  for (size_t i = 0; i < hits.size(); ++i)
  {
    if (cellOf[i] >= 0)
    {
      m_cellHits[fill[cellOf[i]]++] = static_cast<int>(i);
    }
  }
}

template <class F>
void SiHitSeedReco::forEachHitNear(int layer, double zbinCenter, double zHalf,
                                   double phibinCenter, double phiHalf, F&& f) const
{
  const int izlo = std::max(0, (int) std::floor(zbinCenter - zHalf));
  const int izhi = std::min(m_nz - 1, (int) std::floor(zbinCenter + zHalf));
  if (izlo > izhi)
  {
    return;
  }
  int iplo = (int) std::floor(phibinCenter - phiHalf);
  int iphi = (int) std::floor(phibinCenter + phiHalf);
  if (iphi - iplo + 1 >= m_nphi)
  {
    iplo = 0;
    iphi = m_nphi - 1;  // window covers the whole circle: visit each cell once
  }
  for (int iz = izlo; iz <= izhi; ++iz)
  {
    for (int ipr = iplo; ipr <= iphi; ++ipr)
    {
      const int ip = (ipr % m_nphi + m_nphi) % m_nphi;
      const int c = (layer * m_nz + iz) * m_nphi + ip;
      for (int k = m_cellStart[c]; k < m_cellStart[c + 1]; ++k)
      {
        f(m_cellHits[k]);
      }
    }
  }
}

void SiHitSeedReco::seedZWindow(int outerLayer, double zOuter, double& center, double& half) const
{
  const double dz = (m_zmax - m_zmin) / m_nz;  // cm per bin in this event
  if (m_zMode == ZWindowMode::Vertex && std::isfinite(m_eventData->vertex_z))
  {
    const double slope = m_radius[outerLayer - 1] / m_radius[outerLayer] - 1.0;
    center = (slope * (zOuter - m_eventData->vertex_z) + m_zVtxOffsetCm) / dz;
    half = m_zVtxHalfCm / dz;
    return;
  }
  if (m_zMode == ZWindowMode::Band || m_zMode == ZWindowMode::Vertex)  // Vertex without a vertex: band
  {
    const double zscale = std::max(std::abs(m_zmin), std::abs(m_zmax));
    center = m_zBandOffset - (zscale > 0 ? m_zBandShiftAtEdge * zOuter / zscale : 0.0);
    half = m_zBandHalf;
  }
  else if (m_zMode == ZWindowMode::Centimetre)
  {
    center = (m_zInterceptCm + m_zSlopeCm * zOuter) / dz;
    half = m_zHalfCm / dz;
  }
  else
  {
    center = m_zIntercept + m_zSlope * zOuter;
    half = m_zHalf;
  }
}

Eigen::Vector3d SiHitSeedReco::xyz(int layer, double z, double phi) const
{
  return {m_radius[layer] * std::cos(phi), m_radius[layer] * std::sin(phi), z};
}

SiHitSeedReco::Line3D SiHitSeedReco::fitLine(const std::vector<Eigen::Vector3d>& p) const
{
  Line3D l;
  l.c = Eigen::Vector3d::Zero();
  for (const auto& x : p)
  {
    l.c += x;
  }
  l.c /= p.size();
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto& x : p)
  {
    auto d = x - l.c;
    cov += d * d.transpose();
  }
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
  l.d = es.eigenvectors().col(2).normalized();
  return l;
}

bool SiHitSeedReco::propagate(const Line3D& l, int layer, const Eigen::Vector3d& near, double& z, double& phi) const
{
  const double r = m_radius[layer], a = l.d.x() * l.d.x() + l.d.y() * l.d.y();
  if (a < 1e-12)
  {
    return false;
  }
  const double b = 2 * (l.c.x() * l.d.x() + l.c.y() * l.d.y());
  const double q = l.c.x() * l.c.x() + l.c.y() * l.c.y() - r * r;
  const double disc = b * b - 4 * a * q;
  if (disc < 0)
  {
    return false;
  }
  const double s = std::sqrt(disc);
  const Eigen::Vector3d p1 = l.c + ((-b - s) / (2 * a)) * l.d;
  const Eigen::Vector3d p2 = l.c + ((-b + s) / (2 * a)) * l.d;
  const auto& p = ((p1 - near).squaredNorm() < (p2 - near).squaredNorm()) ? p1 : p2;
  z = p.z();
  phi = std::atan2(p.y(), p.x());
  if (phi < 0)
  {
    phi += twopi;
  }
  return true;
}

std::vector<SiHitSeedReco::Candidate> SiHitSeedReco::seedCandidates(const SiBlob& seed, int inner, const std::vector<char>& used) const
{
  std::unordered_map<int, Candidate> best;
  auto better = [](const Candidate& a, const Candidate& b)
  {
    return std::make_tuple(std::abs(a.res_phi), std::abs(a.res_z) / a.z_half, a.blob) <
           std::make_tuple(std::abs(b.res_phi), std::abs(b.res_z) / b.z_half, b.blob);
  };
  for (int oh : seed.hit_ids)
  {
    const auto& o = m_eventData->hits[oh];
    double center = 0, zhalf = 0;
    seedZWindow(seed.layer, o.z, center, zhalf);
    forEachHitNear(inner, o.zbin + center, zhalf, o.phibin, m_seedPhiHalf, [&](int ih)
    {
      const int b = m_hitBlob[ih];
      if (b < 0 || used[b])
      {
        return;
      }
      const auto& h = m_eventData->hits[ih];
      const double dphi = wrapDelta(h.phibin - o.phibin, m_nphi);
      const double dzres = h.zbin - (o.zbin + center);
      if (std::abs(dphi) > m_seedPhiHalf || std::abs(dzres) > zhalf)
      {
        return;
      }
      Candidate c;
      c.blob = b;
      c.hit = ih;
      c.layer = inner;
      c.outer_hit = oh;
      c.z = h.z;
      c.phi = h.phi;
      c.outer_z = o.z;
      c.outer_phi = o.phi;
      c.res_z = dzres;
      c.res_phi = dphi;
      c.z_center = center;
      c.z_half = zhalf;
      c.score = (dphi / m_seedPhiHalf) * (dphi / m_seedPhiHalf) + (dzres / zhalf) * (dzres / zhalf);
      auto jt = best.find(b);
      if (jt == best.end() || better(c, jt->second))
      {
        best[b] = c;
      }
    });
  }
  std::vector<Candidate> out;
  out.reserve(best.size());
  for (auto& [_, c] : best)
  {
    out.push_back(c);
  }
  std::sort(out.begin(), out.end(), better);  // closest in Uphi first
  return out;
}

bool SiHitSeedReco::buildChain(const SiBlob& seed, const Candidate& cand, const std::vector<char>& globallyUsed, SiHitChain& out) const
{
  // Blobs taken by this candidate chain (at most one per layer) - no copy of the global mask.
  std::array<int, nlayers> chainBlobs{};
  int nChainBlobs = 0;
  auto isUsed = [&](int b)
  {
    return globallyUsed[b] || std::find(chainBlobs.begin(), chainBlobs.begin() + nChainBlobs, b) != chainBlobs.begin() + nChainBlobs;
  };
  chainBlobs[nChainBlobs++] = seed.id;
  chainBlobs[nChainBlobs++] = cand.blob;

  struct A
  {
    int blob, hit;
    double z, phi;
  };
  std::map<int, A> acc;
  acc[seed.layer] = {seed.id, cand.outer_hit, cand.outer_z, cand.outer_phi};
  acc[cand.layer] = {cand.blob, cand.hit, cand.z, cand.phi};
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;

  SiChainStep s;
  s.kind = SiChainStep::Seed;
  s.from_layer = seed.layer;
  s.to_layer = cand.layer;
  s.blob_id = cand.blob;
  s.hit_id = cand.hit;
  s.ref_z = cand.outer_z;
  s.ref_phi = cand.outer_phi;
  s.pred_z = cand.outer_z + cand.z_center * dz;
  s.pred_phi = cand.outer_phi;
  s.z = cand.z;
  s.phi = cand.phi;
  s.delta_z_bins = (cand.z - cand.outer_z) / dz;
  s.delta_phi_bins = wrapDelta(cand.phi - cand.outer_phi, twopi) / dp;
  s.res_z_bins = cand.res_z;
  s.res_phi_bins = cand.res_phi;
  s.z_half_window_bins = cand.z_half;
  s.phi_half_window_bins = m_seedPhiHalf;
  s.score = cand.score;
  out.steps.push_back(s);

  auto propagateTo = [&](int kind, int layer, const std::vector<int>& fitlayers) -> bool
  {
    std::vector<Eigen::Vector3d> pts;
    for (int l : fitlayers)
    {
      pts.push_back(xyz(l, acc.at(l).z, acc.at(l).phi));
    }
    if (pts.size() < 2)
    {
      return false;
    }
    const auto line = fitLine(pts);
    int ref = acc.begin()->first;
    double rd = 1e9;
    for (auto& [l, a] : acc)
    {
      const double d = std::abs(m_radius[l] - m_radius[layer]);
      if (d < rd)
      {
        rd = d;
        ref = l;
      }
    }
    double pz, pp;
    if (!propagate(line, layer, xyz(ref, acc[ref].z, acc[ref].phi), pz, pp))
    {
      return false;
    }
    int bestBlob = -1, bestHit = -1;
    double bestScore = 1e99, bz = 0, bp = 0, rz = 0, rp = 0;
    forEachHitNear(layer, (pz - m_zmin) / dz, m_propZHalf, pp / dp, m_propPhiHalf, [&](int hid)
    {
      const int b = m_hitBlob[hid];
      if (b < 0 || isUsed(b))
      {
        return;
      }
      const auto& h = m_eventData->hits[hid];
      const double rrphi = wrapDelta(h.phi - pp, twopi) / dp;
      const double rrz = (h.z - pz) / dz;
      if (std::abs(rrphi) > m_propPhiHalf || std::abs(rrz) > m_propZHalf)
      {
        return;
      }
      const double sc = (rrphi / m_propPhiHalf) * (rrphi / m_propPhiHalf) + (rrz / m_propZHalf) * (rrz / m_propZHalf);
      if (sc < bestScore)
      {
        bestScore = sc;
        bestBlob = b;
        bestHit = hid;
        bz = h.z;
        bp = h.phi;
        rz = rrz;
        rp = rrphi;
      }
    });
    if (bestBlob < 0)
    {
      return false;
    }
    chainBlobs[nChainBlobs++] = bestBlob;
    acc[layer] = {bestBlob, bestHit, bz, bp};
    SiChainStep st;
    st.kind = kind;
    st.from_layer = ref;
    st.to_layer = layer;
    st.blob_id = bestBlob;
    st.hit_id = bestHit;
    st.ref_z = acc[ref].z;
    st.ref_phi = acc[ref].phi;
    st.pred_z = pz;
    st.pred_phi = pp;
    st.z = bz;
    st.phi = bp;
    st.delta_z_bins = (bz - acc[ref].z) / dz;
    st.delta_phi_bins = wrapDelta(bp - acc[ref].phi, twopi) / dp;
    st.res_z_bins = rz;
    st.res_phi_bins = rp;
    st.z_half_window_bins = m_propZHalf;
    st.phi_half_window_bins = m_propPhiHalf;
    st.score = bestScore;
    out.steps.push_back(st);
    return true;
  };

  for (int l = 2; l >= 0; --l)
  {
    if (!acc.count(l))
    {
      std::vector<int> fit;
      for (auto& [x, _] : acc)
      {
        if (x <= 2)
        {
          fit.push_back(x);
        }
      }
      propagateTo(SiChainStep::MvtxPropagation, l, fit);
    }
  }
  std::vector<int> mvtx;
  for (auto& [l, _] : acc)
  {
    if (l <= 2)
    {
      mvtx.push_back(l);
    }
  }
  for (int l = 3; l <= 6; ++l)
  {
    std::vector<int> fit = m_refitWithIntt ? std::vector<int>{} : mvtx;
    if (m_refitWithIntt)
    {
      for (auto& [x, _] : acc)
      {
        fit.push_back(x);
      }
    }
    propagateTo(SiChainStep::InttPropagation, l, fit);
  }

  out.n_mvtx = 0;
  out.n_intt = 0;
  for (auto it = acc.rbegin(); it != acc.rend(); ++it)
  {
    out.layers.push_back(it->first);
    out.blob_ids.push_back(it->second.blob);
    out.point_hit_ids.push_back(it->second.hit);
    out.z.push_back(it->second.z);
    out.phi.push_back(it->second.phi);
    if (it->first <= 2)
    {
      ++out.n_mvtx;
    }
    else
    {
      ++out.n_intt;
    }
    for (int h : m_eventData->blobs[it->second.blob].hit_ids)
    {
      out.hit_ids.push_back(h);
    }
  }
  out.score = 0;
  for (auto& st : out.steps)
  {
    out.score += st.score;
  }
  return (int) out.layers.size() >= m_minLayers && out.n_mvtx >= m_minMvtxLayers && out.n_intt >= m_minInttLayers;
}

void SiHitSeedReco::findChains()
{
  std::vector<char> used(m_eventData->blobs.size(), 0);
  for (auto pair : {std::pair<int, int>{2, 1}, {1, 0}})
  {
    std::vector<int> seeds;
    for (const auto& b : m_eventData->blobs)
    {
      if (b.layer == pair.first)
      {
        seeds.push_back(b.id);
      }
    }
    // Largest blob first; ties by blob id (deterministic, as in the notebook).
    std::sort(seeds.begin(), seeds.end(), [&](int a, int b)
    {
      const auto na = m_eventData->blobs[a].hit_ids.size(), nb = m_eventData->blobs[b].hit_ids.size();
      return na != nb ? na > nb : a < b;
    });
    for (int sid : seeds)
    {
      if (used[sid])
      {
        continue;
      }
      const auto& seed = m_eventData->blobs[sid];
      for (const auto& c : seedCandidates(seed, pair.second, used))
      {
        SiHitChain ch;
        if (!buildChain(seed, c, used, ch))
        {
          continue;
        }
        ch.id = m_eventData->chains.size();
        for (int b : ch.blob_ids)
        {
          used[b] = 1;
        }
        m_eventData->chains.push_back(std::move(ch));
        break;
      }
    }
  }
}

// event->vertex_z_linefit [cm]: straight-line fit of dz [bins] vs outer z [cm] over the seed links,
// vertex_z_linefit = z where the line crosses dz = 0.  Iterated: links further than 2.5 rms from the
// line are dropped and the line refitted (otherwise random combinations accepted inside the
// window pull the line towards dz = 0).  NOTE: it only sees links the seed window accepted,
// so with a fixed window it mostly reproduces the window's own line; use vertex_z instead.
// dz_vs_z_intercept_bins / dz_vs_z_slope_bins_per_cm are the final line parameters.
void SiHitSeedReco::fitVertexFromSeedLinks()
{
  std::vector<double> x, y;
  for (const auto& c : m_eventData->chains)
  {
    for (const auto& s : c.steps)
    {
      if (s.kind == SiChainStep::Seed)
      {
        x.push_back(s.ref_z);
        y.push_back(s.delta_z_bins);
      }
    }
  }
  std::vector<char> keep(x.size(), 1);
  double a = 0, b = 0;
  bool ok = false;
  int nkept = 0;
  for (int it = 0; it < 4; ++it)
  {
    double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t i = 0; i < x.size(); ++i)
    {
      if (keep[i])
      {
        n += 1;
        sx += x[i];
        sy += y[i];
        sxx += x[i] * x[i];
        sxy += x[i] * y[i];
      }
    }
    const double den = n * sxx - sx * sx;
    if (n < 3 || den <= 0)
    {
      break;
    }
    b = (n * sxy - sx * sy) / den;
    a = (sy - b * sx) / n;
    ok = true;
    nkept = static_cast<int>(n);
    double s2 = 0;
    for (size_t i = 0; i < x.size(); ++i)
    {
      if (keep[i])
      {
        const double r = y[i] - (a + b * x[i]);
        s2 += r * r;
      }
    }
    const double rms = std::sqrt(s2 / (n - 2));
    if (rms <= 0)
    {
      break;
    }
    for (size_t i = 0; i < x.size(); ++i)
    {
      keep[i] = std::abs(y[i] - (a + b * x[i])) <= 2.5 * rms;
    }
  }
  m_eventData->n_vertex_links = nkept;
  if (!ok)
  {
    return;
  }
  m_eventData->dz_vs_z_intercept_bins = a;
  m_eventData->dz_vs_z_slope_bins_per_cm = b;
  if (std::abs(b) > 1e-12)
  {
    m_eventData->vertex_z_linefit = -a / b;
  }
}

// event->vertex_z [cm]: window-independent tracklet vertex, computed before chain finding.
// Every MVTX blob pair (outer layer L, inner L-1; L = 2, 1) with |dUphi| <= m_tvPhiHalf bins and
// |dz| <= m_tvMaxDzCm gives the straight-line extrapolation to R = 0:
//   z0 = z_out - (z_in - z_out) * R_out / (R_in - R_out).
// True pairs pile up at the vertex, random pairs are spread out.  The densest window of width
// m_tvPeakCm is found in the z0 histogram and vertex_z = mean z0 inside it (twice refined).
void SiHitSeedReco::findTrackletVertex()
{
  const double dz = (m_zmax - m_zmin) / m_nz;
  const int nb = std::max(1, (int) std::lround((m_tvZmax - m_tvZmin) / m_tvBinCm));
  std::vector<int> hist(nb, 0);
  std::vector<double> z0s;
  std::vector<char> seen(m_eventData->blobs.size(), 0);
  std::vector<int> touched;
  for (const int outerLayer : {2, 1})
  {
    const int innerLayer = outerLayer - 1;
    const double ro = m_radius[outerLayer], ri = m_radius[innerLayer];
    for (const auto& o : m_eventData->blobs)
    {
      if (o.layer != outerLayer)
      {
        continue;
      }
      touched.clear();
      forEachHitNear(innerLayer, o.zbin, m_tvMaxDzCm / dz, o.phibin, m_tvPhiHalf, [&](int ih)
      {
        const int b = m_hitBlob[ih];
        if (b < 0 || seen[b])
        {
          return;
        }
        seen[b] = 1;
        touched.push_back(b);
        const auto& in = m_eventData->blobs[b];
        if (std::abs(wrapDelta(in.phibin - o.phibin, m_nphi)) > m_tvPhiHalf || std::abs(in.z - o.z) > m_tvMaxDzCm)
        {
          return;
        }
        const double z0 = o.z - (in.z - o.z) * ro / (ri - ro);
        const int ib = (int) std::floor((z0 - m_tvZmin) / m_tvBinCm);
        if (ib >= 0 && ib < nb)
        {
          ++hist[ib];
          z0s.push_back(z0);
        }
      });
      for (int b : touched)
      {
        seen[b] = 0;
      }
    }
  }
  m_eventData->vertex_tracklet_npairs = static_cast<int>(z0s.size());
  if (z0s.empty())
  {
    return;
  }
  // densest window of w bins
  const int w = std::max(1, (int) std::lround(m_tvPeakCm / m_tvBinCm));
  int best = 0, bestSum = -1, sum = 0;
  for (int i = 0; i < nb; ++i)
  {
    sum += hist[i];
    if (i >= w)
    {
      sum -= hist[i - w];
    }
    if (i >= w - 1 && sum > bestSum)
    {
      bestSum = sum;
      best = i - w + 1;
    }
  }
  double center = m_tvZmin + (best + 0.5 * w) * m_tvBinCm;
  int npeak = 0;
  for (int pass = 0; pass < 2; ++pass)
  {
    double s = 0;
    npeak = 0;
    for (double z0 : z0s)
    {
      if (std::abs(z0 - center) <= 0.5 * m_tvPeakCm)
      {
        s += z0;
        ++npeak;
      }
    }
    if (npeak == 0)
    {
      return;
    }
    center = s / npeak;
  }
  m_eventData->vertex_tracklet_npeak = npeak;
  // Expected flat background in a window of m_tvPeakCm from the pairs outside the peak.
  const double bkg = (double) (z0s.size() - npeak) * m_tvPeakCm / std::max(1e-9, (m_tvZmax - m_tvZmin) - m_tvPeakCm);
  if (npeak < m_tvMinPeak || npeak < m_tvMinSoB * bkg)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": tracklet vertex rejected (peak " << npeak << " pairs at z0 = " << center
                << " cm, background " << bkg << ") -> no vertex, band window used" << std::endl;
    }
    return;  // vertex_z stays NaN
  }
  m_eventData->vertex_z = center;
}
