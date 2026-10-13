#include "SiHitSeedReco.h"

#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>

#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
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
#include <set>
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

const char* statusName(int s)
{
  switch (s)
  {
  case SiTpcSeedInfo::Chained:
    return "chained";
  case SiTpcSeedInfo::TooFewClusters:
    return "too few TPC clusters";
  case SiTpcSeedInfo::FitFailed:
    return "TPC helix fit failed";
  case SiTpcSeedInfo::NoPrediction:
    return "no prediction on the seed layers";
  case SiTpcSeedInfo::NoAnchor:
    return "no hit in the anchor window";
  case SiTpcSeedInfo::NoChain:
    return "no chain passing the layer requirements";
  case SiTpcSeedInfo::Rejected:
    return "best chain fails the agreement cut";
  default:
    return "?";
  }
}
}  // namespace

SiHitSeedReco::SiHitSeedReco(const std::string& name)
  : SubsysReco(name)
{
  m_map.frame = &m_frame;
  m_map.alignment = &m_beamAlign;
}

int SiHitSeedReco::InitRun(PHCompositeNode* topNode)
{
  m_map.frame = &m_frame;
  m_map.alignment = &m_beamAlign;
  if (Verbosity() > 0)
  {
    const auto& c = m_frame.detectorCenterCm();
    std::cout << Name() << ": TPC-driven chaining, TPC from " << m_tpcTrackNodeName << " + " << m_tpcClusterNodeName
              << " (>= " << m_minTpcClusters << " clusters, z offset " << m_tpcZOffset << " cm)" << std::endl
              << Name() << ":   frames: Uphi rotation " << m_frame.uphiRotation() << " rad, detector centre ("
              << 10 * c[0] << ", " << 10 * c[1] << ", " << 10 * c[2] << ") mm, Si halves on TPC beam axis "
              << (m_beamAlign.enabled() ? "on" : "off") << ", Bz " << m_bz << " T" << std::endl
              << Name() << ":   anchor window |r dphi| < " << m_anchorRphi << " cm, |dz| < " << m_anchorZ
              << " cm (max " << m_maxAnchors << " anchors); seed pairs 2->1" << (m_useInnerPair ? ", 1->0" : "") << std::endl
              << Name() << ":   seed window: z centre tanl_TPC * ds + " << m_seedZOffsetCm << " cm, +- " << m_seedZHalfCm << " cm; phi centre "
              << (m_seedPhiFromTpc ? "TPC curvature" : "0") << ", +- " << m_seedPhiHalf << " bins" << std::endl
              << Name() << ":   propagation (hardware coordinates): phi +- " << m_propPhiHalf << " bins, z +- "
              << m_propZHalf << " bins; chains >= " << m_minLayers << " layers, >= "
              << m_minMvtxLayers << " MVTX, >= " << m_minInttLayers << " INTT" << std::endl
              << Name() << ":   selection: most layers, then chi2/ndf vs TPC (sigma r*phi " << m_sigRphi << ", z "
              << m_sigZ[0] << "/" << m_sigZ[1] << " cm, z offset " << (m_freeZOffset ? "free" : "fixed")
              << ", max chi2/ndf " << (m_maxTpcChi2 > 0 ? std::to_string(m_maxTpcChi2) : std::string("off")) << ")"
              << std::endl;
  }
  if (m_applyAlignment)
  {
    const int ret = loadAlignment();
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return ret;
    }
  }
  return createNodes(topNode);
}

int SiHitSeedReco::loadAlignment()
{
  std::unordered_map<TrkrDefs::hitsetkey, SiAlignmentTransform::Params> params;
  const std::string& source = m_alignmentFile;
  const int nread = source.empty() ? -1 : m_alignment.readFile(source, params, Verbosity());
  if (nread < 0)
  {
    std::cerr << Name() << ": cannot read silicon alignment parameters from '" << source
              << "' (set it with setSiAlignmentParamsFile)" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Ideal sensor centres/axes in the DETECTOR frame, from the same nominal geometry that
  // fillHits uses: centre = local (0, 0) of the sensor, local x along +Uphi (tilted by the
  // MVTX stave tilt), local y along +z.  The axes matter only for the new rotation order.
  m_alignment.clear();
  const double rot = mvtx_phi0[2];
  for (const auto& [key, p] : params)
  {
    const int layer = TrkrDefs::getLayer(key);
    if (layer < 0 || layer >= nlayers)
    {
      continue;
    }
    const bool mvtx = layer <= 2;
    const int stave = mvtx ? MvtxDefs::getStaveId(key) : -1;
    const int chip = mvtx ? MvtxDefs::getChipId(key) : -1;
    const int ladderphi = mvtx ? -1 : InttDefs::getLadderPhiId(key);
    const int ladderz = mvtx ? -1 : InttDefs::getLadderZId(key);
    if ((mvtx && (stave >= nstaves[layer] || chip > 8)) ||
        (!mvtx && (ladderphi >= intt_nladder[layer - 3] || ladderz > 3)))
    {
      continue;
    }

    const double phi = intrinsicPhi(layer, stave, ladderphi, 0.0) + rot;
    const double z = longitudinalZ(layer, chip, ladderz, 0.0);
    const Eigen::Vector3d centre(m_radius[layer] * std::cos(phi), m_radius[layer] * std::sin(phi), z);

    const Eigen::Vector3d rhat(std::cos(phi), std::sin(phi), 0.0);
    const Eigen::Vector3d phihat(-std::sin(phi), std::cos(phi), 0.0);
    const double tilt = mvtx ? mvtx_tilt[layer] : 0.0;
    const Eigen::Vector3d ex = std::cos(tilt) * phihat + std::sin(tilt) * rhat;
    const Eigen::Vector3d ey = Eigen::Vector3d::UnitZ();
    Eigen::Matrix3d ideal;
    ideal.col(0) = ex;
    ideal.col(1) = ey;
    ideal.col(2) = ex.cross(ey);

    m_alignment.addSensor(SiAlignmentTransform::sensorId(key), p, centre, ideal);
  }

  std::cout << Name() << ": silicon alignment ON, " << m_alignment.size() << " sensors from " << source
            << (m_alignment.useNewSiliconRotationOrder() ? " (new" : " (old") << " rotation order)" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void SiHitSeedReco::alignPoint(TrkrDefs::hitsetkey hskey, int layer, double& z, double& phi) const
{
  if (!m_applyAlignment)
  {
    return;
  }
  const double rot = mvtx_phi0[2];
  const double p = phi + rot;
  const SiAlignmentTransform::Point det{m_radius[layer] * std::cos(p), m_radius[layer] * std::sin(p), z};
  const auto a = m_alignment.apply(SiAlignmentTransform::sensorId(hskey), det);
  z = a.z;
  phi = std::fmod(std::atan2(a.y, a.x) - rot, twopi);
  if (phi < 0)
  {
    phi += twopi;
  }
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

  const auto tpcTracks = collectTpcTracks(topNode);
  for (const auto& t : tpcTracks)
  {
    processTpcTrack(t);
  }

  if (Verbosity() > 0)
  {
    std::array<unsigned int, 7> st{};
    for (const auto& s : m_eventData->tpc_seeds)
    {
      ++st[std::clamp(s.status, 0, 6)];
    }
    unsigned int nSel = 0;
    for (const auto& c : m_eventData->chains)
    {
      nSel += c.selected;
    }
    std::cout << Name() << ": event " << m_eventData->event
              << " hits=" << m_eventData->hits.size()
              << " blobs=" << m_eventData->blobs.size()
              << " tpc_tracks=" << m_eventData->tpc_seeds.size()
              << " chains=" << m_eventData->chains.size()
              << " selected=" << nSel
              << " | chained " << st[SiTpcSeedInfo::Chained]
              << ", <clusters " << st[SiTpcSeedInfo::TooFewClusters]
              << ", fit failed " << st[SiTpcSeedInfo::FitFailed]
              << ", no prediction " << st[SiTpcSeedInfo::NoPrediction]
              << ", no anchor " << st[SiTpcSeedInfo::NoAnchor]
              << ", no chain " << st[SiTpcSeedInfo::NoChain]
              << ", rejected " << st[SiTpcSeedInfo::Rejected] << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiHitSeedReco::End(PHCompositeNode*)
{
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": " << m_event << " events, " << m_nChains << " chains stored; TPC tracks:";
    for (int s = 0; s < 7; ++s)
    {
      std::cout << " " << statusName(s) << " " << m_statusCount[s] << (s < 6 ? "," : "");
    }
    std::cout << std::endl;
    if (m_hwDdphiProp0 || m_hwDdphiCheck[0] || m_hwDdphiCheck[1])
    {
      std::cout << Name() << ": hw delta-delta-phi (" << m_hwDdphiNSigma << " sigma): layer-0 propagation windows "
                << m_hwDdphiProp0Used;
      for (int l = 0; l < 2; ++l)
      {
        if (m_hwDdphiCheck[l])
        {
          std::cout << ", layer " << l << " checked " << m_hwDdphiChecked[l]
                    << (m_hwDdphiReject ? " rejected " : " outside window ") << m_hwDdphiRejected[l];
        }
      }
      std::cout << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

bool SiHitSeedReco::hwDdphiTerms(int layer, const std::map<int, double>& uphiByLayer, double& dphi, double& dphiPred) const
{
  const int p = layer + 1;
  const auto iL = uphiByLayer.find(layer);
  const auto iP = uphiByLayer.find(p);
  if (iL == uphiByLayer.end() || iP == uphiByLayer.end())
  {
    return false;
  }
  const auto iPP = uphiByLayer.upper_bound(p);  // next outer layer of the chain
  if (iPP == uphiByLayer.end())
  {
    return false;
  }
  const double drPrev = m_radius[p] - m_radius[iPP->first];
  if (std::abs(drPrev) < 1e-3)
  {
    return false;
  }
  dphi = wrapDelta(iL->second - iP->second, twopi);
  dphiPred = wrapDelta(iP->second - iPP->second, twopi) * (m_radius[layer] - m_radius[p]) / drPrev;
  return true;
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
        alignPoint(hskey, layer, h.z, h.phi);
        m_eventData->hits.push_back(h);
      }
    }
  }
  return true;
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

// ============================================================================================
// TPC tracks
// ============================================================================================
std::vector<SiHitSeedReco::TpcTrackInput> SiHitSeedReco::collectTpcTracks(PHCompositeNode* topNode) const
{
  std::vector<TpcTrackInput> out;
  auto* clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  auto* tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName);
  if (!clusters || !tracks)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_tpcClusterNodeName << " or " << m_tpcTrackNodeName
                << " - no TPC-driven chaining in this event" << std::endl;
    }
    return out;
  }
  // TPC = reference frame: clusters as measured (SiTpcBeamAlignment::tpcToTpcFrame is the identity)
  std::map<unsigned int, std::vector<SiTpcHelixFit::Point>> byTrack;
  for (unsigned int i = 0; i < clusters->size(); ++i)
  {
    const Tpc_PolyCluster* c = clusters->get_cluster(i);
    if (!c || !c->isValid())
    {
      continue;
    }
    const auto p = m_beamAlign.tpcToTpcFrame(c->get_centroid_x(), c->get_centroid_y(), c->get_centroid_z());
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
    {
      continue;
    }
    byTrack[c->get_source_assembled_track_id()].push_back({p.x, p.y, p.z + m_tpcZOffset, 1.0, 1.0, 1});
  }
  SiTpcHelixFit::Config cfg;
  cfg.beamX = m_beamAlign.referenceBeamX();
  cfg.beamY = m_beamAlign.referenceBeamY();
  cfg.bz = m_bz;
  for (unsigned int i = 0; i < tracks->size(); ++i)
  {
    const Tpc_PolyTrack* t = tracks->get_track(i);
    if (!t)
    {
      continue;
    }
    TpcTrackInput in;
    in.index = static_cast<int>(i);
    in.id = t->get_track_id();
    in.assembledId = t->get_source_assembled_track_id();
    auto it = byTrack.find(in.assembledId);
    if (it != byTrack.end())
    {
      auto& pts = it->second;
      in.nClusters = static_cast<unsigned int>(pts.size());
      if (pts.size() >= m_minTpcClusters)
      {
        std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b)
                  { return a.x * a.x + a.y * a.y < b.x * b.x + b.y * b.y; });
        in.fit = SiTpcHelixFit::fit(pts, cfg);
      }
    }
    out.push_back(std::move(in));
  }
  return out;
}

// Blobs of `layer` with at least one hit in the anchor window around the TPC prediction
// (both clamshell translations near a boundary), closest first.
std::vector<int> SiHitSeedReco::anchorBlobs(const SiTpcReference::TpcReference& ref, int layer) const
{
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  std::unordered_map<int, double> metric;
  for (const auto* p : {&ref.layers[layer], &ref.alt[layer]})
  {
    if (!p->valid)
    {
      continue;
    }
    const double R = m_radius[layer];
    forEachHitNear(layer, (p->z - m_zmin) / dz, m_anchorZ / dz, p->u / dp, (m_anchorRphi / R) / dp, [&](int ih)
    {
      const int b = m_hitBlob[ih];
      if (b < 0)
      {
        return;
      }
      const auto& h = m_eventData->hits[ih];
      const double rdphi = R * SiTpcReference::wrapPi(h.phi - p->u);
      const double dzc = h.z - p->z;
      if (std::abs(rdphi) > m_anchorRphi || std::abs(dzc) > m_anchorZ)
      {
        return;
      }
      const double m = (rdphi / m_anchorRphi) * (rdphi / m_anchorRphi) + (dzc / m_anchorZ) * (dzc / m_anchorZ);
      auto it = metric.find(b);
      if (it == metric.end() || m < it->second)
      {
        metric[b] = m;
      }
    });
  }
  std::vector<std::pair<double, int>> order;
  order.reserve(metric.size());
  for (const auto& [b, m] : metric)
  {
    order.emplace_back(m, b);
  }
  std::sort(order.begin(), order.end());
  std::vector<int> out;
  for (const auto& [m, b] : order)
  {
    if (out.size() >= m_maxAnchors)
    {
      break;
    }
    out.push_back(b);
  }
  return out;
}

// Seed window outer hit (layer `outer`, z = zOuter) -> layer `inner`, in bins of this event:
//   z  : centre zCenterBins relative to the outer hit, half width zHalfBins
//   phi: centre phiCenterBins relative to the outer hit (TPC curvature)
void SiHitSeedReco::seedWindow(const SiTpcReference::TpcReference& ref, int outer, int inner, double /*zOuter*/,
                               double& zCenterBins, double& zHalfBins, double& phiCenterBins) const
{
  const double dz = (m_zmax - m_zmin) / m_nz, dp = twopi / m_nphi;
  const auto& pOut = ref.layers[outer];
  const auto& pIn = ref.layers[inner];
  // TPC slope only: dz = tanl * ds, ds = transverse path between the two layers along the TPC
  // circle (< 0 going inward).  The TPC pca z is NOT used, so a TPC z offset / t0 does not bias it.
  const double ds = (pOut.valid && pIn.valid) ? (pIn.s - pOut.s) : (m_radius[inner] - m_radius[outer]) / std::max(1e-3, pOut.cosA);
  const double dzcm = ref.tanl * ds;
  zCenterBins = (dzcm + m_seedZOffsetCm) / dz;
  zHalfBins = m_seedZHalfCm / dz;
  phiCenterBins = (m_seedPhiFromTpc && pOut.valid && pIn.valid) ? SiTpcReference::wrapPi(pIn.u - pOut.u) / dp : 0.0;
}

std::vector<SiHitSeedReco::Candidate> SiHitSeedReco::seedCandidates(const SiBlob& seed, int inner,
                                                                   const SiTpcReference::TpcReference& ref) const
{
  std::unordered_map<int, Candidate> best;
  auto better = [](const Candidate& a, const Candidate& b)
  {
    return std::make_tuple(a.score, a.blob) < std::make_tuple(b.score, b.blob);
  };
  for (int oh : seed.hit_ids)
  {
    const auto& o = m_eventData->hits[oh];
    double zc = 0, zh = 0, pc = 0;
    seedWindow(ref, seed.layer, inner, o.z, zc, zh, pc);
    forEachHitNear(inner, o.zbin + zc, zh, o.phibin + pc, m_seedPhiHalf, [&](int ih)
    {
      const int b = m_hitBlob[ih];
      if (b < 0)
      {
        return;
      }
      const auto& h = m_eventData->hits[ih];
      const double dphi = wrapDelta(h.phibin - (o.phibin + pc), m_nphi);
      const double dzres = h.zbin - (o.zbin + zc);
      if (std::abs(dphi) > m_seedPhiHalf || std::abs(dzres) > zh)
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
      c.z_center = zc;
      c.z_half = zh;
      c.phi_center = pc;
      c.score = (dphi / m_seedPhiHalf) * (dphi / m_seedPhiHalf) + (dzres / zh) * (dzres / zh);
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
  std::sort(out.begin(), out.end(), better);
  return out;
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

// Seed pair -> chain, as before: straight line through the MVTX points, propagated to the missing
// MVTX layer and to the INTT layers.  Everything here is in hardware coordinates: a hit is the
// point (R_layer cos Uphi, R_layer sin Uphi, Vz) of the nominal cylinder in the detector frame
// (rotated by the common Uphi rotation, which a straight line does not care about), the same
// frame for all hits, so no conversion is needed.  The chain goes to physics space only
// afterwards (toPhysics).  No global hit mask.
bool SiHitSeedReco::buildChain(const SiBlob& seed, const Candidate& cand, SiHitChain& out) const
{
  // blobs taken by this chain (at most one per layer)
  std::array<int, nlayers> chainBlobs{};
  int nChainBlobs = 0;
  auto isUsed = [&](int b)
  {
    return std::find(chainBlobs.begin(), chainBlobs.begin() + nChainBlobs, b) != chainBlobs.begin() + nChainBlobs;
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
  s.pred_phi = cand.outer_phi + cand.phi_center * dp;
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
    double phiHalf = m_propPhiHalf;
    if (layer == 0 && m_hwDdphiProp0 && acc.count(1) && acc.count(2))
    {
      // phi window from the step-to-step correlation instead of the straight line (z unchanged)
      // layer-0 Uphi is not known yet: only dUphi_pred(0) is used, the dUphi(0) entry is a dummy
      const std::map<int, double> uphi{{0, acc.at(1).phi}, {1, acc.at(1).phi}, {2, acc.at(2).phi}};
      double dummy = 0, pred = 0;
      if (hwDdphiTerms(0, uphi, dummy, pred))
      {
        pp = acc.at(1).phi + m_hwDdphiOffset[0] + m_hwDdphiSlope[0] * pred;
        phiHalf = m_hwDdphiNSigma * m_hwDdphiSigma[0] / dp;
        ++m_hwDdphiProp0Used;
      }
    }
    int bestBlob = -1, bestHit = -1;
    double bestScore = 1e99, bz = 0, bp = 0, rz = 0, rp = 0;
    forEachHitNear(layer, (pz - m_zmin) / dz, m_propZHalf, pp / dp, phiHalf, [&](int hid)
    {
      const int b = m_hitBlob[hid];
      if (b < 0 || isUsed(b))
      {
        return;
      }
      const auto& h = m_eventData->hits[hid];
      const double rrphi = wrapDelta(h.phi - pp, twopi) / dp;
      const double rrz = (h.z - pz) / dz;
      if (std::abs(rrphi) > phiHalf || std::abs(rrz) > m_propZHalf)
      {
        return;
      }
      const double sc = (rrphi / phiHalf) * (rrphi / phiHalf) + (rrz / m_propZHalf) * (rrz / m_propZHalf);
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
    st.phi_half_window_bins = phiHalf;
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

  // step-to-step delta-delta-phi check on the complete chain (layers 0, 1)
  double ddphiScore = 0;
  if (m_hwDdphiCheck[0] || m_hwDdphiCheck[1])
  {
    std::map<int, double> uphi;
    for (const auto& [l, a] : acc)
    {
      uphi[l] = a.phi;
    }
    for (int l = 0; l < 2; ++l)
    {
      double d = 0, pred = 0;
      if (!m_hwDdphiCheck[l] || !hwDdphiTerms(l, uphi, d, pred) || m_hwDdphiSigma[l] <= 0)
      {
        continue;
      }
      const double pull = (d - m_hwDdphiOffset[l] - m_hwDdphiSlope[l] * pred) / m_hwDdphiSigma[l];
      out.hw_ddphi_pull[l] = pull;
      ++m_hwDdphiChecked[l];
      if (std::abs(pull) > m_hwDdphiNSigma)
      {
        ++m_hwDdphiRejected[l];
        if (m_hwDdphiReject)
        {
          return false;
        }
      }
      ddphiScore += pull * pull;
    }
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
  out.score = ddphiScore;
  for (auto& st : out.steps)
  {
    out.score += st.score;
  }
  return (int) out.layers.size() >= m_minLayers && out.n_mvtx >= m_minMvtxLayers && out.n_intt >= m_minInttLayers;
}

// Chain -> physics space (TPC frame) exactly as SiTrajectoryFitter::layerPoints, then residuals
// to the TPC helix at the radius of every point and the agreement chi2/ndf.
void SiHitSeedReco::toPhysics(const SiTpcHelixFit::Result& tpcFit, SiHitChain& chain) const
{
  chain.phys_x.clear();
  chain.phys_y.clear();
  chain.phys_z.clear();
  chain.res_rdphi.clear();
  chain.res_dz.clear();
  for (size_t i = 0; i < chain.layers.size(); ++i)
  {
    const int layer = chain.layers[i];
    const auto& blob = m_eventData->blobs[chain.blob_ids[i]];
    double sx = 0, sy = 0, sz = 0;
    for (int hid : blob.hit_ids)
    {
      const auto& h = m_eventData->hits[hid];
      const auto det = m_frame.toDetector(layer, h.phi, h.z);
      sx += det.x;
      sy += det.y;
      sz += det.z;
    }
    const double n = std::max<size_t>(1, blob.hit_ids.size());
    const double r = SiDetectorFrame::layerRadius(layer);
    const double phi = std::atan2(sy, sx);
    const SiDetectorFrame::Point det{r * std::cos(phi), r * std::sin(phi), sz / n};
    const auto glob = m_frame.toGlobal(det);
    const auto tpc = m_beamAlign.siHalfToTpcAxis(m_beamAlign.siHalf(phi - m_frame.uphiRotation()), glob);
    chain.phys_x.push_back(tpc.x);
    chain.phys_y.push_back(tpc.y);
    chain.phys_z.push_back(tpc.z);

    double rdphi = std::numeric_limits<double>::quiet_NaN();
    double dzr = std::numeric_limits<double>::quiet_NaN();
    const double rad = std::hypot(tpc.x, tpc.y);
    SiTpcHelixFit::Pos hp;
    double s = 0;
    if (SiTpcHelixFit::atRadius(tpcFit, rad, tpc.z, hp, s))
    {
      rdphi = rad * SiTpcReference::wrapPi(std::atan2(tpc.y, tpc.x) - std::atan2(hp.y, hp.x));
      dzr = tpc.z - hp.z;
    }
    chain.res_rdphi.push_back(rdphi);
    chain.res_dz.push_back(dzr);
  }

  double sw = 0, swz = 0;
  int nphi = 0, nz = 0;
  for (size_t i = 0; i < chain.layers.size(); ++i)
  {
    if (std::isfinite(chain.res_dz[i]))
    {
      const double sig = m_sigZ[chain.layers[i] <= 2 ? 0 : 1];
      sw += 1.0 / (sig * sig);
      swz += chain.res_dz[i] / (sig * sig);
      ++nz;
    }
  }
  chain.tpc_dz_offset = sw > 0 ? swz / sw : 0.0;
  const double off = m_freeZOffset ? chain.tpc_dz_offset : 0.0;
  double chi2 = 0;
  for (size_t i = 0; i < chain.layers.size(); ++i)
  {
    if (std::isfinite(chain.res_rdphi[i]))
    {
      chi2 += (chain.res_rdphi[i] / m_sigRphi) * (chain.res_rdphi[i] / m_sigRphi);
      ++nphi;
    }
    if (std::isfinite(chain.res_dz[i]))
    {
      const double sig = m_sigZ[chain.layers[i] <= 2 ? 0 : 1];
      const double d = (chain.res_dz[i] - off) / sig;
      chi2 += d * d;
    }
  }
  chain.tpc_ndf = nphi + nz - ((m_freeZOffset && nz > 0) ? 1 : 0);
  chain.tpc_chi2 = chain.tpc_ndf > 0 ? chi2 / chain.tpc_ndf : std::numeric_limits<double>::quiet_NaN();
}

void SiHitSeedReco::processTpcTrack(const TpcTrackInput& in)
{
  SiTpcSeedInfo info;
  info.tpc_track_index = in.index;
  info.tpc_track_id = in.id;
  info.tpc_assembled_track_id = in.assembledId;
  info.n_clusters = in.nClusters;
  auto finish = [&](int status)
  {
    info.status = status;
    ++m_statusCount[std::clamp(status, 0, 6)];
    m_eventData->tpc_seeds.push_back(info);
    if (Verbosity() > 1)
    {
      std::cout << Name() << ":   TPC track " << in.index << " (assembled " << in.assembledId << ", "
                << in.nClusters << " clusters, pt " << info.pt << ", z0 " << info.z0 << "): " << statusName(status)
                << ", anchors " << info.n_anchors << ", chains " << info.n_chains << ", selected chain "
                << info.selected_chain << std::endl;
    }
  };

  if (in.nClusters < m_minTpcClusters)
  {
    finish(SiTpcSeedInfo::TooFewClusters);
    return;
  }
  if (in.fit.fitted())
  {
    info.pt = in.fit.pt;
    info.phi = in.fit.phi;
    info.tanl = in.fit.tanl;
    info.z0 = in.fit.z0;
    info.charge = in.fit.charge;
  }
  SiTpcReference::TpcReference ref;
  if (!SiTpcReference::buildTpcReference(in.fit, m_map, ref, m_boundaryMargin))
  {
    finish(SiTpcSeedInfo::FitFailed);
    return;
  }

  // ---- chains of this TPC track (hits may be shared, identical chains once)
  std::vector<SiHitChain> chains;
  std::set<std::vector<int>> signatures;
  bool anyPrediction = false;
  std::vector<std::pair<int, int>> pairs{{2, 1}};
  if (m_useInnerPair)
  {
    pairs.emplace_back(1, 0);
  }
  for (const auto& [outer, inner] : pairs)
  {
    if (chains.size() >= m_maxChainsPerTrack)
    {
      break;
    }
    if (!ref.layers[outer].valid)
    {
      continue;
    }
    anyPrediction = true;
    const auto anchors = anchorBlobs(ref, outer);
    info.n_anchors += static_cast<int>(anchors.size());
    for (int a : anchors)
    {
      if (chains.size() >= m_maxChainsPerTrack)
      {
        break;
      }
      const auto& seed = m_eventData->blobs[a];
      const auto cands = seedCandidates(seed, inner, ref);
      const size_t ncand = std::min<size_t>(cands.size(), m_maxSeedCand);
      for (size_t k = 0; k < ncand && chains.size() < m_maxChainsPerTrack; ++k)
      {
        SiHitChain ch;
        if (!buildChain(seed, cands[k], ch))
        {
          continue;
        }
        std::vector<int> sig = ch.blob_ids;
        std::sort(sig.begin(), sig.end());
        if (!signatures.insert(sig).second)
        {
          continue;
        }
        ch.seed_outer_layer = outer;
        ch.anchor_blob = a;
        chains.push_back(std::move(ch));
      }
    }
  }
  info.n_chains = static_cast<int>(chains.size());
  if (!anyPrediction)
  {
    finish(SiTpcSeedInfo::NoPrediction);
    return;
  }
  if (info.n_anchors == 0)
  {
    finish(SiTpcSeedInfo::NoAnchor);
    return;
  }
  if (chains.empty())
  {
    finish(SiTpcSeedInfo::NoChain);
    return;
  }

  // ---- physics space, agreement with the TPC helix, ranking
  for (auto& ch : chains)
  {
    toPhysics(in.fit, ch);
  }
  std::sort(chains.begin(), chains.end(), [](const SiHitChain& a, const SiHitChain& b)
  {
    if (a.n_layers() != b.n_layers())
    {
      return a.n_layers() > b.n_layers();
    }
    const bool fa = std::isfinite(a.tpc_chi2), fb = std::isfinite(b.tpc_chi2);
    if (fa != fb)
    {
      return fa;
    }
    if (fa && a.tpc_chi2 != b.tpc_chi2)
    {
      return a.tpc_chi2 < b.tpc_chi2;
    }
    return a.score < b.score;
  });
  for (size_t r = 0; r < chains.size(); ++r)
  {
    chains[r].rank = static_cast<int>(r);
  }
  const bool pass = m_maxTpcChi2 <= 0 || (std::isfinite(chains[0].tpc_chi2) && chains[0].tpc_chi2 <= m_maxTpcChi2);
  chains[0].selected = pass;

  // ---- store
  for (auto& ch : chains)
  {
    if (m_storeSelectedOnly && !ch.selected)
    {
      continue;
    }
    ch.id = static_cast<int>(m_eventData->chains.size());
    ch.tpc_track_index = in.index;
    ch.tpc_track_id = in.id;
    ch.tpc_assembled_track_id = in.assembledId;
    if (ch.selected)
    {
      info.selected_chain = ch.id;
    }
    if (Verbosity() > 2)
    {
      std::cout << Name() << ":     chain " << ch.id << " rank " << ch.rank << (ch.selected ? " SELECTED" : "")
                << " layers " << ch.n_layers() << " (MVTX " << ch.n_mvtx << ", INTT " << ch.n_intt << ") chi2/ndf "
                << ch.tpc_chi2 << " (ndf " << ch.tpc_ndf << ") dz offset " << ch.tpc_dz_offset << " cm, seed pair "
                << ch.seed_outer_layer << "->" << ch.seed_outer_layer - 1 << std::endl;
    }
    m_eventData->chains.push_back(std::move(ch));
    ++m_nChains;
  }
  finish(pass ? SiTpcSeedInfo::Chained : SiTpcSeedInfo::Rejected);
}