#include "Full_PolyTrackMatcher.h"

#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrackv1.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyClusterContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>

namespace
{
  float finite_float(double value)
  {
    return std::isfinite(value) ? static_cast<float>(value) : std::numeric_limits<float>::quiet_NaN();
  }

  double square(double value)
  {
    return value * value;
  }
}

Full_PolyTrackMatcher::Full_PolyTrackMatcher(const std::string& name)
  : SubsysReco(name)
{
}

Full_PolyTrackMatcher::~Full_PolyTrackMatcher()
{
  if (m_qaFile)
  {
    m_qaFile->Close();
    delete m_qaFile;
    m_qaFile = nullptr;
  }
}

int Full_PolyTrackMatcher::Init(PHCompositeNode*)
{
  createQaObjects();
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackMatcher::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackMatcher::getNodes(PHCompositeNode* topNode)
{
  m_tpcTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_tpcTrackNodeName);
  if (!m_tpcTracks)
  {
    std::cerr << Name() << "::getNodes - missing " << m_tpcTrackNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tpcClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  if (!m_tpcClusters)
  {
    std::cerr << Name() << "::getNodes - missing " << m_tpcClusterNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  if (!m_crossingDecisions && Verbosity() > 0)
  {
    std::cout << Name() << "::getNodes - optional " << m_crossingDecisionNodeName << " not found" << std::endl;
  }

  m_trkrClusters = findNode::getClass<TrkrClusterContainer>(topNode, m_trkrClusterNodeName);
  if (!m_trkrClusters)
  {
    std::cerr << Name() << "::getNodes - missing " << m_trkrClusterNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_actsGeometry)
  {
    std::cerr << Name() << "::getNodes - missing ActsGeometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackMatcher::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_fullTracks = findNode::getClass<Full_PolyTrackContainer>(topNode, m_outputNodeName);
  if (!m_fullTracks)
  {
    m_fullTracks = new Full_PolyTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_fullTracks, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void Full_PolyTrackMatcher::createQaObjects()
{
  if (!m_writeQA || m_qaFile)
  {
    return;
  }

  m_qaFile = new TFile(m_qaFileName.c_str(), "RECREATE");
  m_candidateTree = new TTree("full_polytrack_candidates", "Full_PolyTrack silicon association candidates");
  m_candidateTree->Branch("event", &m_qaRow.event, "event/i");
  m_candidateTree->Branch("tpc_track_id", &m_qaRow.tpc_track_id, "tpc_track_id/i");
  m_candidateTree->Branch("source_assembled_track_id", &m_qaRow.source_assembled_track_id, "source_assembled_track_id/i");
  m_candidateTree->Branch("crossing", &m_qaRow.crossing, "crossing/I");
  m_candidateTree->Branch("layer", &m_qaRow.layer, "layer/i");
  m_candidateTree->Branch("ladder", &m_qaRow.ladder, "ladder/i");
  m_candidateTree->Branch("sensor", &m_qaRow.sensor, "sensor/i");
  m_candidateTree->Branch("cluster_key", &m_qaRow.cluster_key, "cluster_key/l");
  m_candidateTree->Branch("x", &m_qaRow.x, "x/F");
  m_candidateTree->Branch("y", &m_qaRow.y, "y/F");
  m_candidateTree->Branch("z", &m_qaRow.z, "z/F");
  m_candidateTree->Branch("r", &m_qaRow.r, "r/F");
  m_candidateTree->Branch("phi", &m_qaRow.phi, "phi/F");
  m_candidateTree->Branch("pred_x", &m_qaRow.pred_x, "pred_x/F");
  m_candidateTree->Branch("pred_y", &m_qaRow.pred_y, "pred_y/F");
  m_candidateTree->Branch("pred_z", &m_qaRow.pred_z, "pred_z/F");
  m_candidateTree->Branch("pred_phi", &m_qaRow.pred_phi, "pred_phi/F");
  m_candidateTree->Branch("rdphi", &m_qaRow.rdphi, "rdphi/F");
  m_candidateTree->Branch("dz", &m_qaRow.dz, "dz/F");
  m_candidateTree->Branch("chi2", &m_qaRow.chi2, "chi2/F");
  m_candidateTree->Branch("chain_chi2", &m_qaRow.chain_chi2, "chain_chi2/F");
  m_candidateTree->Branch("chain_ndf", &m_qaRow.chain_ndf, "chain_ndf/F");
  m_candidateTree->Branch("chain_score", &m_qaRow.chain_score, "chain_score/F");
  m_candidateTree->Branch("chain_silicon_clusters", &m_qaRow.chain_silicon_clusters, "chain_silicon_clusters/i");
  m_candidateTree->Branch("chain_missing_layers", &m_qaRow.chain_missing_layers, "chain_missing_layers/i");
  m_candidateTree->Branch("chain_missing_mask", &m_qaRow.chain_missing_mask, "chain_missing_mask/i");
  m_candidateTree->Branch("accepted", &m_qaRow.accepted, "accepted/I");
  m_candidateTree->Branch("selected", &m_qaRow.selected, "selected/I");
}

double Full_PolyTrackMatcher::wrapPhi(double phi) const
{
  while (phi > M_PI) { phi -= 2.0 * M_PI; }
  while (phi <= -M_PI) { phi += 2.0 * M_PI; }
  return phi;
}

double Full_PolyTrackMatcher::unwrapPhiNear(double phi, double reference) const
{
  while (phi - reference > M_PI) { phi -= 2.0 * M_PI; }
  while (phi - reference <= -M_PI) { phi += 2.0 * M_PI; }
  return phi;
}

unsigned int Full_PolyTrackMatcher::layerBit(unsigned int layer) const
{
  return layer < 32U ? (1U << layer) : 0U;
}

bool Full_PolyTrackMatcher::getGlobalClusterPosition(TrkrDefs::cluskey key, TrkrCluster* cluster, SpacePoint& point) const
{
  if (!cluster || !m_actsGeometry)
  {
    return false;
  }

  const auto global = m_actsGeometry->getGlobalPosition(key, cluster);
  point.key = key;
  point.layer = TrkrDefs::getLayer(key);
  point.x = global.x();
  point.y = global.y();
  point.z = global.z();
  point.r = std::hypot(point.x, point.y);
  point.phi = std::atan2(point.y, point.x);

  const auto trkrid = TrkrDefs::getTrkrId(key);
  if (trkrid == TrkrDefs::inttId)
  {
    point.ladder = InttDefs::getLadderPhiId(key);
    point.sensor = InttDefs::getLadderZId(key);
  }
  else if (trkrid == TrkrDefs::mvtxId)
  {
    point.ladder = MvtxDefs::getStaveId(key);
    point.sensor = MvtxDefs::getChipId(key);
  }
  else
  {
    point.ladder = TrkrDefs::getPhiElement(key);
    point.sensor = TrkrDefs::getZElement(key);
  }

  return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) && point.r > 0.0;
}

std::vector<Full_PolyTrackMatcher::SpacePoint> Full_PolyTrackMatcher::collectSiliconClusters() const
{
  std::vector<SpacePoint> points;
  if (!m_trkrClusters)
  {
    return points;
  }

  for (const TrkrDefs::TrkrId trkrid : {TrkrDefs::mvtxId, TrkrDefs::inttId})
  {
    for (const unsigned int layer : m_matchLayers)
    {
      const auto hitset_keys = m_trkrClusters->getHitSetKeys(trkrid, static_cast<uint8_t>(layer));
      for (const TrkrDefs::hitsetkey hitset_key : hitset_keys)
      {
        const auto range = m_trkrClusters->getClusters(hitset_key);
        for (auto iter = range.first; iter != range.second; ++iter)
        {
          SpacePoint point;
          if (getGlobalClusterPosition(iter->first, iter->second, point))
          {
            points.push_back(point);
          }
        }
      }
    }
  }
  std::sort(points.begin(), points.end(), [](const SpacePoint& a, const SpacePoint& b) {
    if (a.layer != b.layer) { return a.layer > b.layer; }
    return a.r > b.r;
  });
  return points;
}

std::map<unsigned int, std::vector<const Tpc_PolyCluster*>> Full_PolyTrackMatcher::collectTpcClustersByTrack() const
{
  std::map<unsigned int, std::vector<const Tpc_PolyCluster*>> clusters_by_track;
  if (!m_tpcClusters)
  {
    return clusters_by_track;
  }

  for (unsigned int i = 0; i < m_tpcClusters->size(); ++i)
  {
    const Tpc_PolyCluster* cluster = m_tpcClusters->get_cluster(i);
    if (cluster)
    {
      clusters_by_track[cluster->get_source_assembled_track_id()].push_back(cluster);
    }
  }
  return clusters_by_track;
}

const TpcCrossingDecision* Full_PolyTrackMatcher::findCrossingDecision(unsigned int source_assembled_track_id) const
{
  if (!m_crossingDecisions)
  {
    return nullptr;
  }
  for (unsigned int i = 0; i < m_crossingDecisions->size(); ++i)
  {
    const TpcCrossingDecision* decision = m_crossingDecisions->get_decision(i);
    if (decision && decision->get_assembled_track_id() == source_assembled_track_id)
    {
      return decision;
    }
  }
  return nullptr;
}

Full_PolyTrackMatcher::TrajectoryState Full_PolyTrackMatcher::fitTrajectory(const std::vector<SpacePoint>& points) const
{
  TrajectoryState state;
  if (points.size() < 2)
  {
    return state;
  }

  double sum_r = 0.0;
  double sum_r2 = 0.0;
  double sum_phi = 0.0;
  double sum_rphi = 0.0;
  double sum_z = 0.0;
  double sum_rz = 0.0;
  const double ref_phi = points.front().phi;
  for (const SpacePoint& point : points)
  {
    const double phi = unwrapPhiNear(point.phi, ref_phi);
    sum_r += point.r;
    sum_r2 += point.r * point.r;
    sum_phi += phi;
    sum_rphi += point.r * phi;
    sum_z += point.z;
    sum_rz += point.r * point.z;
  }

  const double n = static_cast<double>(points.size());
  const double denom = n * sum_r2 - sum_r * sum_r;
  if (std::fabs(denom) < 1.0e-12)
  {
    return state;
  }

  state.phi_slope = (n * sum_rphi - sum_r * sum_phi) / denom;
  state.phi_intercept = (sum_phi - state.phi_slope * sum_r) / n;
  state.z_slope = (n * sum_rz - sum_r * sum_z) / denom;
  state.z_intercept = (sum_z - state.z_slope * sum_r) / n;
  state.valid = std::isfinite(state.phi_slope) && std::isfinite(state.phi_intercept) &&
                std::isfinite(state.z_slope) && std::isfinite(state.z_intercept);
  return state;
}

bool Full_PolyTrackMatcher::predictAtRadius(const TrajectoryState& state, double r,
                                            double& pred_phi, double& pred_z,
                                            double& pred_x, double& pred_y) const
{
  if (!state.valid || !std::isfinite(r) || r <= 0.0)
  {
    return false;
  }
  pred_phi = wrapPhi(state.phi_intercept + state.phi_slope * r);
  pred_z = state.z_intercept + state.z_slope * r;
  pred_x = r * std::cos(pred_phi);
  pred_y = r * std::sin(pred_phi);
  return std::isfinite(pred_phi) && std::isfinite(pred_z) && std::isfinite(pred_x) && std::isfinite(pred_y);
}

Full_PolyTrackMatcher::Chain Full_PolyTrackMatcher::extendWithHit(const Chain& chain, const SpacePoint& point,
                                                                  double pred_phi, double pred_z,
                                                                  double pred_x, double pred_y,
                                                                  double rdphi, double dz, double chi2) const
{
  Chain out = chain;
  ChainHit hit;
  hit.point = point;
  hit.pred_phi = pred_phi;
  hit.pred_z = pred_z;
  hit.pred_x = pred_x;
  hit.pred_y = pred_y;
  hit.rdphi = rdphi;
  hit.dz = dz;
  hit.chi2 = chi2;
  out.hits.push_back(hit);
  out.fit_points.push_back(point);
  out.state = fitTrajectory(out.fit_points);
  out.chi2 += chi2;
  out.ndf += 2.0;
  out.score = out.chi2 + m_missingLayerPenalty * out.n_missing;
  return out;
}

Full_PolyTrackMatcher::Chain Full_PolyTrackMatcher::extendMissing(const Chain& chain, unsigned int layer) const
{
  Chain out = chain;
  out.missing_mask |= layerBit(layer);
  out.n_missing += 1U;
  out.score = out.chi2 + m_missingLayerPenalty * out.n_missing;
  return out;
}

void Full_PolyTrackMatcher::fillQaRow(const QaRow& row)
{
  if (!m_writeQA || !m_candidateTree)
  {
    return;
  }

  m_qaRow = row;
  m_candidateTree->Fill();
  if (row.accepted == 0)
  {
    return;
  }

  const unsigned int hist_key = row.layer * 10000U + row.ladder * 100U + row.sensor;
  auto make_hist = [&](std::map<unsigned int, TH2F*>& histmap, const char* prefix, const char* title,
                       int ny, double ymin, double ymax) -> TH2F* {
    auto found = histmap.find(hist_key);
    if (found != histmap.end())
    {
      return found->second;
    }
    std::ostringstream name;
    name << prefix << "_l" << row.layer << "_lad" << row.ladder << "_sen" << row.sensor;
    std::ostringstream htitle;
    htitle << title << " layer " << row.layer << " ladder " << row.ladder << " sensor " << row.sensor;
    TH2F* hist = new TH2F(name.str().c_str(), htitle.str().c_str(), 120, -15.0, 15.0, ny, ymin, ymax);
    histmap[hist_key] = hist;
    return hist;
  };
  auto make_phi_hist = [&](std::map<unsigned int, TH2F*>& histmap, const char* prefix, const char* title,
                           int ny, double ymin, double ymax) -> TH2F* {
    auto found = histmap.find(hist_key);
    if (found != histmap.end())
    {
      return found->second;
    }
    std::ostringstream name;
    name << prefix << "_l" << row.layer << "_lad" << row.ladder << "_sen" << row.sensor;
    std::ostringstream htitle;
    htitle << title << " layer " << row.layer << " ladder " << row.ladder << " sensor " << row.sensor;
    TH2F* hist = new TH2F(name.str().c_str(), htitle.str().c_str(), 128, -M_PI, M_PI, ny, ymin, ymax);
    histmap[hist_key] = hist;
    return hist;
  };

  make_hist(m_hRdphiVsZ, "h_rdphi_vs_z", "r*dphi vs z", 160, -m_looseRdphiWindow, m_looseRdphiWindow)->Fill(row.z, row.rdphi);
  make_hist(m_hDzVsZ, "h_dz_vs_z", "dz vs z", 160, -m_looseDzWindow, m_looseDzWindow)->Fill(row.z, row.dz);
  make_phi_hist(m_hRdphiVsPhi, "h_rdphi_vs_phi", "r*dphi vs phi", 160, -m_looseRdphiWindow, m_looseRdphiWindow)->Fill(row.phi, row.rdphi);
  make_phi_hist(m_hDzVsPhi, "h_dz_vs_phi", "dz vs phi", 160, -m_looseDzWindow, m_looseDzWindow)->Fill(row.phi, row.dz);
}

std::vector<Full_PolyTrackMatcher::Chain> Full_PolyTrackMatcher::buildChains(
    const Tpc_PolyTrack& track,
    const std::vector<const Tpc_PolyCluster*>& tpc_clusters,
    const std::vector<SpacePoint>& silicon_points)
{
  Chain seed;
  seed.fit_points.reserve(tpc_clusters.size() + m_matchLayers.size());
  for (const Tpc_PolyCluster* cluster : tpc_clusters)
  {
    if (!cluster)
    {
      continue;
    }
    SpacePoint point;
    point.key = cluster->get_trkr_cluster_key();
    point.layer = point.key != TrkrDefs::CLUSKEYMAX ? TrkrDefs::getLayer(point.key) : 0;
    point.x = cluster->get_centroid_x();
    point.y = cluster->get_centroid_y();
    point.z = cluster->get_centroid_z();
    point.r = std::hypot(point.x, point.y);
    point.phi = std::atan2(point.y, point.x);
    if (std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z) && point.r > 0.0)
    {
      seed.fit_points.push_back(point);
    }
  }

  if (seed.fit_points.size() < 2)
  {
    SpacePoint p0;
    p0.x = track.get_x();
    p0.y = track.get_y();
    p0.z = track.get_z();
    p0.r = std::hypot(p0.x, p0.y);
    p0.phi = std::atan2(p0.y, p0.x);
    SpacePoint p1 = p0;
    const double pt = std::hypot(track.get_px(), track.get_py());
    if (pt > 0.0)
    {
      p1.x += 20.0 * track.get_px() / pt;
      p1.y += 20.0 * track.get_py() / pt;
      p1.z += 20.0 * track.get_pz() / std::max(pt, 1.0e-9);
      p1.r = std::hypot(p1.x, p1.y);
      p1.phi = std::atan2(p1.y, p1.x);
    }
    if (p0.r > 0.0 && p1.r > 0.0)
    {
      seed.fit_points.push_back(p0);
      seed.fit_points.push_back(p1);
    }
  }

  seed.state = fitTrajectory(seed.fit_points);
  if (!seed.state.valid)
  {
    return {};
  }
  std::vector<Chain> chains{seed};

  const TpcCrossingDecision* crossing = findCrossingDecision(track.get_source_assembled_track_id());
  const int selected_crossing = crossing ? crossing->get_selected_crossing() : 0;

  for (unsigned int layer : m_matchLayers)
  {
    std::vector<Chain> next_chains;
    for (const Chain& chain : chains)
    {
      std::vector<Chain> layer_extensions;
      for (const SpacePoint& point : silicon_points)
      {
        if (point.layer != layer)
        {
          continue;
        }

        double pred_phi = 0.0;
        double pred_z = 0.0;
        double pred_x = 0.0;
        double pred_y = 0.0;
        if (!predictAtRadius(chain.state, point.r, pred_phi, pred_z, pred_x, pred_y))
        {
          continue;
        }

        const double dphi = wrapPhi(point.phi - pred_phi);
        const double rdphi = point.r * dphi;
        const double dz = point.z - pred_z;
        if (std::fabs(rdphi) > m_looseRdphiWindow || std::fabs(dz) > m_looseDzWindow)
        {
          continue;
        }

        const double chi2 = square(rdphi / std::max(m_sigmaRdphi, 1.0e-9)) +
                            square(dz / std::max(m_sigmaDz, 1.0e-9));
        layer_extensions.push_back(extendWithHit(chain, point, pred_phi, pred_z, pred_x, pred_y, rdphi, dz, chi2));

        QaRow row;
        row.event = m_event;
        row.tpc_track_id = track.get_track_id();
        row.source_assembled_track_id = track.get_source_assembled_track_id();
        row.crossing = selected_crossing;
        row.layer = point.layer;
        row.ladder = point.ladder;
        row.sensor = point.sensor;
        row.cluster_key = static_cast<unsigned long long>(point.key);
        row.x = finite_float(point.x);
        row.y = finite_float(point.y);
        row.z = finite_float(point.z);
        row.r = finite_float(point.r);
        row.phi = finite_float(point.phi);
        row.pred_x = finite_float(pred_x);
        row.pred_y = finite_float(pred_y);
        row.pred_z = finite_float(pred_z);
        row.pred_phi = finite_float(pred_phi);
        row.rdphi = finite_float(rdphi);
        row.dz = finite_float(dz);
        row.chi2 = finite_float(chi2);
        row.chain_chi2 = finite_float(chain.chi2 + chi2);
        row.chain_ndf = finite_float(chain.ndf + 2.0);
        row.chain_score = finite_float(chain.chi2 + chi2 + m_missingLayerPenalty * chain.n_missing);
        row.chain_silicon_clusters = static_cast<unsigned int>(chain.hits.size() + 1U);
        row.chain_missing_layers = chain.n_missing;
        row.chain_missing_mask = chain.missing_mask;
        row.accepted = 1;
        row.selected = 0;
        fillQaRow(row);
      }

      std::sort(layer_extensions.begin(), layer_extensions.end(), [](const Chain& a, const Chain& b) {
        if (a.hits.size() != b.hits.size()) { return a.hits.size() > b.hits.size(); }
        return a.score < b.score;
      });
      if (layer_extensions.size() > m_maxBranchesPerLayer)
      {
        layer_extensions.resize(m_maxBranchesPerLayer);
      }
      next_chains.insert(next_chains.end(), layer_extensions.begin(), layer_extensions.end());
      next_chains.push_back(extendMissing(chain, layer));
    }

    std::sort(next_chains.begin(), next_chains.end(), [](const Chain& a, const Chain& b) {
      if (a.hits.size() != b.hits.size()) { return a.hits.size() > b.hits.size(); }
      return a.score < b.score;
    });
    if (next_chains.size() > m_maxChains)
    {
      next_chains.resize(m_maxChains);
    }
    chains.swap(next_chains);
  }

  return chains;
}

const Full_PolyTrackMatcher::Chain* Full_PolyTrackMatcher::selectBestChain(const std::vector<Chain>& chains) const
{
  const Chain* best = nullptr;
  for (const Chain& chain : chains)
  {
    if (chain.hits.size() < m_minSiliconClusters)
    {
      continue;
    }
    if (!best ||
        chain.hits.size() > best->hits.size() ||
        (chain.hits.size() == best->hits.size() && chain.score < best->score))
    {
      best = &chain;
    }
  }
  return best;
}

void Full_PolyTrackMatcher::fillQaForChain(unsigned int tpc_track_id, unsigned int source_id, int crossing,
                                           const Chain& chain, bool selected)
{
  for (const ChainHit& hit : chain.hits)
  {
    QaRow row;
    row.event = m_event;
    row.tpc_track_id = tpc_track_id;
    row.source_assembled_track_id = source_id;
    row.crossing = crossing;
    row.layer = hit.point.layer;
    row.ladder = hit.point.ladder;
    row.sensor = hit.point.sensor;
    row.cluster_key = static_cast<unsigned long long>(hit.point.key);
    row.x = finite_float(hit.point.x);
    row.y = finite_float(hit.point.y);
    row.z = finite_float(hit.point.z);
    row.r = finite_float(hit.point.r);
    row.phi = finite_float(hit.point.phi);
    row.pred_x = finite_float(hit.pred_x);
    row.pred_y = finite_float(hit.pred_y);
    row.pred_z = finite_float(hit.pred_z);
    row.pred_phi = finite_float(hit.pred_phi);
    row.rdphi = finite_float(hit.rdphi);
    row.dz = finite_float(hit.dz);
    row.chi2 = finite_float(hit.chi2);
    row.chain_chi2 = finite_float(chain.chi2);
    row.chain_ndf = finite_float(chain.ndf);
    row.chain_score = finite_float(chain.score);
    row.chain_silicon_clusters = static_cast<unsigned int>(chain.hits.size());
    row.chain_missing_layers = chain.n_missing;
    row.chain_missing_mask = chain.missing_mask;
    row.accepted = 1;
    row.selected = selected ? 1 : 0;
    fillQaRow(row);
  }
}

void Full_PolyTrackMatcher::fillQaForUnmatchedSeed(const Tpc_PolyTrack& tpc_track, int crossing, const Chain& chain)
{
  QaRow row;
  row.event = m_event;
  row.tpc_track_id = tpc_track.get_track_id();
  row.source_assembled_track_id = tpc_track.get_source_assembled_track_id();
  row.crossing = crossing;
  row.layer = 999U;
  row.ladder = 999U;
  row.sensor = 999U;
  row.cluster_key = static_cast<unsigned long long>(TrkrDefs::CLUSKEYMAX);
  row.x = finite_float(tpc_track.get_seed_x0());
  row.y = finite_float(tpc_track.get_seed_y0());
  row.z = finite_float(tpc_track.get_seed_z0());
  row.r = finite_float(std::hypot(tpc_track.get_seed_x0(), tpc_track.get_seed_y0()));
  row.phi = finite_float(tpc_track.get_seed_phi());
  row.pred_x = row.x;
  row.pred_y = row.y;
  row.pred_z = row.z;
  row.pred_phi = row.phi;
  row.rdphi = 0.0F;
  row.dz = 0.0F;
  row.chi2 = 0.0F;
  row.chain_chi2 = finite_float(chain.chi2);
  row.chain_ndf = finite_float(chain.ndf);
  row.chain_score = finite_float(chain.score);
  row.chain_silicon_clusters = static_cast<unsigned int>(chain.hits.size());
  row.chain_missing_layers = chain.n_missing;
  row.chain_missing_mask = chain.missing_mask;
  row.accepted = 0;
  row.selected = 1;
  fillQaRow(row);
}

void Full_PolyTrackMatcher::fillTrack(const Tpc_PolyTrack& tpc_track, const Chain& chain)
{
  Full_PolyTrackv1* out = new Full_PolyTrackv1();
  out->set_event(m_event);
  out->set_track_id(m_fullTracks->size());
  out->set_tpc_poly_track_id(tpc_track.get_track_id());
  out->set_source_assembled_track_id(tpc_track.get_source_assembled_track_id());
  out->set_n_tpc_clusters(tpc_track.get_nclusters());
  out->set_n_missing_layers(chain.n_missing);
  out->set_missing_layer_mask(chain.missing_mask);
  out->set_fit_status((chain.state.valid || tpc_track.get_fit_status() != 0) ? 1 : 0);
  out->set_chi2(chain.chi2);
  out->set_ndf(chain.ndf);
  out->set_score(chain.score);
  out->set_charge(tpc_track.get_charge());
  out->set_seed_x(tpc_track.get_x());
  out->set_seed_y(tpc_track.get_y());
  out->set_seed_z(tpc_track.get_z());
  out->set_seed_px(tpc_track.get_px());
  out->set_seed_py(tpc_track.get_py());
  out->set_seed_pz(tpc_track.get_pz());

  if (chain.state.valid && !chain.hits.empty())
  {
    const double r0 = chain.hits.back().point.r;
    double phi = 0.0;
    double z = 0.0;
    double x = 0.0;
    double y = 0.0;
    predictAtRadius(chain.state, r0, phi, z, x, y);
    out->set_x(x);
    out->set_y(y);
    out->set_z(z);
    const double next_r = r0 + 1.0;
    double phi2 = 0.0;
    double z2 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    predictAtRadius(chain.state, next_r, phi2, z2, x2, y2);
    out->set_px(x2 - x);
    out->set_py(y2 - y);
    out->set_pz(z2 - z);
  }
  else
  {
    out->set_x(tpc_track.get_x());
    out->set_y(tpc_track.get_y());
    out->set_z(tpc_track.get_z());
    out->set_px(tpc_track.get_px());
    out->set_py(tpc_track.get_py());
    out->set_pz(tpc_track.get_pz());
  }

  for (const ChainHit& hit : chain.hits)
  {
    out->add_cluster_key(hit.point.key);
    out->add_silicon_state(hit.point.layer, hit.point.key,
                           finite_float(hit.point.x), finite_float(hit.point.y), finite_float(hit.point.z),
                           finite_float(hit.pred_x), finite_float(hit.pred_y), finite_float(hit.pred_z),
                           finite_float(hit.rdphi), finite_float(hit.dz), finite_float(hit.chi2));
  }

  m_fullTracks->add_track(out);
}

int Full_PolyTrackMatcher::process_event(PHCompositeNode* topNode)
{
  if (!m_tpcTracks || !m_tpcClusters || !m_trkrClusters || !m_actsGeometry || !m_fullTracks)
  {
    if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  m_fullTracks->Reset();
  const std::vector<SpacePoint> silicon_points = collectSiliconClusters();
  const auto tpc_clusters_by_track = collectTpcClustersByTrack();

  for (unsigned int i = 0; i < m_tpcTracks->size(); ++i)
  {
    const Tpc_PolyTrack* tpc_track = m_tpcTracks->get_track(i);
    if (!tpc_track)
    {
      continue;
    }
    auto cluster_iter = tpc_clusters_by_track.find(tpc_track->get_source_assembled_track_id());
    const std::vector<const Tpc_PolyCluster*> empty_clusters;
    const std::vector<const Tpc_PolyCluster*>& tpc_clusters =
        cluster_iter != tpc_clusters_by_track.end() ? cluster_iter->second : empty_clusters;
    std::vector<Chain> chains = buildChains(*tpc_track, tpc_clusters, silicon_points);
    const Chain* best = selectBestChain(chains);
    Chain unmatched_chain;
    const Chain* selected_chain = best;
    if (!selected_chain)
    {
      selected_chain = chains.empty() ? &unmatched_chain : &chains.front();
    }

    fillTrack(*tpc_track, *selected_chain);
    const TpcCrossingDecision* crossing = findCrossingDecision(tpc_track->get_source_assembled_track_id());
    const int selected_crossing = crossing ? crossing->get_selected_crossing() : 0;
    if (best)
    {
      fillQaForChain(tpc_track->get_track_id(), tpc_track->get_source_assembled_track_id(),
                     selected_crossing, *best, true);
    }
    else
    {
      fillQaForUnmatchedSeed(*tpc_track, selected_crossing, *selected_chain);
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " tpc_tracks=" << m_tpcTracks->size()
              << " silicon_clusters=" << silicon_points.size()
              << " full_tracks=" << m_fullTracks->size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Full_PolyTrackMatcher::End(PHCompositeNode*)
{
  if (m_qaFile)
  {
    m_qaFile->cd();
    if (m_candidateTree)
    {
      m_candidateTree->Write();
    }
    for (auto& item : m_hRdphiVsZ) { item.second->Write(); }
    for (auto& item : m_hDzVsZ) { item.second->Write(); }
    for (auto& item : m_hRdphiVsPhi) { item.second->Write(); }
    for (auto& item : m_hDzVsPhi) { item.second->Write(); }
    m_qaFile->Close();
    delete m_qaFile;
    m_qaFile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
