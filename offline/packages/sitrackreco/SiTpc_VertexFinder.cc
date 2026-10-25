#include "SiTpc_VertexFinder.h"

#include "SiHitSeedData.h"
#include "SiTpcHelixFit.h"
#include "SiTpc_Track.h"
#include "SiTpc_TrackContainer.h"
#include "Si_Trajectory.h"
#include "Si_TrajectoryContainer.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <globalvertex/SvtxVertexMap_v1.h>
#include <globalvertex/SvtxVertex_v3.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

namespace
{
  bool good(double x)
  {
    return std::isfinite(x) && std::fabs(x) < 1.0e30;
  }

  const char* sourceName(SiTpc_VertexFinder::Source s)
  {
    switch (s)
    {
    case SiTpc_VertexFinder::Source::Si:
      return "Si only (SI_TRAJECTORY)";
    case SiTpc_VertexFinder::Source::Tpc:
      return "TPC only (TPC_POLYTRACKS)";
    default:
      return "Si+TPC (SITPC_TRACKS)";
    }
  }
}  // namespace

//____________________________________________________________________________
SiTpc_VertexFinder::SiTpc_VertexFinder(const std::string& name, Source source)
  : SubsysReco(name)
  , m_source(source)
{
  // defaults per source; every setter overrides them
  PairVertexFinder::Config& c = m_finder.config();
  switch (source)
  {
  case Source::Si:
    m_inputNodeName = "SI_TRAJECTORY";
    m_vertexMapName = "SiPolyVertexMap";
    m_minPt = 0.0;  // short lever arm: pt not reliable, straight-line fits allowed
    c.pairDcaCut = 0.1;
    c.outlierPairCut = 0.2;
    c.sigmaRPhi = 0.02;
    c.sigmaZ = 0.05;
    c.msTerm = 0.003;
    c.maxDcaXY = 1.0;
    break;
  case Source::Tpc:
    m_inputNodeName = "TPC_POLYTRACKS";
    m_vertexMapName = "TpcPolyVertexMap";
    m_minPt = 0.1;
    c.pairDcaCut = 0.3;
    c.outlierPairCut = 0.5;
    c.sigmaRPhi = 0.1;
    c.sigmaZ = 0.15;
    c.msTerm = 0.05;
    c.maxDcaXY = 5.0;
    break;
  case Source::SiTpc:
    m_inputNodeName = "SITPC_TRACKS";
    m_vertexMapName = "SiTpcPolyVertexMap";
    m_minPt = 0.1;
    c.pairDcaCut = 0.1;
    c.outlierPairCut = 0.2;
    c.sigmaRPhi = 0.02;
    c.sigmaZ = 0.05;
    c.msTerm = 0.003;
    c.maxDcaXY = 1.0;
    break;
  }
}

//____________________________________________________________________________
int SiTpc_VertexFinder::InitRun(PHCompositeNode* topNode)
{
  if (!createNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  const auto& c = m_finder.config();
  std::cout << Name() << "::InitRun - " << sourceName(m_source) << " " << m_inputNodeName << " -> " << m_vertexMapName
            << (m_groupByCrossing ? ", per bunch crossing (" + m_decisionNodeName + ", tier <= " + std::to_string(m_maxDecisionTier) + ")" : std::string(", all crossings together"))
            << std::endl
            << Name() << "::InitRun - pairs: dca < " << c.pairDcaCut << " cm, outlier cut " << c.outlierPairCut << " cm, >= "
            << c.minTracksPerVertex << " tracks; fit: sigma rphi " << c.sigmaRPhi << " cm, z " << c.sigmaZ << " cm, ms " << c.msTerm
            << " cm GeV (lines: " << m_lineSigmaRPhi << ", " << m_lineSigmaZ << " cm), track chi2 < " << c.maxTrackChi2 << "; beam (" << c.beamX << ", " << c.beamY << ") cm, |dca_xy| < "
            << c.maxDcaXY << " cm, |z| < " << c.maxAbsZ << " cm" << (c.useBeamSpotConstraint ? ", beam-spot constraint" : "")
            << "; output shift (" << m_shift[0] << ", " << m_shift[1] << ", " << m_shift[2] << ") cm"
            << (m_source == Source::SiTpc ? (m_siTpcPositionFromSi ? "; SiTpc position from the Si points" : "; SiTpc position from the Si+TPC fit") : "")
            << (m_source != Source::Tpc ? "; tracks with >= " + std::to_string(m_minMvtx) + " MVTX points" : std::string())
            << (m_source == Source::SiTpc && m_maxMatchDphi > 0 ? ", |match dphi| < " + std::to_string(m_maxMatchDphi) + " rad" : std::string())
            << (m_minPt > 0 ? ", pt > " + std::to_string(m_minPt) + " GeV" : std::string())
            << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
bool SiTpc_VertexFinder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator dstIter(dstNode);
  auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);
  if (!m_vertexMap)
  {
    m_vertexMap = new SvtxVertexMap_v1;
    svtxNode->addNode(new PHIODataNode<PHObject>(m_vertexMap, m_vertexMapName, "PHObject"));
  }
  return true;
}

//____________________________________________________________________________
bool SiTpc_VertexFinder::crossingOf(unsigned int assembledId, short& crossing) const
{
  const TpcCrossingDecision* d = m_decisions ? m_decisions->get_decision(assembledId) : nullptr;
  if (!d || static_cast<int>(d->get_selected_tier()) > m_maxDecisionTier)
  {
    return false;
  }
  crossing = d->get_selected_crossing();
  return true;
}

//____________________________________________________________________________
void SiTpc_VertexFinder::add(PairVertexFinder::TrackModel m, bool hasCrossing, short crossing)
{
  if (m.is_line && m.sigmaRPhi <= 0)
  {
    m.sigmaRPhi = m_lineSigmaRPhi;
    m.sigmaZ = m_lineSigmaZ;
  }
  if (!m_finder.acceptAtBeam(m))
  {
    return;
  }
  if (m_groupByCrossing && !hasCrossing)
  {
    ++m_nNoCrossing;
    if (m_requireCrossing)
    {
      return;
    }
  }
  ++m_nSelected;
  m_groups[m_groupByCrossing && hasCrossing ? crossing : kNoCrossing].push_back(m);
}

//____________________________________________________________________________
void SiTpc_VertexFinder::collectSi(PHCompositeNode* topNode)
{
  auto* trajs = findNode::getClass<Si_TrajectoryContainer>(topNode, m_inputNodeName);
  if (!trajs)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_inputNodeName << std::endl;
    }
    return;
  }
  auto* seeds = findNode::getClass<SiHitSeedEvent>(topNode, m_siSeedNodeName);
  for (unsigned int i = 0; i < trajs->size(); ++i)
  {
    const Si_Trajectory* t = trajs->get_trajectory(i);
    if (!t)
    {
      continue;
    }
    ++m_nInput;
    const int status = t->get_fit_status();
    if (!(status == Si_Trajectory::Ok || (m_useLines && status == Si_Trajectory::StraightLine)))
    {
      continue;
    }
    if (t->get_n_mvtx() + t->get_n_intt() < m_minSiPoints || t->get_n_mvtx() < m_minMvtx)
    {
      continue;
    }
    const double pt = t->get_pt();
    if (good(pt) && pt < m_minPt)
    {
      continue;
    }
    PairVertexFinder::TrackModel m;
    if (!PairVertexFinder::fromPerigee(t->get_id(), t->get_pca_x(), t->get_pca_y(), t->get_z0(), t->get_phi(), t->get_tanl(),
                                       t->get_circle_x(), t->get_circle_y(), t->get_radius(), t->get_helicity(), pt, m))
    {
      continue;
    }
    // crossing: decision of the chain's TPC track, else the chain's own crossing
    bool has = false;
    short crossing = 0;
    const int chainId = t->get_chain_id();
    if (seeds && chainId >= 0)
    {
      const SiHitChain* chain = nullptr;
      if (static_cast<std::size_t>(chainId) < seeds->chains.size() && seeds->chains[static_cast<std::size_t>(chainId)].id == chainId)
      {
        chain = &seeds->chains[static_cast<std::size_t>(chainId)];
      }
      else
      {
        for (const auto& c : seeds->chains)
        {
          if (c.id == chainId)
          {
            chain = &c;
            break;
          }
        }
      }
      if (chain)
      {
        has = crossingOf(chain->tpc_assembled_track_id, crossing);
        if (!has && !m_decisions)
        {
          crossing = chain->tpc_crossing;
          has = true;
        }
      }
    }
    add(m, has, crossing);
  }
}

//____________________________________________________________________________
void SiTpc_VertexFinder::collectTpc(PHCompositeNode* topNode)
{
  auto* tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_inputNodeName);
  if (!tracks)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_inputNodeName << std::endl;
    }
    return;
  }
  for (unsigned int i = 0; i < tracks->size(); ++i)
  {
    const Tpc_PolyTrack* trk = tracks->get_track(i);
    if (!trk)
    {
      continue;
    }
    ++m_nInput;
    if (trk->get_fit_status() <= 0 || trk->get_nclusters() < m_minTpcClusters)
    {
      continue;
    }
    if (trk->get_ndf() > 0.0 && trk->get_chi2() / trk->get_ndf() > m_maxChi2Ndf)
    {
      continue;
    }
    const double x = trk->get_x(), y = trk->get_y(), z = trk->get_z();
    const double px = trk->get_px(), py = trk->get_py(), pz = trk->get_pz(), q = trk->get_charge();
    if (!good(x) || !good(y) || !good(z) || !good(px) || !good(py) || !good(pz) || !good(q))
    {
      continue;
    }
    const double pt = std::hypot(px, py);
    if (pt <= 1.0e-9)
    {
      continue;
    }
    // as Tpc_PolyTrackVertexFinder::buildModel
    PairVertexFinder::TrackModel m;
    m.track_id = trk->get_track_id();
    m.x0 = x;
    m.y0 = y;
    m.z0 = z;
    m.phi0 = std::atan2(py, px);
    m.tanl = pz / pt;
    m.is_line = std::fabs(q) < 0.5 || std::fabs(m_bz) < 1.0e-9;
    if (m.is_line)
    {
      if (!m_useLines)
      {
        continue;
      }
      m.dir = Eigen::Vector3d(px, py, pz).normalized();
      m.p = 0.0;
    }
    else
    {
      if (pt < m_minPt)
      {
        continue;
      }
      m.p = pt * std::sqrt(1.0 + m.tanl * m.tanl);
      const double rMom = pt / (0.003 * std::fabs(m_bz));
      const double sphi = std::sin(m.phi0), cphi = std::cos(m.phi0);
      bool useStored = false;
      const double cx = trk->get_helix_x0(), cy = trk->get_helix_y0();
      if (good(cx) && good(cy) && (std::fabs(cx) + std::fabs(cy)) > 0.0)
      {
        const double dx = cx - x, dy = cy - y;
        const double rg = std::hypot(dx, dy);
        const double along = dx * cphi + dy * sphi;  // vanishes: centre is normal to the tangent
        if (rg > 0.0 && std::fabs(along) < 0.05 * rg && std::fabs(rg - rMom) < 0.2 * rMom)
        {
          useStored = true;
          m.xc = cx;
          m.yc = cy;
          m.radius = rg;
          m.h = (-dx * sphi + dy * cphi) >= 0.0 ? 1.0 : -1.0;
        }
      }
      if (!useStored)
      {
        m.radius = rMom;
        m.h = (q * m_bz > 0.0) ? -1.0 : 1.0;  // positive charge in +Bz bends clockwise
        m.xc = x - m.radius * m.h * sphi;
        m.yc = y + m.radius * m.h * cphi;
      }
      if (!good(m.radius) || m.radius <= 0.0)
      {
        continue;
      }
    }
    if (trk->get_cov(0, 0) > 0.0 && trk->get_cov(1, 1) > 0.0 && trk->get_cov(2, 2) > 0.0)
    {
      m.has_cov = true;
      for (unsigned int a = 0; a < 3; ++a)
      {
        for (unsigned int b = 0; b < 3; ++b)
        {
          m.cov(a, b) = trk->get_cov(a, b);
        }
      }
    }
    short crossing = 0;
    const bool has = crossingOf(trk->get_source_assembled_track_id(), crossing);
    add(m, has, crossing);
  }
}

//____________________________________________________________________________
void SiTpc_VertexFinder::collectSiTpc(PHCompositeNode* topNode)
{
  auto* tracks = findNode::getClass<SiTpc_TrackContainer>(topNode, m_inputNodeName);
  if (!tracks)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_inputNodeName << std::endl;
    }
    return;
  }
  for (unsigned int i = 0; i < tracks->size(); ++i)
  {
    const SiTpc_Track* t = tracks->get_track(i);
    if (!t || !t->isValid())
    {
      continue;
    }
    ++m_nInput;
    const int status = t->get_fit_status();
    if (!(status == SiTpcHelixFit::Ok || (m_useLines && status == SiTpcHelixFit::StraightLine)))
    {
      continue;
    }
    if (t->get_n_si_points() < m_minSiPoints || t->get_n_tpc_points() < m_minTpcClusters)
    {
      continue;
    }
    // quality of the Si seed: MVTX points (2-point chains define position and direction with no redundancy)
    unsigned int nMvtx = 0;
    for (unsigned int k = 0; k < t->size_points(); ++k)
    {
      nMvtx += t->get_point_type(k) == SiTpc_Track::SiPoint && t->get_point_layer(k) >= 0 && t->get_point_layer(k) <= 2;
    }
    if (nMvtx < m_minMvtx)
    {
      ++m_nRejMvtx;
      continue;
    }
    if (m_maxMatchDphi > 0.0 && !(std::fabs(t->get_match_dphi()) < m_maxMatchDphi))
    {
      ++m_nRejMatch;
      continue;
    }
    const double pt = t->get_pt();
    if (good(pt) && pt < m_minPt)
    {
      continue;
    }
    PairVertexFinder::TrackModel m;
    bool ok = false;
    if (m_siTpcPositionFromSi)
    {
      std::vector<Eigen::Vector3d> pts;
      for (unsigned int k = 0; k < t->size_points(); ++k)
      {
        if (t->get_point_type(k) == SiTpc_Track::SiPoint)
        {
          pts.emplace_back(t->get_point_x(k), t->get_point_y(k), t->get_point_z(k));
        }
      }
      ok = PairVertexFinder::fromPointsWithCurvature(t->get_id(), pts, t->get_radius(), t->get_helicity(), t->get_tanl(), pt, m);
      m_nSiPosition += ok;
    }
    if (!ok && !PairVertexFinder::fromPerigee(t->get_id(), t->get_pca_x(), t->get_pca_y(), t->get_z0(), t->get_phi(), t->get_tanl(),
                                              t->get_circle_x(), t->get_circle_y(), t->get_radius(), t->get_helicity(), pt, m))
    {
      continue;
    }
    short crossing = 0;
    const bool has = crossingOf(t->get_tpc_assembled_track_id(), crossing);
    add(m, has, crossing);
  }
}

//____________________________________________________________________________
int SiTpc_VertexFinder::process_event(PHCompositeNode* topNode)
{
  ++m_nEvents;
  if (!m_vertexMap && !createNodes(topNode))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_vertexMap->Reset();
  m_groups.clear();
  m_decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  if (m_groupByCrossing && !m_decisions && m_source != Source::Si && Verbosity() > 0)
  {
    std::cout << Name() << ": no " << m_decisionNodeName << " - tracks have no crossing" << std::endl;
  }

  switch (m_source)
  {
  case Source::Si:
    collectSi(topNode);
    break;
  case Source::Tpc:
    collectTpc(topNode);
    break;
  case Source::SiTpc:
    collectSiTpc(topNode);
    break;
  }

  unsigned int nVtx = 0;
  for (const auto& [crossing, models] : m_groups)
  {
    ++m_nGroups;
    PairVertexFinder::Stats stats;
    const auto vertices = m_finder.find(models, &stats);
    for (const auto& v : vertices)
    {
      auto svtx = std::make_unique<SvtxVertex_v3>();
      svtx->set_x(static_cast<float>(v.position.x() + m_shift[0]));
      svtx->set_y(static_cast<float>(v.position.y() + m_shift[1]));
      svtx->set_z(static_cast<float>(v.position.z() + m_shift[2]));
      svtx->set_t(0.0);
      svtx->set_chisq(static_cast<float>(v.chi2));
      svtx->set_ndof(static_cast<unsigned int>(std::max(0, v.ndf)));
      svtx->set_beam_crossing(crossing);
      for (unsigned int a = 0; a < 3; ++a)
      {
        for (unsigned int b = 0; b < 3; ++b)
        {
          svtx->set_error(a, b, static_cast<float>(v.covariance(a, b)));
        }
      }
      for (const unsigned int it : v.tracks)
      {
        svtx->insert_track(models[it].track_id);
      }
      m_vertexMap->insert(svtx.release());  // the map assigns the vertex id
      ++nVtx;
      ++m_nVertices;
      m_nTracksInVertices += v.tracks.size();
      ++m_multiplicity[v.tracks.size()];
      m_sigmaZ.push_back(std::sqrt(std::max(0.0, v.covariance(2, 2))));
      if (Verbosity() > 1)
      {
        std::cout << Name() << ": crossing " << crossing << " vertex (" << v.position.x() << ", " << v.position.y() << ", "
                  << v.position.z() << ") +- (" << std::sqrt(v.covariance(0, 0)) << ", " << std::sqrt(v.covariance(1, 1)) << ", "
                  << std::sqrt(v.covariance(2, 2)) << ") tracks " << v.tracks.size() << " of " << models.size() << " chi2/ndf "
                  << v.chi2 << "/" << v.ndf << (stats.looseRetry ? " (3x pair cut)" : "") << std::endl;
      }
    }
  }
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << m_nEvents << " selected tracks " << m_nSelected << " (cumulative), crossings "
              << m_groups.size() << ", vertices " << nVtx << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
int SiTpc_VertexFinder::End(PHCompositeNode* /*unused*/)
{
  std::cout << Name() << "::End - " << sourceName(m_source) << ": " << m_nEvents << " events, input tracks " << m_nInput
            << ", selected " << m_nSelected << " (without a usable crossing " << m_nNoCrossing << (m_requireCrossing ? ", skipped" : "")
            << "), crossings with tracks " << m_nGroups << ", vertices " << m_nVertices << ", tracks in vertices "
            << m_nTracksInVertices << std::endl;
  if (m_source == Source::SiTpc)
  {
    std::cout << Name() << "::End - SiTpc track selection: rejected with < " << m_minMvtx << " MVTX points " << m_nRejMvtx
              << ", with |match dphi| >= " << m_maxMatchDphi << " rad " << m_nRejMatch << std::endl;
  }
  if (m_source == Source::SiTpc && m_siTpcPositionFromSi)
  {
    std::cout << Name() << "::End - tracks positioned by their Si points with the Si+TPC curvature " << m_nSiPosition << std::endl;
  }
  std::cout << Name() << "::End - tracks per vertex:";
  for (const auto& [n, count] : m_multiplicity)
  {
    std::cout << " " << n << ":" << count;
  }
  if (!m_sigmaZ.empty())
  {
    std::vector<double> s = m_sigmaZ;
    std::nth_element(s.begin(), s.begin() + static_cast<std::ptrdiff_t>(s.size() / 2), s.end());
    std::cout << std::endl << Name() << "::End - median vertex z error " << s[s.size() / 2] << " cm";
  }
  std::cout << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
