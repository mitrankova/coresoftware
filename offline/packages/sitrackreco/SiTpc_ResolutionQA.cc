#include "SiTpc_ResolutionQA.h"

#include "SiTpc_Track.h"
#include "SiTpc_TrackContainer.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <Eigen/Dense>

#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <vector>

namespace
{
  constexpr float fnan = std::numeric_limits<float>::quiet_NaN();

  enum Source
  {
    Own = 0,
    Svtx = 1
  };

  // a track near the beam, in the frame of its source
  struct QTrack
  {
    unsigned int id = 0;
    short crossing = SHRT_MAX;
    double x = 0, y = 0, z = 0;  // PCA to the beam line [cm]
    double phi = 0, tanl = 0, pt = 0;
    int charge = 0;
    int nmvtx = 0, nintt = 0, ntpc = 0;
    float sigD0 = fnan;  // reported by the source (svtx only)
    bool good = false;   // passes the vertex-track cuts
    double wxy = 0, wz = 0;
    // results
    float dcaBeam = fnan, dcaVtxXY = fnan, dcaVtxZ = fnan, sigVtxXY = fnan, sigVtxZ = fnan;
    int vtxNtrk = 0;
    double eta() const { return std::asinh(tanl); }
  };

  struct QVertex
  {
    short crossing = SHRT_MAX;
    double x = 0, y = 0, z = 0;
    std::array<double, 3> err{0, 0, 0};
    std::vector<int> tracks;  // indices into the QTrack vector
  };

  struct VertexFit
  {
    bool ok = false;
    Eigen::Vector3d v = Eigen::Vector3d::Zero();
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    double chi2 = 0;
    std::vector<int> used;
  };

  Eigen::Vector3d normalXY(const QTrack& t) { return {-std::sin(t.phi), std::cos(t.phi), 0}; }
  Eigen::Vector3d gradZ(const QTrack& t) { return {-t.tanl * std::cos(t.phi), -t.tanl * std::sin(t.phi), 1}; }

  // residuals of a track at the vertex v: xy (signed, Acts convention) and z (track - vertex)
  std::pair<double, double> residuals(const QTrack& t, const Eigen::Vector3d& v)
  {
    const Eigen::Vector3d p(t.x, t.y, t.z);
    const double rxy = normalXY(t).dot(p - v);
    const double s = std::cos(t.phi) * (v.x() - t.x) + std::sin(t.phi) * (v.y() - t.y);
    const double rz = (t.z + t.tanl * s) - v.z();
    return {rxy, rz};
  }

  // straight-line vertex fit, iterative removal of the worst track above chi2Cut
  VertexFit fitVertex(const std::vector<QTrack>& tracks, std::vector<int> idx, double chi2Cut, int minTracks)
  {
    VertexFit fit;
    while (static_cast<int>(idx.size()) >= minTracks && idx.size() >= 2)
    {
      Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
      Eigen::Vector3d b = Eigen::Vector3d::Zero();
      for (const int i : idx)
      {
        const auto& t = tracks[i];
        const Eigen::Vector3d p(t.x, t.y, t.z);
        const Eigen::Vector3d n = normalXY(t), g = gradZ(t);
        const double h = t.z - t.tanl * (std::cos(t.phi) * t.x + std::sin(t.phi) * t.y);
        A += t.wxy * n * n.transpose() + t.wz * g * g.transpose();
        b += t.wxy * n * n.dot(p) + t.wz * g * h;
      }
      const double det = A.determinant();
      if (!(std::abs(det) > 0) || !std::isfinite(det))
      {
        return fit;
      }
      const Eigen::Matrix3d C = A.inverse();
      const Eigen::Vector3d v = C * b;
      double chi2 = 0, worst = -1;
      std::size_t iworst = 0;
      for (std::size_t k = 0; k < idx.size(); ++k)
      {
        const auto& t = tracks[idx[k]];
        const auto [rxy, rz] = residuals(t, v);
        const double c = t.wxy * rxy * rxy + t.wz * rz * rz;
        chi2 += c;
        if (c > worst)
        {
          worst = c;
          iworst = k;
        }
      }
      if (worst > chi2Cut && static_cast<int>(idx.size()) > minTracks && idx.size() > 2)
      {
        idx.erase(idx.begin() + static_cast<long>(iworst));
        continue;
      }
      fit.ok = true;
      fit.v = v;
      fit.cov = C;
      fit.chi2 = chi2;
      fit.used = idx;
      return fit;
    }
    return fit;
  }
}  // namespace

//____________________________________________________________________________________________________
struct SiTpc_ResolutionQA::Rows
{
  struct Track
  {
    int src = 0, event = 0, id = 0, crossing = 0, charge = 0, nmvtx = 0, nintt = 0, ntpc = 0, good = 0, vtx_ntrk = 0;
    float pt = fnan, eta = fnan, phi = fnan, x = fnan, y = fnan, z = fnan;
    float dca_beam = fnan, sig_d0 = fnan, dca_vtx_xy = fnan, dca_vtx_z = fnan, sig_vtx_xy = fnan, sig_vtx_z = fnan;
  } t;
  struct Vertex
  {
    int src = 0, event = 0, crossing = 0, ntrk_finder = 0, ntrk_good = 0, ntrk_fit = 0, nhalf = 0;
    float fx = fnan, fy = fnan, fz = fnan, fex = fnan, fey = fnan, fez = fnan;  // vertex finder
    float qx = fnan, qy = fnan, qz = fnan, qex = fnan, qey = fnan, qez = fnan, qchi2 = fnan;  // QA refit
    float beamx = fnan, beamy = fnan;  // beam line at qz
    float split_dx = fnan, split_dy = fnan, split_dz = fnan;  // (v1 - v2) / sqrt(2)
  } v;
  struct Match
  {
    int event = 0, crossing = 0, own_id = 0, svtx_id = 0, own_good = 0, svtx_good = 0;
    float own_pt = fnan, svtx_pt = fnan, eta = fnan, phi = fnan;
    float own_dca_beam = fnan, svtx_dca_beam = fnan, own_dca_vtx_xy = fnan, svtx_dca_vtx_xy = fnan;
    float own_dca_vtx_z = fnan, svtx_dca_vtx_z = fnan, own_sig_vtx_xy = fnan, svtx_sig_vtx_xy = fnan;
  } m;
};

//____________________________________________________________________________________________________
SiTpc_ResolutionQA::SiTpc_ResolutionQA(const std::string& name, const std::string& fileName)
  : SubsysReco(name)
  , m_rows(std::make_unique<Rows>())
  , m_fileName(fileName)
{
}

SiTpc_ResolutionQA::~SiTpc_ResolutionQA() = default;

//____________________________________________________________________________________________________
int SiTpc_ResolutionQA::InitRun(PHCompositeNode* /*topNode*/)
{
  m_file = new TFile(m_fileName.c_str(), "RECREATE");
  auto& t = m_rows->t;
  m_trackTree = new TTree("tracks", "src 0 own, 1 svtx; dca [cm]");
  m_trackTree->Branch("src", &t.src);
  m_trackTree->Branch("event", &t.event);
  m_trackTree->Branch("id", &t.id);
  m_trackTree->Branch("crossing", &t.crossing);
  m_trackTree->Branch("charge", &t.charge);
  m_trackTree->Branch("nmvtx", &t.nmvtx);
  m_trackTree->Branch("nintt", &t.nintt);
  m_trackTree->Branch("ntpc", &t.ntpc);
  m_trackTree->Branch("good", &t.good);
  m_trackTree->Branch("vtx_ntrk", &t.vtx_ntrk);
  m_trackTree->Branch("pt", &t.pt);
  m_trackTree->Branch("eta", &t.eta);
  m_trackTree->Branch("phi", &t.phi);
  m_trackTree->Branch("x", &t.x);
  m_trackTree->Branch("y", &t.y);
  m_trackTree->Branch("z", &t.z);
  m_trackTree->Branch("dca_beam", &t.dca_beam);
  m_trackTree->Branch("sig_d0", &t.sig_d0);
  m_trackTree->Branch("dca_vtx_xy", &t.dca_vtx_xy);
  m_trackTree->Branch("dca_vtx_z", &t.dca_vtx_z);
  m_trackTree->Branch("sig_vtx_xy", &t.sig_vtx_xy);
  m_trackTree->Branch("sig_vtx_z", &t.sig_vtx_z);

  auto& v = m_rows->v;
  m_vertexTree = new TTree("vertices", "src 0 own, 1 svtx; f*: vertex finder, q*: QA refit, split_*: (v1 - v2) / sqrt(2) [cm]");
  m_vertexTree->Branch("src", &v.src);
  m_vertexTree->Branch("event", &v.event);
  m_vertexTree->Branch("crossing", &v.crossing);
  m_vertexTree->Branch("ntrk_finder", &v.ntrk_finder);
  m_vertexTree->Branch("ntrk_good", &v.ntrk_good);
  m_vertexTree->Branch("ntrk_fit", &v.ntrk_fit);
  m_vertexTree->Branch("nhalf", &v.nhalf);
  m_vertexTree->Branch("fx", &v.fx);
  m_vertexTree->Branch("fy", &v.fy);
  m_vertexTree->Branch("fz", &v.fz);
  m_vertexTree->Branch("fex", &v.fex);
  m_vertexTree->Branch("fey", &v.fey);
  m_vertexTree->Branch("fez", &v.fez);
  m_vertexTree->Branch("qx", &v.qx);
  m_vertexTree->Branch("qy", &v.qy);
  m_vertexTree->Branch("qz", &v.qz);
  m_vertexTree->Branch("qex", &v.qex);
  m_vertexTree->Branch("qey", &v.qey);
  m_vertexTree->Branch("qez", &v.qez);
  m_vertexTree->Branch("qchi2", &v.qchi2);
  m_vertexTree->Branch("beamx", &v.beamx);
  m_vertexTree->Branch("beamy", &v.beamy);
  m_vertexTree->Branch("split_dx", &v.split_dx);
  m_vertexTree->Branch("split_dy", &v.split_dy);
  m_vertexTree->Branch("split_dz", &v.split_dz);

  auto& m = m_rows->m;
  m_matchTree = new TTree("matched", "own and svtx track of the same particle");
  m_matchTree->Branch("event", &m.event);
  m_matchTree->Branch("crossing", &m.crossing);
  m_matchTree->Branch("own_id", &m.own_id);
  m_matchTree->Branch("svtx_id", &m.svtx_id);
  m_matchTree->Branch("own_good", &m.own_good);
  m_matchTree->Branch("svtx_good", &m.svtx_good);
  m_matchTree->Branch("own_pt", &m.own_pt);
  m_matchTree->Branch("svtx_pt", &m.svtx_pt);
  m_matchTree->Branch("eta", &m.eta);
  m_matchTree->Branch("phi", &m.phi);
  m_matchTree->Branch("own_dca_beam", &m.own_dca_beam);
  m_matchTree->Branch("svtx_dca_beam", &m.svtx_dca_beam);
  m_matchTree->Branch("own_dca_vtx_xy", &m.own_dca_vtx_xy);
  m_matchTree->Branch("svtx_dca_vtx_xy", &m.svtx_dca_vtx_xy);
  m_matchTree->Branch("own_dca_vtx_z", &m.own_dca_vtx_z);
  m_matchTree->Branch("svtx_dca_vtx_z", &m.svtx_dca_vtx_z);
  m_matchTree->Branch("own_sig_vtx_xy", &m.own_sig_vtx_xy);
  m_matchTree->Branch("svtx_sig_vtx_xy", &m.svtx_sig_vtx_xy);

  std::cout << Name() << ": own " << m_ownTrackNode << " / " << m_ownVertexNode << " beam (" << m_ownBeam[0] << ", " << m_ownBeam[1]
            << ", slopes " << m_ownBeam[2] << ", " << m_ownBeam[3] << "); svtx " << m_svtxTrackNode << " / " << m_svtxVertexNode
            << " beam (" << m_svtxBeam[0] << ", " << m_svtxBeam[1] << ", slopes " << m_svtxBeam[2] << ", " << m_svtxBeam[3] << ")" << std::endl
            << "  vertex tracks: pt > " << m_minPt << ", MVTX >= " << m_minMvtx << ", TPC >= " << m_minTpc
            << "; errors xy " << m_errModel[0] << " (+) " << m_errModel[1] << "/p, z " << m_errModel[2] << " (+) " << m_errModel[3]
            << "/p cm; chi2 cut " << m_chi2Cut << "; output " << m_fileName << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpc_ResolutionQA::process_event(PHCompositeNode* topNode)
{
  auto* crossings = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingNode);

  std::array<std::vector<QTrack>, 2> tracks;
  std::array<std::vector<QVertex>, 2> vertices;
  std::array<bool, 2> present{false, false};

  // ---- own tracks
  if (auto* own = findNode::getClass<SiTpc_TrackContainer>(topNode, m_ownTrackNode))
  {
    present[Own] = true;
    std::map<unsigned int, int> byId;
    for (unsigned int i = 0; i < own->size(); ++i)
    {
      const SiTpc_Track* s = own->get_track(i);
      if (!s || s->get_fit_status() != 1 || s->get_charge() == 0)  // SiTpcHelixFit::Ok
      {
        continue;
      }
      QTrack q;
      q.id = s->get_id();
      q.x = s->get_pca_x();
      q.y = s->get_pca_y();
      q.z = s->get_z0();
      q.phi = s->get_phi();
      q.tanl = s->get_tanl();
      q.pt = s->get_pt();
      q.charge = s->get_charge();
      for (unsigned int k = 0; k < s->size_points(); ++k)
      {
        const int layer = s->get_point_layer(k);
        if (s->get_point_type(k) == SiTpc_Track::SiPoint)
        {
          ++(layer <= 2 ? q.nmvtx : q.nintt);
        }
        else
        {
          ++q.ntpc;
        }
      }
      if (crossings)
      {
        const TpcCrossingDecision* dec = crossings->get_decision(s->get_tpc_assembled_track_id());
        q.crossing = dec ? dec->get_selected_crossing() : SHRT_MAX;
      }
      byId[q.id] = static_cast<int>(tracks[Own].size());
      tracks[Own].push_back(q);
    }
    if (auto* vmap = findNode::getClass<SvtxVertexMap>(topNode, m_ownVertexNode))
    {
      for (const auto& [key, vtx] : *vmap)
      {
        QVertex qv;
        qv.crossing = vtx->get_beam_crossing();
        qv.x = vtx->get_x();
        qv.y = vtx->get_y();
        qv.z = vtx->get_z();
        qv.err = {std::sqrt(vtx->get_error(0, 0)), std::sqrt(vtx->get_error(1, 1)), std::sqrt(vtx->get_error(2, 2))};
        for (auto it = vtx->begin_tracks(); it != vtx->end_tracks(); ++it)
        {
          const auto f = byId.find(*it);
          if (f != byId.end())
          {
            qv.tracks.push_back(f->second);
          }
        }
        vertices[Own].push_back(qv);
      }
    }
  }

  // ---- svtx tracks
  if (auto* tmap = findNode::getClass<SvtxTrackMap>(topNode, m_svtxTrackNode))
  {
    present[Svtx] = true;
    std::map<unsigned int, int> byId;
    for (const auto& [key, trk] : *tmap)
    {
      if (!trk)
      {
        continue;
      }
      QTrack q;
      q.id = key;
      q.crossing = trk->get_crossing();
      q.x = trk->get_x();
      q.y = trk->get_y();
      q.z = trk->get_z();
      q.phi = std::atan2(trk->get_py(), trk->get_px());
      q.pt = std::hypot(trk->get_px(), trk->get_py());
      q.tanl = q.pt > 0 ? trk->get_pz() / q.pt : 0;
      q.charge = trk->get_charge();
      if (const TrackSeed* si = trk->get_silicon_seed())
      {
        for (auto it = si->begin_cluster_keys(); it != si->end_cluster_keys(); ++it)
        {
          const auto id = TrkrDefs::getTrkrId(*it);
          q.nmvtx += id == TrkrDefs::mvtxId;
          q.nintt += id == TrkrDefs::inttId;
        }
      }
      if (const TrackSeed* tpc = trk->get_tpc_seed())
      {
        for (auto it = tpc->begin_cluster_keys(); it != tpc->end_cluster_keys(); ++it)
        {
          q.ntpc += TrkrDefs::getTrkrId(*it) == TrkrDefs::tpcId;
        }
      }
      const Eigen::Vector3d n = normalXY(q);
      Eigen::Matrix2d cxy;
      cxy << trk->get_error(0, 0), trk->get_error(0, 1), trk->get_error(1, 0), trk->get_error(1, 1);
      const double var = n.head<2>().dot(cxy * n.head<2>());
      q.sigD0 = var > 0 ? static_cast<float>(std::sqrt(var)) : fnan;
      byId[q.id] = static_cast<int>(tracks[Svtx].size());
      tracks[Svtx].push_back(q);
    }
    if (auto* vmap = findNode::getClass<SvtxVertexMap>(topNode, m_svtxVertexNode))
    {
      for (const auto& [key, vtx] : *vmap)
      {
        QVertex qv;
        qv.crossing = vtx->get_beam_crossing();
        qv.x = vtx->get_x();
        qv.y = vtx->get_y();
        qv.z = vtx->get_z();
        qv.err = {std::sqrt(vtx->get_error(0, 0)), std::sqrt(vtx->get_error(1, 1)), std::sqrt(vtx->get_error(2, 2))};
        for (auto it = vtx->begin_tracks(); it != vtx->end_tracks(); ++it)
        {
          const auto f = byId.find(*it);
          if (f != byId.end())
          {
            qv.tracks.push_back(f->second);
          }
        }
        vertices[Svtx].push_back(qv);
      }
    }
  }

  // ---- per source: DCA to the beam, QA vertex refits, unbiased DCA, split vertices
  for (int src : {Own, Svtx})
  {
    if (!present[src])
    {
      continue;
    }
    const auto& beam = src == Own ? m_ownBeam : m_svtxBeam;
    auto& trks = tracks[src];
    for (auto& t : trks)
    {
      const double p = t.pt * std::sqrt(1 + t.tanl * t.tanl);
      t.wxy = 1. / (m_errModel[0] * m_errModel[0] + std::pow(m_errModel[1] / p, 2));
      t.wz = 1. / (m_errModel[2] * m_errModel[2] + std::pow(m_errModel[3] / p, 2));
      t.good = std::isfinite(t.x) && std::isfinite(t.phi) && std::isfinite(t.tanl) && t.pt >= m_minPt && t.nmvtx >= m_minMvtx && t.ntpc >= m_minTpc;
      const Eigen::Vector3d b(beam[0] + beam[2] * t.z, beam[1] + beam[3] * t.z, t.z);
      t.dcaBeam = static_cast<float>(normalXY(t).dot(Eigen::Vector3d(t.x, t.y, t.z) - b));
    }

    std::mt19937 rng(static_cast<unsigned int>(1000003 * m_event + 17 * src + 1));
    for (const auto& qv : vertices[src])
    {
      std::vector<int> good;
      for (const int i : qv.tracks)
      {
        if (trks[i].good)
        {
          good.push_back(i);
        }
      }
      auto& v = m_rows->v;
      v = Rows::Vertex();
      v.src = src;
      v.event = m_event;
      v.crossing = qv.crossing;
      v.ntrk_finder = static_cast<int>(qv.tracks.size());
      v.ntrk_good = static_cast<int>(good.size());
      v.fx = qv.x;
      v.fy = qv.y;
      v.fz = qv.z;
      v.fex = qv.err[0];
      v.fey = qv.err[1];
      v.fez = qv.err[2];

      const VertexFit full = fitVertex(trks, good, m_chi2Cut, 2);
      if (full.ok)
      {
        v.ntrk_fit = static_cast<int>(full.used.size());
        v.qx = full.v.x();
        v.qy = full.v.y();
        v.qz = full.v.z();
        v.qex = std::sqrt(full.cov(0, 0));
        v.qey = std::sqrt(full.cov(1, 1));
        v.qez = std::sqrt(full.cov(2, 2));
        v.qchi2 = full.chi2;
        v.beamx = beam[0] + beam[2] * full.v.z();
        v.beamy = beam[1] + beam[3] * full.v.z();
      }

      // unbiased DCA of every good track of this vertex
      for (const int i : good)
      {
        std::vector<int> others;
        for (const int j : good)
        {
          if (j != i)
          {
            others.push_back(j);
          }
        }
        if (static_cast<int>(others.size()) < m_minOther)
        {
          continue;
        }
        const VertexFit fo = fitVertex(trks, others, m_chi2Cut, m_minOther);
        if (!fo.ok)
        {
          continue;
        }
        auto& t = trks[i];
        const auto [rxy, rz] = residuals(t, fo.v);
        const Eigen::Vector3d n = normalXY(t), g = gradZ(t);
        t.dcaVtxXY = static_cast<float>(rxy);
        t.dcaVtxZ = static_cast<float>(rz);
        t.sigVtxXY = static_cast<float>(std::sqrt(n.dot(fo.cov * n)));
        t.sigVtxZ = static_cast<float>(std::sqrt(g.dot(fo.cov * g)));
        t.vtxNtrk = static_cast<int>(fo.used.size());
      }

      // split vertex
      if (full.ok && static_cast<int>(full.used.size()) >= m_minSplit)
      {
        std::vector<int> shuffled = full.used;
        std::shuffle(shuffled.begin(), shuffled.end(), rng);
        const std::size_t half = shuffled.size() / 2;
        const std::vector<int> a(shuffled.begin(), shuffled.begin() + static_cast<long>(half));
        const std::vector<int> b(shuffled.begin() + static_cast<long>(half), shuffled.begin() + static_cast<long>(2 * half));
        // no outlier removal in the halves: both keep n / 2 tracks
        const VertexFit fa = fitVertex(trks, a, 1e30, 2);
        const VertexFit fb = fitVertex(trks, b, 1e30, 2);
        if (fa.ok && fb.ok)
        {
          const Eigen::Vector3d d = (fa.v - fb.v) / std::sqrt(2.);
          v.nhalf = static_cast<int>(half);
          v.split_dx = d.x();
          v.split_dy = d.y();
          v.split_dz = d.z();
        }
      }
      m_vertexTree->Fill();
    }

    for (const auto& t : trks)
    {
      auto& r = m_rows->t;
      r.src = src;
      r.event = m_event;
      r.id = static_cast<int>(t.id);
      r.crossing = t.crossing;
      r.charge = t.charge;
      r.nmvtx = t.nmvtx;
      r.nintt = t.nintt;
      r.ntpc = t.ntpc;
      r.good = t.good ? 1 : 0;
      r.vtx_ntrk = t.vtxNtrk;
      r.pt = t.pt;
      r.eta = t.eta();
      r.phi = t.phi;
      r.x = t.x;
      r.y = t.y;
      r.z = t.z;
      r.dca_beam = t.dcaBeam;
      r.sig_d0 = t.sigD0;
      r.dca_vtx_xy = t.dcaVtxXY;
      r.dca_vtx_z = t.dcaVtxZ;
      r.sig_vtx_xy = t.sigVtxXY;
      r.sig_vtx_z = t.sigVtxZ;
      m_trackTree->Fill();
    }
  }

  // ---- own <-> svtx matching (same crossing, closest in phi / eta / pt)
  if (present[Own] && present[Svtx])
  {
    std::vector<bool> usedSvtx(tracks[Svtx].size(), false);
    for (const auto& o : tracks[Own])
    {
      int best = -1;
      double bestScore = std::numeric_limits<double>::max();
      for (std::size_t j = 0; j < tracks[Svtx].size(); ++j)
      {
        const auto& s = tracks[Svtx][j];
        if (usedSvtx[j] || s.crossing != o.crossing)
        {
          continue;
        }
        const double dphi = std::remainder(s.phi - o.phi, 2 * M_PI);
        const double deta = s.eta() - o.eta();
        const double dpt = (s.pt - o.pt) / o.pt;
        if (std::abs(dphi) > m_match[0] || std::abs(deta) > m_match[1] || std::abs(dpt) > m_match[2])
        {
          continue;
        }
        const double score = std::pow(dphi / m_match[0], 2) + std::pow(deta / m_match[1], 2) + std::pow(dpt / m_match[2], 2);
        if (score < bestScore)
        {
          bestScore = score;
          best = static_cast<int>(j);
        }
      }
      if (best < 0)
      {
        continue;
      }
      usedSvtx[best] = true;
      const auto& s = tracks[Svtx][best];
      auto& m = m_rows->m;
      m.event = m_event;
      m.crossing = o.crossing;
      m.own_id = static_cast<int>(o.id);
      m.svtx_id = static_cast<int>(s.id);
      m.own_good = o.good;
      m.svtx_good = s.good;
      m.own_pt = o.pt;
      m.svtx_pt = s.pt;
      m.eta = o.eta();
      m.phi = o.phi;
      m.own_dca_beam = o.dcaBeam;
      m.svtx_dca_beam = s.dcaBeam;
      m.own_dca_vtx_xy = o.dcaVtxXY;
      m.svtx_dca_vtx_xy = s.dcaVtxXY;
      m.own_dca_vtx_z = o.dcaVtxZ;
      m.svtx_dca_vtx_z = s.dcaVtxZ;
      m.own_sig_vtx_xy = o.sigVtxXY;
      m.svtx_sig_vtx_xy = s.sigVtxXY;
      m_matchTree->Fill();
    }
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpc_ResolutionQA::End(PHCompositeNode* /*topNode*/)
{
  if (m_file)
  {
    m_file->cd();
    m_trackTree->Write();
    m_vertexTree->Write();
    m_matchTree->Write();
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }
  std::cout << Name() << "::End - " << m_event << " events written to " << m_fileName << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
