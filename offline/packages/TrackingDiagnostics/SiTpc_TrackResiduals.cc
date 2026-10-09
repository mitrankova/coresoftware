#include "SiTpc_TrackResiduals.h"

#include "tpctrackreco/Tpc_PolyCluster.h"
#include "tpctrackreco/Tpc_PolyClusterContainer.h"

#include <sitrackreco/SiTpcHelixFit.h>
#include <sitrackreco/SiTpc_Track.h>
#include <sitrackreco/SiTpc_TrackContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <format>
#include <iostream>
#include <limits>

namespace
{
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

  double wrap_phi(double phi)
  {
    const double pi = std::acos(-1.0);
    while (phi > pi)
    {
      phi -= 2.0 * pi;
    }
    while (phi <= -pi)
    {
      phi += 2.0 * pi;
    }
    return phi;
  }
}  // namespace

SiTpc_TrackResiduals::SiTpc_TrackResiduals(const std::string& name, const std::string& outfilename)
  : SubsysReco(name)
  , m_outfilename(outfilename)
{
}

SiTpc_TrackResiduals::~SiTpc_TrackResiduals()
{
  delete m_outfile;
}

int SiTpc_TrackResiduals::Init(PHCompositeNode* /*unused*/)
{
  // cppcheck-suppress publicAllocationError
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open output file " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("residuals", "Si+TPC track residuals");
  m_tree->Branch("event", &m_event, "event/i");
  m_tree->Branch("track_id", &m_trackId, "track_id/i");
  m_tree->Branch("si_chain_id", &m_siChainId, "si_chain_id/I");
  m_tree->Branch("si_trajectory_index", &m_siTrajIndex, "si_trajectory_index/i");
  m_tree->Branch("tpc_track_id", &m_tpcTrackId, "tpc_track_id/i");
  m_tree->Branch("tpc_assembled_track_id", &m_tpcAssembledId, "tpc_assembled_track_id/i");
  m_tree->Branch("fit_status", &m_fitStatus, "fit_status/I");
  m_tree->Branch("n_si_points", &m_nSi, "n_si_points/i");
  m_tree->Branch("n_tpc_points", &m_nTpc, "n_tpc_points/i");
  m_tree->Branch("match_dphi", &m_matchDphi, "match_dphi/D");
  m_tree->Branch("match_deta", &m_matchDeta, "match_deta/D");
  m_tree->Branch("match_dz0", &m_matchDz0, "match_dz0/D");
  m_tree->Branch("match_chi2", &m_matchChi2, "match_chi2/D");
  m_tree->Branch("si_phi", &m_siPhi, "si_phi/D");
  m_tree->Branch("si_eta", &m_siEta, "si_eta/D");
  m_tree->Branch("si_z0", &m_siZ0, "si_z0/D");
  m_tree->Branch("tpc_phi", &m_tpcPhi, "tpc_phi/D");
  m_tree->Branch("tpc_eta", &m_tpcEta, "tpc_eta/D");
  m_tree->Branch("tpc_z0", &m_tpcZ0, "tpc_z0/D");
  m_tree->Branch("tpc_pt", &m_tpcPt, "tpc_pt/D");
  m_tree->Branch("tpc_charge", &m_tpcCharge, "tpc_charge/I");
  m_tree->Branch("pt", &m_pt, "pt/D");
  m_tree->Branch("eta", &m_eta, "eta/D");
  m_tree->Branch("theta", &m_theta, "theta/D");
  m_tree->Branch("phi", &m_phi, "phi/D");
  m_tree->Branch("charge", &m_charge, "charge/I");
  m_tree->Branch("dca", &m_dca, "dca/D");
  m_tree->Branch("z0", &m_z0, "z0/D");
  m_tree->Branch("tanl", &m_tanl, "tanl/D");
  m_tree->Branch("R", &m_R, "R/D");
  m_tree->Branch("tpc_z_offset", &m_tpcZOffset, "tpc_z_offset/D");
  m_tree->Branch("circle_rms", &m_circleRms, "circle_rms/D");
  m_tree->Branch("z_rms", &m_zRms, "z_rms/D");
  m_tree->Branch("pca_x", &m_pcaX, "pca_x/D");
  m_tree->Branch("pca_y", &m_pcaY, "pca_y/D");
  m_tree->Branch("x0", &m_pcaX, "x0/D");
  m_tree->Branch("y0", &m_pcaY, "y0/D");
  m_tree->Branch("n_mvtx", &m_nMvtx, "n_mvtx/i");
  m_tree->Branch("n_intt", &m_nIntt, "n_intt/i");
  m_tree->Branch("n_tpc", &m_nTpcLayers, "n_tpc/i");
  m_tree->Branch("res_rphi_layer", m_resRPhiLayer, std::format("res_rphi_layer[{}]/D", kNLayers).c_str());
  m_tree->Branch("res_z_layer", m_resZLayer, std::format("res_z_layer[{}]/D", kNLayers).c_str());
  m_tree->Branch("rms_rphi_mvtx", &m_rmsRPhi[0], "rms_rphi_mvtx/D");
  m_tree->Branch("rms_rphi_intt", &m_rmsRPhi[1], "rms_rphi_intt/D");
  m_tree->Branch("rms_rphi_tpc", &m_rmsRPhi[2], "rms_rphi_tpc/D");
  m_tree->Branch("rms_z_mvtx", &m_rmsZ[0], "rms_z_mvtx/D");
  m_tree->Branch("rms_z_intt", &m_rmsZ[1], "rms_z_intt/D");
  m_tree->Branch("rms_z_tpc", &m_rmsZ[2], "rms_z_tpc/D");
  m_tree->Branch("point_type", &m_pointType);
  m_tree->Branch("layer", &m_layer);
  m_tree->Branch("source", &m_source);
  m_tree->Branch("point_x", &m_x);
  m_tree->Branch("point_y", &m_y);
  m_tree->Branch("point_z", &m_z);
  m_tree->Branch("point_r", &m_r);
  m_tree->Branch("point_phi", &m_pointPhi);
  m_tree->Branch("state_x", &m_stateX);
  m_tree->Branch("state_y", &m_stateY);
  m_tree->Branch("state_z", &m_stateZ);
  m_tree->Branch("state_r", &m_stateR);
  m_tree->Branch("state_phi", &m_statePhi);
  m_tree->Branch("delta_phi", &m_deltaPhi);
  m_tree->Branch("residual_rphi", &m_residualRPhi);
  m_tree->Branch("residual_z", &m_residualZ);
  m_tree->Branch("cluster_adc", &m_clusterAdc);
  m_tree->Branch("cluster_pad_size", &m_clusterPadSize);

  if (m_writeMatching)
  {
    m_matchTree = new TTree("matching", "Si-TPC match candidates");
    m_matchTree->Branch("event", &m_event, "event/i");
    m_matchTree->Branch("si_index", &m_cSi, "si_index/i");
    m_matchTree->Branch("tpc_index", &m_cTpc, "tpc_index/i");
    m_matchTree->Branch("dphi", &m_cDphi, "dphi/F");
    m_matchTree->Branch("deta", &m_cDeta, "deta/F");
    m_matchTree->Branch("dz0", &m_cDz0, "dz0/F");
    m_matchTree->Branch("chi2", &m_cChi2, "chi2/F");
    m_matchTree->Branch("accepted", &m_cAccepted, "accepted/I");
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void SiTpc_TrackResiduals::reset_tree_values()
{
  m_event = m_evt;
  m_trackId = 0;
  m_siChainId = -1;
  m_siTrajIndex = 0;
  m_tpcTrackId = 0;
  m_tpcAssembledId = 0;
  m_fitStatus = 0;
  m_nSi = 0;
  m_nTpc = 0;
  m_matchDphi = m_matchDeta = m_matchDz0 = m_matchChi2 = kNaN;
  m_siPhi = m_siEta = m_siZ0 = kNaN;
  m_tpcPhi = m_tpcEta = m_tpcZ0 = m_tpcPt = kNaN;
  m_tpcCharge = 0;
  m_pt = m_eta = m_theta = m_phi = m_dca = m_z0 = m_tanl = kNaN;
  m_charge = 0;
  m_R = m_tpcZOffset = m_circleRms = m_zRms = m_pcaX = m_pcaY = kNaN;
  m_nMvtx = m_nIntt = m_nTpcLayers = 0;
  for (int l = 0; l < kNLayers; ++l)
  {
    m_resRPhiLayer[l] = kNaN;
    m_resZLayer[l] = kNaN;
  }
  for (int d = 0; d < 3; ++d)
  {
    m_rmsRPhi[d] = kNaN;
    m_rmsZ[d] = kNaN;
  }
  for (auto* v : {&m_x, &m_y, &m_z, &m_r, &m_pointPhi, &m_stateX, &m_stateY, &m_stateZ, &m_stateR, &m_statePhi,
                  &m_deltaPhi, &m_residualRPhi, &m_residualZ, &m_clusterAdc})
  {
    v->clear();
  }
  m_pointType.clear();
  m_layer.clear();
  m_source.clear();
  m_clusterPadSize.clear();
}

int SiTpc_TrackResiduals::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  m_tracks = findNode::getClass<SiTpc_TrackContainer>(topNode, m_trackNodeName);
  if (!m_tracks)
  {
    if (Verbosity() > 0)
    {
      std::cerr << Name() << " - missing " << m_trackNodeName << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_clusterNodeName);

  // ---- matching tree
  m_event = m_evt;
  if (m_matchTree)
  {
    for (unsigned int i = 0; i < m_tracks->size_candidates(); ++i)
    {
      m_cSi = m_tracks->get_candidate_si(i);
      m_cTpc = m_tracks->get_candidate_tpc(i);
      m_cDphi = m_tracks->get_candidate_dphi(i);
      m_cDeta = m_tracks->get_candidate_deta(i);
      m_cDz0 = m_tracks->get_candidate_dz0(i);
      m_cChi2 = m_tracks->get_candidate_chi2(i);
      m_cAccepted = m_tracks->get_candidate_accepted(i);
      m_matchTree->Fill();
    }
  }

  // ---- residuals tree
  unsigned int nfilled = 0;
  for (unsigned int it = 0; it < m_tracks->size(); ++it)
  {
    const SiTpc_Track* t = m_tracks->get_track(it);
    if (!t || !t->isValid())
    {
      continue;
    }
    const bool curved = t->get_fit_status() == SiTpcHelixFit::Ok;
    if (curved && (t->get_pt() < m_minPt || t->get_pt() > m_maxPt))
    {
      continue;
    }
    if (t->get_n_tpc_points() < m_minTpcClusters || t->get_n_si_points() < m_minSiPoints)
    {
      continue;
    }

    reset_tree_values();
    m_trackId = t->get_id();
    m_siChainId = t->get_si_chain_id();
    m_siTrajIndex = t->get_si_trajectory_index();
    m_tpcTrackId = t->get_tpc_track_id();
    m_tpcAssembledId = t->get_tpc_assembled_track_id();
    m_fitStatus = t->get_fit_status();
    m_nSi = t->get_n_si_points();
    m_nTpc = t->get_n_tpc_points();
    m_matchDphi = t->get_match_dphi();
    m_matchDeta = t->get_match_deta();
    m_matchDz0 = t->get_match_dz0();
    m_matchChi2 = t->get_match_chi2();
    m_siPhi = t->get_si_phi();
    m_siEta = t->get_si_eta();
    m_siZ0 = t->get_si_z0();
    m_tpcPhi = t->get_tpc_phi();
    m_tpcEta = t->get_tpc_eta();
    m_tpcZ0 = t->get_tpc_z0();
    m_tpcPt = t->get_tpc_pt();
    m_tpcCharge = t->get_tpc_charge();
    m_pt = t->get_pt();
    m_eta = t->get_eta();
    m_theta = t->get_theta();
    m_phi = t->get_phi();
    m_charge = t->get_charge();
    m_dca = t->get_dca();
    m_z0 = t->get_z0();
    m_tanl = t->get_tanl();
    m_R = t->get_radius();
    m_tpcZOffset = t->get_tpc_z_offset();
    m_circleRms = t->get_circle_rms();
    m_zRms = t->get_z_rms();
    m_pcaX = t->get_pca_x();
    m_pcaY = t->get_pca_y();

    const SiTpcHelixFit::Result helix = SiTpcHelixFit::fromTrack(*t);
    double layerSumRPhi[kNLayers] = {}, layerSumZ[kNLayers] = {};
    int layerN[kNLayers] = {};
    double detSumRPhi[3] = {}, detSumZ[3] = {};
    int detN[3] = {};
    for (unsigned int k = 0; k < t->size_points(); ++k)
    {
      const int type = t->get_point_type(k);
      const double x = t->get_point_x(k), y = t->get_point_y(k);
      const double zOff = type == SiTpc_Track::TpcPoint ? m_tpcZOffset : 0.0;
      const double z = t->get_point_z(k) - zOff;  // in the fit's z convention
      const double r = std::hypot(x, y);
      SiTpcHelixFit::Pos st;
      double s = 0;
      const bool ok = SiTpcHelixFit::atRadius(helix, r, z, st, s);

      m_pointType.push_back(type);
      m_layer.push_back(t->get_point_layer(k));
      m_source.push_back(t->get_point_source(k));
      m_x.push_back(x);
      m_y.push_back(y);
      m_z.push_back(t->get_point_z(k));
      m_r.push_back(r);
      const double phiPoint = std::atan2(y, x);
      m_pointPhi.push_back(phiPoint);
      if (ok)
      {
        const double phiState = std::atan2(st.y, st.x);
        const double dphi = wrap_phi(phiPoint - phiState);
        m_stateX.push_back(st.x);
        m_stateY.push_back(st.y);
        m_stateZ.push_back(st.z + zOff);  // in the point's z convention
        m_stateR.push_back(std::hypot(st.x, st.y));
        m_statePhi.push_back(phiState);
        m_deltaPhi.push_back(dphi);
        m_residualRPhi.push_back(r * dphi);
        m_residualZ.push_back(z - st.z);

        const int layer = t->get_point_layer(k);
        const int det = type == SiTpc_Track::TpcPoint ? 2 : (layer <= 2 ? 0 : 1);
        if (layer >= 0 && layer < kNLayers)
        {
          layerSumRPhi[layer] += r * dphi;
          layerSumZ[layer] += z - st.z;
          ++layerN[layer];
        }
        detSumRPhi[det] += r * dphi * r * dphi;
        detSumZ[det] += (z - st.z) * (z - st.z);
        ++detN[det];
      }
      else
      {
        for (auto* v : {&m_stateX, &m_stateY, &m_stateZ, &m_stateR, &m_statePhi, &m_deltaPhi, &m_residualRPhi, &m_residualZ})
        {
          v->push_back(kNaN);
        }
      }

      const Tpc_PolyCluster* c = (type == SiTpc_Track::TpcPoint && m_clusters) ? m_clusters->get_cluster(t->get_point_source(k)) : nullptr;
      m_clusterAdc.push_back(c ? c->get_adc() : kNaN);
      m_clusterPadSize.push_back(c ? c->get_phi_width() : 0);
      ++nfilled;

      // point counts per detector (Si: layer 0-2 MVTX, 3-6 INTT)
      if (type == SiTpc_Track::TpcPoint)
      {
        ++m_nTpcLayers;
      }
      else
      {
        ++(t->get_point_layer(k) <= 2 ? m_nMvtx : m_nIntt);
      }
    }
    for (int l = 0; l < kNLayers; ++l)
    {
      if (layerN[l] > 0)
      {
        m_resRPhiLayer[l] = layerSumRPhi[l] / layerN[l];
        m_resZLayer[l] = layerSumZ[l] / layerN[l];
      }
    }
    for (int d = 0; d < 3; ++d)
    {
      if (detN[d] > 0)
      {
        m_rmsRPhi[d] = std::sqrt(detSumRPhi[d] / detN[d]);
        m_rmsZ[d] = std::sqrt(detSumZ[d] / detN[d]);
      }
    }
    m_tree->Fill();
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_evt << " Si+TPC tracks=" << m_tracks->size()
              << " residuals=" << nfilled << " candidates=" << m_tracks->size_candidates() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpc_TrackResiduals::End(PHCompositeNode* /*unused*/)
{
  if (m_outfile)
  {
    m_outfile->cd();
    if (m_tree)
    {
      m_tree->Write();
    }
    if (m_matchTree)
    {
      m_matchTree->Write();
    }
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
    m_tree = nullptr;
    m_matchTree = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}