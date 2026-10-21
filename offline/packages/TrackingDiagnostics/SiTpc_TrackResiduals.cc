#include "SiTpc_TrackResiduals.h"

#include "tpctrackreco/TpcCrossingDecision.h"
#include "tpctrackreco/TpcCrossingDecisionContainer.h"
#include "tpctrackreco/Tpc_PolyCluster.h"
#include "tpctrackreco/Tpc_PolyClusterContainer.h"

#include <sitrackreco/SiHitSeedData.h>
#include <sitrackreco/SiTpcHelixFit.h>
#include <sitrackreco/SiTpc_Track.h>
#include <sitrackreco/SiTpc_TrackContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/InttDefs.h>

#include <TFile.h>
#include <TH3.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <map>

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

  int charge_index(int charge)
  {
    if (charge > 0)
    {
      return 1;
    }
    if (charge < 0)
    {
      return 0;
    }
    return -1;
  }

  // nominal Si layer radii [cm], as in SiHitSeedReco
  constexpr double kSiRadius[7] = {2.523, 3.336, 4.148, 7.188 - 0.0036, 7.732 - 0.0036, 9.680 - 0.0036, 10.262 - 0.0036};

  constexpr const char* kChargeName[2] = {"neg", "pos"};
  constexpr const char* kChargeTitle[2] = {"q<0", "q>0"};
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
  m_tree->Branch("tpc_crossing", &m_tpcCrossing, "tpc_crossing/I");
  m_tree->Branch("tpc_crossing_tier", &m_tpcCrossingTier, "tpc_crossing_tier/I");
  m_tree->Branch("intt_crossing", &m_inttCrossing, "intt_crossing/I");
  m_tree->Branch("intt_n_hits", &m_inttNHits, "intt_n_hits/I");
  m_tree->Branch("intt_crossing_spread", &m_inttCrossingSpread, "intt_crossing_spread/I");
  m_tree->Branch("crossing_diff", &m_crossingDiff, "crossing_diff/I");
  m_tree->Branch("n_tracks_crossing", &m_nTracksCrossing, "n_tracks_crossing/I");
  m_tree->Branch("n_si_chains_crossing", &m_nSiChainsCrossing, "n_si_chains_crossing/I");
  m_tree->Branch("track_vertex_z", &m_trackVertexZ, "track_vertex_z/D");
  m_tree->Branch("vertex_z", &m_vertexZ, "vertex_z/D");
  m_tree->Branch("vertex_dz", &m_vertexDz, "vertex_dz/D");
  m_tree->Branch("vertex_score", &m_vertexScore, "vertex_score/D");
  m_tree->Branch("vertex_status", &m_vertexStatus, "vertex_status/I");
  m_tree->Branch("vertex_own_z", &m_vertexOwnZ, "vertex_own_z/I");
  m_tree->Branch("vertex_n_si", &m_vertexNSi, "vertex_n_si/I");
  m_tree->Branch("vertex_z_mean", &m_vertexZMean, "vertex_z_mean/D");
  m_tree->Branch("vertex_z_rms", &m_vertexZRms, "vertex_z_rms/D");
  m_tree->Branch("vertex_z_err", &m_vertexZErr, "vertex_z_err/D");
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
  m_tree->Branch("dphi_layer", m_dphiLayer, std::format("dphi_layer[{}]/D", kNLayers).c_str());
  m_tree->Branch("dphi_prev_layer", m_dphiPrevLayer, std::format("dphi_prev_layer[{}]/D", kNLayers).c_str());
  m_tree->Branch("prev_layer", m_prevLayer, std::format("prev_layer[{}]/I", kNLayers).c_str());
  m_tree->Branch("ddphi_layer", m_ddphiLayer, std::format("ddphi_layer[{}]/D", kNLayers).c_str());
  m_tree->Branch("hw_has_chain", &m_hwHasChain, "hw_has_chain/I");
  m_tree->Branch("hw_prev_layer", m_hwPrevLayer, std::format("hw_prev_layer[{}]/I", kNSiLayers).c_str());
  m_tree->Branch("hw_step_kind", m_hwStepKind, std::format("hw_step_kind[{}]/I", kNSiLayers).c_str());
  m_tree->Branch("hw_uphi", m_hwUphi, std::format("hw_uphi[{}]/D", kNSiLayers).c_str());
  m_tree->Branch("hw_vz", m_hwVz, std::format("hw_vz[{}]/D", kNSiLayers).c_str());
  m_tree->Branch("hw_dphi_step", m_hwDphiStep, std::format("hw_dphi_step[{}]/D", kNSiLayers).c_str());
  m_tree->Branch("hw_dphi_pred", m_hwDphiPred, std::format("hw_dphi_pred[{}]/D", kNSiLayers).c_str());
  m_tree->Branch("hw_ddphi", m_hwDdphi, std::format("hw_ddphi[{}]/D", kNSiLayers).c_str());
  m_tree->Branch("hw_dz_step", m_hwDzStep, std::format("hw_dz_step[{}]/D", kNSiLayers).c_str());
  m_tree->Branch("hw_res_phi", m_hwResPhi, std::format("hw_res_phi[{}]/D", kNSiLayers).c_str());
  m_tree->Branch("hw_res_z", m_hwResZ, std::format("hw_res_z[{}]/D", kNSiLayers).c_str());
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

  if (m_writeVertices)
  {
    m_vertexTree = new TTree("vertices", "collision vertices (one per decided crossing)");
    m_vertexTree->Branch("event", &m_event, "event/i");
    m_vertexTree->Branch("crossing", &m_vCrossing, "crossing/I");
    m_vertexTree->Branch("vertex_z", &m_vZ, "vertex_z/D");
    m_vertexTree->Branch("n_tracks", &m_vNTracks, "n_tracks/I");
    m_vertexTree->Branch("n_si", &m_vNSi, "n_si/I");
    m_vertexTree->Branch("n_tier0", &m_vNTier0, "n_tier0/I");
    m_vertexTree->Branch("z_mean", &m_vZMean, "z_mean/D");
    m_vertexTree->Branch("z_rms", &m_vZRms, "z_rms/D");
    m_vertexTree->Branch("z_err", &m_vZErr, "z_err/D");
    m_vertexTree->Branch("track_z", &m_vTrackZ);
    m_vertexTree->Branch("track_own_z", &m_vOwnZ);
    m_vertexTree->Branch("track_tier", &m_vTier);
    m_vertexTree->Branch("track_eta", &m_vEta);
    m_vertexTree->Branch("track_pt", &m_vPt);
    m_vertexTree->Branch("track_phi", &m_vPhi);
    m_vertexTree->Branch("track_charge", &m_vCharge);
  }

  if (m_writeDdphiHists)
  {
    book_ddphi_hists();
  }
  if (m_writeHwHists)
  {
    book_hw_hists();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void SiTpc_TrackResiduals::book_hw_hists()
{
  m_outfile->cd();
  for (int q = 0; q < 2; ++q)
  {
    m_h3HwDphiStepVsLayer[q] = new TH3F(
        std::format("h3_hw_dphi_step_vs_layer_{}", kChargeName[q]).c_str(),
        std::format("hardware {}: #DeltaU_{{#phi}} step into layer L;U_{{#phi}}(L) - U_{{#phi}}(prev) [rad];layer L;p_{{T}}^{{TPC}} [GeV]",
                    kChargeTitle[q])
            .c_str(),
        m_hwStepBins, -m_hwStepMax, m_hwStepMax,
        kNSiLayers, -0.5, kNSiLayers - 0.5,
        m_ptBins, m_ptMin, m_ptMax);
    for (int l = 0; l < kNSiLayers; ++l)
    {
      m_h3HwDphiVsPred[l][q] = new TH3F(
          std::format("h3_hw_dphi_vs_pred_l{}_{}", l, kChargeName[q]).c_str(),
          std::format("hardware layer {} {}: #DeltaU_{{#phi}} step vs prediction from previous step (const d#phi/dr);"
                      "#DeltaU_{{#phi}}^{{pred}} [rad];#DeltaU_{{#phi}}(L{}) [rad];p_{{T}}^{{TPC}} [GeV]",
                      l, kChargeTitle[q], l)
              .c_str(),
          m_hwStepBins, -m_hwStepMax, m_hwStepMax,
          m_hwStepBins, -m_hwStepMax, m_hwStepMax,
          m_ptBins, m_ptMin, m_ptMax);
    }
    m_h3HwDdphiVsLayer[q] = new TH3F(
        std::format("h3_hw_ddphi_vs_layer_{}", kChargeName[q]).c_str(),
        std::format("hardware {}: #Delta#DeltaU_{{#phi}} = #DeltaU_{{#phi}} - #DeltaU_{{#phi}}^{{pred}};#Delta#DeltaU_{{#phi}} [rad];layer L;p_{{T}}^{{TPC}} [GeV]",
                    kChargeTitle[q])
            .c_str(),
        m_hwDdphiBins, -m_hwDdphiMax, m_hwDdphiMax,
        kNSiLayers, -0.5, kNSiLayers - 0.5,
        m_ptBins, m_ptMin, m_ptMax);
    m_h3HwResPhiVsLayer[q] = new TH3F(
        std::format("h3_hw_res_phi_vs_layer_{}", kChargeName[q]).c_str(),
        std::format("hardware {}: chain step residual U_{{#phi}} - pred;U_{{#phi}} - U_{{#phi}}^{{pred}} [rad];layer L;p_{{T}}^{{TPC}} [GeV]",
                    kChargeTitle[q])
            .c_str(),
        m_hwResBins, -m_hwResPhiMax, m_hwResPhiMax,
        kNSiLayers, -0.5, kNSiLayers - 0.5,
        m_ptBins, m_ptMin, m_ptMax);
    m_h3HwResZVsLayer[q] = new TH3F(
        std::format("h3_hw_res_z_vs_layer_{}", kChargeName[q]).c_str(),
        std::format("hardware {}: chain step residual V_{{z}} - pred;V_{{z}} - V_{{z}}^{{pred}} [cm];layer L;p_{{T}}^{{TPC}} [GeV]",
                    kChargeTitle[q])
            .c_str(),
        m_hwResBins, -m_hwResZMax, m_hwResZMax,
        kNSiLayers, -0.5, kNSiLayers - 0.5,
        m_ptBins, m_ptMin, m_ptMax);
  }
}

void SiTpc_TrackResiduals::fill_hw(const SiHitChain& chain, double tpcPt, int tpcCharge)
{
  m_hwHasChain = 1;
  const std::size_t n = std::min({chain.layers.size(), chain.phi.size(), chain.z.size()});

  // per-layer hardware position, chain order = outermost first
  for (std::size_t i = 0; i < n; ++i)
  {
    const int l = chain.layers[i];
    if (l < 0 || l >= kNSiLayers)
    {
      continue;
    }
    m_hwUphi[l] = chain.phi[i];
    m_hwVz[l] = chain.z[i];
  }

  // step-to-step delta-delta-phi
  for (std::size_t i = 1; i < n; ++i)
  {
    const int l = chain.layers[i];
    const int p = chain.layers[i - 1];
    if (l < 0 || l >= kNSiLayers || p < 0 || p >= kNSiLayers)
    {
      continue;
    }
    m_hwPrevLayer[l] = p;
    m_hwDphiStep[l] = wrap_phi(chain.phi[i] - chain.phi[i - 1]);
    m_hwDzStep[l] = chain.z[i] - chain.z[i - 1];
    if (i >= 2)
    {
      const int pp = chain.layers[i - 2];
      const double drPrev = kSiRadius[p] - kSiRadius[pp];
      if (pp >= 0 && pp < kNSiLayers && std::abs(drPrev) > 1.0e-3 && std::isfinite(m_hwDphiStep[p]))
      {
        m_hwDphiPred[l] = m_hwDphiStep[p] * (kSiRadius[l] - kSiRadius[p]) / drPrev;
        m_hwDdphi[l] = m_hwDphiStep[l] - m_hwDphiPred[l];
      }
    }
  }

  // residuals of the chaining steps themselves
  for (const auto& st : chain.steps)
  {
    const int l = st.to_layer;
    if (l < 0 || l >= kNSiLayers)
    {
      continue;
    }
    m_hwStepKind[l] = st.kind;
    m_hwResPhi[l] = wrap_phi(st.phi - st.pred_phi);
    m_hwResZ[l] = st.z - st.pred_z;
  }

  if (!m_writeHwHists || !std::isfinite(tpcPt))
  {
    return;
  }
  const int q = charge_index(tpcCharge);
  if (q < 0)
  {
    return;
  }
  for (int l = 0; l < kNSiLayers; ++l)
  {
    if (std::isfinite(m_hwDphiStep[l]))
    {
      m_h3HwDphiStepVsLayer[q]->Fill(m_hwDphiStep[l], l, tpcPt);
    }
    if (std::isfinite(m_hwDdphi[l]))
    {
      m_h3HwDphiVsPred[l][q]->Fill(m_hwDphiPred[l], m_hwDphiStep[l], tpcPt);
      m_h3HwDdphiVsLayer[q]->Fill(m_hwDdphi[l], l, tpcPt);
    }
    if (std::isfinite(m_hwResPhi[l]))
    {
      m_h3HwResPhiVsLayer[q]->Fill(m_hwResPhi[l], l, tpcPt);
    }
    if (std::isfinite(m_hwResZ[l]))
    {
      m_h3HwResZVsLayer[q]->Fill(m_hwResZ[l], l, tpcPt);
    }
  }
}

void SiTpc_TrackResiduals::book_ddphi_hists()
{
  // booked while m_outfile is the current directory, so the file owns them
  m_outfile->cd();
  for (int q = 0; q < 2; ++q)
  {
    for (int l = 0; l < kNSiLayers; ++l)
    {
      m_h3DphiVsPrev[l][q] = new TH3F(
          std::format("h3_dphi_vs_prev_l{}_{}", l, kChargeName[q]).c_str(),
          std::format("layer {} {}: #Delta#phi vs #Delta#phi of previous (outer) layer;"
                      "#Delta#phi_{{prev}} [rad];#Delta#phi_{{L{}}} [rad];p_{{T}} [GeV]",
                      l, kChargeTitle[q], l)
              .c_str(),
          m_dphiBins, -m_dphiMax, m_dphiMax,
          m_dphiBins, -m_dphiMax, m_dphiMax,
          m_ptBins, m_ptMin, m_ptMax);
    }
    m_h3DdphiVsLayer[q] = new TH3F(
        std::format("h3_ddphi_vs_layer_{}", kChargeName[q]).c_str(),
        std::format("{}: #Delta#Delta#phi = #Delta#phi_{{L}} - #Delta#phi_{{prev}};"
                    "#Delta#Delta#phi [rad];layer;p_{{T}} [GeV]",
                    kChargeTitle[q])
            .c_str(),
        m_ddphiBins, -m_ddphiMax, m_ddphiMax,
        kNLayers, -0.5, kNLayers - 0.5,
        m_ptBins, m_ptMin, m_ptMax);
    m_h3DdphiVsDphiPrev[q] = new TH3F(
        std::format("h3_ddphi_vs_dphi_prev_{}", kChargeName[q]).c_str(),
        std::format("{}: #Delta#Delta#phi vs #Delta#phi_{{prev}} per Si layer;"
                    "#Delta#phi_{{prev}} [rad];#Delta#Delta#phi [rad];layer",
                    kChargeTitle[q])
            .c_str(),
        m_dphiBins, -m_dphiMax, m_dphiMax,
        m_ddphiBins, -m_ddphiMax, m_ddphiMax,
        kNSiLayers, -0.5, kNSiLayers - 0.5);
  }
}

void SiTpc_TrackResiduals::fill_ddphi_hists()
{
  // per-layer chain from the outside in, as in TpcSiliconCrossingMatcher:
  // the "previous" residual of a layer is the one of the next outer layer that has a point
  int prev = -1;
  for (int l = kNLayers - 1; l >= 0; --l)
  {
    if (!std::isfinite(m_dphiLayer[l]))
    {
      continue;
    }
    if (prev >= 0)
    {
      m_prevLayer[l] = prev;
      m_dphiPrevLayer[l] = m_dphiLayer[prev];
      m_ddphiLayer[l] = wrap_phi(m_dphiLayer[l] - m_dphiLayer[prev]);
    }
    prev = l;
  }

  if (!m_writeDdphiHists || m_fitStatus != SiTpcHelixFit::Ok || !std::isfinite(m_pt))
  {
    return;
  }
  const int q = charge_index(m_charge);
  if (q < 0)
  {
    return;
  }
  for (int l = 0; l < kNLayers; ++l)
  {
    if (m_prevLayer[l] < 0 || !std::isfinite(m_ddphiLayer[l]))
    {
      continue;
    }
    m_h3DdphiVsLayer[q]->Fill(m_ddphiLayer[l], l, m_pt);
    if (l < kNSiLayers)
    {
      m_h3DphiVsPrev[l][q]->Fill(m_dphiPrevLayer[l], m_dphiLayer[l], m_pt);
      m_h3DdphiVsDphiPrev[q]->Fill(m_dphiPrevLayer[l], m_ddphiLayer[l], l);
    }
  }
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
  m_tpcCrossing = m_inttCrossing = m_crossingDiff = kNoCrossing;
  m_tpcCrossingTier = -1;
  m_inttNHits = m_inttCrossingSpread = 0;
  m_nTracksCrossing = m_nSiChainsCrossing = 0;
  m_trackVertexZ = m_vertexZ = m_vertexDz = m_vertexScore = kNaN;
  m_vertexZMean = m_vertexZRms = m_vertexZErr = kNaN;
  m_vertexStatus = m_vertexOwnZ = m_vertexNSi = 0;
  for (int l = 0; l < kNLayers; ++l)
  {
    m_resRPhiLayer[l] = kNaN;
    m_resZLayer[l] = kNaN;
    m_dphiLayer[l] = kNaN;
    m_dphiPrevLayer[l] = kNaN;
    m_prevLayer[l] = -1;
    m_ddphiLayer[l] = kNaN;
  }
  m_hwHasChain = 0;
  for (int l = 0; l < kNSiLayers; ++l)
  {
    m_hwPrevLayer[l] = -1;
    m_hwStepKind[l] = -1;
    m_hwUphi[l] = m_hwVz[l] = kNaN;
    m_hwDphiStep[l] = m_hwDphiPred[l] = m_hwDdphi[l] = m_hwDzStep[l] = kNaN;
    m_hwResPhi[l] = m_hwResZ[l] = kNaN;
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
  m_siSeeds = findNode::getClass<SiHitSeedEvent>(topNode, m_siSeedNodeName);
  m_decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingNodeName);

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

  // ---- tracks per crossing in this event: final Si+TPC tracks and selected Si chains
  m_tracksPerCrossing.clear();
  m_siChainsPerCrossing.clear();
  m_vertices.clear();
  for (unsigned int it = 0; it < m_tracks->size(); ++it)
  {
    const SiTpc_Track* t = m_tracks->get_track(it);
    const TpcCrossingDecision* d = (t && t->isValid() && m_decisions) ? m_decisions->get_decision(t->get_tpc_assembled_track_id()) : nullptr;
    if (!d)
    {
      continue;
    }
    ++m_tracksPerCrossing[d->get_selected_crossing()];

    // collision vertex of this crossing (tier 0 / 1 decisions carry a vertex, tier 2 not)
    // collisions: tier 0 / 1 decisions with a collision vertex (untimed lone tracks of
    // TpcCrossingResolver keep their Si association but have no vertex: silicon_vertex_z NaN)
    const unsigned char tier = d->get_selected_tier();
    if (tier > 1 || !std::isfinite(d->get_silicon_vertex_z()))
    {
      continue;
    }
    VertexSummary& v = m_vertices[d->get_selected_crossing()];
    const bool ownZ = std::isfinite(d->get_selected_score());  // TPC-only matches have no score
    v.z = d->get_silicon_vertex_z();
    ++v.nTracks;
    v.nSi += ownZ;
    v.nTier0 += tier == 0;
    v.trackZ.push_back(d->get_tpc_z0() + d->get_delta_z());
    v.ownZ.push_back(ownZ ? 1 : 0);
    v.tier.push_back(tier);
    v.eta.push_back(static_cast<float>(t->get_eta()));
    v.pt.push_back(static_cast<float>(t->get_pt()));
    v.phi.push_back(static_cast<float>(t->get_phi()));
    v.charge.push_back(t->get_charge());
  }
  for (auto& [crossing, v] : m_vertices)
  {
    double sum = 0, sum2 = 0;
    for (std::size_t i = 0; i < v.trackZ.size(); ++i)
    {
      if (v.ownZ[i])
      {
        sum += v.trackZ[i];
        sum2 += v.trackZ[i] * v.trackZ[i];
      }
    }
    v.zMean = v.nSi > 0 ? sum / v.nSi : kNaN;
    v.zRms = v.nSi > 1 ? std::sqrt(std::max(0.0, (sum2 - v.nSi * v.zMean * v.zMean) / (v.nSi - 1))) : kNaN;
    v.zErr = v.nSi > 0 ? m_vertexTrackZSigma / std::sqrt(v.nSi) : kNaN;
    if (m_vertexTree)
    {
      m_vCrossing = crossing;
      m_vZ = v.z;
      m_vNTracks = v.nTracks;
      m_vNSi = v.nSi;
      m_vNTier0 = v.nTier0;
      m_vZMean = v.zMean;
      m_vZRms = v.zRms;
      m_vZErr = v.zErr;
      m_vTrackZ = v.trackZ;
      m_vOwnZ = v.ownZ;
      m_vTier = v.tier;
      m_vEta = v.eta;
      m_vPt = v.pt;
      m_vPhi = v.phi;
      m_vCharge = v.charge;
      m_vertexTree->Fill();
    }
  }
  if (m_siSeeds)
  {
    for (const auto& ch : m_siSeeds->chains)
    {
      if (ch.selected)
      {
        ++m_siChainsPerCrossing[ch.tpc_crossing];
      }
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
    double layerSumRPhi[kNLayers] = {}, layerSumZ[kNLayers] = {}, layerSumDphi[kNLayers] = {};
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
          layerSumDphi[layer] += dphi;
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
        m_dphiLayer[l] = layerSumDphi[l] / layerN[l];
      }
    }
    fill_ddphi_hists();

    // hardware chain of this track
    if (m_siSeeds && m_siChainId >= 0)
    {
      const SiHitChain* chain = nullptr;
      if (static_cast<std::size_t>(m_siChainId) < m_siSeeds->chains.size() &&
          m_siSeeds->chains[m_siChainId].id == m_siChainId)
      {
        chain = &m_siSeeds->chains[m_siChainId];
      }
      else
      {
        for (const auto& c : m_siSeeds->chains)
        {
          if (c.id == m_siChainId)
          {
            chain = &c;
            break;
          }
        }
      }
      if (chain)
      {
        fill_hw(*chain, m_tpcPt, m_tpcCharge);
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
    // ---- crossing: TPC (resolved / decided) and INTT (time buckets of the chain's INTT hits)
    if (const TpcCrossingDecision* d = m_decisions ? m_decisions->get_decision(m_tpcAssembledId) : nullptr)
    {
      m_tpcCrossing = d->get_selected_crossing();
      const unsigned char tier = d->get_selected_tier();
      m_tpcCrossingTier = tier == std::numeric_limits<unsigned char>::max() ? -1 : static_cast<int>(tier);
      const auto itr = m_tracksPerCrossing.find(m_tpcCrossing);
      m_nTracksCrossing = itr != m_tracksPerCrossing.end() ? itr->second : 0;
      const auto isc = m_siChainsPerCrossing.find(m_tpcCrossing);
      m_nSiChainsCrossing = isc != m_siChainsPerCrossing.end() ? isc->second : 0;

      m_vertexStatus = d->get_status();
      if (tier <= 1)
      {
        // own Si z at the beam, also for tracks without a collision vertex
        m_trackVertexZ = d->get_tpc_z0() + d->get_delta_z();
        m_vertexOwnZ = std::isfinite(d->get_selected_score()) ? 1 : 0;
      }
      const auto iv = std::isfinite(d->get_silicon_vertex_z()) ? m_vertices.find(m_tpcCrossing) : m_vertices.end();
      if (iv != m_vertices.end())
      {
        m_vertexZ = d->get_silicon_vertex_z();
        m_vertexDz = d->get_delta_z();
        m_vertexScore = d->get_selected_score();
        m_vertexNSi = iv->second.nSi;
        m_vertexZMean = iv->second.zMean;
        m_vertexZRms = iv->second.zRms;
        m_vertexZErr = iv->second.zErr;
      }
    }
    if (m_siSeeds && m_siChainId >= 0)
    {
      for (const auto& c : m_siSeeds->chains)
      {
        if (c.id != m_siChainId)
        {
          continue;
        }
        std::map<int, int> buckets;  // time bucket -> hits
        for (const int hid : c.hit_ids)
        {
          if (hid >= 0 && hid < static_cast<int>(m_siSeeds->hits.size()) && m_siSeeds->hits[hid].layer >= 3)
          {
            ++buckets[InttDefs::getTimeBucketId(m_siSeeds->hits[hid].hitsetkey)];
            ++m_inttNHits;
          }
        }
        if (!buckets.empty())
        {
          m_inttCrossing = std::max_element(buckets.begin(), buckets.end(), [](const auto& a, const auto& b)
                                            { return a.second < b.second; })
                               ->first;
          m_inttCrossingSpread = buckets.rbegin()->first - buckets.begin()->first;
        }
        break;
      }
    }
    if (m_tpcCrossing != kNoCrossing && m_inttCrossing != kNoCrossing)
    {
      m_crossingDiff = m_inttCrossing - m_tpcCrossing;
    }

    m_tree->Fill();
  }

  // keep the trees readable if the job does not reach End (crash, kill, time limit)
  m_tree->AutoSave("SaveSelf");
  if (m_matchTree)
  {
    m_matchTree->AutoSave("SaveSelf");
  }
  if (m_vertexTree)
  {
    m_vertexTree->AutoSave("SaveSelf");
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_evt << " Si+TPC tracks=" << m_tracks->size()
              << " residuals=" << nfilled << " candidates=" << m_tracks->size_candidates() << " vertices=" << m_vertices.size()
              << std::endl;
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
    if (m_vertexTree)
    {
      m_vertexTree->Write();
    }
    if (m_writeDdphiHists)
    {
      for (int q = 0; q < 2; ++q)
      {
        for (int l = 0; l < kNSiLayers; ++l)
        {
          if (m_h3DphiVsPrev[l][q])
          {
            m_h3DphiVsPrev[l][q]->Write();
          }
        }
        if (m_h3DdphiVsLayer[q])
        {
          m_h3DdphiVsLayer[q]->Write();
        }
        if (m_h3DdphiVsDphiPrev[q])
        {
          m_h3DdphiVsDphiPrev[q]->Write();
        }
      }
    }
    if (m_writeHwHists)
    {
      for (int q = 0; q < 2; ++q)
      {
        for (TH3F* h : {m_h3HwDphiStepVsLayer[q], m_h3HwDdphiVsLayer[q], m_h3HwResPhiVsLayer[q], m_h3HwResZVsLayer[q]})
        {
          if (h)
          {
            h->Write();
          }
        }
        for (int l = 0; l < kNSiLayers; ++l)
        {
          if (m_h3HwDphiVsPred[l][q])
          {
            m_h3HwDphiVsPred[l][q]->Write();
          }
        }
      }
    }
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
    m_tree = nullptr;
    m_matchTree = nullptr;
    m_vertexTree = nullptr;
    // histograms were owned by the file and deleted by Close()
    for (int q = 0; q < 2; ++q)
    {
      for (int l = 0; l < kNSiLayers; ++l)
      {
        m_h3DphiVsPrev[l][q] = nullptr;
      }
      m_h3DdphiVsLayer[q] = nullptr;
      m_h3DdphiVsDphiPrev[q] = nullptr;
      m_h3HwDphiStepVsLayer[q] = m_h3HwDdphiVsLayer[q] = m_h3HwResPhiVsLayer[q] = m_h3HwResZVsLayer[q] = nullptr;
      for (int l = 0; l < kNSiLayers; ++l)
      {
        m_h3HwDphiVsPred[l][q] = nullptr;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}