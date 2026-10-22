#include "SiTpc_VertexQA.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/Tpc_AssembledTrack.h>
#include <tpctrackreco/Tpc_AssembledTrackContainer.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <map>
#include <set>

namespace
{
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
  constexpr const char* kTag[3] = {"si", "tpc", "sitpc"};
  constexpr const char* kCmpName[5] = {"Si - SiTPC", "TPC - SiTPC", "resolver - SiTPC", "Si - resolver", "TPC - Si"};

  double median(std::vector<double> v)
  {
    if (v.empty())
    {
      return kNaN;
    }
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(v.size() / 2), v.end());
    return v[v.size() / 2];
  }
  double robustSigma(const std::vector<double>& v)
  {
    const double m = median(v);
    std::vector<double> d;
    d.reserve(v.size());
    for (const double x : v)
    {
      d.push_back(std::fabs(x - m));
    }
    return 1.4826 * median(d);
  }
}  // namespace

SiTpc_VertexQA::SiTpc_VertexQA(const std::string& name, const std::string& outfile)
  : SubsysReco(name)
  , m_outfileName(outfile)
{
}

SiTpc_VertexQA::~SiTpc_VertexQA()
{
  delete m_file;
}

int SiTpc_VertexQA::Init(PHCompositeNode* /*unused*/)
{
  m_file = new TFile(m_outfileName.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open " << m_outfileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_tree = new TTree("crossings", "collision vertices per crossing: resolver, Si only, TPC only, Si+TPC");
  m_tree->Branch("event", &m_event, "event/i");
  m_tree->Branch("crossing", &m_crossing, "crossing/I");
  m_tree->Branch("res_n", &m_resN, "res_n/I");
  m_tree->Branch("res_z", &m_resZ, "res_z/D");
  for (int k = 0; k < kN; ++k)
  {
    const std::string t = kTag[k];
    m_tree->Branch((t + "_nvtx").c_str(), &m_nvtx[k], (t + "_nvtx/I").c_str());
    m_tree->Branch((t + "_n").c_str(), &m_n[k], (t + "_n/I").c_str());
    m_tree->Branch((t + "_x").c_str(), &m_x[k], (t + "_x/D").c_str());
    m_tree->Branch((t + "_y").c_str(), &m_y[k], (t + "_y/D").c_str());
    m_tree->Branch((t + "_z").c_str(), &m_z[k], (t + "_z/D").c_str());
    m_tree->Branch((t + "_ex").c_str(), &m_ex[k], (t + "_ex/D").c_str());
    m_tree->Branch((t + "_ey").c_str(), &m_ey[k], (t + "_ey/D").c_str());
    m_tree->Branch((t + "_ez").c_str(), &m_ez[k], (t + "_ez/D").c_str());
    m_tree->Branch((t + "_chi2ndf").c_str(), &m_chi2ndf[k], (t + "_chi2ndf/D").c_str());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpc_VertexQA::process_event(PHCompositeNode* topNode)
{
  ++m_event;

  // resolver collisions: crossing -> (sum z, n)
  std::map<int, std::pair<double, int>> resolver;
  // (decisions are keyed by the assembled TPC track id: loop over the assembled tracks)
  auto* decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  auto* assembled = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_assembledNodeName);
  if (decisions && assembled)
  {
    for (unsigned int i = 0; i < assembled->size(); ++i)
    {
      const Tpc_AssembledTrack* t = assembled->get_track(i);
      const TpcCrossingDecision* d = t ? decisions->get_decision(t->get_track_id()) : nullptr;
      if (!d || d->get_selected_tier() > 1 || !std::isfinite(d->get_silicon_vertex_z()))
      {
        continue;
      }
      auto& r = resolver[d->get_selected_crossing()];
      r.first += d->get_silicon_vertex_z();
      ++r.second;
    }
  }

  // vertex maps: crossing -> (best vertex, number of vertices)
  std::array<std::map<int, std::pair<const SvtxVertex*, int>>, kN> best;
  std::set<int> crossings;
  for (const auto& [c, r] : resolver)
  {
    crossings.insert(c);
  }
  for (int k = 0; k < kN; ++k)
  {
    auto* map = findNode::getClass<SvtxVertexMap>(topNode, m_mapNames[k]);
    if (!map)
    {
      continue;
    }
    for (const auto& [key, v] : *map)
    {
      if (!v)
      {
        continue;
      }
      const int c = v->get_beam_crossing();
      auto& b = best[k][c];
      ++b.second;
      if (!b.first || v->size_tracks() > b.first->size_tracks())
      {
        b.first = v;
      }
      crossings.insert(c);
    }
  }

  for (const int c : crossings)
  {
    ++m_nCrossings;
    m_crossing = c;
    const auto ir = resolver.find(c);
    m_resN = ir != resolver.end() ? ir->second.second : 0;
    m_resZ = ir != resolver.end() ? ir->second.first / ir->second.second : kNaN;
    for (int k = 0; k < kN; ++k)
    {
      const auto ib = best[k].find(c);
      const SvtxVertex* v = ib != best[k].end() ? ib->second.first : nullptr;
      m_nvtx[k] = ib != best[k].end() ? ib->second.second : 0;
      m_n[k] = v ? static_cast<int>(v->size_tracks()) : 0;
      m_x[k] = v ? v->get_x() : kNaN;
      m_y[k] = v ? v->get_y() : kNaN;
      m_z[k] = v ? v->get_z() : kNaN;
      m_ex[k] = v ? std::sqrt(std::max(0.0F, v->get_error(0, 0))) : kNaN;
      m_ey[k] = v ? std::sqrt(std::max(0.0F, v->get_error(1, 1))) : kNaN;
      m_ez[k] = v ? std::sqrt(std::max(0.0F, v->get_error(2, 2))) : kNaN;
      m_chi2ndf[k] = (v && v->get_ndof() > 0) ? v->get_chisq() / v->get_ndof() : kNaN;
      m_nCrossingsWithVertex[k] += v ? 1 : 0;
    }
    m_tree->Fill();

    // comparisons (resolver error taken as 0.22 cm / sqrt(n))
    const double resErr = m_resN > 0 ? 0.22 / std::sqrt(m_resN) : kNaN;
    auto cmp = [&](int i, double z1, double e1, double z2, double e2)
    {
      if (std::isfinite(z1) && std::isfinite(z2))
      {
        m_dz[i].push_back(z1 - z2);
        const double e = std::hypot(e1, e2);
        if (std::isfinite(e) && e > 0)
        {
          m_pull[i].push_back((z1 - z2) / e);
        }
      }
    };
    cmp(0, m_z[0], m_ez[0], m_z[2], m_ez[2]);
    cmp(1, m_z[1], m_ez[1], m_z[2], m_ez[2]);
    cmp(2, m_resZ, resErr, m_z[2], m_ez[2]);
    cmp(3, m_z[0], m_ez[0], m_resZ, resErr);
    cmp(4, m_z[1], m_ez[1], m_z[0], m_ez[0]);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpc_VertexQA::End(PHCompositeNode* /*unused*/)
{
  std::cout << Name() << "::End - crossings " << m_nCrossings << ", with a vertex: Si " << m_nCrossingsWithVertex[0] << ", TPC "
            << m_nCrossingsWithVertex[1] << ", Si+TPC " << m_nCrossingsWithVertex[2] << std::endl;
  for (int i = 0; i < kCmp; ++i)
  {
    std::cout << Name() << "::End - z " << kCmpName[i] << ": " << m_dz[i].size() << " crossings, median "
              << median(m_dz[i]) << " cm, robust sigma " << robustSigma(m_dz[i]) << " cm; pull robust sigma "
              << robustSigma(m_pull[i]) << std::endl;
  }
  if (m_file)
  {
    m_file->cd();
    m_tree->Write();
    m_file->Close();
    delete m_file;
    m_file = nullptr;
    m_tree = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
