#include "SiTpcHitChainResolver.h"

#include "SiHitGeometryTable.h"
#include "SiTpcHitMatch.h"
#include "SiTpcHitMatchContainer.h"
#include "Si_ChainContainerv1.h"
#include "Si_Chainv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <trackbase/TrkrDefs.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>

SiTpcHitChainResolver::SiTpcHitChainResolver(const std::string& name)
  : SubsysReco(name)
{
}

int SiTpcHitChainResolver::InitRun(PHCompositeNode* topNode)
{
  m_chains = findNode::getClass<Si_ChainContainer>(topNode, m_outputNodeName);
  if (m_chains)
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
  m_chains = new Si_ChainContainerv1();
  svtxNode->addNode(new PHIODataNode<PHObject>(m_chains, m_outputNodeName, "PHObject"));
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTpcHitChainResolver::process_event(PHCompositeNode* topNode)
{
  const auto begin = std::chrono::steady_clock::now();
  if (!m_chains)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  m_chains->Reset();
  auto* matches = findNode::getClass<SiTpcHitMatchContainer>(topNode, m_inputNodeName);
  if (!matches)
  {
    if (Verbosity() > 0)
    {
      std::cout << Name() << ": missing " << m_inputNodeName << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  std::vector<SiTpcHitMatch*> ordered;
  for (unsigned int i = 0; i < matches->size(); ++i)
  {
    SiTpcHitMatch* m = matches->get_match(i);
    if (!m)
    {
      continue;
    }
    m->set_selected(false);
    if (m->isValid())
    {
      ordered.push_back(m);
    }
  }
  std::sort(ordered.begin(), ordered.end(), [](const SiTpcHitMatch* a, const SiTpcHitMatch* b)
            {
              const unsigned int la = a->get_n_mvtx_layers() + a->get_n_intt_layers();
              const unsigned int lb = b->get_n_mvtx_layers() + b->get_n_intt_layers();
              if (la != lb) return la > lb;
              if (a->get_n_mvtx_layers() != b->get_n_mvtx_layers()) return a->get_n_mvtx_layers() > b->get_n_mvtx_layers();
              if (a->get_chi2() != b->get_chi2()) return a->get_chi2() < b->get_chi2();
              return a->get_tpc_track_index() < b->get_tpc_track_index();
            });

  auto pack = [](const SiTpcHitMatch::HitIndex& h)
  { return (static_cast<uint64_t>(h.first) << 32) | static_cast<uint64_t>(h.second); };
  std::set<unsigned int> usedTracks;
  std::set<uint64_t> usedHits;
  unsigned int nConflicts = 0;
  for (SiTpcHitMatch* m : ordered)
  {
    if (usedTracks.count(m->get_tpc_track_index()))
    {
      continue;
    }
    if (m_uniqueHits)
    {
      bool conflict = false;
      for (const auto& h : m->get_hit_indices())
      {
        if (usedHits.count(pack(h)))
        {
          conflict = true;
          break;
        }
      }
      if (conflict)
      {
        ++nConflicts;
        if (Verbosity() > 1)
        {
          std::cout << Name() << " skip_shared_hit tpc_track_index=" << m->get_tpc_track_index()
                    << " chi2=" << m->get_chi2() << std::endl;
        }
        continue;
      }
    }
    m->set_selected(true);
    usedTracks.insert(m->get_tpc_track_index());
    for (const auto& h : m->get_hit_indices())
    {
      usedHits.insert(pack(h));
    }

    auto* chain = new Si_Chainv1();
    chain->set_chain_id(m->get_id());
    unsigned int nMvtx = 0, nIntt = 0;
    for (const auto& h : m->get_hit_indices())
    {
      chain->add_hit_index(h.first, h.second);
      if (SiHitGeometryTable::isMvtx(TrkrDefs::getLayer(h.first)))
      {
        ++nMvtx;
      }
      else
      {
        ++nIntt;
      }
    }
    chain->set_n_mvtx_hits(nMvtx);
    chain->set_n_intt_hits(nIntt);
    chain->set_vertex(m->get_si_z0());
    chain->set_slope(m->get_tpc_tanl());
    chain->set_phi(m->get_tpc_phi());
    m_chains->add_chain(chain);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": candidates=" << ordered.size() << " selected=" << m_chains->size()
              << " shared_hit_conflicts=" << nConflicts << " seconds="
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
