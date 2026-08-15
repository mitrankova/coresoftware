#include "PHGarfieldRawHitsQA.h"

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <tpctrackreco/IdealPadMap.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TH2F.h>

#include <array>
#include <cassert>
#include <cmath>
#include <format>
#include <iostream>
#include <string>

namespace
{
  int layer_to_module(const unsigned int layer)
  {
    if (layer < 7 || layer > 54)
    {
      return -1;
    }
    const int module = static_cast<int>((layer - 7) / 16);
    return (module >= 0 && module < 3) ? module : -1;
  }

  bool is_good_number(const double value)
  {
    return std::isfinite(value) && std::fabs(value) < 1.0e30;
  }
}  // namespace

PHGarfieldRawHitsQA::PHGarfieldRawHitsQA(const std::string& name)
  : SubsysReco(name)
  , m_idealPadMap(new IdealPadMap())
{
}

PHGarfieldRawHitsQA::~PHGarfieldRawHitsQA()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
}

int PHGarfieldRawHitsQA::InitRun(PHCompositeNode* /*topNode*/)
{
  if (!m_idealPadMap)
  {
    m_idealPadMap = new IdealPadMap();
  }

  if (!m_idealPadMap->is_loaded() && m_idealPadMap->load_from_cdb(Verbosity()) != 0)
  {
    std::cout << PHWHERE << " " << Name() << " failed to load IdealPadMap from CDB" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  createHistos();

  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  bool missing_histogram = false;
  for (int side = 0; side < 2; ++side)
  {
    for (int module = 0; module < 3; ++module)
    {
      for (int sector = 0; sector < 12; ++sector)
      {
        m_hAdcPadVsLayer[side][module][sector] = dynamic_cast<TH2*>(hm->getHisto(
            getHistoName(std::format("h_adc_pad_vs_layer_side{}_sec{:02}_R{}",
                                     side, sector, module + 1))));
        missing_histogram = missing_histogram || !m_hAdcPadVsLayer[side][module][sector];
      }
    }

    m_hAdcPhiVsRadius[side] = dynamic_cast<TH2*>(hm->getHisto(
        getHistoName(std::format("h_adc_phi_vs_radius_side{}", side))));
    m_hHighAdcHitsPhiVsRadius[side] = dynamic_cast<TH2*>(hm->getHisto(
        getHistoName(std::format("h_nhits_adc_gt{}_phi_vs_radius_side{}",
                                 m_highAdcThreshold, side))));
    missing_histogram = missing_histogram || !m_hAdcPhiVsRadius[side] || !m_hHighAdcHitsPhiVsRadius[side];
  }

  if (missing_histogram)
  {
    std::cout << PHWHERE << " " << Name() << " failed to retrieve one or more registered histograms" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHGarfieldRawHitsQA::getNodes(PHCompositeNode* topNode)
{
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, m_hitSetNodeName);
  if (!m_hits)
  {
    std::cout << PHWHERE << " Missing " << m_hitSetNodeName << std::endl;
    return false;
  }

  return true;
}

std::string PHGarfieldRawHitsQA::getHistoName(const std::string& base) const
{
  return std::format("{}_{}", base, Name());
}

double PHGarfieldRawHitsQA::idealRadius(const unsigned int layer) const
{
  if (!m_idealPadMap)
  {
    return 0.0;
  }
  return m_idealPadMap->get_radius(layer);
}

double PHGarfieldRawHitsQA::idealPhi(const int side,
                                     const unsigned int sector,
                                     const unsigned int layer,
                                     const unsigned int pad) const
{
  if (!m_idealPadMap)
  {
    return 0.0;
  }

  const unsigned int pads_per_sector = m_idealPadMap->get_pads_per_sector_for_layer(layer);
  if (pads_per_sector == 0U)
  {
    return 0.0;
  }

  const unsigned int local_pad = pad % pads_per_sector;
  return m_idealPadMap->get_phi(side, sector, layer, local_pad);
}

void PHGarfieldRawHitsQA::createHistos()
{
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  std::array<double, 49> radius_edges{};
  bool have_radius_edges = true;
  std::array<double, 48> radii{};
  for (unsigned int ilayer = 0; ilayer < radii.size(); ++ilayer)
  {
    const unsigned int layer = ilayer + 7;
    radii[ilayer] = idealRadius(layer);
    have_radius_edges = have_radius_edges && is_good_number(radii[ilayer]);
  }
  for (unsigned int i = 0; i < radius_edges.size() && have_radius_edges; ++i)
  {
    if (i == 0)
    {
      radius_edges[i] = radii[i] - 0.5 * (radii[i + 1] - radii[i]);
    }
    else if (i + 1 == radius_edges.size())
    {
      radius_edges[i] = radii[i - 1] + 0.5 * (radii[i - 1] - radii[i - 2]);
    }
    else
    {
      radius_edges[i] = 0.5 * (radii[i - 1] + radii[i]);
    }
    have_radius_edges = have_radius_edges && is_good_number(radius_edges[i]);
  }
  if (!have_radius_edges)
  {
    for (unsigned int i = 0; i < radius_edges.size(); ++i)
    {
      radius_edges[i] = 20.0 + 1.5 * static_cast<double>(i);
    }
  }

  for (int side = 0; side < 2; ++side)
  {
    for (int module = 0; module < 3; ++module)
    {
      const int layer_min = 7 + 16 * module;
      const int layer_max = layer_min + 15;
      const unsigned int pads_per_sector = m_idealPadMap ? m_idealPadMap->get_pads_per_sector(static_cast<unsigned int>(module)) : 0U;
      const unsigned int n_pads = pads_per_sector > 0U ? pads_per_sector : 1U;

      for (int sector = 0; sector < 12; ++sector)
      {
        const int pad_min = static_cast<int>(n_pads) * sector;
        const int pad_max = pad_min + static_cast<int>(n_pads) - 1;
        const std::string name = getHistoName(std::format("h_adc_pad_vs_layer_side{}_sec{:02}_R{}",
                                                          side, sector, module + 1));
        auto* hist = new TH2F(name.c_str(),
                              std::format("TPC raw hit ADC side {} sector {} R{};Pad;Layer;Sum ADC",
                                          side, sector, module + 1)
                                  .c_str(),
                              static_cast<int>(n_pads),
                              static_cast<double>(pad_min) - 0.5,
                              static_cast<double>(pad_max) + 0.5,
                              16,
                              static_cast<double>(layer_min) - 0.5,
                              static_cast<double>(layer_max) + 0.5);
        hist->SetStats(false);
        hm->registerHisto(hist);
      }
    }

    const std::string adc_name = getHistoName(std::format("h_adc_phi_vs_radius_side{}", side));
    auto* adc_hist = new TH2F(adc_name.c_str(),
                              std::format("TPC raw hit ADC side {};#phi;R [cm];Sum ADC", side).c_str(),
                              720, -std::acos(-1.0), std::acos(-1.0),
                              48, radius_edges.data());
    adc_hist->SetStats(false);
    hm->registerHisto(adc_hist);

    const std::string nhits_name = getHistoName(std::format("h_nhits_adc_gt{}_phi_vs_radius_side{}",
                                                            m_highAdcThreshold, side));
    auto* nhits_hist = new TH2F(nhits_name.c_str(),
                                std::format("TPC raw hits ADC > {} side {};#phi;R [cm];Hits",
                                            m_highAdcThreshold, side)
                                    .c_str(),
                                720, -std::acos(-1.0), std::acos(-1.0),
                                48, radius_edges.data());
    nhits_hist->SetStats(false);
    hm->registerHisto(nhits_hist);
  }
}

int PHGarfieldRawHitsQA::process_event(PHCompositeNode* topNode)
{
  if (!m_hits && !getNodes(topNode))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  TrkrHitSetContainer::ConstRange hitset_range = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hsiter = hitset_range.first;
       hsiter != hitset_range.second; ++hsiter)
  {
    TrkrHitSet* hitset = hsiter->second;
    if (!hitset)
    {
      continue;
    }

    const TrkrDefs::hitsetkey hsk = hsiter->first;
    const unsigned int layer = TrkrDefs::getLayer(hsk);
    const int module = layer_to_module(layer);
    if (module < 0)
    {
      continue;
    }

    const int side = static_cast<int>(TpcDefs::getSide(hsk));
    if (side < 0 || side > 1)
    {
      continue;
    }

    const unsigned int sector = TpcDefs::getSectorId(hsk);
    if (sector >= 12U)
    {
      continue;
    }

    const double radius = idealRadius(layer);
    if (!is_good_number(radius))
    {
      continue;
    }

    TrkrHitSet::ConstRange hit_range = hitset->getHits();
    for (TrkrHitSet::ConstIterator hiter = hit_range.first;
         hiter != hit_range.second; ++hiter)
    {
      const TrkrDefs::hitkey hk = hiter->first;
      const TrkrHit* hit = hiter->second;
      if (!hit)
      {
        continue;
      }

      const unsigned int pad = TpcDefs::getPad(hk);
      const unsigned int adc = hit->getAdc();
      const double phi = idealPhi(side, sector, layer, pad);
      if (!is_good_number(phi))
      {
        continue;
      }

      m_hAdcPadVsLayer[side][module][sector]->Fill(static_cast<double>(pad),
                                                   static_cast<double>(layer),
                                                   static_cast<double>(adc));
      m_hAdcPhiVsRadius[side]->Fill(phi, radius, static_cast<double>(adc));
      if (adc > m_highAdcThreshold)
      {
        m_hHighAdcHitsPhiVsRadius[side]->Fill(phi, radius);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
