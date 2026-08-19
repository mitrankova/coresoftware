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

#include <TH1F.h>
#include <TH2F.h>

#include <array>
#include <cassert>
#include <cmath>
#include <format>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace
{
  // ============================================================
  // Raw-channel waveform cleaning
  // ============================================================
  constexpr unsigned int kNTimeBins = 700;
  constexpr unsigned int kMinAcceptedTimeBin = 270;
  constexpr unsigned int kMaxAcceptedTimeBin = 700;

  // Persistent channel noise: channel fires through a sizable fraction
  // of the frame (including the beginning of the frame).
  constexpr double kPersistentOccupancy = 0.35;
  constexpr double kPersistentEarlyOccupancy = 0.30;
  constexpr unsigned int kPersistentEarlyBins = 100;
  constexpr unsigned int kMinPersistentHits = 20;

  // A real pulse on top of a noisy/saturated baseline must exceed both
  // an absolute and a relative threshold.
  constexpr double kSignalAboveBaselineADC = 30.0;
  constexpr double kSignalAboveBaselineFraction = 0.50;

  // Saturated SAMPA: after the onset keep only the first N time bins,
  // then suppress the plateau until a new pulse rises above it.
  constexpr unsigned int kKeepBinsAfterFire = 5;
  constexpr double kSaturatedTailOccupancy = 0.35;
  constexpr unsigned int kMinSaturatedTailBins = 20;

  struct TimeAdc
  {
    unsigned int tbin = 0;
    double adc = 0.0;
  };

  double median(std::vector<double> values)
  {
    if (values.empty())
    {
      return 0.0;
    }
    const auto middle = values.begin() + values.size() / 2;
    std::nth_element(values.begin(), middle, values.end());
    double value = *middle;
    if (values.size() % 2 == 0)
    {
      const auto lower = std::max_element(values.begin(), middle);
      value = 0.5 * (value + *lower);
    }
    return value;
  }

  double signal_threshold(const double baseline)
  {
    return std::max(kSignalAboveBaselineADC, kSignalAboveBaselineFraction * baseline);
  }

  bool is_persistent_noise(const std::vector<TimeAdc>& waveform)
  {
    if (waveform.size() < kMinPersistentHits)
    {
      return false;
    }

    unsigned int early_hits = 0;
    for (const auto& sample : waveform)
    {
      if (sample.tbin < kPersistentEarlyBins)
      {
        ++early_hits;
      }
    }

    const double occupancy = static_cast<double>(waveform.size()) / static_cast<double>(kNTimeBins);
    const double early_occupancy = static_cast<double>(early_hits) / static_cast<double>(kPersistentEarlyBins);
    return occupancy >= kPersistentOccupancy && early_occupancy >= kPersistentEarlyOccupancy;
  }

  double persistent_baseline(const std::vector<TimeAdc>& waveform)
  {
    // Prefer the beginning of the frame, before a physics pulse is likely
    // to dominate. The median is robust against occasional large ADC pulses.
    std::vector<double> values;
    values.reserve(waveform.size());
    for (const auto& sample : waveform)
    {
      if (sample.tbin < kPersistentEarlyBins)
      {
        values.push_back(sample.adc);
      }
    }
    if (values.size() < 5)
    {
      values.clear();
      for (const auto& sample : waveform)
      {
        values.push_back(sample.adc);
      }
    }
    return median(std::move(values));
  }

  int find_saturation_onset(const std::vector<TimeAdc>& waveform)
  {
    // Find the earliest hit for which the remainder of the frame has a
    // persistent occupancy. A normal short pulse therefore does not qualify.
    if (waveform.size() < kMinSaturatedTailBins)
    {
      return -1;
    }

    for (std::size_t i = 0; i < waveform.size(); ++i)
    {
      const unsigned int onset = waveform[i].tbin;
      const unsigned int tail_start = onset + kKeepBinsAfterFire;
      if (tail_start >= kNTimeBins)
      {
        break;
      }

      const auto first_tail = std::lower_bound(
          waveform.begin(), waveform.end(), tail_start,
          [](const TimeAdc& sample, const unsigned int tbin)
          { return sample.tbin < tbin; });

      const unsigned int tail_hits = static_cast<unsigned int>(waveform.end() - first_tail);
      const unsigned int tail_bins = kNTimeBins - tail_start;
      if (tail_hits >= kMinSaturatedTailBins &&
          static_cast<double>(tail_hits) / static_cast<double>(tail_bins) >= kSaturatedTailOccupancy)
      {
        // Do not call a channel saturated if it was already persistently noisy
        // from the beginning of the frame; that is handled separately.
        const unsigned int before_bins = std::max(1U, onset);
        const unsigned int before_hits = static_cast<unsigned int>(i);
        const double before_occupancy =
            static_cast<double>(before_hits) / static_cast<double>(before_bins);
        if (before_occupancy < kPersistentEarlyOccupancy)
        {
          return static_cast<int>(onset);
        }
      }
    }
    return -1;
  }

  std::vector<TimeAdc> clean_waveform(std::vector<TimeAdc> waveform)
  {
    if (waveform.empty())
    {
      return {};
    }

    std::sort(
        waveform.begin(), waveform.end(),
        [](const TimeAdc& lhs, const TimeAdc& rhs)
        { return lhs.tbin < rhs.tbin; });

    std::vector<TimeAdc> cleaned;
    cleaned.reserve(waveform.size());

    // ----------------------------------------------------------
    // Case 1: channel fires throughout the frame.
    // Remove its robust per-channel baseline and keep only a pulse
    // that rises significantly above that baseline.
    // ----------------------------------------------------------
    if (is_persistent_noise(waveform))
    {
      const double baseline = persistent_baseline(waveform);
      const double threshold = signal_threshold(baseline);

      for (const auto& sample : waveform)
      {
        const double excess = sample.adc - baseline;
        if (excess > threshold)
        {
          cleaned.push_back({sample.tbin, excess});
        }
      }
      return cleaned;
    }

    // ----------------------------------------------------------
    // Case 2: quiet channel fires and then stays on (SAMPA saturation).
    // ----------------------------------------------------------
    const int saturation_onset = find_saturation_onset(waveform);
    if (saturation_onset < 0)
    {
      return waveform;
    }

    const unsigned int onset = static_cast<unsigned int>(saturation_onset);
    const unsigned int first_keep_end = std::min(kNTimeBins, onset + kKeepBinsAfterFire);

    // Keep the original first five bins after the onset.
    for (const auto& sample : waveform)
    {
      if (sample.tbin >= onset && sample.tbin < first_keep_end)
      {
        cleaned.push_back(sample);
      }
    }

    // Estimate the saturated plateau after the first kept bins.
    std::vector<double> plateau_values;
    for (const auto& sample : waveform)
    {
      if (sample.tbin >= first_keep_end)
      {
        plateau_values.push_back(sample.adc);
      }
    }
    double plateau = median(std::move(plateau_values));
    double threshold = signal_threshold(plateau);

    // While suppressed, a new pulse above the plateau re-opens a 5-bin
    // acceptance window. Its ADC is baseline-subtracted.
    unsigned int accept_until = first_keep_end;
    bool in_retrigger_window = false;

    for (const auto& sample : waveform)
    {
      if (sample.tbin < first_keep_end)
      {
        continue;
      }

      if (sample.tbin < accept_until && in_retrigger_window)
      {
        const double excess = sample.adc - plateau;
        if (excess > 0.0)
        {
          cleaned.push_back({sample.tbin, excess});
        }
        continue;
      }

      const double excess = sample.adc - plateau;
      if (excess > threshold)
      {
        in_retrigger_window = true;
        accept_until = std::min(kNTimeBins, sample.tbin + kKeepBinsAfterFire);
        cleaned.push_back({sample.tbin, excess});
      }
      else
      {
        in_retrigger_window = false;
      }
    }

    return cleaned;
  }

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
  m_hNEvents = dynamic_cast<TH1*>(hm->getHisto(getHistoName("h_nEvents")));
  missing_histogram = missing_histogram || !m_hNEvents;

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

  const std::string nevents_name = getHistoName("h_nEvents");
  auto* nevents_hist = new TH1F(nevents_name.c_str(), "Number of processed events;;Events", 1, 0.5, 1.5);
  nevents_hist->SetStats(false);
  hm->registerHisto(nevents_hist);

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
  m_hNEvents->Fill(1.0);

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

    // Collect the complete time waveform for every pad/channel first.
    std::map<unsigned int, std::vector<TimeAdc>> channel_waveforms;

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
      const unsigned int tbin = TpcDefs::getTBin(hk);
      if (tbin >= kNTimeBins)
      {
        continue;
      }

      channel_waveforms[pad].push_back({tbin, static_cast<double>(hit->getAdc())});
    }

    // Clean each channel waveform and only then fill the spatial QA maps.
    for (auto& [pad, waveform] : channel_waveforms)
    {
      const double phi = idealPhi(side, sector, layer, pad);
      if (!is_good_number(phi))
      {
        continue;
      }

      const auto cleaned = clean_waveform(std::move(waveform));
      for (const auto& sample : cleaned)
      {
        if (sample.tbin <= kMinAcceptedTimeBin || sample.tbin >= kMaxAcceptedTimeBin)
        {
          continue;
        }

        const double adc = sample.adc;
        if (adc <= 0.0)
        {
          continue;
        }

        m_hAdcPadVsLayer[side][module][sector]->Fill(
            static_cast<double>(pad),
            static_cast<double>(layer),
            adc);
        m_hAdcPhiVsRadius[side]->Fill(phi, radius, adc);

        if (adc > static_cast<double>(m_highAdcThreshold))
        {
          m_hHighAdcHitsPhiVsRadius[side]->Fill(phi, radius);
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
