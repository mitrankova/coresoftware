#include "SingleMvtxPoolInput.h"
#include "Fun4AllStreamingInputManager.h"
#include "mvtx_pool.h"

#include <ffarawobjects/MvtxFeeIdInfov1.h>
#include <ffarawobjects/MvtxRawEvtHeaderv2.h>
#include <ffarawobjects/MvtxRawHitContainerv1.h>
#include <ffarawobjects/MvtxRawHitv1.h>
#include <fun4all/Fun4AllUtils.h>
#include "MvtxRawDefs.h"

#include <frog/FROG.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/sphenix_constants.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

#include <cassert>
#include <cmath>
#include <memory>
#include <set>

SingleMvtxPoolInput::SingleMvtxPoolInput(const std::string &name)
  : SingleStreamingInput(name)
{
  plist = new Packet *[2];
  m_rawHitContainerName = "MVTXRAWHIT";
}

SingleMvtxPoolInput::~SingleMvtxPoolInput()
{
  delete[] plist;
  for (auto &iter : poolmap)
  {
    if (Verbosity() > 2)
    {
      std::cout << "deleting mvtx pool for id " << iter.first << std::endl;
    }
    delete (iter.second);
  }
}

void SingleMvtxPoolInput::FillPool(const uint64_t minBCO)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (GetEventiterator() == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return;
    }
  }

  //  std::set<uint64_t> saved_beamclocks;
  while (GetSomeMoreEvents())
  {
    std::unique_ptr<Event> evt(GetEventiterator()->getNextEvent());
    while (!evt)
    {
      fileclose();
      if (!OpenNextFile())
      {
        AllDone(1);
        return;
      }
      evt.reset(GetEventiterator()->getNextEvent());
    }
    if (Verbosity() > 2)
    {
      std::cout << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }
    RunNumber(evt->getRunNumber());
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      continue;
    }
    int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, 2);

    if (npackets > 2)
    {
      exit(1);
    }
    for (int i = 0; i < npackets; i++)
    {
      // Ignoring packet not from MVTX detector
      if (Verbosity() > 1)
      {
        plist[i]->identify();
      }

      if (poolmap.find(plist[i]->getIdentifier()) == poolmap.end())
      {
        if (Verbosity() > 1)
        {
          std::cout << "starting new mvtx pool for packet " << plist[i]->getIdentifier() << std::endl;
        }
        poolmap[plist[i]->getIdentifier()] = new mvtx_pool();
      }
      poolmap[plist[i]->getIdentifier()]->addPacket(plist[i]);
      delete plist[i];
    }
    for (auto &iter : poolmap)
    {
      mvtx_pool *pool = iter.second;
      int num_feeId = pool->get_feeidSet_size();
      if (Verbosity() > 1)
      {
        std::cout << "Number of feeid in RCDAQ events: " << num_feeId << " for packet "
                  << iter.first << std::endl;
      }
      if (num_feeId > 0)
      {
        for (int i_fee{0}; i_fee < num_feeId; ++i_fee)
        {
          // auto feeId = pool->iValue(i_fee, "FEEID");
          auto feeId = pool->get_feeid(i_fee);
          auto link = MvtxRawDefs::decode_feeid(feeId);

          //          auto hbfSize = plist[i]->iValue(feeId, "NR_HBF");
          // auto num_strobes_old = pool->iValue(feeId, "NR_STROBES");
          // auto num_L1Trgs_old = pool->iValue(feeId, "NR_PHYS_TRG");
          auto num_strobes = pool->get_strbSet_size(feeId);
          auto num_L1Trgs = pool->get_trgSet_size(feeId);
          for (int iL1 = 0; iL1 < num_L1Trgs; ++iL1)
          {
            //  auto l1Trg_bco = pool->lValue(feeId, iL1, "L1_IR_BCO");
            auto l1Trg_bco = pool->get_L1_IR_BCO(feeId, iL1);
            //            auto l1Trg_bc  = plist[i]->iValue(feeId, iL1, "L1_IR_BC");
            m_FeeGTML1BCOMap[feeId].insert(l1Trg_bco);
            gtmL1BcoSet.emplace(l1Trg_bco);
          }
          m_FeeStrobeMap[feeId] += num_strobes;
          for (int i_strb{0}; i_strb < num_strobes; ++i_strb)
          {
            // auto strb_detField = pool->iValue(feeId, i_strb, "TRG_DET_FIELD");
            // uint64_t strb_bco = pool->lValue(feeId, i_strb, "TRG_IR_BCO");
            // auto strb_bc = pool->iValue(feeId, i_strb, "TRG_IR_BC");
            // auto num_hits = pool->iValue(feeId, i_strb, "TRG_NR_HITS");
            auto strb_detField = pool->get_TRG_DET_FIELD(feeId, i_strb);
            uint64_t strb_bco = pool->get_TRG_IR_BCO(feeId, i_strb);
            auto strb_bc = pool->get_TRG_IR_BC(feeId, i_strb);
            auto num_hits = pool->get_TRG_NR_HITS(feeId, i_strb);
            m_BclkStack.insert(strb_bco);
            m_FEEBclkMap[feeId] = strb_bco;

            if (strb_bco < minBCO - m_NegativeBco)
            {
              continue;
            }

            if (Verbosity() > 4)
            {
              std::cout << "evtno: " << EventSequence << ", Fee: " << feeId;
              std::cout << " Layer: " << link.layer << " Stave: " << link.stave;
              std::cout << " GBT: " << link.gbtid << ", bco: 0x" << std::hex << strb_bco << std::dec;
              std::cout << ", n_hits: " << num_hits << std::endl;
            }
            auto hits = pool->get_hits(feeId, i_strb);
            for (auto &&hit : hits)
            {
              auto newhit = std::make_unique<MvtxRawHitv1>();
              newhit->set_bco(strb_bco);
              newhit->set_strobe_bc(strb_bc);
              newhit->set_chip_bc(hit->bunchcounter);
              newhit->set_layer_id(link.layer);
              newhit->set_stave_id(link.stave);
              newhit->set_chip_id(
                  MvtxRawDefs::gbtChipId_to_staveChipId[link.gbtid][hit->chip_id]);
              newhit->set_row(hit->row_pos);
              newhit->set_col(hit->col_pos);
              if (StreamingInputManager())
              {
                StreamingInputManager()->AddMvtxRawHit(strb_bco, newhit.get());
              }
              m_MvtxRawHitMap[strb_bco].push_back(newhit.release());
            }
            if (StreamingInputManager())
            {
              StreamingInputManager()->AddMvtxFeeIdInfo(strb_bco, feeId, strb_detField);
            }
          }
        }
      }
    }
    // Assign L1 trg to Strobe windows data.
    for (auto &lv1Bco : gtmL1BcoSet)
    {
      auto it = m_BclkStack.lower_bound(lv1Bco);
      auto const strb_it = (it == m_BclkStack.begin()) ? (*it == lv1Bco ? it : m_BclkStack.cend()) : --it;
      if (strb_it != m_BclkStack.cend())
      {
        if (StreamingInputManager())
        {
          StreamingInputManager()->AddMvtxL1TrgBco(*strb_it, lv1Bco);
        }
      }
      else if (m_BclkStack.empty())
      {
        continue;
      }
      else
      {
        std::cout << "ERROR: lv1Bco: 0x" << std::hex << lv1Bco << std::dec
                  << " is less than minimun strobe bco 0x" << std::hex
                  << *m_BclkStack.begin() << std::dec << std::endl;
        // assert(0);
      }
    }
    gtmL1BcoSet.clear();
  }
}

void SingleMvtxPoolInput::Print(const std::string &what) const
{
  // TODO: adapt to MVTX case

  if (what == "ALL" || what == "FEEBCLK")
  {
    for (auto bcliter : m_FEEBclkMap)
    {
      std::cout << "FEE" << bcliter.first << " bclk: 0x"
                << std::hex << bcliter.second << std::dec << std::endl;
    }
  }
  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto &bcliter : m_MvtxRawHitMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (const auto &feeiter : bcliter.second)
      {
        std::cout << "fee: " << feeiter->get_stave_id()
                  << " at " << std::hex << feeiter << std::dec << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "STACK")
  {
    for (const auto &iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
  if (what == "ALL" || what == "GET_NR_STROBES")
  {
    for (const auto &iter : m_FeeStrobeMap)
    {
      std::cout << "Total number of strobes for feeid: " << iter.first << ", " << iter.second << std::endl;
    }
  }
}

void SingleMvtxPoolInput::CleanupUsedPackets(const uint64_t bclk)
{
  m_BclkStack.erase(m_BclkStack.begin(), m_BclkStack.upper_bound(bclk));
  for (auto it = m_MvtxRawHitMap.begin(); it != m_MvtxRawHitMap.end() && (it->first <= bclk); it = m_MvtxRawHitMap.erase(it))
  {
    for (const auto &rawhit : it->second)
    {
      delete rawhit;
    }
  }
  m_MvtxRawHitMap.erase(m_MvtxRawHitMap.begin(), m_MvtxRawHitMap.upper_bound(bclk));
  m_FeeStrobeMap.erase(m_FeeStrobeMap.begin(), m_FeeStrobeMap.upper_bound(bclk));
  for (auto &[feeid, gtmbcoset] : m_FeeGTML1BCOMap)
  {
    gtmbcoset.erase(gtmbcoset.begin(), gtmbcoset.upper_bound(bclk));
  }
}

bool SingleMvtxPoolInput::CheckPoolDepth(const uint64_t bclk)
{
  // if (m_FEEBclkMap.size() < 10)
  // {
  //   std::cout << "not all FEEs in map: " << m_FEEBclkMap.size() << std::endl;
  //   return true;
  // }
  for (auto iter : m_FEEBclkMap)
  {
    if (Verbosity() > 2)
    {
      std::cout << iter.first << " my bclk 0x" << std::hex << iter.second
                << " req: 0x" << bclk << std::dec << std::endl;
    }
    // equal case when we have more strobe with same bco
    // due not synchronization
    if (iter.second <= bclk)
    {
      if (Verbosity() > 1)
      {
        std::cout << "FEE " << iter.first << " beamclock 0x" << std::hex << iter.second
                  << " smaller than req bclk: 0x" << bclk << std::dec << std::endl;
      }
      return false;
    }
  }
  return true;
}

void SingleMvtxPoolInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
  //  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  return;
}

bool SingleMvtxPoolInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (m_MvtxRawHitMap.empty())
  {
    return true;
  }
  uint64_t lowest_bclk = m_MvtxRawHitMap.begin()->first;
  //  lowest_bclk += m_BcoRange;
  lowest_bclk += m_BcoRange;
  std::set<int> toerase;
  for (auto bcliter : m_FEEBclkMap)
  {
    if (bcliter.second <= lowest_bclk)
    {
      uint64_t highest_bclk = m_MvtxRawHitMap.rbegin()->first;
      if ((highest_bclk - m_MvtxRawHitMap.begin()->first) < MaxBclkDiff())
      {
        // std::cout << "FEE " << bcliter.first << " bclk: "
        // 		<< std::hex << bcliter.second << ", req: " << lowest_bclk
        // 		 << " low: 0x" <<  m_MvtxRawHitMap.begin()->first << ", high: " << highest_bclk << ", delta: " << std::dec << (highest_bclk-m_MvtxRawHitMap.begin()->first)
        // 		<< std::dec << std::endl;
        return true;
      }
      else
      {
        std::cout << PHWHERE << Name() << ": erasing FEE " << bcliter.first
                  << " with stuck bclk: " << std::hex << bcliter.second
                  << " current bco range: 0x" << m_MvtxRawHitMap.begin()->first
                  << ", to: 0x" << highest_bclk << ", delta: " << std::dec
                  << (highest_bclk - m_MvtxRawHitMap.begin()->first)
                  << std::dec << std::endl;
        toerase.insert(bcliter.first);
      }
    }
  }
  for (auto iter : toerase)
  {
    m_FEEBclkMap.erase(iter);
  }
  return false;

  //  if (CheckPoolDepth(m_MvtxRawHitMap.begin()->first))
  //  {
  //     if (m_MvtxRawHitMap.size() >= 200)
  //     {
  //       return false;
  //     }
  // //  }
  //   return true;
}

void SingleMvtxPoolInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "MVTX"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("MVTX");
    dstNode->addNode(detNode);
  }

  MvtxRawEvtHeader *mvtxEH = findNode::getClass<MvtxRawEvtHeader>(detNode, m_rawEventHeaderName);
  if (!mvtxEH)
  {
    mvtxEH = new MvtxRawEvtHeaderv2();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(mvtxEH, m_rawEventHeaderName, "PHObject");
    detNode->addNode(newNode);
  }

  MvtxRawHitContainer *mvtxhitcont = findNode::getClass<MvtxRawHitContainer>(detNode, m_rawHitContainerName);
  if (!mvtxhitcont)
  {
    mvtxhitcont = new MvtxRawHitContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(mvtxhitcont, m_rawHitContainerName, "PHObject");
    detNode->addNode(newNode);
  }
}

void SingleMvtxPoolInput::ConfigureStreamingInputManager()
{
  auto [runnumber, segment] = Fun4AllUtils::GetRunSegment(*(GetFileList().begin()));

  if (m_readStrWidthFromDB)
  {
    m_strobeWidth = MvtxRawDefs::getStrobeLength(runnumber);
    if (std::isnan(m_strobeWidth))
    {
      std::cout << PHWHERE << "WARNING: Strobe length is not defined for run " << runnumber;
      std::cout << " neither in the OCDB or DAQ DB. Exiting SingleMvtxPoolInput." << std::endl;
      // std::cout << "Defaulting to 89 mus strobe length" << std::endl;
      // m_strobeWidth = 89.;
      exit(1);
    }
  }

  if (!m_mvtx_is_standalone)
  {
    if (m_strobeWidth > 88.)
    {
      m_BcoRange = 1000;
      m_NegativeBco = 1000;
    }
    else if (m_strobeWidth > 9 && m_strobeWidth < 11)
    {
      m_BcoRange = 500;
      m_NegativeBco = 500;
    }
    else if (m_strobeWidth < 1)  // triggered mode
    {
      m_BcoRange = 3;
      m_NegativeBco = 0;
      if (StreamingInputManager())
      {
        StreamingInputManager()->runMvtxTriggered(true);
      }
    }
    else  // catchall for anyting else to set to a range based on the rhic clock
    {
      m_BcoRange = std::ceil(m_strobeWidth * 1000. / sphenix_constants::time_between_crossings);
      m_NegativeBco = std::ceil(m_strobeWidth * 1000. / sphenix_constants::time_between_crossings);
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << "Mvtx strobe length " << m_strobeWidth << std::endl;
    std::cout << "Mvtx BCO range and negative bco range set based on strobe length " << m_BcoRange << ", " << m_NegativeBco << std::endl;
  }

  if (StreamingInputManager())
  {
    StreamingInputManager()->SetMvtxBcoRange(m_BcoRange);
    StreamingInputManager()->SetMvtxNegativeBco(m_NegativeBco);
  }
  return;
}
