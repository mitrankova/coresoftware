/// ===========================================================================
/*! \file   TrksInJetQAInJetFiller.h
 *  \author Derek Anderson
 *  \date   04.03.2024
 *
 *  A submodule for the TrksInJetsQA F4A module to produce
 *  QA histograms for tracks and more in jets
 */
/// ===========================================================================

#ifndef TRKSINJETQAINJETFILLER_H
#define TRKSINJETQAINJETFILLER_H

// module utilities
#include "TrksInJetQABaseFiller.h"
#include "TrksInJetQADefs.h"

// g4eval includes
#include <g4eval/ClusKeyIter.h>

// jet includes
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

// particle flow includes
#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>

// phool includes
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// tracking includes
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

// c+ utilities
#include <cassert>
#include <vector>

// ============================================================================
//! In-jet histogram filler for TrksInJetQA module
// ============================================================================
/*! This histogram filler defines how to fill histograms
 *  for in-jet populations.
 */
class TrksInJetQAInJetFiller : public TrksInJetQABaseFiller
{
 public:
  ///! enumerates additional nodes to grab
  enum Node
  {
    Flow
  };

  // ctor/dtor
  using TrksInJetQABaseFiller::TrksInJetQABaseFiller;
  ~TrksInJetQAInJetFiller() override = default;

  // inherited public methods
  void Fill(PHCompositeNode* topNode) override;

 private:
  // private methods
  void GetNode(const int node, PHCompositeNode* topNode);
  void FillJetAndTrackQAHists(PHCompositeNode* topNode);
  void FillClustAndHitQAHists(SvtxTrack* track);
  void GetCstTracks(Jet* jet, PHCompositeNode* topNode);
  void GetNonCstTracks(Jet* jet);
  static bool IsCstNotRelevant(const uint32_t type);
  bool IsTrkInList(const uint32_t id);
  static double GetTrackJetDist(SvtxTrack* track, Jet* jet);
  TrksInJetQADefs::PFObject* GetPFObject(const uint32_t id, PHCompositeNode* topNode);
  static SvtxTrack* GetTrkFromPFO(TrksInJetQADefs::PFObject* pfo);

  ///! node of particle flow elements
  TrksInJetQADefs::PFObjectStore* m_flowStore {nullptr};

  // for tracks in jet
  std::vector<SvtxTrack*> m_trksInJet;

};  // end TrksInJetQAInJetFiller

#endif

// end ========================================================================
