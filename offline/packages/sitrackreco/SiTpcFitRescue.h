// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCFITRESCUE_H
#define SITRACKRECO_SITPCFITRESCUE_H

// SiTpcFitRescue: keeps the TPC part of Si+TPC tracks whose ACTS fit failed.
//
// PHActsTrkFitter drops a track completely when its fit fails (PropagatorError "maximum number of steps",
// MagneticFieldError "interpolation out of bounds"): with the own Si chains this happens mostly for wrong /
// poor Si associations (large seed dca, MVTX-only chains), and the TPC track, which would fit on its own, is
// lost.  Two instances of this module around a second PHActsTrkFitter recover it as a TPC-only track:
//
//   PHActsTrkFitter                  SvtxTrackSeedContainer          -> SvtxTrackMap
//   SiTpcFitRescue(MakeSeeds)        Si+TPC seeds without a track    -> SvtxTrackSeedContainer_Rescue
//                                    (TPC seed + crossing carrier silicon seed: no clusters, TPC helix, z0 at
//                                    the seed's crossing, as SiTpcTrackSeedFiller::TpcOnlyCrossing::FromDecision)
//   PHActsTrkFitter (2nd instance)   SvtxTrackSeedContainer_Rescue   -> SvtxTrackMap_Rescue
//   SiTpcFitRescue(MergeTracks)      SvtxTrackMap_Rescue appended to SvtxTrackMap (new keys = ids)
//
// Then PHSimpleVertexFinder / TrackResiduals see every TPC track.  The rescued tracks have no MVTX / INTT
// clusters (the failed Si association is dropped); the carrier silicon seeds are appended to
// SiliconTrackSeedContainer.  A track is "fitted" if a track of SvtxTrackMap has its TPC seed.
#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class SiTpcFitRescue : public SubsysReco
{
 public:
  enum class Mode
  {
    MakeSeeds,
    MergeTracks
  };

  explicit SiTpcFitRescue(Mode mode, const std::string& name = "SiTpcFitRescue");
  ~SiTpcFitRescue() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // ---- nodes (defaults as in the macro)
  void setSvtxSeedContainerName(const std::string& s) { m_svtxSeedName = s; }
  void setSiliconSeedContainerName(const std::string& s) { m_siSeedName = s; }
  void setTpcSeedContainerName(const std::string& s) { m_tpcSeedName = s; }
  void setTrackMapName(const std::string& s) { m_trackMapName = s; }
  void setRescueSeedContainerName(const std::string& s) { m_rescueSeedName = s; }
  void setRescueTrackMapName(const std::string& s) { m_rescueTrackMapName = s; }
  // MakeSeeds: also retry failed TPC-only seeds (default false: the same fit would fail again)
  void setRetryTpcOnly(bool b) { m_retryTpcOnly = b; }

 private:
  int makeSeeds(PHCompositeNode* topNode);
  int mergeTracks(PHCompositeNode* topNode);

  Mode m_mode;
  std::string m_svtxSeedName{"SvtxTrackSeedContainer"};
  std::string m_siSeedName{"SiliconTrackSeedContainer"};
  std::string m_tpcSeedName{"TpcTrackSeedContainer"};
  std::string m_trackMapName{"SvtxTrackMap"};
  std::string m_rescueSeedName{"SvtxTrackSeedContainer_Rescue"};
  std::string m_rescueTrackMapName{"SvtxTrackMap_Rescue"};
  bool m_retryTpcOnly{false};

  unsigned long m_nEvents{0};
  unsigned long m_nSeeds{0};         // SvtxTrackSeeds seen
  unsigned long m_nFailed{0};        // ... without a fitted track
  unsigned long m_nFailedTpcOnly{0}; // ... of them TPC-only (not retried unless setRetryTpcOnly)
  unsigned long m_nRescueSeeds{0};   // rescue seeds made
  unsigned long m_nNoCrossing{0};    // failed seeds without crossing / side: not retried
  unsigned long m_nMerged{0};        // rescued tracks appended to the track map
};

#endif
