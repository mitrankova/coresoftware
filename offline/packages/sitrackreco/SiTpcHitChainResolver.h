// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SITPCHITCHAINRESOLVER_H
#define SITRACKRECO_SITPCHITCHAINRESOLVER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class Si_ChainContainer;
class SiTpcHitMatchContainer;

// SiTpcHitChainResolver
//
// Hit-level version of TpcSiliconCrossingResolver for the TPC-driven chains of
// SiTpcHitChainMatcher.  Ranking: number of layers, then number of MVTX layers, then chi2,
// then TPC track index.  Greedy: every TPC track and every silicon HIT is used at most once.
// Sets SiTpcHitMatch::selected and writes the selected chains as Si_Chainv1 to
// SI_TPC_HITCHAINS (DST/SVTX):
//   vertex = silicon chain z0 at the beam, slope = TPC tan(lambda), phi = TPC phi at the pca,
//   hit indices = all hits of all groups.
class SiTpcHitChainResolver : public SubsysReco
{
 public:
  explicit SiTpcHitChainResolver(const std::string& name = "SiTpcHitChainResolver");
  ~SiTpcHitChainResolver() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& s) { m_inputNodeName = s; }
  void setOutputNodeName(const std::string& s) { m_outputNodeName = s; }
  // false: a candidate sharing a hit with an already selected chain is still selected (QA)
  void setRequireUniqueHits(bool v) { m_uniqueHits = v; }

 private:
  std::string m_inputNodeName = "SI_TPC_HITMATCHES";
  std::string m_outputNodeName = "SI_TPC_HITCHAINS";
  bool m_uniqueHits = true;
  Si_ChainContainer* m_chains = nullptr;
};

#endif
