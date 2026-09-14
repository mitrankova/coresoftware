#ifndef TPCTRACKRECO_TPCSILICONCROSSINGRESOLVER_H
#define TPCTRACKRECO_TPCSILICONCROSSINGRESOLVER_H
#include <fun4all/SubsysReco.h>
#include <string>
class PHCompositeNode;
class TpcSiliconMatchCandidateContainer;
class TpcSiliconCrossingResolver : public SubsysReco
{
 public:
  explicit TpcSiliconCrossingResolver(const std::string& name = "TpcSiliconCrossingResolver");
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  void setCandidateNodeName(const std::string& value) { m_candidateNodeName = value; }
 private:
  std::string m_candidateNodeName{"TPC_SILICON_MATCH_CANDIDATES"};
  TpcSiliconMatchCandidateContainer* m_candidates{nullptr};
};
#endif
