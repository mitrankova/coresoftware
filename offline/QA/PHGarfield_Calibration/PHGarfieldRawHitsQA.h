// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHGARFIELDCALIBRATION_PHGARFIELDRAWHITSQA_H
#define PHGARFIELDCALIBRATION_PHGARFIELDRAWHITSQA_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class IdealPadMap;
class PHCompositeNode;
class TH1;
class TH2;
class TrkrHitSetContainer;

class PHGarfieldRawHitsQA : public SubsysReco
{
 public:
  explicit PHGarfieldRawHitsQA(const std::string& name = "PHGarfieldRawHitsQA");
  ~PHGarfieldRawHitsQA() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setHitSetNodeName(const std::string& name) { m_hitSetNodeName = name; }
  void setHighAdcThreshold(unsigned int value) { m_highAdcThreshold = value; }

 private:
  void createHistos();
  bool getNodes(PHCompositeNode* topNode);
  std::string getHistoName(const std::string& base) const;

  double idealRadius(unsigned int layer) const;
  double idealPhi(int side, unsigned int sector, unsigned int layer, unsigned int pad) const;

  std::string m_hitSetNodeName{"TRKR_HITSET"};
  unsigned int m_highAdcThreshold{100};

  TrkrHitSetContainer* m_hits{nullptr};
  IdealPadMap* m_idealPadMap{nullptr};

  TH1* m_hNEvents{nullptr};

  std::array<std::array<std::array<TH2*, 12>, 3>, 2> m_hAdcPadVsLayer{};
  std::array<TH2*, 2> m_hAdcPhiVsRadius{};
  std::array<TH2*, 2> m_hHighAdcHitsPhiVsRadius{};
};

#endif
