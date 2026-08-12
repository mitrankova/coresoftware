// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHGARFIELDCALIBRATION_PHGARFIELDPOLYRESIDUALQA_H
#define PHGARFIELDCALIBRATION_PHGARFIELDPOLYRESIDUALQA_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class PHCompositeNode;
class TH1;
class TH2;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrackContainer;

class PHGarfieldPolyResidualQA : public SubsysReco
{
 public:
  explicit PHGarfieldPolyResidualQA(const std::string& name = "PHGarfieldPolyResidualQA",
                                    double kEffSide0 = 1.0,
                                    double kEffSide1 = 1.0,
                                    double spaceChargeScaleSide0 = 1.0,
                                    double spaceChargeScaleSide1 = 1.0);
  ~PHGarfieldPolyResidualQA() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setClusterNodeName(const std::string& name) { m_clusterNodeName = name; }
  void setTpcPolyTrackNodeName(const std::string& name) { m_trackNodeName = name; }
  void setMagneticFieldTesla(double value) { m_magneticFieldTesla = value; }
  void setUseStraightLineTracks(bool value) { m_useStraightLineTracks = value; }

  void setKEffSide0(double value) { m_kEffSide0 = value; }
  void setKEffSide1(double value) { m_kEffSide1 = value; }
  void setSpaceChargeScaleSide0(double value) { m_spaceChargeScaleSide0 = value; }
  void setSpaceChargeScaleSide1(double value) { m_spaceChargeScaleSide1 = value; }

 private:
  void createHistos();
  bool getNodes(PHCompositeNode* topNode);
  std::string getHistoName(const std::string& base) const;
  std::string coefficientSuffix() const;

  std::string m_clusterNodeName{"TPC_POLYCLUSTERS"};
  std::string m_trackNodeName{"TPC_POLYTRACKS"};

  double m_kEffSide0{1.0};
  double m_kEffSide1{1.0};
  double m_spaceChargeScaleSide0{1.0};
  double m_spaceChargeScaleSide1{1.0};
  double m_magneticFieldTesla{1.4};
  bool m_useStraightLineTracks{false};

  Tpc_PolyClusterContainer* m_clusters{nullptr};
  Tpc_PolyTrackContainer* m_tracks{nullptr};

  std::array<std::array<std::array<TH1*, 3>, 2>, 2> m_hRPhiResidual{};
  std::array<std::array<TH1*, 3>, 2> m_hZResidual{};
  std::array<std::array<TH2*, 2>, 2> m_hRPhiVsLayer{};
  std::array<TH2*, 2> m_hZVsLayer{};
};

#endif
