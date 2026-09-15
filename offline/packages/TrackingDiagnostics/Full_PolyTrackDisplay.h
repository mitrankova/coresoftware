#ifndef TRACKINGDIAGNOSTICS_FULLPOLYTRACKDISPLAY_H
#define TRACKINGDIAGNOSTICS_FULLPOLYTRACKDISPLAY_H

#include <fun4all/SubsysReco.h>

#include <string>

class ActsGeometry;
class Full_PolyTrackContainer;
class PHCompositeNode;
class PHField;
class TFile;
class Tpc_PolyClusterContainer;
class TrkrClusterContainer;

class Full_PolyTrackDisplay : public SubsysReco
{
 public:
  Full_PolyTrackDisplay(const std::string& name = "Full_PolyTrackDisplay",
                        const std::string& outfilename = "full_polytrack_display.root",
                        const std::string& fullTrackNodeName = "FULL_POLYTRACKS",
                        unsigned int maxEventDisplays = 5);
  ~Full_PolyTrackDisplay() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setFullTrackNodeName(const std::string& value) { m_fullTrackNodeName = value; }
  void setTpcClusterNodeName(const std::string& value) { m_tpcClusterNodeName = value; }
  void setTrkrClusterNodeName(const std::string& value) { m_trkrClusterNodeName = value; }
  void setZRange(double minimum, double maximum) { m_zmin = minimum; m_zmax = maximum; }
  void setXYRange(double maximum) { m_xymax = maximum; }
  void setMagneticFieldTesla(double value) { m_magneticFieldTesla = value; }
  void setUseStraightLineTracks(bool value) { m_useStraightLineTracks = value; }
  void setMinTrackPt(double value) { m_minTrackPt = value; }
  void setMinMvtxHits(unsigned int value) { m_minMvtxHits = value; }
  void setMinInttHits(unsigned int value) { m_minInttHits = value; }
  void setDrawMeasurements(bool value) { m_drawMeasurements = value; }
  void setDrawFittedTrajectory(bool value) { m_drawFittedTrajectory = value; }

 private:
  bool getNodes(PHCompositeNode*);

  std::string m_outfilename;
  std::string m_fullTrackNodeName;
  std::string m_tpcClusterNodeName{"TPC_POLYCLUSTERS_CROSSING_CORRECTED"};
  std::string m_trkrClusterNodeName{"TRKR_CLUSTER"};
  unsigned int m_maxEventDisplays{5};
  unsigned int m_event{0};
  unsigned int m_eventsSaved{0};
  unsigned int m_minMvtxHits{0};
  unsigned int m_minInttHits{0};
  double m_zmin{-102.0};
  double m_zmax{102.0};
  double m_xymax{85.0};
  double m_magneticFieldTesla{1.4};
  double m_minTrackPt{0.1};
  bool m_useStraightLineTracks{false};
  bool m_drawMeasurements{true};
  bool m_drawFittedTrajectory{true};

  TFile* m_outfile{nullptr};
  Full_PolyTrackContainer* m_fullTracks{nullptr};
  Tpc_PolyClusterContainer* m_tpcClusters{nullptr};
  TrkrClusterContainer* m_trkrClusters{nullptr};
  ActsGeometry* m_actsGeometry{nullptr};
  const PHField* m_field{nullptr};
};

#endif
