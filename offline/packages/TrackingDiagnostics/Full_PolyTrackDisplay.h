#ifndef FULL_POLYTRACKDISPLAY_H
#define FULL_POLYTRACKDISPLAY_H

#include <fun4all/SubsysReco.h>

#include <string>

class ActsGeometry;
class Full_PolyTrackContainer;
class PHCompositeNode;
class TFile;
class TpcCrossingDecisionContainer;
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

  void setFull_PolyTrackNodeName(const std::string& n) { m_fullTrackNodeName = n; }
  void setTpc_PolyClusterNodeName(const std::string& n) { m_tpcClusterNodeName = n; }
  void setTrkrClusterNodeName(const std::string& n) { m_trkrClusterNodeName = n; }
  void setCrossingDecisionNodeName(const std::string& n) { m_crossingDecisionNodeName = n; }
  void setZRange(double zmin, double zmax)
  {
    m_zmin = zmin;
    m_zmax = zmax;
  }
  void setXYRange(double xymax) { m_xymax = xymax; }
  void setDrawTrackLines(bool v) { m_drawTrackLines = v; }
  void setMagneticFieldTesla(double b) { m_magneticFieldTesla = b; }
  void setUseStraightLineTracks(bool v) { m_useStraightLineTracks = v; }
  void setDrawTpcOnlyFullPolyTracks(bool v) { m_drawTpcOnlyFullPolyTracks = v; }
  void setSkipTpcOnlyFullPolyTracks(bool v) { m_drawTpcOnlyFullPolyTracks = !v; }
  void setDrawUnusedSiliconSeeds(bool v) { m_drawUnusedSiliconSeeds = v; }
  void setUnusedSiliconSeedMarkerSize(double v) { m_unusedSiliconSeedMarkerSize = v; }
  void setMinMvtxHits(unsigned int v) { m_minMvtxHits = v; }
  void setMinInttHits(unsigned int v) { m_minInttHits = v; }
  void setMinSiliconHits(unsigned int mvtx, unsigned int intt)
  {
    m_minMvtxHits = mvtx;
    m_minInttHits = intt;
  }

 private:
  bool get_nodes(PHCompositeNode* topNode);

  std::string m_outfilename;
  std::string m_fullTrackNodeName;
  std::string m_tpcClusterNodeName;
  std::string m_trkrClusterNodeName;
  std::string m_crossingDecisionNodeName;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt;
  unsigned int m_eventsSaved;
  unsigned int m_minMvtxHits;
  unsigned int m_minInttHits;
  double m_zmin;
  double m_zmax;
  double m_xymax;
  double m_magneticFieldTesla;
  double m_unusedSiliconSeedMarkerSize;
  bool m_drawTrackLines;
  bool m_useStraightLineTracks;
  bool m_drawTpcOnlyFullPolyTracks;
  bool m_drawUnusedSiliconSeeds;

  TFile* m_outfile;
  Full_PolyTrackContainer* m_fullTracks;
  Tpc_PolyClusterContainer* m_tpcClusters;
  TrkrClusterContainer* m_trkrClusters;
  ActsGeometry* m_actsGeometry;
  TpcCrossingDecisionContainer* m_crossingDecisions;
};

#endif
