#ifndef TPCTRACKRECO_TPCSILICONCROSSINGMATCHER_H
#define TPCTRACKRECO_TPCSILICONCROSSINGMATCHER_H

#include "BeamFrameTransform.h"
#include <fun4all/SubsysReco.h>
#include <array>
#include <string>

class ActsGeometry;
class PHCompositeNode;
class TpcCrossingTrajectoryContainer;
class TpcSiliconMatchCandidateContainer;
class TrkrClusterContainer;

class TpcSiliconCrossingMatcher : public SubsysReco
{
 public:
  explicit TpcSiliconCrossingMatcher(const std::string& name = "TpcSiliconCrossingMatcher");
  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  void setTrajectoryNodeName(const std::string& value) { m_trajectoryNodeName = value; }
  void setClusterNodeName(const std::string& value) { m_clusterNodeName = value; }
  void setOutputNodeName(const std::string& value) { m_outputNodeName = value; }
  void setMaxDz(float value) { m_maxDz = value; }
  void setMaxDPhi(float value) { m_maxDPhi = value; }
  void setMinMvtx(unsigned int value) { m_minMvtx = value; }
  void setMinIntt(unsigned int value) { m_minIntt = value; }
  void setTpcBeamLine(double x0, double dxdz, double y0, double dydz) { m_beamFrame.setTpcBeamLine({x0, dxdz, y0, dydz}); }
  void setMvtxBeamLine(double x0, double dxdz, double y0, double dydz) { m_beamFrame.setMvtxBeamLine({x0, dxdz, y0, dydz}); }
  void setInttBeamLine(double x0, double dxdz, double y0, double dydz) { m_beamFrame.setInttBeamLine({x0, dxdz, y0, dydz}); }
  void setDynamicPhiMean(unsigned int layer, double offset, double slope) { if (layer < 7) { m_phiOffset[layer] = offset; m_phiSlope[layer] = slope; } }
 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  std::string m_trajectoryNodeName{"TPC_CROSSING_TRAJECTORIES"};
  std::string m_clusterNodeName{"TRKR_CLUSTER"};
  std::string m_outputNodeName{"TPC_SILICON_MATCH_CANDIDATES"};
  float m_maxDz{1.0F};
  float m_maxDPhi{0.03F};
  unsigned int m_minMvtx{2};
  unsigned int m_minIntt{1};
  BeamFrameTransform m_beamFrame;
  std::array<double, 7> m_phiOffset{};
  std::array<double, 7> m_phiSlope{{1., 1., 1., 1., 1., 1., 1.}};
  TpcCrossingTrajectoryContainer* m_trajectories{nullptr};
  TrkrClusterContainer* m_clusters{nullptr};
  ActsGeometry* m_geometry{nullptr};
  TpcSiliconMatchCandidateContainer* m_candidates{nullptr};
};
#endif
