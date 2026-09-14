#ifndef TPCTRACKRECO_TPCCROSSINGFINDER_H
#define TPCTRACKRECO_TPCCROSSINGFINDER_H

#include "TpcCrossingDecision.h"
#include "TpcDriftPolylineLookup.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

class PHCompositeNode;
class SvtxVertexMap;
class Tpc_AssembledTrack;
class Tpc_AssembledTrackContainer;
class TpcCrossingDecisionContainerv1;
class TrkrClusterContainer;
class TrkrHitSetContainer;

class TpcCrossingFinder : public SubsysReco
{
 public:
  explicit TpcCrossingFinder(const std::string& name = "TpcCrossingFinder");
  ~TpcCrossingFinder() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setVertexMapNodeName(const std::string& n) { m_vertexMapNodeName = n; }

  void setRequireSiliconVertex(bool v) { m_requireSiliconVertex = v; }
  void setResolveAmbiguousWithoutVertex(bool v) { m_resolveAmbiguousWithoutVertex = v; }
  void setPreferTriggeredCrossing(bool v) { m_preferTriggeredCrossing = v; }
  void setTriggeredCrossing(short v) { m_triggeredCrossing = v; }
  void setTriggeredMode(bool v) { m_triggeredMode = v; }
  void setCollisionZ(double v) { m_collisionZ = v; }
  void setMaxVertexDz(double v) { m_maxVertexDz = v; }
  void setMaxTier2BeamlineZ(double v) { m_maxTier2BeamlineZ = v; }
  void setMaxCandidateVertexZ(double v) { m_maxCandidateVertexZ = v; }
  void setMinBestSecondSeparation(double v) { m_minBestSecondSeparation = v; }
  void setTpcGeometryTolerance(double radial_cm, double z_cm, double central_membrane_cm)
  {
    m_radialTolerance = radial_cm;
    m_zTolerance = z_cm;
    m_centralMembraneTolerance = central_membrane_cm;
  }
  void setTpcHalfLength(double v) { m_tpcHalfLength = v; }

 private:
  struct Point
  {
    TrkrDefs::hitsetkey hitsetkey{0};
    TrkrDefs::hitkey hitkey{0};
    unsigned int layer{0};
    unsigned int side{0};
    unsigned int pad{0};
    unsigned int tbin{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  struct SiliconVertexHypothesis
  {
    short crossing{0};
    unsigned int vertex_id{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double sigma_z{0.0};
    unsigned int ntracks{0};
  };

  struct ZFitResult
  {
    bool valid{false};
    double slope{0.0};
    double intercept{0.0};
    double chi2{0.0};
    int ndf{-1};
    double s_at_pca{0.0};
    double minimum_radius{0.0};
    double z_at_pca{0.0};
    double z_at_r0{0.0};
    std::vector<Point> points;
    std::vector<float> path_length;
  };

  struct Candidate
  {
    short crossing{0};
    bool tpc_valid{false};
    bool has_silicon_vertex{false};
    bool vertex_compatible{false};
    unsigned int silicon_vertex_id{0};
    double tpc_z0{0.0};
    double silicon_vertex_z{0.0};
    double delta_z{0.0};
    unsigned char rejection_status{0};
    unsigned char confidence_tier{std::numeric_limits<unsigned char>::max()};
    double confidence_score{std::numeric_limits<double>::quiet_NaN()};
    TpcCrossingCandidate qa;
  };

  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool make_xyz_point(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, short crossing, Point& p) const;
  bool find_time_extrema(const Tpc_AssembledTrack* track, TrkrDefs::hitsetkey& min_hsk, TrkrDefs::hitkey& min_hk, TrkrDefs::hitsetkey& max_hsk, TrkrDefs::hitkey& max_hk) const;
  std::set<short> get_available_crossings() const;
  std::set<short> get_intt_crossings() const;
  std::map<short, std::vector<SiliconVertexHypothesis>> get_vertices_by_crossing() const;
  std::vector<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>> select_representatives(const Tpc_AssembledTrack* track, TrkrDefs::hitsetkey min_hsk, TrkrDefs::hitkey min_hk, TrkrDefs::hitsetkey max_hsk, TrkrDefs::hitkey max_hk) const;
  bool estimate_tpc_z0(std::vector<Point>& points, double& z0) const;
  ZFitResult estimate_tpc_z0_diagnostics(std::vector<Point> points) const;
  bool point_in_tpc(const Point& p) const;
  bool point_in_correct_side(const Point& p) const;
  Candidate test_candidate(const Tpc_AssembledTrack* track, short crossing, TrkrDefs::hitsetkey min_hsk, TrkrDefs::hitkey min_hk, TrkrDefs::hitsetkey max_hsk, TrkrDefs::hitkey max_hk, const std::map<short, std::vector<SiliconVertexHypothesis>>& vertices_by_crossing) const;

  std::string m_inputNodeName{"TPC_ASSEMBLEDTRACKS"};
  std::string m_outputNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_vertexMapNodeName{"SvtxVertexMap"};

  Tpc_AssembledTrackContainer* m_assembledTracks{nullptr};
  TpcCrossingDecisionContainerv1* m_decisions{nullptr};
  TrkrHitSetContainer* m_hits{nullptr};
  TrkrClusterContainer* m_clusterMap{nullptr};
  SvtxVertexMap* m_vertexMap{nullptr};
  TpcDriftPolylineLookup* m_driftLookup{nullptr};

  unsigned int m_event{0};
  double m_tpcHalfLength{105.5};
  double m_radialTolerance{2.0};
  double m_zTolerance{1.0};
  double m_centralMembraneTolerance{1.0};
  double m_maxVertexDz{2.0};
  double m_minBestSecondSeparation{0.4};
  double m_collisionZ{0.0};
  double m_maxTier2BeamlineZ{40.0};
  double m_maxCandidateVertexZ{20.0};
  bool m_requireSiliconVertex{false};
  bool m_resolveAmbiguousWithoutVertex{true};
  bool m_preferTriggeredCrossing{false};
  short m_triggeredCrossing{0};
  bool m_triggeredMode{false};
};

#endif  // TPCTRACKRECO_TPCCROSSINGFINDER_H
