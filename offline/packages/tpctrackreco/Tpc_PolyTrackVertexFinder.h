// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCPOLYTRACKVERTEXFINDER_H
#define TPCTRACKRECO_TPCPOLYTRACKVERTEXFINDER_H

/*!
 * \file  Tpc_PolyTrackVertexFinder.h
 * \brief Primary (collision) vertex finder using Tpc_PolyTrackReco tracks.
 *
 * Algorithm (pair finding follows PHSimpleVertexFinder):
 *  1. Select good Tpc_PolyTracks and build a helix (or line) model for each.
 *  2. Linearise every track at its PCA to the beam spot and find all track
 *     pairs whose 3D DCA is below the pair cut (retry with 3x cut if none).
 *  3. Group tracks into vertex candidates as connected components of the
 *     pair graph; reject pairs whose midpoint is far from the candidate
 *     median and regroup (splits merged pile-up candidates).
 *  4. Fit each candidate with a weighted least-squares vertex fit:
 *       minimise sum_i r_i^T W_i r_i, r_i = perpendicular residual of
 *       the vertex to track i, W_i from a (sigma_rphi, sigma_z, MS) model
 *       or from the track covariance. Tracks are re-linearised at their
 *       helix PCA to the current vertex each iteration; the worst track is
 *       dropped while its chi2 exceeds the cut. Optional beam-spot constraint.
 *  5. Optionally absorb unused tracks compatible with a vertex and refit.
 *  6. Write SvtxVertex_v3 objects to an SvtxVertexMap (track ids are the
 *     Tpc_PolyTrack track ids) and, optionally, the collision-vertex fields
 *     of the Tpc_PolyTrackVertexContainer.
 */

#include <fun4all/SubsysReco.h>

#include <Eigen/Dense>

#include <string>
#include <vector>

class PHCompositeNode;
class SvtxVertexMap;
class Tpc_PolyTrack;
class Tpc_PolyTrackContainer;
class Tpc_PolyTrackVertexContainer;

class Tpc_PolyTrackVertexFinder : public SubsysReco
{
 public:
  explicit Tpc_PolyTrackVertexFinder(const std::string& name = "Tpc_PolyTrackVertexFinder");
  ~Tpc_PolyTrackVertexFinder() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // --- I/O ---
  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setVertexMapName(const std::string& n) { m_vertexMapName = n; }
  //! also overwrite collision-vertex fields of Tpc_PolyTrackVertexContainer
  //! (Tpc_PolyTrackVertexer must run before this module)
  void setFillPolyVertexContainer(bool v = true) { m_fillPolyVertexContainer = v; }
  void setPolyVertexContainerName(const std::string& n) { m_polyVertexNodeName = n; }

  // --- field / beam ---
  void setMagneticFieldTesla(double v) { m_magneticFieldTesla = v; }
  void zeroField(bool v = true) { m_zeroField = v; }
  void setBeamSpot(double x, double y) { m_beamX = x; m_beamY = y; }
  void setBeamSpotSigma(double sx, double sy) { m_beamSigmaX = sx; m_beamSigmaY = sy; }
  void setUseBeamSpotConstraint(bool v = true) { m_useBeamSpotConstraint = v; }
  void setBeamCrossing(short int c) { m_beamCrossing = c; }

  // --- track selection ---
  void setMinClusters(unsigned int v) { m_minClusters = v; }
  void setTrackPtCut(double v) { m_minPt = v; }
  void setTrackQualityCut(double v) { m_maxChi2Ndf = v; }
  void setMaxDcaXY(double v) { m_maxDcaXY = v; }
  void setMaxAbsZ0(double v) { m_maxAbsZ0 = v; }

  // --- pair finding (PHSimpleVertexFinder-like) ---
  void setDcaCut(double v) { m_pairDcaCut = v; }
  void setOutlierPairCut(double v) { m_outlierPairCut = v; }
  void setMinTracksPerVertex(unsigned int v) { m_minTracksPerVertex = v; }

  // --- vertex fit ---
  //! constant track position resolution at the vertex [cm]
  void setTrackResolution(double sigma_rphi, double sigma_z)
  {
    m_sigmaRPhi = sigma_rphi;
    m_sigmaZ = sigma_z;
  }
  //! multiple-scattering term: sigma_ms = k / p   [cm * GeV]
  void setMultipleScatteringTerm(double k) { m_msTerm = k; }
  //! use the (x,y,z) block of Tpc_PolyTrack::get_cov instead of the model
  void setUseTrackCovariance(bool v = true) { m_useTrackCovariance = v; }
  void setMaxTrackChi2(double v) { m_maxTrackChi2 = v; }
  void setMaxIterations(unsigned int v) { m_maxIterations = v; }
  void setConvergenceTolerance(double v) { m_convergenceTol = v; }
  void setAbsorbUnusedTracks(bool v = true) { m_absorbUnusedTracks = v; }

  //! geometric model of one track used by the vertex finder
  struct TrackModel
  {
    unsigned int index{0};     //!< position in the input container
    unsigned int track_id{0};  //!< Tpc_PolyTrack::get_track_id()
    unsigned int nclusters{0};
    bool is_line{false};
    // reference point (track PCA to the z axis for helices)
    double x0{0.0};
    double y0{0.0};
    double z0{0.0};
    double phi0{0.0};   //!< transverse direction angle at the reference point
    double tanl{0.0};   //!< pz / pt
    double p{0.0};      //!< total momentum (GeV), 0 for lines
    // helix
    double xc{0.0};
    double yc{0.0};
    double radius{0.0};
    double h{1.0};      //!< +1 counter-clockwise, -1 clockwise
    // line
    Eigen::Vector3d dir{0.0, 0.0, 1.0};
    // optional measured covariance of (x,y,z)
    bool has_cov{false};
    Eigen::Matrix3d cov{Eigen::Matrix3d::Zero()};
  };

  //! track linearised at a point: closest point + unit tangent there
  struct LinearTrack
  {
    bool ok{false};
    Eigen::Vector3d point{Eigen::Vector3d::Zero()};
    Eigen::Vector3d dir{0.0, 0.0, 1.0};
  };

  struct VertexFit
  {
    bool ok{false};
    Eigen::Vector3d position{Eigen::Vector3d::Zero()};
    Eigen::Matrix3d covariance{Eigen::Matrix3d::Zero()};
    double chi2{0.0};
    int ndf{0};
    std::vector<unsigned int> tracks;  //!< indices into m_models
  };

 private:
  struct TrackPair
  {
    unsigned int a{0};
    unsigned int b{0};
    double dca{0.0};
    Eigen::Vector3d pca_a{Eigen::Vector3d::Zero()};
    Eigen::Vector3d pca_b{Eigen::Vector3d::Zero()};
    Eigen::Vector3d midpoint() const { return 0.5 * (pca_a + pca_b); }
  };

  bool getNodes(PHCompositeNode* topNode);
  bool createNodes(PHCompositeNode* topNode);

  bool buildModel(const Tpc_PolyTrack* trk, unsigned int index, TrackModel& model) const;
  LinearTrack linearize(const TrackModel& model, const Eigen::Vector3d& target) const;
  Eigen::Matrix3d weightMatrix(const TrackModel& model, const LinearTrack& lin) const;

  static double dcaTwoLines(const LinearTrack& t1, const LinearTrack& t2,
                            Eigen::Vector3d& pca1, Eigen::Vector3d& pca2);
  std::vector<TrackPair> findPairs(const std::vector<LinearTrack>& lin, double dcacut) const;
  static std::vector<std::vector<unsigned int>> connectedComponents(
      unsigned int nnodes, const std::vector<TrackPair>& pairs);
  std::vector<TrackPair> removeOutlierPairs(const std::vector<TrackPair>& pairs) const;

  VertexFit fitVertex(std::vector<unsigned int> tracks, const Eigen::Vector3d& seed) const;
  double trackChi2(const TrackModel& model, const Eigen::Vector3d& vtx) const;

  void writeVertices(const std::vector<VertexFit>& vertices);

  // nodes
  std::string m_inputNodeName{"TPC_POLYTRACKS"};
  std::string m_vertexMapName{"TpcPolyVertexMap"};
  std::string m_polyVertexNodeName{"TPC_POLYTRACKVERTICES"};
  Tpc_PolyTrackContainer* m_polyTracks{nullptr};
  SvtxVertexMap* m_vertexMap{nullptr};
  Tpc_PolyTrackVertexContainer* m_polyVertices{nullptr};
  bool m_fillPolyVertexContainer{false};

  // field / beam
  double m_magneticFieldTesla{1.4};
  bool m_zeroField{false};
  double m_beamX{0.0};
  double m_beamY{0.0};
  double m_beamSigmaX{0.05};
  double m_beamSigmaY{0.05};
  bool m_useBeamSpotConstraint{false};
  short int m_beamCrossing{0};

  // selection
  unsigned int m_minClusters{10};
  double m_minPt{0.1};
  double m_maxChi2Ndf{1.0e9};
  double m_maxDcaXY{5.0};
  double m_maxAbsZ0{150.0};

  // pair finding (TPC-only defaults, cm)
  double m_pairDcaCut{0.3};
  double m_outlierPairCut{0.5};
  unsigned int m_minTracksPerVertex{2};

  // fit (TPC-only defaults, cm) - tune on data/simulation
  double m_sigmaRPhi{0.1};
  double m_sigmaZ{0.15};
  double m_msTerm{0.05};
  bool m_useTrackCovariance{false};
  double m_maxTrackChi2{12.0};  // ~99.75% for 2 dof
  unsigned int m_maxIterations{30};
  double m_convergenceTol{1.0e-4};
  bool m_absorbUnusedTracks{true};

  // per-event work space
  std::vector<TrackModel> m_models;

  // statistics
  unsigned long m_nEvents{0};
  unsigned long m_nEventsWithVertex{0};
  unsigned long m_nVertices{0};
};

#endif
