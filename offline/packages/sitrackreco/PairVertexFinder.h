// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_PAIRVERTEXFINDER_H
#define SITRACKRECO_PAIRVERTEXFINDER_H

// Track-source independent collision-vertex finder (no Fun4All dependence), the algorithm of
// Tpc_PolyTrackVertexFinder (TpcVertexer branch) moved into one class so that silicon-only,
// TPC-only and Si+TPC tracks use exactly the same code:
//
//  1. every track is a TrackModel: helix (circle centre, radius, helicity, reference point,
//     direction, tanl) or straight line, plus an optional position resolution of its own
//  2. tracks linearised at the beam spot; all pairs with a 3D DCA below the pair cut
//     (pair finding as PHSimpleVertexFinder; helices re-linearised at the pair midpoint;
//     retry with 3x the cut if no pair passes)
//  3. candidates = connected components of the pair graph after removing pairs whose midpoint
//     is farther than the outlier cut from the candidate median (re-examined in later passes:
//     splits pile-up candidates)
//  4. weighted least-squares fit per candidate: minimise sum r_i^T W_i r_i, r_i = vertex - PCA
//     of track i, W_i = B C^-1 B^T in the plane normal to the track, C from (sigma_rphi,
//     sigma_z) + multiple scattering k/p (or the track covariance); re-linearised every
//     iteration; the worst track is dropped while its chi2 exceeds the cut; optional
//     beam-spot constraint in x, y
//  5. optionally unused tracks join their most compatible vertex, which is refitted
//
// Units cm, GeV.  The caller groups the tracks (per bunch crossing, for example) and calls
// find() once per group.

#include <Eigen/Dense>

#include <vector>

class PairVertexFinder
{
 public:
  struct Config
  {
    // beam spot (transverse) used for the linearisation, the selection and the constraint
    double beamX{0.0};
    double beamY{0.0};
    double beamSigmaX{0.05};
    double beamSigmaY{0.05};
    bool useBeamSpotConstraint{false};
    // track and pair acceptance around the beam
    double maxDcaXY{5.0};
    double maxAbsZ{30.0};
    // pair finding
    double pairDcaCut{0.3};
    double outlierPairCut{0.5};
    bool retryLooseDcaCut{true};
    unsigned int minTracksPerVertex{2};
    // fit: default position resolution of a track at the vertex [cm], multiple scattering k/p
    double sigmaRPhi{0.1};
    double sigmaZ{0.15};
    double msTerm{0.05};
    bool useTrackCovariance{false};
    double maxTrackChi2{12.0};  // ~99.75% for 2 dof
    unsigned int maxIterations{30};
    double convergenceTol{1.0e-4};
    bool absorbUnusedTracks{true};
  };

  //! geometric model of one track
  struct TrackModel
  {
    unsigned int track_id{0};  //!< id written to the vertex (caller's convention)
    bool is_line{false};
    // reference point (PCA to the beam for helices) and direction of motion there
    double x0{0.0};
    double y0{0.0};
    double z0{0.0};
    double phi0{0.0};
    double tanl{0.0};
    double p{0.0};  //!< total momentum [GeV], 0: no multiple-scattering term
    // helix: centre, radius, h = +1 counter-clockwise, -1 clockwise (seen from +z)
    double xc{0.0};
    double yc{0.0};
    double radius{0.0};
    double h{1.0};
    // line direction (unit)
    Eigen::Vector3d dir{0.0, 0.0, 1.0};
    // resolution of this track at the vertex [cm]; <= 0: Config value
    double sigmaRPhi{-1.0};
    double sigmaZ{-1.0};
    // optional measured covariance of (x, y, z)
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
    std::vector<unsigned int> tracks;  //!< indices into the models given to find()
  };

  // per-call statistics (find)
  struct Stats
  {
    unsigned int pairs{0};
    unsigned int cleanPairs{0};
    unsigned int candidates{0};
    bool looseRetry{false};
  };

  PairVertexFinder() = default;
  explicit PairVertexFinder(const Config& c)
    : m_cfg(c)
  {
  }

  Config& config() { return m_cfg; }
  const Config& config() const { return m_cfg; }

  //! helix model from the perigee convention of Si_Trajectory / SiTpc_Track: pca (x, y), z0 at the
  //! pca, phi = direction of motion at the pca, tanl, circle centre, radius, helicity (+1 CCW),
  //! pt [GeV] (NaN: no multiple scattering term).  A non-finite or very large radius gives a line.
  static bool fromPerigee(unsigned int id, double pcaX, double pcaY, double z0, double phi, double tanl,
                          double cx, double cy, double radius, int helicity, double pt, TrackModel& m);

  //! helix through the silicon points with a given curvature (from a longer fit, e.g. the Si+TPC
  //! radius and helicity): circle of radius R through the innermost and outermost point (side
  //! from the helicity, short arc), reference point = innermost point, z = z_ref + tanl * s with
  //! z_ref from all points weighted by their z resolution (MVTX r < 6 cm: 1, INTT: 1e-4; arc length s
  //! from the innermost point).  Not finite radius: line
  //! through the two points.  Points: (x, y, z) [cm], at least two at different radius.
  static bool fromPointsWithCurvature(unsigned int id, const std::vector<Eigen::Vector3d>& points, double radius, int helicity,
                                      double tanl, double pt, TrackModel& m);

  //! the track's PCA to the beam is within maxDcaXY and maxAbsZ
  bool acceptAtBeam(const TrackModel& m) const;

  //! vertices of one group of tracks, largest first
  std::vector<VertexFit> find(const std::vector<TrackModel>& models, Stats* stats = nullptr) const;

  //! weighted least-squares fit of ALL given tracks (no pair finding), started at `seed`: the step 4
  //! of find() - re-linearised every iteration, the worst track dropped while its chi2 exceeds
  //! maxTrackChi2 (keeping >= minTracksPerVertex), beam-spot constraint as configured.
  //! VertexFit::tracks are indices into `models`.  Used to refit a vertex without some of its
  //! tracks (e.g. the daughters of a D0 candidate) with the same weights.
  VertexFit fit(const std::vector<TrackModel>& models, const Eigen::Vector3d& seed) const;

  LinearTrack linearize(const TrackModel& m, const Eigen::Vector3d& target) const;
  double trackChi2(const TrackModel& m, const Eigen::Vector3d& vtx) const;

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

  Eigen::Matrix3d weightMatrix(const TrackModel& m, const LinearTrack& lin) const;
  static double dcaTwoLines(const LinearTrack& t1, const LinearTrack& t2, Eigen::Vector3d& pca1, Eigen::Vector3d& pca2);
  std::vector<TrackPair> findPairs(const std::vector<TrackModel>& models, const std::vector<LinearTrack>& lin, double dcacut) const;
  static std::vector<std::vector<unsigned int>> connectedComponents(unsigned int nnodes, const std::vector<TrackPair>& pairs);
  std::vector<TrackPair> removeOutlierPairs(unsigned int nnodes, const std::vector<TrackPair>& pairs) const;
  VertexFit fitVertex(const std::vector<TrackModel>& models, std::vector<unsigned int> tracks, const Eigen::Vector3d& seed) const;

  Config m_cfg;
};

#endif
