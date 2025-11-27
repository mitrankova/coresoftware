// PHSimpleGeomProp.h
#ifndef PHSIMPLEGEOMPROP_H
#define PHSIMPLEGEOMPROP_H


#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <trackbase/ActsGeometry.h>
#include <tpc/TpcGlobalPositionWrapper.h>

#include <trackbase/TrkrClusterIterationMapv1.h>

#include <nanoflann.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <map>
#include <memory>
#include <string>
#include <vector>

// forward declarations
class PHCompositeNode;
class PHG4TpcGeomContainer;
class TrkrCluster;

class PHSimpleGeomProp : public SubsysReco
{
 public:
  PHSimpleGeomProp(const std::string &name = "PHSimpleGeomProp");
  ~PHSimpleGeomProp() override = default;

  // standard Fun4All interface
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // configuration setters
  void set_use_truth_clusters(bool value) { _use_truth_clusters = value; }
  void set_iteration(unsigned int n) { _n_iteration = n; }
  void set_max_seeds(unsigned int n) { _max_seeds = n; }
  void set_max_dist(double d) { _max_dist = d; }
  void set_max_propagation_steps(int n) { _max_propagation_steps = n; }
  void set_num_threads(int n) { m_num_threads = n; }
  void set_do_ghost_rejection(bool v) { m_ghostrejection = v; }
    void set_pp_mode(bool v) { _pp_mode = v; }

  // robust-circle-fit parameters (from S,L,X0/Tukey formalism)
  void set_tukey_C(double c) { m_tukey_C = c; }
  void set_tukey_n(double n) { m_tukey_n = n; }
  void set_circle_max_iter(int n) { m_circle_max_iter = n; }

 private:
  // convenience typedefs
  using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;


  template <typename T>
  struct KDPointCloud
  {
    std::vector<std::vector<T>> pts;

    // number of points
    inline size_t kdtree_get_point_count() const
    {
      return pts.size();
    }

    // Euclidean squared distance between query point p1 and data point idx_p2
    inline T kdtree_distance(const T *p1,
                             const size_t idx_p2,
                             size_t /*size*/) const
    {
      const T d0 = p1[0] - pts[idx_p2][0];
      const T d1 = p1[1] - pts[idx_p2][1];
      const T d2 = p1[2] - pts[idx_p2][2];
      return d0 * d0 + d1 * d1 + d2 * d2;
    }

    // coordinate accessor
    inline T kdtree_get_pt(const size_t idx, int dim) const
    {
      return pts[idx][dim];
    }

    // optional bounding-box: not used
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &) const
    {
      return false;
    }
  };


  using KDTree_t = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<double, KDPointCloud<double>>,
      KDPointCloud<double>, 3>;

  enum class PropagationDirection
  {
    Inward = 0,
    Outward = 1
  };

  // core helpers
  int get_nodes(PHCompositeNode *topNode);

  PositionMap PrepareKDTrees();
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey key, TrkrCluster *cluster) const;

  std::vector<TrkrDefs::cluskey> PropagateTrackSimple(
      TrackSeed *track,
      const std::vector<TrkrDefs::cluskey> &seed_keys,
      PropagationDirection direction,
      const PositionMap &globalPositions) const;

  TrackSeed_v2 makeSeedFromChain(
      const std::vector<TrkrDefs::cluskey> &chain,
      const PositionMap &globalPositions) const;

  // simple filtering/publishing
  std::vector<std::vector<TrkrDefs::cluskey>> RemoveBadClusters(
      const std::vector<std::vector<TrkrDefs::cluskey>> &chains,
      const PositionMap &globalPositions) const;

  void rejectAndPublishSeeds(std::vector<TrackSeed_v2> &seeds,
                             const PositionMap &positions,
                             std::vector<float> &trackChi2);

  void publishSeeds(const std::vector<TrackSeed_v2> &seeds);

 private:
  // node pointers
  ActsGeometry *m_tgeometry = nullptr;
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  TrkrClusterContainer *_cluster_map = nullptr;
  TrackSeedContainer *_track_map = nullptr;
  TrkrClusterIterationMapv1 *_iteration_map = nullptr;

  // geometry
  std::vector<double> radii;

  // KD-tree structures (per layer)
  std::vector<std::shared_ptr<KDPointCloud<double>>> _ptclouds;
  std::vector<std::shared_ptr<KDTree_t>> _kdtrees;

  // configuration
  bool _use_truth_clusters = false;
  unsigned int _n_iteration = 0;
  unsigned int _max_seeds = 0;  // 0 = no limit
  double _max_dist = 5.0;       // cm, geometric matching window
  int _max_propagation_steps = 60;
  int m_num_threads = 1;
  bool m_ghostrejection = true;

  // robust circle fit parameters (from thesis figs)
  double m_tukey_C = 3.0;     // C_T in Eq. (3.8)
  double m_tukey_n = 1.5;     // n in Eq. (3.8)
  int m_circle_max_iter = 10; // max robust-fit iterations

  // use distortion-corrected positions? (same semantics as KF module)
  bool _pp_mode = false;
};

#endif  // PHSIMPLEGEOMPROP_H
