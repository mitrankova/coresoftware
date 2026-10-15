#ifndef SIHITSEEDDATA_H
#define SIHITSEEDDATA_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

// Silicon hit in HARDWARE (intrinsic) coordinates, as decoded and aligned by SiHitSeedReco.
struct SiHitPoint
{
  int layer = -1;
  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::HITSETKEYMAX;
  TrkrDefs::hitkey hitkey = TrkrDefs::HITKEYMAX;
  unsigned int row = 0;
  unsigned int col = 0;
  int stave = -1;
  int chip = -1;
  int ladderphi = -1;
  int ladderz = -1;
  float adc = 0;
  double lx = 0;
  double ly = 0;
  double z = 0;       // nominal Vz [cm]
  double phi = 0;     // intrinsic Uphi [rad]
  double zbin = 0;
  double phibin = 0;
};

struct SiBlob
{
  int id = -1;
  int layer = -1;
  std::vector<int> hit_ids;
  std::vector<std::pair<int, int>> cells;
  double z = 0;
  double phi = 0;
  double zbin = 0;
  double phibin = 0;
  double adc_sum = 0;
};

struct SiChainStep
{
  enum Kind
  {
    Seed = 0,
    MvtxPropagation = 1,
    InttPropagation = 2
  };

  int kind = Seed;
  int from_layer = -1;
  int to_layer = -1;
  int blob_id = -1;
  int hit_id = -1;
  double ref_z = 0;
  double ref_phi = 0;
  double pred_z = 0;
  double pred_phi = 0;
  double z = 0;
  double phi = 0;
  double delta_z_bins = 0;
  double delta_phi_bins = 0;
  double res_z_bins = 0;
  double res_phi_bins = 0;
  double z_half_window_bins = 0;
  double phi_half_window_bins = 0;
  double score = 0;
};

// One silicon chain found for one TPC track.  Hits may be shared between chains (of the same or
// of different TPC tracks): the chaining is TPC-driven and does not mask used hits.
// Per-layer vectors (layers, blob_ids, point_hit_ids, z, phi, phys_*, res_*) are ordered
// outermost layer first.
struct SiHitChain
{
  int id = -1;  // index in SiHitSeedEvent::chains
  std::vector<int> blob_ids;
  std::vector<int> hit_ids;  // all hits of the member blobs
  std::vector<int> layers;
  std::vector<int> point_hit_ids;
  std::vector<double> z;    // hardware Vz of the representative hit [cm]
  std::vector<double> phi;  // hardware Uphi of the representative hit [rad]
  std::vector<SiChainStep> steps;
  int n_mvtx = 0;
  int n_intt = 0;
  double score = 0;  // sum of the step scores (window-normalised)

  // ---- step-to-step delta-delta-phi in hardware coordinates (SiHitSeedReco::setHwDdphi*), for
  // layers 0 and 1: pull = (dUphi(L) - offset_L - slope_L * dUphi_pred(L)) / sigma_L with
  // dUphi(L) = Uphi(L) - Uphi(L+1), dUphi_pred = dUphi(L+1) * (r_L - r_{L+1}) / (r_{L+1} - r_PP),
  // PP = next outer chain layer.  NaN if not evaluated (check off, or L, L+1 or PP missing).
  double hw_ddphi_pull[2] = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  // ---- TPC association
  int tpc_track_index = -1;  // index in TPC_POLYTRACKS (or in the crossing-trajectory container)
  unsigned int tpc_track_id = 0;
  unsigned int tpc_assembled_track_id = 0;
  short tpc_crossing = 0;  // crossing assumed for the TPC trajectory
  int seed_outer_layer = -1;  // outer layer of the seed pair (2: pair 2->1, 1: pair 1->0)
  int anchor_blob = -1;       // blob in the TPC anchor window on the outer seed layer

  // ---- physics space: one point per layer, TPC frame (detector -> global -> clamshell half on
  // the TPC beam axis), blob centroid projected on the nominal layer cylinder
  std::vector<double> phys_x;
  std::vector<double> phys_y;
  std::vector<double> phys_z;
  // residuals of those points to the TPC helix at the same radius [cm]
  std::vector<double> res_rdphi;
  std::vector<double> res_dz;
  std::vector<double> res_dz_fit;  // res_dz minus the fitted z model (offset [+ slope s]) [cm]
  double tpc_dz_offset = 0;  // fitted z offset Si - TPC [cm] (z model Offset / OffsetSlope)
  double tpc_dz_slope = 0;   // fitted dz/ds Si - TPC (z model OffsetSlope), ~ tanl correction
  double tpc_chi2_rphi = std::numeric_limits<double>::quiet_NaN();  // sum (r*dphi/sigma)^2
  double tpc_chi2_z = std::numeric_limits<double>::quiet_NaN();     // sum ((dz - model)/sigma)^2, unweighted
  double tpc_chi2 = std::numeric_limits<double>::quiet_NaN();  // chi2 / ndf of the agreement
  int tpc_ndf = 0;
  int tpc_nres = 0;  // chain layers on which the TPC helix could be evaluated (finite r*dphi)
  int rank = -1;          // 0 = best chain of its TPC track (most layers, then best tpc_chi2)
  bool selected = false;  // rank 0 and passing the agreement cut
  int n_layers() const { return n_mvtx + n_intt; }
};

// One entry per TPC track seen by SiHitSeedReco.
struct SiTpcSeedInfo
{
  enum Status
  {
    Chained = 0,
    TooFewClusters = 1,
    FitFailed = 2,
    NoPrediction = 3,   // TPC helix does not reach the outer seed layers
    NoAnchor = 4,       // no silicon hit in the anchor window
    NoChain = 5,        // anchors found, no chain passing the layer requirements
    Rejected = 6        // chains found, best one fails the agreement cut
  };
  int tpc_track_index = -1;
  unsigned int tpc_track_id = 0;
  unsigned int tpc_assembled_track_id = 0;
  short tpc_crossing = 0;
  unsigned int n_clusters = 0;
  int status = TooFewClusters;
  double pt = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  double tanl = std::numeric_limits<double>::quiet_NaN();
  double z0 = std::numeric_limits<double>::quiet_NaN();  // pca z [cm], TPC frame
  int charge = 0;
  int n_anchors = 0;
  int n_chains = 0;
  int selected_chain = -1;  // SiHitChain::id, -1 if none
};

// Transient Fun4All node payload produced by SiHitSeedReco and consumed by SiTrajectoryFitter,
// SiTpcTrackMatcher and SiTpc_EventDisplay.  It is not written to DST.
class SiHitSeedEvent : public PHObject
{
 public:
  ~SiHitSeedEvent() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    unsigned int nSel = 0;
    for (const auto& c : chains)
    {
      nSel += c.selected;
    }
    os << "SiHitSeedEvent: event=" << event
       << " hits=" << hits.size()
       << " blobs=" << blobs.size()
       << " tpc_seeds=" << tpc_seeds.size()
       << " chains=" << chains.size()
       << " selected=" << nSel
       << std::endl;
  }

  void Reset() override
  {
    hits.clear();
    blobs.clear();
    chains.clear();
    tpc_seeds.clear();
  }

  int isValid() const override { return 1; }

  unsigned long long event = 0;
  std::vector<SiHitPoint> hits;
  std::vector<SiBlob> blobs;
  std::vector<SiHitChain> chains;
  std::vector<SiTpcSeedInfo> tpc_seeds;
};

#endif