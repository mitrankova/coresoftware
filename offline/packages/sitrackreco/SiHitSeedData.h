#ifndef SIHITSEEDDATA_H
#define SIHITSEEDDATA_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

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
  double z = 0;       // notebook nominal Vz [cm]
  double phi = 0;     // notebook intrinsic Uphi [rad]
  double zbin = 0;
  double phibin = 0;
};

struct SiClusterPoint
{
  TrkrDefs::cluskey key = TrkrDefs::CLUSKEYMAX;
  int layer = -1;
  int stave = -1;
  int chip = -1;
  int ladderphi = -1;
  int ladderz = -1;
  float adc = 0;
  double lx = 0;
  double ly = 0;
  double z = 0;
  double phi = 0;
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

struct SiHitChain
{
  int id = -1;
  std::vector<int> blob_ids;
  std::vector<int> hit_ids;
  std::vector<int> layers;
  std::vector<int> point_hit_ids;
  std::vector<double> z;
  std::vector<double> phi;
  std::vector<SiChainStep> steps;
  int n_mvtx = 0;
  int n_intt = 0;
  double score = 0;
};

// Transient Fun4All node payload produced by SiHitSeedReco and consumed by
// SiHitSeedQA / SiHitSeedDisplay.  It is not written to DST.
class SiHitSeedEvent : public PHObject
{
 public:
  ~SiHitSeedEvent() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SiHitSeedEvent: event=" << event
       << " hits=" << hits.size()
       << " clusters=" << clusters.size()
       << " blobs=" << blobs.size()
       << " chains=" << chains.size()
       << std::endl;
  }

  void Reset() override
  {
    hits.clear();
    clusters.clear();
    blobs.clear();
    chains.clear();
    vertex_z = std::numeric_limits<double>::quiet_NaN();
    vertex_z_linefit = std::numeric_limits<double>::quiet_NaN();
    vertex_tracklet_npairs = 0;
    vertex_tracklet_npeak = 0;
    dz_vs_z_intercept_bins = std::numeric_limits<double>::quiet_NaN();
    dz_vs_z_slope_bins_per_cm = std::numeric_limits<double>::quiet_NaN();
    n_vertex_links = 0;
  }

  int isValid() const override { return 1; }

  unsigned long long event = 0;
  std::vector<SiHitPoint> hits;
  std::vector<SiClusterPoint> clusters;
  std::vector<SiBlob> blobs;
  std::vector<SiHitChain> chains;

  // vertex_z [cm, nominal z frame]: tracklet vertex found BEFORE chain finding, independent of
  // the seed z window (peak of the z0 distribution of all MVTX blob pairs, see SiHitSeedReco).
  double vertex_z = std::numeric_limits<double>::quiet_NaN();
  int vertex_tracklet_npairs = 0;  // pairs filled into the z0 histogram
  int vertex_tracklet_npeak = 0;   // pairs within the peak window
  // vertex_z_linefit [cm]: z where the line fit of dz[bins] vs z over the ACCEPTED seed links
  // crosses dz = 0.  Biased by the seed window (it can only see links the window accepted).
  double vertex_z_linefit = std::numeric_limits<double>::quiet_NaN();
  double dz_vs_z_intercept_bins = std::numeric_limits<double>::quiet_NaN();
  double dz_vs_z_slope_bins_per_cm = std::numeric_limits<double>::quiet_NaN();
  int n_vertex_links = 0;
};

#endif
