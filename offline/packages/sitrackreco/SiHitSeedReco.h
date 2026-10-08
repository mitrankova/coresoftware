#ifndef SIHITSEEDRECO_H
#define SIHITSEEDRECO_H

#include "SiHitSeedData.h"

#include <fun4all/SubsysReco.h>
#include <Eigen/Core>

#include <array>
#include <string>
#include <vector>

class PHCompositeNode;

class SiHitSeedReco : public SubsysReco
{
 public:
  explicit SiHitSeedReco(const std::string& name = "SiHitSeedReco");
  ~SiHitSeedReco() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setHitNodeName(const std::string& s) { m_hitNodeName = s; }
  void setClusterNodeName(const std::string& s) { m_clusterNodeName = s; }
  void setOutputNodeName(const std::string& s) { m_outputNodeName = s; }

  void setMinLayers(int v) { m_minLayers = v; }
  void setMinMvtxLayers(int v) { m_minMvtxLayers = v; }
  void setMinInttLayers(int v) { m_minInttLayers = v; }
  void setRefitWithIntt(bool v) { m_refitWithIntt = v; }

  // ---------------------------------------------------------------------------
  // Seed z window (outer MVTX hit -> next inner layer), in z bins:
  //   accept if |dz_bins - center(z_outer)| <= half,  dz_bins = zbin_inner - zbin_outer
  //
  // Band mode (DEFAULT): center(z) = offsetBins - shiftAtEdgeBins * z / Zmax, constant half-width.
  //   setSeedZBand(5.0)            -> flat band, dz in [-5, +5] bins for every z   (default)
  //   setSeedZBand(4.0, 4.0)       -> old notebook window: +-4 bins, center +4 at -Zmax,
  //                                   0 at z=0, -4 at +Zmax
  //   setSeedZBand(6.0, 4.0, 1.0)  -> same tilt, 1 bin higher, +-6 bins wide
  //   Zmax = max(|zmin|, |zmax|) of the event's z binning.
  // Parametrized mode: center(z) = interceptBins + slopeBinsPerCm * z
  //   setSeedZWindow(-0.1896, -0.4511, 2.787)  -> 95% band fit from the notebook
  // Centimetre mode: the same straight line but in cm, so it does not change when the
  // event's bin width changes (z bin = (zmax - zmin)/nz is recomputed every event):
  //   accept if |dz_cm - (interceptCm + slope * z_outer_cm)| <= halfWidthCm
  //   For tracks from (0,0,zv): dz = (R_in/R_out - 1)(z_outer - zv), i.e.
  //   slope = R_in/R_out - 1 (-0.196 for 2->1, -0.244 for 1->0) and interceptCm = -slope * zv.
  //   setSeedZWindowCm(0.0, -0.22, 1.5) -> line through (zv = 0), +-1.5 cm
  // ---------------------------------------------------------------------------
  void setSeedZBand(double halfWidthBins, double shiftAtEdgeBins = 0.0, double offsetBins = 0.0)
  {
    m_zMode = ZWindowMode::Band;
    m_zBandHalf = halfWidthBins;
    m_zBandShiftAtEdge = shiftAtEdgeBins;
    m_zBandOffset = offsetBins;
  }
  // Vertex mode (recommended): the window follows the event vertex found before chain finding
  // (vertex_z, tracklet z0 peak).  For each seed pair the exact line through the vertex is used:
  //   center_cm = (R_in/R_out - 1) * (z_outer - vertex_z) + offsetCm,  half = halfWidthCm.
  //   setSeedZWindowVertex(1.5)        -> +-1.5 cm around the line through the vertex
  //   setSeedZWindowVertex(2.0, 0.3)   -> +-2 cm, line moved up by 0.3 cm
  // If no vertex is found in the event, the band settings (setSeedZBand) are used instead.
  void setSeedZWindowVertex(double halfWidthCm, double offsetCm = 0.0)
  {
    m_zMode = ZWindowMode::Vertex;
    m_zVtxHalfCm = halfWidthCm;
    m_zVtxOffsetCm = offsetCm;
  }
  // Tracklet vertex finder: all MVTX blob pairs (2->1 and 1->0) with |dUphi| <= phiHalfBins
  // and |dz| <= maxDzCm are extrapolated to R = 0; vertex_z = peak of that z0 histogram
  // (binCm wide bins, peak = densest peakWidthCm window, then mean of z0 inside it).
  void setTrackletVertex(double phiHalfBins, double maxDzCm, double binCm = 0.2, double peakWidthCm = 1.0)
  {
    m_tvPhiHalf = phiHalfBins;
    m_tvMaxDzCm = maxDzCm;
    m_tvBinCm = binCm;
    m_tvPeakCm = peakWidthCm;
  }
  void setTrackletVertexRange(double zminCm, double zmaxCm) { m_tvZmin = zminCm; m_tvZmax = zmaxCm; }
  // Quality of the tracklet vertex: the peak must contain at least minPairs pairs AND be at least
  // minSignalOverBkg times the flat combinatorial background expected in the same window.
  // Otherwise vertex_z = NaN (no vertex) and the vertex-mode window falls back to setSeedZBand.
  void setTrackletVertexQuality(int minPairs, double minSignalOverBkg = 3.0)
  {
    m_tvMinPeak = minPairs;
    m_tvMinSoB = minSignalOverBkg;
  }
  void setSeedZWindowCm(double interceptCm, double slope, double halfWidthCm)
  {
    m_zMode = ZWindowMode::Centimetre;
    m_zInterceptCm = interceptCm;
    m_zSlopeCm = slope;
    m_zHalfCm = halfWidthCm;
  }
  void setSeedZWindow(double interceptBins, double slopeBinsPerCm, double halfWidthBins)
  {
    m_zMode = ZWindowMode::Parametrized;
    m_zIntercept = interceptBins;
    m_zSlope = slopeBinsPerCm;
    m_zHalf = halfWidthBins;
  }

  // Phi windows in Uphi bins (1 bin = 2*pi/nphi).
  void setPhiHalfWindowBins(double v) { m_seedPhiHalf = v; m_propPhiHalf = v; }  // both at once
  void setSeedPhiHalfWindowBins(double v) { m_seedPhiHalf = v; }
  void setPropagationPhiHalfWindowBins(double v) { m_propPhiHalf = v; }
  void setPropagationZHalfWindowBins(double v) { m_propZHalf = v; }

  // Number of bins. The z range is recomputed each event from the hits (as in the notebook).
  void setZBins(int n) { m_nz = n; }
  // Kept for existing macros: only the bin count is used (zmin/zmax are recomputed per event).
  void setZBinning(int n, double /*zmin*/, double /*zmax*/) { m_nz = n; }
  void setPhiBins(int n) { m_nphi = n; }

 private:
  enum class ZWindowMode { Band, Parametrized, Centimetre, Vertex };

  struct Line3D { Eigen::Vector3d c; Eigen::Vector3d d; };
  struct Candidate
  {
    int blob = -1;
    int hit = -1;
    int layer = -1;
    int outer_hit = -1;
    double z = 0, phi = 0, outer_z = 0, outer_phi = 0;
    double res_z = 0, res_phi = 0, z_center = 0, z_half = 0, score = 0;
  };

  int createNodes(PHCompositeNode*);
  bool fillHits(PHCompositeNode*);
  void fillClusters(PHCompositeNode*);
  void buildBlobs();
  void buildHitGrid();
  void findChains();
  void findTrackletVertex();      // event->vertex_z [cm], before chain finding
  void fitVertexFromSeedLinks();  // event->vertex_z_linefit [cm], after chain finding

  static double wrapDelta(double x, double period);
  static double circularMean(const std::vector<double>& values);
  std::pair<double, double> localCoordinates(int layer, unsigned int row, unsigned int col, int ladderz) const;
  double intrinsicPhi(int layer, int stave, int ladderphi, double lx) const;
  double longitudinalZ(int layer, int chip, int ladderz, double ly) const;
  Eigen::Vector3d xyz(int layer, double z, double phi) const;
  Line3D fitLine(const std::vector<Eigen::Vector3d>& points) const;
  bool propagate(const Line3D&, int layer, const Eigen::Vector3d& near, double& z, double& phi) const;
  void seedZWindow(int outerLayer, double zOuter, double& center, double& half) const;

  // Calls f(hitIndex) for every hit of `layer` in the (z, phi) bin cells overlapping
  // [zbinCenter +- zHalf] x [phibinCenter +- phiHalf] (phi periodic).
  template <class F>
  void forEachHitNear(int layer, double zbinCenter, double zHalf,
                      double phibinCenter, double phiHalf, F&& f) const;

  std::vector<Candidate> seedCandidates(const SiBlob& seed, int innerLayer, const std::vector<char>& used) const;
  bool buildChain(const SiBlob& seed, const Candidate& cand, const std::vector<char>& globallyUsed,
                  SiHitChain& out) const;

  std::string m_hitNodeName = "TRKR_HITSET";
  std::string m_clusterNodeName = "TRKR_CLUSTER";
  std::string m_outputNodeName = "SI_HIT_SEED_EVENT";
  unsigned long long m_event = 0;
  SiHitSeedEvent* m_eventData = nullptr;

  int m_minLayers = 3;
  int m_minMvtxLayers = 2;
  int m_minInttLayers = 0;
  bool m_refitWithIntt = false;

  // seed windows
  double m_seedPhiHalf = 2.0;
  ZWindowMode m_zMode = ZWindowMode::Band;
  double m_zBandHalf = 5.0;
  double m_zBandShiftAtEdge = 0.0;
  double m_zBandOffset = 0.0;
  double m_zInterceptCm = 0.0;
  double m_zSlopeCm = -0.22;
  double m_zHalfCm = 1.5;
  double m_zVtxHalfCm = 1.5;
  double m_zVtxOffsetCm = 0.0;
  // tracklet vertex finder
  double m_tvPhiHalf = 2.0;
  double m_tvMaxDzCm = 4.0;
  double m_tvBinCm = 0.2;
  double m_tvPeakCm = 1.0;
  int m_tvMinPeak = 15;
  double m_tvMinSoB = 3.0;
  double m_tvZmin = -30.0;
  double m_tvZmax = 30.0;
  double m_zIntercept = -0.18960782937488402;
  double m_zSlope = -0.4510638124504644;
  double m_zHalf = 2.7868024058621024;
  // propagation windows
  double m_propPhiHalf = 2.0;
  double m_propZHalf = 5.0;

  int m_nz = 120;
  int m_nphi = 180;
  double m_zmin = -20.0;
  double m_zmax = 20.0;

  // Per-event hit lookup: hits of each (layer, iz, iphi) cell, CSR layout.
  std::vector<int> m_hitBlob;    // blob id of each hit, -1 if not in a blob
  std::vector<int> m_cellStart;  // size 7*nz*nphi + 1
  std::vector<int> m_cellHits;

  const std::array<double, 7> m_radius{{2.523, 3.336, 4.148, 7.188 - 0.0036, 7.732 - 0.0036, 9.680 - 0.0036, 10.262 - 0.0036}};
};

#endif
