#pragma once

// SiTpc_ResolutionQA
//
// Same resolution QA for the own Si+TPC tracks (SITPC_TRACKS + SiTpcPolyVertexMap) and for SvtxTrack
// based tracks (e.g. the ACTS fit: SvtxTrackMap + SvtxVertexMap), so a run before and a run after the
// ACTS fit (or both sources in one run) can be compared with identical cuts and the identical vertex fit.
// A source whose track node is missing is skipped.
//
// Per track (tree "tracks", src = 0 own, 1 svtx):
//   - DCA_xy to the beam line (x0 + dxdz z, y0 + dydz z), signed: (pca - beam) . (-sin phi, cos phi)
//     (the Acts perigee convention).  Width vs pT = track resolution (+) beam size (~13 um).
//   - unbiased DCA_xy and DCA_z to the vertex: the vertex of the track's collision (vertex finder's track
//     list) refitted WITHOUT the track by the QA vertex fit, + the vertex error projected on the track
//     (sigma_track^2 ~ width^2 - <sigma_vtx^2>).
// Per vertex (tree "vertices"):
//   - QA refit of all its tracks (no beam constraint), the vertex finder's own position and errors
//   - split-vertex resolution: tracks split randomly in two halves, both fitted, (v1 - v2) / sqrt(2)
//     vs tracks per half (same fit for every source: compares the tracks, not the vertex finders).
// Tree "matched": own and svtx track of the same particle in the same event (same crossing, closest in
// phi, eta, pT), for the comparison on one track sample.
//
// QA vertex fit: straight lines at the PCA, minimises sum (d_xy / s_xy)^2 + (d_z / s_z)^2 with the same
// pT-dependent weights for every source (setTrackErrorModel), iterative removal of tracks with chi2 > cut.
//
// Frames: give each source its own beam line in its own frame (own frame: (0, 0) with the Si beam-spot
// correction; ACTS / official frame: the Si beam spot).

#include <fun4all/SubsysReco.h>

#include <array>
#include <memory>
#include <string>

class PHCompositeNode;
class TFile;
class TTree;

class SiTpc_ResolutionQA : public SubsysReco
{
 public:
  explicit SiTpc_ResolutionQA(const std::string& name = "SiTpc_ResolutionQA",
                              const std::string& fileName = "sitpc_resolution_qa.root");
  ~SiTpc_ResolutionQA() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // ---- own reconstruction (SiTpc_TrackContainer + SvtxVertexMap with SiTpc track ids)
  void setOwnNodes(const std::string& tracks, const std::string& vertices)
  {
    m_ownTrackNode = tracks;
    m_ownVertexNode = vertices;
  }
  void setOwnBeamLineCm(double x0, double y0, double dxdz = 0, double dydz = 0) { m_ownBeam = {x0, y0, dxdz, dydz}; }

  // ---- SvtxTrack source (SvtxTrackMap + SvtxVertexMap with SvtxTrack ids), e.g. the ACTS fit
  void setSvtxNodes(const std::string& tracks, const std::string& vertices)
  {
    m_svtxTrackNode = tracks;
    m_svtxVertexNode = vertices;
  }
  void setSvtxBeamLineCm(double x0, double y0, double dxdz = 0, double dydz = 0) { m_svtxBeam = {x0, y0, dxdz, dydz}; }

  // ---- crossings of the own tracks (assembled TPC track id -> crossing)
  void setCrossingDecisionNodeName(const std::string& s) { m_crossingNode = s; }

  // ---- QA vertex fit
  // track errors at the vertex [cm]: s^2 = a^2 + (b / p)^2, p in GeV (same for every source)
  void setTrackErrorModel(double axy, double bxy, double az, double bz) { m_errModel = {axy, bxy, az, bz}; }
  void setVertexChi2Cut(double c) { m_chi2Cut = c; }
  // tracks used in the vertex fits
  void setMinPt(double pt) { m_minPt = pt; }
  void setMinMvtx(int n) { m_minMvtx = n; }
  void setMinTpc(int n) { m_minTpc = n; }
  // unbiased DCA: at least this many OTHER tracks in the refitted vertex
  void setMinOtherTracks(int n) { m_minOther = n; }
  // split vertices: at least this many tracks (two halves of n / 2)
  void setMinSplitTracks(int n) { m_minSplit = n; }

  // ---- own <-> svtx matching
  void setMatchWindow(double dphi, double deta, double relpt) { m_match = {dphi, deta, relpt}; }

 private:
  struct Rows;
  std::unique_ptr<Rows> m_rows;

  std::string m_fileName;
  std::string m_ownTrackNode{"SITPC_TRACKS"};
  std::string m_ownVertexNode{"SiTpcPolyVertexMap"};
  std::string m_svtxTrackNode{"SvtxTrackMap"};
  std::string m_svtxVertexNode{"SvtxVertexMap"};
  std::string m_crossingNode{"TPC_CROSSING_DECISIONS"};

  std::array<double, 4> m_ownBeam{0, 0, 0, 0};
  std::array<double, 4> m_svtxBeam{0, 0, 0, 0};
  std::array<double, 4> m_errModel{0.0020, 0.0060, 0.0040, 0.0100};
  double m_chi2Cut{16};
  double m_minPt{0.3};
  int m_minMvtx{3};
  int m_minTpc{10};
  int m_minOther{2};
  int m_minSplit{4};
  std::array<double, 3> m_match{0.02, 0.02, 0.1};

  int m_event{0};
  TFile* m_file{nullptr};
  TTree* m_trackTree{nullptr};
  TTree* m_vertexTree{nullptr};
  TTree* m_matchTree{nullptr};
};
