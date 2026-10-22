// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKINGDIAGNOSTICS_SITPCVERTEXQA_H
#define TRACKINGDIAGNOSTICS_SITPCVERTEXQA_H

// Compares the collision vertices of one event, crossing by crossing:
//   resolver : TPC_CROSSING_DECISIONS (TpcCrossingResolver, read through TPC_ASSEMBLEDTRACKS),
//              silicon_vertex_z of the tier 0 / 1
//              decisions with a collision (n = number of such tracks)
//   si, tpc, sitpc : SvtxVertexMaps of SiTpc_VertexFinder (Source::Si, Tpc, SiTpc); per crossing
//              the vertex with the most tracks (n_vtx = vertices of the source in that crossing)
//
// Tree "crossings": one entry per (event, crossing) present in any of the inputs
//   event, crossing, res_n, res_z,
//   {si,tpc,sitpc}_nvtx, _n, _x, _y, _z, _ex, _ey, _ez, _chi2ndf   (NaN / 0 if none)
// End: median and robust sigma of the z differences and of the pulls dz / sqrt(e1^2 + e2^2)
// (pull width ~ 1 when the vertex errors are right: tune setTrackResolution with it).

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

class SiTpc_VertexQA : public SubsysReco
{
 public:
  explicit SiTpc_VertexQA(const std::string& name = "SiTpc_VertexQA", const std::string& outfile = "sitpc_vertex_qa.root");
  ~SiTpc_VertexQA() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setVertexMapNames(const std::string& si, const std::string& tpc, const std::string& sitpc)
  {
    m_mapNames = {si, tpc, sitpc};
  }
  void setCrossingDecisionNodeName(const std::string& n) { m_decisionNodeName = n; }
  void setAssembledTrackNodeName(const std::string& n) { m_assembledNodeName = n; }

 private:
  static constexpr int kN = 3;  // si, tpc, sitpc
  std::string m_outfileName;
  std::array<std::string, kN> m_mapNames{{"SiPolyVertexMap", "TpcPolyVertexMap", "SiTpcPolyVertexMap"}};
  std::string m_decisionNodeName{"TPC_CROSSING_DECISIONS"};
  std::string m_assembledNodeName{"TPC_ASSEMBLEDTRACKS"};

  TFile* m_file{nullptr};
  TTree* m_tree{nullptr};
  unsigned int m_event{0};
  int m_crossing{0};
  int m_resN{0};
  double m_resZ{0};
  std::array<int, kN> m_nvtx{};
  std::array<int, kN> m_n{};
  std::array<double, kN> m_x{}, m_y{}, m_z{}, m_ex{}, m_ey{}, m_ez{}, m_chi2ndf{};

  // End summary: z differences and pulls, [comparison]
  //   0 si - sitpc, 1 tpc - sitpc, 2 resolver - sitpc, 3 si - resolver, 4 tpc - si
  static constexpr int kCmp = 5;
  std::array<std::vector<double>, kCmp> m_dz;
  std::array<std::vector<double>, kCmp> m_pull;
  std::array<unsigned long, kN> m_nCrossingsWithVertex{};
  unsigned long m_nCrossings{0};
};

#endif
