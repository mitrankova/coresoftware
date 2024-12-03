/*
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef KFPARTICLESPHENIX_KFPARTICLESPHENIX_H
#define KFPARTICLESPHENIX_KFPARTICLESPHENIX_H

#include "KFParticle_eventReconstruction.h"

// There is something broken in this package. clang format puts
// #include "KFParticle_DST.h" first if there is no space but then
// loading KFParticle makes root barf. I have no time to track this down
// right now, it must be something in KFParticle_eventReconstruction.h and
// the include files it uses (some include guard misfiring?)

#include "KFParticle_DST.h"
#include "KFParticle_nTuple.h"

// sPHENIX stuff
#include <fun4all/SubsysReco.h>

// KFParticle stuff
#include <KFParticle.h>

#include <algorithm>  // for max
#include <memory>     // for allocator_traits<>::valu...
#include <string>
#include <utility>  // for pair
#include <vector>   // for vector

class PHCompositeNode;
class TFile;

class KFParticle_sPHENIX : public SubsysReco, public KFParticle_nTuple, public KFParticle_DST, protected KFParticle_eventReconstruction
{
 public:
  KFParticle_sPHENIX();

  explicit KFParticle_sPHENIX(const std::string &name);

  ~KFParticle_sPHENIX() override = default;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  /**
   * If verbosity is > 0, this will print out all candidate information:
   * masses, momenta and positions for mothers, intermediates and final state tracks,
   * PV position, number of vertices and number of tracks in the event (multiplicity)
   */
  void printParticles(const KFParticle &motherParticle,
                      const KFParticle &chosenVertex,
                      const std::vector<KFParticle> &daughterParticles,
                      const std::vector<KFParticle> &intermediateParticles,
                      const int numPVs, const int numTracks);

  int parseDecayDescriptor();

  /// Parameters for the user to vary

  void setDecayDescriptor(const std::string &decayDescriptor) { m_decayDescriptor = decayDescriptor; }

  static const int max_particles = 99;

  void setMotherName(const std::string &mother_name)
  {
    m_mother_name = mother_name;
    m_mother_name_Tools = mother_name;
  }

  void hasIntermediateStates()
  {
    m_has_intermediates = true;
    m_has_intermediates_nTuple = true;
    m_has_intermediates_sPHENIX = true;
    m_has_intermediates_DST = true;
  }

  void setNumberOfTracks(int num_tracks)
  {
    m_num_tracks = num_tracks;
    m_num_tracks_nTuple = num_tracks;
  }

  void setNumberTracksFromIntermeditateState(const std::vector<int> &num_tracks)
  {
    for (unsigned int i = 0; i < num_tracks.size(); ++i)
    {
      m_num_tracks_from_intermediate.push_back(num_tracks[i]);
      m_num_tracks_from_intermediate_nTuple.push_back(num_tracks[i]);
    }
  }

  void setNumberOfIntermediateStates(int n_intermediates)
  {
    m_num_intermediate_states = n_intermediates;
    m_num_intermediate_states_nTuple = n_intermediates;
  }

  void getChargeConjugate()
  {
    m_get_charge_conjugate_nTuple = true;
    m_get_charge_conjugate = true;
  }

  void setDaughters(std::vector<std::pair<std::string, int> /*unused*/> daughter_list)
  {
    for (unsigned int i = 0; i < daughter_list.size(); ++i)
    {
      m_daughter_name.push_back(daughter_list[i].first);
      m_daughter_charge.push_back(daughter_list[i].second);
    }
  }

  void setIntermediateStates(std::vector<std::pair<std::string, int> /*unused*/> intermediate_list)
  {
    for (unsigned int i = 0; i < intermediate_list.size(); ++i)
    {
      m_intermediate_name_ntuple.push_back(intermediate_list[i].first);
      m_intermediate_name.push_back(intermediate_list[i].first);
      m_intermediate_charge.push_back(intermediate_list[i].second);
    }
  }

  void setMinimumMass(float min_mass) { m_min_mass = min_mass; }

  void setMaximumMass(float max_mass) { m_max_mass = max_mass; }

  void setDecayTimeRange(float min_decayTime, float max_decayTime)
  {
    m_min_decayTime = min_decayTime;
    m_max_decayTime = max_decayTime;
  }

  void setDecayLengthRange(float min_decayLength, float max_decayLength)
  {
    m_min_decayLength = min_decayLength;
    m_max_decayLength = max_decayLength;
  }

  void setMinimumTrackPT(float pt) { m_track_pt = pt; }

  void setMaximumTrackPTchi2(float ptchi2) { m_track_ptchi2 = ptchi2; }

  void setMinimumTrackIP(float ip) { m_track_ip = ip; }

  void setMinimumTrackIPchi2(float ipchi2) { m_track_ipchi2 = ipchi2; }

  void setMaximumTrackchi2nDOF(float trackchi2ndof) { m_track_chi2ndof = trackchi2ndof; }

  void setMinMVTXhits(int nHits) { m_nMVTXHits = nHits; }

  void setMinTPChits(int nHits) { m_nTPCHits = nHits; }

  void setMaximumDaughterDCA(float dca) { m_comb_DCA = dca; }

  void setMaximumVertexchi2nDOF(float vertexchi2nDOF) { m_vertex_chi2ndof = vertexchi2nDOF; }

  void setMinimumRadialSV(float min_rad_sv) { m_min_radial_SV = min_rad_sv; }

  void setFlightDistancechi2(float fdchi2) { m_fdchi2 = fdchi2; }

  void setMinDIRA(float dira_min) { m_dira_min = dira_min; }

  void setMaxDIRA(float dira_max) { m_dira_max = dira_max; }

  void setMotherPT(float mother_pt) { m_mother_pt = mother_pt; }

  void setMotherIPchi2(float mother_ipchi2) { m_mother_ipchi2 = mother_ipchi2; }

  void setMaximumMotherVertexVolume(float vertexvol) { m_mother_vertex_volume = vertexvol; }

  void constrainToPrimaryVertex()
  {
    m_constrain_to_vertex = true;
    m_constrain_to_vertex_nTuple = true;
    m_constrain_to_vertex_sPHENIX = true;
  }

  void useFakePrimaryVertex()
  {
    m_use_fake_pv = true;
    m_use_fake_pv_nTuple = true;
  }

  void allowZeroMassTracks() { m_allowZeroMassTracks = true; }

  void dontExtraolateTracksToSV()
  {
    m_extrapolateTracksToSV = m_extrapolateTracksToSV_nTuple = false;
  }

  void constrainIntermediateMasses() { m_constrain_int_mass = true; }

  void setIntermediateMassRange(std::vector<std::pair<float, float> /*unused*/> intermediate_mass_range)
  {
    for (unsigned int i = 0; i < intermediate_mass_range.size(); ++i) m_intermediate_mass_range.push_back(intermediate_mass_range[i]);
  }

  void setIntermediateMinPT(const std::vector<float> &intermediate_min_pt)
  {
    m_intermediate_min_pt = intermediate_min_pt;
  }

  void setIntermediateMinIP(const std::vector<float> &intermediate_min_IP)
  {
    for (unsigned int i = 0; i < intermediate_min_IP.size(); ++i) m_intermediate_min_ip.push_back(intermediate_min_IP[i]);
  }

  void setIntermediateIPRange(const std::vector<std::pair<float, float> /*unused*/> &intermediate_IP_range)
  {
    for (unsigned int i = 0; i < intermediate_IP_range.size(); ++i)
    {
      m_intermediate_min_ip.push_back(intermediate_IP_range[i].first);
      m_intermediate_max_ip.push_back(intermediate_IP_range[i].second);
    }
  }

  void setIntermediateMinIPchi2(const std::vector<float> &intermediate_min_IPchi2)
  {
    for (unsigned int i = 0; i < intermediate_min_IPchi2.size(); ++i) m_intermediate_min_ipchi2.push_back(intermediate_min_IPchi2[i]);
  }

  void setIntermediateIPchi2Range(const std::vector<std::pair<float, float> /*unused*/> &intermediate_IPchi2_range)
  {
    for (unsigned int i = 0; i < intermediate_IPchi2_range.size(); ++i)
    {
      m_intermediate_min_ipchi2.push_back(intermediate_IPchi2_range[i].first);
      m_intermediate_max_ipchi2.push_back(intermediate_IPchi2_range[i].second);
    }
  }

  void setIntermediateMinDIRA(const std::vector<float> &intermediate_min_DIRA)
  {
    for (unsigned int i = 0; i < intermediate_min_DIRA.size(); ++i) m_intermediate_min_dira.push_back(intermediate_min_DIRA[i]);
  }

  void setIntermediateMinFDchi2(const std::vector<float> &intermediate_min_FDchi2)
  {
    for (unsigned int i = 0; i < intermediate_min_FDchi2.size(); ++i) m_intermediate_min_fdchi2.push_back(intermediate_min_FDchi2[i]);
  }

  void setIntermediateMaxVertexVolume(const std::vector<float> &intermediate_max_vertexvol)
  {
    for (unsigned int i = 0; i < intermediate_max_vertexvol.size(); ++i) m_intermediate_vertex_volume.push_back(intermediate_max_vertexvol[i]);
  }

  void use2Dmatching(){ m_use_2D_matching_tools = true; }

  void useMVA() { m_require_mva = true; }

  void setNumMVAPars(unsigned int nPars) { m_nPars = nPars; }

  void setMVAVarList(std::vector<std::string> mva_variable_list)
  {
    for (unsigned int i = 0; i < mva_variable_list.size(); ++i) m_mva_variable_list.push_back(mva_variable_list[i]);
  }

  void setMVAType(const std::string &mva_type) { m_mva_type = mva_type; }

  void setMVAWeightsPath(const std::string &mva_weights_path) { m_mva_path = mva_weights_path; }

  void setMVACutValue(float cut_value) { m_mva_cut_value = cut_value; }

  void saveDST() { m_save_dst = true; }

  void dontSaveTrackContainer() { m_write_track_container = false; }

  void dontSaveParticleContainer() { m_write_particle_container = false; }

  void setContainerName(const std::string &name) { m_container_name = name; }

  void dontSaveOutput() { m_save_output = false; }

  void setOutputName(const std::string &name) { m_outfile_name = name; }

  void doTruthMatching() { m_truth_matching = true; }

  void getDetectorInfo() { m_detector_info = true; }

  void getCaloInfo() { m_calo_info = true; }

  void getAllPVInfo() { m_get_all_PVs = true; }

  /// Use alternate vertex and track fitters
  void setVertexMapNodeName(const std::string &vtx_map_node_name) { m_vtx_map_node_name = m_vtx_map_node_name_nTuple = vtx_map_node_name; }

  /// Use alternate vertex and track fitters
  void setTrackMapNodeName(const std::string &trk_map_node_name) { m_trk_map_node_name = m_trk_map_node_name_nTuple = m_origin_track_map_node_name = trk_map_node_name; }

  void magFieldFile(const std::string &fname) { m_magField = fname; }

  void getField();

 private:
  bool m_has_intermediates_sPHENIX;
  bool m_constrain_to_vertex_sPHENIX;
  bool m_require_mva;
  bool m_save_dst;
  bool m_save_output;
  std::string m_outfile_name;
  TFile *m_outfile;
  std::string m_decayDescriptor;
  std::string m_magField = "FIELDMAP_TRACKING";
};

#endif  // KFPARTICLESPHENIX_KFPARTICLESPHENIX_H
