#ifndef TPCTRACKRECO_TPCCROSSINGTRAJECTORY_H
#define TPCTRACKRECO_TPCCROSSINGTRAJECTORY_H

#include <phool/PHObject.h>

#include <array>
#include <cmath>
#include <iostream>
#include <limits>

class TpcCrossingTrajectory : public PHObject
{
 public:
  enum StateIndex : unsigned int { X = 0, Y, Z, Phi, Theta, QOverP, StateSize };
  struct LayerState
  {
    unsigned int layer{0};
    float x{0};
    float y{0};
    float z{0};
    float phi{0};
    bool valid{false};
  };

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new TpcCrossingTrajectory(*this); }

  unsigned int get_parent_track_id() const { return m_parentTrackId; }
  void set_parent_track_id(unsigned int value) { m_parentTrackId = value; }
  unsigned int get_source_assembled_track_id() const { return m_sourceAssembledTrackId; }
  void set_source_assembled_track_id(unsigned int value) { m_sourceAssembledTrackId = value; }
  short get_crossing() const { return m_crossing; }
  void set_crossing(short value) { m_crossing = value; }
  short get_reference_crossing() const { return m_referenceCrossing; }
  void set_reference_crossing(short value) { m_referenceCrossing = value; }
  float get_state(unsigned int i) const { return i < StateSize ? m_state[i] : nan(); }
  void set_state(unsigned int i, float value) { if (i < StateSize) m_state[i] = value; }
  float get_delta(unsigned int i) const { return i < StateSize ? m_delta[i] : nan(); }
  void set_delta(unsigned int i, float value) { if (i < StateSize) m_delta[i] = value; }
  float get_covariance(unsigned int i, unsigned int j) const { return i < StateSize && j < StateSize ? m_covariance[i * StateSize + j] : nan(); }
  void set_covariance(unsigned int i, unsigned int j, float value) { if (i < StateSize && j < StateSize) m_covariance[i * StateSize + j] = value; }
  float get_linear_chi2() const { return m_linearChi2; }
  void set_linear_chi2(float value) { m_linearChi2 = value; }
  unsigned int size_layer_states() const { return m_nLayerStates; }
  const LayerState* get_layer_state(unsigned int i) const { return i < m_nLayerStates ? &m_layerStates[i] : nullptr; }
  void add_layer_state(const LayerState& value) { if (m_nLayerStates < m_layerStates.size()) m_layerStates[m_nLayerStates++] = value; }

 private:
  static float nan() { return std::numeric_limits<float>::quiet_NaN(); }
  unsigned int m_parentTrackId{0};
  unsigned int m_sourceAssembledTrackId{0};
  short m_crossing{0};
  short m_referenceCrossing{0};
  std::array<float, StateSize> m_state{};
  std::array<float, StateSize> m_delta{};
  std::array<float, StateSize * StateSize> m_covariance{};
  float m_linearChi2{nan()};
  std::array<LayerState, 7> m_layerStates{};
  unsigned int m_nLayerStates{0};

  ClassDefOverride(TpcCrossingTrajectory, 2)
};

#endif
