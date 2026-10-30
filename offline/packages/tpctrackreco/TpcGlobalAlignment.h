// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCGLOBALALIGNMENT_H
#define TPCTRACKRECO_TPCGLOBALALIGNMENT_H

#include <fun4all/SubsysReco.h>

#include <Acts/Definitions/Algebra.hpp>

#include <string>

class ActsGeometry;
class PHCompositeNode;

// Rigid TPC -> silicon alignment applied to the TPC CLUSTER POSITIONS, not to the ACTS surfaces.
//
// Every TPC cluster position in this chain goes through ActsGeometry::transformTpcEnvelopeToWorld
// (TpcPolyClusterTrkrClusterConverter, Tpc_PolyTripletReco, Tpc_PolyTrackReco, TpcClusterizer,
// SiTpcActsFitter, Tpc_PolyClusterResiduals, and SiHitSeedReco in SiTpcBeamAlignment::TpcFrame::World).
// That transform E (envelope -> world) is built by MakeActsGeometry from the TPC placement / tilt stored
// in TPCGEOMCONTAINER (Tracking_Geometry).  This module replaces it, once, at InitRun, by
//
//     E' = A * E ,     A(p) = t + pivot + R (p - pivot) ,    R = Rz(rz) * Ry(ry) * Rx(rx)
//
// The ACTS surfaces are NOT moved.  The TPC surface lookup (get_tpc_surface_from_coords) uses the
// inverse transform, so a moved cluster is still assigned to the surface of its pad.  The converter
// keeps the in-plane (r*phi) and z part of the moved position; the out-of-plane (radial) part is dropped,
// which does not change a track (to first order it is a shift along the track).
//
// Why not the TPC rows of the alignment file: those move the TPC SURFACES.  The poly converter computes
// the local position of a cluster from its global position on the already moved surface
// (local = surface'^-1 * global), so for a pure translation the r*phi and z parts cancel and only the
// radial part survives: the rows do (almost) nothing to the tracks.  Zero them when using this module.
//
// Small-angle meaning of the parameters (pivot = 0):
//     x' = x + ry z - rz y + tx ,   y' = y - rx z + rz x + ty ,   z' = z - ry x + rx y + tz
// A TPC beam line x = a + b z, y = a' + b' z (measured in the world frame WITHOUT this module) is moved
// onto the silicon line x = c + d z, y = c' + d' z by
//     tx = c - a ,  ty = c' - a' ,  ry = d - b ,  rx = b' - d'      (setFromBeamLines)
//
// Register AFTER TrackingInit() (MakeActsGeometry) and BEFORE every module that uses TPC positions
// (InitRun runs in registration order; the transform is replaced at InitRun).
class TpcGlobalAlignment : public SubsysReco
{
 public:
  explicit TpcGlobalAlignment(const std::string& name = "TpcGlobalAlignment");
  ~TpcGlobalAlignment() override = default;

  int InitRun(PHCompositeNode* topNode) override;

  void setEnabled(bool on) { m_enabled = on; }

  // translation [cm] / [mm]
  void setTranslationCm(double x, double y, double z)
  {
    m_t = Acts::Vector3(x, y, z);
  }
  void setTranslationMm(double x, double y, double z) { setTranslationCm(0.1 * x, 0.1 * y, 0.1 * z); }

  // rotations [rad] about the global x, y, z axes through the pivot, applied as Rz * Ry * Rx
  void setRotation(double rx, double ry, double rz)
  {
    m_rot = Acts::Vector3(rx, ry, rz);
  }

  // rotation centre [cm] (default: origin)
  void setPivotCm(double x, double y, double z) { m_pivot = Acts::Vector3(x, y, z); }

  // Parameters from two beam lines measured in the CURRENT world frame (this module off, TPC rows zero):
  //   TPC-only vertices: x = tpcX0 + tpcDxDz z, y = tpcY0 + tpcDyDz z     [cm, rad]
  //   Si-only  vertices: x = siX0  + siDxDz  z, y = siY0  + siDyDz  z
  // applyTilt = false: translation only.  The pivot is set to the origin, rz = 0, tz = 0.
  void setFromBeamLines(double tpcX0, double tpcY0, double tpcDxDz, double tpcDyDz,
                        double siX0, double siY0, double siDxDz, double siDyDz, bool applyTilt = true);

  // A as an Acts transform (cm)
  Acts::Transform3 transform() const;

 private:
  static Acts::Transform3 envelopeToWorld(const ActsGeometry* geometry);
  void printTests(const ActsGeometry* geometry, const std::string& label) const;

  bool m_enabled = true;
  Acts::Vector3 m_t = Acts::Vector3::Zero();      // cm
  Acts::Vector3 m_rot = Acts::Vector3::Zero();    // rad
  Acts::Vector3 m_pivot = Acts::Vector3::Zero();  // cm

  bool m_applied = false;
  Acts::Transform3 m_lastSet = Acts::Transform3::Identity();  // envelope -> world set by this module
};

#endif
