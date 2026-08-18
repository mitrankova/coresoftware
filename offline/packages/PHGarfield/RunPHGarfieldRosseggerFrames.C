// ROOT macro to create a PHGarfield Rossegger map from insulating-frame surface charge.
//
// Frame charge is represented as a thin sheet on the exposed frame polygons.
// The Rossegger solution supplies the conducting TPC field-cage/endcap boundaries.
//
// Examples:
//   root -b -q 'RunPHGarfieldRosseggerFrames.C()'       // full frame
//   root -b -q 'RunPHGarfieldRosseggerFrames.C("radial.root", "radial_3d.root", 0.0, true, 1, 72, 40, 2, 0, 1, "tpc_frame_geometry.csv", 200000.0, 0.01, 1, 1.0, 0.0)'
//   root -b -q 'RunPHGarfieldRosseggerFrames.C("side.root", "side_3d.root", 0.0, true, 1, 72, 40, 2, 0, 1, "tpc_frame_geometry.csv", 200000.0, 0.01, 2, 1.0, 1.0)'

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libPHGarfield.so)

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phgarfield/PHGarfieldRossegger.h>

#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

int RunPHGarfieldRosseggerFrames(
    const std::string& garfield_output = "frames_side0_2d_rossegger_geom_v17_rad_only.root",
    const std::string& diagnostic_field_3d_output = "frames_side0_3d_geom_v17_rad_only.root",
    const double frame_reference_phi = 0.0,
    const bool write_diagnostic_field_3d = true,
    const unsigned int tpc_side_for_2d_output = 0,
    const unsigned int mode_phi_max = 72,
    const unsigned int n_radial_modes = 40,
    const unsigned int n_longitudinal_modes = 2,
    const unsigned int radial_job_index = 0,
    const unsigned int n_radial_jobs = 1,
    const std::string& frame_geometry_file = "tpc_frame_geometry.csv",
    const double frame_surface_charge_density_nC_per_m2 = 200000.0,
    const double frame_sheet_thickness_cm = 0.01,
    const unsigned int frame_charge_pieces = 0,
    const double radial_rail_weight = 1.0,
    const double side_rail_weight = 0.0)
{
  if (fs::exists(garfield_output) &&
      (!write_diagnostic_field_3d || fs::exists(diagnostic_field_3d_output)))
  {
    std::cout << "Requested frame-only Rossegger outputs already exist. Skipping calculation.\n"
              << "  " << garfield_output << "\n";
    if (write_diagnostic_field_3d) { std::cout << "  " << diagnostic_field_3d_output << "\n"; }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto* se = Fun4AllServer::instance();
  auto* rossegger = new PHGarfieldRossegger("PHGarfieldRosseggerFrames");

  rossegger->setGeometryCm(21.78, 76.28, 102.325);
  rossegger->setSourceRadiusCm(21.78, 76.28);

  rossegger->setPhiModulation(0.0, 0.0, 0.0, 0.0);

  // Frame mode keeps exact frame source radial edges; source Nphi and observation bins remain the active resolution knobs.
//  rossegger->setSourceGrid(18, 72, 16);
//  rossegger->setObservationGrid(32, 36, 16);
 // rossegger->setModeTruncation(24, 12, 12);

  rossegger->setSourceGrid(18, 360, 1);
  rossegger->setObservationGrid(192, 360, 24);
  rossegger->setModeTruncation(mode_phi_max, n_radial_modes, n_longitudinal_modes);

  //rossegger->setSourceGrid(18, 288, 24);
  //rossegger->setObservationGrid(48, 180, 24);
 // rossegger->setModeTruncation(48, 48, 24);
  rossegger->setAutoAxisymmetric(false);

  rossegger->setTpcSide(tpc_side_for_2d_output);
  rossegger->setRadialJob(radial_job_index, n_radial_jobs);

  // Frames only: insulating surface charge, no previous IBF charge model.
  rossegger->setUseRealTpcSourceGeometry(false);
  rossegger->setUseFrameChargeModel(true);
  rossegger->setFrameReferencePhi(frame_reference_phi);
  rossegger->setFrameGeometryFile(frame_geometry_file);
  rossegger->setFrameBoundaryPotential(1.0);  // dimensionless source-strength scale
  rossegger->setFrameSurfaceChargeDensity(frame_surface_charge_density_nC_per_m2);
  rossegger->setFrameSheetThicknessCm(frame_sheet_thickness_cm);
  rossegger->setFrameChargeWeighting(
      PHGarfieldRossegger::FrameChargeWeighting::ProportionalToArea);
  rossegger->setFrameChargePieces(
      static_cast<PHGarfieldRossegger::FrameChargePieces>(frame_charge_pieces));
  rossegger->setFrameRailWeights(radial_rail_weight, side_rail_weight);
  rossegger->setDivideChargeByGain(false);
  rossegger->setNormalizeGainWeightedTotal(false);

  rossegger->setGarfieldOutputFile(garfield_output);
  rossegger->setWriteField3D(write_diagnostic_field_3d);
  rossegger->setField3DOutputFile(diagnostic_field_3d_output);
  rossegger->setWritePHGarfieldField3D(false);
  rossegger->setVerifyOutput(true);

  se->registerSubsystem(rossegger);

  const int status = se->run(1);
  se->End();
  return status;
}
