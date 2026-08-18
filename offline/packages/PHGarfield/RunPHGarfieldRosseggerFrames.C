// ROOT macro to create a PHGarfield Rossegger map from charged frame volumes.
//
// Frame charge uses the ordinary Rossegger volume-source boundary conditions,
// with each frame protruding into z by its own thickness.
//
// Example:
//   root -b -q 'RunPHGarfieldRosseggerFrames.C()'

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
    const std::string& garfield_output = "frames_side1_2d_rossegger_geom_v14.root",
    const std::string& diagnostic_field_3d_output = "frames_side1_3d_geom_v14.root",
    const double frame_reference_phi = 0.0,
    const bool write_diagnostic_field_3d = true,
    const unsigned int tpc_side_for_2d_output = 1,
    const unsigned int mode_phi_max = 24,
    const unsigned int n_radial_modes = 15,
    const unsigned int n_longitudinal_modes = 2,
    const unsigned int radial_job_index = 0,
    const unsigned int n_radial_jobs = 1,
    const std::string& frame_geometry_file = "tpc_frame_geometry.csv")
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

  rossegger->setSourceGrid(18, 144, 24);
  rossegger->setObservationGrid(192, 240, 24);
  rossegger->setModeTruncation(mode_phi_max, n_radial_modes, n_longitudinal_modes);

  //rossegger->setSourceGrid(18, 288, 24);
  //rossegger->setObservationGrid(48, 180, 24);
 // rossegger->setModeTruncation(48, 48, 24);
  rossegger->setAutoAxisymmetric(false);

  rossegger->setTpcSide(tpc_side_for_2d_output);
  rossegger->setRadialJob(radial_job_index, n_radial_jobs);

  // Frames only: charged frame volume, no previous IBF charge model.
  rossegger->setUseRealTpcSourceGeometry(false);
  rossegger->setUseFrameChargeModel(true);
  rossegger->setFrameReferencePhi(frame_reference_phi);
  rossegger->setFrameGeometryFile(frame_geometry_file);
  rossegger->setFrameBoundaryPotential(1.0);  // frame source-strength scale
  rossegger->setFrameChargeWeighting(
      PHGarfieldRossegger::FrameChargeWeighting::ProportionalToArea);
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
