// ROOT macro to create a PHGarfield Rossegger map from frame charge only.
//
// This launcher intentionally does not use the analytic radial source model,
// pad-placement source geometry, or layer/sector gain maps. The charge density
// is constant inside the frame material and zero elsewhere.
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
    const std::string& garfield_output = "frames_only_2d_rossegger.root",
    const std::string& diagnostic_field_3d_output = "frames_only_fields_3d.root",
    const double frame_reference_phi = 0.0,
    const bool write_diagnostic_field_3d = true,
    const unsigned int tpc_side_for_2d_output = 0,
    const unsigned int radial_job_index = 0,
    const unsigned int n_radial_jobs = 1)
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

  // Constant frame charge density: reference density [nC/m^3], k_eff, alpha.
  // Alpha is ignored by the frame-only charge model.
  rossegger->setDensity(20.0, 5.0, 0.0);
  rossegger->setPhiModulation(0.0, 0.0, 0.0, 0.0);

  // Frame mode replaces the source radial grid with exact frame radial edges.
  // Only source Nphi and Nz are used here.
  //rossegger->setSourceGrid(18, 72, 16);
  //rossegger->setObservationGrid(32, 36, 16);
  rossegger->setSourceGrid(18, 216, 48);
  rossegger->setObservationGrid(96, 108, 48);
  rossegger->setModeTruncation(24, 7, 7);
  rossegger->setAutoAxisymmetric(false);

  rossegger->setTpcSide(tpc_side_for_2d_output);
  rossegger->setRadialJob(radial_job_index, n_radial_jobs);

  // This is the important part: frames only, no previous charge model.
  rossegger->setUseRealTpcSourceGeometry(false);
  rossegger->setUseFrameChargeModel(true);
  rossegger->setFrameReferencePhi(frame_reference_phi);
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
