// ROOT macro to create PHGarfield Rossegger electric-field maps.
//
// Example, from a built sPHENIX environment:
//   root -b -q 'RunPHGarfieldRossegger.C()'
//
// Real TPC source geometry mode expects the pad-placement text file and the
// layer/sector gain ROOT file to exist at the paths passed below.

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libPHGarfield.so)

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phgarfield/PHGarfieldRossegger.h>

#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

int RunPHGarfieldRossegger(
    const std::string& garfield_output = "simple_2d_rossegger.root",
    const bool use_real_tpc_source_geometry = true,
    const std::string& pad_placement_file = "input/TPC_pad_placement.txt",
    const std::string& gain_map_file = "input/layer_gain_79513_Mariia_side01.root",
    const bool write_diagnostic_field_3d = true,
    const std::string& diagnostic_field_3d_output = "sphenix_rossegger_fields_3d.root",
    const bool write_phgarfield_field_3d = true,
    const std::string& phgarfield_field_3d_side0 = "sphenix_3d_ibf_field_side0_South_v1.root",
    const std::string& phgarfield_field_3d_side1 = "sphenix_3d_ibf_field_side1_North_v1.root",
    const unsigned int tpc_side_for_2d_output = 0,
    const unsigned int radial_job_index = 0,
    const unsigned int n_radial_jobs = 1)
{
  if (fs::exists(garfield_output) &&
      (!write_diagnostic_field_3d || fs::exists(diagnostic_field_3d_output)) &&
      (!write_phgarfield_field_3d ||
       (fs::exists(phgarfield_field_3d_side0) && fs::exists(phgarfield_field_3d_side1))))
  {
    std::cout << "Requested Rossegger map outputs already exist. Skipping calculation.\n"
              << "  " << garfield_output << "\n";
    if (write_diagnostic_field_3d) { std::cout << "  " << diagnostic_field_3d_output << "\n"; }
    if (write_phgarfield_field_3d)
    {
      std::cout << "  " << phgarfield_field_3d_side0 << "\n"
                << "  " << phgarfield_field_3d_side1 << "\n";
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto* se = Fun4AllServer::instance();
  auto* rossegger = new PHGarfieldRossegger("PHGarfieldRossegger");

  // Geometry in cm. These match the newer direct-Green-function notebook.
  rossegger->setGeometryCm(21.78, 76.28, 102.325);
  rossegger->setSourceRadiusCm(22.8, 75.43);

  // Charge model: reference density [nC/m^3], k_eff, radial power alpha.
  rossegger->setDensity(20.0, 1.0, 1.22);
  rossegger->setPhiModulation(0.0, 0.0, 0.0, 0.0);

  // Notebook development grid. Increase these for production-quality maps.
  rossegger->setSourceGrid(18, 72, 16);
  rossegger->setObservationGrid(32, 36, 16);
  rossegger->setModeTruncation(24, 7, 7);
  rossegger->setAutoAxisymmetric(false);

  rossegger->setTpcSide(tpc_side_for_2d_output);
  rossegger->setRadialJob(radial_job_index, n_radial_jobs);

  rossegger->setUseRealTpcSourceGeometry(use_real_tpc_source_geometry);
  rossegger->setPadPlacementFile(pad_placement_file);
  rossegger->setGainMapFile(gain_map_file);
  rossegger->setGainHistograms("hGainMap_side0_South", "hGainMap_side1_North");
  rossegger->setDivideChargeByGain(true);
  rossegger->setNormalizeGainWeightedTotal(true);

  rossegger->setGarfieldOutputFile(garfield_output);
  rossegger->setWriteField3D(write_diagnostic_field_3d);
  rossegger->setField3DOutputFile(diagnostic_field_3d_output);
  rossegger->setWritePHGarfieldField3D(write_phgarfield_field_3d);
  rossegger->setPHGarfieldField3DOutputFiles(phgarfield_field_3d_side0, phgarfield_field_3d_side1);
  rossegger->setVerifyOutput(true);

  se->registerSubsystem(rossegger);

  const int status = se->run(1);
  se->End();
  return status;
}
