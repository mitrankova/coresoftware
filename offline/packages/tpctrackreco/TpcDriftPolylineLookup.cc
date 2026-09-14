#include "TpcDriftPolylineLookup.h"

#include "IdealPadMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <phgarfield/PHGarfield.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>

#include <TPolyLine3D.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>

namespace
{
  constexpr double PhiConsistencyTolerance = 1.0e-10;

  constexpr const char* CdbEField2D = "Tpc_PolySeeding_EField";
  constexpr const char* CdbEField3DSide0 = "Tpc_PolySeeding_EField3D_Side0";
  constexpr const char* CdbEField3DSide1 = "Tpc_PolySeeding_EField3D_Side1";
  constexpr const char* CdbModuleFrames3DSide0 = "Tpc_PolySeeding_ModuleFrames3D_Side0";
  constexpr const char* CdbModuleFrames3DSide1 = "Tpc_PolySeeding_ModuleFrames3D_Side1";
  constexpr const char* CdbKEff = "Tpc_PolyClusterizer_kEff";

  double wrap_phi(double phi)
  {
    while (phi > M_PI)
    {
      phi -= 2.0 * M_PI;
    }
    while (phi <= -M_PI)
    {
      phi += 2.0 * M_PI;
    }
    return phi;
  }

  double unwrap_phi_near(double phi, const double reference)
  {
    while (phi - reference > M_PI)
    {
      phi -= 2.0 * M_PI;
    }
    while (phi - reference <= -M_PI)
    {
      phi += 2.0 * M_PI;
    }
    return phi;
  }

  double clamp_unit(const double value)
  {
    return std::max(0.0, std::min(1.0, value));
  }

  double phi_sample_fraction(const unsigned int sample)
  {
    if (TpcDriftPolylineLookup::NPhiSamples <= 1U)
    {
      return 0.0;
    }
    return static_cast<double>(sample) / static_cast<double>(TpcDriftPolylineLookup::NPhiSamples - 1U);
  }

  template <class T>
  bool merge_value(T& destination,
                   bool& destinationOverride,
                   const T& source,
                   const bool sourceOverride,
                   const char* label,
                   const std::string& sourceName)
  {
    if (!sourceOverride)
    {
      return true;
    }

    if (destinationOverride && destination != source)
    {
      std::cerr << "TpcDriftPolylineLookup::registerConfiguration - conflicting explicit " << label
                << " from " << sourceName << std::endl;
      return false;
    }

    destination = source;
    destinationOverride = true;
    return true;
  }
}  // namespace

TpcDriftPolylineLookup::TpcDriftPolylineLookup() = default;
TpcDriftPolylineLookup::~TpcDriftPolylineLookup() = default;

TpcDriftPolylineLookup* TpcDriftPolylineLookup::get(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    return nullptr;
  }

  PHNodeIterator iter(topNode);
  auto* runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    return nullptr;
  }

  return findNode::getClass<TpcDriftPolylineLookup>(runNode, NodeName);
}

TpcDriftPolylineLookup* TpcDriftPolylineLookup::getOrCreate(PHCompositeNode* topNode, const int verbosity)
{
  if (!topNode)
  {
    return nullptr;
  }

  if (auto* lookup = get(topNode))
  {
    return lookup;
  }

  PHNodeIterator iter(topNode);
  auto* runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << "TpcDriftPolylineLookup::getOrCreate - RUN node missing" << std::endl;
    return nullptr;
  }

  auto* lookup = new TpcDriftPolylineLookup();
  runNode->addNode(new PHDataNode<TpcDriftPolylineLookup>(lookup, NodeName));
  if (verbosity > 0)
  {
    std::cout << "TpcDriftPolylineLookup::getOrCreate - created transient RUN/" << NodeName << std::endl;
  }
  return lookup;
}

bool TpcDriftPolylineLookup::registerConfiguration(const Config& config,
                                           const std::string& source,
                                           const int verbosity)
{
  if (m_initialized)
  {
    std::cerr << "TpcDriftPolylineLookup::registerConfiguration - configuration from " << source
              << " arrived after lookup initialization" << std::endl;
    return false;
  }

  bool ok = true;
  ok &= merge_value(m_config.t0, m_config.t0Override, config.t0, config.t0Override, "t0", source);
  ok &= merge_value(m_config.tpcAdcClock, m_config.tpcAdcClockOverride, config.tpcAdcClock, config.tpcAdcClockOverride, "TPC ADC clock", source);
  ok &= merge_value(m_config.crossingPeriodNs, m_config.crossingPeriodNsOverride, config.crossingPeriodNs, config.crossingPeriodNsOverride, "crossing period", source);
  ok &= merge_value(m_config.reverseDriftStepNs, m_config.reverseDriftStepNsOverride, config.reverseDriftStepNs, config.reverseDriftStepNsOverride, "reverse-drift step", source);

  if (config.startZOverride)
  {
    if (m_config.startZOverride &&
        (m_config.startZSouth != config.startZSouth || m_config.startZNorth != config.startZNorth))
    {
      std::cerr << "TpcDriftPolylineLookup::registerConfiguration - conflicting explicit start-z values from " << source << std::endl;
      ok = false;
    }
    else
    {
      m_config.startZSouth = config.startZSouth;
      m_config.startZNorth = config.startZNorth;
      m_config.startZOverride = true;
    }
  }

  ok &= merge_value(m_config.kEffSide0, m_config.kEffSide0Override, config.kEffSide0, config.kEffSide0Override, "kEff side 0", source);
  ok &= merge_value(m_config.kEffSide1, m_config.kEffSide1Override, config.kEffSide1, config.kEffSide1Override, "kEff side 1", source);
  ok &= merge_value(m_config.cmVoltageDefault, m_config.cmVoltageDefaultOverride, config.cmVoltageDefault, config.cmVoltageDefaultOverride, "CM voltage", source);
  ok &= merge_value(m_config.frameChargeScale, m_config.frameChargeScaleOverride, config.frameChargeScale, config.frameChargeScaleOverride, "frame charge scale", source);
  ok &= merge_value(m_config.useSurveyGeometry, m_config.useSurveyGeometryOverride, config.useSurveyGeometry, config.useSurveyGeometryOverride, "survey geometry flag", source);
  ok &= merge_value(m_config.use2DElectricFieldMap, m_config.use2DElectricFieldMapOverride, config.use2DElectricFieldMap, config.use2DElectricFieldMapOverride, "2D field-map flag", source);
  ok &= merge_value(m_config.tpcMove, m_config.tpcMoveOverride, config.tpcMove, config.tpcMoveOverride, "TPC translation", source);
  ok &= merge_value(m_config.fieldCageVoltageOffsets, m_config.fieldCageVoltageOffsetsOverride, config.fieldCageVoltageOffsets, config.fieldCageVoltageOffsetsOverride, "field-cage voltage offsets", source);

  for (unsigned int i = 0; i < m_config.tpcRotations.size(); ++i)
  {
    ok &= merge_value(m_config.tpcRotations[i], m_config.tpcRotationOverride[i], config.tpcRotations[i], config.tpcRotationOverride[i], "TPC rotation", source);
  }

  ok &= merge_value(m_config.electricFieldMap, m_config.electricFieldMapOverride, config.electricFieldMap, config.electricFieldMapOverride, "2D electric-field map", source);
  ok &= merge_value(m_config.field3DCoefficientFile, m_config.field3DCoefficientFileOverride, config.field3DCoefficientFile, config.field3DCoefficientFileOverride, "kEff coefficient file", source);
  ok &= merge_value(m_config.field3DSide0, m_config.field3DSide0Override, config.field3DSide0, config.field3DSide0Override, "3D electric-field map side 0", source);
  ok &= merge_value(m_config.field3DSide1, m_config.field3DSide1Override, config.field3DSide1, config.field3DSide1Override, "3D electric-field map side 1", source);
  ok &= merge_value(m_config.framesSide0, m_config.framesSide0Override, config.framesSide0, config.framesSide0Override, "frame electric-field map side 0", source);
  ok &= merge_value(m_config.framesSide1, m_config.framesSide1Override, config.framesSide1, config.framesSide1Override, "frame electric-field map side 1", source);

  if (verbosity > 1 && ok)
  {
    std::cout << "TpcDriftPolylineLookup::registerConfiguration - accepted configuration from " << source << std::endl;
  }
  return ok;
}

bool TpcDriftPolylineLookup::load_cdb_inputs(const int verbosity)
{
  auto resolve_file = [verbosity](const std::string& payload, std::string& filename, const bool manualOverride) -> bool
  {
    if (manualOverride)
    {
      if (verbosity > 0)
      {
        std::cout << "TpcDriftPolylineLookup::load_cdb_inputs - manual override for " << payload << ": " << filename << std::endl;
      }
      if (filename.empty())
      {
        std::cerr << "TpcDriftPolylineLookup::load_cdb_inputs - manual filename is empty for " << payload << std::endl;
        return false;
      }
      return true;
    }

    filename = CDBInterface::instance()->getUrl(payload);
    if (filename.empty())
    {
      std::cerr << "TpcDriftPolylineLookup::load_cdb_inputs - CDB payload not found: " << payload << std::endl;
      return false;
    }
    if (verbosity > 0)
    {
      std::cout << "TpcDriftPolylineLookup::load_cdb_inputs - loaded " << payload << " from CDB: " << filename << std::endl;
    }
    return true;
  };

  bool ok = true;
  m_electricFieldMap = m_config.electricFieldMap;
  m_field3DSide0 = m_config.field3DSide0;
  m_field3DSide1 = m_config.field3DSide1;
  m_framesSide0 = m_config.framesSide0;
  m_framesSide1 = m_config.framesSide1;

  if (m_config.use2DElectricFieldMap)
  {
    if (!resolve_file(CdbEField2D, m_electricFieldMap, m_config.electricFieldMapOverride))
    {
      ok = false;
    }
  }
  else
  {
    if (!resolve_file(CdbEField3DSide0, m_field3DSide0, m_config.field3DSide0Override))
    {
      ok = false;
    }
    if (!resolve_file(CdbEField3DSide1, m_field3DSide1, m_config.field3DSide1Override))
    {
      ok = false;
    }
  }

  if (!resolve_file(CdbModuleFrames3DSide0, m_framesSide0, m_config.framesSide0Override))
  {
    ok = false;
  }
  if (!resolve_file(CdbModuleFrames3DSide1, m_framesSide1, m_config.framesSide1Override))
  {
    ok = false;
  }

  m_kEffSide0 = m_config.kEffSide0;
  m_kEffSide1 = m_config.kEffSide1;
  if (!m_config.kEffSide0Override || !m_config.kEffSide1Override)
  {
    std::string kefffile;
    if (m_config.field3DCoefficientFileOverride)
    {
      kefffile = m_config.field3DCoefficientFile;
      if (verbosity > 0)
      {
        std::cout << "TpcDriftPolylineLookup::load_cdb_inputs - manual kEff coefficient file: " << kefffile << std::endl;
      }
    }
    else
    {
      kefffile = CDBInterface::instance()->getUrl(CdbKEff);
      if (verbosity > 0)
      {
        std::cout << "TpcDriftPolylineLookup::load_cdb_inputs - kEff coefficient file from CDB: " << kefffile << std::endl;
      }
    }

    if (kefffile.empty())
    {
      std::cerr << "TpcDriftPolylineLookup::load_cdb_inputs - kEff coefficient file is empty" << std::endl;
      ok = false;
    }
    else
    {
      auto keffcdbtree = std::make_unique<CDBTTree>(kefffile);
      keffcdbtree->LoadCalibrations();
      if (!m_config.kEffSide0Override)
      {
        m_kEffSide0 = keffcdbtree->GetSingleFloatValue("keffside0");
      }
      if (!m_config.kEffSide1Override)
      {
        m_kEffSide1 = keffcdbtree->GetSingleFloatValue("keffside1");
      }
    }
  }

  if (verbosity > 0)
  {
    std::cout << "TpcDriftPolylineLookup::load_cdb_inputs - final kEff values: side0 = " << m_kEffSide0
              << ", side1 = " << m_kEffSide1 << std::endl;
  }
  return ok;
}

bool TpcDriftPolylineLookup::initialize(PHCompositeNode* topNode, const int verbosity)
{
  if (m_initialized)
  {
    return true;
  }
  if (!topNode)
  {
    return false;
  }

  auto* geomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!geomContainer)
  {
    std::cerr << "TpcDriftPolylineLookup::initialize - missing TPCGEOMCONTAINER" << std::endl;
    return false;
  }

  auto* layergeom = geomContainer->GetLayerCellGeom(20);
  if (!layergeom)
  {
    std::cerr << "TpcDriftPolylineLookup::initialize - missing TPC layer-20 geometry" << std::endl;
    return false;
  }

  m_idealPadMap = std::make_unique<IdealPadMap>();
  if (m_idealPadMap->load_from_cdb(verbosity) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << "TpcDriftPolylineLookup::initialize - failed to load IdealPadMap" << std::endl;
    return false;
  }

  m_tpcMove = m_config.tpcMove;
  m_tpcRotations = m_config.tpcRotations;
  if (m_config.useSurveyGeometry)
  {
    m_tpcMove = {{layergeom->get_place_x(), layergeom->get_place_y(), layergeom->get_place_z()}};
    m_tpcRotations = {{{{layergeom->get_rot_x(), layergeom->get_rot_y(), layergeom->get_rot_z()}}, {{0.0, 0.0, 0.0}}}};
  }

  if (m_config.startZOverride)
  {
    m_startZSouth = m_config.startZSouth;
    m_startZNorth = m_config.startZNorth;
  }
  else
  {
    m_startZSouth = -(layergeom->get_max_driftlength() + layergeom->get_CM_halfwidth());
    m_startZNorth = layergeom->get_max_driftlength() + layergeom->get_CM_halfwidth();
  }

  if (verbosity > 0)
  {
    std::cout << "TpcDriftPolylineLookup::initialize - start z south=" << m_startZSouth
              << " cm north=" << m_startZNorth << " cm" << std::endl;
  }

  if (!load_cdb_inputs(verbosity))
  {
    return false;
  }

  m_garfield = std::make_unique<PHGarfield>("TpcDriftPolylineLookup_PHGarfield", "", m_kEffSide0, m_kEffSide1);
  configure_garfield(m_garfield.get());
  if (m_garfield->InitRun(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "TpcDriftPolylineLookup::initialize - PHGarfield InitRun failed" << std::endl;
    return false;
  }

  if (!build_drift_lookup(verbosity))
  {
    return false;
  }

  m_initialized = true;
  return true;
}

void TpcDriftPolylineLookup::configure_garfield(PHGarfield* garfield) const
{
  if (!garfield)
  {
    return;
  }

  if (m_config.use2DElectricFieldMap)
  {
    garfield->SetElectricFieldMap(m_electricFieldMap);
  }
  else
  {
    garfield->SetElectricFieldMap3D(m_field3DSide0, m_field3DSide1);
  }

  garfield->SetFrameElectricFieldMap3D(m_framesSide0, m_framesSide1);
  garfield->SetFrameChargeScale(m_config.frameChargeScale);
  garfield->SetUseIFCVoltageDistortion(true);
  garfield->SetUseOFCVoltageDistortion(true);
  garfield->SetFieldCageVoltageOffsets(m_config.fieldCageVoltageOffsets[0],
                                       m_config.fieldCageVoltageOffsets[1],
                                       m_config.fieldCageVoltageOffsets[2],
                                       m_config.fieldCageVoltageOffsets[3]);
  garfield->MoveTpc(m_tpcMove[0], m_tpcMove[1], m_tpcMove[2]);
  for (const auto& rotation : m_tpcRotations)
  {
    garfield->RotateTpc(rotation[0], rotation[1], rotation[2]);
  }
  garfield->SetCMVoltageDefault(m_config.cmVoltageDefault);
}

unsigned int TpcDriftPolylineLookup::drift_lookup_index(const unsigned int layer_index,
                                                const unsigned int side,
                                                const unsigned int sector,
                                                const unsigned int sample)
{
  return (((layer_index * NSides + side) * NSectors + sector) * NPhiSamples + sample);
}

const TpcDriftPolylineLookup::DriftPoint& TpcDriftPolylineLookup::driftPoint(const DriftPolyline& polyline,
                                                             const std::size_t index) const
{
  return m_driftPoints[static_cast<std::size_t>(polyline.offset) + index];
}

bool TpcDriftPolylineLookup::build_drift_lookup(const int verbosity)
{
  if (!m_idealPadMap || !m_garfield)
  {
    return false;
  }
  if (m_config.reverseDriftStepNs <= 0.0 || !std::isfinite(m_config.reverseDriftStepNs))
  {
    return false;
  }

  for (auto& polyline : m_driftLookup)
  {
    polyline = DriftPolyline{};
  }
  m_driftPoints.clear();
  m_maxLookupTimeNs = std::numeric_limits<double>::max();

  unsigned int nbuilt = 0;
  for (unsigned int layer = FirstLayer; layer <= LastLayer; ++layer)
  {
    const unsigned int layer_index = layer - FirstLayer;
    const double radius = m_idealPadMap->get_radius(layer);
    const unsigned int pads_per_sector = m_idealPadMap->get_pads_per_sector_for_layer(layer);
    if (!std::isfinite(radius) || pads_per_sector == 0U)
    {
      std::cerr << "TpcDriftPolylineLookup::build_drift_lookup - invalid geometry for layer " << layer << std::endl;
      return false;
    }

    for (unsigned int side = 0; side < NSides; ++side)
    {
      const double z0 = (side == 0U) ? m_startZSouth : m_startZNorth;
      for (unsigned int sector = 0; sector < NSectors; ++sector)
      {
        for (unsigned int sample = 0; sample < NPhiSamples; ++sample)
        {
          const unsigned int local_phibin = static_cast<unsigned int>(std::llround(
              phi_sample_fraction(sample) * static_cast<double>(pads_per_sector - 1U)));
          const unsigned int global_pad = sector * pads_per_sector + local_phibin;
          const double phi_local = m_idealPadMap->get_phi(side, sector, layer, local_phibin);
          const double phi_global = m_idealPadMap->get_phi(side, layer, global_pad);
          if (!std::isfinite(phi_local) || !std::isfinite(phi_global))
          {
            std::cerr << "TpcDriftPolylineLookup::build_drift_lookup - invalid phi for layer " << layer
                      << " side " << side << " sector " << sector << " sample " << sample << std::endl;
            return false;
          }

          if (std::abs(wrap_phi(phi_global - phi_local)) > PhiConsistencyTolerance)
          {
            std::cerr << "TpcDriftPolylineLookup::build_drift_lookup - inconsistent IdealPadMap phi overloads"
                      << " layer " << layer << " side " << side << " sector " << sector
                      << " local_phibin " << local_phibin << " global_pad " << global_pad << std::endl;
            return false;
          }

          const double x0 = radius * std::cos(phi_local);
          const double y0 = radius * std::sin(phi_local);
          std::unique_ptr<TPolyLine3D> drift(m_garfield->ReverseDrift(x0, y0, z0, m_config.reverseDriftStepNs));
          if (!drift || drift->GetN() <= 0)
          {
            std::cerr << "TpcDriftPolylineLookup::build_drift_lookup - ReverseDrift failed for layer " << layer
                      << " side " << side << " sector " << sector << " sample " << sample << std::endl;
            return false;
          }

          const int npoints = drift->GetN();
          const Float_t* xyz = drift->GetP();
          if (!xyz || npoints <= 0)
          {
            std::cerr << "TpcDriftPolylineLookup::build_drift_lookup - empty drift points for layer " << layer
                      << " side " << side << " sector " << sector << " sample " << sample << std::endl;
            return false;
          }

          if (m_driftPoints.capacity() == 0U)
          {
            m_driftPoints.reserve(static_cast<std::size_t>(NPolylines) * static_cast<std::size_t>(npoints));
          }

          DriftPolyline& polyline = m_driftLookup[drift_lookup_index(layer_index, side, sector, sample)];
          if (m_driftPoints.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
          {
            std::cerr << "TpcDriftPolylineLookup::build_drift_lookup - compact lookup offset overflow" << std::endl;
            return false;
          }
          polyline.phi = phi_local;
          polyline.offset = static_cast<std::uint32_t>(m_driftPoints.size());
          polyline.count = static_cast<std::uint32_t>(npoints);

          for (int ipoint = 0; ipoint < npoints; ++ipoint)
          {
            const int idx = 3 * ipoint;
            const double output_r = std::hypot(static_cast<double>(xyz[idx]), static_cast<double>(xyz[idx + 1]));
            const double output_phi = unwrap_phi_near(
                std::atan2(static_cast<double>(xyz[idx + 1]), static_cast<double>(xyz[idx])), phi_local);
            DriftPoint point;
            point.delta_r = static_cast<float>(output_r - radius);
            point.delta_phi = static_cast<float>(output_phi - phi_local);
            point.z = xyz[idx + 2];
            m_driftPoints.push_back(point);
          }

          m_maxLookupTimeNs = std::min(
              m_maxLookupTimeNs,
              static_cast<double>(npoints - 1) * m_config.reverseDriftStepNs);
          ++nbuilt;
        }
      }
    }
  }

  // Preserve the clusterizer's existing phi-symmetry QA.
  for (unsigned int layer = FirstLayer; layer <= LastLayer; ++layer)
  {
    const unsigned int layer_index = layer - FirstLayer;
    const double radius = m_idealPadMap->get_radius(layer);
    for (unsigned int side = 0; side < NSides; ++side)
    {
      const DriftPolyline& reference = m_driftLookup[drift_lookup_index(layer_index, side, 0, 0)];
      bool symmetry_ok = reference.count > 0U && std::isfinite(radius);
      unsigned int length_mismatches = 0;
      double max_delta_r = 0.0;
      double max_delta_phi_arc = 0.0;
      double max_delta_z = 0.0;

      for (unsigned int sector = 0; sector < NSectors; ++sector)
      {
        for (unsigned int sample = 0; sample < NPhiSamples; ++sample)
        {
          const DriftPolyline& polyline = m_driftLookup[drift_lookup_index(layer_index, side, sector, sample)];
          if (polyline.count != reference.count)
          {
            ++length_mismatches;
            symmetry_ok = false;
          }

          const std::size_t npoints = std::min<std::size_t>(reference.count, polyline.count);
          for (std::size_t ipoint = 0; ipoint < npoints; ++ipoint)
          {
            const DriftPoint& ref_point = driftPoint(reference, ipoint);
            const DriftPoint& point = driftPoint(polyline, ipoint);
            max_delta_r = std::max(max_delta_r, std::abs(static_cast<double>(point.delta_r - ref_point.delta_r)));
            max_delta_phi_arc = std::max(
                max_delta_phi_arc,
                std::abs(radius * wrap_phi(static_cast<double>(point.delta_phi - ref_point.delta_phi))));
            max_delta_z = std::max(max_delta_z, std::abs(static_cast<double>(point.z - ref_point.z)));
          }
        }
      }

      if (max_delta_r > 0.001 || max_delta_phi_arc > 0.001 || max_delta_z > 0.001)
      {
        symmetry_ok = false;
      }
      if (verbosity > 0)
      {
        std::cout << "TpcDriftPolylineLookup::build_drift_lookup - phi symmetry "
                  << (symmetry_ok ? "holds" : "broken")
                  << " layer=" << layer << " side=" << side
                  << " max_dr=" << max_delta_r
                  << " max_r_dphi=" << max_delta_phi_arc
                  << " max_dz=" << max_delta_z
                  << " length_mismatches=" << length_mismatches << std::endl;
      }
    }
  }

  if (verbosity > 0)
  {
    std::cout << "TpcDriftPolylineLookup::build_drift_lookup - built " << nbuilt
              << " drift polylines, points=" << m_driftPoints.size()
              << " max_lookup_time_ns=" << m_maxLookupTimeNs << std::endl;
  }

  return nbuilt == NPolylines && std::isfinite(m_maxLookupTimeNs) && m_maxLookupTimeNs > 0.0;
}

double TpcDriftPolylineLookup::crossingToDriftTimeNs(const unsigned int tbin, const short crossing) const
{
  return (static_cast<double>(tbin) - m_config.t0) * m_config.tpcAdcClock
         - static_cast<double>(crossing) * m_config.crossingPeriodNs;
}

bool TpcDriftPolylineLookup::sample_drift_lookup(const unsigned int layer,
                                         const unsigned int side,
                                         const unsigned int pad,
                                         const unsigned int tbin,
                                         const short crossing,
                                         Point& point) const
{
  if (!m_initialized || !m_idealPadMap)
  {
    return false;
  }
  if (layer < FirstLayer || layer > LastLayer || side >= NSides)
  {
    return false;
  }
  if (m_config.reverseDriftStepNs <= 0.0 || !std::isfinite(m_config.reverseDriftStepNs))
  {
    return false;
  }

  const unsigned int pads_per_sector = m_idealPadMap->get_pads_per_sector_for_layer(layer);
  if (pads_per_sector == 0U)
  {
    return false;
  }
  const unsigned int sector = pad / pads_per_sector;
  if (sector >= NSectors)
  {
    return false;
  }

  const double hit_radius = m_idealPadMap->get_radius(layer);
  const double hit_phi = m_idealPadMap->get_phi(side, layer, pad);
  if (!std::isfinite(hit_radius) || !std::isfinite(hit_phi))
  {
    return false;
  }

  const double target_time_ns = crossingToDriftTimeNs(tbin, crossing);
  if (target_time_ns <= 0.0 || !std::isfinite(target_time_ns))
  {
    return false;
  }

  const unsigned int layer_index = layer - FirstLayer;
  std::array<const DriftPolyline*, NPhiSamples> samples{};
  std::array<double, NPhiSamples> sample_phi{};
  for (unsigned int sample = 0; sample < NPhiSamples; ++sample)
  {
    samples[sample] = &m_driftLookup[drift_lookup_index(layer_index, side, sector, sample)];
    if (samples[sample]->count == 0U)
    {
      return false;
    }
    sample_phi[sample] = samples[sample]->phi;
    if (sample > 0U)
    {
      sample_phi[sample] = unwrap_phi_near(sample_phi[sample], sample_phi[sample - 1U]);
    }
  }

  const double unwrapped_hit_phi = unwrap_phi_near(hit_phi, sample_phi[NPhiSamples / 2U]);
  const bool increasing = sample_phi[NPhiSamples - 1U] >= sample_phi[0];

  bool bracket_found = false;
  unsigned int sample0 = 0;
  unsigned int sample1 = 0;
  double phi_fraction = 0.0;
  if ((increasing && unwrapped_hit_phi <= sample_phi[0]) ||
      (!increasing && unwrapped_hit_phi >= sample_phi[0]))
  {
    bracket_found = true;
  }
  else if ((increasing && unwrapped_hit_phi >= sample_phi[NPhiSamples - 1U]) ||
           (!increasing && unwrapped_hit_phi <= sample_phi[NPhiSamples - 1U]))
  {
    sample0 = NPhiSamples - 1U;
    sample1 = NPhiSamples - 1U;
    bracket_found = true;
  }
  else
  {
    for (unsigned int sample = 0; sample + 1U < NPhiSamples; ++sample)
    {
      const bool in_interval = increasing
          ? (unwrapped_hit_phi >= sample_phi[sample] && unwrapped_hit_phi <= sample_phi[sample + 1U])
          : (unwrapped_hit_phi <= sample_phi[sample] && unwrapped_hit_phi >= sample_phi[sample + 1U]);
      if (!in_interval)
      {
        continue;
      }

      sample0 = sample;
      sample1 = sample + 1U;
      const double denom = sample_phi[sample1] - sample_phi[sample0];
      phi_fraction = (denom != 0.0)
          ? clamp_unit((unwrapped_hit_phi - sample_phi[sample0]) / denom)
          : 0.0;
      bracket_found = true;
      break;
    }
  }
  if (!bracket_found)
  {
    return false;
  }

  auto sample_time = [this, target_time_ns](const DriftPolyline& polyline,
                                            double& delta_r,
                                            double& delta_phi,
                                            double& point_z) -> bool
  {
    const int npoints = static_cast<int>(polyline.count);
    if (npoints <= 0)
    {
      return false;
    }

    const double max_time_ns = static_cast<double>(npoints - 1) * m_config.reverseDriftStepNs;
    if (target_time_ns > max_time_ns)
    {
      return false;
    }

    const double fbin = target_time_ns / m_config.reverseDriftStepNs;
    const int i0 = std::min(static_cast<int>(std::floor(fbin)), npoints - 1);
    const int i1 = std::min(i0 + 1, npoints - 1);
    const double frac = fbin - static_cast<double>(i0);

    const DriftPoint& p0 = driftPoint(polyline, static_cast<std::size_t>(i0));
    const DriftPoint& p1 = driftPoint(polyline, static_cast<std::size_t>(i1));
    const double dphi0 = static_cast<double>(p0.delta_phi);
    const double dphi1 = unwrap_phi_near(static_cast<double>(p1.delta_phi), dphi0);
    delta_r = static_cast<double>(p0.delta_r) + frac * static_cast<double>(p1.delta_r - p0.delta_r);
    delta_phi = dphi0 + frac * (dphi1 - dphi0);
    point_z = static_cast<double>(p0.z) + frac * static_cast<double>(p1.z - p0.z);
    return std::isfinite(delta_r) && std::isfinite(delta_phi) && std::isfinite(point_z);
  };

  double delta_r0 = 0.0;
  double delta_phi0 = 0.0;
  double z0 = 0.0;
  const bool valid0 = sample_time(*samples[sample0], delta_r0, delta_phi0, z0);

  double delta_r1 = 0.0;
  double delta_phi1 = 0.0;
  double z1 = 0.0;
  const bool same_sample = sample0 == sample1;
  const bool valid1 = same_sample ? valid0 : sample_time(*samples[sample1], delta_r1, delta_phi1, z1);
  if (same_sample)
  {
    delta_r1 = delta_r0;
    delta_phi1 = delta_phi0;
    z1 = z0;
  }
  if (!valid0 && !valid1)
  {
    return false;
  }

  double delta_r = 0.0;
  double delta_phi = 0.0;
  double z = 0.0;
  if (valid0 && valid1 && !same_sample)
  {
    delta_r = delta_r0 + phi_fraction * (delta_r1 - delta_r0);
    const double unwrapped_delta_phi1 = unwrap_phi_near(delta_phi1, delta_phi0);
    delta_phi = delta_phi0 + phi_fraction * (unwrapped_delta_phi1 - delta_phi0);
    z = z0 + phi_fraction * (z1 - z0);
  }
  else if (valid0)
  {
    delta_r = delta_r0;
    delta_phi = delta_phi0;
    z = z0;
  }
  else
  {
    delta_r = delta_r1;
    delta_phi = delta_phi1;
    z = z1;
  }

  const double radius = hit_radius + delta_r;
  const double output_phi = unwrapped_hit_phi + delta_phi;
  point.x = radius * std::cos(output_phi);
  point.y = radius * std::sin(output_phi);
  point.z = z;
  return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z);
}

bool TpcDriftPolylineLookup::getPosition(const unsigned int layer,
                                 const unsigned int side,
                                 const unsigned int pad,
                                 const unsigned int tbin,
                                 const short crossing,
                                 Point& point) const
{
  return sample_drift_lookup(layer, side, pad, tbin, crossing, point);
}

bool TpcDriftPolylineLookup::getPosition(const TrkrDefs::hitsetkey hitsetkey,
                                 const TrkrDefs::hitkey hitkey,
                                 const short crossing,
                                 Point& point) const
{
  return getPosition(TrkrDefs::getLayer(hitsetkey),
                     TpcDefs::getSide(hitsetkey),
                     TpcDefs::getPad(hitkey),
                     TpcDefs::getTBin(hitkey),
                     crossing,
                     point);
}

bool TpcDriftPolylineLookup::getDelta(const unsigned int layer,
                              const unsigned int side,
                              const unsigned int pad,
                              const unsigned int tbin,
                              const short crossing,
                              const short referenceCrossing,
                              Point& delta) const
{
  Point position;
  Point reference;
  if (!getPosition(layer, side, pad, tbin, crossing, position) ||
      !getPosition(layer, side, pad, tbin, referenceCrossing, reference))
  {
    return false;
  }

  delta.x = position.x - reference.x;
  delta.y = position.y - reference.y;
  delta.z = position.z - reference.z;
  return true;
}

bool TpcDriftPolylineLookup::getDelta(const TrkrDefs::hitsetkey hitsetkey,
                              const TrkrDefs::hitkey hitkey,
                              const short crossing,
                              const short referenceCrossing,
                              Point& delta) const
{
  return getDelta(TrkrDefs::getLayer(hitsetkey),
                  TpcDefs::getSide(hitsetkey),
                  TpcDefs::getPad(hitkey),
                  TpcDefs::getTBin(hitkey),
                  crossing,
                  referenceCrossing,
                  delta);
}

bool TpcDriftPolylineLookup::hasEntry(const unsigned int layer,
                              const unsigned int side,
                              const unsigned int pad,
                              const unsigned int tbin,
                              const short crossing) const
{
  Point point;
  return getPosition(layer, side, pad, tbin, crossing, point);
}
