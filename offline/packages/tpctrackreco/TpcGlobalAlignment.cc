#include "TpcGlobalAlignment.h"

#include <trackbase/ActsGeometry.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <Eigen/Geometry>

#include <cmath>
#include <iomanip>
#include <iostream>

TpcGlobalAlignment::TpcGlobalAlignment(const std::string& name)
  : SubsysReco(name)
{
}

void TpcGlobalAlignment::setFromBeamLines(double tpcX0, double tpcY0, double tpcDxDz, double tpcDyDz,
                                          double siX0, double siY0, double siDxDz, double siDyDz, bool applyTilt)
{
  m_pivot = Acts::Vector3::Zero();
  m_t = Acts::Vector3(siX0 - tpcX0, siY0 - tpcY0, 0.0);
  if (applyTilt)
  {
    // x' = x + ry z  -> slope b + ry = d ;  y' = y - rx z -> slope b' - rx = d'
    m_rot = Acts::Vector3(tpcDyDz - siDyDz, siDxDz - tpcDxDz, 0.0);
  }
  else
  {
    m_rot = Acts::Vector3::Zero();
  }
}

Acts::Transform3 TpcGlobalAlignment::transform() const
{
  const Eigen::AngleAxisd ax(m_rot.x(), Eigen::Vector3d::UnitX());
  const Eigen::AngleAxisd ay(m_rot.y(), Eigen::Vector3d::UnitY());
  const Eigen::AngleAxisd az(m_rot.z(), Eigen::Vector3d::UnitZ());
  const Eigen::Quaternion<double> q = az * ay * ax;

  Acts::Transform3 a = Acts::Transform3::Identity();
  a.linear() = q.matrix();
  a.translation() = m_t + m_pivot - q.matrix() * m_pivot;
  return a;
}

// ActsGeometry exposes only the world -> envelope setter; recover E (envelope -> world) exactly from
// the public mapping (it is affine)
Acts::Transform3 TpcGlobalAlignment::envelopeToWorld(const ActsGeometry* geometry)
{
  const Acts::Vector3 o = geometry->transformTpcEnvelopeToWorld(Acts::Vector3(0, 0, 0));
  Eigen::Matrix3d r;
  for (int c = 0; c < 3; ++c)
  {
    Acts::Vector3 e = Acts::Vector3::Zero();
    e(c) = 1.0;
    r.col(c) = geometry->transformTpcEnvelopeToWorld(e) - o;
  }
  Acts::Transform3 t = Acts::Transform3::Identity();
  t.linear() = r;
  t.translation() = o;
  return t;
}

void TpcGlobalAlignment::printTests(const ActsGeometry* geometry, const std::string& label) const
{
  std::cout << Name() << ": " << label << " envelope -> world (mm):";
  for (double z : {113.025, 0.0, -113.025})
  {
    const Acts::Vector3 w = geometry->transformTpcEnvelopeToWorld(Acts::Vector3(0, 0, z)) * 10.0;
    std::cout << "  (0,0," << 10 * z << ") -> (" << std::fixed << std::setprecision(3) << w.x() << ", " << w.y() << ", "
              << w.z() << ")";
  }
  std::cout << std::defaultfloat << std::endl;
}

int TpcGlobalAlignment::InitRun(PHCompositeNode* topNode)
{
  if (!m_enabled)
  {
    std::cout << Name() << ": disabled - TPC envelope -> world transform unchanged" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto* geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!geometry)
  {
    std::cerr << Name() << ": no ActsGeometry node - register this module after TrackingInit()" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const Acts::Transform3 current = envelopeToWorld(geometry);
  if (m_applied && current.isApprox(m_lastSet, 1e-12))
  {
    // new run, same geometry object: already applied, do not apply twice
    std::cout << Name() << ": already applied to this ActsGeometry" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  printTests(geometry, "before");

  const Acts::Transform3 a = transform();
  const Acts::Transform3 updated = a * current;
  geometry->set_tpc_world_envelope_transform(updated.inverse());
  m_lastSet = updated;
  m_applied = true;

  std::cout << Name() << ": TPC clusters moved by A = T * Rz*Ry*Rx about pivot (" << m_pivot.x() << ", " << m_pivot.y() << ", "
            << m_pivot.z() << ") cm:  t = (" << 10 * m_t.x() << ", " << 10 * m_t.y() << ", " << 10 * m_t.z() << ") mm,  rx = "
            << 1e3 * m_rot.x() << " ry = " << 1e3 * m_rot.y() << " rz = " << 1e3 * m_rot.z()
            << " mrad (ACTS surfaces unchanged)" << std::endl;
  printTests(geometry, "after ");

  // check the inverse that the surface lookup uses
  const Acts::Vector3 p(30.0, 40.0, 50.0);
  const Acts::Vector3 back = geometry->transformTpcWorldToEnvelope(geometry->transformTpcEnvelopeToWorld(p));
  if ((back - p).norm() > 1e-9)
  {
    std::cerr << Name() << ": envelope <-> world round trip failed by " << (back - p).norm() << " cm" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}