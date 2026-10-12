// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SIALIGNMENTTRANSFORM_H
#define SITRACKRECO_SIALIGNMENTTRANSFORM_H

#include "SiDetectorFrame.h"

#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

// Silicon (MVTX + INTT) alignment without ACTS.
//
// Reads the same parameter file format as trackbase/AlignmentTransformation:
//     hitsetkey  alpha beta gamma  dx dy dz  [dgrx dgry dgrz]
// angles in rad, translations in mm (ACTS units).  TPC / TPOT lines are ignored.
//
// AlignmentTransformation::newMakeTransform builds, for silicon,
//   old order (default):  T = Tglob * Rglob * Tacts * Rloc * Racts
//   new order          :  T = Tglob * Rglob * Tacts * Racts * Rloc
// with Tacts*Racts the ideal sensor transform (local -> global).  Acting on an ideal global
// point p = c + Racts * p_local (c = ideal sensor centre) this is
//   old order:  p' = Rglob * (c + Rloc * (p - c)) + t
//   new order:  p' = Rglob * (c + Racts * Rloc * Racts^T * (p - c)) + t
//   Rloc = Rz(gamma) Ry(beta) Rx(alpha),  Rglob = Rz(dgrz) Ry(dgry) Rx(dgrx),  t = (dx,dy,dz)
// so the old order needs only the ideal sensor centre (exact, no ACTS surface needed); the
// new order also needs the ideal sensor axes (approximated here, see SiHitSeedReco).
//
// The transform acts in the DETECTOR frame of SiDetectorFrame (ideal sPHENIX global frame,
// before the optional detector-centre shift).  All lengths here in cm.
//
// Not supported: INTT survey-geometry mode (useInttSurveyGeometry), random perturbations.
class SiAlignmentTransform
{
 public:
  using Point = SiDetectorFrame::Point;

  struct Params
  {
    Eigen::Vector3d localAngles{0, 0, 0};     // alpha, beta, gamma [rad]
    Eigen::Vector3d translationMm{0, 0, 0};   // dx, dy, dz [mm]
    Eigen::Vector3d globalAngles{0, 0, 0};    // dgrx, dgry, dgrz [rad]
  };

  // Sensor identity with the readout-time bits (MVTX strobe, INTT crossing) dropped.
  static uint32_t sensorId(int layer, int a, int b)
  {
    return (static_cast<uint32_t>(layer) << 16) | (static_cast<uint32_t>(a & 0xff) << 8) | static_cast<uint32_t>(b & 0xff);
  }
  // MVTX: (layer, stave, chip);  INTT: (layer, ladderphi, ladderz)
  static uint32_t sensorId(TrkrDefs::hitsetkey key)
  {
    const int layer = TrkrDefs::getLayer(key);
    if (TrkrDefs::getTrkrId(key) == TrkrDefs::mvtxId)
    {
      return sensorId(layer, MvtxDefs::getStaveId(key), MvtxDefs::getChipId(key));
    }
    return sensorId(layer, InttDefs::getLadderPhiId(key), InttDefs::getLadderZId(key));
  }

  static bool isSilicon(TrkrDefs::hitsetkey key)
  {
    const auto id = TrkrDefs::getTrkrId(key);
    return id == TrkrDefs::mvtxId || id == TrkrDefs::inttId;
  }

  // ---- parsing (same rules as AlignmentTransformation::createMap)
  // Returns the number of silicon entries read, -1 if the file cannot be opened.
  int readFile(const std::string& filename, std::unordered_map<TrkrDefs::hitsetkey, Params>& out, int verbosity = 0) const
  {
    std::ifstream in(filename);
    if (!in.is_open())
    {
      return -1;
    }
    int nsi = 0;
    std::string str;
    while (std::getline(in, str))
    {
      str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char ch) { return !std::isspace(ch); }));
      if (str.empty() || str.rfind("//", 0) == 0 || str[0] == '#')
      {
        continue;
      }
      std::stringstream count_ss(str);
      std::string dummy;
      int count = 0;
      while (count_ss >> dummy)
      {
        ++count;
      }

      TrkrDefs::hitsetkey key = 0;
      Params p;
      std::stringstream ss(str);
      ss >> key >> p.localAngles(0) >> p.localAngles(1) >> p.localAngles(2) >> p.translationMm(0) >> p.translationMm(1) >> p.translationMm(2);
      if (count >= 10)
      {
        ss >> p.globalAngles(0) >> p.globalAngles(1) >> p.globalAngles(2);
      }
      if (ss.fail())
      {
        std::cout << "SiAlignmentTransform::readFile - invalid line: " << str << std::endl;
        return -1;
      }
      if (!isSilicon(key))
      {
        continue;
      }
      out[key] = p;
      ++nsi;
      if (verbosity > 1)
      {
        std::cout << "  " << key << "  " << p.localAngles.transpose() << "  " << p.translationMm.transpose()
                  << "  " << p.globalAngles.transpose() << std::endl;
      }
    }
    return nsi;
  }

  // ---- building
  void setUseNewSiliconRotationOrder(bool flag) { m_newRotationOrder = flag; }
  bool useNewSiliconRotationOrder() const { return m_newRotationOrder; }

  void clear() { m_delta.clear(); }
  size_t size() const { return m_delta.size(); }

  // centreCm: ideal sensor centre (detector frame, cm).  idealRotation: columns = ideal sensor
  // local x, y, normal in the detector frame (used only for the new rotation order).
  void addSensor(uint32_t id, const Params& p, const Eigen::Vector3d& centreCm, const Eigen::Matrix3d& idealRotation)
  {
    const Eigen::Matrix3d rloc = rotationZYX(p.localAngles);
    const Eigen::Matrix3d rglob = rotationZYX(p.globalAngles);
    const Eigen::Matrix3d r = m_newRotationOrder ? Eigen::Matrix3d(idealRotation * rloc * idealRotation.transpose()) : rloc;

    Eigen::Affine3d d = Eigen::Affine3d::Identity();
    d.pretranslate(-centreCm);
    d.prerotate(r);
    d.pretranslate(centreCm);
    d.prerotate(rglob);
    d.pretranslate(0.1 * p.translationMm);  // mm -> cm
    m_delta[id] = d;
  }

  // ---- applying
  bool has(uint32_t id) const { return m_delta.count(id) != 0; }

  // Ideal detector-frame point -> aligned detector-frame point.  Unknown sensors unchanged.
  Point apply(uint32_t id, const Point& det) const
  {
    auto it = m_delta.find(id);
    if (it == m_delta.end())
    {
      return det;
    }
    const Eigen::Vector3d q = it->second * Eigen::Vector3d(det.x, det.y, det.z);
    return {q.x(), q.y(), q.z()};
  }

  const Eigen::Affine3d* delta(uint32_t id) const
  {
    auto it = m_delta.find(id);
    return it == m_delta.end() ? nullptr : &it->second;
  }

  // Rz(a2) * Ry(a1) * Rx(a0), as gamma * beta * alpha in AlignmentTransformation
  static Eigen::Matrix3d rotationZYX(const Eigen::Vector3d& a)
  {
    const Eigen::Quaterniond q = Eigen::AngleAxisd(a(2), Eigen::Vector3d::UnitZ()) *
                                 Eigen::AngleAxisd(a(1), Eigen::Vector3d::UnitY()) *
                                 Eigen::AngleAxisd(a(0), Eigen::Vector3d::UnitX());
    return q.toRotationMatrix();
  }

 private:
  bool m_newRotationOrder = false;  // same default as AlignmentTransformation
  std::unordered_map<uint32_t, Eigen::Affine3d> m_delta;
};

#endif
