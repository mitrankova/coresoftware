// Copyright 2026 sPHENIX Collaboration
#ifndef TPCTRACKRECO_BEAMFRAMETRANSFORM_H
#define TPCTRACKRECO_BEAMFRAMETRANSFORM_H

#include <initializer_list>

/** Convert detector-specific global coordinates to one beam-centered frame. */
class BeamFrameTransform
{
 public:
  struct BeamLine
  {
    double x0{0.0};
    double dxdz{0.0};
    double y0{0.0};
    double dydz{0.0};
  };

  struct Position
  {
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };

  BeamFrameTransform()
    //: m_tpcBeamLine{-0.0103, -0.0013, 0.1814, -0.0002}
    //, m_mvtxBeamLine{-0.0407, -0.0015, 0.1645, -0.0001}
    //, m_inttBeamLine{-0.0407, -0.0015, 0.1645, -0.0001}
    : m_tpcBeamLine{0,0, 0,0}
    , m_mvtxBeamLine{0, 0,0, 0}
    , m_inttBeamLine{0, 0, 0, 0}
  {
  }

  void setTpcBeamLine(const BeamLine& line) { m_tpcBeamLine = line; }
  void setMvtxBeamLine(const BeamLine& line) { m_mvtxBeamLine = line; }
  void setInttBeamLine(const BeamLine& line) { m_inttBeamLine = line; }

  const BeamLine& tpcBeamLine() const { return m_tpcBeamLine; }
  const BeamLine& mvtxBeamLine() const { return m_mvtxBeamLine; }
  const BeamLine& inttBeamLine() const { return m_inttBeamLine; }

  Position toBeamFrame(const BeamLine& line, double x, double y, double z) const
  {
    return {x - (line.x0 + line.dxdz * z),
            y - (line.y0 + line.dydz * z), z};
  }

  bool validate(double tolerance = 1.0e-12) const
  {
    const auto valid = [this, tolerance](const BeamLine& line)
    {
      for (const double z : {-100.0, 0.0, 100.0})
      {
        const auto centered = toBeamFrame(
            line, line.x0 + line.dxdz * z, line.y0 + line.dydz * z, z);
        if (centered.x < -tolerance || centered.x > tolerance ||
            centered.y < -tolerance || centered.y > tolerance)
        {
          return false;
        }
      }
      return true;
    };
    return valid(m_tpcBeamLine) && valid(m_mvtxBeamLine) && valid(m_inttBeamLine);
  }

 private:
  BeamLine m_tpcBeamLine;
  BeamLine m_mvtxBeamLine;
  BeamLine m_inttBeamLine;
};

#endif  // TPCTRACKRECO_BEAMFRAMETRANSFORM_H
