// Copyright 2026 sPHENIX Collaboration
#ifndef TPCTRACKRECO_BEAMFRAMETRANSFORM_H
#define TPCTRACKRECO_BEAMFRAMETRANSFORM_H

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
    : m_tpcBeamLine{-0.0103, -0.0013, 0.1814, -0.0002}
    , m_mvtxBeamLine{-0.0407, -0.0015, 0.1645, -0.0001}
    // Temporary default until an independent INTT calibration is available.
    , m_inttBeamLine{-0.0407, -0.0015, 0.1645, -0.0001}
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

 private:
  BeamLine m_tpcBeamLine;
  BeamLine m_mvtxBeamLine;
  BeamLine m_inttBeamLine;
};

#endif
