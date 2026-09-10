// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCCROSSINGCLUSTERCORRECTOR_H
#define TPCTRACKRECO_TPCCROSSINGCLUSTERCORRECTOR_H

#include "Tpc_PolyClusterizer.h"

#include <string>

/** Rebuild TPC poly-clusters from raw pad/time-bin hits at the silicon-refined crossing.
 *
 * This reuses Tpc_PolyClusterizer's IdealPadMap and PHGarfield reverse-drift
 * implementation. The original TPC_POLYCLUSTERS node is left untouched.
 */
class TpcCrossingClusterCorrector : public Tpc_PolyClusterizer
{
 public:
  explicit TpcCrossingClusterCorrector(const std::string& name = "TpcCrossingClusterCorrector")
    : Tpc_PolyClusterizer(name)
  {
    setCrossingDecisionNodeName("TPC_SILICON_CROSSING_DECISIONS");
    setOutputNodeName("TPC_POLYCLUSTERS_CROSSING_CORRECTED");
  }

  ~TpcCrossingClusterCorrector() override = default;
};

#endif
