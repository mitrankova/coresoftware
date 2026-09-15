#ifndef TPCTRACKRECO_TPCCROSSINGCLUSTERPOSITION_H
#define TPCTRACKRECO_TPCCROSSINGCLUSTERPOSITION_H

#include <trackbase/TrkrDefs.h>

#include <array>
#include <vector>

class TpcDriftPolylineLookup;
class Tpc_PolyCluster;

class TpcCrossingClusterPosition
{
 public:
  struct HitPosition
  {
    TrkrDefs::hitsetkey hitsetkey{};
    TrkrDefs::hitkey hitkey{};
    std::array<double, 3> reference{};
    std::array<double, 3> candidate{};
    std::array<double, 3> delta{};
  };

  static bool get(const Tpc_PolyCluster& cluster,
                  const TpcDriftPolylineLookup& lookup,
                  short crossing,
                  short referenceCrossing,
                  std::array<double, 3>& position,
                  std::vector<HitPosition>* hits = nullptr);
};

#endif
