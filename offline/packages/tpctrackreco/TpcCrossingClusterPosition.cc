#include "TpcCrossingClusterPosition.h"

#include "TpcDriftPolylineLookup.h"
#include "Tpc_PolyCluster.h"

#include <cmath>

bool TpcCrossingClusterPosition::get(
    const Tpc_PolyCluster& cluster,
    const TpcDriftPolylineLookup& lookup,
    const short crossing,
    const short referenceCrossing,
    std::array<double, 3>& position,
    std::vector<HitPosition>* hits)
{
  const auto& indices = cluster.get_hit_indices();
  if (indices.empty()) return false;

  if (hits)
  {
    hits->clear();
    hits->reserve(indices.size());
  }

  std::array<double, 3> deltaSum{};
  for (const auto& index : indices)
  {
    TpcDriftPolylineLookup::Point delta;
    if (!lookup.getDelta(index.first, index.second, crossing, referenceCrossing, delta)) return false;
    deltaSum[0] += delta.x;
    deltaSum[1] += delta.y;
    deltaSum[2] += delta.z;

    if (hits)
    {
      TpcDriftPolylineLookup::Point reference;
      TpcDriftPolylineLookup::Point candidate;
      if (!lookup.getPosition(index.first, index.second, referenceCrossing, reference) ||
          !lookup.getPosition(index.first, index.second, crossing, candidate)) return false;
      hits->push_back({index.first, index.second,
                       {reference.x, reference.y, reference.z},
                       {candidate.x, candidate.y, candidate.z},
                       {delta.x, delta.y, delta.z}});
    }
  }

  const double inverseCount = 1.0 / static_cast<double>(indices.size());
  position = {cluster.get_centroid_x() + deltaSum[0] * inverseCount,
              cluster.get_centroid_y() + deltaSum[1] * inverseCount,
              cluster.get_centroid_z() + deltaSum[2] * inverseCount};
  return std::isfinite(position[0]) && std::isfinite(position[1]) && std::isfinite(position[2]);
}
