// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SITRACKRECO_SIHITINDEX_H
#define SITRACKRECO_SIHITINDEX_H

#include "SiHitGeometryTable.h"
#include "SiHitSeedData.h"

#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrDefs.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

// One silicon hit in the per-event index: physics-space position (SiHitGeometryTable) plus the
// readout address needed for the geometric neighbour test.
struct SiIndexedHit
{
  float u = 0;  // intrinsic azimuth [0, 2pi), detector frame
  float z = 0;  // [cm]
  float r = 0;  // [cm]
  int hit = -1;  // index in SiHitSeedEvent::hits
  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::HITSETKEYMAX;
  TrkrDefs::hitkey hitkey = TrkrDefs::HITKEYMAX;
  uint16_t row = 0;
  uint16_t col = 0;
  int16_t phiId = -1;  // MVTX stave / INTT ladder phi
  int16_t zId = -1;    // MVTX chip / INTT ladder z
  int time = 0;        // MVTX strobe / INTT time bucket
};

// SiHitIndex
//
// Per-event, per-layer index of the silicon hits:
//   - azimuth buckets of fixed width in r*phi (setBucketRphiCm), the same physical size on
//     every layer;
//   - inside a bucket the hits are sorted by z, so a z range is two binary searches.
// A window query visits only the buckets overlapping [u - du, u + du] and only the hits with
// z in [zlo, zhi].  No trig anywhere: u, z, r come from the geometry table.
class SiHitIndex
{
 public:
  static constexpr int kNLayers = SiHitGeometryTable::kNLayers;

  void setBucketRphiCm(double w) { m_bucketRphi = w > 0 ? w : 0.1; }

  void build(const std::vector<SiHitPoint>& hits, const SiHitGeometryTable& geo)
  {
    std::array<std::vector<std::pair<int, SiIndexedHit>>, kNLayers> tmp;
    for (int layer = 0; layer < kNLayers; ++layer)
    {
      const double R = SiHitGeometryTable::radius(layer);
      m_nb[layer] = std::max(1, static_cast<int>(std::floor(SiHitGeometryTable::kTwoPi * R / m_bucketRphi)));
      m_bw[layer] = SiHitGeometryTable::kTwoPi / m_nb[layer];
      tmp[layer].clear();
    }
    for (std::size_t i = 0; i < hits.size(); ++i)
    {
      const auto& h = hits[i];
      if (h.layer < 0 || h.layer >= kNLayers)
      {
        continue;
      }
      const auto p = geo.position(h.layer, h.row, h.phi, h.z);
      SiIndexedHit e;
      e.u = p.u;
      e.z = p.z;
      e.r = p.r;
      e.hit = static_cast<int>(i);
      e.hitsetkey = h.hitsetkey;
      e.hitkey = h.hitkey;
      e.row = static_cast<uint16_t>(h.row);
      e.col = static_cast<uint16_t>(h.col);
      if (SiHitGeometryTable::isMvtx(h.layer))
      {
        e.phiId = static_cast<int16_t>(h.stave);
        e.zId = static_cast<int16_t>(h.chip);
        e.time = MvtxDefs::getStrobeId(h.hitsetkey);
      }
      else
      {
        e.phiId = static_cast<int16_t>(h.ladderphi);
        e.zId = static_cast<int16_t>(h.ladderz);
        e.time = InttDefs::getTimeBucketId(h.hitsetkey);
      }
      int b = static_cast<int>(e.u / m_bw[h.layer]);
      b = std::min(std::max(b, 0), m_nb[h.layer] - 1);
      tmp[h.layer].emplace_back(b, e);
    }
    for (int layer = 0; layer < kNLayers; ++layer)
    {
      auto& t = tmp[layer];
      std::sort(t.begin(), t.end(), [](const auto& a, const auto& b)
                { return a.first != b.first ? a.first < b.first : a.second.z < b.second.z; });
      auto& out = m_hits[layer];
      auto& start = m_start[layer];
      out.clear();
      out.reserve(t.size());
      start.assign(m_nb[layer] + 1, 0);
      for (const auto& [b, e] : t)
      {
        ++start[b + 1];
        out.push_back(e);
      }
      for (int b = 0; b < m_nb[layer]; ++b)
      {
        start[b + 1] += start[b];
      }
    }
  }

  const std::vector<SiIndexedHit>& hits(int layer) const { return m_hits[layer]; }
  std::size_t size(int layer) const { return m_hits[layer].size(); }
  std::size_t size() const
  {
    std::size_t n = 0;
    for (const auto& v : m_hits)
    {
      n += v.size();
    }
    return n;
  }

  // f(k) for every hit k of `layer` (index into hits(layer)) in the buckets covering
  // [uc - du, uc + du] (periodic) with z in [zlo, zhi].  The caller does the exact u test.
  template <class F>
  void forEachInWindow(int layer, double uc, double du, double zlo, double zhi, F&& f) const
  {
    const int nb = m_nb[layer];
    const double bw = m_bw[layer];
    const auto& v = m_hits[layer];
    const auto& start = m_start[layer];
    if (v.empty() || zhi < zlo)
    {
      return;
    }
    int b0 = 0, b1 = nb - 1;
    if (du < M_PI)
    {
      b0 = static_cast<int>(std::floor((uc - du) / bw));
      b1 = static_cast<int>(std::floor((uc + du) / bw));
      if (b1 - b0 + 1 >= nb)
      {
        b0 = 0;
        b1 = nb - 1;
      }
    }
    const auto byZ = [](const SiIndexedHit& h, double z) { return h.z < z; };
    for (int b = b0; b <= b1; ++b)
    {
      const int bb = ((b % nb) + nb) % nb;
      const auto first = v.begin() + start[bb];
      const auto last = v.begin() + start[bb + 1];
      for (auto it = std::lower_bound(first, last, zlo, byZ); it != last && it->z <= zhi; ++it)
      {
        f(static_cast<int>(it - v.begin()));
      }
    }
  }

 private:
  double m_bucketRphi = 0.1;  // [cm]
  std::array<int, kNLayers> m_nb{};
  std::array<double, kNLayers> m_bw{};
  std::array<std::vector<SiIndexedHit>, kNLayers> m_hits;
  std::array<std::vector<int>, kNLayers> m_start;
};

#endif
