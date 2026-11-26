#include "Tpc_PolyTrackHits.h"

#include "tpctrackreco/IdealPadMap.h"
#include "tpctrackreco/TpcCrossingDecision.h"
#include "tpctrackreco/TpcCrossingDecisionContainer.h"
#include "tpctrackreco/Tpc_PolyCluster.h"
#include "tpctrackreco/Tpc_PolyClusterContainer.h"
#include "tpctrackreco/Tpc_PolyTrack.h"
#include "tpctrackreco/Tpc_PolyTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

Tpc_PolyTrackHits::Tpc_PolyTrackHits(const std::string& name,
                                     const std::string& outfilename)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_clusterNodeName("TPC_POLYCLUSTERS")
  , m_crossingDecisionNodeName("TPC_CROSSING_DECISIONS")
  , m_finalTrackNodeName("TPC_POLYTRACKS")
  , m_hitSetNodeName("TRKR_HITSET")
{
}

Tpc_PolyTrackHits::~Tpc_PolyTrackHits()
{
  if (m_outfile)
  {
    delete m_outfile;
    m_outfile = nullptr;
  }
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
}

int Tpc_PolyTrackHits::Init(PHCompositeNode* /*unused*/)
{
  // cppcheck-suppress publicAllocationError
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open output file " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("hits", "TPC hits: on polytracks and unassigned (track fields = 32767)");
  m_tree->Branch("event", &m_event, "event/i");
  m_tree->Branch("poly_track_id", &m_polyTrackId, "poly_track_id/i");
  m_tree->Branch("source_assembled_track_id", &m_sourceAssembledTrackId, "source_assembled_track_id/i");
  m_tree->Branch("crossing", &m_crossing, "crossing/S");
  m_tree->Branch("poly_cluster_id", &m_polyClusterId, "poly_cluster_id/i");
  m_tree->Branch("ntpc_clusters", &m_ntpcClusters, "ntpc_clusters/i");
  m_tree->Branch("pt", &m_pt, "pt/D");
  m_tree->Branch("pad", &m_pad, "pad/i");
  m_tree->Branch("timebin", &m_timebin, "timebin/i");
  m_tree->Branch("layer", &m_layer, "layer/i");
  m_tree->Branch("sector", &m_sector, "sector/i");
  m_tree->Branch("side", &m_side, "side/i");
  m_tree->Branch("phi", &m_phi, "phi/D");
  m_tree->Branch("radius", &m_radius, "radius/D");
  m_tree->Branch("poly_cluster_x", &m_polyClusterX, "poly_cluster_x/D");
  m_tree->Branch("poly_cluster_y", &m_polyClusterY, "poly_cluster_y/D");
  m_tree->Branch("poly_cluster_z", &m_polyClusterZ, "poly_cluster_z/D");
  m_tree->Branch("poly_cluster_rms_x", &m_polyClusterRmsX, "poly_cluster_rms_x/D");
  m_tree->Branch("poly_cluster_rms_y", &m_polyClusterRmsY, "poly_cluster_rms_y/D");
  m_tree->Branch("poly_cluster_rms_z", &m_polyClusterRmsZ, "poly_cluster_rms_z/D");
  m_tree->Branch("poly_cluster_rphi_error", &m_polyClusterRPhiError, "poly_cluster_rphi_error/D");
  m_tree->Branch("poly_cluster_z_error", &m_polyClusterZError, "poly_cluster_z_error/D");
  m_tree->Branch("x", &m_x, "x/D");
  m_tree->Branch("y", &m_y, "y/D");
  m_tree->Branch("z", &m_z, "z/D");
  m_tree->Branch("adc", &m_adc, "adc/D");

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyTrackHits::InitRun(PHCompositeNode* /*topNode*/)
{
  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - failed to load IdealPadMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool Tpc_PolyTrackHits::get_nodes(PHCompositeNode* topNode)
{
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_clusterNodeName);
  if (!m_clusters)
  {
    std::cerr << Name() << " - missing " << m_clusterNodeName << std::endl;
    return false;
  }

  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  if (!m_crossingDecisions)
  {
    std::cerr << Name() << " - missing " << m_crossingDecisionNodeName << std::endl;
    return false;
  }

  m_finalTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_finalTrackNodeName);
  if (!m_finalTracks)
  {
    std::cerr << Name() << " - missing " << m_finalTrackNodeName << std::endl;
    return false;
  }

  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, m_hitSetNodeName);
  if (!m_hits)
  {
    std::cerr << Name() << " - missing " << m_hitSetNodeName << std::endl;
    return false;
  }

  return true;
}

void Tpc_PolyTrackHits::reset_tree_values()
{
  // Everything that may be undefined for a hit is set to short max (32767)
  m_event = m_evt;
  m_polyTrackId = kUndefUInt;
  m_sourceAssembledTrackId = kUndefUInt;
  m_crossing = kUndefShort;
  m_polyClusterId = kUndefUInt;
  m_ntpcClusters = kUndefUInt;
  m_pt = kUndefDouble;
  m_pad = kUndefUInt;
  m_timebin = kUndefUInt;
  m_layer = kUndefUInt;
  m_sector = kUndefUInt;
  m_side = kUndefUInt;
  m_phi = kUndefDouble;
  m_radius = kUndefDouble;
  m_polyClusterX = kUndefDouble;
  m_polyClusterY = kUndefDouble;
  m_polyClusterZ = kUndefDouble;
  m_polyClusterRmsX = kUndefDouble;
  m_polyClusterRmsY = kUndefDouble;
  m_polyClusterRmsZ = kUndefDouble;
  m_polyClusterRPhiError = kUndefDouble;
  m_polyClusterZError = kUndefDouble;
  m_x = kUndefDouble;
  m_y = kUndefDouble;
  m_z = kUndefDouble;
  m_adc = kUndefDouble;
}

bool Tpc_PolyTrackHits::fill_hit_fields(uint32_t hitsetkey, uint32_t hitkey, const TrkrHit* hit)
{
  const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  const unsigned int side = TpcDefs::getSide(hitsetkey);
  const unsigned int sector = TpcDefs::getSectorId(hitsetkey);
  const unsigned int pad = TpcDefs::getPad(hitkey);
  const unsigned int timebin = TpcDefs::getTBin(hitkey);
  const double phi = m_idealPadMap->get_phi(side, layer, pad);
  const double radius = m_idealPadMap->get_radius(layer);
  if (!std::isfinite(phi) || !std::isfinite(radius))
  {
    return false;
  }

  m_pad = pad;
  m_timebin = timebin;
  m_layer = layer;
  m_sector = sector;
  m_side = side;
  m_phi = phi;
  m_radius = radius;
  m_adc = hit->getAdc();
  return true;
}

void Tpc_PolyTrackHits::fill_cluster_fields(const Tpc_PolyCluster* cluster, unsigned int ihit)
{
  m_polyClusterId = cluster->get_cluster_id();
  m_polyClusterX = cluster->get_centroid_x();
  m_polyClusterY = cluster->get_centroid_y();
  m_polyClusterZ = cluster->get_centroid_z();
  m_polyClusterRmsX = cluster->get_rms_x();
  m_polyClusterRmsY = cluster->get_rms_y();
  m_polyClusterRmsZ = cluster->get_rms_z();
  m_polyClusterRPhiError = std::hypot(cluster->get_rms_x(), cluster->get_rms_y());
  m_polyClusterZError = std::fabs(cluster->get_rms_z());
  m_x = cluster->get_hit_x(ihit);
  m_y = cluster->get_hit_y(ihit);
  m_z = cluster->get_hit_z(ihit);
}

int Tpc_PolyTrackHits::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  if (!get_nodes(topNode) || !m_tree || !m_idealPadMap)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // clusters grouped by assembled track, and hit -> (first) cluster lookup for all valid clusters
  std::map<unsigned int, std::vector<const Tpc_PolyCluster*> > clusters_by_source_assembled_track_id;
  std::unordered_map<uint64_t, std::pair<const Tpc_PolyCluster*, unsigned int> > cluster_by_hit;
  for (unsigned int icluster = 0; icluster < m_clusters->size(); ++icluster)
  {
    const Tpc_PolyCluster* cluster = m_clusters->get_cluster(icluster);
    if (!cluster || !cluster->isValid())
    {
      continue;
    }
    clusters_by_source_assembled_track_id[cluster->get_source_assembled_track_id()].push_back(cluster);

    for (unsigned int ihit = 0; ihit < cluster->size_hits(); ++ihit)
    {
      const Tpc_PolyCluster::HitIndex hit_index = cluster->get_hit_index(ihit);
      cluster_by_hit.emplace(make_hit_id(hit_index.first, hit_index.second), std::make_pair(cluster, ihit));
    }
  }

  std::map<unsigned int, const TpcCrossingDecision*> crossing_by_source_assembled_track_id;
  for (unsigned int idecision = 0; idecision < m_crossingDecisions->size(); ++idecision)
  {
    const TpcCrossingDecision* decision = m_crossingDecisions->get_decision(idecision);
    if (!decision || !decision->isValid())
    {
      continue;
    }
    crossing_by_source_assembled_track_id[decision->get_assembled_track_id()] = decision;
  }

  // hits that belong to a found track (written on-track, or skipped by an optional track cut)
  std::unordered_set<uint64_t> track_hits;

  unsigned int nfilled_track = 0;
  unsigned int nfilled_free = 0;
  unsigned int ntracks_written = 0;
  const unsigned int npoly_tracks = m_finalTracks->size();

  // ---------------- pass 1: all found tracks ----------------
  for (unsigned int itrack = 0; itrack < npoly_tracks; ++itrack)
  {
    const Tpc_PolyTrack* poly_track = m_finalTracks->get_track(itrack);
    if (!poly_track || !poly_track->isValid())
    {
      continue;
    }

    const unsigned int source_assembled_track_id = poly_track->get_source_assembled_track_id();
    const auto cluster_iter = clusters_by_source_assembled_track_id.find(source_assembled_track_id);
    if (cluster_iter == clusters_by_source_assembled_track_id.end())
    {
      continue;  // no clusters -> nothing to write
    }
    const std::vector<const Tpc_PolyCluster*>& track_clusters = cluster_iter->second;
    const unsigned int ntpc_clusters = track_clusters.size();

    const double pt = std::hypot(poly_track->get_px(), poly_track->get_py());

    const auto crossing_iter = crossing_by_source_assembled_track_id.find(source_assembled_track_id);
    const short crossing = (crossing_iter != crossing_by_source_assembled_track_id.end())
                               ? crossing_iter->second->get_selected_crossing()
                               : kUndefShort;

    // optional cuts (disabled by default). Rejected tracks still claim their hits.
    const bool pass_pt = (m_minPt < 0) || (std::isfinite(pt) && pt >= m_minPt);
    const bool pass_cuts = pass_pt && ntpc_clusters >= m_minTpcClusters;

    bool wrote_any = false;
    for (const Tpc_PolyCluster* cluster : track_clusters)
    {
      for (unsigned int ihit = 0; ihit < cluster->size_hits(); ++ihit)
      {
        const Tpc_PolyCluster::HitIndex hit_index = cluster->get_hit_index(ihit);
        track_hits.insert(make_hit_id(hit_index.first, hit_index.second));
        if (!pass_cuts)
        {
          continue;
        }

        TrkrHitSet* hitset = m_hits->findHitSet(hit_index.first);
        if (!hitset)
        {
          continue;
        }
        TrkrHit* hit = hitset->getHit(hit_index.second);
        if (!hit)
        {
          continue;
        }

        reset_tree_values();
        if (!fill_hit_fields(hit_index.first, hit_index.second, hit))
        {
          continue;
        }
        m_polyTrackId = poly_track->get_track_id();
        m_sourceAssembledTrackId = source_assembled_track_id;
        m_crossing = crossing;
        m_ntpcClusters = ntpc_clusters;
        m_pt = std::isfinite(pt) ? pt : kUndefDouble;
        fill_cluster_fields(cluster, ihit);
        m_tree->Fill();
        ++nfilled_track;
        wrote_any = true;
      }
    }
    if (wrote_any)
    {
      ++ntracks_written;
    }
  }

  // ---------------- pass 2: TPC hits not on any track ----------------
  if (m_saveUnassignedHits)
  {
    const TrkrHitSetContainer::ConstRange hitset_range = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
    for (auto hitset_iter = hitset_range.first; hitset_iter != hitset_range.second; ++hitset_iter)
    {
      const TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
      TrkrHitSet* hitset = hitset_iter->second;
      if (!hitset)
      {
        continue;
      }

      const TrkrHitSet::ConstRange hit_range = hitset->getHits();
      for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
      {
        const TrkrDefs::hitkey hitkey = hit_iter->first;
        const TrkrHit* hit = hit_iter->second;
        if (!hit)
        {
          continue;
        }

        const uint64_t hit_id = make_hit_id(hitsetkey, hitkey);
        if (track_hits.count(hit_id))
        {
          continue;
        }

        reset_tree_values();  // track fields stay at short max
        if (!fill_hit_fields(hitsetkey, hitkey, hit))
        {
          continue;
        }

        // hit may still sit in a cluster that is not attached to any track
        const auto cl_iter = cluster_by_hit.find(hit_id);
        if (cl_iter != cluster_by_hit.end())
        {
          fill_cluster_fields(cl_iter->second.first, cl_iter->second.second);
        }

        m_tree->Fill();
        ++nfilled_free;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_evt
              << " poly_tracks=" << npoly_tracks
              << " tracks_written=" << ntracks_written
              << " track_hits=" << nfilled_track
              << " unassigned_hits=" << nfilled_free << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyTrackHits::End(PHCompositeNode* /*unused*/)
{
  if (m_outfile)
  {
    m_outfile->cd();
    if (m_tree)
    {
      m_tree->Write();
    }
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
    m_tree = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}