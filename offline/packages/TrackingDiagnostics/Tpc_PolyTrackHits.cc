#include "Tpc_PolyTrackHits.h"

#include "tpctrackreco/IdealPadMap.h"
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
#include <vector>

Tpc_PolyTrackHits::Tpc_PolyTrackHits(const std::string& name,
                                     const std::string& outfilename)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_clusterNodeName("TPC_POLYCLUSTERS")
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

  m_tree = new TTree("hits", "TPC hits on selected polytracks");
  m_tree->Branch("event", &m_event, "event/i");
  m_tree->Branch("poly_track_id", &m_polyTrackId, "poly_track_id/i");
  m_tree->Branch("source_assembled_track_id", &m_sourceAssembledTrackId, "source_assembled_track_id/i");
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
  m_event = m_evt;
  m_polyTrackId = 0;
  m_sourceAssembledTrackId = 0;
  m_polyClusterId = 0;
  m_ntpcClusters = 0;
  m_pt = std::numeric_limits<double>::quiet_NaN();
  m_pad = 0;
  m_timebin = 0;
  m_layer = 0;
  m_sector = 0;
  m_side = 0;
  m_phi = std::numeric_limits<double>::quiet_NaN();
  m_radius = std::numeric_limits<double>::quiet_NaN();
  m_x = std::numeric_limits<double>::quiet_NaN();
  m_y = std::numeric_limits<double>::quiet_NaN();
  m_z = std::numeric_limits<double>::quiet_NaN();
  m_adc = std::numeric_limits<double>::quiet_NaN();
}

int Tpc_PolyTrackHits::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  if (!get_nodes(topNode) || !m_tree || !m_idealPadMap)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  std::map<unsigned int, std::vector<const Tpc_PolyCluster*> > clusters_by_source_assembled_track_id;
  for (unsigned int icluster = 0; icluster < m_clusters->size(); ++icluster)
  {
    const Tpc_PolyCluster* cluster = m_clusters->get_cluster(icluster);
    if (!cluster || !cluster->isValid())
    {
      continue;
    }
    clusters_by_source_assembled_track_id[cluster->get_source_assembled_track_id()].push_back(cluster);
  }

  unsigned int nfilled = 0;
  const unsigned int npoly_tracks = m_finalTracks->size();
  for (unsigned int itrack = 0; itrack < npoly_tracks; ++itrack)
  {
    const Tpc_PolyTrack* poly_track = m_finalTracks->get_track(itrack);
    if (!poly_track || !poly_track->isValid())
    {
      continue;
    }

    const double px = poly_track->get_px();
    const double py = poly_track->get_py();
    const double pt = std::hypot(px, py);
    if (!std::isfinite(pt) || pt <= m_minPt)
    {
      continue;
    }

    const auto cluster_iter = clusters_by_source_assembled_track_id.find(poly_track->get_source_assembled_track_id());
    if (cluster_iter == clusters_by_source_assembled_track_id.end())
    {
      continue;
    }

    const std::vector<const Tpc_PolyCluster*>& track_clusters = cluster_iter->second;
    const unsigned int ntpc_clusters = track_clusters.size();
    if (ntpc_clusters <= m_minTpcClusters)
    {
      continue;
    }

    for (const Tpc_PolyCluster* cluster : track_clusters)
    {
      if (!cluster || !cluster->isValid())
      {
        continue;
      }

      for (unsigned int ihit = 0; ihit < cluster->size_hits(); ++ihit)
      {
        const Tpc_PolyCluster::HitIndex hit_index = cluster->get_hit_index(ihit);
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

        const unsigned int layer = TrkrDefs::getLayer(hit_index.first);
        const unsigned int side = TpcDefs::getSide(hit_index.first);
        const unsigned int sector = TpcDefs::getSectorId(hit_index.first);
        const unsigned int pad = TpcDefs::getPad(hit_index.second);
        const unsigned int timebin = TpcDefs::getTBin(hit_index.second);
        const double phi = m_idealPadMap->get_phi(side, layer, pad);
        const double radius = m_idealPadMap->get_radius(layer);
        if (!std::isfinite(phi) || !std::isfinite(radius))
        {
          continue;
        }

        reset_tree_values();
        m_event = m_evt;
        m_polyTrackId = poly_track->get_track_id();
        m_sourceAssembledTrackId = poly_track->get_source_assembled_track_id();
        m_polyClusterId = cluster->get_cluster_id();
        m_ntpcClusters = ntpc_clusters;
        m_pt = pt;
        m_pad = pad;
        m_timebin = timebin;
        m_layer = layer;
        m_sector = sector;
        m_side = side;
        m_phi = phi;
        m_radius = radius;
        m_x = cluster->get_hit_x(ihit);
        m_y = cluster->get_hit_y(ihit);
        m_z = cluster->get_hit_z(ihit);
        m_adc = hit->getAdc();
        m_tree->Fill();
        ++nfilled;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_evt
              << " poly_tracks=" << npoly_tracks
              << " hits=" << nfilled << std::endl;
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
