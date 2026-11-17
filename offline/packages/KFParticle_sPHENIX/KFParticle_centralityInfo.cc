#include "KFParticle_centralityInfo.h"

#include <calotrigger/MinimumBiasInfo.h>
#include <centrality/CentralityInfo.h>

#include <phool/getClass.h>

#include <TTree.h>

#include <cmath>

bool KFParticle_centralityInfo::buildCentralityBranches(PHCompositeNode *topNode, TTree *m_tree)
{
  CentralityInfo *centrality = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!centrality)
  {
    return false;
  }

  m_tree->Branch("centrality_mbd_bin", &m_centrality_mbd_bin, "centrality_mbd_bin/I");
  m_tree->Branch("is_min_bias", &m_is_min_bias, "is_min_bias/O");

  return true;
}

void KFParticle_centralityInfo::fillCentralityBranches(PHCompositeNode *topNode, bool crossingIsZero)
{
  MinimumBiasInfo *mbinfo = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
  if (mbinfo)
  {
    m_is_min_bias = mbinfo->isAuAuMinimumBias();
  }

  if (!crossingIsZero)
  {
    return;
  }

  CentralityInfo *centrality = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!centrality)
  {
    return;
  }

  //! CentralityReco only sets this for events passing the minimum bias
  //! selection, so has_centile is false for the rest and the defaults stand.
  //! That is the normal case for a non minimum bias event, not an error.
  if (!centrality->has_centile(CentralityInfo::PROP::mbd_NS))
  {
    return;
  }

  const float centile = centrality->get_centile(CentralityInfo::PROP::mbd_NS);

  //! CentralityReco stores the centile as (bin + 1) / NDIVS and the stored float
  //! sits just below the true value, so round rather than truncate. The -1
  //! converts the 1 indexed CentralityReco bin to the 0 to 99 convention.
  m_centrality_mbd_bin = (int) std::lround(centile * 100.) - 1;
}

void KFParticle_centralityInfo::resetCentralityBranches()
{
  m_centrality_mbd_bin = -999;
  m_is_min_bias = false;
}
