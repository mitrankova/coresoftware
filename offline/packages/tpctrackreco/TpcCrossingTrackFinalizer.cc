#include "TpcCrossingTrackFinalizer.h"
#include "Full_PolyTrackContainerv1.h"
#include "Full_PolyTrackv1.h"
#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
#include "TpcSiliconMatchCandidate.h"
#include "TpcSiliconMatchCandidateContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <iostream>

TpcCrossingTrackFinalizer::TpcCrossingTrackFinalizer(const std::string& name) : SubsysReco(name) {}
int TpcCrossingTrackFinalizer::getNodes(PHCompositeNode* topNode)
{
  m_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_trackNodeName);
  m_trajectories = findNode::getClass<TpcCrossingTrajectoryContainer>(topNode, m_trajectoryNodeName);
  m_candidates = findNode::getClass<TpcSiliconMatchCandidateContainer>(topNode, m_candidateNodeName);
  if (!m_tracks || !m_trajectories || !m_candidates)
  {
    std::cerr << Name() << "::getNodes - missing full-track input node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrackFinalizer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dst = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst) { dst = new PHCompositeNode("DST"); topNode->addNode(dst); }
  m_output = findNode::getClass<Full_PolyTrackContainer>(topNode, m_outputNodeName);
  if (!m_output)
  {
    m_output = new Full_PolyTrackContainerv1;
    dst->addNode(new PHIODataNode<PHObject>(m_output, m_outputNodeName, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int TpcCrossingTrackFinalizer::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  m_event = 0;
  return createNodes(topNode);
}
int TpcCrossingTrackFinalizer::process_event(PHCompositeNode*)
{
  ++m_event;
  m_output->Reset();
  unsigned int outputId = 0;
  for (unsigned int i = 0; i < m_candidates->size(); ++i)
  {
    const auto* candidate = m_candidates->get(i);
    if (!candidate || !candidate->get_selected()) continue;
    const Tpc_PolyTrack* parent = nullptr;
    for (unsigned int j = 0; j < m_tracks->size(); ++j)
    {
      const auto* value = m_tracks->get_track(j);
      if (value && value->get_track_id() == candidate->get_parent_track_id()) { parent = value; break; }
    }
    if (!parent) continue;
    const TpcCrossingTrajectory* trajectory = nullptr;
    for (unsigned int j = 0; j < m_trajectories->size(); ++j)
    {
      const auto* value = m_trajectories->get(j);
      if (value && value->get_parent_track_id() == candidate->get_parent_track_id() && value->get_crossing() == candidate->get_crossing()) { trajectory = value; break; }
    }
    if (!trajectory) continue;

    auto* full = new Full_PolyTrackv1;
    full->set_event(m_event);
    full->set_track_id(outputId++);
    full->set_parent_tpc_track_id(parent->get_track_id());
    full->set_source_assembled_track_id(parent->get_source_assembled_track_id());
    full->set_crossing(candidate->get_crossing());
    full->set_status(1);
    full->set_n_mvtx(candidate->get_n_mvtx());
    full->set_n_intt(candidate->get_n_intt());
    full->set_score(candidate->get_score());
    full->set_max_abs_dz(candidate->get_max_abs_dz());
    full->set_max_abs_ddphi(candidate->get_max_abs_ddphi());
    full->set_fit_status(parent->get_fit_status());
    full->set_chi2(parent->get_chi2());
    full->set_ndf(parent->get_ndf());
    full->set_x(trajectory->get_state(TpcCrossingTrajectory::X));
    full->set_y(trajectory->get_state(TpcCrossingTrajectory::Y));
    full->set_z(trajectory->get_state(TpcCrossingTrajectory::Z));
    full->set_px(parent->get_px()); full->set_py(parent->get_py()); full->set_pz(parent->get_pz()); full->set_charge(parent->get_charge());
    for (unsigned int row = 0; row < 6; ++row) for (unsigned int col = 0; col < 6; ++col) full->set_cov(row, col, parent->get_cov(row, col));
    for (const auto key : parent->get_cluster_keys()) full->add_tpc_cluster_key(key);
    for (const auto key : candidate->get_silicon_cluster_keys()) full->add_silicon_cluster_key(key);
    m_output->add_track(full);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
