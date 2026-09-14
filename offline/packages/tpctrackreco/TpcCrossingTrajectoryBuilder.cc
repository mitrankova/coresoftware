#include "TpcCrossingTrajectoryBuilder.h"

#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"
#include "TpcCrossingTrajectory.h"
#include "TpcCrossingTrajectoryContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <set>

TpcCrossingTrajectoryBuilder::TpcCrossingTrajectoryBuilder(const std::string& name)
  : SubsysReco(name)
{
}

int TpcCrossingTrajectoryBuilder::getNodes(PHCompositeNode* topNode)
{
  m_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_trackNodeName);
  m_decisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_decisionNodeName);
  if (!m_tracks || !m_decisions)
  {
    std::cerr << Name() << "::getNodes - missing "
              << (!m_tracks ? m_trackNodeName : m_decisionNodeName) << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCrossingTrajectoryBuilder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  m_trajectories = findNode::getClass<TpcCrossingTrajectoryContainer>(topNode, m_outputNodeName);
  if (!m_trajectories)
  {
    m_trajectories = new TpcCrossingTrajectoryContainer;
    dstNode->addNode(new PHIODataNode<PHObject>(m_trajectories, m_outputNodeName, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCrossingTrajectoryBuilder::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  return createNodes(topNode);
}

int TpcCrossingTrajectoryBuilder::process_event(PHCompositeNode*)
{
  m_trajectories->Reset();
  for (unsigned int i = 0; i < m_tracks->size(); ++i)
  {
    const auto* track = m_tracks->get_track(i);
    if (!track || track->get_fit_status() == 0) continue;
    const auto* decision = m_decisions->get_decision(track->get_source_assembled_track_id());
    if (!decision) continue;

    const short referenceCrossing = decision->get_selected_crossing();
    float referenceZ0 = std::numeric_limits<float>::quiet_NaN();
    std::set<short> crossings;
    for (unsigned int j = 0; j < decision->get_number_of_candidates(); ++j)
    {
      const auto* candidate = decision->get_candidate(j);
      if (!candidate || !candidate->tpc_valid) continue;
      crossings.insert(candidate->crossing);
      if (candidate->crossing == referenceCrossing) referenceZ0 = candidate->tpc_z_at_r0;
    }
    if (crossings.empty()) crossings.insert(referenceCrossing);

    const double px = track->get_px();
    const double py = track->get_py();
    const double pz = track->get_pz();
    const double pt = std::hypot(px, py);
    const double p = std::hypot(pt, pz);
    if (!(p > 0.) || !(pt > 0.)) continue;

    for (const short crossing : crossings)
    {
      float zShift = 0.F;
      float linearChi2 = 0.F;
      for (unsigned int j = 0; j < decision->get_number_of_candidates(); ++j)
      {
        const auto* candidate = decision->get_candidate(j);
        if (!candidate || candidate->crossing != crossing) continue;
        if (std::isfinite(candidate->tpc_z_at_r0) && std::isfinite(referenceZ0))
          zShift = candidate->tpc_z_at_r0 - referenceZ0;
        linearChi2 = candidate->z_fit_chi2;
        break;
      }

      auto* trajectory = new TpcCrossingTrajectory;
      trajectory->set_parent_track_id(track->get_track_id());
      trajectory->set_crossing(crossing);
      trajectory->set_reference_crossing(referenceCrossing);
      trajectory->set_state(TpcCrossingTrajectory::X, track->get_x());
      trajectory->set_state(TpcCrossingTrajectory::Y, track->get_y());
      trajectory->set_state(TpcCrossingTrajectory::Z, track->get_z() + zShift);
      trajectory->set_state(TpcCrossingTrajectory::Phi, std::atan2(py, px));
      trajectory->set_state(TpcCrossingTrajectory::Theta, std::atan2(pt, pz));
      trajectory->set_state(TpcCrossingTrajectory::QOverP, track->get_charge() / p);
      trajectory->set_delta(TpcCrossingTrajectory::Z, zShift);
      trajectory->set_linear_chi2(linearChi2);

      const double x0 = track->get_x();
      const double y0 = track->get_y();
      const double z0 = track->get_z() + zShift;
      for (unsigned int layer = 0; layer < m_siliconRadii.size(); ++layer)
      {
        const double radius = m_siliconRadii[layer];
        const double b = x0 * px + y0 * py;
        const double c = x0 * x0 + y0 * y0 - radius * radius;
        const double disc = b * b - pt * pt * c;
        TpcCrossingTrajectory::LayerState state;
        state.layer = layer;
        if (disc >= 0.)
        {
          const double root = std::sqrt(disc);
          const double s1 = (-b + root) / (pt * pt);
          const double s2 = (-b - root) / (pt * pt);
          const double s = std::fabs(s1) < std::fabs(s2) ? s1 : s2;
          state.x = x0 + s * px;
          state.y = y0 + s * py;
          state.z = z0 + s * pz;
          state.phi = std::atan2(state.y, state.x);
          state.valid = std::isfinite(state.z);
        }
        trajectory->add_layer_state(state);
      }
      m_trajectories->add(trajectory);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
