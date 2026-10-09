#include "SiTrajectoryFitter.h"

#include "SiHitSeedData.h"
#include "SiTpcHelixFit.h"
#include "Si_Trajectoryv1.h"
#include "Si_TrajectoryContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>

SiTrajectoryFitter::SiTrajectoryFitter(const std::string& name)
  : SubsysReco(name)
{
}

int SiTrajectoryFitter::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    const auto& c = m_frame.detectorCenterCm();
    std::cout << Name() << ": Uphi rotation " << m_frame.uphiRotation() << " rad"
              << ", detector centre (" << 10 * c[0] << ", " << 10 * c[1] << ", " << 10 * c[2] << ") mm"
              << ", beam (" << m_beamX << ", " << m_beamY << ") cm"
              << ", Bz " << m_bz << " T"
              << ", use INTT " << m_useIntt << ", min points " << m_minPoints
              << ", straight line if radial lever arm < " << m_minCurvatureLeverArm << " cm"
              << ", beam constraint " << (m_beamConstraint ? "on" : "off")
              << ", z weights MVTX/INTT " << m_zWeightMvtx << "/" << m_zWeightIntt << std::endl;
    if (m_alignment.enabled())
    {
      const char* names[] = {"TPC", "Si half A", "Si half B"};
      for (int part = SiTpcBeamAlignment::SiHalfA; part <= SiTpcBeamAlignment::SiHalfB; ++part)
      {
        const auto& l = m_alignment.beamLine(part);
        std::cout << Name() << ":   beam-axis alignment " << names[part] << ": x0 " << 10 * l.x0 << " mm, y0 "
                  << 10 * l.y0 << " mm, dx/dz " << 1e3 * l.dxdz << " mrad, dy/dz " << 1e3 * l.dydz << " mrad" << std::endl;
      }
    }
    else
    {
      std::cout << Name() << ":   beam-axis alignment off" << std::endl;
    }
  }
  return createNodes(topNode);
}

int SiTrajectoryFitter::createNodes(PHCompositeNode* topNode)
{
  m_container = findNode::getClass<Si_TrajectoryContainer>(topNode, m_outputNodeName);
  if (m_container)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << ": DST node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator dstIter(dstNode);
  auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_container = new Si_TrajectoryContainerv1();
  svtxNode->addNode(new PHIODataNode<PHObject>(m_container, m_outputNodeName, "PHObject"));

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": created DST/SVTX/" << m_outputNodeName << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int SiTrajectoryFitter::process_event(PHCompositeNode* topNode)
{
  auto* event = findNode::getClass<SiHitSeedEvent>(topNode, m_inputNodeName);
  if (!event)
  {
    std::cerr << Name() << ": missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!m_container)
  {
    m_container = findNode::getClass<Si_TrajectoryContainer>(topNode, m_outputNodeName);
    if (!m_container)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  m_container->Reset();

  std::array<unsigned int, 5> nStatus{};  // NotFitted, Ok, TooFewPoints, Degenerate, StraightLine
  unsigned int nLowPt = 0;
  for (const auto& chain : event->chains)
  {
    auto* trajectory = new Si_Trajectoryv1();
    trajectory->set_id(m_container->size());
    trajectory->set_chain_id(chain.id);

    for (int hid : chain.hit_ids)
    {
      if (hid < 0 || hid >= static_cast<int>(event->hits.size()))
      {
        continue;
      }
      const auto& h = event->hits[hid];
      if (!m_useIntt && h.layer > 2)
      {
        continue;
      }
      trajectory->add_hit_index(h.hitsetkey, h.hitkey);
    }

    const auto points = layerPoints(*event, chain.hit_ids);
    unsigned int nMvtx = 0;
    unsigned int nIntt = 0;
    for (const auto& p : points)
    {
      trajectory->add_point(p.layer, p.global.x, p.global.y, p.global.z, p.nhits);
      (p.layer <= 2 ? nMvtx : nIntt)++;
    }
    trajectory->set_n_mvtx(nMvtx);
    trajectory->set_n_intt(nIntt);

    fit(points, *trajectory);
    const int status = trajectory->get_fit_status();
    ++nStatus[std::clamp(status, 0, 4)];
    const bool lowPt = status == Si_Trajectory::Ok && trajectory->get_pt() < m_reportPt;
    nLowPt += lowPt;
    // Verbosity 2: every chain that is not fitted or is fitted with a low pt
    if (Verbosity() > 2 || (Verbosity() > 1 && (!trajectory->is_fitted() || lowPt)))
    {
      std::cout << Name() << ":   chain " << chain.id << " status " << status
                << " pt " << trajectory->get_pt() << " circle_rms " << trajectory->get_circle_rms() << " cm, layers:";
      for (const auto& p : points)
      {
        std::cout << " " << p.layer << "(" << p.nhits << " hits, r=" << std::hypot(p.global.x, p.global.y) << ")";
      }
      std::cout << std::endl;
      if (Verbosity() > 2)
      {
        trajectory->identify();
      }
    }
    m_container->add_trajectory(trajectory);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": event " << event->event << " chains=" << event->chains.size()
              << " circle=" << nStatus[Si_Trajectory::Ok] << " (pt<" << m_reportPt << " GeV: " << nLowPt << ")"
              << " straight=" << nStatus[Si_Trajectory::StraightLine]
              << " too_few_points=" << nStatus[Si_Trajectory::TooFewPoints]
              << " degenerate=" << nStatus[Si_Trajectory::Degenerate] << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

// hits -> detector frame -> global frame, then one centroid per layer (on the layer cylinder
// in the detector frame, so the shift is applied to a point that lies on the nominal layer).
std::vector<SiTrajectoryFitter::LayerPoint> SiTrajectoryFitter::layerPoints(
    const SiHitSeedEvent& event, const std::vector<int>& hitIds) const
{
  struct Sum
  {
    double x = 0, y = 0, z = 0;
    unsigned int n = 0;
  };
  std::map<int, Sum> sums;  // ordered by layer = ordered outward
  for (int hid : hitIds)
  {
    if (hid < 0 || hid >= static_cast<int>(event.hits.size()))
    {
      continue;
    }
    const auto& h = event.hits[hid];
    if (h.layer < 0 || h.layer >= SiDetectorFrame::kNLayers || (!m_useIntt && h.layer > 2))
    {
      continue;
    }
    // step 1: detector frame
    const auto det = m_frame.toDetector(h.layer, h.phi, h.z);
    auto& s = sums[h.layer];
    s.x += det.x;
    s.y += det.y;
    s.z += det.z;
    ++s.n;
  }

  std::vector<LayerPoint> out;
  out.reserve(sums.size());
  for (const auto& [layer, s] : sums)
  {
    const double r = SiDetectorFrame::layerRadius(layer);
    const double phi = std::atan2(s.y, s.x);
    const SiDetectorFrame::Point det{r * std::cos(phi), r * std::sin(phi), s.z / s.n};
    LayerPoint p;
    p.layer = layer;
    // step 2: detector shift (global frame); step 3: beam-axis alignment of its clamshell half
    p.global = m_alignment.toBeamAxis(m_alignment.siHalf(phi - m_frame.uphiRotation()), m_frame.toGlobal(det));
    p.nhits = s.n;
    out.push_back(p);
  }
  return out;
}

void SiTrajectoryFitter::fit(const std::vector<LayerPoint>& points, Si_Trajectory& t) const
{
  std::vector<SiTpcHelixFit::Point> pts;
  pts.reserve(points.size());
  for (const auto& p : points)
  {
    pts.push_back({p.global.x, p.global.y, p.global.z, 1.0, p.layer <= 2 ? m_zWeightMvtx : m_zWeightIntt, 0});
  }
  SiTpcHelixFit::Config cfg;
  cfg.beamX = m_beamX;
  cfg.beamY = m_beamY;
  cfg.bz = m_bz;
  cfg.minLeverArm = m_minCurvatureLeverArm;
  cfg.beamWeight = m_beamConstraint ? m_beamWeight : 0.0;
  const auto r = SiTpcHelixFit::fit(pts, cfg, m_minPoints);
  t.set_fit_status(r.status);
  if (!r.fitted())
  {
    return;
  }
  t.set_circle_x(r.cx);
  t.set_circle_y(r.cy);
  t.set_radius(r.R);
  t.set_helicity(r.helicity);
  t.set_circle_rms(r.circleRms);
  t.set_pca_x(r.pcaX);
  t.set_pca_y(r.pcaY);
  t.set_phi(r.phi);
  t.set_dca(r.dca);
  t.set_pt(r.pt);
  t.set_charge(r.charge);
  t.set_z0(r.z0);
  t.set_tanl(r.tanl);
  t.set_z_rms(r.zRms);
}

SiTrajectoryFitter::Circle SiTrajectoryFitter::fitCircleTaubin(const std::vector<double>& x,
                                                                const std::vector<double>& y)
{
  return fitCircleTaubin(x, y, std::vector<double>(std::min(x.size(), y.size()), 1.0));
}

SiTrajectoryFitter::Circle SiTrajectoryFitter::fitCircleTaubin(const std::vector<double>& x,
                                                                const std::vector<double>& y,
                                                                const std::vector<double>& w)
{
  const auto c = SiTpcHelixFit::fitCircleTaubin(x, y, w);
  Circle out;
  out.x = c.x;
  out.y = c.y;
  out.r = c.r;
  out.ok = c.ok;
  return out;
}
