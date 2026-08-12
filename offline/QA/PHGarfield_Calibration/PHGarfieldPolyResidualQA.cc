#include "PHGarfieldPolyResidualQA.h"

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <TH1D.h>
#include <TH2F.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cctype>
#include <format>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    const double pi = std::acos(-1.0);
    while (phi > pi)
    {
      phi -= 2.0 * pi;
    }
    while (phi <= -pi)
    {
      phi += 2.0 * pi;
    }
    return phi;
  }

  std::string encode_value(const double value)
  {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(4) << value;
    std::string out = stream.str();
    while (out.size() > 1 && out.back() == '0')
    {
      out.pop_back();
    }
    if (!out.empty() && out.back() == '.')
    {
      out.pop_back();
    }
    for (char& c : out)
    {
      if (c == '-')
      {
        c = 'm';
      }
      else if (c == '.')
      {
        c = 'p';
      }
      else if (!std::isalnum(static_cast<unsigned char>(c)))
      {
        c = '_';
      }
    }
    return out;
  }

  int region_from_layer(const unsigned int layer)
  {
    if (layer >= 7 && layer <= 22)
    {
      return 0;
    }
    if (layer >= 23 && layer <= 38)
    {
      return 1;
    }
    if (layer >= 39 && layer <= 54)
    {
      return 2;
    }
    return -1;
  }

  unsigned int cluster_layer(const Tpc_PolyCluster* cluster)
  {
    if (!cluster || cluster->size_hits() == 0)
    {
      return 0xffffffffU;
    }
    const Tpc_PolyCluster::HitIndex hit_index = cluster->get_hit_index(0);
    return TrkrDefs::getLayer(hit_index.first);
  }

  bool line_xy_at_z(const Tpc_PolyTrack* trk,
                    const double z,
                    const double arc_direction,
                    double& x_state,
                    double& y_state)
  {
    if (!trk || trk->get_fit_status() == 0 || !std::isfinite(z))
    {
      return false;
    }

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        std::fabs(pz) < 1.0e-12)
    {
      return false;
    }

    const double dz = z - z0;
    x_state = x0 + arc_direction * px / pz * dz;
    y_state = y0 + arc_direction * py / pz * dz;
    return std::isfinite(x_state) && std::isfinite(y_state);
  }

  bool project_track_to_z(const Tpc_PolyTrack* trk,
                          const double z,
                          const double magnetic_field_tesla,
                          const double arc_direction,
                          const bool use_straight_line,
                          double& x_state,
                          double& y_state)
  {
    if (!trk || trk->get_fit_status() == 0)
    {
      return false;
    }

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    const double charge = trk->get_charge();

    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(z))
    {
      return false;
    }

    if (use_straight_line || std::fabs(charge * magnetic_field_tesla) < 1.0e-12)
    {
      return line_xy_at_z(trk, z, arc_direction, x_state, y_state);
    }

    const double pt = std::hypot(px, py);
    if (pt <= 0.0 || std::fabs(pz) < 1.0e-12)
    {
      return false;
    }

    const double signed_radius = pt / (0.003 * charge * magnetic_field_tesla);
    const double radius = std::fabs(signed_radius);
    if (radius <= 0.0 || !std::isfinite(radius))
    {
      return false;
    }

    const double tx = px / pt;
    const double ty = py / pt;
    const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
    const double xc = x0 + sign * radius * ty;
    const double yc = y0 - sign * radius * tx;
    const double phi0 = std::atan2(y0 - yc, x0 - xc);
    const double dzds = pz / pt;
    if (std::fabs(dzds) < 1.0e-12)
    {
      return false;
    }

    const double arc = arc_direction * (z - z0) / dzds;
    const double phi = phi0 - sign * arc / radius;
    x_state = xc + radius * std::cos(phi);
    y_state = yc + radius * std::sin(phi);
    return std::isfinite(x_state) && std::isfinite(y_state);
  }

  bool line_z_at_radius(const Tpc_PolyTrack* trk,
                        const double target_r,
                        const double reference_z,
                        const double arc_direction,
                        double& z_state)
  {
    if (!trk || trk->get_fit_status() == 0 || !std::isfinite(target_r) || target_r <= 0.0)
    {
      return false;
    }

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        std::fabs(pz) < 1.0e-12)
    {
      return false;
    }

    const double ax = arc_direction * px / pz;
    const double ay = arc_direction * py / pz;
    const double a = ax * ax + ay * ay;
    const double b = 2.0 * (x0 * ax + y0 * ay);
    const double c = x0 * x0 + y0 * y0 - target_r * target_r;
    if (a < 1.0e-20)
    {
      return false;
    }

    const double disc = b * b - 4.0 * a * c;
    if (disc < -1.0e-8)
    {
      return false;
    }
    const double root = std::sqrt(std::max(0.0, disc));
    const double dz1 = (-b - root) / (2.0 * a);
    const double dz2 = (-b + root) / (2.0 * a);
    const double z1 = z0 + dz1;
    const double z2 = z0 + dz2;
    z_state = std::fabs(z1 - reference_z) <= std::fabs(z2 - reference_z) ? z1 : z2;
    return std::isfinite(z_state);
  }

  bool helix_z_at_radius(const Tpc_PolyTrack* trk,
                         const double target_r,
                         const double reference_z,
                         const double magnetic_field_tesla,
                         const double arc_direction,
                         double& z_state)
  {
    if (!trk || trk->get_fit_status() == 0)
    {
      return false;
    }

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    const double charge = trk->get_charge();
    const double pt = std::hypot(px, py);
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(charge) || pt <= 0.0 ||
        std::fabs(charge * magnetic_field_tesla) < 1.0e-12)
    {
      return false;
    }

    const double signed_radius = pt / (0.003 * charge * magnetic_field_tesla);
    const double radius = std::fabs(signed_radius);
    const double tx = px / pt;
    const double ty = py / pt;
    const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
    const double xc = x0 + sign * radius * ty;
    const double yc = y0 - sign * radius * tx;
    const double phi0 = std::atan2(y0 - yc, x0 - xc);
    const double dzds = pz / pt;
    const double center_r = std::hypot(xc, yc);
    if (!std::isfinite(target_r) || target_r <= 0.0 || center_r <= 1.0e-12)
    {
      return false;
    }

    if (center_r > target_r + radius || center_r < std::fabs(target_r - radius))
    {
      return false;
    }

    const double a = (target_r * target_r - radius * radius + center_r * center_r) / (2.0 * center_r);
    const double h2 = target_r * target_r - a * a;
    if (h2 < -1.0e-8)
    {
      return false;
    }

    const double h = std::sqrt(std::max(0.0, h2));
    const double ux = xc / center_r;
    const double uy = yc / center_r;

    double best_z = 0.0;
    double best_abs_dz = std::numeric_limits<double>::max();
    const double pi = std::acos(-1.0);
    for (int isign = -1; isign <= 1; isign += 2)
    {
      const double x = a * ux - static_cast<double>(isign) * h * uy;
      const double y = a * uy + static_cast<double>(isign) * h * ux;
      const double phi_on_circle = std::atan2(y - yc, x - xc);
      for (int k = -4; k <= 4; ++k)
      {
        const double dphi = phi_on_circle - phi0 + 2.0 * pi * static_cast<double>(k);
        const double arc = -sign * radius * dphi;
        const double z = z0 + arc_direction * dzds * arc;
        const double abs_dz = std::fabs(z - reference_z);
        if (abs_dz < best_abs_dz)
        {
          best_abs_dz = abs_dz;
          best_z = z;
        }
      }
    }

    z_state = best_z;
    return best_abs_dz != std::numeric_limits<double>::max() && std::isfinite(z_state);
  }

  double cluster_line_residual2(const Tpc_PolyTrack* poly_track,
                                const std::vector<const Tpc_PolyCluster*>& clusters,
                                const double magnetic_field_tesla,
                                const double arc_direction,
                                const bool use_straight_line)
  {
    if (!poly_track || clusters.empty())
    {
      return std::numeric_limits<double>::max();
    }

    double sum = 0.0;
    unsigned int n = 0;
    for (const Tpc_PolyCluster* cluster : clusters)
    {
      if (!cluster || !cluster->isValid())
      {
        continue;
      }
      const double cx = cluster->get_centroid_x();
      const double cy = cluster->get_centroid_y();
      const double cz = cluster->get_centroid_z();
      if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(cz))
      {
        continue;
      }

      double x = 0.0;
      double y = 0.0;
      if (!project_track_to_z(poly_track, cz, magnetic_field_tesla, arc_direction, use_straight_line, x, y))
      {
        continue;
      }

      const double dx = x - cx;
      const double dy = y - cy;
      sum += dx * dx + dy * dy;
      ++n;
    }

    return n > 0 ? sum / static_cast<double>(n) : std::numeric_limits<double>::max();
  }
}  // namespace

PHGarfieldPolyResidualQA::PHGarfieldPolyResidualQA(const std::string& name,
                                                   const double kEffSide0,
                                                   const double kEffSide1,
                                                   const double spaceChargeScaleSide0,
                                                   const double spaceChargeScaleSide1)
  : SubsysReco(name)
  , m_kEffSide0(kEffSide0)
  , m_kEffSide1(kEffSide1)
  , m_spaceChargeScaleSide0(spaceChargeScaleSide0)
  , m_spaceChargeScaleSide1(spaceChargeScaleSide1)
{
}

int PHGarfieldPolyResidualQA::InitRun(PHCompositeNode* /*topNode*/)
{
  createHistos();

  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  static constexpr std::array<const char*, 2> charge_name = {"qneg", "qpos"};

  for (int side = 0; side < 2; ++side)
  {
    for (int charge_bin = 0; charge_bin < 2; ++charge_bin)
    {
      for (int region = 0; region < 3; ++region)
      {
        m_hRPhiResidual[side][charge_bin][region] = dynamic_cast<TH1*>(hm->getHisto(
            getHistoName(std::format("h_rphi_residual_side{}_{}_R{}",
                                     side, charge_name[charge_bin], region + 1))));
      }

      m_hRPhiVsLayer[side][charge_bin] = dynamic_cast<TH2*>(hm->getHisto(
          getHistoName(std::format("h_rphi_vs_layer_side{}_{}", side, charge_name[charge_bin]))));
    }

    for (int region = 0; region < 3; ++region)
    {
      m_hZResidual[side][region] = dynamic_cast<TH1*>(hm->getHisto(
          getHistoName(std::format("h_z_residual_side{}_R{}", side, region + 1))));
    }

    m_hZVsLayer[side] = dynamic_cast<TH2*>(hm->getHisto(
        getHistoName(std::format("h_z_vs_layer_side{}", side))));
  }

  bool missing_histogram = false;
  for (int side = 0; side < 2; ++side)
  {
    for (int charge_bin = 0; charge_bin < 2; ++charge_bin)
    {
      for (int region = 0; region < 3; ++region)
      {
        missing_histogram = missing_histogram || !m_hRPhiResidual[side][charge_bin][region];
      }
      missing_histogram = missing_histogram || !m_hRPhiVsLayer[side][charge_bin];
    }
    for (int region = 0; region < 3; ++region)
    {
      missing_histogram = missing_histogram || !m_hZResidual[side][region];
    }
    missing_histogram = missing_histogram || !m_hZVsLayer[side];
  }
  if (missing_histogram)
  {
    std::cout << PHWHERE << " " << Name() << " failed to retrieve one or more registered histograms" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHGarfieldPolyResidualQA::getNodes(PHCompositeNode* topNode)
{
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_clusterNodeName);
  if (!m_clusters)
  {
    std::cout << PHWHERE << " Missing " << m_clusterNodeName << std::endl;
    return false;
  }

  m_tracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_trackNodeName);
  if (!m_tracks)
  {
    std::cout << PHWHERE << " Missing " << m_trackNodeName << std::endl;
    return false;
  }

  return true;
}

std::string PHGarfieldPolyResidualQA::coefficientSuffix() const
{
  return std::format("KEffSide0_{}_KEffSide1_{}_SpaceChargeScaleSide0_{}_SpaceChargeScaleSide1_{}",
                     encode_value(m_kEffSide0),
                     encode_value(m_kEffSide1),
                     encode_value(m_spaceChargeScaleSide0),
                     encode_value(m_spaceChargeScaleSide1));
}

std::string PHGarfieldPolyResidualQA::getHistoName(const std::string& base) const
{
  return std::format("{}_{}", base, coefficientSuffix());
}

void PHGarfieldPolyResidualQA::createHistos()
{
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  static constexpr std::array<const char*, 2> charge_name = {"qneg", "qpos"};

  for (int side = 0; side < 2; ++side)
  {
    for (int charge_bin = 0; charge_bin < 2; ++charge_bin)
    {
      for (int region = 0; region < 3; ++region)
      {
        const std::string name = getHistoName(std::format("h_rphi_residual_side{}_{}_R{}",
                                                          side, charge_name[charge_bin], region + 1));
        auto* hist = new TH1D(name.c_str(),
                              std::format("TPC polycluster r#phi residual side {} {} R{};r#Delta#phi [cm];Entries",
                                          side, charge_name[charge_bin], region + 1)
                                  .c_str(),
                              200, -1.0, 1.0);
        hm->registerHisto(hist);
        std::cout<<"PHGarfieldPolyResidualQA:: registerHisto "<<hist<<std::endl;
      }

      const std::string name = getHistoName(std::format("h_rphi_vs_layer_side{}_{}", side, charge_name[charge_bin]));
      auto* hist = new TH2F(name.c_str(),
                            std::format("TPC polycluster r#phi residual vs layer side {} {};Layer;r#Delta#phi [cm]",
                                        side, charge_name[charge_bin])
                                .c_str(),
                            48, 6.5, 54.5, 200, -1.0, 1.0);
      hm->registerHisto(hist);
       std::cout<<"PHGarfieldPolyResidualQA:: registerHisto "<<hist<<std::endl;
    }

    for (int region = 0; region < 3; ++region)
    {
      const std::string name = getHistoName(std::format("h_z_residual_side{}_R{}", side, region + 1));
      auto* hist = new TH1D(name.c_str(),
                            std::format("TPC polycluster z residual side {} R{};#Deltaz [cm];Entries",
                                        side, region + 1)
                                .c_str(),
                            200, -2.0, 2.0);
      hm->registerHisto(hist);
       std::cout<<"PHGarfieldPolyResidualQA:: registerHisto "<<hist<<std::endl;
    }

    const std::string name = getHistoName(std::format("h_z_vs_layer_side{}", side));
    auto* hist = new TH2F(name.c_str(),
                          std::format("TPC polycluster z residual vs layer side {};Layer;#Deltaz [cm]", side).c_str(),
                          48, 6.5, 54.5, 200, -2.0, 2.0);
    hm->registerHisto(hist);
     std::cout<<"PHGarfieldPolyResidualQA:: registerHisto "<<hist<<std::endl;
  }
}

int PHGarfieldPolyResidualQA::process_event(PHCompositeNode* topNode)
{
  if (!m_clusters || !m_tracks)
  {
    if (!getNodes(topNode))
    {

      return Fun4AllReturnCodes::EVENT_OK;
    }
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

  for (unsigned int itrack = 0; itrack < m_tracks->size(); ++itrack)
  {
    const Tpc_PolyTrack* poly_track = m_tracks->get_track(itrack);
    if (!poly_track || !poly_track->isValid() || poly_track->get_fit_status() == 0)
    {
      continue;
    }

    const auto cluster_iter = clusters_by_source_assembled_track_id.find(poly_track->get_source_assembled_track_id());
    if (cluster_iter == clusters_by_source_assembled_track_id.end())
    {
      continue;
    }

    const std::vector<const Tpc_PolyCluster*>& track_clusters = cluster_iter->second;
    const double charge = poly_track->get_charge();
    if (!std::isfinite(charge) || charge == 0.0)
    {
      continue;
    }
    const int charge_bin = charge > 0.0 ? 1 : 0;
    const bool use_straight_line = m_useStraightLineTracks || std::fabs(charge * m_magneticFieldTesla) < 1.0e-12;

    const double forward_residual2 = cluster_line_residual2(poly_track, track_clusters, m_magneticFieldTesla, 1.0, use_straight_line);
    const double reverse_residual2 = cluster_line_residual2(poly_track, track_clusters, m_magneticFieldTesla, -1.0, use_straight_line);
    const double arc_direction = forward_residual2 <= reverse_residual2 ? 1.0 : -1.0;

    for (const Tpc_PolyCluster* cluster : track_clusters)
    {
      if (!cluster || !cluster->isValid())
      {
        continue;
      }

      const int side = cluster->get_side();
      if (side < 0 || side > 1)
      {
        continue;
      }

      const unsigned int layer = cluster_layer(cluster);
      const int region = region_from_layer(layer);
      if (region < 0)
      {
        continue;
      }

      const double cluster_x = cluster->get_centroid_x();
      const double cluster_y = cluster->get_centroid_y();
      const double cluster_z = cluster->get_centroid_z();
      if (!std::isfinite(cluster_x) || !std::isfinite(cluster_y) || !std::isfinite(cluster_z))
      {
        continue;
      }

      double state_x = 0.0;
      double state_y = 0.0;
      if (!project_track_to_z(poly_track, cluster_z, m_magneticFieldTesla, arc_direction, use_straight_line, state_x, state_y))
      {
        continue;
      }

      const double cluster_r = std::hypot(cluster_x, cluster_y);
      const double cluster_phi = std::atan2(cluster_y, cluster_x);
      const double state_phi = std::atan2(state_y, state_x);
      const double residual_rphi = cluster_r * wrap_phi(cluster_phi - state_phi);

      m_hRPhiResidual[side][charge_bin][region]->Fill(residual_rphi);
      m_hRPhiVsLayer[side][charge_bin]->Fill(layer, residual_rphi);

      double state_z = 0.0;
      const bool have_state_z = use_straight_line ? line_z_at_radius(poly_track, cluster_r, cluster_z, arc_direction, state_z) : helix_z_at_radius(poly_track, cluster_r, cluster_z, m_magneticFieldTesla, arc_direction, state_z);
      if (have_state_z)
      {
        const double residual_z = cluster_z - state_z;
        m_hZResidual[side][region]->Fill(residual_z);
        m_hZVsLayer[side]->Fill(layer, residual_z);
      }
    }
  }


  return Fun4AllReturnCodes::EVENT_OK;
}
