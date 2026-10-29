#include "SiTpcActsFitter.h"

#include "SiTpc_Track.h"
#include "SiTpc_TrackContainer.h"
#include "Si_Trajectory.h"
#include "Si_TrajectoryContainer.h"

#include <tpctrackreco/TpcCrossingDecision.h>
#include <tpctrackreco/TpcCrossingDecisionContainer.h>
#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSourceLink.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/Calibrator.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/alignmentTransformationContainer.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrackState_v3.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Material/ISurfaceMaterial.hpp>
#include <Acts/Material/MaterialSlab.hpp>
#include <Acts/Propagator/DirectNavigator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/TrackHelpers.hpp>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  using Stepper = Acts::EigenStepper<>;
  using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
  using DirectFitter = Acts::KalmanFitter<DirectPropagator, Acts::VectorMultiTrajectory>;
  using SurfacePtr = std::shared_ptr<const Acts::Surface>;
  using TrackContainer = ActsTrackFittingAlgorithm::TrackContainer;
  using MeasurementContainer = ActsTrackFittingAlgorithm::MeasurementContainer;

  constexpr double cm = Acts::UnitConstants::cm;  // 10 (Acts length unit = mm)
  constexpr double twopi = 2. * M_PI;

  // detector index used in the QA trees
  constexpr int kMvtx = 0;
  constexpr int kIntt = 1;
  constexpr int kTpc = 2;

  // geometry id volume of the artificial TPC planes (not used by the sPHENIX tracking geometry)
  constexpr Acts::GeometryIdentifier::Value kArtificialVolume = 250;

  // tie margin [cm]: material surfaces at the same path length as a measurement go after it
  // (the direct navigator silently skips surfaces behind the track, see PHActsTrkFitter)
  constexpr double kMaterialTieMargin = 0.01;

  // window of the position match chain point - official cluster [cm] (local 0, local 1)
  constexpr std::array<double, 2> kMvtxMatchWindow{0.03, 0.03};
  constexpr std::array<double, 2> kInttMatchWindow{0.05, 1.5};

  template <class T>
  T square(T x)
  {
    return x * x;
  }

  int tpcModule(unsigned int layer)
  {
    if (layer < 7)
    {
      return 0;
    }
    return std::min(2, static_cast<int>((layer - 7) / 16));
  }

  // one measurement: Si on an official surface, TPC on an artificial plane
  struct Measurement
  {
    int det = kTpc;
    unsigned int layer = 0;
    TrkrDefs::cluskey key = TrkrDefs::CLUSKEYMAX;
    SurfacePtr surface;
    Acts::ActsVector<2> local = Acts::ActsVector<2>::Zero();          // [mm]
    Acts::ActsSquareMatrix<2> cov = Acts::ActsSquareMatrix<2>::Zero();  // [mm^2]
    Acts::Vector3 global = Acts::Vector3::Zero();                       // [cm], fit frame
    double s = 0;                                                       // path length in xy along the guess [cm]
    bool official = true;     // Si: official TRKR_CLUSTER cluster (false: projected chain point)
    double crossingAngle = 0;  // TPC: track direction vs plane normal [rad]
  };

  // circle + straight line in (s, z) in the fit frame, built from the own Si+TPC fit
  struct HelixGuess
  {
    double cx = 0, cy = 0, R = 0;
    double omega = 1;  // +1: counter-clockwise seen from +z
    double x0 = 0, y0 = 0, z0 = 0, tanl = 0;

    double angleOf(double x, double y) const { return std::atan2(y - cy, x - cx); }

    // path length in xy from the start point, along the direction of motion [cm]
    double arcLength(double x, double y) const
    {
      double d = std::remainder(omega * (angleOf(x, y) - angleOf(x0, y0)), twopi);
      if (d < -0.1)
      {
        d += twopi;
      }
      return R * d;
    }

    // first crossing of the circle with a cylinder of radius r around (ox, oy); NaN if none
    double arcLengthToRadius(double r, double ox, double oy) const
    {
      const double dx = cx - ox, dy = cy - oy;
      const double d = std::hypot(dx, dy);
      if (d <= 0 || d > r + R || d < std::abs(r - R))
      {
        return std::numeric_limits<double>::quiet_NaN();
      }
      const double a = (r * r - R * R + d * d) / (2 * d);
      const double h = std::sqrt(std::max(0., r * r - a * a));
      const double ux = dx / d, uy = dy / d;
      double best = std::numeric_limits<double>::quiet_NaN();
      for (const double sign : {-1., 1.})
      {
        const double px = ox + a * ux - sign * h * uy;
        const double py = oy + a * uy + sign * h * ux;
        const double s = arcLength(px, py);
        if (s >= 0 && (std::isnan(best) || s < best))
        {
          best = s;
        }
      }
      return best;
    }

    double zAt(double s) const { return z0 + tanl * s; }

    // unit direction of motion at a point on (or near) the circle
    Acts::Vector3 direction(double x, double y) const
    {
      const double rx = (x - cx) / R, ry = (y - cy) / R;
      Acts::Vector3 d(-omega * ry, omega * rx, tanl);
      return d.normalized();
    }
  };

  // KF source link -> surface (the artificial planes are not in the tracking geometry)
  // KF outlier rejection: a measurement whose predicted chi2 (before the update) exceeds the cut of its
  // detector is not used (kept as an outlier state).  Cut <= 0: off.
  struct ChiOutlierFinder
  {
    double tpcCut = 0;
    double siCut = 0;
    bool operator()(Acts::VectorMultiTrajectory::ConstTrackStateProxy state) const
    {
      if (!state.hasCalibrated() || !state.hasPredicted())
      {
        return false;
      }
      const auto key = state.getUncalibratedSourceLink().template get<ActsSourceLink>().cluskey();
      const double cut = TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId ? tpcCut : siCut;
      if (cut <= 0)
      {
        return false;
      }
      return Acts::calculatePredictedChi2(state) > cut;
    }
  };

  struct IndexSurfaceAccessor
  {
    const std::vector<const Acts::Surface*>* surfaces = nullptr;
    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const
    {
      const auto& sl = sourceLink.get<ActsSourceLink>();
      return surfaces->at(sl.index());
    }
  };

  // one fit: measurements, source links, surface sequence and the result
  struct FitJob
  {
    MeasurementContainer measurements;
    std::vector<Acts::SourceLink> sourceLinks;
    std::vector<const Acts::Surface*> slSurfaces;  // by source link index
    std::vector<const Measurement*> slMeasurements;
    std::vector<const Acts::Surface*> sequence;
    std::shared_ptr<Acts::VectorTrackContainer> trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    std::shared_ptr<Acts::VectorMultiTrajectory> trackStates = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tracks{trackContainer, trackStates};
    std::optional<TrackContainer::TrackProxy> track;
    int nMaterial = 0;
    std::vector<std::pair<double, const Acts::Surface*>> ordered;  // navigation sequence with path length (debug)
    std::string failReason;
  };

  // perigee parameters of a fit, for the QA tree
  struct PerigeeSummary
  {
    int ok = 0;
    float d0 = NAN, z0 = NAN, phi = NAN, theta = NAN, pt = NAN;
    float d0err = NAN, z0err = NAN, chi2 = NAN;
    int ndf = 0;
  };
}  // namespace

//____________________________________________________________________________________________________
struct SiTpcActsFitter::Impl
{
  // nodes
  ActsGeometry* geometry = nullptr;
  alignmentTransformationContainer* transientMap = nullptr;
  SiTpc_TrackContainer* siTpcTracks = nullptr;
  Si_TrajectoryContainer* siTrajectories = nullptr;
  Tpc_PolyClusterContainer* tpcClusters = nullptr;
  TpcCrossingDecisionContainer* crossings = nullptr;
  TrkrClusterContainer* siClusters = nullptr;
  TrkrClusterHitAssoc* hitAssoc = nullptr;
  SvtxTrackMap* trackMap = nullptr;
  TrackSeedContainer* siSeeds = nullptr;
  TrackSeedContainer* tpcSeeds = nullptr;

  // fit machinery
  Acts::GeometryContext gctx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  std::unique_ptr<const Acts::Logger> logger;
  std::unique_ptr<DirectFitter> fitter;
  std::unique_ptr<const Acts::Logger> verboseLogger;
  std::unique_ptr<DirectFitter> verboseFitter;  // debug of failed fits (setDebugFailedFits)
  int nDebugged = 0;
  unsigned long nTpcErrStd = 0, nTpcErrFallback = 0;
  Acts::GainMatrixUpdater updater;
  Acts::GainMatrixSmoother smoother;
  Calibrator calibrator;
  std::vector<const Acts::Surface*> materialSurfaces;  // non-sensitive official material surfaces
  double bzSign = 1;  // sign of B_z at the centre: charge = -omega * sign(B_z)

  // INTT surfaces moved for this event: id -> transient transform before the move
  std::map<Acts::GeometryIdentifier, Acts::Transform3> movedIntt;
  // hitsetkey -> (hitkey -> official cluster key), filled on demand per event
  std::unordered_map<TrkrDefs::hitsetkey, std::unordered_map<TrkrDefs::hitkey, TrkrDefs::cluskey>> hitToCluster;

  // counters
  unsigned long nEvents = 0, nInput = 0, nSelected = 0, nFitOk = 0, nFitFail = 0;
  unsigned long nSiOfficialAssoc = 0, nSiOfficialPosition = 0, nSiProjected = 0, nSiMissing = 0, nSiOtherStrobe = 0;
  unsigned long nTpc = 0, nTpcAngleCut = 0, nTpcMaterial = 0, nChargeFlip = 0;
  unsigned long nRecoveredNoMaterial = 0, nOutlierTpc = 0, nOutlierSi = 0;
  std::map<std::string, unsigned long> rejected;
  // failed fits: reason -> counts per pT bin (< 0.5, 0.5 - 1, 1 - 2, > 2 GeV)
  std::map<std::string, std::array<unsigned long, 4>> failReasons;

  // QA
  TFile* qaFile = nullptr;
  TTree* trackTree = nullptr;
  TTree* stateTree = nullptr;
  struct TrackRow
  {
    int event = 0, track = 0, crossing = 0, charge = 0;
    int nmvtx = 0, nintt = 0, ntpc = 0, nsiProjected = 0, nmaterial = 0, noutliers = 0, retry_no_material = 0;
    float own_pt = NAN, own_phi = NAN, own_tanl = NAN, own_z0 = NAN, own_tpcZOffset = NAN;
    PerigeeSummary full, tpc, si;
    float beamX = NAN, beamY = NAN;
  } trackRow;
  struct StateRow
  {
    int event = 0, track = 0, det = 0, layer = 0, official = 0;
    float gx = NAN, gy = NAN, gz = NAN;  // measurement [cm]
    float fx = NAN, fy = NAN, fz = NAN;  // smoothed track [cm]
    float res0 = NAN, res1 = NAN;        // measurement - smoothed, local [cm]
    float sigMeas0 = NAN, sigMeas1 = NAN, sigFit0 = NAN, sigFit1 = NAN;  // [cm]
    float pt = NAN, s = NAN, crossingAngle = NAN;
  } stateRow;

  void moveInttSurface(const Acts::Surface& surface, const std::array<double, 2>& offset);
  void restoreIntt();
  TrkrDefs::cluskey clusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk);
  // official cluster of a chain hit.  MVTX: MvtxClusterPruner removes the copy of a cluster duplicated in the
  // next strobe (identical pixels) from TRKR_CLUSTER without updating TRKR_CLUSTERHITASSOC, so a chain hit
  // can point to a removed cluster: then the cluster with the same pixel in the neighbouring strobe is used.
  TrkrDefs::cluskey officialClusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, bool& otherStrobe);
};

//____________________________________________________________________________________________________
void SiTpcActsFitter::Impl::moveInttSurface(const Acts::Surface& surface, const std::array<double, 2>& offset)
{
  const auto id = surface.geometryId();
  if (movedIntt.contains(id))
  {
    return;
  }
  const Acts::Transform3 before = transientMap->getTransform(id);
  const Acts::Transform3 shift = geometry->makeAffineTransform(Acts::Vector3(0, 0, 0), Acts::Vector3(offset[0] * cm, offset[1] * cm, 0));
  transientMap->replaceTransform(id, shift * before);
  movedIntt.emplace(id, before);
}

//____________________________________________________________________________________________________
void SiTpcActsFitter::Impl::restoreIntt()
{
  for (const auto& [id, transform] : movedIntt)
  {
    transientMap->replaceTransform(id, transform);
  }
  movedIntt.clear();
}

//____________________________________________________________________________________________________
TrkrDefs::cluskey SiTpcActsFitter::Impl::clusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk)
{
  if (!hitAssoc)
  {
    return TrkrDefs::CLUSKEYMAX;
  }
  auto it = hitToCluster.find(hsk);
  if (it == hitToCluster.end())
  {
    std::unordered_map<TrkrDefs::hitkey, TrkrDefs::cluskey> m;
    if (auto* map = hitAssoc->getClusterMap(hsk))
    {
      for (const auto& [ckey, hitkey] : *map)
      {
        m.emplace(hitkey, ckey);
      }
    }
    it = hitToCluster.emplace(hsk, std::move(m)).first;
  }
  const auto c = it->second.find(hk);
  return c == it->second.end() ? TrkrDefs::CLUSKEYMAX : c->second;
}

//____________________________________________________________________________________________________
TrkrDefs::cluskey SiTpcActsFitter::Impl::officialClusterOfHit(TrkrDefs::hitsetkey hsk, TrkrDefs::hitkey hk, bool& otherStrobe)
{
  otherStrobe = false;
  const auto ckey = clusterOfHit(hsk, hk);
  if (ckey != TrkrDefs::CLUSKEYMAX && siClusters->findCluster(ckey))
  {
    return ckey;
  }
  if (TrkrDefs::getTrkrId(hsk) != TrkrDefs::mvtxId)
  {
    return TrkrDefs::CLUSKEYMAX;
  }
  const auto layer = TrkrDefs::getLayer(hsk);
  const auto stave = MvtxDefs::getStaveId(hsk);
  const auto chip = MvtxDefs::getChipId(hsk);
  const int strobe = MvtxDefs::getStrobeId(hsk);
  for (const int ds : {-1, 1})
  {
    const auto other = clusterOfHit(MvtxDefs::genHitSetKey(layer, stave, chip, strobe + ds), hk);
    if (other != TrkrDefs::CLUSKEYMAX && siClusters->findCluster(other))
    {
      otherStrobe = true;
      return other;
    }
  }
  return TrkrDefs::CLUSKEYMAX;
}

//____________________________________________________________________________________________________
SiTpcActsFitter::SiTpcActsFitter(const std::string& name)
  : SubsysReco(name)
  , m_impl(std::make_unique<Impl>())
{
}

SiTpcActsFitter::~SiTpcActsFitter() = default;

//____________________________________________________________________________________________________
int SiTpcActsFitter::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto& d = *m_impl;
  d.gctx = Acts::GeometryContext(d.transientMap);

  // direct-navigation Kalman fitter with our own source link -> surface accessor
  const auto level = Verbosity() > 5 ? Acts::Logging::VERBOSE : Acts::Logging::FATAL;
  d.logger = Acts::getDefaultLogger(Name(), level);
  Stepper stepper(d.geometry->geometry().magField);
  Acts::DirectNavigator navigator{d.logger->cloneWithSuffix("DirectNavigator")};
  DirectPropagator propagator(std::move(stepper), std::move(navigator), d.logger->cloneWithSuffix("DirectPropagator"));
  d.fitter = std::make_unique<DirectFitter>(std::move(propagator), d.logger->cloneWithSuffix("DirectFitter"));
  if (m_debugFailedFits > 0)
  {
    d.verboseLogger = Acts::getDefaultLogger(Name() + "Debug", Acts::Logging::VERBOSE);
    Stepper vstepper(d.geometry->geometry().magField);
    Acts::DirectNavigator vnavigator{d.verboseLogger->cloneWithSuffix("DirectNavigator")};
    DirectPropagator vpropagator(std::move(vstepper), std::move(vnavigator), d.verboseLogger->cloneWithSuffix("DirectPropagator"));
    d.verboseFitter = std::make_unique<DirectFitter>(std::move(vpropagator), d.verboseLogger->cloneWithSuffix("DirectFitter"));
  }

  // field direction: the charge given to the fit follows the bending of the own track
  {
    const auto& field = *d.geometry->geometry().magField;
    auto cache = field.makeCache(d.geometry->geometry().magFieldContext);
    const auto b = field.getField(Acts::Vector3(0, 0, 0), cache);
    d.bzSign = (b.ok() && b.value().z() < 0) ? -1. : 1.;
    std::cout << Name() << "::InitRun - B_z at the centre " << (b.ok() ? b.value().z() / Acts::UnitConstants::T : NAN) << " T" << std::endl;
  }

  // non-sensitive official material surfaces (beam pipe, layer / approach surfaces, field cages, ...).
  // Material on sensitive Si surfaces comes with the measurement surfaces, on the TPC with the planes.
  d.materialSurfaces.clear();
  if (m_useMaterialSurfaces)
  {
    MaterialSurfaceSelector selector;
    d.geometry->geometry().tGeometry->visitSurfaces(selector, false);
    for (const auto* surface : selector.surfaces)
    {
      if (surface->geometryId().sensitive() == 0)
      {
        d.materialSurfaces.push_back(surface);
      }
    }
  }

  if (m_inttOffset[0] != 0 || m_inttOffset[1] != 0)
  {
    if (!alignmentTransformationContainer::use_alignment)
    {
      std::cout << Name() << "::InitRun - WARNING alignment transforms are not used by the geometry: the INTT offset is"
                << " applied in the sensor plane only (no out-of-plane component)" << std::endl;
    }
  }

  if (!m_qaFileName.empty())
  {
    d.qaFile = new TFile(m_qaFileName.c_str(), "RECREATE");
    d.trackTree = new TTree("tracks", "SiTpcActsFitter tracks (d0, z0 relative to the beam axis at (beamX, beamY))");
    auto& t = d.trackRow;
    d.trackTree->Branch("event", &t.event);
    d.trackTree->Branch("track", &t.track);
    d.trackTree->Branch("crossing", &t.crossing);
    d.trackTree->Branch("charge", &t.charge);
    d.trackTree->Branch("nmvtx", &t.nmvtx);
    d.trackTree->Branch("nintt", &t.nintt);
    d.trackTree->Branch("ntpc", &t.ntpc);
    d.trackTree->Branch("nsi_projected", &t.nsiProjected);
    d.trackTree->Branch("nmaterial", &t.nmaterial);
    d.trackTree->Branch("noutliers", &t.noutliers);
    d.trackTree->Branch("retry_no_material", &t.retry_no_material);
    d.trackTree->Branch("own_pt", &t.own_pt);
    d.trackTree->Branch("own_phi", &t.own_phi);
    d.trackTree->Branch("own_tanl", &t.own_tanl);
    d.trackTree->Branch("own_z0", &t.own_z0);
    d.trackTree->Branch("own_tpc_z_offset", &t.own_tpcZOffset);
    d.trackTree->Branch("beamX", &t.beamX);
    d.trackTree->Branch("beamY", &t.beamY);
    for (auto& [prefix, p] : std::array<std::pair<std::string, PerigeeSummary*>, 3>{{{"", &t.full}, {"tpc_", &t.tpc}, {"si_", &t.si}}})
    {
      d.trackTree->Branch((prefix + "ok").c_str(), &p->ok);
      d.trackTree->Branch((prefix + "d0").c_str(), &p->d0);
      d.trackTree->Branch((prefix + "z0").c_str(), &p->z0);
      d.trackTree->Branch((prefix + "phi").c_str(), &p->phi);
      d.trackTree->Branch((prefix + "theta").c_str(), &p->theta);
      d.trackTree->Branch((prefix + "pt").c_str(), &p->pt);
      d.trackTree->Branch((prefix + "d0err").c_str(), &p->d0err);
      d.trackTree->Branch((prefix + "z0err").c_str(), &p->z0err);
      d.trackTree->Branch((prefix + "chi2").c_str(), &p->chi2);
      d.trackTree->Branch((prefix + "ndf").c_str(), &p->ndf);
    }

    d.stateTree = new TTree("states", "SiTpcActsFitter measurements of the full fit (smoothed, biased residuals)");
    auto& s = d.stateRow;
    d.stateTree->Branch("event", &s.event);
    d.stateTree->Branch("track", &s.track);
    d.stateTree->Branch("det", &s.det);
    d.stateTree->Branch("layer", &s.layer);
    d.stateTree->Branch("official", &s.official);
    d.stateTree->Branch("gx", &s.gx);
    d.stateTree->Branch("gy", &s.gy);
    d.stateTree->Branch("gz", &s.gz);
    d.stateTree->Branch("fx", &s.fx);
    d.stateTree->Branch("fy", &s.fy);
    d.stateTree->Branch("fz", &s.fz);
    d.stateTree->Branch("res0", &s.res0);
    d.stateTree->Branch("res1", &s.res1);
    d.stateTree->Branch("sig_meas0", &s.sigMeas0);
    d.stateTree->Branch("sig_meas1", &s.sigMeas1);
    d.stateTree->Branch("sig_fit0", &s.sigFit0);
    d.stateTree->Branch("sig_fit1", &s.sigFit1);
    d.stateTree->Branch("pt", &s.pt);
    d.stateTree->Branch("s", &s.s);
    d.stateTree->Branch("crossing_angle", &s.crossingAngle);
  }

  std::cout << Name() << "::InitRun - configuration" << std::endl
            << "  Si beam spot (" << m_siBeam[0] << ", " << m_siBeam[1] << ") cm, TPC beam line (" << m_tpcLine[0] << ", "
            << m_tpcLine[1] << ") cm, Si acts z shift " << m_siActsZShift << " cm" << std::endl
            << "  TPC frame " << (m_tpcFrame == TpcFrame::BeamAxis ? "BeamAxis" : "OfficialWorld")
            << ", TPC extra shift (" << m_tpcExtraShift[0] << ", " << m_tpcExtraShift[1] << ", " << m_tpcExtraShift[2] << ") cm"
            << ", per-track TPC z offset " << (m_applyTrackTpcZOffset ? "applied" : "not applied") << std::endl
            << "  INTT offset (" << m_inttOffset[0] << ", " << m_inttOffset[1] << ") cm, Si measurements "
            << (m_siMeasurement == SiMeasurement::OfficialClusters ? "official clusters" : "projected chain points") << std::endl
            << "  TPC errors rphi/z [cm]: " << m_tpcErr[0][0] << "/" << m_tpcErr[0][1] << ", " << m_tpcErr[1][0] << "/" << m_tpcErr[1][1]
            << ", " << m_tpcErr[2][0] << "/" << m_tpcErr[2][1] << ", max crossing angle " << m_maxTpcCrossingAngle << " rad" << std::endl
            << "  material surfaces " << d.materialSurfaces.size() << ", copy TPC material " << m_copyTpcMaterial
            << ", multiple scattering " << m_multipleScattering << ", energy loss " << m_energyLoss << std::endl
            << "  output " << m_trackMapName << ", " << m_siSeedName << ", " << m_tpcSeedName
            << (m_qaFileName.empty() ? "" : ", QA " + m_qaFileName) << (m_qaSubFits ? " (with TPC-only and Si-only fits)" : "")
            << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcActsFitter::getNodes(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
  d.geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  d.transientMap = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainerTransient");
  d.siTpcTracks = findNode::getClass<SiTpc_TrackContainer>(topNode, m_siTpcTrackNodeName);
  d.siTrajectories = findNode::getClass<Si_TrajectoryContainer>(topNode, m_siTrajectoryNodeName);
  d.tpcClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_tpcClusterNodeName);
  d.crossings = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingNodeName);
  d.siClusters = findNode::getClass<TrkrClusterContainer>(topNode, m_siClusterNodeName);
  d.hitAssoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  bool ok = true;
  auto require = [&](const void* p, const std::string& what)
  {
    if (!p)
    {
      std::cout << PHWHERE << " " << Name() << ": missing " << what << std::endl;
      ok = false;
    }
  };
  require(d.geometry, "ActsGeometry");
  require(d.transientMap, "alignmentTransformationContainerTransient");
  require(d.siTpcTracks, m_siTpcTrackNodeName);
  require(d.siTrajectories, m_siTrajectoryNodeName);
  require(d.tpcClusters, m_tpcClusterNodeName);
  require(d.siClusters, m_siClusterNodeName);
  if (!d.crossings)
  {
    std::cout << Name() << ": no " << m_crossingNodeName << " - crossing 0 for every track" << std::endl;
  }
  if (!d.hitAssoc && m_siMeasurement == SiMeasurement::OfficialClusters)
  {
    std::cout << Name() << ": no TRKR_CLUSTERHITASSOC - Si clusters matched by position only" << std::endl;
  }
  return ok ? Fun4AllReturnCodes::EVENT_OK : Fun4AllReturnCodes::ABORTRUN;
}

//____________________________________________________________________________________________________
int SiTpcActsFitter::createNodes(PHCompositeNode* topNode)
{
  auto& d = *m_impl;
  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " " << Name() << ": DST node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator dstIter(dstNode);
  auto* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  d.trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!d.trackMap)
  {
    d.trackMap = new SvtxTrackMap_v2;
    svtxNode->addNode(new PHIODataNode<PHObject>(d.trackMap, m_trackMapName, "PHObject"));
  }
  d.siSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_siSeedName);
  if (!d.siSeeds)
  {
    d.siSeeds = new TrackSeedContainer_v1;
    svtxNode->addNode(new PHIODataNode<PHObject>(d.siSeeds, m_siSeedName, "PHObject"));
  }
  d.tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_tpcSeedName);
  if (!d.tpcSeeds)
  {
    d.tpcSeeds = new TrackSeedContainer_v1;
    svtxNode->addNode(new PHIODataNode<PHObject>(d.tpcSeeds, m_tpcSeedName, "PHObject"));
  }

  // Downstream modules (PHSimpleVertexFinder, ...) require seed containers with the standard names in InitRun,
  // although for SvtxTrack_v4 they use the seed pointers stored on the track (our containers above).
  // Create empty placeholders only if nothing else (e.g. the standard seeding) made them.
  for (const std::string& name : {std::string("SiliconTrackSeedContainer"), std::string("TpcTrackSeedContainer")})
  {
    if (name != m_siSeedName && name != m_tpcSeedName && !findNode::getClass<TrackSeedContainer>(topNode, name))
    {
      svtxNode->addNode(new PHIODataNode<PHObject>(new TrackSeedContainer_v1, name, "PHObject"));
      std::cout << Name() << "::createNodes - created an empty " << name << " (placeholder for downstream modules)" << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcActsFitter::process_event(PHCompositeNode* /*topNode*/)
{
  auto& d = *m_impl;
  ++d.nEvents;
  d.trackMap->Reset();
  d.siSeeds->Reset();
  d.tpcSeeds->Reset();
  d.hitToCluster.clear();
  d.restoreIntt();

  const auto& magContext = d.geometry->geometry().magFieldContext;
  const auto& calibContext = d.geometry->geometry().calibContext;
  const bool shiftIntt = m_inttOffset[0] != 0 || m_inttOffset[1] != 0;
  const bool moveInttSurfaces = shiftIntt && alignmentTransformationContainer::use_alignment;
  const Acts::Vector3 inttShiftMm(m_inttOffset[0] * cm, m_inttOffset[1] * cm, 0);

  auto reject = [&](const std::string& why)
  {
    ++d.rejected[why];
  };

  // ------------------------------------------------------------------------------------------
  // fit of a measurement subset along the guess; the surface sequence is built here
  auto runFit = [&](FitJob& job, const std::vector<const Measurement*>& meas, const HelixGuess& h,
                    const ActsTrackFittingAlgorithm::TrackParameters& seed, const Acts::Surface* perigee,
                    bool withMaterial = true, const DirectFitter* fitter = nullptr) -> bool
  {
    if (meas.empty() || meas.size() > 250)
    {
      job.failReason = "number of measurements 0 or > 250";
      return false;
    }
    std::vector<std::pair<double, const Acts::Surface*>> ordered;
    std::set<const Acts::Surface*> used;
    double sMax = 0;
    for (const auto* m : meas)
    {
      if (!used.insert(m->surface.get()).second)
      {
        continue;  // one measurement per surface
      }
      const auto index = static_cast<ActsSourceLink::Index>(job.measurements.size());
      const std::array<Acts::BoundIndices, 2> indices{Acts::eBoundLoc0, Acts::eBoundLoc1};
      job.measurements.emplaceMeasurement<2>(m->surface->geometryId(), indices, m->local, m->cov);
      ActsSourceLink sl(m->surface->geometryId(), index, m->key);
      job.sourceLinks.emplace_back(sl);
      job.slSurfaces.push_back(m->surface.get());
      job.slMeasurements.push_back(m);
      ordered.emplace_back(m->s, m->surface.get());
      sMax = std::max(sMax, m->s);
    }
    for (const auto* surface : d.materialSurfaces)
    {
      if (!withMaterial)
      {
        break;
      }
      double r = 0, ox = 0, oy = 0;
      const Acts::Vector3 center = surface->center(d.gctx) / cm;
      if (surface->type() == Acts::Surface::Cylinder)
      {
        const auto& bounds = static_cast<const Acts::CylinderBounds&>(surface->bounds());
        r = bounds.get(Acts::CylinderBounds::eR) / cm;
        ox = center.x();
        oy = center.y();
      }
      else
      {
        r = std::hypot(center.x(), center.y());
      }
      const double s = h.arcLengthToRadius(r, ox, oy);
      if (!std::isfinite(s) || s > sMax)
      {
        continue;
      }
      if (surface->type() == Acts::Surface::Cylinder)
      {
        const auto& bounds = static_cast<const Acts::CylinderBounds&>(surface->bounds());
        if (std::abs(h.zAt(s) - center.z()) > bounds.get(Acts::CylinderBounds::eHalfLengthZ) / cm)
        {
          continue;
        }
      }
      ordered.emplace_back(s + kMaterialTieMargin, surface);
      ++job.nMaterial;
    }
    std::stable_sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b)
                     { return a.first < b.first; });
    job.ordered = ordered;
    for (const auto& [s, surface] : ordered)
    {
      job.sequence.push_back(surface);
    }

    CalibratorAdapter calibratorAdapter{d.calibrator, job.measurements};
    IndexSurfaceAccessor accessor{&job.slSurfaces};

    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
    extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&d.updater);
    extensions.smoother.connect<&Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(&d.smoother);
    extensions.calibrator.connect<&CalibratorAdapter::calibrate>(&calibratorAdapter);
    extensions.surfaceAccessor.connect<&IndexSurfaceAccessor::operator()>(&accessor);
    ChiOutlierFinder outlierFinder{m_tpcOutlierChi2, m_siOutlierChi2};
    if (m_tpcOutlierChi2 > 0 || m_siOutlierChi2 > 0)
    {
      extensions.outlierFinder.connect<&ChiOutlierFinder::operator()>(&outlierFinder);
    }

    Acts::PropagatorPlainOptions plainOptions(d.gctx, magContext);
    plainOptions.maxSteps = m_maxSteps;
    Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(
        d.gctx, magContext, calibContext, extensions, plainOptions, perigee, m_multipleScattering, m_energyLoss);

    auto result = (fitter ? fitter : d.fitter.get())->fit(job.sourceLinks.begin(), job.sourceLinks.end(), seed, kfOptions, job.sequence, job.tracks);
    if (!result.ok())
    {
      job.failReason = result.error().message();
      if (Verbosity() > 1)
      {
        std::cout << Name() << ": fit failed - " << job.failReason << std::endl;
      }
      return false;
    }
    if (!result.value().hasReferenceSurface())
    {
      job.failReason = "no parameters at the perigee (extrapolation to the beam failed)";
      return false;
    }
    job.track.emplace(result.value());
    return true;
  };

  auto summarize = [&](const FitJob& job, PerigeeSummary& p)
  {
    p = PerigeeSummary();
    if (!job.track)
    {
      return;
    }
    const auto& t = *job.track;
    const auto& par = t.parameters();
    const auto& cov = t.covariance();
    p.ok = 1;
    p.d0 = par[Acts::eBoundLoc0] / cm;
    p.z0 = par[Acts::eBoundLoc1] / cm;
    p.phi = par[Acts::eBoundPhi];
    p.theta = par[Acts::eBoundTheta];
    p.pt = std::abs(1. / par[Acts::eBoundQOverP]) * std::sin(par[Acts::eBoundTheta]);
    p.d0err = std::sqrt(cov(Acts::eBoundLoc0, Acts::eBoundLoc0)) / cm;
    p.z0err = std::sqrt(cov(Acts::eBoundLoc1, Acts::eBoundLoc1)) / cm;
    const auto ts = Acts::MultiTrajectoryHelpers::trajectoryState(job.tracks.trackStateContainer(), t.tipIndex());
    p.chi2 = ts.chi2Sum;
    p.ndf = static_cast<int>(ts.NDF);
  };

  // ------------------------------------------------------------------------------------------
  for (unsigned int itrack = 0; itrack < d.siTpcTracks->size(); ++itrack)
  {
    const SiTpc_Track* track = d.siTpcTracks->get_track(itrack);
    if (!track)
    {
      continue;
    }
    ++d.nInput;

    // ---- selection
    if (track->get_fit_status() != 1)  // SiTpcHelixFit::Ok (StraightLine has no curvature)
    {
      reject("own fit not ok");
      continue;
    }
    const int ownCharge = track->get_charge();
    const double pt = track->get_pt();
    if (ownCharge == 0 || !std::isfinite(pt) || pt < m_minPt || !(track->get_radius() > 0))
    {
      reject("charge / pt");
      continue;
    }
    if (track->get_si_trajectory_index() >= d.siTrajectories->size())
    {
      reject("no Si trajectory");
      continue;
    }
    const Si_Trajectory* traj = d.siTrajectories->get_trajectory(track->get_si_trajectory_index());
    if (!traj)
    {
      reject("no Si trajectory");
      continue;
    }

    short crossing = 0;
    if (d.crossings)
    {
      const TpcCrossingDecision* decision = d.crossings->get_decision(track->get_tpc_assembled_track_id());
      crossing = decision ? decision->get_selected_crossing() : SHRT_MAX;
    }

    // ---- guess in the fit frame: own frame = fit frame - Si beam spot
    HelixGuess h;
    h.cx = track->get_circle_x() + m_siBeam[0];
    h.cy = track->get_circle_y() + m_siBeam[1];
    h.R = track->get_radius();
    {
      // start = point of the circle closest to the beam axis
      const double vx = m_siBeam[0] - h.cx, vy = m_siBeam[1] - h.cy;
      const double v = std::hypot(vx, vy);
      if (!(v > 0))
      {
        reject("degenerate circle");
        continue;
      }
      h.x0 = h.cx + h.R * vx / v;
      h.y0 = h.cy + h.R * vy / v;
      // sense of rotation from the direction at the own PCA (translation invariant)
      const double tx = std::cos(track->get_phi()), ty = std::sin(track->get_phi());
      h.omega = ((h.x0 - h.cx) * ty - (h.y0 - h.cy) * tx) >= 0 ? 1. : -1.;
      h.tanl = track->get_tanl();
      h.z0 = track->get_z0() - m_siActsZShift;
    }
    if (!std::isfinite(h.cx) || !std::isfinite(h.cy) || !std::isfinite(h.tanl) || !std::isfinite(h.z0))
    {
      reject("invalid guess");
      continue;
    }
    // charge from the sense of rotation and the field (positive charge, B_z > 0: clockwise)
    const int charge = -h.omega * d.bzSign > 0 ? 1 : -1;
    if (charge != ownCharge)
    {
      ++d.nChargeFlip;
    }

    // ---- Si measurements: chain points and hits per layer
    std::map<int, Acts::Vector3> chainPoint;  // own frame [cm]
    for (unsigned int k = 0; k < traj->size_points(); ++k)
    {
      const int layer = traj->get_point_layer(k);
      if (layer >= 0 && layer < 7 && !chainPoint.contains(layer))
      {
        chainPoint[layer] = Acts::Vector3(traj->get_point_x(k), traj->get_point_y(k), traj->get_point_z(k));
      }
    }
    std::map<int, std::vector<Si_Trajectory::HitIndex>> hitsByLayer;
    for (const auto& hit : traj->get_hit_indices())
    {
      const int layer = TrkrDefs::getLayer(hit.first);
      if (layer < 7)
      {
        hitsByLayer[layer].push_back(hit);
      }
    }

    std::vector<Measurement> meas;
    meas.reserve(64);
    unsigned int nProjected = 0;
    for (const auto& [layer, hits] : hitsByLayer)
    {
      const bool isIntt = layer >= 3;
      const auto cp = chainPoint.find(layer);
      // chain point in the fit frame (on the moved INTT surfaces: the own INTT points contain the offset)
      std::optional<Acts::Vector3> target;
      if (cp != chainPoint.end())
      {
        target = Acts::Vector3(cp->second.x() + m_siBeam[0], cp->second.y() + m_siBeam[1], cp->second.z() - m_siActsZShift);
      }

      // majority sensor of the chain hits in this layer
      std::map<TrkrDefs::hitsetkey, int> nPerSensor;
      for (const auto& hit : hits)
      {
        ++nPerSensor[hit.first];
      }
      const auto sensor = std::max_element(nPerSensor.begin(), nPerSensor.end(), [](const auto& a, const auto& b)
                                           { return a.second < b.second; })
                              ->first;

      Measurement m;
      m.det = isIntt ? kIntt : kMvtx;
      m.layer = layer;
      bool done = false;

      if (m_siMeasurement == SiMeasurement::OfficialClusters)
      {
        // 1) official cluster containing most chain hits
        TrkrDefs::cluskey best = TrkrDefs::CLUSKEYMAX;
        std::set<TrkrDefs::cluskey> otherStrobes;  // clusters found only in the neighbouring MVTX strobe
        {
          std::map<TrkrDefs::cluskey, int> nPerCluster;
          for (const auto& hit : hits)
          {
            bool otherStrobe = false;
            const auto ckey = d.officialClusterOfHit(hit.first, hit.second, otherStrobe);
            if (ckey != TrkrDefs::CLUSKEYMAX)
            {
              ++nPerCluster[ckey];
              otherStrobes.insert(ckey);
              if (!otherStrobe)
              {
                otherStrobes.erase(ckey);
              }
            }
          }
          int nbest = 0;
          for (const auto& [ckey, n] : nPerCluster)
          {
            if (n > nbest)
            {
              nbest = n;
              best = ckey;
            }
          }
        }
        bool byPosition = false;
        // 2) closest official cluster of the majority sensor
        if (best == TrkrDefs::CLUSKEYMAX && target)
        {
          const auto surface = d.geometry->maps().getSiliconSurface(sensor);
          if (surface)
          {
            if (moveInttSurfaces && isIntt)
            {
              d.moveInttSurface(*surface, m_inttOffset);
            }
            const Acts::Vector3 lt = surface->localToGlobalTransform(d.gctx).inverse() * (*target * cm);
            const auto& win = isIntt ? kInttMatchWindow : kMvtxMatchWindow;
            double bestDist = std::numeric_limits<double>::max();
            const auto range = d.siClusters->getClusters(sensor);
            for (auto it = range.first; it != range.second; ++it)
            {
              const Acts::Vector2 lc = d.geometry->getLocalCoords(it->first, it->second) * cm;
              const double d0 = std::abs(lc.x() - lt.x()) / cm, d1 = std::abs(lc.y() - lt.y()) / cm;
              if (d0 > win[0] || d1 > win[1])
              {
                continue;
              }
              const double dist = square(d0 / win[0]) + square(d1 / win[1]);
              if (dist < bestDist)
              {
                bestDist = dist;
                best = it->first;
              }
            }
            byPosition = best != TrkrDefs::CLUSKEYMAX;
          }
        }
        if (best != TrkrDefs::CLUSKEYMAX)
        {
          TrkrCluster* cluster = d.siClusters->findCluster(best);
          const auto surface = cluster ? d.geometry->maps().getSurface(best, cluster) : nullptr;
          if (surface)
          {
            if (moveInttSurfaces && isIntt)
            {
              d.moveInttSurface(*surface, m_inttOffset);
            }
            m.key = best;
            m.surface = surface;
            m.local = d.geometry->getLocalCoords(best, cluster) * cm;
            if (shiftIntt && isIntt && !moveInttSurfaces)
            {
              // geometry without alignment transforms: shift in the sensor plane
              const Acts::Vector3 dl = surface->localToGlobalTransform(d.gctx).linear().transpose() * inttShiftMm;
              m.local += dl.head<2>();
            }
            const Acts::Vector3 g = surface->localToGlobal(d.gctx, m.local, Acts::Vector3::UnitZ()) / cm;
            const auto err = ClusterErrorPara::get_clusterv5_modified_error(cluster, std::hypot(g.x(), g.y()), best);
            m.cov(0, 0) = err.first * cm * cm;
            m.cov(1, 1) = err.second * cm * cm;
            m.global = g;
            m.official = true;
            done = true;
            ++(byPosition ? d.nSiOfficialPosition : d.nSiOfficialAssoc);
            d.nSiOtherStrobe += otherStrobes.contains(best);
          }
        }
      }

      // 3) chain point projected on the official surface of its sensor
      if (!done && target)
      {
        const auto surface = d.geometry->maps().getSiliconSurface(sensor);
        if (surface)
        {
          if (moveInttSurfaces && isIntt)
          {
            d.moveInttSurface(*surface, m_inttOffset);
          }
          const Acts::Vector3 lt = surface->localToGlobalTransform(d.gctx).inverse() * (*target * cm);
          m.key = TrkrDefs::genClusKey(sensor, std::numeric_limits<uint32_t>::max() - layer);
          m.surface = surface;
          m.local = lt.head<2>();
          const auto& err = isIntt ? m_inttChainErr : m_mvtxChainErr;
          m.cov(0, 0) = square(err[0] * cm);
          m.cov(1, 1) = square(err[1] * cm);
          m.global = surface->localToGlobal(d.gctx, m.local, Acts::Vector3::UnitZ()) / cm;
          m.official = false;
          done = true;
          ++nProjected;
          ++d.nSiProjected;
        }
      }
      if (!done)
      {
        ++d.nSiMissing;
        continue;
      }
      m.s = h.arcLength(m.global.x(), m.global.y());
      meas.push_back(std::move(m));
    }

    const auto nmvtx = std::count_if(meas.begin(), meas.end(), [](const auto& m)
                                     { return m.det == kMvtx; });
    const auto nsi = meas.size();
    if (nmvtx < static_cast<long>(m_minMvtx) || nsi < m_minSi)
    {
      reject("too few Si measurements");
      continue;
    }

    // ---- TPC: one artificial plane per cluster of the own Si+TPC fit
    const double zOffset = m_applyTrackTpcZOffset ? track->get_tpc_z_offset() : 0.;
    // TPC axis in the fit frame (radial direction of the pad rows)
    const double tpcAxisX = m_tpcFrame == TpcFrame::BeamAxis ? m_siBeam[0] - m_tpcLine[0] + m_tpcExtraShift[0] : m_tpcExtraShift[0];
    const double tpcAxisY = m_tpcFrame == TpcFrame::BeamAxis ? m_siBeam[1] - m_tpcLine[1] + m_tpcExtraShift[1] : m_tpcExtraShift[1];
    unsigned int planeIndex = 0;
    for (unsigned int ip = 0; ip < track->size_points(); ++ip)
    {
      if (track->get_point_type(ip) != SiTpc_Track::TpcPoint)
      {
        continue;
      }
      const unsigned int icluster = track->get_point_source(ip);
      const Tpc_PolyCluster* cluster = icluster < d.tpcClusters->size() ? d.tpcClusters->get_cluster(icluster) : nullptr;
      if (!cluster || !cluster->isValid() || cluster->size_hits() == 0)
      {
        continue;
      }
      const TrkrDefs::hitsetkey hitsetkey = cluster->get_hit_index(0).first;
      const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
      const Acts::Vector3 raw(cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z());
      Acts::Vector3 pos;
      if (m_tpcFrame == TpcFrame::BeamAxis)
      {
        pos = Acts::Vector3(raw.x() - m_tpcLine[0] + m_siBeam[0], raw.y() - m_tpcLine[1] + m_siBeam[1], raw.z() - m_siActsZShift);
      }
      else
      {
        pos = d.geometry->transformTpcEnvelopeToWorld(raw);
      }
      pos += Acts::Vector3(m_tpcExtraShift[0], m_tpcExtraShift[1], m_tpcExtraShift[2] - zOffset);
      if (!pos.allFinite())
      {
        continue;
      }
      ++d.nTpc;

      // plane: normal radial, local x = r*phi direction, local y = z
      const double phi = std::atan2(pos.y() - tpcAxisY, pos.x() - tpcAxisX);
      const Acts::Vector3 nrm(std::cos(phi), std::sin(phi), 0);
      const Acts::Vector3 dir = h.direction(pos.x(), pos.y());
      const Acts::Vector3 dirXY = Acts::Vector3(dir.x(), dir.y(), 0).normalized();
      const double angle = std::acos(std::min(1., std::abs(dirXY.dot(nrm))));
      if (angle > m_maxTpcCrossingAngle)
      {
        ++d.nTpcAngleCut;
        continue;
      }
      Acts::Transform3 transform = Acts::Transform3::Identity();
      Acts::RotationMatrix3 rotation;
      rotation.col(0) = Acts::Vector3(-std::sin(phi), std::cos(phi), 0);
      rotation.col(1) = Acts::Vector3(0, 0, 1);
      rotation.col(2) = nrm;
      transform.linear() = rotation;
      transform.translation() = pos * cm;
      auto plane = Acts::Surface::makeShared<Acts::PlaneSurface>(transform, std::make_shared<const Acts::RectangleBounds>(20 * cm, 120 * cm));
      plane->assignGeometryId(Acts::GeometryIdentifier().withVolume(kArtificialVolume).withLayer(layer).withSensitive(++planeIndex));

      if (m_copyTpcMaterial)
      {
        // material of the official TPC surface at this position (homogeneous slab on the plane)
        TrkrDefs::subsurfkey subsurfkey = 0;
        const auto official = d.geometry->get_tpc_surface_from_coords(hitsetkey, pos, subsurfkey);
        if (official && official->surfaceMaterial())
        {
          const Acts::MaterialSlab slab = official->surfaceMaterial()->materialSlab(Acts::Vector3(pos * cm));
          if (slab.thickness() > 0)
          {
            plane->assignSurfaceMaterial(std::make_shared<const Acts::HomogeneousSurfaceMaterial>(slab));
            ++d.nTpcMaterial;
          }
        }
      }

      TrkrDefs::cluskey key = cluster->get_trkr_cluster_key();
      if (key == TrkrDefs::CLUSKEYMAX || TrkrDefs::getTrkrId(key) != TrkrDefs::tpcId)
      {
        key = TrkrDefs::genClusKey(hitsetkey, cluster->get_cluster_id());
      }

      Measurement m;
      m.det = kTpc;
      m.layer = layer;
      m.key = key;
      m.surface = plane;
      m.local = Acts::ActsVector<2>::Zero();
      // errors: the standard sPHENIX parametrisation (ClusterErrorPara, as PHActsTrkFitter) on top of the
      // TpcClusterizer-type base errors stored by Tpc_PolyClusterizer (rms_x: r*phi, rms_z: z) and the cluster
      // sizes (pads, time bins); or constants per module
      bool stdErrors = false;
      if (m_tpcErrorMode == TpcErrors::Standard)
      {
        TrkrClusterv5 tmp;
        tmp.setPhiError(static_cast<float>(std::hypot(cluster->get_rms_x(), cluster->get_rms_y())));
        tmp.setZError(static_cast<float>(std::abs(cluster->get_rms_z())));
        tmp.setPhiSize(static_cast<char>(std::min(cluster->get_phi_width(), 127U)));
        tmp.setZSize(static_cast<char>(std::min(cluster->get_time_width(), 127U)));
        tmp.setEdge(0);
        tmp.setOverlap(0);
        const auto e = ClusterErrorPara::get_clusterv5_modified_error(&tmp, std::hypot(pos.x(), pos.y()), key);
        if (std::isfinite(e.first) && std::isfinite(e.second) && e.first > 0 && e.second > 0)
        {
          m.cov(0, 0) = e.first * cm * cm * square(m_tpcErrorScale[0]);
          m.cov(1, 1) = e.second * cm * cm * square(m_tpcErrorScale[1]);
          stdErrors = true;
          ++d.nTpcErrStd;
        }
        else
        {
          ++d.nTpcErrFallback;
        }
      }
      if (!stdErrors)
      {
        const auto& err = m_tpcErr[tpcModule(layer)];
        m.cov(0, 0) = square(err[0] * cm);
        m.cov(1, 1) = square(err[1] * cm);
      }
      m.global = pos;
      m.s = h.arcLength(pos.x(), pos.y());
      m.official = false;
      m.crossingAngle = angle;
      meas.push_back(std::move(m));
    }

    const auto ntpc = meas.size() - nsi;
    if (ntpc < m_minTpc)
    {
      reject("too few TPC measurements");
      continue;
    }
    ++d.nSelected;

    std::stable_sort(meas.begin(), meas.end(), [](const auto& a, const auto& b)
                     { return a.s < b.s; });

    // ---- seed at the perigee = line along z through the beam axis
    const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(m_siBeam[0] * cm, m_siBeam[1] * cm, 0));
    const Acts::Vector3 seedDir = h.direction(h.x0, h.y0);
    const double p = pt * std::sqrt(1 + h.tanl * h.tanl);
    auto makeSeed = [&](double sigPos, double sigAngle, double relQop) -> std::optional<ActsTrackFittingAlgorithm::TrackParameters>
    {
      Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = square(sigPos * cm);
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = square(sigPos * cm);
      cov(Acts::eBoundPhi, Acts::eBoundPhi) = square(sigAngle);
      cov(Acts::eBoundTheta, Acts::eBoundTheta) = square(sigAngle);
      cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = square(relQop / p);
      cov(Acts::eBoundTime, Acts::eBoundTime) = square(1. * Acts::UnitConstants::ns);
      const Acts::Vector4 pos4(h.x0 * cm, h.y0 * cm, h.z0 * cm, 10 * Acts::UnitConstants::ns);
      auto r = ActsTrackFittingAlgorithm::TrackParameters::create(d.gctx, perigee, pos4, seedDir, charge / p, cov,
                                                                  Acts::ParticleHypothesis::pion());
      if (!r.ok())
      {
        return std::nullopt;
      }
      return r.value();
    };
    const auto seed = makeSeed(m_seedSigmaPos, m_seedSigmaAngle, m_seedSigmaRelQop);
    if (!seed)
    {
      reject("seed not on perigee");
      continue;
    }

    // ---- full fit
    std::vector<const Measurement*> all;
    for (const auto& m : meas)
    {
      all.push_back(&m);
    }
    FitJob jobMaterial;
    bool ok = runFit(jobMaterial, all, h, *seed, perigee.get());
    bool retriedNoMaterial = false;
    FitJob jobPlain;
    if (!ok && m_retryWithoutMaterial && m_useMaterialSurfaces)
    {
      // diagnosis / recovery: the same fit with the measurement surfaces only
      ok = runFit(jobPlain, all, h, *seed, perigee.get(), false);
      retriedNoMaterial = true;
      if (ok)
      {
        ++d.nRecoveredNoMaterial;
        ++d.failReasons["(recovered without material) " + jobMaterial.failReason][pt < 0.5 ? 0 : (pt < 1.0 ? 1 : (pt < 2.0 ? 2 : 3))];
      }
    }
    FitJob& job = retriedNoMaterial ? jobPlain : jobMaterial;
    if (!ok && d.verboseFitter && d.nDebugged < m_debugFailedFits)
    {
      ++d.nDebugged;
      std::cout << "==== " << Name() << " debug of failed fit " << d.nDebugged << ": track " << itrack << " crossing " << crossing
                << " charge " << charge << " own pt " << pt << " tanl " << h.tanl << " reason: " << jobMaterial.failReason << std::endl
                << "  guess circle (" << h.cx << ", " << h.cy << ") R " << h.R << " omega " << h.omega << ", start (" << h.x0 << ", "
                << h.y0 << ", " << h.z0 << ") cm, seed direction " << seedDir.transpose() << std::endl;
      for (const auto& m : meas)
      {
        std::cout << "  meas det " << m.det << " layer " << m.layer << " s " << m.s << " cm  global (" << m.global.x() << ", "
                  << m.global.y() << ", " << m.global.z() << ") r " << std::hypot(m.global.x(), m.global.y()) << "  sigma ("
                  << std::sqrt(m.cov(0, 0)) / cm << ", " << std::sqrt(m.cov(1, 1)) / cm << ") cm" << (m.official ? "" : " projected")
                  << std::endl;
      }
      for (const auto& [sOrd, surface] : jobMaterial.ordered)
      {
        const bool isMeas = std::any_of(meas.begin(), meas.end(), [surface](const auto& m)
                                        { return m.surface.get() == surface; });
        const Acts::Vector3 c = surface->center(d.gctx) / cm;
        std::cout << "  seq s " << sOrd << " " << (isMeas ? "meas" : "mat ") << " " << surface->geometryId() << " type "
                  << static_cast<int>(surface->type()) << " centre (" << c.x() << ", " << c.y() << ", " << c.z() << ")" << std::endl;
      }
      std::cout << "  ---- Acts verbose refit" << std::endl;
      FitJob debugJob;
      runFit(debugJob, all, h, *seed, perigee.get(), true, d.verboseFitter.get());
      std::cout << "  ---- end of debug (" << (debugJob.track ? "ok" : debugJob.failReason) << ")" << std::endl;
    }
    if (!ok)
    {
      ++d.nFitFail;
      const int ptBin = pt < 0.5 ? 0 : (pt < 1.0 ? 1 : (pt < 2.0 ? 2 : 3));
      ++d.failReasons[job.failReason][ptBin];
      if (d.trackTree)
      {
        // failed fits are kept in the QA tree (ok = 0) to study them against pT, eta, hits
        auto& r = d.trackRow;
        r.event = static_cast<int>(d.nEvents - 1);
        r.track = -1;
        r.crossing = crossing;
        r.charge = charge;
        r.nmvtx = static_cast<int>(nmvtx);
        r.nintt = static_cast<int>(nsi - nmvtx);
        r.ntpc = static_cast<int>(ntpc);
        r.nsiProjected = static_cast<int>(nProjected);
        r.nmaterial = job.nMaterial;
        r.noutliers = 0;
        r.retry_no_material = retriedNoMaterial ? 1 : 0;
        r.own_pt = pt;
        r.own_phi = track->get_phi();
        r.own_tanl = track->get_tanl();
        r.own_z0 = track->get_z0();
        r.own_tpcZOffset = track->get_tpc_z_offset();
        r.beamX = m_siBeam[0];
        r.beamY = m_siBeam[1];
        r.full = PerigeeSummary();
        r.tpc = PerigeeSummary();
        r.si = PerigeeSummary();
        d.trackTree->Fill();
      }
      continue;
    }
    ++d.nFitOk;

    // ---- seeds + SvtxTrack
    auto makeTrackSeed = [&](int det) -> TrackSeed*
    {
      TrackSeed_v2 s;
      s.set_X0(h.cx);
      s.set_Y0(h.cy);
      s.set_qOverR(charge / h.R);
      s.set_Z0(h.z0);
      s.set_slope(h.tanl);
      s.set_phi(std::atan2(seedDir.y(), seedDir.x()));
      s.set_crossing(crossing);
      for (const auto& m : meas)
      {
        if ((det == kTpc) == (m.det == kTpc) && (m.det == kTpc || m.official))
        {
          s.insert_cluster_key(m.key);
        }
      }
      return (det == kTpc ? d.tpcSeeds : d.siSeeds)->insert(&s);
    };

    const unsigned int id = d.trackMap->size();
    SvtxTrack_v4 svtx;
    svtx.set_id(id);
    svtx.set_crossing(crossing);
    svtx.set_silicon_seed(makeTrackSeed(kMvtx));
    svtx.set_tpc_seed(makeTrackSeed(kTpc));
    {
      const auto& t = *job.track;
      const Acts::BoundTrackParameters params(t.referenceSurface().getSharedPtr(), t.parameters(), t.covariance(), t.particleHypothesis());
      SvtxTrackState_v3 out(0.0);
      out.set_localX(0.0);
      out.set_localY(0.0);
      out.set_x(0.0);
      out.set_y(0.0);
      out.set_z(0.0);
      svtx.insert_state(&out);
      const Acts::Vector3 position = params.position(d.gctx) / cm;
      svtx.set_x(position.x());
      svtx.set_y(position.y());
      svtx.set_z(position.z());
      svtx.set_px(params.momentum().x());
      svtx.set_py(params.momentum().y());
      svtx.set_pz(params.momentum().z());
      svtx.set_charge(params.charge() > 0 ? 1 : -1);
      const auto ts = Acts::MultiTrajectoryHelpers::trajectoryState(job.tracks.trackStateContainer(), t.tipIndex());
      svtx.set_chisq(ts.chi2Sum);
      svtx.set_ndf(ts.NDF);
      ActsTransformations transformer;
      transformer.setVerbosity(Verbosity());
      if (params.covariance())
      {
        const auto rotated = transformer.rotateActsCovToSvtxTrack(params);
        for (int i = 0; i < 6; ++i)
        {
          for (int j = 0; j < 6; ++j)
          {
            svtx.set_error(i, j, rotated(i, j));
          }
        }
      }
      if (m_fillStates)
      {
        transformer.fillSvtxTrackStates(job.tracks.trackStateContainer(), t.tipIndex(), &svtx, d.gctx);
      }
    }
    d.trackMap->insertWithKey(&svtx, id);

    int nOutliers = 0;
    job.tracks.trackStateContainer().visitBackwards(job.track->tipIndex(), [&](const auto& state)
                                                    {
      if (state.typeFlags().isOutlier())
      {
        ++nOutliers;
        const auto key = state.getUncalibratedSourceLink().template get<ActsSourceLink>().cluskey();
        ++(TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId ? d.nOutlierTpc : d.nOutlierSi);
      }
      return true; });

    if (Verbosity() > 1)
    {
      std::cout << Name() << ": track " << itrack << " crossing " << crossing << " Si " << nsi << " (projected " << nProjected
                << ") TPC " << ntpc << " material " << job.nMaterial << " own pt " << pt << " -> fit pt " << svtx.get_pt()
                << " chi2/ndf " << svtx.get_chisq() << "/" << svtx.get_ndf() << std::endl;
    }

    // ---- QA
    if (d.trackTree)
    {
      auto& r = d.trackRow;
      r.event = static_cast<int>(d.nEvents - 1);
      r.track = static_cast<int>(id);
      r.crossing = crossing;
      r.charge = charge;
      r.nmvtx = static_cast<int>(nmvtx);
      r.nintt = static_cast<int>(nsi - nmvtx);
      r.ntpc = static_cast<int>(ntpc);
      r.nsiProjected = static_cast<int>(nProjected);
      r.nmaterial = job.nMaterial;
      r.noutliers = nOutliers;
      r.retry_no_material = retriedNoMaterial ? 1 : 0;
      r.own_pt = pt;
      r.own_phi = track->get_phi();
      r.own_tanl = track->get_tanl();
      r.own_z0 = track->get_z0();
      r.own_tpcZOffset = track->get_tpc_z_offset();
      r.beamX = m_siBeam[0];
      r.beamY = m_siBeam[1];
      summarize(job, r.full);
      r.tpc = PerigeeSummary();
      r.si = PerigeeSummary();
      if (m_qaSubFits)
      {
        // loose seeds: the sub-fits must not be pulled to the full-fit perigee
        const auto looseSeed = makeSeed(1.0, 0.1, 0.5);
        if (looseSeed)
        {
          std::vector<const Measurement*> tpcOnly, siOnly;
          for (const auto* m : all)
          {
            (m->det == kTpc ? tpcOnly : siOnly).push_back(m);
          }
          FitJob tpcJob;
          if (runFit(tpcJob, tpcOnly, h, *looseSeed, perigee.get()))
          {
            summarize(tpcJob, r.tpc);
          }
          if (siOnly.size() >= 3)
          {
            FitJob siJob;
            if (runFit(siJob, siOnly, h, *looseSeed, perigee.get()))
            {
              summarize(siJob, r.si);
            }
          }
        }
      }
      d.trackTree->Fill();

      const auto& mj = job.tracks.trackStateContainer();
      mj.visitBackwards(job.track->tipIndex(), [&](const auto& state)
                        {
        if (!state.typeFlags().isMeasurement() || !state.hasSmoothed())
        {
          return true;
        }
        const auto sl = state.getUncalibratedSourceLink().template get<ActsSourceLink>();
        if (sl.index() >= job.slMeasurements.size())
        {
          return true;
        }
        const Measurement* m = job.slMeasurements[sl.index()];
        const auto& sm = state.smoothed();
        const auto& sc = state.smoothedCovariance();
        const Acts::BoundTrackParameters sp(state.referenceSurface().getSharedPtr(), sm, sc, Acts::ParticleHypothesis::pion());
        const Acts::Vector3 f = sp.position(d.gctx) / cm;
        auto& s = d.stateRow;
        s.event = static_cast<int>(d.nEvents - 1);
        s.track = static_cast<int>(id);
        s.det = m->det;
        s.layer = static_cast<int>(m->layer);
        s.official = m->official ? 1 : 0;
        s.gx = m->global.x();
        s.gy = m->global.y();
        s.gz = m->global.z();
        s.fx = f.x();
        s.fy = f.y();
        s.fz = f.z();
        s.res0 = (m->local[0] - sm[Acts::eBoundLoc0]) / cm;
        s.res1 = (m->local[1] - sm[Acts::eBoundLoc1]) / cm;
        s.sigMeas0 = std::sqrt(m->cov(0, 0)) / cm;
        s.sigMeas1 = std::sqrt(m->cov(1, 1)) / cm;
        s.sigFit0 = std::sqrt(sc(Acts::eBoundLoc0, Acts::eBoundLoc0)) / cm;
        s.sigFit1 = std::sqrt(sc(Acts::eBoundLoc1, Acts::eBoundLoc1)) / cm;
        s.pt = svtx.get_pt();
        s.s = m->s;
        s.crossingAngle = m->crossingAngle;
        d.stateTree->Fill();
        return true; });
    }
  }

  d.restoreIntt();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________________________
int SiTpcActsFitter::End(PHCompositeNode* /*topNode*/)
{
  auto& d = *m_impl;
  std::cout << Name() << "::End - events " << d.nEvents << ", input tracks " << d.nInput << ", selected " << d.nSelected
            << ", fit ok " << d.nFitOk << ", fit failed " << d.nFitFail << std::endl
            << "  Si measurements: official cluster by hit association " << d.nSiOfficialAssoc
            << " (of them the copy in the neighbouring MVTX strobe, pruned duplicate " << d.nSiOtherStrobe << "), by position "
            << d.nSiOfficialPosition << ", projected chain points " << d.nSiProjected << ", layers without measurement "
            << d.nSiMissing << std::endl
            << "  TPC clusters " << d.nTpc << ", dropped by the crossing angle " << d.nTpcAngleCut << ", planes with copied material "
            << d.nTpcMaterial << std::endl
            << "  tracks with a fit charge different from the own charge (sense of rotation vs field) " << d.nChargeFlip << std::endl
            << "  fits recovered by the retry without the material surfaces " << d.nRecoveredNoMaterial << std::endl
            << "  TPC errors: standard parametrisation " << d.nTpcErrStd << ", constant (fallback or mode) " << d.nTpc - d.nTpcErrStd << std::endl
            << "  outlier measurements (chi2 cut TPC " << m_tpcOutlierChi2 << ", Si " << m_siOutlierChi2 << "): TPC " << d.nOutlierTpc
            << ", Si " << d.nOutlierSi << std::endl;
  for (const auto& [why, n] : d.rejected)
  {
    std::cout << "  rejected (" << why << "): " << n << std::endl;
  }
  if (!d.failReasons.empty())
  {
    std::cout << "  failed fits by reason (pT < 0.5 / 0.5-1 / 1-2 / > 2 GeV), max propagation steps " << m_maxSteps << ":" << std::endl;
    for (const auto& [why, n] : d.failReasons)
    {
      std::cout << "    " << why << ": " << n[0] << " / " << n[1] << " / " << n[2] << " / " << n[3] << std::endl;
    }
  }
  if (d.qaFile)
  {
    d.qaFile->cd();
    d.trackTree->Write();
    d.stateTree->Write();
    d.qaFile->Close();
    delete d.qaFile;
    d.qaFile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
