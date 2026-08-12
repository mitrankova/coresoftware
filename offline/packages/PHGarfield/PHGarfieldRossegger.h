// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHGARFIELD_PHGARFIELDROSSEGGER_H
#define PHGARFIELD_PHGARFIELDROSSEGGER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <utility>
#include <vector>

class PHCompositeNode;

class PHGarfieldRossegger : public SubsysReco
{
 public:
  explicit PHGarfieldRossegger(const std::string& name = "PHGarfieldRossegger");
  ~PHGarfieldRossegger() override = default;

  int Init(PHCompositeNode*) override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setGeometryCm(double inner_radius, double outer_radius, double half_length);
  void setSourceRadiusCm(double min_radius, double max_radius);
  void setDensity(double reference_density_nC_per_m3, double k_eff, double radial_power_alpha);
  void setPhiModulation(double m1_amplitude, double m1_phase, double m12_amplitude, double m12_phase);
  void setSourceGrid(unsigned int nr, unsigned int nphi, unsigned int nz);
  void setObservationGrid(unsigned int nr, unsigned int nphi, unsigned int nz);
  void setModeTruncation(unsigned int m_phi_max, unsigned int n_radial_modes, unsigned int n_longitudinal_modes);
  void setAutoAxisymmetric(bool value) { m_autoAxisymmetric = value; }

  // Split observation radial bins across condor jobs. Histograms keep full
  // binning, so PART files can be merged with hadd.
  void setRadialJob(unsigned int job_index, unsigned int n_jobs);

  void setGarfieldOutputFile(const std::string& filename) { m_garfieldOutputFile = filename; }
  void setField3DOutputFile(const std::string& filename) { m_field3DOutputFile = filename; }
  void setWriteField3D(bool value) { m_writeField3D = value; }
  void setVerifyOutput(bool value) { m_verifyOutput = value; }

 private:
  int calculate();
  bool validateConfig() const;
  unsigned int effectivePhiMax() const;
  std::pair<unsigned int, unsigned int> radialRange() const;

  double radialBoundaryFunction(double kval, unsigned int mode_m) const;
  double radialBasis(double kval, unsigned int mode_m, double radius_m) const;
  double radialBasisDerivative(double kval, unsigned int mode_m, double radius_m) const;
  std::vector<double> findRadialRoots(unsigned int mode_m, unsigned int n_roots) const;
  std::vector<double> legendreRootsAndWeights(unsigned int n_points, std::vector<double>& weights) const;

  void writeGarfieldRootFile(const std::vector<double>& r_edges_m,
                             const std::vector<double>& z_edges_m,
                             const std::vector<double>& er,
                             const std::vector<double>& ez,
                             unsigned int r_begin,
                             unsigned int r_end) const;

  void writeField3DRootFile(const std::vector<double>& r_source_edges_m,
                            const std::vector<double>& r_obs_edges_m,
                            const std::vector<double>& phi_source_edges,
                            const std::vector<double>& phi_obs_edges,
                            const std::vector<double>& z_source_edges_m,
                            const std::vector<double>& z_obs_edges_m,
                            const std::vector<double>& rho,
                            const std::vector<double>& potential,
                            const std::vector<double>& er,
                            const std::vector<double>& ephi,
                            const std::vector<double>& ez,
                            unsigned int r_begin,
                            unsigned int r_end) const;

  bool verifyOutput() const;

  double m_aCm{21.6};
  double m_bCm{76.4};
  double m_lCm{102.325};
  double m_sourceRMinCm{22.8};
  double m_sourceRMaxCm{75.43};

  double m_rhoReferenceNCPerM3{20.0};
  double m_kEff{1.0};
  double m_radialPowerAlpha{1.4};
  double m_m1Amplitude{0.0};
  double m_m1Phase{0.0};
  double m_m12Amplitude{0.0};
  double m_m12Phase{0.0};

  unsigned int m_nrSource{18};
  unsigned int m_nphiSource{24};
  unsigned int m_nzSource{20};
  unsigned int m_nrObs{24};
  unsigned int m_nphiObs{24};
  unsigned int m_nzObs{28};
  unsigned int m_mPhiMax{4};
  unsigned int m_nRadialModes{18};
  unsigned int m_nLongitudinalModes{22};
  bool m_autoAxisymmetric{true};

  unsigned int m_jobIndex{0};
  unsigned int m_nJobs{1};

  std::string m_garfieldOutputFile{"sphenix_rossegger_garfield_field_1p4.root"};
  std::string m_field3DOutputFile{"sphenix_rossegger_fields_3d.root"};
  bool m_writeField3D{false};
  bool m_verifyOutput{true};
  bool m_done{false};
};

#endif
