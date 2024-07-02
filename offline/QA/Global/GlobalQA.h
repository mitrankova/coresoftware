#ifndef GLOBALQA_GLOBALQA_H
#define GLOBALQA_GLOBALQA_H

#include <fun4all/SubsysReco.h>

#include <string>

// Forward declarations
class PHCompositeNode;
class TH1;
class TH2;

class GlobalQA : public SubsysReco
{
 public:
  //! constructor
  GlobalQA(const std::string& name = "GlobalQA");  // const std::string &filename = "testQA.root"); //int nevents = 100);

  //! destructor
  ~GlobalQA() override;

  //! full initialization
  int Init(PHCompositeNode*) override;

  //! event processing method
  int process_event(PHCompositeNode*) override;

  //! end of run method
  int End(PHCompositeNode*) override;

  int process_g4hits(PHCompositeNode*);
  int process_g4cells(PHCompositeNode*);
  int process_towers(PHCompositeNode*);
  int process_clusters(PHCompositeNode*);

  void Detector(const std::string& name) { detector = name; }
  void set_timing_cut_width(const int& t) { _range = t; }

  void set_debug(bool debug) { m_debug = debug; }
  TH2* LogYHist2D(const std::string& name, const std::string& title, int, double, double, int, double, double);

 private:
  int evtcount = 0;
  int Getpeaktime(TH1* h);
  void createHistos();

  TH1* h_mbd_zvtx = nullptr;
  TH1* h_mbd_zvtx_wide = nullptr;
  TH1* h_calc_zvtx = nullptr;
  TH1* h_calc_zvtx_wide = nullptr;
  TH1* h_mbd_charge_s = nullptr;
  TH1* h_mbd_charge_n = nullptr;
  TH1* h_mbd_nhit_s = nullptr;
  TH1* h_mbd_nhit_n = nullptr;
  
  TH1* h_zdc_zvtx = nullptr;
  TH1* h_zdc_energy_s = nullptr;
  TH1* h_zdc_energy_n = nullptr;

  int _eventcounter{0};
  int _range{1};

  bool m_debug{false};

  std::string detector = "";
  std::string m_outputFileName = "";
  std::string OutputFileName = "";
};

#endif
