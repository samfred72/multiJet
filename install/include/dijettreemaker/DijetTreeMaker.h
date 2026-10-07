#ifndef DIJETTREEMAKER_H
#define DIJETTREEMAKER_H



#include <fun4all/SubsysReco.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <string>
#include <vector>
#include <globalvertex/GlobalVertex.h>
#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1Packetv1.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <calotrigger/LL1Out.h>
#include <calotrigger/LL1Outv1.h>
#include <calotrigger/TriggerPrimitive.h>
#include <calotrigger/TriggerPrimitivev1.h>
#include <calotrigger/TriggerPrimitiveContainer.h>
#include <calotrigger/TriggerPrimitiveContainerv1.h>
#include <calotrigger/TriggerDefs.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfov2.h>
#include <calobase/TowerInfov3.h>
#include <calobase/RawTowerGeom.h>// jet includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jetv2.h>

#include "TTree.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TVector3.h"

#include <phparameter/PHParameters.h>

class CentralityInfo;
class MinimumBiasInfo;
class TVector3;
class PHCompositeNode;
class PHParameters;
class TowerInfoContainer;
class JetContainer;
class PHG4TruthInfoContainer;
class TowerBackground;

class DijetTreeMaker : public SubsysReco
{
 public:

  DijetTreeMaker(const std::string &name = "DijetTreeMaker", const std::string &outfilename = "trees_caloemulator.root");

  virtual ~DijetTreeMaker();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int process_jets(int cone_index, PHCompositeNode *topNode);

  int process_truth_jets(int cone_index, PHCompositeNode *topNode);
  
  void GetNodes (PHCompositeNode *topNode);

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void UseEmulator(bool use) { m_useEmulator = use; }//  float get_eccentricity(std::vector<TVector3> hcaltowers, float oh_sum);
  void UsePsi(bool use) { m_use_psi = use; }
  void CalibrateJets(bool calib) { m_calibrate_jets = calib; }
  void SetVerbosity(int verbo) ;
  void SetTrigger(const std::string trigger) {_trigger = trigger;}
  void SetIsSim(bool use) {isSim = use;}
  void KeepTruth(bool use) {keep_truth = use;}
  void KeepUE(bool use) {keep_ue = use;}
  void SetPtCut(float pt) {pt_cut = pt;}
  void SetSubleadingPt(float pt) {subleading_pt_cut = pt;}
  void SetMaxPtCut(float pt) {max_pt_cut = pt;}
  void SetDijetCut(int b) {m_dijet_cut = b;}
  void SetPtCutTruth(float pt) {pt_cut_truth = pt;}
  void SaveAllEvents(bool all = true) {m_allevents = all;}
  Double_t getDPHI(Double_t phi1, Double_t phi2);
  // Truth-anchored JER smear: mean is always pt_calib (reco), but the resolution
  // width is looked up at the matched truth jet's pt when available, falling back
  // to pt_calib when unmatched (matches gammajet/treemaking/src/CaloAna.cc's
  // smear_pt() fix - see its comment there for why the fallback matters: without
  // it, unmatched jets got zero smearing instead of a reasonable one).
  // z is the jet's standard-normal deviate (shared by every variation); pt_ref is the
  // matched truth pt, or pt_calib when unmatched. sign: 0 = nominal, +1 = high, -1 = low.
  Double_t smear_pt(float pt_calib, float pt_ref, float z, int sign);
  // Seed for this event's smearing: run and event number, plus the truth vertex in MC
  // (MC event numbers restart in every DST segment). The same event always gets the same
  // smearing; different events and jobs get independent ones.
  UInt_t eventSeed(PHCompositeNode *topNode);
  Double_t DeltaR(float x1, float x2, float y1, float y2);
  void saveCalo(bool save){save_calo = save;}
  void isAuAu(bool auau){m_is_auau = auau;}
  void skim(bool save){m_skim = save;}
  void deepJet(bool save){m_deep_jet = save;}
  void setSingleJet(bool single){singlejet = single;}
  void doTimingCut(bool do_time) {m_do_timing_cut = do_time;}
  // Apply the same event-level dijet-candidate selection that
  // multiJet's offline skimmer/makeSkimmedTrees.C::check_dijet_reco used to
  // apply as a *second* pass over hadded trees. With this on, that selection
  // (+ the |mbd_vertex_z|<=60 cut) runs inline in process_event, so events
  // that would have been dropped by the offline skim are never written here
  // either - the single Fun4All pass IS the skimmed output. Overrides
  // SaveAllEvents()/the m_dijet_cut path when true.
  void ApplySkimCuts(bool use) { m_apply_skim_cuts = use; }
  /*
    cone size map
    0 - unsub
    1 - sub
    2 - all truth
    3 - embed truth 0
    4 - embed truth 1
    5 - embed truth 2
   */
  void addRecoConeSize(int conesize, int bkg = 0) { m_reco_cone_sizes[m_n_reco_cone_sizes++] = std::make_pair(conesize, bkg); }
  void addTruthConeSize(int conesize, int bkg = 0) { m_truth_cone_sizes[m_n_truth_cone_sizes++] = std::make_pair(conesize, bkg); }
  void setConepT(int conesize, float pT_threshold) { m_reco_cone_size_threshold[conesize] = pT_threshold; }

  void setRunnumber(int run){b_runnumber = run;}
  void setSegment(int seg){b_segment = seg;}

 private:
  int count = 0;

  bool m_do_timing_cut{true};
  bool m_deep_jet{false};
  bool m_use_psi{true};
  bool m_is_auau{false};
  bool m_skim{false};
  bool keep_ue{0};
  int m_n_reco_cone_sizes = 0;
  int m_n_truth_cone_sizes = 0;

  bool m_bad_event{false};
  std::vector<GlobalVertex::VTXTYPE> m_vtxtypes{};
  std::map<int, std::pair<int, float>> m_reco_cone_sizes{};
  std::map<int, float> m_reco_cone_size_threshold{};
  std::map<int, std::pair<int, int>> m_truth_cone_sizes{};
  bool singlejet{false};
  bool m_useEmulator{false};

  PHParameters _cutParams{"TimingCutParams"}; //variable name is arbitrary`

  TriggerAnalyzer *triggeranalyzer{nullptr};

  void reset_tree_vars();

  // Skim-mode storage trim: keeps only the top 3 jets by pt (calibrated pt
  // when available) for the given cone_index, at indices 0/1/2 = leading/
  // subleading/subsubleading, and drops the rest - applied to every per-jet
  // vector for that cone (pt, calib pt, smear variants, kinematics, etc.),
  // all reordered/truncated together so they stay index-aligned. Called
  // once at the end of process_jets(), only when m_skim is set.
  // A jet is kept if it is in the top n_keep of ANY pt flavor (calib or one of
  // the six JER smears), so a jet that only becomes leading after smearing survives.
  void keepLeadingJets(int cone_index, size_t n_keep = 4);
  // pt vectors of every flavor for a cone: calib, then (MC) the six JER smears
  std::vector<std::vector<float>*> ptFlavors(int cone_index);

  // --- combined single-pass "offline skim" selection ---
  // Mirrors multiJet's old skimmer/makeSkimmedTrees.C::check_dijet_reco, run
  // per reco cone size (radius) on the jet_pt_calib/eta/phi/t vectors this
  // module already fills, OR'd across radii, exactly like that macro's
  // per-radius anypass loop. See DijetTreeMaker::passesOfflineSkimCuts().
  bool m_apply_skim_cuts{false};
  // Skim thresholds sit below the analysis cuts (leading 20, jets 2 and 3 at 7 GeV), so the
  // analysis can re-smear or tighten without reprocessing.
  static constexpr float m_skim_leading_pt_cut = 15.;
  static constexpr float m_skim_subjet_pt_cut = 5.;
  static constexpr float m_skim_dphicut_loose = TMath::Pi() / 2.;  // makeSkimmedTrees.C dphicutloose
  static constexpr float m_skim_vertex_cut = 60.;         // makeSkimmedTrees.C vertex_cut
  // Data skim thresholds are multiplied by the lowest in-situ JES (p_a) the scans
  // try, so events that pass after jet_pt/p_a are kept (as unfolder.cc does).
  static constexpr float m_skim_jes_floor = 0.9;
  // Jets are stored if the raw pt >= pt_cut or any pt flavor >= this
  float m_store_pt_cut = 4;
  bool passesOfflineSkimCuts();

  bool m_allevents = false;
  int _verbosity;
  float m_higheta = 0;
  float m_loweta = 0;
  float m_higheta_truth = 0;
  float m_loweta_truth = 0;

  float dijetcut = 10.;
  bool m_no_cuts = true;
  bool m_dijet_cut = false;
  float max_pt_cut = 10.;
  float pt_cut = 4;
  float pt_cutCalib = 7;
  float subleading_pt_cut = 4;
  float pt_cut_truth = 4;
  bool save_calo = false;
  TFile *_f;
  TTree *_tree;
  TTree *_tree_end;
  TH1D *h_centrality;
  TH1D *h_jet_spectra[100];
  TH1D *h_jet_spectra_etacut[100];
  std::string _trigger;
  std::string _foutname;
  std::string _nodename;
  std::string m_calo_nodename;
  int _i_event;
  bool isSim{0};
  bool keep_truth{1};
  float m_energy_cut = 0.1;  

  float dphicut = 3.*TMath::Pi()/4.;
  TowerInfo *_tower;
  TowerInfoContainer* _towers;

  CentralityInfo *m_central = nullptr;
  MinimumBiasInfo *m_minimumbiasinfo = nullptr;
  TowerBackground *m_towBack1 = nullptr;
  TowerBackground *m_towBack2 = nullptr;

  bool m_has_tower_background = false;
  float m_background_v2 = 0;
  int m_background_strips = 0;
  float m_background_psi2 = 0;
  int m_background_fail = 0;
  float m_background_v2_sub1 = 0;
  int m_background_strips_sub1 = 0;
  float m_background_psi2_sub1 = 0;
  int m_background_fail_sub1 = 0;

  double b_calib_lead_time = -999;
  double b_calib_delta_time = -999;
  double b_calib_mbd_time = -999;
  
  int b_centrality = 0;
  float b_psi = 0;
  int b_ismb = 0;
  int m_calibrate_jets = 0;
  RawClusterContainer *clusters{nullptr};

  RawTowerGeomContainer *tower_geomIH = nullptr;
  RawTowerGeomContainer *tower_geomEM = nullptr;
  RawTowerGeomContainer *tower_geomOH = nullptr;

  TowerInfoContainer *hcalin_sub_towers = nullptr;
  TowerInfoContainer *hcalout_sub_towers = nullptr;
  TowerInfoContainer *emcalre_sub_towers = nullptr;
  TowerInfoContainer *hcalin_towers = nullptr;
  TowerInfoContainer *hcalout_towers = nullptr;
  TowerInfoContainer *emcalre_towers = nullptr;
  TowerInfoContainer *emcal_towers = nullptr;

  int b_runnumber = 0;
  int b_segment = 0;
  ULong64_t first_event;
  ULong64_t last_event;
  ULong64_t first_scaled[64];
  ULong64_t first_live[64];
  ULong64_t first_raw[64];
  ULong64_t last_scaled[64];
  ULong64_t last_live[64];
  ULong64_t last_raw[64];
  ULong64_t b_scaled_scalers[64];
  ULong64_t b_live_scalers[64];
  ULong64_t b_raw_scalers[64];
  
  float vtx_z = 0;
  float b_vertex_z;

  float b_truth_vertex_x;  
  float b_truth_vertex_y;
  float b_truth_vertex_z;
  double b_x1;
  double b_x2;
  int b_id1;
  int b_id2;
  double b_Q;
  float b_time_zero;
  int b_mbd_hit;
  float b_mbd_charge_sum;
  float b_emcal_energy_sum;
  float b_hcalin_energy_sum;
  float b_hcalout_energy_sum;
  int b_mbd_npmt;

  std::vector<float> b_mbd_charge;
  std::vector<float> b_mbd_time;
  std::vector<short> b_mbd_ipmt;
  std::vector<short> b_mbd_side;

  std::vector<short> b_emcal_good;
  std::vector<float> b_emcal_energy;
  std::vector<float> b_emcal_time;
  std::vector<float> b_emcal_etabin;
  std::vector<float> b_emcal_phibin;
  std::vector<float> b_emcal_eta;
  std::vector<float> b_emcal_phi;

  std::vector<short> b_hcalin_good;
  std::vector<float> b_hcalin_energy;
  std::vector<float> b_hcalin_time;
  std::vector<float> b_hcalin_etabin;
  std::vector<float> b_hcalin_phibin;
  std::vector<float> b_hcalin_eta;
  std::vector<float> b_hcalin_phi;

  std::vector<short> b_hcalout_good;
  std::vector<float> b_hcalout_energy;
  std::vector<float> b_hcalout_time;
  std::vector<float> b_hcalout_etabin;
  std::vector<float> b_hcalout_phibin;
  std::vector<float> b_hcalout_eta;
  std::vector<float> b_hcalout_phi;

  int b_truth_particle_n = 0;
  std::vector<int> b_truth_particle_pid = {};
  std::vector<float> b_truth_particle_pt = {};
  std::vector<float>  b_truth_particle_eta = {};
  std::vector<float>  b_truth_particle_phi = {};


  std::vector<float> *b_ue_1[3];
  std::vector<float> *b_ue_2[3];
  std::vector<int> b_njet{};
  std::vector<std::vector<float>*> b_jet_pt;
  std::vector<std::vector<float>*> b_jet_pt_calib;
  // JER-smeared pt, two families (see smear_pt() above): "_reco" is smeared
  // entirely from the reco/calib jet (mean and width both from pt_calib);
  // "_truth" keeps the reco mean but evaluates the resolution width at the
  // matched truth jet's pt when available. Each has nominal/high/low.
  std::vector<std::vector<float>*> b_jet_pt_smear_reco;
  std::vector<std::vector<float>*> b_jet_pt_smear_high_reco;
  std::vector<std::vector<float>*> b_jet_pt_smear_low_reco;
  std::vector<std::vector<float>*> b_jet_pt_smear_truth;
  std::vector<std::vector<float>*> b_jet_pt_smear_high_truth;
  std::vector<std::vector<float>*> b_jet_pt_smear_low_truth;
  // MC: the jet's standard-normal smearing deviate (every smeared variation is
  // pt_calib + z * pt_ref * width_variation(pt_ref)), and the matched truth jet pt (-1 if none)
  std::vector<std::vector<float>*> b_jet_smear_z;
  std::vector<std::vector<float>*> b_jet_truth_pt;
  // MC: one standard-normal deviate per event and radius, for smearing a summed recoil
  float b_recoil_smear_z[16] = {0};
  std::vector<std::vector<float>*> b_jet_et;
  std::vector<std::vector<float>*> b_jet_t;
  std::vector<std::vector<float>*> b_jet_e;
  std::vector<std::vector<float>*> b_jet_eta;
  std::vector<std::vector<float>*> b_jet_eta_det;
  std::vector<std::vector<float>*> b_jet_phi;

  std::vector<std::vector<float>*> b_jet_pt_unsub;
  std::vector<std::vector<float>*> b_jet_e_unsub;
  std::vector<std::vector<float>*> b_jet_eta_unsub;
  std::vector<std::vector<float>*> b_jet_phi_unsub;

  std::vector<std::vector<float>*> b_jet_emcal;
  std::vector<std::vector<float>*> b_jet_hcalin;
  std::vector<std::vector<float>*> b_jet_hcalout;
  std::vector<std::vector<float>*> b_jet_e_emcal;
  std::vector<std::vector<float>*> b_jet_e_hcalin;
  std::vector<std::vector<float>*> b_jet_e_hcalout;
  std::vector<std::vector<float>*> b_jet_rho;
  std::vector<std::vector<float>*> b_jet_rho_emcal;
  std::vector<std::vector<float>*> b_jet_rho_hcalin;
  std::vector<std::vector<float>*> b_jet_rho_hcalout;
  std::vector<std::vector<float>*> b_jet_rhophi;
  std::vector<std::vector<float>*> b_jet_rhophi_emcal;
  std::vector<std::vector<float>*> b_jet_rhophi_hcalin;
  std::vector<std::vector<float>*> b_jet_rhophi_hcalout;
  std::vector<std::vector<float>*> b_jet_rhoeta;
  std::vector<std::vector<float>*> b_jet_rhoeta_emcal;
  std::vector<std::vector<float>*> b_jet_rhoeta_hcalin;
  std::vector<std::vector<float>*> b_jet_rhoeta_hcalout;


  std::vector<int> b_truth_njet;
  std::vector<std::vector<float>*> b_truth_jet_pt;
  std::vector<std::vector<float>*> b_truth_jet_eta;
  std::vector<std::vector<float>*> b_truth_jet_phi;
  std::vector<std::vector<int>*> b_truth_jet_flavor;
  std::vector<std::vector<float>*> b_truth_jet_particle_pt;

  int b_ncluster = 0;
  std::vector<float> b_cluster_e = {};
  std::vector<float> b_cluster_ecore = {};
  std::vector<float> b_cluster_eta = {};
  std::vector<float> b_cluster_phi = {};
  std::vector<float> b_cluster_pt = {};
  std::vector<float> b_cluster_chi2 = {};
  std::vector<float> b_cluster_prob = {};

  int b_seed_njet;
  std::vector<float>* b_seed_jet_pt;
  std::vector<float>* b_seed_jet_eta;
  std::vector<float>* b_seed_jet_phi;
  std::vector<float>* b_seed_jet_D;
  std::vector<float>* b_seed_jet_max_eT;
  std::vector<float>* b_seed_jet_mean_eT;

  int b_seed_njet_sub;
  std::vector<float>* b_seed_jet_pt_sub;
  std::vector<float>* b_seed_jet_eta_sub;
  std::vector<float>* b_seed_jet_phi_sub;

  ULong64_t b_gl1_scaled{0};
  ULong64_t b_gl1_live{0};
  ULong64_t b_gl1_raw{0};

  std::array<int, 4> m_photon_emu_triggernames;
  std::array<int, 4> m_jet_emu_triggernames;


  double b_chargesum{0};
  int m_hitcut{2};
  double b_prodsigma{0};  
  double b_avgsigma{0};
  double b_maxsigma{0};
  double b_proddelta{0};
  double b_avgdelta{0};
  double b_maxdelta{0};
  double b_pileup{0};

  //smearing parameters
  TRandom3 rand;
  // JER smearing templates (fractional resolution vs pt), ported from
  // gammajet/treemaking/src/CaloAna.cc. That file only has templates derived
  // for r04 (jerband_smearing_templates.root has no other radius); by
  // decision, the same r04-derived curves are reused for every radius here
  // rather than falling back to a different methodology per radius.
  TH1D *h_jer_smear_nominal{nullptr};
  TH1D *h_jer_smear_up{nullptr};
  TH1D *h_jer_smear_down{nullptr};

};

#endif 
