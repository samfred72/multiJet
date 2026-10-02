#include "TMath.h"
#include "TLorentzVector.h"
#include "DijetTreeMaker.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <HepMC/PdfInfo.h>

#pragma GCC diagnostic pop

//#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <centrality/CentralityInfov2.h>
#include <calotrigger/MinimumBiasInfo.h>
//for emc clusters
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>
#include <TMath.h>
#include <calotrigger/TriggerRunInfov1.h>
#include <eventplaneinfo/EventplaneinfoMap.h>
#include <eventplaneinfo/Eventplaneinfo.h>
#include <calobase/RawTowerGeom.h>
#include <jetbackground/TowerBackground.h>
#include <jetbase/Jet.h>
#include <jetbase/Jetv1.h>
#include <jetbase/Jetv2.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>
#include <HepMC/SimpleVector.h> 
//for vetex information

#include <globalvertex/GlobalVertexMap.h>
#include <vector>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
// G4Cells includes

#include <iostream>

#include <map>

//____________________________________________________________________________..
DijetTreeMaker::DijetTreeMaker(const std::string &name, const std::string &outfilename):
  SubsysReco(name)

{
  isSim = false;
  _foutname = outfilename;  
  m_calo_nodename = "TOWERINFO_CALIB";
  m_jet_emu_triggernames = {20, 21, 22, 23};

  m_photon_emu_triggernames = {28, 29, 30, 31};

  m_vtxtypes.push_back(GlobalVertex::VTXTYPE::MBD);
}

//____________________________________________________________________________..
DijetTreeMaker::~DijetTreeMaker()
{

}

//____________________________________________________________________________..
int DijetTreeMaker::Init(PHCompositeNode *topNode)
{

  // JER smearing templates - see the comment on h_jer_smear_nominal in the
  // header for why the same r04-derived curves are reused for every radius.
  {
    TFile fjer("/sphenix/user/samfred/projects/gammajet/treemaking/macros/jerband_smearing_templates.root");
    h_jer_smear_nominal = dynamic_cast<TH1D*>(fjer.Get("h_jer_smear_r04_pileup_EMfracJES_nominal"));
    h_jer_smear_up      = dynamic_cast<TH1D*>(fjer.Get("h_jer_smear_r04_pileup_EMfracJES_sysup"));
    h_jer_smear_down    = dynamic_cast<TH1D*>(fjer.Get("h_jer_smear_r04_pileup_EMfracJES_sysdown"));
    if (!h_jer_smear_nominal || !h_jer_smear_up || !h_jer_smear_down)
    {
      std::cout << "could not load JER smearing templates from jerband_smearing_templates.root" << std::endl;
      exit(-1);
    }
    h_jer_smear_nominal->SetDirectory(0);
    h_jer_smear_up     ->SetDirectory(0);
    h_jer_smear_down   ->SetDirectory(0);
  }

  h_centrality = new TH1D("h_centrality",";Centrality Bin;Events", 100, -0.5, 99.5);
  for (int i = 0; i < 100; i++)
  {
    h_jet_spectra[i] = new TH1D(Form("h_jet_spectra_%d", i),";p_{T};Counts", 100, 0, 50);
    h_jet_spectra_etacut[i] = new TH1D(Form("h_jet_spectra_etacut_%d",i),";p_{T};Counts", 100, 0, 50);
  }

  _f = new TFile( _foutname.c_str(), "RECREATE");

  std::cout << " making a file = " <<  _foutname.c_str() << " , _f = " << _f << std::endl;
  triggeranalyzer = new TriggerAnalyzer();
  triggeranalyzer->UseEmulator(m_useEmulator);

  if (!isSim)
  {
    _tree_end = new TTree("ttree_end","run stats");
    _tree_end->Branch("runnumber", &b_runnumber);
    _tree_end->Branch("segment", &b_segment);     
    _tree_end->Branch("scaled_scalers", &b_scaled_scalers, "scaled_scalers[64]/l");
    _tree_end->Branch("live_scalers", &b_live_scalers, "live_scalers[64]/l");
    _tree_end->Branch("raw_scalers", &b_raw_scalers, "raw_scalers[64]/l");
    _tree_end->Branch("scaled_last", &last_scaled, "scaled_last[64]/l");
    _tree_end->Branch("live_last", &last_live, "live_last[64]/l");
    _tree_end->Branch("raw_last", &last_raw, "raw_last[64]/l");
    _tree_end->Branch("scaled_first", &first_scaled, "scaled_first[64]/l");
    _tree_end->Branch("live_first", &first_live, "live_first[64]/l");
    _tree_end->Branch("raw_first", &first_raw, "raw_first[64]/l");
    _tree_end->Branch("first_event", &first_event, "first_event/l");
    _tree_end->Branch("last_event", &last_event, "last_event/l");
  }
  _tree = new TTree("ttree","a persevering date tree");
  if (isSim && keep_truth)
  {
    std::cout << "Initing truth jets" << std::endl;
    _tree->Branch("truth_vertex_x", &b_truth_vertex_x);
    _tree->Branch("truth_vertex_y", &b_truth_vertex_y);
    _tree->Branch("truth_vertex_z", &b_truth_vertex_z);
    _tree->Branch("gen_x1", &b_x1);
    _tree->Branch("gen_x2", &b_x2);
    _tree->Branch("gen_id1", &b_id1);
    _tree->Branch("gen_id2", &b_id2);
    _tree->Branch("gen_Q", &b_Q);
    for (auto cone_size : m_truth_cone_sizes)
    {
      int coneindex = cone_size.first;
      int cone = cone_size.second.first;
      int bkg = cone_size.second.second;
      std::cout << " Cone size : " << cone << " " << bkg << std::endl;
      b_truth_njet.push_back(0);
      b_truth_jet_pt.push_back(new std::vector<float>());
      b_truth_jet_flavor.push_back(new std::vector<int>());
      b_truth_jet_particle_pt.push_back(new std::vector<float>());
      b_truth_jet_eta.push_back(new std::vector<float>());
      b_truth_jet_phi.push_back(new std::vector<float>());
      if (bkg == 0)
      {
        _tree->Branch(Form("truth_jet_pt_%d", cone), b_truth_jet_pt.at(coneindex));
        _tree->Branch(Form("truth_jet_eta_%d", cone), b_truth_jet_eta.at(coneindex));
        _tree->Branch(Form("truth_jet_phi_%d", cone), b_truth_jet_phi.at(coneindex));
        //if (m_is_auau) _tree->Branch(Form("truth_jet_flavor_%d", cone), b_truth_jet_flavor.at(coneindex));
        //if (m_is_auau) _tree->Branch(Form("truth_jet_particle_pt_%d", cone), b_truth_jet_particle_pt.at(coneindex));
      }
      else
      {
        _tree->Branch(Form("truth_e%d_jet_pt_%d", bkg - 1, cone), b_truth_jet_pt.at(coneindex));
        _tree->Branch(Form("truth_e%d_jet_eta_%d", bkg - 1, cone), b_truth_jet_eta.at(coneindex));
        _tree->Branch(Form("truth_e%d_jet_phi_%d", bkg - 1,cone), b_truth_jet_phi.at(coneindex));
        if (m_is_auau)  _tree->Branch(Form("truth_e%d_jet_flavor_%d", bkg - 1, cone), b_truth_jet_flavor.at(coneindex));
        if (m_is_auau) _tree->Branch(Form("truth_e%d_jet_particle_pt_%d", bkg - 1, cone), b_truth_jet_particle_pt.at(coneindex));

      }
    }
  }
  else if (!isSim)
  {
    _tree->Branch("runnumber", &b_runnumber);
    if (!m_skim)  // not kept by makeSkimmedTrees.C's output tree
    {
      _tree->Branch("gl1_live", &b_gl1_live, "gl1_live/l");
      _tree->Branch("gl1_min_bias", &b_gl1_raw, "gl1_min_bias/l");
    }
  }
  _tree->Branch("gl1_scaled", &b_gl1_scaled, "gl1_scaled/l");

  _i_event = 0;

  if (save_calo)
  {
    _tree->Branch("emcal_good",&b_emcal_good);
    _tree->Branch("emcal_energy",&b_emcal_energy);
    _tree->Branch("emcal_time",&b_emcal_time);
    _tree->Branch("emcal_phibin",&b_emcal_phibin);
    _tree->Branch("emcal_etabin",&b_emcal_etabin);
    _tree->Branch("emcal_phi",&b_emcal_phi);
    _tree->Branch("emcal_eta",&b_emcal_eta);

    _tree->Branch("hcalin_good",&b_hcalin_good);
    _tree->Branch("hcalin_energy",&b_hcalin_energy);
    _tree->Branch("hcalin_time",&b_hcalin_time);
    _tree->Branch("hcalin_phibin",&b_hcalin_phibin);
    _tree->Branch("hcalin_etabin",&b_hcalin_etabin);
    _tree->Branch("hcalin_phi",&b_hcalin_phi);
    _tree->Branch("hcalin_eta",&b_hcalin_eta);

    _tree->Branch("hcalout_good",&b_hcalout_good);
    _tree->Branch("hcalout_energy",&b_hcalout_energy);
    _tree->Branch("hcalout_time",&b_hcalout_time);
    _tree->Branch("hcalout_phibin",&b_hcalout_phibin);
    _tree->Branch("hcalout_etabin",&b_hcalout_etabin);
    _tree->Branch("hcalout_phi",&b_hcalout_phi);
    _tree->Branch("hcalout_eta",&b_hcalout_eta);
  }

  if (m_is_auau)
  {
    _tree->Branch("centrality",&b_centrality);
    _tree->Branch("psi",&b_psi);
    _tree->Branch("minbias",&b_ismb);

    _tree->Branch("jet_isflowsub", &m_background_fail);
    _tree->Branch("jet_psi", &m_background_psi2);
    _tree->Branch("jet_v2", &m_background_v2);
    _tree->Branch("jet_strips", &m_background_strips);

    if (keep_ue)
    {	  
      //_tree->Branch("jet_isflowsub_sub1", &m_background_fail_sub1);
      //_tree->Branch("jet_psi_sub1", &m_background_psi2_sub1);
      //_tree->Branch("jet_v2_sub1", &m_background_v2_sub1);
      //_tree->Branch("jet_strips_sub1", &m_background_strips_sub1);

      for (int i = 0; i < 3; i++)
      {
        b_ue_1[i] = new std::vector<float>();
        b_ue_2[i] = new std::vector<float>();
        //_tree->Branch(Form("ue_1_lay_%d", i), b_ue_1[i]);
        _tree->Branch(Form("ue_2_lay_%d", i), b_ue_2[i]);
      }
    }
  }
  if (!m_skim)
  {
    _tree->Branch("mbd_npt", &b_mbd_npmt);
    _tree->Branch("mbd_time", &b_mbd_time);
    _tree->Branch("mbd_charge", &b_mbd_charge);
    _tree->Branch("mbd_side", &b_mbd_side);
    _tree->Branch("mbd_ipmt", &b_mbd_ipmt);
  }

  _tree->Branch("mbd_prodsigma", &b_prodsigma);
  _tree->Branch("mbd_avgsigma", &b_avgsigma);
  _tree->Branch("mbd_maxsigma", &b_maxsigma);
  _tree->Branch("mbd_proddelta", &b_proddelta);
  _tree->Branch("mbd_avgdelta", &b_avgdelta);
  _tree->Branch("mbd_maxdelta", &b_maxdelta);

  _tree->Branch("calib_lead_time", &b_calib_lead_time);
  _tree->Branch("calib_delta_time", &b_calib_delta_time);
  _tree->Branch("calib_mbd_time", &b_calib_mbd_time);

  if (m_is_auau)
  {
    _tree->Branch("mbd_charge_sum", &b_mbd_charge_sum);
    _tree->Branch("emcal_energy_sum", &b_emcal_energy_sum);
    _tree->Branch("hcalin_energy_sum", &b_hcalin_energy_sum);
    _tree->Branch("hcalout_energy_sum", &b_hcalout_energy_sum);
  }
  _tree->Branch("mbd_vertex_z", &b_vertex_z, "mbd_vertex_z/F");
  _tree->Branch("time_zero", &b_time_zero, "time_zero/F");
  _tree->Branch("mbd_hit", &b_mbd_hit, "mbd_hit/I");

  if (keep_ue && m_is_auau)
  {
    b_seed_njet = 0;

    b_seed_jet_pt = new std::vector<float>();
    b_seed_jet_eta = new std::vector<float>();
    b_seed_jet_phi = new std::vector<float>();
    b_seed_jet_D = new std::vector<float>();
    b_seed_jet_mean_eT = new std::vector<float>();
    b_seed_jet_max_eT = new std::vector<float>();


    _tree->Branch("jet_seed_pt", b_seed_jet_pt);
    _tree->Branch("jet_seed_eta", b_seed_jet_eta);
    _tree->Branch("jet_seed_phi", b_seed_jet_phi);
    _tree->Branch("jet_seed_D", b_seed_jet_D);
    _tree->Branch("jet_seed_mean_eT", b_seed_jet_mean_eT);
    _tree->Branch("jet_seed_max_eT", b_seed_jet_max_eT);

    b_seed_njet_sub = 0;

    b_seed_jet_pt_sub = new std::vector<float>();
    b_seed_jet_eta_sub = new std::vector<float>();
    b_seed_jet_phi_sub = new std::vector<float>();


    _tree->Branch("jet_seed_pt_sub", b_seed_jet_pt_sub);
    _tree->Branch("jet_seed_eta_sub", b_seed_jet_eta_sub);
    _tree->Branch("jet_seed_phi_sub", b_seed_jet_phi_sub);

  }

  for (auto cone_size : m_reco_cone_sizes){
    int coneindex = cone_size.first;
    int cone = cone_size.second.first;
    int bkg = cone_size.second.second;

    std::cout << " Cone size : " << cone << std::endl;
    b_njet.push_back(0);
    b_jet_pt.push_back(new std::vector<float>());
    b_jet_pt_calib.push_back(new std::vector<float>());
    b_jet_pt_smear_reco.push_back(new std::vector<float>());
    b_jet_pt_smear_high_reco.push_back(new std::vector<float>());
    b_jet_pt_smear_low_reco.push_back(new std::vector<float>());
    b_jet_pt_smear_truth.push_back(new std::vector<float>());
    b_jet_pt_smear_high_truth.push_back(new std::vector<float>());
    b_jet_pt_smear_low_truth.push_back(new std::vector<float>());
    b_jet_et.push_back(new std::vector<float>());
    b_jet_t.push_back(new std::vector<float>());
    b_jet_e.push_back(new std::vector<float>());
    b_jet_eta.push_back(new std::vector<float>());
    b_jet_eta_det.push_back(new std::vector<float>());
    b_jet_phi.push_back(new std::vector<float>());
    if (m_is_auau)
    {
      b_jet_pt_unsub.push_back(new std::vector<float>());
      b_jet_e_unsub.push_back(new std::vector<float>());
      b_jet_eta_unsub.push_back(new std::vector<float>());
      b_jet_phi_unsub.push_back(new std::vector<float>());
    }
    b_jet_emcal.push_back(new std::vector<float>());
    b_jet_hcalin.push_back(new std::vector<float>());
    b_jet_hcalout.push_back(new std::vector<float>());
    if (m_deep_jet)
    {
      b_jet_e_emcal.push_back(new std::vector<float>());
      b_jet_e_hcalin.push_back(new std::vector<float>());
      b_jet_e_hcalout.push_back(new std::vector<float>());

      b_jet_rho.push_back(new std::vector<float>());
      b_jet_rho_emcal.push_back(new std::vector<float>());
      b_jet_rho_hcalin.push_back(new std::vector<float>());
      b_jet_rho_hcalout.push_back(new std::vector<float>());
      b_jet_rhophi.push_back(new std::vector<float>());
      b_jet_rhophi_emcal.push_back(new std::vector<float>());
      b_jet_rhophi_hcalin.push_back(new std::vector<float>());
      b_jet_rhophi_hcalout.push_back(new std::vector<float>());
      b_jet_rhoeta.push_back(new std::vector<float>());
      b_jet_rhoeta_emcal.push_back(new std::vector<float>());
      b_jet_rhoeta_hcalin.push_back(new std::vector<float>());
      b_jet_rhoeta_hcalout.push_back(new std::vector<float>());
    }
    _tree->Branch(Form("jet_pt_%d%s", cone, (bkg?"_sub":"")), b_jet_pt.at(coneindex));
    if (m_calibrate_jets){
      _tree->Branch(Form("jet_pt_calib_%d%s", cone, (bkg?"_sub":"")), b_jet_pt_calib.at(coneindex));
      // For real data these are just copies of jet_pt_calib (smear_pt() is a
      // no-op when !isSim - see process_jets()); makeSkimmedTrees.C never
      // even read them for data, so skip branching the duplicates there too.
      if (isSim || !m_skim)
      {
        _tree->Branch(Form("jet_pt_smear_reco_%d%s", cone, (bkg?"_sub":"")), b_jet_pt_smear_reco.at(coneindex));
        _tree->Branch(Form("jet_pt_smear_high_reco_%d%s", cone, (bkg?"_sub":"")), b_jet_pt_smear_high_reco.at(coneindex));
        _tree->Branch(Form("jet_pt_smear_low_reco_%d%s", cone, (bkg?"_sub":"")), b_jet_pt_smear_low_reco.at(coneindex));
        _tree->Branch(Form("jet_pt_smear_truth_%d%s", cone, (bkg?"_sub":"")), b_jet_pt_smear_truth.at(coneindex));
        _tree->Branch(Form("jet_pt_smear_high_truth_%d%s", cone, (bkg?"_sub":"")), b_jet_pt_smear_high_truth.at(coneindex));
        _tree->Branch(Form("jet_pt_smear_low_truth_%d%s", cone, (bkg?"_sub":"")), b_jet_pt_smear_low_truth.at(coneindex));
      }
    }
    _tree->Branch(Form("jet_t_%d%s", cone, (bkg?"_sub":"")), b_jet_t.at(coneindex));
    _tree->Branch(Form("jet_e_%d%s", cone, (bkg?"_sub":"")), b_jet_e.at(coneindex));      
    _tree->Branch(Form("jet_eta_%d%s", cone, (bkg?"_sub":"")), b_jet_eta.at(coneindex));
    _tree->Branch(Form("jet_eta_det_%d%s", cone, (bkg?"_sub":"")), b_jet_eta_det.at(coneindex));
    _tree->Branch(Form("jet_phi_%d%s", cone, (bkg?"_sub":"")), b_jet_phi.at(coneindex));
    if (m_is_auau)
    {
      _tree->Branch(Form("jet_pt_unsub_%d%s", cone, (bkg?"_sub":"")), b_jet_pt_unsub.at(coneindex));
      _tree->Branch(Form("jet_e_unsub_%d%s", cone, (bkg?"_sub":"")), b_jet_e_unsub.at(coneindex));
      _tree->Branch(Form("jet_eta_unsub_%d%s", cone, (bkg?"_sub":"")), b_jet_eta_unsub.at(coneindex));
      _tree->Branch(Form("jet_phi_unsub_%d%s", cone, (bkg?"_sub":"")), b_jet_phi_unsub.at(coneindex));
    }
    _tree->Branch(Form("jet_emcal_%d%s", cone, (bkg?"_sub":"")), b_jet_emcal.at(coneindex));
    if (!m_skim)  // not kept by makeSkimmedTrees.C's output tree
    {
      _tree->Branch(Form("jet_hcalout_%d%s", cone, (bkg?"_sub":"")), b_jet_hcalout.at(coneindex));
    }

    if (m_deep_jet)
    {
      _tree->Branch(Form("jet_e_emcal_%d%s", cone, (bkg?"_sub":"")), b_jet_e_emcal.at(coneindex));
      _tree->Branch(Form("jet_e_hcalout_%d%s", cone, (bkg?"_sub":"")), b_jet_e_hcalout.at(coneindex));
      _tree->Branch(Form("jet_e_hcalin_%d%s", cone, (bkg?"_sub":"")), b_jet_e_hcalin.at(coneindex));

      _tree->Branch(Form("jet_rho_%d%s", cone, (bkg?"_sub":"")), b_jet_rho.at(coneindex));
      _tree->Branch(Form("jet_rhophi_%d%s", cone, (bkg?"_sub":"")), b_jet_rhophi.at(coneindex));
      _tree->Branch(Form("jet_rhoeta_%d%s", cone, (bkg?"_sub":"")), b_jet_rhoeta.at(coneindex));
      _tree->Branch(Form("jet_rho_emcal_%d%s", cone, (bkg?"_sub":"")), b_jet_rho_emcal.at(coneindex));
      _tree->Branch(Form("jet_rho_hcalin_%d%s", cone, (bkg?"_sub":"")), b_jet_rho_hcalin.at(coneindex));
      _tree->Branch(Form("jet_rho_hcalout_%d%s", cone, (bkg?"_sub":"")), b_jet_rho_hcalout.at(coneindex));
      //_tree->Branch(Form("jet_rhophi_%d%s", cone, (bkg?"_sub":"")), b_jet_rhophi.at(coneindex));
      _tree->Branch(Form("jet_rhophi_emcal_%d%s", cone, (bkg?"_sub":"")), b_jet_rhophi_emcal.at(coneindex));
      _tree->Branch(Form("jet_rhophi_hcalin_%d%s", cone, (bkg?"_sub":"")), b_jet_rhophi_hcalin.at(coneindex));
      _tree->Branch(Form("jet_rhophi_hcalout_%d%s", cone, (bkg?"_sub":"")), b_jet_rhophi_hcalout.at(coneindex));
      _tree->Branch(Form("jet_rhoeta_emcal_%d%s", cone, (bkg?"_sub":"")), b_jet_rhoeta_emcal.at(coneindex));
      _tree->Branch(Form("jet_rhoeta_hcalin_%d%s", cone, (bkg?"_sub":"")), b_jet_rhoeta_hcalin.at(coneindex));
      _tree->Branch(Form("jet_rhoeta_hcalout_%d%s", cone, (bkg?"_sub":"")), b_jet_rhoeta_hcalout.at(coneindex));
    }
  }
  if (!m_skim)
  {
    _tree->Branch("ncluster", &b_ncluster);
    _tree->Branch("cluster_ecore", &b_cluster_ecore);
    _tree->Branch("cluster_e", &b_cluster_e);
    _tree->Branch("cluster_pt", &b_cluster_pt);
    _tree->Branch("cluster_chi2", &b_cluster_chi2);
    _tree->Branch("cluster_prob", &b_cluster_prob);
    _tree->Branch("cluster_eta", &b_cluster_eta);
    _tree->Branch("cluster_phi", &b_cluster_phi);
  }
  std::cout << "Done initing the treemaker"<<std::endl;  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetTreeMaker::InitRun(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

void DijetTreeMaker::SetVerbosity(int verbo){
  _verbosity = verbo;
  return;
}

void DijetTreeMaker::reset_tree_vars()
{

  b_ismb = 0;
  b_centrality = -999;
  b_mbd_hit = 0;
  b_mbd_npmt = 0;
  b_mbd_charge_sum = 0;
  b_emcal_energy_sum = 0;
  b_hcalin_energy_sum = 0;
  b_hcalout_energy_sum = 0;
  b_mbd_charge.clear();
  b_mbd_time.clear();
  b_mbd_ipmt.clear();
  b_mbd_side.clear();

  b_ncluster = 0;
  b_cluster_e.clear();
  b_cluster_ecore.clear();
  b_cluster_eta.clear();
  b_cluster_phi.clear();
  b_cluster_pt.clear();
  b_cluster_chi2.clear();
  b_cluster_prob.clear();
  b_calib_lead_time = -999;
  b_calib_delta_time = -999;
  b_calib_mbd_time = -999;
  if (keep_ue && m_is_auau){
    for (int i = 0; i < 3; i++){
      b_ue_1[i]->clear();
      b_ue_2[i]->clear();
    }

    b_seed_njet = 0;
    b_seed_jet_pt->clear();
    b_seed_jet_phi->clear();
    b_seed_jet_eta->clear();
    b_seed_jet_D->clear();
    b_seed_jet_mean_eT->clear();
    b_seed_jet_max_eT->clear();

    b_seed_njet_sub = 0;
    b_seed_jet_pt_sub->clear();
    b_seed_jet_phi_sub->clear();
    b_seed_jet_eta_sub->clear();
  }

  for (int i = 0; i < m_n_reco_cone_sizes; i++)
  {
    b_njet.at(i) = 0;
    b_jet_pt.at(i)->clear();
    b_jet_pt_calib.at(i)->clear();
    b_jet_pt_smear_reco.at(i)->clear();
    b_jet_pt_smear_high_reco.at(i)->clear();
    b_jet_pt_smear_low_reco.at(i)->clear();
    b_jet_pt_smear_truth.at(i)->clear();
    b_jet_pt_smear_high_truth.at(i)->clear();
    b_jet_pt_smear_low_truth.at(i)->clear();
    b_jet_et.at(i)->clear();
    b_jet_t.at(i)->clear();
    b_jet_e.at(i)->clear();
    b_jet_eta.at(i)->clear();
    b_jet_eta_det.at(i)->clear();
    b_jet_phi.at(i)->clear();
    if (m_is_auau)
    {
      b_jet_pt_unsub.at(i)->clear();
      b_jet_e_unsub.at(i)->clear();
      b_jet_eta_unsub.at(i)->clear();
      b_jet_phi_unsub.at(i)->clear();
    }
    b_jet_emcal.at(i)->clear();
    b_jet_hcalin.at(i)->clear();
    b_jet_hcalout.at(i)->clear();
    if (m_deep_jet)
    {
      b_jet_e_emcal.at(i)->clear();
      b_jet_e_hcalin.at(i)->clear();
      b_jet_e_hcalout.at(i)->clear();

      b_jet_rho.at(i)->clear();
      b_jet_rho_emcal.at(i)->clear();
      b_jet_rho_hcalin.at(i)->clear();
      b_jet_rho_hcalout.at(i)->clear();
      b_jet_rhophi.at(i)->clear();
      b_jet_rhophi_emcal.at(i)->clear();
      b_jet_rhophi_hcalin.at(i)->clear();
      b_jet_rhophi_hcalout.at(i)->clear();
      b_jet_rhoeta.at(i)->clear();
      b_jet_rhoeta_emcal.at(i)->clear();
      b_jet_rhoeta_hcalin.at(i)->clear();
      b_jet_rhoeta_hcalout.at(i)->clear();
    }
  }
  if (save_calo)
  {
    b_emcal_good.clear();
    b_emcal_energy.clear();
    b_emcal_time.clear();
    b_emcal_phibin.clear();
    b_emcal_etabin.clear();
    b_emcal_phi.clear();
    b_emcal_eta.clear();

    b_hcalin_good.clear();
    b_hcalin_energy.clear();
    b_hcalin_time.clear();
    b_hcalin_phibin.clear();
    b_hcalin_etabin.clear();
    b_hcalin_phi.clear();
    b_hcalin_eta.clear();

    b_hcalout_good.clear();  
    b_hcalout_energy.clear();
    b_hcalout_time.clear();
    b_hcalout_phibin.clear();
    b_hcalout_etabin.clear();
    b_hcalout_phi.clear();
    b_hcalout_eta.clear();
  }
  b_gl1_scaled = 0x0;
  b_gl1_live = 0x0;
  b_gl1_raw = 0x0;

  if (isSim && keep_truth)
  {
    b_truth_vertex_x = 0;
    b_truth_vertex_y = 0;
    b_truth_vertex_z = 0;
    for (int i = 0; i < m_n_truth_cone_sizes; i++)
    {
      b_truth_njet.at(i) = 0;
      b_truth_jet_pt.at(i)->clear();
      b_truth_jet_flavor.at(i)->clear();
      b_truth_jet_particle_pt.at(i)->clear();
      b_truth_jet_eta.at(i)->clear();
      b_truth_jet_phi.at(i)->clear();
    }
    b_x1 = 0;//pdf->x1();
    b_x2 = 0;//pdf->x2();

    b_id1 = 0;//pdf->id1();  // incoming parton PDG
    b_id2 = 0;//pdf->id2();

    b_Q = 0;//pdf->scalePDF();
  }

  return;
}

//____________________________________________________________________________..
// Keeps only the top n_keep jets (by calibrated pt, falling back to raw pt
// if calibration wasn't applied for this cone) for the given cone_index -
// indices 0/1/2 after this call are leading/subleading/subsubleading.
// Every per-jet vector for that cone is reordered/truncated together so
// they stay index-aligned with each other.
void DijetTreeMaker::keepLeadingJets(int cone_index, size_t n_keep)
{
  std::vector<float> *pt = b_jet_pt.at(cone_index);
  size_t n = pt->size();
  if (n == 0) return;

  std::vector<float> *rankPt = (b_jet_pt_calib.at(cone_index)->size() == n) ? b_jet_pt_calib.at(cone_index) : pt;

  std::vector<size_t> order(n);
  for (size_t k = 0; k < n; k++) order[k] = k;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return rankPt->at(a) > rankPt->at(b); });
  // Even when n <= n_keep there's nothing to truncate, but the sort above
  // still has to run and be applied below: FastJet's own cluster-sequence
  // ordering isn't pt-sorted, so without this pass index 0/1/2 wouldn't
  // reliably be leading/subleading/subsubleading.
  if (order.size() > n_keep) order.resize(n_keep);

  auto reorder = [&order, n](std::vector<float> *v)
  {
    if (!v || v->size() != n) return; // not filled for this cone (e.g. calib/smear vectors when disabled)
    std::vector<float> kept;
    kept.reserve(order.size());
    for (size_t idx : order) kept.push_back(v->at(idx));
    *v = kept;
  };

  reorder(b_jet_pt.at(cone_index));
  reorder(b_jet_pt_calib.at(cone_index));
  reorder(b_jet_pt_smear_reco.at(cone_index));
  reorder(b_jet_pt_smear_high_reco.at(cone_index));
  reorder(b_jet_pt_smear_low_reco.at(cone_index));
  reorder(b_jet_pt_smear_truth.at(cone_index));
  reorder(b_jet_pt_smear_high_truth.at(cone_index));
  reorder(b_jet_pt_smear_low_truth.at(cone_index));
  reorder(b_jet_et.at(cone_index));
  reorder(b_jet_t.at(cone_index));
  reorder(b_jet_e.at(cone_index));
  reorder(b_jet_eta.at(cone_index));
  reorder(b_jet_eta_det.at(cone_index));
  reorder(b_jet_phi.at(cone_index));
  reorder(b_jet_emcal.at(cone_index));
  reorder(b_jet_hcalin.at(cone_index));
  reorder(b_jet_hcalout.at(cone_index));

  if (m_is_auau)
  {
    reorder(b_jet_pt_unsub.at(cone_index));
    reorder(b_jet_e_unsub.at(cone_index));
    reorder(b_jet_eta_unsub.at(cone_index));
    reorder(b_jet_phi_unsub.at(cone_index));
  }

  if (m_deep_jet)
  {
    reorder(b_jet_e_emcal.at(cone_index));
    reorder(b_jet_e_hcalin.at(cone_index));
    reorder(b_jet_e_hcalout.at(cone_index));
    reorder(b_jet_rho.at(cone_index));
    reorder(b_jet_rho_emcal.at(cone_index));
    reorder(b_jet_rho_hcalin.at(cone_index));
    reorder(b_jet_rho_hcalout.at(cone_index));
    reorder(b_jet_rhoeta.at(cone_index));
    reorder(b_jet_rhoeta_emcal.at(cone_index));
    reorder(b_jet_rhoeta_hcalin.at(cone_index));
    reorder(b_jet_rhoeta_hcalout.at(cone_index));
    reorder(b_jet_rhophi.at(cone_index));
    reorder(b_jet_rhophi_emcal.at(cone_index));
    reorder(b_jet_rhophi_hcalin.at(cone_index));
    reorder(b_jet_rhophi_hcalout.at(cone_index));
  }
}

//____________________________________________________________________________..
// Combined single-pass "offline skim" selection - faithful port of
// multiJet's old skimmer/makeSkimmedTrees.C::check_dijet_reco (run per
// reco cone size / radius, OR'd across radii, same as that macro's
// anypass loop), plus its |mbd_vertex_z|<=60 event cut. Runs on the
// jet_pt_calib/eta/phi/t vectors this module has just filled for the
// current event via process_jets(), so it needs no second pass over the
// data. Note: the old skimmer also computed a GL1 trigger-bit-22 flag
// ("passes_Trigger") but never actually used it to cut events - that is
// intentionally NOT enforced here either, to keep this selection
// identical to what the offline skim actually did (not what it computed).
bool DijetTreeMaker::passesOfflineSkimCuts()
{
  if (std::isnan(b_vertex_z) || fabs(b_vertex_z) > m_skim_vertex_cut)
  {
    return false;
  }

  struct SkimJet
  {
    float pt, eta, phi, t;
  };

  for (auto cone_size : m_reco_cone_sizes)
  {
    int coneindex = cone_size.first;
    int bkg = cone_size.second.second;
    if (bkg != 0) continue;  // makeSkimmedTrees.C only ever looked at the unsub jet_* branches

    std::vector<SkimJet> jets;
    jets.reserve(b_jet_pt_calib.at(coneindex)->size());
    for (size_t j = 0; j < b_jet_pt_calib.at(coneindex)->size(); j++)
    {
      float pt = b_jet_pt_calib.at(coneindex)->at(j);
      if (pt < pt_cutCalib) continue;  // makeSkimmedTrees.C calib_cut prefilter
      SkimJet sj;
      sj.pt = pt;
      sj.eta = b_jet_eta.at(coneindex)->at(j);
      sj.phi = b_jet_phi.at(coneindex)->at(j);
      sj.t = isSim ? 0.f : b_jet_t.at(coneindex)->at(j);  // makeSkimmedTrees.C: isMC ? 0 : jet_t
      jets.push_back(sj);
    }

    if (jets.size() < 3) continue;

    std::sort(jets.begin(), jets.end(), [](const SkimJet &a, const SkimJet &b) { return a.pt > b.pt; });

    const SkimJet &leading = jets.at(0);
    const SkimJet &subleading = jets.at(1);
    const SkimJet &subsubleading = jets.at(2);

    float dphir = getDPHI(leading.phi, subleading.phi);
    float dphirr = getDPHI(leading.phi, subsubleading.phi);

    if (!(leading.pt >= m_skim_leading_pt_cut &&
          subleading.pt >= pt_cutCalib &&
          subsubleading.pt >= pt_cutCalib &&
          dphir >= dphicut &&
          dphirr >= m_skim_dphicut_loose))
    {
      continue;
    }

    double jetdeltatime = 17.6 * (leading.t - subleading.t);
    double jetleadtime = 17.6 * (leading.t);
    bool passleadtime = (TMath::Abs(jetleadtime + 2.0) < 6.0);
    bool passdijettime = (TMath::Abs(jetdeltatime) < 3.0);

    if (!(passleadtime && passdijettime)) continue;

    return true;  // matches makeSkimmedTrees.C: anypass = true; break;
  }

  return false;
}

//____________________________________________________________________________..
int DijetTreeMaker::process_event(PHCompositeNode *topNode)
{

  if (count % 10 == 0) std::cout << "Event: " << count << std::endl;
  count++;
  PHNodeIterator itNode(topNode);
  PHCompositeNode* parNode = dynamic_cast<PHCompositeNode*>(itNode.findFirst("PHCompositeNode","PAR"));
  PdbParameterMap* flagNode; 
  flagNode = findNode::getClass<PdbParameterMap>(parNode, "TimingCutParams"); 

  m_no_cuts = false;

  if(flagNode)
  {
    _cutParams.FillFrom(flagNode);
  }
  else
  {
    m_no_cuts = true;
  }


  bool dijet_candidate = false;
  m_bad_event = false;

  reset_tree_vars();
  if (!isSim)
  {
    Gl1Packet* gl1p = findNode::getClass<Gl1Packet>(topNode, 14001);
    if (!gl1p)
    {
      std::cout << "no gl1 packet " << std::endl;
      exit(-1);
    }
    //Gl1Packet* gl1p = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    b_gl1_scaled = gl1p->getScaledVector();
    b_gl1_live = gl1p->getLiveVector();
    b_gl1_raw = gl1p->lValue(10, 1);

    if (_i_event == 0)
    {
      first_event = gl1p->getEvtSequence();
      for (int i = 0; i < 64; i++)
      {
        first_raw[i] = gl1p->lValue(i, 0);
        first_live[i] = gl1p->lValue(i, 1);
        first_scaled[i] = gl1p->lValue(i, 2);
      }
    }
    else
    {
      last_event = gl1p->getEvtSequence();
      for (int i = 0; i < 64; i++)
      {
        last_raw[i] = gl1p->lValue(i, 0);
        last_live[i] = gl1p->lValue(i, 1);
        last_scaled[i] = gl1p->lValue(i, 2);
      }
    }
  }

  if (m_is_auau)
  {
    // if (!isSim)
    // 	{
    // 	  if (( (b_gl1_scaled >> 10) & 0x1U ) != 0x1U)
    // 	    {
    // 	      return Fun4AllReturnCodes::EVENT_OK;
    // 	    }
    // 	}
    m_minimumbiasinfo = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
    m_central = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");  

    EventplaneinfoMap *m_evpmap = findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
    b_psi = 0;
    if (!m_evpmap)
    {
      m_use_psi = false;
    }
    if (m_use_psi)
    {
      if(!(m_evpmap->empty()))
      {
        auto EPDNS = m_evpmap->get(EventplaneinfoMap::sEPDNS);
        b_psi = EPDNS->get_shifted_psi(2);
      }
    }
    if (m_central)
      b_centrality = (m_central->has_centile(CentralityInfo::PROP::mbd_NS)?m_central->get_centrality_bin(CentralityInfo::PROP::mbd_NS) : -99);
    if (m_minimumbiasinfo)
      b_ismb = (m_minimumbiasinfo->isAuAuMinimumBias()? 1 : 0);


    m_towBack1 = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub1");

    m_towBack2 = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");


    // if (!b_ismb)
    // 	return Fun4AllReturnCodes::EVENT_OK;

    //h_centrality->Fill(b_centrality);

    m_background_fail = 0;
    m_background_v2 = -99;
    m_background_strips = 0;
    m_background_psi2 = -99;

    m_background_fail_sub1 = 0;
    m_background_v2_sub1 = -99;
    m_background_strips_sub1 = 0;
    m_background_psi2_sub1 = -99;

    if (m_towBack2 && m_towBack1)
    {
      m_has_tower_background = true;
      m_background_v2 = m_towBack2->get_v2();
      m_background_v2_sub1 = m_towBack1->get_v2();
      m_background_psi2 = m_towBack2->get_Psi2();
      m_background_psi2_sub1 = m_towBack1->get_Psi2();
      m_background_strips = m_towBack2->get_nStripsUsedForFlow();
      //m_background_fail = (m_towBack2->get_flow_failure_flag()? 1: 0);
      m_background_strips_sub1 = m_towBack1->get_nStripsUsedForFlow();
      //m_background_fail_sub1 = (m_towBack1->get_flow_failure_flag()? 1: 0);
      if (keep_ue)
      {
        for (int i = 0; i < 3; i++)
        {
          std::vector<float> layerUE1 = m_towBack1->get_UE(i);
          for (float ue : layerUE1)
          {
            b_ue_1[i]->push_back(ue);
          }
          std::vector<float> layerUE2 = m_towBack2->get_UE(i);
          for (float ue : layerUE2)
          {
            b_ue_2[i]->push_back(ue);
          }
        }
      }
    }
  }

  _i_event++;

  if (isSim)
  {
    triggeranalyzer->decodeTriggers(topNode);

    for (int j = 0; j < 4; j++)
    {
      if (triggeranalyzer->didTriggerFire(m_photon_emu_triggernames[j]))
      {
        unsigned int bit = 28 + j;
        b_gl1_scaled |= 0x1 << bit;
      }
      if (triggeranalyzer->didTriggerFire(m_jet_emu_triggernames[j]))
      {
        unsigned int bit = 20 + j;
        b_gl1_scaled |= 0x1 << bit;
      }
    }
    b_gl1_scaled &= 0x00000000ffffffff;
  }

  clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");

  tower_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");


  vtx_z = 0;
  b_vertex_z = -999;

  std::vector<GlobalVertex*> vertices = vertexmap->get_gvtxs_with_type(m_vtxtypes);
  if(!vertices.empty())
  {
    if(vertices.at(0))
    {
      vtx_z = vertices.at(0)->get_z();
      b_vertex_z = vtx_z;
      b_time_zero = vertices.at(0)->get_t();
    }
    if(vertices.size() > 1 && Verbosity() > 0)
    {
      std::cout << "TowerJetInput::WARNING!! More than one vertex of selected type!" << std::endl;
    }
  }

  if (std::isnan(b_vertex_z))
  {
    b_vertex_z = -999;
    vtx_z = 0;
  }

  //std::cout << "MBD Vertex: " << b_vertex_z << std::endl;

  if (Verbosity() > 5) std::cout << "Mbd Vertex = " << b_vertex_z << std::endl;


  hcalin_sub_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_HCALIN");
  hcalout_sub_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_HCALOUT");
  emcalre_sub_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_CEMC_RETOWER");

  if (m_is_auau)
  {
    hcalin_sub_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_HCALIN_SUB1");
    hcalout_sub_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_HCALOUT_SUB1");
    emcalre_sub_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_CEMC_RETOWER_SUB1");
    hcalin_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_HCALIN");
    hcalout_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_HCALOUT");
    emcalre_towers = findNode::getClass<TowerInfoContainer>(topNode, m_calo_nodename + "_CEMC_RETOWER");
    if (!hcalin_towers)
    {
      std::cout << "no hcalin towers" << std::endl;
      exit(-1);
    }
    if (!hcalout_towers)
    {
      std::cout << "no hcalout towers" << std::endl;
      exit(-1);
    }
    if (!emcalre_towers)
    {
      std::cout << "no emcalre towers" << std::endl;
      exit(-1);
    }
    if (!hcalin_sub_towers)
    {
      std::cout << "no hcalin sub towers" << std::endl;
      exit(-1);
    }
    if (!hcalout_sub_towers)
    {
      std::cout << "no hcalout sub towers" << std::endl;
      exit(-1);
    }
    if (!emcalre_sub_towers)
    {
      std::cout << "no emcalre sub towers" << std::endl;
      exit(-1);
    }

  }

  //const RawTowerDefs::keytype eekey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, 0);
  //const RawTowerDefs::keytype ohkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, 0, 0);
  // float em_r = tower_geomEM->get_tower_geometry(eekey)->get_center_radius();
  // float oh_r = tower_geomOH->get_tower_geometry(ohkey)->get_center_radius();

  // float em_boundary = 1.1;		   
  // float em_theta = 2*atan(exp(-1*em_boundary));
  // float em_z = em_r/tan(em_theta);
  // float pos_theta = atan(em_r/(em_z - b_vertex_z));
  // float neg_theta = atan(em_r/(em_z + b_vertex_z));

  // float oh_z = oh_r/tan(em_theta);
  // float pos_oh_theta = atan(oh_r/(oh_z - b_vertex_z));
  // float neg_oh_theta = atan(oh_r/(oh_z + b_vertex_z));

  // float new_em_low_boundary = log(tan(neg_theta/2));
  // float new_em_high_boundary = -1*log(tan(pos_theta/2));
  // float new_oh_low_boundary = log(tan(neg_oh_theta/2));
  // float new_oh_high_boundary = -1*log(tan(pos_oh_theta/2));

  m_loweta = -1.1;//std::max(new_em_low_boundary, new_oh_low_boundary);

  m_higheta = 1.1;//std::max(new_em_high_boundary, new_oh_high_boundary);

  if (isSim && keep_truth)
  {

    PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfo)
    {   
      PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
      b_truth_vertex_z = gvertex->get_z();
      b_truth_vertex_x = gvertex->get_x();
      b_truth_vertex_y = gvertex->get_y();


      // float tpos_em_theta = atan(em_r/(em_z - b_truth_vertex_z));
      // float tneg_em_theta = atan(em_r/(em_z + b_truth_vertex_z));
      // float tpos_oh_theta = atan(oh_r/(oh_z - b_truth_vertex_z));
      // float tneg_oh_theta = atan(oh_r/(oh_z + b_truth_vertex_z));

      // float tnew_em_low_boundary = log(tan(tneg_em_theta/2));
      // float tnew_em_high_boundary = -1*log(tan(tpos_em_theta/2));
      // float tnew_oh_low_boundary = log(tan(tneg_oh_theta/2));
      // float tnew_oh_high_boundary = -1*log(tan(tpos_oh_theta/2));

      m_loweta_truth = -1.1;//std::max(tnew_em_low_boundary, tnew_oh_low_boundary);

      m_higheta_truth = 1.1;//std::max(tnew_em_high_boundary, tnew_oh_high_boundary);

    }  
    for (int i = 0; i < m_n_truth_cone_sizes; i++)
    {
      dijet_candidate |= process_truth_jets(i, topNode);
    }      
  }

  for (int i = 0; i < m_n_reco_cone_sizes; i++)
  {
    dijet_candidate |= process_jets(i, topNode);
  }

  if (m_apply_skim_cuts)
  {
    // Combined single-pass mode: this IS the offline skim's event
    // selection, so it alone decides whether to keep the event -
    // SaveAllEvents()/m_dijet_cut/m_bad_event from the old paths above
    // are not consulted here.
    if (!passesOfflineSkimCuts())
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }
  else if ((!dijet_candidate || m_bad_event) && !m_allevents)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (keep_ue)
  {
    std::string recoSeed_JetName = "AntiKt_TowerInfo_HIRecoSeedsRaw_r02";
    JetContainer *seed_jets = findNode::getClass<JetContainer>(topNode, recoSeed_JetName);
    if (seed_jets)
    {
      for (auto jet : *seed_jets)
      {

        if (jet->get_pt() < 5) continue;

        b_seed_njet++;

        b_seed_jet_pt->push_back(jet->get_pt());
        b_seed_jet_eta->push_back(jet->get_eta());
        b_seed_jet_phi->push_back(jet->get_phi());

        std::map<int, double> constituent_ETsum;

        for (const auto &comp : jet->get_comp_vec())
        {
          int comp_ieta = -1;
          int comp_iphi = -1;
          float comp_ET = 0;
          int comp_isBad = -99;
          TowerInfo *towerinfo;
          RawTowerGeom *tower_geom;
          if (comp.first == 5 || comp.first == 26)
          {
            towerinfo = hcalin_towers->get_tower_at_channel(comp.second);
            unsigned int towerkey = hcalin_towers->encode_key(comp.second);
            comp_ieta = hcalin_towers->getTowerEtaBin(towerkey);
            comp_iphi = hcalin_towers->getTowerPhiBin(towerkey);
            const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, comp_ieta, comp_iphi);
            tower_geom = tower_geomIH->get_tower_geometry(key);
            comp_ET = towerinfo->get_energy() / cosh(tower_geom->get_eta());
            comp_isBad = towerinfo->get_isHot() || towerinfo->get_isNoCalib() || towerinfo->get_isNotInstr() || towerinfo->get_isBadChi2();
          }
          else if (comp.first == 7 || comp.first == 27)
          {
            towerinfo = hcalout_towers->get_tower_at_channel(comp.second);
            unsigned int towerkey = hcalout_towers->encode_key(comp.second);
            comp_ieta = hcalout_towers->getTowerEtaBin(towerkey);
            comp_iphi = hcalout_towers->getTowerPhiBin(towerkey);
            const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, comp_ieta, comp_iphi);
            tower_geom = tower_geomOH->get_tower_geometry(key);
            comp_ET = towerinfo->get_energy() / cosh(tower_geom->get_eta());
            comp_isBad = towerinfo->get_isHot() || towerinfo->get_isNoCalib() || towerinfo->get_isNotInstr() || towerinfo->get_isBadChi2();
          }
          else if (comp.first == 13 || comp.first == 28)
          {
            towerinfo = emcalre_towers->get_tower_at_channel(comp.second);
            unsigned int towerkey = emcalre_towers->encode_key(comp.second);
            comp_ieta = emcalre_towers->getTowerEtaBin(towerkey);
            comp_iphi = emcalre_towers->getTowerPhiBin(towerkey);
            const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, comp_ieta, comp_iphi);
            tower_geom = tower_geomIH->get_tower_geometry(key);
            comp_ET = towerinfo->get_energy() / cosh(tower_geom->get_eta());
            comp_isBad = towerinfo->get_isHot() || towerinfo->get_isNoCalib() || towerinfo->get_isNotInstr() || towerinfo->get_isBadChi2();
          }
          if (comp_isBad)
          {
            continue;
          }
          int comp_ikey = (1000 * comp_ieta) + comp_iphi;


          constituent_ETsum[comp_ikey] += comp_ET;

        }

        float constituent_max_ET = 0;
        float constituent_sum_ET = 0;
        int nconstituents = 0;
        for (auto &map_iter : constituent_ETsum)
        {
          nconstituents++;
          constituent_sum_ET += map_iter.second;
          constituent_max_ET = std::max<double>(map_iter.second, constituent_max_ET);

        }


        float mean_constituent_ET = constituent_sum_ET / nconstituents;
        float seed_D = constituent_max_ET / mean_constituent_ET;
        //float rho = b_rho_val;//* static_cast<float>(nconstituents);
        b_seed_jet_D->push_back(seed_D);
        b_seed_jet_mean_eT->push_back(mean_constituent_ET);
        b_seed_jet_max_eT->push_back(constituent_max_ET);
      }
    }

    //TowerBackground *m_towBack2 = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub2");

    // float m_background_fail = 0;
    // float m_background_v2 = -99;
    // float m_background_psi2 = -99;
    // bool m_has_tower_background = false;
    // if (m_towBack2)
    //   {
    //     m_has_tower_background = true;
    //     m_background_v2 = m_towBack2->get_v2();
    //     m_background_psi2 = m_towBack2->get_Psi2();
    //     m_background_fail = (m_towBack2->get_flow_failure_flag()? 1: 0);
    //   }

    recoSeed_JetName = "AntiKt_TowerInfo_HIRecoSeedsSub_r02";
    JetContainer *seed_jets2 = findNode::getClass<JetContainer>(topNode, recoSeed_JetName);
    if (seed_jets2)
    {
      for (auto jet : *seed_jets2)
      {

        if (jet->get_pt() < 5) continue;

        b_seed_jet_pt_sub->push_back(jet->get_pt());
        b_seed_jet_eta_sub->push_back(jet->get_eta());
        b_seed_jet_phi_sub->push_back(jet->get_phi());
      }
    }

  }

  int size;
  if (emcal_towers && save_calo)
  {
    size = emcal_towers->size(); //online towers should be the same!
    for (int channel = 0; channel < size;channel++)
    {
      _tower = emcal_towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      //if (energy < m_energy_cut) continue;
      float time = _tower->get_time();
      unsigned int towerkey = emcal_towers->encode_key(channel);
      int ieta = emcal_towers->getTowerEtaBin(towerkey);
      int iphi = emcal_towers->getTowerPhiBin(towerkey);
      short good = (_tower->get_isGood() ? 1:0);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
      float tower_phi = tower_geomEM->get_tower_geometry(key)->get_phi();
      float tower_r = tower_geomEM->get_tower_geometry(key)->get_center_radius();
      float tower_z = tower_geomEM->get_tower_geometry(key)->get_center_z();
      float tower_theta = atan(tower_r/(tower_z - b_vertex_z));
      if (tower_theta < 0) tower_theta = TMath::Pi() + tower_theta;

      float tower_eta = -1*log(tan(tower_theta/2));


      // float tower_eta = -1*log(tan(atan(tower_r/(tower_z - b_vertex_z))/2));
      // float tower_eta = -1*log(tan(atan(tower_r/(tower_z - b_vertex_z))/2));

      b_emcal_good.push_back(good);
      b_emcal_energy.push_back(energy);
      b_emcal_time.push_back(time);
      b_emcal_etabin.push_back(ieta);
      b_emcal_phibin.push_back(iphi);
      b_emcal_eta.push_back(tower_eta);
      b_emcal_phi.push_back(tower_phi);
      if (good)
        b_emcal_energy_sum += energy;

    }
  }

  if (hcalin_towers && save_calo)
  {

    size = hcalin_towers->size(); //online towers should be the same!
    for (int channel = 0; channel < size;channel++)
    {
      _tower = hcalin_towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      //	  if (energy < m_energy_cut) continue;
      float time = _tower->get_time();	  
      short good = (_tower->get_isGood() ? 1:0);
      unsigned int towerkey = hcalin_towers->encode_key(channel);
      int ieta = hcalin_towers->getTowerEtaBin(towerkey);
      int iphi = hcalin_towers->getTowerPhiBin(towerkey);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
      float tower_phi = tower_geomIH->get_tower_geometry(key)->get_phi();
      float tower_r = tower_geomIH->get_tower_geometry(key)->get_center_radius();
      float tower_z = tower_geomIH->get_tower_geometry(key)->get_center_z();

      float tower_theta = atan(tower_r/(tower_z - b_vertex_z));
      if (tower_theta < 0) tower_theta = TMath::Pi() + tower_theta;

      float tower_eta = -1*log(tan(tower_theta/2));
      // float tower_eta = -1*log(tan(atan(tower_r/(tower_z - b_vertex_z))/2));
      // float tower_eta = tower_geomIH->get_tower_geometry(key)->get_eta();
      b_hcalin_good.push_back(good);
      b_hcalin_energy.push_back(energy);
      b_hcalin_time.push_back(time);
      b_hcalin_etabin.push_back(ieta);
      b_hcalin_phibin.push_back(iphi);
      b_hcalin_eta.push_back(tower_eta);
      b_hcalin_phi.push_back(tower_phi);
      if (good)
        b_hcalin_energy_sum += energy;
    }
  }
  if (hcalout_towers && save_calo)
  {

    size = hcalout_towers->size(); //online towers should be the same!
    for (int channel = 0; channel < size;channel++)
    {
      _tower = hcalout_towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      //if (energy < m_energy_cut) continue;
      float time = _tower->get_time();	  
      unsigned int towerkey = hcalout_towers->encode_key(channel);
      int ieta = hcalout_towers->getTowerEtaBin(towerkey);
      int iphi = hcalout_towers->getTowerPhiBin(towerkey);
      short good = (_tower->get_isGood() ? 1:0);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
      float tower_phi = tower_geomOH->get_tower_geometry(key)->get_phi();
      float tower_r = tower_geomOH->get_tower_geometry(key)->get_center_radius();
      float tower_z = tower_geomOH->get_tower_geometry(key)->get_center_z();

      float tower_theta = atan(tower_r/(tower_z - b_vertex_z));
      if (tower_theta < 0) tower_theta = TMath::Pi() + tower_theta;

      float tower_eta = -1*log(tan(tower_theta/2));
      // float tower_eta = -1*log(tan(atan(tower_r/(tower_z - b_vertex_z))/2));
      // float tower_eta = tower_geomOH->get_tower_geometry(key)->get_eta();


      b_hcalout_good.push_back(good);
      b_hcalout_energy.push_back(energy);
      b_hcalout_time.push_back(time);
      b_hcalout_etabin.push_back(ieta);
      b_hcalout_phibin.push_back(iphi);
      b_hcalout_eta.push_back(tower_eta);
      b_hcalout_phi.push_back(tower_phi);
      if (good)
        b_hcalout_energy_sum += energy;
    }
  }


  if (clusters && !m_skim)
  {
    CLHEP::Hep3Vector vertex(0, 0, b_vertex_z);

    auto clusterrange = clusters->getClusters();
    for (auto iclus = clusterrange.first; iclus != clusterrange.second; ++iclus)
    {

      RawCluster *cluster = iclus->second;
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      float energy = E_vec_cluster.mag();
      float pt = E_vec_cluster.perp();
      float clus_eta = E_vec_cluster.pseudoRapidity();
      float clus_phi = E_vec_cluster.phi();
      float chi2 = cluster->get_chi2();
      float prob = cluster->get_prob();
      float ecore = cluster->get_ecore();
      if (energy < 1) continue;
      b_ncluster++;
      b_cluster_e.push_back(energy);
      b_cluster_ecore.push_back(ecore);
      b_cluster_eta.push_back(clus_eta);
      b_cluster_phi.push_back(clus_phi);
      b_cluster_pt.push_back(pt);
      b_cluster_chi2.push_back(chi2);
      b_cluster_prob.push_back(prob);
    }
  }

  MbdPmtContainer *pmts_mbd = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");

  b_chargesum = 0;
  b_prodsigma = 0;
  b_avgsigma = 0;
  b_maxsigma = 0;
  b_proddelta = 0;
  b_avgdelta = 0;
  b_maxdelta = 0;

  int npmt[2] = {0};
  if (pmts_mbd)
  {
    int hit_north = 0;
    int hit_south = 0;

    double rmstime[2] = {0};
    double maxtime[2] = {0};
    double mintime[2] = {0};
    double sumtime[2] = {0};
    for (int i = 0 ; i < 128; i++)
    {
      MbdPmtHit *tmp_pmt = pmts_mbd->get_pmt(i);

      float charge = tmp_pmt->get_q();
      float time = tmp_pmt->get_time();

      if (fabs(time) < 25. && charge > 0.4)
      {
        b_chargesum += charge;
        int side = i/64;
        npmt[side]++;
        if (time > maxtime[side])
          maxtime[side] = time;
        if (time < mintime[side])
          mintime[side] = time;

        rmstime[side] += time*time;
        sumtime[side] += time;

        if (i < 64) hit_south = 1;
        else  hit_north = 1;
        b_mbd_charge_sum += charge;
        b_mbd_charge.push_back(charge);
        b_mbd_time.push_back(time);
        b_mbd_ipmt.push_back(i);
        b_mbd_side.push_back(i/64);
        b_mbd_npmt++;
      }
    }

    if (npmt[0] >= m_hitcut)
    {
      rmstime[0]/=(float)npmt[0];
      sumtime[0]/=(float)npmt[0];
      rmstime[0] = sqrt(rmstime[0] -sumtime[0]*sumtime[0]);
    }
    else
    {
      maxtime[0] = 0;
      mintime[0] = 0;
      rmstime[0] = 0;
      sumtime[0] = 0;
    }

    if (npmt[1] >= m_hitcut)
    {
      rmstime[1]/=(float)npmt[1];
      sumtime[1]/=(float)npmt[1];
      rmstime[1] = sqrt(rmstime[1] -sumtime[1]*sumtime[1]);
    }
    else
    {
      maxtime[1] = 0;
      mintime[1] = 0;
      rmstime[1] = 0;
      sumtime[1] = 0;
    }

    b_prodsigma = rmstime[0]*rmstime[1];
    b_avgsigma = (rmstime[0] + rmstime[1])/2.;
    b_maxsigma = std::max(rmstime[0], rmstime[1]);

    b_proddelta = (maxtime[1] - mintime[1])*(maxtime[0] - mintime[0]);
    b_avgdelta = ((maxtime[1] - mintime[1]) + (maxtime[0] - mintime[0]))/2.;
    b_maxdelta = std::max(maxtime[1] - mintime[1], maxtime[0] - mintime[0]);

    b_mbd_hit = (hit_south && hit_north ? 1 : 0);
  }


  // UE for each layer and each eta slice for both iterations

  // each layer
  _tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}
int DijetTreeMaker::process_truth_jets(int cone_index, PHCompositeNode* topNode)
{
  float max_secondmax_pt[2]={0};

  bool dijet_candidate = false;
  int isbkg = m_truth_cone_sizes[cone_index].second;
  int cone_size = m_truth_cone_sizes[cone_index].first;
  float jet_radius = cone_size*0.1;

  std::string truthJetName = Form("AntiKt_Truth_r%02d", cone_size);

  if (isbkg > 0)
  {
    truthJetName = Form("AntiKt_Truth_e%d_r%02d", isbkg - 1, cone_size);
  }

  JetContainer *jetstruth = findNode::getClass<JetContainer>(topNode, truthJetName);

  PHHepMCGenEventMap* geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  HepMC::GenEvent* hepmc_event = nullptr;
  if (geneventmap && false)
  {
    PHHepMCGenEvent* genevt = geneventmap->get(2);
    if (genevt) hepmc_event = genevt->getEvent();

    const HepMC::PdfInfo* pdf = hepmc_event->pdf_info();
    if (pdf) {
      b_x1 = pdf->x1();
      b_x2 = pdf->x2();

      b_id1 = pdf->id1();  // incoming parton PDG
      b_id2 = pdf->id2();

      b_Q = pdf->scalePDF();

    }
  }


  if (jetstruth)
  {
    // zero out counters

    for (auto jet : *jetstruth)
    {

      if (jet->get_pt() < pt_cut_truth) continue;
      //if (jet->get_eta() > (m_higheta_truth - cone_size*0.1) || jet->get_eta() < (m_loweta_truth + cone_size*0.1)) continue;
      b_truth_jet_pt.at(cone_index)->push_back(jet->get_pt());
      b_truth_jet_eta.at(cone_index)->push_back(jet->get_eta());
      b_truth_jet_phi.at(cone_index)->push_back(jet->get_phi());

      if (Verbosity() >= 2)
      {
        std::cout << "pt/eta/phi = " << jet->get_pt() << " / " << jet->get_eta() << " / " << jet->get_phi() << std::endl; 
      }

      for (int im = 0; im < 2; im++)
      {
        if (jet->get_pt() > max_secondmax_pt[im])
        {
          if (im == 0) max_secondmax_pt[1] = max_secondmax_pt[0];
          max_secondmax_pt[im] = jet->get_pt();
          break;
        }
      }

      float jet_eta = jet->get_eta();
      float jet_phi = jet->get_phi();
      int flavor = 0;
      float max_pt = 0;
      if (m_is_auau && false)
      {
        for (HepMC::GenEvent::particle_const_iterator p = hepmc_event->particles_begin();
            p != hepmc_event->particles_end(); ++p)
        {
          HepMC::GenParticle *particle = *p;
          if (!particle) continue;

          int pid = abs(particle->pdg_id());
          int status = particle->status();
          // Select outgoing partons from hard scattering (status 23) or
          // partons before hadronization (status 21, 22)
          if (status != 23 && status != 21 && status != 22) continue;

          // Only consider quarks (1-6) and gluons (21)
          if (!(pid >= 1 && pid <= 6) && pid != 21) continue;

          HepMC::FourVector momentum = particle->momentum();
          float part_pt = sqrt(momentum.px() * momentum.px() + momentum.py() * momentum.py());
          float part_eta = momentum.eta();
          float part_phi = momentum.phi();

          // Require minimum pT for parton matching
          if (part_pt < 3.0) continue;

          // Calculate angular distance between parton and jet

          if (part_phi > TMath::Pi())
          {
            part_phi = 2*TMath::Pi() - part_phi;
          }
          float deta = fabs(jet_eta - part_eta);
          float dphi = fabs(jet_phi - part_phi);

          if (dphi > TMath::Pi())
          {
            dphi = 2*TMath::Pi() - dphi;
          }
          float dr = sqrt(deta*deta + dphi*dphi);

          // Match parton to jet if within jet radius and has highest pT
          if (dr < jet_radius && part_pt > max_pt)
          {
            max_pt = part_pt;
            flavor = pid;
          }	   	  
        }
      }
      b_truth_jet_flavor.at(cone_index)->push_back(flavor);
      b_truth_jet_particle_pt.at(cone_index)->push_back(max_pt);
    }

    if (max_secondmax_pt[0] > pt_cut_truth) dijet_candidate = true;

    for (int im = 0; im < 2; im++)
    {
      max_secondmax_pt[im]=0;
    }
  }

  return dijet_candidate;
}
int DijetTreeMaker::process_jets(int cone_index, PHCompositeNode* topNode)
{
  float max_secondmax_pt[2]={0};
  float max_secondmax_phi[2]={0};
  float max_secondmax_t[2]={0};
  int isbkg = m_reco_cone_sizes[cone_index].second;
  int cone_size = m_reco_cone_sizes[cone_index].first;
  float max_pt_cone = m_reco_cone_size_threshold[cone_size];
  bool dijet_candidate = false;

  //if (m_is_auau) && (b_centrality < 0 || b_centrality > 99)) return 0;

  if (!m_no_cuts)
  {
    b_calib_lead_time = _cutParams.get_double_param("corrMaxJett");
    b_calib_delta_time = _cutParams.get_double_param("corrSubJett");
    b_calib_mbd_time = _cutParams.get_double_param("mbd_time");
  }
  if (isbkg > 2) return 0;

  std::string recoJetName = Form("AntiKt_TowerInfo_r%02d", cone_size);
  std::string calibJetName = Form("AntiKt_TowerInfo_r%02d_calib", cone_size);

  if (isbkg == 1)
  {
    recoJetName = Form("AntiKt_TowerInfo_r%02d_Sub1", cone_size);
  }


  JetContainer *jets = findNode::getClass<JetContainer>(topNode, recoJetName);
  JetContainer *jets_calib = findNode::getClass<JetContainer>(topNode, calibJetName);
  if (!jets)
  {
    std::cout << " no reco_jets" << std::endl;
    exit(-1);
  }
  if (!jets_calib && m_calibrate_jets)
  {
    std::cout << " no calib jets" << std::endl;
    exit(-1);
  }
  if (jets)
  {
    if (Verbosity() >= 2){
      std::cout << "Going through " << recoJetName << std::endl;
    }

    // Truth matching
    std::string truthJetName = Form("AntiKt_Truth_r%02d", cone_size);
    JetContainer * truthjets = findNode::getClass<JetContainer>(topNode, truthJetName);
    std::vector<std::pair<int,float>> matched_jets;
    
    if (isSim) {
      struct Pair { float pt_r, pt_t, dr; bool isvalid = 1; int i_r, i_t; };
      std::vector<std::vector<Pair>> pairs(jets->size(), std::vector<Pair>(truthjets->size()));
      std::vector<Pair> pair_list;

      int ireco = 0;
      for (auto jet : *jets) {
        
        Jet *jet_r = jets_calib->get_jet(ireco);
        
        int itruth = 0;
        if (jet_r->get_pt() < pt_cutCalib) {ireco++; continue; }
        for (auto jet_t : *truthjets) {
          if (jet_t->get_pt() < 5) {itruth++; continue; }
          Pair temp;
          temp.pt_r = jet_r->get_pt();
          temp.pt_t = jet_t->get_pt();
          temp.dr = DeltaR(jet_r->get_eta(), jet_r->get_phi(), jet_t->get_eta(), jet_t->get_phi());
          temp.i_r = ireco;
          temp.i_t = itruth;
          temp.isvalid = 1;
          pairs[ireco][itruth] = temp;
          pair_list.push_back(temp);
          //std::cout << std::fixed << std::setprecision(2) << pairs[ireco][itruth].dr << "," << pairs[ireco][itruth].pt_r << "," << pairs[ireco][itruth].pt_t << " ";
          itruth++;
        }
        //std::cout << std::endl;
        ireco++;
      }
      //std::cout << endl;
      std::sort(pair_list.begin(), pair_list.end(),
          [](const Pair& a, const Pair& b) {
          return a.dr < b.dr;
          });

      for (int ip = 0; ip < pair_list.size(); ip++) {
        Pair p = pair_list.at(ip);
        float dr = p.dr;
        if (dr >= cone_size / 10.0 * 3.0/4.0) break;

        if (!pairs[p.i_r][p.i_t].isvalid) continue;

        matched_jets.push_back(std::make_pair(p.i_r, p.pt_t));
        for (int ireco = 0; ireco < jets->size(); ireco++) {
          pairs[ireco][p.i_t].isvalid = 0;
        }
        for (int itruth = 0; itruth < truthjets->size(); itruth++) {
          pairs[p.i_r][itruth].isvalid = 0;
        }
      }

      //for (int a = 0; a < matched_jets.size(); a++) {
      //  cout << matched_jets.at(a).first << " " << matched_jets.at(a).second << endl;
      //}
      //std::cout << endl;
    }

    // zero out counters
    int ijet = 0;
    for (auto jet : *jets)
    {

      float calibpt = 0;
      if (m_calibrate_jets && !isbkg){	      
        Jet *jet_calib = jets_calib->get_jet(ijet);
        if (!jet_calib){
          std::cout << "no calib jet" << std::endl;
          exit(-1);
        }
        calibpt = jet_calib->get_pt();
      }

      // "_reco": mean and width both from the calibrated reco jet (matches
      // CaloAna.cc's pt_smear_reco/pt_smear_high_reco/pt_smear_low_reco).
      float smearpT_reco       = rand.Gaus(calibpt, calibpt * h_jer_smear_nominal->Interpolate(calibpt));
      float smearpT_high_reco  = rand.Gaus(calibpt, calibpt * h_jer_smear_up->Interpolate(calibpt));
      float smearpT_low_reco   = rand.Gaus(calibpt, calibpt * h_jer_smear_down->Interpolate(calibpt));
      // "_truth": reco mean, but the resolution width is looked up at the
      // matched truth jet's pt (falls back to calibpt if unmatched) - see
      // smear_pt()'s header comment.
      float smearpT_truth      = smear_pt(calibpt, ijet, matched_jets,  0);
      float smearpT_high_truth = smear_pt(calibpt, ijet, matched_jets, +1);
      float smearpT_low_truth  = smear_pt(calibpt, ijet, matched_jets, -1);

      if (!isSim) {
        smearpT_reco       = calibpt;
        smearpT_high_reco  = calibpt;
        smearpT_low_reco   = calibpt;
        smearpT_truth      = calibpt;
        smearpT_high_truth = calibpt;
        smearpT_low_truth  = calibpt;
      }

      ijet++;
      if (jet->get_pt() < pt_cut && calibpt < pt_cutCalib &&
          smearpT_reco < pt_cutCalib && smearpT_high_reco < pt_cutCalib && smearpT_low_reco < pt_cutCalib &&
          smearpT_truth < pt_cutCalib && smearpT_high_truth < pt_cutCalib && smearpT_low_truth < pt_cutCalib){
        continue;
      }
      //if (jet->get_pt() < 3) std::cout << cone_size << " " << jet->get_pt() << "*******************************" << std::endl;


      //if (jet->get_eta() > (m_higheta-0.1*cone_size) || jet->get_eta() < (m_loweta + 0.1*cone_size)) continue;
      // Total constituents

      int n_comp_total = 0;
      int n_comp_ihcal = 0;
      int n_comp_ohcal = 0;
      int n_comp_emcal = 0;

      // energy without event subtraction

      float jet_total_eT = 0;
      float jet_total_eX = 0;
      float jet_total_eY = 0;
      float jet_total_eZ = 0;
      float jet_total_e = 0;

      float jet_total_eT_unsub = 0;
      float jet_total_eX_unsub = 0;
      float jet_total_eY_unsub = 0;
      float jet_total_eZ_unsub = 0;
      float jet_total_e_unsub = 0;

      float eFrac_ihcal = 0;
      float eFrac_ohcal = 0;
      float eFrac_emcal = 0;

      float eTFrac_ihcal = 0;
      float eTFrac_ohcal = 0;
      float eTFrac_emcal = 0;

      float energy_hcalin = 0;
      float energy_hcalout = 0;
      float energy_emcal = 0;

      float rho_emcal = 0;
      float rho_hcalin = 0;
      float rho_hcalout = 0;

      float rhoeta_emcal = 0;
      float rhoeta_hcalin = 0;
      float rhoeta_hcalout = 0;
      float rhophi_emcal = 0;
      float rhophi_hcalin = 0;
      float rhophi_hcalout = 0;

      float jet_eta = jet->get_eta();
      float jet_phi = jet->get_phi();

      const RawTowerDefs::keytype tkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0,0);
      const RawTowerDefs::keytype tikey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, 0,0);
      const RawTowerDefs::keytype tokey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, 0,0);
      double r = tower_geomEM->get_tower_geometry(tkey)->get_center_radius();
      double z0 = sinh(jet_eta) * r;
      double z = z0 + b_vertex_z;
      double det_eta = asinh(z / r);  // eta after shift from vertex
      double emcal_r = r;
      double hcalin_r = tower_geomIH->get_tower_geometry(tikey)->get_center_radius();
      double hcalout_r = tower_geomOH->get_tower_geometry(tokey)->get_center_radius();

      double jet_t_avg = 0;
      double jet_t_Esum = 0;

      for (auto comp : jet->get_comp_vec())
      {

        unsigned int channel = comp.second;
        TowerInfo *tower;
        float tower_eT = 0;
        float tower_e = 0;
        if (comp.first == 26 || comp.first == 30)
        {  // IHcal

          tower = hcalin_sub_towers->get_tower_at_channel(channel);

          if (!tower || !tower_geomIH)
          {
            continue;
          }

          if (!tower->get_isGood()) continue;

          unsigned int calokey = hcalin_sub_towers->encode_key(channel);

          int ieta = hcalin_sub_towers->getTowerEtaBin(calokey);
          int iphi = hcalin_sub_towers->getTowerPhiBin(calokey);

          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
          float tower_phi = tower_geomIH->get_tower_geometry(key)->get_phi();
          float old_tower_eta = tower_geomIH->get_tower_geometry(key)->get_eta();

          float tower_z = hcalin_r*sinh(old_tower_eta);
          float tower_eta = asinh((tower_z - vtx_z)/hcalin_r);
          // tower energy and tower_eT

          tower_e = tower->get_energy();
          tower_eT = tower->get_energy() / std::cosh(tower_eta);

          jet_total_e += tower_e;
          jet_total_eT += tower_eT;
          jet_total_eX += tower_eT * cos(tower_phi);
          jet_total_eY += tower_eT * sin(tower_phi);
          jet_total_eZ += tower_eT * sinh(tower_eta);

          eTFrac_ihcal += tower_eT;
          eFrac_ihcal += tower_e;

          float deta = fabs(jet_eta - tower_eta);
          float dphi = fabs(jet_phi - tower_phi);

          if (dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;

          rho_hcalin += sqrt(dphi*dphi + deta*deta) * ( tower_e );
          rhophi_hcalin += dphi * ( tower_e );
          rhoeta_hcalin += deta * ( tower_e );
          energy_hcalin += tower_e;

          // get rid of UE subtraction if there
          if(tower_e > 0.1){
            jet_t_Esum += tower_e;
            jet_t_avg += tower_e*tower->get_time();
          }

          if (isbkg && m_is_auau)
          {
            float UE = m_towBack2->get_UE(1).at(ieta);
            if (m_background_v2 != -99)
              UE = UE * (1 + 2 * m_background_v2 * cos(2 * (tower_phi - m_background_psi2)));

            tower_eT = (tower->get_energy() + UE) / cosh(tower_eta);
            tower_e = (tower->get_energy() + UE);
            jet_total_e_unsub += tower_e;
            jet_total_eT_unsub += tower_eT;
            jet_total_eX_unsub += tower_eT * cos(tower_phi);
            jet_total_eY_unsub += tower_eT * sin(tower_phi);
            jet_total_eZ_unsub += tower_eT * sinh(tower_eta);

          }


          n_comp_ihcal++;
          n_comp_total++;
        }
        else if (comp.first == 27 || comp.first == 31)
        {  // OHCAL

          tower = hcalout_sub_towers->get_tower_at_channel(channel);
          if (!tower || !tower_geomOH)
          {
            continue;
          }
          if (!tower->get_isGood()) continue;
          unsigned int calokey = hcalout_sub_towers->encode_key(channel);

          int ieta = hcalout_sub_towers->getTowerEtaBin(calokey);
          int iphi = hcalout_sub_towers->getTowerPhiBin(calokey);

          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
          float tower_phi = tower_geomOH->get_tower_geometry(key)->get_phi();
          float old_tower_eta = tower_geomOH->get_tower_geometry(key)->get_eta();

          float tower_z = hcalout_r*sinh(old_tower_eta);
          float tower_eta = asinh((tower_z - vtx_z)/hcalout_r);
          // tower energy and tower_eT
          tower_e = tower->get_energy();
          tower_eT = tower->get_energy() / std::cosh(tower_eta);

          jet_total_e += tower_e;
          jet_total_eT += tower_eT;
          jet_total_eX += tower_eT * cos(tower_phi);
          jet_total_eY += tower_eT * sin(tower_phi);
          jet_total_eZ += tower_eT * sinh(tower_eta);


          float deta = fabs(jet_eta - tower_eta);
          float dphi = fabs(jet_phi - tower_phi);

          if (dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;

          rho_hcalout += sqrt(dphi*dphi + deta*deta) * ( tower_e );
          rhophi_hcalout += dphi * ( tower_e );
          rhoeta_hcalout += deta * ( tower_e );

          energy_hcalout += tower_e;

          eFrac_ohcal += tower_e;
          eTFrac_ohcal += tower_eT;

          // get rid of UE subtraction if there
          if(tower_e > 0.1)
          {
            jet_t_Esum += tower_e;
            jet_t_avg += tower_e*tower->get_time();
          }

          if (isbkg && m_is_auau)
          {
            float UE = m_towBack2->get_UE(2).at(ieta);
            if (m_background_v2 != -99)
              UE = UE * (1 + 2 * m_background_v2 * cos(2 * (tower_phi - m_background_psi2)));
            tower_eT = (tower->get_energy() + UE) / cosh(tower_eta);
            tower_e = (tower->get_energy() + UE);
            jet_total_e_unsub += tower_e;
            jet_total_eT_unsub += tower_eT;
            jet_total_eX_unsub += tower_eT * cos(tower_phi);
            jet_total_eY_unsub += tower_eT * sin(tower_phi);
            jet_total_eZ_unsub += tower_eT * sinh(tower_eta);

          }

          n_comp_ohcal++;
          n_comp_total++;
        }	       
        else if (comp.first == 28 || comp.first == 29)
        {  // EMCAL

          tower = emcalre_sub_towers->get_tower_at_channel(channel);

          if (!tower || !tower_geomIH)
          {
            continue;
          }
          if (!tower->get_isGood()) continue;

          unsigned int calokey = emcalre_sub_towers->encode_key(channel);

          int ieta = emcalre_sub_towers->getTowerEtaBin(calokey);
          int iphi = emcalre_sub_towers->getTowerPhiBin(calokey);

          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
          float tower_phi = tower_geomIH->get_tower_geometry(key)->get_phi();
          float old_tower_eta = tower_geomIH->get_tower_geometry(key)->get_eta();

          float tower_z = emcal_r*sinh(old_tower_eta);
          float tower_eta = asinh((tower_z - vtx_z)/emcal_r);
          // tower energy and tower_eT
          tower_e = tower->get_energy();
          tower_eT = tower->get_energy() / std::cosh(tower_eta);

          jet_total_e += tower_e;
          jet_total_eT += tower_eT;
          jet_total_eX += tower_eT * cos(tower_phi);
          jet_total_eY += tower_eT * sin(tower_phi);
          jet_total_eZ += tower_eT * sinh(tower_eta);

          float deta = fabs(jet_eta - tower_eta);
          float dphi = fabs(jet_phi - tower_phi);

          if (dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;

          rho_emcal += sqrt(dphi*dphi + deta*deta) * ( tower_e );
          rhophi_emcal += dphi * ( tower_e );
          rhoeta_emcal += deta * ( tower_e );

          energy_emcal += tower_e;

          eFrac_emcal += tower_e;
          eTFrac_emcal += tower_eT;

          // get rid of UE subtraction if there
          if(tower_e > 0.1)
          {
            jet_t_Esum += tower_e;
            jet_t_avg += tower_e*tower->get_time();
          }

          if (isbkg && m_is_auau)
          {
            float UE = m_towBack2->get_UE(0).at(ieta);
            if (m_background_v2 != -99)
              UE = UE * (1 + 2 * m_background_v2 * cos(2 * (tower_phi - m_background_psi2)));
            tower_eT = (tower->get_energy() + UE) / cosh(tower_eta);
            tower_e = (tower->get_energy() + UE);
            jet_total_e_unsub += tower_e;
            jet_total_eT_unsub += tower_eT;
            jet_total_eX_unsub += tower_eT * cos(tower_phi);
            jet_total_eY_unsub += tower_eT * sin(tower_phi);
            jet_total_eZ_unsub += tower_eT * sinh(tower_eta);

          }

          n_comp_emcal++;
          n_comp_total++;
        }	       

      }
      jet_t_avg /= jet_t_Esum;

      if (Verbosity() >= 2)
      {
        std::cout << "Reco " << cone_size << " JET: " << ijet << std::endl;
        std::cout << "Energy = " << jet_total_eT << " ( " << eFrac_emcal << "// " << eFrac_ihcal << "// " << eFrac_ohcal << std::endl;
        std::cout << "pt/eta/phi = " << jet->get_pt() << " / " << jet_eta << " / " << jet_phi << std::endl; 
      }

      eTFrac_ohcal /= jet_total_eT;
      eTFrac_emcal /= jet_total_eT;
      eTFrac_ihcal /= jet_total_eT;

      eFrac_ohcal /= jet_total_e;
      eFrac_emcal /= jet_total_e;
      eFrac_ihcal /= jet_total_e;

      if (jet_total_e < 0) continue;


      if (isbkg && m_is_auau)
      {

        Jet* unsubjet = new Jetv2();
        unsubjet->set_px(jet_total_eX_unsub);
        unsubjet->set_py(jet_total_eY_unsub);
        unsubjet->set_pz(jet_total_eZ_unsub);
        unsubjet->set_e(jet_total_e_unsub);
        if (jet_total_e_unsub < 0) continue;


        // if (fabs(jet->get_eta()) < 1.1)
        // 	{
        // 	  h_jet_spectra[b_centrality]->Fill(jet->get_pt());
        // 	}
        // if (fabs(jet->get_eta()) < 0.8)
        // 	{
        // 	  h_jet_spectra_etacut[b_centrality]->Fill(jet->get_pt());
        // 	}

        b_jet_pt_unsub.at(cone_index)->push_back(unsubjet->get_pt());
        b_jet_e_unsub.at(cone_index)->push_back(jet_total_e_unsub);
        b_jet_eta_unsub.at(cone_index)->push_back(unsubjet->get_eta());
        b_jet_phi_unsub.at(cone_index)->push_back(unsubjet->get_phi());


        delete unsubjet;

      }

      //if (jet->get_pt() < 3) std::cout << jet->get_pt() << std::endl;
      b_jet_pt.at(cone_index)->push_back(jet->get_pt());
      if (!isbkg && m_calibrate_jets){ 
        b_jet_pt_calib.at(cone_index)->push_back(calibpt);

        b_jet_pt_smear_reco.at(cone_index)->push_back(smearpT_reco);
        b_jet_pt_smear_high_reco.at(cone_index)->push_back(smearpT_high_reco);
        b_jet_pt_smear_low_reco.at(cone_index)->push_back(smearpT_low_reco);
        b_jet_pt_smear_truth.at(cone_index)->push_back(smearpT_truth);
        b_jet_pt_smear_high_truth.at(cone_index)->push_back(smearpT_high_truth);
        b_jet_pt_smear_low_truth.at(cone_index)->push_back(smearpT_low_truth);
      }
      b_jet_eta.at(cone_index)->push_back(jet->get_eta());
      b_jet_eta_det.at(cone_index)->push_back(det_eta);
      b_jet_phi.at(cone_index)->push_back(jet->get_phi());
      b_jet_e.at(cone_index)->push_back(jet_total_e);
      b_jet_t.at(cone_index)->push_back(jet_t_avg);
      b_jet_et.at(cone_index)->push_back(jet_total_eT);
      b_jet_emcal.at(cone_index)->push_back(eTFrac_emcal);
      b_jet_hcalin.at(cone_index)->push_back(eTFrac_ihcal);
      b_jet_hcalout.at(cone_index)->push_back(eTFrac_ohcal);

      if (m_deep_jet)
      {
        float rho_jet = rho_emcal + rho_hcalin + rho_hcalout;
        float rhoeta_jet = rhoeta_emcal + rhoeta_hcalin + rhoeta_hcalout;
        float rhophi_jet = rhophi_emcal + rhophi_hcalin + rhophi_hcalout;

        rho_jet/=jet_total_e;
        rho_hcalin/=energy_hcalin;
        rho_hcalout/=energy_hcalout;
        rho_emcal/=energy_emcal;
        rhoeta_jet/=jet_total_e;
        rhoeta_hcalin/=energy_hcalin;
        rhoeta_hcalout/=energy_hcalout;
        rhoeta_emcal/=energy_emcal;
        rhophi_jet/=jet_total_e;
        rhophi_hcalin/=energy_hcalin;
        rhophi_hcalout/=energy_hcalout;
        rhophi_emcal/=energy_emcal;

        b_jet_e_emcal.at(cone_index)->push_back(eFrac_emcal);
        b_jet_e_hcalin.at(cone_index)->push_back(eFrac_ihcal);
        b_jet_e_hcalout.at(cone_index)->push_back(eFrac_ohcal);

        b_jet_rho.at(cone_index)->push_back(rho_jet);
        b_jet_rho_emcal.at(cone_index)->push_back(rho_emcal);
        b_jet_rho_hcalin.at(cone_index)->push_back(rho_hcalin);
        b_jet_rho_hcalout.at(cone_index)->push_back(rho_hcalout);
        b_jet_rhoeta.at(cone_index)->push_back(rhoeta_jet);
        b_jet_rhoeta_emcal.at(cone_index)->push_back(rhoeta_emcal);
        b_jet_rhoeta_hcalin.at(cone_index)->push_back(rhoeta_hcalin);
        b_jet_rhoeta_hcalout.at(cone_index)->push_back(rhoeta_hcalout);
        b_jet_rhophi.at(cone_index)->push_back(rhophi_jet);
        b_jet_rhophi_emcal.at(cone_index)->push_back(rhophi_emcal);
        b_jet_rhophi_hcalin.at(cone_index)->push_back(rhophi_hcalin);
        b_jet_rhophi_hcalout.at(cone_index)->push_back(rhophi_hcalout);

      }



      for (int im = 0; im < 2; im++)
      {
        if (jet->get_pt() > max_secondmax_pt[im])
        {
          if (im == 0) max_secondmax_pt[1] = max_secondmax_pt[0];
          max_secondmax_pt[im] = jet->get_pt();
          if (im == 0) max_secondmax_t[1] = max_secondmax_t[0];
          max_secondmax_t[im] = jet_t_avg;
          if (im == 0) max_secondmax_phi[1] = max_secondmax_phi[0];
          max_secondmax_phi[im] = jet->get_phi();
          break;
        }
      }
    }
  }

  dijet_candidate = m_allevents;

  if (!m_is_auau && !m_dijet_cut && max_secondmax_pt[0] > max_pt_cone) dijet_candidate = true;

  if (!m_is_auau && m_dijet_cut)
  {
    bool ptcutpass = (max_secondmax_pt[0] > max_pt_cone && max_secondmax_pt[1] > subleading_pt_cut);
    float dphi = max_secondmax_phi[0] - max_secondmax_phi[1];
    if (dphi < -1*TMath::Pi())
    {
      dphi += 2*TMath::Pi();
    }
    if (dphi > TMath::Pi())
    {
      dphi -= 2*TMath::Pi();
    }

    double jetdeltatime = 17.6*(max_secondmax_t[1] - max_secondmax_t[0]);
    double jetleadtime = 17.6*(max_secondmax_t[0]);

    bool passleadtime = ( TMath::Abs(jetleadtime + 2.0) < 6.0 );
    bool passdijettime = (TMath::Abs(jetdeltatime) < 3.0);	  

    bool passbothtime = (passdijettime) && (passleadtime);
    if (!m_do_timing_cut) passbothtime = true;

    dijet_candidate = ptcutpass && (fabs(dphi) >= dphicut);// && (passbothtime);

    if (!isSim && cone_size == 4 && !passbothtime) m_bad_event = true;

    if (!dijet_candidate)
    {
      b_jet_pt.at(cone_index)->clear();//(jet->get_pt());
      if (!isbkg && m_calibrate_jets){ 
        b_jet_pt_calib.at(cone_index)->clear();
        b_jet_pt_smear_reco.at(cone_index)->clear();
        b_jet_pt_smear_high_reco.at(cone_index)->clear();
        b_jet_pt_smear_low_reco.at(cone_index)->clear();
        b_jet_pt_smear_truth.at(cone_index)->clear();
        b_jet_pt_smear_high_truth.at(cone_index)->clear();
        b_jet_pt_smear_low_truth.at(cone_index)->clear();
      }
      b_jet_eta.at(cone_index)->clear();//(jet->get_eta());
      b_jet_eta_det.at(cone_index)->clear();//(det_eta);
      b_jet_phi.at(cone_index)->clear();//(jet->get_phi());
      b_jet_e.at(cone_index)->clear();//(jet_total_e);
      b_jet_t.at(cone_index)->clear();//(jet_t_avg);
      b_jet_et.at(cone_index)->clear();//(jet_total_eT);
      b_jet_emcal.at(cone_index)->clear();//(eTFrac_emcal);
      b_jet_hcalin.at(cone_index)->clear();//(eTFrac_ihcal);
      b_jet_hcalout.at(cone_index)->clear();//(eTFrac_ohcal);

      if (m_deep_jet)
      {

        b_jet_e_emcal.at(cone_index)->clear();//(eFrac_emcal);
        b_jet_e_hcalin.at(cone_index)->clear();//(eFrac_ihcal);
        b_jet_e_hcalout.at(cone_index)->clear();//(eFrac_ohcal);

        b_jet_rho.at(cone_index)->clear();//(rho_jet);
        b_jet_rho_emcal.at(cone_index)->clear();//(rho_emcal);
        b_jet_rho_hcalin.at(cone_index)->clear();//(rho_hcalin);
        b_jet_rho_hcalout.at(cone_index)->clear();//(rho_hcalout);
        b_jet_rhoeta.at(cone_index)->clear();//(rhoeta_jet);
        b_jet_rhoeta_emcal.at(cone_index)->clear();//(rhoeta_emcal);
        b_jet_rhoeta_hcalin.at(cone_index)->clear();//(rhoeta_hcalin);
        b_jet_rhoeta_hcalout.at(cone_index)->clear();//(rhoeta_hcalout);
        b_jet_rhophi.at(cone_index)->clear();//(rhophi_jet);
        b_jet_rhophi_emcal.at(cone_index)->clear();//(rhophi_emcal);
        b_jet_rhophi_hcalin.at(cone_index)->clear();//(rhophi_hcalin);
        b_jet_rhophi_hcalout.at(cone_index)->clear();//(rhophi_hcalout);
      }
    }
  }
  if (m_is_auau && max_secondmax_pt[0] > max_pt_cut && isbkg) dijet_candidate = true;

  for (int im = 0; im < 2; im++)
  {
    max_secondmax_pt[im]=0;
    max_secondmax_phi[im]=0;
  }

  if (m_skim) keepLeadingJets(cone_index);

  return (dijet_candidate ? 1 : 0);

}

void DijetTreeMaker::GetNodes(PHCompositeNode* topNode)
{


}

int DijetTreeMaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "DijetTreeMaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  }


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetTreeMaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
  {
    std::cout << "DijetTreeMaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetTreeMaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "DijetTreeMaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  std::cout<<"Total events: "<<_i_event<<std::endl;

  for (int i = 0; i < 64; i++)
  {
    b_scaled_scalers[i] = (last_scaled[i] - first_scaled[i]);
    b_live_scalers[i] = (last_live[i] - first_live[i]);
    b_raw_scalers[i] = (last_raw[i] - first_raw[i]);
  }

  if (!isSim)
  {
    _tree_end->Fill();
  }

  if (m_is_auau && false)
  {
    h_centrality->Write();

    for (int i = 0 ; i < 100; i++)
    {
      h_jet_spectra[i]->Write();
      h_jet_spectra_etacut[i]->Write();
    }
  }
  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int DijetTreeMaker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "DijetTreeMaker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

Double_t DijetTreeMaker::getDPHI(Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;

  //3.141592653589
  if ( dphi > TMath::Pi() )
    dphi = dphi - 2. * TMath::Pi();
  if ( dphi <= -TMath::Pi() )
    dphi = dphi + 2. * TMath::Pi();

  if ( TMath::Abs(dphi) > TMath::Pi() ) {
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than TMath::Pi() " << std::endl;
    std::cout << " " << phi1 << ", " << phi2 << ", " << dphi << std::endl;
  }

  return dphi;

}

Double_t DijetTreeMaker::smear_pt(float pt_calib, int ijet, const std::vector<std::pair<int,float>> &matched_jets, int sign) {
  bool matched = false;
  float truthval = 0;
  for (size_t i = 0; i < matched_jets.size(); i++) {
    if (matched_jets[i].first == ijet) {
      truthval = matched_jets[i].second;
      matched = true;
      break;
    }
  }
  // Unmatched jets fall back to pt_calib as the resolution-lookup reference
  // instead of leaving them with zero smearing - ported from
  // gammajet/treemaking/src/CaloAna.cc::smear_pt(), see its comment there.
  float pt_ref = matched ? truthval : pt_calib;
  const TH1D *h_width = (sign > 0) ? h_jer_smear_up : (sign < 0) ? h_jer_smear_down : h_jer_smear_nominal;
  float width = h_width->Interpolate(pt_ref);
  return rand.Gaus(pt_calib, pt_ref*width);
}
Double_t DijetTreeMaker::DeltaR(float x1, float y1, float x2, float y2) {
  float deta = std::abs(x1-x2);
  float dphi = std::abs(y1-y2);
  if (dphi > M_PI) dphi -= 2*M_PI;

  double dr = TMath::Sqrt(deta*deta + dphi*dphi);
  return dr;
}
