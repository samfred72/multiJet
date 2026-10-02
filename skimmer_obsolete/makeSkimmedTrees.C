#include <iostream>

using std::cout;
using std::endl;

bool default1 = false;

struct jet
{
  int id;
  float pt = 0;
  float t = 0;
  float eta = 0;
  float eta_det = 0;
  float phi = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;

  void print()
  {
    std::cout << "Jet " << id << std::endl;
    std::cout << "    pt/eta/phi : " << pt << " / " << eta << " / " << phi << std::endl;
  };

};

const float truth_cut = 3;
const float calib_cut = 7;
const float leading_cut = 20;
const float dphicut = 3*TMath::Pi()/4.;
const float dphicutloose = TMath::Pi()/2.;
const float vertex_cut = 60;

float getDPHI(float phi1, float phi2);
bool check_dijet_reco(std::vector<struct jet> myrecojets);
bool check_dijet_truth(std::vector<struct jet> mytruthjets);
double getDR(struct jet j1, struct jet j2)
{ 
  double dphi = fabs(j1.phi - j2.phi);
  if (dphi > TMath::Pi())
  {
    dphi = 2*TMath::Pi() - dphi;
  }

  double dR = sqrt(TMath::Power(j1.eta - j2.eta, 2) + TMath::Power(dphi, 2));
  return dR;
}
void makeSkimmedTrees(int runn = 0, const char * sample = "Jet12", const char * simname = "pythia")
{
  bool isMC = (runn == 0);

  std::string infile =   Form("/sphenix/tg/tg01/jets/samfred/multiJet_hadded/multijet_%s_%s.root", simname, sample);
  const char * outfile = Form("/sphenix/tg/tg01/jets/samfred/multiJet_skimmed/TREE_MULTIJET_SKIM_%s_%s.root", simname, sample);
  if (!isMC) {
    infile = Form("/sphenix/tg/tg01/jets/samfred/multiJet_hadded/multijet_%i.root", runn);
    outfile = Form("/sphenix/tg/tg01/jets/samfred/multiJet_skimmed/TREE_MULTIJET_SKIM_%i.root", runn);
  }

  TFile *f = new TFile(infile.c_str(),"r");

  TTree *ttree = (TTree *) f->Get("ttree");

  int njetR = 7;
  ULong64_t gl1_scaled;
  std::vector<std::vector<float>*> truth_jet_pt(njetR);
  std::vector<std::vector<float>*> truth_jet_eta(njetR);
  std::vector<std::vector<float>*> truth_jet_phi(njetR);

  std::vector<std::vector<float>*> reco_jet_pt(njetR);
  std::vector<std::vector<float>*> reco_jet_pt_calib(njetR);
  std::vector<std::vector<float>*> reco_jet_pt_smear(njetR);
  std::vector<std::vector<float>*> reco_jet_pt_smearRECO(njetR);
  std::vector<std::vector<float>*> reco_jet_pt_smearHIGH(njetR);
  std::vector<std::vector<float>*> reco_jet_pt_smearLOW(njetR);
  std::vector<std::vector<float>*> reco_jet_emcal(njetR);
  std::vector<std::vector<float>*> reco_jet_e(njetR);
  std::vector<std::vector<float>*> reco_jet_eta(njetR);
  std::vector<std::vector<float>*> reco_jet_eta_det(njetR);
  std::vector<std::vector<float>*> reco_jet_phi(njetR);
  std::vector<std::vector<float>*> reco_jet_t(njetR);
  float truth_vertex_z;

  double calib_lead_time = 0;
  double calib_delta_time = 0;
  double calib_mbd_time = 0;
  int runnumber;

  float time_zero;  
  float mbd_vertex_z;
  int mbd_hit;
  double mbd_maxdelta;
  double mbd_avgdelta;
  double mbd_proddelta;
  double mbd_maxsigma;
  double mbd_avgsigma;
  double mbd_prodsigma;

  cout << "Setting up ttrees..." << endl;
  for (int i = 0; i < njetR; i++) {
    int conesize = i + 2;
    if (isMC) {
      ttree->SetBranchAddress(Form("truth_jet_pt_%d", conesize), &truth_jet_pt.at(i));
      ttree->SetBranchAddress(Form("truth_jet_eta_%d", conesize), &truth_jet_eta.at(i));
      ttree->SetBranchAddress(Form("truth_jet_phi_%d", conesize), &truth_jet_phi.at(i));
      ttree->SetBranchAddress(Form("jet_pt_smear_%d", conesize), &reco_jet_pt_smear.at(i));
      ttree->SetBranchAddress(Form("jet_pt_smearRECO_%d", conesize), &reco_jet_pt_smearRECO.at(i));
      ttree->SetBranchAddress(Form("jet_pt_smearHIGH_%d", conesize), &reco_jet_pt_smearHIGH.at(i));
      ttree->SetBranchAddress(Form("jet_pt_smearLOW_%d", conesize), &reco_jet_pt_smearLOW.at(i));
    }
    else {
      ttree->SetBranchAddress(Form("jet_t_%d", conesize), &reco_jet_t.at(i));
    }
    ttree->SetBranchAddress(Form("jet_pt_%d", conesize), &reco_jet_pt.at(i));
    ttree->SetBranchAddress(Form("jet_pt_calib_%d", conesize), &reco_jet_pt_calib.at(i));
    ttree->SetBranchAddress(Form("jet_emcal_%d", conesize), &reco_jet_emcal.at(i));
    ttree->SetBranchAddress(Form("jet_e_%d", conesize), &reco_jet_e.at(i));
    ttree->SetBranchAddress(Form("jet_eta_%d", conesize), &reco_jet_eta.at(i));
    ttree->SetBranchAddress(Form("jet_eta_det_%d", conesize), &reco_jet_eta_det.at(i));
    ttree->SetBranchAddress(Form("jet_phi_%d", conesize), &reco_jet_phi.at(i));
  }
  if (isMC) {
    ttree->SetBranchAddress("truth_vertex_z", &truth_vertex_z);
  }
  else {
    ttree->SetBranchAddress("runnumber", &runnumber);
    ttree->SetBranchAddress("calib_lead_time", &calib_lead_time);
    ttree->SetBranchAddress("calib_delta_time", &calib_delta_time);
    ttree->SetBranchAddress("calib_mbd_time", &calib_mbd_time);
  }
  ttree->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);
  ttree->SetBranchAddress("mbd_hit", &mbd_hit);
  ttree->SetBranchAddress("time_zero", &time_zero);
  ttree->SetBranchAddress("gl1_scaled", &gl1_scaled);
  ttree->SetBranchAddress("mbd_avgsigma", &mbd_avgsigma);
  ttree->SetBranchAddress("mbd_prodsigma", &mbd_prodsigma);
  ttree->SetBranchAddress("mbd_maxsigma", &mbd_maxsigma);
  ttree->SetBranchAddress("mbd_avgdelta", &mbd_avgdelta);
  ttree->SetBranchAddress("mbd_proddelta", &mbd_proddelta);
  ttree->SetBranchAddress("mbd_maxdelta", &mbd_maxdelta);

  int nevents = ttree->GetEntries();

  TFile *fout = TFile::Open(outfile, "RECREATE");
  TTree *tout = new TTree("ttree", "skimmed tree");

  for (int i = 0; i < njetR; i++) {
    int conesize = i + 2;
    if (isMC) {
      tout->Branch(Form("truth_jet_pt_%d", conesize), &truth_jet_pt.at(i));
      tout->Branch(Form("truth_jet_eta_%d", conesize), &truth_jet_eta.at(i));
      tout->Branch(Form("truth_jet_phi_%d", conesize), &truth_jet_phi.at(i));
      tout->Branch(Form("jet_pt_smear_%d", conesize), &reco_jet_pt_smear.at(i));
      tout->Branch(Form("jet_pt_smearRECO_%d", conesize), &reco_jet_pt_smearRECO.at(i));
      tout->Branch(Form("jet_pt_smearHIGH_%d", conesize), &reco_jet_pt_smearHIGH.at(i));
      tout->Branch(Form("jet_pt_smearLOW_%d", conesize), &reco_jet_pt_smearLOW.at(i));
    }
    else {
      tout->Branch(Form("jet_t_%d", conesize), &reco_jet_t.at(i));
    }
    tout->Branch(Form("jet_pt_%d", conesize), &reco_jet_pt.at(i));
    tout->Branch(Form("jet_pt_calib_%d", conesize), &reco_jet_pt_calib.at(i));
    tout->Branch(Form("jet_emcal_%d", conesize), &reco_jet_emcal.at(i));
    tout->Branch(Form("jet_e_%d", conesize), &reco_jet_e.at(i));
    tout->Branch(Form("jet_eta_%d", conesize), &reco_jet_eta.at(i));
    tout->Branch(Form("jet_eta_det_%d", conesize), &reco_jet_eta_det.at(i));
    tout->Branch(Form("jet_phi_%d", conesize), &reco_jet_phi.at(i));
  }
  if (isMC) {
    tout->Branch("truth_vertex_z", &truth_vertex_z);
  }
  else {
    tout->Branch("runnumber", &runnumber);
    tout->Branch("calib_lead_time", &calib_lead_time);
    tout->Branch("calib_delta_time", &calib_delta_time);
    tout->Branch("calib_mbd_time", &calib_mbd_time);
  }

  tout->Branch("mbd_vertex_z", &mbd_vertex_z);
  tout->Branch("mbd_avgsigma", &mbd_avgsigma);
  tout->Branch("mbd_prodsigma", &mbd_prodsigma);
  tout->Branch("mbd_maxsigma", &mbd_maxsigma);
  tout->Branch("mbd_avgdelta", &mbd_avgdelta);
  tout->Branch("mbd_proddelta", &mbd_proddelta);
  tout->Branch("mbd_maxdelta", &mbd_maxdelta);
  tout->Branch("mbd_hit", &mbd_hit);
  tout->Branch("time_zero", &time_zero);
  tout->Branch("gl1_scaled", &gl1_scaled);

  std::vector<std::vector<struct jet>> mytruthjets(njetR);
  std::vector<std::vector<struct jet>> myrecojets(njetR);
  std::vector<std::vector<struct jet>> myrecojets_smear(njetR);
  std::vector<std::vector<struct jet>> myrecojets_smearRECO(njetR);
  std::vector<std::vector<struct jet>> myrecojets_smearHIGH(njetR);
  std::vector<std::vector<struct jet>> myrecojets_smearLOW(njetR);


  cout << "Now run..." << endl;
  for (int i = 0; i < nevents; i++)
  {
    ttree->GetEntry(i);
    if (i % 10000 == 0) std::cout << "Event " << i << " \r "<< std::flush;
    if (std::isnan(mbd_vertex_z) || fabs(mbd_vertex_z) > vertex_cut) continue;
    //if (!mbd_hit) continue;
    
    // check emulated triggertrigger bits
    bool passes_Trigger = ((gl1_scaled >> 22) & 0x1);
    bool anypass = false;

    for (int ir = 0; ir < njetR; ir++) {
      int ntruthjets = (isMC ? truth_jet_pt.at(ir)->size() : 0);      
      int nrecojets = reco_jet_pt.at(ir)->size();

      mytruthjets.at(ir).clear();
      myrecojets.at(ir).clear();
      myrecojets_smear.at(ir).clear();
      myrecojets_smearRECO.at(ir).clear();
      myrecojets_smearHIGH.at(ir).clear();
      myrecojets_smearLOW.at(ir).clear();

      for (int j = 0; j < ntruthjets;j++)
      {

        if (truth_jet_pt.at(ir)->at(j) < truth_cut) continue;

        struct jet tempjet;

        tempjet.pt = truth_jet_pt.at(ir)->at(j);
        tempjet.eta = truth_jet_eta.at(ir)->at(j);
        tempjet.phi = truth_jet_phi.at(ir)->at(j);
        tempjet.id = j;

        mytruthjets.at(ir).push_back(tempjet);	  	  

      }

      std::sort(mytruthjets.at(ir).begin(), mytruthjets.at(ir).end(), [] (auto a, auto b) { return a.pt > b.pt; });

      for (int j = 0; j < nrecojets;j++){
        if (reco_jet_pt_calib.at(ir)->at(j) < calib_cut && 
            (!isMC || 
             (reco_jet_pt_smear.at(ir)->at(j) < calib_cut &&
              reco_jet_pt_smearRECO.at(ir)->at(j) < calib_cut &&
              reco_jet_pt_smearHIGH.at(ir)->at(j) < calib_cut &&
              reco_jet_pt_smearLOW.at(ir)->at(j) < calib_cut
             )
            )
           ) continue;
        if (reco_jet_e.at(ir)->at(j) < 0) continue;

        struct jet tempjet;
        struct jet tempjet_smear;
        struct jet tempjet_smearRECO;
        struct jet tempjet_smearHIGH;
        struct jet tempjet_smearLOW;

        tempjet.pt = reco_jet_pt_calib.at(ir)->at(j);
        tempjet.emcal = reco_jet_emcal.at(ir)->at(j);
        tempjet.eta = reco_jet_eta.at(ir)->at(j);
        tempjet.eta_det = reco_jet_eta_det.at(ir)->at(j);
        tempjet.phi = reco_jet_phi.at(ir)->at(j);
        tempjet.id = j;
        tempjet.t = (isMC ? 0 : reco_jet_t.at(ir)->at(j));

        myrecojets.at(ir).push_back(tempjet);

        if (isMC) {
          tempjet_smear.pt = reco_jet_pt_smear.at(ir)->at(j);
          tempjet_smear.emcal = reco_jet_emcal.at(ir)->at(j);
          tempjet_smear.eta = reco_jet_eta.at(ir)->at(j);
          tempjet_smear.eta_det = reco_jet_eta_det.at(ir)->at(j);
          tempjet_smear.phi = reco_jet_phi.at(ir)->at(j);
          tempjet_smear.id = j;
          tempjet_smear.t = 0;
          
          tempjet_smearRECO.pt = reco_jet_pt_smearRECO.at(ir)->at(j);
          tempjet_smearRECO.emcal = reco_jet_emcal.at(ir)->at(j);
          tempjet_smearRECO.eta = reco_jet_eta.at(ir)->at(j);
          tempjet_smearRECO.eta_det = reco_jet_eta_det.at(ir)->at(j);
          tempjet_smearRECO.phi = reco_jet_phi.at(ir)->at(j);
          tempjet_smearRECO.id = j;
          tempjet_smearRECO.t = 0;

          tempjet_smearHIGH.pt = reco_jet_pt_smearHIGH.at(ir)->at(j);
          tempjet_smearHIGH.emcal = reco_jet_emcal.at(ir)->at(j);
          tempjet_smearHIGH.eta = reco_jet_eta.at(ir)->at(j);
          tempjet_smearHIGH.eta_det = reco_jet_eta_det.at(ir)->at(j);
          tempjet_smearHIGH.phi = reco_jet_phi.at(ir)->at(j);
          tempjet_smearHIGH.id = j;
          tempjet_smearHIGH.t = 0;

          tempjet_smearLOW.pt = reco_jet_pt_smearLOW.at(ir)->at(j);
          tempjet_smearLOW.emcal = reco_jet_emcal.at(ir)->at(j);
          tempjet_smearLOW.eta = reco_jet_eta.at(ir)->at(j);
          tempjet_smearLOW.eta_det = reco_jet_eta_det.at(ir)->at(j);
          tempjet_smearLOW.phi = reco_jet_phi.at(ir)->at(j);
          tempjet_smearLOW.id = j;
          tempjet_smearLOW.t = 0;

          myrecojets_smear.at(ir).push_back(tempjet_smear);
          myrecojets_smearRECO.at(ir).push_back(tempjet_smearRECO);
          myrecojets_smearHIGH.at(ir).push_back(tempjet_smearHIGH);
          myrecojets_smearLOW.at(ir).push_back(tempjet_smearLOW);
        }
      }


      std::sort(  myrecojets.at(ir).begin(),           myrecojets.at(ir).end(),           [] (auto a, auto b) { return a.pt > b.pt; });
      if (isMC) {
        std::sort(myrecojets_smear.at(ir).begin(),     myrecojets_smear.at(ir).end(),     [] (auto a, auto b) { return a.pt > b.pt; });
        std::sort(myrecojets_smearRECO.at(ir).begin(), myrecojets_smearRECO.at(ir).end(), [] (auto a, auto b) { return a.pt > b.pt; });
        std::sort(myrecojets_smearHIGH.at(ir).begin(), myrecojets_smearHIGH.at(ir).end(), [] (auto a, auto b) { return a.pt > b.pt; });
        std::sort(myrecojets_smearLOW.at(ir).begin(),  myrecojets_smearLOW.at(ir).end(),  [] (auto a, auto b) { return a.pt > b.pt; });
      }

      if (check_dijet_reco(myrecojets.at(ir)) ||
          (isMC && 
           (
            check_dijet_reco(myrecojets_smear.at(ir)) ||
            check_dijet_reco(myrecojets_smearRECO.at(ir)) ||
            check_dijet_reco(myrecojets_smearHIGH.at(ir)) ||
            check_dijet_reco(myrecojets_smearLOW.at(ir)) //||
            //check_dijet_truth(mytruthjets))
          )
         )
        ) {
          anypass = true;
          break;
        }
    }
    if (anypass) tout->Fill();
  }
  std::cout << "Writing file to " << fout->GetName() << std::endl;
  fout->Write();
  fout->Close();
}

bool check_dijet_reco(std::vector<struct jet> myrecojets)
{

  if (myrecojets.size() < 3) return false;

  auto leading_iter = myrecojets.begin();

  auto subleading_iter = myrecojets.begin() + 1;
  auto subsubleading_iter = myrecojets.begin() + 2;

  float dphir = getDPHI(leading_iter->phi, subleading_iter->phi);
  float dphirr = getDPHI(leading_iter->phi, subsubleading_iter->phi);

  if (!(leading_iter->pt >= leading_cut && 
        subleading_iter->pt >= calib_cut && 
        subsubleading_iter->pt >= calib_cut && 
        dphir >= dphicut && 
        dphirr >= dphicutloose
        )
      ) return false;

  double jetdeltatime = 17.6*(leading_iter->t - subleading_iter->t);
  double jetleadtime = 17.6*(leading_iter->t);
  bool passleadtime = ( TMath::Abs(jetleadtime +2.0) < 6.0 );
  bool passdijettime = (TMath::Abs(jetdeltatime) < 3.0);	  

  bool passbothtime = (passdijettime) && (passleadtime);

  if (!passbothtime) return false;

  return true;

}
/*
bool check_dijet_truth(std::vector<struct jet> mytruthjets)
{

  if (mytruthjets.size() < 2) return false;

  auto leading_iter = mytruthjets.begin();

  auto subleading_iter = mytruthjets.begin() + 1;

  float dphir = getDPHI(leading_iter->phi, subleading_iter->phi);

  if (!(leading_iter->pt >= 10 && subleading_iter->pt >= 5 && dphir >= dphicut)) return false;

  return true;


}
*/

float getDPHI(float phi1, float phi2)
{
  float dphi = std::fabs(phi1 - phi2);
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return dphi;
}
