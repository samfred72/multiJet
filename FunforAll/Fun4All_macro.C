#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <fstream>
#include <filesystem>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>

#include <G4_Global.C>
#include <GlobalVariables.C>
#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>
#include <globalvertex/GlobalVertexReco.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloWaveformProcessing.h>
//
//#include <calotowerbuilder/CaloTowerBuilder.h>
//
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

//#include <caloana/ClusterIso.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/DeadHotMapLoader.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterDeadHotMask.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <caloreco/TowerInfoDeadHotMask.h>
#include <caloreco/PhotonClusterBuilder.h>
#include <clusteriso/ClusterIso.h>

#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetCalib.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <jetbackground/TimingCut.h>

// #include <runtowerinfo/RunTowerInfo.h>
#include <dijettreemaker/DijetTreeMaker.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexReco.h>
#include <mbd/MbdReco.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4centrality/PHG4CentralityReco.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <calotrigger/CaloTriggerEmulator.h>

#include <centrality/CentralityReco.h>
#include <calotrigger/MinimumBiasClassifier.h>

#include "HIJetReco.C"
#include <Calo_Calib.C>

#include <sstream>
#include <fstream>
#include <string>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libdijettreemaker.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4vertex.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libzdcinfo.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libclusteriso.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libFROG.so)

#endif

// Real GL1 prescale lookup - meaningless for simulation, only called when !isMC.
void get_scaledowns(int runnumber, int scaledowns[]){

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","phnxro","");
  if (db) {
    printf("Server info: %s\n", db->ServerInfo());
  }
  else {
    printf("bad\n");
  }
  TSQLRow *row;
  TSQLResult *res;
  TString cmd = "";
  for (int is = 0; is < 64; is++) {
    const char * sql = Form("select scaledown%02d from gl1_scaledown where runnumber = %d;", is, runnumber);
    res = db->Query(sql);
    int nrows = res->GetRowCount();
    int nfields = res->GetFieldCount();
    for (int i = 0; i < nrows; i++) {
      row = res->Next();
      for (int j = 0; j < nfields; j++) {
        scaledowns[is] = stoi(row->GetField(j));
      }
      delete row;
    }
    delete res;
  }
  delete db;
}

// Combined data/MC treemaking driver (formerly Fun4All_macro.C + MCFun4All_macro.C).
// isMC selects the file-handling path: data's `infile` is a flat list of one
// DST_TOWERS path per line; MC's `infile` is a list of one line per event
// chunk, each line holding 5 whitespace-separated DST paths in the fixed
// order cluster/global/truth/jet/mbd, matching the split-DST MC production.
// trigger/sim only matter for isMC (they name the output file).
void Fun4All_macro(const char* infile="/sphenix/user/samfred/projects/filelists/goldenruns_pp24_calojet_split15/queue_52403_v0011.list",
                    bool test=true,
                    bool isMC=false,
                    const char* trigger="MC",
                    const char* sim="pythia"){

    //=====================
    // Filename management
    //=====================

    string outdir = isMC ? "/sphenix/tg/tg01/jets/samfred/multiJet_MC/" : "/sphenix/tg/tg01/jets/samfred/multiJet";
    void * dirf = gSystem->OpenDirectory(outdir.c_str());
    if(dirf) gSystem->FreeDirectory(dirf);
    else {gSystem->mkdir(outdir.c_str(), kTRUE);}

    std::cout << "infile : " << infile << " isMC : " << isMC << std::endl;

    // Parsed input file paths, by DST type. Only the branch matching isMC is filled.
    std::vector<std::string> infile_towers;                                      // data: DST_TOWERS
    std::vector<std::string> infile_cluster, infile_global, infile_truth, infile_jet, infile_mbd;  // MC: split DSTs
    std::string firstfilename;

    if (isMC) {
      if (!std::filesystem::exists(infile)) {
        std::cout << "File '" << infile << "' does not exist." << std::endl;
        return;
      }
      std::ifstream file(infile);
      std::string line;
      while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string dst1, dst2, dst3, dst4, dst5;
        if (iss >> dst1 >> dst2 >> dst3 >> dst4 >> dst5) {
          infile_cluster.push_back(dst1);
          infile_global.push_back(dst2);
          infile_truth.push_back(dst3);
          infile_jet.push_back(dst4);
          infile_mbd.push_back(dst5);
        }
      }
      if (infile_global.empty()) {
        std::cerr << "Error reading first filename or file is empty!" << std::endl;
        return;
      }
      firstfilename = infile_global.at(0);
    }
    else {
      std::ifstream infilestream(infile);
      if (!infilestream.is_open()) {
        std::cerr << "Error opening file: " << infile << std::endl;
        return;
      }
      std::string fn;
      while (infilestream >> fn) {
        infile_towers.push_back(fn);
      }
      if (infile_towers.empty()) {
        std::cerr << "Error reading first filename or file is empty!" << std::endl;
        return;
      }
      firstfilename = infile_towers.front();
      std::cout << "First filename: " << firstfilename << std::endl;
    }

    std::filesystem::path p(infile);
    std::string fname = p.filename().stem().string();
    string outfile = "";
    if (test) {
      // Always write test output to the current directory - avoids
      // hardcoding a path into either project's tree.
      outfile = "test.root";
    }
    else if (isMC) {
      outfile = Form("outtree_%s_%s_%s.root", fname.c_str(), sim, trigger);
    }
    else {
      outfile = Form("outtree_%s.root", fname.c_str());
    }

    Fun4AllServer *se = Fun4AllServer::instance();
    int verbosity = 0;

    se->Verbosity(verbosity);
    recoConsts *rc = recoConsts::instance();

    pair<int, int> runseg = Fun4AllUtils::GetRunSegment(firstfilename.c_str());
    int runnumber = runseg.first;
    int segment = runseg.second;

    //=====================
    // conditions DB flags
    //=====================

    // global tag
    rc->set_StringFlag("CDB_GLOBALTAG", isMC ? "MDC2" : "ProdA_2024");
    rc->set_uint64Flag("TIMESTAMP",runnumber);

    if (!isMC) {
      // scaledowns (real GL1 prescales; not meaningful for simulation)
      int m_scaledowns[64];
      get_scaledowns(runnumber,m_scaledowns);
    }

    //=====================
    // Global reco
    //=====================
    MbdReco *mbdreco = new MbdReco();
    se->registerSubsystem(mbdreco);
    GlobalVertexReco *gvertex = new GlobalVertexReco();
    se->registerSubsystem(gvertex);

    if (isMC) {
      //=====================
      // Trigger subsystem
      //=====================
      // MC has no real GL1/trigger-primitive data, so the trigger has to be
      // emulated before DijetTreeMaker's TriggerAnalyzer (UseEmulator(true),
      // set below) can read it. Not needed for data, which has real GL1 info.
      TriggerRunInfoReco *trig = new TriggerRunInfoReco();
      trig->UseEmulator(true); // for sim
      se->registerSubsystem(trig);
      CaloTriggerEmulator *te = new CaloTriggerEmulator("CALOTRIGGEREMULATOR");
      te->setNSamples(12);
      te->setTriggerSample(6);
      te->setTriggerDelay(5);
      te->SetIsData(false); //for sim
      te->setEmcalLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/emcal_ll1_lut_0.50tr_new.root");
      te->setHcalinLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/hcalin_ll1_lut_0.50tr_new.root");
      te->setHcaloutLUTFile("/sphenix/user/dlis/Projects/macros/CDBTest/hcalout_ll1_lut_0.50tr_new.root");
      se->registerSubsystem(te);
    }

    //====================
    // Calo Calib
    //====================
    Process_Calo_Calib();

    //====================
    // Jet reco
    //====================
    std::vector<float> doUnsubJet_radius = {0.2,0.3,0.4,0.5,0.6,0.7,0.8};

    // retowering
    std::string jetreco_input_prefix = "TOWERINFO_CALIB";
    RetowerCEMC *_retowerCEMC;
    _retowerCEMC = new RetowerCEMC();
    _retowerCEMC->Verbosity(verbosity);
    _retowerCEMC->set_towerinfo(true);
    _retowerCEMC->set_frac_cut(0.5); //fraction of retower that must be masked to mask the full retower
    _retowerCEMC->set_towerNodePrefix(jetreco_input_prefix);
    se->registerSubsystem(_retowerCEMC);

    // Jet reco
    JetReco *_jetRecoUnsub = new JetReco();
    _jetRecoUnsub->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER, jetreco_input_prefix));
    _jetRecoUnsub->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO, jetreco_input_prefix));
    _jetRecoUnsub->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO, jetreco_input_prefix));
    for (int ir = 0; ir < doUnsubJet_radius.size(); ++ir) {
      _jetRecoUnsub->add_algo(new FastJetAlgoSub(Jet::ANTIKT, doUnsubJet_radius[ir]), "AntiKt_TowerInfo_r0" + std::to_string((int)(10*doUnsubJet_radius[ir])));
    }
    _jetRecoUnsub->set_algo_node("ANTIKT");
    _jetRecoUnsub->set_input_node("TOWER");
    _jetRecoUnsub->Verbosity(verbosity);
    se->registerSubsystem(_jetRecoUnsub);

    for (int ir = 0; ir < doUnsubJet_radius.size(); ++ir) {
      cout << "Input: " << Form("AntiKt_TowerInfo_r0%i",(int)(doUnsubJet_radius[ir]*10)) << endl;
      cout << "Output: " << Form("AntiKt_TowerInfo_r0%i_calib",(int)(doUnsubJet_radius[ir]*10)) << endl << endl;
      JetCalib * jetCalib = new JetCalib(Form("JetCalib0%i",(int)(doUnsubJet_radius[ir]*10)));
      jetCalib->set_InputNode(Form("AntiKt_TowerInfo_r0%i",(int)(doUnsubJet_radius[ir]*10)));
      jetCalib->set_OutputNode(Form("AntiKt_TowerInfo_r0%i_calib",(int)(doUnsubJet_radius[ir]*10)));
      jetCalib->set_JetRadius(doUnsubJet_radius[ir]);
      // set_ZvrtxNode() does not exist in the current JetCalib API (vertex
      // lookup is internal now); set_ApplyZvrtxDependentCalib/EtaDependentCalib
      // only take effect under the legacy method - see JetCalib.h and
      // gammajet/treemaking/macros/MCFun4All_macro.C's jetCalibOld block.
      jetCalib->set_UseEMfracCalib(false);
      jetCalib->set_ApplyZvrtxDependentCalib(true);
      jetCalib->set_ApplyEtaDependentCalib(true);
      se->registerSubsystem(jetCalib);
    }

    //======================
    // File inputs
    //======================

    if (isMC) {
      Fun4AllInputManager *incluster = new Fun4AllDstInputManager("DST_CLUSTER");
      for (auto &fn : infile_cluster) { incluster->AddFile(fn); cout << fn << endl; }
      se->registerInputManager(incluster);

      Fun4AllInputManager *inglobal = new Fun4AllDstInputManager("DST_GLOBAL");
      for (auto &fn : infile_global) inglobal->AddFile(fn);
      se->registerInputManager(inglobal);

      Fun4AllInputManager *injet = new Fun4AllDstInputManager("DST_JET");
      for (auto &fn : infile_jet) injet->AddFile(fn);
      se->registerInputManager(injet);

      Fun4AllInputManager *inmbd = new Fun4AllDstInputManager("DST_MBD");
      for (auto &fn : infile_mbd) inmbd->AddFile(fn);
      se->registerInputManager(inmbd);

      Fun4AllInputManager *intruth = new Fun4AllDstInputManager("DST_TRUTH");
      for (auto &fn : infile_truth) intruth->AddFile(fn);
      se->registerInputManager(intruth);
    }
    else {
      Fun4AllInputManager *intower = new Fun4AllDstInputManager("DST_TOWERS");
      for (auto &fn : infile_towers) {
        intower->AddFile(fn);
        std::cout << "input file : " << fn << std::endl;
      }
      se->registerInputManager(intower);
    }

    //======================
    // Calo Ana
    //======================

    if (!isMC) {
      // Timing cut on real HCALOUT tower timing - not meaningful for
      // simulation, so only registered for data.
      // arguments (only first required): 1. jet node name to use, 2. name
      // for fun4all module, 3. doAbort, 4. ohTowerName. The third tells the
      // module whether or not to abort events that fail the cuts.
      TimingCut* tc = new TimingCut("AntiKt_TowerInfo_r04","TimingCutModule",false,"TOWERINFO_CALIB_HCALOUT");
      se->registerSubsystem(tc);
    }

    DijetTreeMaker *tt1 = new DijetTreeMaker("DijetTreemaker",outfile.c_str());
    tt1->CalibrateJets(true);
    tt1->SetIsSim(isMC);
    tt1->UseEmulator(true);
    tt1->SaveAllEvents(false); // combined single-pass mode: selection now applied inline, see ApplySkimCuts()
    tt1->SetPtCutTruth(3);
    tt1->SetPtCut(3);
    tt1->skim(true);
    tt1->ApplySkimCuts(true); // fold makeSkimmedTrees.C selection into this one Fun4All pass
    tt1->deepJet(false);
    tt1->addRecoConeSize(2,0);
    tt1->addTruthConeSize(2,0);
    tt1->addRecoConeSize(3,0);
    tt1->addTruthConeSize(3,0);
    tt1->addRecoConeSize(4,0);
    tt1->addTruthConeSize(4,0);
    tt1->addRecoConeSize(5,0);
    tt1->addTruthConeSize(5,0);
    tt1->addRecoConeSize(6,0);
    tt1->addTruthConeSize(6,0);
    tt1->addRecoConeSize(7,0);
    tt1->addTruthConeSize(7,0);
    tt1->addRecoConeSize(8,0);
    tt1->addTruthConeSize(8,0);
    tt1->Verbosity(0);
    se->registerSubsystem(tt1);

    std::cout << "now run..." << std::endl;
    if (test) {
      se->run(isMC ? 100 : 1000);
    }
    else {
      se->run();
    }
    se->End();
    std::cout << "ok done.. " << std::endl;
    std::cout << "Written to " << outfile << std::endl;

    delete se;
}
