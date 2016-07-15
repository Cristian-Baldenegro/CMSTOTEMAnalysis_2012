
//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

//OUR OWN CLASSES TO READ THE TREE
#include "MassParticles.h"
#include "MyBaseJet.h"
#include "MyBeamSpot.h"
#include "MyCaloJet.h"
#include "MyCastorDigi.h"
#include "MyCastorJet.h"
#include "MyCastorRecHit.h"
#include "MyDiJet.h"
#include "MyElectron.h"
#include "MyEvtId.h"
#include "MyFwdGap.h"
#include "MyGenJet.h"
#include "MyGenKin.h"
#include "MyGenMet.h"
#include "MyGenPart.h"
#include "MyHLTrig.h"
#include "MyJet.h"
#include "MyL1Trig.h"
#include "MyL1TrigOld.h"
//#include "MyMITEvtSel.h"
#include "MyMet.h"
#include "MyMuon.h"
#include "MyPFCand.h"
#include "MyPFJet.h"
#include "MyPUSumInfo.h"
#include "MyPart.h"
#include "MySimVertex.h"
#include "MyTracks.h"
#include "MyVertex.h"
#include "MyFSCHit.h"
#include "MyFSCDigi.h"

// TOTEM data formats
#include "T1Event.h"
#include "T2Event.h"
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpReconstructedProtonPair.h"
#include "RPRootDumpTrackInfo.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPatternInfo.h"

#include "/storage/lhuertas/uerj-1/CMSTOTEM/data/analysis_tools.h"
#include "/storage/lhuertas/uerj-storage1/antoniov/Workspace/Analysis/CMSSW_5_3_4_patch2/src/CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "/storage/lhuertas/uerj-storage1/antoniov/Workspace/Analysis/CMSSW_5_3_4_patch2/src/CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//STANDARD C++ INCLUDES
#include <iostream>
//#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#define PI 3.141592653589793
using namespace std;

Double_t fFermiLike(Double_t *x, Double_t *par) {
  Double_t result = 0;
  result = 1.0/(TMath::Exp((par[0]-TMath::Sqrt(x[0]))/par[1]) + 1);
  return result;
}


void data_SDdijet_jets_systematics_CMSTOTEM_ntuple(bool jet = true, bool sys_up = false,  const Int_t nevt_max = -1){
//void data_SDdijet_CMSTOTEM(string const& outputFileName = "/storage/lhuertas/uerj-1/CMSTOTEM/data/root_files/data_SDdijet_CMSTOTEM.root", const Int_t nevt_max = -1){

  TString obj = (jet) ? "jet": "pf";
  TString sys = (sys_up) ? "up": "dw";
  TString file_name = "data_SDdijet_CMSTOTEM_" + obj + "_" + sys + ".root";
  TString outputFileName = "/storage/lhuertas/uerj-1/CMSTOTEM/data/root_files/" + file_name;
  cout<<outputFileName<<endl;
  
  bool verbose = false;
  string treeName = "cms_totem";
  string jetCollName = "ak5PFJets";
  string jetCorrName = "ak5PFL2L3Residual";
  double ptJetMin = 15.0;
  double etaJetMax = 4.0;
  double etaMaxThreshold = 2.0;
  
  bool selectBunchCrossing = false;
  bool selectVertex = true;
  bool selectJets = false;
  bool selectEtaMax = false;
  bool selectEtaMin = false;
  bool selectZeroHitsT2Plus = false;
  bool selectZeroHitsT2Minus = false;
  bool selectSingleArmRecProton = true;
  bool selectDoubleArmRecProton = false;
  bool selectElastic = false;
  bool selectNonElastic = false;

  ThresholdsPerRegion thresholdsPFlow;
  thresholdsPFlow[Barrel] = ThresholdsPerType(); 
  thresholdsPFlow[Endcap] = ThresholdsPerType(); 
  thresholdsPFlow[Transition] = ThresholdsPerType(); 
  thresholdsPFlow[Endcap] = ThresholdsPerType(); 
  resetPFThresholds(thresholdsPFlow[Barrel]);
  resetPFThresholds(thresholdsPFlow[Endcap]);
  resetPFThresholds(thresholdsPFlow[Transition]);
  resetPFThresholds(thresholdsPFlow[Forward]);

  thresholdsPFlow[Barrel][MyPFCand::h0]            = make_pair(-1.,1.4);
  thresholdsPFlow[Barrel][MyPFCand::gamma]         = make_pair(-1.,0.9);
  thresholdsPFlow[Endcap][MyPFCand::h0]            = make_pair(-1.,2.7);
  thresholdsPFlow[Endcap][MyPFCand::gamma]         = make_pair(-1.,2.5);
  thresholdsPFlow[Transition][MyPFCand::h0]        = make_pair(-1.,3.8);
  thresholdsPFlow[Transition][MyPFCand::gamma]     = make_pair(-1.,2.5);
  thresholdsPFlow[Transition][MyPFCand::h_HF]      = make_pair(-1.,4.0);
  thresholdsPFlow[Transition][MyPFCand::egamma_HF] = make_pair(-1.,3.5);
  thresholdsPFlow[Forward][MyPFCand::h_HF]         = make_pair(-1.,4.0);
  thresholdsPFlow[Forward][MyPFCand::egamma_HF]    = make_pair(-1.,3.5);

  ThresholdsPerType::const_iterator pfThreshold = thresholdsPFlow[Barrel].begin();
  ThresholdsPerType::const_iterator pfThresholds_end = thresholdsPFlow[Barrel].end(); 
  ostringstream oss;
  oss << "Using the following PF thresholds:\n";
  for(; pfThreshold != pfThresholds_end; ++pfThreshold){
     int key = pfThreshold->first;    
     oss << "  " << key << ": "
                 << "(" << thresholdsPFlow[Barrel][key].first
                 << "," << thresholdsPFlow[Barrel][key].second << ")  "
                 << "(" << thresholdsPFlow[Endcap][key].first
                 << "," << thresholdsPFlow[Endcap][key].second << ")  "
                 << "(" << thresholdsPFlow[Transition][key].first
                 << "," << thresholdsPFlow[Transition][key].second << ")  "
                 << "(" << thresholdsPFlow[Forward][key].first
                 << "," << thresholdsPFlow[Forward][key].second << ")\n";   
  }
  cout << oss.str();
  //==============================


  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;

  vector<string> hltPathNames;
  hltPathNames.push_back("HLT_L1DoubleEG3_FwdVeto_v1");
  hltPathNames.push_back("HLT_L1DoubleMu0_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20_RomanPotsOR_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20part1_v1");
  hltPathNames.push_back("HLT_L1DoubleJet24_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20part2_v1");
  hltPathNames.push_back("HLT_L1Tech40_BPTXAND_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_1_v1");
  hltPathNames.push_back("HLT_L1Tech_HF9OR10_v1");
  hltPathNames.push_back("HLT_T1minbias_Tech55_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_2_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_3_v1");
  hltPathNames.push_back("HLT_RomanPots_Tech52_v1");
  hltPathNames.push_back("HLT_L1Tech54_ZeroBias_v1");
  hltPathNames.push_back("HLT_ZeroBias_v7");

  // Declaration of histograms
  map<string,TH1F*> histosTH1F;
  
  vector<string> selections;
  selections.push_back("All");
  selections.push_back("BunchCrossing");
  selections.push_back("HLT");
  selections.push_back("Vertex");
  selections.push_back("Jet");
  selections.push_back("EtaMax");
  selections.push_back("EtaMin");
  selections.push_back("ZeroHitsT2Plus");
  selections.push_back("ZeroHitsT2Minus");
  selections.push_back("SingleArmRP");
  selections.push_back("DoubleArmRP");
  selections.push_back("Elastic");
  selections.push_back("NonElastic");


  double xi_cms_minus, xi_cms_plus, x_right, x_left;
  double xi_totem_right, xi_totem_left, t_totem_right, t_totem_left, xi_cms_minus_totem, xi_cms_plus_totem, xi_totem_xicmscut_right, xi_totem_xicmscut_left;
  double beta_proton_right, beta_proton_left, phi_proton_right, phi_proton_left, thetax_proton_right, thetax_proton_left, thetay_proton_right, thetay_proton_left;
  double eff_trigger, jet1_pt, jet1_eta, jet1_phi, jet2_pt, jet2_eta, jet2_phi;
  bool  valid_proton_right, valid_proton_left, rp_right_top, rp_right_bottom, rp_left_top, rp_left_bottom;
  double vtx_x, vtx_y, vtx_z, x_pos_024, x_pos_124, x_pos_025, x_pos_125, y_pos_025, y_pos_125, y_pos_024, y_pos_124;
  TTree* small_tree = new TTree("small_tree","");
  small_tree->Branch("xi_cms_minus",&xi_cms_minus,"xi_cms_minus/D");
  small_tree->Branch("xi_cms_plus",&xi_cms_plus,"xi_cms_plus/D");
  small_tree->Branch("x_right",&x_right,"x_right/D");
  small_tree->Branch("x_left",&x_left,"x_left/D");
  small_tree->Branch("xi_totem_right",&xi_totem_right,"xi_totem_right/D");
  small_tree->Branch("xi_totem_left",&xi_totem_left,"xi_totem_left/D");
  small_tree->Branch("t_totem_right",&t_totem_right,"t_totem_right/D");
  small_tree->Branch("t_totem_left",&t_totem_left,"t_totem_left/D");
  small_tree->Branch("xi_cms_minus_totem",&xi_cms_minus_totem,"xi_cms_minus_totem/D");
  small_tree->Branch("xi_cms_plus_totem",&xi_cms_plus_totem,"xi_cms_plus_totem/D");
  small_tree->Branch("xi_totem_xicmscut_right",&xi_totem_xicmscut_right,"xi_totem_xicmscut_right/D");
  small_tree->Branch("xi_totem_xicmscut_left",&xi_totem_xicmscut_left,"xi_totem_xicmscut_left/D");
  small_tree->Branch("valid_proton_right",&valid_proton_right,"valid_proton_right/O");
  small_tree->Branch("valid_proton_left",&valid_proton_left,"valid_proton_left/O");
  small_tree->Branch("rp_right_top",&rp_right_top,"rp_right_top/O");
  small_tree->Branch("rp_right_bottom",&rp_right_bottom,"rp_right_bottom/O");
  small_tree->Branch("rp_left_top",&rp_left_top,"rp_left_top/O");
  small_tree->Branch("rp_left_bottom",&rp_left_bottom,"rp_left_bottom/O");
  small_tree->Branch("beta_proton_right",&beta_proton_right,"beta_proton_right/D");
  small_tree->Branch("beta_proton_left",&beta_proton_left,"beta_proton_left/D");
  small_tree->Branch("phi_proton_right",&phi_proton_right,"phi_proton_right/D");
  small_tree->Branch("phi_proton_left",&phi_proton_left,"phi_proton_left/D");
  small_tree->Branch("thetax_proton_right",&thetax_proton_right,"thetax_proton_right/D");
  small_tree->Branch("thetax_proton_left",&thetax_proton_left,"thetax_proton_left/D");
  small_tree->Branch("thetay_proton_right",&thetay_proton_right,"thetay_proton_right/D");
  small_tree->Branch("thetay_proton_left",&thetay_proton_left,"thetay_proton_left/D");
  small_tree->Branch("eff_trigger",&eff_trigger,"eff_trigger/D");
  small_tree->Branch("jet1_pt",&jet1_pt,"jet1_pt/D");
  small_tree->Branch("jet2_pt",&jet2_pt,"jet2_pt/D");
  small_tree->Branch("jet1_eta",&jet1_eta,"jet1_eta/D");
  small_tree->Branch("jet2_eta",&jet2_eta,"jet2_eta/D");
  small_tree->Branch("jet1_phi",&jet1_phi,"jet1_phi/D");
  small_tree->Branch("jet2_phi",&jet2_phi,"jet2_phi/D");
  small_tree->Branch("vtx_x",&vtx_x,"vtx_x/D");
  small_tree->Branch("vtx_y",&vtx_y,"vtx_y/D");
  small_tree->Branch("vtx_z",&vtx_z,"vtx_z/D");
  small_tree->Branch("x_pos_024",&x_pos_024,"x_pos_024/D");
  small_tree->Branch("x_pos_025",&x_pos_025,"x_pos_025/D");
  small_tree->Branch("x_pos_124",&x_pos_124,"x_pos_124/D");
  small_tree->Branch("x_pos_125",&x_pos_125,"x_pos_125/D");
  small_tree->Branch("y_pos_024",&y_pos_024,"y_pos_024/D");
  small_tree->Branch("y_pos_124",&y_pos_124,"y_pos_124/D");
  small_tree->Branch("y_pos_025",&y_pos_025,"y_pos_025/D");
  small_tree->Branch("y_pos_125",&y_pos_125,"y_pos_125/D");

  small_tree->SetDirectory(0);


TH1F* jet1_hist = new TH1F("jet1_hist","",100,0,300);
TH1F* jet1_up_hist = new TH1F("jet1_up_hist","",100,0,300);
TH1F* jet1_dw_hist = new TH1F("jet1_dw_hist","",100,0,300);
  //===================
  //===================
  // Instantiate JES uncertainty sources
  //const int nsrc = 19;//16;
  //const char* srcnames[nsrc] =
  //{"Absolute", "HighPtExtra", /*"SinglePion",*/ "SinglePionECAL"/*new*/, "SinglePionHCAL"/*new*/,
  // "Flavor", "Time",
  // "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
  // "RelativePtEC1"/*new*/, "RelativePtEC2"/*new*/, "RelativePtHF"/*new*/,
  // "RelativeStatEC2", "RelativeStatHF", /*"RelativeFSR",*/ /*"RelativeSample",*/ /*new*/
  // "PileUpDataMC", /*"PileUpOOT",*/ "PileUpBias", /*"PileUpJetRate"*/
  // /*"PileUpPt",*/ "PileUpPtBB"/*new*/, "PileUpPtEC"/*new*/, "PileUpPtHF"/*new*/};

  const int n_jec_unc_src = 12;
  const char* jec_unc_srcnames[n_jec_unc_src] = {"Absolute", "HighPtExtra", "Flavor", "Time",
                                         "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
                                         "RelativePtEC1"/*new*/, "RelativePtEC2"/*new*/, "RelativePtHF"/*new*/,
                                         "RelativeStatEC2", "RelativeStatHF"};
  std::vector<JetCorrectionUncertainty*> jec_unc_src(n_jec_unc_src);

  for (int isrc = 0; isrc < n_jec_unc_src; isrc++) {

     const char *name = jec_unc_srcnames[isrc];
     JetCorrectorParameters *p = new JetCorrectorParameters("jes_uncertainty/txt/Fall12_V7_DATA_UncertaintySources_AK5PF.txt", name);
     JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
     jec_unc_src[isrc] = unc;
  } // for isrc

  // Total uncertainty for reference
  JetCorrectionUncertainty *jec_unc_total = new JetCorrectionUncertainty(*(new JetCorrectorParameters("jes_uncertainty/txt/Fall12_V7_DATA_UncertaintySources_AK5PF.txt", "Total")));
  //===================
  //===================

  //===================
  int i_tot = 0 , nevt_tot = 0;


   const char *ext=".root";
 
   vector<TString>* vdirs = new vector<TString>; 
   vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198902-8369_8371-V00-02-00/Jets1/");
   vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198902-8369_8371-V00-02-00/Jets2/");
   vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/Jets1/");
   vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/Jets2/");
   
   vector<TString>* vfiles = new vector<TString>;
   for(vector<TString>::iterator itdirs = vdirs->begin(); itdirs != vdirs->end(); ++itdirs){
      TString& dirname = *itdirs;
       //vector<TString>* vfiles = new vector<TString>; 
      TSystemDirectory dir(dirname, dirname);
      TList *files = dir.GetListOfFiles();
      if (files) {
         TSystemFile *file;
         TString fname;
         TIter next(files);
         while ((file=(TSystemFile*)next())) {
             fname = file->GetName();
             if (!file->IsDirectory() && fname.EndsWith(ext)) {
                 TString root_file = dirname + string(fname.Data());
	         vfiles->push_back(root_file); cout<<root_file<<endl;      
             }
         }   
      } 
   }

  //std::ofstream ofs ("tracks.txt");
  //ofs << "  event  " << "  eta   " << "  phi  "<<endl;
 
  //Declaration of tree and its branches variables
//   TTree* tree = new TTree(treeName.c_str(),"");
  TTree* tree = NULL;
  MyEvtId*           evtId        = NULL;
  MyL1TrigOld*       l1Trig       = NULL;  
  MyHLTrig*          hltTrig      = NULL;
  //vector<MyGenPart>* genPart      = NULL;
  vector<MyTracks>*  track_coll   = NULL;
  vector<MyVertex>*  vertex_coll  = NULL;
  vector<MyPFJet>*   pfJet_coll   = NULL;
  vector<MyPFCand>*  pFlow_coll   = NULL;
  vector<MyFSCHit>*  fscHits_coll = NULL;
  vector<MyFSCDigi>* fscDigis_coll = NULL;
  //===================
  T2Event* t2_event = NULL;
  RPRootDumpReconstructedProton* rec_proton_left  = NULL;
  RPRootDumpReconstructedProton* rec_proton_right = NULL;
  RPRootDumpReconstructedProtonPair* rec_proton_pair  = NULL;
  map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
  map<unsigned int, RPRootDumpDigiInfo*> rp_digi_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_par_patterns_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_nonpar_patterns_info;
  map<unsigned int, std::vector<RPRootDumpTrackInfo>*> rp_multi_track_info;
  //===================  
  int n_evt_total = 0; 
  int n_evt_trigger = 0; 
  int n_evt_vtx = 0; 
  int n_evt_jets = 0; 
  int n_evt_PF = 0; 
  int n_evt_etapf = 0; 
  int n_evt_pass_threshold = 0; 
  int n_evt_proton_right = 0; 
  int n_evt_proton_right_cut = 0; 
  int n_evt_proton_left = 0; 
  int n_evt_proton_left_cut = 0; 

  //starting Loop over files, stops at end of list of files or when reached nevt_max
  for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){
  
    TFile* file = TFile::Open(*itfiles,"READ");
    
    //getting the tree form the current file
    tree = (TTree*) file->Get( treeName.c_str() );

    //ofstream ofs;
    //ofs.Open ("eff_ptJet2.txt");

    //Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

    //adding branches to the tree ----------------------------------------------------------------------
    tree->SetBranchAddress("cmsEvtUA",&evtId);
    tree->SetBranchAddress("cmsTrigUA",&l1Trig);
    tree->SetBranchAddress("cmsHLTTrigUA",&hltTrig);
    tree->SetBranchAddress("cmsTracksUA",&track_coll);
    tree->SetBranchAddress("cmsVerticesUA",&vertex_coll);
//    tree->SetBranchAddress("cmsPFJetsUA",&pfJet_coll);
     tree->SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll);
    tree->SetBranchAddress("cmsParticleFlowUA",&pFlow_coll);
    //tree->SetBranchAddress("cmsFSCHitsUA",&fscHits_coll);
    //tree->SetBranchAddress("cmsFSCDigisUA",&fscDigis_coll);
    tree->SetBranchAddress("branchT2EV.",&t2_event);
    tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
    tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
    tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);
//     tree->SetBranchAddress("Evt",&evtId);
//     tree->SetBranchAddress("L1TrigOld",&l1Trig);
//     tree->SetBranchAddress("HLTrig",&hltTrig);
//     tree->SetBranchAddress("generalTracks",&track_coll); 
//     tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
//     tree->SetBranchAddress("ak5PFJets",&pfJet_coll);
//     tree->SetBranchAddress("particleFlow",&pFlow_coll);
    //if(isMC) tree->SetBranchAddress("genPart",&genPart);
    std::vector<unsigned int> rp_list;
    rp_list.push_back(20); rp_list.push_back(21); rp_list.push_back(22); rp_list.push_back(23); rp_list.push_back(24); rp_list.push_back(25);
    rp_list.push_back(120); rp_list.push_back(121); rp_list.push_back(122); rp_list.push_back(123); rp_list.push_back(124); rp_list.push_back(125);
    char br_name[200];
    for (unsigned int a = 0; a < 2; ++a) {
       int s = 2;
       for (unsigned int r = 0; r < 6; r++) {
          unsigned int id = 100 * a + 10 * s + r;
          if( std::find(rp_list.begin(), rp_list.end(), id) == rp_list.end() ) continue;

          sprintf(br_name, "track_rp_%u.", id);
          std::cout << br_name << std::endl;
          tree->SetBranchAddress(br_name, &rp_track_info[id]);
       }
    } 
    
    /*//Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;*/
 
    int n_evt_Jet = 0;
    //double eff_sum = 0; 

 
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){

      //printing the % of events done every 10k evts
      if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i_evt);

      ++n_evt_total;
      bool passedHLT_jets = false;
      bool passedvtx = false;
      bool jet1_selected = false;
      bool jet2_selected = false;
      bool PF_eta_max = false;
      bool PF_eta_min = false;
      bool passed_threshold = false;
      
      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
      double event_weight = 1.;
      //double event_weight_eff;
      //double event_weight_averagept_eff;

      map<string,bool>::iterator it_hlt = (*hltTrig).HLTmap.begin();
      map<string,bool>::iterator it_hlt_end = (*hltTrig).HLTmap.end();
      for(; it_hlt != it_hlt_end; ++it_hlt){
         string const& hltName = it_hlt->first;
         vector<string>::const_iterator it_pos = find(hltPathNames.begin(),hltPathNames.end(),hltName);
         if(it_pos != hltPathNames.end()){
            size_t idx = it_pos - hltPathNames.begin();//cout <<hltName<<endl;
            if( hltName == "HLT_L1DoubleJet20part1_v1" || hltName == "HLT_L1DoubleJet20part2_v1"){
	      if( it_hlt->second == true ){
	         passedHLT_jets = true; 
//		 histosTH1F["hltTrigFired"]->Fill( idx, event_weight );
	      }
	    }
	 }   
	         /*for(int ibin = 1; ibin <= histosTH1F["hltTrigFired"]->GetNbinsX(); ++ibin){
            if( hltName.c_str() != histosTH1F["hltTrigFired"]->GetXaxis()->GetBinLabel(ibin) ) continue;
            
            if( it_hlt->second ) 
               histosTH1F["hltTrigFired"]->Fill( histosTH1F["hltTrigFired"]->GetBinCenter( ibin ) );
         }*/ 
      }
      if(!passedHLT_jets) continue;
      ++n_evt_trigger;      
      //-------------------------------------------------------------------------------------------------
      //filling pt distribution for the generated particles
      //ie those from pythia generator, without reconstruction
      /*if(isMC){
        for(vector<MyGenPart>::iterator p=genPart->begin() ; p!=genPart->end() ; p++ )
          pt_gen->Fill(p->Pt());
      }*/
      
      //-------------------------------------------------------------------------------------------------
 	    
     // Vertices
      /*for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
        if (it_vtx!=vertex_coll->begin()) continue;
         int idx_vtx = it_vtx - vertex_coll->begin();
         if( it_vtx->fake ) continue;
         if( !it_vtx->validity ) continue;
         if( it_vtx->ndof<4 ) continue;
      }
*/  
      //if(!passedvtx) continue;
      MyVertex& primaryVertex = vertex_coll->at(0); 
      double prim_vtx_r = sqrt( primaryVertex.x*primaryVertex.x + primaryVertex.y*primaryVertex.y );
      bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity && primaryVertex.ndof > 4);// && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);
      if (!select_Vertex) continue;
      ++n_evt_vtx;
      
      // Tracks
/*      int n_tracks_selected = 0;
      for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
         histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
         histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
         histosTH1F["track_phi"]->Fill( it_trk->Phi(), event_weight );
         if( it_trk->Pt() < 0.5 ) continue;
         if( fabs( it_trk->Eta() ) > 2.5 ) continue;
         if( ( it_trk->dz / it_trk->edz ) > 5. ) continue;
         if( ( it_trk->d0 / it_trk->ed0 ) > 5. ) continue;
   
         if( !it_trk->quality[2] ) continue;

         ++n_tracks_selected;
//        ofs << i_evt << it_trk->Eta() << it_trk->Phi() << endl;
      }
      histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );
*/
      
 
     //Jets with pt>30Gev and !eta!<2
      Double_t Jet1_pt; 
      Double_t Jet1_pz, Jet1_px, Jet1_py; 
      Double_t Jet1_E; 
      Double_t Jet2_pt, Jet2_px, Jet2_py; 
      Double_t Jet2_E; 
      Double_t Jet2_pz; 
      Double_t Jet1_eta; 
      Double_t Jet2_eta; 
      Double_t Jet1_phi; 
      Double_t Jet2_phi;
      //Double_t eff; 
      //Double_t averagept_eff; 

      ///Fit 
     // TF1* func = new TF1("func", fFermiLike, 0., 20., 2);
     // func->SetParameter(0,5.525);
      //func->SetParameter(1,0.529);      

       
      for(vector<MyPFJet>::iterator it_jet = pfJet_coll->begin() ; it_jet != pfJet_coll->end() ; ++it_jet){
         map<string,MyBaseJet>::iterator it_map = it_jet->mapjet.begin();
         for(; it_map != it_jet->mapjet.end(); ++it_map)
            if(verbose) cout << it_map->first << endl;

         MyBaseJet const& basejet = it_jet->mapjet[jetCorrName];
     //    histosTH1F["jet_pt"]->Fill( basejet.Pt(), event_weight  );
        // if(basejet.Pt() > 0.) histosTH1F["jet_eta"]->Fill( basejet.Eta(), event_weight  );
       //  histosTH1F["jet_phi"]->Fill( basejet.Phi(), event_weight  );
      }
      if( pfJet_coll->size() > 0 ){
	 MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[jetCorrName];
	 Jet1_E = leadingJet.E(); 
	 Jet1_px = leadingJet.Px(); 
	 Jet1_py = leadingJet.Py(); 
	 Jet1_pz = leadingJet.Pz(); 
	 Jet1_pt = leadingJet.Pt(); 
	 Jet1_eta = leadingJet.Eta(); 
	 Jet1_phi = leadingJet.Phi(); 


         double sum2_up(0), sum2_dw(0);
         for (int isrc = 0; isrc < n_jec_unc_src; isrc++) {
            //cout << isrc << endl; 
            JetCorrectionUncertainty *unc = jec_unc_src[isrc];
            unc->setJetPt(Jet1_pt);
            unc->setJetEta(Jet1_eta);
            double sup = unc->getUncertainty(true); // up variation
            //cout << sup << endl;
            unc->setJetPt(Jet1_pt);
            unc->setJetEta(Jet1_eta);
            double sdw = unc->getUncertainty(false); // down variation
            //cout << sdw << endl;

            sum2_up += pow(max(sup,sdw),2);
            sum2_dw += pow(min(sup,sdw),2);
         } // for isrc
         jet1_hist->Fill(Jet1_pt, 1.);
         jet1_up_hist->Fill(Jet1_pt*(1+sqrt(sum2_up)), 1.);
         jet1_dw_hist->Fill(Jet1_pt*(1-sqrt(sum2_dw)), 1.);

         //cout << sqrt(sum2_up) << endl;
         //cout << sqrt(sum2_dw) << endl;
         jec_unc_total->setJetPt(Jet1_pt);
         jec_unc_total->setJetEta(Jet1_eta);
         double uncert_total = jec_unc_total->getUncertainty(true);
         //cout << uncert_total << endl;

         if(verbose) std::cout << "Leading jet pT, uncertainties total, individual sources up/dw: "
                               << Jet1_pt << ", " << uncert_total
                               << ", " << sqrt(sum2_up) << ", " << sqrt(sum2_dw) << std::endl;
         Jet1_pt = (sys_up == true) ? Jet1_pt*(1. + sqrt(sum2_up)) : Jet1_pt*(1. - sqrt(sum2_dw));
	 
	 if(Jet1_pt > 30. && fabs(Jet1_eta)<4.4 ) jet1_selected = true;
	 
      }
      //if(!jet1_selected) continue;
      
      if( pfJet_coll->size() > 1 ){
	 MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet[jetCorrName];
         Jet2_E = secondJet.E(); 
         Jet2_px = secondJet.Px(); 
         Jet2_py = secondJet.Py(); 
         Jet2_pz = secondJet.Pz(); 
         Jet2_pt = secondJet.Pt(); 
	 Jet2_eta = secondJet.Eta(); 
	 Jet2_phi = secondJet.Phi(); 

         double sum2_up_jet2(0), sum2_dw_jet2(0);
         for (int isrc = 0; isrc < n_jec_unc_src; isrc++) {
            //cout << isrc << endl; 
            JetCorrectionUncertainty *unc = jec_unc_src[isrc];
            unc->setJetPt(Jet2_pt);
            unc->setJetEta(Jet2_eta);
            double sup = unc->getUncertainty(true); // up variation
            //cout << sup << endl;
            unc->setJetPt(Jet2_pt);
            unc->setJetEta(Jet2_eta);
            double sdw = unc->getUncertainty(false); // down variation
            //cout << sdw << endl;

            sum2_up_jet2 += pow(max(sup,sdw),2);
            sum2_dw_jet2 += pow(min(sup,sdw),2);
         } // for isrc
         //cout << sqrt(sum2_up) << endl;
         //cout << sqrt(sum2_dw) << endl;
         jec_unc_total->setJetPt(Jet2_pt);
         jec_unc_total->setJetEta(Jet2_eta);
         double uncert_total = jec_unc_total->getUncertainty(true);
         //cout << uncert_total << endl;
         if(verbose) std::cout << "Second jet pT, uncertainties total, individual sources up/dw: "
                               << Jet2_pt << ", " << uncert_total
                               << ", " << sqrt(sum2_up_jet2) << ", " << sqrt(sum2_dw_jet2) << std::endl;

         Jet2_pt = (sys_up == true) ? Jet2_pt*(1. + sqrt(sum2_up_jet2)) : Jet2_pt*(1. - sqrt(sum2_dw_jet2));

	 if(Jet2_pt > 30. && fabs(Jet2_eta)<4.4 )  jet2_selected = true;
	 
      }
      //if(!jet2_selected) continue;
      
      //double average_pt = (Jet1_pt+Jet2_pt)/2;
      //eff = 1/func->Eval(Jet2_pt);
      //eff_sum += func->Eval(Jet2_pt);
      //averagept_eff = 1/func->Eval(average_pt);
    
      //double protons_correction = (1.067+1.066+1.060+1.063)/4; 

      //event_weight_eff = protons_correction/eff;
      //event_weight_averagept_eff = 1/averagept_eff;
      if(jet1_selected && jet2_selected) ++n_evt_jets;

      double x_plus = ((Jet1_E+Jet1_pz)+(Jet2_E+Jet2_pz))/8000;
      double x_minus = ((Jet1_E-Jet1_pz)+(Jet2_E-Jet2_pz))/8000;
      double x_minus_sel = (x_minus<x_plus) ? x_minus : x_plus;
      double x_plus_sel = (x_minus>x_plus) ? x_plus : x_minus;
      double mass_jets= sqrt(pow(Jet1_E+Jet2_E,2)-pow(Jet1_px+Jet2_px,2)-pow(Jet1_py+Jet2_py,2)-pow(Jet1_pz+Jet2_pz,2));

 
      // Particle-flow
      double soma1 = 0;
      double soma2 = 0;
      double eta_max=-999.;
      double eta_min=999;
      double cm = 8000;
      
      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
	 int partType = it_pfcand->particleId;
	 double eta = it_pfcand->Eta();
	 double energy = it_pfcand->Energy();
	 double pz = it_pfcand->Pz();
	 
         // Apply thresholds
         if( !pflowThreshold(*it_pfcand,thresholdsPFlow) ) continue;

	     soma1 += (energy + pz);
	     soma2 += (energy - pz);

             if (eta > eta_max) {eta_max = eta; PF_eta_max = true;} 
	     if (eta < eta_min) {eta_min = eta; PF_eta_min = true;}

      //   }  
      }
      //if(!PF_eta_max) continue;
      //if(!PF_eta_min) continue;
      double xi_plus_Reco = soma1/cm;
      double xi_minus_Reco = soma2/cm;
      double delta_eta_maxmin = eta_max - eta_min;
      double correction = 1/0.8;


      //if(jet1_selected && jet2_selected && PF_eta_max && PF_eta_min && passed_threshold) ++n_evt_pass_threshold;


      // TOTEM T2
      vector<double> const& t2_trk_entryX = t2_event->TrkEntryX;
      vector<double> const& t2_trk_entryY = t2_event->TrkEntryY;
      vector<double> const& t2_trk_entryZ =  t2_event->TrkEntryZ;
      vector<double> const& t2_trk_chiProb =  t2_event->TrkChiProb;

      int n_t2_tracks_selected = 0;
      int n_t2_tracks_selected_zplus = 0;
      int n_t2_tracks_selected_zminus = 0;
      size_t n_t2_tracks = t2_trk_chiProb.size();
      for(size_t i_t2_trk = 0; i_t2_trk < n_t2_tracks; ++i_t2_trk){
         double trk_entryZ = t2_trk_entryZ[i_t2_trk];
         int zside = ( trk_entryZ >= 0. ) ? 1 : -1;
         // Select tracks
         if( t2_trk_chiProb[i_t2_trk] < 0.2 ) continue;

         ++n_t2_tracks_selected;
         if( zside > 0 ) ++n_t2_tracks_selected_zplus;
         else            ++n_t2_tracks_selected_zminus;

      }
      if( selectZeroHitsT2Plus && (n_t2_tracks_selected_zplus > 0) ) continue;
      //histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Plus", event_weight );

      if( selectZeroHitsT2Minus && (n_t2_tracks_selected_zminus > 0) ) continue;
      //histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Minus", event_weight );

      bool proton_right_valid = rec_proton_right->valid;
      bool proton_left_valid = rec_proton_left->valid;
      if( selectSingleArmRecProton && (proton_right_valid && proton_left_valid) ) continue;
     // histosTH1F["EventSelection"]->Fill( "SingleArmRP", event_weight );

      if( selectDoubleArmRecProton && !(proton_right_valid && proton_left_valid) ) continue;
      //histosTH1F["EventSelection"]->Fill( "DoubleArmRP", event_weight );

      bool tag_elastic_top45_bot56 = elastic_top45_bot56(rp_track_info);      
      bool tag_elastic_bot45_top56 = elastic_bot45_top56(rp_track_info);      
      if( selectElastic && !(tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
      //histosTH1F["EventSelection"]->Fill( "Elastic", event_weight );

      if( selectNonElastic && (tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
      //histosTH1F["EventSelection"]->Fill( "NonElastic", event_weight );

 

      //fiducial cuts
      bool rp_track_valid_120 = rp_track_info[120]->valid; 
      bool rp_track_valid_121 = rp_track_info[121]->valid;
      bool rp_track_valid_122 = rp_track_info[122]->valid;
      bool rp_track_valid_123 = rp_track_info[123]->valid; 
      bool rp_track_valid_124 = rp_track_info[124]->valid;
      bool rp_track_valid_125 = rp_track_info[125]->valid;
      bool rp_track_valid_020 = rp_track_info[20]->valid;
      bool rp_track_valid_021 = rp_track_info[21]->valid;
      bool rp_track_valid_022 = rp_track_info[22]->valid;
      bool rp_track_valid_023 = rp_track_info[23]->valid;
      bool rp_track_valid_024 = rp_track_info[24]->valid;
      bool rp_track_valid_025 = rp_track_info[25]->valid;
      double rp_x_024 = rp_track_info[24]->x;
      double rp_y_024 = rp_track_info[24]->y;
      double rp_x_025 = rp_track_info[25]->x;
      double rp_y_025 = rp_track_info[25]->y;
      double rp_x_124 = rp_track_info[124]->x;
      double rp_y_124 = rp_track_info[124]->y;
      double rp_x_125 = rp_track_info[125]->x;
      double rp_y_125 = rp_track_info[125]->y;

      bool cut_rp_024 =  rp_x_024>0 && rp_x_024<6 && rp_y_024>8.4 && rp_y_024<29 ;
      bool cut_rp_025 =  rp_x_025>0 && rp_x_025<6 && rp_y_025<-8.4 && rp_y_025>-29 ;
      bool cut_rp_124 =  rp_x_124>0 && rp_x_124<6 && rp_y_124>8.4 && rp_y_124<27 ;
      bool cut_rp_125 =  rp_x_125>0 && rp_x_125<6 && rp_y_125<-8.4 && rp_y_125>-27 ;

      int rp_hits_120 = rp_track_info[120]->entries;
      int rp_hits_121 = rp_track_info[121]->entries;
      int rp_hits_122 = rp_track_info[122]->entries;
      int rp_hits_123 = rp_track_info[123]->entries;
      int rp_hits_124 = rp_track_info[124]->entries;
      int rp_hits_125 = rp_track_info[125]->entries;

      bool rp_track_right_top = rp_track_valid_120 && rp_track_valid_124 && cut_rp_124;
      bool rp_track_right_bottom = rp_track_valid_121 && rp_track_valid_125 && cut_rp_125;
      bool rp_track_accept_right = rp_track_right_top || rp_track_right_bottom;// || ( rp_track_valid_122 && rp_track_valid_123 );
      bool rp_track_left_top = rp_track_valid_020 && rp_track_valid_024 && cut_rp_024;
      bool rp_track_left_bottom = rp_track_valid_021 && rp_track_valid_025 && cut_rp_025;
      bool rp_track_accept_left = rp_track_left_top  || rp_track_left_bottom;// || ( rp_track_valid_022 && rp_track_valid_023 );

      

      // RP protons
      double xi_totem;
      double chi2_proton_right = rec_proton_right->chi2;
      double chindf_proton_right = rec_proton_right->chindf;
      double xi_proton_right = rec_proton_right->xi;//cout<<xi_proton_right<<endl;
      double t_proton_right = rec_proton_right->t;
      bool xi_region_right = -xi_proton_right>=0.03 && -xi_proton_right<=0.1;
      bool t_region_right = fabs(t_proton_right)>=0.03 && fabs(t_proton_right)<=1;
      bool good_proton_right = proton_right_valid; // && (xi_proton_right < 0.);
      double chi2_proton_left = rec_proton_left->chi2;
      double chindf_proton_left = rec_proton_left->chindf;
      double xi_proton_left = rec_proton_left->xi;//cout<<xi_proton_left<<endl;
      double t_proton_left = rec_proton_left->t;
      bool xi_region_left = -xi_proton_left>=0.03 && -xi_proton_left<=0.1;
      bool t_region_left = -t_proton_left>=0.03 && -t_proton_left<=1;
      bool good_proton_left = proton_left_valid;// && (xi_proton_left < 0.);
      double mass_x = sqrt(xi_proton_right*8000);
      double xi_cms_totem_minus_cut = xi_minus_Reco+xi_proton_right<0; 
      double xi_cms_totem_plus_cut = xi_plus_Reco+xi_proton_left<0; 

      xi_cms_minus = xi_minus_Reco;
      xi_cms_plus = xi_plus_Reco;
      jet1_pt = Jet1_pt;
      jet2_pt = Jet2_pt;
      jet1_eta = Jet1_eta;
      jet2_eta = Jet2_eta;
      jet1_phi = Jet1_phi;
      jet2_phi = Jet2_phi;

      valid_proton_right = proton_right_valid;
      valid_proton_left = proton_left_valid;
      rp_right_top = rp_track_right_top;
      rp_right_bottom = rp_track_right_bottom;
      rp_left_top = rp_track_left_top;
      rp_left_bottom = rp_track_left_bottom;
      xi_totem_right = -xi_proton_right; 
      t_totem_right = -t_proton_right;
      xi_cms_minus_totem = xi_minus_Reco+xi_proton_right;
      xi_totem_left = -xi_proton_left;
      t_totem_left = -t_proton_left;
      xi_cms_plus_totem = xi_plus_Reco+xi_proton_left;
      beta_proton_right = x_minus/-xi_proton_right;
      beta_proton_left = x_plus/-xi_proton_left;
      thetax_proton_right = rec_proton_right->thx;
      thetay_proton_right = rec_proton_right->thy;
      phi_proton_right = rec_proton_right->phi;
      thetax_proton_left = rec_proton_left->thx;
      thetay_proton_left = rec_proton_left->thy;
      phi_proton_left = rec_proton_left->phi;
      x_right = x_minus;
      x_left = x_plus;
      vtx_x = primaryVertex.x;
      vtx_y = primaryVertex.y;
      vtx_z = primaryVertex.z;
      x_pos_024 = rp_x_024;
      x_pos_025 = rp_x_025;
      y_pos_024 = rp_y_024;
      y_pos_025 = rp_y_025;
      x_pos_124 = rp_x_124;
      x_pos_125 = rp_x_125;
      y_pos_124 = rp_y_124;
      y_pos_125 = rp_y_125;

      small_tree->Fill();

    }//end loop for events
    file->Close();
//ofs.close();
  
  }//end of loop over files

  //output file
  TFile* output = new TFile(outputFileName,"RECREATE");
  output->cd();


  small_tree->Write();
  jet1_hist->Write();
  jet1_up_hist->Write();
  jet1_dw_hist->Write();
 
  output->Close();

}
