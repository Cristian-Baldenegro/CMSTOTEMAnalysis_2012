
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
#include <TRandom.h>

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

#include "analysis_tools.h"

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



//void zb_ntuple(string const& outputFileName = "/storage/lhuertas/uerj-1/CMSTOTEM/data/root_files/zb_ntuple.root", const Int_t nevt_max = -1){
void zb_ntuple(string const& outputFileName = "/afs/cern.ch/user/l/lhuertas/zb_ntuple.root", const Int_t nevt_max = -1){
  
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


  double xi_cms_minus_st, xi_cms_plus, x_right, x_left, beta_proton_right, beta_proton_left, jet1_pt, jet1_eta, jet1_phi, jet2_pt, jet2_eta, jet2_phi;
  double xi_totem_right_st, xi_totem_left, t_totem_right, t_totem_left, xi_cms_minus_totem, xi_cms_plus_totem, xi_totem_xicmscut_right, xi_totem_xicmscut_left;
  bool  valid_proton_right, valid_proton_left, rp_right, rp_left, valid_vtx;
  bool rp_right_top_near, rp_right_top_far, rp_right_bottom_near, rp_right_bottom_far, rp_left_top_near, rp_left_top_far, rp_left_bottom_near, rp_left_bottom_far;
  double rp_x_020, rp_y_020, rp_x_120, rp_y_120, rp_x_021, rp_y_021, rp_x_024, rp_y_024, rp_x_025, rp_y_025, rp_x_121, rp_y_121, rp_x_124, rp_y_124, rp_x_125, rp_y_125;
  TTree* small_tree = new TTree("small_tree","");
  small_tree->Branch("xi_cms_minus",&xi_cms_minus_st,"xi_cms_minus/D");
  small_tree->Branch("xi_cms_plus",&xi_cms_plus,"xi_cms_plus/D");
  small_tree->Branch("xi_totem_right",&xi_totem_right_st,"xi_totem_right/D");
  small_tree->Branch("xi_totem_left",&xi_totem_left,"xi_totem_left/D");
  small_tree->Branch("t_totem_right",&t_totem_right,"t_totem_right/D");
  small_tree->Branch("t_totem_left",&t_totem_left,"t_totem_left/D");
  small_tree->Branch("xi_cms_minus_totem",&xi_cms_minus_totem,"xi_cms_minus_totem/D");
  small_tree->Branch("xi_cms_plus_totem",&xi_cms_plus_totem,"xi_cms_plus_totem/D");
  small_tree->Branch("xi_totem_xicmscut_right",&xi_totem_xicmscut_right,"xi_totem_xicmscut_right/D");
  small_tree->Branch("xi_totem_xicmscut_left",&xi_totem_xicmscut_left,"xi_totem_xicmscut_left/D");
  small_tree->Branch("valid_proton_right",&valid_proton_right,"valid_proton_right/O");
  small_tree->Branch("valid_proton_left",&valid_proton_left,"valid_proton_left/O");
  //small_tree->Branch("rp_right",&rp_right,"rp_right/O");
  small_tree->Branch("rp_right_top_near",&rp_right_top_near,"rp_right_top_near/O");
  small_tree->Branch("rp_right_top_far",&rp_right_top_far,"rp_right_top_far/O");
  small_tree->Branch("rp_right_bottom_near",&rp_right_bottom_near,"rp_right_bottom_near/O");
  small_tree->Branch("rp_right_bottom_far",&rp_right_bottom_far,"rp_right_bottom_far/O");
  //small_tree->Branch("rp_left",&rp_left,"rp_left/O");
  small_tree->Branch("rp_left_top_near",&rp_left_top_near,"rp_left_top_near/O");
  small_tree->Branch("rp_left_top_far",&rp_left_top_far,"rp_left_top_far/O");
  small_tree->Branch("rp_left_bottom_near",&rp_left_bottom_near,"rp_left_bottom_near/O");
  small_tree->Branch("rp_left_bottom_far",&rp_left_bottom_far,"rp_left_bottom_far/O");
  small_tree->Branch("valid_vtx",&valid_vtx,"valid_vtx/O");
  small_tree->Branch("x_right",&x_right,"x_right/D");
  small_tree->Branch("x_left",&x_left,"x_left/D");
  small_tree->Branch("beta_proton_right",&beta_proton_right,"beta_proton_right/D");
  small_tree->Branch("beta_proton_left",&beta_proton_left,"beta_proton_left/D");
  small_tree->Branch("jet1_pt",&jet1_pt,"jet1_pt/D");
  small_tree->Branch("jet2_pt",&jet2_pt,"jet2_pt/D");
  small_tree->Branch("jet1_eta",&jet1_eta,"jet1_eta/D");
  small_tree->Branch("jet2_eta",&jet2_eta,"jet2_eta/D");
  small_tree->Branch("jet1_phi",&jet1_phi,"jet1_phi/D");
  small_tree->Branch("jet2_phi",&jet2_phi,"jet2_phi/D");
  small_tree->Branch("rp_x_020",&rp_x_020,"rp_x_020/D");
  small_tree->Branch("rp_y_020",&rp_y_020,"rp_y_020/D");
  small_tree->Branch("rp_x_021",&rp_x_021,"rp_x_021/D");
  small_tree->Branch("rp_y_021",&rp_y_021,"rp_y_021/D");
  small_tree->Branch("rp_x_024",&rp_x_024,"rp_x_024/D");
  small_tree->Branch("rp_y_024",&rp_y_024,"rp_y_024/D");
  small_tree->Branch("rp_x_025",&rp_x_025,"rp_x_025/D");
  small_tree->Branch("rp_y_025",&rp_y_025,"rp_y_025/D");
  /*small_tree->Branch("rp_x_120",&rp_x_120,"rp_x_120/D");
  small_tree->Branch("rp_y_120",&rp_y_120,"rp_y_120/D");
  small_tree->Branch("rp_x_121",&rp_x_121,"rp_x_121/D");
  small_tree->Branch("rp_y_121",&rp_y_121,"rp_y_121/D");
  small_tree->Branch("rp_x_124",&rp_x_124,"rp_x_124/D");
  small_tree->Branch("rp_y_124",&rp_y_124,"rp_y_124/D");
  small_tree->Branch("rp_x_125",&rp_x_125,"rp_x_125/D");
  small_tree->Branch("rp_y_125",&rp_y_125,"rp_y_125/D");
*/
  small_tree->SetDirectory(0);
TH1F* test = new TH1F("test","",50,0,0.2);
  map<string,TH2F*> histosTH2F;
  histosTH2F["proton_y_vs_x_rp_024_025"] = new TH2F("proton_y_vs_x_rp_024_025","proton_y_vs_x_rp_024_025",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_024_025_accept"] = new TH2F("proton_y_vs_x_rp_024_025_accept","proton_y_vs_x_rp_024_025",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_124_125"] = new TH2F("proton_y_vs_x_rp_124_125","proton_y_vs_x_rp_124_125",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_124_125_accept"] = new TH2F("proton_y_vs_x_rp_124_125_accept","proton_y_vs_x_rp_124_125",100,-10,10,100,-40,40);

  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
  //===================
  int i_tot = 0 , nevt_tot = 0;


   const char *ext=".root";
 
   vector<TString>* vdirs = new vector<TString>; 
   //vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198902-8369_8371-V00-02-00/ZeroBias/");
   //vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/ZeroBias/");
   //vdirs->push_back("/afs/cern.ch/user/l/lhuertas/work/uerj1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198902-8369_8371-V00-02-00/ZeroBias/");
   //vdirs->push_back("/afs/cern.ch/user/l/lhuertas/work/uerj1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/ZeroBias/");
   vdirs->push_back("root://eoscms.cern.ch//store/group/phys_diffraction/CMSTOTEM_2012/MergedNtuples/HighBeta/reReco/198902-8369_8371-V00-02-00_new/ZeroBias/");
   vdirs->push_back("root://eoscms.cern.ch//store/group/phys_diffraction/CMSTOTEM_2012/MergedNtuples/HighBeta/reReco/198903-8372-V00-02-00_new/ZeroBias/");
   
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
                 if (root_file == "/afs/cern.ch/user/l/lhuertas/work/uerj1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/ZeroBias/val8372ea_totem_ntuple_UABaseTree_CMS-TOTEM_17_1_Fwu.root") continue;
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
    if (file->IsZombie()) continue;

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


     xi_cms_minus_st = -999.; xi_cms_plus = 999.; xi_totem_right_st = 0; xi_totem_left= 999; t_totem_right = 0; t_totem_left = -999.; xi_cms_minus_totem = -999.; xi_cms_plus_totem = 999.;
     xi_totem_xicmscut_right = -999.; xi_totem_xicmscut_left= 999.;
     valid_proton_right = false; valid_proton_left = false; rp_right = false; rp_left = false;

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

      
      //-------------------------------------------------------------------------------------------------
 	    

      MyVertex& primaryVertex = vertex_coll->at(0); 
      valid_vtx = ( !primaryVertex.fake && primaryVertex.validity && primaryVertex.ndof > 4);// && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);
      
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
      for(vector<MyPFJet>::iterator it_jet = pfJet_coll->begin() ; it_jet != pfJet_coll->end() ; ++it_jet){
         map<string,MyBaseJet>::iterator it_map = it_jet->mapjet.begin();
         for(; it_map != it_jet->mapjet.end(); ++it_map)
            if(verbose) cout << it_map->first << endl;

         MyBaseJet const& basejet = it_jet->mapjet[jetCorrName];
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
	 
	 if(Jet1_pt > 30. && Jet1_eta<4.4 ) jet1_selected = true;
	 
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

	 if(Jet2_pt > 30. && Jet2_eta<4.4 )  jet2_selected = true;
	 
      }
      //if(!jet2_selected) continue;
      
      if(jet1_selected && jet2_selected) ++n_evt_jets;

      double x_plus = ((Jet1_E+Jet1_pz)+(Jet2_E+Jet2_pz))/8000;
      double x_minus = ((Jet1_E-Jet1_pz)+(Jet2_E-Jet2_pz))/8000;
 
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
         //if( t2_trk_chiProb[i_t2_trk] < 0.2 ) continue;

         ++n_t2_tracks_selected;
         if( zside > 0 ) ++n_t2_tracks_selected_zplus;
         else            ++n_t2_tracks_selected_zminus;

      }
      if( selectZeroHitsT2Plus && (n_t2_tracks_selected_zplus > 0) ) continue;

      if( selectZeroHitsT2Minus && (n_t2_tracks_selected_zminus > 0) ) continue;


      bool proton_right_valid = rec_proton_right->valid;
      bool proton_left_valid = rec_proton_left->valid;

      if( proton_right_valid && proton_left_valid ) ++n_events_DoubleArm;

      bool tag_elastic_top45_bot56 = elastic_top45_bot56(rp_track_info);      
      bool tag_elastic_bot45_top56 = elastic_bot45_top56(rp_track_info);      
      //if( tag_elastic_top45_bot56 || tag_elastic_bot45_top56  )  ++n_events_Elastic;

      // Select single-arm events (inclusive)
      if( selectSingleArmRecProton && !(proton_right_valid || proton_left_valid) ) continue;

      if( selectDoubleArmRecProton && !(proton_right_valid && proton_left_valid) ) continue;
      if( selectElastic && !(tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
      // Veto elastic-tagged events
      if( selectNonElastic && (tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
 

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
      rp_x_020 = rp_track_info[20]->x;
      rp_y_020 = rp_track_info[20]->y;
      rp_x_021 = rp_track_info[21]->x;
      rp_y_021 = rp_track_info[21]->y;
      rp_x_120 = rp_track_info[120]->x;
      rp_y_120 = rp_track_info[120]->y;
      rp_x_121 = rp_track_info[121]->x;
      rp_y_121 = rp_track_info[121]->y;
      rp_x_024 = rp_track_info[24]->x;
      rp_y_024 = rp_track_info[24]->y;
      rp_x_025 = rp_track_info[25]->x;
      rp_y_025 = rp_track_info[25]->y;
      rp_x_124 = rp_track_info[124]->x;
      rp_y_124 = rp_track_info[124]->y;
      rp_x_125 = rp_track_info[125]->x;
      rp_y_125 = rp_track_info[125]->y;

      histosTH2F["proton_y_vs_x_rp_024_025"]->Fill( rp_x_024, rp_y_024, event_weight );
      histosTH2F["proton_y_vs_x_rp_024_025"]->Fill( rp_x_025, rp_y_025, event_weight );
      histosTH2F["proton_y_vs_x_rp_124_125"]->Fill( rp_x_124, rp_y_124, event_weight );
      histosTH2F["proton_y_vs_x_rp_124_125"]->Fill( rp_x_125, rp_y_125, event_weight );

      bool cut_rp_024 =  rp_x_024>0 && rp_x_024<6 && rp_y_024>8.4 && rp_y_024<29 ;
      bool cut_rp_025 =  rp_x_025>0 && rp_x_025<6 && rp_y_025<-8.4 && rp_y_024>-29 ;
      bool cut_rp_124 =  rp_x_124>0 && rp_x_124<6 && rp_y_124>8.4 && rp_y_124<27 ;
      bool cut_rp_125 =  rp_x_125>0 && rp_x_125<6 && rp_y_125<-8.4 && rp_y_124>-27 ;

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
      double xi_proton_right = rec_proton_right->xi;//cout<<xi_proton_right<<endl;
      double t_proton_right = rec_proton_right->t;
      bool proton_right_valid = rec_proton_right->valid; // && (xi_proton_right < 0.);
      bool proton_left_valid = rec_proton_left->valid; // && (xi_proton_left < 0.);
      double xi_proton_left = rec_proton_left->xi;//cout<<xi_proton_left<<endl;
      double t_proton_left = rec_proton_left->t;
      bool good_proton_left = proton_left_valid;// && (xi_proton_left < 0.);
    
      test->Fill(-xi_proton_right); 

      jet1_pt = Jet1_pt;
      jet2_pt = Jet2_pt;
      jet1_eta = Jet1_eta;
      jet2_eta = Jet2_eta;
      jet1_phi = Jet1_phi;
      jet2_phi = Jet2_phi;

      valid_proton_right = proton_right_valid;
      valid_proton_left = proton_left_valid;
      rp_right = rp_track_accept_right;
      rp_left = rp_track_accept_left;
      rp_right_top_near = rp_track_valid_120; 
      rp_right_top_far = rp_track_valid_124; 
      rp_right_bottom_near = rp_track_valid_121; 
      rp_right_bottom_far = rp_track_valid_125; 
      rp_left_top_near = rp_track_valid_020;
      rp_left_top_far = rp_track_valid_024; 
      rp_left_bottom_near = rp_track_valid_021;
      rp_left_bottom_far = rp_track_valid_025;

      //if(good_proton_right && rp_track_accept_right){
         xi_totem_right_st = -xi_proton_right; 
         t_totem_right = -t_proton_right;
        // if (xi_minus_Reco>0.12) xi_totem_xicmscut_right = -xi_proton_right; 
      //}



      //if (good_proton_left && rp_track_accept_left){
          xi_totem_left = -xi_proton_left;
         t_totem_left = -t_proton_left;
        // if (xi_plus_Reco>0.12) xi_totem_xicmscut_left = -xi_proton_left; 
      //}

      x_right = x_minus;
      x_left = x_plus;
      beta_proton_right = x_minus/-xi_proton_right;
      beta_proton_left = x_plus/-xi_proton_left;

    small_tree->Fill(); 


    }//end loop for events
    file->Close();
//ofs.close();
  
  }//end of loop over files

  //output file
  TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  output->cd();



  small_tree->Write();

test->Write();


  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
                                  it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();



 
  output->Close();

}
