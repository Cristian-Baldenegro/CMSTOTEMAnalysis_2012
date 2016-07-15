
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

//ROOUNFOLD CLASSES
#include "/storage/lhuertas/uerj-1/CMSTOTEM/mc/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/storage/lhuertas/uerj-1/CMSTOTEM/mc/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/storage/lhuertas/uerj-1/CMSTOTEM/mc/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"



using namespace std;



void pomwig_fullsim(string const& outputFileName = "/storage/lhuertas/uerj-1/CMSTOTEM/mc/Workspace/root_files/pomwig_fullsim.root", const Int_t nevt_max = 1000000){
  
  string treeName = "TotemNtuple";

  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;

  // Declaration of histograms
  map<string,TH1F*> histosTH1F;
  
  Float_t tbins_4[9] = {0.03, 0.08, 0.13, 0.21, 0.31, 0.41,  0.55, 0.75, 1.};
  Float_t tbins_2[11] = {0.03, 0.07, 0.11, 0.15, 0.20, 0.25, 0.32, 0.42, 0.52, 0.65, 1.};

  float xibins[12] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2};

  histosTH1F["proton_right_t_rec"] = new TH1F("proton_right_t_rec", "-t" , 8, tbins_4);
  histosTH1F["proton_right_t_rec_bin2"] = new TH1F("proton_right_t_rec_bin2", "-t" , 10, tbins_2);
  histosTH1F["proton_right_t_rec_kin"] = new TH1F("proton_right_t_rec_kin", "-t" , 8, tbins_4);
  histosTH1F["proton_right_t_rec_kin_bin2"] = new TH1F("proton_right_t_rec_kin_bin2", "-t" , 10, tbins_2);
  histosTH1F["proton_right_t_gen"] = new TH1F("proton_right_t_gen", "-t" , 8, tbins_4);
  histosTH1F["proton_right_t_gen_kin"] = new TH1F("proton_right_t_gen_kin", "-t" , 8, tbins_4);
  histosTH1F["proton_right_xi_rec"] = new TH1F("proton_right_xi_rec", "#xi Right RPs" , 11, xibins);
  histosTH1F["proton_right_xi_rec_kin"] = new TH1F("proton_right_xi_rec_kin", "#xi Right RPs" , 11, xibins);
  histosTH1F["proton_right_xi_gen"] = new TH1F("proton_right_xi_gen", "#xi Right RPs" , 11, xibins);
  histosTH1F["proton_right_xi_gen_kin"] = new TH1F("proton_right_xi_gen_kin", "#xi Right RPs" , 11, xibins);

  histosTH1F["proton_left_t_rec"] = new TH1F("proton_left_t_rec", "-t" , 8, tbins_4);
  histosTH1F["proton_left_t_rec_bin2"] = new TH1F("proton_left_t_rec_bin2", "-t" , 10, tbins_2);
  histosTH1F["proton_left_t_rec_kin"] = new TH1F("proton_left_t_rec_kin", "-t" , 8, tbins_4);
  histosTH1F["proton_left_t_rec_kin_bin2"] = new TH1F("proton_left_t_rec_kin_bin2", "-t" , 10, tbins_2);
  histosTH1F["proton_left_t_gen"] = new TH1F("proton_left_t_gen", "-t" , 8, tbins_4);
  histosTH1F["proton_left_t_gen_kin"] = new TH1F("proton_left_t_gen_kin", "-t" , 8, tbins_4);
  histosTH1F["proton_left_xi_rec"] = new TH1F("proton_left_xi_rec", "#xi Right RPs" , 11, xibins);
  histosTH1F["proton_left_xi_rec_kin"] = new TH1F("proton_left_xi_rec_kin", "#xi Right RPs" , 11, xibins);
  histosTH1F["proton_left_xi_gen"] = new TH1F("proton_left_xi_gen", "#xi Right RPs" , 11, xibins);
  histosTH1F["proton_left_xi_gen_kin"] = new TH1F("proton_left_xi_gen_kin", "#xi Right RPs" , 11, xibins);

  map<string,TH2F*> histosTH2F;
  histosTH2F["proton_y_vs_x_rp_024_025"] = new TH2F("proton_y_vs_x_rp_024_025","proton_y_vs_x_rp_024_025",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_024_025_accept"] = new TH2F("proton_y_vs_x_rp_024_025_accept","proton_y_vs_x_rp_024_025",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_124_125"] = new TH2F("proton_y_vs_x_rp_124_125","proton_y_vs_x_rp_124_125",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_124_125_accept"] = new TH2F("proton_y_vs_x_rp_124_125_accept","proton_y_vs_x_rp_124_125",100,-10,10,100,-40,40);
  histosTH2F["t_xi_gen_minus"] = new TH2F("t_xi_gen_minus", "t_xi" , 26, 0, 1, 26, 0, 0.1);
  histosTH2F["t_xi_gen_minus_xicut"] = new TH2F("t_xi_gen_minus_xicut", "t_xi" , 15, 0, 1., 15, 0, 0.1);
  histosTH2F["t_xi_gen_minus_bin"] = new TH2F("t_xi_gen_minus_bin", "t_xi" , 15, 0,1, 15, 0, 0.1);
  histosTH2F["t_xi_rec_minus"] = new TH2F("t_xi_rec_minus", "t_xi" , 26, 0, 1, 26, 0, 0.1);
  histosTH2F["t_xi_rec_minus_xicut"] = new TH2F("t_xi_rec_minus_xicut", "t_xi" , 15, 0, 1., 15, 0, 0.1);
  histosTH2F["t_xi_rec_minus_bin"] = new TH2F("t_xi_rec_minus_bin", "t_xi" , 15, 0, 1., 15, 0, 0.1);

  RooUnfoldResponse response_xi_right_test (histosTH1F["proton_right_xi_rec_kin"],histosTH1F["proton_right_xi_gen_kin"],"unfolded_xi_right_test","unfolded_xi_right_test");
  RooUnfoldResponse response_xi_right (histosTH1F["proton_right_xi_rec_kin"],histosTH1F["proton_right_xi_gen_kin"], "unfolded_xi_right","unfolded_xi_right");
  histosTH2F["response_xi_right"] = (TH2F*) response_xi_right.Hresponse();
  //RooUnfoldResponse response_t_right (histosTH1F["proton_right_t_rec_kin"],histosTH1F["proton_right_t_gen_kin"],"unfolded_t_right","unfolded_t_right");
  //histosTH2F["response_t_right"] = (TH2F*) response_t_right.Hresponse();
  //RooUnfoldResponse response_t_right_test (histosTH1F["proton_right_t_rec_kin"],histosTH1F["proton_right_t_gen_kin"],"unfolded_t_right_test","unfolded_t_right_test");
   histosTH2F["response_t_right"] = new TH2F("response_t_right","response_t_right",8, tbins_4,8, tbins_4);
   TH2F* response_t_right_test = new TH2F("response_t_right_test","response_t_right_test",8, tbins_4,8, tbins_4);

  for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
      it->second->Sumw2();
  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
  //===================
  int i_tot = 0 , nevt_tot = 0;

 
   vector<TString>* vfiles = new vector<TString>;
   vfiles->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/mc/RPPomwigSDbeta90energy4p0TeV_July2012_totem_ntuple.root/");
   
 
  //Declaration of tree and its branches variables
//   TTree* tree = new TTree(treeName.c_str(),"");
  TTree* tree = NULL;
  MyEvtId*           evtId        = NULL;
  MyL1TrigOld*       l1Trig       = NULL;  
  MyHLTrig*          hltTrig      = NULL;
  //vector<MyGenPart>* genPart      = NULL;
  //vector<MyTracks>*  track_coll   = NULL;
  vector<MyVertex>*  vertex_coll  = NULL;
  vector<MyPFJet>*   pfJet_coll   = NULL;
  vector<MyPFCand>*  pFlow_coll   = NULL;
  vector<MyFSCHit>*  fscHits_coll = NULL;
  vector<MyFSCDigi>* fscDigis_coll = NULL;
  //===================
  T2Event* t2_event = NULL;
  RPRootDumpReconstructedProton* rec_proton_left  = NULL;
  RPRootDumpReconstructedProton* rec_proton_right = NULL;
  RPRootDumpReconstructedProton* gen_proton_left  = NULL;
  RPRootDumpReconstructedProton* gen_proton_right = NULL;
  RPRootDumpReconstructedProtonPair* rec_proton_pair  = NULL;
  map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
  map<unsigned int, RPRootDumpDigiInfo*> rp_digi_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_par_patterns_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_nonpar_patterns_info;
  map<unsigned int, std::vector<RPRootDumpTrackInfo>*> rp_multi_track_info;
  //===================  

  
  TFile* file = new TFile("/storage/lhuertas/uerj-1/CMSTOTEM/samples/mc/RPPomwigSDbeta90energy4p0TeV_July2012_totem_ntuple.root");
    
  //getting the tree form the current file
  tree = (TTree*) file->Get( treeName.c_str() );

  //Getting number of events
  int nev = int(tree->GetEntriesFast());
  nevt_tot += nev;
  cout <<"The current file has " << nev << " entries  "<< endl;

  //adding branches to the tree ----------------------------------------------------------------------
  //tree->SetBranchAddress("cmsFSCHitsUA",&fscHits_coll);
  //tree->SetBranchAddress("cmsFSCDigisUA",&fscDigis_coll);
  tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
  tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
  tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);
  tree->SetBranchAddress("sim_prot_left.",&gen_proton_left);
  tree->SetBranchAddress("sim_prot_right.",&gen_proton_right);
//     tree->SetBranchAddress("Evt",&evtId);
//     tree->SetBranchAddress("L1TrigOld",&l1Trig);
//     tree->SetBranchAddress("HLTrig",&hltTrig);
//     tree->SetBranchAddress("generalTracks",&track_coll); 
//     tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
//     tree->SetBranchAddress("ak5PFJets",&pfJet_coll);
//     tree->SetBranchAddress("particleFlow",&pFlow_coll);
  //if(isMC) tree->SetBranchAddress("genPart",&genPart);
  std::vector<unsigned int> rp_list;
  rp_list.push_back(20); rp_list.push_back(21); rp_list.push_back(24); rp_list.push_back(25);
  rp_list.push_back(120); rp_list.push_back(121); rp_list.push_back(124); rp_list.push_back(125);
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
 

  //starting loop over events, stops when reached end of file or nevt_max
  for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){

      //printing the % of events done every 10k evts
      if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i_evt);

      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
      double event_weight = 1.;
      //double event_weight_eff;
      //double event_weight_averagept_eff;

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

      bool rp_track_accept_right = ( rp_track_valid_120 && rp_track_valid_124 && cut_rp_124 )  || ( rp_track_valid_121 && rp_track_valid_125 && cut_rp_125 );// || ( track_valid_122 && track_valid_123 );
      bool rp_track_accept_left = ( rp_track_valid_020 && rp_track_valid_024 && cut_rp_024 )  || ( rp_track_valid_021 && rp_track_valid_025 && cut_rp_025 );// || ( track_valid_122 && track_valid_123 );

      // RP protons
      double xi_rec_proton_right = rec_proton_right->xi;//cout<<xi_proton_right<<endl;
      double t_rec_proton_right = rec_proton_right->t;
      bool xi_rec_region_right = -xi_rec_proton_right<=0.1;
      bool t_rec_region_right = fabs(t_rec_proton_right)>=0.03 && fabs(t_rec_proton_right)<=1;
      bool good_rec_proton_right = rec_proton_right->valid; // && (xi_proton_right < 0.);

      double xi_rec_proton_left = rec_proton_left->xi;//cout<<xi_proton_left<<endl;
      double t_rec_proton_left = rec_proton_left->t;
      bool xi_rec_region_left = -xi_rec_proton_left<=0.1;
      bool t_rec_region_left = fabs(t_rec_proton_left)>=0.03 && fabs(t_rec_proton_left)<=1;
      bool good_rec_proton_left = rec_proton_left->valid; // && (xi_proton_left < 0.);

      double xi_gen_proton_right = gen_proton_right->xi;//cout<<xi_proton_right<<endl;
      double t_gen_proton_right = gen_proton_right->t;
      bool xi_gen_region_right = -xi_gen_proton_right<=0.1;
      bool t_gen_region_right = fabs(t_gen_proton_right)>=0.03 && fabs(t_gen_proton_right)<=1;
      bool good_gen_proton_right = gen_proton_right->valid; // && (xi_proton_right < 0.);

      double xi_gen_proton_left = gen_proton_left->xi;//cout<<xi_proton_left<<endl;
      double t_gen_proton_left = gen_proton_left->t;
      bool xi_gen_region_left = -xi_gen_proton_left<=0.1;
      bool t_gen_region_left = fabs(t_gen_proton_left)>=0.03 && fabs(t_gen_proton_left)<=1;
      bool good_gen_proton_left = gen_proton_left->valid; // && (xi_proton_left < 0.);

      if(good_rec_proton_right && rp_track_accept_right){
         histosTH1F["proton_right_xi_rec"]->Fill(-xi_rec_proton_right, event_weight );
         histosTH1F["proton_right_t_rec"]->Fill( fabs(t_rec_proton_right), event_weight );
         histosTH1F["proton_right_t_rec_bin2"]->Fill( fabs(t_rec_proton_right), event_weight );
         if (t_rec_region_right && xi_rec_region_right){
                //response_t_right_test->Fill(fabs(t_rec_proton_right),fabs(t_gen_proton_right),event_weight);
             //if (!(t_gen_region_right && xi_gen_region_right)) response_t_right_test.Fake(fabs(t_rec_proton_right),event_weight);
             //if (t_rec_region_right && !t_gen_region_right) response_t_right.Fake(fabs(t_rec_proton_right),event_weight);
             if (-xi_rec_proton_right && !-xi_gen_proton_right) response_xi_right.Fake(-xi_rec_proton_right,event_weight);
             if (-xi_rec_proton_right && !-xi_gen_proton_right) response_xi_right_test.Fake(-xi_rec_proton_right,event_weight);
             histosTH1F["proton_right_xi_rec_kin"]->Fill( -xi_rec_proton_right, event_weight );
             histosTH1F["proton_right_t_rec_kin"]->Fill( fabs(t_rec_proton_right), event_weight );
             histosTH1F["proton_right_t_rec_kin_bin2"]->Fill( fabs(t_rec_proton_right), event_weight );
             histosTH2F["t_xi_rec_minus"]->Fill( fabs(t_rec_proton_right), -xi_rec_proton_right, event_weight);
             histosTH2F["t_xi_rec_minus_bin"]->Fill( fabs(t_rec_proton_right), -xi_rec_proton_right, event_weight );
             if (-xi_rec_proton_right>=0.03) histosTH2F["t_xi_rec_minus_xicut"]->Fill( fabs(t_rec_proton_right), -xi_rec_proton_right, event_weight);
             if (t_gen_region_right && xi_gen_region_right){
                histosTH2F["response_t_right"]->Fill( fabs(t_rec_proton_right), fabs(t_gen_proton_right), event_weight );
                //response_t_right.Fill(fabs(t_rec_proton_right),fabs(t_gen_proton_right),event_weight);
                response_xi_right.Fill(-xi_rec_proton_right,-xi_gen_proton_right,event_weight);
                response_xi_right_test.Fill(-xi_rec_proton_right,-xi_gen_proton_right,event_weight);
             } 
         }
      }
      if(good_rec_proton_left && rp_track_accept_left){
         histosTH1F["proton_left_xi_rec"]->Fill(-xi_rec_proton_left, event_weight );
         histosTH1F["proton_left_t_rec"]->Fill( fabs(t_rec_proton_left), event_weight );
         histosTH1F["proton_left_t_rec_bin2"]->Fill( fabs(t_rec_proton_left), event_weight );
         if (t_rec_region_left && xi_rec_region_left){
             histosTH1F["proton_left_xi_rec_kin"]->Fill( -xi_rec_proton_left, event_weight );
             histosTH1F["proton_left_t_rec_kin"]->Fill( fabs(t_rec_proton_left), event_weight );
             histosTH1F["proton_left_t_rec_kin_bin2"]->Fill( fabs(t_rec_proton_left), event_weight );
         }
      }


      if(good_gen_proton_right){
         histosTH1F["proton_right_xi_gen"]->Fill(-xi_gen_proton_right, event_weight );
         histosTH1F["proton_right_t_gen"]->Fill( fabs(t_gen_proton_right), event_weight );
         if (t_gen_region_right && xi_gen_region_right){
             //if (!(t_rec_region_right && xi_rec_region_right)) response_t_right_test.Miss(fabs(t_gen_proton_right),event_weight);
             //if (!t_rec_region_right && t_gen_region_right) response_t_right.Miss(fabs(t_gen_proton_right),event_weight);
             if (!-xi_rec_proton_right && -xi_gen_proton_right) response_xi_right.Miss(-xi_gen_proton_right,event_weight);
             if (!-xi_rec_proton_right) response_xi_right_test.Miss(-xi_gen_proton_right,event_weight);
             histosTH1F["proton_right_xi_gen_kin"]->Fill( -xi_gen_proton_right, event_weight );
             histosTH1F["proton_right_t_gen_kin"]->Fill( fabs(t_gen_proton_right), event_weight );
             histosTH2F["t_xi_gen_minus"]->Fill( fabs(t_gen_proton_right), -xi_gen_proton_right, event_weight);
             histosTH2F["t_xi_gen_minus_bin"]->Fill( fabs(t_gen_proton_right), -xi_gen_proton_right, event_weight );
             if (-xi_gen_proton_right>=0.03) histosTH2F["t_xi_gen_minus_xicut"]->Fill( fabs(t_gen_proton_right), -xi_gen_proton_right, event_weight);
         }
      }
      if(good_gen_proton_left){
         histosTH1F["proton_left_xi_gen"]->Fill(-xi_gen_proton_left, event_weight );
         histosTH1F["proton_left_t_gen"]->Fill( fabs(t_gen_proton_left), event_weight );
         if (t_gen_region_left && xi_gen_region_left){
             histosTH1F["proton_left_xi_gen_kin"]->Fill( -xi_gen_proton_left, event_weight );
             histosTH1F["proton_left_t_gen_kin"]->Fill( fabs(t_gen_proton_left), event_weight );
         }
      }


  }//end loop for events
  file->Close();
  
  //output file
  TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  output->cd();

  RooUnfoldResponse matrix (histosTH1F["proton_right_t_rec_kin"],histosTH1F["proton_right_t_gen_kin"],response_t_right_test);
  RooUnfoldBinByBin t_jjp_minus_unfold_test (&matrix, histosTH1F["proton_right_t_rec_kin"]);
  TH1F* t_jjp_minus_data_unfolded_test = (TH1F*) t_jjp_minus_unfold_test.Hreco();
  t_jjp_minus_data_unfolded_test->Write();

  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();
                                  it_histo != histosTH1F.end(); ++it_histo)
     (*it_histo).second->Write();
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
                                  it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();

 
  output->Close();

}
