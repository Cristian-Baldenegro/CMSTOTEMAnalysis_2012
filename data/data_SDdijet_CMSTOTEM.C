
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
using namespace std;

Double_t fFermiLike(Double_t *x, Double_t *par) {
  Double_t result = 0;
  result = 1.0/(TMath::Exp((par[0]-TMath::Sqrt(x[0]))/par[1]) + 1);
  return result;
}


void data_SDdijet_CMSTOTEM(string const& outputFileName = "/storage/lhuertas/uerj-1/CMSTOTEM/data/root_files/data_SDdijet_CMSTOTEM.root", const Int_t nevt_max = -1){
  
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


  double weight_st, xi_cms_st, xi_totem_st, xi_totem_sel, t_right_totem_sel, t_left_totem_sel, xi_cms_minus_totem_st;
  TTree* small_tree = new TTree("small_tree","");
  small_tree->Branch("weight",&weight_st,"weight/D");
  small_tree->Branch("xi_cms",&xi_cms_st,"xi_cms/D");
  small_tree->Branch("xi_totem",&xi_totem_st,"xi_totem/D");
  small_tree->Branch("xi_totem_sel",&xi_totem_sel,"xi_totem_sel/D");
  small_tree->Branch("t_right_totem_sel",&t_right_totem_sel,"t_right_totem_sel/D");
  small_tree->Branch("t_left_totem_sel",&t_left_totem_sel,"t_left_totem_sel/D");
  small_tree->Branch("xi_cms_minus_totem",&xi_cms_minus_totem_st,"xi_cms_minus_totem/D");

  float xbins = 50;
  float bin2 = 15;
  float bin_logx = 15;
//  Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.22, 0.30, 0.40, 0.50, 0.65, 1.};
  Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.16, 0.20, 0.25, 0.32, 0.42, 0.52, 0.65, 1.};
  Float_t tbins_2[11] = {0.03, 0.07, 0.11, 0.15, 0.20, 0.25, 0.32, 0.42, 0.52, 0.65, 1.};
  Float_t tbins_3[10] = {0.03, 0.07, 0.11, 0.18, 0.27, 0.36,  0.45, 0.56, 0.75, 1.};
  Float_t tbins_4[9] = {0.03, 0.08, 0.13, 0.21, 0.31, 0.41,  0.55, 0.75, 1.};
//Float_t x[13] = { -0.05, -0.03, -0.01, 0.005, 0.01, 0.03, 0.05, 0.07, 0.09, 0.12, 0.15, 0.17, 0.2};
  Float_t xi_bins_2[8] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.2}; 
  Float_t xi_bins_old[9] = {-0.05, 0, 0.015, 0.03, 0.047, 0.064, 0.08, 0.1, 0.2}; 
  float_t bin_sasha[9] ={0.0003, 0.002, 0.0045, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1};

  TH1F* event_selection = new TH1F("event_selection", "event_selection", 10, 0, 10);

  int nBinsEventSelection = selections.size();
  histosTH1F["EventSelection"] = new TH1F("EventSelection","EventSelection",nBinsEventSelection,0,nBinsEventSelection);
  for(size_t k = 0; k < selections.size(); ++k)
     histosTH1F["EventSelection"]->GetXaxis()->SetBinLabel( (k + 1), selections[k].c_str() );

  histosTH1F["bunchCrossingNumber"] = new TH1F("bunchCrossingNumber", "bunchCrossingNumber" , 3900 , 0 , 3900);

  histosTH1F["decisionPhysTrig"] = new TH1F("decisionPhysTrig", "decisionPhysTrig" , 128 , 0 , 128);
  histosTH1F["decisionTechTrig"] = new TH1F("decisionTechTrig", "decisionTechTrig" , 64 , 0 , 64);

  int nBinsHLT = hltPathNames.size(); 
  histosTH1F["hltTrigFired"] = new TH1F("hltTrigFired", "hltTrigFired" , nBinsHLT , 0 , nBinsHLT);
  for(size_t k = 0; k < nBinsHLT; ++k) 
     histosTH1F["hltTrigFired"]->GetXaxis()->SetBinLabel( (k + 1) , hltPathNames[k].c_str() );

  histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
  histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);
  histosTH1F["vertex_multiplicity"] = new TH1F("vertex_multiplicity", "n vertex" , 100 , 0 , 100);
  histosTH1F["proton_right_vtx_zpos_jetcut"] = new TH1F("proton_right_vtx_zpos_jetcut", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["proton_right_vtx_zpos_protontag"] = new TH1F("proton_right_vtx_zpos_protontag", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["proton_right_vtx_zpos_protonkin"] = new TH1F("proton_right_vtx_zpos_protonkin", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["proton_right_vtx_zpos_backgcut"] = new TH1F("proton_right_vtx_zpos_backgcut", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["proton_left_vtx_zpos_protontag"] = new TH1F("proton_left_vtx_zpos_protontag", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["proton_left_vtx_zpos_protonkin"] = new TH1F("proton_left_vtx_zpos_protonkin", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["proton_left_vtx_zpos_backgcut"] = new TH1F("proton_left_vtx_zpos_backgcut", "z(vtx)" , 150 , -30. , 30.);

  //histosTH1F["pt_gen"] = new TH1F("pt_gen" , "pt_gen;pt;nTracks" , 120 , 0 , 6);
  histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 150 , 0. , 15.);
  histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 200 , -5.2 , 5.2);
  histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 200 , -M_PI , M_PI);
  histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);
  
  histosTH1F["jet_pt"] = new TH1F("jet_pt", "p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["jet_eta"] = new TH1F("jet_eta", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["jet_phi"] = new TH1F("jet_phi", "#phi(jet)" , 20 , -M_PI , M_PI);

  histosTH1F["proton_right_ptjet1_jetcut"] = new TH1F("proton_right_ptjet1_jetcut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_right_ptjet1_protontag"] = new TH1F("proton_right_ptjet1_protontag", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_right_ptjet1_protonkin"] = new TH1F("proton_right_ptjet1_protonkin", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_right_ptjet1_backgcut"] = new TH1F("proton_right_ptjet1_backgcut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_right_etajet1_jetcut"] = new TH1F("proton_right_etajet1_jetcut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_right_etajet1_protontag"] = new TH1F("proton_right_etajet1_protontag", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_right_etajet1_protonkin"] = new TH1F("proton_right_etajet1_protonkin", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_right_etajet1_backgcut"] = new TH1F("proton_right_etajet1_backgcut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_left_ptjet1_protontag"] = new TH1F("proton_left_ptjet1_protontag", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_left_ptjet1_protonkin"] = new TH1F("proton_left_ptjet1_protonkin", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_left_ptjet1_backgcut"] = new TH1F("proton_left_ptjet1_backgcut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_left_etajet1_protontag"] = new TH1F("proton_left_etajet1_protontag", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_left_etajet1_protonkin"] = new TH1F("proton_left_etajet1_protonkin", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_left_etajet1_backgcut"] = new TH1F("proton_left_etajet1_backgcut", "#eta(jet)" , 20 , -5.2 , 5.2);

  histosTH1F["leadingJet_pt"] = new TH1F("leadingJet_pt", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_eta"] = new TH1F("leadingJet_eta", "eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["leadingJet_pt_right_signal"] = new TH1F("leadingJet_pt_right_signal", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_right_signal_kint_kinxi"] = new TH1F("leadingJet_pt_right_signal_kint_kinxi", "p_{T}(jet)" , bin2, 0. , 200.);
  histosTH1F["leadingJet_eta_right_signal_kint_kinxi"] = new TH1F("leadingJet_eta_right_signal_kint_kinxi", "eta(jet)" , bin2, 0. , 200.);
  histosTH1F["leadingJet_eta_right_signal"] = new TH1F("leadingJet_eta_right_signal", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_phi_right_signal"] = new TH1F("leadingJet_phi_right_signal", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["leadingJet_eta_right_signal_kint_kinxi"] = new TH1F("leadingJet_eta_right_signal_kint_kinxi", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_phi_right_signal_kint_kinxi"] = new TH1F("leadingJet_phi_right_signal_kint_kinxi", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["leadingJet_pt_left_signal"] = new TH1F("leadingJet_pt_left_signal", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_left_signal_kint_kinxi"] = new TH1F("leadingJet_pt_left_signal_kint_kinxi", "p_{T}(jet)" , bin2, 0. , 200.);

  histosTH1F["proton_right_ptjet2_jetcut"] = new TH1F("proton_right_ptjet2_jetcut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_right_ptjet2_protontag"] = new TH1F("proton_right_ptjet2_protontag", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_right_ptjet2_protonkin"] = new TH1F("proton_right_ptjet2_protonkin", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_right_ptjet2_backgcut"] = new TH1F("proton_right_ptjet2_backgcut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_right_etajet2_jetcut"] = new TH1F("proton_right_etajet2_jetcut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_right_etajet2_protontag"] = new TH1F("proton_right_etajet2_protontag", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_right_etajet2_protonkin"] = new TH1F("proton_right_etajet2_protonkin", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_right_etajet2_backgcut"] = new TH1F("proton_right_etajet2_backgcut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_left_ptjet2_protontag"] = new TH1F("proton_left_ptjet2_protontag", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_left_ptjet2_protonkin"] = new TH1F("proton_left_ptjet2_protonkin", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_left_ptjet2_backgcut"] = new TH1F("proton_left_ptjet2_backgcut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["proton_left_etajet2_protontag"] = new TH1F("proton_left_etajet2_protontag", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_left_etajet2_protonkin"] = new TH1F("proton_left_etajet2_protonkin", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["proton_left_etajet2_backgcut"] = new TH1F("proton_left_etajet2_backgcut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_pt"] = new TH1F("secondJet_pt", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_eta"] = new TH1F("secondJet_eta", "eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["secondJet_pt_right_signal"] = new TH1F("secondJet_pt_right_signal", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_right_signal_kint_kinxi"] = new TH1F("secondJet_pt_right_signal_kint_kinxi", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_eta_right_signal_kint_kinxi"] = new TH1F("secondJet_eta_right_signal_kint_kinxi", "eta(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_eta_right_signal"] = new TH1F("secondJet_eta_right_signal", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_phi_right_signal"] = new TH1F("secondJet_phi_right_signal", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["secondJet_eta_right_signal_kint_kinxi"] = new TH1F("secondJet_eta_right_signal_kint_kinxi", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_phi_right_signal_kint_kinxi"] = new TH1F("secondJet_phi_right_signal_kint_kinxi", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["secondJet_pt_left_signal"] = new TH1F("secondJet_pt_left_signal", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_left_signal_kint_kinxi"] = new TH1F("secondJet_pt_left_signal_kint_kinxi", "p_{T}(jet)" , bin2 , 0. , 200.);

  histosTH1F["DeltaPtJet_right_signal"] = new TH1F("Delta_pt_Jet_right_signal", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_right_signal"] = new TH1F("Delta_eta_Jet_right_signal", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_right_signal"] = new TH1F("Delta_phi_Jet_right_signal", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["DeltaPtJet_right_signal_kint_kinxi"] = new TH1F("Delta_pt_Jet_right_signal_kint_kinxi", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_right_signal_kint_kinxi"] = new TH1F("Delta_eta_Jet_right_signal_kint_kinxi", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_right_signal_kint_kinxi"] = new TH1F("Delta_phi_Jet_right_signal_kint_kinxi", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["EtaJet_average_right_signal"] = new TH1F("eta_Jet_average_right_signal", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_right_signal_kint_kinxi"] = new TH1F("eta_Jet_average_right_signal_kint_kinxi", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["DeltaPtJet_left_signal"] = new TH1F("Delta_pt_Jet_left_signal", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_left_signal"] = new TH1F("Delta_eta_Jet_left_signal", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_left_signal"] = new TH1F("Delta_phi_Jet_left_signal", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["DeltaPtJet_left_signal_kint_kinxi"] = new TH1F("Delta_pt_Jet_left_signal_kint_kinxi", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_left_signal_kint_kinxi"] = new TH1F("Delta_eta_Jet_left_signal_kint_kinxi", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_left_signal_kint_kinxi"] = new TH1F("Delta_phi_Jet_left_signal_kint_kinxi", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["EtaJet_average_left_signal"] = new TH1F("eta_Jet_average_left_signal", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_signal_kint_kinxi"] = new TH1F("eta_Jet_average_left_signal_kint_kinxi", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["massJets"] = new TH1F("massJet", "M_{jet}" , bin2 , 0 , 450);
  histosTH1F["massJet_right_signal_kint_kinxi"] = new TH1F("massJet_right_signal_kint_kinxi", "M_{jet}" , bin2 , 0 , 450);
  histosTH1F["mass_x_right_signal_kint_kinxi"] = new TH1F("mass_x_right_signal_kint_kinxi", "M_{X}" , bin2 , 0 , 8000);

  histosTH1F["Eta_max"] = new TH1F("Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Eta_min"] = new TH1F("Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 20 , 0,11);
  histosTH1F["xi_plus_Reco"] = new TH1F("xi_plus", "#xi^{+}" , 7, xi_bins_2);
  histosTH1F["xi_plus_Reco_bin"] = new TH1F("xi_plus_bin", "#xi^{+}" , 8, bin_sasha);
  histosTH1F["log_xi_plus_Reco"] = new TH1F("log_xi_plus", "Log #xi^{+}" , 20 , -3,0.5);
  histosTH1F["log_x_minus_before"] = new TH1F("log_x_minus_before", "Log x^{-}" , bin_logx , -4, 0);
  histosTH1F["log_x_plus_before"] = new TH1F("log_x_plus_before", "Log x^{+}" , bin_logx , -4, 0);
  histosTH1F["log_x_minus"] = new TH1F("log_x_minus", "Log x^{-}" , bin_logx , -4, 0);
  histosTH1F["log_x_minus_kint_kinxi"] = new TH1F("log_x_minus_kint_kinxi", "Log x^{-}" , bin_logx , -4, 0);
  histosTH1F["log_x_plus"] = new TH1F("log_x_plus", "Log x^{+}" , bin_logx , -4, 0);
  histosTH1F["log_x_plus_kint_kinxi"] = new TH1F("log_x_plus_kint_kinxi", "Log x^{+}" , bin_logx , -4, 0);
  histosTH1F["log_x_minus_jj"] = new TH1F("log_x_minus_jj", "Log x^{-}" , bin_logx , -4, 0);
  histosTH1F["log_x_minus_sel_jj"] = new TH1F("log_x_minus_sel_jj", "Log x^{-}" , bin_logx , -4, 0);
  histosTH1F["log_x_plus_jj"] = new TH1F("log_x_plus_jj", "Log x^{+}" , bin_logx , -4, 0);
  histosTH1F["log_x_plus_sel_jj"] = new TH1F("log_x_plus_sel_jj", "Log x^{+}" , bin_logx , -4, 0);
    
  histosTH1F["energy_pfcand"] = new TH1F("energy_pfcand", "energy" , 50, 0, 8000);
  histosTH1F["pz_pfcand"] = new TH1F("pz_pfcand", "pz" , 50, -4000, 4000);
  histosTH1F["energy_pfcand_thresholds"] = new TH1F("energy_pfcand_thresholds", "energy" , 50, 0, 8000);
  histosTH1F["pz_pfcand_thresholds"] = new TH1F("pz_pfcand_thresholds", "pz" , 50, -4000, 4000);

  histosTH1F["fscHit_energy"] = new TH1F("fscHit_energy", "FSC hit energy" , 150 , -100. , 200.);
  histosTH1F["fscHit_time"] = new TH1F("fscHit_time", "FSC hit time" , 150 , 0. , 300.);

  histosTH1F["t2_track_chi2Prob_zplus"] = new TH1F("t2_track_chi2Prob_zplus", "#chi^{2}" , 100 , 0. , 1.);
  histosTH1F["t2_track_entryX_zplus"] = new TH1F("t2_track_entryX_zplus", "x_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_entryY_zplus"] = new TH1F("t2_track_entryY_zplus", "y_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_multiplicity_zplus"] = new TH1F("t2_track_multiplicity_zplus", "n tracks" , 100 , 0 , 100);
  histosTH1F["t2_track_chi2Prob_zminus"] = new TH1F("t2_track_chi2Prob_zminus", "#chi^{2}" , 100 , 0. , 1.);
  histosTH1F["t2_track_entryX_zminus"] = new TH1F("t2_track_entryX_zminus", "x_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_entryY_zminus"] = new TH1F("t2_track_entryY_zminus", "y_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_multiplicity_zminus"] = new TH1F("t2_track_multiplicity_zminus", "n tracks" , 100 , 0 , 100);

  histosTH1F["proton_right_xi_binsasha"] = new TH1F("proton_right_xi_binsasha", "#xi^{-}" , 8, bin_sasha);
  histosTH1F["xi_minus_Reco"] = new TH1F("xi_minus", "#xi^{-}" , 7, xi_bins_2);
  histosTH1F["xi_minus_Reco_bin"] = new TH1F("xi_minus_bin", "#xi^{-}" , 8, bin_sasha);
  histosTH1F["xi_minus_Reco_gap"] = new TH1F("xi_minus_gap", "#xi^{-}" , 8, bin_sasha);
  histosTH1F["proton_right_t"] = new TH1F("proton_right_t", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_sel"] = new TH1F("proton_right_t_sel", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_cut"] = new TH1F("proton_right_t_cut", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_signal"] = new TH1F("proton_right_t_signal", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_totem"] = new TH1F("proton_right_t_totem", "-t" , 10, tbins_2);
  histosTH1F["proton_right_t_signal_kint_kinxi"] = new TH1F("proton_right_t_signal_kint_kinxi", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_signal_kint_kinxi_bin"] = new TH1F("proton_right_t_signal_kint_kinxi_bin", "-t" , 10, tbins_2);
  histosTH1F["proton_right_t_signal_kint_kinxi_bin3"] = new TH1F("proton_right_t_signal_kint_kinxi_bin3", "-t" , 9, tbins_3);
  histosTH1F["proton_right_t_signal_kint_kinxi_bin4"] = new TH1F("proton_right_t_signal_kint_kinxi_bin4", "-t" , 8, tbins_4);
  histosTH1F["proton_right_t_jetcut"] = new TH1F("proton_right_t_jetcut", "-t" , 10, tbins_2);
  histosTH1F["proton_right_t_protontag"] = new TH1F("proton_right_t_protontag", "-t" , 10, tbins_2);
  histosTH1F["proton_right_t_protonkin"] = new TH1F("proton_right_t_protonkin", "-t" , 10, tbins_2);
  histosTH1F["proton_right_t_backgcut"] = new TH1F("proton_right_t_backgcut", "-t" , 10, tbins_2);
  histosTH1F["proton_right_t_true"] = new TH1F("proton_right_t_true", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_true_constbin"] = new TH1F("proton_right_t_true_constbin", "-t" , 20 , 0, 1);
  histosTH1F["proton_right_t_signal_bin"] = new TH1F("proton_right_t_sigal_bin", "-t" , xbins , 0, 1);
  histosTH1F["proton_right_t_halo"] = new TH1F("proton_right_t_halo", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_halo_constbin"] = new TH1F("proton_right_t_halo_constbin", "-t" , 20 , 0, 1);
  histosTH1F["halo_right_pt30"] = new TH1F("halo_right_pt30", "-t halo" , 11 , tbins);
  histosTH1F["halo_right_pt30_constbin"] = new TH1F("halo_right_pt30_constbin", "-t halo" , 20 , 0, 1);
  histosTH1F["proton_right_chi2"] = new TH1F("proton_right_chi2", "#chi^{2}" , 20 , 0. , 100.);
  histosTH1F["proton_right_xi_signal"] = new TH1F("proton_right_xi_signal", "#xi Right RPs" , bin2 , -0.05, 0.2);
  histosTH1F["proton_right_xi_signal_kinxi"] = new TH1F("proton_right_xi_signal_kinxi", "#xi Right RPs" , 8, bin_sasha);
  histosTH1F["proton_right_xi_signal_kinxi_cut"] = new TH1F("proton_right_xi_signal_kinxi_cut", "#xi Right RPs" , 8, bin_sasha);
  histosTH1F["proton_right_xi_signal_kint_kinxi"] = new TH1F("proton_right_xi_signal_kint_kinxi", "#xi Right RPs" , 8, xi_bins_old);
  histosTH1F["proton_right_xi_signal_kint_kinxi_bin"] = new TH1F("proton_right_xi_signal_kint_kinxi_bin", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_right_xi_jetcut"] = new TH1F("proton_right_xi_jetcut", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_right_xi_protontag"] = new TH1F("proton_right_xi_protontag", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_right_xi_protonkin"] = new TH1F("proton_right_xi_protonkin", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_right_xi_backgcut"] = new TH1F("proton_right_xi_backgcut", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_right_xi_totem"] = new TH1F("proton_right_xi_totem", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi Right RPs" , xbins ,  -0.05, 0.2);
  histosTH1F["proton_right_xi_kin"] = new TH1F("proton_right_xi_kin", "#xi Right RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_right_xi_kint"] = new TH1F("proton_right_xi_kint", "#xi Right RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_right_xi_kint_kinxi"] = new TH1F("proton_right_xi_kint_kinxi", "#xi Right RPs" , 25 , -0.05, 0.2);
  histosTH1F["proton_minus_xi_kint_sel"] = new TH1F("proton_minus_xi_kint_sel", "#xi Right RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_minus_xi_kint"] = new TH1F("proton_minus_xi_kint", "#xi Right RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_right_xi_kint_bin25"] = new TH1F("proton_right_xi_kint_bin25", "#xi Right RPs" , 25 , -0.05, 0.2);
  histosTH1F["proton_right_xi_sel"] = new TH1F("proton_right_xi_sel", "#xi Right RPs" , bin2 , -0.05, 0.2);
  histosTH1F["proton_right_xi_halo"] = new TH1F("proton_right_xi_halo", "#xi Right RPS" , xbins , -0.05, 0.2);
  histosTH1F["proton_right_beta"] = new TH1F("proton_right_beta", "#beta Right RPs" , bin2 , 0, 0.9);
  histosTH1F["proton_right_beta_kint_kinxi"] = new TH1F("proton_right_beta_kint_kinxi", "#beta Right RPs" , bin2 , 0, 0.9);
  histosTH1F["proton_right_xi_cut"] = new TH1F("proton_right_xi_tcut", "#xi Right RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_right_logXi"] = new TH1F("proton_right_logXi","log(#xi)",20,-5.,0.);

  histosTH1F["proton_right_thx_signal_kint_kinxi"] = new TH1F("proton_right_thx_signal_kint_kinxi", "theta_x", 50, -1, 1);
  histosTH1F["proton_right_thy_signal_kint_kinxi"] = new TH1F("proton_right_thy_signal_kint_kinxi", "theta_y", 50, -1, 1);

  histosTH1F["proton_left_t_jetcut"] = new TH1F("proton_left_t_jetcut", "-t" , 10, tbins_2);
  histosTH1F["proton_left_t_protontag"] = new TH1F("proton_left_t_protontag", "-t" , 10, tbins_2);
  histosTH1F["proton_left_t_protonkin"] = new TH1F("proton_left_t_protonkin", "-t" , 10, tbins_2);
  histosTH1F["proton_left_t_backgcut"] = new TH1F("proton_left_t_backgcut", "-t" , 10, tbins_2);
  histosTH1F["proton_left_t_signal_kint_kinxi_bin3"] = new TH1F("proton_left_t_signal_kint_kinxi_bin3", "-t" , 9, tbins_3);
  histosTH1F["proton_left_t_signal_kint_kinxi_bin4"] = new TH1F("proton_left_t_signal_kint_kinxi_bin4", "-t" , 8, tbins_4);
  histosTH1F["proton_left_t"] = new TH1F("proton_left_t", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_sel"] = new TH1F("proton_left_t_sel", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_cut"] = new TH1F("proton_left_t_cut", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_signal"] = new TH1F("proton_left_t_signal", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_signal_kint_kinxi"] = new TH1F("proton_left_t_signal_kint_kinxi", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_signal_kint_kinxi_bin"] = new TH1F("proton_left_t_signal_kint_kinxi_bin", "-t" , 10, tbins_2);
  histosTH1F["proton_left_t_true"] = new TH1F("proton_left_t_true", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_true_constbin"] = new TH1F("proton_left_t_true_constbin", "-t" , 20 , 0, 1);
  histosTH1F["proton_left_t_signal_bin"] = new TH1F("proton_left_t_sigal_bin", "-t" , xbins , 0, 1); 
  histosTH1F["proton_left_t_halo"] = new TH1F("proton_left_t_halo", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_halo_constbin"] = new TH1F("proton_left_t_halo_constbin", "-t" , 20 , 0, 1);
  histosTH1F["halo_left_pt30"] = new TH1F("halo_left_pt30", "-t halo" , 11 , tbins);
  histosTH1F["halo_left_pt30_constbin"] = new TH1F("halo_left_pt30_constbin", "-t halo" , 20 , 0, 1);
  histosTH1F["proton_left_chi2"] = new TH1F("proton_left_chi2", "#chi^{2}" , 20 , 0. , 100.);
  histosTH1F["proton_left_xi_signal"] = new TH1F("proton_left_xi_signal", "#xi Left RPs" , bin2 , -0.05, 0.2);
  histosTH1F["proton_left_xi_signal_kint_kinxi"] = new TH1F("proton_left_xi_signal_kint_kinxi", "#xi Left RPs" , bin2 , -0.05, 0.2);
  histosTH1F["proton_left_xi_signal_kint_kinxi_bin"] = new TH1F("proton_left_xi_signal_kint_kinxi_bin", "#xi Left RPs" , 7, xi_bins_2);
  histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi Left RPs" , xbins ,  -0.05, 0.2);
  histosTH1F["proton_left_xi_kin"] = new TH1F("proton_left_xi_kin", "#xi Left RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_left_xi_kint"] = new TH1F("proton_left_xi_kint", "#xi Left RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_left_xi_kint_kinxi"] = new TH1F("proton_left_xi_kint_kinxi", "#xi Left RPs" , 25, -0.05, 0.2);
  histosTH1F["proton_left_xi_kint_bin25"] = new TH1F("proton_left_xi_kint_bin25", "#xi Left RPs" , 25 , -0.05, 0.2);
  histosTH1F["proton_left_xi_sel"] = new TH1F("proton_left_xi_sel", "#xi Left RPs" , bin2 , -0.05, 0.2);
  histosTH1F["proton_left_xi_halo"] = new TH1F("proton_left_xi_halo", "#xi Left RPS" , xbins , -0.05, 0.2);
  histosTH1F["proton_left_beta"] = new TH1F("proton_left_beta", "#beta Left RPs" , bin2 , 0, 0.9);
  histosTH1F["proton_left_beta_kint_kinxi"] = new TH1F("proton_left_beta_kint_kinxi", "#beta Left RPs" , bin2 , 0, 0.9);
  histosTH1F["proton_left_xi_cut"] = new TH1F("proton_left_xi_tcut", "#xi Left RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_left_logXi"] = new TH1F("proton_left_logXi","log(#xi)",20,-5.,0.);
  histosTH1F["proton_plus_xi_kint_kinxi"] = new TH1F("proton_plus_xi_kint_kinxi", "#xi Left RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_plus_xi_kint_kinxi_sel"] = new TH1F("proton_plus_xi_kint_kinxi_sel", "#xi Left RPs" , xbins , -0.05, 0.2);
  histosTH1F["proton_left_xi_jetcut"] = new TH1F("proton_left_xi_jetcut", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_left_xi_protontag"] = new TH1F("proton_left_xi_protontag", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_left_xi_protonkin"] = new TH1F("proton_left_xi_protonkin", "#xi Right RPs" , 7, xi_bins_2);
  histosTH1F["proton_left_xi_backgcut"] = new TH1F("proton_left_xi_backgcut", "#xi Right RPs" , 7, xi_bins_2);

  histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
  histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);
  histosTH1F["pf_xiMinus_plus_proton_left_xi"] = new TH1F("pf_xiMinus_plus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
  
  Float_t bin[16] = {-0.4, -0.112, -0.096, -0.08, -0.064, -0.048, -0.032, -0.016, 0, 0.048, 0.112, 0.176, 0.24, 0.304, 0.368, 0.4};
  histosTH1F["xitotem_xicms_rightRPs"] = new TH1F("xitotem_xicms_rightRPs", "Right RPs" , xbins , -0.4 , 0.4);
  histosTH1F["xitotem_xicms_rightRPs_kin"] = new TH1F("xitotem_xicms_rightRPs_kin", "Right RPs" , xbins, -0.4 , 0.4);
  histosTH1F["xitotem_xicms_rightRPs_kint_kinxi"] = new TH1F("xitotem_xicms_rightRPs_kint_kinxi", "Right RPs" , 15, bin);
  histosTH1F["xitotem_xicms_rightRPs_tcut"] = new TH1F("xitotem_xicms_rightRPs_tcut", "Right RPs" , xbins, -0.4 , 0.4);
  histosTH1F["xitotem_xicms_rightRPs_cut"] = new TH1F("xitotem_xicms_rightRPs_cut", "Right RPs" ,xbins , -0.4 , 0.4);
  histosTH1F["xi_cms_totem_background_simulated_right"] = new TH1F("xitotem_xicms_rightRPs_simulated_right", "Right RPs" , 25 , -0.4 , 0.4);
  histosTH1F["xi_cms_totem_background_simulated_left"] = new TH1F("xitotem_xicms_rightRPs_simulated_left", "Left RPs" , 25, -0.4 , 0.4);
  histosTH1F["xitotem_xicms_leftRPs"] = new TH1F("xitotem_xicms_leftRPs", "Left RPs" , xbins , -0.4 , 0.4);
  histosTH1F["xitotem_xicms_leftRPs_kin"] = new TH1F("xitotem_xicms_leftRPs_kin", "Left RPs" , xbins, -0.4 , 0.4);
  histosTH1F["xitotem_xicms_leftRPs_kint_kinxi"] = new TH1F("xitotem_xicms_leftRPs_kint_kinxi", "Left RPs" , 15, bin);

  
  map<string,TH2F*> histosTH2F;
  histosTH2F["proton_y_vs_x_rp_024_025"] = new TH2F("proton_y_vs_x_rp_024_025","proton_y_vs_x_rp_024_025",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_024_025_accept"] = new TH2F("proton_y_vs_x_rp_024_025_accept","proton_y_vs_x_rp_024_025",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_124_125"] = new TH2F("proton_y_vs_x_rp_124_125","proton_y_vs_x_rp_124_125",100,-10,10,100,-40,40);
  histosTH2F["proton_y_vs_x_rp_124_125_accept"] = new TH2F("proton_y_vs_x_rp_124_125_accept","proton_y_vs_x_rp_124_125",100,-10,10,100,-40,40);
  histosTH2F["t2_track_multiplicity_vs_track_multiplicity"] = new TH2F("t2_track_multiplicity_vs_track_multiplicity","t2_track_multiplicity_vs_track_multiplicity", 100 , 0 , 100, 100 , 0 , 100);
  histosTH2F["t2_track_multiplicity_vs_leadingJet_pt"] = new TH2F("t2_track_multiplicity_vs_leadingJet_pt","t2_track_multiplicity_vs_leadingJet_pt", 150 , 0. , 150., 100 , 0 , 100);
  histosTH2F["t2_track_entryY_vs_entryX_zplus"] = new TH2F("t2_track_entryY_vs_entryX_zplus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);
  histosTH2F["t2_track_entryY_vs_entryX_zminus"] = new TH2F("t2_track_entryY_vs_entryX_zminus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);

  histosTH2F["proton_right_logXi_vs_pf_logXiPlus"] = new TH2F("proton_right_logXi_vs_pf_logXiPlus","proton_right_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_right_logXi_vs_pf_logXiMinus"] = new TH2F("proton_right_logXi_vs_pf_logXiMinus","proton_right_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_right_logXi_vs_t"] = new TH2F("proton_right_logXi_vs_t","proton_right_logXi_vs_t", 200, 0., 5., 200, -5., 0.);
  histosTH2F["proton_right_xi_vs_pf_xiMinus"] = new TH2F("proton_right_xi_vs_pf_xiMinus","proton_right_xi_vs_pf_xiMinus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_right_xi_vs_pf_xiPlus"] = new TH2F("proton_right_xi_vs_pf_xiPlus","proton_right_xi_vs_pf_xiPlus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_right_t_vs_leadingJet_pt"] = new TH2F("proton_right_t_vs_leadingJet_pt","proton_right_t_vs_leadingJet_pt", 150 , 0. , 150., 200 , 0. , 5.);

  histosTH2F["proton_left_logXi_vs_pf_logXiPlus"] = new TH2F("proton_left_logXi_vs_pf_logXiPlus","proton_left_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_left_logXi_vs_pf_logXiMinus"] = new TH2F("proton_left_logXi_vs_pf_logXiMinus","proton_left_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_left_logXi_vs_t"] = new TH2F("proton_left_logXi_vs_t","proton_left_logXi_vs_t", 200, 0., 5., 200, -5., 0.);
  histosTH2F["proton_left_xi_vs_pf_xiMinus"] = new TH2F("proton_left_xi_vs_pf_xiMinus","proton_left_xi_vs_pf_xiMinus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_left_xi_vs_pf_xiPlus"] = new TH2F("proton_left_xi_vs_pf_xiPlus","proton_left_xi_vs_pf_xiPlus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_left_t_vs_leadingJet_pt"] = new TH2F("proton_left_t_vs_leadingJet_pt","proton_left_t_vs_leadingJet_pt", 150 , 0. , 150., 200 , 0. , 5.);
  
  double energyMin = -10.;
  double energyMax = 190.;
  int nBinsEnergy = 1000;
  histosTH2F["energyVsEtaAllTypes"] = new TH2F("energyVsEtaAllTypes","energy Vs Eta AllTypes",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaUndefined"] = new TH2F("energyVsEtaUndefined","energy Vs Eta Undefined",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energy Vs Eta Charged Hadron"] = new TH2F("energyVsEtaChargedHadron","energy Vs Eta Charged Hadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaElectron"] = new TH2F("energyVsEtaElectron","energy Vs Eta Electron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaMuon"] = new TH2F("energyVsEtaMuon","energy Vs Eta Muon",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaGamma"] = new TH2F("energyVsEtaGamma","energy Vs Eta Gamma",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaNeutralHadron"] = new TH2F("energyVsEtaNeutralHadron","energy Vs Eta Neutral Hadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHF"] = new TH2F("energyVsEtaHadronHF","energy Vs Eta HadronHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFEcalEnergy"] = new TH2F("energyVsEtaHadronHFEcalEnergy","energyVsEtaHadronHFEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFNoEcalEnergy"] = new TH2F("energyVsEtaHadronHFNoEcalEnergy","energyVsEtaHadronHFNoEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaEGammaHF"] = new TH2F("energyVsEtaEGammaHF","energy Vs Eta GammaHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["xi+Vseta_max"] = new TH2F("xi+Vseta_max","#xi^{+] Vs #eta_{max}",200,0,5.2,nBinsEnergy,0,0.5);

  for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
      it->second->Sumw2();
  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
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

     weight_st = -1.; 
     xi_cms_st = -999.; xi_totem_st = -999.; xi_totem_sel = -999.; t_right_totem_sel = -999.; t_left_totem_sel = -999.; xi_cms_minus_totem_st = -999.;
 
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
      weight_st = event_weight;

      for (int itrig = 0 ; itrig < 128 ; ++itrig){
         if( l1Trig->PhysTrigWord[itrig] == 1) 
            histosTH1F["decisionPhysTrig"]->Fill( itrig, event_weight );
      }
        
      for (int itrig = 0 ; itrig < 64 ; ++itrig){
         if( l1Trig->TechTrigWord[itrig] == 1 )
            histosTH1F["decisionTechTrig"]->Fill( itrig, event_weight );
      }

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
		 histosTH1F["hltTrigFired"]->Fill( idx, event_weight );
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
      histosTH1F["vtx_xpos"]->Fill(primaryVertex.x, event_weight);
      histosTH1F["vtx_ypos"]->Fill(primaryVertex.y, event_weight);
      histosTH1F["vtx_zpos"]->Fill(primaryVertex.z, event_weight);
      histosTH1F["vertex_multiplicity"]->Fill( n_evt_vtx, event_weight );
      
      // Tracks
      int n_tracks_selected = 0;
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
         histosTH1F["jet_pt"]->Fill( basejet.Pt(), event_weight  );
         if(basejet.Pt() > 0.) histosTH1F["jet_eta"]->Fill( basejet.Eta(), event_weight  );
         histosTH1F["jet_phi"]->Fill( basejet.Phi(), event_weight  );
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

      double x_plus = ((Jet1_E+Jet1_pz)+(Jet2_E+Jet2_pz))/8000;
      double x_minus = ((Jet1_E-Jet1_pz)+(Jet2_E-Jet2_pz))/8000;
      double x_minus_sel = (x_minus<x_plus) ? x_minus : x_plus;
      double x_plus_sel = (x_minus>x_plus) ? x_plus : x_minus;
      double mass_jets= sqrt(pow(Jet1_E+Jet2_E,2)-pow(Jet1_px+Jet2_px,2)-pow(Jet1_py+Jet2_py,2)-pow(Jet1_pz+Jet2_pz,2));
      if (jet1_selected && jet2_selected){
         ++n_evt_jets;
         histosTH1F["leadingJet_pt"]->Fill( Jet1_pt, 1. );
         histosTH1F["secondJet_pt"]->Fill( Jet2_pt, 1. );
         histosTH1F["leadingJet_eta"]->Fill( Jet1_eta, 1. );
         histosTH1F["secondJet_eta"]->Fill( Jet2_eta, 1. );
         histosTH1F["massJets"]->Fill( mass_jets, 1. );
         histosTH1F["log_x_minus_jj"]->Fill( log10(x_minus), 1. );
         histosTH1F["log_x_minus_sel_jj"]->Fill( log10(x_minus_sel), 1. );
         histosTH1F["log_x_plus_jj"]->Fill( log10(x_plus), 1. );
         histosTH1F["log_x_plus_sel_jj"]->Fill( log10(x_plus_sel), 1. );
      }

 
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
         histosTH1F["energy_pfcand"]->Fill( energy, 1. );
         histosTH1F["pz_pfcand"]->Fill( pz, 1. );
	 
         // Apply thresholds
         if( !pflowThreshold(*it_pfcand,thresholdsPFlow) ) continue;
         histosTH1F["energy_pfcand_thresholds"]->Fill( energy, 1. );
         histosTH1F["pz_pfcand_thresholds"]->Fill( pz, 1. );

/*
         // HF eta rings 29, 30, 40, 41
         if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || (fabs(eta) >= 4.730) ) continue;
      //   if( ( (fabs(eta) <= 2.866) && (fabs(eta) > 3.152) ) || (fabs(eta) <= 4.730) ){
	 
	     bool Barrel_region = false;
             bool Endcap_region = false;
             bool Transition_region = false;
             bool Forward_region = false;
	     bool passed_threshold = false;
	 
             if( eta>=-1.4 && eta<=1.4 ) Barrel_region = true;
             else if( (eta>=-2.6 && eta<-1.4) || (eta>1.4 && eta<=2.6) ) Endcap_region = true;
             else if( (eta>=-3.2 && eta<-2.6) || (eta>2.6 && eta<=3.2) ) Transition_region = true;
             else if( eta<-3.2 || eta>3.2 ) Forward_region = true;
             else cout << "ERROR!!!!!!!!!" << endl;

	     //Applying threshold
	     if (Barrel_region == true){
	       //if(partType == MyPFCand::h0 && energy>=1.4) Treshold_Barrel1==true;
	       if(partType == MyPFCand::h0 && energy > 1.4) passed_threshold = true; 
	       if(partType == MyPFCand::gamma && energy>0.9) passed_threshold = true; 
	     }  
	     if (Endcap_region == true){
	       if(partType == MyPFCand::h0 && energy>2.7) passed_threshold = true; 
	       if(partType == MyPFCand::gamma && energy>2.5) passed_threshold = true; 
	     }  
	     if (Transition_region == true){
	       if(partType == MyPFCand::h0 && energy>3.8) passed_threshold = true; 
	       if(partType == MyPFCand::gamma && energy>2.5) passed_threshold = true; 
	       if(partType == MyPFCand::h_HF && energy>4) passed_threshold = true; 
	       if(partType == MyPFCand::egamma_HF && energy>3.5) passed_threshold = true; 
	     }  
	     if (Forward_region == true){
	       if(partType == MyPFCand::h_HF && energy>4) passed_threshold = true; 
	       if(partType == MyPFCand::egamma_HF && energy>3.5) passed_threshold = true; 
	     }  

	 if(jet1_selected && jet2_selected && PF_eta_max && PF_eta_min && passed_threshold ) ++n_evt_pass_threshold;//continue;
	 if(!passed_threshold ) continue;
*/	 
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


      if(jet1_selected && jet2_selected && PF_eta_max && PF_eta_min){
         ++n_evt_etapf;   
         histosTH1F["Eta_max"]->Fill( eta_max, event_weight  );
         histosTH1F["Eta_min"]->Fill( eta_min, event_weight  );
         histosTH1F["Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );
         histosTH1F["xi_plus_Reco"]->Fill( xi_plus_Reco, event_weight  );
         histosTH1F["xi_plus_Reco_bin"]->Fill( xi_plus_Reco, event_weight  );

         if (eta_min>-3) histosTH1F["xi_minus_Reco_gap"]->Fill( xi_minus_Reco, event_weight  );
         histosTH1F["xi_minus_Reco"]->Fill( xi_minus_Reco, event_weight  );
         histosTH1F["xi_minus_Reco_bin"]->Fill( xi_minus_Reco, event_weight  );
         histosTH1F["log_xi_plus_Reco"]->Fill( log10(xi_plus_Reco), event_weight  );
         histosTH2F["xi+Vseta_max"]->Fill( eta_max, xi_plus_Reco, event_weight );
      }
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
         if( zside > 0 )
            histosTH1F["t2_track_chi2Prob_zplus"]->Fill( t2_trk_chiProb[i_t2_trk], event_weight );
         else
            histosTH1F["t2_track_chi2Prob_zminus"]->Fill( t2_trk_chiProb[i_t2_trk], event_weight );

         // Select tracks
         if( t2_trk_chiProb[i_t2_trk] < 0.2 ) continue;

         ++n_t2_tracks_selected;
         if( zside > 0 ) ++n_t2_tracks_selected_zplus;
         else            ++n_t2_tracks_selected_zminus;

         if( zside > 0 ){
            histosTH1F["t2_track_entryX_zplus"]->Fill( t2_trk_entryX[i_t2_trk], event_weight );
            histosTH1F["t2_track_entryY_zplus"]->Fill( t2_trk_entryY[i_t2_trk], event_weight );
            histosTH2F["t2_track_entryY_vs_entryX_zplus"]->Fill( t2_trk_entryX[i_t2_trk], t2_trk_entryY[i_t2_trk], event_weight );
         } else{
            histosTH1F["t2_track_entryX_zminus"]->Fill( t2_trk_entryX[i_t2_trk], event_weight );
            histosTH1F["t2_track_entryY_zminus"]->Fill( t2_trk_entryY[i_t2_trk], event_weight );
            histosTH2F["t2_track_entryY_vs_entryX_zminus"]->Fill( t2_trk_entryX[i_t2_trk], t2_trk_entryY[i_t2_trk], event_weight );
         }
      }
      if( selectZeroHitsT2Plus && (n_t2_tracks_selected_zplus > 0) ) continue;
      histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Plus", event_weight );

      if( selectZeroHitsT2Minus && (n_t2_tracks_selected_zminus > 0) ) continue;
      histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Minus", event_weight );

      bool proton_right_valid = rec_proton_right->valid;
      bool proton_left_valid = rec_proton_left->valid;
      if( selectSingleArmRecProton && (proton_right_valid && proton_left_valid) ) continue;
      histosTH1F["EventSelection"]->Fill( "SingleArmRP", event_weight );

      if( selectDoubleArmRecProton && !(proton_right_valid && proton_left_valid) ) continue;
      histosTH1F["EventSelection"]->Fill( "DoubleArmRP", event_weight );

      bool tag_elastic_top45_bot56 = elastic_top45_bot56(rp_track_info);      
      bool tag_elastic_bot45_top56 = elastic_bot45_top56(rp_track_info);      
      if( selectElastic && !(tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
      histosTH1F["EventSelection"]->Fill( "Elastic", event_weight );

      if( selectNonElastic && (tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
      histosTH1F["EventSelection"]->Fill( "NonElastic", event_weight );

 

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

      bool rp_track_accept_right = ( rp_track_valid_120 && rp_track_valid_124 && cut_rp_124 )  || ( rp_track_valid_121 && rp_track_valid_125 && cut_rp_125 );// || ( rp_track_valid_122 && rp_track_valid_123 );
      bool rp_track_accept_left = ( rp_track_valid_020 && rp_track_valid_024 && cut_rp_024 )  || ( rp_track_valid_021 && rp_track_valid_025 && cut_rp_025 );// || ( rp_track_valid_022 && rp_track_valid_023 );

      // RP protons
      double xi_totem;
      double chi2_proton_right = rec_proton_right->chi2;
      double chindf_proton_right = rec_proton_right->chindf;
      double xi_proton_right = rec_proton_right->xi;//cout<<xi_proton_right<<endl;
      double t_proton_right = rec_proton_right->t;
      double thetax_proton_right = rec_proton_right->thx;
      double thetay_proton_right = rec_proton_right->thy;
      double beta_proton_right = x_minus/-xi_proton_right;
      bool xi_region_right = -xi_proton_right>=0.03 && -xi_proton_right<=0.1;
      bool t_region_right = fabs(t_proton_right)>=0.03 && fabs(t_proton_right)<=1;
      bool good_proton_right = proton_right_valid; // && (xi_proton_right < 0.);
      double chi2_proton_left = rec_proton_left->chi2;
      double chindf_proton_left = rec_proton_left->chindf;
      double xi_proton_left = rec_proton_left->xi;//cout<<xi_proton_left<<endl;
      double t_proton_left = rec_proton_left->t;
      double thetax_proton_left = rec_proton_left->thx;
      double thetay_proton_left = rec_proton_left->thy;
      double beta_proton_left = x_plus/-xi_proton_left;
      bool xi_region_left = -xi_proton_left>=0.03 && -xi_proton_left<=0.1;
      bool t_region_left = -t_proton_left>=0.03 && -t_proton_left<=1;
      bool good_proton_left = proton_left_valid;// && (xi_proton_left < 0.);
      double mass_x = sqrt(xi_proton_right*8000);

      Double_t deltapt = abs(Jet1_pt - Jet2_pt);
      Double_t deltaeta = abs(Jet1_eta - Jet2_eta);
      Double_t deltaphi = abs(Jet1_phi - Jet2_phi);
      Double_t eta_average = 0.5*(Jet1_eta + Jet2_eta);

      if (jet1_selected && jet2_selected){
          histosTH1F["proton_right_xi_jetcut"]->Fill(-xi_proton_right, event_weight );
          histosTH1F["proton_right_t_jetcut"]->Fill(fabs(t_proton_right), event_weight );
          histosTH1F["proton_left_xi_jetcut"]->Fill(-xi_proton_left, event_weight );
          histosTH1F["proton_left_t_jetcut"]->Fill(fabs(t_proton_left), event_weight );
          histosTH1F["proton_right_ptjet1_jetcut"]->Fill(Jet1_pt, event_weight );
          histosTH1F["proton_right_ptjet2_jetcut"]->Fill(Jet2_pt, event_weight );
          histosTH1F["proton_right_etajet1_jetcut"]->Fill(Jet1_eta, event_weight );
          histosTH1F["proton_right_etajet2_jetcut"]->Fill(Jet2_eta, event_weight );
          histosTH1F["proton_right_vtx_zpos_jetcut"]->Fill(primaryVertex.z, event_weight);
      } 
      if (jet1_selected && jet2_selected && good_proton_right && rp_track_accept_right){
          histosTH1F["proton_right_xi_protontag"]->Fill(-xi_proton_right, event_weight );
          histosTH1F["proton_right_t_protontag"]->Fill(fabs(t_proton_right), event_weight );
          histosTH1F["proton_right_ptjet1_protontag"]->Fill(Jet1_pt, event_weight );
          histosTH1F["proton_right_ptjet2_protontag"]->Fill(Jet2_pt, event_weight );
          histosTH1F["proton_right_etajet1_protontag"]->Fill(Jet1_eta, event_weight );
          histosTH1F["proton_right_etajet2_protontag"]->Fill(Jet2_eta, event_weight );
          histosTH1F["proton_right_vtx_zpos_protontag"]->Fill(primaryVertex.z, event_weight);
      } 
      if (jet1_selected && jet2_selected && good_proton_right && rp_track_accept_right && t_region_right && -xi_proton_right<=0.1){
          histosTH1F["proton_right_xi_protonkin"]->Fill(-xi_proton_right, event_weight );
          histosTH1F["proton_right_t_protonkin"]->Fill(fabs(t_proton_right), event_weight );
          histosTH1F["proton_right_ptjet1_protonkin"]->Fill(Jet1_pt, event_weight );
          histosTH1F["proton_right_ptjet2_protonkin"]->Fill(Jet2_pt, event_weight );
          histosTH1F["proton_right_etajet1_protonkin"]->Fill(Jet1_eta, event_weight );
          histosTH1F["proton_right_etajet2_protonkin"]->Fill(Jet2_eta, event_weight );
          histosTH1F["proton_right_vtx_zpos_protonkin"]->Fill(primaryVertex.z, event_weight);
      } 
      if (jet1_selected && jet2_selected && good_proton_right && rp_track_accept_right && t_region_right && -xi_proton_right<=0.1 && xi_minus_Reco+xi_proton_right<0){
          histosTH1F["proton_right_xi_backgcut"]->Fill(-xi_proton_right, event_weight );
          histosTH1F["proton_right_t_backgcut"]->Fill(fabs(t_proton_right), event_weight );
          histosTH1F["proton_right_ptjet1_backgcut"]->Fill(Jet1_pt, event_weight );
          histosTH1F["proton_right_ptjet2_backgcut"]->Fill(Jet2_pt, event_weight );
          histosTH1F["proton_right_etajet1_backgcut"]->Fill(Jet1_eta, event_weight );
          histosTH1F["proton_right_etajet2_backgcut"]->Fill(Jet2_eta, event_weight );
          histosTH1F["proton_right_vtx_zpos_backgcut"]->Fill(primaryVertex.z, event_weight);
      } 

      if(jet1_selected && jet2_selected && good_proton_right ){
          histosTH1F["proton_right_xi_binsasha"]->Fill(-xi_proton_right, event_weight );
      }
 
      if (good_proton_right && rp_track_accept_right && t_region_right && -xi_proton_right<=0.1){ 
          histosTH1F["proton_right_xi_totem"]->Fill(-xi_proton_right, event_weight );
          histosTH1F["proton_right_t_totem"]->Fill( fabs(t_proton_right), event_weight );
      }
      if(jet1_selected && jet2_selected && good_proton_right && rp_track_accept_right){
         ++n_evt_proton_right;
         histosTH1F["proton_right_chi2"]->Fill( chi2_proton_right, event_weight );
         histosTH1F["proton_right_xi"]->Fill(-xi_proton_right, event_weight );
         histosTH1F["xitotem_xicms_rightRPs"]->Fill( xi_minus_Reco+xi_proton_right, event_weight );
         histosTH1F["proton_right_t"]->Fill( fabs(t_proton_right), event_weight );
         histosTH1F["log_x_minus_before"]->Fill( log10(x_minus), event_weight  );
         if (t_region_right){
             histosTH1F["proton_right_xi_kint"]->Fill( -xi_proton_right, event_weight );
             histosTH1F["proton_right_xi_kint_bin25"]->Fill( -xi_proton_right, event_weight );}
         if (t_region_right && xi_minus_Reco+xi_proton_right<0){ 
            histosTH1F["proton_right_xi_sel"]->Fill( -xi_proton_right, event_weight );
            histosTH1F["proton_right_t_sel"]->Fill( fabs(t_proton_right), event_weight );
         }

         if (-xi_proton_right<=0.1 && t_region_right){ 
             histosTH1F["proton_right_xi_kint_kinxi"]->Fill( -xi_proton_right, event_weight );
             histosTH1F["proton_minus_xi_kint"]->Fill( xi_minus_Reco, event_weight );
             histosTH1F["xitotem_xicms_rightRPs_kint_kinxi"]->Fill( xi_minus_Reco+xi_proton_right, event_weight );
             xi_cms_st = xi_minus_Reco;
             xi_totem_st = -xi_proton_right;
             xi_cms_minus_totem_st = xi_minus_Reco+xi_proton_right;
             //small_tree->Fill();
         }
         histosTH1F["proton_right_xi_signal_kinxi"]->Fill( -xi_proton_right, event_weight );
         if (xi_minus_Reco+xi_proton_right<0 ) histosTH1F["proton_right_xi_signal_kinxi_cut"]->Fill( -xi_proton_right, event_weight );
         if (-xi_proton_right<=0.1 && t_region_right && xi_minus_Reco+xi_proton_right<0){
             t_right_totem_sel = fabs(t_proton_right);
             xi_totem_sel = -xi_proton_right;
             small_tree->Fill();
             ++n_evt_proton_right_cut;
             histosTH1F["proton_minus_xi_kint_sel"]->Fill( xi_minus_Reco, event_weight );
             histosTH1F["proton_right_thx_signal_kint_kinxi"]->Fill( thetax_proton_right, event_weight );     
             histosTH1F["proton_right_thy_signal_kint_kinxi"]->Fill( thetay_proton_right, event_weight );     
             histosTH1F["proton_right_t_signal_kint_kinxi"]->Fill( fabs(t_proton_right), event_weight );     
             histosTH1F["proton_right_t_signal_kint_kinxi_bin"]->Fill( fabs(t_proton_right), event_weight );     
             histosTH1F["proton_right_t_signal_kint_kinxi_bin3"]->Fill( fabs(t_proton_right), event_weight );     
             histosTH1F["proton_right_xi_signal_kint_kinxi"]->Fill(-xi_proton_right, event_weight );
             histosTH1F["proton_right_xi_signal_kint_kinxi_bin"]->Fill(-xi_proton_right, event_weight );
             histosTH1F["proton_right_beta_kint_kinxi"]->Fill(beta_proton_right, event_weight );
             histosTH1F["log_x_minus_kint_kinxi"]->Fill( log10(x_minus), event_weight  );
             histosTH1F["leadingJet_pt_right_signal_kint_kinxi"]->Fill(Jet1_pt, event_weight  );
             histosTH1F["secondJet_pt_right_signal_kint_kinxi"]->Fill(Jet2_pt, event_weight  );
             histosTH1F["leadingJet_eta_right_signal_kint_kinxi"]->Fill(Jet1_eta, event_weight  );
             histosTH1F["secondJet_eta_right_signal_kint_kinxi"]->Fill(Jet2_eta, event_weight  );
             histosTH1F["EtaJet_average_right_signal_kint_kinxi"]->Fill( eta_average, event_weight  );
             histosTH1F["DeltaPtJet_right_signal_kint_kinxi"]->Fill( deltapt, event_weight  );
             histosTH1F["DeltaEtaJet_right_signal_kint_kinxi"]->Fill( deltaeta, event_weight  );
             histosTH1F["DeltaPhiJet_right_signal_kint_kinxi"]->Fill( deltaphi, event_weight  );
             histosTH1F["massJet_right_signal_kint_kinxi"]->Fill( mass_jets, event_weight  );
             histosTH1F["mass_x_right_signal_kint_kinxi"]->Fill( mass_x, event_weight  );
         }
         if (xi_region_right && t_region_right) {
             histosTH2F["proton_y_vs_x_rp_124_125_accept"]->Fill( rp_x_124, rp_y_124, event_weight );
             histosTH2F["proton_y_vs_x_rp_124_125_accept"]->Fill( rp_x_125, rp_y_125, event_weight );
             histosTH1F["xitotem_xicms_rightRPs_kin"]->Fill( xi_minus_Reco+xi_proton_right, event_weight );
             histosTH1F["proton_right_xi_kin"]->Fill(-xi_proton_right, event_weight );

             if (xi_minus_Reco>0.12) histosTH1F["proton_right_xi_cut"]->Fill(-xi_proton_right, event_weight );

             if(xi_minus_Reco+xi_proton_right>0){
               histosTH1F["proton_right_xi_halo"]->Fill(-xi_proton_right, event_weight );
               histosTH1F["proton_right_t_halo"]->Fill(fabs(t_proton_right), event_weight );
             }

             histosTH1F["proton_right_logXi"]->Fill( log10(-xi_proton_right), event_weight );
             histosTH1F["pf_xiMinus_minus_proton_right_xi"]->Fill( (xi_minus_Reco + xi_proton_right), event_weight );
             histosTH2F["proton_right_logXi_vs_pf_logXiPlus"]->Fill( log10(xi_plus_Reco),log10(-xi_proton_right), event_weight );
             histosTH2F["proton_right_logXi_vs_pf_logXiMinus"]->Fill( log10(xi_minus_Reco),log10(-xi_proton_right), event_weight );
             histosTH2F["proton_right_xi_vs_pf_xiPlus"]->Fill( xi_plus_Reco, -xi_proton_right, event_weight );
             histosTH2F["proton_right_xi_vs_pf_xiMinus"]->Fill( xi_minus_Reco, -xi_proton_right, event_weight );
             histosTH2F["proton_right_logXi_vs_t"]->Fill( fabs(t_proton_right), log10(-xi_proton_right), event_weight );
             histosTH1F["proton_right_t_cut"]->Fill( fabs(t_proton_right), event_weight );
             histosTH1F["proton_right_logXi"]->Fill( log10(-xi_proton_right), event_weight );

             if (xi_minus_Reco+xi_proton_right<0){
                histosTH1F["proton_right_t_signal"]->Fill( fabs(t_proton_right), event_weight );
                histosTH1F["proton_right_xi_signal"]->Fill(-xi_proton_right, event_weight );
                histosTH1F["proton_right_beta"]->Fill(beta_proton_right, event_weight );
                histosTH1F["log_x_minus"]->Fill( log10(x_minus), event_weight  );
                histosTH1F["leadingJet_pt_right_signal"]->Fill(Jet1_pt, event_weight  );
                histosTH1F["secondJet_pt_right_signal"]->Fill(Jet2_pt, event_weight  );
                histosTH1F["EtaJet_average_right_signal"]->Fill( eta_average, event_weight  );                    
                histosTH1F["DeltaPtJet_right_signal"]->Fill( deltapt, event_weight  );
                histosTH1F["DeltaEtaJet_right_signal"]->Fill( deltaeta, event_weight  );
                histosTH1F["DeltaPhiJet_right_signal"]->Fill( deltaphi, event_weight  );
             } 
         
             if( pfJet_coll->size() > 0 ){
               histosTH2F["proton_right_t_vs_leadingJet_pt"]->Fill( Jet1_pt, fabs(t_proton_right), event_weight );
             }
         }
      }

      if (jet1_selected && jet2_selected && good_proton_left && rp_track_accept_left){
          histosTH1F["proton_left_xi_protontag"]->Fill(-xi_proton_left, event_weight );
          histosTH1F["proton_left_t_protontag"]->Fill(fabs(t_proton_left), event_weight );
          histosTH1F["proton_left_ptjet1_protontag"]->Fill(Jet1_pt, event_weight );
          histosTH1F["proton_left_ptjet2_protontag"]->Fill(Jet2_pt, event_weight );
          histosTH1F["proton_left_etajet1_protontag"]->Fill(Jet1_eta, event_weight );
          histosTH1F["proton_left_etajet2_protontag"]->Fill(Jet2_eta, event_weight );
          histosTH1F["proton_left_vtx_zpos_protontag"]->Fill(primaryVertex.z, event_weight);
      } 
      if (jet1_selected && jet2_selected && good_proton_left && rp_track_accept_left && t_region_left && -xi_proton_left<=0.1){
          histosTH1F["proton_left_xi_protonkin"]->Fill(-xi_proton_left, event_weight );
          histosTH1F["proton_left_t_protonkin"]->Fill(fabs(t_proton_left), event_weight );
          histosTH1F["proton_left_ptjet1_protonkin"]->Fill(Jet1_pt, event_weight );
          histosTH1F["proton_left_ptjet2_protonkin"]->Fill(Jet2_pt, event_weight );
          histosTH1F["proton_left_etajet1_protonkin"]->Fill(Jet1_eta, event_weight );
          histosTH1F["proton_left_etajet2_protonkin"]->Fill(Jet2_eta, event_weight );
          histosTH1F["proton_left_vtx_zpos_protonkin"]->Fill(primaryVertex.z, event_weight);
      } 
      if (jet1_selected && jet2_selected && good_proton_left && rp_track_accept_left && t_region_left && -xi_proton_left<=0.1 && xi_plus_Reco+xi_proton_left<0){
          histosTH1F["proton_left_xi_backgcut"]->Fill(-xi_proton_left, event_weight );
          histosTH1F["proton_left_t_backgcut"]->Fill(fabs(t_proton_left), event_weight );
          histosTH1F["proton_left_ptjet1_backgcut"]->Fill(Jet1_pt, event_weight );
          histosTH1F["proton_left_ptjet2_backgcut"]->Fill(Jet2_pt, event_weight );
          histosTH1F["proton_left_etajet1_backgcut"]->Fill(Jet1_eta, event_weight );
          histosTH1F["proton_left_etajet2_backgcut"]->Fill(Jet2_eta, event_weight );
          histosTH1F["proton_left_vtx_zpos_backgcut"]->Fill(primaryVertex.z, event_weight);
      } 

      if(jet1_selected && jet2_selected && good_proton_left && rp_track_accept_left){
         ++n_evt_proton_left;
         histosTH1F["proton_left_chi2"]->Fill( chi2_proton_left, event_weight );
         histosTH1F["proton_left_xi"]->Fill(-xi_proton_left, event_weight );
         histosTH1F["xitotem_xicms_leftRPs"]->Fill( xi_plus_Reco+xi_proton_left, event_weight );
         histosTH1F["proton_left_t"]->Fill( fabs(t_proton_left), event_weight );
         histosTH1F["log_x_plus_before"]->Fill( log10(x_plus), event_weight  );
         if (t_region_left){
            histosTH1F["proton_left_xi_kint"]->Fill( -xi_proton_left, event_weight );
            histosTH1F["proton_left_xi_kint_bin25"]->Fill( -xi_proton_left, event_weight );}
         if (t_region_left && xi_plus_Reco+xi_proton_left<0){
            histosTH1F["proton_left_xi_sel"]->Fill( -xi_proton_left, event_weight );
            histosTH1F["proton_left_t_sel"]->Fill( fabs(t_proton_left), event_weight );
         }

         if (-xi_proton_left<=0.1 && t_region_left){
             histosTH1F["xitotem_xicms_leftRPs_kint_kinxi"]->Fill( xi_plus_Reco+xi_proton_left, event_weight );
             histosTH1F["proton_left_xi_kint_kinxi"]->Fill( -xi_proton_left, event_weight );
             histosTH1F["proton_plus_xi_kint_kinxi"]->Fill( xi_plus_Reco, event_weight );
         }
         if (-xi_proton_left<=0.1 && t_region_left && xi_plus_Reco+xi_proton_left<0){
             ++n_evt_proton_left_cut;
             t_left_totem_sel = fabs(t_proton_left);
             small_tree->Fill();
             histosTH1F["proton_plus_xi_kint_kinxi_sel"]->Fill( xi_plus_Reco, event_weight );
             histosTH1F["proton_left_t_signal_kint_kinxi"]->Fill( fabs(t_proton_left), event_weight );
             histosTH1F["proton_left_t_signal_kint_kinxi_bin"]->Fill( fabs(t_proton_left), event_weight );
             histosTH1F["proton_left_t_signal_kint_kinxi_bin3"]->Fill( fabs(t_proton_left), event_weight );     
             histosTH1F["proton_left_t_signal_kint_kinxi_bin4"]->Fill( fabs(t_proton_left), event_weight );     
             histosTH1F["proton_left_xi_signal_kint_kinxi"]->Fill(-xi_proton_left, event_weight );
             histosTH1F["proton_left_xi_signal_kint_kinxi_bin"]->Fill(-xi_proton_left, event_weight );
             histosTH1F["proton_left_beta_kint_kinxi"]->Fill(beta_proton_left, event_weight );
             histosTH1F["log_x_plus_kint_kinxi"]->Fill( log10(x_plus), event_weight  );
             histosTH1F["leadingJet_pt_left_signal_kint_kinxi"]->Fill(Jet1_pt, event_weight  );
             histosTH1F["secondJet_pt_left_signal_kint_kinxi"]->Fill(Jet2_pt, event_weight  );
             histosTH1F["EtaJet_average_left_signal_kint_kinxi"]->Fill( eta_average, event_weight  );
             histosTH1F["DeltaPtJet_left_signal_kint_kinxi"]->Fill( deltapt, event_weight  );
             histosTH1F["DeltaEtaJet_left_signal_kint_kinxi"]->Fill( deltaeta, event_weight  );
             histosTH1F["DeltaPhiJet_left_signal_kint_kinxi"]->Fill( deltaphi, event_weight  );
         }
         if (xi_region_left && t_region_left) {
             histosTH2F["proton_y_vs_x_rp_024_025_accept"]->Fill( rp_x_024, rp_y_024, event_weight );
             histosTH2F["proton_y_vs_x_rp_024_025_accept"]->Fill( rp_x_025, rp_y_025, event_weight );
             histosTH1F["xitotem_xicms_leftRPs_kin"]->Fill( xi_plus_Reco+xi_proton_left, event_weight );
             histosTH1F["proton_left_xi_kin"]->Fill(-xi_proton_left, event_weight );

             if (xi_plus_Reco>0.12) histosTH1F["proton_left_xi_cut"]->Fill(-xi_proton_left, event_weight );

             if(xi_plus_Reco+xi_proton_left>0){
               histosTH1F["proton_left_xi_halo"]->Fill(-xi_proton_left, event_weight );
               histosTH1F["proton_left_t_halo"]->Fill(fabs(t_proton_left), event_weight );
             }

             histosTH1F["pf_xiMinus_plus_proton_left_xi"]->Fill( (xi_plus_Reco + xi_proton_left), event_weight );
             histosTH2F["proton_left_logXi_vs_pf_logXiPlus"]->Fill( log10(xi_plus_Reco),log10(-xi_proton_left), event_weight );
             histosTH2F["proton_left_logXi_vs_pf_logXiMinus"]->Fill( log10(xi_plus_Reco),log10(-xi_proton_left), event_weight );
             histosTH2F["proton_left_xi_vs_pf_xiPlus"]->Fill( xi_plus_Reco, -xi_proton_left, event_weight );
             histosTH2F["proton_left_xi_vs_pf_xiMinus"]->Fill( xi_plus_Reco, -xi_proton_left, event_weight );
             histosTH2F["proton_left_logXi_vs_t"]->Fill( fabs(t_proton_left), log10(-xi_proton_left), event_weight );
             histosTH1F["proton_left_t_cut"]->Fill( fabs(t_proton_left), event_weight );
             histosTH1F["proton_left_logXi"]->Fill( log10(-xi_proton_left), event_weight );

             if (xi_plus_Reco+xi_proton_left<0){
                histosTH1F["proton_left_t_signal"]->Fill( fabs(t_proton_left), event_weight );
                histosTH1F["proton_left_xi_signal"]->Fill(-xi_proton_left, event_weight );
                histosTH1F["proton_left_beta"]->Fill(beta_proton_left, event_weight );
                histosTH1F["log_x_plus"]->Fill( log10(x_plus), event_weight  );
                histosTH1F["leadingJet_pt_left_signal"]->Fill(Jet1_pt, event_weight  );
                histosTH1F["secondJet_pt_left_signal"]->Fill(Jet2_pt, event_weight  );
                histosTH1F["EtaJet_average_left_signal"]->Fill( eta_average, event_weight  );
                histosTH1F["DeltaPtJet_left_signal"]->Fill( deltapt, event_weight  );
                histosTH1F["DeltaEtaJet_left_signal"]->Fill( deltaeta, event_weight  );
                histosTH1F["DeltaPhiJet_left_signal"]->Fill( deltaphi, event_weight  );
             }

             if( pfJet_coll->size() > 0 ){
               histosTH2F["proton_left_t_vs_leadingJet_pt"]->Fill( Jet1_pt, fabs(t_proton_left), event_weight );
             }
         }
      }
 

    }//end loop for events
    file->Close();
//ofs.close();
  
  }//end of loop over files

  //output file
  TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  output->cd();


  event_selection->SetBinContent(1, n_evt_total);
  event_selection->SetBinContent(2, n_evt_trigger);
  event_selection->SetBinContent(3, n_evt_vtx);
  event_selection->SetBinContent(4, n_evt_jets);
  event_selection->SetBinContent(5, n_evt_etapf);
  event_selection->SetBinContent(6, n_evt_pass_threshold);
  event_selection->SetBinContent(7, n_evt_proton_right);
  event_selection->SetBinContent(8, n_evt_proton_right_cut);
  event_selection->SetBinContent(9, n_evt_proton_left);
  event_selection->SetBinContent(10, n_evt_proton_left_cut);
  event_selection->GetXaxis()->SetBinLabel(1, "total");
  event_selection->GetXaxis()->SetBinLabel(2, "trigger");
  event_selection->GetXaxis()->SetBinLabel(3, "vtx");
  event_selection->GetXaxis()->SetBinLabel(4, "jets");
  event_selection->GetXaxis()->SetBinLabel(5, "etapf");
  event_selection->GetXaxis()->SetBinLabel(6, "threshold");
  event_selection->GetXaxis()->SetBinLabel(7, "proton_right");
  event_selection->GetXaxis()->SetBinLabel(8, "proton_right_cut");
  event_selection->GetXaxis()->SetBinLabel(9, "proton_left");
  event_selection->GetXaxis()->SetBinLabel(10, "proton_left_cut");
  event_selection->Write();
  small_tree->Write();

  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();
                                  it_histo != histosTH1F.end(); ++it_histo)
     (*it_histo).second->Write();
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
                                  it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();

 
  output->Close();

}
