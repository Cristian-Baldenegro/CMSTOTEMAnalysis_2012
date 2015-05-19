//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
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

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

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
#include "MyCaloTower.h"
#include "MyCaloJet.h"
#include "MyPFJet.h"
#include "MyPUSumInfo.h"
#include "MyPart.h"
#include "MySimVertex.h"
#include "MyTracks.h"
#include "MyVertex.h"
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpTrackInfo.h"

#include "rp_aperture_config.h"
#include "analysis_tools.h"
#include "beam_vtx_smearing.h"


//ROOUNFOLD CLASSES
#include "/storage/lhuertas/uerj-1/CMSTOTEM/mc/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/storage/lhuertas/uerj-1/CMSTOTEM/mc/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/storage/lhuertas/uerj-1/CMSTOTEM/mc/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#define PI 3.141592653589793
using namespace std;

Double_t beta_fit(Double_t *x, Double_t *par ){
  Double_t result = 0;
  //result = par[0]+par[1]*x[0]+ par[2]*x[0]*x[0];
  //result = par[0]+par[1]*x[0]+ par[2]*pow(x[0],2) + par[3]*pow(x[0],3);
  result = par[0]+par[1]*x[0]+ par[2]*pow(x[0],2) + par[3]*pow(x[0],3) +par[4]*pow(x[0],4);
  //result = par[0]+par[1]*x[0]+ par[2]*pow(x[0],2) + par[3]*pow(x[0],3) +par[4]*pow(x[0],4) +par[5]*pow(x[0],5) +par[6]*pow(x[0],6);
  return result;
}
Double_t fFermiLike(Double_t *x, Double_t *par) {
  Double_t result = 0;
  result = 1.0/(TMath::Exp((par[0]-TMath::Sqrt(x[0]))/par[1]) + 1);
  return result;
}

class MyZeroBiasData {
   public:
      MyZeroBiasData() {}
      ~MyZeroBiasData() {}

      bool   proton_rec_left_valid;
      double proton_rec_left_t;
      double proton_rec_left_xi;
      bool   proton_rec_right_valid;
      double proton_rec_right_t;
      double proton_rec_right_xi;
      bool   vtx_valid;
      double vtx_ndof;
      double vtx_x;
      double vtx_y;
      double vtx_z;
      double leadingJet_pt;
      double leadingJet_eta;
      double secondJet_pt;
      double secondJet_eta;
      bool rp_track_valid_120;
      bool rp_track_valid_121;
      bool rp_track_valid_122;
      bool rp_track_valid_123;
      bool rp_track_valid_124;
      bool rp_track_valid_125;
      bool rp_track_valid_020;
      bool rp_track_valid_021;
      bool rp_track_valid_022;
      bool rp_track_valid_023;
      bool rp_track_valid_024;
      bool rp_track_valid_025;
      double rp_x_024;
      double rp_y_024;
      double rp_x_025;
      double rp_y_025;
      double rp_x_124;
      double rp_y_124;
      double rp_x_125;
      double rp_y_125;
      double xi_cms_plus;
      double xi_cms_minus;
};

void diffractive(string const& mc = "pomwig", bool reggeon = false, bool side_minus = true, bool side_plus = false, bool reweigth = false, const Int_t nevt_max = -1){
  
  TString file_name, side;
  if (side_minus && !side_plus) side = "minus";
  if (!side_minus && side_plus) side = "plus";
  if (side_minus && side_plus) side = "bothsides";
  TString rew = (reweigth) ? "_rew" : "";
  TString regg = (reggeon) ? "_reggeon" : ""; 
  file_name = mc + regg + "_" + side + rew + ".root";
  TString outputFileName = "/storage/lhuertas/uerj-1/CMSTOTEM/mc/Workspace/root_files/" + file_name; 
  cout<<outputFileName<<endl;

  //TFile *data = new TFile("/afs/cern.ch/user/l/lhuertas/data_SDdijet_CMSTOTEM_vtx.root");
  TFile *data = new TFile("/afs/cern.ch/user/l/lhuertas/data_SDdijet_CMSTOTEM.root");
  //TFile *data = new TFile("/storage/lhuertas/uerj-1/CMSTOTEM/data/root_files/data_SDdijet_CMSTOTEM.root");

  bool verbose = false;
  string treeName = "evt";//"cms_totem";
  string jetCollName = "ak5PFJets";
//  string jetCorrName = "ak5PFL2L3Residual";
  //string jetCorrName = "ak5PFL2L3"; 
  string jetCorrName = "ak5PFJets"; 

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


  //beta fit - reweigth
  TF1* func_right = new TF1("func_right", beta_fit, 0., 1., 5);
  //<s>=0.10-test      
  func_right->SetParameter(0, 2.02798);
  func_right->SetParameter(1, -16.7293);
  func_right->SetParameter(2, 82.8748);
  func_right->SetParameter(3, -164.542);
  func_right->SetParameter(4, 114.142);

  TF1* func_left = new TF1("func_left", beta_fit, 0., 1., 5);
  //<s>=0.10-test      
  func_left->SetParameter(0, 1.45284);
  func_left->SetParameter(1, -9.0417);
  func_left->SetParameter(2, 32.8926);
  func_left->SetParameter(3, -46.2235);
  func_left->SetParameter(4, 22.5018);

  //beta fit - reweigth - pythia8
  TF1* func_right_pythia = new TF1("func_right_pythia", beta_fit, 0., 1., 5);
  //<s>=0.10-test      
  func_right_pythia->SetParameter(0, 1.5213);
  func_right_pythia->SetParameter(1, -10.2655);
  func_right_pythia->SetParameter(2, 53.7009);
  func_right_pythia->SetParameter(3, -117.164);
  func_right_pythia->SetParameter(4, 89.4077);

  TF1* func_left_pythia = new TF1("func_left_pythia", beta_fit, 0., 1., 5);
  //<s>=0.10-test      
  func_left_pythia->SetParameter(0, 1.15756);
  func_left_pythia->SetParameter(1, -4.80457);
  func_left_pythia->SetParameter(2, 13.8974);
  func_left_pythia->SetParameter(3, -15.0937);
  func_left_pythia->SetParameter(4, 5.66318);

  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;

  // Declaration of histograms
  map<string,TH1F*> histosTH1F;
//   histosTH1F["decisionPhysTrig"] = new TH1F("decisionPhysTrig", "decisionPhysTrig" , 128 , 0 , 128);
//   histosTH1F["decisionTechTrig"] = new TH1F("decisionTechTrig", "decisionTechTrig" , 64 , 0 , 64);
  float xbins = 50;
  float bin2 = 15;

  double weight_st, xi_cms_st, xi_totem_st, xi_totem_sel, xi_cms_minus_totem_st;
  TTree* small_tree = new TTree("small_tree","");
  small_tree->Branch("weight",&weight_st,"weight/D");
  //small_tree->Branch("xi_cms",&xi_cms_st,"xi_cms/D");
  //small_tree->Branch("xi_totem",&xi_totem_st,"xi_totem/D");
  small_tree->Branch("xi_totem_sel",&xi_totem_sel,"xi_totem_sel/D");
  //small_tree->Branch("xi_cms_minus_totem",&xi_cms_minus_totem_st,"xi_cms_minus_totem/D");

  histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
  histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);

  //histosTH1F["pt_gen"] = new TH1F("pt_gen" , "pt_gen;pt;nTracks" , 120 , 0 , 6);
  histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 150 , 0. , 15.);
  histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 200 , -5.2 , 5.2);
  histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 200 , -M_PI , M_PI);
  histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);
  
  histosTH1F["jet_pt"] = new TH1F("jet_pt", "p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["jet_eta"] = new TH1F("jet_eta", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["jet_phi"] = new TH1F("jet_phi", "#phi(jet)" , 20 , -M_PI , M_PI);

  histosTH1F["leadingJet_pt_rec_signal_kin_cut"] = new TH1F("leadingJet_pt_rec_signal_kin_cut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_gen_signal_kin_cut_recsel"] = new TH1F("leadingJet_pt_gen_signal_kin_cut_recsel", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_phi_rec_signal_kin_cut"] = new TH1F("leadingJet_phi_rec_signal_kin_cut", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["leadingJet_phi_gen_signal_kin_cut_recsel"] = new TH1F("leadingJet_phi_rec_signal_kin_cut_recsel", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["leadingJet_eta_rec_signal_kin_cut"] = new TH1F("leadingJet_eta_rec_signal_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_eta_gen_signal_kin_cut_recsel"] = new TH1F("leadingJet_eta_rec_signal_kin_cut_recsel", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_pt_rec_signal_kint_kinxi_cut"] = new TH1F("leadingJet_pt_rec_signal_kint_kinxi_cut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_rec_signal_kint_kinxi_cut_triggereff"] = new TH1F("leadingJet_pt_rec_signal_kint_kinxi_cut_triggereff", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_rec_signal_kint_kinxi_cut_protoneff"] = new TH1F("leadingJet_pt_rec_signal_kint_kinxi_cut_protoneff", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_rec_signal_kint_kinxi_cut_zvtx"] = new TH1F("leadingJet_pt_rec_signal_kint_kinxi_cut_zvtx", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_eta_rec_signal_kint_kinxi_cut"] = new TH1F("leadingJet_eta_rec_signal_kint_kinxi_cut", "eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["leadingJet_pt_gen_signal_kint_kinxi_cut_recsel"] = new TH1F("leadingJet_pt_gen_signal_kint_kinxi_cut_recsel", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_gen_signal_kint_kinxi_cut"] = new TH1F("leadingJet_pt_gen_signal_kint_kinxi_cut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_eta_gen_signal_kint_kinxi_cut"] = new TH1F("leadingJet_eta_gen_signal_kint_kinxi_cut", "eta(jet)" , bin2 , -5.2 , 5.2);

  histosTH1F["secondJet_eta_rec_signal_kin_cut"] = new TH1F("secondJet_eta_rec_signal_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_eta_gen_signal_kin_cut_recsel"] = new TH1F("secondJet_eta_gen_signal_kin_cut_recsel", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_phi_rec_signal_kin_cut"] = new TH1F("secondJet_phi_rec_signal_kin_cut", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["secondJet_phi_gen_signal_kin_cut_recsel"] = new TH1F("secondJet_phi_gen_signal_kin_cut_recsel", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["secondJet_pt_rec_signal_kin_cut"] = new TH1F("secondJet_pt_rec_signal_kin_cut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_gen_signal_kin_cut_recsel"] = new TH1F("secondJet_pt_gen_signal_kin_cut_recsel", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_eta_signal_kin_cut"] = new TH1F("secondJet_eta_signal_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_pt_rec_signal_kint_kinxi_cut"] = new TH1F("secondJet_pt_rec_signal_kint_kinxi_cut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_rec_signal_kint_kinxi_cut_triggereff"] = new TH1F("secondJet_pt_rec_signal_kint_kinxi_cut_triggereff", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_rec_signal_kint_kinxi_cut_protoneff"] = new TH1F("secondJet_pt_rec_signal_kint_kinxi_cut_protoneff", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_rec_signal_kint_kinxi_cut_zvtx"] = new TH1F("secondJet_pt_rec_signal_kint_kinxi_cut_zvtx", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_eta_rec_signal_kint_kinxi_cut"] = new TH1F("secondJet_eta_rec_signal_kint_kinxi_cut", "eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["secondJet_pt_gen_signal_kint_kinxi_cut_recsel"] = new TH1F("secondJet_pt_gen_signal_kint_kinxi_cut_recsel", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_gen_signal_kint_kinxi_cut"] = new TH1F("secondJet_pt_gen_signal_kint_kinxi_cut", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_eta_gen_signal_kint_kinxi_cut"] = new TH1F("secondJet_eta_gen_signal_kint_kinxi_cut", "eta(jet)" , bin2 , -5.2 , 5.2);

  histosTH1F["DeltaPtJet_rec_signal_kin_cut"] = new TH1F("Delta_pt_Jet_rec_signal_kin_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_rec_signal_kin_cut"] = new TH1F("Delta_eta_Jet_rec_signal_kin_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_rec_signal_kin_cut"] = new TH1F("Delta_phi_Jet_rec_signal_kin_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_rec_signal_kin_cut"] = new TH1F("mass_Jet_rec_signal_kin_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["DeltaPtJet_gen_signal_kin_cut_recsel"] = new TH1F("Delta_pt_Jet_gen_signal_kin_cut_recsel", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_gen_signal_kin_cut_recsel"] = new TH1F("Delta_eta_Jet_gen_signal_kin_cut_recsel", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_gen_signal_kin_cut_recsel"] = new TH1F("Delta_phi_Jet_gen_signal_kin_cut_recsel", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_gen_signal_kin_cut_recsel"] = new TH1F("mass_Jet_gen_signal_kin_cut_recsel", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["EtaJet_average_rec_signal_kin_cut"] = new TH1F("eta_Jet_average_rec_signal_kin_cut", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_gen_signal_kin_cut_recsel"] = new TH1F("eta_Jet_average_gen_signal_kin_cut_recsel", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_rec_signal_kint_kinxi_cut"] = new TH1F("eta_Jet_average_rec_signal_kint_kinxi_cut", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_rec_signal_kint_kinxi_cut_triggereff"] = new TH1F("eta_Jet_average_rec_signal_kint_kinxi_cut_triggereff", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_rec_signal_kint_kinxi_cut_protoneff"] = new TH1F("eta_Jet_average_rec_signal_kint_kinxi_cut_protoneff", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_rec_signal_kint_kinxi_cut_zvtx"] = new TH1F("eta_Jet_average_rec_signal_kint_kinxi_cut_zvtx", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_gen_signal_kint_kinxi_cut_recsel"] = new TH1F("eta_Jet_average_gen_signal_kint_kinxi_cut_recsel", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_gen_signal_kint_kinxi_cut"] = new TH1F("eta_Jet_average_gen_signal_kint_kinxi_cut", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["massJet_rec_signal_kint_kinxi_cut"] = new TH1F("mass_Jet_rec_signal_kint_kinxi_cut", "Mass(jet)" , bin2 , 0 , 450);
  histosTH1F["massJet_gen_signal_kint_kinxi_cut"] = new TH1F("mass_Jet_gen_signal_kint_kinxi_cut", "Mass(jet)" , bin2 , 0 , 450);
  histosTH1F["mass_x_rec_signal_kint_kinxi_cut"] = new TH1F("mass_x_rec_signal_kint_kinxi_cut", "M_{X}" , bin2 , 0 , 450);
  histosTH1F["mass_x_gen_signal_kint_kinxi_cut"] = new TH1F("mass_x_gen_signal_kint_kinxi_cut", "M_{X}" , bin2 , 0 , 450);

  histosTH1F["Eta_max"] = new TH1F("Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Eta_min"] = new TH1F("Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 20 , 0,11);
  histosTH1F["xi_plus_Reco"] = new TH1F("xi+", "#xi^{+}" , 20 , 0,4);
  histosTH1F["xi_minus_Reco"] = new TH1F("xi-", "#xi^{-}" , 20 , 0,4);
  histosTH1F["logxi_plus"] = new TH1F("logxi+", "Log #xi^{+}" , 20 , -3,0.5);
  histosTH1F["logxi_plus_gen"] = new TH1F("logxi+_gen", "Log #xi_{+}^{gen}" , 20 , -3,0.5);
  histosTH1F["logxi_minus_gen"] = new TH1F("logxi-_gen", "Log #xi_{-}^{gen}" , 20 , -3,0.5);
  histosTH1F["correction"] = new TH1F("correction", "Correction factor" , 20 , 0,2);
  histosTH1F["resolution_after"] = new TH1F("resolution_after", "Resolution" , 20 , -2,2);
  histosTH1F["resolution_before"] = new TH1F("resolution_before", "Resolution" , 20 , -2,2);
  histosTH1F["log_x_minus"] = new TH1F("log_x_minus", "Log x^{-}" , 20 , -4, 0);
  histosTH1F["log_x_minus_accepted"] = new TH1F("log_x_minus_accepted", "Log x^{-}" , 20 , -4, 0);

  histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);

//  Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.16, 0.20, 0.25, 0.30, 0.40, 0.50, 0.65, 1.};
  Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.16, 0.20, 0.25, 0.32, 0.42, 0.52, 0.65, 1.};
  Float_t tbins_2[11] = {0.03, 0.07, 0.11, 0.15, 0.20, 0.25, 0.32, 0.42, 0.52, 0.65, 1.};
  Float_t tbins_3[10] = {0.03, 0.07, 0.11, 0.18, 0.24, 0.32,  0.43, 0.53, 0.65, 1.};
  Float_t tbins_4[9] = {0.03, 0.08, 0.13, 0.21, 0.31, 0.41,  0.55, 0.75, 1.};
  float xi_bins[10] = {-0.05, 0, 0.015, 0.030, 0.047, 0.064, 0.08, 0.1, 0.12, 0.2};
  float xi_bins_2[8] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.2};
  //float bin_sasha[4] = {0.0003, 0.002, 0.0045, 0.01};
  float bin_sasha[9] ={0.0003, 0.002, 0.0045, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1};
  histosTH1F["xi_rec_proton_signal"] = new TH1F("xi_rec_proton_signal", "xi_proton" , xbins, -0.05, 0.2);
  histosTH1F["xi_gen_proton_signal_recsel"] = new TH1F("xi_gen_proton_signal_recsel", "xi_proton" , xbins, -0.05, 0.2);
  histosTH1F["xi_rec_proton_signal_kint"] = new TH1F("xi_rec_proton_signal_kint", "xi_proton" , xbins, -0.05, 0.2);
  histosTH1F["xi_rec_proton_signal_kinxi_kint"] = new TH1F("xi_rec_proton_signal_kinxi_kint", "xi_proton" , 25, -0.05, 0.2);
  histosTH1F["xi_cms_rec_proton_signal_kint"] = new TH1F("xi_cms_rec_proton_signal_kint", "xi_proton" , xbins, -0.05, 0.2);
  histosTH1F["xi_cms_rec_proton_signal_kint_cut"] = new TH1F("xi_cms_rec_proton_signal_kint_cut", "xi_proton" , xbins, -0.05, 0.2);
  histosTH1F["xi_rec_proton_signal_kint_bin25"] = new TH1F("xi_rec_proton_signal_kint_bin25", "xi_proton" , 25, -0.05, 0.2);
  histosTH1F["xi_gen_proton_signal_kint_recsel"] = new TH1F("xi_gen_proton_signal_kint_recsel", "xi_proton" , xbins, -0.05, 0.2);
  histosTH1F["xi_rec_proton_signal_kint_cut"] = new TH1F("xi_rec_proton_signal_kint_cut", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_gen_proton_signal_kint_cut_recsel"] = new TH1F("xi_gen_proton_signal_kint_cut_recsel", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_rec_proton_signal_kint_kinxi_cut"] = new TH1F("xi_rec_proton_signal_kint_kinxi_cut", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_rec_proton_signal_kint_kinxi_cut_bin15"] = new TH1F("xi_rec_proton_signal_kint_kinxi_cut_bin15", "xi_proton" , 15, -0.05, 0.2);
  histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin"] = new TH1F("xi_rec_proton_minus_signal_kint_kinxi_cut_bin", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_triggereff"] = new TH1F("xi_rec_proton_minus_signal_kint_kinxi_cut_bin_triggereff", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_protoneff"] = new TH1F("xi_rec_proton_minus_signal_kint_kinxi_cut_bin_protoneff", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_zvtx"] = new TH1F("xi_rec_proton_minus_signal_kint_kinxi_cut_bin_zvtx", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_smear"] = new TH1F("xi_rec_proton_minus_signal_kint_kinxi_cut_bin_smear", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_fidcut"] = new TH1F("xi_rec_proton_minus_signal_kint_kinxi_cut_bin_fidcut", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_slopeup"] = new TH1F("xi_rec_proton_minus_signal_kint_kinxi_cut_bin_slopeup", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_slopedw"] = new TH1F("xi_rec_proton_minus_signal_kint_kinxi_cut_bin_slopedw", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_gen_proton_minus_signal_kint_kinxi_cut_bin_slopeup"] = new TH1F("xi_gen_proton_minus_signal_kint_kinxi_cut_bin_slopeup", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_gen_proton_minus_signal_kint_kinxi_cut_bin_slopedw"] = new TH1F("xi_gen_proton_minus_signal_kint_kinxi_cut_bin_slopedw", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_plus_signal_kint_kinxi_cut_bin"] = new TH1F("xi_rec_proton_plus_signal_kint_kinxi_cut_bin", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_gen_proton_signal_kint_kinxi_cut_recsel"] = new TH1F("xi_gen_proton_signal_kint_kinxi_cut_recsel", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_gen_proton_signal_kint_kinxi_cut"] = new TH1F("xi_gen_proton_signal_kint_kinxi_cut", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_gen_proton_minus_signal_kint_kinxi_cut_bin"] = new TH1F("xi_gen_proton_minus_signal_kint_kinxi_cut_bin", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_gen_proton_plus_signal_kint_kinxi_cut_bin"] = new TH1F("xi_gen_proton_plus_signal_kint_kinxi_cut_bin", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_rec_proton_kint_kinxi_bin_nojet"] = new TH1F("xi_rec_proton_kint_kinxi_bin_nojet", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_gen_proton_kint_kinxi_bin_nojet"] = new TH1F("xi_gen_proton_kint_kinxi_bin_nojet", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_cms_rec_pt20"] = new TH1F("xi_cms_rec_pt20", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_cms_gen_pt20"] = new TH1F("xi_cms_gen_pt20", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_gen_pt20_proton_signal_kinxi_cut"] = new TH1F("xi_gen_pt20_proton_signal_kinxi_cut", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_gen_pt20_proton_signal_kinxi_cut_bin"] = new TH1F("xi_gen_pt20_proton_signal_kinxi_cut_bin", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_rec_proton_signal_kin_cut"] = new TH1F("xi_rec_proton_signal_kin_cut", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_gen_proton_signal_kin_cut_recsel"] = new TH1F("xi_gen_proton_signal_kin_cut_recsel", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_gen_proton_signal_kin_cut"] = new TH1F("xi_gen_proton_signal_kin_cut", "xi_proton" , 9, xi_bins);
  histosTH1F["xi_cms_rec"] = new TH1F("xi_cms_rec", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_cms_gen"] = new TH1F("xi_cms_gen", "xi_proton" , 7, xi_bins_2);
  histosTH1F["xi_cms_rec_bin"] = new TH1F("xi_cms_rec_bin", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_cms_rec_gap_bin"] = new TH1F("xi_cms_rec_gap_bin", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_rec_proton_minus_binsasha"] = new TH1F("xi_rec_proton_minus_binsasha", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_gen_proton_minus_binsasha"] = new TH1F("xi_gen_proton_minus_binsasha", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_cms_gen_bin"] = new TH1F("xi_cms_gen_bin", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_cms_gen_gap_bin"] = new TH1F("xi_cms_gen_gap_bin", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_rec_proton_signal_kinxi_cut"] = new TH1F("xi_rec_proton_signal_kinxi_cut", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_rec_proton_signal_kinxi"] = new TH1F("xi_rec_proton_signal_kinxi", "xi_proton" , 8, bin_sasha);
  histosTH1F["xi_gen_proton_signal_kinxi"] = new TH1F("xi_gen_proton_signal_kinxi", "xi_proton" , 8, bin_sasha);

  histosTH1F["mass_gen_proton_signal"] = new TH1F("mass_gen_proton_signal","mass", 50, -0.00001, 0.00001);
  histosTH1F["energy_gen_proton_signal"] = new TH1F("energy_gen_proton_signal","energy", 50, 2000,5000);
  histosTH1F["px_gen_proton_signal"] = new TH1F("px_gen_proton_signal","px", 50, -4000,4000);
  histosTH1F["py_gen_proton_signal"] = new TH1F("py_gen_proton_signal","py", 50, -4000,4000);
  histosTH1F["pt_gen_proton_signal"] = new TH1F("pt_gen_proton_signal","pt", 50, -1,1);
  histosTH1F["pt2_gen_proton_signal"] = new TH1F("pt2_gen_proton_signal","pt2", 50, -1,1);
  
  histosTH1F["t_rec_proton_signal"] = new TH1F("t_rec_proton_signal", "t_proton" , 11, tbins);
  histosTH1F["t_gen_proton_signal_recsel"] = new TH1F("t_gen_proton_signal_recsel", "t_proton" , 11, tbins);
  histosTH1F["t_rec_proton_signal_kin"] = new TH1F("t_rec_proton_signal_kin", "t_proton" , 11, tbins);
  histosTH1F["t_rec_proton_signal_kin_cut"] = new TH1F("t_rec_proton_signal_kin_cut", "t_proton" , 11, tbins);
  histosTH1F["t_gen_proton_signal_kin_cut_recsel"] = new TH1F("t_gen_proton_signal_kin_cut_recsel", "t_proton" , 11, tbins);
  histosTH1F["t_rec_proton_signal_kint_kinxi_cut"] = new TH1F("t_rec_proton_signal_kint_kinxi_cut", "t_proton" , 11, tbins);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_slopeup"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin_slopeup", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_slopedw"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin_slopedw", "t_proton" , 10, tbins_2);
  histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin_slopeup"] = new TH1F("t_gen_proton_minus_signal_kint_kinxi_cut_bin_slopeup", "t_proton" , 10, tbins_2);
  histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin_slopedw"] = new TH1F("t_gen_proton_minus_signal_kint_kinxi_cut_bin_slopedw", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_triggereff"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin_triggereff", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_protoneff"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin_protoneff", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_zvtx"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin_zvtx", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin3"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin3", "t_proton" , 9, tbins_3);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin4"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin4", "t_proton" , 8, tbins_4);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_smear"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_smear", "t_proton" , 11, tbins);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_smear"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin_smear", "t_proton" , 10, tbins_2);
  histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_smear"] = new TH1F("t_gen_proton_minus_signal_kint_kinxi_cut_smear", "t_proton" , 11, tbins);
  histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin_smear"] = new TH1F("t_gen_proton_minus_signal_kint_kinxi_cut_bin_smear", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_fidcut"] = new TH1F("t_rec_proton_minus_signal_kint_kinxi_cut_bin_fidcut", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_plus_signal_kint_kinxi_cut_bin"] = new TH1F("t_rec_proton_plus_signal_kint_kinxi_cut_bin", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_plus_signal_kint_kinxi_cut_bin4"] = new TH1F("t_rec_proton_plus_signal_kint_kinxi_cut_bin4", "t_proton" , 8, tbins_4);
  histosTH1F["t_gen_proton_signal_kint_kinxi_cut_recsel"] = new TH1F("t_gen_proton_signal_kint_kinxi_cut_recsel", "t_proton" , 11, tbins);
  histosTH1F["t_gen_proton_signal_kint_kinxi_cut"] = new TH1F("t_gen_proton_signal_kint_kinxi_cut", "t_proton" , 11, tbins);
  histosTH1F["t_gen_proton_kint_kinxi_bin_nojet"] = new TH1F("t_gen_proton_kint_kinxi_bin_nojet", "t_proton" , 10, tbins_2);
  histosTH1F["t_rec_proton_kint_kinxi_bin_nojet"] = new TH1F("t_rec_proton_kint_kinxi_bin_nojet", "t_proton" , 10, tbins_2);
  histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin"] = new TH1F("t_gen_proton_minus_signal_kint_kinxi_cut_bin", "t_proton" , 10, tbins_2);
  histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin3"] = new TH1F("t_gen_proton_minus_signal_kint_kinxi_cut_bin3", "t_proton" , 9, tbins_3);
  histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin4"] = new TH1F("t_gen_proton_minus_signal_kint_kinxi_cut_bin4", "t_proton" , 8, tbins_4);
  histosTH1F["t_gen_proton_plus_signal_kint_kinxi_cut_bin"] = new TH1F("t_gen_proton_plus_signal_kint_kinxi_cut_bin", "t_proton" , 10, tbins_2);
  histosTH1F["t_gen_proton_plus_signal_kint_kinxi_cut_bin4"] = new TH1F("t_gen_proton_plus_signal_kint_kinxi_cut_bin4", "t_proton" , 8, tbins_4);
  histosTH1F["t_gen_proton_signal_kint_kinxi_cut_bin_recsel"] = new TH1F("t_gen_proton_signal_kint_kinxi_cut_bin_recsel", "t_proton" , xbins, 0, 1);
  histosTH1F["t_rec_proton_signal_acep_kint_kinxi"] = new TH1F("t_rec_proton_signal_acep_kint_kinxi", "t_proton" , 11, tbins);
  histosTH1F["t_rec_proton_signal_acep_kint_kinxi_bin"] = new TH1F("t_rec_proton_signal_acep_kint_kinxi_bin", "t_proton" , xbins, 0, 1);

  histosTH1F["thx_proton_kint_kinxi_cut"] = new TH1F("thx_proton_kint_kinxi_cut", "thx_proton" , 20, -5e-4, 5e-4);
  histosTH1F["thy_proton_kint_kinxi_cut"] = new TH1F("thy_proton_kint_kinxi_cut", "thy_proton" , 20, -5e-4, 5e-4);

  histosTH1F["beta_rec_proton_signal"] = new TH1F("beta_rec_proton_signal", "beta_proton" , xbins, 0, 0.9 );
  histosTH1F["beta_gen_proton_signal_recsel"] = new TH1F("beta_gen_proton_signal_recsel", "beta_proton" , xbins, 0, 0.9 );
  histosTH1F["beta_rec_proton_signal_kin_cut"] = new TH1F("beta_rec_proton_signal_kin_cut", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_gen_proton_signal_kin_cut_recsel"] = new TH1F("beta_gen_proton_signal_kin_cut_recsel", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_rec_proton_minus_signal_kint_kinxi_cut"] = new TH1F("beta_rec_proton_minus_signal_kint_kinxi_cut", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_rec_proton_minus_signal_kint_kinxi_cut_triggereff"] = new TH1F("beta_rec_proton_minus_signal_kint_kinxi_cut_triggereff", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_rec_proton_minus_signal_kint_kinxi_cut_protoneff"] = new TH1F("beta_rec_proton_minus_signal_kint_kinxi_cut_protoneff", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_rec_proton_minus_signal_kint_kinxi_cut_zvtx"] = new TH1F("beta_rec_proton_minus_signal_kint_kinxi_cut_zvtx", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_rec_proton_plus_signal_kint_kinxi_cut"] = new TH1F("beta_rec_proton_plus_signal_kint_kinxi_cut", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_gen_proton_signal_kint_kinxi_cut_recsel"] = new TH1F("beta_gen_proton_signal_kint_kinxi_cut_recsel", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_gen_proton_minus_signal_kint_kinxi_cut"] = new TH1F("beta_gen_proton_minus_signal_kint_kinxi_cut", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["beta_gen_proton_plus_signal_kint_kinxi_cut"] = new TH1F("beta_gen_proton_plus_signal_kint_kinxi_cut", "beta_proton" , bin2, 0, 0.9 );
  histosTH1F["t_rec_gen_minus"] = new TH1F("t_rec_gen_minus", "t" , 50, -1, 1 );
  histosTH1F["xi_rec_gen_minus"] = new TH1F("xi_rec_gen_minus", "xi" , 50, -0.2, 0.2 );
  histosTH1F["t_rec_gen_plus"] = new TH1F("t_rec_gen_plus", "t" , 50, -1, 1 );
  histosTH1F["xi_rec_gen_plus"] = new TH1F("xi_rec_gen_plus", "xi" , 50, -0.2, 0.2 );
  histosTH1F["log_x_parton_rec_signal_kin_cut"] = new TH1F("log_x_parton_rec_signal_kin_cut", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_gen_signal_kin_cut_recsel"] = new TH1F("log_x_parton_gen_signal_kin_cut_recsel", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec"] = new TH1F("log_x_parton_minus_rec", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_rec"] = new TH1F("log_x_parton_plus_rec", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_gen"] = new TH1F("log_x_parton_minus_gen", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_gen"] = new TH1F("log_x_parton_plus_gen", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_gen_recsel"] = new TH1F("log_x_parton_gen_recsel", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut"] = new TH1F("log_x_parton_minus_rec_signal_kint_kinxi_cut", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_triggereff"] = new TH1F("log_x_parton_minus_rec_signal_kint_kinxi_cut_triggereff", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_protoneff"] = new TH1F("log_x_parton_minus_rec_signal_kint_kinxi_cut_protoneff", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_zvtx"] = new TH1F("log_x_parton_minus_rec_signal_kint_kinxi_cut_zvtx", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_smear"] = new TH1F("log_x_parton_minus_rec_signal_kint_kinxi_cut_smear", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_fidcut"] = new TH1F("log_x_parton_minus_rec_signal_kint_kinxi_cut_fidcut", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_slopeup"] = new TH1F("log_x_parton_minus_rec_signal_kint_kinxi_cut_slopeup", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_slopedw"] = new TH1F("log_x_parton_minus_rec_signal_kint_kinxi_cut_slopedw", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_gen_signal_kint_kinxi_cut_slopeup"] = new TH1F("log_x_parton_minus_gen_signal_kint_kinxi_cut_slopeup", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_gen_signal_kint_kinxi_cut_slopedw"] = new TH1F("log_x_parton_minus_gen_signal_kint_kinxi_cut_slopedw", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_rec_signal_kint_kinxi_cut"] = new TH1F("log_x_parton_plus_rec_signal_kint_kinxi_cut", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_gen_signal_kint_kinxi_cut_recsel"] = new TH1F("log_x_parton_gen_signal_kint_kinxi_cut_recsel", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_gen_signal_kint_kinxi_cut"] = new TH1F("log_x_parton_minus_gen_signal_kint_kinxi_cut", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_gen_signal_kint_kinxi_cut"] = new TH1F("log_x_parton_plus_gen_signal_kint_kinxi_cut", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_rec_gen"] = new TH1F("log_x_parton_minus_rec_gen_resol", "log_x_parton_rec_gen", 50, -4, 0);
  histosTH1F["log_x_parton_plus_rec_gen"] = new TH1F("log_x_parton_plus_rec_gen_resol", "log_x_parton_rec_gen", 50, -4, 0);

  Float_t bin[16] = {-0.4, -0.112, -0.096, -0.08, -0.064, -0.048, -0.032, -0.016, 0, 0.048, 0.112, 0.176, 0.24, 0.304, 0.368, 0.4};
  histosTH1F["xi_cms_totem_rec_signal"] = new TH1F("xi_cms_totem_rec_signal", "xi_cms_totem" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_gen_signal_recsel"] = new TH1F("xi_cms_totem_gen_signal_recsel", "xi_cms_totem" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_rec_signal_kin"] = new TH1F("xi_cms_totem_rec_signal_kin", "xi_cms_totem" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_rec_signal_kin_bin"] = new TH1F("xi_cms_totem_rec_signal_kin_bin", "xi_cms_totem" , 50, -0.4, 0.4);
  histosTH1F["xi_cms_totem_rec_signal_kinxi_kint"] = new TH1F("xi_cms_totem_rec_signal_kinxi_kint", "xi_cms_totem" , 15, bin);

  histosTH1F["thx_proton"] = new TH1F("thx_proton", "p_proton_smear", 50, -5e-4, 5e-4);
  histosTH1F["thx_proton_smear"] = new TH1F("thx_proton_smear", "p_proton_smear", 50, -5e-4, 5e-4);
  histosTH1F["thy_proton"] = new TH1F("thy_proton", "p_proton_smear", 50, -5e-4, 5e-4);
  histosTH1F["thy_proton_smear"] = new TH1F("thy_proton_smear", "p_proton_smear", 50, -5e-4, 5e-4);
  histosTH1F["t_proton"] = new TH1F("t_proton", "p_proton_smear", 50, 0, 1);
  histosTH1F["t_proton_smear"] = new TH1F("t_proton_smear", "p_proton_smear", 50, 0, 1);

  map<string,TH2F*> histosTH2F;
  double energyMin = -10.;
  double energyMax = 190.;
  int nBinsEnergy = 1000;
  histosTH2F["energyVsEtaAllTypes"] = new TH2F("energyVsEtaAllTypes","energyVsEtaAllTypes",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaUndefined"] = new TH2F("energyVsEtaUndefined","energyVsEtaUndefined",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaChargedHadron"] = new TH2F("energyVsEtaChargedHadron","energyVsEtaChargedHadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaElectron"] = new TH2F("energyVsEtaElectron","energyVsEtaElectron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaMuon"] = new TH2F("energyVsEtaMuon","energyVsEtaMuon",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaGamma"] = new TH2F("energyVsEtaGamma","energyVsEtaGamma",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaNeutralHadron"] = new TH2F("energyVsEtaNeutralHadron","energyVsEtaNeutralHadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHF"] = new TH2F("energyVsEtaHadronHF","energyVsEtaHadronHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFEcalEnergy"] = new TH2F("energyVsEtaHadronHFEcalEnergy","energyVsEtaHadronHFEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFNoEcalEnergy"] = new TH2F("energyVsEtaHadronHFNoEcalEnergy","energyVsEtaHadronHFNoEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaEGammaHF"] = new TH2F("energyVsEtaEGammaHF","energyVsEtaEGammaHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["xi_plus_reco_gen"] = new TH2F("xi+_reco_gen","xi+",82,0,0.5,82,0,0.5);
  histosTH2F["xi_minus_reco_gen"] = new TH2F("xi-_reco_gen","xi-",82,0,0.5,82,0,0.5);
  histosTH2F["logxi_plus_reco_gen"] = new TH2F("logxi+_reco_gen","xi+",82,-3,0.5,82,-3,0.5);
  histosTH2F["logxi_minus_reco_gen"] = new TH2F("logxi-","xi-",82,-3,0.5,82,-3,0.5);

  histosTH2F["rp_track_pos_y_vs_x_020"] = new TH2F("rp_track_pos_y_vs_x_020", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
  histosTH2F["rp_track_pos_y_vs_x_120"] = new TH2F("rp_track_pos_y_vs_x_120", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);

  histosTH2F["proton_plus_xi_vs_t"] = new TH2F("proton_plus_xi_vs_t","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
  histosTH2F["proton_minus_xi_vs_t"] = new TH2F("proton_minus_xi_vs_t","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

  histosTH2F["proton_plus_xi_vs_t_accepted"] = new TH2F("proton_plus_xi_vs_t_accepted","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
  histosTH2F["proton_minus_xi_vs_t_accepted"] = new TH2F("proton_minus_xi_vs_t_accepted","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

  histosTH2F["pos_y_vs_x_proton_plus_accepted_020"] = new TH2F("pos_y_vs_x_proton_plus_accepted_020", "pos_y_vs_x_proton_plus" , 200, -0.05, 0.05, 200, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_minus_accepted_120"] = new TH2F("pos_y_vs_x_proton_minus_accepted_120", "pos_y_vs_x_proton_minus" , 200, -0.05, 0.05, 200, -0.05, 0.05);

  histosTH2F["pos_y_vs_x_proton_plus_024_025"] = new TH2F("pos_y_vs_x_proton_plus_024_025", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_plus_024_025_accept"] = new TH2F("pos_y_vs_x_proton_plus_024_025_accept", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_minus_124_125"] = new TH2F("pos_y_vs_x_proton_minus_124_125", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_minus_124_125_accept"] = new TH2F("pos_y_vs_x_proton_minus_124_125_accept", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);

  histosTH2F["delta_xi_thx_minus"] = new TH2F("delta_xi_thx_minus", "delta_xi_thx_minus" , 50, -0.0001, 0.0001, 50, -0.04, 0.04);
  histosTH2F["delta_xi_thx_plus"] = new TH2F("delta_xi_thx_plus", "delta_xi_thx_plus" , 50, -0.0001, 0.0001, 50, -0.04, 0.04);
  histosTH2F["log_x_parton_resol"] = new TH2F("log_x_parton_resol", "log_x_parton" , 50, -4, 0, 50, -4, 4);
  histosTH2F["log_x_parton_rec_gen"] = new TH2F("log_x_parton_rec_gen", "log_x_parton_rec_gen", 50, -4, 0, 50, -4, 0);
  histosTH2F["beta_rec_gen"] = new TH2F("beta_rec_gen", "beta_rec_gen", bin2, 0, 0.9, bin2, 0, 0.9);
  histosTH2F["xi_totem_rec_gen"] = new TH2F("xi_totem_rec_gen", "xi_proton" , xbins, -0.05, 0.2, xbins, -0.05, 0.2);

  RooUnfoldResponse response_beta_minus (50, 0, 0.9,"unfolded_beta_minus","unfolded_beta_minus");
  RooUnfoldResponse response_beta_plus (50, 0, 0.9,"unfolded_beta_plus","unfolded_beta_plus");
  //histosTH2F["Response_beta"] = (TH2F*) response_beta.Hresponse();
  RooUnfoldResponse response_x_minus (bin2, -4, 0,"unfolded_x_minus","unfolded_x_minus");
  histosTH2F["response_x_minus"] = (TH2F*) response_x_minus.Hresponse();

  RooUnfoldResponse response_x_minus_test (bin2, -4, 0,"unfolded_x_minus_test","unfolded_x_minus_test");
  histosTH2F["response_x_minus_test"] = (TH2F*) response_x_minus_test.Hresponse();

  RooUnfoldResponse response_x_plus (bin2, -4, 0,"unfolded_x_plus","unfolded_x_plus");
  histosTH2F["response_x_plus"] = (TH2F*) response_x_plus.Hresponse();

  RooUnfoldResponse response_t_minus (histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin4"],histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin4"],"unfolded_t_minus","unfolded_t_minus");
  histosTH2F["response_t_minus"] = (TH2F*) response_t_minus.Hresponse();

  RooUnfoldResponse response_t_minus_gauss (histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_smear"],histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_smear"],"unfolded_t_minus_gauss","unfolded_t_minus_gauss");
  histosTH2F["response_t_minus_gauss"] = (TH2F*) response_t_minus_gauss.Hresponse();

  RooUnfoldResponse response_t_minus_test (histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin4"],histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin4"],"unfolded_t_minus_closuretest","unfolded_t_minus_closuertest");
  histosTH2F["response_t_minus_test"] = (TH2F*) response_t_minus_test.Hresponse();

  RooUnfoldResponse response_t_plus (histosTH1F["t_rec_proton_plus_signal_kint_kinxi_cut_bin4"],histosTH1F["t_gen_proton_plus_signal_kint_kinxi_cut_bin4"],"unfolded_t_plus","unfolded_t_plus");
  histosTH2F["response_t_plus"] = (TH2F*) response_t_plus.Hresponse();

  RooUnfoldResponse response_xi_minus (histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin"],histosTH1F["xi_gen_proton_minus_signal_kint_kinxi_cut_bin"],"unfolded_xi_minus","unfolded_xi_minus");
  histosTH2F["response_xi_minus"] = (TH2F*) response_xi_minus.Hresponse();

  RooUnfoldResponse response_xi_minus_test (histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin"],histosTH1F["xi_gen_proton_minus_signal_kint_kinxi_cut_bin"],"unfolded_xi_minus_test","unfolded_xi_minus_test");
  histosTH2F["response_xi_minus_test"] = (TH2F*) response_xi_minus_test.Hresponse();

  RooUnfoldResponse response_xi_plus (histosTH1F["xi_rec_proton_plus_signal_kint_kinxi_cut_bin"],histosTH1F["xi_gen_proton_plus_signal_kint_kinxi_cut_bin"],"unfolded_xi_plus","unfolded_xi_plus");
  histosTH2F["response_xi_plus"] = (TH2F*) response_xi_plus.Hresponse();

  TH1F* event_selection = new TH1F("event_selection", "event_selection", 9, 0, 9);

  for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it){
      it->second->Sumw2();
      }
  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
  gStyle->SetPalette(1);

  
  //===================
  int i_tot = 0 , nevt_tot = 0;
  
  // MC files
  const char *ext=".root";
  vector<TString>* vdirs = new vector<TString>;
  if (mc == "pomwig"){
     if (side_minus){ 
        if (!reggeon) vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/mc/Pomwig_SDDijetsMinus_8TeV/");
        else vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/mc/Pomwig_Reggeon_SDDijetsMinus_8TeV/");
     }
     if (side_plus){
        if (!reggeon) vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/mc/Pomwig_plus_test/");//SDDijetsPlus_8TeV/");
        else vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/mc/Pomwig_Reggeon_SDDijetsPlus_8TeV/");//SDDijetsPlus_8TeV/");
     }
  }
  if (mc == "pythia8_diff") vdirs->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/mc/Pythia8_SD_DD_Dijets_8Tev_Pt20/");

  vector<TString>* vfiles = new vector<TString>;
  for(vector<TString>::iterator itdirs = vdirs->begin(); itdirs != vdirs->end(); ++itdirs){
      TString& dirname = *itdirs;
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
  delete vdirs;

  //Declaration of tree and its branches variables
  TTree* tree = new TTree(treeName.c_str(),"");
  MyEvtId*           evtId        = NULL;
//   MyL1TrigOld*       l1Trig       = NULL;  
//   MyHLTrig*          hltTrig      = NULL;
  vector<MyGenPart>* genPart      = NULL;
  vector<MyTracks>*  track_coll   = NULL;
  vector<MyVertex>*  vertex_coll  = NULL;
  vector<MyCaloJet>*   caloJet_coll   = NULL;
  vector<MyPFJet>*   pfJet_coll   = NULL;
  vector<MyGenJet>*   genJet_coll   = NULL;
  vector<MyPFCand>*  pFlow_coll   = NULL;
  vector<MyCaloTower>*  caloTowers_coll   = NULL;
  MyGenKin*  genKin   = NULL;
  //=================================================

  //ZeroBias Files
  string treeNameZB = "cms_totem";
  TChain treeZB("cms_totem");
  vector<TString>* vdirs_zb = new vector<TString>;
  //vdirs_zb->push_back("/storage/lhuertas/uerj-1/CMSTotem/samples/data/MergedNtuples/HighBeta/198902-8369_8371-V00-02-00/ZeroBias/");
  vdirs_zb->push_back("/storage/lhuertas/uerj-1/CMSTOTEM/samples/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/ZeroBias/");
  vector<TString>* vfiles_zb = new vector<TString>;
  for(vector<TString>::iterator itdirs_zb = vdirs_zb->begin(); itdirs_zb != vdirs_zb->end(); ++itdirs_zb){
      TString& dirname_zb = *itdirs_zb;
      TSystemDirectory dir(dirname_zb, dirname_zb);
      TList *files = dir.GetListOfFiles();
      if (files) {
         TSystemFile *file;
         TString fname;
         TIter next(files);
         while ((file=(TSystemFile*)next())) {
             fname = file->GetName();
             if (!file->IsDirectory() && fname.EndsWith(ext)) {
                 TString root_file_zb = dirname_zb + string(fname.Data());
                 vfiles_zb->push_back(root_file_zb); cout<<root_file_zb<<endl;
                 treeZB.Add(root_file_zb);
             }
         }
      }
   }//cout<<"n_events_total ZB "<<treeZB.GetEntries()<<endl;
   delete vdirs_zb;
   delete vfiles_zb;


  //TTree* treeZB = new TTree(treeNameZB.c_str(),"");
  MyEvtId*           evtIdZB        = NULL;
  RPRootDumpReconstructedProton* rec_proton_left  = NULL;
  RPRootDumpReconstructedProton* rec_proton_right = NULL;
  map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
  vector<MyVertex>*  vertex_coll_ZB  = NULL;
  vector<MyPFCand>*  pFlow_coll_ZB   = NULL;
  vector<MyPFJet>*   pfJet_coll_ZB   = NULL;
//  TFile* fileZB = TFile::Open(fileNameZB.c_str(),"READ");

//  treeZB = (TTree*)fileZB->Get( treeNameZB.c_str() );
  int nevZB = int(treeZB.GetEntries());

  treeZB.SetBranchAddress("cmsEvtUA",&evtIdZB);
  treeZB.SetBranchAddress("rec_prot_left.",&rec_proton_left);
  treeZB.SetBranchAddress("rec_prot_right.",&rec_proton_right);
  treeZB.SetBranchAddress("cmsVerticesUA",&vertex_coll_ZB);
  treeZB.SetBranchAddress("cmsParticleFlowUA",&pFlow_coll_ZB);
  treeZB.SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll_ZB);
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
          treeZB.SetBranchAddress(br_name, &rp_track_info[id]);
       }
  }

  std::vector<MyZeroBiasData> zeroBias;

  for(int i_evt = 0; i_evt < nevZB && i_evt < nevt_max_corr; ++i_evt){

     if( ((i_evt+1) % 10000) == 0) cout <<int(double(i_evt+1)/1000)<<"k done"<<endl;

     //Filling the variables defined setting branches
     treeZB.GetEntry(i_evt);

     MyZeroBiasData mydata;
     mydata.proton_rec_left_valid = rec_proton_left->valid;
     mydata.proton_rec_left_t = rec_proton_left->t;
     mydata.proton_rec_left_xi = rec_proton_left->xi;
     mydata.proton_rec_right_valid = rec_proton_right->valid;
     mydata.proton_rec_right_t = rec_proton_right->t;
     mydata.proton_rec_right_xi = rec_proton_right->xi;

     MyVertex& primaryVertex = vertex_coll_ZB->at(0);
     /*for(vector<MyVertex>::iterator it_vtx_zb = vertex_coll_ZB->begin() ; it_vtx_zb != vertex_coll_ZB->end() ; ++it_vtx_zb){
        mydata.vtx_ndof = it_vtx_zb->ndof;    
        mydata.vtx_valid = it_vtx_zb->validity;
     }*/
     mydata.vtx_ndof = primaryVertex.ndof;
     mydata.vtx_valid = primaryVertex.validity;
     mydata.vtx_x = primaryVertex.x;
     mydata.vtx_y = primaryVertex.y;
     mydata.vtx_z = primaryVertex.z;

     if( pfJet_coll_ZB->size() > 0 ){
         MyBaseJet const& leadingJet = ( pfJet_coll_ZB->at(0) ).mapjet["ak5PFL2L3Residual"];
         mydata.leadingJet_pt = leadingJet.Pt();
         mydata.leadingJet_eta = leadingJet.Eta();
     }
     if( pfJet_coll_ZB->size() > 1 ){
         MyBaseJet const& secondJet = ( pfJet_coll_ZB->at(1) ).mapjet["ak5PFL2L3Residual"];
         mydata.secondJet_pt = secondJet.Pt();
         mydata.secondJet_eta = secondJet.Eta();
     }

     mydata.rp_track_valid_020 = rp_track_info[20]->valid;
     mydata.rp_track_valid_021 = rp_track_info[21]->valid;
     mydata.rp_track_valid_024 = rp_track_info[24]->valid;
     mydata.rp_track_valid_025 = rp_track_info[25]->valid;
     mydata.rp_track_valid_120 = rp_track_info[120]->valid;
     mydata.rp_track_valid_121 = rp_track_info[121]->valid;
     mydata.rp_track_valid_124 = rp_track_info[124]->valid;
     mydata.rp_track_valid_125 = rp_track_info[125]->valid;
     mydata.rp_x_024 = rp_track_info[24]->x;
     mydata.rp_y_024 = rp_track_info[24]->y;
     mydata.rp_x_025 = rp_track_info[25]->x;
     mydata.rp_y_025 = rp_track_info[25]->y;
     mydata.rp_x_124 = rp_track_info[124]->x;
     mydata.rp_y_124 = rp_track_info[124]->y;
     mydata.rp_x_125 = rp_track_info[125]->x;
     mydata.rp_y_125 = rp_track_info[125]->y;

     double sum_plus = 0;
     double sum_minus = 0;
     for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll_ZB->begin(); it_pfcand != pFlow_coll_ZB->end(); ++it_pfcand){
         int partType = it_pfcand->particleId;
         double eta = it_pfcand->Eta();
         double energy = it_pfcand->Energy();
         double pz = it_pfcand->Pz();
         sum_plus += (energy + pz);
         sum_minus += (energy - pz);
     }
     mydata.xi_cms_plus = sum_plus/8000;
     mydata.xi_cms_minus = sum_minus/8000;

     zeroBias.push_back(mydata);
  }
  //fileZB->Close();

  cout << zeroBias.size() << " events analyzed" << endl;

 // for(vector<MyTOTEMData>::const_iterator it_zb = zeroBias.begin(); it_zb != zeroBias.end(); it_zb++){
 //    cout<<it_zb->proton_rec_right_t<<endl;  
 // }

  //================================================



  rp_aperture_config();
  gRandom->SetSeed(12345);
  
  double nevents_vtx = 0; 
  double nevents_jet_rec = 0; 
  double nevents_proton_rec = 0; 
  double nevents_proton_gen = 0; 
  double nevents_pf = 0; 
  double nevents_jet_gen = 0; 
  double nevents_total = 0; 
  double nevents_jetrec_protkintkinxi = 0;
  double nevents_jetrec_protkintkinxi_rp = 0;
  double nevents_jetrec_protkintkinxi_rp_cut = 0;

  double nweight_total = 0; 
  double weight_total_leadingJet = 0; 
  double weight_total_secondJet = 0; 
  double weight_total_leadingJet_selected = 0; 
  double weight_total_secondJet_selected = 0; 
  double weight_total_Jet_selected = 0; 
  double weight_total_PF_selected = 0; 
  //starting Loop over files, stops at end of list of files or when reached nevt_max
  for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){
  
    TFile* file = TFile::Open(*itfiles,"READ");
    
    //getting the tree form the current file
    tree = (TTree*) file->Get( treeName.c_str() );

    //Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

    //adding branches to the tree ----------------------------------------------------------------------
    tree->SetBranchAddress("evtId",&evtId);
    tree->SetBranchAddress("generalTracks",&track_coll); 
    tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
    tree->SetBranchAddress("ak5CaloJets",&caloJet_coll);
    tree->SetBranchAddress("ak5PFJets",&pfJet_coll);
    tree->SetBranchAddress("ak5GenJets",&genJet_coll);
    tree->SetBranchAddress("particleFlow",&pFlow_coll);
    tree->SetBranchAddress("genKin",&genKin);
    tree->SetBranchAddress("genPart",&genPart);
    tree->SetBranchAddress("caloTowers",&caloTowers_coll);
  
    /*//Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;*/

     weight_st = -1.; 
     xi_cms_st = -999.; xi_totem_st = -999.; xi_totem_sel=-999.; xi_cms_minus_totem_st = -999.;
  
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){
    
    //printing the % of events done every 10k evts
    if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i_evt);

      bool passedHLT = false;
      bool passedvtx = false;
      bool jet1_rec_selected = false;
      bool jet2_rec_selected = false;
      bool jet1_gen_selected = false;
      bool jet2_gen_selected = false;
      bool jet1_pt20_rec_selected = false;
      bool jet2_pt20_rec_selected = false;
      bool jet1_pt20_gen_selected = false;
      bool jet2_pt20_gen_selected = false;
      bool pz_proton_max = false;
      bool PF_eta_max = false;
      bool PF_eta_min = false;
      bool xi_negat_gen = false;
      bool xi_posit_gen = false;
      
      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
      double weight = genKin->genWeight; 
      double event_weight = (mc == "pomwig") ? 1.0 : weight;
      nweight_total += weight; 
      ++nevents_total;
  
      bool sd_minus_pythia = false;
      bool sd_plus_pythia = false;
      bool dd_pythia = false;

      if (mc == "pythia8_diff"){
         int process_id = genKin->MCProcId;
         if (process_id == 103) sd_minus_pythia = true;
         if (process_id == 104) sd_plus_pythia = true;
         if (process_id == 105) dd_pythia = true;
      }

 
 	    
     // Vertices
      for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
//        if (it_vtx!=vertex_coll->begin()) continue;
//	if( it_vtx->ndof>4 ) passedvtx = true;   
            //histosTH1F["vtx_zpos"]->Fill( it_vtx->z, event_weight );
            //histosTH1F["vtx_xpos"]->Fill( it_vtx->x, event_weight );
            //histosTH1F["vtx_ypos"]->Fill( it_vtx->y, event_weight );}
      }
//      if(!passedvtx) continue;
      MyVertex& primaryVertex = vertex_coll->at(0);
      bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity && primaryVertex.ndof > 4);// && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);
      //if (!select_Vertex) continue;
      if (select_Vertex) ++nevents_vtx;
      histosTH1F["vtx_xpos"]->Fill( primaryVertex.x, event_weight );
      histosTH1F["vtx_ypos"]->Fill( primaryVertex.y, event_weight );
      histosTH1F["vtx_zpos"]->Fill( primaryVertex.z, event_weight );
      //double zvtx = gRandom->Gaus(-0.35840,5.8481)/gRandom->Gaus(0.414212,6.14986);  
//      TF1* zvtx_corr_data = new TF1("zvtx_corr_data", "gaus", -30,30);
//      zvtx_corr_data->SetParameter(0,-0.35840);
//      zvtx_corr_data->SetParameter(1,5.8481);
//      TF1* zvtx_corr_mc = new TF1("zvtx_corr_mc", "gaus", -30,30);
//      zvtx_corr_mc->SetParameter(0,0.414212);
//      zvtx_corr_mc->SetParameter(1, 6.14986);
      double zvtx = 1;//zvtx_corr_data->Eval(primaryVertex.z)/zvtx_corr_mc->Eval(primaryVertex.z);
// cout<<zvtx<<endl;

      // Tracks
      int n_tracks_selected = 0;
      for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
         histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
         histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
         histosTH1F["track_phi"]->Fill( it_trk->Phi(), event_weight );

         ++n_tracks_selected;
      }
      histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );



      //Jets with pt>30Gev and !eta!<2
      Double_t Jet1_pt_rec; 
      Double_t Jet2_pt_rec; 
      Double_t Jet1_eta_rec; 
      Double_t Jet2_eta_rec; 
      Double_t Jet1_phi_rec, Jet1_px_rec, Jet1_py_rec, Jet1_pz_rec, Jet1_energy_rec; 
      Double_t Jet2_phi_rec, Jet2_px_rec, Jet2_py_rec, Jet2_pz_rec, Jet2_energy_rec;
      
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
         weight_total_leadingJet += event_weight;
	 Jet1_pt_rec = leadingJet.Pt(); 
	 Jet1_eta_rec = leadingJet.Eta(); 
	 Jet1_phi_rec = leadingJet.Phi(); 
	 Jet1_px_rec = leadingJet.Px(); 
	 Jet1_py_rec = leadingJet.Py(); 
	 Jet1_pz_rec = leadingJet.Pz(); 
	 Jet1_energy_rec = leadingJet.E(); 
      }
      //if(!jet1_rec_selected) continue;
      
      if( pfJet_coll->size() > 1 ){
	 MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet[jetCorrName];
         weight_total_secondJet += event_weight;
         Jet2_pt_rec = secondJet.Pt(); 
	 Jet2_eta_rec = secondJet.Eta(); 
	 Jet2_phi_rec = secondJet.Phi(); 
         Jet2_px_rec = secondJet.Px();
         Jet2_py_rec = secondJet.Py();
         Jet2_pz_rec = secondJet.Pz();
         Jet2_energy_rec = secondJet.E();
      }
      //if(!jet2_rec_selected) continue;
      weight_total_Jet_selected += event_weight;
//       histosTH1F["DeltaPhiJet"]->Fill( deltaphi, event_weight  );	
      if(Jet1_pt_rec > 20. && fabs(Jet1_eta_rec)<4.4) jet1_pt20_rec_selected = true;
      if(Jet1_pt_rec > 30. && fabs(Jet1_eta_rec)<4.4) jet1_rec_selected = true;
      if(Jet2_pt_rec > 20. && fabs(Jet2_eta_rec)<4.4 ) jet2_pt20_rec_selected = true;
      if(Jet2_pt_rec > 30. && fabs(Jet2_eta_rec)<4.4 ) jet2_rec_selected = true; 
      double mass_jets_rec= sqrt(pow(Jet1_energy_rec+Jet2_energy_rec,2)-pow(Jet1_px_rec+Jet2_px_rec,2)-pow(Jet1_py_rec+Jet2_py_rec,2)-pow(Jet1_pz_rec+Jet2_pz_rec,2));
      double x_minus_rec = ((Jet1_energy_rec-Jet1_pz_rec)+(Jet2_energy_rec-Jet2_pz_rec))/8000; 
      double x_plus_rec = ((Jet1_energy_rec+Jet1_pz_rec)+(Jet2_energy_rec+Jet2_pz_rec))/8000; 

      ///Fit 
      TF1* func_trigger = new TF1("func_trigger", fFermiLike, 0., 20., 2);
      func_trigger->SetParameter(0,5.525);
      func_trigger->SetParameter(1,0.529);
      double eff_trigger = func_trigger->Eval(Jet2_pt_rec);



      //Jet generated level         
      double leadingJet_pt_gen = -999;
      double secondJet_pt_gen = -999;
      double Jet1_energy_gen;
      double Jet1_px_gen;
      double Jet1_py_gen;
      double Jet1_pz_gen;
      double Jet2_energy_gen;
      double Jet2_px_gen;
      double Jet2_py_gen;
      double Jet2_pz_gen;
      double Jet1_eta_gen;
      double Jet1_phi_gen;
      double Jet2_eta_gen;
      double Jet2_phi_gen;


      for(vector<MyGenJet>::iterator it_genjet = genJet_coll->begin(); it_genjet != genJet_coll->end(); ++it_genjet){
         double jet_pt_gen = it_genjet->Pt();
         double jet_eta_gen = it_genjet->Eta();
         double jet_ene_gen = it_genjet->E();
         double jet_px_gen = it_genjet->Px();
         double jet_py_gen = it_genjet->Py();
         double jet_pz_gen = it_genjet->Pz();
         double jet_phi_gen = it_genjet->Phi();

        // if (fabs(jet_eta_gen)>4.4) continue;

         if (jet_pt_gen>leadingJet_pt_gen){
             leadingJet_pt_gen = jet_pt_gen;
             Jet1_energy_gen = jet_ene_gen;
             Jet1_px_gen = jet_px_gen;
             Jet1_py_gen = jet_py_gen;
             Jet1_pz_gen = jet_pz_gen;
             Jet1_eta_gen = jet_eta_gen;
             Jet1_phi_gen = jet_phi_gen;
         }
         if (jet_pt_gen>secondJet_pt_gen && jet_pt_gen<leadingJet_pt_gen){
             secondJet_pt_gen = jet_pt_gen;
             Jet2_energy_gen = jet_ene_gen;
             Jet2_px_gen = jet_px_gen;
             Jet2_py_gen = jet_py_gen;
             Jet2_pz_gen = jet_pz_gen;
             Jet2_eta_gen = jet_eta_gen;
             Jet2_phi_gen = jet_phi_gen;
         }
      }
      if(leadingJet_pt_gen>30. && fabs(Jet1_eta_gen)<4.4) jet1_gen_selected = true;
      if(secondJet_pt_gen>30. && fabs(Jet2_eta_gen)<4.4) jet2_gen_selected = true;
      if(leadingJet_pt_gen>20. && fabs(Jet1_eta_gen)<4.4) jet1_pt20_gen_selected = true;
      if(secondJet_pt_gen>20. && fabs(Jet2_eta_gen)<4.4) jet2_pt20_gen_selected = true;

      double mass_jets_gen= sqrt(pow(Jet1_energy_gen+Jet2_energy_gen,2)-pow(Jet1_px_gen+Jet2_px_gen,2)-pow(Jet1_py_gen+Jet2_py_gen,2)-pow(Jet1_pz_gen+Jet2_pz_gen,2));
      double x_minus_gen = ((Jet1_energy_gen-Jet1_pz_gen)+(Jet2_energy_gen-Jet2_pz_gen))/8000;
      double x_plus_gen = ((Jet1_energy_gen+Jet1_pz_gen)+(Jet2_energy_gen+Jet2_pz_gen))/8000;


      // Particle-flow
      double soma1 = 0;
      double soma2 = 0;
      double eta_max=-999.;
      double eta_min=999.;
      double cm = 8000;
      bool pf = false;

      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
         int partType = it_pfcand->particleId;
         double eta = it_pfcand->Eta();
         double energy = it_pfcand->Energy();
         double pz = it_pfcand->Pz();

         // HF eta rings 29, 30, 40, 41
         //if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || (fabs(eta) >= 4.730) ) continue;

         // Apply thresholds
         if( !pflowThreshold(*it_pfcand,thresholdsPFlow) ) continue;

         soma1 += (energy + pz);
         soma2 += (energy - pz);

         if (eta > eta_max) {eta_max = eta; PF_eta_max = true;}
         if (eta < eta_min) {eta_min = eta; PF_eta_min = true;}


         histosTH2F["energyVsEtaAllTypes"]->Fill( eta, energy, event_weight );

         if(partType == MyPFCand::X)
            histosTH2F["energyVsEtaUndefined"]->Fill( eta, energy, event_weight );
         else if(partType == MyPFCand::h)
            histosTH2F["energyVsEtaChargedHadron"]->Fill( eta, energy, event_weight );
         else if(partType == MyPFCand::e)
            histosTH2F["energyVsEtaElectron"]->Fill( eta, energy, event_weight );
         else if(partType == MyPFCand::mu)
            histosTH2F["energyVsEtaMuon"]->Fill( eta, energy, event_weight );
         else if(partType == MyPFCand::gamma)
            histosTH2F["energyVsEtaGamma"]->Fill( eta, energy, event_weight );
         else if(partType == MyPFCand::h0)
            histosTH2F["energyVsEtaNeutralHadron"]->Fill( eta, energy, event_weight );
         else if(partType == MyPFCand::h_HF){
            histosTH2F["energyVsEtaHadronHF"]->Fill( eta, energy, event_weight );}
         else if(partType == MyPFCand::egamma_HF)
            histosTH2F["energyVsEtaEGammaHF"]->Fill( eta, energy, event_weight );

       }
       if(!PF_eta_max) continue;
       if(!PF_eta_min) continue;
       weight_total_PF_selected += event_weight;

       ++nevents_pf;
      // if(PF_eta_max || PF_eta_min) ++nevents_pf;

       double xi_plus_Reco = soma1/cm;
       double xi_minus_Reco = soma2/cm;
       double delta_eta_maxmin = eta_max - eta_min;


 
      //GenPart
      double genEPlusPz = 0;
      double genEMinusPz = 0;
     // double cm = 8000;
      Double_t proton_pi = 4000;
      Double_t proton_pz_plus=-999;
      Double_t proton_px_plus = -999.;
      Double_t proton_py_plus = -999.;
      Double_t proton_energy_plus = 0.;
      Double_t proton_mass_minus=999;
      Double_t proton_pz_minus=999;
      Double_t proton_px_minus = 999.;
      Double_t proton_py_minus = 999.;
      Double_t proton_pt_minus = 999.;
      Double_t proton_eta_minus = 999.;
      Double_t proton_energy_minus = 0.;
      Double_t px_gen, pt_gen, mass_gen;
      Double_t py_gen;
      Double_t pz_gen;
      Double_t energy_gen;
      Double_t proton_pf;
      Double_t eta_gen;
      
      for(vector<MyGenPart>::iterator it_genpart = genPart->begin(); it_genpart != genPart->end(); ++it_genpart){
 
	 //double eta_gen = it_genpart->Eta();
         int status = it_genpart->status;
         int id = it_genpart->pdgId;
	 
	 if (status == 1) {
            energy_gen = it_genpart->Energy();
            px_gen = it_genpart->Px();
            py_gen = it_genpart->Py();
            pz_gen = it_genpart->Pz();
            pt_gen = it_genpart->Pt();
            eta_gen = it_genpart->Eta();
            mass_gen = it_genpart->M();
	    proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);  
	    if (id != 2212) {
   	       genEPlusPz += (energy_gen + pz_gen);
	       genEMinusPz += (energy_gen - pz_gen);
            }
	    if (id == 2212) {
             double pz_cut = 0.7*proton_pi;
             if (fabs(pz_gen) > pz_cut){

	        if (pz_gen > proton_pz_plus) {
                   if (mc == "pythia8_diff" && !sd_plus_pythia) continue;
                    proton_pz_plus = pz_gen; proton_energy_plus = energy_gen;
                    proton_px_plus = px_gen; proton_py_plus = py_gen;       
                }
                if (pz_gen < proton_pz_minus) {
                   if (mc == "pythia8_diff" && !sd_minus_pythia) continue;
                   proton_pz_minus = pz_gen; proton_energy_minus = energy_gen;
                   proton_px_minus = px_gen; proton_py_minus = py_gen;
                   proton_pt_minus = pt_gen; proton_mass_minus = mass_gen;
                   proton_eta_minus = eta_gen; 
                }
             }
            }
	 }
      }
      histosTH1F["energy_gen_proton_signal"]->Fill( proton_energy_minus, 1 );
      histosTH1F["px_gen_proton_signal"]->Fill( proton_px_minus, 1 );
      histosTH1F["py_gen_proton_signal"]->Fill( proton_py_minus, 1 );
      histosTH1F["pt_gen_proton_signal"]->Fill( proton_pt_minus, 1 );
      histosTH1F["pt2_gen_proton_signal"]->Fill( proton_pt_minus*proton_pt_minus, 1 );
      histosTH1F["mass_gen_proton_signal"]->Fill( proton_mass_minus, 1 );
//cout<<proton_mass_minus<<endl;

      
      double xi_plus_gen = genEPlusPz/cm; //cout<<xi1_gen<<endl;
      double xi_minus_gen = genEMinusPz/cm;
      double xi_minus_proton_gen = -1.;
      double xi_minus_proton_smear_gen = -1.;
      double t_minus_proton_gen = 0.;
      double p_minus_proton_gen = 0.;
      double p_minus_proton_smear_gen = 0.;
      double t_minus_proton_smear_gen = 0.;
      double thx_minus_proton = 0.;
      double thy_minus_proton = 0.;
      double thx_minus_proton_smear = 0.;
      double thy_minus_proton_smear = 0.;
      //double thx_proton_minus = 0.;
      //double thy_proton_minus = 0.;

      double xi_plus_proton_gen = -1.;
      double xi_plus_proton_smear_gen = -1.;
      double t_plus_proton_gen = 0.;
      double p_plus_proton_gen = 0.;
      double p_plus_proton_smear_gen = 0.;
      double t_plus_proton_smear_gen = 0.;
      double thx_plus_proton = 0.;
      double thy_plus_proton = 0.;
      double thx_plus_proton_smear = 0.;
      double thy_plus_proton_smear = 0.;

      double proton_px_minus_smear = 0;
      double proton_py_minus_smear = 0;
      double proton_pz_minus_smear = 0;
      double proton_energy_minus_smear = 0;
      double proton_px_plus_smear = 0;
      double proton_py_plus_smear = 0;
      double proton_pz_plus_smear = 0;
      double proton_energy_plus_smear = 0;
      double v_x = 0;
      double v_y = 0;
      double v_z = 0;

     // generate vertex smearing
      vtx_smearing(v_x, v_y, v_z);

      //rp parametrization
      bool proton_minus_rp_accept_120 = false;
      bool proton_minus_rp_accept_121 = false;
      bool proton_minus_rp_accept_122 = false;
      bool proton_minus_rp_accept_123 = false;
      bool proton_minus_rp_accept_124 = false;
      bool proton_minus_rp_accept_125 = false;
      bool proton_minus_rp_accept_020 = false;

      bool proton_plus_rp_accept_020 = false;
      bool proton_plus_rp_accept_021 = false;
      bool proton_plus_rp_accept_022 = false;
      bool proton_plus_rp_accept_023 = false;
      bool proton_plus_rp_accept_024 = false;
      bool proton_plus_rp_accept_025 = false;
      bool proton_plus_rp_accept_120 = false;

      bool fiducial_cut_rp_024=false;
      bool fiducial_cut_rp_025=false;
      bool fiducial_cut_rp_124=false;
      bool fiducial_cut_rp_125=false;
      bool fiducial_cut_rp_124_smear=false;
      bool fiducial_cut_rp_125_smear=false;

      std::map<int,std::vector<double> > proton_plus_pars;
      std::map<int,std::vector<double> > proton_minus_pars;

      if( proton_pz_plus > 0.){
         //beam smearing
         beam_smearing(proton_px_plus, proton_py_plus, proton_pz_plus, proton_energy_plus, proton_px_plus_smear, proton_py_plus_smear, proton_pz_plus_smear, proton_energy_plus_smear);

         xi_plus_proton_gen =  ( 1 - (proton_pz_plus/proton_pi) );
         //xi_proton_smear_gen =  ( 1 - (proton_pz_smear/proton_pi) );
         //t_proton_gen = -2*( (proton_pi*proton_energy_plus) - (proton_pi*proton_pz_plus) );
         //t_proton_smear_gen = -2*( (proton_pi*proton_energy_smear) - (proton_pi*proton_pz_smear) );
         //thx_proton = atan(-proton_px_smear/proton_pi);
         //thy_proton = atan(proton_py_smear/proton_pi);

         p_plus_proton_gen = sqrt(proton_px_plus*proton_px_plus+proton_py_plus*proton_py_plus+proton_pz_plus*proton_pz_plus);
         p_plus_proton_smear_gen = sqrt(proton_px_plus_smear*proton_px_plus_smear+proton_py_plus_smear*proton_py_plus_smear+proton_pz_plus_smear*proton_pz_plus_smear);
         thx_plus_proton = atan(-proton_px_plus/proton_pi);
         thy_plus_proton = atan(proton_py_plus/proton_pi);
         thx_plus_proton_smear = atan(-proton_px_plus_smear/proton_pi) + gRandom->Gaus(0,25.10e-6);
         thy_plus_proton_smear = atan(proton_py_plus_smear/proton_pi) + gRandom->Gaus(0,2.42e-6);
         t_plus_proton_gen = -p_plus_proton_gen*p_plus_proton_gen*((thx_plus_proton*thx_plus_proton)+(thy_plus_proton*thy_plus_proton)); 
         t_plus_proton_smear_gen = -p_plus_proton_smear_gen*p_plus_proton_smear_gen*((thx_plus_proton_smear*thx_plus_proton_smear)+(thy_plus_proton_smear*thy_plus_proton_smear)); 
         //t_proton_smear_gen = -2*( (proton_pi*proton_energy_smear) + (proton_pi*proton_pz_smear) ); 

         double delta_thx_plus = thx_plus_proton_smear-thx_plus_proton;
         xi_plus_proton_smear_gen = (400*delta_thx_plus) + xi_plus_proton_gen; 
         histosTH2F["delta_xi_thx_plus"]->Fill(delta_thx_plus, xi_plus_proton_smear_gen-xi_plus_proton_gen, event_weight);

         //FIXME
         double out_x, out_thx, out_y, out_thy, out_xi;
         proton_plus_rp_accept_020 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 20, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[20] = std::vector<double>(5,0.);
         proton_plus_pars[20][0] = out_x; proton_plus_pars[20][1] = out_y;
         proton_plus_pars[20][2] = out_thx; proton_plus_pars[20][3] = out_thy;
         proton_plus_pars[20][4] = out_xi;

         proton_plus_rp_accept_024 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 24, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[24] = std::vector<double>(5,0.);
         proton_plus_pars[24][0] = out_x; proton_plus_pars[24][1] = out_y;
         proton_plus_pars[24][2] = out_thx; proton_plus_pars[24][3] = out_thy;
         proton_plus_pars[24][4] = out_xi;
         histosTH2F["pos_y_vs_x_proton_plus_024_025"]->Fill( proton_plus_pars[24][0], proton_plus_pars[24][1] , event_weight );
         fiducial_cut_rp_024 = proton_plus_pars[24][0]>0 && proton_plus_pars[24][0]<0.006 && proton_plus_pars[24][1]>0.0084 && proton_plus_pars[24][1]<0.029;
         //fiducial_cut_rp_024 = proton_plus_pars[24][0]>0 && proton_plus_pars[24][0]<0.006 && proton_plus_pars[24][1]>0.0082 && proton_plus_pars[24][1]<0.029;
      
         proton_plus_rp_accept_025 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 25, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[25] = std::vector<double>(5,0.);
         proton_plus_pars[25][0] = out_x; proton_plus_pars[25][1] = out_y;
         proton_plus_pars[25][2] = out_thx; proton_plus_pars[25][3] = out_thy;
         proton_plus_pars[25][4] = out_xi;
         histosTH2F["pos_y_vs_x_proton_plus_024_025"]->Fill( proton_plus_pars[25][0], proton_plus_pars[25][1] , event_weight );
         fiducial_cut_rp_025 = proton_plus_pars[25][0]>0 && proton_plus_pars[25][0]<0.006 && proton_plus_pars[25][1]<-0.0084 && proton_plus_pars[25][1]>-0.029;
         //fiducial_cut_rp_025 = proton_plus_pars[25][0]>0 && proton_plus_pars[25][0]<0.006 && proton_plus_pars[25][1]<-0.0082 && proton_plus_pars[25][1]>-0.029;

         proton_plus_rp_accept_021 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 21);
         proton_plus_rp_accept_022 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 22);
         proton_plus_rp_accept_023 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 23);
         proton_plus_rp_accept_024 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 24);
         proton_plus_rp_accept_025 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 25);
         proton_plus_rp_accept_020 = protonRPDetected(0., thx_plus_proton, 0., thy_plus_proton, -xi_plus_proton_gen, 20);

         if (fiducial_cut_rp_024 || fiducial_cut_rp_025){
             histosTH2F["pos_y_vs_x_proton_plus_024_025_accept"]->Fill( proton_plus_pars[24][0], proton_plus_pars[24][1] , event_weight );
             histosTH2F["pos_y_vs_x_proton_plus_024_025_accept"]->Fill( proton_plus_pars[25][0], proton_plus_pars[25][1] , event_weight );
         }
      }
 
      double t_minus_proton_gen_gauss;
      if( proton_pz_minus < 0.){
         //beam smearing
         beam_smearing(proton_px_minus, proton_py_minus, proton_pz_minus, proton_energy_minus, proton_px_minus_smear, proton_py_minus_smear, proton_pz_minus_smear, proton_energy_minus_smear);
         xi_minus_proton_gen = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
         //xi_proton_smear_gen = (proton_pz_smear < 0.) ? ( 1 + (proton_pz_smear/proton_pi) ) : -1.;
         t_minus_proton_gen_gauss = -2*( (proton_pi*proton_energy_minus) + (proton_pi*proton_pz_minus) ); 

         p_minus_proton_gen = sqrt(proton_px_minus*proton_px_minus+proton_py_minus*proton_py_minus+proton_pz_minus*proton_pz_minus);
         p_minus_proton_smear_gen = sqrt(proton_px_minus_smear*proton_px_minus_smear+proton_py_minus_smear*proton_py_minus_smear+proton_pz_minus_smear*proton_pz_minus_smear);
         thx_minus_proton = atan(-proton_px_minus/proton_pi);
         thy_minus_proton = atan(proton_py_minus/proton_pi);
         thx_minus_proton_smear = atan(-proton_px_minus/proton_pi) + gRandom->Gaus(0,25.10e-6);
         thy_minus_proton_smear = atan(proton_py_minus_smear/proton_pi) + gRandom->Gaus(0,2.42e-6);
         t_minus_proton_gen = -p_minus_proton_gen*p_minus_proton_gen*((thx_minus_proton*thx_minus_proton)+(thy_minus_proton*thy_minus_proton)); 
         t_minus_proton_smear_gen = -p_minus_proton_smear_gen*p_minus_proton_smear_gen*((thx_minus_proton_smear*thx_minus_proton_smear)+(thy_minus_proton_smear*thy_minus_proton_smear)); 
         //t_proton_smear_gen = -2*( (proton_pi*proton_energy_smear) + (proton_pi*proton_pz_smear) ); 

         double delta_thx_minus = thx_minus_proton_smear-thx_minus_proton;
         xi_minus_proton_smear_gen = (400*delta_thx_minus) + xi_minus_proton_gen; 
         histosTH2F["delta_xi_thx_minus"]->Fill(delta_thx_minus, xi_minus_proton_smear_gen-xi_minus_proton_gen, event_weight);
         histosTH1F["thx_proton"]->Fill(thx_minus_proton, event_weight);
         histosTH1F["thx_proton_smear"]->Fill(thx_minus_proton_smear, event_weight);
         histosTH1F["thy_proton"]->Fill(thy_minus_proton, event_weight);
         histosTH1F["thy_proton_smear"]->Fill(thy_minus_proton_smear, event_weight);

         double out_x, out_thx, out_y, out_thy, out_xi;
         proton_minus_rp_accept_120 = protonRPDetected(v_x, thx_minus_proton, v_y, thy_minus_proton, -xi_minus_proton_smear_gen, 120, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[120] = std::vector<double>(5,0.);
         proton_minus_pars[120][0] = out_x; proton_minus_pars[120][1] = out_y;
         proton_minus_pars[120][2] = out_thx; proton_minus_pars[120][3] = out_thy;
         proton_minus_pars[120][4] = out_xi;

         proton_minus_rp_accept_124 = protonRPDetected(v_x, thx_minus_proton, v_y, thy_minus_proton, -xi_minus_proton_smear_gen, 124, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[124] = std::vector<double>(5,0.);
         proton_minus_pars[124][0] = out_x; proton_minus_pars[124][1] = out_y;
         proton_minus_pars[124][2] = out_thx; proton_minus_pars[124][3] = out_thy;
         proton_minus_pars[124][4] = out_xi;
         histosTH2F["pos_y_vs_x_proton_minus_124_125"]->Fill( proton_minus_pars[124][0], proton_minus_pars[124][1] , event_weight );
         fiducial_cut_rp_124 = proton_minus_pars[124][0]>0 && proton_minus_pars[124][0]<0.006 && proton_minus_pars[124][1]>0.0084 && proton_minus_pars[124][1]<0.027;
         fiducial_cut_rp_124_smear = proton_minus_pars[124][0]>0 && proton_minus_pars[124][0]<0.006 && proton_minus_pars[124][1]>0.0082 && proton_minus_pars[124][1]<0.027;//lower y bound by 200microm

         proton_minus_rp_accept_125 = protonRPDetected(v_x, thx_minus_proton, v_y, thy_minus_proton, -xi_minus_proton_smear_gen, 125, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[125] = std::vector<double>(5,0.);
         proton_minus_pars[125][0] = out_x; proton_minus_pars[125][1] = out_y;
         proton_minus_pars[125][2] = out_thx; proton_minus_pars[125][3] = out_thy;
         proton_minus_pars[125][4] = out_xi;
         histosTH2F["pos_y_vs_x_proton_minus_124_125"]->Fill( proton_minus_pars[125][0], proton_minus_pars[125][1] , event_weight );
         fiducial_cut_rp_125 = proton_minus_pars[125][0]>0 && proton_minus_pars[125][0]<0.006 && proton_minus_pars[125][1]<-0.0084 && proton_minus_pars[125][1]>-0.027;
         fiducial_cut_rp_125_smear = proton_minus_pars[125][0]>0 && proton_minus_pars[125][0]<0.006 && proton_minus_pars[125][1]<-0.0082 && proton_minus_pars[125][1]>-0.027;

         if (fiducial_cut_rp_124 || fiducial_cut_rp_125){ 
             histosTH2F["pos_y_vs_x_proton_minus_124_125_accept"]->Fill( proton_minus_pars[124][0], proton_minus_pars[124][1] , event_weight );
             histosTH2F["pos_y_vs_x_proton_minus_124_125_accept"]->Fill( proton_minus_pars[125][0], proton_minus_pars[125][1] , event_weight );
         }

         proton_minus_rp_accept_121 = protonRPDetected(v_x, thx_minus_proton, v_y, thy_minus_proton, -xi_minus_proton_smear_gen, 121);
         proton_minus_rp_accept_122 = protonRPDetected(v_x, thx_minus_proton, v_y, thy_minus_proton, -xi_minus_proton_smear_gen, 122);
         proton_minus_rp_accept_123 = protonRPDetected(v_x, thx_minus_proton, v_y, thy_minus_proton, -xi_minus_proton_smear_gen, 123);
         //proton_minus_rp_accept_124 = protonRPDetected(v_x, thx_proton_minus, v_y, thy_proton_minus, -xi_proton_minus_smear, 124);
         //proton_minus_rp_accept_125 = protonRPDetected(v_x, thx_proton_minus, v_y, thy_proton_minus, -xi_proton_minus_smear, 125);
         //proton_minus_rp_accept_120 = protonRPDetected(v_x, thx_proton_minus, v_y, thy_proton_minus, -xi_proton_minus_smear, 120);
      }

      //Access Zero-bias data
      //gRandom->SetSeed(12345);
      int i_evt_ZB = 0 + gRandom->Rndm()*(zeroBias.size());
      MyZeroBiasData const & zeroBiasData = zeroBias.at(i_evt_ZB);
      bool valid_vtx = zeroBiasData.vtx_valid;
      double ndof_vtx = zeroBiasData.vtx_ndof;
      //if (valid_vtx && ndof_vtx>=4) continue;
      double xi_cms_plus_zb = zeroBiasData.xi_cms_plus;
      double xi_cms_minus_zb = zeroBiasData.xi_cms_minus;
      double xi_minus_cms = xi_cms_minus_zb+xi_minus_Reco;
      double xi_plus_cms = xi_cms_plus_zb+xi_plus_Reco;
 

      //totem proton reconstructed
      double xi_proton_rec, xi_proton_smear_rec;
      double t_proton_rec, t_proton_smear_rec;
      double proton_beta_rec, proton_beta_gen;

      //gauss smearing
      float sigma_xi56=0.00720615 - 0.0418783*xi_minus_proton_gen + 0.0999515*xi_minus_proton_gen*xi_minus_proton_gen; // sigma56 vs xi from Hubert
      double xi_minus_proton_rec_gauss = xi_minus_proton_gen + gRandom->Gaus(0,sigma_xi56);
      double sigma_t56=0.233365*t_minus_proton_gen - 0.0975751*t_minus_proton_gen*t_minus_proton_gen;  // sigma_t56 vs t from Hubert
      double t_minus_proton_rec_gauss = t_minus_proton_gen_gauss + gRandom->Gaus(0,sigma_t56);

      //theta smearing
      double xi_minus_proton_smear_rec = xi_minus_proton_smear_gen;// + gRandom->Gaus(0,sigma_xi56_smear);
      double t_minus_proton_smear_rec = t_minus_proton_smear_gen;// + gRandom->Gaus(0,sigma_t56_smear);
      double proton_beta_minus_rec = x_minus_rec/xi_minus_proton_smear_rec;
      double proton_beta_minus_gen = x_minus_gen/xi_minus_proton_gen;
      //float sigma_xi45 = 0.00714986 - 0.0408903*xi_proton_gen + 0.0965813*xi_proton_gen*xi_proton_gen;
      //float sigma_xi45_smear = 0.00714986 - 0.0408903*xi_proton_smear_gen + 0.0965813*xi_proton_smear_gen*xi_proton_smear_gen;
      //xi_proton_rec = xi_proton_gen + gRandom->Gaus(0,sigma_xi45);
      double xi_plus_proton_smear_rec = xi_plus_proton_smear_gen;// + gRandom->Gaus(0,sigma_xi45_smear);
      //double sigma_t45 = 0.233365*t_proton_gen - 0.0975751*t_proton_gen*t_proton_gen;
      //double sigma_t45_smear = 0.233365*t_proton_smear_gen - 0.0975751*t_proton_smear_gen*t_proton_smear_gen;
      //t_proton_rec = t_proton_gen + gRandom->Gaus(0,sigma_t45);
      double t_plus_proton_smear_rec = t_plus_proton_smear_gen;// + gRandom->Gaus(0,sigma_t45_smear);
      double proton_beta_plus_rec = x_plus_rec/xi_plus_proton_smear_rec;
      double proton_beta_plus_gen = x_plus_gen/xi_plus_proton_gen;


      //rp_accept
      bool proton_minus_rp_accept_mc = ( proton_minus_rp_accept_120 && proton_minus_rp_accept_124 && fiducial_cut_rp_124 ) || ( proton_minus_rp_accept_121 && proton_minus_rp_accept_125 && fiducial_cut_rp_125);// || ( proton_minus_rp_accept_122 && proton_minus_rp_accept_123 );
      bool proton_plus_rp_accept_mc = ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 && fiducial_cut_rp_024 ) || ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 && fiducial_cut_rp_025 );// || ( proton_plus_rp_accept_022 && proton_plus_rp_accept_023 );
      bool proton_minus_rp_accept_mc_smear = ( proton_minus_rp_accept_120 && proton_minus_rp_accept_124 && fiducial_cut_rp_124_smear ) || ( proton_minus_rp_accept_121 && proton_minus_rp_accept_125 && fiducial_cut_rp_125_smear);// || ( proton_minus_rp_accept_122 && proton_minus_rp_accept_123 );

      bool rp_minus_sel =  !proton_plus_rp_accept_mc && proton_minus_rp_accept_mc;
      bool rp_plus_sel =  proton_plus_rp_accept_mc && !proton_minus_rp_accept_mc;

      bool jet_rec = jet1_rec_selected && jet2_rec_selected;
      bool jet_gen = jet1_gen_selected && jet2_gen_selected;
      bool jet_pt20_rec = jet1_pt20_rec_selected && jet2_pt20_rec_selected;
      bool jet_pt20_gen = jet1_pt20_gen_selected && jet2_pt20_gen_selected;
      bool proton_kinec_accept_t_minus_rec = fabs(t_minus_proton_smear_rec)>=0.03 && fabs(t_minus_proton_smear_rec)<=1;
      bool proton_kinec_accept_xi_minus_rec = xi_minus_proton_smear_rec<=0.1;
      bool proton_kinec_accept_t_minus_gen = fabs(t_minus_proton_gen)>=0.03 && fabs(t_minus_proton_gen)<=1;
      bool proton_kinec_accept_xi_minus_gen = xi_minus_proton_gen<=0.1;
      bool xi_cms_totem_minus_cut =  xi_minus_cms-xi_minus_proton_smear_rec<=0;
      bool proton_kinec_accept_t_plus_rec = fabs(t_plus_proton_smear_rec)>=0.03 && fabs(t_plus_proton_smear_rec)<=1;
      bool proton_kinec_accept_xi_plus_rec = xi_plus_proton_smear_rec<=0.1;
      bool proton_kinec_accept_t_plus_gen = fabs(t_plus_proton_gen)>=0.03 && fabs(t_plus_proton_gen)<=1;
      bool proton_kinec_accept_xi_plus_gen = xi_plus_proton_gen<=0.1;
      bool xi_cms_totem_plus_cut =  xi_plus_cms-xi_plus_proton_smear_rec<=0;

      //...
      Double_t deltapt_rec = abs(Jet1_pt_rec - Jet2_pt_rec);
      Double_t deltapt_gen = abs(leadingJet_pt_gen - secondJet_pt_gen);
      Double_t deltaeta_rec = abs(Jet1_eta_rec - Jet2_eta_rec);
      Double_t deltaeta_gen = abs(Jet1_eta_gen - Jet2_eta_gen);
      Double_t deltaphi_rec = abs(Jet1_phi_rec - Jet2_phi_rec);
      Double_t deltaphi_gen = abs(Jet1_phi_gen - Jet2_phi_gen);
      Double_t eta_average_rec = 0.5*(Jet1_eta_rec + Jet2_eta_rec);
      Double_t eta_average_gen = 0.5*(Jet1_eta_gen + Jet2_eta_gen);
      double mass_x_minus_gen = sqrt(xi_minus_proton_gen*8000);
      double mass_x_minus_rec = sqrt(xi_minus_proton_smear_rec*8000);

      double reweigth_beta_minus;
      double reweigth_beta_plus;
      if (mc == "pomwig"){
          reweigth_beta_minus = (reweigth && proton_beta_minus_gen<=0.7) ? func_right->Eval(proton_beta_minus_gen) : 1; 
          reweigth_beta_plus = (reweigth && proton_beta_plus_gen<=0.7) ? func_left->Eval(proton_beta_plus_gen) : 1; 
      }      
      if (mc == "pythia8_diff"){
          reweigth_beta_minus = (reweigth && proton_beta_minus_gen<=0.7) ? func_right_pythia->Eval(proton_beta_minus_gen) : 1; 
          reweigth_beta_plus = (reweigth && proton_beta_plus_gen<=0.7) ? func_left_pythia->Eval(proton_beta_plus_gen) : 1; 
      } 
     
      double eff_proton = 0.94;//proton efficiency
      double event_weight_rec_minus = event_weight*reweigth_beta_minus*eff_trigger*eff_proton;
      double event_weight_rec_plus = event_weight*reweigth_beta_plus*eff_trigger*eff_proton;
      bool rec_selection_minus = false;
      bool gen_selection_minus = false;
      bool rec_selection_plus = false;
      bool gen_selection_plus = false;

      //reweight t-slope
      double t_slope_up = 4.50509+0.221904;
      double t_slope_down = 4.50509-0.221904;
      double reweight_t_slope_up = exp(-t_slope_up*fabs(t_minus_proton_gen))/exp(-4.5050*fabs(t_minus_proton_gen));
      double reweight_t_slope_dw = exp(-t_slope_down*fabs(t_minus_proton_gen))/exp(-4.5050*fabs(t_minus_proton_gen));
      
      if (jet_gen) {
          ++nevents_jet_gen; 
          histosTH1F["log_x_parton_minus_gen"]->Fill( log10(x_minus_gen), reweigth_beta_minus  );
          histosTH1F["log_x_parton_plus_gen"]->Fill( log10(x_plus_gen), reweigth_beta_plus  );
          histosTH1F["xi_cms_gen"]->Fill( xi_minus_gen, 1.);
          histosTH1F["xi_cms_gen_bin"]->Fill( xi_minus_gen, 1.);
          histosTH1F["xi_gen_proton_minus_binsasha"]->Fill( xi_minus_proton_gen, 1. );
          if (proton_eta_minus>-3) histosTH1F["xi_cms_gen_gap_bin"]->Fill( xi_minus_gen, 1.);
      }

      if (jet_rec && select_Vertex) {
         ++nevents_jet_rec;
         histosTH1F["log_x_parton_minus_rec"]->Fill( log10(x_minus_rec), event_weight_rec_minus  );
         histosTH1F["log_x_parton_plus_rec"]->Fill( log10(x_plus_rec), event_weight_rec_plus  );
         histosTH1F["log_x_parton_gen_recsel"]->Fill( log10(x_minus_gen), event_weight_rec_minus  );
         histosTH1F["xi_cms_rec"]->Fill( xi_minus_Reco, event_weight*eff_trigger );
         histosTH1F["xi_cms_rec_bin"]->Fill( xi_minus_Reco, event_weight*eff_trigger );
         histosTH1F["xi_rec_proton_minus_binsasha"]->Fill( xi_minus_proton_smear_rec , event_weight*eff_trigger );
         if (eta_min>-3) histosTH1F["xi_cms_rec_gap_bin"]->Fill( xi_minus_Reco, event_weight*eff_trigger );
      }
      
      if (jet_gen && proton_kinec_accept_t_minus_gen && proton_kinec_accept_xi_minus_gen ){
          ++nevents_proton_gen;
          gen_selection_minus = true;
          //histosTH1F["t_gen_proton_signal_right_kint_kinxi_cut_slopetplus"]->Fill( t_random_slopeplus, reweigth_beta );
          //histosTH1F["t_gen_proton_signal_right_kint_kinxi_cut_slopetminus"]->Fill( t_random_slopeminus, reweigth_beta );
          histosTH1F["t_gen_proton_signal_kint_kinxi_cut"]->Fill( fabs(t_minus_proton_gen), reweigth_beta_minus );
          histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin"]->Fill( fabs(t_minus_proton_gen), reweigth_beta_minus );
          histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin3"]->Fill( fabs(t_minus_proton_gen), reweigth_beta_minus );
          histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin4"]->Fill( fabs(t_minus_proton_gen), reweigth_beta_minus );
          histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin_slopeup"]->Fill( fabs(t_minus_proton_gen), reweigth_beta_minus*reweight_t_slope_up );
          histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin_slopedw"]->Fill( fabs(t_minus_proton_gen), reweigth_beta_minus*reweight_t_slope_dw );
          histosTH1F["leadingJet_pt_gen_signal_kint_kinxi_cut"]->Fill( leadingJet_pt_gen, reweigth_beta_minus  );
          histosTH1F["secondJet_pt_gen_signal_kint_kinxi_cut"]->Fill( secondJet_pt_gen, reweigth_beta_minus );
          histosTH1F["leadingJet_eta_gen_signal_kint_kinxi_cut"]->Fill( Jet1_eta_gen, reweigth_beta_minus  );
          histosTH1F["secondJet_eta_gen_signal_kint_kinxi_cut"]->Fill( Jet2_eta_gen, reweigth_beta_minus );
          histosTH1F["EtaJet_average_gen_signal_kint_kinxi_cut"]->Fill( eta_average_gen, reweigth_beta_minus );
          histosTH1F["massJet_gen_signal_kint_kinxi_cut"]->Fill( mass_jets_gen, reweigth_beta_minus );
          histosTH1F["mass_x_gen_signal_kint_kinxi_cut"]->Fill( mass_x_minus_gen, reweigth_beta_minus );
          histosTH1F["xi_gen_proton_signal_kint_kinxi_cut"]->Fill( xi_minus_proton_gen , reweigth_beta_minus );
          histosTH1F["xi_gen_proton_minus_signal_kint_kinxi_cut_bin"]->Fill( xi_minus_proton_gen , reweigth_beta_minus );
          histosTH1F["xi_gen_proton_minus_signal_kint_kinxi_cut_bin_slopeup"]->Fill( xi_minus_proton_gen , reweigth_beta_minus*reweight_t_slope_up );
          histosTH1F["xi_gen_proton_minus_signal_kint_kinxi_cut_bin_slopedw"]->Fill( xi_minus_proton_gen , reweigth_beta_minus*reweight_t_slope_dw );
          histosTH1F["beta_gen_proton_minus_signal_kint_kinxi_cut"]->Fill( proton_beta_minus_gen , reweigth_beta_minus );
          histosTH1F["log_x_parton_minus_gen_signal_kint_kinxi_cut"]->Fill( log10(x_minus_gen) , reweigth_beta_minus );
          histosTH1F["log_x_parton_minus_gen_signal_kint_kinxi_cut_slopeup"]->Fill( log10(x_minus_gen) , reweigth_beta_minus*reweight_t_slope_up );
          histosTH1F["log_x_parton_minus_gen_signal_kint_kinxi_cut_slopedw"]->Fill( log10(x_minus_gen) , reweigth_beta_minus*reweight_t_slope_dw );
          if ( !xi_minus_proton_smear_rec && xi_minus_proton_gen ) response_xi_minus.Miss (xi_minus_proton_gen, reweigth_beta_minus);
          if ( !xi_minus_proton_smear_rec && xi_minus_proton_gen ) response_xi_minus_test.Miss (xi_minus_proton_gen, reweigth_beta_minus);
          if ( !fabs(t_minus_proton_smear_rec) && fabs(t_minus_proton_gen) ) response_t_minus.Miss (fabs(t_minus_proton_gen), reweigth_beta_minus);
          if ( !proton_beta_minus_rec && proton_beta_minus_gen ) response_beta_minus.Miss (proton_beta_minus_gen, reweigth_beta_minus);
          if ( !log10(x_minus_rec) && log10(x_minus_gen) ) response_x_minus.Miss (log10(x_minus_gen), reweigth_beta_minus);
          if ( !log10(x_minus_rec) && log10(x_minus_gen) ) response_x_minus_test.Miss (log10(x_minus_gen), reweigth_beta_minus);
          if ( !fabs(t_minus_proton_smear_rec) && fabs(t_minus_proton_gen) ) response_t_minus_test.Miss (fabs(t_minus_proton_gen), reweigth_beta_minus);
      }

      if (jet_gen && proton_kinec_accept_t_plus_gen && proton_kinec_accept_xi_plus_gen ){
          histosTH1F["t_gen_proton_plus_signal_kint_kinxi_cut_bin"]->Fill( fabs(t_plus_proton_gen), reweigth_beta_plus );
          histosTH1F["t_gen_proton_plus_signal_kint_kinxi_cut_bin4"]->Fill( fabs(t_plus_proton_gen), reweigth_beta_plus );
          histosTH1F["xi_gen_proton_plus_signal_kint_kinxi_cut_bin"]->Fill( xi_plus_proton_gen , reweigth_beta_plus );
          histosTH1F["beta_gen_proton_plus_signal_kint_kinxi_cut"]->Fill( proton_beta_plus_gen , reweigth_beta_plus );
          histosTH1F["log_x_parton_plus_gen_signal_kint_kinxi_cut"]->Fill( log10(x_plus_gen) , reweigth_beta_plus );
          if ( !xi_plus_proton_smear_rec && xi_plus_proton_gen ) response_xi_plus.Miss (xi_plus_proton_gen, reweigth_beta_plus);
          if ( !fabs(t_plus_proton_smear_rec) && fabs(t_plus_proton_gen) ) response_t_plus.Miss (fabs(t_plus_proton_gen), reweigth_beta_plus);
          if ( !proton_beta_plus_rec && proton_beta_plus_gen ) response_beta_plus.Miss (proton_beta_plus_gen, reweigth_beta_plus);
          if ( !log10(x_plus_rec) && log10(x_plus_gen) ) response_x_plus.Miss (log10(x_plus_gen), reweigth_beta_plus);
      } 
      if (proton_kinec_accept_t_minus_gen && proton_kinec_accept_xi_minus_gen ){
          histosTH1F["t_gen_proton_kint_kinxi_bin_nojet"]->Fill( fabs(t_minus_proton_gen), reweigth_beta_minus );
          histosTH1F["xi_gen_proton_kint_kinxi_bin_nojet"]->Fill( xi_minus_proton_gen , reweigth_beta_minus );
      }
      if (jet_pt20_gen ){
          histosTH1F["xi_cms_gen_pt20"]->Fill( xi_minus_gen , event_weight );
          histosTH1F["xi_gen_pt20_proton_signal_kinxi_cut"]->Fill( xi_minus_proton_gen , 1. );
          histosTH1F["xi_gen_pt20_proton_signal_kinxi_cut_bin"]->Fill( xi_minus_proton_gen , 1. );
      }
      if (jet_pt20_rec && select_Vertex )histosTH1F["xi_cms_rec_pt20"]->Fill( xi_minus_Reco , event_weight_rec_minus );
      
      
      if (jet_rec && select_Vertex && proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec)++nevents_jetrec_protkintkinxi;
      if (jet_rec && select_Vertex && proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec && rp_minus_sel)++nevents_jetrec_protkintkinxi_rp;
      if (jet_rec && select_Vertex && proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec && rp_minus_sel && xi_cms_totem_minus_cut)++nevents_jetrec_protkintkinxi_rp_cut;

      if (rp_minus_sel && proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec){
          histosTH1F["t_rec_proton_kint_kinxi_bin_nojet"]->Fill( fabs(t_minus_proton_smear_rec), event_weight*reweigth_beta_minus*eff_proton );
          histosTH1F["xi_rec_proton_kint_kinxi_bin_nojet"]->Fill( xi_minus_proton_smear_rec , event_weight*reweigth_beta_minus*eff_proton );
      }

      if (!proton_plus_rp_accept_mc && proton_minus_rp_accept_mc_smear && jet_rec && select_Vertex && proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec && xi_cms_totem_minus_cut){
          histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_fidcut"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
          histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_fidcut"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus );
          histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_fidcut"]->Fill( log10(x_minus_rec) , event_weight_rec_minus );
      }


      if (jet_gen && xi_minus_proton_gen<=0.1 && fabs(t_minus_proton_gen_gauss)>=0.03 && fabs(t_minus_proton_gen_gauss)<=1. ){
         histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_smear"]->Fill( fabs(t_minus_proton_gen_gauss) , 1.);//reweigth_beta_minus);
         histosTH1F["t_gen_proton_minus_signal_kint_kinxi_cut_bin_smear"]->Fill( fabs(t_minus_proton_gen_gauss) , reweigth_beta_minus);
         if ( !fabs(t_minus_proton_rec_gauss) && fabs(t_minus_proton_gen_gauss) ) response_t_minus_gauss.Miss (fabs(t_minus_proton_gen_gauss), 1.);//reweigth_beta_minus);
      } 

      if (rp_minus_sel && jet_rec && select_Vertex && xi_minus_proton_rec_gauss<=0.1 && fabs(t_minus_proton_rec_gauss)>=0.03 && fabs(t_minus_proton_rec_gauss)<=1. && xi_minus_cms-xi_minus_proton_rec_gauss<=0){
         histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_smear"]->Fill( xi_minus_proton_rec_gauss , event_weight_rec_minus );
         histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_smear"]->Fill( log10(x_minus_rec) , event_weight_rec_minus );
         histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_smear"]->Fill( fabs(t_minus_proton_rec_gauss) , event_weight_rec_minus );
         histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_smear"]->Fill( fabs(t_minus_proton_rec_gauss) , 1.);//event_weight_rec_minus );
         if (fabs(t_minus_proton_rec_gauss) && !fabs(t_minus_proton_gen_gauss) ) response_t_minus_gauss.Fake (fabs(t_minus_proton_rec_gauss), 1.);//event_weight_rec_minus);
         if (jet_gen && xi_minus_proton_gen<=0.1 && fabs(t_minus_proton_gen_gauss)>=0.03 && fabs(t_minus_proton_gen_gauss)<=1. )
             response_t_minus_gauss.Fill ( fabs(t_minus_proton_rec_gauss), fabs(t_minus_proton_gen_gauss), 1.);
      }



 
       // reconstructed proton from MC in the RP acceptance
      if (rp_minus_sel && jet_rec && select_Vertex){
          histosTH2F["beta_rec_gen"]->Fill( proton_beta_minus_rec, proton_beta_minus_gen, event_weight_rec_minus );
          histosTH2F["log_x_parton_rec_gen"]->Fill( log10(x_minus_rec), log10(x_minus_gen) , event_weight_rec_minus );
          histosTH2F["xi_totem_rec_gen"]->Fill( xi_minus_proton_smear_rec, xi_minus_proton_gen , event_weight_rec_minus );
          histosTH1F["t_rec_proton_signal"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus );
          histosTH1F["beta_rec_proton_signal"]->Fill( proton_beta_minus_rec , event_weight_rec_minus );
          histosTH1F["xi_rec_proton_signal"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
          histosTH1F["xi_cms_totem_rec_signal"]->Fill( xi_minus_cms-xi_minus_proton_smear_rec , event_weight_rec_minus );
          histosTH1F["t_gen_proton_signal_recsel"]->Fill( fabs(t_minus_proton_gen) , event_weight_rec_minus );
          histosTH1F["beta_gen_proton_signal_recsel"]->Fill( proton_beta_minus_gen , event_weight_rec_minus );
          histosTH1F["xi_gen_proton_signal_recsel"]->Fill( xi_minus_proton_gen , event_weight_rec_minus );
          ///kinematic region
          if (proton_kinec_accept_t_minus_rec){
              histosTH1F["xi_rec_proton_signal_kint"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
              histosTH1F["xi_rec_proton_signal_kint_bin25"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
              histosTH1F["xi_gen_proton_signal_kint_recsel"]->Fill( xi_minus_proton_gen , event_weight_rec_minus );}
          if (proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec && xi_minus_proton_smear_rec>=0.03) histosTH1F["xi_cms_totem_rec_signal_kin"]->Fill( xi_minus_cms-xi_minus_proton_smear_rec , event_weight_rec_minus );
          if (proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec){
              histosTH1F["xi_rec_proton_signal_kinxi_kint"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
              histosTH1F["xi_cms_rec_proton_signal_kint"]->Fill( xi_minus_cms , event_weight_rec_minus );
              histosTH1F["xi_cms_totem_rec_signal_kinxi_kint"]->Fill( xi_minus_cms-xi_minus_proton_smear_rec , event_weight_rec_minus );
              //xi_cms_st = xi_cms_minus;
              //xi_totem_st = xi_proton_minus_rec;
              //xi_cms_minus_totem_st = xi_cms_minus-xi_proton_minus_rec;
             // small_tree->Fill();
          }
          ///xi_cms -xi_totem cut 
          if (proton_kinec_accept_t_minus_rec && xi_cms_totem_minus_cut){ 
              histosTH1F["xi_rec_proton_signal_kint_cut"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
              histosTH1F["xi_gen_proton_signal_kint_cut_recsel"]->Fill( xi_minus_proton_gen , event_weight_rec_minus);}

          if (proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec && xi_minus_proton_smear_rec>=0.03 && xi_cms_totem_minus_cut){
              histosTH1F["xi_rec_proton_signal_kin_cut"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
              histosTH1F["t_rec_proton_signal_kin_cut"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus );
              histosTH1F["beta_rec_proton_signal_kin_cut"]->Fill( proton_beta_minus_rec , event_weight_rec_minus );
              histosTH1F["leadingJet_pt_rec_signal_kin_cut"]->Fill( Jet1_pt_rec, event_weight_rec_minus  );
              histosTH1F["secondJet_pt_rec_signal_kin_cut"]->Fill( Jet2_pt_rec, event_weight_rec_minus  );
              histosTH1F["log_x_parton_rec_signal_kin_cut"]->Fill( log10(x_minus_rec) , event_weight_rec_minus );
              histosTH1F["EtaJet_average_rec_signal_kin_cut"]->Fill( eta_average_rec, event_weight_rec_minus  );
              histosTH1F["DeltaPtJet_rec_signal_kin_cut"]->Fill( deltapt_rec, event_weight_rec_minus  );
              histosTH1F["DeltaEtaJet_rec_signal_kin_cut"]->Fill( deltaeta_rec, event_weight_rec_minus  );
              histosTH1F["DeltaPhiJet_rec_signal_kin_cut"]->Fill( deltaphi_rec, event_weight_rec_minus  );

              //gen
              if (jet_gen && proton_kinec_accept_t_minus_gen && proton_kinec_accept_xi_minus_gen && xi_minus_proton_gen>=0.03){
                 histosTH1F["xi_gen_proton_signal_kin_cut_recsel"]->Fill( xi_minus_proton_gen , event_weight_rec_minus );
                 histosTH1F["t_gen_proton_signal_kin_cut_recsel"]->Fill( fabs(t_minus_proton_gen) , event_weight_rec_minus );
                 histosTH1F["beta_gen_proton_signal_kin_cut_recsel"]->Fill( proton_beta_minus_gen , event_weight_rec_minus );
                 histosTH1F["log_x_parton_gen_signal_kin_cut_recsel"]->Fill( log10(x_minus_gen) , event_weight_rec_minus );
                 histosTH1F["leadingJet_pt_gen_signal_kin_cut_recsel"]->Fill( leadingJet_pt_gen, event_weight_rec_minus  );
                 histosTH1F["secondJet_pt_gen_signal_kin_cut_recsel"]->Fill( secondJet_pt_gen, event_weight_rec_minus );
                 histosTH1F["EtaJet_average_gen_signal_kin_cut_recsel"]->Fill( eta_average_gen, event_weight_rec_minus );
                 histosTH1F["DeltaPtJet_gen_signal_kin_cut_recsel"]->Fill( deltapt_gen, event_weight_rec_minus  );
                 histosTH1F["DeltaEtaJet_gen_signal_kin_cut_recsel"]->Fill( deltaeta_gen, event_weight_rec_minus  );
                 histosTH1F["DeltaPhiJet_gen_signal_kin_cut_recsel"]->Fill( deltaphi_gen, event_weight_rec_minus  );
              }
          }
          histosTH1F["xi_rec_proton_signal_kinxi"]->Fill( xi_minus_proton_smear_rec, event_weight);
          if (xi_cms_totem_minus_cut) histosTH1F["xi_rec_proton_signal_kinxi_cut"]->Fill( xi_minus_proton_smear_rec, event_weight_rec_minus);
          if (proton_kinec_accept_t_minus_rec && proton_kinec_accept_xi_minus_rec && xi_cms_totem_minus_cut){
              ++nevents_proton_rec;
              rec_selection_minus = true;
         histosTH1F["t_proton"]->Fill(fabs(t_minus_proton_gen), event_weight);
         histosTH1F["t_proton_smear"]->Fill(fabs(t_minus_proton_smear_gen), event_weight);
              histosTH1F["xi_cms_rec_proton_signal_kint_cut"]->Fill( xi_minus_cms , event_weight_rec_minus );
              histosTH1F["xi_rec_proton_signal_kint_kinxi_cut"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
              histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
              histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_triggereff"]->Fill( xi_minus_proton_smear_rec , event_weight*eff_proton );
              histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_protoneff"]->Fill( xi_minus_proton_smear_rec , event_weight*eff_proton );
              histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_zvtx"]->Fill( xi_minus_proton_smear_rec , event_weight*zvtx );
              histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_slopeup"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus*reweight_t_slope_up );
              histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin_slopedw"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus*reweight_t_slope_dw );
              histosTH1F["xi_rec_proton_signal_kint_kinxi_cut_bin15"]->Fill( xi_minus_proton_smear_rec , event_weight_rec_minus );
              histosTH1F["t_rec_proton_signal_kint_kinxi_cut"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus );
              histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus );
              histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_triggereff"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight*eff_proton );
              histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_protoneff"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight*eff_proton );
              histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_zvtx"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight*zvtx );
              histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin3"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus );
              histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin4"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus );
              histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_slopeup"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus*reweight_t_slope_up );
              histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin_slopedw"]->Fill( fabs(t_minus_proton_smear_rec) , event_weight_rec_minus*reweight_t_slope_dw );
              histosTH1F["beta_rec_proton_minus_signal_kint_kinxi_cut"]->Fill( proton_beta_minus_rec , event_weight_rec_minus );
              histosTH1F["beta_rec_proton_minus_signal_kint_kinxi_cut_triggereff"]->Fill( proton_beta_minus_rec, event_weight*eff_proton );
              histosTH1F["beta_rec_proton_minus_signal_kint_kinxi_cut_protoneff"]->Fill( proton_beta_minus_rec , event_weight*eff_proton );
              histosTH1F["beta_rec_proton_minus_signal_kint_kinxi_cut_zvtx"]->Fill( proton_beta_minus_rec , event_weight*zvtx );
              histosTH1F["leadingJet_pt_rec_signal_kint_kinxi_cut"]->Fill( Jet1_pt_rec, event_weight_rec_minus  );
              histosTH1F["leadingJet_pt_rec_signal_kint_kinxi_cut_triggereff"]->Fill( Jet1_pt_rec, event_weight*eff_proton );
              histosTH1F["leadingJet_pt_rec_signal_kint_kinxi_cut_protoneff"]->Fill( Jet1_pt_rec , event_weight*eff_proton );
              histosTH1F["leadingJet_pt_rec_signal_kint_kinxi_cut_zvtx"]->Fill( Jet1_pt_rec , event_weight*zvtx );
              histosTH1F["secondJet_pt_rec_signal_kint_kinxi_cut"]->Fill( Jet2_pt_rec, event_weight_rec_minus  );
              histosTH1F["secondJet_pt_rec_signal_kint_kinxi_cut_triggereff"]->Fill( Jet2_pt_rec, event_weight*eff_proton );
              histosTH1F["secondJet_pt_rec_signal_kint_kinxi_cut_protoneff"]->Fill( Jet2_pt_rec , event_weight*eff_proton );
              histosTH1F["secondJet_pt_rec_signal_kint_kinxi_cut_zvtx"]->Fill( Jet2_pt_rec , event_weight*zvtx );
              histosTH1F["leadingJet_eta_rec_signal_kint_kinxi_cut"]->Fill( Jet1_eta_rec, event_weight_rec_minus  );
              histosTH1F["secondJet_eta_rec_signal_kint_kinxi_cut"]->Fill( Jet2_eta_rec, event_weight_rec_minus  );
              histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut"]->Fill( log10(x_minus_rec) , event_weight_rec_minus );
              histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_triggereff"]->Fill( log10(x_minus_rec), event_weight*eff_proton );
              histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_protoneff"]->Fill( log10(x_minus_rec) , event_weight*eff_proton );
              histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_zvtx"]->Fill( log10(x_minus_rec) , event_weight*zvtx );
              histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_slopeup"]->Fill( log10(x_minus_rec) , event_weight_rec_minus*reweight_t_slope_up );
              histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut_slopedw"]->Fill( log10(x_minus_rec) , event_weight_rec_minus*reweight_t_slope_dw );
              histosTH1F["EtaJet_average_rec_signal_kint_kinxi_cut"]->Fill( eta_average_rec, event_weight_rec_minus  );
              histosTH1F["EtaJet_average_rec_signal_kint_kinxi_cut_triggereff"]->Fill( eta_average_rec, event_weight*eff_proton );
              histosTH1F["EtaJet_average_rec_signal_kint_kinxi_cut_protoneff"]->Fill( eta_average_rec, event_weight*eff_proton );
              histosTH1F["EtaJet_average_rec_signal_kint_kinxi_cut_zvtx"]->Fill( eta_average_rec, event_weight*zvtx );
              histosTH1F["massJet_rec_signal_kint_kinxi_cut"]->Fill( mass_jets_rec, event_weight_rec_minus  );
              histosTH1F["mass_x_rec_signal_kint_kinxi_cut"]->Fill( mass_x_minus_rec, reweigth_beta_minus );
              histosTH1F["t_rec_gen_minus"]->Fill( fabs(t_minus_proton_smear_rec)-fabs(t_minus_proton_gen) , event_weight_rec_minus );
              histosTH1F["xi_rec_gen_minus"]->Fill( xi_minus_proton_smear_rec-xi_minus_proton_gen , event_weight_rec_minus );
              histosTH1F["log_x_parton_minus_rec_gen"]->Fill( log10(x_minus_rec)-log10(x_minus_gen), event_weight_rec_minus);
              histosTH2F["log_x_parton_resol"]->Fill( log10(x_minus_rec), log10(x_minus_rec)-log10(x_minus_gen), 1.);
              histosTH1F["thx_proton_kint_kinxi_cut"]->Fill( thx_minus_proton, event_weight_rec_minus );
              histosTH1F["thy_proton_kint_kinxi_cut"]->Fill( thy_minus_proton, event_weight_rec_minus );
              if ( xi_minus_proton_smear_rec && !xi_minus_proton_gen ) response_xi_minus.Fake (xi_minus_proton_smear_rec, event_weight_rec_minus);
              if ( xi_minus_proton_smear_rec && !xi_minus_proton_gen ) response_xi_minus_test.Fake (xi_minus_proton_smear_rec, event_weight_rec_minus);
              if ( fabs(t_minus_proton_smear_rec) && !fabs(t_minus_proton_gen) ) response_t_minus.Fake (fabs(t_minus_proton_smear_rec), event_weight_rec_minus);
              if ( proton_beta_minus_rec && !proton_beta_minus_gen ) response_beta_minus.Fake (proton_beta_minus_rec, event_weight_rec_minus);
              if ( log10(x_minus_rec) && !log10(x_minus_gen) ) response_x_minus.Fake (log10(x_minus_rec), event_weight_rec_minus);
              if ( log10(x_minus_rec) && !log10(x_minus_gen) ) response_x_minus_test.Fake (log10(x_minus_rec), event_weight_rec_minus);
              if ( fabs(t_minus_proton_smear_rec) && !fabs(t_minus_proton_gen) ) response_t_minus_test.Fake (fabs(t_minus_proton_smear_rec), event_weight_rec_minus);

              //gen
              if (jet_gen && proton_kinec_accept_t_minus_gen && proton_kinec_accept_xi_minus_gen ){
                 response_xi_minus.Fill ( xi_minus_proton_smear_rec, xi_minus_proton_gen, 1.);
                 response_xi_minus_test.Fill ( xi_minus_proton_smear_rec, xi_minus_proton_gen, 1.);
                 response_t_minus.Fill ( fabs(t_minus_proton_smear_rec), fabs(t_minus_proton_gen), 1.);
                 response_beta_minus.Fill (proton_beta_minus_rec, proton_beta_minus_gen, 1.);
                 response_x_minus.Fill (log10(x_minus_rec), log10(x_minus_gen), 1.);
                 response_x_minus_test.Fill (log10(x_minus_rec), log10(x_minus_gen), 1.);
                 response_t_minus_test.Fill ( fabs(t_minus_proton_smear_rec), fabs(t_minus_proton_gen), 1. );
                 histosTH1F["xi_gen_proton_signal_kint_kinxi_cut_recsel"]->Fill( xi_minus_proton_gen , event_weight_rec_minus );
                 histosTH1F["t_gen_proton_signal_kint_kinxi_cut_recsel"]->Fill( fabs(t_minus_proton_gen) , event_weight_rec_minus );
                 histosTH1F["beta_gen_proton_signal_kint_kinxi_cut_recsel"]->Fill( proton_beta_minus_gen , event_weight_rec_minus );
                 histosTH1F["log_x_parton_gen_signal_kint_kinxi_cut_recsel"]->Fill( log10(x_minus_gen) , event_weight_rec_minus );
                 histosTH1F["leadingJet_pt_gen_signal_kint_kinxi_cut_recsel"]->Fill( leadingJet_pt_gen, event_weight_rec_minus  );
                 histosTH1F["secondJet_pt_gen_signal_kint_kinxi_cut_recsel"]->Fill( secondJet_pt_gen, event_weight_rec_minus );
                 histosTH1F["EtaJet_average_gen_signal_kint_kinxi_cut_recsel"]->Fill( eta_average_gen, event_weight_rec_minus );
              }
          }
      }
      if (rp_plus_sel && jet_rec && select_Vertex){
          if (proton_kinec_accept_t_plus_rec && proton_kinec_accept_xi_plus_rec && xi_cms_totem_plus_cut){
              rec_selection_plus = true;
              histosTH1F["xi_rec_proton_plus_signal_kint_kinxi_cut_bin"]->Fill( xi_plus_proton_smear_rec , event_weight_rec_plus );
              histosTH1F["t_rec_proton_plus_signal_kint_kinxi_cut_bin"]->Fill( fabs(t_plus_proton_smear_rec) , event_weight_rec_plus );
              histosTH1F["t_rec_proton_plus_signal_kint_kinxi_cut_bin4"]->Fill( fabs(t_plus_proton_smear_rec) , event_weight_rec_plus );
              histosTH1F["beta_rec_proton_plus_signal_kint_kinxi_cut"]->Fill( proton_beta_plus_rec , event_weight_rec_plus );
              histosTH1F["log_x_parton_plus_rec_signal_kint_kinxi_cut"]->Fill( log10(x_plus_rec) , event_weight_rec_plus );
              histosTH1F["t_rec_gen_plus"]->Fill( fabs(t_plus_proton_smear_rec)-fabs(t_plus_proton_gen) , event_weight_rec_plus );
              histosTH1F["xi_rec_gen_plus"]->Fill( xi_plus_proton_smear_rec-xi_plus_proton_gen , event_weight_rec_plus );
              histosTH1F["log_x_parton_plus_rec_gen"]->Fill( log10(x_plus_rec)-log10(x_plus_gen) , event_weight_rec_plus );
              if ( xi_plus_proton_smear_rec && !xi_plus_proton_gen ) response_xi_plus.Fake (xi_plus_proton_smear_rec, event_weight_rec_plus);
              if ( fabs(t_plus_proton_smear_rec) && !fabs(t_plus_proton_gen) ) response_t_plus.Fake (fabs(t_plus_proton_smear_rec), event_weight_rec_plus);
              if ( proton_beta_plus_rec && !proton_beta_plus_gen ) response_beta_plus.Fake (proton_beta_plus_rec, event_weight_rec_plus);
              if ( log10(x_plus_rec) && !log10(x_plus_gen) ) response_x_plus.Fake (log10(x_plus_rec), event_weight_rec_plus);

              if (jet_gen && proton_kinec_accept_t_plus_gen && proton_kinec_accept_xi_plus_gen ){
                 response_xi_plus.Fill ( xi_plus_proton_smear_rec, xi_plus_proton_gen, 1.);
                 response_t_plus.Fill ( fabs(t_plus_proton_smear_rec), fabs(t_plus_proton_gen), 1.);
                 response_beta_plus.Fill (proton_beta_plus_rec, proton_beta_plus_gen, 1.);
                 response_x_plus.Fill (log10(x_plus_rec), log10(x_plus_gen), 1.);
              }
          }
      }


    }//end loop for events
   // cout <<"After the jet selection " << nevents_jets << " events  "<< endl;
   // cout <<"After GenPart selection " << nevents_gen << " events "<< endl;
   // cout <<"After PF selection " << nevents_pf << " events "<< endl;
   // cout <<"  "<< endl;
   file->Close();
 
  }//end of loop over files
//   cout <<"After selection " << nevents_accepted << endl;

     
  //output file
  TFile* output = new TFile(outputFileName,"RECREATE");
  output->cd();

  event_selection->SetBinContent(1, nevents_total); 
  event_selection->SetBinContent(2, nevents_vtx); 
  event_selection->SetBinContent(3, nevents_pf); 
  event_selection->SetBinContent(4, nevents_jet_rec); 
  event_selection->SetBinContent(5, nevents_jet_gen); 
  event_selection->SetBinContent(6, nevents_proton_gen); 
  event_selection->SetBinContent(7, nevents_jetrec_protkintkinxi); 
  event_selection->SetBinContent(8, nevents_jetrec_protkintkinxi_rp); 
  event_selection->SetBinContent(9, nevents_jetrec_protkintkinxi_rp_cut); 
  event_selection->GetXaxis()->SetBinLabel(1, "total");
  event_selection->GetXaxis()->SetBinLabel(2, "vtx");
  event_selection->GetXaxis()->SetBinLabel(3, "pf");
  event_selection->GetXaxis()->SetBinLabel(4, "jetrec");
  event_selection->GetXaxis()->SetBinLabel(5, "jetgen");
  event_selection->GetXaxis()->SetBinLabel(6, "jetgen-protongen");
  event_selection->GetXaxis()->SetBinLabel(7, "jetrec-protonrec");
  event_selection->GetXaxis()->SetBinLabel(8, "jetrec-protonrec-rp");
  event_selection->GetXaxis()->SetBinLabel(9, "jetrec-protonrec-rp-cut");
  event_selection->Write();


  //unfolding Bayes
  TH1F* log_x_right_proton_kint_kinxi_data = (TH1F*)data->Get("log_x_minus_kint_kinxi");
  TH1F* log_x_left_proton_kint_kinxi_data = (TH1F*)data->Get("log_x_plus_kint_kinxi");
  TH1F* xi_right_proton_kint_kinxi_data = (TH1F*)data->Get("proton_right_xi_signal_kint_kinxi_bin");
  TH1F* xi_left_proton_kint_kinxi_data = (TH1F*)data->Get("proton_left_xi_signal_kint_kinxi_bin");
  TH1F* t_right_proton_kint_kinxi_data = (TH1F*)data->Get("proton_right_t_signal_kint_kinxi_bin");
  TH1F* t_left_proton_kint_kinxi_data = (TH1F*)data->Get("proton_left_t_signal_kint_kinxi_bin");


  //if (side_minus){
//     RooUnfoldBayes x_jjp_minus_unfold_bayes (&response_x_minus, log_x_right_proton_kint_kinxi_data, 4);
//     RooUnfoldBayes x_jjp_minus_unfold_test (&response_x_minus_test, histosTH1F["log_x_parton_minus_rec_signal_kint_kinxi_cut"], 4);
//     RooUnfoldBayes xi_jjp_minus_unfold_bayes (&response_xi_minus, xi_right_proton_kint_kinxi_data, 4);
//     TH1F* x_jjp_minus_data_unfolded_bayes = (TH1F*) x_jjp_minus_unfold_bayes.Hreco();
//     TH1F* x_jjp_minus_data_unfolded_test = (TH1F*) x_jjp_minus_unfold_test.Hreco();
//     TH1F* xi_jjp_minus_data_unfolded_bayes = (TH1F*) xi_jjp_minus_unfold_bayes.Hreco();
//     x_jjp_minus_data_unfolded_bayes->Write();
//     x_jjp_minus_data_unfolded_test->Write();
//     xi_jjp_minus_data_unfolded_bayes->Write();
  //}
/*  if (side_plus){    
     RooUnfoldBayes x_jjp_plus_unfold_bayes (&response_x_plus, log_x_left_proton_kint_kinxi_data, 4);
     RooUnfoldBayes xi_jjp_plus_unfold_bayes (&response_xi_plus, xi_left_proton_kint_kinxi_data, 4);
     RooUnfoldBayes t_jjp_plus_unfold_bayes (&response_t_plus, t_left_proton_kint_kinxi_data, 4);
     TH1F* x_jjp_plus_data_unfolded_bayes = (TH1F*) x_jjp_plus_unfold_bayes.Hreco();
     TH1F* xi_jjp_plus_data_unfolded_bayes = (TH1F*) xi_jjp_plus_unfold_bayes.Hreco();
     TH1F* t_jjp_plus_data_unfolded_bayes = (TH1F*) t_jjp_plus_unfold_bayes.Hreco();
     x_jjp_plus_data_unfolded_bayes->Write();
     xi_jjp_plus_data_unfolded_bayes->Write();
     t_jjp_plus_data_unfolded_bayes->Write();
  }
*/
  float cross_section;
  if (mc == "pomwig" && !reggeon && side_minus && side_plus) cross_section = 2*2.1189e7; //pb cross section for pomwig-pomeron
  if (mc == "pomwig" && !reggeon && ((side_minus && !side_plus)||(!side_minus && side_plus))) cross_section = 2.1189e7; //pb cross section for pomwig-pomeron
  if (mc == "pomwig" && reggeon && side_minus && side_plus) cross_section = 2*5.9695e6;//pb cross section for pomwig-reggeon
  if (mc == "pomwig" && reggeon && ((side_minus && !side_plus)||(!side_minus && side_plus))) cross_section = 5.9695e6;//pb cross section for pomwig-reggeon
  if (mc == "pythia8_diff") cross_section = 2.073e+10; //pb cross section pythia8

  float luminity_HLT_L1Jet1_198902 = 0.015879;//pb ---- luminity for LP_Jets1_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet2_198902 = 0.015879;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet1_198903 = 0.008698;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet2_198903 = 0.008698;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity = luminity_HLT_L1Jet1_198902 + luminity_HLT_L1Jet2_198902 + luminity_HLT_L1Jet1_198903 + luminity_HLT_L1Jet2_198903;
  
  float n_events = luminity*cross_section;
  float f1 = (float) nevents_total;
  // float f1 = (float) nweight_total;
  Double_t scale = n_events/f1;
  cout<<"eventos  "<<nevents_total<<"   pesos "<<nweight_total<<"   cross section: "<<cross_section<< "  scale "<< scale<< endl; 

  small_tree->Write();


  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin(); it_histo != histosTH1F.end(); ++it_histo){
     (*it_histo).second->Scale(scale);
     (*it_histo).second->Write();
  }
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin(); it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();


//     RooUnfoldBayes t_jjp_minus_unfold_bayes (&response_t_minus_gauss, histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_smear"], 4);
//     TH1F* t_jjp_minus_data_unfolded_bayes = (TH1F*) t_jjp_minus_unfold_bayes.Hreco();
//     t_jjp_minus_data_unfolded_bayes->Write();

  /*   RooUnfoldBayes t_jjp_minus_unfold_test (&response_t_minus_test, histosTH1F["t_rec_proton_minus_signal_kint_kinxi_cut_bin"], 4);
     TH1F* t_jjp_minus_data_unfolded_test = (TH1F*) t_jjp_minus_unfold_test.Hreco();
     t_jjp_minus_data_unfolded_test->Write();

     RooUnfoldBayes xi_jjp_minus_unfold_test (&response_xi_minus_test, histosTH1F["xi_rec_proton_minus_signal_kint_kinxi_cut_bin"], 4);
     TH1F* xi_jjp_minus_data_unfolded_test = (TH1F*) xi_jjp_minus_unfold_test.Hreco();
     xi_jjp_minus_data_unfolded_test->Write();

*/
  output->Close();
}
