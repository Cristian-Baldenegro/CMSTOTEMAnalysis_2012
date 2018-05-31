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
#include <exception>

//ROOUNFOLD CLASSES
//#include "/home/lhuertas/afs-work/CMSTOTEM_analysis/MC/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
//#include "/home/lhuertas/afs-work/CMSTOTEM_analysis/MC/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
//#include "/home/lhuertas/afs-work/CMSTOTEM_analysis/MC/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetResolution.h"
#include "CondFormats/JetMETObjects/src/SimpleJetResolution.cc"
#include "smearedJet.h"

#define PI 3.141592653589793
using namespace std;

bool sortByPt(MyBaseJet const& jet1, MyBaseJet const& jet2){
   return ( jet1.Pt() > jet2.Pt() );
}

JetCorrectorParameters *L2Relative, *L3Absolute;
vector<JetCorrectorParameters> vecL2Relative, vecL3Absolute;



void non_diffractive_ntuple(string const& mc = "pythia6", string const& pt_range = "QCD_pt_15_30", const Int_t nevt_max = -1){

  //mc = "pythia6" || "pythia8_nondiff" || "herwig"
  

  bool pt_15_30 = false;
  bool pt_30_50 = false;
  bool pt_50_80 = false;
  bool pt_80_120 = false;
  bool pt_120_170 = false;
  bool pt_170_300 = false;
  bool pt_300_470 = false;

  if (pt_range == "QCD_pt_15_30") pt_15_30 = true;
  if (pt_range == "QCD_pt_30_50") pt_30_50 = true;
  if (pt_range == "QCD_pt_50_80") pt_50_80 = true;
  if (pt_range == "QCD_pt_80_120") pt_80_120 = true;
  if (pt_range == "QCD_pt_120_170") pt_120_170 = true;
  if (pt_range == "QCD_pt_170_300") pt_170_300 = true;
  if (pt_range == "QCD_pt_300_470") pt_300_470 = true;

  int variation = -1;// scale factor variation -> 0=nominal; 1=up, -1=down
  TString sf_var;  
  if (variation == 1) sf_var = "up";
  if (variation == -1) sf_var = "dw";
  if (variation == 0) sf_var = "";
  
  //TString outputFileName = "/storage/lhuertas/uerj-1/mc/Workspace/root_files/" + mc + "_" + pt_range + "_ntuple.root";
  TString outputFileName = (mc == "herwig") ? "/afs/cern.ch/user/l/lhuertas/" + mc + "_" + pt_range + "_ntuple_jetid_jer" + sf_var +  ".root" :
          "/afs/cern.ch/user/l/lhuertas/" + mc + "_ntuple_jetid_jer" + sf_var +  ".root";

  cout<<outputFileName<<endl;

  
  bool verbose = false;
  string treeName = "evt";//"cms_totem";
  string jetCollName = "ak5PFJets";
  string jetCorrName = "ak5PFL2L3Residual";
  jetCorrName = "ak5PFL2L3"; 

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
  
  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;

  TFile* output = new TFile(outputFileName,"RECREATE");
  output->cd();

  double weight, xi_rec_cms_right, xi_rec_cms_left, xi_rec_totem_right, xi_rec_totem_left, t_rec_totem_right, t_rec_totem_left, beta_rec_right, beta_rec_left;
  double xi_gen_cms_right, xi_gen_cms_left, xi_gen_totem_right, xi_gen_totem_left, t_gen_totem_right, t_gen_totem_left, beta_gen_right, beta_gen_left;
  double jet1_rec_pt, jet1_rec_eta, jet1_rec_phi, jet2_rec_pt, jet2_rec_eta, jet2_rec_phi, x_rec_right, x_rec_left; 
  double jet1_gen_pt, jet1_gen_eta, jet1_gen_phi, jet2_gen_pt, jet2_gen_eta, jet2_gen_phi, x_gen_right, x_gen_left; 
  bool rp_right, rp_left, rp_right_top_near, rp_right_top_far, rp_right_bottom_near, rp_right_bottom_far, rp_left_top_near, rp_left_top_far, rp_left_bottom_near, rp_left_bottom_far;
  double xi_true;
  double rp_xpos_024, rp_ypos_024, rp_xpos_025, rp_ypos_025, rp_xpos_124, rp_ypos_124, rp_xpos_125, rp_ypos_125;
  TTree* small_tree = new TTree("small_tree","");
  small_tree->Branch("weight",&weight,"weight/D");
  small_tree->Branch("xi_rec_cms_right",&xi_rec_cms_right,"xi_rec_cms_right/D");
  small_tree->Branch("xi_rec_cms_left",&xi_rec_cms_left,"xi_rec_cms_left/D");
  small_tree->Branch("xi_gen_cms_right",&xi_gen_cms_right,"xi_gen_cms_right/D");
  small_tree->Branch("xi_gen_cms_left",&xi_gen_cms_left,"xi_gen_cms_left/D");
  small_tree->Branch("xi_rec_totem_right",&xi_rec_totem_right,"xi_rec_totem_right/D");
  small_tree->Branch("xi_rec_totem_left",&xi_rec_totem_left,"xi_rec_totem_left/D");
  small_tree->Branch("xi_gen_totem_right",&xi_gen_totem_right,"xi_gen_totem_right/D");
  small_tree->Branch("xi_gen_totem_left",&xi_gen_totem_left,"xi_gen_totem_left/D");
  small_tree->Branch("t_rec_totem_right",&t_rec_totem_right,"t_rec_totem_right/D");
  small_tree->Branch("t_rec_totem_left",&t_rec_totem_left,"t_rec_totem_left/D");
  small_tree->Branch("t_gen_totem_right",&t_gen_totem_right,"t_gen_totem_right/D");
  small_tree->Branch("t_gen_totem_left",&t_gen_totem_left,"t_gen_totem_left/D");
  small_tree->Branch("x_rec_right",&x_rec_right,"x_rec_right/D");
  small_tree->Branch("x_rec_left",&x_rec_left,"x_rec_left/D");
  small_tree->Branch("x_gen_right",&x_gen_right,"x_gen_right/D");
  small_tree->Branch("x_gen_left",&x_gen_left,"x_gen_left/D");
  small_tree->Branch("beta_rec_right",&beta_rec_right,"beta_rec_right/D");
  small_tree->Branch("beta_rec_left",&beta_rec_left,"beta_rec_left/D");
  small_tree->Branch("beta_gen_right",&beta_gen_right,"beta_gen_right/D");
  small_tree->Branch("beta_gen_left",&beta_gen_left,"beta_gen_left/D");
  small_tree->Branch("jet1_rec_pt",&jet1_rec_pt,"jet1_rec_pt/D");
  small_tree->Branch("jet1_rec_eta",&jet1_rec_eta,"jet1_rec_eta/D");
  small_tree->Branch("jet1_rec_phi",&jet1_rec_phi,"jet1_rec_phi/D");
  small_tree->Branch("jet1_gen_pt",&jet1_gen_pt,"jet1_gen_pt/D");
  small_tree->Branch("jet1_gen_eta",&jet1_gen_eta,"jet1_gen_eta/D");
  small_tree->Branch("jet1_gen_phi",&jet1_gen_phi,"jet1_gen_phi/D");
  small_tree->Branch("jet2_rec_pt",&jet2_rec_pt,"jet2_rec_pt/D");
  small_tree->Branch("jet2_rec_eta",&jet2_rec_eta,"jet2_rec_eta/D");
  small_tree->Branch("jet2_rec_phi",&jet2_rec_phi,"jet2_rec_phi/D");
  small_tree->Branch("jet2_gen_pt",&jet2_gen_pt,"jet2_gen_pt/D");
  small_tree->Branch("jet2_gen_eta",&jet2_gen_eta,"jet2_gen_eta/D");
  small_tree->Branch("jet2_gen_phi",&jet2_gen_phi,"jet2_gen_phi/D");
  small_tree->Branch("rp_xpos_024",&rp_xpos_024,"rp_xpos_024/D");
  small_tree->Branch("rp_ypos_024",&rp_ypos_024,"rp_ypos_024/D");
  small_tree->Branch("rp_xpos_025",&rp_xpos_025,"rp_xpos_025/D");
  small_tree->Branch("rp_ypos_025",&rp_ypos_025,"rp_ypos_025/D");
  small_tree->Branch("rp_xpos_124",&rp_xpos_124,"rp_xpos_124/D");
  small_tree->Branch("rp_ypos_124",&rp_ypos_124,"rp_ypos_124/D");
  small_tree->Branch("rp_xpos_125",&rp_xpos_125,"rp_xpos_125/D");
  small_tree->Branch("rp_ypos_125",&rp_ypos_125,"rp_ypos_125/D");
  small_tree->Branch("rp_right",&rp_right,"rp_right/O");
  small_tree->Branch("rp_right_top_near",&rp_right_top_near,"rp_right_top_near/O");
  small_tree->Branch("rp_right_top_far",&rp_right_top_far,"rp_right_top_far/O");
  small_tree->Branch("rp_right_bottom_near",&rp_right_bottom_near,"rp_right_bottom_near/O");
  small_tree->Branch("rp_right_bottom_far",&rp_right_bottom_far,"rp_right_bottom_far/O");
  small_tree->Branch("rp_left",&rp_left,"rp_left/O");
  small_tree->Branch("rp_left_top_near",&rp_left_top_near,"rp_left_top_near/O");
  small_tree->Branch("rp_left_top_far",&rp_left_top_far,"rp_left_top_far/O");
  small_tree->Branch("rp_left_bottom_near",&rp_left_bottom_near,"rp_left_bottom_near/O");
  small_tree->Branch("rp_left_bottom_far",&rp_left_bottom_far,"rp_left_bottom_far/O");
  small_tree->Branch("xi_true",&xi_true,"xi_true/D");


  gStyle->SetPalette(1);

  //===================
  int i_tot = 0 , nevt_tot = 0;

  //MC Files
  TString dirname;
  if (mc == "herwig") dirname = "/storage1/lhuertas/CMSTOTEM/samples/mc/QCD_pt_15_470_Herwig6/" + pt_range + "/";
  if (mc == "pythia6") dirname = "/storage1/lhuertas/CMSTOTEM/samples/mc/Pythia6_QCD/";
  if (mc == "pythia8_nondiff") dirname = "/storage1/lhuertas/CMSTOTEM/samples/mc/Pythia8_QCD_Pt-15to3000_Tune4C_Flat_8TeV/";
  if (mc == "pythia8_nondiff_CUETP8S1") dirname = "/storage1/lhuertas/CMSTOTEM/samples/mc/Pythia8_QCD_Pt-15to3000_CUETP8S1_Flat_8TeV/";
  if (mc == "pythia8_nondiff_CUETP8M1") dirname = "/storage1/lhuertas/CMSTOTEM/samples/mc/Pythia8_QCD_Pt-15to3000_CUETP8M1_Flat_8TeV/";
  TString /*const char*/ ext=".root";
  vector<TString>* vfiles = new vector<TString>; 
  
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

  //Declaration of tree and its branches variables
  TTree* tree = new TTree(treeName.c_str(),"");
  MyEvtId*           evtId        = NULL;
//   MyL1TrigOld*       l1Trig       = NULL;  
//   MyHLTrig*          hltTrig      = NULL;
  vector<MyGenPart>* genPart      = NULL;
  vector<MyTracks>*  track_coll   = NULL;
  vector<MyVertex>*  vertex_coll  = NULL;
  vector<MyPFJet>*   pfJet_coll   = NULL;
  vector<MyPFCand>*  pFlow_coll   = NULL;
  vector<MyGenJet>*   genJet_coll   = NULL;
  MyGenKin*  genKin   = NULL;
  //===============================
  // Jet energy corrections
  //L2Relative = new JetCorrectorParameters("/storage/lhuertas/uerj-1/mc/Workspace/Winter14_V8/Winter14_V8_MC_L2Relative_AK5PF.txt");
  L2Relative = new JetCorrectorParameters("/storage1/lhuertas/CMSTOTEM/mc/Workspace/Winter14_V8/Winter14_V8_MC_L2Relative_AK5PF.txt");
  //L3Absolute = new JetCorrectorParameters("/storage/lhuertas/uerj-1/mc/Workspace/Winter14_V8/Winter14_V8_MC_L3Absolute_AK5PF.txt");
  L3Absolute = new JetCorrectorParameters("/storage1/lhuertas/CMSTOTEM/mc/Workspace/Winter14_V8/Winter14_V8_MC_L3Absolute_AK5PF.txt");
  vecL2Relative.push_back(*L2Relative);
  vecL3Absolute.push_back(*L3Absolute);
  FactorizedJetCorrector *jecL2Relative   = new FactorizedJetCorrector(vecL2Relative);
  FactorizedJetCorrector *jecL3Absolute   = new FactorizedJetCorrector(vecL3Absolute);

  // Jet energy resolution
  string ak5Tag = "/storage1/lhuertas/CMSTOTEM/mc/Workspace/JetResolutionInputAK5PF.txt";
  JetCorrectorParameters *AK5PFPar    = new JetCorrectorParameters(ak5Tag);
  SimpleJetResolution *ak5PFResolution =  new SimpleJetResolution(*AK5PFPar);
  //cout<<"AK5PFPar isValid:"<<AK5PFPar->isValid()<<endl;
  //cout<<"AK5PFPar size:"<<AK5PFPar->size()<<endl;



  double nevents_total =0;
  rp_aperture_config();

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
    tree->SetBranchAddress("ak5PFJets",&pfJet_coll);
    tree->SetBranchAddress("particleFlow",&pFlow_coll);
    tree->SetBranchAddress("genKin",&genKin);
    tree->SetBranchAddress("genPart",&genPart);
    tree->SetBranchAddress("ak5GenJets",&genJet_coll);
  
    weight= -1.; 
    xi_rec_cms_right = -999.; xi_rec_cms_left = 999.; xi_rec_totem_right = -999.; xi_rec_totem_left = 999.; t_rec_totem_right = -999.; t_rec_totem_left=999.;
    rp_right = false; rp_left = false;

    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){
    
    //printing the % of events done every 10k evts
    if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i_evt);

      bool passedHLT = false;
      bool passedvtx = false;
      bool jet1_selected = false;
      bool jet2_selected = false;
      bool pz_proton_max = false;
      bool PF_eta_max = false;
      bool PF_eta_min = false;
      bool xi_negat_gen = false;
      bool xi_posit_gen = false;

      
      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
      double event_weight = genKin->genWeight; 
      
      weight = event_weight;
      ++nevents_total;
      //-------------------------------------------------------------------------------------------------
 	    
     // Vertices
      for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
      //  if (it_vtx!=vertex_coll->begin()) continue;
	if( it_vtx->ndof>4 ) passedvtx = true;   
      }
      //if(!passedvtx) continue;
      MyVertex& primaryVertex = vertex_coll->at(0);
      double prim_vtx_r = sqrt( primaryVertex.x*primaryVertex.x + primaryVertex.y*primaryVertex.y );
      bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity && primaryVertex.ndof > 4);// && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);
      


      //Jet generated level         
      double leadingJet_pt_gen = -999;
      double secondJet_pt_gen = -999;
      double thirdJet_pt_gen = -999;
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
      double Jet3_energy_gen, Jet3_pz_gen;
      double jet_pt_gen;
      double jet_eta_gen;
      double jet_ene_gen;
      double jet_px_gen;
      double jet_py_gen;
      double jet_pz_gen;
      double jet_phi_gen;

         
      for(vector<MyGenJet>::iterator it_genjet = genJet_coll->begin(); it_genjet != genJet_coll->end(); ++it_genjet){
         jet_pt_gen = it_genjet->Pt();
         jet_eta_gen = it_genjet->Eta();
         jet_ene_gen = it_genjet->E();
         jet_px_gen = it_genjet->Px();
         jet_py_gen = it_genjet->Py();
         jet_pz_gen = it_genjet->Pz();
         jet_phi_gen = it_genjet->Phi();
         //cout<<i_evt<<" "<<pfJet_coll->size();
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
         if (jet_pt_gen>thirdJet_pt_gen && jet_pt_gen<secondJet_pt_gen){
             thirdJet_pt_gen = jet_pt_gen;
             Jet3_energy_gen = jet_ene_gen;
             Jet3_pz_gen = jet_pz_gen;
         }

      }
      //if(leadingJet_pt_gen < 30.) continue;
      //if(secondJet_pt_gen < 30.) continue;
      bool jet1_selected_gen = leadingJet_pt_gen > 30 && fabs(Jet1_eta_gen)<4.4;
      bool jet2_selected_gen = secondJet_pt_gen > 30 && fabs(Jet2_eta_gen)<4.4;
      bool selectJet_gen = jet1_selected_gen && jet2_selected_gen;
      bool selectJet20_gen = leadingJet_pt_gen > 20 && fabs(Jet1_eta_gen)<4.4 && secondJet_pt_gen > 20 && fabs(Jet2_eta_gen)<4.4;

      double mass_jets_gen= sqrt(pow(Jet1_energy_gen+Jet2_energy_gen,2)-pow(Jet1_px_gen+Jet2_px_gen,2)-pow(Jet1_py_gen+Jet2_py_gen,2)-pow(Jet1_pz_gen+Jet2_pz_gen,2));
      double x_plus_gen = ((Jet1_energy_gen+Jet1_pz_gen)+(Jet2_energy_gen+Jet2_pz_gen))/8000;
      double x_minus_gen = ((Jet1_energy_gen-Jet1_pz_gen)+(Jet2_energy_gen-Jet2_pz_gen))/8000;
      double x_gen_selec = (x_minus_gen<x_plus_gen) ? x_minus_gen : x_plus_gen;
      double x_plus_gen_3j = ((Jet1_energy_gen+Jet1_pz_gen)+(Jet2_energy_gen+Jet2_pz_gen)+(Jet3_energy_gen+Jet3_pz_gen))/8000;
      double x_minus_gen_3j = ((Jet1_energy_gen-Jet1_pz_gen)+(Jet2_energy_gen-Jet2_pz_gen)+(Jet3_energy_gen-Jet3_pz_gen))/8000;
 

      //Jets with pt>30Gev and !eta!<2
      Double_t Jet1_pt; 
      Double_t Jet2_pt; 
      Double_t Jet1_eta; 
      Double_t Jet2_eta; 
      Double_t Jet1_phi, Jet1_px, Jet1_py, Jet1_pz, Jet1_energy; 
      Double_t Jet2_phi, Jet2_px, Jet2_py, Jet2_pz, Jet2_energy;
      Double_t Jet3_pt, Jet3_pz, Jet3_energy, Jet3_eta;
      
      vector<MyBaseJet> JetVectorCorrected;
      JetVectorCorrected.resize( pfJet_coll->size() );
      size_t idx_jet = 0;
      for(vector<MyPFJet>::iterator it_jet = pfJet_coll->begin(); it_jet != pfJet_coll->end() ; ++it_jet,++idx_jet){
         map<string,MyBaseJet>::iterator it_map = it_jet->mapjet.begin();
         for(; it_map != it_jet->mapjet.end(); ++it_map)
            if(verbose) cout << it_map->first << endl;

         if( !loosePFJetID(*it_jet,jetCorrName) ) continue;

         MyBaseJet const& basejet = it_jet->mapjet[jetCorrName];

         TLorentzVector oldJet;
         oldJet.SetPxPyPzE(basejet.Px(), basejet.Py(), basejet.Pz(), basejet.E());
         TLorentzVector UnCorrectedJet = oldJet*(1/basejet.jec);
         //cout<<Jet1_px<<"  "<<UnCorrectedJet1.Px()<<endl;   
         // ---- Evaluating the L2Relative correction factor ---- //
         jecL2Relative->setJetPt(UnCorrectedJet.Pt());
         jecL2Relative->setJetEta(UnCorrectedJet.Eta());
         double corFactorL2Relative = jecL2Relative->getCorrection();
               //cout<<"L2Relative Cor Factor"<<corFactorL2Relative<<endl;
         TLorentzVector JetL2Relative = UnCorrectedJet*corFactorL2Relative;

         // ---- Evaluating the L3Absolute correction factor ---- //
         jecL3Absolute->setJetPt(JetL2Relative.Pt());
         jecL3Absolute->setJetEta(JetL2Relative.Eta());
         double corFactorL3Absolute = jecL3Absolute->getCorrection();
         TLorentzVector JetL2RelativeL3Absolute = JetL2Relative*corFactorL3Absolute;

         double CoorFactor = JetL2RelativeL3Absolute.Pt()/UnCorrectedJet.Pt();

         JetVectorCorrected[idx_jet].SetPxPyPzE( JetL2RelativeL3Absolute.Px(),
                                                 JetL2RelativeL3Absolute.Py(),
                                                 JetL2RelativeL3Absolute.Pz(),
                                                 JetL2RelativeL3Absolute.E() );
         JetVectorCorrected[idx_jet].jec = CoorFactor;


      }

      std::stable_sort(JetVectorCorrected.begin(),JetVectorCorrected.end(),sortByPt);


      //Jet energy resolution

      std::vector<float> fx, fY;
      double jetResolution;
      double jetScaleFactor;
      double smearFactor;
      //vector<MyGenJet> genJetMatched;
      MyGenJet const* genJetMatched;
      TRandom rdm_number(12345);
      vector<MyBaseJet> smearedJets;
      smearedJets.resize( pfJet_coll->size() );
      idx_jet = 0;

      for(vector<MyBaseJet>::iterator it_jet = JetVectorCorrected.begin(); it_jet != JetVectorCorrected.end() ; ++it_jet,++idx_jet){

          if (it_jet->Pt()==0 && it_jet->Pz()==0)continue; // Warning in <TVector3::PseudoRapidity>: transvers momentum = 0! return +/- 10e10 
          //cout << it_jet->Pt() << endl;
          double etaJet = it_jet->Eta();
          //if (it_jet->Pt() == 0 && it_jet->Pz() == 0) etaJet = 0;
          //else etaJet = it_jet->Eta();
	  fx.push_back(etaJet);  // Jet Eta
          fY.push_back(it_jet->Pt()); // Jet PT
          fY.push_back(0); // Number of truth pileup
          jetResolution = ak5PFResolution->resolution(fx,fY);
          jetScaleFactor = getScaleFactor(etaJet,variation);
          genJetMatched = getMatch(*it_jet, *genJet_coll, jetResolution);      

          //for(vector<MyGenJet>::iterator it_genjet = genJetMatched.begin(); it_genjet != genJetMatched.end() ; ++it_genjet){
	  if (genJetMatched){
             double dPt = it_jet->Pt() - genJetMatched->Pt();
             //cout << dPt << endl;
             smearFactor = 1 + (jetScaleFactor - 1.)*dPt/it_jet->Pt();

	  }
          else{
             double sigma = jetResolution * std::sqrt(jetScaleFactor * jetScaleFactor - 1);
             smearFactor = 1. + rdm_number.Gaus(0, sigma); 
          }

          fx.clear();	
          fY.clear();
	
          TLorentzVector JetCorrected;
          JetCorrected.SetPxPyPzE( it_jet->Px(), it_jet->Py(), it_jet->Pz(), it_jet->E());
          TLorentzVector JetSmeared = JetCorrected*smearFactor;
          smearedJets[idx_jet].SetPxPyPzE( JetSmeared.Px(), JetSmeared.Py(), JetSmeared.Pz(), JetSmeared.E());
     }

      std::stable_sort(smearedJets.begin(),smearedJets.end(),sortByPt);


/*      if( pfJet_coll->size() > 0 ){
	 //MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet["ak5PFJets"];
	 MyBaseJet const& leadingJet = ( smearedJets.at(0) );
	 Jet1_pt = leadingJet.Pt(); 
	 Jet1_phi = leadingJet.Phi(); 
	 Jet1_px = leadingJet.Px(); 
	 Jet1_py = leadingJet.Py(); 
	 Jet1_pz = leadingJet.Pz(); 
	 Jet1_energy = leadingJet.E(); 
	 Jet1_eta = (Jet1_pz == 0 && Jet1_pt == 0) ? 0 : leadingJet.Eta(); 
	 
 	 if(leadingJet.Pt() > 30. && fabs(leadingJet.Eta())<4.4 ){ 
           jet1_selected = true;}
      }
      
      if( pfJet_coll->size() > 1 ){
	 //MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet["ak5PFJets"];
	 MyBaseJet const& secondJet = ( smearedJets.at(1) );
         Jet2_pt = secondJet.Pt(); 
	 Jet2_phi = secondJet.Phi(); 
         Jet2_px = secondJet.Px();
         Jet2_py = secondJet.Py();
         Jet2_pz = secondJet.Pz();
         Jet2_energy = secondJet.E();
	 Jet2_eta = (Jet2_pz == 0 && Jet2_pt == 0) ? 0 : secondJet.Eta(); 
	 
 	 if(secondJet.Pt() > 30. && fabs(secondJet.Eta())<4.4 ){  
            jet2_selected = true;}
      }
      if( pfJet_coll->size() > 2 ){
	 //MyBaseJet const& thirdJet = ( pfJet_coll->at(2) ).mapjet["ak5PFJets"];
	 MyBaseJet const& thirdJet = ( smearedJets.at(2) );
         Jet3_pt = thirdJet.Pt(); 
         Jet3_pz = thirdJet.Pz();
         Jet3_energy = thirdJet.E();
	 Jet3_eta = (Jet3_pz == 0 && Jet3_pt == 0) ? 0 : thirdJet.Eta(); 
	 
      }
      bool selectJet_rec = jet1_selected && jet2_selected;
      //if( !selectJet) continue;
      double mass_jets= sqrt(pow(Jet1_energy+Jet2_energy,2)-pow(Jet1_px+Jet2_px,2)-pow(Jet1_py+Jet2_py,2)-pow(Jet1_pz+Jet2_pz,2));
      double x_plus = ((Jet1_energy+Jet1_pz)+(Jet2_energy+Jet2_pz))/8000;
      double x_minus = ((Jet1_energy-Jet1_pz)+(Jet2_energy-Jet2_pz))/8000;
      double x_selec = (x_minus<x_plus) ?  x_minus : x_plus;
      double x_plus_3j = ((Jet1_energy+Jet1_pz)+(Jet2_energy+Jet2_pz)+(Jet3_energy+Jet3_pz))/8000;
      double x_minus_3j = ((Jet1_energy-Jet1_pz)+(Jet2_energy-Jet2_pz)+(Jet3_energy-Jet3_pz))/8000;




      // Particle-flow
      double soma1 = 0;
      double soma2 = 0;
      double eta_max=-999.;
      double eta_min=999.;
      double cm = 8000;
      bool pf_thresholds = false;     
 
      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
	 int partType = it_pfcand->particleId;
	 double eta = it_pfcand->Eta();
	 double energy = it_pfcand->Energy();
	 double pz = it_pfcand->Pz();

	 // HF eta rings 29, 30, 40, 41
//         if( ( (fabs(eta) <= 2.866) && (fabs(eta) > 3.152) ) || (fabs(eta) <= 4.730) ){

         // Apply thresholds
         if( !pflowThreshold(*it_pfcand,thresholdsPFlow) ) continue;
         if( pflowThreshold(*it_pfcand,thresholdsPFlow) ) pf_thresholds = true;

	    soma1 += (energy + pz);
	    soma2 += (energy - pz);

            if (eta > eta_max) { eta_max = eta; PF_eta_max = true;} 
	    if (eta < eta_min) { eta_min = eta; PF_eta_min = true;}

  
         
       }


       double xi_plus_Reco = soma1/cm;
       double xi_minus_Reco = soma2/cm;
       double delta_eta_maxmin = eta_max - eta_min;  
       
       //double resolution_before = (xi_plus_gen-xi_plus_Reco)/xi_plus_gen;
       //double xi_reconst = xi_plus_Reco/0.8;
       //double resolution_after = (xi_plus_gen-xi_reconst)/xi_plus_gen;



      //GenPart
      double genEPlusPz = 0;
      double genEMinusPz = 0;
      double proton_pi = 4000;
      double proton_pz_plus=-999;
      double proton_px_plus = -999.;
      double proton_py_plus = -999.;
      double proton_energy_plus = 0.;
      double proton_pz_minus=999;
      double proton_px_minus = 999.;
      double proton_pt_minus = 999.;
      double proton_py_minus = 999.;
      double proton_energy_minus = 0.;
      double px_gen, pt_gen;
      double py_gen;
      double pz_gen;
      double energy_gen;
      double proton_pf;
      double mass_gen; 
      double mx_gen = 0; 
      
      for(vector<MyGenPart>::iterator it_genpart = genPart->begin(); it_genpart != genPart->end(); ++it_genpart){
 
	 double eta_gen = it_genpart->Eta();
         int status = it_genpart->status;
         int id = it_genpart->pdgId;

	 if (status == 1) {
            energy_gen = it_genpart->Energy();
            px_gen = it_genpart->Px();
            py_gen = it_genpart->Py();
            pz_gen = it_genpart->Pz();
            pt_gen = it_genpart->Pt();
            mass_gen = it_genpart->M();
	    proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);  
	    if (id != 2212) {
   	       genEPlusPz += (energy_gen + pz_gen);
	       genEMinusPz += (energy_gen - pz_gen);
               mx_gen += mass_gen;
            }
	    if (id == 2212) {
             double pz_cut = 0.7*proton_pi;
             if (fabs(pz_gen) > pz_cut){

	        if (pz_gen > proton_pz_plus) {
                    proton_pz_plus = pz_gen; proton_energy_plus = energy_gen;
                    proton_px_plus = px_gen; proton_py_plus = py_gen;       
                }
                if (pz_gen < proton_pz_minus) {
                   proton_pz_minus = pz_gen; proton_energy_minus = energy_gen;
                   proton_px_minus = px_gen; proton_py_minus = py_gen;
                   proton_pt_minus = pt_gen;
                }
             }
            }
	 }
      }
	 

      double xi_minus_gen = genEMinusPz/cm;
      double xi_plus_gen = genEPlusPz/cm;
      double correction = xi_plus_Reco/xi_plus_gen;
      xi_true = mx_gen*mx_gen/(cm);

      double xi_proton_plus = -1.;
      double xi_proton_minus = -1.;
      double t_proton_plus = 0.;
      double t_proton_minus = 0.;
      double thx_proton_plus = 0.;
      double thy_proton_plus = 0.;
      double thx_proton_minus = 0.;
      double thy_proton_minus = 0.;


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

      std::map<int,std::vector<double> > proton_plus_pars;
      std::map<int,std::vector<double> > proton_minus_pars;

      bool fiducial_cut_rp_024= false; 
      bool fiducial_cut_rp_025 = false;
      if(proton_pz_plus > 0.){
         xi_proton_plus =  ( 1 - (proton_pz_plus/proton_pi) );
         //t_proton_plus = -2*( (proton_pi*proton_energy_plus) - (proton_pi*proton_pz_plus) );
         TLorentzVector vec_pi(0.,0.,proton_pi,proton_pi);
         TLorentzVector vec_pf(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus);
         TLorentzVector vec_t = (vec_pf - vec_pi);
         t_proton_plus = vec_t.Mag2();
         thx_proton_plus = atan(-proton_px_plus/proton_pi);
         thy_proton_plus = atan(proton_py_plus/proton_pi);

         //FIXME
         double out_x, out_thx, out_y, out_thy, out_xi;
         proton_plus_rp_accept_020 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 20, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[20] = std::vector<double>(5,0.);
         proton_plus_pars[20][0] = out_x; proton_plus_pars[20][1] = out_y;
         proton_plus_pars[20][2] = out_thx; proton_plus_pars[20][3] = out_thy;
         proton_plus_pars[20][4] = out_xi;

         proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[24] = std::vector<double>(5,0.);
         proton_plus_pars[24][0] = out_x; proton_plus_pars[24][1] = out_y;
         proton_plus_pars[24][2] = out_thx; proton_plus_pars[24][3] = out_thy;
         proton_plus_pars[24][4] = out_xi;
         fiducial_cut_rp_024 = proton_plus_pars[24][0]>0 && proton_plus_pars[24][0]<0.006 && proton_plus_pars[24][1]>0.0084 && proton_plus_pars[24][1]<0.029;

         proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[25] = std::vector<double>(5,0.);
         proton_plus_pars[25][0] = out_x; proton_plus_pars[25][1] = out_y;
         proton_plus_pars[25][2] = out_thx; proton_plus_pars[25][3] = out_thy;
         proton_plus_pars[25][4] = out_xi;
         fiducial_cut_rp_025 = proton_plus_pars[25][0]>0 && proton_plus_pars[25][0]<0.006 && proton_plus_pars[25][1]<-0.0084 && proton_plus_pars[25][1]>-0.029;

         rp_xpos_024 = proton_plus_pars[24][0];
         rp_ypos_024 = proton_plus_pars[24][1];
         rp_xpos_025 = proton_plus_pars[25][0];
         rp_ypos_025 = proton_plus_pars[25][1];

         proton_plus_rp_accept_021 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 21);
         proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22);
         proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23);
         proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24);
         proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25);
         proton_plus_rp_accept_120 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 120);
      }

      bool fiducial_cut_rp_124= false; 
      bool fiducial_cut_rp_125 = false;
      if(proton_pz_minus < 0.){
         xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
         // t_proton_minus = -2*( (proton_pi*proton_energy_minus) + (proton_pi*proton_pz_minus) ); 
         TLorentzVector vec_pi(0.,0.,-proton_pi,proton_pi);
         TLorentzVector vec_pf(proton_px_minus,proton_py_minus,proton_pz_minus,proton_energy_minus);
         TLorentzVector vec_t = (vec_pf - vec_pi);
         t_proton_minus = vec_t.Mag2();
 
         thx_proton_minus = atan(-proton_px_minus/proton_pi);
         thy_proton_minus = atan(proton_py_minus/proton_pi);

         double out_x, out_thx, out_y, out_thy, out_xi;
         proton_minus_rp_accept_120 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 120, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[120] = std::vector<double>(5,0.);
         proton_minus_pars[120][0] = out_x; proton_minus_pars[120][1] = out_y;
         proton_minus_pars[120][2] = out_thx; proton_minus_pars[120][3] = out_thy;
         proton_minus_pars[120][4] = out_xi;

         proton_minus_rp_accept_124 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 124, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[124] = std::vector<double>(5,0.);
         proton_minus_pars[124][0] = out_x; proton_minus_pars[124][1] = out_y;
         proton_minus_pars[124][2] = out_thx; proton_minus_pars[124][3] = out_thy;
         proton_minus_pars[124][4] = out_xi;
         fiducial_cut_rp_124 = proton_minus_pars[124][0]>0 && proton_minus_pars[124][0]<0.006 && proton_minus_pars[124][1]>0.0084 && proton_minus_pars[124][1]<0.027;

         proton_minus_rp_accept_125 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 125, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[125] = std::vector<double>(5,0.);
         proton_minus_pars[125][0] = out_x; proton_minus_pars[125][1] = out_y;
         proton_minus_pars[125][2] = out_thx; proton_minus_pars[125][3] = out_thy;
         proton_minus_pars[125][4] = out_xi;
         fiducial_cut_rp_125 = proton_minus_pars[125][0]>0 && proton_minus_pars[125][0]<0.006 && proton_minus_pars[125][1]<-0.0084 && proton_minus_pars[125][1]>-0.027;
         rp_xpos_124 = proton_minus_pars[124][0]; 
         rp_ypos_124 = proton_minus_pars[124][1]; 
         rp_xpos_125 = proton_minus_pars[125][0]; 
         rp_ypos_125 = proton_minus_pars[125][1]; 

         proton_minus_rp_accept_121 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 121);
         proton_minus_rp_accept_122 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 122);
         proton_minus_rp_accept_123 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 123);
         proton_minus_rp_accept_124 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 124);
         proton_minus_rp_accept_125 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 125);
         proton_minus_rp_accept_020 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 20);
      }



     //totem proton reconstructed
      float sigma_xi56 = 0.00720615 - 0.0418783*xi_proton_minus + 0.0999515*xi_proton_minus*xi_proton_minus; // sigma56 vs xi from Hubert
      float xi_proton_minus_rec = xi_proton_minus + gRandom->Gaus(0,sigma_xi56);
      double sigma_t56 = 0.233365*t_proton_minus - 0.0975751*t_proton_minus*t_proton_minus;  // sigma_t56 vs t from Hubert
      double t_proton_minus_rec = t_proton_minus + gRandom->Gaus(0,sigma_t56);
      float sigma_xi45=0.00714986 - 0.0408903*xi_proton_plus + 0.0965813*xi_proton_plus*xi_proton_plus;   // sigma45 vs xi from Hubert
      float xi_proton_plus_rec = xi_proton_plus + gRandom->Gaus(0,sigma_xi45);
      double sigma_t45=0.233365*t_proton_plus - 0.0975751*t_proton_plus*t_proton_plus;     // sigma_t45 vs t from Hubert
      double t_proton_plus_rec = t_proton_plus + gRandom->Gaus(0,sigma_t45);
      double proton_minus_beta = (Jet3_pt>20) ? x_minus_3j/xi_proton_minus_rec : x_minus/xi_proton_minus_rec;
      double proton_plus_beta = (Jet3_pt>20) ? x_plus_3j/xi_proton_plus_rec : x_plus/xi_proton_plus_rec;



      //rp_accept
      bool proton_minus_rp_accept_mc = ( proton_minus_rp_accept_120 && proton_minus_rp_accept_124 && fiducial_cut_rp_124) || ( proton_minus_rp_accept_121 && proton_minus_rp_accept_125 && fiducial_cut_rp_125 );// || ( proton_minus_rp_accept_122 && proton_minus_rp_accept_123 );
      bool proton_plus_rp_accept_mc = ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 && fiducial_cut_rp_024 ) || ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 && fiducial_cut_rp_025 );// || ( proton_plus_rp_accept_022 && proton_plus_rp_accept_023 );


      if (!(select_Vertex && pf_thresholds)) continue;
      
      rp_right = proton_minus_rp_accept_mc;
      rp_right_top_near =  proton_minus_rp_accept_120;
      rp_right_top_far =  proton_minus_rp_accept_124;
      rp_right_bottom_near =  proton_minus_rp_accept_121;
      rp_right_bottom_far =  proton_minus_rp_accept_124;
      rp_left = proton_plus_rp_accept_mc;
      rp_left_top_near =  proton_plus_rp_accept_020;
      rp_left_top_far =  proton_plus_rp_accept_024;
      rp_left_bottom_near =  proton_plus_rp_accept_021;
      rp_left_bottom_far =  proton_plus_rp_accept_025;
      xi_rec_cms_right = xi_minus_Reco;
      xi_rec_cms_left = xi_plus_Reco;
      xi_rec_totem_right =  xi_proton_minus_rec;    
      xi_rec_totem_left =  xi_proton_plus_rec;    
      t_rec_totem_right = fabs(t_proton_minus_rec);    
      t_rec_totem_left = fabs(t_proton_plus_rec);    
      xi_gen_cms_right = xi_minus_gen;
      xi_gen_cms_left = xi_plus_gen;
      xi_gen_totem_right =  xi_proton_minus;    
      xi_gen_totem_left =  xi_proton_plus;    
      t_gen_totem_right = fabs(t_proton_minus);    
      t_gen_totem_left = fabs(t_proton_plus);    
      beta_rec_right =  proton_minus_beta;
      beta_rec_left = proton_plus_beta;
      beta_gen_right = (thirdJet_pt_gen>20) ? x_minus_gen_3j/xi_proton_minus : x_minus_gen/xi_proton_minus;
      beta_gen_left = (thirdJet_pt_gen>20) ? x_plus_gen_3j/xi_proton_plus : x_plus_gen/xi_proton_plus;
      jet1_rec_pt = Jet1_pt;
      jet1_rec_eta = Jet1_eta;
      jet1_rec_phi = Jet1_phi;
      jet2_rec_pt = Jet2_pt;
      jet2_rec_eta = Jet2_eta;
      jet2_rec_phi = Jet2_phi;
      x_rec_right = (Jet3_pt>20) ? x_minus_3j : x_minus;
      x_rec_left = (Jet3_pt>20) ? x_plus_3j : x_plus;
      jet1_gen_pt = leadingJet_pt_gen;
      jet1_gen_eta = Jet1_eta_gen;
      jet1_gen_phi = Jet1_phi_gen;
      jet2_gen_pt = secondJet_pt_gen;
      jet2_gen_eta = Jet2_eta_gen;
      jet2_gen_phi = Jet2_phi_gen;
      x_gen_right = (thirdJet_pt_gen > 20) ? x_minus_gen_3j : x_minus_gen;
      x_gen_left = (thirdJet_pt_gen>20) ? x_plus_gen_3j : x_plus_gen;
*/
      small_tree->Fill();    
    }//end loop for events

   file->Close();
 
  }//end of loop over files
  

  //output file
  //TFile* output = new TFile(outputFileName,"RECREATE");
  output->cd();
 
  float cross_section; //2.213e10; //pb cross section
  if (pt_15_30) cross_section = 7.9e08;//7.59e08;
  if (pt_30_50) cross_section = 5.32e07;
  if (pt_50_80) cross_section = 6.55e06;//6.48e06;
  if (pt_80_120) cross_section = 8.36e05;//8.54e05;
  if (pt_120_170) cross_section = 1.27e05;
  if (pt_170_300) cross_section = 2.79e04;//2.83e04;
  if (pt_300_470) cross_section = 1.40e03;
  if (mc == "pythia6") cross_section = 2.998e+10; //2.213e10; //pb cross section
  if (mc == "pythia8_nondiff") cross_section = 1.17e+09; //2.213e10; //pb cross section
  if (mc == "pythia8_nondiff_CUETP8S1") cross_section = 1.279e+09; //2.213e10; //pb cross section
  if (mc == "pythia8_nondiff_CUETP8M1") cross_section = 1.262e+09; //2.213e10; //pb cross section

  float luminity_HLT_L1Jet1_198902 = 0.015879;//pb ---- luminity for LP_Jets1_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet1_198903 = 0.008698;//pb ---- luminity for LP_Jets1_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198903
  float luminity_HLT_L1Jet2_198902 = 0.015879;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet2_198903 = 0.008698;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198903
  float luminity = luminity_HLT_L1Jet1_198902 + luminity_HLT_L1Jet1_198903 + luminity_HLT_L1Jet2_198902 + luminity_HLT_L1Jet2_198903;
 
  float n_events = luminity*cross_section;
  
  float f1 = (float) nevents_total;
  //float f1 = (float) nweight_total;
  Double_t scale = n_events/f1;
  cout << "cross section:  "<< cross_section << "   Total events:   " << nevt_tot << "   scale  "<< scale << endl;   
  
//   float f3 = (float) weight_total_PF_selected ; 
//   Double_t scale_PF = 1.0/f3;
  
  //histosTH2F["xi_plus_reco_gen"]->SetOption("colz");
  //histosTH2F["xi_minus_reco_gen"]->SetOption("colz");
  //histosTH2F["logxi_plus_reco_gen"]->SetOption("colz");
  //histosTH2F["logxi_minus_reco_gen"]->SetOption("colz");

  //rdm_evt_zb->Write();

  small_tree->SetDirectory(0);

  small_tree->Write();

  output->Close();
}
