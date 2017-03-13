//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>
#include <TMath.h>
#include <TVectorT.h>
#include <THStack.h>
#include <TLatex.h>
#include <TRandom.h>
#include <TKey.h>
//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;


class MCInclusiveAnalysis
{
public:
	MCInclusiveAnalysis();
	~MCInclusiveAnalysis();
	void getMCHistos( string const& mc = "pythia6");
  void addHerwigHistos ();

private:
	TFile * outfile;
};

MCInclusiveAnalysis::MCInclusiveAnalysis (void){
}

MCInclusiveAnalysis::~MCInclusiveAnalysis (void){
}


void MCInclusiveAnalysis::getMCHistos ( string const& mc){

    TFile* pythia = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pythia6_QCD_pt_15_3000_ntuple.root","READ");
    TFile* pythia8_QCD_Tune4C = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pythia8_nondiff_QCD_pt_15_3000_ntuple.root","READ");
    TFile* pythia8_QCD_CUETP8M1 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pythia8_nondiff_CUETP8M1_QCD_pt_15to3000_ntuple.root","READ");
    TFile* pythia8_QCD_CUETP8S1 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pythia8_nondiff_CUETP8S1_QCD_15to3000_ntuple.root","READ");
    TFile* herwig_15_30 = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/herwig_QCD_pt_15_30_ntuple.root","READ");
    TFile* herwig_30_50 = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/herwig_QCD_pt_30_50_ntuple.root","READ");
    TFile* herwig_50_80 = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/herwig_QCD_pt_50_80_ntuple.root","READ");
    TFile* herwig_80_120 = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/herwig_QCD_pt_80_120_ntuple.root","READ");
    TFile* herwig_120_170 = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/herwig_QCD_pt_120_170_ntuple.root","READ");
    TFile* herwig_170_300 = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/herwig_QCD_pt_170_300_ntuple.root","READ");
    TFile* mc_file;

    float bin_sasha[9] ={0.0003, 0.002, 0.0045, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1};
    float xi_bins[12] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2};

    map<string,TH1F*> TH1F_histos;
    map<string,TH2F*> TH2F_histos;

    TH1F_histos["pt_jet1"] = new TH1F("pt_jet1","", 15, 0, 200);
    TH1F_histos["pt_jet2"] = new TH1F("pt_jet2","", 15, 0, 200);
    TH1F_histos["eta_jet1"] = new TH1F("eta_jet1","", 20, -5.2, 5.2);
    TH1F_histos["eta_jet2"] = new TH1F("eta_jet2","", 20, -5.2, 5.2);
    TH1F_histos["log_x_rec_minus"] = new TH1F("log_x_rec_minus","",15, -4, 0);
    TH1F_histos["log_x_gen_minus"] = new TH1F("log_x_gen_minus","",15, -4, 0);
    TH1F_histos["log_x_rec_plus"] = new TH1F("log_x_rec_plus","",15, -4, 0);
    TH1F_histos["log_x_gen_plus"] = new TH1F("log_x_gen_plus","",15, -4, 0);
    TH1F_histos["xi_cms_rec_right"] = new TH1F("xi_cms_rec_right","",11,xi_bins);
    TH1F_histos["xi_cms_gen_right"] = new TH1F("xi_cms_gen_right","",11,xi_bins);
    TH1F_histos["xi_cms_rec_right_sasha"] = new TH1F("xi_cms_rec_right_sasha","",8, bin_sasha);
    TH1F_histos["xi_cms_gen_right_sasha"] = new TH1F("xi_cms_gen_right_sasha","",8, bin_sasha);
    TH1F_histos["reco_logx_minus"] = new TH1F("reco_logx_minus","",15, -4, 0);
    TH1F_histos["truth_logx_minus"] = new TH1F("truth_logx_minus","",15, -4, 0);
    TH2F_histos["response_logx_minus"] = new TH2F("response_logx_minus","",15, -4, 0,15, -4, 0);
    TH1F_histos["reco_logx_plus"] = new TH1F("reco_logx_plus","",15, -4, 0);
    TH1F_histos["truth_logx_plus"] = new TH1F("truth_logx_plus","",15, -4, 0);
    TH2F_histos["response_logx_plus"] = new TH2F("response_logx_plus","",15, -4, 0,15, -4, 0);
    
    for(map<string,TH1F*>::const_iterator histos = TH1F_histos.begin(); histos != TH1F_histos.end(); ++histos){
       histos->second->Sumw2();
    }

    double norm; 
    string treeName = "small_tree"; 
    string fileName;
    if (mc == "pythia6") {mc_file = pythia; norm = 1488.52; fileName = "histos_pythia6.root"; }
    if (mc == "pythia8_tune4C") {mc_file = pythia8_QCD_Tune4C; norm = 57.5102; fileName = "histos_pythia8_4c_nondiff.root"; }
    if (mc == "pythia8_CUETP8M1") {mc_file = pythia8_QCD_CUETP8M1; norm =(0.049154*1.027e+09)/2203000; fileName = "histos_pythia8_cuetp8m1_nondiff.root"; }// 28.1581;}
    if (mc == "pythia8_CUETP8S1") {mc_file = pythia8_QCD_CUETP8S1; norm = (0.049154*1.24e+09)/2390000; fileName = "histos_pythia8_cuetp8s1_nondiff.root"; }//26.3046;}
    if (mc == "herwig_15_30") {mc_file = herwig_15_30; norm = 0.000236165; fileName = "histos_herwig_pt_15_30.root"; }
    if (mc == "herwig_30_50") {mc_file = herwig_30_50; norm = 0.000246105; fileName = "histos_herwig_pt_30_50.root"; }
    if (mc == "herwig_50_80") {mc_file = herwig_50_80; norm = 0.000243068; fileName = "histos_herwig_pt_50_80.root"; }
    if (mc == "herwig_80_120") {mc_file = herwig_80_120; norm = 0.000251334; fileName = "histos_herwig_pt_80_120.root"; }
    if (mc == "herwig_120_170") {mc_file = herwig_120_170; norm = 0.000246304; fileName = "histos_herwig_pt_120_170.root"; }
    if (mc == "herwig_170_300") {mc_file = herwig_170_300; norm = 0.000249124; fileName = "histos_herwig_pt_170_300.root"; }

    outfile = new TFile(fileName.c_str(), "RECREATE");
    cout << "output: " << fileName <<endl;

    TTree* tree= (TTree*) mc_file->Get( treeName.c_str() );
    int nev = int(tree->GetEntriesFast());
    cout << mc <<" non_diffractive file has " << nev << " entries  " << endl;
 
    double jet1_rec_pt, jet1_rec_eta, jet1_rec_phi, jet2_rec_pt, jet2_rec_eta, jet2_rec_phi;
    double jet1_gen_pt, jet1_gen_eta, jet1_gen_phi, jet2_gen_pt, jet2_gen_eta, jet2_gen_phi;
    double x_rec_right, x_rec_left, xi_rec_cms_right, xi_gen_cms_right, x_gen_right, x_gen_left;
    double weight;
    tree->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt);
    tree->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt);
    tree->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta);
    tree->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta);
    tree->SetBranchAddress("jet1_rec_phi",&jet1_rec_phi);
    tree->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt);
    tree->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt);
    tree->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta);
    tree->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta);
    tree->SetBranchAddress("jet2_rec_phi",&jet2_rec_phi);
    tree->SetBranchAddress("x_rec_right",&x_rec_right);
    tree->SetBranchAddress("x_gen_right",&x_gen_right);
    tree->SetBranchAddress("x_rec_left",&x_rec_left);
    tree->SetBranchAddress("x_gen_left",&x_gen_left);
    tree->SetBranchAddress("xi_rec_cms_right",&xi_rec_cms_right);
    tree->SetBranchAddress("xi_gen_cms_right",&xi_gen_cms_right);
    tree->SetBranchAddress("weight",&weight);

    double pt_threshold = 40.;

    for(int i_evt = 0; i_evt < nev; ++i_evt){
        tree->GetEntry(i_evt);

        bool jet_rec_sel = jet1_rec_pt>pt_threshold && jet2_rec_pt>pt_threshold && fabs(jet1_rec_eta)<4.4 && fabs(jet2_rec_eta)<4.4;
        bool jet_gen_sel = jet1_gen_pt>pt_threshold && jet2_gen_pt>pt_threshold && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4;
           
        if (jet_rec_sel){
            TH1F_histos["pt_jet1"]->Fill(jet1_rec_pt, weight);
            TH1F_histos["pt_jet2"]->Fill(jet2_rec_pt, weight);
            TH1F_histos["eta_jet1"]->Fill(jet1_rec_eta, weight);
            TH1F_histos["eta_jet2"]->Fill(jet2_rec_eta, weight);
            TH1F_histos["log_x_rec_minus"]->Fill(log10(x_rec_right),weight);
            TH1F_histos["log_x_rec_plus"]->Fill(log10(x_rec_left),weight);
            TH1F_histos["xi_cms_rec_right"]->Fill(xi_rec_cms_right, weight);
            TH1F_histos["xi_cms_rec_right_sasha"]->Fill(xi_rec_cms_right, weight);
            TH1F_histos["reco_logx_minus"]->Fill(log10(x_rec_right),weight);
            TH1F_histos["reco_logx_plus"]->Fill(log10(x_rec_left),weight);
            if (jet_gen_sel){
                TH2F_histos["response_logx_minus"]->Fill(log10(x_rec_right), log10(x_gen_right), weight*norm);
                TH2F_histos["response_logx_plus"]->Fill(log10(x_rec_left), log10(x_gen_left), weight*norm);
            }
        } 
        if (jet_gen_sel){
            TH1F_histos["xi_cms_gen_right"]->Fill(xi_gen_cms_right, weight);
            TH1F_histos["log_x_gen_minus"]->Fill(log10(x_gen_right),weight);
            TH1F_histos["log_x_gen_plus"]->Fill(log10(x_gen_left),weight);
            TH1F_histos["truth_logx_minus"]->Fill(log10(x_gen_right),weight);
            TH1F_histos["truth_logx_plus"]->Fill(log10(x_gen_left),weight);
        }
        if(jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4) TH1F_histos["xi_cms_gen_right_sasha"]->Fill(xi_gen_cms_right, weight);
 
    }

    for(map<string,TH1F*>::const_iterator histos = TH1F_histos.begin(); histos != TH1F_histos.end(); ++histos){
       histos->second->Scale(norm);
       histos->second->Write();
    }
    for(map<string,TH2F*>::const_iterator histos = TH2F_histos.begin(); histos != TH2F_histos.end(); ++histos){
       histos->second->Write();
    }
}
void MCInclusiveAnalysis::addHerwigHistos (void){

  TFile* herwig1 = TFile::Open("histos_herwig_pt_15_30.root","READ");
  TFile* herwig2 = TFile::Open("histos_herwig_pt_30_50.root","READ");
  TFile* herwig3 = TFile::Open("histos_herwig_pt_50_80.root","READ");
  TFile* herwig4 = TFile::Open("histos_herwig_pt_80_120.root","READ");
  TFile* herwig5 = TFile::Open("histos_herwig_pt_120_170.root","READ");
  TFile* herwig6 = TFile::Open("histos_herwig_pt_170_300.root","READ");
  TFile* outfile = new TFile("histos_herwig_nondiff.root", "RECREATE");


  TIter next_herwig1 (herwig1->GetListOfKeys());
  TKey* key_herwig1;
  while ((key_herwig1 = (TKey*)next_herwig1())){
    TClass *cl_herwig1 = gROOT->GetClass(key_herwig1->GetClassName());
    if (!cl_herwig1->InheritsFrom("TH1")) continue;
    std::string histname_herwig1 = key_herwig1->GetName();
    TH1F* histo_herwig1 = (TH1F*)key_herwig1->ReadObj();

    TIter next_herwig2 (herwig2->GetListOfKeys());
    TKey* key_herwig2;
    while ((key_herwig2 = (TKey*)next_herwig2())){
        TClass *cl_herwig2 = gROOT->GetClass(key_herwig2->GetClassName());
        if (!cl_herwig2->InheritsFrom("TH1")) continue;
        std::string histname_herwig2 = key_herwig2->GetName();
        TH1F* histos_herwig2 = (TH1F*)key_herwig2->ReadObj();
        if (histname_herwig1 == histname_herwig2) histo_herwig1->Add(histos_herwig2); 
    }

    TIter next_herwig3 (herwig3->GetListOfKeys());
    TKey* key_herwig3;
    while ((key_herwig3 = (TKey*)next_herwig3())){
        TClass *cl_herwig3 = gROOT->GetClass(key_herwig3->GetClassName());
        if (!cl_herwig3->InheritsFrom("TH1")) continue;
        std::string histname_herwig3 = key_herwig3->GetName();
        TH1F* histos_herwig3 = (TH1F*)key_herwig3->ReadObj();
        if (histname_herwig1 == histname_herwig3) histo_herwig1->Add(histos_herwig3); 
    }
    TIter next_herwig4 (herwig4->GetListOfKeys());
    TKey* key_herwig4;
    while ((key_herwig4 = (TKey*)next_herwig4())){
        TClass *cl_herwig4 = gROOT->GetClass(key_herwig4->GetClassName());
        if (!cl_herwig4->InheritsFrom("TH1")) continue;
        std::string histname_herwig4 = key_herwig4->GetName();
        TH1F* histos_herwig4 = (TH1F*)key_herwig4->ReadObj();
        if (histname_herwig1 == histname_herwig4) histo_herwig1->Add(histos_herwig4); 
    }
    TIter next_herwig5 (herwig5->GetListOfKeys());
    TKey* key_herwig5;
    while ((key_herwig5 = (TKey*)next_herwig5())){
        TClass *cl_herwig5 = gROOT->GetClass(key_herwig5->GetClassName());
        if (!cl_herwig5->InheritsFrom("TH1")) continue;
        std::string histname_herwig5 = key_herwig5->GetName();
        TH1F* histos_herwig5 = (TH1F*)key_herwig5->ReadObj();
        if (histname_herwig1 == histname_herwig5) histo_herwig1->Add(histos_herwig5); 
    }
    TIter next_herwig6 (herwig6->GetListOfKeys());
    TKey* key_herwig6;
    while ((key_herwig6 = (TKey*)next_herwig6())){
        TClass *cl_herwig6 = gROOT->GetClass(key_herwig6->GetClassName());
        if (!cl_herwig6->InheritsFrom("TH1")) continue;
        std::string histname_herwig6 = key_herwig6->GetName();
        TH1F* histos_herwig6 = (TH1F*)key_herwig6->ReadObj();
        if (histname_herwig1 == histname_herwig6) histo_herwig1->Add(histos_herwig6); 
    }

    histo_herwig1->Write();
  }
  outfile->Close();    
}

// int main(void)
// {
//  MCInclusiveAnalysis pythia;
//  // pythia.getMCHistos("herwig_15_30");
//  // pythia.getMCHistos("herwig_30_50");
//  // pythia.getMCHistos("herwig_50_80");
//  // pythia.getMCHistos("herwig_80_120");
//  // pythia.getMCHistos("herwig_120_170");
//  // pythia.getMCHistos("herwig_170_300");
//  pythia.addHerwigHistos();
//  return 0;
// }