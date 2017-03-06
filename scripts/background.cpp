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


class RdmDataEvent {
   public:
      RdmDataEvent() {}
      ~RdmDataEvent() {}

      double xi_totem_right_rdm;
      double xi_totem_left_rdm;
      double t_totem_right_rdm;
      double t_totem_left_rdm;
      double x_right_rdm;
      double x_left_rdm;
      double beta_right_rdm;
      double beta_left_rdm;
      double jet1_pt_right_rdm;
      double jet1_pt_left_rdm;
      double jet2_pt_right_rdm;
      double jet2_pt_left_rdm;
      bool rp_right_rdm;
      bool rp_left_rdm;
      bool valid_vtx_rdm;
      bool valid_proton_right_rdm;
      bool valid_proton_left_rdm;
};      

class BackgroundZB
{
public:
	BackgroundZB();
	~BackgroundZB();
	void getBackgroundHistos(bool full = false, bool vtx_sel = false);
  void BackgNormalisation(TH1F* data, TH1F* zb, double &norm );

};

BackgroundZB::BackgroundZB(void){}

BackgroundZB::~BackgroundZB(void){}

void BackgroundZB::getBackgroundHistos(bool full, bool vtx_sel){

	string treeName = "small_tree";  

    float tbins[9] = {0.03, 0.07, 0.11, 0.21, 0.31, 0.45, 0.58, 0.72, 1.};
    float xi_bins[12] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2};
    float bin[16] = {-0.4, -0.112, -0.096, -0.08, -0.064, -0.048, -0.032, -0.016, 0, 0.048, 0.112, 0.176, 0.24, 0.304, 0.368, 0.4};

    TFile* mc = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pythia6_QCD_pt_15_3000_ntuple.root","READ");
    // TFile* zb = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/rereco/zb_ntuple_novtxsel.root","READ");
    TFile* zb = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/zb_ntuple_noveto.root","READ");

    TFile* outfile = new TFile("histos_zb_background.root", "RECREATE");
    map<string,TH1F*> TH1F_histos;
    TH1F_histos["xi_cms_minus_totem_right_backg"] = new TH1F("xi_cms_minus_totem_right_backg","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_left_backg"] = new TH1F("xi_cms_minus_totem_left_backg","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_right_backg_bin"] = new TH1F("xi_cms_minus_totem_right_backg_bin","", 15, bin);
    TH1F_histos["xi_cms_minus_totem_left_backg_bin"] = new TH1F("xi_cms_minus_totem_left_backg_bin","", 15, bin);
    TH1F_histos["xi_right_backg"] = new TH1F("xi_right_backg","",50,-0.04,0.2);
    TH1F_histos["xi_left_backg"] = new TH1F("xi_left_backg","",50,-0.04,0.2);
    TH1F_histos["xi_right_cut_backg"] = new TH1F("xi_right_cut_backg","",11, xi_bins);
    TH1F_histos["xi_left_cut_backg"] = new TH1F("xi_left_cut_backg","", 11, xi_bins);
    TH1F_histos["t_right_backg"] = new TH1F("t_right_backg","",8, tbins);
    TH1F_histos["t_left_backg"] = new TH1F("t_left_backg","", 8, tbins);
    TH1F_histos["t_right_cut_backg"] = new TH1F("t_right_cut_backg","",8, tbins);
    TH1F_histos["t_left_cut_backg"] = new TH1F("t_left_cut_backg","",8, tbins);
    TH1F_histos["log_x_right_cut_backg"] = new TH1F("log_x_right_cut_backg","",15, -4, 0);
    TH1F_histos["log_x_left_cut_backg"] = new TH1F("log_x_left_cut_backg","",15, -4, 0);
    TH1F_histos["beta_right_cut_backg"] = new TH1F("beta_right_cut_backg","",15,0,1);
    TH1F_histos["beta_left_cut_backg"] = new TH1F("beta_left_cut_backg","",15,0,1);

    TTree* tree_zb;
    tree_zb = (TTree*) zb->Get( treeName.c_str() );
    int nev_zb = int(tree_zb->GetEntriesFast());
    cout <<"The zb file has " << nev_zb << " entries : " << endl;
 
    double xi_proton_right_zb, xi_proton_left_zb, x_right_zb, x_left_zb, beta_right_zb, beta_left_zb, jet1_pt_zb, jet2_pt_zb;
    double t_proton_right_zb, t_proton_left_zb;
    bool rp_right_zb, rp_left_zb, valid_vtx_zb, valid_proton_right_zb, valid_proton_left_zb;
    tree_zb->SetBranchAddress("xi_totem_right",&xi_proton_right_zb);
    tree_zb->SetBranchAddress("xi_totem_left",&xi_proton_left_zb);
    tree_zb->SetBranchAddress("x_right",&x_right_zb);
    tree_zb->SetBranchAddress("x_left",&x_left_zb);
    tree_zb->SetBranchAddress("beta_proton_right",&beta_right_zb);
    tree_zb->SetBranchAddress("beta_proton_left",&beta_left_zb);
    tree_zb->SetBranchAddress("t_totem_right",&t_proton_right_zb);
    tree_zb->SetBranchAddress("t_totem_left",&t_proton_left_zb);
    tree_zb->SetBranchAddress("rp_right",&rp_right_zb);
    tree_zb->SetBranchAddress("rp_left",&rp_left_zb);
    tree_zb->SetBranchAddress("valid_vtx",&valid_vtx_zb);
    tree_zb->SetBranchAddress("valid_proton_right",&valid_proton_right_zb);
    tree_zb->SetBranchAddress("valid_proton_left",&valid_proton_left_zb);
    tree_zb->SetBranchAddress("jet1_pt",&jet1_pt_zb);
    tree_zb->SetBranchAddress("jet2_pt",&jet2_pt_zb);

    std::vector<RdmDataEvent> pu;
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev_zb; ++i_evt){
        tree_zb->GetEntry(i_evt);

        RdmDataEvent rdm;
        if (full==true){
           rdm.valid_vtx_rdm = valid_vtx_zb;
           rdm.valid_proton_right_rdm = valid_proton_right_zb;
           rdm.valid_proton_left_rdm = valid_proton_left_zb;
           rdm.rp_right_rdm = rp_right_zb;
           rdm.rp_left_rdm = rp_left_zb;
           rdm.xi_totem_right_rdm = xi_proton_right_zb;
           rdm.xi_totem_left_rdm = xi_proton_left_zb;
           rdm.t_totem_right_rdm = t_proton_right_zb;
           rdm.t_totem_left_rdm = t_proton_left_zb;
           rdm.x_right_rdm = x_right_zb;
           rdm.x_left_rdm = x_left_zb;
           rdm.beta_right_rdm = beta_right_zb;
           rdm.beta_left_rdm = beta_left_zb;
        }

        else {
          if (vtx_sel && valid_vtx_zb) continue;

          if (/*!valid_vtx_zb &&*/ valid_proton_right_zb && rp_right_zb && !valid_proton_left_zb && !rp_left_zb && t_proton_right_zb>0.03 && t_proton_right_zb<1.){
              rdm.xi_totem_right_rdm = xi_proton_right_zb;
              rdm.t_totem_right_rdm = t_proton_right_zb;
              rdm.valid_vtx_rdm = valid_vtx_zb;
              rdm.valid_proton_right_rdm = valid_proton_right_zb;
              rdm.rp_right_rdm = rp_right_zb;
              rdm.x_right_rdm = x_right_zb;
              rdm.beta_right_rdm = beta_right_zb;
          } 
          if (/*!valid_vtx_zb &&*/ valid_proton_left_zb && rp_left_zb && !valid_proton_right_zb && !rp_right_zb && t_proton_left_zb>0.03 && t_proton_left_zb<1.){
              rdm.xi_totem_left_rdm = xi_proton_left_zb;
              rdm.t_totem_left_rdm = t_proton_left_zb;
              rdm.valid_proton_left_rdm = valid_proton_left_zb;
              rdm.rp_left_rdm = rp_left_zb;
              rdm.x_left_rdm = x_left_zb;
              rdm.beta_left_rdm = beta_left_zb;
           }
        }
           
        pu.push_back(rdm);

    }

 	TRandom* rdm = new TRandom;
  	rdm->SetSeed(12345);

    TTree* tree_non_diff = (TTree*) mc->Get( treeName.c_str() );
    int nev = int(tree_non_diff->GetEntriesFast());
    cout <<"The non diff. file has " << nev << " entries : " << endl;
    int nevent_backg_right = 0;
    int nevent_backg_left = 0;

    double xi_rec_proton_right, xi_rec_proton_left, weight, xi_rec_cms_minus, xi_rec_cms_plus, x_rec_right,x_rec_left,beta_rec_right, beta_rec_left;
    double t_rec_proton_right, t_rec_proton_left;
    double jet1_rec_pt, jet1_rec_eta, jet2_rec_pt, jet2_rec_eta;
    bool rp_right, rp_left;
    tree_non_diff->SetBranchAddress("xi_rec_cms_right",&xi_rec_cms_minus);
    tree_non_diff->SetBranchAddress("xi_rec_cms_left",&xi_rec_cms_plus);
    tree_non_diff->SetBranchAddress("x_rec_right",&x_rec_right);
    tree_non_diff->SetBranchAddress("x_rec_left",&x_rec_left);
    tree_non_diff->SetBranchAddress("beta_rec_right",&beta_rec_right);
    tree_non_diff->SetBranchAddress("beta_rec_left",&beta_rec_left);
    tree_non_diff->SetBranchAddress("xi_rec_totem_right",&xi_rec_proton_right);
    tree_non_diff->SetBranchAddress("xi_rec_totem_left",&xi_rec_proton_left);
    tree_non_diff->SetBranchAddress("t_rec_totem_right",&t_rec_proton_right);
    tree_non_diff->SetBranchAddress("t_rec_totem_left",&t_rec_proton_left);
    tree_non_diff->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt);
    tree_non_diff->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta);
    tree_non_diff->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt);
    tree_non_diff->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta);
    tree_non_diff->SetBranchAddress("rp_right",&rp_right);
    tree_non_diff->SetBranchAddress("rp_left",&rp_left);
    tree_non_diff->SetBranchAddress("weight",&weight);

    for(int i_evt = 0; i_evt < nev; ++i_evt){
        tree_non_diff->GetEntry(i_evt);
      
        bool signal_right = rp_right && !rp_left; 
        bool signal_left = rp_left && !rp_right;
        bool jet_rec_sel = jet1_rec_pt>40 && jet2_rec_pt>40 && fabs(jet1_rec_eta)<4.4 && fabs(jet2_rec_eta)<4.4;

        int i_evt_backg = 0 + rdm->Rndm()*(pu.size());
        RdmDataEvent const & pu_zb = pu.at(i_evt_backg);
        double xi_right_totem_pu = pu_zb.xi_totem_right_rdm;
        double t_right_totem_pu = pu_zb.t_totem_right_rdm;
        double xi_left_totem_pu = pu_zb.xi_totem_left_rdm;
        double t_left_totem_pu = pu_zb.t_totem_left_rdm;
        bool rp_right_pu = pu_zb.rp_right_rdm;
        bool rp_left_pu = pu_zb.rp_left_rdm;
        bool valid_vtx_pu = pu_zb.valid_vtx_rdm;
        bool valid_proton_right_pu = pu_zb.valid_proton_right_rdm;
        bool valid_proton_left_pu = pu_zb.valid_proton_left_rdm;
        double beta_right_pu = pu_zb.beta_right_rdm;
        double beta_left_pu = pu_zb.beta_left_rdm;
        double x_right_pu = pu_zb.x_right_rdm;
        double x_left_pu = pu_zb.x_left_rdm;

        bool full_right = /*!valid_vtx_zb && */valid_proton_right_pu && rp_right_pu && !valid_proton_left_pu && !rp_left_pu;
        bool full_left = /*!valid_vtx_zb && */valid_proton_left_pu && rp_left_pu && !valid_proton_right_pu && !rp_right_pu;
	
	
        if (!signal_right && jet_rec_sel && xi_right_totem_pu>0 && xi_right_totem_pu<0.1 && t_right_totem_pu>0.03 && t_right_totem_pu<1. ){
            if (full==true && full_right){
              if (vtx_sel && valid_vtx_zb) continue;
               TH1F_histos["xi_cms_minus_totem_right_backg"]->Fill(xi_rec_cms_minus - xi_right_totem_pu, weight);
               TH1F_histos["xi_cms_minus_totem_right_backg_bin"]->Fill(xi_rec_cms_minus - xi_right_totem_pu, weight);
               TH1F_histos["xi_right_backg"]->Fill(xi_right_totem_pu, weight);
               TH1F_histos["t_right_backg"]->Fill(t_right_totem_pu, weight);
               if(xi_rec_cms_minus - xi_right_totem_pu<0){
                 TH1F_histos["xi_right_cut_backg"]->Fill(xi_right_totem_pu, weight);
                 TH1F_histos["t_right_cut_backg"]->Fill(t_right_totem_pu, weight);
                 TH1F_histos["log_x_right_cut_backg"]->Fill(log10(x_right_pu),weight);
                 TH1F_histos["beta_right_cut_backg"]->Fill(beta_right_pu, weight);
               }
            }
            if (full==false){
               TH1F_histos["xi_cms_minus_totem_right_backg"]->Fill(xi_rec_cms_minus - xi_right_totem_pu, weight);
               TH1F_histos["xi_cms_minus_totem_right_backg_bin"]->Fill(xi_rec_cms_minus - xi_right_totem_pu, weight);
               TH1F_histos["xi_right_backg"]->Fill(xi_right_totem_pu, weight);
               TH1F_histos["t_right_backg"]->Fill(t_right_totem_pu, weight);
               if(xi_rec_cms_minus - xi_right_totem_pu<0){++nevent_backg_right;
                 TH1F_histos["xi_right_cut_backg"]->Fill(xi_right_totem_pu, weight);
                 TH1F_histos["t_right_cut_backg"]->Fill(t_right_totem_pu, weight); 
                 TH1F_histos["log_x_right_cut_backg"]->Fill(log10(x_right_pu), weight);
                 TH1F_histos["beta_right_cut_backg"]->Fill(beta_right_pu, weight);
               }
            }      
        }
        if (!signal_left && jet_rec_sel && xi_left_totem_pu>0 && xi_left_totem_pu<0.1 && t_left_totem_pu>0.03 && t_left_totem_pu<1.){
            if (full==true && full_left){
              if (vtx_sel && valid_vtx_zb) continue;
               TH1F_histos["xi_cms_minus_totem_left_backg"]->Fill(xi_rec_cms_plus - xi_left_totem_pu, weight);
               TH1F_histos["xi_cms_minus_totem_left_backg_bin"]->Fill(xi_rec_cms_plus - xi_left_totem_pu, weight);
               TH1F_histos["xi_left_backg"]->Fill(xi_left_totem_pu, weight);
               TH1F_histos["t_left_backg"]->Fill(t_left_totem_pu, weight);
               if(xi_rec_cms_plus - xi_left_totem_pu<0){
                 TH1F_histos["xi_left_cut_backg"]->Fill(xi_left_totem_pu, weight);
                 TH1F_histos["t_left_cut_backg"]->Fill(t_left_totem_pu, weight);
                 TH1F_histos["log_x_left_cut_backg"]->Fill(log10(x_left_pu), weight);
                 TH1F_histos["beta_left_cut_backg"]->Fill(beta_left_pu, weight);
               }
            }
            if (full==false){
               TH1F_histos["xi_cms_minus_totem_left_backg"]->Fill(xi_rec_cms_plus - xi_left_totem_pu, weight);
               TH1F_histos["xi_cms_minus_totem_left_backg_bin"]->Fill(xi_rec_cms_plus - xi_left_totem_pu, weight);
               TH1F_histos["xi_left_backg"]->Fill(xi_left_totem_pu, weight);
               TH1F_histos["t_left_backg"]->Fill(t_left_totem_pu, weight);
               if(xi_rec_cms_plus - xi_left_totem_pu<0){++nevent_backg_left;
                 TH1F_histos["xi_left_cut_backg"]->Fill(xi_left_totem_pu, weight);
                 TH1F_histos["t_left_cut_backg"]->Fill(t_left_totem_pu, weight);
                 TH1F_histos["log_x_left_cut_backg"]->Fill(log10(x_left_pu), weight);
                 TH1F_histos["beta_left_cut_backg"]->Fill(beta_left_pu, weight);
               }
            }
               
        }
    }

    for(map<string,TH1F*>::const_iterator histos = TH1F_histos.begin(); histos != TH1F_histos.end(); ++histos){
        histos->second->Write();
    } 

    outfile->Close();  

}

void BackgroundZB::BackgNormalisation(TH1F* data, TH1F* zb, double &norm){

  //input the xi_cms - xi_totem histos
  TAxis *axis = data->GetXaxis();
  int bmin = axis->FindBin(0.0); 
  int bmax = axis->FindBin(0.4); 

  double area_data = data->Integral(bmin, bmax);
  double area_pu = zb->Integral(bmin, bmax);
  norm = area_data/area_pu;

}

// int main(void)
// {
// 	BackgroundZB zb_background;
// 	zb_background.getBackgroundHistos(false, false);
// 	return 0;
// }