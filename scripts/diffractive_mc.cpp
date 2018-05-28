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
#include <TObject.h>
//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;

///ROOUNFOLD
#include "/Users/lina/cernbox/doctorado/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/lina/cernbox/doctorado/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/lina/cernbox/doctorado/RooUnfold/src/RooUnfoldBinByBin.h"
#include "/Users/lina/cernbox/doctorado/RooUnfold/src/RooUnfoldInvert.h"

Double_t beta_fit(Double_t *x, Double_t *par ){
  Double_t result = 0;
  result = par[0]+par[1]*x[0]+ par[2]*x[0]*x[0];
  return result;
}
Double_t fFermiLike(Double_t *x, Double_t *par) {
  Double_t result = 0;
  result = 1.0/(TMath::Exp((par[0]-TMath::Sqrt(x[0]))/par[1]) + 1);
  return result;
}


double corr_thx_vs_xi [5][18]; //(xi, thx) 
void xi_thx_corr(double xi_proton, double thx_proton, double &corr){

    float xi_bins[12] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2};
    float thetax_bins[21] = {-0.0003, -0.00027, -0.00024, -0.00021, -0.00018, -0.00015, -0.00012, -0.00009, -0.00006, -0.00003, 0, 0.00003,
                            0.00006, 0.00009, 0.00012, 0.00015, 0.00018, 0.00021, 0.00024, 0.00027, 0.0003};

    for (int i = 2; i <= 6; ++i){
         for (int j = 2; j <= 19; ++j){
              if (xi_proton<xi_bins[i] && xi_proton>xi_bins[i-1] && thx_proton<thetax_bins[j] && thx_proton>thetax_bins[j-1]) corr = corr_thx_vs_xi[i-2][j-2];
              if (xi_proton<xi_bins[i] && xi_proton>xi_bins[i-1] && thx_proton>thetax_bins[19]) corr = corr_thx_vs_xi[i-2][17];
              if (xi_proton>xi_bins[6] && thx_proton<thetax_bins[j] && thx_proton>thetax_bins[j-1]) corr = corr_thx_vs_xi[4][j-2];
        }
    }
}


class MCAnalysis
{
public:
	MCAnalysis();
	~MCAnalysis();
	void getMCHistos( string const& mc = "pomwig_pom", bool reweight = false, bool reweight_slope = false, bool unc_rp_y_up = false, 
		bool unc_rp_y_dw = false, bool unc_gauss = false, bool unc_rp_x = false, bool trigger_up = false, bool trigger_dw = false);
	// void getMCHistos( TH1F* &t_rec_right, string const& mc = "pomwig_regg", bool reweight = false, bool reweight_slope = false, bool unc_rp_y_up = false, 
	// 	bool unc_rp_y_dw = false, bool unc_gauss = false, bool unc_rp_x = false);
	// void getMCHistos( TH1F* &t_rec_right, string const& mc = "pythia8_4c", bool reweight = false, bool reweight_slope = false, bool unc_rp_y_up = false, 
	// 	bool unc_rp_y_dw = false, bool unc_gauss = false, bool unc_rp_x = false);
	// void getMCHistos( TH1F* &t_rec_right, string const& mc = "pythia8_cuetp8m1", bool reweight = false, bool reweight_slope = false, bool unc_rp_y_up = false, 
	// 	bool unc_rp_y_dw = false, bool unc_gauss = false, bool unc_rp_x = false);

    void addHistosPomwigPomRegg ();
    void getUnfoldedHistos (TH1F* reco, TH1F* truth, TH2F* response, TH1F* data, TH1F* &histo_unfolded, int n_iter = 4);
	void getUnfoldingTest (TH1F* reco, TH1F* truth, TH2F* response, TH1F* data, TH2F* &chi2_unfolded_vs_iter, TH2F* &delta_chi2_unfolded_vs_iter,
	double &chi2_smeared);

private:
	TFile * outfile;
};

MCAnalysis::MCAnalysis (void){
	cout << "File with histos from MC is created" << endl;
}

MCAnalysis::~MCAnalysis (void){
}


void MCAnalysis::getMCHistos ( string const& mc, bool reweight, bool reweight_slope, bool unc_rp_y_up, bool unc_rp_y_dw, 
		bool unc_gauss,  bool unc_rp_x, bool trigger_up, bool trigger_dw){

    string treeName = "small_tree"; 
    string fileName, type, rew;

    // correction due proton reconstruction
    bool rec_corr = true; double max = -10;
    if (mc == "pomwig_regg") rec_corr = false; 
    TFile* accep_thx_vs_xi = TFile::Open("~/cernbox/doctorado/note/scripts/final_analysis/accep_param_vs_fullsim_gen.root","READ");
    TH2F* accep_thx_vs_xi_nojet_right = (TH2F*)accep_thx_vs_xi->Get("theta_x_vs_xi_gen_right_rp_nojet");
    double sum = 0;
    for (int i = 2; i <= 6; ++i)
    { 
        for (int j = 2; j <= 19; ++j)
        {
            // cout<<i<<" "<<j<<" "<<accep_thx_vs_xi_nojet_right->GetBinContent(i,j)<<endl;
            corr_thx_vs_xi[i-2][j-2] = accep_thx_vs_xi_nojet_right->GetBinContent(i,j);
            sum += accep_thx_vs_xi_nojet_right->GetBinContent(i,j);
            if (corr_thx_vs_xi[i-2][j-2] > max) max = corr_thx_vs_xi[i-2][j-2];
        }   
    }accep_thx_vs_xi->Close();
    cout << "corr: " << sum/90 << "  "<< max<<endl;
    //...........................

    float tbins[9] = {0.03, 0.07, 0.11, 0.21, 0.31, 0.45,  0.58, 0.72, 1.};
    float xi_bins[12] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2};
    float bin_sasha[9] ={0.0003, 0.002, 0.0045, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1};
    float bin[16] = {-0.4, -0.112, -0.096, -0.08, -0.064, -0.048, -0.032, -0.016, 0, 0.048, 0.112, 0.176, 0.24, 0.304, 0.368, 0.4};
    float thetax_bins[21] = {-0.0003, -0.00027, -0.00024, -0.00021, -0.00018, -0.00015, -0.00012, -0.00009, -0.00006, -0.00003, 0, 0.00003,
                            0.00006, 0.00009, 0.00012, 0.00015, 0.00018, 0.00021, 0.00024, 0.00027, 0.0003};

    map<string,TH1F*> TH1F_right_histos;
    map<string,TH1F*> TH1F_left_histos;
    map<string,TH2F*> TH2F_right_histos;
    map<string,TH2F*> TH2F_left_histos;
    TF1* func = new TF1("func", beta_fit, 0., 1., 3);
    double t_slope_data = 5.12750;
    double t_slope_mc;
    ////// normalization factors
    double norm_pomwig_minus_pom = (0.049154*4.04e+06)/1.90724e+06;
    double norm_pomwig_minus_regg = 0.458476;
    double norm_pomwig_plus_pom = (0.049154*4.04e+06)/1.83134e+06;
    double norm_pomwig_plus_regg = 0.303177;
    double norm_pythia_4c = 2109.82*1.7e-4; 
    double norm_pythia_cuetp8m1 = 15790.5*4.3e-5;
    double scale_rew_right, scale_rew_left, scale_rew_slope_right, scale_rew_slope_left;

    double eff_proton = 0.94; // proton reconstruction efficiency

    if (mc == "pomwig_pom" || mc == "pomwig_regg") {
    	func->SetParameter(0, 1.47798);
        func->SetParameter(1, -4.08005);
        func->SetParameter(2, 4.59568);
		/////////// reweight slope
		t_slope_mc = 5.31985;
        scale_rew_right = 1.06364355;
        scale_rew_left = 1.0633095;
        scale_rew_slope_right = 0.997521;
        scale_rew_slope_left = 0.997532;
    }
 
	if (mc == "pythia8_4c"){  
	    func->SetParameter(0, 1.40943);//1.39861);//1.6199);
	    func->SetParameter(1, -3.39277);//-3.46778);//-4.38166);
	    func->SetParameter(2, 3.57715);//3.70415);//5.45443);
	    scale_rew_right = 1.032153795;
	    scale_rew_left = 1.03091912;
	    t_slope_mc = 6.08193;
	    scale_rew_slope_right = 0.978077;
	    scale_rew_slope_left = 0.9740678;
	} 
	if (mc == "pythia8_cuetp8m1" ){  
	    func->SetParameter(0, 1.36136);//1.39861);//1.6199);
	    func->SetParameter(1, -2.99797);//-3.46778);//-4.38166);
	    func->SetParameter(2, 3.14615);//3.70415);//5.45443);
	    scale_rew_right = 1.0333083;//0.8854668;
	    scale_rew_left = 1.0352818;//0.89104998;
	    t_slope_mc = 6.07664;
	    scale_rew_slope_right = 0.972148;
	    scale_rew_slope_left = 0.9770;
	}

    TH1F_right_histos["correction"] = new TH1F("correction","",20,-1, 2);
	TH1F_right_histos["xi_cms_minus_totem_rec_right"] = new TH1F("xi_cms_minus_totem_rec_right","",50,-0.4,0.4);
	TH1F_right_histos["xi_cms_minus_totem_rec_right_top"] = new TH1F("xi_cms_minus_totem_rec_right_top","",50,-0.4,0.4);
	TH1F_right_histos["xi_cms_minus_totem_rec_right_bottom"] = new TH1F("xi_cms_minus_totem_rec_right_bottom","",50,-0.4,0.4);
	TH1F_right_histos["xi_cms_minus_totem_rec_right_bin"] = new TH1F("xi_cms_minus_totem_rec_right_bin","", 15, bin);
	TH1F_right_histos["xi_rec_right"] = new TH1F("xi_rec_right","",50,-0.04,0.2);
	TH1F_right_histos["xi_rec_right_cut"] = new TH1F("xi_rec_right_cut","",11,xi_bins);
    TH1F_right_histos["xi_gen_right_sasha"] = new TH1F("xi_gen_right_sasha","",8, bin_sasha);
	TH1F_right_histos["xi_rec_right_cut_top"] = new TH1F("xi_rec_right_cut_top","",11,xi_bins);
	TH1F_right_histos["xi_rec_right_cut_bottom"] = new TH1F("xi_rec_right_cut_bottom","",11,xi_bins);
	TH1F_right_histos["xi_rec_right_cut_nojet"] = new TH1F("xi_rec_right_cut_nojet","",11,xi_bins);
    TH1F_right_histos["xi_gen_right_cut_nojet"] = new TH1F("xi_gen_right_cut_nojet","",11,xi_bins);
    TH1F_right_histos["xi_gen_right_cut_nojet_top"] = new TH1F("xi_gen_right_cut_nojet_top","",11,xi_bins);
	TH1F_right_histos["xi_gen_right_cut_nojet_bottom"] = new TH1F("xi_gen_right_cut_nojet_bottom","",11,xi_bins);
    TH1F_right_histos["xi_gen_rp_right_cut_nojet"] = new TH1F("xi_gen_rp_right_cut_nojet","",11,xi_bins);
    TH1F_right_histos["xi_gen_rp_right_cut_nojet_top"] = new TH1F("xi_gen_rp_right_cut_nojet_top","",11,xi_bins);
    TH1F_right_histos["xi_gen_rp_right_cut_nojet_bottom"] = new TH1F("xi_gen_rp_right_cut_nojet_bottom","",11,xi_bins);
	TH1F_right_histos["xi_gen_rp_right_cut_nojet_aftercorr"] = new TH1F("xi_gen_rp_right_cut_nojet_aftercorr","",11,xi_bins);
	TH1F_right_histos["xi_gen_right"] = new TH1F("xi_gen_right_cut","",11,xi_bins);
	TH1F_right_histos["t_rec_right_cut"] = new TH1F("t_rec_right_cut","", 8, tbins);
	TH1F_right_histos["t_rec_right_cut_top"] = new TH1F("t_rec_right_cut_top","", 8, tbins);
	TH1F_right_histos["t_rec_right_cut_bottom"] = new TH1F("t_rec_right_cut_bottom","", 8, tbins);
    TH1F_right_histos["t_rec_right_cut_nojet"] = new TH1F("t_rec_right_cut_nojet","", 8, tbins);
    TH1F_right_histos["t_gen_right_cut_nojet"] = new TH1F("t_gen_right_cut_nojet","", 8, tbins);
    TH1F_right_histos["t_gen_right_cut_nojet_top"] = new TH1F("t_gen_right_cut_nojet_top","", 8, tbins);
    TH1F_right_histos["t_gen_right_cut_nojet_bottom"] = new TH1F("t_gen_right_cut_nojet_bottom","", 8, tbins);
    TH1F_right_histos["t_gen_rp_right_cut_nojet"] = new TH1F("t_gen_rp_right_cut_nojet","", 8, tbins);
    TH1F_right_histos["t_gen_rp_right_cut_nojet_top"] = new TH1F("t_gen_rp_right_cut_nojet_top","", 8, tbins);
    TH1F_right_histos["t_gen_rp_right_cut_nojet_bottom"] = new TH1F("t_gen_rp_right_cut_nojet_bottom","", 8, tbins);
    TH1F_right_histos["t_gen_rp_right_cut_nojet_aftercorr"] = new TH1F("t_gen_rp_right_cut_nojet_aftercorr","", 8, tbins);
    TH1F_right_histos["t_gen_right_cut"] = new TH1F("t_gen_right_cut","", 8, tbins);
    TH1F_right_histos["beta_right_cut"] = new TH1F("beta_right_cut","",15,0,1);
    TH1F_right_histos["xi_cms_rec_right"] = new TH1F("xi_cms_rec_right","",11,xi_bins);
    TH1F_right_histos["xi_cms_gen_right"] = new TH1F("xi_cms_gen_right","",11,xi_bins);
    TH1F_right_histos["xi_cms_rec_right_sasha"] = new TH1F("xi_cms_rec_right_sasha","",8, bin_sasha);
    TH1F_right_histos["xi_cms_gen_right_sasha"] = new TH1F("xi_cms_gen_right_sasha","",8, bin_sasha);
    TH1F_right_histos["log_x_rec_right_cut"] = new TH1F("log_x_rec_right_cut","",15, -4, 0);
    TH1F_right_histos["log_x_rec_right_cut_top"] = new TH1F("log_x_rec_right_cut_top","",15, -4, 0);
    TH1F_right_histos["log_x_rec_right_cut_bottom"] = new TH1F("log_x_rec_right_cut_bottom","",15, -4, 0);
    TH1F_right_histos["log_x_gen_right_cut"] = new TH1F("log_x_gen_right_cut","",15, -4, 0);
    TH1F_right_histos["pt_jet1_right_cut"] = new TH1F("pt_jet1_right_cut","", 15, 0, 200);
    TH1F_right_histos["pt_jet2_right_cut"] = new TH1F("pt_jet2_right_cut","", 15, 0, 200);
    TH1F_right_histos["eta_jet1_right_cut"] = new TH1F("eta_jet1_right_cut","", 20, -5.2, 5.2);
    TH1F_right_histos["eta_jet2_right_cut"] = new TH1F("eta_jet2_right_cut","", 20, -5.2, 5.2);
    TH1F_right_histos["xi_rec_right_cut_sasha"] = new TH1F("xi__rec_right_sasha","",8, bin_sasha);
    TH1F_right_histos["xi_rec_cms_minus_sasha"] = new TH1F("xi_rec_cms_minus_sasha","",8, bin_sasha);
    TH1F_right_histos["xi_gen_right_cut_sasha"] = new TH1F("xi__gen_right_sasha","",8, bin_sasha);
    TH1F_right_histos["xi_cms_rec_right_cut_sasha"] = new TH1F("xi_cms_rec_right_cut_sasha","",8, bin_sasha);
    TH1F_right_histos["xi_cms_gen_right_cut_sasha"] = new TH1F("xi_cms_gen_right_cut_sasha","",8, bin_sasha);
    TH1F_right_histos["x_rec_gen_right"] = new TH1F("x_rec_gen_right","",20,-0.4,0.4);
    TH1F_right_histos["th_x_rec_right"] = new TH1F("th_x_rec_right","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_rec_right"] = new TH1F("th_y_rec_right","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_x_gen_right"] = new TH1F("th_x_gen_right","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_gen_right"] = new TH1F("th_y_gen_right","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_x_gen_right_nojet"] = new TH1F("th_x_gen_right_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_x_gen_right_nojet_top"] = new TH1F("th_x_gen_right_nojet_top","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_x_gen_right_nojet_bottom"] = new TH1F("th_x_gen_right_nojet_bottom","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_x_gen_rp_right_nojet"] = new TH1F("th_x_gen_rp_right_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_x_gen_rp_right_nojet_top"] = new TH1F("th_x_gen_rp_right_nojet_top","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_x_gen_rp_right_nojet_bottom"] = new TH1F("th_x_gen_rp_right_nojet_bottom","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_x_gen_rp_right_nojet_aftercorr"] = new TH1F("th_x_gen_rp_right_nojet_aftercorr","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_gen_right_nojet"] = new TH1F("th_y_gen_right_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_gen_right_nojet_top"] = new TH1F("th_y_gen_right_nojet_top","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_gen_right_nojet_bottom"] = new TH1F("th_y_gen_right_nojet_bottom","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_gen_rp_right_nojet"] = new TH1F("th_y_gen_rp_right_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_gen_rp_right_nojet_top"] = new TH1F("th_y_gen_rp_right_nojet_top","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_gen_rp_right_nojet_bottom"] = new TH1F("th_y_gen_rp_right_nojet_bottom","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["th_y_gen_rp_right_nojet_aftercorr"] = new TH1F("th_y_gen_rp_right_nojet_aftercorr","",20, -0.4e-3, 0.4e-3);
    TH1F_right_histos["delta_eta_jets_right_cut"] = new TH1F("delta_eta_jets_right_cut","", 40, -5.2, 5.2);
    TH1F_right_histos["delta_phi_jets_right_cut"] = new TH1F("delta_phi_jets_right_cut","", 40, -5.2, 5.2);
    TH1F_right_histos["mass_x_rec_right"] = new TH1F("mass_x_right","",20, 0, 1800);
    TH1F_right_histos["mass_jj_rec_right"] = new TH1F("mass_jj_right","",20, 0, 1000);
    TH1F_right_histos["r_jj_rec_right"] = new TH1F("r_jj_right","",20, 0, 1);
    TH1F_right_histos["truth_t_right"] = new TH1F("truth_t_right","", 8, tbins);
    TH1F_right_histos["reco_t_right"] = new TH1F("reco_t_right","", 8, tbins);
    TH1F_right_histos["truth_xi_right"] = new TH1F("truth_xi_right","", 11,xi_bins);
    TH1F_right_histos["reco_xi_right"] = new TH1F("reco_xi_right","", 11,xi_bins);
    TH1F_right_histos["truth_logx_right"] = new TH1F("truth_x_right","", 15, -4, 0);
    TH1F_right_histos["reco_logx_right"] = new TH1F("reco_x_right","", 15, -4, 0);

    TH2F_right_histos["response_t_right"] = new TH2F("response_t_right", "", 8, tbins, 8, tbins);
    TH2F_right_histos["response_xi_right"] = new TH2F("response_xi_right", "", 11, xi_bins, 11, xi_bins);
    TH2F_right_histos["response_logx_right"] = new TH2F("response_logx_right", "", 15, -4, 0, 15, -4, 0);
    TH2F_right_histos["rp_pos_right"] = new TH2F("rp_pos_right","",100,-1,10,100,-40,40);
    TH2F_right_histos["pt2_xi_gen_minus"] = new TH2F("pt2_xi_gen_minus","",20,0,1,20,0,0.1);
    TH2F_right_histos["pt2_xi_gen_minus_rp"] = new TH2F("pt2_xi_gen_minus_rp","",20,0,1,20,0,0.1);
    TH2F_right_histos["t_gen_vs_xi_gen_minus"] = new TH2F("t_gen_vs_xi_gen_minus","", 8, tbins, 11, xi_bins);
    TH2F_right_histos["t_rec_vs_xi_rec_minus"] = new TH2F("t_rec_vs_xi_rec_minus","", 8, tbins, 11, xi_bins);
    TH2F_right_histos["t_minus_th2"] = new TH2F("t_minus_th2","", 8, tbins,8,tbins);
    TH2F_right_histos["xi_minus_th2"] = new TH2F("xi_minus_th2","",11,xi_bins,11,xi_bins);
    TH2F_right_histos["logx_minus_th2"] = new TH2F("logx_minus_th2","", 15, -4, 0,15, -4, 0);
    TH2F_right_histos["log_t_vs_log_xi_right"] = new TH2F("log_t_vs_log_xi_right","",100,-4,1,100,-4,-0.5);
    TH2F_right_histos["delta_thx_vs_delta_xi_right"] = new TH2F("delta_thx_vs_delta_xi_right", "", 20, -0.1e-3, 0.1e-3, 20, -0.04, 0.04);
    TH2F_right_histos["theta_x_vs_xi_gen_nojet"] = new TH2F("theta_x_vs_xi_gen_right_nojet","", 11,xi_bins, 20, thetax_bins);
    TH2F_right_histos["theta_x_vs_xi_gen_rp_nojet"] = new TH2F("theta_x_vs_xi_gen_right_rp_nojet","", 11,xi_bins, 20, thetax_bins);
    TH2F_right_histos["theta_x_vs_xi_gen_rp_nojet_aftercorr"] = new TH2F("theta_x_vs_xi_gen_right_rp_nojet_aftercorr","", 11,xi_bins, 20, thetax_bins);

    TH1F_left_histos["xi_cms_minus_totem_rec_left"] = new TH1F("xi_cms_minus_totem_rec_left","",50,-0.4,0.4);
    TH1F_left_histos["xi_cms_minus_totem_rec_left_top"] = new TH1F("xi_cms_minus_totem_rec_left_top","",50,-0.4,0.4);
    TH1F_left_histos["xi_cms_minus_totem_rec_left_bottom"] = new TH1F("xi_cms_minus_totem_rec_left_bottom","",50,-0.4,0.4);
    TH1F_left_histos["xi_cms_minus_totem_rec_left_bin"] = new TH1F("xi_cms_minus_totem_rec_left_bin","", 15, bin);
    TH1F_left_histos["xi_rec_left"] = new TH1F("xi_rec_left","",50,-0.04,0.2);
    TH1F_left_histos["xi_rec_left_cut"] = new TH1F("xi_rec_left_cut","",11,xi_bins);
    TH1F_left_histos["xi_rec_left_cut_top"] = new TH1F("xi_rec_left_cut_top","",11,xi_bins);
    TH1F_left_histos["xi_rec_left_cut_bottom"] = new TH1F("xi_rec_left_cut_bottom","",11,xi_bins);
    TH1F_left_histos["xi_rec_left_cut_nojet"] = new TH1F("xi_rec_left_cut_nojet","",11,xi_bins);
    TH1F_left_histos["xi_gen_left_cut_nojet"] = new TH1F("xi_gen_left_cut_nojet","",11,xi_bins);
    TH1F_left_histos["xi_gen_left_cut_nojet_top"] = new TH1F("xi_gen_left_cut_nojet_top","",11,xi_bins);
    TH1F_left_histos["xi_gen_left_cut_nojet_bottom"] = new TH1F("xi_gen_left_cut_nojet_bottom","",11,xi_bins);
    TH1F_left_histos["xi_gen_rp_left_cut_nojet"] = new TH1F("xi_gen_rp_left_cut_nojet","",11,xi_bins);
    TH1F_left_histos["xi_gen_rp_left_cut_nojet_top"] = new TH1F("xi_gen_rp_left_cut_nojet_top","",11,xi_bins);
    TH1F_left_histos["xi_gen_rp_left_cut_nojet_bottom"] = new TH1F("xi_gen_rp_left_cut_nojet_bottom","",11,xi_bins);
    TH1F_left_histos["xi_gen_rp_left_cut_nojet_aftercorr"] = new TH1F("xi_gen_rp_left_cut_nojet_aftercorr","",11,xi_bins);
    TH1F_left_histos["xi_gen_left"] = new TH1F("xi_gen_left_cut","",11,xi_bins);
    TH1F_left_histos["t_rec_left_cut"] = new TH1F("t_rec_left_cut","", 8, tbins);
    TH1F_left_histos["t_rec_left_cut_top"] = new TH1F("t_rec_left_cut_top","", 8, tbins);
    TH1F_left_histos["t_rec_left_cut_bottom"] = new TH1F("t_rec_left_cut_bottom","", 8, tbins);
    TH1F_left_histos["t_rec_left_cut_nojet"] = new TH1F("t_rec_left_cut_nojet","", 8, tbins);
    TH1F_left_histos["t_gen_left_cut_nojet"] = new TH1F("t_gen_left_cut_nojet","", 8, tbins);
    TH1F_left_histos["t_gen_left_cut_nojet_top"] = new TH1F("t_gen_left_cut_nojet_top","", 8, tbins);
    TH1F_left_histos["t_gen_left_cut_nojet_bottom"] = new TH1F("t_gen_left_cut_nojet_bottom","", 8, tbins);
    TH1F_left_histos["t_gen_rp_left_cut_nojet"] = new TH1F("t_gen_rp_left_cut_nojet","", 8, tbins);
    TH1F_left_histos["t_gen_rp_left_cut_nojet_top"] = new TH1F("t_gen_rp_left_cut_nojet_top","", 8, tbins);
    TH1F_left_histos["t_gen_rp_left_cut_nojet_bottom"] = new TH1F("t_gen_rp_left_cut_nojet_bottom","", 8, tbins);
    TH1F_left_histos["t_gen_rp_left_cut_nojet_aftercorr"] = new TH1F("t_gen_rp_left_cut_nojet_aftercorr","", 8, tbins);
    TH1F_left_histos["t_gen_left_cut"] = new TH1F("t_gen_left_cut","", 8, tbins);
    TH1F_left_histos["beta_left_cut"] = new TH1F("beta_left_cut","",15,0,1);
    TH1F_left_histos["xi_cms_rec_left"] = new TH1F("xi_cms_rec_left","",11,xi_bins);
    TH1F_left_histos["xi_cms_gen_left"] = new TH1F("xi_cms_gen_left","",11,xi_bins);
    TH1F_left_histos["xi_cms_rec_left_sasha"] = new TH1F("xi_cms_rec_left_sasha","",8, bin_sasha);
    TH1F_left_histos["xi_cms_gen_left_sasha"] = new TH1F("xi_cms_gen_left_sasha","",8, bin_sasha);
    TH1F_left_histos["log_x_rec_left_cut"] = new TH1F("log_x_rec_left_cut","",15, -4, 0);
    TH1F_left_histos["log_x_rec_left_cut_top"] = new TH1F("log_x_rec_left_cut_top","",15, -4, 0);
    TH1F_left_histos["log_x_rec_left_cut_bottom"] = new TH1F("log_x_rec_left_cut_bottom","",15, -4, 0);
    TH1F_left_histos["log_x_gen_left_cut"] = new TH1F("log_x_gen_left_cut","",15, -4, 0);
    TH1F_left_histos["pt_jet1_left_cut"] = new TH1F("pt_jet1_left_cut","", 15, 0, 200);
    TH1F_left_histos["pt_jet2_left_cut"] = new TH1F("pt_jet2_left_cut","", 15, 0, 200);
    TH1F_left_histos["eta_jet1_left_cut"] = new TH1F("eta_jet1_left_cut","", 20, -5.2, 5.2);
    TH1F_left_histos["eta_jet2_left_cut"] = new TH1F("eta_jet2_left_cut","", 20, -5.2, 5.2);
    TH1F_left_histos["xi_rec_left_cut_sasha"] = new TH1F("xi__rec_left_sasha","",11,xi_bins);
    TH1F_left_histos["xi_gen_left_cut_sasha"] = new TH1F("xi__gen_left_sasha","",11,xi_bins);
    TH1F_left_histos["xi_cms_rec_left_cut_sasha"] = new TH1F("xi_cms_rec_left_cut_sasha","",8, bin_sasha);
    TH1F_left_histos["xi_cms_gen_left_cut_sasha"] = new TH1F("xi_cms_gen_left_cut_sasha","",8, bin_sasha);
    TH1F_left_histos["x_rec_gen_left"] = new TH1F("x_rec_gen_left","",20,-0.4,0.4);
    TH1F_left_histos["th_x_rec_left"] = new TH1F("th_x_rec_left","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_rec_left"] = new TH1F("th_y_rec_left","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_x_gen_left"] = new TH1F("th_x_gen_left","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_gen_left"] = new TH1F("th_y_gen_left","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_x_gen_left_nojet"] = new TH1F("th_x_gen_left_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_x_gen_left_nojet_top"] = new TH1F("th_x_gen_left_nojet_top","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_x_gen_left_nojet_bottom"] = new TH1F("th_x_gen_left_nojet_bottom","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_x_gen_rp_left_nojet"] = new TH1F("th_x_gen_rp_left_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_x_gen_rp_left_nojet_top"] = new TH1F("th_x_gen_rp_left_nojet_top","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_x_gen_rp_left_nojet_bottom"] = new TH1F("th_x_gen_rp_left_nojet_bottom","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_x_gen_rp_left_nojet_aftercorr"] = new TH1F("th_x_gen_rp_left_nojet_aftercorr","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_gen_left_nojet"] = new TH1F("th_y_gen_left_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_gen_left_nojet_top"] = new TH1F("th_y_gen_left_nojet_top","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_gen_left_nojet_bottom"] = new TH1F("th_y_gen_left_nojet_bottom","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_gen_rp_left_nojet"] = new TH1F("th_y_gen_rp_left_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_gen_rp_left_nojet_top"] = new TH1F("th_y_gen_rp_left_nojet_top","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_gen_rp_left_nojet_bottom"] = new TH1F("th_y_gen_rp_left_nojet_bottom","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["th_y_gen_rp_left_nojet_aftercorr"] = new TH1F("th_y_gen_rp_left_nojet_aftercorr","",20, -0.4e-3, 0.4e-3);
    TH1F_left_histos["delta_eta_jets_left_cut"] = new TH1F("delta_eta_jets_left_cut","", 40, -5.2, 5.2);
    TH1F_left_histos["delta_phi_jets_left_cut"] = new TH1F("delta_phi_jets_left_cut","", 40, -5.2, 5.2);
    TH1F_left_histos["mass_x_rec_left"] = new TH1F("mass_x_left","",20, 0, 1800);
    TH1F_left_histos["mass_jj_rec_left"] = new TH1F("mass_jj_left","",20, 0, 1000);
    TH1F_left_histos["r_jj_rec_left"] = new TH1F("r_jj_left","",20, 0, 1);
    TH1F_left_histos["truth_t_left"] = new TH1F("truth_t_left","", 8, tbins);
    TH1F_left_histos["reco_t_left"] = new TH1F("reco_t_left","", 8, tbins);
    TH1F_left_histos["truth_xi_left"] = new TH1F("truth_xi_left","", 11,xi_bins);
    TH1F_left_histos["reco_xi_left"] = new TH1F("reco_xi_left","", 11,xi_bins);
    TH1F_left_histos["truth_logx_left"] = new TH1F("truth_x_left","", 15, -4, 0);
    TH1F_left_histos["reco_logx_left"] = new TH1F("reco_x_left","", 15, -4, 0);
    TH2F_left_histos["theta_x_vs_xi_gen_nojet"] = new TH2F("theta_x_vs_xi_gen_left_nojet","", 11,xi_bins, 20, thetax_bins);
    TH2F_left_histos["theta_x_vs_xi_gen_rp_nojet"] = new TH2F("theta_x_vs_xi_gen_rp_left_nojet","", 11,xi_bins, 20, thetax_bins);
    TH2F_left_histos["theta_x_vs_xi_gen_rp_nojet_aftercorr"] = new TH2F("theta_x_vs_xi_gen_rp_left_nojet_aftercorr","", 11,xi_bins, 20, thetax_bins);

    TH2F_left_histos["response_t_left"] = new TH2F("response_t_left", "", 8, tbins, 8, tbins);
    TH2F_left_histos["response_xi_left"] = new TH2F("response_xi_left", "", 11, xi_bins, 11, xi_bins);
    TH2F_left_histos["response_logx_left"] = new TH2F("response_logx_left", "", 15, -4, 0, 15, -4, 0);
    TH2F_left_histos["rp_pos_left"] = new TH2F("rp_pos_left","",100,-1,10,100,-40,40);
    TH2F_left_histos["pt2_xi_gen_plus"] = new TH2F("pt2_xi_gen_plus","",20,0,1,20,0,0.1);
    TH2F_left_histos["pt2_xi_gen_plus_rp"] = new TH2F("pt2_xi_gen_plus_rp","",20,0,1,20,0,0.1);
    TH2F_left_histos["t_gen_vs_xi_gen_plus"] = new TH2F("t_gen_vs_xi_gen_plus","", 8, tbins, 11, xi_bins);
    TH2F_left_histos["t_rec_vs_xi_rec_plus"] = new TH2F("t_rec_vs_xi_rec_plus","", 8, tbins, 11, xi_bins);
    TH2F_left_histos["t_plus_th2"] = new TH2F("t_plus_th2","", 8, tbins,8,tbins);
    TH2F_left_histos["xi_plus_th2"] = new TH2F("xi_plus_th2","",11,xi_bins,11,xi_bins);
    TH2F_left_histos["logx_plus_th2"] = new TH2F("logx_plus_th2","", 15, -4, 0,15, -4, 0);
    TH2F_left_histos["log_t_vs_log_xi_left"] = new TH2F("log_t_vs_log_xi_left","",100,-4,1,100,-4,-0.5);
    TH2F_left_histos["delta_thx_vs_delta_xi_left"] = new TH2F("delta_thx_vs_delta_xi_left", "", 20, -0.1e-3, 0.1e-3, 20, -0.04, 0.04);
    
    for(map<string,TH1F*>::const_iterator histos = TH1F_right_histos.begin(); histos != TH1F_right_histos.end(); ++histos){
       histos->second->Sumw2();
    }
    for(map<string,TH1F*>::const_iterator histos = TH1F_left_histos.begin(); histos != TH1F_left_histos.end(); ++histos){
       histos->second->Sumw2();
    }

    //Open ntuples
    TFile* pomwig_pom = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pomwig_minus_ntuple_beam_smearing.root","READ");
    TFile* pomwig_regg = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pomwig_reggeon_minus_ntuple_beam_smearing.root","READ");
    TFile* pomwig_pom_plus = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pomwig_plus_ntuple_beam_smearing.root","READ");
    TFile* pomwig_regg_plus = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pomwig_reggeon_plus_ntuple_beam_smearing.root","READ");
    TFile* pythia8 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pythia8_diff_ntuple_beam_smearing.root","READ");
    TFile* pythia8_CUETP8M1 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pythia8_diff_CUETP8M1_ntuple_beam_smearing.root","READ");

    ///trigger efficiency
    TF1* func_trigger = new TF1("func_trigger", fFermiLike, 0., 200, 2);
    func_trigger->SetParameter(0,4.992);
    func_trigger->SetParameter(1,0.4885);

    if (trigger_up){
        func_trigger->SetParameter(0, 4.992 + 0.06335);
        func_trigger->SetParameter(1, 0.4885 + 0.04156);
    }
    if (trigger_dw){
        func_trigger->SetParameter(0, 4.992 - 0.06335);
        func_trigger->SetParameter(1, 0.4885 - 0.04156);
    }

    //output file
    if (unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_gauss && !trigger_up && !trigger_dw) type = "unc_rp_y_up";
    else if (!unc_rp_y_up && unc_rp_y_dw && !unc_rp_x && !unc_gauss && !trigger_up && !trigger_dw) type = "unc_rp_y_dw";
    else if (!unc_rp_y_up && !unc_rp_y_dw && unc_rp_x && !unc_gauss && !trigger_up && !trigger_dw) type = "unc_rp_x";
    else if (!unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && unc_gauss && !trigger_up && !trigger_dw) type = "unc_gauss";
    else if (!unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_gauss && trigger_up && !trigger_dw) type = "unc_trigger_up";
    else if (!unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_gauss && !trigger_up && trigger_dw) type = "unc_trigger_dw";
    else type = "nominal";

    if (reweight && !reweight_slope) rew = "betarew";
    else if (!reweight && reweight_slope) rew = "sloperew";
    else if (reweight && reweight_slope) rew = "betarew_sloperew";
    else rew = "norew";

	if (mc == "pomwig_pom") fileName = "histos_pomwig_pom_" + rew + "_" + type + ".root";	
	else if (mc == "pomwig_regg") fileName = "histos_pomwig_regg_" + rew + "_" + type + ".root";	
	else if (mc == "pythia8_4c") fileName = "histos_pythia8_4c_" + rew + "_" + type + ".root";
	else if(mc == "pythia8_cuetp8m1") fileName = "histos_pythia8_cuetp8m1_" + rew + "_" + type + ".root";
	cout << "output: " << fileName << endl; 

    outfile = new TFile(fileName.c_str(), "RECREATE");
    ///.............


	TFile* file_minus; double norm_minus;
	if (mc == "pomwig_pom"){ file_minus = pomwig_pom; norm_minus = norm_pomwig_minus_pom;}
	if (mc == "pomwig_regg"){ file_minus = pomwig_regg; norm_minus = norm_pomwig_minus_regg;}
	if (mc == "pythia8_4c"){ file_minus = pythia8; norm_minus = norm_pythia_4c;}
	if (mc == "pythia8_cuetp8m1"){ file_minus = pythia8_CUETP8M1; norm_minus = norm_pythia_cuetp8m1;}

	TTree* tree_minus = (TTree*) file_minus->Get( treeName.c_str() );
    int nev = int(tree_minus->GetEntriesFast());
    cout << mc <<" minus file has " << nev << " entries : " << endl;
 
    double jet1_rec_pt, jet1_rec_eta, jet1_rec_phi, jet2_rec_pt, jet2_rec_eta, jet2_rec_phi;
    double jet1_gen_pt, jet1_gen_eta, jet1_gen_phi, jet2_gen_pt, jet2_gen_eta, jet2_gen_phi;
    double xi_rec_cms_minus, xi_rec_proton_right, x_rec_right, x_gen_right;
    double xi_gen_cms_minus, xi_gen_proton_right, xi_rec_proton_right_gauss, xi_rec_proton_right_theta;
    double t_rec_proton_right, t_rec_proton_right_theta, t_rec_proton_right_gauss, beta_rec_proton_right;
    double t_gen_proton_right, beta_gen_proton_right, rp_xpos_124, rp_xpos_125, rp_ypos_124, rp_ypos_125, t_rec_proton_right_old;
    bool rp_right, rp_right_accep_top, rp_right_accep_bottom;
    double theta_x_plus, theta_x_minus, theta_y_plus, theta_y_minus;
    double theta_x_plus_smear, theta_x_minus_smear, theta_y_plus_smear, theta_y_minus_smear, pz_minus_smear, pz_plus_smear, e_minus_smear, e_plus_smear;
    double px_minus, py_minus, pz_minus, e_minus, mass_minus, mjj2_rec_minus;
    int nVtx;
    tree_minus->SetBranchAddress("xi_rec_cms_right",&xi_rec_cms_minus);
    tree_minus->SetBranchAddress("xi_gen_cms_right",&xi_gen_cms_minus);
    tree_minus->SetBranchAddress("xi_rec_proton_right",&xi_rec_proton_right);
    tree_minus->SetBranchAddress("xi_gen_proton_right",&xi_gen_proton_right);
    tree_minus->SetBranchAddress("t_rec_proton_right",&t_rec_proton_right);
    tree_minus->SetBranchAddress("t_gen_proton_right",&t_gen_proton_right);
    tree_minus->SetBranchAddress("beta_rec_right",&beta_rec_proton_right);
    tree_minus->SetBranchAddress("beta_gen_right",&beta_gen_proton_right);
    tree_minus->SetBranchAddress("rp_right",&rp_right);
    tree_minus->SetBranchAddress("rp_right_accep_bottom",&rp_right_accep_bottom);
    tree_minus->SetBranchAddress("rp_right_accep_top",&rp_right_accep_top);
    tree_minus->SetBranchAddress("rp_xpos_124",&rp_xpos_124);
    tree_minus->SetBranchAddress("rp_ypos_124",&rp_ypos_124);
    tree_minus->SetBranchAddress("rp_xpos_125",&rp_xpos_125);
    tree_minus->SetBranchAddress("rp_ypos_125",&rp_ypos_125);
    tree_minus->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt);
    tree_minus->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt);
    tree_minus->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta);
    tree_minus->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta);
    tree_minus->SetBranchAddress("jet1_rec_phi",&jet1_rec_phi);
    tree_minus->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt);
    tree_minus->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt);
    tree_minus->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta);
    tree_minus->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta);
    tree_minus->SetBranchAddress("jet2_rec_phi",&jet2_rec_phi);
    tree_minus->SetBranchAddress("x_rec_right",&x_rec_right);
    tree_minus->SetBranchAddress("x_gen_right",&x_gen_right);
    tree_minus->SetBranchAddress("xi_rec_proton_right_gauss",&xi_rec_proton_right_gauss);
    tree_minus->SetBranchAddress("t_rec_proton_right_gauss",&t_rec_proton_right_gauss);
    tree_minus->SetBranchAddress("theta_x_minus",&theta_x_minus);
    tree_minus->SetBranchAddress("theta_y_minus",&theta_y_minus);
    tree_minus->SetBranchAddress("theta_x_minus_smear",&theta_x_minus_smear);
    tree_minus->SetBranchAddress("theta_y_minus_smear",&theta_y_minus_smear);
    tree_minus->SetBranchAddress("px_proton_right",&px_minus);
    tree_minus->SetBranchAddress("py_proton_right",&py_minus);
    tree_minus->SetBranchAddress("pz_proton_right",&pz_minus);
    tree_minus->SetBranchAddress("mass_proton_right",&mass_minus);
    tree_minus->SetBranchAddress("nVtx",&nVtx);
    tree_minus->SetBranchAddress("mjj2_rec",&mjj2_rec_minus);

double regg = 0;
    double pt_threshold = 40.;
    for(int i_evt = 0; i_evt < nev; ++i_evt){
        tree_minus->GetEntry(i_evt);
           
        // if (single_vertex && nVtx!=1) continue;
        // if (!single_vertex && nVtx<1) continue;
        if (nVtx<1) continue;

if (jet1_gen_pt>30 && jet2_gen_pt>30) ++regg;

        double eff_trigger = func_trigger->Eval(jet2_rec_pt);

        xi_rec_proton_right = (unc_gauss == false) ? xi_rec_proton_right : xi_rec_proton_right_gauss;
        t_rec_proton_right = (unc_gauss == false) ? t_rec_proton_right : t_rec_proton_right_gauss; 

        double reweigth_beta_minus = (beta_gen_proton_right<=0.7) ? func->Eval(beta_gen_proton_right) : 1;
        double reweight_slope_minus = (fabs(t_gen_proton_right<=0.45)) ? scale_rew_slope_right*t_slope_data*exp(-t_slope_data*fabs(t_gen_proton_right))/(t_slope_mc*exp(-t_slope_mc*fabs(t_gen_proton_right))) : 1;//0.95 normalization factor
        double event_weight_minus;
        if (reweight && !reweight_slope) event_weight_minus = scale_rew_right*reweigth_beta_minus; 
        if (reweight && reweight_slope) event_weight_minus = scale_rew_right*reweigth_beta_minus*reweight_slope_minus;
        if (!reweight && reweight_slope) event_weight_minus = reweight_slope_minus;
        if (!reweight && !reweight_slope) event_weight_minus = 1;

      	bool jet_rec_sel = jet1_rec_pt>pt_threshold && jet2_rec_pt>pt_threshold && fabs(jet1_rec_eta)<4.4 && fabs(jet2_rec_eta)<4.4;
        bool jet_gen_sel = jet1_gen_pt>pt_threshold && jet2_gen_pt>pt_threshold && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4;
        bool proton_rec_sel =  xi_rec_proton_right>0 && xi_rec_proton_right<0.1 && fabs(t_rec_proton_right)>0.03 && fabs(t_rec_proton_right)<1;
        bool proton_gen_sel =  xi_gen_proton_right>0 && xi_gen_proton_right<0.1 && fabs(t_gen_proton_right)>0.03 && fabs(t_gen_proton_right)<1;

        bool fid_cuts_nom_right_top = rp_xpos_124>0 && rp_xpos_124<0.007 && rp_ypos_124 >0.0084 && rp_ypos_124<0.027 ;
        bool fid_cuts_nom_right_bottom = rp_xpos_125>0 && rp_xpos_125<0.007 && rp_ypos_125 <-0.0084 && rp_ypos_125>-0.027 ;
        bool fid_cuts_unc_y_up_right_top = rp_xpos_124>0 && rp_xpos_124<0.007 && rp_ypos_124 >0.0084 && rp_ypos_124<0.0272 ;
        bool fid_cuts_unc_y_up_right_bottom = rp_xpos_125>0 && rp_xpos_125<0.007 && rp_ypos_125 <-0.0084 && rp_ypos_125>-0.0272 ;
        bool fid_cuts_unc_y_dw_right_top = rp_xpos_124>0 && rp_xpos_124<0.007 && rp_ypos_124 >0.0082 && rp_ypos_124<0.027 ;
        bool fid_cuts_unc_y_dw_right_bottom = rp_xpos_125>0 && rp_xpos_125<0.007 && rp_ypos_125 <-0.0082 && rp_ypos_125>-0.027 ;
        bool fid_cuts_unc_x_right_top = rp_xpos_124>0.00 && rp_xpos_124<0.006 && rp_ypos_124 >0.0084 && rp_ypos_124<0.027 ;
        bool fid_cuts_unc_x_right_bottom = rp_xpos_125>0.00 && rp_xpos_125<0.006 && rp_ypos_125 <-0.0084 && rp_ypos_125>-0.027 ;

        bool rp_right_sel = false; 
        if(unc_rp_y_up == false && unc_rp_y_dw == false && unc_rp_x == false) rp_right_sel = (rp_right_accep_top && fid_cuts_nom_right_top) || (rp_right_accep_bottom && fid_cuts_nom_right_bottom);
        if(unc_rp_y_up == true && unc_rp_y_dw == false && unc_rp_x == false) rp_right_sel = (rp_right_accep_top && fid_cuts_unc_y_up_right_top) || (rp_right_accep_bottom && fid_cuts_unc_y_up_right_bottom);
        if(unc_rp_y_up == false && unc_rp_y_dw == true && unc_rp_x == false) rp_right_sel = (rp_right_accep_top && fid_cuts_unc_y_dw_right_top) || (rp_right_accep_bottom && fid_cuts_unc_y_dw_right_bottom);
        if(unc_rp_y_up == false && unc_rp_y_dw == false && unc_rp_x == true) rp_right_sel = (rp_right_accep_top && fid_cuts_unc_x_right_top) || (rp_right_accep_bottom && fid_cuts_unc_x_right_bottom);

        double corr_xi = 1;//(xi_rec_proton_right_pom > 0.08) ? 1.32864 : 1.;
        if (rec_corr) xi_thx_corr(xi_gen_proton_right, theta_x_minus, corr_xi);

        TH1F_right_histos["correction"]->Fill(xi_rec_cms_minus/xi_gen_cms_minus, 1.);

        if ((rp_right_accep_top && fid_cuts_nom_right_top) && !(rp_right_accep_bottom && fid_cuts_nom_right_bottom)) TH2F_right_histos["rp_pos_right"]->Fill(rp_xpos_124*1000, rp_ypos_124*1000, 1.);
        if (!(rp_right_accep_top && fid_cuts_nom_right_top) && (rp_right_accep_bottom && fid_cuts_nom_right_bottom)) TH2F_right_histos["rp_pos_right"]->Fill(rp_xpos_125*1000, rp_ypos_125*1000, 1.);

        TH2F_right_histos["delta_thx_vs_delta_xi_right"]->Fill(theta_x_minus_smear-theta_x_minus, xi_rec_proton_right-xi_gen_proton_right, 1.);

        if(proton_gen_sel){
             TH1F_right_histos["t_gen_right_cut_nojet"]->Fill(fabs(t_gen_proton_right), 1.);
             TH1F_right_histos["xi_gen_right_cut_nojet"]->Fill(xi_gen_proton_right, 1.);
             TH2F_right_histos["t_gen_vs_xi_gen_minus"]->Fill( fabs(t_gen_proton_right), xi_gen_proton_right, 1.);
             TH1F_right_histos["t_rec_right_cut_nojet"]->Fill(fabs(t_rec_proton_right), 1.);
             TH1F_right_histos["xi_rec_right_cut_nojet"]->Fill(xi_rec_proton_right, 1.);
             TH1F_right_histos["th_x_gen_right_nojet"]->Fill(theta_x_minus, 1.);
             TH1F_right_histos["th_y_gen_right_nojet"]->Fill(theta_y_minus, 1.);
             TH2F_right_histos["theta_x_vs_xi_gen_nojet"]->Fill(xi_gen_proton_right, theta_x_minus, 1.);
             if(rp_right_sel){
               TH1F_right_histos["t_gen_rp_right_cut_nojet"]->Fill(fabs(t_gen_proton_right), 1.);
               TH1F_right_histos["xi_gen_rp_right_cut_nojet"]->Fill(xi_gen_proton_right, 1.);
               TH1F_right_histos["th_x_gen_rp_right_nojet"]->Fill(theta_x_minus, 1.);
               TH1F_right_histos["th_y_gen_rp_right_nojet"]->Fill(theta_y_minus, 1.);
               TH2F_right_histos["theta_x_vs_xi_gen_rp_nojet"]->Fill(xi_gen_proton_right, theta_x_minus, 1.);
               TH1F_right_histos["t_gen_rp_right_cut_nojet_aftercorr"]->Fill(fabs(t_gen_proton_right), 1./corr_xi);
               TH1F_right_histos["xi_gen_rp_right_cut_nojet_aftercorr"]->Fill(xi_gen_proton_right, 1./corr_xi);
               TH1F_right_histos["th_x_gen_rp_right_nojet_aftercorr"]->Fill(theta_x_minus, 1./corr_xi);
               TH1F_right_histos["th_y_gen_rp_right_nojet_aftercorr"]->Fill(theta_y_minus, 1./corr_xi);
               TH2F_right_histos["theta_x_vs_xi_gen_rp_nojet_aftercorr"]->Fill(xi_gen_proton_right, theta_x_minus, 1./corr_xi);
             }
            if(theta_y_minus >= 0){
               TH1F_right_histos["t_gen_right_cut_nojet_top"]->Fill(fabs(t_gen_proton_right), 1.);
               TH1F_right_histos["xi_gen_right_cut_nojet_top"]->Fill(xi_gen_proton_right, 1.);
               TH1F_right_histos["th_x_gen_right_nojet_top"]->Fill(theta_x_minus, 1.);
               TH1F_right_histos["th_y_gen_right_nojet_top"]->Fill(theta_y_minus, 1.);
               if (rp_right_accep_top && fid_cuts_nom_right_top){ 
                  TH1F_right_histos["t_gen_rp_right_cut_nojet_top"]->Fill(fabs(t_gen_proton_right), 1./corr_xi);
                  TH1F_right_histos["xi_gen_rp_right_cut_nojet_top"]->Fill(xi_gen_proton_right, 1./corr_xi);
                  TH1F_right_histos["th_x_gen_rp_right_nojet_top"]->Fill(theta_x_minus, 1./corr_xi);
                  TH1F_right_histos["th_y_gen_rp_right_nojet_top"]->Fill(theta_y_minus, 1./corr_xi);
                } 
            }     
            if(theta_y_minus < 0){
               TH1F_right_histos["t_gen_right_cut_nojet_bottom"]->Fill(fabs(t_gen_proton_right), 1.);
               TH1F_right_histos["xi_gen_right_cut_nojet_bottom"]->Fill(xi_gen_proton_right, 1.);
               TH1F_right_histos["th_x_gen_right_nojet_bottom"]->Fill(theta_x_minus, 1.);
               TH1F_right_histos["th_y_gen_right_nojet_bottom"]->Fill(theta_y_minus, 1.);
               if (rp_right_accep_bottom && fid_cuts_nom_right_bottom){ 
                  TH1F_right_histos["t_gen_rp_right_cut_nojet_bottom"]->Fill(fabs(t_gen_proton_right), 1./corr_xi);
                  TH1F_right_histos["xi_gen_rp_right_cut_nojet_bottom"]->Fill(xi_gen_proton_right, 1./corr_xi);
                  TH1F_right_histos["th_x_gen_rp_right_nojet_bottom"]->Fill(theta_x_minus, 1./corr_xi);
                  TH1F_right_histos["th_y_gen_rp_right_nojet_bottom"]->Fill(theta_y_minus, 1./corr_xi);
                } 
            }
        }   

        if (rp_right && jet_gen_sel){
             TH2F_right_histos["pt2_xi_gen_minus_rp"]->Fill( fabs(t_gen_proton_right)*(1 - xi_gen_proton_right), xi_gen_proton_right, 1);
        }
      
        if (jet_rec_sel){
             TH1F_right_histos["xi_cms_rec_right_sasha"]->Fill(xi_rec_cms_minus,1);
             TH1F_right_histos["xi_cms_rec_right"]->Fill(xi_rec_cms_minus, 1);
        }

        if (jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4){
              TH1F_right_histos["xi_cms_gen_right_sasha"]->Fill(xi_gen_cms_minus, event_weight_minus);
              TH1F_right_histos["xi_gen_right_sasha"]->Fill(xi_gen_proton_right, event_weight_minus);
        }

        if (jet_gen_sel){
              TH1F_right_histos["xi_cms_gen_right"]->Fill(xi_gen_cms_minus, event_weight_minus);
              TH2F_right_histos["pt2_xi_gen_minus"]->Fill( fabs(t_gen_proton_right)/*(1 + xi_gen_proton_right)*/, xi_gen_proton_right, event_weight_minus);
              if (rp_right) TH2F_right_histos["pt2_xi_gen_minus_rp"]->Fill( fabs(t_gen_proton_right)/**(1 + xi_gen_proton_right)*/, xi_gen_proton_right, event_weight_minus);
        }   

        TH2F_right_histos["log_t_vs_log_xi_right"]->Fill(log10(fabs(t_gen_proton_right)), log10(xi_gen_proton_right));

        if (jet_rec_sel && (rp_right_accep_top || rp_right_accep_bottom)){
            TH1F_right_histos["xi_rec_right_cut_sasha"]->Fill(xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
            TH1F_right_histos["xi_rec_cms_minus_sasha"]->Fill(1.25*xi_rec_cms_minus, event_weight_minus*eff_trigger*eff_proton/corr_xi);
        }    

        if (jet_rec_sel && proton_rec_sel && rp_right_sel){
              TH1F_right_histos["xi_cms_minus_totem_rec_right"]->Fill(xi_rec_cms_minus - xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
              TH1F_right_histos["xi_cms_minus_totem_rec_right_bin"]->Fill(xi_rec_cms_minus - xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
              TH1F_right_histos["xi_rec_right"]->Fill(xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
              if (xi_rec_cms_minus - xi_rec_proton_right<0){
                  TH1F_right_histos["xi_rec_right_cut"]->Fill(xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["t_rec_right_cut"]->Fill(fabs(t_rec_proton_right), event_weight_minus*eff_trigger*eff_proton/corr_xi); 
                  TH1F_right_histos["beta_right_cut"]->Fill(beta_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["log_x_rec_right_cut"]->Fill(log10(x_rec_right), event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["pt_jet1_right_cut"]->Fill(jet1_rec_pt, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["pt_jet2_right_cut"]->Fill(jet2_rec_pt, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["eta_jet1_right_cut"]->Fill(jet1_rec_eta, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["eta_jet2_right_cut"]->Fill(jet2_rec_eta, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["x_rec_gen_right"]->Fill(log10(x_rec_right)-log10(x_gen_right),1);
                  TH1F_right_histos["th_x_rec_right"]->Fill(theta_x_minus_smear, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["th_y_rec_right"]->Fill(theta_y_minus_smear, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["delta_eta_jets_right_cut"]->Fill(jet1_rec_eta - jet2_rec_eta, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["delta_phi_jets_right_cut"]->Fill(jet1_rec_phi - jet2_rec_phi, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["xi_cms_rec_right_cut_sasha"]->Fill(xi_rec_cms_minus, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["mass_jj_rec_right"]->Fill(sqrt(mjj2_rec_minus), event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["mass_x_rec_right"]->Fill(4000*sqrt(xi_rec_proton_right), event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["r_jj_rec_right"]->Fill(sqrt(mjj2_rec_minus)/(4000*sqrt(xi_rec_proton_right)), event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["reco_t_right"]->Fill(fabs(t_rec_proton_right), event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["reco_xi_right"]->Fill(xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                  TH1F_right_histos["reco_logx_right"]->Fill(log10(x_rec_right), event_weight_minus*eff_trigger*eff_proton/corr_xi);

                  // if (!(jet_gen_sel && proton_gen_sel)){
                  //    t_minus_response.Fake(fabs(t_rec_proton_right), event_weight_minus*norm_minus*eff_proton/corr_xi);
                  //    xi_minus_response.Fake(xi_rec_proton_right, event_weight_minus*norm_minus*eff_proton/corr_xi);
                  //    logx_minus_response.Fake(log10(x_rec_right), event_weight_minus*norm_minus*eff_proton/corr_xi);
                  //  }
                  // if (!(jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4 && proton_gen_sel))   
                  //     xi_cms_minus_response.Fake(xi_rec_cms_minus, event_weight_minus*norm_minus*eff_proton/corr_xi);
                  if (jet_gen_sel && proton_gen_sel){
                     // t_minus_response.Fill(fabs(t_rec_proton_right), fabs(t_gen_proton_right), event_weight_minus*norm_minus*eff_proton/corr_xi);
                     // xi_minus_response.Fill(xi_rec_proton_right, xi_gen_proton_right, event_weight_minus*norm_minus*eff_proton/corr_xi);
                     // logx_minus_response.Fill(log10(x_rec_right),log10(x_gen_right), event_weight_minus*norm_minus*eff_proton/corr_xi);
                     TH2F_right_histos["xi_minus_th2"]->Fill(xi_rec_proton_right, xi_gen_proton_right, event_weight_minus*norm_minus*eff_trigger*eff_proton/corr_xi);
                     TH2F_right_histos["logx_minus_th2"]->Fill(log10(x_rec_right),log10(x_gen_right), event_weight_minus*norm_minus*eff_trigger*eff_proton/corr_xi);
                     TH2F_right_histos["t_minus_th2"]->Fill(fabs(t_rec_proton_right), fabs(t_gen_proton_right), event_weight_minus*norm_minus*eff_trigger*eff_proton/corr_xi);
                     TH2F_right_histos["response_t_right"]->Fill(fabs(t_rec_proton_right), fabs(t_gen_proton_right), event_weight_minus*eff_proton*norm_minus*eff_trigger/corr_xi);
                     TH2F_right_histos["response_xi_right"]->Fill(xi_rec_proton_right, xi_gen_proton_right, event_weight_minus*eff_proton*norm_minus*eff_trigger/corr_xi);
                     TH2F_right_histos["response_logx_right"]->Fill(log10(x_rec_right), log10(x_gen_right), event_weight_minus*eff_proton*norm_minus*eff_trigger/corr_xi);
                   }  
// 	              //     if(jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4 && proton_gen_sel)
// 	              //        xi_cms_minus_response.Fill(xi_rec_cms_minus, xi_gen_cms_minus, event_weight_minus*norm_minus*eff_proton/corr_xi);
              }
        }

        if (jet_rec_sel && proton_rec_sel){
        	if (rp_right_accep_top && fid_cuts_nom_right_top) TH1F_right_histos["xi_cms_minus_totem_rec_right_top"]->Fill(xi_rec_cms_minus - xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
        	if (rp_right_accep_bottom && fid_cuts_nom_right_bottom) TH1F_right_histos["xi_cms_minus_totem_rec_right_bottom"]->Fill(xi_rec_cms_minus - xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
            if (xi_rec_cms_minus - xi_rec_proton_right<0){
            	if (rp_right_accep_top && fid_cuts_nom_right_top){
	                TH1F_right_histos["xi_rec_right_cut_top"]->Fill(xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                    TH1F_right_histos["t_rec_right_cut_top"]->Fill(fabs(t_rec_proton_right),event_weight_minus*eff_trigger*eff_proton/corr_xi); 
                    TH1F_right_histos["log_x_rec_right_cut_top"]->Fill(log10(x_rec_right), event_weight_minus*eff_trigger*eff_proton/corr_xi);
            	}
            	if (rp_right_accep_bottom && fid_cuts_nom_right_bottom){
	                    TH1F_right_histos["xi_rec_right_cut_bottom"]->Fill(xi_rec_proton_right, event_weight_minus*eff_trigger*eff_proton/corr_xi);
                    TH1F_right_histos["t_rec_right_cut_bottom"]->Fill(fabs(t_rec_proton_right), event_weight_minus*eff_trigger*eff_proton/corr_xi); 
                    TH1F_right_histos["log_x_rec_right_cut_bottom"]->Fill(log10(x_rec_right), event_weight_minus*eff_trigger*eff_proton/corr_xi);
            	}
            }
        }

        if (jet_gen_sel && proton_gen_sel){
             TH1F_right_histos["t_gen_right_cut"]->Fill(fabs(t_gen_proton_right), event_weight_minus); 
             TH1F_right_histos["xi_gen_right"]->Fill(xi_gen_proton_right, event_weight_minus);
             TH1F_right_histos["xi_gen_right_cut_sasha"]->Fill(xi_gen_proton_right, event_weight_minus);
             TH1F_right_histos["log_x_gen_right_cut"]->Fill(log10(x_gen_right), event_weight_minus);
             TH1F_right_histos["th_x_gen_right"]->Fill(theta_x_minus, event_weight_minus);
             TH1F_right_histos["th_y_gen_right"]->Fill(theta_y_minus, event_weight_minus);
             TH1F_right_histos["truth_t_right"]->Fill(fabs(t_gen_proton_right), event_weight_minus);
             TH1F_right_histos["truth_xi_right"]->Fill(xi_gen_proton_right, event_weight_minus);
             TH1F_right_histos["truth_logx_right"]->Fill(log10(x_gen_right), event_weight_minus);
// 	             // if (!(jet_rec_sel && proton_rec_sel && rp_right_sel && xi_rec_cms_minus - xi_rec_proton_right<0)){
// 	             //    t_minus_response.Miss(fabs(t_gen_proton_right), event_weight_minus*norm_minus);
// 	             //    xi_minus_response.Miss(xi_gen_proton_right, event_weight_minus*norm_minus);
// 	             //    logx_minus_response.Miss(log10(x_gen_right), event_weight_minus*norm_minus);
// 	             // }   
        } 

// 	        if(jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4 && proton_gen_sel){
// 	             TH1F_right_histos["xi_cms_gen_right_cut_sasha"]->Fill(xi_gen_cms_minus, event_weight_minus);
// 	             // if (!(jet_rec_sel && proton_rec_sel && rp_right_sel && xi_rec_cms_minus - xi_rec_proton_right<0)){
// 	             //    xi_cms_minus_response.Miss(xi_gen_cms_minus, event_weight_minus*norm_minus);
// 	             //  }  
// 	        }    
    }

cout<<regg<<endl;
    TFile* file_plus; double norm_plus = 0;
    if (mc == "pomwig_pom"){ file_plus = pomwig_pom_plus; norm_plus = norm_pomwig_plus_pom;}
    if (mc == "pomwig_regg"){ file_plus = pomwig_regg_plus; norm_plus = norm_pomwig_plus_regg;}
	if (mc == "pythia8_4c"){ file_plus = pythia8; norm_plus = norm_pythia_4c;}
	if (mc == "pythia8_cuetp8m1"){ file_plus = pythia8_CUETP8M1; norm_plus = norm_pythia_cuetp8m1;}

    TTree* tree_plus = (TTree*) file_plus->Get( treeName.c_str() );
    int nev_plus = int(tree_plus->GetEntriesFast());
    cout << mc <<" plus file has " << nev_plus << " entries : " << endl;

    double jet1_rec_pt_left, jet1_rec_eta_left, jet1_rec_phi_left, jet2_rec_pt_left, jet2_rec_eta_left, jet2_rec_phi_left;
    double jet1_gen_pt_left, jet1_gen_eta_left, jet1_gen_phi_left, jet2_gen_pt_left, jet2_gen_eta_left, jet2_gen_phi_left;
    double xi_rec_cms_plus, xi_rec_proton_left, x_rec_left, x_gen_left;
    double xi_gen_cms_plus, xi_gen_proton_left, xi_rec_proton_left_gauss, xi_rec_proton_left_theta;
    double t_rec_proton_left, t_rec_proton_left_theta, t_rec_proton_left_gauss, beta_rec_proton_left;
    double t_gen_proton_left, beta_gen_proton_left, rp_xpos_024, rp_xpos_025, rp_ypos_024, rp_ypos_025, t_rec_proton_left_old;
    bool rp_left, rp_left_accep_top, rp_left_accep_bottom;
    double px_plus, py_plus, pz_plus, e_plus, mass_plus, mjj2_rec_plus;
    int nVtx_left;
    tree_plus->SetBranchAddress("xi_rec_cms_left",&xi_rec_cms_plus);
    tree_plus->SetBranchAddress("xi_gen_cms_left",&xi_gen_cms_plus);
    tree_plus->SetBranchAddress("xi_rec_proton_left",&xi_rec_proton_left);
    tree_plus->SetBranchAddress("xi_gen_proton_left",&xi_gen_proton_left);
    tree_plus->SetBranchAddress("t_rec_proton_left",&t_rec_proton_left);
    tree_plus->SetBranchAddress("t_gen_proton_left",&t_gen_proton_left);
    tree_plus->SetBranchAddress("beta_rec_left",&beta_rec_proton_left);
    tree_plus->SetBranchAddress("beta_gen_left",&beta_gen_proton_left);
    tree_plus->SetBranchAddress("rp_left",&rp_left);
    tree_plus->SetBranchAddress("rp_left_accep_bottom",&rp_left_accep_bottom);
    tree_plus->SetBranchAddress("rp_left_accep_top",&rp_left_accep_top);
    tree_plus->SetBranchAddress("rp_xpos_24",&rp_xpos_024);
    tree_plus->SetBranchAddress("rp_ypos_24",&rp_ypos_024);
    tree_plus->SetBranchAddress("rp_xpos_25",&rp_xpos_025);
    tree_plus->SetBranchAddress("rp_ypos_25",&rp_ypos_025);
    tree_plus->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt_left);
    tree_plus->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt_left);
    tree_plus->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta_left);
    tree_plus->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta_left);
    tree_plus->SetBranchAddress("jet1_rec_phi",&jet1_rec_phi_left);
    tree_plus->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt_left);
    tree_plus->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt_left);
    tree_plus->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta_left);
    tree_plus->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta_left);
    tree_plus->SetBranchAddress("jet2_rec_phi",&jet2_rec_phi_left);
    tree_plus->SetBranchAddress("x_rec_left",&x_rec_left);
    tree_plus->SetBranchAddress("x_gen_left",&x_gen_left);
    tree_plus->SetBranchAddress("xi_rec_proton_left_gauss",&xi_rec_proton_left_gauss);
    tree_plus->SetBranchAddress("t_rec_proton_left_gauss",&t_rec_proton_left_gauss);
    tree_plus->SetBranchAddress("theta_x_plus",&theta_x_plus);
    tree_plus->SetBranchAddress("theta_y_plus",&theta_y_plus);
    tree_plus->SetBranchAddress("theta_x_plus_smear",&theta_x_plus_smear);
    tree_plus->SetBranchAddress("theta_y_plus_smear",&theta_y_plus_smear);
    tree_plus->SetBranchAddress("px_proton_left",&px_plus);
    tree_plus->SetBranchAddress("py_proton_left",&py_plus);
    tree_plus->SetBranchAddress("pz_proton_left",&pz_plus);
    tree_plus->SetBranchAddress("mass_proton_left",&mass_plus);
    tree_plus->SetBranchAddress("nVtx",&nVtx_left);
    tree_plus->SetBranchAddress("mjj2_rec",&mjj2_rec_plus);

    for(int i_evt = 0; i_evt < nev_plus; ++i_evt){
  	    tree_plus->GetEntry(i_evt);
       
        // if (single_vertex && nVtx_left!=1) continue;
        if (nVtx_left<1) continue;

        ///trigger efficiency
        double eff_trigger = func_trigger->Eval(jet2_rec_pt_left);

        xi_rec_proton_left = (unc_gauss == false) ? xi_rec_proton_left : xi_rec_proton_left_gauss;
        t_rec_proton_left = (unc_gauss == false) ? t_rec_proton_left : t_rec_proton_left_gauss; 

        double reweigth_beta_plus = (beta_gen_proton_left<=0.7) ? func->Eval(beta_gen_proton_left) : 1;
        double reweight_slope_plus = (fabs(t_gen_proton_left<=0.45)) ? scale_rew_slope_left*t_slope_data*exp(-t_slope_data*fabs(t_gen_proton_left))/(t_slope_mc*exp(-t_slope_mc*fabs(t_gen_proton_left))) : 1;//0.95 normalization factor
        double event_weight_plus;
        if (reweight && !reweight_slope) event_weight_plus = scale_rew_left*reweigth_beta_plus; 
        if (reweight && reweight_slope) event_weight_plus = scale_rew_left*reweigth_beta_plus*reweight_slope_plus;
        if (!reweight && reweight_slope) event_weight_plus = reweight_slope_plus;
        if (!reweight && !reweight_slope) event_weight_plus = 1;

        bool jet_rec_sel = jet1_rec_pt_left>pt_threshold && jet2_rec_pt_left>pt_threshold && fabs(jet1_rec_eta_left)<4.4 && fabs(jet2_rec_eta_left)<4.4;
        bool jet_gen_sel = jet1_gen_pt_left>pt_threshold && jet2_gen_pt_left>pt_threshold && fabs(jet1_gen_eta_left)<4.4 && fabs(jet2_gen_eta_left)<4.4;
        bool proton_rec_sel =  xi_rec_proton_left>0 && xi_rec_proton_left<0.1 && fabs(t_rec_proton_left)>0.03 && fabs(t_rec_proton_left)<1;
        bool proton_gen_sel =  xi_gen_proton_left>0 && xi_gen_proton_left<0.1 && fabs(t_gen_proton_left)>0.03 && fabs(t_gen_proton_left)<1;

        bool fid_cuts_nom_left_top = rp_xpos_024>0 && rp_xpos_024<0.007 && rp_ypos_024 >0.0084 && rp_ypos_024<0.027 ;
        bool fid_cuts_nom_left_bottom = rp_xpos_025>0 && rp_xpos_025<0.007 && rp_ypos_025 <-0.0084 && rp_ypos_025>-0.027 ;
        bool fid_cuts_unc_y_up_left_top = rp_xpos_024>0 && rp_xpos_024<0.007 && rp_ypos_024 >0.0084 && rp_ypos_024<0.0272 ;
        bool fid_cuts_unc_y_up_left_bottom = rp_xpos_025>0 && rp_xpos_025<0.007 && rp_ypos_025 <-0.0084 && rp_ypos_025>-0.0272 ;
        bool fid_cuts_unc_y_dw_left_top = rp_xpos_024>0 && rp_xpos_024<0.007 && rp_ypos_024 >0.0082 && rp_ypos_024<0.027 ;
        bool fid_cuts_unc_y_dw_left_bottom = rp_xpos_025>0 && rp_xpos_025<0.007 && rp_ypos_025 <-0.0082 && rp_ypos_025>-0.027 ;
        bool fid_cuts_unc_x_left_top = rp_xpos_024>0.00 && rp_xpos_024<0.006 && rp_ypos_024 >0.0084 && rp_ypos_024<0.027 ;
        bool fid_cuts_unc_x_left_bottom = rp_xpos_025>0.00 && rp_xpos_025<0.006 && rp_ypos_025 <-0.0084 && rp_ypos_025>-0.027 ;

      	bool rp_left_sel = false; 
      	if(unc_rp_y_up == false && unc_rp_y_dw == false && unc_rp_x == false) rp_left_sel = (rp_left_accep_top && fid_cuts_nom_left_top) || (rp_left_accep_bottom && fid_cuts_nom_left_bottom);
      	if(unc_rp_y_up == true && unc_rp_y_dw == false && unc_rp_x == false) rp_left_sel = (rp_left_accep_top && fid_cuts_unc_y_up_left_top) || (rp_left_accep_bottom && fid_cuts_unc_y_up_left_bottom);
      	if(unc_rp_y_up == false && unc_rp_y_dw == true && unc_rp_x == false) rp_left_sel = (rp_left_accep_top && fid_cuts_unc_y_dw_left_top) || (rp_left_accep_bottom && fid_cuts_unc_y_dw_left_bottom);
      	if(unc_rp_y_up == false && unc_rp_y_dw == false && unc_rp_x == true) rp_left_sel = (rp_left_accep_top && fid_cuts_unc_x_left_top) || (rp_left_accep_bottom && fid_cuts_unc_x_left_bottom);

        double corr_xi = 1;//(xi_rec_proton_right_pom > 0.08) ? 1.32864 : 1.;
        if (rec_corr) xi_thx_corr(xi_gen_proton_left, theta_x_plus, corr_xi);

      	if ((rp_left_accep_top && fid_cuts_nom_left_top) && !(rp_left_accep_bottom && fid_cuts_nom_left_bottom)) TH2F_left_histos["rp_pos_left"]->Fill(rp_xpos_024*1000, rp_ypos_024*1000, 1.);
        if (!(rp_left_accep_top && fid_cuts_nom_left_top) && (rp_left_accep_bottom && fid_cuts_nom_left_bottom)) TH2F_left_histos["rp_pos_left"]->Fill(rp_xpos_025*1000, rp_ypos_025*1000, 1.);

	    TH2F_left_histos["delta_thx_vs_delta_xi_left"]->Fill(theta_x_plus_smear-theta_x_plus, xi_rec_proton_left-xi_gen_proton_left, 1.);

        if(proton_gen_sel){
             TH1F_left_histos["t_gen_left_cut_nojet"]->Fill(fabs(t_gen_proton_left), 1.);
             TH1F_left_histos["xi_gen_left_cut_nojet"]->Fill(xi_gen_proton_left, 1.);
             TH2F_left_histos["t_gen_vs_xi_gen_plus"]->Fill( fabs(t_gen_proton_left), xi_gen_proton_left, 1.);
             TH1F_left_histos["t_rec_left_cut_nojet"]->Fill(fabs(t_rec_proton_left), 1.);
             TH1F_left_histos["xi_rec_left_cut_nojet"]->Fill(xi_rec_proton_left, 1.);
             TH1F_left_histos["th_x_gen_left_nojet"]->Fill(theta_x_plus, 1.);
             TH1F_left_histos["th_y_gen_left_nojet"]->Fill(theta_y_plus, 1.);
             TH2F_left_histos["theta_x_vs_xi_gen_nojet"]->Fill(xi_gen_proton_left, theta_x_plus, 1.);
             if(rp_left_sel){
               TH1F_left_histos["t_gen_rp_left_cut_nojet"]->Fill(fabs(t_gen_proton_left), 1.);
               TH1F_left_histos["xi_gen_rp_left_cut_nojet"]->Fill(xi_gen_proton_left, 1.);
               TH1F_left_histos["th_x_gen_rp_left_nojet"]->Fill(theta_x_plus, 1.);
               TH1F_left_histos["th_y_gen_rp_left_nojet"]->Fill(theta_y_plus, 1.);
               TH2F_left_histos["theta_x_vs_xi_gen_rp_nojet"]->Fill(xi_gen_proton_left, theta_x_plus, 1.);
               TH1F_left_histos["t_gen_rp_left_cut_nojet_aftercorr"]->Fill(fabs(t_gen_proton_left), 1./corr_xi);
               TH1F_left_histos["xi_gen_rp_left_cut_nojet_aftercorr"]->Fill(xi_gen_proton_left, 1./corr_xi);
               TH1F_left_histos["th_x_gen_rp_left_nojet_aftercorr"]->Fill(theta_x_plus, 1./corr_xi);
               TH1F_left_histos["th_y_gen_rp_left_nojet_aftercorr"]->Fill(theta_y_plus, 1./corr_xi);
               TH2F_left_histos["theta_x_vs_xi_gen_rp_nojet_aftercorr"]->Fill(xi_gen_proton_left, theta_x_plus, 1./corr_xi);
             }
            if(theta_y_plus >= 0){
               TH1F_left_histos["t_gen_left_cut_nojet_top"]->Fill(fabs(t_gen_proton_left), 1.);
               TH1F_left_histos["xi_gen_left_cut_nojet_top"]->Fill(xi_gen_proton_left, 1.);
               TH1F_left_histos["th_x_gen_left_nojet_top"]->Fill(theta_x_plus, 1.);
               TH1F_left_histos["th_y_gen_left_nojet_top"]->Fill(theta_y_plus, 1.);
               if (rp_left_accep_top && fid_cuts_nom_left_top){ 
                  TH1F_left_histos["t_gen_rp_left_cut_nojet_top"]->Fill(fabs(t_gen_proton_left), 1./corr_xi);
                  TH1F_left_histos["xi_gen_rp_left_cut_nojet_top"]->Fill(xi_gen_proton_left, 1./corr_xi);
                  TH1F_left_histos["th_x_gen_rp_left_nojet_top"]->Fill(theta_x_plus, 1./corr_xi);
                  TH1F_left_histos["th_y_gen_rp_left_nojet_top"]->Fill(theta_y_plus, 1./corr_xi);
                } 
            }     
            if(theta_y_plus < 0){
               TH1F_left_histos["t_gen_left_cut_nojet_bottom"]->Fill(fabs(t_gen_proton_left), 1.);
               TH1F_left_histos["xi_gen_left_cut_nojet_bottom"]->Fill(xi_gen_proton_left, 1.);
               TH1F_left_histos["th_x_gen_left_nojet_bottom"]->Fill(theta_x_plus, 1.);
               TH1F_left_histos["th_y_gen_left_nojet_bottom"]->Fill(theta_y_plus, 1.);
               if (rp_left_accep_bottom && fid_cuts_nom_left_bottom){ 
                  TH1F_left_histos["t_gen_rp_left_cut_nojet_bottom"]->Fill(fabs(t_gen_proton_left), 1./corr_xi);
                  TH1F_left_histos["xi_gen_rp_left_cut_nojet_bottom"]->Fill(xi_gen_proton_left, 1./corr_xi);
                  TH1F_left_histos["th_x_gen_rp_left_nojet_bottom"]->Fill(theta_x_plus, 1./corr_xi);
                  TH1F_left_histos["th_y_gen_rp_left_nojet_bottom"]->Fill(theta_y_plus, 1./corr_xi);
                } 
            }
        }   


    	if(proton_gen_sel){
        	TH1F_left_histos["t_gen_left_cut_nojet"]->Fill(fabs(t_gen_proton_left), 1.);
         	TH1F_left_histos["xi_gen_left_cut_nojet"]->Fill(xi_gen_proton_left, 1.);
         	TH2F_left_histos["t_gen_vs_xi_gen_plus"]->Fill( fabs(t_gen_proton_left), xi_gen_proton_left, 1.);
        }   

        if(proton_rec_sel && rp_left_sel){
        	TH1F_left_histos["t_rec_left_cut_nojet"]->Fill(fabs(t_rec_proton_left), 0.94);
         	TH1F_left_histos["xi_rec_left_cut_nojet"]->Fill(xi_rec_proton_left, 0.94);
        	TH2F_left_histos["t_rec_vs_xi_rec_plus"]->Fill( fabs(t_rec_proton_left), xi_rec_proton_left, 0.94);
        }   

	        if (rp_left && jet_gen_sel){
	         TH2F_left_histos["pt2_xi_gen_plus_rp"]->Fill( fabs(t_gen_proton_left)*(1 - xi_gen_proton_left), xi_gen_proton_left, 1);
        }
  
    	if (jet_rec_sel){
           TH1F_left_histos["xi_cms_rec_left_sasha"]->Fill(xi_rec_cms_plus,1);
           TH1F_left_histos["xi_cms_rec_left"]->Fill(xi_rec_cms_plus, 1);
       	}

        if (jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4){
          TH1F_left_histos["xi_cms_gen_left_sasha"]->Fill(xi_gen_cms_plus,1.);
        }

        if (jet_gen_sel){
          TH1F_left_histos["xi_cms_gen_left"]->Fill(xi_gen_cms_plus, event_weight_plus);
          TH2F_left_histos["pt2_xi_gen_plus"]->Fill( fabs(t_gen_proton_left)/*(1 + xi_gen_proton_left)*/, xi_gen_proton_left, event_weight_plus);
          if (rp_left) TH2F_left_histos["pt2_xi_gen_plus_rp"]->Fill( fabs(t_gen_proton_left)/**(1 + xi_gen_proton_left)*/, xi_gen_proton_left, event_weight_plus);
        }   

        TH2F_left_histos["log_t_vs_log_xi_left"]->Fill(log10(fabs(t_gen_proton_left)), log10(xi_gen_proton_left));

        if (jet_rec_sel && proton_rec_sel && rp_left_sel){
          TH1F_left_histos["xi_cms_minus_totem_rec_left"]->Fill(xi_rec_cms_plus - xi_rec_proton_left, event_weight_plus);
          TH1F_left_histos["xi_cms_minus_totem_rec_left_bin"]->Fill(xi_rec_cms_plus - xi_rec_proton_left, event_weight_plus);
          TH1F_left_histos["xi_rec_left"]->Fill(xi_rec_proton_left, event_weight_plus);
          if (xi_rec_cms_plus - xi_rec_proton_left<0){
              TH1F_left_histos["xi_rec_left_cut"]->Fill(xi_rec_proton_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["xi_rec_left_cut_sasha"]->Fill(xi_rec_proton_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["t_rec_left_cut"]->Fill(fabs(t_rec_proton_left), event_weight_plus*eff_trigger*eff_proton/corr_xi); 
              TH1F_left_histos["beta_left_cut"]->Fill(beta_rec_proton_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["log_x_rec_left_cut"]->Fill(log10(x_rec_left), event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["pt_jet1_left_cut"]->Fill(jet1_rec_pt_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["pt_jet2_left_cut"]->Fill(jet2_rec_pt_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["eta_jet1_left_cut"]->Fill(jet1_rec_eta_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["eta_jet2_left_cut"]->Fill(jet2_rec_eta_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["x_rec_gen_left"]->Fill(log10(x_rec_left)-log10(x_gen_left),1);
              TH1F_left_histos["th_x_rec_left"]->Fill(theta_x_plus_smear, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["th_y_rec_left"]->Fill(theta_y_plus_smear, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["delta_eta_jets_left_cut"]->Fill(jet1_rec_eta_left - jet2_rec_eta_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["delta_phi_jets_left_cut"]->Fill(jet1_rec_phi_left - jet2_rec_phi_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["xi_cms_rec_left_cut_sasha"]->Fill(xi_rec_cms_plus, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["mass_jj_rec_left"]->Fill(sqrt(mjj2_rec_plus), event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["mass_x_rec_left"]->Fill(4000*sqrt(xi_rec_proton_left), event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["r_jj_rec_left"]->Fill(sqrt(mjj2_rec_plus)/(4000*sqrt(xi_rec_proton_left)), event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["reco_t_left"]->Fill(fabs(t_rec_proton_left), event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["reco_xi_left"]->Fill(xi_rec_proton_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
              TH1F_left_histos["reco_logx_left"]->Fill(log10(x_rec_left), event_weight_plus*eff_trigger*eff_proton/corr_xi);

        //       // if (!(jet_gen_sel && proton_gen_sel)){
        //       //    t_plus_response.Fake(fabs(t_rec_proton_left), event_weight_plus*norm_plus*eff_proton/corr_xi);
        //       //    xi_plus_response.Fake(xi_rec_proton_left, event_weight_plus*norm_plus*eff_proton/corr_xi);
        //       //    logx_plus_response.Fake(log10(x_rec_left), event_weight_plus*norm_plus*eff_proton/corr_xi);
        //       //  }
        //       // if (!(jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4 && proton_gen_sel))   
        //       //     xi_cms_plus_response.Fake(xi_rec_cms_plus, event_weight_plus*norm_plus*eff_proton/corr_xi);
              if (jet_gen_sel && proton_gen_sel){
        //          // t_plus_response.Fill(fabs(t_rec_proton_left), fabs(t_gen_proton_left), event_weight_plus*norm_plus*eff_proton/corr_xi);
        //          // xi_plus_response.Fill(xi_rec_proton_left, xi_gen_proton_left, event_weight_plus*norm_plus*eff_proton/corr_xi);
        //          // logx_plus_response.Fill(log10(x_rec_left),log10(x_gen_left), event_weight_plus*norm_plus*eff_proton/corr_xi);
                 TH2F_left_histos["xi_plus_th2"]->Fill(xi_rec_proton_left, xi_gen_proton_left, event_weight_plus*norm_plus*eff_proton*eff_trigger/corr_xi);
                 TH2F_left_histos["logx_plus_th2"]->Fill(log10(x_rec_left),log10(x_gen_left), event_weight_plus*norm_plus*eff_proton*eff_trigger/corr_xi);
                 TH2F_left_histos["t_plus_th2"]->Fill(fabs(t_rec_proton_left), fabs(t_gen_proton_left), event_weight_plus*norm_plus*eff_proton*eff_trigger/corr_xi);
                 TH2F_left_histos["response_t_left"]->Fill(fabs(t_rec_proton_left), fabs(t_gen_proton_left), event_weight_plus*eff_proton*norm_plus*eff_trigger/corr_xi);
                 TH2F_left_histos["response_xi_left"]->Fill(xi_rec_proton_left, xi_gen_proton_left, event_weight_plus*eff_proton*norm_plus*eff_trigger/corr_xi);
                 TH2F_left_histos["response_logx_left"]->Fill(log10(x_rec_left), log10(x_gen_left), event_weight_plus*eff_proton*norm_plus*eff_trigger/corr_xi);
                }  
        //         // if(jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4 && proton_gen_sel)
        //         //  xi_cms_plus_response.Fill(xi_rec_cms_plus, xi_gen_cms_plus, event_weight_plus*norm_plus*eff_proton/corr_xi);
            }
        }     
      
        if (jet_rec_sel && proton_rec_sel){
        	if (rp_left_accep_top && fid_cuts_nom_left_top) TH1F_left_histos["xi_cms_minus_totem_rec_left_top"]->Fill(xi_rec_cms_plus - xi_rec_proton_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
        	if (rp_left_accep_bottom && fid_cuts_nom_left_bottom) TH1F_left_histos["xi_cms_minus_totem_rec_left_bottom"]->Fill(xi_rec_cms_plus - xi_rec_proton_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
            if (xi_rec_cms_plus - xi_rec_proton_left<0){
            	if (rp_left_accep_top && fid_cuts_nom_left_top){
	                TH1F_left_histos["xi_rec_left_cut_top"]->Fill(xi_rec_proton_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
                    TH1F_left_histos["t_rec_left_cut_top"]->Fill(fabs(t_rec_proton_left), event_weight_plus*eff_trigger*eff_proton/corr_xi); 
                    TH1F_left_histos["log_x_rec_left_cut_top"]->Fill(log10(x_rec_left), event_weight_plus*eff_trigger*eff_proton/corr_xi);
            	}
            	if (rp_left_accep_bottom && fid_cuts_nom_left_bottom){
	                TH1F_left_histos["xi_rec_left_cut_bottom"]->Fill(xi_rec_proton_left, event_weight_plus*eff_trigger*eff_proton/corr_xi);
                    TH1F_left_histos["t_rec_left_cut_bottom"]->Fill(fabs(t_rec_proton_left), event_weight_plus*eff_trigger*eff_proton/corr_xi); 
                    TH1F_left_histos["log_x_rec_left_cut_bottom"]->Fill(log10(x_rec_left), event_weight_plus*eff_trigger*eff_proton/corr_xi);
            	}
            }
        }

        if (jet_gen_sel && proton_gen_sel){
            TH1F_left_histos["t_gen_left_cut"]->Fill(fabs(t_gen_proton_left), event_weight_plus); 
         	TH1F_left_histos["xi_gen_left"]->Fill(xi_gen_proton_left, event_weight_plus);
         	TH1F_left_histos["xi_gen_left_cut_sasha"]->Fill(xi_gen_proton_left, event_weight_plus);
         	TH1F_left_histos["log_x_gen_left_cut"]->Fill(log10(x_gen_left), event_weight_plus);
         	TH1F_left_histos["th_x_gen_left"]->Fill(theta_x_plus, event_weight_plus);
         	TH1F_left_histos["th_y_gen_left"]->Fill(theta_y_plus, event_weight_plus);
            TH1F_left_histos["truth_t_left"]->Fill(fabs(t_gen_proton_left), event_weight_plus);
            TH1F_left_histos["truth_xi_left"]->Fill(xi_gen_proton_left, event_weight_plus);
            TH1F_left_histos["truth_logx_left"]->Fill(log10(x_gen_left), event_weight_plus);
         	// if (!(jet_rec_sel && proton_rec_sel && rp_left_sel && xi_rec_cms_plus - xi_rec_proton_left<0)){
          //   	t_plus_response.Miss(fabs(t_gen_proton_left), event_weight_plus*norm_plus);
          //   	xi_plus_response.Miss(xi_gen_proton_left, event_weight_plus*norm_plus);
          //   	logx_plus_response.Miss(log10(x_gen_left), event_weight_plus*norm_plus);
         	// }   
        } 

        // if(jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4 && proton_gen_sel){
        //    TH1F_left_histos["xi_cms_gen_left_cut_sasha"]->Fill(xi_gen_cms_plus, event_weight_plus);
        //    // if (!(jet_rec_sel && proton_rec_sel && rp_left_sel && xi_rec_cms_plus - xi_rec_proton_left<0)){
        //    //    xi_cms_plus_response.Miss(xi_gen_cms_plus, event_weight_plus*norm_plus);
        //    //  }  
        // }    
    } 

    for(map<string,TH1F*>::const_iterator histos = TH1F_right_histos.begin(); histos != TH1F_right_histos.end(); ++histos){
       histos->second->Scale(norm_minus);
       histos->second->Write();
    }   
    for(map<string,TH1F*>::const_iterator histos = TH1F_left_histos.begin(); histos != TH1F_left_histos.end(); ++histos){
       histos->second->Scale(norm_plus);
       histos->second->Write();
    }
    for(map<string,TH2F*>::const_iterator histos = TH2F_right_histos.begin(); histos != TH2F_right_histos.end(); ++histos){
       histos->second->Write();
    }   
    for(map<string,TH2F*>::const_iterator histos = TH2F_left_histos.begin(); histos != TH2F_left_histos.end(); ++histos){
       histos->second->Write();
    }

    outfile->Close();   
	return;
}


void MCAnalysis::addHistosPomwigPomRegg (void){

	TFile* histos_pomwig_pom = TFile::Open("histos_pomwig_pom_betarew_sloperew_nominal.root","READ");
	TFile* histos_pomwig_regg = TFile::Open("histos_pomwig_regg_betarew_sloperew_nominal.root","READ");
    TFile* outfile = new TFile("pomwig_pom_plus_regg_betarew_sloperew_nominal.root", "RECREATE");


	if (!histos_pomwig_pom || !histos_pomwig_regg) { cout << "Pomwig files do not exist" << endl; return; }

	TIter next_pom (histos_pomwig_pom->GetListOfKeys());
	TKey* key_pom;
	while ((key_pom = (TKey*)next_pom())){
 		TClass *cl_pom = gROOT->GetClass(key_pom->GetClassName());
    	if (!cl_pom->InheritsFrom("TH1")) continue;
      	std::string histname_pom = key_pom->GetName();
      	TH1F* histo_pomeron = (TH1F*)key_pom->ReadObj();

  		TIter next_regg (histos_pomwig_regg->GetListOfKeys());
		TKey* key_regg;
		while ((key_regg = (TKey*)next_regg())){
			TClass *cl_regg = gROOT->GetClass(key_regg->GetClassName());
      		if (!cl_regg->InheritsFrom("TH1")) continue;
      		std::string histname_regg = key_regg->GetName();
      		TH1F* histo_reggeon = (TH1F*)key_regg->ReadObj();
      		 if (histname_pom == histname_regg){ histo_pomeron->Add(histo_reggeon); histo_pomeron->Write();}
		}	
	}
	outfile->Close();    
}

void MCAnalysis::getUnfoldedHistos (TH1F* reco, TH1F* truth, TH2F* response, TH1F* data, TH1F* &histo_unfolded, int n_iter){

	//Unfolding data
	RooUnfoldResponse response_matrix (reco, truth, response, "unfold", "unfold");
	RooUnfoldBayes bayes (&response_matrix, data, n_iter);
	histo_unfolded = (TH1F*) bayes.Hreco();
}

void MCAnalysis::getUnfoldingTest (TH1F* reco, TH1F* truth, TH2F* response, TH1F* data, TH2F* &chi2_unfolded_vs_iter, TH2F* &delta_chi2_unfolded_vs_iter,
	double &chi2_smeared){

	//response matrix setup
	RooUnfoldResponse response_matrix (reco, truth, response, "unfold", "unfold");

	//Bottom line test
	double chi2_unfolded[40];
    chi2_unfolded_vs_iter = new TH2F("chi2_unfolded_vs_iter","",100,0,40,100,0.0001,10);
    delta_chi2_unfolded_vs_iter = new TH2F("delta_chi2_unfolded_vs_iter","",2000,0,40,2000,0.0001,1);

	for (int i = 1; i < 40; ++i){
		RooUnfoldBayes bayes_per_iter (&response_matrix, data, i);
		TH1F* histo_unfolded_per_iter = (TH1F*) bayes_per_iter.Hreco();

	    chi2_unfolded[i] = 0; 
        for (int j = 1; j <= histo_unfolded_per_iter->GetNbinsX()-1; ++j){
       		 if(histo_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded[i] += pow(histo_unfolded_per_iter->GetBinContent(j) - truth->GetBinContent(j), 2)/(pow(histo_unfolded_per_iter->GetBinError(j),2));
    	} 
	    chi2_unfolded_vs_iter->Fill(i, log10(chi2_unfolded[i])); 
	}

    for (int j = 1; j <= data->GetNbinsX()-1; ++j){
        if (data->GetBinError(j)>0) chi2_smeared += pow(data->GetBinContent(j) - reco->GetBinContent(j), 2)/(pow(data->GetBinError(j),2));
    }	

    //chi2 variation
    for (int i = 2; i<40; ++i){
    	if (chi2_unfolded[i]!=0) delta_chi2_unfolded_vs_iter->Fill(i, fabs(chi2_unfolded[i]-chi2_unfolded[i-1])/chi2_unfolded[i]); 
	}


}

// int main(void)
// {
// 	MCAnalysis pomwig, regg;
//  //    pomwig.getMCHistos("pomwig_regg", true,true,false,false,false,false,false, false);
//  //    // pomwig.getMCHistos("pomwig_pom", true,true,false,false,false,false, false, true);
//  // //    pomwig.getMCHistos("pomwig_regg", false,false,false,false,true,false);
// 	// pomwig.getMCHistos("pythia8_4c", true,true,false,false,false,false,false,false);
// 	// pomwig.getMCHistos("pythia8_cuetp8m1", true,true,false,false,false,false,false,false);
//     // pomwig.getMCHistos("pomwig_regg", true,false,false,false,false,false);
//     // pomwig.getMCHistos("pomwig_regg", false,true,false,false,false,false);
//     // pomwig.getMCHistos("pomwig_regg", true,true,false,false,false,false);
//     // pomwig.getMCHistos("pomwig_regg", true,true,true,false,false,false);
//     // pomwig.getMCHistos("pomwig_regg", true,true,false,true,false,false);
//     // pomwig.getMCHistos("pomwig_regg", true,true,false,false,true,false);
//     // pomwig.getMCHistos("pomwig_regg", true,true,false,false,false,true);
// 	// rec->Draw();
// 	pomwig.addHistosPomwigPomRegg();
// 	// pomwig.getUnfoldedHistos(rec, gen, response);
// 	return 0;
// }