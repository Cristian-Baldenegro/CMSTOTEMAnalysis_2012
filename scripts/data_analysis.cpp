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
#include <TAxis.h>
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
#include <TFractionFitter.h>
#include "Math/MinimizerOptions.h"
#include "TMinuit.h"


//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
using namespace std;

float tbins[9] = {0.03, 0.07, 0.11, 0.21, 0.31, 0.45,  0.58, 0.72, 1.};
float xi_bins[12] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2};
float bin_sasha[9] ={0.0003, 0.002, 0.0045, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1};
float bin[16] = {-0.4, -0.112, -0.096, -0.08, -0.064, -0.048, -0.032, -0.016, 0, 0.048, 0.112, 0.176, 0.24, 0.304, 0.368, 0.4};


class BHRdmDataEvent {
   public:
      BHRdmDataEvent() {}
      ~BHRdmDataEvent() {}

      double xi_totem_right_rdm;
      double xi_totem_left_rdm;
      double t_totem_right_rdm;
      double t_totem_left_rdm;
      bool rp_right_rdm;
      bool rp_left_rdm;
      bool valid_vtx_rdm;
      bool valid_proton_right_rdm;
      bool valid_proton_left_rdm;
      double xi_cms_right;
      double xi_cms_left;
      double x_right;
      double x_left;
      double xi_totem_xicmscut_right;
      double t_totem_xicmscut_right;
      double xi_totem_xicmscut_left;
      double t_totem_xicmscut_left;
}; 

class DataAnalysis
{
public:
	DataAnalysis();
	void getHistos( bool unc_jec_up = false, bool unc_jec_dw = false, bool unc_pf_up = false, bool unc_pf_dw = false, bool unc_rp_y_up = false, bool unc_rp_y_dw = false,
			bool unc_rp_x = false, bool unc_xi_up = false, bool unc_xi_dw = false);
    void AssymError (TH1F* histo_unfolded_nominal, TH1F* histo_jes_up_unfolded, TH1F* histo_jes_dw_unfolded, TH1F* histo_pf_up_unfolded, 
    TH1F* histo_pf_dw_unfolded, TH1F* histo_xi_up_unfolded, TH1F* histo_xi_dw_unfolded, TH1F* histo_norew_unfolded, TH1F* histo_betarew_unfolded, 
    TH1F* histo_rp_up_unfolded, TH1F* histo_rp_dw_unfolded, TH1F* histo_rp_x_unfolded, TH1F* histo_gauss_unfolded, TH1F* histo_hera_backg_unfolded,
    TH1F* histo_iter_unfolded, TH1F* histo_right_unfolded, TH1F* histo_left_unfolded, TH1F* histo_unfolded_mc1, TH1F* histo_unfolded_m2, 
    TH1F* histo_unfolded_mc3, TH1F* histo_mc1_unfolded_mc2, TH1F* histo_mc1_unfolded_mc3, TH1F* histo_mc2_unfolded_mc1, TH1F* histo_mc2_unfolded_mc3,
    TH1F* histo_mc3_unfolded_mc1, TH1F* histo_mc3_unfolded_mc2, TH1F* histo_gen_mc1, TH1F* histo_gen_mc2, TH1F* histo_gen_mc3, 
    TH1F* histo_trigger_up, TH1F* histo_trigger_dw, TH1F* y_ref, TGraphAsymmErrors* &error,
    string const& result = "t", bool histos = false);
    void AbsoluteSigma (TH1F* data, TH1F* mc_gen, TH1F* mc_nom_rec, double &sigma_data, double &error_data);
    void AssymErrorInclusive (TH1F* histo_unfolded_nominal, TH1F* histo_jes_up_unfolded, TH1F* histo_jes_dw_unfolded, TH1F* histo_pf_up_unfolded, 
    TH1F* histo_pf_dw_unfolded, TH1F* histo_iter_unfolded, TH1F* histo_right_unfolded, TH1F* histo_left_unfolded, TH1F* histo_unfolded_mc1, TH1F* histo_unfolded_mc2, 
    TH1F* histo_unfolded_mc3, TH1F* histo_mc1_unfolded_mc2, TH1F* histo_mc1_unfolded_mc3, TH1F* histo_mc2_unfolded_mc1, TH1F* histo_mc2_unfolded_mc3,
    TH1F* histo_mc3_unfolded_mc1, TH1F* histo_mc3_unfolded_mc2, TH1F* histo_gen_mc1, TH1F* histo_gen_mc2, TH1F* histo_gen_mc3, TH1F* y_ref, TGraphAsymmErrors* &error);
	~DataAnalysis();

private:
	TFile * outfile;
};

DataAnalysis::DataAnalysis (void){
}

void DataAnalysis::getHistos (bool unc_jec_up, bool unc_jec_dw, bool unc_pf_up, bool unc_pf_dw, bool unc_rp_y_up, bool unc_rp_y_dw, bool unc_rp_x, bool unc_xi_up, 
					bool unc_xi_dw) {

    // TFile* data_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_jecwinter.root","READ");
	TFile* data_file;
    if (unc_jec_up) data_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_jet_up_noveto.root","READ");
    else if (unc_jec_dw) data_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_jet_dw_noveto.root","READ");
    else if (unc_pf_up) data_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_pf_up_noveto.root","READ");
    else if (unc_pf_dw) data_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_pf_dw_noveto.root","READ");
    else data_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_noveto.root","READ");

    string treeName = "small_tree"; 
    string fileName;

    if (unc_jec_up && !unc_jec_dw && !unc_pf_up && !unc_pf_dw && !unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_xi_up && !unc_xi_dw) 
    	fileName = "data_unc_jec_up.root";
    else if (!unc_jec_up && unc_jec_dw && !unc_pf_up && !unc_pf_dw && !unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_xi_up && !unc_xi_dw) 
    	fileName = "data_unc_jec_dw.root";
    else if (!unc_jec_up && !unc_jec_dw && unc_pf_up && !unc_pf_dw && !unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_xi_up && !unc_xi_dw) 
    	fileName = "data_unc_pf_up.root";
    else if (!unc_jec_up && !unc_jec_dw && !unc_pf_up && unc_pf_dw && !unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_xi_up && !unc_xi_dw) 
    	fileName = "data_unc_pf_dw.root";
    else if (!unc_jec_up && !unc_jec_dw && !unc_pf_up && !unc_pf_dw && !unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_xi_up && !unc_xi_dw) 
    	fileName = "data_nominal.root";
    else if (!unc_jec_up && !unc_jec_dw && !unc_pf_up && !unc_pf_dw && unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_xi_up && !unc_xi_dw) 
    	fileName = "data_unc_rp_y_up.root";
    else if (!unc_jec_up && !unc_jec_dw && !unc_pf_up && !unc_pf_dw && !unc_rp_y_up && unc_rp_y_dw && !unc_rp_x && !unc_xi_up && !unc_xi_dw) 
        fileName = "data_unc_rp_y_dw.root";
    else if (!unc_jec_up && !unc_jec_dw && !unc_pf_up && !unc_pf_dw && !unc_rp_y_up && !unc_rp_y_dw && unc_rp_x && !unc_xi_up && !unc_xi_dw) 
    	fileName = "data_unc_rp_x.root";
    else if (!unc_jec_up && !unc_jec_dw && !unc_pf_up && !unc_pf_dw && !unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && unc_xi_up && !unc_xi_dw) 
    	fileName = "data_unc_xi_up.root";
    else if (!unc_jec_up && !unc_jec_dw && !unc_pf_up && !unc_pf_dw && !unc_rp_y_up && !unc_rp_y_dw && !unc_rp_x && !unc_xi_up && unc_xi_dw) 
    	fileName = "data_unc_xi_dw.root";
    else { cout << "Wrong initial values, set them again..." << endl; return;}

    cout << "output: " << fileName << endl;

    outfile = new TFile(fileName.c_str(), "RECREATE");

    map<string,TH1F*> TH1F_histos;
    map<string,TH2F*> TH2F_histos;
    TH1F_histos["xi_cms_minus_totem_right"] = new TH1F("xi_cms_minus_totem_right","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_right_top"] = new TH1F("xi_cms_minus_totem_right_top","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_right_bottom"] = new TH1F("xi_cms_minus_totem_right_bottom","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_right_bin"] = new TH1F("xi_cms_minus_totem_right_bin","", 15, bin);
    TH1F_histos["xi_cms_minus_totem_left"] = new TH1F("xi_cms_minus_totem_left","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_left_top"] = new TH1F("xi_cms_minus_totem_left_top","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_left_bottom"] = new TH1F("xi_cms_minus_totem_left_bottom","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_left_bin"] = new TH1F("xi_cms_minus_totem_left_bin","", 15, bin);
    TH1F_histos["xi_cms_minus_totem_right_bh"] = new TH1F("xi_cms_minus_totem_right_bh","",50,-0.4,0.4);
    TH1F_histos["xi_cms_minus_totem_left_bh"] = new TH1F("xi_cms_minus_totem_left_bh","",50,-0.4,0.4);
    TH1F_histos["xi_right"] = new TH1F("xi_right","",50,-0.04,0.2);
    TH1F_histos["xi_left"] = new TH1F("xi_left","",50,-0.04,0.2);
    TH1F_histos["xi_right_cut"] = new TH1F("xi_right_cut","",11,xi_bins);
    TH1F_histos["xi_right_cut_top"] = new TH1F("xi_right_cut_top","",11,xi_bins);
    TH1F_histos["xi_right_cut_bottom"] = new TH1F("xi_right_cut_bottom","",11,xi_bins);
    TH1F_histos["xi_right_bh"] = new TH1F("xi_right_bh","",50,-0.04,0.2);
    TH1F_histos["xi_right_bh_cut"] = new TH1F("xi_right_bh_cut","",11,xi_bins);
    TH1F_histos["xi_right_cut_noeff"] = new TH1F("xi_right_cut_noeff","",11,xi_bins);
    TH1F_histos["xi_left_cut"] = new TH1F("xi_left_cut","",11,xi_bins);
    TH1F_histos["xi_left_cut_top"] = new TH1F("xi_left_cut_top","",11,xi_bins);
    TH1F_histos["xi_left_cut_bottom"] = new TH1F("xi_left_cut_bottom","",11,xi_bins);
    TH1F_histos["xi_left_bh"] = new TH1F("xi_left_bh","",50,-0.04,0.2);
    TH1F_histos["xi_left_bh_cut"] = new TH1F("xi_left_bh_cut","",11,xi_bins);
    TH1F_histos["xi_left_cut_noeff"] = new TH1F("xi_left_cut_noeff","",11,xi_bins);
    TH1F_histos["t_right_cut"] = new TH1F("t_right_cut","", 8, tbins);
    TH1F_histos["t_right_cut_top"] = new TH1F("t_right_cut_top","", 8, tbins);
    TH1F_histos["t_right_cut_bottom"] = new TH1F("t_right_cut_bottom","", 8, tbins);
    TH1F_histos["t_right_bh_cut"] = new TH1F("t_right_bh_cut","", 8, tbins);
    TH1F_histos["t_right_cut_zb"] = new TH1F("t_right_cut_zb","", 8, tbins);
    TH1F_histos["t_right_cut_full"] = new TH1F("t_right_cut_full","", 8, tbins);
    TH1F_histos["t_right_cut_bh"] = new TH1F("t_right_cut_bh","", 8, tbins);
    TH1F_histos["t_right_cut_noeff"] = new TH1F("t_right_cut_noeff","", 8, tbins);
    TH1F_histos["t_left_cut"] = new TH1F("t_left_cut","", 8, tbins);
    TH1F_histos["t_left_cut_top"] = new TH1F("t_left_cut_top","", 8, tbins);
    TH1F_histos["t_left_cut_bottom"] = new TH1F("t_left_cut_bottom","", 8, tbins);
    TH1F_histos["t_left_bh_cut"] = new TH1F("t_left_bh_cut","", 8, tbins);
    TH1F_histos["t_left_cut_zb"] = new TH1F("t_left_cut_zb","", 8, tbins);
    TH1F_histos["t_left_cut_full"] = new TH1F("t_left_cut_full","", 8, tbins);
    TH1F_histos["t_left_cut_bh"] = new TH1F("t_left_cut_bh","", 8, tbins);
    TH1F_histos["t_left_cut_noeff"] = new TH1F("t_left_cut_noeff","", 8, tbins);
    TH1F_histos["beta_right_cut"] = new TH1F("beta_right_cut","",15,0,1);
    TH1F_histos["beta_left_cut"] = new TH1F("beta_left_cut","",15,0,1);
    TH1F_histos["pt_jet1"] = new TH1F("pt_jet1","", 15, 0, 200);
    TH1F_histos["pt_jet1_right_cut"] = new TH1F("pt_jet1_right_cut","", 15, 0, 200);
    TH1F_histos["pt_jet1_left_cut"] = new TH1F("pt_jet1_left_cut","", 15, 0, 200);
    TH1F_histos["pt_jet2"] = new TH1F("pt_jet2","", 15, 0, 200);
    TH1F_histos["pt_jet2_right_cut"] = new TH1F("pt_jet2_right_cut","", 15, 0, 200);
    TH1F_histos["pt_jet2_left_cut"] = new TH1F("pt_jet2_left_cut","", 15, 0, 200);
    TH1F_histos["eta_jet1"] = new TH1F("eta_jet1","", 20, -5.2, 5.2);
    TH1F_histos["eta_jet1_right_cut"] = new TH1F("eta_jet1_right_cut","", 20, -5.2, 5.2);
    TH1F_histos["eta_jet1_left_cut"] = new TH1F("eta_jet1_left_cut","", 20, -5.2, 5.2);
    TH1F_histos["eta_jet2"] = new TH1F("eta_jet2","", 20, -5.2, 5.2);
    TH1F_histos["eta_jet2_right_cut"] = new TH1F("eta_jet2_right_cut","", 20, -5.2, 5.2);
    TH1F_histos["eta_jet2_left_cut"] = new TH1F("eta_jet2_left_cut","", 20, -5.2, 5.2);
    TH1F_histos["delta_eta_jets_right_cut"] = new TH1F("delta_eta_jets_right_cut","", 40, -5.2, 5.2);
    TH1F_histos["delta_phi_jets_right_cut"] = new TH1F("delta_phi_jets_right_cut","", 40, -5.2, 5.2);
    TH1F_histos["delta_eta_jets_left_cut"] = new TH1F("delta_eta_jets_left_cut","", 40, -5.2, 5.2);
    TH1F_histos["delta_phi_jets_left_cut"] = new TH1F("delta_phi_jets_left_cut","", 40, -5.2, 5.2);
    TH1F_histos["log_x_right"] = new TH1F("log_x_right","",15, -4, 0);
    TH1F_histos["log_x_right_noeff"] = new TH1F("log_x_right_noeff","",15, -4, 0);
    TH1F_histos["log_x_right_cut"] = new TH1F("log_x_right_cut","",15, -4, 0);
    TH1F_histos["log_x_right_cut_top"] = new TH1F("log_x_right_cut_top","",15, -4, 0);
    TH1F_histos["log_x_right_cut_bottom"] = new TH1F("log_x_right_cut_bottom","",15, -4, 0);
    TH1F_histos["log_x_right_bh_cut"] = new TH1F("log_x_right_bh_cut","",15, -4, 0);
    TH1F_histos["log_x_right_cut_noeff"] = new TH1F("log_x_right_cut_noeff","",15, -4, 0);
    TH1F_histos["log_x_left"] = new TH1F("log_x_left","",15, -4, 0);
    TH1F_histos["log_x_left_noeff"] = new TH1F("log_x_left_noeff","",15, -4, 0);
    TH1F_histos["log_x_left_cut"] = new TH1F("log_x_left_cut","",15, -4, 0);
    TH1F_histos["log_x_left_cut_top"] = new TH1F("log_x_left_cut_top","",15, -4, 0);
    TH1F_histos["log_x_left_cut_bottom"] = new TH1F("log_x_left_cut_bottom","",15, -4, 0);
    TH1F_histos["log_x_left_bh_cut"] = new TH1F("log_x_left_bh_cut","",15, -4, 0);
    TH1F_histos["log_x_left_cut_noeff"] = new TH1F("log_x_left_cut_noeff","",15, -4, 0);
    TH1F_histos["xi_right_cut_sasha"] = new TH1F("xi_right_cut_sasha","",8, bin_sasha);
    TH1F_histos["xi_cms_minus_sasha"] = new TH1F("xi_cms_minus_sasha","",8, bin_sasha);
    TH1F_histos["xi_left_cut_sasha"] = new TH1F("xi_left_cut_sasha","",8, bin_sasha);
    TH1F_histos["xi_cms_right_cut_sasha"] = new TH1F("xi_cms_right_sasha_data","",8, bin_sasha);
    TH1F_histos["xi_cms_left_cut_sasha"] = new TH1F("xi_cms_left_sasha_data","",8, bin_sasha);
    TH1F_histos["sigma_xi_cms_right_sasha"] = new TH1F("sigma_xi_cms_right_sasha","",8, bin_sasha);
    TH1F_histos["sigma_xi_cms_left_sasha"] = new TH1F("sigma_xi_cms_left_sasha","",8, bin_sasha);
    TH1F_histos["vtx_z"] = new TH1F("vtx_z","",150,-30, 30);
    TH1F_histos["th_x_right"] = new TH1F("th_x_right","",20, -0.4e-3, 0.4e-3);
    TH1F_histos["th_y_right"] = new TH1F("th_y_right","",20, -0.4e-3, 0.4e-3);
    TH1F_histos["th_x_left"] = new TH1F("th_x_left","",20, -0.4e-3, 0.4e-3);
    TH1F_histos["th_y_left"] = new TH1F("th_y_left","",20, -0.4e-3, 0.4e-3);
    TH1F_histos["mass_x_right"] = new TH1F("mass_x_right","",20, 0, 1800);
    TH1F_histos["mass_x_left"] = new TH1F("mass_x_left","",20, 0, 1800);
    TH1F_histos["mass_jj_right"] = new TH1F("mass_jj_right","",20, 0, 1000);
    TH1F_histos["mass_jj_left"] = new TH1F("mass_jj_left","",20, 0, 1000);
    TH1F_histos["r_jj_right"] = new TH1F("r_jj_right","",20, 0, 1);
    TH1F_histos["r_jj_left"] = new TH1F("r_jj_left","",20, 0, 1);
    TH1F_histos["x_pos_right_top"] = new TH1F("x_pos_right_top", "", 20, -1, 10);
    TH1F_histos["y_pos_right_top"] = new TH1F("y_pos_right_top", "", 20, 0, 35);
    TH1F_histos["x_pos_right_bottom"] = new TH1F("x_pos_right_bottom", "", 20, -1, 10);
    TH1F_histos["y_pos_right_bottom"] = new TH1F("y_pos_right_bottom", "", 20, -35, 0);
    TH1F_histos["x_pos_left_top"] = new TH1F("x_pos_left_top", "", 20, -1, 10);
    TH1F_histos["y_pos_left_top"] = new TH1F("y_pos_left_top", "", 20, 0, 35);
    TH1F_histos["x_pos_left_bottom"] = new TH1F("x_pos_left_bottom", "", 20, -1, 10);
    TH1F_histos["y_pos_left_bottom"] = new TH1F("y_pos_left_bottom", "", 20, -35, 0);
    TH1F_histos["x_pos_right_cut_top"] = new TH1F("x_pos_right_cut_top", "", 20, -1, 10);
    TH1F_histos["y_pos_right_cut_top"] = new TH1F("y_pos_right_cut_top", "", 20, 0, 35);
    TH1F_histos["x_pos_right_cut_bottom"] = new TH1F("x_pos_right_cut_bottom", "", 20, -1, 10);
    TH1F_histos["y_pos_right_cut_bottom"] = new TH1F("y_pos_right_cut_bottom", "", 20, -35, 0);
    TH1F_histos["x_pos_left_cut_top"] = new TH1F("x_pos_left_cut_top", "", 20, -1, 10);
    TH1F_histos["y_pos_left_cut_top"] = new TH1F("y_pos_left_cut_top", "", 20, 0, 35);
    TH1F_histos["x_pos_left_cut_bottom"] = new TH1F("x_pos_left_cut_bottom", "", 20, -1, 10);
    TH1F_histos["y_pos_left_cut_bottom"] = new TH1F("y_pos_left_cut_bottom", "", 20, -35, 0);
    TH2F_histos["proton_y_vs_x_rp_024_025_accept"] = new TH2F("proton_y_vs_x_rp_024_025_accept","proton_y_vs_x_rp_024_025",100,-10,10,100,-40,40);
    TH2F_histos["proton_y_vs_x_rp_120_121_accept"] = new TH2F("proton_y_vs_x_rp_120_121_accept","proton_y_vs_x_rp_120_121",100,-10,10,100,-40,40);
    TH2F_histos["proton_y_vs_x_rp_124_125_accept"] = new TH2F("proton_y_vs_x_rp_124_125_accept","proton_y_vs_x_rp_124_125",100,-10,10,100,-40,40);

    for(map<string,TH1F*>::const_iterator histos = TH1F_histos.begin(); histos != TH1F_histos.end(); ++histos){
       histos->second->Sumw2();
    }

	TTree* tree_data;
    tree_data = (TTree*) data_file->Get( treeName.c_str() );
    int nev_data = int(tree_data->GetEntriesFast());
    cout <<"The data file has " << nev_data << " entries : " << endl;

    double jet1_pt, jet1_eta, jet1_phi, jet2_pt, jet2_eta, jet2_phi, mjj2;
    double xi_cms_minus_data, xi_proton_right_data, xi_cms_plus_data, xi_proton_left_data, x_right_data, x_left_data;
    double t_proton_right_data, t_proton_left_data, beta_proton_right_data, beta_proton_left_data, eff_trigger;
    bool valid_proton_right_data, valid_proton_left_data;
    bool rp_right_top_data, rp_right_bottom_data, rp_left_top_data, rp_left_bottom_data; 
    bool primVtx = false;
    int nVtx;
    double x_pos_024, x_pos_025, x_pos_124, x_pos_125, y_pos_024, y_pos_025, y_pos_124, y_pos_125, vtx_x, vtx_y, vtx_z_pos, x_pos_120, x_pos_121, y_pos_120, y_pos_121, lumisec;
    double thx_proton_right, thx_proton_left, thy_proton_right, thy_proton_left;
    tree_data->SetBranchAddress("eff_trigger",&eff_trigger);
    tree_data->SetBranchAddress("xi_cms_minus",&xi_cms_minus_data);
    tree_data->SetBranchAddress("xi_cms_plus",&xi_cms_plus_data);
    tree_data->SetBranchAddress("xi_totem_right",&xi_proton_right_data);
    tree_data->SetBranchAddress("xi_totem_left",&xi_proton_left_data);
    tree_data->SetBranchAddress("t_totem_right",&t_proton_right_data);
    tree_data->SetBranchAddress("t_totem_left",&t_proton_left_data);
    tree_data->SetBranchAddress("beta_proton_right",&beta_proton_right_data);
    tree_data->SetBranchAddress("beta_proton_left",&beta_proton_left_data);
    tree_data->SetBranchAddress("x_right",&x_right_data);
    tree_data->SetBranchAddress("x_left",&x_left_data);
    tree_data->SetBranchAddress("valid_proton_right",&valid_proton_right_data);
    tree_data->SetBranchAddress("valid_proton_left",&valid_proton_left_data);
    tree_data->SetBranchAddress("rp_right_top",&rp_right_top_data);
    tree_data->SetBranchAddress("rp_right_bottom",&rp_right_bottom_data);
    tree_data->SetBranchAddress("rp_left_top",&rp_left_top_data);
    tree_data->SetBranchAddress("rp_left_bottom",&rp_left_bottom_data);
    tree_data->SetBranchAddress("jet1_pt",&jet1_pt);
    tree_data->SetBranchAddress("jet1_eta",&jet1_eta);
    tree_data->SetBranchAddress("jet1_phi",&jet1_phi);
    tree_data->SetBranchAddress("jet2_pt",&jet2_pt);
    tree_data->SetBranchAddress("jet2_eta",&jet2_eta);
    tree_data->SetBranchAddress("jet2_phi",&jet2_phi);
    tree_data->SetBranchAddress("x_pos_024",&x_pos_024);
    tree_data->SetBranchAddress("y_pos_024",&y_pos_024);
    tree_data->SetBranchAddress("x_pos_025",&x_pos_025);
    tree_data->SetBranchAddress("y_pos_025",&y_pos_025);
    tree_data->SetBranchAddress("x_pos_120",&x_pos_120);
    tree_data->SetBranchAddress("y_pos_120",&y_pos_120);
    tree_data->SetBranchAddress("x_pos_121",&x_pos_121);
    tree_data->SetBranchAddress("y_pos_121",&y_pos_121);
    tree_data->SetBranchAddress("x_pos_124",&x_pos_124);
    tree_data->SetBranchAddress("y_pos_124",&y_pos_124);
    tree_data->SetBranchAddress("x_pos_125",&x_pos_125);
    tree_data->SetBranchAddress("y_pos_125",&y_pos_125);
    tree_data->SetBranchAddress("vtx_x",&vtx_x);
    tree_data->SetBranchAddress("vtx_y",&vtx_y);
    tree_data->SetBranchAddress("vtx_z",&vtx_z_pos);
    tree_data->SetBranchAddress("select_Vertex",&primVtx);
    tree_data->SetBranchAddress("nVtx",&nVtx);
    tree_data->SetBranchAddress("thetax_proton_right",&thx_proton_right);
    tree_data->SetBranchAddress("thetay_proton_right",&thy_proton_right);
    tree_data->SetBranchAddress("thetax_proton_left",&thx_proton_left);
    tree_data->SetBranchAddress("thetay_proton_left",&thy_proton_left);
    tree_data->SetBranchAddress("mjj2",&mjj2);

    ///-------- Beam halo estimate--------
    TRandom* rdm = new TRandom;
    rdm->SetSeed(12345);
    std::vector<BHRdmDataEvent> beamhalo;
    for(int i_evt = 0; i_evt < nev_data; ++i_evt){
        tree_data->GetEntry(i_evt);
        bool rp_right_data = rp_right_bottom_data || rp_right_top_data;
        bool rp_left_data = rp_left_bottom_data || rp_left_top_data;

        BHRdmDataEvent bh;
        bh.xi_cms_right = xi_cms_minus_data;
        bh.xi_cms_left = xi_cms_plus_data;
        if (xi_cms_minus_data>0.12 && valid_proton_right_data && rp_right_data && xi_proton_right_data>0 && xi_proton_right_data<0.1 && fabs(t_proton_right_data)>0.03 && fabs(t_proton_right_data)<1){
            bh.xi_totem_xicmscut_right = xi_proton_right_data;
            bh.t_totem_xicmscut_right = t_proton_right_data;
            bh.x_right = x_right_data;
        }    
        if (xi_cms_plus_data>0.12 && valid_proton_left_data && rp_left_data && xi_proton_left_data>0 && xi_proton_left_data<0.1 && t_proton_left_data>0.03 && t_proton_left_data<1){
            bh.xi_totem_xicmscut_left = xi_proton_left_data;
            bh.t_totem_xicmscut_left = t_proton_left_data;
            bh.x_left = x_left_data;
        }    
        beamhalo.push_back(bh);

    }
    ///-------- Beam halo estimate --------


    ///-------- Analysis of data ----------- 
    double eff_proton = 0.94;
    // double eff_vtx;
    // if (single_vertex) eff_vtx = 0.95;
    // else eff_vtx = 1;
    double pt_threshold = 40.;
    int nevents_hera_right=0;
    int nevents_right=0;
    int nevents_hera_left=0;
    int nevents_left=0;
    int nevents_jets = 0;
    int nevents_rp_right = 0;
    int nevents_rp_left = 0;
    int nevents_kin_proton_right = 0;
    int nevents_kin_proton_left = 0;
    int nevents_single_arm = 0;

    for(int i_evt = 0; i_evt < nev_data; ++i_evt){

        tree_data->GetEntry(i_evt);
        eff_trigger = 1.;

        // if (single_vertex && nVtx!=1) continue;
        // if (!single_vertex && nVtx<1) continue;
        if (nVtx<1) continue;

        if (unc_xi_up) {
           xi_proton_right_data = xi_proton_right_data + xi_proton_right_data*0.1;
           xi_proton_left_data = xi_proton_left_data + xi_proton_left_data*0.1;
        }   
        if (unc_xi_dw){
           xi_proton_right_data = xi_proton_right_data - xi_proton_right_data*0.1;
           xi_proton_left_data = xi_proton_left_data - xi_proton_left_data*0.1;
        }
        
        bool jet_sel = jet1_pt>pt_threshold && jet2_pt>pt_threshold && fabs(jet1_eta)<4.4 && fabs(jet2_eta)<4.4;
        bool proton_right_kin_sel =  xi_proton_right_data>0 && xi_proton_right_data<0.1 && t_proton_right_data>0.03 && t_proton_right_data<1;
        bool proton_left_kin_sel = xi_proton_left_data>0 && xi_proton_left_data<0.1 && t_proton_left_data>0.03 && t_proton_left_data<1;

      	bool fid_cut_right_top = x_pos_124>0 && x_pos_124<7 && y_pos_124 >8.4 && y_pos_124<27;
	    bool fid_cut_right_bottom = x_pos_125>0 && x_pos_125<7 && y_pos_125<-8.4 && y_pos_125>-27;
        bool fid_cut_unc_y_up_right_top = x_pos_124>0 && x_pos_124<7 && y_pos_124 >8.4 && y_pos_124<27.2 ;
        bool fid_cut_unc_y_up_right_bottom = x_pos_125>0 && x_pos_125<7 && y_pos_125 <-8.4 && y_pos_125>-27.2 ;
        bool fid_cut_unc_y_dw_right_top = x_pos_124>0 && x_pos_124<7 && y_pos_124 >8.2 && y_pos_124<27;
        bool fid_cut_unc_y_dw_right_bottom = x_pos_125>0 && x_pos_125<7 && y_pos_125<-8.2 && y_pos_125>-27;
        bool fid_cut_unc_x_right_top = x_pos_124>0 && x_pos_124<6 && y_pos_124 >8.4 && y_pos_124<27;
        bool fid_cut_unc_x_right_bottom = x_pos_125>0 && x_pos_125<6 && y_pos_125<-8.4 && y_pos_125>-27;

        bool fid_cut_left_top = x_pos_024>0 && x_pos_024<7 && y_pos_024 >8.4 && y_pos_024<27;
        bool fid_cut_left_bottom = x_pos_025>0 && x_pos_025<7 && y_pos_025<-8.4 && y_pos_025>-27;
        bool fid_cut_unc_y_up_left_top = x_pos_024>0 && x_pos_024<7 && y_pos_024 >8.4 && y_pos_024<27.2 ;
        bool fid_cut_unc_y_up_left_bottom = x_pos_025>0 && x_pos_025<7 && y_pos_025 <-8.4 && y_pos_025>-27.2 ;
        bool fid_cut_unc_y_dw_left_top = x_pos_024>0 && x_pos_024<7 && y_pos_024 >8.2 && y_pos_024<27 ;
        bool fid_cut_unc_y_dw_left_bottom = x_pos_025>0 && x_pos_025<7 && y_pos_025 <-8.2 && y_pos_025>-27 ;
        bool fid_cut_unc_x_left_top = x_pos_024>0 && x_pos_024<6 && y_pos_024 >8.4 && y_pos_024<27;
        bool fid_cut_unc_x_left_bottom = x_pos_025>0 && x_pos_025<6 && y_pos_025<-8.4 && y_pos_025>-27;

        bool rp_left_data = false;
        if (unc_rp_y_up == false && unc_rp_y_dw == false && unc_rp_x == false) rp_left_data = (rp_left_top_data && fid_cut_left_top) || (rp_left_bottom_data && fid_cut_left_bottom);
        if (unc_rp_y_up == true && unc_rp_y_dw == false && unc_rp_x == false) rp_left_data = (rp_left_top_data && fid_cut_unc_y_up_left_top) || (rp_left_bottom_data && fid_cut_unc_y_up_left_bottom);
        if (unc_rp_y_up == false && unc_rp_y_dw == true && unc_rp_x == false) rp_left_data = (rp_left_top_data && fid_cut_unc_y_dw_left_top) || (rp_left_bottom_data && fid_cut_unc_y_dw_left_bottom);
        if (unc_rp_y_up == false && unc_rp_y_dw == false && unc_rp_x == true) rp_left_data = (rp_left_top_data && fid_cut_unc_x_left_top) || (rp_left_bottom_data && fid_cut_unc_x_left_bottom);

        bool rp_right_data = false;
        if (unc_rp_y_up == false && unc_rp_y_dw == false && unc_rp_x == false) rp_right_data = (rp_right_top_data && fid_cut_right_top) || (rp_right_bottom_data && fid_cut_right_bottom);
        if (unc_rp_y_up == true && unc_rp_y_dw == false && unc_rp_x == false) rp_right_data = (rp_right_top_data && fid_cut_unc_y_up_right_top) || (rp_right_bottom_data && fid_cut_unc_y_up_right_bottom);
        if (unc_rp_y_up == false && unc_rp_y_dw == true && unc_rp_x == false) rp_right_data = (rp_right_top_data && fid_cut_unc_y_dw_right_top) || (rp_right_bottom_data && fid_cut_unc_y_dw_right_bottom);
        if (unc_rp_y_up == false && unc_rp_y_dw == false && unc_rp_x == true) rp_right_data = (rp_right_top_data && fid_cut_unc_x_right_top) || (rp_right_bottom_data && fid_cut_unc_x_right_bottom);

        int i_evt_backg_1 = 0 + rdm->Rndm()*(beamhalo.size());
        int i_evt_backg_2 = 0 + rdm->Rndm()*(beamhalo.size());
        BHRdmDataEvent const & beamhalo_data_1 = beamhalo.at(i_evt_backg_1);
        BHRdmDataEvent const & beamhalo_data_2 = beamhalo.at(i_evt_backg_2);
        double xi_right_cms_bh = beamhalo_data_1.xi_cms_right;
        double xi_right_totem_bh = beamhalo_data_2.xi_totem_xicmscut_right;
        double x_right_bh = beamhalo_data_2.x_right;
        double x_left_bh = beamhalo_data_2.x_left;
        double t_right_totem_bh = beamhalo_data_2.t_totem_xicmscut_right;
        double xi_left_cms_bh = beamhalo_data_1.xi_cms_left;
        double xi_left_totem_bh = beamhalo_data_2.xi_totem_xicmscut_left;
        double t_left_totem_bh = beamhalo_data_2.t_totem_xicmscut_left;
        TH1F_histos["vtx_z"]->Fill(vtx_z_pos, 1);
	
        // TF1* func_trigger = new TF1("func_trigger", fFermiLike, 0., 200, 2);
        // func_trigger->SetParameter(0,5);
        // func_trigger->SetParameter(1,0.5);
        // cout<< "eff_trigger: "<<func_trigger->Eval(40)<<endl;

        if (jet_sel){++nevents_jets;
            TH1F_histos["pt_jet1"]->Fill(jet1_pt, 1/(eff_trigger));
            TH1F_histos["pt_jet2"]->Fill(jet2_pt, 1/(eff_trigger));
            TH1F_histos["eta_jet1"]->Fill(jet1_eta, 1/(eff_trigger));
            TH1F_histos["eta_jet2"]->Fill(jet2_eta, 1/(eff_trigger));
            TH1F_histos["log_x_right"]->Fill(log10(x_right_data),1/(eff_trigger));
            TH1F_histos["log_x_right_noeff"]->Fill(log10(x_right_data),1);
            TH1F_histos["log_x_left"]->Fill(log10(x_left_data),1/(eff_trigger));
            TH1F_histos["log_x_left_noeff"]->Fill(log10(x_left_data),1);
            TH1F_histos["sigma_xi_cms_right_sasha"]->Fill(xi_cms_minus_data, 1/(eff_trigger));
            TH1F_histos["sigma_xi_cms_left_sasha"]->Fill(xi_cms_plus_data, 1/(eff_trigger));
        }

        // if (valid_proton_right_data && valid_proton_left_data) continue; // double arm
        if (!(valid_proton_right_data || valid_proton_left_data)) continue;
        if (jet_sel) ++nevents_single_arm;
        if (jet_sel && rp_right_data) ++nevents_rp_right;
        if (jet_sel && rp_left_data) ++nevents_rp_left;

        // if (jet_sel && valid_proton_right_data && rp_right_data && xi_right_totem_bh<0.1 &&  t_right_totem_bh>0.03 && t_right_totem_bh<1 && xi_right_cms_bh - xi_right_totem_bh<0)t_right_bh_cut->Fill(t_right_totem_bh,0.8/(eff_trigger)); 
        if (jet_sel && valid_proton_right_data && (rp_right_top_data || rp_right_bottom_data)){
           TH1F_histos["xi_right_cut_sasha"]->Fill(xi_proton_right_data, 1/(eff_trigger));
           TH1F_histos["xi_cms_minus_sasha"]->Fill(xi_cms_minus_data, 1/(eff_trigger));
        }   
        if (jet_sel && proton_right_kin_sel && valid_proton_right_data && rp_right_data){ ++nevents_kin_proton_right;
            TH1F_histos["xi_cms_minus_totem_right"]->Fill(xi_cms_minus_data - xi_proton_right_data, 1.);
            TH1F_histos["xi_cms_minus_totem_right_bin"]->Fill(xi_cms_minus_data - xi_proton_right_data, 1.);
            TH1F_histos["xi_cms_minus_totem_right_bh"]->Fill(xi_right_cms_bh - xi_right_totem_bh, 1.); 
            TH1F_histos["xi_right"]->Fill(xi_proton_right_data, 1./(eff_trigger));
            TH1F_histos["xi_right_bh"]->Fill(xi_right_totem_bh, 1.); 
            if (xi_right_totem_bh>0 && xi_right_totem_bh<0.1 &&  t_right_totem_bh>0.03 && t_right_totem_bh<1 && xi_right_cms_bh - xi_right_totem_bh<0){
               TH1F_histos["t_right_cut_bh"]->Fill(t_right_totem_bh, 1);
                    TH1F_histos["t_right_bh_cut"]->Fill(t_right_totem_bh,1); 
                    ++nevents_hera_right;
                    TH1F_histos["xi_right_bh_cut"]->Fill(xi_right_totem_bh, 1.); 
                    TH1F_histos["log_x_right_bh_cut"]->Fill(log10(x_right_bh), 1.); 
            }   
            if (xi_cms_minus_data - xi_proton_right_data<0){++nevents_right;
                TH1F_histos["xi_right_cut"]->Fill(xi_proton_right_data, 1/(eff_trigger));
                TH1F_histos["xi_right_cut_noeff"]->Fill(xi_proton_right_data, 1);
                TH1F_histos["t_right_cut"]->Fill(t_proton_right_data, 1/(eff_trigger));
                TH1F_histos["t_right_cut_zb"]->Fill(t_proton_right_data, 1/(eff_trigger));
                TH1F_histos["t_right_cut_full"]->Fill(t_proton_right_data, 1/(eff_trigger));
                TH1F_histos["t_right_cut_noeff"]->Fill(t_proton_right_data, 1);
                TH1F_histos["beta_right_cut"]->Fill(beta_proton_right_data, 1/(eff_trigger));
                TH1F_histos["log_x_right_cut"]->Fill(log10(x_right_data),1./(eff_trigger));
                TH1F_histos["log_x_right_cut_noeff"]->Fill(log10(x_right_data),1.);
                TH1F_histos["pt_jet1_right_cut"]->Fill(jet1_pt, 1/(eff_trigger));
                TH1F_histos["eta_jet1_right_cut"]->Fill(jet1_eta,1./(eff_trigger));
                TH1F_histos["pt_jet2_right_cut"]->Fill(jet2_pt,1./(eff_trigger));
                TH1F_histos["eta_jet2_right_cut"]->Fill(jet2_eta,1./(eff_trigger));
                TH1F_histos["th_x_right"]->Fill(thx_proton_right, 1./(eff_trigger));
                TH1F_histos["th_y_right"]->Fill(thy_proton_right, 1./(eff_trigger));
                TH1F_histos["delta_eta_jets_right_cut"]->Fill(jet1_eta-jet2_eta, 1./(eff_trigger)); 
                TH1F_histos["delta_phi_jets_right_cut"]->Fill(jet1_phi-jet2_phi, 1./(eff_trigger));
                TH1F_histos["xi_cms_right_cut_sasha"]->Fill(xi_cms_minus_data, 1/(eff_trigger));
                TH1F_histos["mass_jj_right"]->Fill(sqrt(mjj2), 1/(eff_trigger));
                TH1F_histos["mass_x_right"]->Fill(4000*sqrt(xi_proton_right_data), 1/(eff_trigger));
                TH1F_histos["r_jj_right"]->Fill(sqrt(mjj2)/(4000*sqrt(xi_proton_right_data)), 1/(eff_trigger));

                // if(xi_right_cms_bh - xi_right_totem_bh<0){
                //     TH1F_histos["t_right_bh_cut"]->Fill(t_right_totem_bh,1); 
                //     ++nevents_hera_right;
                //     TH1F_histos["xi_right_bh_cut"]->Fill(xi_right_totem_bh, 1.); 
                //     TH1F_histos["log_x_right_bh_cut"]->Fill(log10(x_right_bh), 1.); 
                // }
            }
        }
        if (jet_sel && proton_right_kin_sel && valid_proton_right_data){
            if (rp_right_top_data && fid_cut_right_top){
                TH1F_histos["xi_cms_minus_totem_right_top"]->Fill(xi_cms_minus_data - xi_proton_right_data, 1.);
                TH1F_histos["x_pos_right_top"]->Fill(x_pos_124, 1./(eff_trigger));
                TH1F_histos["y_pos_right_top"]->Fill(y_pos_124, 1./(eff_trigger));
            }    
            if (rp_right_bottom_data && fid_cut_right_bottom){
                TH1F_histos["xi_cms_minus_totem_right_bottom"]->Fill(xi_cms_minus_data - xi_proton_right_data, 1.);
                TH1F_histos["x_pos_right_bottom"]->Fill(x_pos_125, 1./(eff_trigger));
                TH1F_histos["y_pos_right_bottom"]->Fill(y_pos_125, 1./(eff_trigger));
            }
            if (xi_cms_minus_data - xi_proton_right_data<0){
                if (rp_right_top_data && fid_cut_right_top){
                    TH1F_histos["xi_right_cut_top"]->Fill(xi_proton_right_data, 1/(eff_trigger));
                    TH1F_histos["t_right_cut_top"]->Fill(t_proton_right_data, 1/(eff_trigger));
                    TH1F_histos["log_x_right_cut_top"]->Fill(log10(x_right_data),1./(eff_trigger));
                    TH1F_histos["x_pos_right_cut_top"]->Fill(x_pos_124, 1./(eff_trigger));
                    TH1F_histos["y_pos_right_cut_top"]->Fill(y_pos_124, 1./(eff_trigger));
                }   
                if (rp_right_bottom_data && fid_cut_right_bottom){
                    TH1F_histos["xi_right_cut_bottom"]->Fill(xi_proton_right_data, 1/(eff_trigger));
                    TH1F_histos["t_right_cut_bottom"]->Fill(t_proton_right_data, 1/(eff_trigger));
                    TH1F_histos["log_x_right_cut_bottom"]->Fill(log10(x_right_data),1./(eff_trigger));
                    TH1F_histos["x_pos_right_cut_bottom"]->Fill(x_pos_125, 1./(eff_trigger));
                    TH1F_histos["y_pos_right_cut_bottom"]->Fill(y_pos_125, 1./(eff_trigger));
                }   
            }    
        }    

        if (jet_sel && proton_left_kin_sel && valid_proton_left_data && rp_left_data){++nevents_kin_proton_left;
            TH1F_histos["xi_cms_minus_totem_left"]->Fill(xi_cms_plus_data - xi_proton_left_data, 1.); 
            TH1F_histos["xi_cms_minus_totem_left_bin"]->Fill(xi_cms_plus_data - xi_proton_left_data, 1.);
            TH1F_histos["xi_cms_minus_totem_left_bh"]->Fill(xi_left_cms_bh - xi_left_totem_bh, 1.); 
            TH1F_histos["xi_left"]->Fill(xi_proton_left_data, 1./(eff_trigger));
            TH1F_histos["xi_left_bh"]->Fill(xi_left_totem_bh, 1.); 
            if (xi_left_totem_bh>0 && xi_left_totem_bh<0.1 &&  t_left_totem_bh>0.03 && t_left_totem_bh<1 && xi_left_cms_bh - xi_left_totem_bh<0){
               TH1F_histos["t_left_cut_bh"]->Fill(t_left_totem_bh, 1);
                  ++nevents_hera_left;
                  TH1F_histos["xi_left_bh_cut"]->Fill(xi_left_totem_bh, 1.); 
                  TH1F_histos["t_left_bh_cut"]->Fill(t_left_totem_bh, 1.); 
                  TH1F_histos["log_x_left_bh_cut"]->Fill(log10(x_left_bh), 1.); 
            }   
            if (xi_cms_plus_data - xi_proton_left_data<0){++nevents_left;
                TH1F_histos["xi_left_cut"]->Fill(xi_proton_left_data, 1/(eff_trigger));
                TH1F_histos["xi_left_cut_noeff"]->Fill(xi_proton_left_data, 1);
                TH1F_histos["xi_left_cut_sasha"]->Fill(xi_proton_left_data, 1/(eff_trigger));
                TH1F_histos["t_left_cut"]->Fill(t_proton_left_data, 1/(eff_trigger));
                TH1F_histos["t_left_cut_zb"]->Fill(t_proton_left_data, 1/(eff_trigger));
                TH1F_histos["t_left_cut_full"]->Fill(t_proton_left_data, 1/(eff_trigger));
                TH1F_histos["t_left_cut_noeff"]->Fill(t_proton_left_data, 1);
                TH1F_histos["beta_left_cut"]->Fill(beta_proton_left_data, 1/(eff_trigger));
                TH1F_histos["log_x_left_cut"]->Fill(log10(x_left_data),1./(eff_trigger));
                TH1F_histos["log_x_left_cut_noeff"]->Fill(log10(x_left_data),1.);
                TH1F_histos["pt_jet1_left_cut"]->Fill(jet1_pt,1./(eff_trigger));
                TH1F_histos["eta_jet1_left_cut"]->Fill(jet1_eta,1./(eff_trigger));
                TH1F_histos["pt_jet2_left_cut"]->Fill(jet2_pt,1./(eff_trigger));
                TH1F_histos["eta_jet2_left_cut"]->Fill(jet2_eta,1./(eff_trigger));
                TH1F_histos["th_x_left"]->Fill(thx_proton_left, 1./(eff_trigger)); 
                TH1F_histos["th_y_left"]->Fill(thy_proton_left, 1./(eff_trigger)); 
                TH1F_histos["delta_eta_jets_left_cut"]->Fill(jet1_eta-jet2_eta, 1./(eff_trigger)); 
                TH1F_histos["delta_phi_jets_left_cut"]->Fill(jet1_phi-jet2_phi, 1./(eff_trigger));
                TH1F_histos["xi_cms_left_cut_sasha"]->Fill(xi_cms_plus_data, 1/(eff_trigger));
                TH1F_histos["mass_jj_left"]->Fill(sqrt(mjj2), 1/(eff_trigger));
                TH1F_histos["mass_x_left"]->Fill(4000*sqrt(xi_proton_left_data), 1/(eff_trigger));
                TH1F_histos["r_jj_left"]->Fill(sqrt(mjj2)/(4000*sqrt(xi_proton_left_data)), 1/(eff_trigger));

                // if(xi_left_cms_bh - xi_left_totem_bh<0){
                //   ++nevents_hera_left;
                //   TH1F_histos["xi_left_bh_cut"]->Fill(xi_left_totem_bh, 1.); 
                //   TH1F_histos["t_left_bh_cut"]->Fill(t_left_totem_bh, 1.); 
                //   TH1F_histos["log_x_left_bh_cut"]->Fill(log10(x_left_bh), 1.); 
                // }   
            }
        }


        if (jet_sel && proton_left_kin_sel && valid_proton_left_data){
            if (rp_left_top_data && fid_cut_left_top){
                TH1F_histos["xi_cms_minus_totem_left_top"]->Fill(xi_cms_plus_data - xi_proton_left_data, 1.);
                TH1F_histos["x_pos_left_top"]->Fill(x_pos_024, 1./(eff_trigger));
                TH1F_histos["y_pos_left_top"]->Fill(y_pos_024, 1./(eff_trigger));
            }    
            if (rp_left_bottom_data && fid_cut_left_bottom){
                TH1F_histos["xi_cms_minus_totem_left_bottom"]->Fill(xi_cms_plus_data - xi_proton_left_data, 1.);
                TH1F_histos["x_pos_left_bottom"]->Fill(x_pos_025, 1./(eff_trigger));
                TH1F_histos["y_pos_left_bottom"]->Fill(y_pos_025, 1./(eff_trigger));
            }
            if (xi_cms_plus_data - xi_proton_left_data<0){
                if (rp_left_top_data && fid_cut_left_top){
                    TH1F_histos["xi_left_cut_top"]->Fill(xi_proton_left_data, 1/(eff_trigger));
                    TH1F_histos["t_left_cut_top"]->Fill(t_proton_left_data, 1/(eff_trigger));
                    TH1F_histos["log_x_left_cut_top"]->Fill(log10(x_left_data),1./(eff_trigger));
                    TH1F_histos["x_pos_left_cut_top"]->Fill(x_pos_024, 1./(eff_trigger));
                    TH1F_histos["y_pos_left_cut_top"]->Fill(y_pos_024, 1./(eff_trigger));
                }   
                if (rp_left_bottom_data && fid_cut_left_bottom){
                    TH1F_histos["xi_left_cut_bottom"]->Fill(xi_proton_left_data, 1/(eff_trigger));
                    TH1F_histos["t_left_cut_bottom"]->Fill(t_proton_left_data, 1/(eff_trigger));
                    TH1F_histos["log_x_left_cut_bottom"]->Fill(log10(x_left_data),1./(eff_trigger));
                    TH1F_histos["x_pos_left_cut_bottom"]->Fill(x_pos_025, 1./(eff_trigger));
                    TH1F_histos["y_pos_left_cut_bottom"]->Fill(y_pos_025, 1./(eff_trigger));
                }   
            }    
        }    

        //rp_pos
        if (jet_sel && proton_right_kin_sel && valid_proton_right_data && rp_right_top_data && !rp_right_bottom_data &&  fid_cut_right_top){
            TH2F_histos["proton_y_vs_x_rp_120_121_accept"]->Fill(x_pos_120, y_pos_120, 1.);
            TH2F_histos["proton_y_vs_x_rp_124_125_accept"]->Fill(x_pos_124, y_pos_124, 1.);
        }
        if (jet_sel && proton_right_kin_sel && valid_proton_right_data && !rp_right_top_data && rp_right_bottom_data && fid_cut_right_bottom){
            TH2F_histos["proton_y_vs_x_rp_120_121_accept"]->Fill(x_pos_121, y_pos_121, 1.);
            TH2F_histos["proton_y_vs_x_rp_124_125_accept"]->Fill(x_pos_125, y_pos_125, 1.);
        }
        if (jet_sel && proton_left_kin_sel && valid_proton_left_data && rp_left_top_data && !rp_left_bottom_data && fid_cut_left_top){
            TH2F_histos["proton_y_vs_x_rp_024_025_accept"]->Fill(x_pos_024, y_pos_024, 1.);
            // TH2F_histos["proton_y_vs_x_rp_024_025_accept"]->Fill(x_pos_025, y_pos_025, 1.);
          }
        if (jet_sel && proton_left_kin_sel && valid_proton_left_data && !rp_left_top_data && rp_left_bottom_data && fid_cut_left_bottom){
            // TH2F_histos["proton_y_vs_x_rp_024_025_accept"]->Fill(x_pos_024, y_pos_024, 1.);
            TH2F_histos["proton_y_vs_x_rp_024_025_accept"]->Fill(x_pos_025, y_pos_025, 1.);
        }


    }

    cout << "Jets: "<< nevents_jets << endl;
    cout << "Single Arm: "<< nevents_single_arm << endl; 
    cout << "Right rp: "<< nevents_rp_right << "  Left rp: "<< nevents_rp_left << endl; 
    cout << "Right kin: "<< nevents_kin_proton_right << "  Left kin: "<< nevents_kin_proton_left << endl; 
    cout << "Right events: "<< nevents_right << "  Left events: "<< nevents_left << endl; 

    for(map<string,TH1F*>::const_iterator histos = TH1F_histos.begin(); histos != TH1F_histos.end(); ++histos){
       histos->second->SetMarkerSize(1);
       histos->second->SetMarkerStyle(20);
       histos->second->SetMarkerColor(1);
       histos->second->SetLineColor(1);
       histos->second->Write();
    }   
    for(map<string,TH2F*>::const_iterator histos = TH2F_histos.begin(); histos != TH2F_histos.end(); ++histos){
       histos->second->Write();
    }   

    outfile->Close();

    return;
}

void DataAnalysis::AssymError (TH1F* histo_unfolded_nominal, TH1F* histo_jes_up_unfolded, TH1F* histo_jes_dw_unfolded, TH1F* histo_pf_up_unfolded, 
    TH1F* histo_pf_dw_unfolded, TH1F* histo_xi_up_unfolded, TH1F* histo_xi_dw_unfolded, TH1F* histo_norew_unfolded, TH1F* histo_betarew_unfolded, 
    TH1F* histo_rp_up_unfolded, TH1F* histo_rp_dw_unfolded, TH1F* histo_rp_x_unfolded, TH1F* histo_gauss_unfolded, TH1F* histo_hera_backg_unfolded,
    TH1F* histo_iter_unfolded, TH1F* histo_right_unfolded, TH1F* histo_left_unfolded, TH1F* histo_unfolded_mc1, TH1F* histo_unfolded_mc2, 
    TH1F* histo_unfolded_mc3, TH1F* histo_mc1_unfolded_mc2, TH1F* histo_mc1_unfolded_mc3, TH1F* histo_mc2_unfolded_mc1, TH1F* histo_mc2_unfolded_mc3,
    TH1F* histo_mc3_unfolded_mc1, TH1F* histo_mc3_unfolded_mc2, TH1F* histo_gen_mc1, TH1F* histo_gen_mc2, TH1F* histo_gen_mc3, 
    TH1F* histo_trigger_up, TH1F* histo_trigger_dw, TH1F* y_ref, TGraphAsymmErrors* &error,
    string const& result, bool histos){
    //result = "t" || "xi" || "x" || "ratio"
    //systematic error band

    int Nbins = histo_unfolded_nominal->GetNbinsX(); 
    float bins[Nbins]; float bins_histo[Nbins]; double bin = -4;
    float t_bins[10] = {0, 0.03, 0.07, 0.11, 0.21, 0.31, 0.45,  0.58, 0.72, 1.};
    float xibins[6] = {0, 0.028, 0.048, 0.065, 0.08, 0.1};
    for (int i = 0; i <= Nbins+1; ++i){
        if (result == "t") bins[i] = t_bins[i]; 
        if (result == "xi") bins[i] = xibins[i];
        if (result == "x" || result == "ratio") {bins[i] = bin; bin = bin + (4./15.);}
        // cout<<bins[i]<<endl;
    }  
    bin = -4;  
    for (int i = 0; i <= Nbins; ++i){
        if (result == "t") bins_histo[i] = tbins[i]; 
        if (result == "xi") bins_histo[i] = xi_bins[i];
        if (result == "x" || result == "ratio") {bins_histo[i] = bin; bin = bin + (4./15.);}
        // cout<<bins_histo[i]<<endl;
    }    
    TH1F* unc_jesup = new TH1F("unc_jesup", "", Nbins, bins_histo);
    TH1F* unc_jesdw = new TH1F("unc_jesdw", "", Nbins, bins_histo);
    TH1F* unc_pfup = new TH1F("unc_pfup", "", Nbins, bins_histo);
    TH1F* unc_pfdw = new TH1F("unc_pfdw", "", Nbins, bins_histo);
    TH1F* unc_xiup = new TH1F("unc_xiup", "", Nbins, bins_histo);
    TH1F* unc_xidw = new TH1F("unc_xidw", "", Nbins, bins_histo);
    TH1F* unc_rew_up = new TH1F("unc_rew_up", "", Nbins, bins_histo);
    TH1F* unc_rew_dw = new TH1F("unc_rew_dw", "", Nbins, bins_histo);
    TH1F* unc_rp_up = new TH1F("unc_rp_up", "", Nbins, bins_histo);
    TH1F* unc_rp_dw = new TH1F("unc_rp_dw", "", Nbins, bins_histo);
    TH1F* unc_slope_up = new TH1F("unc_slope_up", "", Nbins, bins_histo);
    TH1F* unc_slope_dw = new TH1F("unc_slope_dw", "", Nbins, bins_histo);
    TH1F* unc_gauss_up = new TH1F("unc_gauss_up", "", Nbins, bins_histo);
    TH1F* unc_gauss_dw = new TH1F("unc_gauss_dw", "", Nbins, bins_histo);
    TH1F* unc_backg_up = new TH1F("unc_backg_up", "", Nbins, bins_histo);
    TH1F* unc_backg_dw = new TH1F("unc_backg_dw", "", Nbins, bins_histo);
    TH1F* unc_iter_up = new TH1F("unc_iter_up", "", Nbins, bins_histo);
    TH1F* unc_iter_dw = new TH1F("unc_iter_dw", "", Nbins, bins_histo);
    TH1F* unc_sector_up = new TH1F("unc_sector_up", "", Nbins, bins_histo);
    TH1F* unc_sector_dw = new TH1F("unc_sector_dw", "", Nbins, bins_histo);
    TH1F* unc_accep_up = new TH1F("unc_accep_up", "", Nbins, bins_histo);
    TH1F* unc_accep_dw = new TH1F("unc_accep_dw", "", Nbins, bins_histo);
    TH1F* unc_bias_up = new TH1F("unc_bias_up", "", Nbins, bins_histo);
    TH1F* unc_bias_dw = new TH1F("unc_bias_dw", "", Nbins, bins_histo);
    TH1F* unc_trigger_up = new TH1F("unc_trigger_up", "", Nbins, bins_histo);
    TH1F* unc_trigger_dw = new TH1F("unc_trigger_dw", "", Nbins, bins_histo);
    TH1F* bias_mc1_unfolded_mc2 = new TH1F("bias_mc1_unfolded_mc2", "", Nbins, bins_histo);
    TH1F* bias_mc1_unfolded_mc3 = new TH1F("bias_mc1_unfolded_mc3", "", Nbins, bins_histo);
    TH1F* bias_mc2_unfolded_mc1 = new TH1F("bias_mc2_unfolded_mc1", "", Nbins, bins_histo);
    TH1F* bias_mc2_unfolded_mc3 = new TH1F("bias_mc2_unfolded_mc3", "", Nbins, bins_histo);
    TH1F* bias_mc3_unfolded_mc1 = new TH1F("bias_mc3_unfolded_mc1", "", Nbins, bins_histo);
    TH1F* bias_mc3_unfolded_mc2 = new TH1F("bias_mc3_unfolded_mc2", "", Nbins, bins_histo);
    TH1F* bias_mc1_unfolded_mc2_dw = new TH1F("bias_mc1_unfolded_mc2", "", Nbins, bins_histo);
    TH1F* bias_mc1_unfolded_mc3_dw = new TH1F("bias_mc1_unfolded_mc3", "", Nbins, bins_histo);
    TH1F* bias_mc2_unfolded_mc1_dw = new TH1F("bias_mc2_unfolded_mc1", "", Nbins, bins_histo);
    TH1F* bias_mc2_unfolded_mc3_dw = new TH1F("bias_mc2_unfolded_mc3", "", Nbins, bins_histo);
    TH1F* bias_mc3_unfolded_mc1_dw = new TH1F("bias_mc3_unfolded_mc1", "", Nbins, bins_histo);
    TH1F* bias_mc3_unfolded_mc2_dw = new TH1F("bias_mc3_unfolded_mc2", "", Nbins, bins_histo);

    double y[Nbins+1], y_nom[Nbins+1], y_up[Nbins+1], y_dw[Nbins+1], x[Nbins+1], deltaX[Nbins+1], jes_up[Nbins+1], jes_dw[Nbins+1], jes_up_dw[Nbins+1], jes_dw_up[Nbins+1];
    double xi_up[Nbins+1], xi_dw[Nbins+1], xi_up_dw[Nbins+1], xi_dw_up[Nbins+1], pf_up[Nbins+1], pf_dw[Nbins+1], pf_up_dw[Nbins+1], pf_dw_up[Nbins+1];
    double rew_beta[Nbins+1], rp_x[Nbins+1], rp_y_up[Nbins+1], rp_y_dw[Nbins+1], rp[Nbins+1], gauss[Nbins+1], b_slope[Nbins+1], backg[Nbins+1], iteration[Nbins+1], sector[Nbins+1];
    double bias[6], max_diff_bias[Nbins+1], max_diff_accep[Nbins+1], min_diff_accep[Nbins+1], trigger[Nbins+1];

    double av_jesup = 0; double av_jesdw = 0; double av_pfup = 0; double av_pfdw = 0; double av_xiup = 0; double av_xidw = 0; double av_backg = 0;
    double av_rew = 0; double av_gauss = 0; double av_slope = 0; double av_rp = 0; double av_iter = 0; double av_sector = 0; double av_accep = 0; double av_bias = 0;
    double max_jesup = -100; double max_jesdw = -100; double max_pfup = -100; double max_pfdw = -100; double max_xiup = -100; double max_xidw = -100; double max_backg = -100;
    double max_rew = -100; double max_gauss = -100; double max_slope = -100; double max_rp = -100; double max_iter = -100; double max_sector = -100; double max_accep = -100; double max_bias = -100;
    double av_trigger = 0; double max_trigger = -100;

    TH1F* accep[3]; 
    accep[0] = histo_unfolded_mc1; 
    accep[1] = histo_unfolded_mc2;
    accep[2] = histo_unfolded_mc3;

    int initial = 1;
    if (result == "xi") {Nbins = 6; initial = 2;}; 
    if (result == "x" || result == "ratio") {Nbins = 10; initial = 4; bin = -4+(3*4./15.);}; 


    for(int i = initial; i <= Nbins; ++i){
    // for(int i = 4; i <= 4; ++i){
    // for(int i = 2; i <= 2; ++i){
    // for(int i = 6; i <= 6; ++i){
    // for(int i = 8; i <= 8; ++i){
 
        y[i] = y_ref->GetBinContent(i); 
        y_nom[i] = histo_unfolded_nominal->GetBinContent(i);
        if (result == "t"){
            x[i] = (bins[i+1] + bins[i])/2.;
            deltaX[i] = x[i] - bins[i]; 
        }      
        if (result == "xi"){
            x[i] = (bins[i-1] + bins[i-2])/2.;
            deltaX[i] = x[i] - bins[i-2];   
        }
        if (result == "x" || result == "ratio") { 
            deltaX[i] = 2./15.; 
            x[i] = bin + (2./15.);//cout<<i<<"  "<<bin<<"  "<<x_logx[i]<<endl;
            bin = bin + (4./15.); 
        }    

        jes_up[i] = (histo_jes_up_unfolded->GetBinContent(i) > y_nom[i]) ? histo_jes_up_unfolded->GetBinContent(i)/y_nom[i]: 1 ;
        jes_up_dw[i] = (histo_jes_up_unfolded->GetBinContent(i) < y_nom[i]) ? y_nom[i]/histo_jes_up_unfolded->GetBinContent(i) : 1;
        jes_dw[i] = (y_nom[i] > histo_jes_dw_unfolded->GetBinContent(i)) ? y_nom[i]/histo_jes_dw_unfolded->GetBinContent(i) : 1 ;
        jes_dw_up[i] = (y_nom[i] < histo_jes_dw_unfolded->GetBinContent(i)) ? histo_jes_dw_unfolded->GetBinContent(i)/y_nom[i] : 1;
        if (jes_up[i] > jes_dw_up[i]) jes_up[i] = jes_up[i];
        else jes_up[i] = jes_dw_up[i];
        if (jes_dw[i] > jes_up_dw[i]) jes_dw[i] = jes_dw[i];
        else jes_dw[i] = jes_up_dw[i];

        pf_up[i] = (histo_pf_up_unfolded->GetBinContent(i) > y_nom[i]) ? histo_pf_up_unfolded->GetBinContent(i)/y_nom[i]: 1 ;
        pf_up_dw[i] = (histo_pf_up_unfolded->GetBinContent(i) < y_nom[i]) ? y_nom[i]/histo_pf_up_unfolded->GetBinContent(i) : 1;
        pf_dw[i] = (y_nom[i] > histo_pf_dw_unfolded->GetBinContent(i)) ? y_nom[i]/histo_pf_dw_unfolded->GetBinContent(i) : 1 ;
        pf_dw_up[i] = (y_nom[i] < histo_pf_dw_unfolded->GetBinContent(i)) ? histo_pf_dw_unfolded->GetBinContent(i)/y_nom[i] : 1;
        if (pf_up[i] > pf_dw_up[i]) pf_up[i] = pf_up[i];
        else pf_up[i] = pf_dw_up[i];
        if (pf_dw[i] > pf_up_dw[i]) pf_dw[i] = pf_dw[i];
        else pf_dw[i] = pf_up_dw[i];

        xi_up[i] = (histo_xi_up_unfolded->GetBinContent(i) > y_nom[i]) ? histo_xi_up_unfolded->GetBinContent(i)/y_nom[i] : 1;
        xi_up_dw[i] = (histo_xi_up_unfolded->GetBinContent(i) < y_nom[i]) ? y_nom[i]/histo_xi_up_unfolded->GetBinContent(i) : 1;
        xi_dw[i] = (y_nom[i] > histo_xi_dw_unfolded->GetBinContent(i)) ? y_nom[i]/histo_xi_dw_unfolded->GetBinContent(i) : 1;
        xi_dw_up[i] = (y_nom[i] < histo_xi_dw_unfolded->GetBinContent(i)) ? histo_xi_dw_unfolded->GetBinContent(i)/y_nom[i] : 1;
        if (xi_up[i] > xi_dw_up[i]) xi_up[i] = xi_up[i];
        else xi_up[i] = xi_dw_up[i];
        if (xi_dw[i] > xi_up_dw[i]) xi_dw[i] = xi_dw[i];
        else xi_dw[i] = xi_up_dw[i];
        if (xi_dw[i]>2) xi_dw[i] = 1;

        rew_beta[i] = fabs(histo_norew_unfolded->GetBinContent(i) - histo_betarew_unfolded->GetBinContent(i))/histo_betarew_unfolded->GetBinContent(i);
        rp_y_up[i] = fabs(histo_rp_up_unfolded->GetBinContent(i) - y_nom[i])/y_nom[i];
        rp_y_dw[i] = fabs(histo_rp_dw_unfolded->GetBinContent(i) - y_nom[i])/y_nom[i];
        rp_x[i] = fabs(histo_rp_x_unfolded->GetBinContent(i) - y_nom[i])/y_nom[i];
        rp[i] = sqrt(pow(rp_y_up[i], 2) + pow(rp_y_dw[i], 2) + pow(rp_x[i], 2));
        gauss[i] = fabs(histo_gauss_unfolded->GetBinContent(i) - y_nom[i])/y_nom[i];
        b_slope[i] =  fabs(histo_betarew_unfolded->GetBinContent(i) - y_nom[i])/y_nom[i];
        backg[i] = fabs(histo_hera_backg_unfolded->GetBinContent(i) - y_nom[i])/y_nom[i];
        iteration[i] = fabs(histo_iter_unfolded->GetBinContent(i) - y_nom[i])/y_nom[i];
        sector[i] = fabs(histo_right_unfolded->GetBinContent(i) - histo_left_unfolded->GetBinContent(i))/y_nom[i];
        trigger[i] = fabs(histo_trigger_up->GetBinContent(i) - histo_trigger_dw->GetBinContent(i))/y_nom[i];

        bias[0] = fabs(histo_mc1_unfolded_mc2->GetBinContent(i)-histo_gen_mc1->GetBinContent(i))/histo_gen_mc1->GetBinContent(i);
        bias[1] = fabs(histo_mc1_unfolded_mc3->GetBinContent(i)-histo_gen_mc1->GetBinContent(i))/histo_gen_mc1->GetBinContent(i);
        bias[2] = fabs(histo_mc2_unfolded_mc1->GetBinContent(i)-histo_gen_mc2->GetBinContent(i))/histo_gen_mc2->GetBinContent(i);
        bias[3] = fabs(histo_mc2_unfolded_mc3->GetBinContent(i)-histo_gen_mc2->GetBinContent(i))/histo_gen_mc2->GetBinContent(i);
        bias[4] = fabs(histo_mc3_unfolded_mc1->GetBinContent(i)-histo_gen_mc3->GetBinContent(i))/histo_gen_mc3->GetBinContent(i);
        bias[5] = fabs(histo_mc3_unfolded_mc2->GetBinContent(i)-histo_gen_mc3->GetBinContent(i))/histo_gen_mc3->GetBinContent(i);
        max_diff_bias[i] = -10;           
        for (int j = 0; j<6; ++j){
            if (bias[j] > max_diff_bias[i]) max_diff_bias[i] = bias[j];
        }

        max_diff_accep[i] = -10;
        min_diff_accep[i] = 10000;
        for (int j = 0; j<3; ++j){
            if (accep[j]->GetBinContent(i) > max_diff_accep[i]) max_diff_accep[i] = accep[j]->GetBinContent(i);
            if (accep[j]->GetBinContent(i) < min_diff_accep[i]) min_diff_accep[i] = accep[j]->GetBinContent(i);
        }

        y_up[i] = y[i]*sqrt(pow(jes_up[i]-1, 2) + pow(pf_up[i]-1, 2) + pow(xi_up[i]-1, 2) + pow(fabs(backg[i])/2, 2) + pow(fabs(rew_beta[i])/2, 2) +  
              pow(fabs(rp[i])/2, 2) + pow(fabs(gauss[i])/2, 2) + pow(fabs(sector[i])/2, 2) + pow(fabs(b_slope[i])/2, 2) + pow(fabs(trigger[i])/2, 2) + 
              pow(fabs(iteration[i])/2, 2) + + pow(max_diff_bias[i]/2,2) + pow((max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2, 2));
        y_dw[i] = y[i]*sqrt(pow(jes_dw[i]-1, 2) + pow(pf_dw[i]-1, 2) + pow(xi_dw[i]-1, 2) + pow(fabs(backg[i])/2, 2) + pow(fabs(rew_beta[i])/2, 2) +
              pow(fabs(rp[i])/2, 2) + pow(fabs(gauss[i])/2, 2) + pow(fabs(sector[i])/2, 2) + pow(fabs(b_slope[i])/2, 2) + pow(fabs(iteration[i])/2, 2) + 
              pow(max_diff_bias[i]/2, 2) + + pow((max_diff_accep[i] - min_diff_accep[i])/y_nom[i]/2, 2) + pow(fabs(trigger[i])/2, 2));
 
        unc_jesup->SetBinContent(i, jes_up[i]-1);
        unc_jesdw->SetBinContent(i, jes_dw[i]-1);
        unc_pfup->SetBinContent(i, pf_up[i]-1);
        unc_pfdw->SetBinContent(i, pf_dw[i]-1);
        unc_xiup->SetBinContent(i, xi_up[i]-1);
        unc_xidw->SetBinContent(i, xi_dw[i]-1);
        unc_backg_up->SetBinContent(i, fabs(backg[i])/2);
        unc_backg_dw->SetBinContent(i, -fabs(backg[i])/2);
        unc_rew_up->SetBinContent(i, fabs(rew_beta[i])/2);
        unc_rew_dw->SetBinContent(i, -fabs(rew_beta[i])/2);
        unc_slope_up->SetBinContent(i, fabs(b_slope[i])/2);
        unc_slope_dw->SetBinContent(i, -fabs(b_slope[i])/2);
        unc_rp_up->SetBinContent(i, fabs(rp[i])/2);
        unc_rp_dw->SetBinContent(i, -fabs(rp[i])/2);
        unc_gauss_up->SetBinContent(i, fabs(gauss[i])/2);
        unc_gauss_dw->SetBinContent(i, -fabs(gauss[i])/2);
        unc_sector_up->SetBinContent(i, fabs(sector[i])/2);
        unc_sector_dw->SetBinContent(i, -fabs(sector[i])/2);
        unc_iter_up->SetBinContent(i, fabs(iteration[i])/2);
        unc_iter_dw->SetBinContent(i, -fabs(iteration[i])/2);
        unc_accep_up->SetBinContent(i, (max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2);
        unc_accep_dw->SetBinContent(i, -(max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2);
        unc_bias_up->SetBinContent(i, max_diff_bias[i]/2);
        unc_bias_dw->SetBinContent(i, -max_diff_bias[i]/2);
        unc_trigger_up->SetBinContent(i, fabs(trigger[i])/2);
        unc_trigger_dw->SetBinContent(i, -fabs(trigger[i])/2);

        bias_mc1_unfolded_mc2->SetBinContent(i, bias[0]);
        bias_mc1_unfolded_mc3->SetBinContent(i, bias[1]);
        bias_mc2_unfolded_mc1->SetBinContent(i, bias[2]);
        bias_mc2_unfolded_mc3->SetBinContent(i, bias[3]);
        bias_mc3_unfolded_mc1->SetBinContent(i, bias[4]);
        bias_mc3_unfolded_mc2->SetBinContent(i, bias[5]);
        bias_mc1_unfolded_mc2_dw->SetBinContent(i, -bias[0]);
        bias_mc1_unfolded_mc3_dw->SetBinContent(i, -bias[1]);
        bias_mc2_unfolded_mc1_dw->SetBinContent(i, -bias[2]);
        bias_mc2_unfolded_mc3_dw->SetBinContent(i, -bias[3]);
        bias_mc3_unfolded_mc1_dw->SetBinContent(i, -bias[4]);
        bias_mc3_unfolded_mc2_dw->SetBinContent(i, -bias[5]);

        av_jesup += fabs(jes_up[i]-1);
        av_jesdw += fabs(jes_dw[i]-1);
        av_pfup += fabs(pf_up[i]-1);
        av_pfdw += fabs(pf_dw[i]-1);
        av_xiup += fabs(xi_up[i]-1);
        av_xidw += fabs(xi_dw[i]-1);
        av_backg += fabs(backg[i])/2;
        av_rew += fabs(rew_beta[i])/2;
        av_slope += fabs(b_slope[i])/2;
        av_rp += fabs(rp[i])/2;
        av_gauss += fabs(gauss[i])/2;
        av_sector += fabs(sector[i])/2;
        av_iter += fabs(iteration[i])/2;
        av_accep += (max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2;
        av_bias += max_diff_bias[i]/2;
        av_trigger += fabs(trigger[i])/2;
        if (fabs(jes_up[i]-1) > max_jesup) max_jesup = fabs(jes_up[i]-1);
        if (fabs(jes_dw[i]-1) > max_jesdw) max_jesdw = fabs(jes_dw[i]-1);
        if (fabs(pf_up[i]-1) > max_pfup) max_pfup = fabs(pf_up[i]-1);
        if (fabs(pf_dw[i]-1) > max_pfdw) max_pfdw = fabs(pf_dw[i]-1);
        if (fabs(xi_up[i]-1) > max_xiup) max_xiup = fabs(xi_up[i]-1);
        if (fabs(xi_dw[i]-1) > max_xidw) max_xidw = fabs(xi_dw[i]-1);
        if (fabs(backg[i])/2 > max_backg) max_backg = fabs(backg[i])/2;
        if (fabs(rew_beta[i])/2 > max_rew) max_rew = fabs(rew_beta[i])/2;
        if (fabs(b_slope[i])/2 > max_slope) max_slope = fabs(b_slope[i])/2;
        if (fabs(rp[i])/2 > max_rp) max_rp = fabs(rp[i])/2;
        if (fabs(gauss[i])/2 > max_gauss) max_gauss = fabs(gauss[i])/2;
        if (fabs(sector[i])/2 > max_sector) max_sector = fabs(sector[i])/2;
        if (fabs(iteration[i])/2 > max_iter) max_iter = fabs(iteration[i])/2;
        if ((max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2 > max_accep) max_accep = (max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2;
        if (max_diff_bias[i]/2 > max_bias) max_bias = max_diff_bias[i]/2;
        if (fabs(trigger[i])/2 > max_trigger) max_trigger = fabs(trigger[i])/2;

    }

    error = new TGraphAsymmErrors(Nbins+1, x, y, deltaX, deltaX, y_dw, y_up);

    int n_bins;
    if (result == "t") n_bins = 8; 
    if (result == "xi") n_bins = 5;
    if (result == "x" || result == "ratio") n_bins = 7;
    // n_bins = 1;

    cout <<  " " << endl;
    cout <<  "------------------ Max and Min systematic errors -------------------" << endl;
    cout << "JES:   av: " << av_jesup/n_bins << " / " << av_jesdw/n_bins << "   max: " << max_jesup << " / " << max_jesdw << endl;
    cout << "Trigger:   av: +/- " << av_trigger/n_bins << "   max: " << max_trigger << endl;
    cout << "PF:    av: " << av_pfup/n_bins << " / " << av_pfdw/n_bins << "   max: " << max_pfup << " / " << max_pfdw << endl;
    cout << "Background:   av: +/- " << av_backg/n_bins << "   max: " << max_backg << endl;
    cout << "RP:   av: +/- " << av_rp/n_bins << "   max: " << max_rp << endl;
    cout << "Beam divergence:   av: +/- " << av_gauss/n_bins << "   max: " << max_gauss << endl;
    cout << "Horizontal dispersion:  " << av_xiup/n_bins << " / " << av_xidw/n_bins << "   max: " << max_xiup << " / " << max_xidw << endl;
    cout << "t slope:   av: +/- " << av_slope/n_bins << "   max: " << max_slope << endl;
    cout << "Beta reweight:   av: +/- " << av_rew/n_bins << "   max: " << max_rew << endl;
    cout << "Acceptance & unfolding:   av: +/- " << av_accep/n_bins << "   max: " << max_accep << endl;
    cout << "Unfolding regularisation:   av: +/- " << av_iter/n_bins << "   max: " << max_iter << endl;
    cout << "Unfolding Bias:   av: +/- " << av_bias/n_bins << "   max: " << max_bias << endl;
    cout << "Sector:   av: +/- " << av_sector/n_bins << "   max: " << max_sector << endl;

    double total_av_up = sqrt(pow(av_jesup/n_bins, 2) + pow(av_pfup/n_bins, 2) + pow(av_backg/n_bins, 2) + pow(av_rp/n_bins, 2) + pow(av_gauss/n_bins, 2) +
        pow(av_xiup/n_bins, 2) + pow(av_slope/n_bins, 2) + pow(av_rew/n_bins, 2) + pow(av_accep/n_bins, 2) + pow(av_iter/n_bins, 2) + pow(av_bias/n_bins, 2) +
        pow(av_sector/n_bins, 2) + pow(av_trigger/n_bins, 2));
    double total_av_dw = sqrt(pow(av_jesdw/n_bins, 2) + pow(av_pfdw/n_bins, 2) + pow(av_backg/n_bins, 2) + pow(av_rp/n_bins, 2) + pow(av_gauss/n_bins, 2) +
        pow(av_xidw/n_bins, 2) + pow(av_slope/n_bins, 2) + pow(av_rew/n_bins, 2) + pow(av_accep/n_bins, 2) + pow(av_iter/n_bins, 2) + pow(av_bias/n_bins, 2) +
        pow(av_sector/n_bins, 2) + pow(av_trigger/n_bins, 2));

    double total_max_up = sqrt(pow(max_jesup, 2) + pow(max_pfup, 2) + pow(max_backg, 2) + pow(max_rp, 2) + pow(max_gauss, 2) +
        pow(max_xiup, 2) + pow(max_slope, 2) + pow(max_rew, 2) + pow(max_accep, 2) + pow(max_iter, 2) + pow(max_bias, 2) +
        pow(max_sector, 2) + pow(max_trigger, 2));
    double total_max_dw = sqrt(pow(max_jesdw, 2) + pow(max_pfdw, 2) + pow(max_backg, 2) + pow(max_rp, 2) + pow(max_gauss, 2) +
        pow(max_xidw, 2) + pow(max_slope, 2) + pow(max_rew, 2) + pow(max_accep, 2) + pow(max_iter, 2) + pow(max_bias, 2) +
        pow(max_sector, 2) + pow(max_trigger, 2));
    cout << "--------------------------------------------------------------------- " << endl;
    cout << "Total:   av: "<< total_av_up << " / " << total_av_dw << "   max: " << total_max_up << " / " << total_max_dw << endl;
    cout <<  " " << endl;
    cout <<  " " << endl;


    if (histos){
       TCanvas *c1 = new TCanvas("c1","unc");
        unc_jesup->SetLineColor(2);    
        unc_jesdw->SetLineColor(2);    
        unc_pfup->SetLineColor(4);    
        unc_pfdw->SetLineColor(4);    
        unc_xiup->SetLineColor(3);    
        unc_xidw->SetLineColor(3);    
        unc_backg_up->SetLineColor(6);    
        unc_backg_dw->SetLineColor(6);    
        unc_rew_up->SetLineColor(7);    
        unc_rew_dw->SetLineColor(7);    
        unc_rp_up->SetLineColor(41);    
        unc_rp_dw->SetLineColor(41);    
        unc_gauss_up->SetLineColor(29);    
        unc_gauss_dw->SetLineColor(29);    
        unc_sector_up->SetLineColor(13);    
        unc_sector_dw->SetLineColor(13);    
        unc_accep_up->SetLineColor(45);    
        unc_accep_dw->SetLineColor(45);    
        unc_iter_up->SetLineColor(21);    
        unc_iter_dw->SetLineColor(21);    
        unc_bias_up->SetLineColor(38);    
        unc_bias_dw->SetLineColor(38); 
        unc_trigger_up->SetLineColor(5);    
        unc_trigger_dw->SetLineColor(5);    
           
        if (result == "t"){ unc_jesup->SetXTitle("|t| (GeV)^{2}"); unc_jesup->SetYTitle("#Delta(|t|)"); }   
        if (result == "xi"){ unc_jesup->SetXTitle("#xi");   unc_jesup->SetYTitle("#Delta(#xi)"); }    
        if (result == "x"){ unc_jesup->SetXTitle("log_{10} x"); unc_jesup->SetYTitle("#Delta(log_{10} x)"); }      
        unc_jesup->Draw("hist");    
        unc_jesdw->Draw("histsame");    
        unc_pfup->Draw("histsame");    
        unc_pfdw->Draw("histsame");    
        unc_xiup->Draw("histsame");    
        unc_xidw->Draw("histsame");    
        unc_backg_up->Draw("histsame");    
        unc_backg_dw->Draw("histsame");    
        unc_rew_up->Draw("histsame");    
        unc_rew_dw->Draw("histsame");    
        unc_rp_up->Draw("histsame");    
        unc_rp_dw->Draw("histsame");    
        unc_gauss_up->Draw("histsame");    
        unc_gauss_dw->Draw("histsame");    
        // unc_trigger_up->Draw("histsame");    
        // unc_trigger_dw->Draw("histsame");    
        unc_sector_up->Draw("histsame");    
        unc_sector_dw->Draw("histsame");    
        unc_accep_up->Draw("histsame");    
        unc_accep_dw->Draw("histsame");    
        unc_iter_up->Draw("histsame");    
        unc_iter_dw->Draw("histsame");    
        unc_bias_up->Draw("histsame");    
        unc_bias_dw->Draw("histsame");  
       TLegend *leg1 = new TLegend(0.2,0.75,0.48,0.9);
       leg1->AddEntry(unc_jesup,"Jet Energy Scale","l");
       // leg1->AddEntry(unc_trigger_up,"Trigger efficiency","l");
       leg1->AddEntry(unc_pfup,"Calorimeter Energy Scale","l");
       leg1->AddEntry(unc_backg_up,"Background","l");
       leg1->AddEntry(unc_rp_up,"RPs acceptance","l");
       leg1->AddEntry(unc_gauss_up,"Beam divergence","l");
       leg1->AddEntry(unc_xiup,"Horizontal dispersion","l");
       leg1->AddEntry(unc_slope_up,"t-slope","l");
       leg1->AddEntry(unc_rew_up,"#beta-reweighting","l");
       leg1->AddEntry(unc_accep_up,"Acceptance and unfolding","l");
       leg1->AddEntry(unc_iter_up,"Unfolding regularisation","l");
       leg1->AddEntry(unc_bias_up,"Unfolding bias","l");
       leg1->AddEntry(unc_sector_up,"Sector","l");
       leg1->SetFillColor(0);
       leg1->SetLineColor(0);
       leg1->SetShadowColor(0);
       leg1->Draw();     

       TCanvas *c2 = new TCanvas("c2","bias");
        bias_mc1_unfolded_mc2->SetLineColor(2);
        bias_mc1_unfolded_mc3->SetLineColor(3);  
        bias_mc2_unfolded_mc1->SetLineColor(4);  
        bias_mc2_unfolded_mc3->SetLineColor(6);  
        bias_mc3_unfolded_mc1->SetLineColor(7);  
        bias_mc3_unfolded_mc2->SetLineColor(41);  
        bias_mc1_unfolded_mc2_dw->SetLineColor(2);
        bias_mc1_unfolded_mc3_dw->SetLineColor(3);  
        bias_mc2_unfolded_mc1_dw->SetLineColor(4);  
        bias_mc2_unfolded_mc3_dw->SetLineColor(6);  
        bias_mc3_unfolded_mc1_dw->SetLineColor(7);  
        bias_mc3_unfolded_mc2_dw->SetLineColor(41);  
            
        if (result == "t"){ bias_mc1_unfolded_mc2->SetXTitle("|t| (GeV)^{2}"); bias_mc1_unfolded_mc2->SetYTitle("#Delta(|t|)");}   
        if (result == "xi"){ bias_mc1_unfolded_mc2->SetXTitle("#xi"); bias_mc1_unfolded_mc2->SetYTitle("#Delta(#xi)");}      
        if (result == "x"){ bias_mc1_unfolded_mc2->SetXTitle("log_{10} x"); bias_mc1_unfolded_mc2->SetYTitle("#Delta(log_{10} x)");}       
        bias_mc1_unfolded_mc2->Draw("hist");
        bias_mc1_unfolded_mc3->Draw("histsame");  
        bias_mc2_unfolded_mc1->Draw("histsame");  
        bias_mc2_unfolded_mc3->Draw("histsame");  
        bias_mc3_unfolded_mc1->Draw("histsame");  
        bias_mc3_unfolded_mc2->Draw("histsame");  
        bias_mc1_unfolded_mc2_dw->Draw("histsame");
        bias_mc1_unfolded_mc3_dw->Draw("histsame");  
        bias_mc2_unfolded_mc1_dw->Draw("histsame");  
        bias_mc2_unfolded_mc3_dw->Draw("histsame");  
        bias_mc3_unfolded_mc1_dw->Draw("histsame");  
        bias_mc3_unfolded_mc2_dw->Draw("histsame");  
       TLegend *leg2 = new TLegend(0.2,0.75,0.48,0.9);
       leg2->AddEntry(bias_mc1_unfolded_mc2,"POMWIG unfolded with PYTHIA8 4C","l");
       leg2->AddEntry(bias_mc1_unfolded_mc3,"POMWIG unfolded with PYTHIA8 CUETP8M1","l");
       leg2->AddEntry(bias_mc2_unfolded_mc1,"PYTHIA8 4C unfolded with POMWIG","l");
       leg2->AddEntry(bias_mc2_unfolded_mc3,"PYTHIA8 4C unfolded with PYTHIA8 CUETP8M1","l");
       leg2->AddEntry(bias_mc3_unfolded_mc1,"PYTHIA8 CUETP8M1 unfolded with POMWIG","l");
       leg2->AddEntry(bias_mc3_unfolded_mc2,"PYTHIA8 CUETP8M1 unfolded with PYTHIA8 4C","l");
       leg2->SetFillColor(0);
       leg2->SetLineColor(0);
       leg2->SetShadowColor(0);
       leg2->Draw();     
    }
  
}

void DataAnalysis::AssymErrorInclusive (TH1F* histo_unfolded_nominal, TH1F* histo_jes_up_unfolded, TH1F* histo_jes_dw_unfolded, TH1F* histo_pf_up_unfolded, 
    TH1F* histo_pf_dw_unfolded, TH1F* histo_iter_unfolded, TH1F* histo_right_unfolded, TH1F* histo_left_unfolded, TH1F* histo_unfolded_mc1, TH1F* histo_unfolded_mc2, 
    TH1F* histo_unfolded_mc3, TH1F* histo_mc1_unfolded_mc2, TH1F* histo_mc1_unfolded_mc3, TH1F* histo_mc2_unfolded_mc1, TH1F* histo_mc2_unfolded_mc3,
    TH1F* histo_mc3_unfolded_mc1, TH1F* histo_mc3_unfolded_mc2, TH1F* histo_gen_mc1, TH1F* histo_gen_mc2, TH1F* histo_gen_mc3, TH1F* y_ref, TGraphAsymmErrors* &error){
    //systematic error band

    int Nbins = histo_unfolded_nominal->GetNbinsX(); 
    float bins[Nbins]; float bins_histo[Nbins]; double bin =-4+(3*4./15.);

    double y[Nbins+1], y_nom[Nbins+1], y_up[Nbins+1], y_dw[Nbins+1], x[Nbins+1], deltaX[Nbins+1], jes_up[Nbins+1], jes_dw[Nbins+1], jes_up_dw[Nbins+1], jes_dw_up[Nbins+1];
    double pf_up[Nbins+1], pf_dw[Nbins+1], pf_up_dw[Nbins+1], pf_dw_up[Nbins+1];
    double iteration[Nbins+1], sector[Nbins+1];
    double bias[6], max_diff_bias[Nbins+1], max_diff_accep[Nbins+1], min_diff_accep[Nbins+1];

    double av_jesup = 0; double av_jesdw = 0; double av_pfup = 0; double av_pfdw = 0; 
    double av_iter = 0; double av_sector = 0; double av_accep = 0; double av_bias = 0;
    double max_jesup = -100; double max_jesdw = -100; double max_pfup = -100; double max_pfdw = -100;
    double max_iter = -100; double max_sector = -100; double max_accep = -100; double max_bias = -100;

    TH1F* accep[3]; 
    accep[0] = histo_unfolded_mc1; 
    accep[1] = histo_unfolded_mc2;
    accep[2] = histo_unfolded_mc3;

    for(int i = 4; i <= 10; ++i){
 
        y[i] = y_ref->GetBinContent(i); 
        y_nom[i] = histo_unfolded_nominal->GetBinContent(i);

        deltaX[i] = 2./15.; 
        x[i] = bin + (2./15.);//cout<<i<<"  "<<bin<<"  "<<x_logx[i]<<endl;
        bin = bin + (4./15.); 

        jes_up[i] = (histo_jes_up_unfolded->GetBinContent(i) > y_nom[i]) ? histo_jes_up_unfolded->GetBinContent(i)/y_nom[i]: 1 ;
        jes_up_dw[i] = (histo_jes_up_unfolded->GetBinContent(i) < y_nom[i]) ? y_nom[i]/histo_jes_up_unfolded->GetBinContent(i) : 1;
        jes_dw[i] = (y_nom[i] > histo_jes_dw_unfolded->GetBinContent(i)) ? y_nom[i]/histo_jes_dw_unfolded->GetBinContent(i) : 1 ;
        jes_dw_up[i] = (y_nom[i] < histo_jes_dw_unfolded->GetBinContent(i)) ? histo_jes_dw_unfolded->GetBinContent(i)/y_nom[i] : 1;
        if (jes_up[i] > jes_dw_up[i]) jes_up[i] = jes_up[i];
        else jes_up[i] = jes_dw_up[i];
        if (jes_dw[i] > jes_up_dw[i]) jes_dw[i] = jes_dw[i];
        else jes_dw[i] = jes_up_dw[i];

        pf_up[i] = (histo_pf_up_unfolded->GetBinContent(i) > y_nom[i]) ? histo_pf_up_unfolded->GetBinContent(i)/y_nom[i]: 1 ;
        pf_up_dw[i] = (histo_pf_up_unfolded->GetBinContent(i) < y_nom[i]) ? y_nom[i]/histo_pf_up_unfolded->GetBinContent(i) : 1;
        pf_dw[i] = (y_nom[i] > histo_pf_dw_unfolded->GetBinContent(i)) ? y_nom[i]/histo_pf_dw_unfolded->GetBinContent(i) : 1 ;
        pf_dw_up[i] = (y_nom[i] < histo_pf_dw_unfolded->GetBinContent(i)) ? histo_pf_dw_unfolded->GetBinContent(i)/y_nom[i] : 1;
        if (pf_up[i] > pf_dw_up[i]) pf_up[i] = pf_up[i];
        else pf_up[i] = pf_dw_up[i];
        if (pf_dw[i] > pf_up_dw[i]) pf_dw[i] = pf_dw[i];
        else pf_dw[i] = pf_up_dw[i];

        iteration[i] = fabs(histo_iter_unfolded->GetBinContent(i) - y_nom[i])/y_nom[i];
        sector[i] = fabs(histo_right_unfolded->GetBinContent(i) - histo_left_unfolded->GetBinContent(i))/y_nom[i];

        bias[0] = fabs(histo_mc1_unfolded_mc2->GetBinContent(i)-histo_gen_mc1->GetBinContent(i))/histo_gen_mc1->GetBinContent(i);
        bias[1] = fabs(histo_mc1_unfolded_mc3->GetBinContent(i)-histo_gen_mc1->GetBinContent(i))/histo_gen_mc1->GetBinContent(i);
        bias[2] = fabs(histo_mc2_unfolded_mc1->GetBinContent(i)-histo_gen_mc2->GetBinContent(i))/histo_gen_mc2->GetBinContent(i);
        bias[3] = fabs(histo_mc2_unfolded_mc3->GetBinContent(i)-histo_gen_mc2->GetBinContent(i))/histo_gen_mc2->GetBinContent(i);
        bias[4] = fabs(histo_mc3_unfolded_mc1->GetBinContent(i)-histo_gen_mc3->GetBinContent(i))/histo_gen_mc3->GetBinContent(i);
        bias[5] = fabs(histo_mc3_unfolded_mc2->GetBinContent(i)-histo_gen_mc3->GetBinContent(i))/histo_gen_mc3->GetBinContent(i);
        max_diff_bias[i] = -10;           
        for (int j = 0; j<6; ++j){
            if (bias[j] > max_diff_bias[i]) max_diff_bias[i] = bias[j];
        }

        max_diff_accep[i] = -10;
        min_diff_accep[i] = 1000000;
        for (int j = 0; j<3; ++j){
            if (accep[j]->GetBinContent(i) > max_diff_accep[i]) max_diff_accep[i] = accep[j]->GetBinContent(i);
            if (accep[j]->GetBinContent(i) < min_diff_accep[i]) min_diff_accep[i] = accep[j]->GetBinContent(i);
            cout<<i <<j << "  "<< max_diff_accep[i] << "  "<< min_diff_accep[i]<<endl;
        }

        y_up[i] = y[i]*sqrt(pow(jes_up[i]-1, 2) + pow(pf_up[i]-1, 2) + pow(fabs(sector[i])/2, 2) + pow(fabs(iteration[i])/2, 2) 
          + pow(max_diff_bias[i]/2,2) + pow((max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2, 2));
        y_dw[i] = y[i]*sqrt(pow(jes_dw[i]-1, 2) + pow(pf_dw[i]-1, 2) + pow(fabs(sector[i])/2, 2) + pow(fabs(iteration[i])/2, 2) + 
              pow(max_diff_bias[i]/2, 2) + pow((max_diff_accep[i] - min_diff_accep[i])/y_nom[i]/2, 2));
 

        av_jesup += fabs(jes_up[i]-1);
        av_jesdw += fabs(jes_dw[i]-1);
        av_pfup += fabs(pf_up[i]-1);
        av_pfdw += fabs(pf_dw[i]-1);
        av_sector += fabs(sector[i])/2;
        av_iter += fabs(iteration[i])/2;
        av_accep += (max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2;
        av_bias += max_diff_bias[i]/2;
        if (fabs(jes_up[i]-1) > max_jesup) max_jesup = fabs(jes_up[i]-1);
        if (fabs(jes_dw[i]-1) > max_jesdw) max_jesdw = fabs(jes_dw[i]-1);
        if (fabs(pf_up[i]-1) > max_pfup) max_pfup = fabs(pf_up[i]-1);
        if (fabs(pf_dw[i]-1) > max_pfdw) max_pfdw = fabs(pf_dw[i]-1);
        if (fabs(sector[i])/2 > max_sector) max_sector = fabs(sector[i])/2;
        if (fabs(iteration[i])/2 > max_iter) max_iter = fabs(iteration[i])/2;
        if ((max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2 > max_accep) max_accep = (max_diff_accep[i]-min_diff_accep[i])/y_nom[i]/2;
        if (max_diff_bias[i]/2 > max_bias) max_bias = max_diff_bias[i]/2;

    }

    error = new TGraphAsymmErrors(Nbins+1, x, y, deltaX, deltaX, y_dw, y_up);

    int n_bins = 7;

    cout <<  " " << endl;
    cout <<  "------------------ Max and Min systematic errors -------------------" << endl;
    cout << "JES:   av: " << av_jesup/n_bins << " / " << av_jesdw/n_bins << "   max: " << max_jesup << " / " << max_jesdw << endl;
    cout << "PF:    av: " << av_pfup/n_bins << " / " << av_pfdw/n_bins << "   max: " << max_pfup << " / " << max_pfdw << endl;
    cout << "Acceptance & unfolding:   av: +/- " << av_accep/n_bins << "   max: " << max_accep << endl;
    cout << "Unfolding regularisation:   av: +/- " << av_iter/n_bins << "   max: " << max_iter << endl;
    cout << "Unfolding Bias:   av: +/- " << av_bias/n_bins << "   max: " << max_bias << endl;
    // cout << "Sector:   av: +/- " << av_sector/n_bins << "   max: " << max_sector << endl;

    double total_av_up = sqrt(pow(av_jesup/n_bins, 2) + pow(av_pfup/n_bins, 2) + pow(av_accep/n_bins, 2) + pow(av_iter/n_bins, 2) + pow(av_bias/n_bins, 2));
    double total_av_dw = sqrt(pow(av_jesdw/n_bins, 2) + pow(av_pfdw/n_bins, 2) + pow(av_accep/n_bins, 2) + pow(av_iter/n_bins, 2) + pow(av_bias/n_bins, 2));

    double total_max_up = sqrt(pow(max_jesup, 2) + pow(max_pfup, 2) + pow(max_accep, 2) + pow(max_iter, 2) + pow(max_bias, 2));
    double total_max_dw = sqrt(pow(max_jesdw, 2) + pow(max_pfdw, 2) + pow(max_accep, 2) + pow(max_iter, 2) + pow(max_bias, 2));
    cout << "--------------------------------------------------------------------- " << endl;
    cout << "Total:   av: "<< total_av_up << " / " << total_av_dw << "   max: " << total_max_up << " / " << total_max_dw << endl;
    cout <<  " " << endl;
    cout <<  " " << endl;
}

void DataAnalysis::AbsoluteSigma (TH1F* data, TH1F* mc_gen, TH1F* mc_rec, double &sigma_data, double &error_data){
    //Absolute cross section

    double sum_data_nominal = 0; double sum_error_data = 0; double sum_rec = 0; double sum_gen = 0; 
    double luminosity = 37.5;
    int nbins = data->GetNbinsX();
    int initial = 1;
    // if (result == "x" || result == "ratio") {Nbins = 10; initial = 4; bin = -4+(3*4./15.);}; 

    // for(int i = 4; i <= 10; ++i){
    for(int i = 1; i <= nbins; ++i){
        sum_data_nominal += data->GetBinContent(i); sum_error_data += pow(data->GetBinError(i),2);
        sum_rec += mc_rec->GetBinContent(i);
        sum_gen += mc_gen->GetBinContent(i);
    }

    sigma_data = sum_data_nominal*sum_gen/sum_rec/luminosity;
    error_data = sqrt(sum_error_data)*sum_gen/sum_rec/luminosity;

}

// }
DataAnalysis::~DataAnalysis (void) {}

int main() {
	DataAnalysis data;
    data.getHistos(true,false,false,false,false,false,false,false,false);
    // data.getHistos(false,false,false,false,false,false,false,true,false);
	// data.getHistos(false,false,false,false,false,false,false,false,true);
    DataAnalysis data2;
    // data2.getHistos(false,false,false,true,false,false,false,false,false);
}