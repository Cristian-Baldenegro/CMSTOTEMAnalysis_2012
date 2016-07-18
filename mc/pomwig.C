
#include <pomwig.h>
// #include "data.C"

using namespace std;

///ROOUNFOLD
#include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldBinByBin.h"
#include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldInvert.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldInvert.h"
 
   double scale_right = 1;//1.29963;
   double scale_rew = 0.8456880899762572;//0.84442877844358-promptreco;


void pomwig(TH1F* &xi_cms_minus_totem_rec_right, TH1F* &xi_cms_minus_totem_rec_left, TH1F* &xi_cms_minus_totem_rec_right_bin, TH1F* &xi_cms_minus_totem_rec_left_bin, TH1F* &xi_rec_right, TH1F* &xi_rec_left, 
            TH1F* &xi_rec_right_cut, TH1F* &xi_rec_right_cut_sasha, TH1F* &xi_gen_right_cut_sasha, TH1F* &xi_rec_left_cut_sasha,TH1F* &xi_rec_left_cut, TH1F* &xi_gen_right, TH1F* &xi_gen_left,TH1F* &t_rec_right_cut, TH1F* &t_gen_right_cut, 
            TH1F* &t_rec_left_cut, TH1F* &t_gen_left_cut, TH1F* &beta_right_cut, TH1F* &beta_left_cut, TH1F* &xi_cms_rec_right, TH1F* &xi_cms_gen_right,TH1F* &xi_cms_rec_left, TH1F* &xi_cms_gen_left,
            TH1F* &xi_cms_rec_right_sasha, TH1F* &xi_cms_gen_right_sasha,TH1F* &xi_cms_rec_left_sasha, TH1F* &xi_cms_gen_left_sasha, 
            TH1F* &log_x_rec_right_cut, TH1F* &log_x_gen_right_cut, TH1F* &log_x_rec_left_cut,  TH1F* &log_x_gen_left_cut, 
            TH1F* &pt_jet1_right_cut, TH1F* &pt_jet1_left_cut, TH1F* &pt_jet2_right_cut, TH1F* &pt_jet2_left_cut, TH1F* &eta_jet1_right_cut, TH1F* &eta_jet1_left_cut, TH1F* &eta_jet2_right_cut, 
            TH1F* &eta_jet2_left_cut,  TH1F* &t_right_unfolded,  TH1F* &t_left_unfolded, TH1F* &xi_right_unfolded,  TH1F* &xi_left_unfolded,
            TH1F* &logx_right_unfolded,  TH1F* &logx_left_unfolded, TH1F* &t_right_unfolded_noeff,  TH1F* &t_left_unfolded_noeff, TH1F* &xi_right_unfolded_noeff,  TH1F* &xi_left_unfolded_noeff,
            TH1F* &logx_right_cut_unfolded_noeff,  TH1F* &logx_left_cut_unfolded_noeff, TH1F* &logx_right_unfolded_noeff,  TH1F* &logx_left_unfolded_noeff,  TH1F* &t_right_unfolded_binbybin,  
            TH1F* &t_left_unfolded_binbybin, TH1F* &xi_right_unfolded_binbybin,  TH1F* &xi_left_unfolded_binbybin,
            TH1F* &logx_right_unfolded_binbybin,  TH1F* &logx_left_unfolded_binbybin,  TH1F* &t_right_unfolded_up,  TH1F* &t_right_unfolded_dw, TH1F* &t_left_unfolded_up,  TH1F* &t_left_unfolded_dw, TH1F* &xi_right_unfolded_up,  TH1F* &xi_right_unfolded_dw,
            TH1F* &xi_left_unfolded_up,  TH1F* &xi_left_unfolded_dw, TH1F* &log_x_right_unfolded_up,  TH1F* &log_x_right_unfolded_dw, TH1F* &log_x_left_unfolded_up,  TH1F* &log_x_left_unfolded_dw,  
            TH1F* &t_right_unfolded_up_pf,  TH1F* &t_right_unfolded_dw_pf, TH1F* &t_left_unfolded_up_pf,  TH1F* &t_left_unfolded_dw_pf, TH1F* &xi_right_unfolded_up_pf,  TH1F* &xi_right_unfolded_dw_pf,
            TH1F* &xi_left_unfolded_up_pf,  TH1F* &xi_left_unfolded_dw_pf, TH1F* &log_x_right_unfolded_up_pf,  TH1F* &log_x_right_unfolded_dw_pf, TH1F* &log_x_left_unfolded_up_pf,  TH1F* &log_x_left_unfolded_dw_pf, 
            TH1F* &t_right_unfolded_herabackg, TH1F* &t_left_unfolded_herabackg, TH1F* &xi_right_unfolded_herabackg, TH1F* &xi_left_unfolded_herabackg, TH1F* &logx_right_unfolded_herabackg, TH1F* &logx_left_unfolded_herabackg,
            bool reweight, bool reweight_slope, bool unc_rp, bool unc_gauss, bool sys_xi_up, bool sys_xi_dw, bool rp_x_min, bool rp_x_max){
 

 if (pt_threshold==30){
   func->SetParameter(0, 1.65945);
   func->SetParameter(1, -4.64686);
   func->SetParameter(2, 5.73927);
//   func_right->SetParameter(3, -52.6806);
//   func_right->SetParameter(4, 35.0543);

 }

 if (pt_threshold==35){
   func->SetParameter(0, 1.61473);
   func->SetParameter(1, -4.03202);
   func->SetParameter(2, 5.14115);
//   func_right->SetParameter(3, -52.6806);
//   func_right->SetParameter(4, 35.0543);

 }

 if (pt_threshold==40){
   // func->SetParameter(0, 1.52592);
   // func->SetParameter(1, -3.96672);
   // func->SetParameter(2, 4.33963);
   //rereco 
   func->SetParameter(0, 1.18178);//;/*1.81434);//*/1.67629);
   func->SetParameter(1, -2.7493);//*-5.09231);//*/-5.02845);
   func->SetParameter(2, 3.26971);//*6.11125);//*/6.07907);

 }
    xi_cms_minus_totem_rec_right = new TH1F("xi_cms_minus_totem_rec_right_pom","",50,-0.4,0.4);
    xi_cms_minus_totem_rec_right_bin = new TH1F("xi_cms_minus_totem_rec_right_pom_bin","", 15, bin);
    xi_cms_minus_totem_rec_left = new TH1F("xi_cms_minus_totem_rec_left_pom","",50,-0.4,0.4);
    xi_cms_minus_totem_rec_left_bin = new TH1F("xi_cms_minus_totem_rec_left_pom_bin","", 15, bin);
    xi_rec_right = new TH1F("xi_rec_right_pom","",50,-0.04,0.2);
    xi_rec_left = new TH1F("xi_rec_left_pom","",50,-0.04,0.2);
    xi_rec_right_cut = new TH1F("xi_rec_right_cut_pom","",11,xi_bins);
    xi_gen_right = new TH1F("xi_gen_right_cut_pom","",11,xi_bins);
    xi_rec_left_cut = new TH1F("xi_rec_left_cut_pom","",11,xi_bins);
    xi_gen_left = new TH1F("xi_gen_left_cut_pom","",11,xi_bins);
    t_rec_right_cut = new TH1F("t_rec_right_cut_pom","", 8, tbins);
    t_gen_right_cut = new TH1F("t_gen_right_cut_pom","", 8, tbins);
    t_rec_left_cut = new TH1F("t_rec_left_cut_pom","", 8, tbins);
    t_gen_left_cut = new TH1F("t_gen_left_cut_pom","", 8, tbins);
    beta_right_cut = new TH1F("beta_right_cut_pom","",15,0,1);
    beta_left_cut = new TH1F("beta_left_cut_pom","",15,0,1);
    TH2F* pt2_xi_gen_minus = new TH2F("pt2_xi_gen_minus_pom","",20,0,1,20,0,0.1);
    TH2F* pt2_xi_gen_minus_rp = new TH2F("pt2_xi_gen_minus_rp_pom","",20,0,1,20,0,0.1);
    xi_cms_rec_right = new TH1F("xi_cms_rec_right_pom","",11,xi_bins);
    xi_cms_rec_left = new TH1F("xi_cms_rec_left_pom","",11,xi_bins);
    xi_cms_gen_right = new TH1F("xi_cms_gen_right_pom","",11,xi_bins);
    xi_cms_gen_left = new TH1F("xi_cms_gen_left_pom","",11,xi_bins);
    xi_cms_rec_right_sasha = new TH1F("xi_cms_rec_right_sasha_pom","",8, bin_sasha);
    xi_cms_rec_left_sasha = new TH1F("xi_cms_rec_left_sasha_pom","",8, bin_sasha);
    xi_cms_gen_right_sasha = new TH1F("xi_cms_gen_right_sasha_pom","",8, bin_sasha);
    xi_cms_gen_left_sasha = new TH1F("xi_cms_gen_left_sasha_pom","",8, bin_sasha);
    log_x_rec_right_cut = new TH1F("log_x_rec_right_cut_pom","",15, -4, 0);
    log_x_gen_right_cut = new TH1F("log_x_gen_right_cut_pom","",15, -4, 0);
    log_x_rec_left_cut = new TH1F("log_x_rec_left_cut_pom","",15, -4, 0);
    log_x_gen_left_cut = new TH1F("log_x_gen_left_cut_pom","",15, -4, 0);
    pt_jet1_right_cut = new TH1F("pt_jet1_right_cut_pomwig","", 15, 0, 200);
    pt_jet2_right_cut = new TH1F("pt_jet2_right_cut_pomwig","", 15, 0, 200);
    eta_jet1_right_cut = new TH1F("eta_jet1_right_cut_pomwig","", 20, -5.2, 5.2);
    eta_jet2_right_cut = new TH1F("eta_jet2_right_cut_pomwig","", 20, -5.2, 5.2);
    pt_jet1_left_cut = new TH1F("pt_jet1_left_cut_pomwig","", 15, 0, 200);
    pt_jet2_left_cut = new TH1F("pt_jet2_left_cut_pomwig","", 15, 0, 200);
    eta_jet1_left_cut = new TH1F("eta_jet1_left_cut_pomwig","", 20, -5.2, 5.2);
    eta_jet2_left_cut = new TH1F("eta_jet2_left_cut_pomwig","", 20, -5.2, 5.2);
    TH2F* rp_pos_left = new TH2F("rp_pos_left","",100,-1,10,100,-40,40);
    TH2F* rp_pos_rigth = new TH2F("rp_pos_rigth","",100,-1,10,100,-40,40);
    TH2F* rp_pos_rigth_xmin = new TH2F("rp_pos_rigth","",100,-4,10,100,-40,40);
    TH2F* rp_pos_left_xmin = new TH2F("rp_pos_lefh","",100,-4,10,100,-40,40);
    TH2F* rp_pos_rigth_xmax = new TH2F("rp_pos_rigth","",100,-4,10,100,-40,40);
    TH2F* rp_pos_left_xmax = new TH2F("rp_pos_lefh","",100,-4,10,100,-40,40);
    xi_rec_right_cut_sasha = new TH1F("xi__rec_right_sasha_pom","",11,xi_bins);
    xi_gen_right_cut_sasha = new TH1F("xi__gen_right_sasha_pom","",11,xi_bins);
    xi_rec_left_cut_sasha = new TH1F("xi__rec_left_sasha_pom","",8, bin_sasha);
    TH1F* x_rec_gen = new TH1F("","",20,-0.4,0.4);
    TH1F* resolution = new TH1F("","",20,-0.2,0.2);

    TH2F* pt2_xi_gen_minus_pomwig_rp= new TH2F("pt2_xi_gen_minus_rp_pomwig","",20,0,1,20,0,0.1);
    TH2F* pt2_xi_gen_minus_pomwig= new TH2F("pt2_xi_gen_minus_rp_pomwig","",20,0,1,20,0,0.1);
    TH2F* t_minus_th2 = new TH2F("t_minus_th2","", 8, tbins,8,tbins);
    TH2F* xi_minus_th2 = new TH2F("xi_minus_th2","",11,xi_bins,11,xi_bins);
    TH2F* logx_minus_th2 = new TH2F("logx_minus_th2","", 15, -4, 0,15, -4, 0);

    TH1F* corr_xi_cms_minus =  new TH1F("corr_xi_cms_minus","", 20, 0, 2);
    
    RooUnfoldResponse t_minus_response (t_rec_right_cut, t_gen_right_cut, "t_minus_unfolded", "t_minus_unfolded");
    RooUnfoldResponse t_minus_response_old (t_rec_right_cut, t_gen_right_cut, "t_minus_unfolded", "t_minus_unfolded");
    RooUnfoldResponse t_plus_response (t_rec_left_cut, t_gen_left_cut, "t_plus_unfolded", "t_plus_unfolded");
    RooUnfoldResponse xi_minus_response (xi_rec_right_cut, xi_gen_right, "xi_minus_unfolded", "xi_minus_unfolded");
    RooUnfoldResponse xi_plus_response (xi_rec_left_cut, xi_gen_left, "xi_plus_unfolded", "xi_plus_unfolded");
    RooUnfoldResponse logx_minus_response (log_x_rec_right_cut, log_x_gen_right_cut, "logx_minus_unfolded", "logx_minus_unfolded");
    RooUnfoldResponse logx_plus_response (log_x_rec_left_cut, log_x_gen_left_cut, "logx_plus_unfolded", "logx_plus_unfolded");
    RooUnfoldResponse t_minus_response_back (t_gen_right_cut, t_rec_right_cut, "t_minus_unfolded", "t_minus_unfolded");
    RooUnfoldResponse xi_minus_response_back (xi_gen_right, xi_rec_right_cut, "xi_minus_unfolded", "xi_minus_unfolded");
    RooUnfoldResponse t_plus_response_back (t_rec_left_cut, t_gen_left_cut, "t_plus_unfolded", "t_plus_unfolded");
    RooUnfoldResponse xi_plus_response_back (xi_gen_left, xi_rec_left_cut, "xi_plus_unfolded", "xi_plus_unfolded");
    RooUnfoldResponse logx_minus_response_back (log_x_gen_right_cut, log_x_rec_right_cut, "logx_minus_unfolded", "logx_minus_unfolded");
    

    TTree* tree_pom = (TTree*) pomwig_pom->Get( treeName.c_str() );
    int nev_pom = int(tree_pom->GetEntriesFast());
    cout <<"The pomwig-pomeron file has " << nev_pom << " entries : " << endl;
 
    double jet1_rec_pt_pom, jet1_rec_eta_pom, jet1_rec_phi_pom, jet2_rec_pt_pom, jet2_rec_eta_pom, jet2_rec_phi_pom;
    double jet1_gen_pt_pom, jet1_gen_eta_pom, jet1_gen_phi_pom, jet2_gen_pt_pom, jet2_gen_eta_pom, jet2_gen_phi_pom;
    double xi_rec_cms_minus_pom, xi_rec_proton_right_pom, x_rec_right_pom, x_gen_right_pom;
    double xi_gen_cms_minus_pom, xi_gen_proton_right_pom, xi_rec_proton_right_gauss_pom, xi_rec_proton_right_theta_pom;
    double t_rec_proton_right_pom, t_rec_proton_right_theta_pom, t_rec_proton_right_gauss_pom, beta_rec_proton_right_pom;
    double t_gen_proton_right_pom, beta_gen_proton_right_pom, rp_xpos_124_pom, rp_xpos_125_pom, rp_ypos_124_pom, rp_ypos_125_pom, t_rec_proton_right_old_pom;
    bool rp_right_pom, rp_right_accep_top_pom, rp_right_accep_bottom_pom;
    double theta_x_plus_smear_pom, theta_x_minus_smear_pom, theta_y_plus_smear_pom, theta_y_minus_smear_pom, pz_minus_smear_pom, pz_plus_smear_pom, e_minus_smear_pom, e_plus_smear_pom;
    double px_minus, py_minus, pz_minus, e_minus, mass_minus;
    tree_pom->SetBranchAddress("xi_rec_cms_right",&xi_rec_cms_minus_pom);
    tree_pom->SetBranchAddress("xi_gen_cms_right",&xi_gen_cms_minus_pom);
    tree_pom->SetBranchAddress("xi_rec_proton_right",&xi_rec_proton_right_pom);
    tree_pom->SetBranchAddress("xi_gen_proton_right",&xi_gen_proton_right_pom);
    tree_pom->SetBranchAddress("t_rec_proton_right",&t_rec_proton_right_pom);
    tree_pom->SetBranchAddress("t_gen_proton_right",&t_gen_proton_right_pom);
    tree_pom->SetBranchAddress("beta_rec_right",&beta_rec_proton_right_pom);
    tree_pom->SetBranchAddress("beta_gen_right",&beta_gen_proton_right_pom);
    tree_pom->SetBranchAddress("rp_right",&rp_right_pom);
    tree_pom->SetBranchAddress("rp_right_accep_bottom",&rp_right_accep_bottom_pom);
    tree_pom->SetBranchAddress("rp_right_accep_top",&rp_right_accep_top_pom);
    tree_pom->SetBranchAddress("rp_xpos_124",&rp_xpos_124_pom);
    tree_pom->SetBranchAddress("rp_ypos_124",&rp_ypos_124_pom);
    tree_pom->SetBranchAddress("rp_xpos_125",&rp_xpos_125_pom);
    tree_pom->SetBranchAddress("rp_ypos_125",&rp_ypos_125_pom);
    tree_pom->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt_pom);
    tree_pom->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt_pom);
    tree_pom->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta_pom);
    tree_pom->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta_pom);
    tree_pom->SetBranchAddress("jet1_rec_phi",&jet1_rec_phi_pom);
    tree_pom->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt_pom);
    tree_pom->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt_pom);
    tree_pom->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta_pom);
    tree_pom->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta_pom);
    tree_pom->SetBranchAddress("jet2_rec_phi",&jet2_rec_phi_pom);
    tree_pom->SetBranchAddress("x_rec_right",&x_rec_right_pom);
    tree_pom->SetBranchAddress("x_gen_right",&x_gen_right_pom);
    tree_pom->SetBranchAddress("xi_rec_proton_right_gauss",&xi_rec_proton_right_gauss_pom);
    tree_pom->SetBranchAddress("t_rec_proton_right_gauss",&t_rec_proton_right_gauss_pom);
    // tree_pom->SetBranchAddress("theta_x_minus_smear",&theta_x_minus_smear_pom);
    // tree_pom->SetBranchAddress("theta_y_minus_smear",&theta_y_minus_smear_pom);
    // tree_pom->SetBranchAddress("pz_proton_right_smear",&pz_minus_smear_pom);
    // tree_pom->SetBranchAddress("e_proton_right_smear",&e_minus_smear_pom);
    tree_pom->SetBranchAddress("px_proton_right",&px_minus);
    tree_pom->SetBranchAddress("py_proton_right",&py_minus);
    tree_pom->SetBranchAddress("pz_proton_right",&pz_minus);
    tree_pom->SetBranchAddress("e_proton_right",&e_minus);
    tree_pom->SetBranchAddress("mass_proton_right",&mass_minus);

    TH1F* p_right_pomwig = new TH1F("pomwig","",100,0,6000);
    TH1F* p_right_pythia = new TH1F("pythia","",100,0,6000);

    t_gen_right_cut->Sumw2();
    t_gen_left_cut->Sumw2();
    t_rec_right_cut->Sumw2();
    t_rec_left_cut->Sumw2();

    double pi = 4000;
    TH1F* t_old = new TH1F("","",8, tbins);

    for(int i_evt = 0; i_evt < nev_pom; ++i_evt){
        tree_pom->GetEntry(i_evt);
         
        // double px_minus_smear = -pi*tan(theta_x_minus_smear_pom);
        // double py_minus_smear = pi*tan(theta_y_minus_smear_pom);
        // TLorentzVector p_beam_minus_pom (0, 0, -pi, pi);
        // TLorentzVector p_scatt_minus_pom (px_minus_smear, py_minus_smear, pz_minus_smear_pom, e_minus_smear_pom);
        // TLorentzVector t_vec_minus_pom = (p_beam_minus_pom - p_scatt_minus_pom);
        // t_rec_proton_right_theta_pom = t_vec_minus_pom.Mag2(); 

        xi_rec_proton_right_pom = (unc_gauss == false) ? xi_rec_proton_right_pom : xi_rec_proton_right_gauss_pom;
        t_rec_proton_right_pom = (unc_gauss == false) ? t_rec_proton_right_pom : t_rec_proton_right_gauss_pom; 

        double reweigth_beta_pom_minus = (beta_gen_proton_right_pom<=0.7) ? func->Eval(beta_gen_proton_right_pom) : 1;
	    double reweight_slope_pom_minus = 0.95*t_slope_right_data*exp(-t_slope_right_data*fabs(t_gen_proton_right_pom))/(t_slope_right*exp(-t_slope_right*fabs(t_gen_proton_right_pom)));//0.95 normalization factor
        double event_weight_pom_minus;
	    if (reweight && !reweight_slope) event_weight_pom_minus = scale_rew*reweigth_beta_pom_minus; //0.84 normalization factor
        if (reweight && reweight_slope) event_weight_pom_minus = scale_rew*reweigth_beta_pom_minus*reweight_slope_pom_minus;
	    if (!reweight && reweight_slope) event_weight_pom_minus = reweight_slope_pom_minus;
	    if (!reweight && !reweight_slope) event_weight_pom_minus = 1;

        bool jet_rec_sel_pom = jet1_rec_pt_pom>pt_threshold && jet2_rec_pt_pom>pt_threshold && fabs(jet1_rec_eta_pom)<4.4 && fabs(jet2_rec_eta_pom)<4.4;
        bool jet_gen_sel_pom = jet1_gen_pt_pom>pt_threshold && jet2_gen_pt_pom>pt_threshold && fabs(jet1_gen_eta_pom)<4.4 && fabs(jet2_gen_eta_pom)<4.4;
        bool proton_rec_sel_pom =  xi_rec_proton_right_pom>0 && xi_rec_proton_right_pom<0.1 && fabs(t_rec_proton_right_pom)>0.03 && fabs(t_rec_proton_right_pom)<1;
        bool proton_gen_sel_pom =  xi_gen_proton_right_pom>0 && xi_gen_proton_right_pom<0.1 && fabs(t_gen_proton_right_pom)>0.03 && fabs(t_gen_proton_right_pom)<1;

        bool fid_cuts_nom_right_top_pom = rp_xpos_124_pom>0 && rp_xpos_124_pom<0.007 && rp_ypos_124_pom >0.0084 && rp_ypos_124_pom<0.029 ;
        bool fid_cuts_nom_right_bottom_pom = rp_xpos_125_pom>0 && rp_xpos_125_pom<0.007 && rp_ypos_125_pom <-0.0084 && rp_ypos_125_pom>-0.029 ;
        bool fid_cuts_unc_right_top_pom = rp_xpos_124_pom>0 && rp_xpos_124_pom<0.007 && rp_ypos_124_pom >0.0082 && rp_ypos_124_pom<0.027 ;
        bool fid_cuts_unc_right_bottom_pom = rp_xpos_125_pom>0 && rp_xpos_125_pom<0.007 && rp_ypos_125_pom <-0.0082 && rp_ypos_125_pom>-0.027 ;
        bool rp_right_unc_pom = (rp_right_accep_top_pom && fid_cuts_unc_right_top_pom) || (rp_right_accep_bottom_pom && fid_cuts_unc_right_bottom_pom);

        bool fid_cuts_unc_right_top_xmax_pom = rp_xpos_124_pom>-0.0005 && rp_xpos_124_pom<0.008 && rp_ypos_124_pom >0.0084 && rp_ypos_124_pom<0.029 ;
        bool fid_cuts_unc_right_bottom_xmax_pom = rp_xpos_125_pom>-0.0005 && rp_xpos_125_pom<0.008 && rp_ypos_125_pom <-0.0084 && rp_ypos_125_pom>-0.029 ;
        bool fid_cuts_unc_right_top_xmin_pom = rp_xpos_124_pom>0.00 && rp_xpos_124_pom<0.006 && rp_ypos_124_pom >0.0084 && rp_ypos_124_pom<0.029 ;
        bool fid_cuts_unc_right_bottom_xmin_pom = rp_xpos_125_pom>0.00 && rp_xpos_125_pom<0.006 && rp_ypos_125_pom <-0.0084 && rp_ypos_125_pom>-0.029 ;

        // bool rp_right_sel_pom = (unc_rp == false) ? rp_right_pom : rp_right_unc_pom;
        bool rp_right_sel_pom = false; 
        if(unc_rp == false && rp_x_min == false && rp_x_max == false) rp_right_sel_pom = (rp_right_accep_top_pom && fid_cuts_nom_right_top_pom) || (rp_right_accep_bottom_pom && fid_cuts_nom_right_bottom_pom);
        if(unc_rp == true && rp_x_min == false && rp_x_max == false) rp_right_sel_pom = rp_right_unc_pom ;
        if(unc_rp == false && rp_x_min == true && rp_x_max == false) rp_right_sel_pom = (rp_right_accep_top_pom && fid_cuts_unc_right_top_xmin_pom) || (rp_right_accep_bottom_pom && fid_cuts_unc_right_bottom_xmin_pom);
        if(unc_rp == false && rp_x_min == false && rp_x_max == true) rp_right_sel_pom = (rp_right_accep_top_pom && fid_cuts_unc_right_top_xmax_pom) || (rp_right_accep_bottom_pom && fid_cuts_unc_right_bottom_xmax_pom);

        if ((rp_right_accep_top_pom && fid_cuts_unc_right_top_pom) && !(rp_right_accep_bottom_pom && fid_cuts_unc_right_bottom_pom)) rp_pos_rigth->Fill(rp_xpos_124_pom*1000, rp_ypos_124_pom*1000, 1.);
        if (!(rp_right_accep_top_pom && fid_cuts_unc_right_top_pom) && (rp_right_accep_bottom_pom && fid_cuts_unc_right_bottom_pom)) rp_pos_rigth->Fill(rp_xpos_125_pom*1000, rp_ypos_125_pom*1000, 1.);

        if ((rp_right_accep_top_pom && fid_cuts_unc_right_top_xmax_pom) && !(rp_right_accep_bottom_pom && fid_cuts_unc_right_bottom_xmax_pom)) rp_pos_rigth_xmax->Fill(rp_xpos_124_pom*1000, rp_ypos_124_pom*1000, 1.);
        if (!(rp_right_accep_top_pom && fid_cuts_unc_right_top_xmax_pom) && (rp_right_accep_bottom_pom && fid_cuts_unc_right_bottom_xmax_pom)) rp_pos_rigth_xmax->Fill(rp_xpos_125_pom*1000, rp_ypos_125_pom*1000, 1.);

        // if ( rp_right_sel_pom ) rp_pos_rigth->Fill(rp_xpos_124_pom*1000, rp_ypos_124_pom*1000, 1.);
        // if ( rp_right_sel_pom ) rp_pos_rigth->Fill(rp_xpos_125_pom*1000, rp_ypos_125_pom*1000, 1.);

        corr_xi_cms_minus->Fill(xi_rec_cms_minus_pom/xi_gen_cms_minus_pom, 1.);
        if(jet_gen_sel_pom)pt2_xi_gen_minus_pomwig->Fill( fabs(t_gen_proton_right_pom)*(1 - xi_gen_proton_right_pom), xi_gen_proton_right_pom, 1);

        if (rp_right_pom && jet_gen_sel_pom){
           pt2_xi_gen_minus_pomwig_rp->Fill( fabs(t_gen_proton_right_pom)*(1 - xi_gen_proton_right_pom), xi_gen_proton_right_pom, 1);
   	    }
   	
        if (jet_rec_sel_pom){
           xi_cms_rec_right_sasha->Fill(xi_rec_cms_minus_pom,1);
           xi_cms_rec_right->Fill(xi_rec_cms_minus_pom, 1);
        }

        if (jet1_gen_pt_pom>20 && jet2_gen_pt_pom>20 && fabs(jet1_gen_eta_pom)<4.4 && fabs(jet2_gen_eta_pom)<4.4){
            xi_cms_gen_right_sasha->Fill(xi_gen_cms_minus_pom,1.);
        }

        if (jet_gen_sel_pom){
            xi_cms_gen_right->Fill(xi_gen_cms_minus_pom, event_weight_pom_minus);
            pt2_xi_gen_minus->Fill( fabs(t_gen_proton_right_pom)/*(1 + xi_gen_proton_right_pom)*/, xi_gen_proton_right_pom, event_weight_pom_minus);
            if (rp_right_pom) pt2_xi_gen_minus_rp->Fill( fabs(t_gen_proton_right_pom)/**(1 + xi_gen_proton_right_pom)*/, xi_gen_proton_right_pom, event_weight_pom_minus);
        }   
         
        // if (jet_rec_sel_pom && xi_rec_proton_right_pom>0 && xi_rec_proton_right_pom<0.1 && fabs(t_rec_proton_right_old_pom)>0.03 && fabs(t_rec_proton_right_old_pom)<1 && rp_right_sel_pom && xi_rec_cms_minus_pom - xi_rec_proton_right_pom<0){
        //         if (!(jet_gen_sel_pom && proton_gen_sel_pom))
        //            t_minus_response_old.Fake(fabs(t_rec_proton_right_old_pom), event_weight_pom_minus*scale_right*norm_pom);
        //         if (jet_gen_sel_pom && proton_gen_sel_pom)
        //            t_minus_response_old.Fill(fabs(t_rec_proton_right_old_pom), fabs(t_gen_proton_right_pom), event_weight_pom_minus*norm_pom*scale_right);
        // }
        // if (jet_gen_sel_pom && proton_gen_sel_pom && !(jet_rec_sel_pom && xi_rec_proton_right_pom>0 && xi_rec_proton_right_pom<0.1 && fabs(t_rec_proton_right_old_pom)>0.03 && fabs(t_rec_proton_right_old_pom)<1 && rp_right_sel_pom && xi_rec_cms_minus_pom - xi_rec_proton_right_pom<0))
        //            t_minus_response_old.Miss(fabs(t_gen_proton_right_pom), event_weight_pom_minus*scale_right*norm_pom);

        if (jet_rec_sel_pom && proton_rec_sel_pom && rp_right_sel_pom){
            xi_cms_minus_totem_rec_right->Fill(xi_rec_cms_minus_pom - xi_rec_proton_right_pom, event_weight_pom_minus);
            xi_cms_minus_totem_rec_right_bin->Fill(xi_rec_cms_minus_pom - xi_rec_proton_right_pom, event_weight_pom_minus);
            xi_rec_right->Fill(xi_rec_proton_right_pom, event_weight_pom_minus);
            if (xi_rec_cms_minus_pom - xi_rec_proton_right_pom<0){
                xi_rec_right_cut->Fill(xi_rec_proton_right_pom, event_weight_pom_minus);
                xi_rec_right_cut_sasha->Fill(xi_rec_proton_right_pom, event_weight_pom_minus);
                t_rec_right_cut->Fill(fabs(t_rec_proton_right_pom), event_weight_pom_minus); 
                beta_right_cut->Fill(beta_rec_proton_right_pom, event_weight_pom_minus);
                log_x_rec_right_cut->Fill(log10(x_rec_right_pom), event_weight_pom_minus);
		        pt_jet1_right_cut->Fill(jet1_rec_pt_pom, event_weight_pom_minus);
		        pt_jet2_right_cut->Fill(jet2_rec_pt_pom, event_weight_pom_minus);
		        eta_jet1_right_cut->Fill(jet1_rec_eta_pom, event_weight_pom_minus);
		        eta_jet2_right_cut->Fill(jet2_rec_eta_pom, event_weight_pom_minus);
                x_rec_gen->Fill(log10(x_rec_right_pom)-log10(x_gen_right_pom),1);
                if (!(jet_gen_sel_pom && proton_gen_sel_pom)){
                   t_minus_response.Fake(fabs(t_rec_proton_right_pom), event_weight_pom_minus*scale_right*norm_pom);
                   xi_minus_response.Fake(xi_rec_proton_right_pom, event_weight_pom_minus*scale_right*norm_pom);
                   logx_minus_response.Fake(log10(x_rec_right_pom), event_weight_pom_minus*scale_right*norm_pom);
                   t_minus_response_back.Miss(fabs(t_rec_proton_right_pom), event_weight_pom_minus*scale_right*norm_pom);
                   xi_minus_response_back.Miss(xi_rec_proton_right_pom, event_weight_pom_minus*scale_right*norm_pom);
                   logx_minus_response_back.Miss(log10(x_rec_right_pom), event_weight_pom_minus*scale_right*norm_pom);
                 }  
                if (jet_gen_sel_pom && proton_gen_sel_pom){
                   t_minus_response.Fill(fabs(t_rec_proton_right_pom), fabs(t_gen_proton_right_pom), event_weight_pom_minus*norm_pom*scale_right);
                   t_minus_response_back.Fill(fabs(t_gen_proton_right_pom), fabs(t_rec_proton_right_pom), event_weight_pom_minus*norm_pom*scale_right);
                   xi_minus_response.Fill(xi_rec_proton_right_pom, xi_gen_proton_right_pom, event_weight_pom_minus*norm_pom*scale_right);
                   xi_minus_response_back.Fill(xi_gen_proton_right_pom, xi_rec_proton_right_pom, event_weight_pom_minus*norm_pom*scale_right);
                   logx_minus_response.Fill(log10(x_rec_right_pom),log10(x_gen_right_pom), event_weight_pom_minus*norm_pom*scale_right);
                   logx_minus_response_back.Fill(log10(x_gen_right_pom),log10(x_rec_right_pom), event_weight_pom_minus*norm_pom*scale_right);
                   xi_minus_th2->Fill(xi_rec_proton_right_pom, xi_gen_proton_right_pom, event_weight_pom_minus*norm_pom*scale_right);
                   logx_minus_th2->Fill(log10(x_rec_right_pom),log10(x_gen_right_pom), event_weight_pom_minus*norm_pom*scale_right);
                   t_minus_th2->Fill(fabs(t_rec_proton_right_pom), fabs(t_gen_proton_right_pom), event_weight_pom_minus*norm_pom*scale_right);
		        }  
            }
        }

        if (jet_gen_sel_pom && proton_gen_sel_pom){
           t_gen_right_cut->Fill(fabs(t_gen_proton_right_pom), event_weight_pom_minus); 
           xi_gen_right->Fill(xi_gen_proton_right_pom, event_weight_pom_minus);
           xi_gen_right_cut_sasha->Fill(xi_gen_proton_right_pom, event_weight_pom_minus);
           log_x_gen_right_cut->Fill(log10(x_gen_right_pom), event_weight_pom_minus);
           if (!(jet_rec_sel_pom && proton_rec_sel_pom && rp_right_sel_pom && xi_rec_cms_minus_pom - xi_rec_proton_right_pom<0)){
              t_minus_response.Miss(fabs(t_gen_proton_right_pom), event_weight_pom_minus*norm_pom);
              t_minus_response_back.Fake(fabs(t_gen_proton_right_pom), event_weight_pom_minus*norm_pom);
              xi_minus_response.Miss(xi_gen_proton_right_pom, event_weight_pom_minus*norm_pom);
              xi_minus_response_back.Fake(xi_gen_proton_right_pom, event_weight_pom_minus*norm_pom);
              logx_minus_response.Miss(log10(x_gen_right_pom), event_weight_pom_minus*norm_pom);
              logx_minus_response_back.Fake(log10(x_gen_right_pom), event_weight_pom_minus*norm_pom);
           }   
        } 

        // if(jet1_gen_pt_pom>20 && jet2_gen_pt_pom>20 && fabs(jet1_gen_eta_pom)<4.4 && fabs(jet2_gen_eta_pom)<4.4 && proton_gen_sel_pom) xi_gen_right_sasha->Fill(xi_gen_proton_right_pom, event_weight_pom_minus);


    }

    TTree* tree_pom_plus = (TTree*) pomwig_pom_plus->Get( treeName.c_str() );
    int nev_pom_plus = int(tree_pom_plus->GetEntriesFast());
    cout <<"The pomwig-pomeron-plus file has " << nev_pom_plus << " entries : " << endl;
 
    double jet1_rec_pt_pom_plus, jet1_rec_eta_pom_plus, jet1_rec_phi_pom_plus, jet2_rec_pt_pom_plus, jet2_rec_eta_pom_plus, jet2_rec_phi_pom_plus;
    double jet1_gen_pt_pom_plus, jet1_gen_eta_pom_plus, jet1_gen_phi_pom_plus, jet2_gen_pt_pom_plus, jet2_gen_eta_pom_plus, jet2_gen_phi_pom_plus;
    double xi_rec_cms_plus_pom, xi_rec_proton_left_pom, x_rec_left_pom, x_gen_left_pom, xi_rec_proton_left_theta_pom, xi_rec_proton_left_gauss_pom;
    double xi_gen_cms_plus_pom, xi_gen_proton_left_pom;
    double t_rec_proton_left_pom, t_rec_proton_left_gauss_pom, t_rec_proton_left_theta_pom, beta_rec_proton_left_pom,t_rec_proton_left_old_pom;
    double t_gen_proton_left_pom, beta_gen_proton_left_pom, rp_xpos_24_pom, rp_xpos_25_pom, rp_ypos_24_pom, rp_ypos_25_pom;
    double px_plus, py_plus, pz_plus, mass_plus, e_plus;
    bool rp_left_pom, rp_left_accep_bottom_pom, rp_left_accep_top_pom;
    tree_pom_plus->SetBranchAddress("xi_rec_cms_left",&xi_rec_cms_plus_pom);
    tree_pom_plus->SetBranchAddress("xi_gen_cms_left",&xi_gen_cms_plus_pom);
    tree_pom_plus->SetBranchAddress("xi_rec_proton_left",&xi_rec_proton_left_pom);
    tree_pom_plus->SetBranchAddress("xi_gen_proton_left",&xi_gen_proton_left_pom);
    tree_pom_plus->SetBranchAddress("t_rec_proton_left",&t_rec_proton_left_pom);
    tree_pom_plus->SetBranchAddress("t_gen_proton_left",&t_gen_proton_left_pom);
    tree_pom_plus->SetBranchAddress("beta_rec_left",&beta_rec_proton_left_pom);
    tree_pom_plus->SetBranchAddress("beta_gen_left",&beta_gen_proton_left_pom);
    tree_pom_plus->SetBranchAddress("rp_left",&rp_left_pom);
    tree_pom_plus->SetBranchAddress("rp_left_accep_bottom",&rp_left_accep_bottom_pom);
    tree_pom_plus->SetBranchAddress("rp_left_accep_top",&rp_left_accep_top_pom);
    tree_pom_plus->SetBranchAddress("rp_xpos_24",&rp_xpos_24_pom);
    tree_pom_plus->SetBranchAddress("rp_ypos_24",&rp_ypos_24_pom);
    tree_pom_plus->SetBranchAddress("rp_xpos_25",&rp_xpos_25_pom);
    tree_pom_plus->SetBranchAddress("rp_ypos_25",&rp_ypos_25_pom);
    tree_pom_plus->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt_pom_plus);
    tree_pom_plus->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta_pom_plus);
    tree_pom_plus->SetBranchAddress("jet1_rec_phi",&jet1_rec_phi_pom_plus);
    tree_pom_plus->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt_pom_plus);
    tree_pom_plus->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta_pom_plus);
    tree_pom_plus->SetBranchAddress("jet2_rec_phi",&jet2_rec_phi_pom_plus);
    tree_pom_plus->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt_pom_plus);
    tree_pom_plus->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta_pom_plus);
    tree_pom_plus->SetBranchAddress("jet1_gen_phi",&jet1_gen_phi_pom_plus);
    tree_pom_plus->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt_pom_plus);
    tree_pom_plus->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta_pom_plus);
    tree_pom_plus->SetBranchAddress("jet2_gen_phi",&jet2_gen_phi_pom_plus);
    tree_pom_plus->SetBranchAddress("x_rec_left",&x_rec_left_pom);
    tree_pom_plus->SetBranchAddress("x_gen_left",&x_gen_left_pom);
    tree_pom_plus->SetBranchAddress("xi_rec_proton_left_gauss",&xi_rec_proton_left_gauss_pom);
    tree_pom_plus->SetBranchAddress("t_rec_proton_left_gauss",&t_rec_proton_left_gauss_pom);
    // tree_pom_plus->SetBranchAddress("theta_x_plus_smear",&theta_x_plus_smear_pom);
    // tree_pom_plus->SetBranchAddress("theta_y_plus_smear",&theta_y_plus_smear_pom);
    // tree_pom_plus->SetBranchAddress("pz_proton_left_smear",&pz_plus_smear_pom);
    // tree_pom_plus->SetBranchAddress("e_proton_left_smear",&e_plus_smear_pom);
    tree_pom_plus->SetBranchAddress("px_proton_left",&px_plus);
    tree_pom_plus->SetBranchAddress("py_proton_left",&py_plus);
    tree_pom_plus->SetBranchAddress("pz_proton_left",&pz_plus);
    tree_pom_plus->SetBranchAddress("e_proton_left",&e_plus);
    tree_pom_plus->SetBranchAddress("mass_proton_left",&mass_plus);
    TH1F* t_gen_left_pt20 = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH1F* t_gen_left_pt20_nomass = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH1F* t_gen_left_pt30 = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH1F* t_gen_left_pt40 = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH1F* t_gen_left_pt50 = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH2F* log_t_vs_xi = new TH2F("","",100,-3.,1,100,-3.5,1);
    TH2F* log_t_vs_xi_nomass = new TH2F("","",100,-3.,1,100,-3.5,1);

    TH1F* t_old_left = new TH1F("","",8,tbins);
    TH1F* t_old_left_regg = new TH1F("","",8,tbins);

    for(int i_evt = 0; i_evt < nev_pom_plus; ++i_evt){
        tree_pom_plus->GetEntry(i_evt);

        // double px_plus_smear = -pi*tan(theta_x_plus_smear_pom);
        // double py_plus_smear = pi*tan(theta_y_plus_smear_pom);
        // TLorentzVector p_beam_plus_pom (0, 0, pi, pi);
        // TLorentzVector p_scatt_plus_pom (px_plus_smear, py_plus_smear, pz_plus_smear_pom, e_plus_smear_pom);
        // TLorentzVector t_vec_plus_pom = (p_beam_plus_pom - p_scatt_plus_pom);
        // t_rec_proton_left_theta_pom = t_vec_plus_pom.Mag2(); 

        xi_rec_proton_left_pom = (unc_gauss == false) ? xi_rec_proton_left_pom : xi_rec_proton_left_gauss_pom;
        t_rec_proton_left_pom = (unc_gauss == false) ? t_rec_proton_left_pom : t_rec_proton_left_gauss_pom;

        double reweigth_beta_pom_plus = (beta_gen_proton_left_pom<=0.7) ? func->Eval(beta_gen_proton_left_pom) : 1;
	    double reweight_slope_pom_plus = 0.95*t_slope_left_data*exp(-t_slope_left_data*fabs(t_gen_proton_left_pom))/(t_slope_left*exp(-t_slope_left*fabs(t_gen_proton_left_pom)));
        double event_weight_pom_plus;
	    if (reweight && !reweight_slope) event_weight_pom_plus = scale_rew*reweigth_beta_pom_plus;
        if (reweight && reweight_slope) event_weight_pom_plus = scale_rew*reweigth_beta_pom_plus*reweight_slope_pom_plus;
	    if (!reweight && reweight_slope) event_weight_pom_plus = reweight_slope_pom_plus;
	    if (!reweight && !reweight_slope) event_weight_pom_plus = 1;

        bool jet_rec_sel_pom_plus = jet1_rec_pt_pom_plus>pt_threshold && jet2_rec_pt_pom_plus>pt_threshold && fabs(jet1_rec_eta_pom_plus)<4.4 && fabs(jet2_rec_eta_pom_plus)<4.4;
        bool jet_gen_sel_pom_plus = jet1_gen_pt_pom_plus>pt_threshold && jet2_gen_pt_pom_plus>pt_threshold && fabs(jet1_gen_eta_pom_plus)<4.4 && fabs(jet2_gen_eta_pom_plus)<4.4;
        bool proton_rec_sel_pom_plus =  xi_rec_proton_left_pom>0 && xi_rec_proton_left_pom<0.1 && fabs(t_rec_proton_left_pom)>0.03 && fabs(t_rec_proton_left_pom)<1;
        bool proton_gen_sel_pom_plus = xi_gen_proton_left_pom>0 && xi_gen_proton_left_pom<0.1 && fabs(t_gen_proton_left_pom)>0.03 && fabs(t_gen_proton_left_pom)<1;

        bool fid_cuts_nom_left_top_pom = rp_xpos_24_pom>0 && rp_xpos_24_pom<0.007 && rp_ypos_24_pom >0.0084 && rp_ypos_24_pom<0.027 ;
        bool fid_cuts_nom_left_bottom_pom = rp_xpos_25_pom>0 && rp_xpos_25_pom<0.007 && rp_ypos_25_pom <-0.0084 && rp_ypos_25_pom>-0.027 ;
        bool fid_cuts_unc_left_top_pom = rp_xpos_24_pom>0 && rp_xpos_24_pom<0.007 && rp_ypos_24_pom >0.0082 && rp_ypos_24_pom<0.029 ;
        bool fid_cuts_unc_left_bottom_pom = rp_xpos_25_pom>0 && rp_xpos_25_pom<0.007 && rp_ypos_25_pom <-0.0082 && rp_ypos_25_pom>-0.029 ;
        bool fid_cuts_unc_left_top_xmin_pom = rp_xpos_24_pom>0.000 && rp_xpos_24_pom<0.006 && rp_ypos_24_pom >0.0084 && rp_ypos_24_pom<0.027 ;
        bool fid_cuts_unc_left_top_xmax_pom = rp_xpos_24_pom>-0.0005 && rp_xpos_24_pom<0.008 && rp_ypos_24_pom >0.0084 && rp_ypos_24_pom<0.027 ;
        bool fid_cuts_unc_left_bottom_xmin_pom = rp_xpos_25_pom>0.000 && rp_xpos_25_pom<0.006 && rp_ypos_25_pom <-0.0084 && rp_ypos_25_pom>-0.027 ;
        bool fid_cuts_unc_left_bottom_xmax_pom = rp_xpos_25_pom>-0.0005 && rp_xpos_25_pom<0.008 && rp_ypos_25_pom <-0.0084 && rp_ypos_25_pom>-0.027 ;

        bool rp_left_unc_pom = (rp_left_accep_top_pom && fid_cuts_unc_left_top_pom) || (rp_left_accep_bottom_pom && fid_cuts_unc_left_bottom_pom);
        // bool rp_left_sel_pom = (unc_rp == false) ? rp_left_pom : rp_left_unc_pom;
        bool rp_left_sel_pom = false; 
        if(unc_rp == false && rp_x_min == false && rp_x_max == false) rp_left_sel_pom = (rp_left_accep_top_pom && fid_cuts_nom_left_top_pom) || (rp_left_accep_bottom_pom && fid_cuts_nom_left_bottom_pom);
        if(unc_rp == true && rp_x_min == false && rp_x_max == false) rp_left_sel_pom = rp_left_unc_pom ;
        if(unc_rp == false && rp_x_min == true && rp_x_max == false) rp_left_sel_pom = (rp_left_accep_top_pom && fid_cuts_unc_left_top_xmin_pom) || (rp_left_accep_bottom_pom && fid_cuts_unc_left_bottom_xmin_pom);
        if(unc_rp == false && rp_x_min == false && rp_x_max == true) rp_left_sel_pom = (rp_left_accep_top_pom && fid_cuts_unc_left_top_xmax_pom) || (rp_left_accep_bottom_pom && fid_cuts_unc_left_bottom_xmax_pom);

 
    if ((rp_left_accep_top_pom && fid_cuts_unc_left_top_pom) && !(rp_left_accep_bottom_pom && fid_cuts_unc_left_bottom_pom)) rp_pos_left->Fill(rp_xpos_24_pom*1000, rp_ypos_24_pom*1000, 1.);
    if (!(rp_left_accep_top_pom && fid_cuts_unc_left_top_pom) && (rp_left_accep_bottom_pom && fid_cuts_unc_left_bottom_pom)) rp_pos_left->Fill(rp_xpos_25_pom*1000, rp_ypos_25_pom*1000, 1.);
	
    if ((rp_left_accep_top_pom && fid_cuts_unc_left_top_xmax_pom) && !(rp_left_accep_bottom_pom && fid_cuts_unc_left_bottom_xmax_pom)) rp_pos_left_xmax->Fill(rp_xpos_24_pom*1000, rp_ypos_24_pom*1000, 1.);
    if (!(rp_left_accep_top_pom && fid_cuts_unc_left_top_xmax_pom) && (rp_left_accep_bottom_pom && fid_cuts_unc_left_bottom_xmax_pom)) rp_pos_left_xmax->Fill(rp_xpos_25_pom*1000, rp_ypos_25_pom*1000, 1.);

        if (jet_rec_sel_pom_plus){
           xi_cms_rec_left->Fill(xi_rec_cms_plus_pom,1.);
           xi_cms_rec_left_sasha->Fill(xi_rec_cms_plus_pom,1.);
        }
        if (jet1_gen_pt_pom_plus>20 && jet2_gen_pt_pom_plus>20 && fabs(jet1_gen_eta_pom_plus)<4.4 && fabs(jet2_gen_eta_pom_plus)<4.4){
            xi_cms_gen_left->Fill(xi_gen_cms_plus_pom,1.);
            xi_cms_gen_left_sasha->Fill(xi_gen_cms_plus_pom,1.);
        }   
        if (rp_left_pom && jet_gen_sel_pom_plus){
   	}
    
        if (jet_rec_sel_pom_plus && xi_rec_proton_left_pom>0 && xi_rec_proton_left_pom<0.1 && fabs(t_rec_proton_left_old_pom)>0.03 && fabs(t_rec_proton_left_old_pom)<1 && rp_left_sel_pom && xi_rec_cms_plus_pom - xi_rec_proton_left_pom<0)
            t_old_left->Fill(fabs(t_rec_proton_left_old_pom), event_weight_pom_plus);

        if (jet_rec_sel_pom_plus && proton_rec_sel_pom_plus && rp_left_sel_pom){
            xi_cms_minus_totem_rec_left->Fill(xi_rec_cms_plus_pom - xi_rec_proton_left_pom, event_weight_pom_plus);
            xi_cms_minus_totem_rec_left_bin->Fill(xi_rec_cms_plus_pom - xi_rec_proton_left_pom, event_weight_pom_plus);
            xi_rec_left->Fill(xi_rec_proton_left_pom, event_weight_pom_plus);
            if (xi_rec_cms_plus_pom - xi_rec_proton_left_pom<0){
                xi_rec_left_cut->Fill(xi_rec_proton_left_pom, event_weight_pom_plus);
 // resolution->Fill(xi_rec_proton_left_pom - xi_gen_proton_left_pom, 1);
 // resolution->Fill(fabs(t_rec_proton_left_pom) - fabs(t_gen_proton_left_pom), 1);
 // resolution->Fill(log10(x_rec_left_pom) - log10(x_gen_left_pom), 1);
                xi_rec_left_cut_sasha->Fill(xi_rec_proton_left_pom, event_weight_pom_plus);
                t_rec_left_cut->Fill(fabs(t_rec_proton_left_pom), event_weight_pom_plus);
                beta_left_cut->Fill(beta_rec_proton_left_pom, event_weight_pom_plus);
                log_x_rec_left_cut->Fill(log10(x_rec_left_pom), event_weight_pom_plus);
	            pt_jet1_left_cut->Fill(jet1_rec_pt_pom_plus, event_weight_pom_plus);
	            pt_jet2_left_cut->Fill(jet2_rec_pt_pom_plus, event_weight_pom_plus);
	            eta_jet1_left_cut->Fill(jet1_rec_eta_pom_plus, event_weight_pom_plus);
	            eta_jet2_left_cut->Fill(jet2_rec_eta_pom_plus, event_weight_pom_plus);
                if (!(jet_gen_sel_pom_plus && proton_gen_sel_pom_plus)){
                   t_plus_response.Fake(fabs(t_rec_proton_left_pom), event_weight_pom_plus*norm_pom_plus);
                   xi_plus_response.Fake(xi_rec_proton_left_pom, event_weight_pom_plus*norm_pom_plus);
                   logx_plus_response.Fake(log10(x_rec_left_pom), event_weight_pom_plus*norm_pom_plus);
                   t_plus_response_back.Miss(fabs(t_rec_proton_left_pom), event_weight_pom_plus*norm_pom_plus);
                   xi_plus_response_back.Miss(xi_rec_proton_left_pom, event_weight_pom_plus*norm_pom_plus);
                }   
		        if (jet_gen_sel_pom_plus && proton_gen_sel_pom_plus){
                   t_plus_response.Fill(fabs(t_rec_proton_left_pom), fabs(t_gen_proton_left_pom), event_weight_pom_plus*norm_pom_plus);
                   xi_plus_response.Fill(xi_rec_proton_left_pom, xi_gen_proton_left_pom, event_weight_pom_plus*norm_pom_plus);
                   logx_plus_response.Fill(log10(x_rec_left_pom),log10(x_gen_left_pom), event_weight_pom_plus*norm_pom_plus);
                   t_plus_response_back.Fill(fabs(t_gen_proton_left_pom), fabs(t_rec_proton_left_pom), event_weight_pom_plus*norm_pom_plus);
                   xi_plus_response_back.Fill(xi_gen_proton_left_pom, xi_rec_proton_left_pom, event_weight_pom_plus*norm_pom_plus);
		        }   
            }
        }
        if (jet_gen_sel_pom_plus && proton_gen_sel_pom_plus){
           t_gen_left_cut->Fill(fabs(t_gen_proton_left_pom), event_weight_pom_plus);
           xi_gen_left->Fill(xi_gen_proton_left_pom, event_weight_pom_plus);
           log_x_gen_left_cut->Fill(log10(x_gen_left_pom), event_weight_pom_plus);
           if (!(jet_rec_sel_pom_plus && proton_rec_sel_pom_plus && rp_left_sel_pom && xi_rec_cms_plus_pom - xi_rec_proton_left_pom<0)){
               t_plus_response.Miss(fabs(t_gen_proton_left_pom), event_weight_pom_plus*norm_pom_plus);
              xi_plus_response.Miss(xi_gen_proton_left_pom, event_weight_pom_plus*norm_pom_plus);
              logx_plus_response.Miss(log10(x_gen_left_pom), event_weight_pom_plus*norm_pom_plus);}
               t_plus_response_back.Fake(fabs(t_gen_proton_left_pom), event_weight_pom_plus*norm_pom_plus);
              xi_plus_response_back.Fake(xi_gen_proton_left_pom, event_weight_pom_plus*norm_pom_plus);
        }

        // if(jet1_gen_pt_pom_plus>20 && jet2_gen_pt_pom_plus>20 && fabs(jet1_gen_eta_pom_plus)<4.4 && fabs(jet2_gen_eta_pom_plus)<4.4 && proton_gen_sel_pom_plus) xi_gen_left_sasha->Fill(xi_gen_proton_left_pom, event_weight_pom_plus);


        // double pi_mass = sqrt(pi*pi - M_P*M_P);
        // TLorentzVector p_beam_plus (0, 0, pi, pi);
        // double e_gen_plus_from_p_nomass = sqrt(px_plus*px_plus + py_plus*py_plus + pz_plus*pz_plus);// + M_P*M_P);
        // double e_gen_plus_from_p = sqrt(px_plus*px_plus + py_plus*py_plus + pz_plus*pz_plus + M_P*M_P);
        // TLorentzVector p_scatt_plus_nomass (px_plus, py_plus, pz_plus, e_gen_plus_from_p_nomass);
        // TLorentzVector p_scatt_plus (px_plus, py_plus, pz_plus, e_gen_plus_from_p);
        // TLorentzVector t_vec_plus_nomass = (p_beam_plus - p_scatt_plus_nomass);
        // TLorentzVector t_vec_plus = (p_beam_plus - p_scatt_plus);
        // double t_rec_proton_left_theta_nomass = t_vec_plus_nomass.Mag2(); 
        // double t_rec_proton_left_theta = t_vec_plus.Mag2(); 
        // /*if(jet1_gen_pt_pom_plus>20 && fabs(jet1_gen_eta_pom_plus)<5)*/ t_gen_left_pt20_nomass->Fill(fabs(t_rec_proton_left_theta_nomass), event_weight_pom_plus);
        // /*if(jet1_gen_pt_pom_plus>20 && fabs(jet1_gen_eta_pom_plus)<5)*/ t_gen_left_pt20->Fill(fabs(t_rec_proton_left_theta), event_weight_pom_plus);
        log_t_vs_xi->Fill(log10(fabs(t_rec_proton_left_pom)), log10(xi_gen_proton_left_pom),1.);
        // log_t_vs_xi->Fill(log10(fabs(t_rec_proton_left_theta)), log10(xi_gen_proton_left_pom),1.);
        // log_t_vs_xi_nomass->Fill(log10(fabs(t_rec_proton_left_theta_nomass)), log10(xi_gen_proton_left_pom),1.);
        // if(jet1_gen_pt_pom_plus>30 && fabs(jet1_gen_eta_pom_plus)<5) t_gen_left_pt30->Fill(fabs(t_rec_proton_left_theta), event_weight_pom_plus);
        // if(jet1_gen_pt_pom_plus>40 && fabs(jet1_gen_eta_pom_plus)<5) t_gen_left_pt40->Fill(fabs(t_rec_proton_left_theta), event_weight_pom_plus);
        // if(jet1_gen_pt_pom_plus>50 && fabs(jet1_gen_eta_pom_plus)<5) t_gen_left_pt50->Fill(fabs(t_rec_proton_left_theta), event_weight_pom_plus);


    }

    //reggeon contribution
    reggeon_histos["xi_cms_minus_totem_rec_right"] = new TH1F("xi_cms_minus_totem_rec_right_regg","",50,-0.4,0.4);
    reggeon_histos["xi_cms_minus_totem_rec_right_bin"] = new TH1F("xi_cms_minus_totem_rec_right_regg_bin","", 15, bin);
    reggeon_histos["xi_cms_minus_totem_rec_left"] = new TH1F("xi_cms_minus_totem_rec_left_regg","",50,-0.4,0.4);
    reggeon_histos["xi_cms_minus_totem_rec_left_bin"] = new TH1F("xi_cms_minus_totem_rec_left_regg_bin","", 15, bin);
    reggeon_histos["xi_rec_right"] = new TH1F("xi_rec_right_regg","",50,-0.04,0.2);
    reggeon_histos["xi_gen_right"] = new TH1F("xi_gen_right_regg","",50,-0.04,0.2);
    reggeon_histos["xi_rec_left"] = new TH1F("xi_rec_left_regg","",50,-0.04,0.2);
    reggeon_histos["xi_gen_left"] = new TH1F("xi_gen_left_regg","",50,-0.04,0.2);
    reggeon_histos["xi_rec_right_cut"] = new TH1F("xi_rec_right_cut_regg","",11,xi_bins);
    reggeon_histos["xi_rec_right_cut_sasha"] = new TH1F("xi_rec_right_cut_regg","",11,xi_bins);
    reggeon_histos["xi_rec_left_cut_sasha"] = new TH1F("xi_rec_left_cut_regg","",8, bin_sasha);
    reggeon_histos["xi_gen_right_cut"] = new TH1F("xi_gen_right_cut_regg","",11,xi_bins);
    reggeon_histos["xi_gen_right_cut_sasha"] = new TH1F("xi_gen_right_cut_regg","",11,xi_bins);
    reggeon_histos["xi_rec_left_cut"] = new TH1F("xi_rec_left_cut_regg","",11,xi_bins);
    reggeon_histos["xi_gen_left_cut"] = new TH1F("xi_gen_left_cut_regg","",11,xi_bins);
    reggeon_histos["xi_rec_right_response"] = new TH1F("xi_rec_right_response_regg","",11,xi_bins);
    reggeon_histos["xi_rec_left_response"] = new TH1F("xi_rec_left_response_regg","",11,xi_bins);
    reggeon_histos["xi_gen_right_response"] = new TH1F("xi_gen_right_response_regg","",11,xi_bins);
    reggeon_histos["xi_gen_left_response"] = new TH1F("xi_gen_left_response_regg","",11,xi_bins);
    reggeon_histos["t_rec_right_cut"] = new TH1F("t_rec_right_cut_regg","", 8, tbins);
    reggeon_histos["t_rec_right_response"] = new TH1F("t_rec_right_response_regg","", 8, tbins);
    reggeon_histos["t_gen_right_cut"] = new TH1F("t_gen_right_cut_regg","", 8, tbins);
    reggeon_histos["t_gen_right_response"] = new TH1F("t_gen_right_response_regg","", 8, tbins);
    reggeon_histos["t_rec_left_cut"] = new TH1F("t_rec_left_cut_regg","", 8, tbins);
    reggeon_histos["t_rec_left_response"] = new TH1F("t_rec_left_response_regg","", 8, tbins);
    reggeon_histos["t_gen_left_cut"] = new TH1F("t_gen_left_cut_regg","", 8, tbins);
    reggeon_histos["t_gen_left_response"] = new TH1F("t_gen_left_response_regg","", 8, tbins);
    reggeon_histos["beta_right_cut"] = new TH1F("beta_right_cut_regg","",15,0,1);
    reggeon_histos["beta_left_cut"] = new TH1F("beta_left_cut_regg","",15,0,1);
    reggeon_histos["xi_cms_rec_right"] = new TH1F("xi_cms_rec_right_regg","",5,0,0.1);
    reggeon_histos["xi_cms_gen_right"] = new TH1F("xi_cms_gen_right_regg","",5,0,0.1);
    reggeon_histos["xi_cms_rec_left"] = new TH1F("xi_cms_rec_left_regg","",11,xi_bins);
    reggeon_histos["xi_cms_gen_left"] = new TH1F("xi_cms_gen_left_regg","",11,xi_bins);
    reggeon_histos["xi_cms_rec_right_sasha"] = new TH1F("xi_cms_rec_right_sasha_regg","",8, bin_sasha);
    reggeon_histos["xi_cms_gen_right_sasha"] = new TH1F("xi_cms_gen_right_sasha_regg","",8, bin_sasha);
    reggeon_histos["xi_cms_rec_left_sasha"] = new TH1F("xi_cms_rec_left_sasha_regg","",8, bin_sasha);
    reggeon_histos["xi_cms_gen_left_sasha"] = new TH1F("xi_cms_gen_left_sasha_regg","",8, bin_sasha);
    reggeon_histos["log_x_rec_right_cut"] = new TH1F("log_x_rec_right_cut_regg","",15, -4, 0);
    reggeon_histos["log_x_gen_right_cut"] = new TH1F("log_x_gen_right_cut_regg","",15, -4, 0);
    reggeon_histos["log_x_rec_left_cut"] = new TH1F("log_x_rec_left_cut_regg","",15, -4, 0);
    reggeon_histos["log_x_gen_left_cut"] = new TH1F("log_x_gen_left_cut_regg","",15, -4, 0);
    reggeon_histos["log_x_rec_right_response"] = new TH1F("log_x_rec_right_response_regg","",15, -4, 0);
    reggeon_histos["log_x_rec_left_response"] = new TH1F("log_x_rec_left_response_regg","",15, -4, 0);
    reggeon_histos["log_x_gen_right_response"] = new TH1F("log_x_gen_right_response_regg","",15, -4, 0);
    reggeon_histos["log_x_gen_left_response"] = new TH1F("log_x_gen_left_response_regg","",15, -4, 0);
    reggeon_histos["pt_jet1_right_cut"] = new TH1F("pt_jet1_right_cut_regg","", 15, 0, 200);
    reggeon_histos["pt_jet2_right_cut"] = new TH1F("pt_jet2_right_cut_regg","", 15, 0, 200);
    reggeon_histos["eta_jet1_right_cut"] = new TH1F("eta_jet1_right_cut_regg","", 20, -5.2, 5.2);
    reggeon_histos["eta_jet2_right_cut"] = new TH1F("eta_jet2_right_cut_regg","", 20, -5.2, 5.2);
    reggeon_histos["pt_jet1_left_cut"] = new TH1F("pt_jet1_left_cut_regg","", 15, 0, 200);
    reggeon_histos["pt_jet2_left_cut"] = new TH1F("pt_jet2_left_cut_regg","", 15, 0, 200);
    reggeon_histos["eta_jet1_left_cut"] = new TH1F("eta_jet1_left_cut_regg","", 20, -5.2, 5.2);
    reggeon_histos["eta_jet2_left_cut"] = new TH1F("eta_jet2_left_cut_regg","", 20, -5.2, 5.2);
    // reggeon_histos["xi_gen_right_sasha"] = new TH1F("xi_gen_right_sasha_regg","",8, bin_sasha);
    // reggeon_histos["xi_gen_left_sasha"] = new TH1F("xi_gen_left_sasha_regg","",8, bin_sasha);
 
    TTree* tree_regg = (TTree*) pomwig_regg->Get( treeName.c_str() );
    int nev_regg = int(tree_regg->GetEntriesFast());
    cout <<"The pomwig-reggeon file has " << nev_regg << " entries : " << endl;
   reggeon_histos["t_gen_right_cut"]->Sumw2();
   reggeon_histos["t_gen_left_cut"]->Sumw2();
   reggeon_histos["t_rec_right_cut"]->Sumw2();
   reggeon_histos["t_rec_left_cut"]->Sumw2();

    double jet1_rec_pt_regg, jet1_rec_eta_regg, jet1_rec_phi_regg, jet2_rec_pt_regg, jet2_rec_eta_regg, jet2_rec_phi_regg;
    double jet1_gen_pt_regg, jet1_gen_eta_regg, jet1_gen_phi_regg, jet2_gen_pt_regg, jet2_gen_eta_regg, jet2_gen_phi_regg;
    double xi_rec_cms_minus_regg, xi_rec_proton_right_regg, xi_rec_proton_right_gauss_regg, xi_rec_proton_right_theta_regg, x_rec_right_regg, x_gen_right_regg;
    double xi_gen_cms_minus_regg, xi_gen_proton_right_regg;
    double t_rec_proton_right_regg, t_rec_proton_right_theta_regg, t_rec_proton_right_gauss_regg, beta_rec_proton_right_regg, rp_xpos_124_regg, rp_xpos_125_regg, rp_ypos_124_regg, rp_ypos_125_regg;
    double t_gen_proton_right_regg, beta_gen_proton_right_regg, t_rec_proton_right_old_regg;
    bool rp_right_regg, rp_right_accep_bottom_regg, rp_right_accep_top_regg;
    double theta_x_plus_smear_regg, theta_x_minus_smear_regg, theta_y_plus_smear_regg, theta_y_minus_smear_regg, pz_minus_smear_regg, pz_plus_smear_regg, e_minus_smear_regg, e_plus_smear_regg;
    tree_regg->SetBranchAddress("xi_rec_cms_right",&xi_rec_cms_minus_regg);
    tree_regg->SetBranchAddress("xi_gen_cms_right",&xi_gen_cms_minus_regg);
    tree_regg->SetBranchAddress("xi_rec_proton_right",&xi_rec_proton_right_regg);
    tree_regg->SetBranchAddress("xi_gen_proton_right",&xi_gen_proton_right_regg);
    tree_regg->SetBranchAddress("t_rec_proton_right",&t_rec_proton_right_regg);
    tree_regg->SetBranchAddress("t_gen_proton_right",&t_gen_proton_right_regg);
    tree_regg->SetBranchAddress("beta_rec_right",&beta_rec_proton_right_regg);
    tree_regg->SetBranchAddress("beta_gen_right",&beta_gen_proton_right_regg);
    tree_regg->SetBranchAddress("rp_right",&rp_right_regg);
    tree_regg->SetBranchAddress("rp_right_accep_bottom",&rp_right_accep_bottom_regg);
    tree_regg->SetBranchAddress("rp_right_accep_top",&rp_right_accep_top_regg);
    tree_regg->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt_regg);
    tree_regg->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta_regg);
    tree_regg->SetBranchAddress("jet1_rec_phi",&jet1_rec_phi_regg);
    tree_regg->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt_regg);
    tree_regg->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta_regg);
    tree_regg->SetBranchAddress("jet2_rec_phi",&jet2_rec_phi_regg);
    tree_regg->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt_regg);
    tree_regg->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta_regg);
    tree_regg->SetBranchAddress("jet1_gen_phi",&jet1_gen_phi_regg);
    tree_regg->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt_regg);
    tree_regg->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta_regg);
    tree_regg->SetBranchAddress("jet2_gen_phi",&jet2_gen_phi_regg);
    tree_regg->SetBranchAddress("x_rec_right",&x_rec_right_regg);
    tree_regg->SetBranchAddress("x_gen_right",&x_gen_right_regg);
    tree_regg->SetBranchAddress("rp_xpos_124",&rp_xpos_124_regg);
    tree_regg->SetBranchAddress("rp_ypos_124",&rp_ypos_124_regg);
    tree_regg->SetBranchAddress("rp_xpos_125",&rp_xpos_125_regg);
    tree_regg->SetBranchAddress("rp_ypos_125",&rp_ypos_125_regg);
    tree_regg->SetBranchAddress("xi_rec_proton_right_gauss",&xi_rec_proton_right_gauss_regg);
    tree_regg->SetBranchAddress("t_rec_proton_right_gauss",&t_rec_proton_right_gauss_regg);
    // tree_regg->SetBranchAddress("theta_x_minus_smear",&theta_x_minus_smear_regg);
    // tree_regg->SetBranchAddress("theta_y_minus_smear",&theta_y_minus_smear_regg);
    // tree_regg->SetBranchAddress("pz_proton_right_smear",&pz_minus_smear_regg);
    // tree_regg->SetBranchAddress("e_proton_right_smear",&e_minus_smear_regg);

    TH1F* t_old_regg = new TH1F("","", 8, tbins);

    for(int i_evt = 0; i_evt < nev_regg; ++i_evt){
        tree_regg->GetEntry(i_evt);

        // double px_minus_smear_regg = -pi*tan(theta_x_minus_smear_regg);
        // double py_minus_smear_regg = pi*tan(theta_y_minus_smear_regg);
        // TLorentzVector p_beam_minus_regg (0, 0, -pi, pi);
        // TLorentzVector p_scatt_minus_regg (px_minus_smear_regg, py_minus_smear_regg, pz_minus_smear_regg, e_minus_smear_regg);
        // TLorentzVector t_vec_minus_regg = (p_beam_minus_regg - p_scatt_minus_regg);
        // t_rec_proton_right_theta_regg = t_vec_minus_regg.Mag2(); 

        xi_rec_proton_right_regg = (unc_gauss == false) ? xi_rec_proton_right_regg : xi_rec_proton_right_gauss_regg;
        t_rec_proton_right_regg = (unc_gauss == false) ? t_rec_proton_right_regg : t_rec_proton_right_gauss_regg;
      
        double reweigth_beta_regg_minus = (beta_gen_proton_right_regg<=0.7) ? func->Eval(beta_gen_proton_right_regg) : 1;
   	    double reweight_slope_regg_minus = 0.95*t_slope_right_data*exp(-t_slope_right_data*fabs(t_gen_proton_right_regg))/(t_slope_right*exp(-t_slope_right*fabs(t_gen_proton_right_regg)));
        double event_weight_regg_minus;
	    if (reweight && !reweight_slope) event_weight_regg_minus = scale_rew*reweigth_beta_regg_minus;
        if (reweight && reweight_slope) event_weight_regg_minus = scale_rew*reweigth_beta_regg_minus*reweight_slope_regg_minus;
	    if (!reweight && reweight_slope) event_weight_regg_minus = reweight_slope_regg_minus;
	   if (!reweight && !reweight_slope) event_weight_regg_minus = 1;

        bool jet_rec_sel_regg = jet1_rec_pt_regg>pt_threshold && jet2_rec_pt_regg>pt_threshold && fabs(jet1_rec_eta_regg)<4.4 && fabs(jet2_rec_eta_regg)<4.4;
        bool jet_gen_sel_regg = jet1_gen_pt_regg>pt_threshold && jet2_gen_pt_regg>pt_threshold && fabs(jet1_gen_eta_regg)<4.4 && fabs(jet2_gen_eta_regg)<4.4;
        bool proton_rec_sel_regg =  xi_rec_proton_right_regg>0 && xi_rec_proton_right_regg<0.1 && fabs(t_rec_proton_right_regg)>0.03 && fabs(t_rec_proton_right_regg)<1;
        bool proton_gen_sel_regg =  xi_gen_proton_right_regg>0 && xi_gen_proton_right_regg<0.1 && fabs(t_gen_proton_right_regg)>0.03 && fabs(t_gen_proton_right_regg)<1;

        bool fid_cuts_nom_right_top_regg = rp_xpos_124_regg>0 && rp_xpos_124_regg<0.007 && rp_ypos_124_regg >0.0084 && rp_ypos_124_regg<0.029 ;
        bool fid_cuts_nom_right_bottom_regg = rp_xpos_125_regg>0 && rp_xpos_125_regg<0.007 && rp_ypos_125_regg <-0.0084 && rp_ypos_125_regg>-0.029 ;
        bool fid_cuts_unc_right_top_regg = rp_xpos_124_regg>0 && rp_xpos_124_regg<0.007 && rp_ypos_124_regg >0.0082 && rp_ypos_124_regg<0.027 ;
        bool fid_cuts_unc_right_bottom_regg = rp_xpos_125_regg>0 && rp_xpos_125_regg<0.007 && rp_ypos_125_regg <-0.0082 && rp_ypos_125_regg>-0.027 ;
        bool rp_right_unc_regg = (rp_right_accep_top_regg && fid_cuts_unc_right_top_regg) || (rp_right_accep_bottom_regg && fid_cuts_unc_right_bottom_regg);

        bool fid_cuts_unc_right_top_xmax_regg = rp_xpos_124_regg>-0.0005 && rp_xpos_124_regg<0.008 && rp_ypos_124_regg >0.0084 && rp_ypos_124_regg<0.029 ;
        bool fid_cuts_unc_right_bottom_xmax_regg = rp_xpos_125_regg>-0.0005 && rp_xpos_125_regg<0.008 && rp_ypos_125_regg <-0.0084 && rp_ypos_125_regg>-0.029 ;
        bool fid_cuts_unc_right_top_xmin_regg = rp_xpos_124_regg>0.00 && rp_xpos_124_regg<0.006 && rp_ypos_124_regg >0.0084 && rp_ypos_124_regg<0.029 ;
        bool fid_cuts_unc_right_bottom_xmin_regg = rp_xpos_125_regg>0.00 && rp_xpos_125_regg<0.006 && rp_ypos_125_regg <-0.0084 && rp_ypos_125_regg>-0.029 ;

        // bool rp_right_sel_regg = (unc_rp == false) ? rp_right_regg : rp_right_unc_regg;

        bool rp_right_sel_regg = false; 
        if(unc_rp == false && rp_x_min == false && rp_x_max == false) rp_right_sel_regg = (rp_right_accep_top_regg && fid_cuts_nom_right_top_regg) || (rp_right_accep_bottom_regg && fid_cuts_nom_right_bottom_regg);
        if(unc_rp == true && rp_x_min == false && rp_x_max == false) rp_right_sel_regg = rp_right_unc_regg ;
        if(unc_rp == false && rp_x_min == true && rp_x_max == false) rp_right_sel_regg = (rp_right_accep_top_regg && fid_cuts_unc_right_top_xmin_regg) || (rp_right_accep_bottom_regg && fid_cuts_unc_right_bottom_xmin_regg);
        if(unc_rp == false && rp_x_min == false && rp_x_max == true) rp_right_sel_regg = (rp_right_accep_top_regg && fid_cuts_unc_right_top_xmax_regg) || (rp_right_accep_bottom_regg && fid_cuts_unc_right_bottom_xmax_regg);

        if (jet_rec_sel_regg){
           reggeon_histos["xi_cms_rec_right_sasha"]->Fill(xi_rec_cms_minus_regg,1);
           reggeon_histos["xi_cms_rec_right"]->Fill(xi_rec_cms_minus_regg, 1);
        }   
            
        if (jet_gen_sel_regg) reggeon_histos["xi_cms_gen_right"]->Fill(xi_gen_cms_minus_regg, event_weight_regg_minus);
        
        if (jet1_gen_pt_regg>20 && jet2_gen_pt_regg>20 && fabs(jet1_gen_eta_regg)<4.4 && fabs(jet2_gen_eta_regg)<4.4){
           reggeon_histos["xi_cms_gen_right_sasha"]->Fill(xi_gen_cms_minus_regg,1.);
           pt2_xi_gen_minus->Fill( fabs(t_gen_proton_right_regg)/*(1 + xi_gen_proton_right_pom)*/, xi_gen_proton_right_regg, event_weight_regg_minus);
           if (rp_right_regg) pt2_xi_gen_minus_rp->Fill( fabs(t_gen_proton_right_regg)/**(1 + xi_gen_proton_right_pom)*/, xi_gen_proton_right_regg, event_weight_regg_minus);
        }   
           
        // if (jet_rec_sel_regg && xi_rec_proton_right_regg>0 && xi_rec_proton_right_regg<0.1 && fabs(t_rec_proton_right_old_regg)>0.03 && fabs(t_rec_proton_right_old_regg)<1 && rp_right_sel_regg && xi_rec_cms_minus_regg - xi_rec_proton_right_regg<0){
        //  t_old_regg->Fill(fabs(t_rec_proton_right_old_regg), event_weight_regg_minus);
        //         if (!(jet_gen_sel_regg && proton_gen_sel_regg))
        //            t_minus_response_old.Fake(fabs(t_rec_proton_right_old_regg), event_weight_regg_minus*scale_right*norm_regg);
        // }
        // if (jet_gen_sel_regg && proton_gen_sel_regg && !(jet_rec_sel_regg && xi_rec_proton_right_regg>0 && xi_rec_proton_right_regg<0.1 && fabs(t_rec_proton_right_old_regg)>0.03 && fabs(t_rec_proton_right_old_regg)<1 && rp_right_sel_regg && xi_rec_cms_minus_regg - xi_rec_proton_right_regg<0))
        //       t_minus_response_old.Miss(fabs(t_gen_proton_right_regg), event_weight_regg_minus*norm_regg);

        if (jet_rec_sel_regg && proton_rec_sel_regg && rp_right_sel_regg){
            reggeon_histos["xi_cms_minus_totem_rec_right"]->Fill(xi_rec_cms_minus_regg - xi_rec_proton_right_regg, event_weight_regg_minus);
            reggeon_histos["xi_cms_minus_totem_rec_right_bin"]->Fill(xi_rec_cms_minus_regg - xi_rec_proton_right_regg, event_weight_regg_minus);
            reggeon_histos["xi_rec_right"]->Fill(xi_rec_proton_right_regg, event_weight_regg_minus);
            if (xi_rec_cms_minus_regg - xi_rec_proton_right_regg<0){
                reggeon_histos["xi_rec_right_cut"]->Fill(xi_rec_proton_right_regg, event_weight_regg_minus);
                reggeon_histos["xi_rec_right_cut_sasha"]->Fill(xi_rec_proton_right_regg, event_weight_regg_minus);
                reggeon_histos["t_rec_right_cut"]->Fill(fabs(t_rec_proton_right_regg), event_weight_regg_minus);
                reggeon_histos["beta_right_cut"]->Fill(beta_rec_proton_right_regg, event_weight_regg_minus);
                reggeon_histos["log_x_rec_right_cut"]->Fill(log10(x_rec_right_regg), event_weight_regg_minus);
                reggeon_histos["pt_jet1_right_cut"]->Fill(jet1_rec_pt_regg, event_weight_regg_minus);
                reggeon_histos["pt_jet2_right_cut"]->Fill(jet2_rec_pt_regg, event_weight_regg_minus);
                reggeon_histos["eta_jet1_right_cut"]->Fill(jet1_rec_eta_regg, event_weight_regg_minus);
                reggeon_histos["eta_jet2_right_cut"]->Fill(jet2_rec_eta_regg, event_weight_regg_minus);
                if (!(jet_gen_sel_regg && proton_gen_sel_regg)){
                   t_minus_response.Fake(fabs(t_rec_proton_right_regg), event_weight_regg_minus*scale_right*norm_regg);
                   t_minus_response_back.Miss(fabs(t_rec_proton_right_regg), event_weight_regg_minus*scale_right*norm_regg);
                   xi_minus_response.Fake(xi_rec_proton_right_regg, event_weight_regg_minus*scale_right*norm_regg);
                   xi_minus_response_back.Miss(xi_rec_proton_right_regg, event_weight_regg_minus*scale_right*norm_regg);
                   logx_minus_response.Fake(log10(x_rec_right_regg), event_weight_regg_minus*scale_right*norm_regg);
                   logx_minus_response_back.Miss(log10(x_rec_right_regg), event_weight_regg_minus*scale_right*norm_regg);
                 }  
		       if (jet_gen_sel_regg && proton_gen_sel_regg){
                   t_minus_response.Fill(fabs(t_rec_proton_right_regg), fabs(t_gen_proton_right_regg), event_weight_regg_minus*norm_regg*scale_right);
                   t_minus_response_back.Fill(fabs(t_gen_proton_right_regg), fabs(t_rec_proton_right_regg), event_weight_regg_minus*norm_regg*scale_right);
                   xi_minus_response.Fill(xi_rec_proton_right_regg, xi_gen_proton_right_regg, event_weight_regg_minus*norm_regg*scale_right);
                   xi_minus_response_back.Fill(xi_gen_proton_right_regg, xi_rec_proton_right_regg, event_weight_regg_minus*norm_regg*scale_right);
                   logx_minus_response.Fill(log10(x_rec_right_regg),log10(x_gen_right_regg), event_weight_regg_minus*norm_regg*scale_right);
                   logx_minus_response_back.Fill(log10(x_gen_right_regg),log10(x_rec_right_regg), event_weight_regg_minus*norm_regg*scale_right);
                   t_minus_th2->Fill(fabs(t_rec_proton_right_regg), fabs(t_gen_proton_right_regg), event_weight_regg_minus*norm_regg*scale_right);
                   xi_minus_th2->Fill(xi_rec_proton_right_regg, xi_gen_proton_right_regg, event_weight_regg_minus*norm_regg*scale_right);
                   logx_minus_th2->Fill(log10(x_rec_right_regg),log10(x_gen_right_regg), event_weight_regg_minus*norm_regg*scale_right);
		        }	  
	        }

        }
        if (jet_gen_sel_regg && proton_gen_sel_regg){
            reggeon_histos["t_gen_right_cut"]->Fill(fabs(t_gen_proton_right_regg), event_weight_regg_minus);
            reggeon_histos["xi_gen_right_cut"]->Fill(xi_gen_proton_right_regg, event_weight_regg_minus);
            reggeon_histos["xi_gen_right_cut_sasha"]->Fill(xi_gen_proton_right_regg, event_weight_regg_minus);
            reggeon_histos["log_x_gen_right_cut"]->Fill(log10(x_gen_right_regg), event_weight_regg_minus);
           if (!(jet_rec_sel_regg && proton_rec_sel_regg && rp_right_sel_regg && xi_rec_cms_minus_regg - xi_rec_proton_right_regg<0)){
              t_minus_response.Miss(fabs(t_gen_proton_right_regg), event_weight_regg_minus*norm_regg);
              t_minus_response_back.Fake(fabs(t_gen_proton_right_regg), event_weight_regg_minus*norm_regg);
              xi_minus_response.Miss(xi_gen_proton_right_regg, event_weight_regg_minus*norm_regg);
              xi_minus_response_back.Fake(xi_gen_proton_right_regg, event_weight_regg_minus*norm_regg);
              logx_minus_response.Miss(log10(x_gen_right_regg), event_weight_regg_minus*norm_regg);
              logx_minus_response_back.Fake(log10(x_gen_right_regg), event_weight_regg_minus*norm_regg);
           }   
        }    
        // if(jet1_gen_pt_regg>20 && jet2_gen_pt_regg>20 && fabs(jet1_gen_eta_regg)<4.4 && fabs(jet2_gen_eta_regg)<4.4 && proton_gen_sel_regg) reggeon_histos["xi_gen_right_sasha"]->Fill(xi_gen_proton_right_regg, event_weight_regg_minus);
    }

    TTree* tree_regg_plus = (TTree*) pomwig_regg_plus->Get( treeName.c_str() );
    int nev_regg_plus = int(tree_regg_plus->GetEntriesFast());
    cout <<"The pomwig-reggeon-plus file has " << nev_regg_plus << " entries : " << endl;
 
    double jet1_rec_pt_regg_plus, jet1_rec_eta_regg_plus, jet1_rec_phi_regg_plus, jet2_rec_pt_regg_plus, jet2_rec_eta_regg_plus, jet2_rec_phi_regg_plus;
    double jet1_gen_pt_regg_plus, jet1_gen_eta_regg_plus, jet1_gen_phi_regg_plus, jet2_gen_pt_regg_plus, jet2_gen_eta_regg_plus, jet2_gen_phi_regg_plus;
    double xi_rec_cms_plus_regg, xi_rec_proton_left_regg, xi_rec_proton_left_theta_regg, xi_rec_proton_left_gauss_regg, x_rec_left_regg, x_gen_left_regg;
    double xi_gen_cms_plus_regg, xi_gen_proton_left_regg, rp_xpos_24_regg_plus, rp_xpos_25_regg_plus, rp_ypos_24_regg_plus, rp_ypos_25_regg_plus;
    double t_rec_proton_left_regg, t_rec_proton_left_theta_regg, t_rec_proton_left_gauss_regg, beta_rec_proton_left_regg,t_rec_proton_left_old_regg;
    double t_gen_proton_left_regg, beta_gen_proton_left_regg;
    bool rp_left_regg, rp_left_accep_bottom_regg, rp_left_accep_top_regg;
    tree_regg_plus->SetBranchAddress("xi_rec_cms_left",&xi_rec_cms_plus_regg);
    tree_regg_plus->SetBranchAddress("xi_gen_cms_left",&xi_gen_cms_plus_regg);
    tree_regg_plus->SetBranchAddress("xi_rec_proton_left",&xi_rec_proton_left_regg);
    tree_regg_plus->SetBranchAddress("xi_gen_proton_left",&xi_gen_proton_left_regg);
    tree_regg_plus->SetBranchAddress("t_rec_proton_left",&t_rec_proton_left_regg);
    tree_regg_plus->SetBranchAddress("t_gen_proton_left",&t_gen_proton_left_regg);
    tree_regg_plus->SetBranchAddress("beta_rec_left",&beta_rec_proton_left_regg);
    tree_regg_plus->SetBranchAddress("beta_gen_left",&beta_gen_proton_left_regg);
    tree_regg_plus->SetBranchAddress("rp_left",&rp_left_regg);
    tree_regg_plus->SetBranchAddress("rp_left_accep_bottom",&rp_left_accep_bottom_regg);
    tree_regg_plus->SetBranchAddress("rp_left_accep_top",&rp_left_accep_top_regg);
    tree_regg_plus->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt_regg_plus);
    tree_regg_plus->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta_regg_plus);
    tree_regg_plus->SetBranchAddress("jet1_rec_phi",&jet1_rec_phi_regg_plus);
    tree_regg_plus->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt_regg_plus);
    tree_regg_plus->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta_regg_plus);
    tree_regg_plus->SetBranchAddress("jet2_rec_phi",&jet2_rec_phi_regg_plus);
    tree_regg_plus->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt_regg_plus);
    tree_regg_plus->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta_regg_plus);
    tree_regg_plus->SetBranchAddress("jet1_gen_phi",&jet1_gen_phi_regg_plus);
    tree_regg_plus->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt_regg_plus);
    tree_regg_plus->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta_regg_plus);
    tree_regg_plus->SetBranchAddress("jet2_gen_phi",&jet2_gen_phi_regg_plus);
    tree_regg_plus->SetBranchAddress("x_rec_left",&x_rec_left_regg);
    tree_regg_plus->SetBranchAddress("x_gen_left",&x_gen_left_regg);
    tree_regg_plus->SetBranchAddress("rp_xpos_24",&rp_xpos_24_regg_plus);
    tree_regg_plus->SetBranchAddress("rp_ypos_24",&rp_ypos_24_regg_plus);
    tree_regg_plus->SetBranchAddress("rp_xpos_25",&rp_xpos_25_regg_plus);
    tree_regg_plus->SetBranchAddress("rp_ypos_25",&rp_ypos_25_regg_plus);
    tree_regg_plus->SetBranchAddress("xi_rec_proton_left_gauss",&xi_rec_proton_left_gauss_regg);
    tree_regg_plus->SetBranchAddress("t_rec_proton_left_gauss",&t_rec_proton_left_gauss_regg);
    // tree_regg_plus->SetBranchAddress("theta_x_plus_smear",&theta_x_plus_smear_regg);
    // tree_regg_plus->SetBranchAddress("theta_y_plus_smear",&theta_y_plus_smear_regg);
    // tree_regg_plus->SetBranchAddress("pz_proton_left_smear",&pz_plus_smear_regg);
    // tree_regg_plus->SetBranchAddress("e_proton_left_smear",&e_plus_smear_regg);

    for(int i_evt = 0; i_evt < nev_regg_plus; ++i_evt){
        tree_regg_plus->GetEntry(i_evt);

        // double px_plus_smear_regg = -pi*tan(theta_x_plus_smear_regg);
        // double py_plus_smear_regg = pi*tan(theta_y_plus_smear_regg);
        // TLorentzVector p_beam_plus_regg (0, 0, pi, pi);
        // TLorentzVector p_scatt_plus_regg (px_plus_smear_regg, py_plus_smear_regg, pz_plus_smear_regg, e_plus_smear_regg);
        // TLorentzVector t_vec_plus_regg = (p_beam_plus_regg - p_scatt_plus_regg);
        // t_rec_proton_left_theta_regg = t_vec_plus_regg.Mag2(); 

        xi_rec_proton_left_regg = (unc_gauss == false) ? xi_rec_proton_left_regg : xi_rec_proton_left_gauss_regg;
        t_rec_proton_left_regg = (unc_gauss == false) ? t_rec_proton_left_regg : t_rec_proton_left_gauss_regg;

        double reweigth_beta_regg_plus = (beta_gen_proton_left_regg<=0.7) ? func->Eval(beta_gen_proton_left_regg) : 1;
  	    double reweight_slope_regg_plus = 0.95*t_slope_left_data*exp(-t_slope_left_data*fabs(t_gen_proton_left_regg))/(t_slope_left*exp(-t_slope_left*fabs(t_gen_proton_left_regg)));//0.95 normalization factor
        double event_weight_regg_plus;
	    if (reweight && !reweight_slope) event_weight_regg_plus = scale_rew*reweigth_beta_regg_plus;
        if (reweight && reweight_slope) event_weight_regg_plus = scale_rew*reweigth_beta_regg_plus*reweight_slope_regg_plus;
	    if (!reweight && reweight_slope) event_weight_regg_plus = reweight_slope_regg_plus;
	    if (!reweight && !reweight_slope) event_weight_regg_plus = 1;

        bool jet_rec_sel_regg_plus = jet1_rec_pt_regg_plus>pt_threshold && jet2_rec_pt_regg_plus>pt_threshold && fabs(jet1_rec_eta_regg_plus)<4.4 && fabs(jet2_rec_eta_regg_plus)<4.4;
        bool jet_gen_sel_regg_plus = jet1_gen_pt_regg_plus>pt_threshold && jet2_gen_pt_regg_plus>pt_threshold && fabs(jet1_gen_eta_regg_plus)<4.4 && fabs(jet2_gen_eta_regg_plus)<4.4;
        bool proton_rec_sel_regg_plus =  xi_rec_proton_left_regg>0 && xi_rec_proton_left_regg<0.1 && fabs(t_rec_proton_left_regg)>0.03 && fabs(t_rec_proton_left_regg)<1;
        bool proton_gen_sel_regg_plus = xi_gen_proton_left_regg>0 && xi_gen_proton_left_regg<0.1 && fabs(t_gen_proton_left_regg)>0.03 && fabs(t_gen_proton_left_regg)<1;

        bool fid_cuts_nom_left_top_regg = rp_xpos_24_regg_plus>0 && rp_xpos_24_regg_plus<0.007 && rp_ypos_24_regg_plus >0.0084 && rp_ypos_24_regg_plus<0.027 ;
        bool fid_cuts_nom_left_bottom_regg = rp_xpos_25_regg_plus>0 && rp_xpos_25_regg_plus<0.007 && rp_ypos_25_regg_plus <-0.0084 && rp_ypos_25_regg_plus>-0.027 ;
        bool fid_cuts_unc_left_top_regg = rp_xpos_24_regg_plus>0 && rp_xpos_24_regg_plus<0.007 && rp_ypos_24_regg_plus >0.0082 && rp_ypos_24_regg_plus<0.029 ;
        bool fid_cuts_unc_left_bottom_regg = rp_xpos_25_regg_plus>0 && rp_xpos_25_regg_plus<0.007 && rp_ypos_25_regg_plus <-0.0082 && rp_ypos_25_regg_plus>-0.029 ;
        bool fid_cuts_unc_left_top_xmin_regg = rp_xpos_24_regg_plus>0.00 && rp_xpos_24_regg_plus<0.006 && rp_ypos_24_regg_plus >0.0084 && rp_ypos_24_regg_plus<0.027 ;
        bool fid_cuts_unc_left_top_xmax_regg = rp_xpos_24_regg_plus>-0.0005 && rp_xpos_24_regg_plus<0.008 && rp_ypos_24_regg_plus >0.0084 && rp_ypos_24_regg_plus<0.027 ;
        bool fid_cuts_unc_left_bottom_xmin_regg = rp_xpos_25_regg_plus>0.00 && rp_xpos_25_regg_plus<0.006 && rp_ypos_25_regg_plus <-0.0084 && rp_ypos_25_regg_plus>-0.027 ;
        bool fid_cuts_unc_left_bottom_xmax_regg = rp_xpos_25_regg_plus>-0.0005 && rp_xpos_25_regg_plus<0.008 && rp_ypos_25_regg_plus <-0.0084 && rp_ypos_25_regg_plus>-0.027 ;
        bool rp_left_unc_regg = (rp_left_accep_top_regg && fid_cuts_unc_left_top_regg) || (rp_left_accep_bottom_regg && fid_cuts_unc_left_bottom_regg);
        // bool rp_left_sel_regg = (unc_rp == false) ? rp_left_regg : rp_left_unc_regg;
	
        bool rp_left_sel_regg = false; 
        if(unc_rp == false && rp_x_min == false && rp_x_max == false) rp_left_sel_regg = (rp_left_accep_top_regg && fid_cuts_nom_left_top_regg) || (rp_left_accep_bottom_regg && fid_cuts_nom_left_bottom_regg);
        if(unc_rp == true && rp_x_min == false && rp_x_max == false) rp_left_sel_regg = rp_left_unc_regg ;
        if(unc_rp == false && rp_x_min == true && rp_x_max == false) rp_left_sel_regg = (rp_left_accep_top_regg && fid_cuts_unc_left_top_xmin_regg) || (rp_left_accep_bottom_regg && fid_cuts_unc_left_bottom_xmin_regg);
        if(unc_rp == false && rp_x_min == false && rp_x_max == true) rp_left_sel_regg = (rp_left_accep_top_regg && fid_cuts_unc_left_top_xmax_regg) || (rp_left_accep_bottom_regg && fid_cuts_unc_left_bottom_xmax_regg);

         if (jet_rec_sel_regg_plus){
           reggeon_histos["xi_cms_rec_left"]->Fill(xi_rec_cms_plus_regg,1.);
           reggeon_histos["xi_cms_rec_left_sasha"]->Fill(xi_rec_cms_plus_regg,1.);
        }   
        if (jet1_gen_pt_regg_plus>20 && jet2_gen_pt_regg_plus>20 && fabs(jet1_gen_eta_regg_plus)<4.4 && fabs(jet2_gen_eta_regg_plus)<4.4){
           reggeon_histos["xi_cms_gen_left"]->Fill(xi_gen_cms_plus_regg,1.);
           reggeon_histos["xi_cms_gen_left_sasha"]->Fill(xi_gen_cms_plus_regg,1.);
        }   

        if (jet_rec_sel_regg_plus && xi_rec_proton_left_regg>0 && xi_rec_proton_left_regg<0.1 && fabs(t_rec_proton_left_old_regg)>0.03 && fabs(t_rec_proton_left_old_regg)<1 && rp_left_sel_regg && xi_rec_cms_plus_regg - xi_rec_proton_left_regg<0)
            t_old_left_regg->Fill(fabs(t_rec_proton_left_old_regg), event_weight_regg_plus);

       if (jet_rec_sel_regg_plus && proton_rec_sel_regg_plus && rp_left_sel_regg){
            reggeon_histos["xi_cms_minus_totem_rec_left"]->Fill(xi_rec_cms_plus_regg - xi_rec_proton_left_regg, event_weight_regg_plus);
            reggeon_histos["xi_cms_minus_totem_rec_left_bin"]->Fill(xi_rec_cms_plus_regg - xi_rec_proton_left_regg, event_weight_regg_plus);
            reggeon_histos["xi_rec_left"]->Fill(xi_rec_proton_left_regg, event_weight_regg_plus);
            if (xi_rec_cms_plus_regg - xi_rec_proton_left_regg<0){
                reggeon_histos["t_rec_left_cut"]->Fill(fabs(t_rec_proton_left_regg), event_weight_regg_plus);
                reggeon_histos["xi_rec_left_cut"]->Fill(xi_rec_proton_left_regg, event_weight_regg_plus);
                reggeon_histos["xi_rec_left_cut_sasha"]->Fill(xi_rec_proton_left_regg, event_weight_regg_plus);
                reggeon_histos["beta_left_cut"]->Fill(beta_rec_proton_left_regg, event_weight_regg_plus);
                reggeon_histos["log_x_rec_left_cut"]->Fill(log10(x_rec_left_regg), event_weight_regg_plus);
                reggeon_histos["pt_jet1_left_cut"]->Fill(jet1_rec_pt_regg_plus, event_weight_regg_plus);
                reggeon_histos["pt_jet2_left_cut"]->Fill(jet2_rec_pt_regg_plus, event_weight_regg_plus);
                reggeon_histos["eta_jet1_left_cut"]->Fill(jet1_rec_eta_regg_plus, event_weight_regg_plus);
                reggeon_histos["eta_jet2_left_cut"]->Fill(jet2_rec_eta_regg_plus, event_weight_regg_plus);
                if (!(jet_gen_sel_regg_plus && proton_gen_sel_regg_plus)){
                   t_plus_response.Fake(fabs(t_rec_proton_left_regg), event_weight_regg_plus*norm_regg_plus);
                   xi_plus_response.Fake(xi_rec_proton_left_regg, event_weight_regg_plus*norm_regg_plus);
                   logx_plus_response.Fake(log10(x_rec_left_regg), event_weight_regg_plus*norm_regg_plus);
                   t_plus_response_back.Miss(fabs(t_rec_proton_left_regg), event_weight_regg_plus*norm_regg_plus);
                   xi_plus_response_back.Miss(xi_rec_proton_left_regg, event_weight_regg_plus*norm_regg_plus);
                }   
		        if (jet_gen_sel_regg_plus && proton_gen_sel_regg_plus){
                   t_plus_response.Fill(fabs(t_rec_proton_left_regg), fabs(t_gen_proton_left_regg), event_weight_regg_plus*norm_regg_plus);
                   xi_plus_response.Fill(xi_rec_proton_left_regg, xi_gen_proton_left_regg, event_weight_regg_plus*norm_regg_plus);
                   logx_plus_response.Fill(log10(x_rec_left_regg),log10(x_gen_left_regg), event_weight_regg_plus*norm_regg_plus);
                   t_plus_response_back.Fill(fabs(t_gen_proton_left_regg), fabs(t_rec_proton_left_regg), event_weight_regg_plus*norm_regg_plus);
                   xi_plus_response_back.Fill(xi_gen_proton_left_regg, xi_rec_proton_left_regg, event_weight_regg_plus*norm_regg_plus);
		        }   
            }
        }
        if (jet_gen_sel_regg_plus && proton_gen_sel_regg_plus){
           reggeon_histos["t_gen_left_cut"]->Fill(fabs(t_gen_proton_left_regg), event_weight_regg_plus);
           reggeon_histos["xi_gen_left_cut"]->Fill(xi_gen_proton_left_regg, event_weight_regg_plus);
           reggeon_histos["log_x_gen_left_cut"]->Fill(log10(x_gen_left_regg), event_weight_regg_plus);
           if (!(jet_rec_sel_regg_plus && proton_rec_sel_regg_plus && rp_left_sel_regg && xi_rec_cms_plus_regg - xi_rec_proton_left_regg<0)){
               t_plus_response.Miss(fabs(t_gen_proton_left_regg), event_weight_regg_plus*norm_regg_plus);
              xi_plus_response.Miss(xi_gen_proton_left_regg, event_weight_regg_plus*norm_regg_plus);
              logx_plus_response.Miss(log10(x_gen_left_regg), event_weight_regg_plus*norm_regg_plus);
               t_plus_response_back.Fake(fabs(t_gen_proton_left_regg), event_weight_regg_plus*norm_regg_plus);
              xi_plus_response_back.Fake(xi_gen_proton_left_regg, event_weight_regg_plus*norm_regg_plus);
           }
        }   
        // if(jet1_gen_pt_regg_plus>20 && jet2_gen_pt_regg_plus>20 && fabs(jet1_gen_eta_regg_plus)<4.4 && fabs(jet2_gen_eta_regg_plus)<4.4 && proton_gen_sel_regg_plus) reggeon_histos["xi_gen_left_sasha"]->Fill(xi_gen_proton_left_regg, event_weight_regg_plus);
    }


   
    xi_cms_minus_totem_rec_right->Scale(norm_pom);
    xi_cms_minus_totem_rec_right_bin->Scale(norm_pom);
    xi_cms_minus_totem_rec_left->Scale(norm_pom_plus);
    xi_cms_minus_totem_rec_left_bin->Scale(norm_pom_plus);
    xi_rec_right->Scale(norm_pom);
    xi_rec_left->Scale(norm_pom_plus);
    xi_rec_right_cut->Scale(norm_pom);
    xi_rec_right_cut_sasha->Scale(norm_pom);
    // xi_gen_right_cut_sasha->Scale(norm_pom);
    xi_rec_left_cut_sasha->Scale(norm_pom_plus);
    xi_gen_right->Scale(norm_pom);
    // xi_gen_right_sasha->Scale(norm_pom);
    // xi_gen_left_sasha->Scale(norm_pom_plus);
    xi_rec_left_cut->Scale(norm_pom_plus);
    xi_gen_left->Scale(norm_pom_plus);
    t_rec_right_cut->Scale(norm_pom);
    t_gen_right_cut->Scale(norm_pom);
    t_rec_left_cut->Scale(norm_pom_plus);
    t_gen_left_cut->Scale(norm_pom_plus);
    beta_right_cut->Scale(norm_pom);  
    beta_left_cut->Scale(norm_pom_plus);  
    xi_cms_rec_right->Scale(norm_pom);  
    xi_cms_gen_right->Scale(norm_pom);  
    xi_cms_rec_left->Scale(norm_pom);  
    xi_cms_gen_left->Scale(norm_pom);  
    xi_cms_rec_right_sasha->Scale(norm_pom);  
    xi_cms_gen_right_sasha->Scale(norm_pom);  
    xi_cms_rec_left_sasha->Scale(norm_pom);  
    xi_cms_gen_left_sasha->Scale(norm_pom);  
    log_x_rec_right_cut->Scale(norm_pom);
    log_x_gen_right_cut->Scale(norm_pom);
    log_x_rec_left_cut->Scale(norm_pom_plus);
    log_x_gen_left_cut->Scale(norm_pom_plus);
    pt_jet1_right_cut->Scale(norm_pom);
    pt_jet1_left_cut->Scale(norm_pom_plus);
    pt_jet2_right_cut->Scale(norm_pom);
    pt_jet2_left_cut->Scale(norm_pom_plus);
    eta_jet1_right_cut->Scale(norm_pom);
    eta_jet1_left_cut->Scale(norm_pom_plus);
    eta_jet2_right_cut->Scale(norm_pom);
    eta_jet2_left_cut->Scale(norm_pom_plus);
    t_old->Scale(norm_pom);
    t_old_regg->Scale(norm_regg);
    t_old_left->Scale(norm_pom_plus);
    t_old_left_regg->Scale(norm_regg_plus);
    t_old->Add(t_old_regg);
    t_old->Scale(survival_prob*scale_right);
    t_old_left->Add(t_old_left_regg);
    t_old_left->Scale(survival_prob);
    t_gen_left_pt20->Scale(norm_pom_plus*survival_prob);
    t_gen_left_pt30->Scale(norm_pom_plus*survival_prob);
    t_gen_left_pt40->Scale(norm_pom_plus*survival_prob);
    t_gen_left_pt50->Scale(norm_pom_plus*survival_prob);
   
   
    reggeon_histos["xi_cms_minus_totem_rec_right"]->Scale(norm_regg);
    reggeon_histos["xi_cms_minus_totem_rec_right_bin"]->Scale(norm_regg);
    reggeon_histos["xi_cms_minus_totem_rec_left"]->Scale(norm_regg_plus);
    reggeon_histos["xi_cms_minus_totem_rec_left_bin"]->Scale(norm_regg_plus);
    reggeon_histos["xi_rec_right"]->Scale(norm_regg);
    reggeon_histos["xi_rec_left"]->Scale(norm_regg_plus);
    reggeon_histos["xi_rec_right_cut"]->Scale(norm_regg);
    reggeon_histos["xi_rec_right_cut_sasha"]->Scale(norm_regg);
    reggeon_histos["xi_rec_left_cut_sasha"]->Scale(norm_regg_plus);
    reggeon_histos["xi_gen_right_cut"]->Scale(norm_regg);
    // reggeon_histos["xi_gen_right_cut_sasha"]->Scale(norm_regg);
    reggeon_histos["xi_rec_left_cut"]->Scale(norm_regg_plus);
    reggeon_histos["xi_gen_left_cut"]->Scale(norm_regg_plus);
    reggeon_histos["t_rec_right_cut"]->Scale(norm_regg);
    reggeon_histos["t_gen_right_cut"]->Scale(norm_regg);
    reggeon_histos["t_rec_left_cut"]->Scale(norm_regg_plus);    
    reggeon_histos["t_gen_left_cut"]->Scale(norm_regg_plus);    
    reggeon_histos["beta_right_cut"]->Scale(norm_regg);  
    reggeon_histos["beta_left_cut"]->Scale(norm_regg_plus);  
    reggeon_histos["log_x_rec_right_cut"]->Scale(norm_regg);
    reggeon_histos["log_x_gen_right_cut"]->Scale(norm_regg);
    reggeon_histos["log_x_rec_left_cut"]->Scale(norm_regg_plus);
    reggeon_histos["log_x_gen_left_cut"]->Scale(norm_regg_plus);
    reggeon_histos["pt_jet1_right_cut"]->Scale(norm_regg);
    reggeon_histos["pt_jet1_left_cut"]->Scale(norm_regg_plus);
    reggeon_histos["pt_jet2_right_cut"]->Scale(norm_regg);
    reggeon_histos["pt_jet2_left_cut"]->Scale(norm_regg_plus);
    reggeon_histos["eta_jet1_right_cut"]->Scale(norm_regg);
    reggeon_histos["eta_jet1_left_cut"]->Scale(norm_regg_plus);
    reggeon_histos["eta_jet2_right_cut"]->Scale(norm_regg);
    reggeon_histos["eta_jet2_left_cut"]->Scale(norm_regg_plus);
    // reggeon_histos["xi_gen_right_sasha"]->Scale(norm_regg);
    // reggeon_histos["xi_gen_left_sasha"]->Scale(norm_regg_plus);

    reggeon_histos["xi_cms_rec_right"]->Scale(norm_regg);  
    reggeon_histos["xi_cms_gen_right"]->Scale(norm_regg);  
    reggeon_histos["xi_cms_rec_left"]->Scale(norm_regg_plus);  
    reggeon_histos["xi_cms_gen_left"]->Scale(norm_regg_plus);  
    reggeon_histos["xi_cms_rec_right_sasha"]->Scale(norm_regg);  
    reggeon_histos["xi_cms_gen_right_sasha"]->Scale(norm_regg);  
    reggeon_histos["xi_cms_rec_left_sasha"]->Scale(norm_regg_plus);  
    reggeon_histos["xi_cms_gen_left_sasha"]->Scale(norm_regg_plus);  

    xi_cms_minus_totem_rec_right->Add(reggeon_histos["xi_cms_minus_totem_rec_right"]);
    xi_cms_minus_totem_rec_right_bin->Add(reggeon_histos["xi_cms_minus_totem_rec_right_bin"]);
    xi_cms_minus_totem_rec_left->Add(reggeon_histos["xi_cms_minus_totem_rec_left"]);
    xi_cms_minus_totem_rec_left_bin->Add(reggeon_histos["xi_cms_minus_totem_rec_left_bin"]);
    xi_rec_right->Add(reggeon_histos["xi_rec_right"]);
    xi_rec_left->Add(reggeon_histos["xi_rec_left"]);
    xi_rec_right_cut->Add(reggeon_histos["xi_rec_right_cut"]);
    xi_rec_right_cut_sasha->Add(reggeon_histos["xi_rec_right_cut_sasha"]);
    xi_gen_right_cut_sasha->Add(reggeon_histos["xi_gen_right_cut_sasha"]);
    xi_rec_left_cut_sasha->Add(reggeon_histos["xi_rec_left_cut_sasha"]);
    xi_gen_right->Add(reggeon_histos["xi_gen_right_cut"]);
    // xi_gen_right_sasha->Add(reggeon_histos["xi_gen_right_sasha"]);
    // xi_gen_left_sasha->Add(reggeon_histos["xi_gen_left_sasha"]);
    xi_rec_left_cut->Add(reggeon_histos["xi_rec_left_cut"]);
    xi_gen_left->Add(reggeon_histos["xi_gen_left_cut"]);
    t_rec_right_cut->Add(reggeon_histos["t_rec_right_cut"]);
    t_gen_right_cut->Add(reggeon_histos["t_gen_right_cut"]);
    // t_rec_left_cut->Add(reggeon_histos["t_rec_left_cut"]);
    t_gen_left_cut->Add(reggeon_histos["t_gen_left_cut"]);
    beta_right_cut->Add(reggeon_histos["beta_right_cut"]);  
    beta_left_cut->Add(reggeon_histos["beta_left_cut"]);  
    xi_cms_rec_right->Add(reggeon_histos["xi_cms_rec_right"]);
    xi_cms_gen_right->Add(reggeon_histos["xi_cms_gen_right"]);
    xi_cms_rec_left->Add(reggeon_histos["xi_cms_rec_left"]);
    xi_cms_gen_left->Add(reggeon_histos["xi_cms_gen_left"]);
    xi_cms_rec_right_sasha->Add(reggeon_histos["xi_cms_rec_right_sasha"]);
    xi_cms_gen_right_sasha->Add(reggeon_histos["xi_cms_gen_right_sasha"]);
    xi_cms_rec_left_sasha->Add(reggeon_histos["xi_cms_rec_left_sasha"]);
    xi_cms_gen_left_sasha->Add(reggeon_histos["xi_cms_gen_left_sasha"]);
    log_x_rec_right_cut->Add(reggeon_histos["log_x_rec_right_cut"]);
    log_x_gen_right_cut->Add(reggeon_histos["log_x_gen_right_cut"]);
    log_x_rec_left_cut->Add(reggeon_histos["log_x_rec_left_cut"]);
    log_x_gen_left_cut->Add(reggeon_histos["log_x_gen_left_cut"]);
    pt_jet1_right_cut->Add(reggeon_histos["pt_jet1_right_cut"]);
    pt_jet2_right_cut->Add(reggeon_histos["pt_jet2_right_cut"]);
    pt_jet1_left_cut->Add(reggeon_histos["pt_jet1_left_cut"]);
    pt_jet2_left_cut->Add(reggeon_histos["pt_jet2_left_cut"]);
    eta_jet1_right_cut->Add(reggeon_histos["eta_jet1_right_cut"]);
    eta_jet2_right_cut->Add(reggeon_histos["eta_jet2_right_cut"]);
    eta_jet1_left_cut->Add(reggeon_histos["eta_jet1_left_cut"]);
    eta_jet2_left_cut->Add(reggeon_histos["eta_jet2_left_cut"]);


    xi_cms_minus_totem_rec_right->Scale(survival_prob*scale_right);
    xi_cms_minus_totem_rec_right_bin->Scale(survival_prob*scale_right);
    xi_rec_right->Scale(survival_prob*scale_right);
    xi_rec_right_cut->Scale(survival_prob*scale_right);
    xi_rec_right_cut_sasha->Scale(survival_prob*scale_right);
    xi_gen_right_cut_sasha->Scale(survival_prob);
    xi_gen_right->Scale(survival_prob/*scale_right*/);
    // xi_gen_right_sasha->Scale(survival_prob*scale_right);
    t_rec_right_cut->Scale(survival_prob*scale_right);
    t_gen_right_cut->Scale(survival_prob/*scale_right*/);
    beta_right_cut->Scale(survival_prob*scale_right);  
    log_x_rec_right_cut->Scale(survival_prob*scale_right);
    log_x_gen_right_cut->Scale(survival_prob/*scale_right*/);
    pt_jet1_right_cut->Scale(survival_prob*scale_right);
    pt_jet2_right_cut->Scale(survival_prob*scale_right);
    eta_jet1_right_cut->Scale(survival_prob*scale_right);
    eta_jet2_right_cut->Scale(survival_prob*scale_right);
    xi_cms_rec_right_sasha->Scale(scale_right);
    xi_cms_gen_right_sasha->Scale(survival_prob);

    xi_rec_left_cut_sasha->Scale(survival_prob);
    // xi_gen_left_sasha->Scale(survival_prob);
    xi_cms_minus_totem_rec_left->Scale(survival_prob);
    xi_cms_minus_totem_rec_left_bin->Scale(survival_prob);
    xi_rec_left->Scale(survival_prob);
    xi_rec_left_cut->Scale(survival_prob);
    xi_gen_left->Scale(survival_prob);
    t_rec_left_cut->Scale(survival_prob);
    t_gen_left_cut->Scale(survival_prob);
    beta_left_cut->Scale(survival_prob); 
    log_x_rec_left_cut->Scale(survival_prob);
    log_x_gen_left_cut->Scale(survival_prob);
    pt_jet1_left_cut->Scale(survival_prob);
    pt_jet2_left_cut->Scale(survival_prob);
    eta_jet1_left_cut->Scale(survival_prob);
    eta_jet2_left_cut->Scale(survival_prob);
  

    TH1F* t_right_cut = new TH1F("t_right_cut","", 8, tbins);
    TH1F* t_right_cut_bh = new TH1F("t_right_cut","", 8, tbins);
    TH1F* t_right_cut_noeff = new TH1F("t_right_cut_noeff","", 8, tbins);
    TH1F* t_left_cut = new TH1F("t_right_cut","", 8, tbins);
    TH1F* t_left_cut_bh = new TH1F("t_right_cut","", 8, tbins);
    TH1F* t_left_cut_noeff = new TH1F("t_right_cut_noeff","", 8, tbins);
    TH1F* xi_right_cut = new TH1F("xi_right_cut","",11,xi_bins);
    TH1F* xi_right_cut_bh = new TH1F("xi_right_cut","",11,xi_bins);
    TH1F* xi_right_cut_noeff = new TH1F("xi_right_cut_noeff","",11,xi_bins);
    TH1F* xi_left_cut = new TH1F("xi_left_cut","",11,xi_bins);
    TH1F* xi_left_cut_bh = new TH1F("xi_left_cut","",11,xi_bins);
    TH1F* xi_left_cut_noeff = new TH1F("xi_left_cut_noeff","",11,xi_bins);
    TH1F* log_x_right_cut = new TH1F("log_x_right_cut","",15, -4, 0);
    TH1F* log_x_right_cut_bh = new TH1F("log_x_right_cut","",15, -4, 0);
    TH1F* log_x_right_cut_noeff = new TH1F("log_x_right_cut_noeff","",15, -4, 0);
    TH1F* log_x_right_noeff = new TH1F("log_x_right_noeff","",15, -4, 0);
    TH1F* log_x_left_cut = new TH1F("log_x_left_cut","",15, -4, 0);
    TH1F* log_x_left_cut_bh = new TH1F("log_x_left_cut","",15, -4, 0);
    TH1F* log_x_left_cut_noeff = new TH1F("log_x_left_cut_noeff","",15, -4, 0);
    TH1F* log_x_left_noeff = new TH1F("log_x_left_noeff","",15, -4, 0);
    TH1F* h1;
    TH1F* h2; 
    TH1F* h3; TH1F* h4; TH1F* h5; TH1F* h6; TH1F* h7; TH1F* h8; TH1F* h9; TH1F* h10; TH1F* h11; TH1F* h12; TH1F* h13; TH1F* h14; TH1F* h15; TH1F* h16; TH1F* h17; TH1F* h18; TH1F* h19; TH1F* h20; 
    TH1F* h21; TH1F* h22; TH1F* h23; TH1F* h24; TH1F* h25; TH1F* h26; TH1F* h27; TH1F* h28; TH1F* h29; TH1F* h30; TH1F* h31; TH1F* h32; TH1F* h33; TH1F* h34; TH1F* h35; TH1F* h36; TH1F* h37; TH1F* h38;
     TH1F* h39; TH1F* h40; TH1F* h41; 
    TH1F* h42; TH1F* h43; TH1F* h44; TH1F* h45; TH1F* h46; TH1F* h47; TH1F* h48; TH1F* h49; TH1F* h50; TH1F* h51; TH1F* h52; TH1F* h53; TH1F* h54; TH1F* h55; TH1F* h56;
    TH1F* sasha = new TH1F("sigma_xi_cms_right_sasha","",8, bin_sasha);

    if (unc_rp && !rp_x_min && !rp_x_max){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, h53, sasha, h55, h56, true, true,false,false, false, false);}

    if (unc_rp == false && rp_x_min == true && rp_x_max == false){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, h53, sasha, h55, h56, false, true,true,false, false, false);}

    if (unc_rp == false && rp_x_min == false && rp_x_max == true){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, h53, sasha, h55, h56, false, true,false,true, false, false);}

    if (!unc_rp && !rp_x_min && !rp_x_max){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, h53, sasha, h55, h56,  false, true, false, false, false, false);}

    if (sys_xi_up == true){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, h53, sasha, h55, h56,  false, true, false, false, true, false);}

    if (sys_xi_dw == true){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, h53, sasha, h55, h56,  false, true, false, false, false, true);}
        

    TH1F* t_right_cut_herabackg = (TH1F*) t_right_cut->Clone();
    TH1F* t_left_cut_herabackg = (TH1F*) t_left_cut->Clone();
    TH1F* xi_right_cut_herabackg = (TH1F*) xi_right_cut->Clone();
    TH1F* xi_left_cut_herabackg = (TH1F*) xi_left_cut->Clone();
    TH1F* log_x_right_cut_herabackg = (TH1F*) log_x_right_cut->Clone();
    TH1F* log_x_left_cut_herabackg = (TH1F*) log_x_left_cut->Clone();
    TH1F* t_right_cut_backg;
    TH1F* t_left_cut_backg;
    TH1F* xi_right_cut_backg;
    TH1F* xi_left_cut_backg;
    TH1F* log_x_right_cut_backg;
    TH1F* log_x_left_cut_backg;

    pu_zb (h1, h2, h3, h4, h5, h6, xi_right_cut_backg, xi_left_cut_backg, h7, h8, t_right_cut_backg, t_left_cut_backg, log_x_right_cut_backg, log_x_left_cut_backg,h9, h10, false, true);

    t_right_cut_herabackg->Add(t_right_cut_bh, -1); 
    t_left_cut_herabackg->Add(t_left_cut_bh, -1); 
    xi_right_cut_herabackg->Add(xi_right_cut_bh, -1); 
    xi_left_cut_herabackg->Add(xi_right_cut_bh, -1); 
    log_x_right_cut_herabackg->Add(log_x_right_cut_bh, -1); 
    log_x_left_cut_herabackg->Add(log_x_left_cut_bh, -1); 

    // t_right_cut->Add(t_right_cut_backg, -1); 
    // t_left_cut->Add(t_left_cut_backg, -1); 
    // xi_right_cut->Add(xi_right_cut_backg, -1); 
    // xi_left_cut->Add(xi_left_cut_backg, -1); 
    // log_x_right_cut->Add(log_x_right_cut_backg, -1); 
    // log_x_left_cut->Add(log_x_left_cut_backg, -1); 

    RooUnfoldBayes unfold_right_t (&t_minus_response, t_right_cut, 5);
    t_right_unfolded = (TH1F*)unfold_right_t.Hreco();
    RooUnfoldBayes unfold_right_t_noeff (&t_minus_response, t_right_cut_noeff, 4);
    t_right_unfolded_noeff = (TH1F*)unfold_right_t_noeff.Hreco();
    RooUnfoldBinByBin unfold_right_t_binbybin (&t_minus_response, t_right_cut);
    t_right_unfolded_binbybin = (TH1F*)unfold_right_t_binbybin.Hreco();
    RooUnfoldBayes unfold_right_t_herabackg (&t_minus_response, t_right_cut_herabackg, 4);
    t_right_unfolded_herabackg = (TH1F*)unfold_right_t_herabackg.Hreco();
     
    RooUnfoldBayes unfold_left_t (&t_plus_response, t_left_cut, 4);
    t_left_unfolded = (TH1F*)unfold_left_t.Hreco();
    RooUnfoldBayes unfold_left_t_noeff (&t_plus_response, t_left_cut_noeff, 4);
    t_left_unfolded_noeff = (TH1F*)unfold_left_t_noeff.Hreco();
    RooUnfoldBinByBin unfold_left_t_binbybin (&t_plus_response, t_left_cut);
    t_left_unfolded_binbybin = (TH1F*)unfold_left_t_binbybin.Hreco();
    RooUnfoldBayes unfold_left_t_herabackg (&t_plus_response, t_left_cut_herabackg, 4);
    t_left_unfolded_herabackg = (TH1F*)unfold_left_t_herabackg.Hreco();

    RooUnfoldBayes unfold_right_xi (&xi_minus_response, xi_right_cut, 4);
    xi_right_unfolded = (TH1F*)unfold_right_xi.Hreco();
    RooUnfoldBayes unfold_right_xi_noeff (&xi_minus_response, xi_right_cut_noeff, 4);
    xi_right_unfolded_noeff = (TH1F*)unfold_right_xi_noeff.Hreco();
    RooUnfoldBinByBin unfold_right_xi_binbybin (&xi_minus_response, xi_right_cut);
    xi_right_unfolded_binbybin = (TH1F*)unfold_right_xi_binbybin.Hreco();
    RooUnfoldBayes unfold_right_xi_herabackg (&xi_minus_response, xi_right_cut_herabackg, 4);
    xi_right_unfolded_herabackg = (TH1F*)unfold_right_xi_herabackg.Hreco();

    RooUnfoldBayes unfold_left_xi (&xi_plus_response, xi_left_cut, 4);
    xi_left_unfolded = (TH1F*)unfold_left_xi.Hreco();
    RooUnfoldBayes unfold_left_xi_noeff (&xi_plus_response, xi_left_cut_noeff, 4);
    xi_left_unfolded_noeff = (TH1F*)unfold_left_xi_noeff.Hreco();
    RooUnfoldBinByBin unfold_left_xi_binbybin (&xi_plus_response, xi_left_cut);
    xi_left_unfolded_binbybin = (TH1F*)unfold_left_xi_binbybin.Hreco();
    RooUnfoldBayes unfold_left_xi_herabackg (&xi_plus_response, xi_left_cut_herabackg, 4);
    xi_left_unfolded_herabackg = (TH1F*)unfold_left_xi_herabackg.Hreco();

    RooUnfoldBayes unfold_right_logx (&logx_minus_response, log_x_right_cut, 4);
    logx_right_unfolded = (TH1F*)unfold_right_logx.Hreco();
    RooUnfoldBayes unfold_right_logx_cut_noeff (&logx_minus_response, log_x_right_cut_noeff, 4);
    logx_right_cut_unfolded_noeff = (TH1F*)unfold_right_logx_cut_noeff.Hreco();
    RooUnfoldBayes unfold_right_logx_noeff (&logx_minus_response, log_x_right_noeff, 4);
    logx_right_unfolded_noeff = (TH1F*)unfold_right_logx_noeff.Hreco();
    RooUnfoldBinByBin unfold_right_logx_binbybin (&logx_minus_response, log_x_right_cut);
    logx_right_unfolded_binbybin = (TH1F*)unfold_right_logx_binbybin.Hreco();
    RooUnfoldBayes unfold_right_logx_herabackg (&logx_minus_response, log_x_right_cut_herabackg, 4);
    logx_right_unfolded_herabackg = (TH1F*)unfold_right_logx_herabackg.Hreco();

    RooUnfoldBayes unfold_left_logx (&logx_plus_response, log_x_left_cut, 4);
    logx_left_unfolded = (TH1F*)unfold_left_logx.Hreco();
    RooUnfoldBayes unfold_left_logx_cut_noeff (&logx_plus_response, log_x_left_cut_noeff, 4);
    logx_left_cut_unfolded_noeff = (TH1F*)unfold_left_logx_cut_noeff.Hreco();
    RooUnfoldBayes unfold_left_logx_noeff (&logx_plus_response, log_x_left_noeff, 4);
    logx_left_unfolded_noeff = (TH1F*)unfold_left_logx_noeff.Hreco();
    RooUnfoldBinByBin unfold_left_logx_binbybin (&logx_plus_response, log_x_left_cut);
    logx_left_unfolded_binbybin = (TH1F*)unfold_left_logx_binbybin.Hreco();
    RooUnfoldBayes unfold_left_logx_herabackg (&logx_plus_response, log_x_left_cut_herabackg, 4);
    logx_left_unfolded_herabackg = (TH1F*)unfold_left_logx_herabackg.Hreco();

TH1F* t_right_cut_up; TH1F* t_left_cut_up; TH1F* xi_right_cut_up; TH1F* xi_left_cut_up; TH1F* log_x_minus_up; TH1F* log_x_plus_up; TH1F* log_x_right_cut_up; TH1F* log_x_left_cut_up;
TH1F* t_right_cut_dw; TH1F* t_left_cut_dw; TH1F* xi_right_cut_dw; TH1F* xi_left_cut_dw; TH1F* log_x_minus_dw; TH1F* log_x_plus_dw; TH1F* log_x_right_cut_dw; TH1F* log_x_left_cut_dw;
TH1F* t_right_cut_up_pf; TH1F* t_left_cut_up_pf; TH1F* xi_right_cut_up_pf; TH1F* xi_left_cut_up_pf; TH1F* log_x_minus_up_pf; TH1F* log_x_plus_up_pf; TH1F* log_x_right_cut_up_pf; TH1F* log_x_left_cut_up_pf;
TH1F* t_right_cut_dw_pf; TH1F* t_left_cut_dw_pf; TH1F* xi_right_cut_dw_pf; TH1F* xi_left_cut_dw_pf; TH1F* log_x_minus_dw_pf; TH1F* log_x_plus_dw_pf; TH1F* log_x_right_cut_dw_pf; TH1F* log_x_left_cut_dw_pf;
unc(t_right_cut_up, t_left_cut_up, xi_right_cut_up, xi_left_cut_up, log_x_minus_up, log_x_plus_up, log_x_right_cut_up, log_x_left_cut_up, true, true, true);
unc(t_right_cut_dw, t_left_cut_dw, xi_right_cut_dw, xi_left_cut_dw, log_x_minus_dw, log_x_plus_dw, log_x_right_cut_dw, log_x_left_cut_dw, true, false, true);
unc(t_right_cut_up_pf, t_left_cut_up_pf, xi_right_cut_up_pf, xi_left_cut_up_pf, log_x_minus_up_pf, log_x_plus_up_pf, log_x_right_cut_up_pf, log_x_left_cut_up_pf, false, true, true);
unc(t_right_cut_dw_pf, t_left_cut_dw_pf, xi_right_cut_dw_pf, xi_left_cut_dw_pf, log_x_minus_dw_pf, log_x_plus_dw_pf, log_x_right_cut_dw_pf, log_x_left_cut_dw_pf, false, false, true);

    RooUnfoldBayes unfold_right_t_up (&t_minus_response, t_right_cut_up, 4);
    t_right_unfolded_up = (TH1F*)unfold_right_t_up.Hreco();
    RooUnfoldBayes unfold_right_t_dw (&t_minus_response, t_right_cut_dw, 4);
    t_right_unfolded_dw = (TH1F*)unfold_right_t_dw.Hreco();

    RooUnfoldBayes unfold_left_t_up (&t_plus_response, t_left_cut_up, 4);
    t_left_unfolded_up = (TH1F*)unfold_left_t_up.Hreco();
    RooUnfoldBayes unfold_left_t_dw (&t_plus_response, t_left_cut_dw, 4);
    t_left_unfolded_dw = (TH1F*)unfold_left_t_dw.Hreco();

    RooUnfoldBayes unfold_right_xi_up (&xi_minus_response, xi_right_cut_up, 4);
    xi_right_unfolded_up = (TH1F*)unfold_right_xi_up.Hreco();
    RooUnfoldBayes unfold_right_xi_dw (&xi_minus_response, xi_right_cut_dw, 4);
    xi_right_unfolded_dw = (TH1F*)unfold_right_xi_dw.Hreco();

    RooUnfoldBayes unfold_left_xi_up (&xi_plus_response, xi_left_cut_up, 4);
    xi_left_unfolded_up = (TH1F*)unfold_left_xi_up.Hreco();
    RooUnfoldBayes unfold_left_xi_dw (&xi_plus_response, xi_left_cut_dw, 4);
    xi_left_unfolded_dw = (TH1F*)unfold_left_xi_dw.Hreco();

    RooUnfoldBayes unfold_right_log_x_up (&logx_minus_response, log_x_right_cut_up, 4);
    log_x_right_unfolded_up = (TH1F*)unfold_right_log_x_up.Hreco();
    RooUnfoldBayes unfold_right_log_x_dw (&logx_minus_response, log_x_right_cut_dw, 4);
    log_x_right_unfolded_dw = (TH1F*)unfold_right_log_x_dw.Hreco();

    RooUnfoldBayes unfold_left_log_x_up (&logx_plus_response, log_x_left_cut_up, 4);
    log_x_left_unfolded_up = (TH1F*)unfold_left_log_x_up.Hreco();
    RooUnfoldBayes unfold_left_log_x_dw (&logx_plus_response, log_x_left_cut_dw, 4);
    log_x_left_unfolded_dw = (TH1F*)unfold_left_log_x_dw.Hreco();

    //pf
    RooUnfoldBayes unfold_right_t_up_pf (&t_minus_response, t_right_cut_up_pf, 4);
    t_right_unfolded_up_pf = (TH1F*)unfold_right_t_up_pf.Hreco();
    RooUnfoldBayes unfold_right_t_dw_pf (&t_minus_response, t_right_cut_dw_pf, 4);
    t_right_unfolded_dw_pf = (TH1F*)unfold_right_t_dw_pf.Hreco();

    RooUnfoldBayes unfold_left_t_up_pf (&t_plus_response, t_left_cut_up_pf, 4);
    t_left_unfolded_up_pf = (TH1F*)unfold_left_t_up_pf.Hreco();
    RooUnfoldBayes unfold_left_t_dw_pf (&t_plus_response, t_left_cut_dw_pf, 4);
    t_left_unfolded_dw_pf = (TH1F*)unfold_left_t_dw_pf.Hreco();

    RooUnfoldBayes unfold_right_xi_up_pf (&xi_minus_response, xi_right_cut_up_pf, 4);
    xi_right_unfolded_up_pf = (TH1F*)unfold_right_xi_up_pf.Hreco();
    RooUnfoldBayes unfold_right_xi_dw_pf (&xi_minus_response, xi_right_cut_dw_pf, 4);
    xi_right_unfolded_dw_pf = (TH1F*)unfold_right_xi_dw_pf.Hreco();

    RooUnfoldBayes unfold_left_xi_up_pf (&xi_plus_response, xi_left_cut_up_pf, 4);
    xi_left_unfolded_up_pf = (TH1F*)unfold_left_xi_up_pf.Hreco();
    RooUnfoldBayes unfold_left_xi_dw_pf (&xi_plus_response, xi_left_cut_dw_pf, 4);
    xi_left_unfolded_dw_pf = (TH1F*)unfold_left_xi_dw_pf.Hreco();

    RooUnfoldBayes unfold_right_log_x_up_pf (&logx_minus_response, log_x_right_cut_up_pf, 4);
    log_x_right_unfolded_up_pf = (TH1F*)unfold_right_log_x_up_pf.Hreco();
    RooUnfoldBayes unfold_right_log_x_dw_pf (&logx_minus_response, log_x_right_cut_dw_pf, 4);
    log_x_right_unfolded_dw_pf = (TH1F*)unfold_right_log_x_dw_pf.Hreco();

    RooUnfoldBayes unfold_left_log_x_up_pf (&logx_plus_response, log_x_left_cut_up_pf, 4);
    log_x_left_unfolded_up_pf = (TH1F*)unfold_left_log_x_up_pf.Hreco();
    RooUnfoldBayes unfold_left_log_x_dw_pf (&logx_plus_response, log_x_left_cut_dw_pf, 4);
    log_x_left_unfolded_dw_pf = (TH1F*)unfold_left_log_x_dw_pf.Hreco();

    //manual binbybin
    // t_rec_right_cut->Divide(t_gen_right_cut);
    // t_rec_left_cut->Divide(t_gen_left_cut);
    // xi_rec_right_cut->Divide(xi_gen_right);
    // xi_rec_left_cut->Divide(xi_gen_left);
    // log_x_rec_right_cut->Divide(log_x_gen_right_cut);
    // log_x_rec_left_cut->Divide(log_x_gen_left_cut);

    // t_right_unfolded_up=(TH1F*)t_right_cut_up->Clone();
    // t_right_unfolded_dw=(TH1F*)t_right_cut_dw->Clone();
    // xi_right_unfolded_up=(TH1F*)xi_right_cut_up->Clone();
    // xi_right_unfolded_dw=(TH1F*)xi_right_cut_dw->Clone();
    // log_x_right_unfolded_up=(TH1F*)log_x_right_cut_up->Clone();
    // log_x_right_unfolded_dw=(TH1F*)log_x_right_cut_dw->Clone();
    // t_left_unfolded_up=(TH1F*)t_left_cut_up->Clone();
    // t_left_unfolded_dw=(TH1F*)t_left_cut_dw->Clone();
    // xi_left_unfolded_up=(TH1F*)xi_left_cut_up->Clone();
    // xi_left_unfolded_dw=(TH1F*)xi_left_cut_dw->Clone();
    // log_x_left_unfolded_up=(TH1F*)log_x_left_cut_up->Clone();
    // log_x_left_unfolded_dw=(TH1F*)log_x_left_cut_dw->Clone();

    //  t_right_unfolded_up->Divide(t_rec_right_cut);
    //  t_right_unfolded_dw->Divide(t_rec_right_cut);
    //  xi_right_unfolded_up->Divide(xi_rec_right_cut);
    //  xi_right_unfolded_dw->Divide(xi_rec_right_cut);
    //  log_x_right_unfolded_up->Divide(log_x_rec_right_cut);
    //  log_x_right_unfolded_dw->Divide(log_x_rec_right_cut);
    //  t_left_unfolded_up->Divide(t_rec_left_cut);
    //  t_left_unfolded_dw->Divide(t_rec_left_cut);
    //  xi_left_unfolded_up->Divide(xi_rec_left_cut);
    //  xi_left_unfolded_dw->Divide(xi_rec_left_cut);
    //  log_x_left_unfolded_up->Divide(log_x_rec_left_cut);
    //  log_x_left_unfolded_dw->Divide(log_x_rec_left_cut);


//     TCanvas *c2 = new TCanvas("","right");
// logx_right_unfolded->Draw();
// log_x_gen_right_cut->Draw("histsame");

    // sasha->Scale(1,"width");
    // sasha->SetBinContent(1,sasha->GetBinContent(1)-9133);
    xi_cms_rec_right_sasha->Scale(0.17,"width");
    // xi_cms_gen_right_sasha->Scale(1,"width");
    // xi_cms_rec_right_sasha->Divide(xi_cms_gen_right_sasha);
    // sasha->Divide(xi_cms_rec_right_sasha);
    // xi_cms_gen_right_sasha->Scale(0.001/luminosity);
    // xi_cms_gen_right_sasha->Draw("E1");
    // // xi_cms_rec_right_sasha->Draw("histsame");



// t_rec_right_cut->Divide(t_gen_right_cut); 
// t_rec_right_cut
// t_rec_left_cut->Divide(t_gen_left_cut); 
// t_rec_left_cut->SetLineColor(2);   
// t_rec_left_cut->Draw("E1same");   
cout<<"POMWIG ENTRIES GEN:    right: "<<t_gen_right_cut->Integral()<<"   left: "<<t_gen_right_cut->Integral()<<endl;
// log_x_right_unfolded_up->Draw("hist");
// log_x_right_unfolded_dw->Draw("histsame");
// logx_right_unfolded->Draw("E1same");

    // RooUnfoldBayes unfold_right_t_test (&t_minus_response, t_rec_right_cut, 4);
    // TH1F* t_right_unfolded_rec = (TH1F*)unfold_right_t_test.Hreco(); 
    // RooUnfoldBayes unfold_right_t_back (&t_minus_response_back, t_right_unfolded_rec, 4);
    // TH1F* t_right_unfolded_back = (TH1F*)unfold_right_t_back.Hreco(); 
    // RooUnfoldBayes unfold_right_t_back_data (&t_minus_response_back, t_right_unfolded, 4);
    // TH1F* t_right_unfolded_back_data = (TH1F*)unfold_right_t_back_data.Hreco(); 

    // RooUnfoldBayes unfold_left_t_test (&t_plus_response, t_rec_left_cut, 4);
    // TH1F* t_left_unfolded_rec = (TH1F*)unfold_left_t_test.Hreco(); 
    // RooUnfoldBayes unfold_left_t_back (&t_plus_response_back, t_left_unfolded_rec, 4);
    // TH1F* t_left_unfolded_back = (TH1F*)unfold_left_t_back.Hreco(); 
    // RooUnfoldBayes unfold_left_t_back_data (&t_plus_response_back, t_left_unfolded, 4);
    // TH1F* t_left_unfolded_back_data = (TH1F*)unfold_left_t_back_data.Hreco(); 
    // t_left_unfolded_back->Draw("hist");
    // t_rec_left_cut->Draw("E1same");


    // RooUnfoldBayes unfold_left_xi_test (&xi_plus_response, xi_rec_left_cut, 4);
    // TH1F* xi_left_unfolded_rec = (TH1F*)unfold_left_xi_test.Hreco(); 
    // RooUnfoldBayes unfold_left_xi_back (&xi_plus_response_back, xi_left_unfolded_rec, 4);
    // TH1F* xi_left_unfolded_back = (TH1F*)unfold_left_xi_back.Hreco(); 
    // xi_left_unfolded_back->Draw();
    // xi_rec_left_cut->Draw("histsame");

    // RooUnfoldBayes unfold_right_xi_test (&xi_minus_response, xi_rec_right_cut, 4);
    // TH1F* xi_right_unfolded_rec = (TH1F*)unfold_right_xi_test.Hreco(); 
    // RooUnfoldBayes unfold_right_xi_back (&xi_minus_response_back, xi_right_unfolded_rec, 4);
    // TH1F* xi_right_unfolded_back = (TH1F*)unfold_right_xi_back.Hreco(); 
    // xi_right_unfolded_back->Draw();
    // xi_rec_right_cut->Draw("histsame");


    // RooUnfoldBayes unfold_right_logx_test (&logx_minus_response, log_x_rec_right_cut, 4);
    // TH1F* logx_right_unfolded_rec = (TH1F*)unfold_right_logx_test.Hreco();

    // RooUnfoldBayes unfold_right_logx_back (&logx_minus_response_back, logx_right_unfolded, 4);
    // TH1F* logx_right_unfolded_back = (TH1F*)unfold_right_logx_back.Hreco();
    // logx_right_unfolded_rec_back->Draw();
    // log_x_rec_right_cut->Draw("histsame");

    // TFile *check = new TFile("~/Dropbox/doctorado/note/scripts/root_files/unfolding_checks_backunfolded.root","UPDATE");
    // // check->WriteTObject(t_right_unfolded,"t_right_unfolded"); 
    // check->WriteTObject(t_right_unfolded_back,"t_right_backunfolded"); 
    // check->WriteTObject(t_right_cut,"t_right_cut"); 
    // check->WriteTObject(xi_right_unfolded_back,"xi_right_backunfolded"); 
    // check->WriteTObject(xi_right_cut,"xi_right_cut"); 
    // check->WriteTObject(logx_right_unfolded_back,"logx_right_backunfolded"); 
    // check->WriteTObject(log_x_right_cut,"logx_right_cut"); 
        
  //   TMatrix matrix_t_minus(8+2,8+2,t_minus_th2->GetArray(),"F");
  //   TDecompSVD svd(matrix_t_minus);
  // //   svd.Print();
  // // cout << "decomposed   " << svd.kDecomposed << endl;
  // // cout << "condition nb " << svd.GetCondition() << endl;
  // // cout << "tolerance    " << svd.GetTol() << endl;
  // // cout << "singular     " << svd.kSingular << endl;
  // const TVectorD singularVal = svd.GetSig();
  // singularVal.Print();


// unfold_right_t.PrintTable(cout,t_gen_right_cut);
// TMatrixD cov_right_t = (TMatrixD)unfold_right_t.Ereco();
// cov_right_t.Draw();
// cout<<unfold_right_t.Chi2(t_gen_right_cut)/*<<" "<<t_right_cut.Chi2(t_gen_right_cut)*/<<endl;

// logx_right_unfolded->Draw();
// logx_left_unfolded->Draw("same");

    double jet1_rec_pt_pythia, jet1_rec_eta_pythia, jet2_rec_pt_pythia, jet2_rec_eta_pythia;
    double jet1_gen_pt_pythia, jet1_gen_eta_pythia, jet2_gen_pt_pythia, jet2_gen_eta_pythia;
    double xi_rec_cms_minus_pythia, xi_rec_proton_right_pythia;
    double xi_gen_cms_minus_pythia, xi_gen_proton_right_pythia;
    double t_rec_proton_right_pythia;
    double t_gen_proton_right_pythia;
    double weight_pythia, x_rec_right_pythia, x_gen_right_pythia;
    bool rp_right_pythia, rp_left_pythia;
    double theta_x_minus_smear_pythia, theta_y_minus_smear_pythia, theta_x_plus_smear_pythia, theta_y_plus_smear_pythia, pz_minus_smear_pythia, pz_plus_smear_pythia, e_minus_smear_pythia, e_plus_smear_pythia, pz_minus_pythia, px_minus_pythia, py_minus_pythia, e_minus_pythia, mass_plus_pythia, mass_minus_pythia;
    double px_plus_pythia, py_plus_pythia, pz_plus_pythia;
        TTree* tree_pythia= (TTree*) pythia8->Get( treeName.c_str() );
    int nev = int(tree_pythia->GetEntriesFast());
    tree_pythia->SetBranchAddress("xi_rec_proton_right",&xi_rec_proton_right_pythia);
    tree_pythia->SetBranchAddress("xi_gen_proton_right",&xi_gen_proton_right_pythia);
    tree_pythia->SetBranchAddress("t_rec_proton_right",&t_rec_proton_right_pythia);
    tree_pythia->SetBranchAddress("t_gen_proton_right",&t_gen_proton_right_pythia);
    tree_pythia->SetBranchAddress("rp_right",&rp_right_pythia);
    tree_pythia->SetBranchAddress("rp_left",&rp_left_pythia);
    tree_pythia->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt_pythia);
    tree_pythia->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt_pythia);
    tree_pythia->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta_pythia);
    tree_pythia->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta_pythia);
    tree_pythia->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt_pythia);
    tree_pythia->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt_pythia);
    tree_pythia->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta_pythia);
    tree_pythia->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta_pythia);
    tree_pythia->SetBranchAddress("weight_mc",&weight_pythia);
    tree_pythia->SetBranchAddress("xi_rec_cms_right",&xi_rec_cms_minus_pythia);
    tree_pythia->SetBranchAddress("xi_gen_cms_right",&xi_gen_cms_minus_pythia);
    tree_pythia->SetBranchAddress("x_rec_right",&x_rec_right_pythia);
    tree_pythia->SetBranchAddress("x_gen_right",&x_gen_right_pythia);
    tree_pythia->SetBranchAddress("px_proton_right",&px_minus_pythia);
    tree_pythia->SetBranchAddress("py_proton_right",&py_minus_pythia);
    tree_pythia->SetBranchAddress("pz_proton_right",&pz_minus_pythia);
    tree_pythia->SetBranchAddress("e_proton_right",&e_minus_pythia);
    tree_pythia->SetBranchAddress("mass_proton_left",&mass_plus_pythia);
    tree_pythia->SetBranchAddress("mass_proton_right",&mass_minus_pythia);
    tree_pythia->SetBranchAddress("px_proton_left",&px_plus_pythia);
    tree_pythia->SetBranchAddress("py_proton_left",&py_plus_pythia);
    tree_pythia->SetBranchAddress("pz_proton_left",&pz_plus_pythia);
    double norm_pythia = 2109.82;
    double eff_pythia = 0.00133*survival_prob*1.22393;
    TH1F* t_rec_right_cut_pythia = new TH1F("t_rec_right_cut_pom_pythia","", 8, tbins);
    TH1F* t_gen_right_cut_pythia = new TH1F("t_gen_right_cut_pom_pythia","", 8, tbins);
    TH1F* xi_rec_right_cut_pythia = new TH1F("xi_rec_right_cut_pom_pythia","", 11, xi_bins);
    TH1F* xi_gen_right_cut_pythia = new TH1F("xi_gen_right_cut_pom_pythia","", 11, xi_bins);
    TH1F* logx_rec_right_cut_pythia = new TH1F("logx_rec_right_cut_pom_pythia","", 15, -4, 0);
    TH1F* logx_gen_right_cut_pythia = new TH1F("logx_gen_right_cut_pom_pythia","", 15, -4, 0);
    RooUnfoldResponse t_minus_response_pythia (t_rec_right_cut_pythia, t_gen_right_cut_pythia, "t_minus_unfolded_pythia", "t_minus_unfolded_pythia");
    RooUnfoldResponse xi_minus_response_pythia (xi_rec_right_cut_pythia, xi_gen_right_cut_pythia, "xi_minus_unfolded_pythia", "xi_minus_unfolded_pythia");
    RooUnfoldResponse logx_minus_response_pythia (logx_rec_right_cut_pythia, logx_gen_right_cut_pythia, "logx_minus_unfolded_pythia", "logx_minus_unfolded_pythia");
    t_rec_right_cut_pythia->Sumw2(); 
    t_gen_right_cut_pythia->Sumw2(); 
    xi_rec_right_cut_pythia->Sumw2(); 
    xi_gen_right_cut_pythia->Sumw2(); 
    logx_rec_right_cut_pythia->Sumw2(); 
    logx_gen_right_cut_pythia->Sumw2(); 
    TH1F* t_gen_left_pythia_jet20 = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH1F* t_gen_left_pythia_jet20_nomass = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH1F* t_gen_left_pythia_jet30 = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH1F* t_gen_left_pythia_jet40 = new TH1F("t_gen_left_cut_pom","", 100, 0,1);
    TH1F* t_gen_left_pythia_jet50 = new TH1F("t_gen_left_cut_pom","", 100, 0,1);

        for(int i_evt = 0; i_evt < nev; ++i_evt){
        tree_pythia->GetEntry(i_evt);
        double event_weight_minus_pythia = weight_pythia;

        bool jet_rec_sel_pythia = jet1_rec_pt_pythia>pt_threshold && jet2_rec_pt_pythia>pt_threshold && fabs(jet1_rec_eta_pythia)<4.4 && fabs(jet2_rec_eta_pythia)<4.4;
        bool jet_gen_sel_pythia = jet1_gen_pt_pythia>pt_threshold && jet2_gen_pt_pythia>pt_threshold && fabs(jet1_gen_eta_pythia)<4.4 && fabs(jet2_gen_eta_pythia)<4.4;
        bool proton_right_rec_sel_pythia =  xi_rec_proton_right_pythia<0.1 && fabs(t_rec_proton_right_pythia)>0.03 && fabs(t_rec_proton_right_pythia)<1;
        bool proton_right_gen_sel_pythia = xi_gen_proton_right_pythia<0.1 && fabs(t_gen_proton_right_pythia)>0.03 && fabs(t_gen_proton_right_pythia)<1;
p_right_pythia->Fill(sqrt(pz_minus_pythia*pz_minus_pythia ), event_weight_minus_pythia);
           
        //right  
        if (jet_rec_sel_pythia && proton_right_rec_sel_pythia && rp_right_pythia){
            if (xi_rec_cms_minus_pythia - xi_rec_proton_right_pythia<0){
                t_rec_right_cut_pythia->Fill(fabs(t_rec_proton_right_pythia), event_weight_minus_pythia);
                xi_rec_right_cut_pythia->Fill(xi_rec_proton_right_pythia, event_weight_minus_pythia);
                logx_rec_right_cut_pythia->Fill(log10(x_rec_right_pythia), event_weight_minus_pythia);
                if(!(jet_gen_sel_pythia && proton_right_gen_sel_pythia)){
                    t_minus_response_pythia.Fake(fabs(t_rec_proton_right_pythia),event_weight_minus_pythia*eff_pythia*norm_pythia);
                    xi_minus_response_pythia.Fake(xi_rec_proton_right_pythia,event_weight_minus_pythia*eff_pythia*norm_pythia);
                    logx_minus_response_pythia.Fake(log10(x_rec_right_pythia),event_weight_minus_pythia*eff_pythia*norm_pythia);

                }
                if (jet_gen_sel_pythia && proton_right_gen_sel_pythia){
                    t_minus_response_pythia.Fill(fabs(t_rec_proton_right_pythia),fabs(t_gen_proton_right_pythia),event_weight_minus_pythia*norm_pythia*eff_pythia);
                    xi_minus_response_pythia.Fill(xi_rec_proton_right_pythia,xi_gen_proton_right_pythia,event_weight_minus_pythia*eff_pythia*norm_pythia);
                    logx_minus_response_pythia.Fill(log10(x_rec_right_pythia),log10(x_gen_right_pythia),event_weight_minus_pythia*eff_pythia*norm_pythia);
                }   
            }

        }
        if (jet_gen_sel_pythia && proton_right_gen_sel_pythia){
           t_gen_right_cut_pythia->Fill(fabs(t_gen_proton_right_pythia), event_weight_minus_pythia); 
           xi_gen_right_cut_pythia->Fill(xi_gen_proton_right_pythia, event_weight_minus_pythia); 
           logx_gen_right_cut_pythia->Fill(log10(x_gen_right_pythia), event_weight_minus_pythia);
           if(!(jet_rec_sel_pythia && proton_right_rec_sel_pythia && rp_right_pythia && xi_rec_cms_minus_pythia - xi_rec_proton_right_pythia<0)){
                    t_minus_response_pythia.Miss(fabs(t_gen_proton_right_pythia),event_weight_minus_pythia*eff_pythia*norm_pythia);
                    xi_minus_response_pythia.Miss(xi_gen_proton_right_pythia,event_weight_minus_pythia*eff_pythia*norm_pythia);
                    logx_minus_response_pythia.Miss(log10(x_gen_right_pythia),event_weight_minus_pythia*eff_pythia*norm_pythia);
           }
        } 

        double pi_mass_pythia = sqrt(pi*pi - M_P*M_P);
        TLorentzVector p_beam_plus_pythia (0, 0, pi, pi);
        double e_gen_plus_from_p_pythia_nomass = sqrt(px_plus_pythia*px_plus_pythia + py_plus_pythia*py_plus_pythia + pz_plus_pythia*pz_plus_pythia);// + M_P*M_P);
        double e_gen_plus_from_p_pythia = sqrt(px_plus_pythia*px_plus_pythia + py_plus_pythia*py_plus_pythia + pz_plus_pythia*pz_plus_pythia + M_P*M_P);
        TLorentzVector p_scatt_plus_pythia_nomass (px_plus_pythia, py_plus_pythia, pz_plus_pythia, e_gen_plus_from_p_pythia_nomass);
        TLorentzVector p_scatt_plus_pythia (px_plus_pythia, py_plus_pythia, pz_plus_pythia, e_gen_plus_from_p_pythia);
        TLorentzVector t_vec_plus_pythia_nomass = (p_beam_plus_pythia - p_scatt_plus_pythia_nomass);
        TLorentzVector t_vec_plus_pythia = (p_beam_plus_pythia - p_scatt_plus_pythia);
        double t_rec_proton_left_theta_pythia_nomass = t_vec_plus_pythia_nomass.Mag2(); 
        double t_rec_proton_left_theta_pythia = t_vec_plus_pythia.Mag2(); 

        /*if (jet1_gen_pt_pythia>20 && fabs(jet1_rec_eta_pythia)<5)*/ t_gen_left_pythia_jet20_nomass->Fill(fabs(t_rec_proton_left_theta_pythia_nomass), event_weight_minus_pythia); 
        /*if (jet1_gen_pt_pythia>20 && fabs(jet1_rec_eta_pythia)<5)*/ t_gen_left_pythia_jet20->Fill(fabs(t_rec_proton_left_theta_pythia), event_weight_minus_pythia); 
        if (jet1_gen_pt_pythia>30 && fabs(jet1_rec_eta_pythia)<5) t_gen_left_pythia_jet30->Fill(fabs(t_rec_proton_left_theta_pythia), event_weight_minus_pythia); 
        if (jet1_gen_pt_pythia>40 && fabs(jet1_rec_eta_pythia)<5) t_gen_left_pythia_jet40->Fill(fabs(t_rec_proton_left_theta_pythia), event_weight_minus_pythia); 
        if (jet1_gen_pt_pythia>50 && fabs(jet1_rec_eta_pythia)<5) t_gen_left_pythia_jet50->Fill(fabs(t_rec_proton_left_theta_pythia), event_weight_minus_pythia); 



    }
    t_rec_right_cut_pythia->Scale(norm_pythia*eff_pythia);
    t_gen_right_cut_pythia->Scale(norm_pythia*eff_pythia);
    xi_rec_right_cut_pythia->Scale(norm_pythia*eff_pythia);
    xi_gen_right_cut_pythia->Scale(norm_pythia*eff_pythia);
    logx_rec_right_cut_pythia->Scale(norm_pythia*eff_pythia);
    logx_gen_right_cut_pythia->Scale(norm_pythia*eff_pythia);

t_gen_left_pythia_jet20->Scale(norm_pythia*eff_pythia);
t_gen_left_pythia_jet30->Scale(norm_pythia*eff_pythia);
t_gen_left_pythia_jet40->Scale(norm_pythia*eff_pythia);
t_gen_left_pythia_jet50->Scale(norm_pythia*eff_pythia);

// /////////////////////////////MC1 and MC2 ////////////////////////////////////////////////////////
// TH2F* chi2_unfolded_vs_iter_t = new TH2F("","",100,0,100,100,-0.5,3);
// TH2F* chi2_unfolded_vs_iter_t_pythia = new TH2F("","",100,0,100,100,-0.5,3);

// TH2F* chi2_unfolded_vs_iter_xi = new TH2F("","",100,0,100,100,-0.5,3);
// TH2F* chi2_unfolded_vs_iter_xi_pythia = new TH2F("","",100,0,100,100,-0.5,3);
// TH2F* chi2_unfolded_vs_iter_logx = new TH2F("","",100,0,100,100,-0.5,3);
// TH2F* chi2_unfolded_vs_iter_logx_pythia = new TH2F("","",100,0,100,100,-0.5,3);
//    // TFile * outfile = new TFile("sigma_xi_ratio_unfold_error.root", "RECREATE");

// for (int i = 1; i<50; ++i){
//     RooUnfoldBayes logx_unfold_per_iter (&logx_minus_response_pythia, log_x_rec_right_cut, i);
//     TH1F* logx_right_unfolded_per_iter = (TH1F*)logx_unfold_per_iter.Hreco();
//     RooUnfoldBayes logx_unfold_per_iter_pythia (&logx_minus_response, logx_rec_right_cut_pythia, i);
//     TH1F* logx_right_unfolded_per_iter_pythia = (TH1F*)logx_unfold_per_iter_pythia.Hreco();
//     double chi2_unfolded = 0; 
//     double chi2_unfolded_pythia = 0; 
//     for (int j = 1; j<=logx_right_unfolded->GetNbinsX(); ++j){
//          if(logx_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded += pow(logx_right_unfolded_per_iter->GetBinContent(j) - log_x_gen_right_cut->GetBinContent(j), 2)/(pow(logx_right_unfolded_per_iter->GetBinError(j),2));
//          if(logx_right_unfolded_per_iter_pythia->GetBinError(j)>0) chi2_unfolded_pythia += pow(logx_right_unfolded_per_iter_pythia->GetBinContent(j) - logx_gen_right_cut_pythia->GetBinContent(j), 2)/(pow(logx_right_unfolded_per_iter_pythia->GetBinError(j),2));
//     }
//     chi2_unfolded_vs_iter_logx->Fill(i, log10(chi2_unfolded)); 
//     chi2_unfolded_vs_iter_logx_pythia->Fill(i, log10(chi2_unfolded_pythia)); 

//     RooUnfoldBayes xi_unfold_per_iter (&xi_minus_response_pythia, xi_rec_right_cut, i);
//     TH1F* xi_right_unfolded_per_iter = (TH1F*)xi_unfold_per_iter.Hreco();
//     RooUnfoldBayes xi_unfold_per_iter_pythia (&xi_minus_response, xi_rec_right_cut_pythia, i);
//     TH1F* xi_right_unfolded_per_iter_pythia = (TH1F*)xi_unfold_per_iter_pythia.Hreco();
//     double chi2_unfolded_xi = 0; 
//     double chi2_unfolded_xi_pythia = 0; 
//     for (int j = 1; j<=xi_right_unfolded->GetNbinsX(); ++j){
//          if(xi_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_xi += pow(xi_right_unfolded_per_iter->GetBinContent(j) - xi_gen_right->GetBinContent(j), 2)/(pow(xi_right_unfolded_per_iter->GetBinError(j),2));
//          if(xi_right_unfolded_per_iter_pythia->GetBinError(j)>0) chi2_unfolded_xi_pythia += pow(xi_right_unfolded_per_iter_pythia->GetBinContent(j) - xi_gen_right_cut_pythia->GetBinContent(j), 2)/(pow(xi_right_unfolded_per_iter_pythia->GetBinError(j),2));
//     }
//     chi2_unfolded_vs_iter_xi->Fill(i, log10(chi2_unfolded_xi)); 
//     chi2_unfolded_vs_iter_xi_pythia->Fill(i, log10(chi2_unfolded_xi_pythia)); 

//     RooUnfoldBayes t_unfold_per_iter (&t_minus_response_pythia, t_rec_right_cut, i);
//     TH1F* t_right_unfolded_per_iter = (TH1F*)t_unfold_per_iter.Hreco();
//     RooUnfoldBayes t_unfold_per_iter_pythia (&t_minus_response, t_rec_right_cut_pythia, i);
//     TH1F* t_right_unfolded_per_iter_pythia = (TH1F*)t_unfold_per_iter_pythia.Hreco();
//     double chi2_unfolded_t = 0; 
//     double chi2_unfolded_t_pythia = 0; 
//     for (int j = 1; j<=t_right_unfolded->GetNbinsX(); ++j){
//          if(t_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_t += pow(t_right_unfolded_per_iter->GetBinContent(j) - t_gen_right_cut->GetBinContent(j), 2)/(pow(t_right_unfolded_per_iter->GetBinError(j),2));
//          if(t_right_unfolded_per_iter_pythia->GetBinError(j)>0) chi2_unfolded_t_pythia += pow(t_right_unfolded_per_iter_pythia->GetBinContent(j) - t_gen_right_cut_pythia->GetBinContent(j), 2)/(pow(t_right_unfolded_per_iter_pythia->GetBinError(j),2));
//     }
//     chi2_unfolded_vs_iter_t->Fill(i, log10(chi2_unfolded_t)); 
//     chi2_unfolded_vs_iter_t_pythia->Fill(i, log10(chi2_unfolded_t_pythia)); 
// }



// TCanvas *c1 = new TCanvas("c1","xi_right");
//     // RooUnfoldBayes unfold_right_xi_test (&xi_minus_response_pythia, xi_rec_right_cut, 4);
//     // TH1F* xi_unfolded_test = (TH1F*)unfold_right_xi_test.Hreco();
//     // xi_unfolded_test->SetXTitle("#xi");;
//     // xi_unfolded_test->SetYTitle("#frac{d#sigma}{d#xi} (nb})");;
//     // xi_unfolded_test->Scale(1/luminosity,"width");
//     // xi_unfolded_test->Draw();
//     // xi_gen_right->Scale(1/luminosity,"width");
//     // xi_gen_right->Draw("histsame");
//    chi2_unfolded_vs_iter_xi->SetXTitle("N_{iter}");
//    chi2_unfolded_vs_iter_xi->SetYTitle("log_{10} #chi^{2}_{unfolded}");
//    chi2_unfolded_vs_iter_xi->SetMarkerColor(1);
//    chi2_unfolded_vs_iter_xi->SetMarkerSize(0.3);
//    chi2_unfolded_vs_iter_xi->SetMarkerStyle(20);
//    chi2_unfolded_vs_iter_xi->Draw();

// TCanvas *c2 = new TCanvas("c2","xi_pythia");
//     // RooUnfoldBayes unfold_right_xi_test_pythia (&xi_minus_response, xi_rec_right_cut_pythia, 4);
//     // TH1F* xi_unfolded_test_pythia = (TH1F*)unfold_right_xi_test_pythia.Hreco();
//     // xi_unfolded_test_pythia->SetXTitle("#xi");;
//     // xi_unfolded_test_pythia->SetYTitle("#frac{d#sigma}{d#xi} (nb})");;
//     // xi_unfolded_test_pythia->Scale(1/luminosity,"width");
//     // xi_unfolded_test_pythia->Draw();
//     // xi_gen_right_cut_pythia->Scale(1/luminosity,"width");
//     // xi_gen_right_cut_pythia->Draw("histsame");
//    chi2_unfolded_vs_iter_xi_pythia->SetXTitle("N_{iter}");
//    chi2_unfolded_vs_iter_xi_pythia->SetYTitle("log_{10} #chi^{2}_{unfolded}");
//    chi2_unfolded_vs_iter_xi_pythia->SetMarkerColor(1);
//    chi2_unfolded_vs_iter_xi_pythia->SetMarkerSize(0.3);
//    chi2_unfolded_vs_iter_xi_pythia->SetMarkerStyle(20);
//    chi2_unfolded_vs_iter_xi_pythia->Draw();

// TCanvas *c3 = new TCanvas("c3","logx_right");
//     // RooUnfoldBayes unfold_right_logx_test (&logx_minus_response_pythia, log_x_rec_right_cut, 4);
//     // TH1F* logx_unfolded_test = (TH1F*)unfold_right_logx_test.Hreco();
//     // logx_unfolded_test->SetXTitle("log_{10} x");;
//     // logx_unfolded_test->SetYTitle("#frac{d#sigma}{dx} (nb})");;
//     // logx_unfolded_test->Scale(1/luminosity,"width");
//     // logx_unfolded_test->Draw();
//     // log_x_gen_right_cut->Scale(1/luminosity,"width");
//     // log_x_gen_right_cut->Draw("histsame");
//    chi2_unfolded_vs_iter_logx->SetXTitle("N_{iter}");
//    chi2_unfolded_vs_iter_logx->SetYTitle("log_{10} #chi^{2}_{unfolded}");
//    chi2_unfolded_vs_iter_logx->SetMarkerColor(1);
//    chi2_unfolded_vs_iter_logx->SetMarkerSize(0.3);
//    chi2_unfolded_vs_iter_logx->SetMarkerStyle(20);
//    chi2_unfolded_vs_iter_logx->Draw();

// TCanvas *c4 = new TCanvas("c4","logx_pythia");
//     // RooUnfoldBayes unfold_right_logx_test_pythia (&logx_minus_response, logx_rec_right_cut_pythia, 4);
//     // TH1F* logx_unfolded_test_pythia = (TH1F*)unfold_right_logx_test_pythia.Hreco();
//     // logx_unfolded_test_pythia->SetXTitle("log_{10} x");;
//     // logx_unfolded_test_pythia->SetYTitle("#frac{d#sigma}{dx} (nb})");;
//     // logx_unfolded_test_pythia->Scale(1/luminosity,"width");
//     // logx_unfolded_test_pythia->Draw();
//     // logx_gen_right_cut_pythia->Scale(1/luminosity,"width");
//     // logx_gen_right_cut_pythia->Draw("histsame");
//    chi2_unfolded_vs_iter_logx_pythia->SetXTitle("N_{iter}");
//    chi2_unfolded_vs_iter_logx_pythia->SetYTitle("log_{10} #chi^{2}_{unfolded}");
//    chi2_unfolded_vs_iter_logx_pythia->SetMarkerColor(1);
//    chi2_unfolded_vs_iter_logx_pythia->SetMarkerSize(0.3);
//    chi2_unfolded_vs_iter_logx_pythia->SetMarkerStyle(20);
//    chi2_unfolded_vs_iter_logx_pythia->Draw();

// TCanvas *c4 = new TCanvas("c4","t_right");
//    chi2_unfolded_vs_iter_t->SetXTitle("N_{iter}");
//    chi2_unfolded_vs_iter_t->SetYTitle("log_{10} #chi^{2}_{unfolded}");
//    chi2_unfolded_vs_iter_t->SetMarkerColor(1);
//    chi2_unfolded_vs_iter_t->SetMarkerSize(0.3);
//    chi2_unfolded_vs_iter_t->SetMarkerStyle(20);
//    chi2_unfolded_vs_iter_t->Draw();

// TCanvas *c5 = new TCanvas("c5","t_pythia");
//    chi2_unfolded_vs_iter_t_pythia->SetXTitle("N_{iter}");
//    chi2_unfolded_vs_iter_t_pythia->SetYTitle("log_{10} #chi^{2}_{unfolded}");
//    chi2_unfolded_vs_iter_t_pythia->SetMarkerColor(1);
//    chi2_unfolded_vs_iter_t_pythia->SetMarkerSize(0.3);
//    chi2_unfolded_vs_iter_t_pythia->SetMarkerStyle(20);
//    chi2_unfolded_vs_iter_t_pythia->Draw();

//     // RooUnfoldBayes unfold_right_t_test (&t_minus_response, t_rec_right_cut_pythia, 100);
//     // TH1F* t_right_unfolded_rec = (TH1F*)unfold_right_t_test.Hreco(); 
//     // t_right_unfolded_rec->SetXTitle("|t| (GeV^{2})");
//     // t_right_unfolded_rec->SetYTitle("#frac{d#sigma}{d|t|} (nb/GeV^{2})");
//     // t_right_unfolded_rec->Scale(luminosity,"width");
//     // // t_right_unfolded_rec->Draw("E1");
//     // t_gen_right_cut_pythia->Scale(luminosity,"width");
//     // // t_gen_right_cut_pythia->Draw("E1same");


////////////////////////////////// data chi2 - bottom line test /////////////////////////////////////////



// TH2F* chi2_unfolded_vs_iter_t = new TH2F("","",100,0,40,100,0,3);
// TH2F* chi2_unfolded_vs_iter_t_left = new TH2F("","",100,0,40,100,0,3);
// TH2F* chi2_unfolded_vs_iter_xi = new TH2F("","",100,0,40,100,0,3);
// TH2F* chi2_unfolded_vs_iter_xi_left = new TH2F("","",100,0,40,100,0,3);
// TH2F* chi2_unfolded_vs_iter_logx = new TH2F("","",100,0,40,100,0,3);
// TH2F* chi2_unfolded_vs_iter_logx_left = new TH2F("","",100,0,40,100,0,3);
// TH2F* delta_chi2_unfolded_vs_iter_t = new TH2F("","",100,0,40,100,0,3);
// TH2F* delta_chi2_unfolded_vs_iter_t_left = new TH2F("","",100,0,40,100,0,3);
// TH2F* delta_chi2_unfolded_vs_iter_xi = new TH2F("","",100,0,40,100,0,3);
// TH2F* delta_chi2_unfolded_vs_iter_xi_left = new TH2F("","",100,0,40,100,0,3);
// TH2F* delta_chi2_unfolded_vs_iter_logx = new TH2F("","",100,0,40,100,0,3);
// TH2F* delta_chi2_unfolded_vs_iter_logx_left = new TH2F("","",100,0,40,100,0,3);
// double chi2_unfolded_t [40]; 
// double chi2_unfolded_t_left[40]; 
// double chi2_unfolded_xi[40]; 
// double chi2_unfolded_xi_left[40]; 
// double chi2_unfolded_logx[40]; 
// double chi2_unfolded_logx_left[40]; 

// for (int i = 1; i<40; ++i){
//     RooUnfoldBayes t_right_unfold_per_iter (&t_minus_response, t_right_cut, i);
//     TH1F* t_right_unfolded_per_iter = (TH1F*)t_right_unfold_per_iter.Hreco();
//     RooUnfoldBayes t_left_unfold_per_iter (&t_plus_response, t_left_cut, i);
//     TH1F* t_left_unfolded_per_iter = (TH1F*)t_left_unfold_per_iter.Hreco();

//     RooUnfoldBayes xi_right_unfold_per_iter (&xi_minus_response, xi_right_cut, i);
//     TH1F* xi_right_unfolded_per_iter = (TH1F*)xi_right_unfold_per_iter.Hreco();
//     RooUnfoldBayes xi_left_unfold_per_iter (&xi_plus_response, xi_left_cut, i);
//     TH1F* xi_left_unfolded_per_iter = (TH1F*)xi_left_unfold_per_iter.Hreco();

//     RooUnfoldBayes logx_right_unfold_per_iter (&logx_minus_response, log_x_right_cut, i);
//     TH1F* logx_right_unfolded_per_iter = (TH1F*)logx_right_unfold_per_iter.Hreco();
//     RooUnfoldBayes logx_left_unfold_per_iter (&logx_plus_response, log_x_left_cut, i);
//     TH1F* logx_left_unfolded_per_iter = (TH1F*)logx_left_unfold_per_iter.Hreco();

//     chi2_unfolded_t[i] = 0; 
//     chi2_unfolded_t_left[i] = 0; 
//     chi2_unfolded_xi[i] = 0; 
//     chi2_unfolded_xi_left[i] = 0; 
//     chi2_unfolded_logx[i] = 0; 
//     chi2_unfolded_logx_left[i] = 0; 
//     for (int j = 1; j<=t_left_unfolded->GetNbinsX(); ++j){
//          if(t_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_t[i] += pow(t_right_unfolded_per_iter->GetBinContent(j) - t_gen_right_cut->GetBinContent(j), 2)/(pow(t_right_unfolded_per_iter->GetBinError(j),2));
//          if(t_left_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_t_left[i] += pow(t_left_unfolded_per_iter->GetBinContent(j) - t_gen_left_cut->GetBinContent(j), 2)/(pow(t_left_unfolded_per_iter->GetBinError(j),2));
//     } 
//     for (int j = 1; j<=xi_left_unfolded->GetNbinsX(); ++j){
//          if(xi_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_xi[i] += pow(xi_right_unfolded_per_iter->GetBinContent(j) - xi_gen_right->GetBinContent(j), 2)/(pow(xi_right_unfolded_per_iter->GetBinError(j),2));
//          if(xi_left_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_xi_left[i] += pow(xi_left_unfolded_per_iter->GetBinContent(j) - xi_gen_left->GetBinContent(j), 2)/(pow(xi_left_unfolded_per_iter->GetBinError(j),2));
//     }
//     for (int j = 1; j<=logx_left_unfolded->GetNbinsX(); ++j){
//          if(logx_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_logx[i] += pow(logx_right_unfolded_per_iter->GetBinContent(j) - log_x_gen_right_cut->GetBinContent(j), 2)/(pow(logx_right_unfolded_per_iter->GetBinError(j),2));
//          if(logx_left_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_logx_left[i] += pow(logx_left_unfolded_per_iter->GetBinContent(j) - log_x_gen_left_cut->GetBinContent(j), 2)/(pow(logx_left_unfolded_per_iter->GetBinError(j),2));
//     }
//     chi2_unfolded_vs_iter_t->Fill(i, log10(chi2_unfolded_t[i])); 
//     chi2_unfolded_vs_iter_t_left->Fill(i, log10(chi2_unfolded_t_left[i])); 
//     chi2_unfolded_vs_iter_xi->Fill(i, log10(chi2_unfolded_xi[i])); 
//     chi2_unfolded_vs_iter_xi_left->Fill(i, log10(chi2_unfolded_xi_left[i])); 
//     chi2_unfolded_vs_iter_logx->Fill(i, log10(chi2_unfolded_logx[i])); 
//     chi2_unfolded_vs_iter_logx_left->Fill(i, log10(chi2_unfolded_logx_left[i])); 
// }

// double chi2_smeared_t = 0;
// double chi2_smeared_t_left = 0;
// double chi2_smeared_xi = 0;
// double chi2_smeared_xi_left = 0;
// double chi2_smeared_logx = 0;
// double chi2_smeared_logx_left = 0;
//     for (int j = 1; j<=t_left_cut->GetNbinsX(); ++j){
//          if(t_right_cut->GetBinError(j)>0) chi2_smeared_t += pow(t_right_cut->GetBinContent(j) - t_rec_right_cut->GetBinContent(j), 2)/(pow(t_right_cut->GetBinError(j),2));
//          if(t_left_cut->GetBinError(j)>0) chi2_smeared_t_left += pow(t_left_cut->GetBinContent(j) - t_rec_left_cut->GetBinContent(j), 2)/(pow(t_left_cut->GetBinError(j),2));
//     } 
//     for (int j = 1; j<=xi_left_cut->GetNbinsX(); ++j){
//          if(xi_right_cut->GetBinError(j)>0) chi2_smeared_xi += pow(xi_right_cut->GetBinContent(j) - xi_rec_right_cut->GetBinContent(j), 2)/(pow(xi_right_cut->GetBinError(j),2));
//          if(xi_left_cut->GetBinError(j)>0) chi2_smeared_xi_left += pow(xi_left_cut->GetBinContent(j) - xi_rec_left_cut->GetBinContent(j), 2)/(pow(xi_left_cut->GetBinError(j),2));
//     }
//     for (int j = 1; j<=log_x_left_cut->GetNbinsX(); ++j){
//          if(log_x_right_cut->GetBinError(j)>0) chi2_smeared_logx += pow(log_x_right_cut->GetBinContent(j) - log_x_rec_right_cut->GetBinContent(j), 2)/(pow(log_x_right_cut->GetBinError(j),2));
//          if(log_x_left_cut->GetBinError(j)>0) chi2_smeared_logx_left += pow(log_x_left_cut->GetBinContent(j) - log_x_rec_left_cut->GetBinContent(j), 2)/(pow(log_x_left_cut->GetBinError(j),2));
//     }



// TCanvas *c1 = new TCanvas("c1","t_right");
// chi2_unfolded_vs_iter_t->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_t->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_t->SetMarkerColor(1);
// chi2_unfolded_vs_iter_t->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_t->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_t->Draw();
// TLine *line_t = new TLine(0,log10(chi2_smeared_t),40,log10(chi2_smeared_t));
// line_t->SetLineColor(kRed);
// line_t->Draw();
// TPad *text_line_t = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_t->SetFillColor(0);
//  // text_lixne_t->Draw();

// TCanvas *c2 = new TCanvas("c2","t_left");
// chi2_unfolded_vs_iter_t_left->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_t_left->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_t_left->SetMarkerColor(1);
// chi2_unfolded_vs_iter_t_left->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_t_left->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_t_left->Draw();
// TLine *line_t_left = new TLine(0,log10(chi2_smeared_t_left),40,log10(chi2_smeared_t_left));
// line_t_left->SetLineColor(kRed);
// line_t_left->Draw();
// TPad *text_line_t_left = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_t_left->SetFillColor(0);
//  // text_line_t_left->Draw();


// TCanvas *c3 = new TCanvas("c3","xi_right");
// chi2_unfolded_vs_iter_xi->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_xi->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_xi->SetMarkerColor(1);
// chi2_unfolded_vs_iter_xi->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_xi->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_xi->Draw();
// TLine *line_xi = new TLine(0,log10(chi2_smeared_xi),40,log10(chi2_smeared_xi));
// line_xi->SetLineColor(kRed);
// line_xi->Draw();
// TPad *text_line_xi = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_xi->SetFillColor(0);
//  // text_line_t->Draw();

// TCanvas *c4 = new TCanvas("c4","xi_left");
// chi2_unfolded_vs_iter_xi_left->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_xi_left->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_xi_left->SetMarkerColor(1);
// chi2_unfolded_vs_iter_xi_left->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_xi_left->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_xi_left->Draw();
// TLine *line_xi_left = new TLine(0,log10(chi2_smeared_xi_left),40,log10(chi2_smeared_xi_left));
// line_xi_left->SetLineColor(kRed);
// line_xi_left->Draw();
// TPad *text_line_xi_left = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_xi_left->SetFillColor(0);
//  // text_line_t_left->Draw();


// TCanvas *c5 = new TCanvas("c5","logx_right");
// chi2_unfolded_vs_iter_logx->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_logx->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_logx->SetMarkerColor(1);
// chi2_unfolded_vs_iter_logx->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_logx->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_logx->Draw();
// TLine *line_logx = new TLine(0,log10(chi2_smeared_logx),40,log10(chi2_smeared_logx));
// line_logx->SetLineColor(kRed);
// line_logx->Draw();
// TPad *text_line_logx = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_logx->SetFillColor(0);
//  // text_line_t->Draw();

// TCanvas *c6 = new TCanvas("c6","logx_left");
// chi2_unfolded_vs_iter_logx_left->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_logx_left->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_logx_left->SetMarkerColor(1);
// chi2_unfolded_vs_iter_logx_left->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_logx_left->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_logx_left->Draw();
// TLine *line_logx_left = new TLine(0,log10(chi2_smeared_logx_left),40,log10(chi2_smeared_logx_left));
// line_logx_left->SetLineColor(kRed);
// line_logx_left->Draw();
// TPad *text_line_logx_left = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_logx_left->SetFillColor(0);
//  // text_line_t_left->Draw();



// for (int i = 2; i<40; ++i){
//     if(chi2_unfolded_t[i]!=0) delta_chi2_unfolded_vs_iter_t->Fill(i, (log10(chi2_unfolded_t[i])-log10(chi2_unfolded_t[i-1]))/log10(chi2_unfolded_t[i-1])); cout<<log10(chi2_unfolded_t[i-1])<<"  "<<log10(chi2_unfolded_t[i])<<"  "<<(chi2_unfolded_t[i]-chi2_unfolded_t[i-1])/chi2_unfolded_t[i-1]<<endl;
//     // if(chi2_unfolded_t_left[i]!=0) delta_chi2_unfolded_vs_iter_t_left->Fill(i, (chi2_unfolded_t_left[i+1]-chi2_unfolded_t_left[i])/chi2_unfolded_t_left[i]); 
//     // if(chi2_unfolded_xi[i]!=0) delta_chi2_unfolded_vs_iter_xi->Fill(i, (chi2_unfolded_xi[i+1]-chi2_unfolded_xi[i])/chi2_unfolded_xi[i]);
//     // if(chi2_unfolded_xi_left[i]!=0) delta_chi2_unfolded_vs_iter_xi_left->Fill(i, (chi2_unfolded_xi_left[i+1]-chi2_unfolded_xi_left[i])/chi2_unfolded_xi_left[i]);
//     // if(chi2_unfolded_logx[i]!=0) delta_chi2_unfolded_vs_iter_logx->Fill(i, (chi2_unfolded_logx[i+1]-chi2_unfolded_logx[i])/chi2_unfolded_logx[i]);
//     // if(chi2_unfolded_logx_left[i]!=0) delta_chi2_unfolded_vs_iter_logx_left->Fill(i, (chi2_unfolded_logx_left[i+1]-chi2_unfolded_logx_left[i])/chi2_unfolded_logx_left[i]);
// }

// TCanvas *c7 = new TCanvas("c7","t_right");
// delta_chi2_unfolded_vs_iter_t->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_t->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_t->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_t->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_t->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_t->Draw();

// TCanvas *c8 = new TCanvas("c8","t_left");
// delta_chi2_unfolded_vs_iter_t_left->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_t_left->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_t_left->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_t_left->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_t_left->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_t_left->Draw();


// TCanvas *c9 = new TCanvas("c9","xi_right");
// delta_chi2_unfolded_vs_iter_xi->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_xi->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_xi->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_xi->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_xi->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_xi->Draw();

// TCanvas *c10 = new TCanvas("c10","xi_left");
// delta_chi2_unfolded_vs_iter_xi_left->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_xi_left->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_xi_left->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_xi_left->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_xi_left->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_xi_left->Draw();


// TCanvas *c11 = new TCanvas("c11","logx_right");
// delta_chi2_unfolded_vs_iter_logx->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_logx->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_logx->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_logx->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_logx->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_logx->Draw();

// TCanvas *c12 = new TCanvas("c12","logx_left");
// delta_chi2_unfolded_vs_iter_logx_left->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_logx_left->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_logx_left->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_logx_left->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_logx_left->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_logx_left->Draw();




// // if(reweight_slope){
//    TFile * outfile = new TFile("unfold_vs_iter_data_pomwig.root", "RECREATE");
//    TH1F* xi_right_unfolded_per_iter[21];
//    TH1F* xi_left_unfolded_per_iter[21];
//    TH1F* logx_right_unfolded_per_iter[21];
//    TH1F* logx_left_unfolded_per_iter[21];
//    TH1F* t_right_unfolded_per_iter[21];
//    TH1F* t_left_unfolded_per_iter[21];

// // //    // for (int j = 1; j<=t_left_unfolded->GetNbinsX(); ++j){
// // //    // t_error_per_iter_left[j] = new TH1F("","",100,0,100);
// // //    //     for (int i = 1; i<50; ++i){
// // //    //          RooUnfoldBayes t_left_unfold_per_iter (&t_plus_response, t_left_cut, i);
// // //    //          TH1F* t_left_unfolded_per_iter = (TH1F*)t_left_unfold_per_iter.Hreco();
// // //    //          t_error_per_iter_left[j]->Fill(i, t_left_unfolded_per_iter->GetBinContent(j));
// // //    //      }
// // //    //  }  
// // TH2F* chi2_unfolded_vs_iter = new TH2F("","",100,0,50,100,0,3);
// // TH2F* chi2_unfolded_vs_iter_left = new TH2F("","",100,0,50,100,0,3);
//    for (int i = 1; i<22; ++i){
//         RooUnfoldBayes xi_right_unfold_per_iter (&xi_minus_response, xi_right_cut, i);
//         xi_right_unfolded_per_iter[i] = (TH1F*) xi_right_unfold_per_iter.Hreco();
//         xi_right_unfolded_per_iter[i]->Scale(1/luminosity,"width");
//         xi_right_unfolded_per_iter[i]->GetXaxis()->SetTitle("#xi");
//         xi_right_unfolded_per_iter[i]->GetYaxis()->SetTitle("#frac{d#sigma}{d#xi}(nb)");

//         RooUnfoldBayes xi_left_unfold_per_iter (&xi_plus_response, xi_left_cut, i);
//         xi_left_unfolded_per_iter[i] = (TH1F*) xi_left_unfold_per_iter.Hreco();
//         xi_left_unfolded_per_iter[i]->Scale(1/luminosity,"width");
//         xi_left_unfolded_per_iter[i]->GetXaxis()->SetTitle("#xi");
//         xi_left_unfolded_per_iter[i]->GetYaxis()->SetTitle("#frac{d#sigma}{d#xi}(nb)");

//         RooUnfoldBayes logx_right_unfold_per_iter (&logx_minus_response, log_x_right_cut, i);
//         logx_right_unfolded_per_iter[i] = (TH1F*) logx_right_unfold_per_iter.Hreco();
//         logx_right_unfolded_per_iter[i]->Scale(1/luminosity,"width");
//         logx_right_unfolded_per_iter[i]->GetXaxis()->SetTitle("log_{10} x");
//         logx_right_unfolded_per_iter[i]->GetYaxis()->SetTitle("#frac{d#sigma}{dx}(nb)");

//         RooUnfoldBayes logx_left_unfold_per_iter (&logx_plus_response, log_x_left_cut, i);
//         logx_left_unfolded_per_iter[i] = (TH1F*) logx_left_unfold_per_iter.Hreco();
//         logx_left_unfolded_per_iter[i]->Scale(1/luminosity,"width");
//         logx_left_unfolded_per_iter[i]->GetXaxis()->SetTitle("log_{10} x");
//         logx_left_unfolded_per_iter[i]->GetYaxis()->SetTitle("#frac{d#sigma}{dx}(nb)");

//            RooUnfoldBayes t_right_unfold_per_iter (&t_minus_response, t_right_cut, i);
//            t_right_unfolded_per_iter[i] = (TH1F*) t_right_unfold_per_iter.Hreco();
//            t_right_unfolded_per_iter[i]->Scale(1/luminosity,"width");
//            t_right_unfolded_per_iter[i]->GetXaxis()->SetTitle("|t| (GeV^{2})");
//            t_right_unfolded_per_iter[i]->GetYaxis()->SetTitle("#frac{d#sigma}{d|t|}(nb/GeV^{2})");

//            RooUnfoldBayes t_left_unfold_per_iter (&t_plus_response, t_left_cut, i);
//            t_left_unfolded_per_iter[i] = (TH1F*) t_left_unfold_per_iter.Hreco();
//            t_left_unfolded_per_iter[i]->Scale(1/luminosity,"width");
//            t_left_unfolded_per_iter[i]->GetXaxis()->SetTitle("|t| (GeV^{2})");
//            t_left_unfolded_per_iter[i]->GetYaxis()->SetTitle("#frac{d#sigma}{d|t|}(nb/GeV^{2})");
          
// //            // RooUnfoldBayes xi_right_unfold_per_iter (&t_minus_response_pythia, t_rec_right_cut, i);
// //            // xi_right_unfolded_per_iter[i] = (TH1F*) t_right_unfold_per_iter.Hreco();
// //            // xi_right_unfolded_per_iter[i]->Scale(1/luminosity,"width");
// //            // xi_right_unfolded_per_iter[i]->GetXaxis()->SetTitle("|t| (GeV^{2})");
// //            // xi_right_unfolded_per_iter[i]->GetYaxis()->SetTitle("#frac{d#sigma}{d|t|}(nb/GeV^{2})");
// //            // RooUnfoldBayes t_right_unfold_per_iter_pythia (&t_minus_response, t_rec_right_cut_pythia, i);
// //            // t_right_unfolded_per_iter_pythia[i] = (TH1F*) t_right_unfold_per_iter_pythia.Hreco();
// //            // t_right_unfolded_per_iter_pythia[i]->Scale(1/luminosity,"width");
// //            // t_right_unfolded_per_iter_pythia[i]->GetXaxis()->SetTitle("|t| (GeV^{2})");
// //            // t_right_unfolded_per_iter_pythia[i]->GetYaxis()->SetTitle("#frac{d#sigma}{d|t|}(nb/GeV^{2})");
//    }
//  outfile->Write();
// // }
// // // //     t_right_unfolded_per_iter[i] = (TH1F*)t_right_unfold_per_iter.Hreco();
// // // // // //     TH1F* t_right_unfolded_per_iter_pythia = (TH1F*)t_right_unfold_per_iter_pythia.Hreco();
// // // // // //     // double_t chi2_unfolded = t_gen_right_cut->Chi2Test(t_right_unfolded_per_iter,"UW NORM P CHI2",0); 
// // // // // //     // chi2_unfolded_vs_iter->Fill(i, log10(chi2_unfolded));
// // //     RooUnfoldBayes logx_left_unfold_per_iter (&logx_plus_response, log_x_left_cut, i);
// // //     TH1F* logx_left_unfolded_per_iter = (TH1F*)logx_left_unfold_per_iter.Hreco();

// // // //          chi2_unfolded_pythia += pow(t_right_unfolded_per_iter_pythia->GetBinContent(j) - t_gen_right_cut_pythia->GetBinContent(j), 2)/pow(t_right_unfolded_per_iter_pythia->GetBinError(j),2);
// //     // t_right_unfolded_per_iter[i] = (TH1F*)t_right_unfold_per_iter.Hreco();
// //     t_right_unfolded_per_iter[i]->Scale(1/luminosity,"width");
// //     t_right_unfolded_per_iter[i]->GetXaxis()->SetTitle("|t|(GeV^{2})");
// //     t_right_unfolded_per_iter[i]->GetYaxis()->SetTitle("#frac{d#sigma}{d|t|}(nb/GeV^{2})");


// // TCanvas *c3 = new TCanvas("c3","xi_left");
 
// // t_gen_right_cut->Scale(1,"width");
// // t_gen_right_cut->Draw("hist");
// // t_right_unfolded->Scale(1,"width");
// // t_right_unfolded->Draw("E1same");
// // // t_right_unfolded_per_iter_2 ->Scale(1,"width");
// // // t_right_unfolded_per_iter_2 ->Draw("E1same");

// // TCanvas *c4 = new TCanvas("c4","xi_left");
// // t_gen_left_cut->Scale(1,"width");
// // t_gen_left_cut->Draw("E1same");
// // t_left_unfolded->Scale(1,"width");
// // t_left_unfolded->Draw("E1same");

// }


// if (!reweight_slope){
// TCanvas *c3 = new TCanvas("c3","xi_left");
//     RooUnfoldBayes t_right_unfold_per_iter_2 (&t_plus_response, t_left_cut, 4);
//     TH1F* t_right_unfolded_per_iter_2 = (TH1F*)t_right_unfold_per_iter_2.Hreco();

// t_gen_right_cut->Add(t_gen_left_cut);
// t_gen_right_cut->Scale(0.5,"width");
// t_gen_right_cut->Draw("hist");
// t_gen_left_cut->Scale(1,"width");
// t_gen_left_cut->Draw("E1same");
//}

// t_left_unfolded->Scale(1,"width");
// t_right_unfolded->Add(t_left_unfolded);
// t_right_unfolded->Scale(0.5);
// t_right_unfolded->Draw("E1");
// // t_left_unfolded->Draw("E1same");


/////////////////////////////// error vs iter ////////////////////////////////////////
   //      for (int j = 1; j<=t_right_unfolded->GetNbinsX(); ++j){
   //      for (int i = 1; i<21; ++i){
   //           RooUnfoldBayes t_right_unfold_per_iter (&t_minus_response_pythia, t_rec_right_cut, i);
   //           TH1F* t_right_unfolded_per_iter = (TH1F*)t_right_unfold_per_iter.Hreco();
   //           error_vs_iter->Fill(i, t_right_unfolded_per_iter->GetBinError(j));
   //      }
   // }





/////////////////////////////// closure test chi2 //////////////////////////////////////////
   // TH2F* chi2_unfolded_vs_iter = new TH2F("","",100,0,100,100,-2,3);
   // double chi2_unfolded[21]; 
   // for (int i = 1; i<21; ++i){
   //      RooUnfoldBayes t_right_unfold_per_iter (&t_minus_response, t_rec_right_cut, i);
   //      TH1F* t_right_unfolded_per_iter = (TH1F*)t_right_unfold_per_iter.Hreco();
   //      chi2_unfolded[i] = 0;
   //      for (int j = 1; j<=t_right_unfolded->GetNbinsX(); ++j){
   //           chi2_unfolded[i] += pow(t_right_unfolded_per_iter->GetBinContent(j) - t_gen_right_cut->GetBinContent(j), 2)/pow(t_right_unfolded_per_iter->GetBinError(j),2);
   //      }
   //  }
   // for (int i = 1; i<21; ++i){
   
   //      chi2_unfolded_vs_iter->Fill(i, (chi2_unfolded[i]-chi2_unfolded[i-1])/chi2_unfolded[i]);
   //  }
   //     chi2_unfolded_vs_iter->SetXTitle("N_{iter}");
   //     chi2_unfolded_vs_iter->SetYTitle("#Delta#chi^{2}/#chi^{2}");
   //     chi2_unfolded_vs_iter->SetMarkerColor(1);
   //     chi2_unfolded_vs_iter->SetMarkerSize(0.3);
   //     chi2_unfolded_vs_iter->SetMarkerStyle(20);
   //     chi2_unfolded_vs_iter->Draw();






//     t_right_unfolded_per_iter->GetXaxis()->SetTitle("|t|(GeV^{2})");
//     t_right_unfolded_per_iter->GetYaxis()->SetTitle("#frac{d#sigma}{d|t|}(nb/GeV^{2})");
//     t_right_unfolded_per_iter->Scale(1/luminosity,"width");
//     t_right_unfolded_per_iter->Draw("E1");
//     t_gen_right_cut_pythia->Scale(1/luminosity,"width");
//     t_gen_right_cut_pythia->Draw("histsame");
// TLegend *leg1_1 = new TLegend(0.2,0.75,0.48,0.9);
//         leg1_1->AddEntry(t_right_unfolded_per_iter,"PYTHIA8 unfolded with POMWIG","lp");
//         leg1_1->AddEntry(t_gen_right_cut_pythia,"PYTHIA8","lp");
//         leg1_1->SetFillColor(0);
//         leg1_1->SetLineColor(0);
//         leg1_1->Draw();  


//     RooUnfoldBayes t_right_unfold_per_iter (&t_minus_response_pythia, t_rec_right_cut, 4);
//     TH1F* t_right_unfolded_per_iter = (TH1F*)t_right_unfold_per_iter.Hreco();
//     t_right_unfolded_per_iter->Draw("E1");    
//     t_gen_right_cut->Draw("E1same");    
// TLegend *leg1_1 = new TLegend(0.2,0.75,0.48,0.9);
//         leg1_1->AddEntry(t_right_unfolded_per_iter,"POMWIG unfolded with PYTHIA8","lp");
//         leg1_1->AddEntry(t_gen_right_cut,"POMWIG","lp");
//         leg1_1->SetFillColor(0);
//         leg1_1->SetLineColor(0);
//         leg1_1->Draw();  

// chi2_unfolded_vs_iter_pythia->SetMarkerColor(2);
// chi2_unfolded_vs_iter_pythia->SetMarkerSize(0.3);
// chi2_unfolded_vs_iter_pythia->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_pythia->Draw("same");
// TLegend *leg1_1 = new TLegend(0.2,0.75,0.48,0.9);
//         leg1_1->AddEntry(chi2_unfolded_vs_iter,"POMWIG","p");
//         leg1_1->AddEntry(chi2_unfolded_vs_iter_pythia,"PYTHIA8","p");
//         leg1_1->SetFillColor(0);
//         leg1_1->SetLineColor(0);
//         leg1_1->Draw();  

// // // // double_t chi2_unfolded = t_right_unfolded->Chi2Test(t_gen_right_cut,"UW NORM P",0); 
// double_t chi2_smeared = t_rec_right_cut->Chi2Test(t_right_cut,"UW P CHI2",0); cout << "CHI2 smeared" << log10(chi2_smeared) << endl;
// double_t chi2_unfolded = t_right_unfolded->Chi2Test(t_gen_right_cut,"UW P CHI2",0); cout << "CHI2 unfolded" << log10(chi2_unfolded) << endl;
// t_right_cut->Scale(1,"width");
// t_right_cut->Draw("E1");
// t_rec_right_cut->Scale(1,"width");
// t_rec_right_cut->Draw("histsame");
// t_right_unfolded->Scale(1,"width");
// t_right_unfolded->Draw("E1same");
// t_gen_right_cut->Scale(1,"width");
// t_gen_right_cut->Draw("histsame");


 
/////////////////////////two sectors comparison /////////////////////////////////////
//    TCanvas *c1 = new TCanvas("c1","sigma_t_averaged",0,0,600,500);
//    c1->Range(0,0,1,1);

//   logx_left_unfolded->Scale(1/luminosity,"width");
//   logx_right_unfolded->Scale(1/luminosity,"width");
   
// // Bottom plot
//   TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.99,0.32);
//   c1_1->Draw();
//   c1_1->cd();
//   c1_1->SetTopMargin(0.01);
//   c1_1->SetBottomMargin(0.3);
//   c1_1->SetRightMargin(0.1);
//   c1_1->SetFillStyle(0);
//   TH1F* sigma_t_45_vs_56 = new TH1F("","",15, -4,0);
//   double sector_diff[11];
//   double err_diff[11];

//   for (int i = 1; i<=logx_right_unfolded->GetNbinsX(); ++i){
//       sector_diff[i] = fabs(logx_left_unfolded->GetBinContent(i) - logx_right_unfolded->GetBinContent(i));
//       err_diff[i] = pow(logx_left_unfolded->GetBinError(i),2) + pow(logx_right_unfolded->GetBinError(i),2);
//       if(err_diff[i]>0) sigma_t_45_vs_56->SetBinContent(i, sector_diff[i]/sqrt(err_diff[i]));

//   }
//   sigma_t_45_vs_56->GetYaxis()->SetTitle("#frac{#left|#Delta#frac{d#sigma}{dx}#right|}{#sqrt{#sigma_{45}^{2}+#sigma_{56}^{2}}}");
//   sigma_t_45_vs_56->GetXaxis()->SetTitle("#xi");
//   sigma_t_45_vs_56->Draw("E1");

//   c1->cd();
//   TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.99,0.99);
//   c1_2->Draw(); 
//   c1_2->cd();
//   c1_2->SetTopMargin(0.1);
//   c1_2->SetBottomMargin(0.01);
//   c1_2->SetRightMargin(0.1);
//   c1_1->SetFillStyle(0);
//   c1_1->SetCanvasSize(600,600);
//   logx_left_unfolded->GetYaxis()->SetTitle("#frac{d#sigma}{dx}(nb)");
//   logx_left_unfolded->Draw("E1");
//   logx_right_unfolded->Draw("E1same");
//   TLegend *leg1_2 = new TLegend(0.2,0.75,0.48,0.9);
//   leg1_2->AddEntry(logx_left_unfolded,"Data - Sector 45","lp");
//   leg1_2->AddEntry(logx_right_unfolded,"Data - Sector 56","lp");
//   leg1_2->SetFillColor(0);
//   leg1_2->SetLineColor(0);
//   leg1_2->Draw();  

//   c1_2->SetLogy(1);


///////////////////////matrices///////////////////

  // TCanvas *c1_1 = new TCanvas("c1_1", "newpad");

   // TH2F* stat_cov_t_right = new TH2F("","",8,tbins,8,tbins); 
   // for (int i = 1; i<=t_right_cut->GetNbinsX(); ++i){
   //    stat_cov_t_right->SetBinContent(i,i, pow(t_right_cut->GetBinError(i),2));
   // }
   // stat_cov_t_right->GetXaxis()->SetTitle("|t| (GeV^{2})");
   // stat_cov_t_right->GetYaxis()->SetTitle("|t| (GeV^{2})");
   // stat_cov_t_right->Draw();

  // // // t_minus_th2->Scale(1/t_minus_th2->GetSumOfWeights());
//     TMatrix t_right_covariance = (TMatrix)unfold_right_t.Ereco();
//     TH2F* t_right_covariance_matrix_unfold = new TH2F(t_right_covariance);
//     TH2F* t_right_covariance_matrix = new TH2F("","",8,tbins,8,tbins);

//     for (int i = 1; i<=t_right_cut->GetNbinsX(); ++i){
//          for (int j = 1; j<=t_right_cut->GetNbinsX(); ++j){
//               t_right_covariance_matrix->SetBinContent(i,j,t_right_covariance_matrix_unfold->GetBinContent(i,j)); 
//         }
//     }    

// t_minus_th2->GetXaxis()->SetTitle("|t|_{rec} (GeV^{2})");
// t_minus_th2->GetYaxis()->SetTitle("|t|_{gen} (GeV^{2})");
   // t_right_covariance_matrix->GetXaxis()->SetTitle("|t| (GeV^{2})");
   // t_right_covariance_matrix->GetYaxis()->SetTitle("|t| (GeV^{2})");
   // t_right_covariance_matrix->GetZaxis()->SetTitle("|t|_{gen} (GeV^{-2})");
   // t_right_covariance_matrix->Draw();
// std::vector<double> sum_weights_t; sum_weights_t.clear(); sum_weights_t.resize( t_minus_th2->GetYaxis()->GetNbins() );
// for(int i=1; i <= t_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= t_minus_th2->GetYaxis()->GetNbins(); ++j){
//       sum_weights_t[j] += t_minus_th2->GetBinContent(i,j); 
//    } 
// }

// for(int i=1; i <= t_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= t_minus_th2->GetYaxis()->GetNbins(); ++j){
//       double val = t_minus_th2->GetBinContent(i,j); 
//       cout << "val: " << val << " - sum weights t: " << sum_weights_t[j] << endl;   
//       if( sum_weights_t[j] > 0. ) t_minus_th2->SetBinContent( i, j, val/sum_weights_t[j] ); 
//    } 
// }

 t_minus_th2->Draw("colz");

//   TCanvas *c1_2 = new TCanvas("c1_2", "newpad");
// // // xi_minus_th2->Scale(1/xi_minus_th2->GetSumOfWeights());
// // xi_minus_th2->GetXaxis()->SetTitle("#xi_{rec}");
// // xi_minus_th2->GetYaxis()->SetTitle("#xi_{gen}");
// // //xi_minus_th2->Draw("colz");
//     TMatrix xi_right_covariance = (TMatrix)unfold_right_xi.Ereco();
//     TH2F* xi_right_covariance_matrix_unfold = new TH2F(xi_right_covariance);
//     TH2F* xi_right_covariance_matrix = new TH2F("","",11,xi_bins,11,xi_bins);
//     for (int i = 1; i<=xi_right_cut->GetNbinsX(); ++i){
//          for (int j = 1; j<=xi_right_cut->GetNbinsX(); ++j){
//               xi_right_covariance_matrix->SetBinContent(i,j,xi_right_covariance_matrix_unfold->GetBinContent(i,j)); 
//         }
//     }    
//    xi_right_covariance_matrix->GetXaxis()->SetTitle("#xi");
//    xi_right_covariance_matrix->GetYaxis()->SetTitle("#xi");
//    xi_right_covariance_matrix->Draw();

// std::vector<double> sum_weights_xi; sum_weights_xi.clear(); sum_weights_xi.resize( xi_minus_th2->GetYaxis()->GetNbins() );
// for(int i=1; i <= xi_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= xi_minus_th2->GetYaxis()->GetNbins(); ++j){
//       sum_weights_xi[j] += xi_minus_th2->GetBinContent(i,j); 
//    } 
// }

// for(int i=1; i <= xi_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= xi_minus_th2->GetYaxis()->GetNbins(); ++j){
//       double val = xi_minus_th2->GetBinContent(i,j); 
//       if( sum_weights_xi[j] > 0. ) xi_minus_th2->SetBinContent( i, j, val/sum_weights_xi[j] ); 
//    } 
// }

// xi_minus_th2->Draw("colz");

//   TCanvas *c1_3 = new TCanvas("c1_3", "newpad");
// // // logx_minus_th2->Scale(1/logx_minus_th2->GetSumOfWeights());
// // logx_minus_th2->GetXaxis()->SetTitle("log_{10} x_{rec}");
// // logx_minus_th2->GetYaxis()->SetTitle("log_{10} x_{gen}");
//     TMatrix logx_right_covariance = (TMatrix)unfold_right_logx.Ereco();
//     TH2F* logx_right_covariance_matrix_unfold = new TH2F(logx_right_covariance);
//     TH2F* logx_right_covariance_matrix = new TH2F("","",15,-4,0,15,-4,0);
//     for (int i = 1; i<=log_x_right_cut->GetNbinsX(); ++i){
//          for (int j = 1; j<=log_x_right_cut->GetNbinsX(); ++j){
//               logx_right_covariance_matrix->SetBinContent(i,j,logx_right_covariance_matrix_unfold->GetBinContent(i,j)); 
//         }
//     }    
//    logx_right_covariance_matrix->GetXaxis()->SetTitle("log_{10} x");
//    logx_right_covariance_matrix->GetYaxis()->SetTitle("log_{10} x");
//    logx_right_covariance_matrix->Draw();

// std::vector<double> sum_weights_x; sum_weights_x.clear(); sum_weights_x.resize( logx_minus_th2->GetYaxis()->GetNbins() );
// for(int i=1; i <= logx_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= logx_minus_th2->GetYaxis()->GetNbins(); ++j){
//       sum_weights_x[j] += logx_minus_th2->GetBinContent(i,j); 
//    } 
// }

// for(int i=1; i <= logx_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= logx_minus_th2->GetYaxis()->GetNbins(); ++j){
//       double val = logx_minus_th2->GetBinContent(i,j); 
//       if( sum_weights_x[j] > 0. ) logx_minus_th2->SetBinContent( i, j, val/sum_weights_x[j] ); 
//    } 
// }
// logx_minus_th2->Draw("colz");


//   TCanvas *c1_3 = new TCanvas("c1_3", "newpad");
// rp_pos_rigth_xmax->Draw();
//   TCanvas *c1_1 = new TCanvas("c1_1", "newpad");
// rp_pos_left_xmax->Draw();

// t_right_unfolded_old->GetXaxis()->SetTitle("|t| (GeV^{2})");
// t_right_unfolded_old->GetYaxis()->SetTitle("d#sigma/d|t| (nb/GeV^{2})");
// t_right_unfolded_old->Scale(1/luminosity,"width");
// t_right_unfolded_old->Draw();
// t_right_unfolded->SetLineColor(2);
// t_right_unfolded->Scale(1/luminosity,"width");
// t_right_unfolded->Draw("E1same");
// // t_right_cut->Scale(1,"width");
// // t_right_cut->Draw("E1same");
//   TLegend *leg1_2 = new TLegend(0.2,0.75,0.48,0.9);
//   // leg1_2->AddEntry(t_right_cut,"Data","lp");
//   leg1_2->AddEntry(t_right_unfolded_old,"Data (t = -p^{2}(#theta_{x}^{2} + #theta_{y}^{2})))","lp");
//   leg1_2->AddEntry(t_right_unfolded,"Data (t = (p_{beam} - p_{scatt})^{2})","lp");
//   leg1_2->SetFillColor(0);
//   leg1_2->SetLineColor(0);
//   leg1_2->Draw();  


  //  t_left_cut->GetXaxis()->SetTitle("|t| (GeV^{2})");
  //  t_left_cut->GetYaxis()->SetTitle("d#sigma/d|t| (nb/GeV^{2})");
  //  t_left_cut->Scale(1,"width"); 
  //  t_left_cut->Draw("E1"); 
  //  t_rec_left_cut->SetLineColor(4); 
  //  t_rec_left_cut->Scale(1,"width"); 
  //  t_rec_left_cut->Draw("histsame"); 
  //  t_old_left->Scale(1,"width"); 
  //  t_old_left->Draw("histsame"); 
  // TLegend *leg1_2 = new TLegend(0.2,0.75,0.48,0.9);
  // leg1_2->AddEntry(t_left_cut,"Data","lp");
  // leg1_2->AddEntry(t_old_left,"POMWIG (t = -p^{2}(#theta_{x}^{2} + #theta_{y}^{2})))","lp");
  // leg1_2->AddEntry(t_rec_left_cut,"POMWIG (t = (p_{beam} - p_{scatt})^{2})","lp");
  // leg1_2->SetFillColor(0);
  // leg1_2->SetLineColor(0);
  // leg1_2->Draw();  


// p_right_pomwig->GetXaxis()->SetTitle("p (GeV)");
// p_right_pomwig->Scale(1/p_right_pomwig->GetSumOfWeights());
// p_right_pomwig->Draw();
// p_right_pythia->Scale(1/p_right_pythia->GetSumOfWeights());
// p_right_pythia->Draw("same");
// t_right_unfolded_herabackg->Scale(1,"width");
// t_right_unfolded_herabackg->Draw();
// t_right_unfolded->Scale(1,"width");
// t_right_unfolded->Draw("histsame");
// t_gen_right_cut->Add(t_gen_left_cut);
// t_gen_right_cut->Scale(0.5,"width");

//   TCanvas *c1 = new TCanvas("c1", "newpad");
// t_gen_left_pt20->SetLineColor(4);
// t_gen_left_pt20->Scale(1/t_gen_left_pt20->GetSumOfWeights());
// t_gen_left_pt20->GetXaxis()->SetTitle("|t| (GeV^{2})");
// t_gen_left_pt20->GetYaxis()->SetTitle("Events");
// t_gen_left_pt20->Draw();
// t_gen_left_pythia_jet20->Scale(1/t_gen_left_pythia_jet20->GetSumOfWeights());
// t_gen_left_pythia_jet20->Draw("same");
// log_t_vs_xi->GetXaxis()->SetTitle("log_{10}(|t|)");
// log_t_vs_xi->GetYaxis()->SetTitle("log_{10}(#xi)");
// log_t_vs_xi->Draw();

//   TCanvas *c2 = new TCanvas("c2", "newpad");
// t_gen_left_pt20_nomass->SetLineColor(4);
// t_gen_left_pt20_nomass->Scale(1/t_gen_left_pt20_nomass->GetSumOfWeights());
// t_gen_left_pt20_nomass->GetXaxis()->SetTitle("|t| (GeV^{2})");
// t_gen_left_pt20_nomass->GetYaxis()->SetTitle("Events");
// t_gen_left_pt20_nomass->Draw();
// t_gen_left_pythia_jet20_nomass->Scale(1/t_gen_left_pythia_jet20_nomass->GetSumOfWeights());
// t_gen_left_pythia_jet20_nomass->Draw("same");
// log_t_vs_xi_nomass->GetXaxis()->SetTitle("log_{10}(|t|)");
// log_t_vs_xi_nomass->GetYaxis()->SetTitle("log_{10}(#xi)");
// log_t_vs_xi_nomass->Draw();

//   TCanvas *c3 = new TCanvas("c3", "newpad");
// t_gen_left_pt40->SetLineColor(4);
// t_gen_left_pt40->Scale(1/t_gen_left_pt40->GetSumOfWeights());
// t_gen_left_pt40->Draw();
// t_gen_left_pythia_jet40->Scale(1/t_gen_left_pythia_jet40->GetSumOfWeights());
// t_gen_left_pythia_jet40->Draw("same");

//   TCanvas *c4 = new TCanvas("c4", "newpad");
// t_gen_left_pt50->SetLineColor(4);
// t_gen_left_pt50->Scale(1/t_gen_left_pt50->GetSumOfWeights());
// t_gen_left_pt50->Draw();
// t_gen_left_pythia_jet50->Scale(1/t_gen_left_pythia_jet50->GetSumOfWeights());
// t_gen_left_pythia_jet50->Draw("same");

// t_left_unfolded->Scale(1,"width");
// t_left_unfolded->Draw("E1");
// t_gen_right_cut->Scale(1,"width");
// t_gen_right_cut->Draw("histsame");
// log_t_vs_xi->Draw("colz");

}




