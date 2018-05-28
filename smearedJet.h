#include "Riostream.h"
#include <sstream> 

MyGenJet const* getMatch(MyBaseJet const& recoJet, vector<MyGenJet>& genJet, double resolution){
         double dRMin = 999;
         double dRMax = 0.25;
         double m_dPt_max_factor = 3;
         //double resolution = 1.;

         TLorentzVector v_genJet;
   	 TLorentzVector v_recoJet;
   	 v_recoJet.SetPxPyPzE(recoJet.Px(), recoJet.Py(), recoJet.Pz(), recoJet.E());
   	 MyGenJet const* matched_genJet = nullptr;
   	 //vector<MyGenJet> matched_genJet;
         //matched_genJet.resize( genJet.size() ); 
         size_t idx_jet = 0;

         for (vector<MyGenJet>::iterator it_genjet = genJet.begin(); it_genjet != genJet.end(); ++it_genjet, ++idx_jet){
   
	    v_genJet.SetPxPyPzE(it_genjet->Px(), it_genjet->Py(), it_genjet->Pz(), it_genjet->E());
   	    double dR = v_genJet.DeltaR(v_recoJet);
       	    double dPt = std::abs(it_genjet->Pt() - recoJet.Pt());

            if (dR > dRMin) continue;
   	    if (dR < dRMax) {
       	       if (dPt > m_dPt_max_factor * resolution * recoJet.Pt()) continue;
       	       dRMin = dR;
               //matched_genJet[idx_jet].SetPxPyPzE(it_genjet->Px(), it_genjet->Py(), it_genjet->Pz(), it_genjet->E());
               matched_genJet = &*it_genjet;
            }

         }
         
         return matched_genJet;
}



double getScaleFactor(double etaJet){

    std::ifstream ifile("/storage1/lhuertas/CMSTOTEM/mc/Workspace/scale_factors.txt",std::ifstream::in);
    float etaMin, etaMax, sf, sf_up, sf_dw, val;
    double scaleFactor = 1.;
    std::string line;  
    std::getline(ifile,line); //skip the first line;

    while(ifile.is_open() && !ifile.eof()) { 
	 ifile >> etaMin >> etaMax >> val >> sf >> sf_dw >> sf_up; 
         //cout<<etaMin << " " << etaMax << " " << sf << "  " << sf_dw << " " << sf_up << endl; 
         if (etaJet<=etaMax && etaJet >etaMin) scaleFactor = sf;
    } 
    ifile.close(); 

    return scaleFactor;

} 
