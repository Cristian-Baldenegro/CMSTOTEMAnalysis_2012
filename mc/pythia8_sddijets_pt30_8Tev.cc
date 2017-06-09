

#include "Pythia8/Pythia.h"

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"

#include "Pythia8Plugins/FastJet3.h"
#define M_P 0.938272046

using namespace Pythia8;

int main() {

  // Number of events, generated and listed ones (for jets).
  int nEvent    = 2000000;

  // Set up the ROOT TFile and TTree.
  double Jet1_pz, Jet1_pt, Jet1_energy, Jet1_phi, Jet1_eta;
  double Jet2_pz, Jet2_pt, Jet2_energy, Jet2_phi, Jet2_eta;
  double xP1, xP2, tP1, tP2, xiP1, xiP2;
  double x1, x2, sigma;
  double proton_pz_plus, proton_px_plus, proton_py_plus, proton_energy_plus, proton_energy_corrected_plus, proton_mass_plus;
  double proton_pz_minus, proton_px_minus, proton_py_minus, proton_energy_minus, proton_energy_corrected_minus, proton_mass_minus;

  TFile *file = TFile::Open("pythia82_sddijets_NLO_MPI-unchecked.root","recreate");

  TTree *Tree = new TTree("Tree","ev1 Tree");
  Tree->Branch("Jet1_pz",&Jet1_pz,"Jet1_pz/D");
  Tree->Branch("Jet1_pt",&Jet1_pt,"Jet1_pt/D");
  Tree->Branch("Jet1_phi",&Jet1_phi,"Jet1_phi/D");
  Tree->Branch("Jet1_eta",&Jet1_eta,"Jet1_eta/D");
  Tree->Branch("Jet1_energy",&Jet1_energy,"Jet1_energy/D");
  Tree->Branch("Jet2_pz",&Jet2_pz,"Jet2_pz/D");
  Tree->Branch("Jet2_pt",&Jet2_pt,"Jet2_pt/D");
  Tree->Branch("Jet2_phi",&Jet2_phi,"Jet2_phi/D");
  Tree->Branch("Jet2_eta",&Jet2_eta,"Jet2_eta/D");
  Tree->Branch("Jet2_energy",&Jet2_energy,"Jet2_energy/D");
  Tree->Branch("xP1",&xP1,"xP1/D");
  Tree->Branch("xP2",&xP2,"xP2/D");
  Tree->Branch("tP1",&tP1,"tP1/D");
  Tree->Branch("tP2",&tP2,"tP2/D");
  Tree->Branch("x1",&x1,"x1/D");
  Tree->Branch("x2",&x2,"x2/D");
  Tree->Branch("sigma",&sigma,"sigma/D");
  Tree->Branch("proton_pz_plus",&proton_pz_plus,"proton_pz_plus/D");
  Tree->Branch("proton_px_plus",&proton_px_plus,"proton_px_plus/D");
  Tree->Branch("proton_py_plus",&proton_py_plus,"proton_py_plus/D");
  Tree->Branch("proton_mass_plus",&proton_mass_plus,"proton_mass_plus/D");
  Tree->Branch("proton_energy_plus",&proton_energy_plus,"proton_energy_plus/D");
  Tree->Branch("proton_energy_corrected_plus",&proton_energy_corrected_plus,"proton_energy_corrected_plus/D");
  Tree->Branch("proton_pz_minus",&proton_pz_minus,"proton_pz_minus/D");
  Tree->Branch("proton_px_minus",&proton_px_minus,"proton_px_minus/D");
  Tree->Branch("proton_py_minus",&proton_py_minus,"proton_py_minus/D");
  Tree->Branch("proton_mass_minus",&proton_mass_minus,"proton_mass_minus/D");
  Tree->Branch("proton_energy_minus",&proton_energy_minus,"proton_energy_minus/D");
  Tree->Branch("proton_energy_corrected_minus",&proton_energy_corrected_minus,"proton_energy_corrected_minus/D");

  // Select common parameters for SlowJet and FastJet analyses.
  int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT.
  double R       = 0.5;    // Jet size.
  double pTMin   = 0.5;    // Min jet pT.
  //double etaMax  = 5.0;    // Pseudorapidity range of detector.
  //int    select  = 2;      // Which particles are included?
  //int    massSet = 2;      // Which mass are they assumed to have?

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;
  Info&  info  = pythia.info;

  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 30.");
  pythia.readString("Tune:pp = 18");

  // Setup of diffractive framework.
  pythia.readString("Diffraction:doHard = on");
  //pythia.readString("Diffraction:sampleType = 4"); //MPI-checked
  pythia.readString("Diffraction:sampleType = 3"); //MPI-unchecked
  pythia.readString("Diffraction:PomFlux = 7");
  pythia.readString("PDF:PomSet = 4"); // H1 2006 Fit B NLO
  //pythia.readString("PDF:PomSet = 6"); // H1 2006 Fit B LO

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // LHC initialization.
  pythia.readString("Beams:eCM = 8000.");
  pythia.init();


  // Set up FastJet jet finder.
  //   one can use either explicitly use antikt, cambridge, etc., or
  //   just use genkt_algorithm with specification of power
  //fastjet::JetAlgorithm algorithm;
  //if (power == -1)      algorithm = fastjet::antikt_algorithm;
  //if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
  //if (power ==  1)      algorithm = fastjet::kt_algorithm;
  //fastjet::JetDefinition jetDef(algorithm, R);
  // there's no need for a pointer to the jetDef (it's a fairly small object)
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);
  std::vector <fastjet::PseudoJet> fjInputs;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    clock_t befGen = clock();
    if (!pythia.next()) continue;
    clock_t aftGen = clock();


    // Begin FastJet analysis: extract particles from event record.
    clock_t befFast = clock();
    fjInputs.resize(0);

    proton_pz_plus = -999;
    proton_pz_minus = 999;
    proton_px_plus = 0;
    proton_px_minus = 0;
    proton_py_plus = 0;
    proton_py_minus = 0;
    proton_energy_plus = 0;
    proton_energy_minus = 0;
    proton_mass_plus = 0;
    proton_mass_minus = 0;
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

      // Create a PseudoJet from the complete Pythia particle.
      fastjet::PseudoJet particleTemp = event[i];

      // Store acceptable particles as input to Fastjet.
      // Conversion to PseudoJet is performed automatically
      // with the help of the code in FastJet3.h.
      fjInputs.push_back( particleTemp);

      // proton information
      if (event[i].isFinal()){
         if (event[i].id() == 2212){
            double pz_min = 0.7*4000;
	    if (fabs(event[i].pz()) < pz_min) continue;
	    if (event[i].pz() > proton_pz_plus){ 
	       proton_pz_plus = event[i].pz();
	       proton_px_plus = event[i].px();
	       proton_py_plus = event[i].py();
	       proton_mass_plus = event[i].m();
	       proton_energy_plus = event[i].e();
               proton_energy_corrected_plus = sqrt(proton_px_plus*proton_px_plus + proton_py_plus*proton_py_plus + proton_pz_plus*proton_pz_plus + M_P*M_P);
	    }
 
            if (event[i].pz() < proton_pz_minus){
               proton_pz_minus = event[i].pz();
               proton_px_minus = event[i].px();
               proton_py_minus = event[i].py();
               proton_mass_minus = event[i].m();
               proton_energy_minus = event[i].e();
               proton_energy_corrected_minus = sqrt(proton_px_minus*proton_px_minus + proton_py_minus*proton_py_minus + proton_pz_minus*proton_pz_minus + M_P*M_P);
            }
         }
      }
    }

    // Run Fastjet algorithm and sort jets in pT order.
    vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
    fastjet::ClusterSequence clustSeq(fjInputs, jetDef);
    inclusiveJets = clustSeq.inclusive_jets(pTMin);
    sortedJets    = sorted_by_pt(inclusiveJets);
    clock_t aftFast = clock();

    // Information of the two leading jets.
    Jet1_pz = sortedJets[0].pz();
    Jet1_pt = sortedJets[0].pt();
    Jet1_eta = sortedJets[0].eta();
    Jet1_phi = sortedJets[0].phi();
    Jet1_energy = sortedJets[0].E();
    Jet2_pz = sortedJets[1].pz();
    Jet2_pt = sortedJets[1].pt();
    Jet2_eta = sortedJets[1].eta();
    Jet2_phi = sortedJets[1].phi();
    Jet2_energy = sortedJets[1].E();
 
    /*if (iEvent < nListJets) {
      cout << "\n --------  FastJet jets, p = " << setw(2) << power
           << "  --------------------------------------------------\n\n "
           << "  i         pT        y      phi  mult chgmult photons"
           << "      hardest  pT in neutral " << endl
           << "                                                       "
           << "  constituent        hadrons " << endl;
      for (int i = 0; i < int(sortedJets.size()); ++i) {
        vector<fastjet::PseudoJet> constituents
          = sortedJets[i].constituents();
        fastjet::PseudoJet hardest
          = fastjet::SelectorNHardest(1)(constituents)[0];
        vector<fastjet::PseudoJet> neutral_hadrons
          = ( fastjet::SelectorIsHadron()
           && fastjet::SelectorIsNeutral())(constituents);
        double neutral_hadrons_pt = join(neutral_hadrons).perp();
        cout << setw(4) << i << fixed << setprecision(3) << setw(11)
             << sortedJets[i].perp() << setw(9)  << sortedJets[i].rap()
             << setw(9) << sortedJets[i].phi_std()
             << setw(6) << constituents.size()
             << setw(8) << fastjet::SelectorIsCharged().count(constituents)
             << setw(8) << fastjet::SelectorId(22).count(constituents)
             << setw(13) << hardest.user_info<Particle>().name()
             << "     " << setw(10) << neutral_hadrons_pt << endl;
      }
      cout << "\n --------  End FastJet Listing  ------------------"
           << "---------------------------------" << endl;
    }*/

    //partons in the hard process
    x1 = info.x1pdf();	
    x2 = info.x2pdf();	

    // x_Pomeron and t.
    if ( info.isHardDiffractiveA() ) {
        xP1 = info.xPomeronB();
        tP1 = info.tPomeronB();
    }
    if ( info.isHardDiffractiveB() ) {
        xP2 = info.xPomeronA();
        tP2 = info.tPomeronA();
    }

    //cross section info
    sigma = info.sigmaGen(0);	

    Tree->Fill();

  // End of event loop.
  }
  pythia.stat();

  Tree->Write();
  // Done.
  return 0;
}
