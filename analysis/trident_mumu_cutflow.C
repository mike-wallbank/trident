#include <iostream>
#include <iomanip>

#include "TStyle.h"

// Style
const int kSignalColor = 1;
const int kBackgroundColor = 2;
const int kSignalFill = 3003;
const int kBackgroundFill = 3003;

struct Cut {
  TString name;
  TString cut;
  Cut(const char* n, const char* c) {
    name = TString(n);
    cut = TString(c);
  }
};

std::map<int,std::pair<std::string,int> > PDGMap = {{13,{"#mu",0}},{211,{"#pi",1}},{11,{"e",2}},{2212,{"p",3}},{321,{"K",4}}};

void trident_mumu_cutflow() {

  // Plot style
  // gStyle->SetPadTickX(1);
  // gStyle->SetPadTickY(1);
  // gStyle->SetOptStat("e");
  // gStyle->SetStatX(0.9);
  // gStyle->SetStatY(0.937);
  // gStyle->SetStatW(0.3);
  // gStyle->SetStatH(0.15);
  // gStyle->SetHistLineColor(1);
  // TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);

  // Open input file
  TFile* fileSig = new TFile("/pnfs/dune/persistent/users/wallbank/trident/data/tree_10cmVertex/TridentMuMuOut_SIG.root", "READ");
  TTree* treeSig = (TTree*)fileSig->Get("TridentMuMuSignal");
  TFile* fileBg = new TFile("/pnfs/dune/persistent/users/wallbank/trident/data/tree_10cmVertex/TridentMuMuOut_BG.root", "READ");
  TTree* treeBg = (TTree*)fileBg->Get("TridentMuMuBackground");

  // Output file
  TFile* fileOut = new TFile("TridentMuMuPlots.root", "RECREATE");
  TCanvas* canv = new TCanvas("canv","",800,600);
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);

  // Define the cuts
  std::vector<Cut> cuts;
  cuts.push_back(Cut("NoCut", ""));
  cuts.push_back(Cut("TwoMuCandidates", "NumTracks==2&&Angle_s>0.5"));
  cuts.push_back(Cut("Angle", "Angle_s<5.5"));
  cuts.push_back(Cut("ShortMuonLength", "ShortMuonLength>3000"));
  cuts.push_back(Cut("DiffLengths", "LongMuonLength-ShortMuonLength<150"));
  cuts.push_back(Cut("VertexEnergy", "VertexActivity<70"));

  // Define scalings
  float expectedPOT = 1.46E21;
  float bgPOTEquiv = 1e20;
  float expectedSig = 150;
  float simSig = 1e4;

  // Add all cuts
  std::stringstream allcuts;
  for (unsigned int cutIt = 0; cutIt < cuts.size(); ++cutIt) {
    if (cuts[cutIt].name.Contains("NoCut"))
      continue;
    allcuts << cuts[cutIt].cut;
    if (cutIt < cuts.size()-1)
      allcuts << "&&";
  }
  cuts.push_back(Cut("AllCuts", allcuts.str().c_str()));

  std::cout << "-------------------------------------------------------------------"
	    << "-------------------------------------------------------------------" << std::endl
	    << setw(50) << "" << setw(31) << "1 year" << std::endl
	    << setw(50) << "" << setw(28) << "      ---------------------------------------------" << std::endl
	    << setw(20) << "Cut" << setw(15) << "Signal" << setw(15) << "Background" 
	    << setw(25) << "Signal" << setw(25) << "Background" << std::endl
	    << "-------------------------------------------------------------------"
	    << "-------------------------------------------------------------------" << std::endl;

  // Look at event numbers and distributions for each of the cuts
  for (std::vector<Cut>::const_iterator cutIt = cuts.begin(); cutIt != cuts.end(); ++cutIt) {

    TEventList signalEvents("SignalEvents"), backgroundEvents("BackgroundEvents");
    treeSig->Draw(">>SignalEvents", cutIt->cut);
    treeBg->Draw(">>BackgroundEvents", cutIt->cut);
    std::cout << setw(20) << cutIt->name
    	      << setw(15) << signalEvents.GetN()
	      << setw(15) << backgroundEvents.GetN()
    	      << setw(25) << signalEvents.GetN() * expectedSig/simSig
	      << setw(25) << backgroundEvents.GetN() * expectedPOT/bgPOTEquiv
    	      << std::endl;

  } // cut

  std::cout << "-------------------------------------------------------------------"
	    << "-------------------------------------------------------------------" << std::endl;

}
