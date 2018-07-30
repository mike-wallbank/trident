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

void trident_mumu_plot() {

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
  TFile* fileIn = new TFile("TridentMuMuOut.root", "READ");
  //TFile* fileIn = new TFile("TridentMuMuOut_Current.root", "READ");
  //TFile* fileIn = new TFile("TridentMuMuOut_Small.root", "READ");
  TTree* treeSig = (TTree*)fileIn->Get("TridentMuMuSignal");
  TTree* treeBg = (TTree*)fileIn->Get("TridentMuMuBackground");

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
  // cuts.push_back(Cut("TwoMuCandidates", "NumTracks==2"));
  // cuts.push_back(Cut("MinAngle", "TMath::Abs(Angle)>0.5"));
  // cuts.push_back(Cut("MaxAngle", "TMath::Abs(Angle)<0.096*(180/TMath::Pi())"));//&&TMath::Abs(Angle)>0.5"));
  // cuts.push_back(Cut("ShortMuonLength", "ShortMuonLength>3000")); // mm
  // cuts.push_back(Cut("DiffLengths", "LongMuonLength-ShortMuonLength<150")); // mm

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

//   double vAcSig, vAcBg, eEnSig, eEnBg;
//   treeSig->SetBranchAddress("VertexActivity", &vAcSig);
//   treeSig->SetBranchAddress("EventEnergy", &eEnSig);
//   treeBg->SetBranchAddress("VertexActivity", &vAcBg);
//   treeBg->SetBranchAddress("EventEnergy", &eEnBg);
//   TH1F* hVertexActivityFracSig = new TH1F("VertexActivityFracSig",";Fraction of Energy Deposited within 5 cm of Vertex;",100,0,1);
//   TH1F* hVertexActivityFracBg = new TH1F("VertexActivityFracBg",";Fraction of Energy Deposited within 5 cm of Vertex;",100,0,1);
//   for (size_t sigIt = 0; sigIt < treeSig->GetEntriesFast(); ++sigIt) {
//     treeSig->GetEntry(sigIt);
//     if (eEnSig != 0)
// 	hVertexActivityFracSig->Fill(vAcSig/eEnSig);
//   }
//   for (size_t bgIt = 0; bgIt < treeBg->GetEntriesFast(); ++bgIt) {
//     treeBg->GetEntry(bgIt);
//     if (eEnBg != 0)
// 	hVertexActivityFracBg->Fill(vAcBg/eEnBg);
//   }
//   canv->cd();
//   canv->Clear();
//   leg->Clear();
//   hVertexActivityFracSig->SetLineColor(kSignalColor);
//   hVertexActivityFracSig->SetFillColor(kSignalColor);
//   hVertexActivityFracSig->SetFillStyle(kSignalFill);
//   hVertexActivityFracSig->Scale(1./hVertexActivityFracSig->Integral());
//   leg->AddEntry(hVertexActivityFracSig, "Signal", "f");
//   hVertexActivityFracBg->SetLineColor(kBackgroundColor);
//   hVertexActivityFracBg->SetFillColor(kBackgroundColor);
//   hVertexActivityFracBg->SetFillStyle(kBackgroundFill);
//   hVertexActivityFracBg->Scale(1./hVertexActivityFracBg->Integral());
//   leg->AddEntry(hVertexActivityFracBg, "Background", "f");
//   if (hVertexActivityFracSig->GetMaximum() > hVertexActivityFracBg->GetMaximum()) {
//     hVertexActivityFracSig->Draw("hist");
//     hVertexActivityFracBg->Draw("hist same");
//   } else {
//     hVertexActivityFracBg->Draw("hist");
//     hVertexActivityFracSig->Draw("hist same");
//   }
//   leg->Draw();
//   fileOut->cd();
//   canv->Write("VertexActivityFrac");
//   hVertexActivityFracSig->Write();
//   hVertexActivityFracBg->Write();    

// }

  // Look at event numbers and distributions for each of the cuts
  for (std::vector<Cut>::const_iterator cutIt = cuts.begin(); cutIt != cuts.end(); ++cutIt) {

    TEventList signalEvents("SignalEvents"), backgroundEvents("BackgroundEvents");
    treeSig->Draw(">>SignalEvents", cutIt->cut);
    // treeBg->Draw(">>BackgroundEvents", cutIt->cut);
    std::cout << setw(20) << cutIt->name
    	      << setw(15) << signalEvents.GetN()
      //<< setw(15) << backgroundEvents.GetN()
    	      << setw(25) << signalEvents.GetN() * expectedSig/simSig
      //      << setw(25) << backgroundEvents.GetN() * expectedPOT/bgPOTEquiv
    	      << std::endl;
    // std::ofstream outFile(Form("%sCutEvents.txt", cutIt->name.Data()));
    // std::cout << "Cut " << cutIt->name << std::endl;
    // for (unsigned int event = 0; event < signalEvents.GetN(); ++event) {
    //   outFile << signalEvents.GetEntry(event) << std::endl;
    //   std::cout << "  Event " << signalEvents.GetEntry(event) << std::endl;
    // }

//     // Track length long
//     TH1F* hTrackLengthLongSig = new TH1F("TrackLengthLongSig"+cutIt->name,";Track length (long) (cm);",100,0,600);
//     treeSig->Draw("LongMuonLength>>TrackLengthLongSig"+cutIt->name, cutIt->cut);
//     TH1F* hTrackLengthLongBg = new TH1F("TrackLengthLongBg"+cutIt->name,";Track length (long) (cm);",100,0,600);
//     treeBg->Draw("LongMuonLength>>TrackLengthLongBg"+cutIt->name, cutIt->cut);
//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hTrackLengthLongSig->SetLineColor(kSignalColor);
//     hTrackLengthLongSig->SetFillColor(kSignalColor);
//     hTrackLengthLongSig->SetFillStyle(kSignalFill);
//     leg->AddEntry(hTrackLengthLongSig, "Signal", "f");
//     hTrackLengthLongBg->SetLineColor(kBackgroundColor);
//     hTrackLengthLongBg->SetFillColor(kBackgroundColor);
//     hTrackLengthLongBg->SetFillStyle(kBackgroundFill);
//     leg->AddEntry(hTrackLengthLongBg, "Background", "f");
//     if (hTrackLengthLongSig->GetMaximum() > hTrackLengthLongBg->GetMaximum()) {
//       hTrackLengthLongSig->Draw("hist");
//       hTrackLengthLongBg->Draw("hist same");
//     } else {
//       hTrackLengthLongBg->Draw("hist");
//       hTrackLengthLongSig->Draw("hist same");
//     }
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("TrackLengthLong"+cutIt->name);
//     hTrackLengthLongSig->Write();
//     hTrackLengthLongBg->Write();    

//     // Track length short
//     TH1F* hTrackLengthShortSig = new TH1F("TrackLengthShortSig"+cutIt->name,";Track length (short) (cm);",100,0,600);
//     treeSig->Draw("ShortMuonLength>>TrackLengthShortSig"+cutIt->name, cutIt->cut);
//     TH1F* hTrackLengthShortBg = new TH1F("TrackLengthShortBg"+cutIt->name,";Track length (short) (cm);",100,0,600);
//     treeBg->Draw("ShortMuonLength>>TrackLengthShortBg"+cutIt->name, cutIt->cut);
//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hTrackLengthShortSig->SetLineColor(kSignalColor);
//     hTrackLengthShortSig->SetFillColor(kSignalColor);
//     hTrackLengthShortSig->SetFillStyle(kSignalFill);
//     leg->AddEntry(hTrackLengthShortSig, "Signal", "f");
//     hTrackLengthShortBg->SetLineColor(kBackgroundColor);
//     hTrackLengthShortBg->SetFillColor(kBackgroundColor);
//     hTrackLengthShortBg->SetFillStyle(kBackgroundFill);
//     leg->AddEntry(hTrackLengthShortBg, "Background", "f");
//     if (hTrackLengthShortSig->GetMaximum() > hTrackLengthShortBg->GetMaximum()) {
//       hTrackLengthShortSig->Draw("hist");
//       hTrackLengthShortBg->Draw("hist same");
//     } else {
//       hTrackLengthShortBg->Draw("hist");
//       hTrackLengthShortSig->Draw("hist same");
//     }
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("TrackLengthShort"+cutIt->name);
//     hTrackLengthShortSig->Write();
//     hTrackLengthShortBg->Write();    

//     // Energy long
//     TH1F* hEnergyLongSig = new TH1F("EnergyLongSig"+cutIt->name,";Energy (long) (MeV);",100,0,1200);
//     treeSig->Draw("LongMuonEnergy>>EnergyLongSig"+cutIt->name, cutIt->cut);
//     TH1F* hEnergyLongBg = new TH1F("EnergyLongBg"+cutIt->name,";Energy (long) (MeV);",100,0,1200);
//     treeBg->Draw("LongMuonEnergy>>EnergyLongBg"+cutIt->name, cutIt->cut);
//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hEnergyLongSig->SetLineColor(kSignalColor);
//     hEnergyLongSig->SetFillColor(kSignalColor);
//     hEnergyLongSig->SetFillStyle(kSignalFill);
//     leg->AddEntry(hEnergyLongSig, "Signal", "f");
//     hEnergyLongBg->SetLineColor(kBackgroundColor);
//     hEnergyLongBg->SetFillColor(kBackgroundColor);
//     hEnergyLongBg->SetFillStyle(kBackgroundFill);
//     leg->AddEntry(hEnergyLongBg, "Background", "f");
//     if (hEnergyLongSig->GetMaximum() > hEnergyLongBg->GetMaximum()) {
//       hEnergyLongSig->Draw("hist");
//       hEnergyLongBg->Draw("hist same");
//     } else {
//       hEnergyLongBg->Draw("hist");
//       hEnergyLongSig->Draw("hist same");
//     }
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("EnergyLong"+cutIt->name);
//     hEnergyLongSig->Write();
//     hEnergyLongBg->Write();    

//     // Energy short
//     TH1F* hEnergyShortSig = new TH1F("EnergyShortSig"+cutIt->name,";Energy (short) (MeV);",100,0,1200);
//     treeSig->Draw("ShortMuonEnergy>>EnergyShortSig"+cutIt->name, cutIt->cut);
//     TH1F* hEnergyShortBg = new TH1F("EnergyShortBg"+cutIt->name,";Energy (short) (MeV);",100,0,1200);
//     treeBg->Draw("ShortMuonEnergy>>EnergyShortBg"+cutIt->name, cutIt->cut);
//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hEnergyShortSig->SetLineColor(kSignalColor);
//     hEnergyShortSig->SetFillColor(kSignalColor);
//     hEnergyShortSig->SetFillStyle(kSignalFill);
//     leg->AddEntry(hEnergyShortSig, "Signal", "f");
//     hEnergyShortBg->SetLineColor(kBackgroundColor);
//     hEnergyShortBg->SetFillColor(kBackgroundColor);
//     hEnergyShortBg->SetFillStyle(kBackgroundFill);
//     leg->AddEntry(hEnergyShortBg, "Background", "f");
//     if (hEnergyShortSig->GetMaximum() > hEnergyShortBg->GetMaximum()) {
//       hEnergyShortSig->Draw("hist");
//       hEnergyShortBg->Draw("hist same");
//     } else {
//       hEnergyShortBg->Draw("hist");
//       hEnergyShortSig->Draw("hist same");
//     }
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("EnergyShort"+cutIt->name);
//     hEnergyShortSig->Write();
//     hEnergyShortBg->Write();    

//     // Angle
//     TH1F* hAngleSig = new TH1F("AngleSig"+cutIt->name,";Angle (deg);",100,0,180);
//     treeSig->Draw("Angle>>AngleSig"+cutIt->name, cutIt->cut);
//     TH1F* hAngleBg = new TH1F("AngleBg"+cutIt->name,";Angle (deg);",100,0,180);
//     treeBg->Draw("Angle>>AngleBg"+cutIt->name, cutIt->cut);
//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hAngleSig->SetLineColor(kSignalColor);
//     hAngleSig->SetFillColor(kSignalColor);
//     hAngleSig->SetFillStyle(kSignalFill);
//     leg->AddEntry(hAngleSig, "Signal", "f");
//     hAngleBg->SetLineColor(kBackgroundColor);
//     hAngleBg->SetFillColor(kBackgroundColor);
//     hAngleBg->SetFillStyle(kBackgroundFill);
//     leg->AddEntry(hAngleBg, "Background", "f");
//     if (hAngleSig->GetMaximum() > hAngleBg->GetMaximum()) {
//       hAngleSig->Draw("hist");
//       hAngleBg->Draw("hist same");
//     } else {
//       hAngleBg->Draw("hist");
//       hAngleSig->Draw("hist same");
//     }
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("Angle"+cutIt->name);
//     hAngleSig->Write();
//     hAngleBg->Write();    

//     // Front separation
//     TH1F* hFrontSeparationSig = new TH1F("FrontSeparationSig"+cutIt->name,";FrontSeparation (cm);",100,0,500);
//     treeSig->Draw("FrontSeparation>>FrontSeparationSig"+cutIt->name, cutIt->cut);
//     TH1F* hFrontSeparationBg = new TH1F("FrontSeparationBg"+cutIt->name,";FrontSeparation (cm);",100,0,500);
//     treeBg->Draw("FrontSeparation>>FrontSeparationBg"+cutIt->name, cutIt->cut);
//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hFrontSeparationSig->SetLineColor(kSignalColor);
//     hFrontSeparationSig->SetFillColor(kSignalColor);
//     hFrontSeparationSig->SetFillStyle(kSignalFill);
//     leg->AddEntry(hFrontSeparationSig, "Signal", "f");
//     hFrontSeparationBg->SetLineColor(kBackgroundColor);
//     hFrontSeparationBg->SetFillColor(kBackgroundColor);
//     hFrontSeparationBg->SetFillStyle(kBackgroundFill);
//     leg->AddEntry(hFrontSeparationBg, "Background", "f");
//     if (hFrontSeparationSig->GetMaximum() > hFrontSeparationBg->GetMaximum()) {
//       hFrontSeparationSig->Draw("hist");
//       hFrontSeparationBg->Draw("hist same");
//     } else {
//       hFrontSeparationBg->Draw("hist");
//       hFrontSeparationSig->Draw("hist same");
//     }
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("FrontSeparation"+cutIt->name);
//     hFrontSeparationSig->Write();
//     hFrontSeparationBg->Write();    

//     // Back separation
//     TH1F* hBackSeparationSig = new TH1F("BackSeparationSig"+cutIt->name,";BackSeparation (cm);",100,0,500);
//     treeSig->Draw("BackSeparation>>BackSeparationSig"+cutIt->name, cutIt->cut);
//     TH1F* hBackSeparationBg = new TH1F("BackSeparationBg"+cutIt->name,";BackSeparation (cm);",100,0,500);
//     treeBg->Draw("BackSeparation>>BackSeparationBg"+cutIt->name, cutIt->cut);
//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hBackSeparationSig->SetLineColor(kSignalColor);
//     hBackSeparationSig->SetFillColor(kSignalColor);
//     hBackSeparationSig->SetFillStyle(kSignalFill);
//     leg->AddEntry(hBackSeparationSig, "Signal", "f");
//     hBackSeparationBg->SetLineColor(kBackgroundColor);
//     hBackSeparationBg->SetFillColor(kBackgroundColor);
//     hBackSeparationBg->SetFillStyle(kBackgroundFill);
//     leg->AddEntry(hBackSeparationBg, "Background", "f");
//     if (hBackSeparationSig->GetMaximum() > hBackSeparationBg->GetMaximum()) {
//       hBackSeparationSig->Draw("hist");
//       hBackSeparationBg->Draw("hist same");
//     } else {
//       hBackSeparationBg->Draw("hist");
//       hBackSeparationSig->Draw("hist same");
//     }
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("BackSeparation"+cutIt->name);
//     hBackSeparationSig->Write();
//     hBackSeparationBg->Write();

//     // Vertex activity
//     TH1F* hVertexActivitySig = new TH1F("VertexActivitySig"+cutIt->name,";Energy Deposited within 5 cm of Vertex (MeV);",100,0,600);
//     treeSig->Draw("VertexActivity>>VertexActivitySig"+cutIt->name, cutIt->cut);
//     TH1F* hVertexActivityBg = new TH1F("VertexActivityBg"+cutIt->name,";Energy Deposited within 5 cm of Vertex (MeV);",100,0,600);
//     treeBg->Draw("VertexActivity>>VertexActivityBg"+cutIt->name, cutIt->cut);
//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hVertexActivitySig->SetLineColor(kSignalColor);
//     hVertexActivitySig->SetFillColor(kSignalColor);
//     hVertexActivitySig->SetFillStyle(kSignalFill);
//     hVertexActivitySig->Scale(1./hVertexActivitySig->Integral());
//     leg->AddEntry(hVertexActivitySig, "Signal", "f");
//     hVertexActivityBg->SetLineColor(kBackgroundColor);
//     hVertexActivityBg->SetFillColor(kBackgroundColor);
//     hVertexActivityBg->SetFillStyle(kBackgroundFill);
//     hVertexActivityBg->Scale(1./hVertexActivityBg->Integral());
//     leg->AddEntry(hVertexActivityBg, "Background", "f");
//     if (hVertexActivitySig->GetMaximum() > hVertexActivityBg->GetMaximum()) {
//       hVertexActivitySig->Draw("hist");
//       hVertexActivityBg->Draw("hist same");
//     } else {
//       hVertexActivityBg->Draw("hist");
//       hVertexActivitySig->Draw("hist same");
//     }
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("VertexActivity"+cutIt->name);
//     hVertexActivitySig->Write();
//     hVertexActivityBg->Write();    

//     // // Vertex activity
//     // TH1F* hVertexActivityFracSig = new TH1F("VertexActivityFracSig"+cutIt->name,";Fraction of Energy Deposited within 5 cm of Vertex;",100,0,1);
//     // treeSig->Draw("VertexActivityFrac/EventEnergy>>VertexActivityFracSig"+cutIt->name, cutIt->cut);
//     // TH1F* hVertexActivityFracBg = new TH1F("VertexActivityFracBg"+cutIt->name,";Fraction of Energy Deposited within 5 cm of Vertex;",100,0,1);
//     // treeBg->Draw("VertexActivityFrac/EventEnergy>>VertexActivityFracBg"+cutIt->name, cutIt->cut);
//     // canv->cd();
//     // canv->Clear();
//     // leg->Clear();
//     // hVertexActivityFracSig->SetLineColor(kSignalColor);
//     // hVertexActivityFracSig->SetFillColor(kSignalColor);
//     // hVertexActivityFracSig->SetFillStyle(kSignalFill);
//     // hVertexActivityFracSig->Scale(1./hVertexActivityFracSig->Integral());
//     // leg->AddEntry(hVertexActivityFracSig, "Signal", "f");
//     // hVertexActivityFracBg->SetLineColor(kBackgroundColor);
//     // hVertexActivityFracBg->SetFillColor(kBackgroundColor);
//     // hVertexActivityFracBg->SetFillStyle(kBackgroundFill);
//     // hVertexActivityFracBg->Scale(1./hVertexActivityFracBg->Integral());
//     // leg->AddEntry(hVertexActivityFracBg, "Background", "f");
//     // if (hVertexActivityFracSig->GetMaximum() > hVertexActivityFracBg->GetMaximum()) {
//     //   hVertexActivityFracSig->Draw("hist");
//     //   hVertexActivityFracBg->Draw("hist same");
//     // } else {
//     //   hVertexActivityFracBg->Draw("hist");
//     //   hVertexActivityFracSig->Draw("hist same");
//     // }
//     // leg->Draw();
//     // fileOut->cd();
//     // canv->Write("VertexActivityFrac"+cutIt->name);
//     // hVertexActivityFracSig->Write();
//     // hVertexActivityFracBg->Write();    

//     // dEdx
//     TH1F* hdEdx = new TH1F("dEdx"+cutIt->name,";dEdx (MeV/cm);",100,0,5);
//     TLegend* dEdxStackLeg = new TLegend(0.68,0.58,0.88,0.88);
//     std::map<int,TH1F*> hdEdxBg;
//     for (std::map<int,std::pair<std::string,int> >::iterator pdgIt = PDGMap.begin(); pdgIt != PDGMap.end(); ++pdgIt) {
//       //      hIncorrectTracks->GetXaxis()->SetBinLabel(pdgIt->second.second+1,pdgIt->second.first.c_str());
//       hdEdxBg[pdgIt->second.second] = new TH1F(Form("dEdxBg%d"+cutIt->name,pdgIt->first),";dEdx (MeV/cm);",100,0,5);
//     }
//     treeSig->Draw("LongMuondEdx>>dEdx"+cutIt->name, cutIt->cut);
//     treeSig->Draw("ShortMuondEdx>>dEdx"+cutIt->name, cutIt->cut);
//     //treeBg->Draw(Form("LongMuondEdx>>dEdxBg%d",LongMuonPDG), cutIt->cut);
//     //treeBg->Draw("Form(\"ShortMuondEdx>>dEdxBg%d"+cutIt->name+",TMath::Abs(ShortMuonPDG))", cutIt->cut);

//     treeBg->Draw("LongMuondEdx>>dEdxBg13"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==13");
//     treeBg->Draw("LongMuondEdx>>dEdxBg211"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==211");
//     treeBg->Draw("LongMuondEdx>>dEdxBg11"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==11");
//     treeBg->Draw("LongMuondEdx>>dEdxBg2212"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==2212");
//     treeBg->Draw("LongMuondEdx>>dEdxBg321"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==321");
//     treeBg->Draw("ShortMuondEdx>>dEdxBg13"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==13");
//     treeBg->Draw("ShortMuondEdx>>dEdxBg211"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==211");
//     treeBg->Draw("ShortMuondEdx>>dEdxBg11"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==11");
//     treeBg->Draw("ShortMuondEdx>>dEdxBg2212"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==2212");
//     treeBg->Draw("ShortMuondEdx>>dEdxBg321"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==321");

//     canv->cd();
//     canv->Clear();
//     leg->Clear();
//     hdEdx->SetLineColor(kSignalColor);
//     hdEdx->SetFillColor(kSignalColor);
//     hdEdx->SetFillStyle(kSignalFill);
//     dEdxStackLeg->AddEntry(hdEdx, "Signal", "l");
//     THStack* hdEdxStack = new THStack("dEdxStack",";dEdx (MeV/cm);");
//     for (std::map<int,std::pair<std::string,int> >::iterator pdgIt = PDGMap.begin(); pdgIt != PDGMap.end(); ++pdgIt) {
//       hdEdxBg[pdgIt->second.second]->SetLineColor(pdgIt->second.second+1);
//       hdEdxBg[pdgIt->second.second]->SetFillColor(pdgIt->second.second+1);
//       hdEdxStack->Add(hdEdxBg[pdgIt->second.second]);
//       dEdxStackLeg->AddEntry(hdEdxBg[pdgIt->second.second],pdgIt->second.first.c_str(),"f");
//     }
//     hdEdx->Draw("hist");
//     hdEdxStack->Draw("hist same");
//     dEdxStackLeg->Draw();
//     fileOut->cd();
//     canv->Write("dEdx"+cutIt->name);


// //     // Correct tracks
// //     TProfile* hCorrectTracks = new TProfile("CorrectTracksSig",";Nu Energy (MeV);",60,0,5000,0,1);
// //     TProfile* hCorrectTracks = new TProfile("CorrectTracksBg",";Nu Energy (MeV);",60,0,5000,0,1);

// //     // Incorrect tracks
// //     TH1F* hIncorrectTracks = new TH1F("IncorrectTracksSig",";Particle;",PDGMap.size(),0,PDGMap.size());
// //     TH1F* hIncorrectTracks = new TH1F("IncorrectTracksBg",";Particle;",PDGMap.size(),0,PDGMap.size());

  } // cut

  std::cout << "-------------------------------------------------------------------"
	    << "-------------------------------------------------------------------" << std::endl;

  // Efficiency/purity

//   // Neutrino energy
//   TH1F* hNuEnergyAll = new TH1F("NuEnergyAll",";Nu Energy (MeV);",100,0,5000);
//   TH1F* hNuEnergyCut = new TH1F("NuEnergyCut",";Nu Energy Cut (MeV);",100,0,5000);

  // Write out and overlay histograms made by trident_mumu.C
  // TH1D* hTrackLengthSig	        = (TH1D*)fileIn->Get("TrackLengthSig");
  // TH1D* hTrackEnergySig		= (TH1D*)fileIn->Get("TrackEnergySig");
  // TH1D* hLongTrackLengthSig	= (TH1D*)fileIn->Get("LongTrackLengthSig");
  // TH1D* hLongTrackEnergySig	= (TH1D*)fileIn->Get("LongTrackEnergySig");
  // TH1D* hShortTrackLengthSig	= (TH1D*)fileIn->Get("ShortTrackLengthSig");
  // TH1D* hShortTrackEnergySig	= (TH1D*)fileIn->Get("ShortTrackEnergySig");
  // TH1D* hTrackLengthBg		= (TH1D*)fileIn->Get("TrackLengthBg");
  // TH1D* hTrackEnergyBg		= (TH1D*)fileIn->Get("TrackEnergyBg");
  // TH1D* hLongTrackLengthBg	= (TH1D*)fileIn->Get("LongTrackLengthBg");
  // TH1D* hLongTrackEnergyBg	= (TH1D*)fileIn->Get("LongTrackEnergyBg");
  // TH1D* hShortTrackLengthBg	= (TH1D*)fileIn->Get("ShortTrackLengthBg");
  // TH1D* hShortTrackEnergyBg	= (TH1D*)fileIn->Get("ShortTrackEnergyBg");
  // TH2D* hVertexActivitySig	= (TH2D*)fileIn->Get("VertexActivitySig");
  // TH2D* hVertexActivityBg       = (TH2D*)fileIn->Get("VertexActivityBg");
  // TGraph* gVertexActivitySig	= (TGraph*)fileIn->Get("gVertexActivitySig");
  // TGraph* gVertexActivityBg	= (TGraph*)fileIn->Get("gVertexActivityBg");
  TH1D* hVertexEnergySig        = (TH1D*)fileIn->Get("VertexEnergySig");
  TH1D* hVertexEnergyBg         = (TH1D*)fileIn->Get("VertexEnergyBg");

  // canv->cd();
  // canv->Clear();
  // leg->Clear();
  // hTrackLengthSig->SetLineColor(1);
  // hTrackLengthBg->SetLineColor(2);
  // hTrackLengthBg->Draw();
  // hTrackLengthSig->Draw("same");
  // leg->AddEntry(hTrackLengthSig, "Signal", "l");
  // leg->AddEntry(hTrackLengthBg, "Background", "l");
  // leg->Draw();
  // fileOut->cd();
  // canv->Write("TrackLength.eps");

  // canv->cd();
  // canv->Clear();
  // leg->Clear();
  // hLongTrackLengthSig->SetLineColor(1);
  // hLongTrackLengthBg->SetLineColor(2);
  // hLongTrackLengthBg->Draw();
  // hLongTrackLengthSig->Draw("same");
  // leg->AddEntry(hLongTrackLengthSig, "Signal", "l");
  // leg->AddEntry(hLongTrackLengthBg, "Background", "l");
  // leg->Draw();
  // fileOut->cd();
  // canv->Write("LongTrackLength.eps");

  // canv->cd();
  // canv->Clear();
  // leg->Clear();
  // hShortTrackLengthSig->SetLineColor(1);
  // hShortTrackLengthBg->SetLineColor(2);
  // hShortTrackLengthBg->Draw();
  // hShortTrackLengthSig->Draw("same");
  // leg->AddEntry(hShortTrackLengthSig, "Signal", "l");
  // leg->AddEntry(hShortTrackLengthBg, "Background", "l");
  // leg->Draw();
  // fileOut->cd();
  // canv->Write("ShortTrackLength.eps");

  // canv->cd();
  // canv->Clear();
  // leg->Clear();
  // hTrackEnergySig->SetLineColor(1);
  // hTrackEnergyBg->SetLineColor(2);
  // hTrackEnergyBg->Draw();
  // hTrackEnergySig->Draw("same");
  // leg->AddEntry(hTrackEnergySig, "Signal", "l");
  // leg->AddEntry(hTrackEnergyBg, "Background", "l");
  // leg->Draw();
  // fileOut->cd();
  // canv->Write("TrackEnergy.eps");

  // canv->cd();
  // canv->Clear();
  // leg->Clear();
  // hLongTrackEnergySig->SetLineColor(1);
  // hLongTrackEnergyBg->SetLineColor(2);
  // hLongTrackEnergyBg->Draw();
  // hLongTrackEnergySig->Draw("same");
  // leg->AddEntry(hLongTrackEnergySig, "Signal", "l");
  // leg->AddEntry(hLongTrackEnergyBg, "Background", "l");
  // leg->Draw();
  // fileOut->cd();
  // canv->Write("LongTrackEnergy.eps");

  // canv->cd();
  // canv->Clear();
  // leg->Clear();
  // hShortTrackEnergySig->SetLineColor(1);
  // hShortTrackEnergyBg->SetLineColor(2);
  // hShortTrackEnergyBg->Draw();
  // hShortTrackEnergySig->Draw("same");
  // leg->AddEntry(hShortTrackEnergySig, "Signal", "l");
  // leg->AddEntry(hShortTrackEnergyBg, "Background", "l");
  // leg->Draw();
  // fileOut->cd();
  // canv->Write("ShortTrackEnergy.eps");

  // canv->cd();
  // canv->Clear();
  // leg->Clear();
  // TMultiGraph* mg = new TMultiGraph();
  // gVertexActivitySig->SetMarkerColor(1);
  // gVertexActivitySig->SetMarkerStyle(8);
  // gVertexActivitySig->SetMarkerSize(0.5);
  // gVertexActivityBg->SetMarkerColor(2);
  // gVertexActivityBg->SetMarkerStyle(8);
  // gVertexActivityBg->SetMarkerSize(0.5);
  // mg->Add(gVertexActivitySig);
  // mg->Add(gVertexActivityBg);
  // mg->Draw("AP");
  // mg->SetNameTitle("VertexActivity",";Distance from vertex (cm);Energy deposited (MeV);");
  // leg->AddEntry(gVertexActivitySig, "Signal", "p");
  // leg->AddEntry(gVertexActivityBg, "Background", "p");
  // leg->Draw();
  // fileOut->cd();
  // canv->Write("VertexActivity.eps");

  canv->cd();
  canv->Clear();
  leg->Clear();
  hVertexEnergySig->SetLineColor(1);
  hVertexEnergyBg->SetLineColor(2);
  hVertexEnergySig->Scale(1./hVertexEnergySig->Integral());
  hVertexEnergyBg->Scale(1./hVertexEnergyBg->Integral());
  hVertexEnergyBg->Draw("hist");
  hVertexEnergySig->Draw("hist same");
  leg->AddEntry(hVertexEnergySig, "Signal", "l");
  leg->AddEntry(hVertexEnergyBg, "Background", "l");
  leg->Draw();
  fileOut->cd();
  canv->Write("VertexEnergy.eps");  

}
