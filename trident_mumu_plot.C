#include <iostream>

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

  // Open input file
  TFile* fileIn = new TFile("TridentMuMuOut.root", "READ");
  TTree* treeSig = (TTree*)fileIn->Get("TridentMuMuSignal");
  TTree* treeBg = (TTree*)fileIn->Get("TridentMuMuBackground");

  // Output file
  TFile* fileOut = new TFile("TridentMuMuPlots.root", "RECREATE");
  TCanvas* canv = new TCanvas("canv","",800,600);
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);

  // Define the cuts
  std::vector<Cut> cuts;
  cuts.push_back(Cut("NoCut", ""));
//   cuts.push_back(Cut("Angle", "TMath::Abs(TVector3(LongMuonMomentumX,LongMuonMomentumY,LongMuonMomentumZ).Angle(
//                                           TVector3(ShortMuonMomentumX,ShortMuonMomentumY,ShortMuonMomentumZ))*180/TMath::Pi()<3.5"));
  cuts.push_back(Cut("ShortMuonLength", "ShortMuonLength>180"));
  cuts.push_back(Cut("DiffLengths", "LongMuonLength-ShortMuonLength<30"));

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

  // Look at variables after each of the cuts
  for (std::vector<Cut>::const_iterator cutIt = cuts.begin(); cutIt != cuts.end(); ++cutIt) {

    // Track length long
    TH1F* hTrackLengthLongSig = new TH1F("TrackLengthLongSig"+cutIt->name,";Track length (long) (cm);",100,0,600);
    treeSig->Draw("LongMuonLength>>TrackLengthLongSig"+cutIt->name, cutIt->cut);
    TH1F* hTrackLengthLongBg = new TH1F("TrackLengthLongBg"+cutIt->name,";Track length (long) (cm);",100,0,600);
    treeBg->Draw("LongMuonLength>>TrackLengthLongBg"+cutIt->name, cutIt->cut);
    canv->cd();
    canv->Clear();
    leg->Clear();
    hTrackLengthLongSig->SetLineColor(kSignalColor);
    hTrackLengthLongSig->SetFillColor(kSignalColor);
    hTrackLengthLongSig->SetFillStyle(kSignalFill);
    leg->AddEntry(hTrackLengthLongSig, "Signal", "f");
    hTrackLengthLongBg->SetLineColor(kBackgroundColor);
    hTrackLengthLongBg->SetFillColor(kBackgroundColor);
    hTrackLengthLongBg->SetFillStyle(kBackgroundFill);
    leg->AddEntry(hTrackLengthLongBg, "Background", "f");
    hTrackLengthLongSig->Draw("hist");
    hTrackLengthLongBg->Draw("hist same");
    leg->Draw();
    fileOut->cd();
    canv->Write("TrackLengthLong"+cutIt->name);
    hTrackLengthLongSig->Write();
    hTrackLengthLongBg->Write();    

    // Track length short
    TH1F* hTrackLengthShortSig = new TH1F("TrackLengthShortSig"+cutIt->name,";Track length (short) (cm);",100,0,600);
    treeSig->Draw("ShortMuonLength>>TrackLengthShortSig"+cutIt->name, cutIt->cut);
    TH1F* hTrackLengthShortBg = new TH1F("TrackLengthShortBg"+cutIt->name,";Track length (short) (cm);",100,0,600);
    treeBg->Draw("ShortMuonLength>>TrackLengthShortBg"+cutIt->name, cutIt->cut);
    canv->cd();
    canv->Clear();
    leg->Clear();
    hTrackLengthShortSig->SetLineColor(kSignalColor);
    hTrackLengthShortSig->SetFillColor(kSignalColor);
    hTrackLengthShortSig->SetFillStyle(kSignalFill);
    leg->AddEntry(hTrackLengthShortSig, "Signal", "f");
    hTrackLengthShortBg->SetLineColor(kBackgroundColor);
    hTrackLengthShortBg->SetFillColor(kBackgroundColor);
    hTrackLengthShortBg->SetFillStyle(kBackgroundFill);
    leg->AddEntry(hTrackLengthShortBg, "Background", "f");
    hTrackLengthShortSig->Draw("hist");
    hTrackLengthShortBg->Draw("hist same");
    leg->Draw();
    fileOut->cd();
    canv->Write("TrackLengthShort"+cutIt->name);
    hTrackLengthShortSig->Write();
    hTrackLengthShortBg->Write();    

    // Energy long
    TH1F* hEnergyLongSig = new TH1F("EnergyLongSig"+cutIt->name,";Energy (long) (MeV);",100,0,1200);
    treeSig->Draw("LongMuonEnergy>>EnergyLongSig"+cutIt->name, cutIt->cut);
    TH1F* hEnergyLongBg = new TH1F("EnergyLongBg"+cutIt->name,";Energy (long) (MeV);",100,0,1200);
    treeBg->Draw("LongMuonEnergy>>EnergyLongBg"+cutIt->name, cutIt->cut);
    canv->cd();
    canv->Clear();
    leg->Clear();
    hEnergyLongSig->SetLineColor(kSignalColor);
    hEnergyLongSig->SetFillColor(kSignalColor);
    hEnergyLongSig->SetFillStyle(kSignalFill);
    leg->AddEntry(hEnergyLongSig, "Signal", "f");
    hEnergyLongBg->SetLineColor(kBackgroundColor);
    hEnergyLongBg->SetFillColor(kBackgroundColor);
    hEnergyLongBg->SetFillStyle(kBackgroundFill);
    leg->AddEntry(hEnergyLongBg, "Background", "f");
    hEnergyLongSig->Draw("hist");
    hEnergyLongBg->Draw("hist same");
    leg->Draw();
    fileOut->cd();
    canv->Write("EnergyLong"+cutIt->name);
    hEnergyLongSig->Write();
    hEnergyLongBg->Write();    

    // Energy short
    TH1F* hEnergyShortSig = new TH1F("EnergyShortSig"+cutIt->name,";Energy (short) (MeV);",100,0,1200);
    treeSig->Draw("ShortMuonEnergy>>EnergyShortSig"+cutIt->name, cutIt->cut);
    TH1F* hEnergyShortBg = new TH1F("EnergyShortBg"+cutIt->name,";Energy (short) (MeV);",100,0,1200);
    treeBg->Draw("ShortMuonEnergy>>EnergyShortBg"+cutIt->name, cutIt->cut);
    canv->cd();
    canv->Clear();
    leg->Clear();
    hEnergyShortSig->SetLineColor(kSignalColor);
    hEnergyShortSig->SetFillColor(kSignalColor);
    hEnergyShortSig->SetFillStyle(kSignalFill);
    leg->AddEntry(hEnergyShortSig, "Signal", "f");
    hEnergyShortBg->SetLineColor(kBackgroundColor);
    hEnergyShortBg->SetFillColor(kBackgroundColor);
    hEnergyShortBg->SetFillStyle(kBackgroundFill);
    leg->AddEntry(hEnergyShortBg, "Background", "f");
    hEnergyShortSig->Draw("hist");
    hEnergyShortBg->Draw("hist same");
    leg->Draw();
    fileOut->cd();
    canv->Write("EnergyShort"+cutIt->name);
    hEnergyShortSig->Write();
    hEnergyShortBg->Write();    

//     // Angle
//     TH1F* hAngleSig = new TH1F("AngleSig"+cutIt->name,";Angle (deg);",100,0,180);
//     treeSig->Draw("TVector3(LongMuonMomentumX,LongMuonMomentumY,LongMuonMomentumZ).Angle(
//                    TVector3(ShortMuonMomentumX,ShortMuonMomentumY,ShortMuonMomentumZ))*180/TMath::Pi()>>AngleSig"+cutIt->name, cutIt->cut);
//     TH1F* hAngleBg = new TH1F("AngleBg"+cutIt->name,";Angle (deg);",100,0,180);
//     treeBg->Draw("TVector3(LongMuonMomentumX,LongMuonMomentumY,LongMuonMomentumZ).Angle(
//                   TVector3(ShortMuonMomentumX,ShortMuonMomentumY,ShortMuonMomentumZ))*180/TMath::Pi()>>AngleBg"+cutIt->name, cutIt->cut);
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
//     hAngleSig->Draw("hist");
//     hAngleBg->Draw("hist same");
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("Angle"+cutIt->name);
//     hAngleSig->Write();
//     hAngleBg->Write();    

//     // Front separation
//     TH1F* hFrontSeparationSig = new TH1F("FrontSeparationSig"+cutIt->name,";FrontSeparation (cm);",100,0,500);
//     treeSig->Draw("(TVector3(LongMuonStartX,LongMuonStartY,LongMuonStartZ)-(
//                     TVector3(ShortMuonStartX,ShortMuonStartY,ShortMuonStartZ))).Mag()>>FrontSeparationSig"+cutIt->name, cutIt->cut);
//     TH1F* hFrontSeparationBg = new TH1F("FrontSeparationBg"+cutIt->name,";FrontSeparation (cm);",100,0,500);
//     treeBg->Draw("(TVector3(LongMuonEndX,LongMuonEndY,LongMuonEndZ)-(
//                    TVector3(ShortMuonEndX,ShortMuonEndY,ShortMuonEndZ))).Mag()>>FrontSeparationBg"+cutIt->name, cutIt->cut);
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
//     hFrontSeparationSig->Draw("hist");
//     hFrontSeparationBg->Draw("hist same");
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("FrontSeparation"+cutIt->name);
//     hFrontSeparationSig->Write();
//     hFrontSeparationBg->Write();    

//     // Back separation
//     TH1F* hBackSeparationSig = new TH1F("BackSeparationSig"+cutIt->name,";BackSeparation (cm);",100,0,500);
//     treeSig->Draw("(TVector3(LongMuonStartX,LongMuonStartY,LongMuonStartZ)-(
//                     TVector3(ShortMuonStartX,ShortMuonStartY,ShortMuonStartZ))).Mag()>>BackSeparationSig"+cutIt->name, cutIt->cut);
//     TH1F* hBackSeparationBg = new TH1F("BackSeparationBg"+cutIt->name,";BackSeparation (cm);",100,0,500);
//     treeBg->Draw("(TVector3(LongMuonEndX,LongMuonEndY,LongMuonEndZ)-(
//                    TVector3(ShortMuonEndX,ShortMuonEndY,ShortMuonEndZ))).Mag()>>BackSeparationBg"+cutIt->name, cutIt->cut);
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
//     hBackSeparationSig->Draw("hist");
//     hBackSeparationBg->Draw("hist same");
//     leg->Draw();
//     fileOut->cd();
//     canv->Write("BackSeparation"+cutIt->name);
//     hBackSeparationSig->Write();
//     hBackSeparationBg->Write();    

    // dEdx
    TH1F* hdEdx = new TH1F("dEdx"+cutIt->name,";dEdx (MeV/cm);",100,0,5);
    TLegend* dEdxStackLeg = new TLegend(0.68,0.58,0.88,0.88);
    std::map<int,TH1F*> hdEdxBg;
    for (std::map<int,std::pair<std::string,int> >::iterator pdgIt = PDGMap.begin(); pdgIt != PDGMap.end(); ++pdgIt) {
      //      hIncorrectTracks->GetXaxis()->SetBinLabel(pdgIt->second.second+1,pdgIt->second.first.c_str());
      hdEdxBg[pdgIt->second.second] = new TH1F(Form("dEdxBg%d"+cutIt->name,pdgIt->first),";dEdx (MeV/cm);",100,0,5);
    }
    treeSig->Draw("LongMuondEdx>>dEdx"+cutIt->name, cutIt->cut);
    treeSig->Draw("ShortMuondEdx>>dEdx"+cutIt->name, cutIt->cut);
    treeBg->Draw("LongMuondEdx>>dEdxBg13"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==13");
    treeBg->Draw("LongMuondEdx>>dEdxBg211"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==211");
    treeBg->Draw("LongMuondEdx>>dEdxBg11"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==11");
    treeBg->Draw("LongMuondEdx>>dEdxBg2212"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==2212");
    treeBg->Draw("LongMuondEdx>>dEdxBg321"+cutIt->name, cutIt->cut+"&&TMath::Abs(LongMuonPDG)==321");
    treeBg->Draw("ShortMuondEdx>>dEdxBg13"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==13");
    treeBg->Draw("ShortMuondEdx>>dEdxBg211"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==211");
    treeBg->Draw("ShortMuondEdx>>dEdxBg11"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==11");
    treeBg->Draw("ShortMuondEdx>>dEdxBg2212"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==2212");
    treeBg->Draw("ShortMuondEdx>>dEdxBg321"+cutIt->name, cutIt->cut+"&&TMath::Abs(ShortMuonPDG)==321");
    canv->cd();
    canv->Clear();
    leg->Clear();
    hdEdx->SetLineColor(kSignalColor);
    hdEdx->SetFillColor(kSignalColor);
    hdEdx->SetFillStyle(kSignalFill);
    dEdxStackLeg->AddEntry(hdEdx, "Signal", "l");
    THStack* hdEdxStack = new THStack("dEdxStack",";dEdx (MeV/cm);");
    for (std::map<int,std::pair<std::string,int> >::iterator pdgIt = PDGMap.begin(); pdgIt != PDGMap.end(); ++pdgIt) {
      hdEdxBg[pdgIt->second.second]->SetLineColor(pdgIt->second.second+1);
      hdEdxBg[pdgIt->second.second]->SetFillColor(pdgIt->second.second+1);
      hdEdxStack->Add(hdEdxBg[pdgIt->second.second]);
      dEdxStackLeg->AddEntry(hdEdxBg[pdgIt->second.second],pdgIt->second.first.c_str(),"f");
    }
    hdEdx->Draw("hist");
    hdEdxStack->Draw("hist same");
    dEdxStackLeg->Draw();
    fileOut->cd();
    canv->Write("dEdx"+cutIt->name);


//     // Correct tracks
//     TProfile* hCorrectTracks = new TProfile("CorrectTracksSig",";Nu Energy (MeV);",60,0,5000,0,1);
//     TProfile* hCorrectTracks = new TProfile("CorrectTracksBg",";Nu Energy (MeV);",60,0,5000,0,1);

//     // Incorrect tracks
//     TH1F* hIncorrectTracks = new TH1F("IncorrectTracksSig",";Particle;",PDGMap.size(),0,PDGMap.size());
//     TH1F* hIncorrectTracks = new TH1F("IncorrectTracksBg",";Particle;",PDGMap.size(),0,PDGMap.size());

  } // cut

  // Efficiency/purity

//   // Neutrino energy
//   TH1F* hNuEnergyAll = new TH1F("NuEnergyAll",";Nu Energy (MeV);",100,0,5000);
//   TH1F* hNuEnergyCut = new TH1F("NuEnergyCut",";Nu Energy Cut (MeV);",100,0,5000);

}


// void FillHists(AnalysisHists& hists, TridentMuMu tridentMuMu) {

//   std::vector<Particle> mus = tridentMuMu.GetMus();

//   // In the instance of multiple identified muons, need to pick two
//   // Right now, take the longest track
//   std::map<float,Particle> muLength;
//   for (std::vector<Particle>::iterator muIt = mus.begin(); muIt != mus.end(); ++muIt)
//     muLength[muIt->Length()] = *muIt;
//   std::map<float,Particle>::reverse_iterator muIt = muLength.rbegin();
//   Particle muLong = muIt->second;
//   Particle muShort = std::next(muIt)->second;

//   // Find if the two selected muons are the correct two!
//   bool correctSelection = false;
//   if ((mus.at(0).ID() == tridentMuMu.GetMuPlus().ID() and mus.at(1).ID() == tridentMuMu.GetMuMinus().ID()) or
//       (mus.at(1).ID() == tridentMuMu.GetMuPlus().ID() and mus.at(0).ID() == tridentMuMu.GetMuMinus().ID()))
//     correctSelection = true;

//   hists.hTrackLengthLong->Fill(muLong.Length());
//   hists.hTrackLengthShort->Fill(muShort.Length());
//   hists.hEnergyLong->Fill(muLong.Energy());
//   hists.hEnergyShort->Fill(muShort.Energy());
//   hists.hdEdx->Fill(muLong.dEdx());
//   hists.hdEdx->Fill(muShort.dEdx());
//   hists.hdEdxBG[PDGMap[abs(muLong.PDG())].second]->Fill(muLong.dEdx());
//   hists.hdEdxBG[PDGMap[abs(muShort.PDG())].second]->Fill(muShort.dEdx());
//   hists.hAngle->Fill(muLong.Direction().Angle(muShort.Direction())*180/TMath::Pi());
//   hists.hFrontSeparation->Fill((muLong.Start()-muShort.Start()).Mag());
//   hists.hBackSeparation->Fill((muLong.End()-muShort.End()).Mag());
//   //hists.hSeparationTrack->Fill((muLong.StartTrack()-muShort.StartTrack()).Mag());
//   hists.hCorrectTracks->Fill(tridentMuMu.GetEventEnergy(), (int)correctSelection);
//   for (std::vector<Particle>::iterator muIt = mus.begin(); muIt != mus.end(); ++muIt)
//     if (PDGMap.count(abs(muIt->PDG())))
//       hists.hIncorrectTracks->Fill(PDGMap[abs(muIt->PDG())].first.c_str(),1);

//   return;
    
// }

// void PlotHists(const AnalysisHists& sig_hists, const AnalysisHists& bg_hists) {

//   // Open the out file
//   TFile outFile("mumu_hists.root","RECREATE");
  
//   sig_hists.hNuEnergyAll->Write();
//   sig_hists.hNuEnergyCut->Write();
//   sig_hists.hTrackLengthLong->Write();
//   sig_hists.hTrackLengthShort->Write();
//   sig_hists.hEnergyLong->Write();
//   sig_hists.hEnergyShort->Write();
//   sig_hists.hdEdx->Write();
//   for (std::map<int,std::pair<std::string,int> >::iterator pdgIt = PDGMap.begin(); pdgIt != PDGMap.end(); ++pdgIt)
//     sig_hists.hdEdxBG[pdgIt->second.second]->Write();
//   sig_hists.hAngle->Write();
//   sig_hists.hFrontSeparation->Write();
//   sig_hists.hBackSeparation->Write();
//   //sig_hists.hSeparationTrack->Write();
//   sig_hists.hCorrectTracks->Write();
//   sig_hists.hIncorrectTracks->Write();

//   bg_hists.hNuEnergyAll->Write();
//   bg_hists.hNuEnergyCut->Write();
//   bg_hists.hTrackLengthLong->Write();
//   bg_hists.hTrackLengthShort->Write();
//   bg_hists.hEnergyLong->Write();
//   bg_hists.hEnergyShort->Write();
//   bg_hists.hdEdx->Write();
//   for (std::map<int,std::pair<std::string,int> >::iterator pdgIt = PDGMap.begin(); pdgIt != PDGMap.end(); ++pdgIt)
//     bg_hists.hdEdxBG[pdgIt->second.second]->Write();
//   bg_hists.hAngle->Write();
//   bg_hists.hFrontSeparation->Write();
//   bg_hists.hBackSeparation->Write();
//   //bg_hists.hSeparationTrack->Write();
//   bg_hists.hCorrectTracks->Write();
//   bg_hists.hIncorrectTracks->Write();

//   // Get a canvas to whack everything on
//   TCanvas* canv = new TCanvas("canv","",800,600);

//   canv->cd();
//   canv->Clear();
//   TEfficiency* hEfficiency = new TEfficiency(*sig_hists.hNuEnergyCut,*sig_hists.hNuEnergyAll);
//   hEfficiency->SetLineColor(1);
//   hEfficiency->SetFillColor(1);
//   hEfficiency->SetFillStyle(3003);
//   TH1D* hAllSelected = (TH1D*)sig_hists.hNuEnergyCut->Clone("AllSelected");
//   hAllSelected->Add(bg_hists.hNuEnergyCut);
//   TEfficiency* hPurity = new TEfficiency(*sig_hists.hNuEnergyCut,*hAllSelected);
//   hPurity->SetLineColor(2);
//   hPurity->SetFillColor(2);
//   hPurity->SetFillStyle(3003);
//   hEfficiency->Draw();
//   hPurity->Draw("same");
//   canv->Write("EffPur");

//   canv->cd();
//   canv->Clear();
//   sig_hists.hCorrectTracks->SetMarkerColor(1);
//   sig_hists.hCorrectTracks->Draw();
//   outFile.cd();
//   canv->Write("CorrectSelection");

//   canv->cd();
//   canv->Clear();
//   bg_hists.hIncorrectTracks->SetLineColor(1);
//   bg_hists.hIncorrectTracks->SetFillColor(1);
//   bg_hists.hIncorrectTracks->SetFillStyle(3003);
//   bg_hists.hIncorrectTracks->Draw();
//   outFile.cd();
//   canv->Write("BackgroundTracks");

//   canv->cd();
//   canv->Clear();
//   sig_hists.hTrackLengthLong->SetLineColor(1);
//   sig_hists.hTrackLengthLong->SetFillColor(1);
//   sig_hists.hTrackLengthLong->SetFillStyle(3003);
//   bg_hists.hTrackLengthLong->SetLineColor(2);
//   bg_hists.hTrackLengthLong->SetFillColor(2);
//   bg_hists.hTrackLengthLong->SetFillStyle(3003);
//   sig_hists.hTrackLengthLong->Draw();
//   bg_hists.hTrackLengthLong->Draw("same");
//   outFile.cd();
//   canv->Write("TrackLengthLong");  

//   canv->cd();
//   canv->Clear();
//   sig_hists.hTrackLengthShort->SetLineColor(1);
//   sig_hists.hTrackLengthShort->SetFillColor(1);
//   sig_hists.hTrackLengthShort->SetFillStyle(3003);
//   bg_hists.hTrackLengthShort->SetLineColor(2);
//   bg_hists.hTrackLengthShort->SetFillColor(2);
//   bg_hists.hTrackLengthShort->SetFillStyle(3003);
//   bg_hists.hTrackLengthShort->Draw();
//   sig_hists.hTrackLengthShort->Draw("same");
//   outFile.cd();
//   canv->Write("TrackLengthShort");
  
//   canv->cd();
//   canv->Clear();
//   sig_hists.hEnergyLong->SetLineColor(1);
//   sig_hists.hEnergyLong->SetFillColor(1);
//   sig_hists.hEnergyLong->SetFillStyle(3003);
//   bg_hists.hEnergyLong->SetLineColor(2);
//   bg_hists.hEnergyLong->SetFillColor(2);
//   bg_hists.hEnergyLong->SetFillStyle(3003);
//   sig_hists.hEnergyLong->Draw();
//   bg_hists.hEnergyLong->Draw("same");
//   outFile.cd();
//   canv->Write("EnergyLong");

//   canv->cd();
//   canv->Clear();
//   sig_hists.hEnergyShort->SetLineColor(1);
//   sig_hists.hEnergyShort->SetFillColor(1);
//   sig_hists.hEnergyShort->SetFillStyle(3003);
//   bg_hists.hEnergyShort->SetLineColor(2);
//   bg_hists.hEnergyShort->SetFillColor(2);
//   bg_hists.hEnergyShort->SetFillStyle(3003);
//   sig_hists.hEnergyShort->Draw();
//   bg_hists.hEnergyShort->Draw("same");
//   outFile.cd();
//   canv->Write("EnergyShort");

//   canv->cd();
//   canv->Clear();
//   sig_hists.hdEdx->SetLineColor(1);
//   sig_hists.hdEdx->SetFillColor(1);
//   sig_hists.hdEdx->SetFillStyle(3003);
//   THStack* hdEdxStack = new THStack("dEdxStack",";dEdx (MeV/cm);");
//   TLegend* dEdxStackLeg = new TLegend(0.68,0.58,0.88,0.88);
//   for (std::map<int,std::pair<std::string,int> >::iterator pdgIt = PDGMap.begin(); pdgIt != PDGMap.end(); ++pdgIt) {
//     bg_hists.hdEdxBG[pdgIt->second.second]->SetLineColor(pdgIt->second.second+1);
//     bg_hists.hdEdxBG[pdgIt->second.second]->SetFillColor(pdgIt->second.second+1);
//     hdEdxStack->Add(bg_hists.hdEdxBG[pdgIt->second.second]);
//     dEdxStackLeg->AddEntry(bg_hists.hdEdxBG[pdgIt->second.second],pdgIt->second.first.c_str(),"f");
//   }
//   sig_hists.hdEdx->Draw();
//   hdEdxStack->Draw("same");
//   dEdxStackLeg->Draw();
//   outFile.cd();
//   canv->Write("dEdx");

//   canv->cd();
//   canv->Clear();
//   sig_hists.hAngle->SetLineColor(1);
//   sig_hists.hAngle->SetFillColor(1);
//   sig_hists.hAngle->SetFillStyle(3003);
//   bg_hists.hAngle->SetLineColor(2);
//   bg_hists.hAngle->SetFillColor(2);
//   bg_hists.hAngle->SetFillStyle(3003);
//   sig_hists.hAngle->Draw();
//   bg_hists.hAngle->Draw("same");
//   outFile.cd();
//   canv->Write("Angle");

//   canv->cd();
//   canv->Clear();
//   sig_hists.hFrontSeparation->SetLineColor(1);
//   sig_hists.hFrontSeparation->SetFillColor(1);
//   sig_hists.hFrontSeparation->SetFillStyle(3003);
//   bg_hists.hFrontSeparation->SetLineColor(2);
//   bg_hists.hFrontSeparation->SetFillColor(2);
//   bg_hists.hFrontSeparation->SetFillStyle(3003);
//   sig_hists.hFrontSeparation->Draw();
//   bg_hists.hFrontSeparation->Draw("same");
//   outFile.cd();
//   canv->Write("FrontSeparation");

//   canv->cd();
//   canv->Clear();
//   sig_hists.hBackSeparation->SetLineColor(1);
//   sig_hists.hBackSeparation->SetFillColor(1);
//   sig_hists.hBackSeparation->SetFillStyle(3003);
//   bg_hists.hBackSeparation->SetLineColor(2);
//   bg_hists.hBackSeparation->SetFillColor(2);
//   bg_hists.hBackSeparation->SetFillStyle(3003);
//   sig_hists.hBackSeparation->Draw();
//   bg_hists.hBackSeparation->Draw("same");
//   outFile.cd();
//   canv->Write("BackSeparation");

// }
