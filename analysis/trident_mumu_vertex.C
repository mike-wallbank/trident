#include <iostream>
#include <iomanip>

#include "TStyle.h"

// Style
const int kSignalColor = 1;
const int kBackgroundColor = 2;
const int kSignalFill = 3003;
const int kBackgroundFill = 3003;

void trident_mumu_vertex() {

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
  TFile* outFile = new TFile("TridentMuMuVertex.root", "RECREATE");
  TCanvas* canv = new TCanvas("canv","",800,600);
  TLegend* leg = new TLegend(0.5,0.6,0.85,0.85);

  // Energy around vertex
  TH1F* EnergyVertexSig = new TH1F("EnergyVertexSig",";Energy within 10cm;",100,0,1000);
  treeSig->Draw("VertexActivity>>EnergyVertexSig");
  TH1F* EnergyVertexBg = new TH1F("EnergyVertexBg",";Energy within 10cm;",100,0,1000);
  treeBg->Draw("VertexActivity>>EnergyVertexBg");
  canv->cd();
  canv->Clear();
  leg->Clear();
  EnergyVertexSig->SetLineColor(1);
  leg->AddEntry(EnergyVertexSig, "Signal", "l");
  EnergyVertexBg->SetLineColor(2);
  leg->AddEntry(EnergyVertexBg, "Background", "l");
  EnergyVertexBg->Draw("hist");
  EnergyVertexSig->Draw("hist same");
  leg->Draw();
  outFile->cd();
  canv->Write("EnergyVertex");

  // Energy fraction around vertex
  TH1F* EnergyVertexFracSig = new TH1F("EnergyFracVertexSig",";Fraction of Energy within 10cm;",100,0,1);
  treeSig->Draw("VertexFraction>>EnergyFracVertexSig");
  TH1F* EnergyVertexFracBg = new TH1F("EnergyFracVertexBg",";Fraction of Energy within 10cm;",100,0,1);
  treeBg->Draw("VertexFraction>>EnergyFracVertexBg");
  canv->cd();
  canv->Clear();
  leg->Clear();
  EnergyVertexFracSig->SetLineColor(1);
  leg->AddEntry(EnergyVertexFracSig, "Signal", "l");
  EnergyVertexFracBg->SetLineColor(2);
  leg->AddEntry(EnergyVertexFracBg, "Background", "l");
  EnergyVertexFracBg->Draw("hist");
  EnergyVertexFracSig->Draw("hist same");
  leg->Draw();
  outFile->cd();
  canv->Write("EnergyFracVertex");

  // Long muon energy fraction around vertex
  TH1F* LongMuonEnergyVertexFracSig = new TH1F("LongMuonEnergyVertexFracSig",";Long Identified Muon Deposited Energy within 10cm;",100,0,1);
  treeSig->Draw("LongMuonVertexFrac>>LongMuonEnergyVertexFracSig");
  TH1F* LongMuonEnergyVertexFracBg = new TH1F("LongMuonEnergyVertexFracBg",";Long Identified Muon Deposited Energy within 10cm;",100,0,1);
  treeBg->Draw("LongMuonVertexFrac>>LongMuonEnergyVertexFracBg");
  canv->cd();
  canv->Clear();
  leg->Clear();
  LongMuonEnergyVertexFracSig->SetLineColor(1);
  leg->AddEntry(LongMuonEnergyVertexFracSig, "Signal", "l");
  LongMuonEnergyVertexFracBg->SetLineColor(2);
  leg->AddEntry(LongMuonEnergyVertexFracBg, "Background", "l");
  LongMuonEnergyVertexFracBg->Draw("hist");
  LongMuonEnergyVertexFracSig->Draw("hist same");
  leg->Draw();
  outFile->cd();
  canv->Write("LongMuonEnergyVertexFrac");

  // Short muon energy fraction around vertex
  TH1F* ShortMuonEnergyVertexFracSig = new TH1F("ShortMuonEnergyVertexFracSig",";Short Identified Muon Deposited Energy within 10cm;",100,0,1);
  treeSig->Draw("ShortMuonVertexFrac>>ShortMuonEnergyVertexFracSig");
  TH1F* ShortMuonEnergyVertexFracBg = new TH1F("ShortMuonEnergyVertexFracBg",";Short Identified Muon Deposited Energy within 10cm;",100,0,1);
  treeBg->Draw("ShortMuonVertexFrac>>ShortMuonEnergyVertexFracBg");
  canv->cd();
  canv->Clear();
  leg->Clear();
  ShortMuonEnergyVertexFracSig->SetLineColor(1);
  leg->AddEntry(ShortMuonEnergyVertexFracSig, "Signal", "l");
  ShortMuonEnergyVertexFracBg->SetLineColor(2);
  leg->AddEntry(ShortMuonEnergyVertexFracBg, "Background", "l");
  ShortMuonEnergyVertexFracBg->Draw("hist");
  ShortMuonEnergyVertexFracSig->Draw("hist same");
  leg->Draw();
  outFile->cd();
  canv->Write("ShortMuonEnergyVertexFrac");

  // // Energy deposited fraction
  // canv->cd();
  // canv->Clear();
  // leg->Clear();
  // TProfile* energyDepSig = (TProfile*)fileSig->Get("VertexEnergyFracSig");
  // TProfile* energyDepBg = (TProfile*)fileBg->Get("VertexEnergyFracBg");
  // energyDepSig->SetLineColor(1);
  // energyDepSig->Draw("hist");
  // leg->AddEntry(energyDepSig, "Signal", "l");
  // energyDepBg->SetLineColor(2);
  // energyDepBg->Draw("hist same");
  // leg->AddEntry(energyDepBg, "Background", "l");
  // leg->Draw();
  // outFile->cd();
  // canv->Write("EnergyDepositedFrac");

  // outFile->Close();
  // delete outFile;

}
