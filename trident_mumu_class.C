////////////////////////////////////////////////////////////////////////////////////////////
// trident_mumu.C, Mike Wallbank, January 2018
//
// Analyses signal and background trident events and saves useful information in
// tree for further evaluation.
//
// Depends on Particle and Event classes.  Event in Justo's separate repository, compile
// using cmake:
//   cd dunend/build/; cmake ..; cmake --build . --target Event;
// Particle may be compiled using root compiler:
//   root[] .L Particle.cxx+
// Both libraries should be included.
//
// Can be run interactively `root "trident_mumu.C(nevents)"`
// or compiled:
//   root[] gInterpreter->AddIncludePath("/dune/app/users/wallbank/trident");
//   root[] gSystem->Load("/dune/app/users/wallbank/trident/Particle_cxx");
//   root[] gSystem->Load("/dune/app/users/wallbank/trident/dunend/build/src/libEvent");
//   root[] gROOT->ProcessLine(".L trident_mumu.C+");
////////////////////////////////////////////////////////////////////////////////////////////

// Load dynamic library containing the definitions of the event classes
R__LOAD_LIBRARY(dunend/build/src/libEvent);
R__LOAD_LIBRARY(Particle_cxx.so);

// C++
#include <iostream>
#include <algorithm>
#include <memory>

// ROOT
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TObjArray.g"

// framework
#include "Particle.h"
#include "dunend/src/MCHit.h"
#include "dunend/src/MCTrack.h"
#include "dunend/src/MCParticle.h"
#include "dunend/src/Event.h"

const int kMaxParticles = 1000;
const bool kMakeMonitorHists = true;

class TridentMuMu {
public:

  TridentMuMu();

  void AddNu(Particle nu);
  void AddMuPlus(Particle mu);
  void AddMuMinus(Particle mu);
  void AddIdentifiedMu(Particle mu);
  double GetEventEnergy();
  double GetVertexActivity();
  Particle GetNu();
  Particle GetMuPlus();
  Particle GetMuMinus();
  Particle GetMu(int mu);
  std::vector<Particle> GetIdentifiedMus();
  std::pair<Particle, Particle> GetCandidateMus();
  double MuonAngle();
  double MuonStartSep();
  double MuonEndSep();
  double MuonParticleSep();
  int NumMuPlus();
  int NumMuMinus();
  int NumIdentifiedMus();
  void SetEventEnergy(double energy);
  void SetVertexActivity(double energy);

private:

  double fEventEnergy;
  double fVertexActivity;
  Particle fNu;
  Particle fMuPlus;
  Particle fMuMinus;
  std::vector<Particle> fIdentifiedMus;
  std::map<double,Particle> fMuMap;

};

TridentMuMu::TridentMuMu() {
  fEventEnergy = 0;
  fVertexActivity = 0;
  fNu = Particle();
  fMuPlus = Particle();
  fMuMinus = Particle();
}

void TridentMuMu::AddNu(Particle nu) {
  fNu = nu;
}

void TridentMuMu::AddMuPlus(Particle mu) {
  fMuPlus = mu;
}

void TridentMuMu::AddMuMinus(Particle mu) {
  fMuMinus = mu;
}

void TridentMuMu::AddIdentifiedMu(Particle mu) {
  fIdentifiedMus.push_back(mu);
  fMuMap[mu.Length()] = mu;
}

double TridentMuMu::GetEventEnergy() {
  return fEventEnergy;
}

Particle TridentMuMu::GetNu() {
  return fNu;
}

Particle TridentMuMu::GetMuPlus() {
  return fMuPlus;
}

Particle TridentMuMu::GetMuMinus() {
  return fMuMinus;
}

Particle TridentMuMu::GetMu(int mu) {
  return fIdentifiedMus.at(mu);
}

std::pair<Particle, Particle> TridentMuMu::GetCandidateMus() {
  std::pair<Particle, Particle> candidateMus;
  if (this->NumIdentifiedMus() < 2) {
    std::cout << "Warning: trying to obtain candidate muons from fewer than two identified muons" << std::endl;
    return candidateMus;
  }
  candidateMus.first = fMuMap.rbegin()->second;
  candidateMus.second = std::next(fMuMap.rbegin())->second;
  return candidateMus;
}

std::vector<Particle> TridentMuMu::GetIdentifiedMus() {
  return fIdentifiedMus;
}

double TridentMuMu::MuonAngle() {
  std::pair<Particle, Particle> candidateMus = this->GetCandidateMus();
  return candidateMus.first.Direction().Angle(candidateMus.second.Direction()) * 180/TMath::Pi();
}

double TridentMuMu::MuonStartSep() {
  std::pair<Particle, Particle> candidateMus = this->GetCandidateMus();
  return (candidateMus.first.Start()-candidateMus.second.Start()).Mag();
}

double TridentMuMu::MuonEndSep() {
  std::pair<Particle, Particle> candidateMus = this->GetCandidateMus();
  return (candidateMus.first.End()-candidateMus.second.End()).Mag();
}

double TridentMuMu::MuonParticleSep() {
  std::pair<Particle, Particle> candidateMus = this->GetCandidateMus();
  return (candidateMus.first.ParticleStart()-candidateMus.second.ParticleStart()).Mag();
}

int TridentMuMu::NumMuPlus() {
  return fMuPlus.ID() == -1 ? 0 : 1;
}

int TridentMuMu::NumMuMinus() {
  return fMuMinus.ID() == -1 ? 0 : 1;
}

int TridentMuMu::NumIdentifiedMus() {
  return fIdentifiedMus.size();
}

void TridentMuMu::SetEventEnergy(double energy) {
  fEventEnergy = energy;
}

void TridentMuMu::SetVertexActivity(double energy) {
  fVertexActivity = energy;
}

std::map<int,std::pair<std::string,int> > PDGMap = {{13,{"#mu",0}},{211,{"#pi",1}},{11,{"e",2}},{2212,{"p",3}},{321,{"K",4}}};

/// Define all data products which may be made to monitor more specific details
TObjArray* DefineDataProducts();

/// Process the events in one of the data samples
std::vector<std::unique_ptr<TridentMuMu> > ProcessEvents(TTree* tree,
							 unsigned int n_events,
							 unsigned int nskip,
							 std::string label,
							 TObjArray* monitorHists);

/// Returns the length of the track
double TrackLength(const MCTrack& track);

/// Returns the deposited energy of the track
double TrackDepositedEnergy(const MCTrack& track);

/// Returns the dE/dx from the first `trackLenght` cm of the track
double TrackdEdx(const MCTrack& track, double trackLength);

/// Returns the direction of the track
TVector3 TrackDirection(const MCTrack& track);

/// Returns the start of the track
TVector3 TrackStart(const MCTrack& track);

/// Returns the end of the track
TVector3 TrackEnd(const MCTrack& track);

/// Write out the event into a TTree for external analysis
void WriteData(std::string outFile,
	       const std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > >& tridents,
	       TObjArray* monitorHists);

/// Run main function
void trident_mumu(int n_events = -1, unsigned int nskip = 0, std::string outFile = "TridentMuMuOut.root") {

  gStyle->SetOptStat(0);

  // gROOT->ProcessLine(".L Particle.cxx+");
  // gSystem->AddIncludePath("dunend/src/");
  // gSystem->Load("dunend/build/src/libEvent");
  // gSystem->Load("Particle_cxx");

  // Load file and tree
  TFile file_sig("/pnfs/dune/persistent/users/jmalbos/Trident/data/sim/mumu/mumu.coh.sgn.00.g4.root","READ");
  TTree* tree_sig = (TTree*)file_sig.Get("Event");
  TChain* tree_bg = new TChain("Event");
  for (unsigned int bg_file = 0; bg_file < 100; ++bg_file)
    tree_bg->Add(Form("/pnfs/dune/persistent/users/jmalbos/Trident/data/sim/mumu/mumu.bkg.%02d.g4.root",bg_file));
  // TFile file_bg("/pnfs/dune/persistent/users/jmalbos/Trident/data/sim/mumu/mumu.bkg.00.root");
  // TTree* tree_bg = (TTree*)file_bg.Get("Event");

  unsigned int sig_events_to_process = n_events == -1 ? tree_sig->GetEntries()
    : TMath::Min((unsigned int)nskip+n_events, (unsigned int)tree_sig->GetEntries());
  unsigned int bg_events_to_process = n_events == -1 ? tree_bg->GetEntries()
    : TMath::Min((unsigned int)nskip+n_events, (unsigned int)tree_bg->GetEntries());

  std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > > allTridents;
  // Any histograms to fill
  TObjArray* monitorHists;
  if (kMakeMonitorHists)
    monitorHists = DefineDataProducts();

  allTridents["Signal"]     = std::move(ProcessEvents(tree_sig, sig_events_to_process, nskip, "Signal", monitorHists));
  allTridents["Background"] = std::move(ProcessEvents(tree_bg,  bg_events_to_process,  nskip, "Background", monitorHists));

  WriteData(outFile, allTridents, monitorHists);

  file_sig.Close();
  //file_bg.Close();

}

TObjArray* DefineDataProducts() {

  TObjArray* objs = new TObjArray();
  TH1D *hTrackLengthSig		= new TH1D("TrackLengthSig",";Track length (cm);",100,0,500);
  TH1D *hTrackEnergySig		= new TH1D("TrackEnergySig",";Track energy (MeV);",100,0,1000);
  TH1D *hLongTrackLengthSig	= new TH1D("LongTrackLengthSig",";Track length (cm);",100,0,500);
  TH1D *hLongTrackEnergySig	= new TH1D("LongTrackEnergySig",";Track energy (MeV);",100,0,1000);
  TH1D *hShortTrackLengthSig	= new TH1D("ShortTrackLengthSig",";Track length (cm);",100,0,500);
  TH1D *hShortTrackEnergySig	= new TH1D("ShortTrackEnergySig",";Track energy (MeV);",100,0,1000);
  TH1D *hTrackLengthBg		= new TH1D("TrackLengthBg",";Track length (cm);",100,0,500);
  TH1D *hTrackEnergyBg		= new TH1D("TrackEnergyBg",";Track energy (MeV);",100,0,1000);
  TH1D *hLongTrackLengthBg	= new TH1D("LongTrackLengthBg",";Track length (cm);",100,0,500);
  TH1D *hLongTrackEnergyBg	= new TH1D("LongTrackEnergyBg",";Track energy (MeV);",100,0,1000);
  TH1D *hShortTrackLengthBg	= new TH1D("ShortTrackLengthBg",";Track length (cm);",100,0,500);
  TH1D *hShortTrackEnergyBg	= new TH1D("ShortTrackEnergyBg",";Track energy (MeV);",100,0,1000);
  TH2D* hVertexActivitySig	= new TH2D("VertexActivitySig",";Distance from vertex (cm);Deposited energy (MeV);",40,0,20,40,0,400);
  TH2D* hVertexActivityBg       = new TH2D("VertexActivityBg",";Distance from vertex (cm);Deposited energy (MeV);",40,0,20,40,0,400);
  TGraph* gVertexActivitySig	= new TGraph();
  gVertexActivitySig->SetName("VertexActivitySig");
  TGraph* gVertexActivityBg	= new TGraph();
  gVertexActivityBg->SetName("VertexActivityBg");

  objs->Add(hTrackLengthSig);
  objs->Add(hTrackEnergySig);
  objs->Add(hLongTrackLengthSig);
  objs->Add(hLongTrackEnergySig);
  objs->Add(hShortTrackLengthSig);
  objs->Add(hShortTrackEnergySig);
  objs->Add(hTrackLengthBg);
  objs->Add(hTrackEnergyBg);
  objs->Add(hLongTrackLengthBg);
  objs->Add(hLongTrackEnergyBg);
  objs->Add(hShortTrackLengthBg);
  objs->Add(hShortTrackEnergyBg);
  objs->Add(hVertexActivitySig);
  objs->Add(hVertexActivityBg);
  objs->Add(gVertexActivitySig);
  objs->Add(gVertexActivityBg);

}

std::vector<std::unique_ptr<TridentMuMu> > ProcessEvents(TTree* tree,
							 unsigned int n_events,
							 unsigned int nskip,
							 std::string label,
							 TObjArray* monitorHists) {

  // Sync tree with empty event
  Event* evt = 0;
  tree->SetBranchAddress("Event", &evt);

  std::vector<std::unique_ptr<TridentMuMu> > tridentMuMus;

  for (unsigned int event = nskip; event < n_events; ++event) {

    if (event % 100 == 0)
      std::cout << "\r" << label << ": processing event " << event << "/" << n_events
		<< " (" << (event)*100/(double)n_events << "%)" << std::flush;

    tree->GetEntry(event);

    // Save event information
    std::map<int,Particle*> particleMap;
    std::vector<MCHit> hits;
    double eventEnergy = 0.;
    TVector3 eventVertex = TVector3(0,0,0);//evt->Vertex();

    for (std::vector<MCParticle>::iterator particleIt = evt->MCParticleContainer.begin();
	 particleIt != evt->MCParticleContainer.end();
	 ++particleIt) {
      //particleMap[particleIt->GetMCID()] = Particle(particleIt->GetMCID());
      Particle* particle = new Particle(particleIt->GetMCID());
      if (eventVertex == TVector3(0,0,0))
	eventVertex = particleIt->GetInitialPositionAndTime().Vect();
      particle->SetParticleStart(particleIt->GetInitialPositionAndTime().Vect());
      particle->SetPDG(particleIt->GetPDGCode());
      particle->SetEnergy(10); //particleIt->GetEnergy();
      particle->SetAssociatedTrack(false);
      particleMap[particleIt->GetMCID()] = particle;
    }
    for (std::vector<MCTrack>::iterator trackIt = evt->MCTrackContainer.begin();
	 trackIt != evt->MCTrackContainer.end();
	 ++trackIt) {
      particleMap[trackIt->GetMCID()]->SetAssociatedTrack(true);
      particleMap[trackIt->GetMCID()]->SetDepositedEnergy(TrackDepositedEnergy(*trackIt));
      particleMap[trackIt->GetMCID()]->SetdEdx(TrackdEdx(*trackIt, 100.));
      particleMap[trackIt->GetMCID()]->SetLength(TrackLength(*trackIt));
      particleMap[trackIt->GetMCID()]->SetDirection(TrackDirection(*trackIt));
      particleMap[trackIt->GetMCID()]->SetStart(TrackStart(*trackIt));
      particleMap[trackIt->GetMCID()]->SetEnd(TrackEnd(*trackIt));
      eventEnergy += TrackDepositedEnergy(*trackIt);
      std::vector<MCHit> trackHits = track.MCHitContainer;
      for (std::vector<MCHit>::const_iterator trackHitIt = trackHits.begin();
	   trackHitIt != trackHits.end();
	   ++trackHitIt)
	trackHits.push_back(*trackHitIt);
    }

    // Save all information about the trident candidate event
    TridentMuMu tridentMuMu;
    tridentMuMu.SetEventEnergy(eventEnergy);
    tridentMuMu.SetEventVertex(eventVertex);

    // Look through hits
    double vertexActivity = 0.;
    for (std::vector<MCHit>::const_iterator hitIt = hits.begin(); hitIt != hits.end; ++hitIt) {
      TVector3 pos = hitIt->GetInitialPositionAndTime().Vect();
      if ((pos-eventVertex).Mag() < 5.)
	vertexActivity += hitIt->GetEnergyDeposit();
      if (kMakeMonitorHists) {
	if (label == "Signal") {
	  ((TH2D*)monitorHists->FindObject("VertexActivitySig"))->Fill((pos-eventVertex).Mag(), hitIt->GetEnergyDeposit());
	  ((TGraph*)monitorHists->FindObject("VertexActivitySig"))
	    ->SetPoint(monitorHists->FindObject("VertexActivitySig")->GetN(), (pos-eventVertex).Mag(), hitIt->GetEnergyDeposit());
	} else {
	  ((TH2D*)monitorHists->FindObject("VertexActivityBg:"))->Fill((pos-eventVertex).Mag(), hitIt->GetEnergyDeposit());
	  ((TGraph*)monitorHists->FindObject("VertexActivityBg"))
	    ->SetPoint(monitorHists->FindObject("VertexActivityBg")->GetN(), (pos-eventVertex).Mag(), hitIt->GetEnergyDeposit());
	}
      }
    }
    tridentMuMu->SetVertexActivity(vertexActivity);

    // Look through particles
    std::map<double, Particle*> particleLength;
    for (std::map<int,Particle*>::const_iterator particleIt = particleMap.begin();
	 particleIt != particleMap.end();
	 ++particleIt) {

      // Signal info
      if (label == "Signal") {
	if (particleIt->second->PDG() == 14)
	  tridentMuMu.AddNu(*(particleIt->second));
	if (particleIt->second->PDG() == 13 and tridentMuMu.NumMuPlus() == 0)
	  tridentMuMu.AddMuPlus(*(particleIt->second));
	if (particleIt->second->PDG() == -13 and tridentMuMu.NumMuMinus() == 0)
	  tridentMuMu.AddMuMinus(*(particleIt->second));
      }

      particleLength[particleIt->second->Length()] = particleIt->second;

      // Selection info
      if (particleIt->second->AssociatedTrack() and particleIt->second->Length() > 50 and particleIt->second->Energy() > 100)
	tridentMuMu.AddIdentifiedMu(*(particleIt->second));

    }

    if (kMakeMonitorHists and particleLength.size() >= 2) {
      Particle* longestParticle = particleLength.rbegin()->second;
      Particle* nextLongestParticle = std::next(particleLength.rbegin())->second;
      if (label == "Signal") {
	((TH1D*)monitorHists->FindObject("TrackLengthSig"))->Fill(longestParticle->Length());
	((TH1D*)monitorHists->FindObject("TrackLengthSig"))->Fill(nextLongestParticle->Length());
	((TH1D*)monitorHists->FindObject("TrackEnergySig"))->Fill(longestParticle->DepositedEnergy());
	((TH1D*)monitorHists->FindObject("TrackEnergySig"))->Fill(nextLongestParticle->DepositedEnergy());
	((TH1D*)monitorHists->FindObject("LongTrackLengthSig"))->Fill(longestParticle->Length());
	((TH1D*)monitorHists->FindObject("LongTrackEnergySig"))->Fill(longestParticle->DepositedEnergy());
	((TH1D*)monitorHists->FindObject("ShortTrackLengthSig"))->Fill(nextLongestParticle->Length());
	((TH1D*)monitorHists->FindObject("ShortTrackEnergySig"))->Fill(nextLongestParticle->DepositedEnergy());
      } else {
	((TH1D*)monitorHists->FindObject("TrackLengthBg"))->Fill(longestParticle->Length());
	((TH1D*)monitorHists->FindObject("TrackLengthBg"))->Fill(nextLongestParticle->Length());
	((TH1D*)monitorHists->FindObject("TrackEnergyBg"))->Fill(longestParticle->DepositedEnergy());
	((TH1D*)monitorHists->FindObject("TrackEnergyBg"))->Fill(nextLongestParticle->DepositedEnergy());
	((TH1D*)monitorHists->FindObject("LongTrackLengthBg"))->Fill(longestParticle->Length());
	((TH1D*)monitorHists->FindObject("LongTrackEnergyBg"))->Fill(longestParticle->DepositedEnergy());
	((TH1D*)monitorHists->FindObject("ShortTrackLengthBg"))->Fill(nextLongestParticle->Length());
	((TH1D*)monitorHists->FindObject("ShortTrackEnergyBg"))->Fill(nextLongestParticle->DepositedEnergy());
      }
    }

    tridentMuMus.push_back(std::make_unique<TridentMuMu>(tridentMuMu));

    for (std::map<int,Particle*>::iterator particleIt = particleMap.begin();
	 particleIt != particleMap.end();
	 ++particleIt)
      delete particleIt->second;

  } // events

  std::cout << std::endl;
  return tridentMuMus;

}

double TrackLength(const MCTrack& track) {

  std::vector<MCHit> trackHits = track.MCHitContainer;
  MCHit firstHit = *trackHits.begin();
  MCHit lastHit = *trackHits.rbegin();

  double trackLength = (firstHit.GetStartPositionAndTime().Vect() - lastHit.GetStopPositionAndTime().Vect()).Mag();

  return trackLength/10.;

}

double TrackDepositedEnergy(const MCTrack& track) {

  std::vector<MCHit> trackHits = track.MCHitContainer;

  double energy = 0;
  for (std::vector<MCHit>::iterator hitIt = trackHits.begin(); hitIt != trackHits.end(); ++hitIt)
    energy += hitIt->GetEnergyDeposit();

  return energy;

}

double TrackdEdx(const MCTrack& track, double dEdxlength) {

  std::vector<MCHit> trackHits = track.MCHitContainer;

  if (trackHits.size() < 2)
    return trackHits.size() ? trackHits.begin()->GetEnergyDeposit() : 0;

  double dE = 0., dx = 0.;
  TVector3 lastHitPos = trackHits.begin()->GetStartPositionAndTime().Vect();
  int nhit = 0;
  for (std::vector<MCHit>::iterator hitIt = trackHits.begin(); hitIt != trackHits.end() and dx < dEdxlength; ++hitIt) {
    TVector3 hitPos = hitIt->GetStartPositionAndTime().Vect();
    dx += (hitPos - lastHitPos).Mag();
    dE += hitIt->GetEnergyDeposit();
    lastHitPos = hitPos;
    ++nhit;
  }

  double dEdx = 0.;

  if (dx)
    dEdx = dE/dx;

  //std::cout << "dE/dx = " << dEdx << " found from " << nhit << " hits" << std::endl;

  return dEdx;

}

TVector3 TrackStart(const MCTrack& track) {

  std::vector<MCHit> trackHits = track.MCHitContainer;

  return trackHits.begin()->GetStartPositionAndTime().Vect();

}

TVector3 TrackEnd(const MCTrack& track) {

  std::vector<MCHit> trackHits = track.MCHitContainer;

  return trackHits.rbegin()->GetStopPositionAndTime().Vect();

}

TVector3 TrackDirection(const MCTrack& track) {

  TVector3 direction;
  
  std::vector<MCHit> trackHits = track.MCHitContainer;
  
  if (trackHits.size() < 2)
    direction = TVector3(-999.,-999.,-999.);

  else if (trackHits.size() == 2) {
    TVector3 hit1 = trackHits[0].GetStartPositionAndTime().Vect();
    TVector3 hit2 = trackHits[1].GetStartPositionAndTime().Vect();
    direction = (hit2-hit1).Unit();
  }

  else if (trackHits.size() > 2) {
    TVector3 hit1 = trackHits[0].GetStartPositionAndTime().Vect();
    TVector3 hit2 = trackHits[trackHits.size()-1].GetStartPositionAndTime().Vect();
    direction = (hit2-hit1).Unit();
  }

  return direction;

}

void WriteData(std::string outFile,
	       const std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > >& tridents,
	       TObjArray* monitorHists) {

  TFile* outF = new TFile(outFile.c_str(), "RECREATE");

  for (std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > >::const_iterator
	 allTridentIt = tridents.begin();
       allTridentIt != tridents.end();
       ++allTridentIt) {

    TTree* outT = new TTree(("TridentMuMu"+allTridentIt->first).c_str(), "");//,"TridentMuMu"+tridentIt->first);

    // Tree variables
    TridentMuMu tTridentMuMu;
    double tEventEnergy;
    int tNumIdentifiedMus;
    int tNumParticles;
    // Nu
    int tNuID, tNuPDG;
    double tNuStartX, tNuStartY, tNuStartZ;
    double tNuDirectionX, tNuDirectionY, tNuDirectionZ;
    // Particles
    int tID[kMaxParticles], tPDG[kMaxParticles];
    double tEnergy[kMaxParticles], tdEdx[kMaxParticles],
      tLength[kMaxParticles],
      tDirectionX[kMaxParticles], tDirectionY[kMaxParticles], tDirectionZ[kMaxParticles],
      tStartX[kMaxParticles], tStartY[kMaxParticles], tStartZ[kMaxParticles],
      tEndX[kMaxParticles], tEndY[kMaxParticles], tEndZ[kMaxParticles];
    // Mu plus (true)
    int tMuPlusID, tMuPlusPDG;
    double tMuPlusEnergy, tMuPlusdEdx,
      tMuPlusLength,
      tMuPlusDirectionX, tMuPlusDirectionY, tMuPlusDirectionZ,
      tMuPlusStartX, tMuPlusStartY, tMuPlusStartZ,
      tMuPlusEndX, tMuPlusEndY, tMuPlusEndZ;
    // Mu minus (true)
    int tMuMinusID, tMuMinusPDG;
    double tMuMinusEnergy, tMuMinusdEdx,
      tMuMinusLength,
      tMuMinusDirectionX, tMuMinusDirectionY, tMuMinusDirectionZ,
      tMuMinusStartX, tMuMinusStartY, tMuMinusStartZ,
      tMuMinusEndX, tMuMinusEndY, tMuMinusEndZ;
    // Short muon (selected)
    int tShortMuonID, tShortMuonPDG;
    double tShortMuonEnergy, tShortMuondEdx,
      tShortMuonLength,
      tShortMuonDirectionX, tShortMuonDirectionY, tShortMuonDirectionZ,
      tShortMuonStartX, tShortMuonStartY, tShortMuonStartZ,
      tShortMuonEndX, tShortMuonEndY, tShortMuonEndZ;
    // Long muon (selected)
    int tLongMuonID, tLongMuonPDG;
    double tLongMuonEnergy, tLongMuondEdx,
      tLongMuonLength,
      tLongMuonDirectionX, tLongMuonDirectionY, tLongMuonDirectionZ,
      tLongMuonStartX, tLongMuonStartY, tLongMuonStartZ,
      tLongMuonEndX, tLongMuonEndY, tLongMuonEndZ;
    // Selected muons
    double tAngle, tFrontSeparation, tBackSeparation;

    // Set branch vars
    // Event
    //outT->Branch("TridentMuMu",        tTridentMuMu);
    outT->Branch("EventEnergy",      &tEventEnergy);
    outT->Branch("NumIdentifiedMus", &tNumIdentifiedMus);
    outT->Branch("NumParticles",     &tNumParticles);
    // Nu
    outT->Branch("NuID",         &tNuID);
    outT->Branch("NuPDG",        &tNuPDG);
    outT->Branch("NuStartX",     &tNuStartX);
    outT->Branch("NuStartY",     &tNuStartY);
    outT->Branch("NuStartZ",     &tNuStartZ);
    outT->Branch("NuDirectionX", &tNuDirectionX);
    outT->Branch("NuDirectionY", &tNuDirectionY);
    outT->Branch("NuDirectionZ", &tNuDirectionZ);
    // All particles
    outT->Branch("ID",         &tID,         "ID[NumParticles]/I");
    outT->Branch("PDG",        &tPDG,        "PDG[NumParticles]/I");
    outT->Branch("Energy",     &tEnergy,     "PDG[NumParticles]/F");
    outT->Branch("dEdx",       &tdEdx,       "dEdx[NumParticles]/F");
    outT->Branch("Length",     &tLength,     "Length[NumParticles]/F");
    outT->Branch("DirectionX", &tDirectionX, "MuPlusDirectionX[NumParticles]/F");
    outT->Branch("DirectionY", &tDirectionY, "MuPlusDirectionY[NumParticles]/F");
    outT->Branch("DirectionZ", &tDirectionZ, "MuPlusDirectionZ[NumParticles]/F");
    outT->Branch("StartX",     &tStartX,     "MuPlusStartX[NumParticles]/F");
    outT->Branch("StartY",     &tStartY,     "MuPlusStartY[NumParticles]/F");
    outT->Branch("StartZ",     &tStartZ,     "MuPlusStartZ[NumParticles]/F");
    outT->Branch("EndX",       &tEndX,       "MuPlusEndX[NumParticles]/F");
    outT->Branch("EndY",       &tEndY,       "MuPlusEndY[NumParticles]/F");
    outT->Branch("EndZ",       &tEndZ,       "MuPlusEndZ[NumParticles]/F");
    // Mu plus (true)
    outT->Branch("MuPlusID",         &tMuPlusID);
    outT->Branch("MuPlusPDG",        &tMuPlusPDG);
    outT->Branch("MuPlusEnergy",     &tMuPlusEnergy);
    outT->Branch("MuPlusdEdx",       &tMuPlusdEdx);
    outT->Branch("MuPlusLength",     &tMuPlusLength);
    outT->Branch("MuPlusDirectionX", &tMuPlusDirectionX);
    outT->Branch("MuPlusDirectionY", &tMuPlusDirectionY);
    outT->Branch("MuPlusDirectionZ", &tMuPlusDirectionZ);
    outT->Branch("MuPlusStartX",     &tMuPlusStartX);
    outT->Branch("MuPlusStartY",     &tMuPlusStartY);
    outT->Branch("MuPlusStartZ",     &tMuPlusStartZ);
    outT->Branch("MuPlusEndX",       &tMuPlusEndX);
    outT->Branch("MuPlusEndY",       &tMuPlusEndY);
    outT->Branch("MuPlusEndZ",       &tMuPlusEndZ);
    // Mu minus (true)
    outT->Branch("MuMinusID",         &tMuMinusID);
    outT->Branch("MuMinusPDG",        &tMuMinusPDG);
    outT->Branch("MuMinusEnergy",     &tMuMinusEnergy);
    outT->Branch("MuMinusdEdx",       &tMuMinusdEdx);
    outT->Branch("MuMinusLength",     &tMuMinusLength);
    outT->Branch("MuMinusDirectionX", &tMuMinusDirectionX);
    outT->Branch("MuMinusDirectionY", &tMuMinusDirectionY);
    outT->Branch("MuMinusDirectionZ", &tMuMinusDirectionZ);
    outT->Branch("MuMinusStartX",     &tMuMinusStartX);
    outT->Branch("MuMinusStartY",     &tMuMinusStartY);
    outT->Branch("MuMinusStartZ",     &tMuMinusStartZ);
    outT->Branch("MuMinusEndX",       &tMuMinusEndX);
    outT->Branch("MuMinusEndY",       &tMuMinusEndY);
    outT->Branch("MuMinusEndZ",       &tMuMinusEndZ);
    // Short muon (selected )
    outT->Branch("ShortMuonID",         &tShortMuonID);
    outT->Branch("ShortMuonPDG",        &tShortMuonPDG);
    outT->Branch("ShortMuonEnergy",     &tShortMuonEnergy);
    outT->Branch("ShortMuondEdx",       &tShortMuondEdx);
    outT->Branch("ShortMuonLength",     &tShortMuonLength);
    outT->Branch("ShortMuonDirectionX", &tShortMuonDirectionX);
    outT->Branch("ShortMuonDirectionY", &tShortMuonDirectionY);
    outT->Branch("ShortMuonDirectionZ", &tShortMuonDirectionZ);
    outT->Branch("ShortMuonStartX",     &tShortMuonStartX);
    outT->Branch("ShortMuonStartY",     &tShortMuonStartY);
    outT->Branch("ShortMuonStartZ",     &tShortMuonStartZ);
    outT->Branch("ShortMuonEndX",       &tShortMuonEndX);
    outT->Branch("ShortMuonEndY",       &tShortMuonEndY);
    outT->Branch("ShortMuonEndZ",       &tShortMuonEndZ);
    // Long muon (selected)
    outT->Branch("LongMuonID",         &tLongMuonID);
    outT->Branch("LongMuonPDG",        &tLongMuonPDG);
    outT->Branch("LongMuonEnergy",     &tLongMuonEnergy);
    outT->Branch("LongMuondEdx",       &tLongMuondEdx);
    outT->Branch("LongMuonLength",     &tLongMuonLength);
    outT->Branch("LongMuonDirectionX", &tLongMuonDirectionX);
    outT->Branch("LongMuonDirectionY", &tLongMuonDirectionY);
    outT->Branch("LongMuonDirectionZ", &tLongMuonDirectionZ);
    outT->Branch("LongMuonStartX",     &tLongMuonStartX);
    outT->Branch("LongMuonStartY",     &tLongMuonStartY);
    outT->Branch("LongMuonStartZ",     &tLongMuonStartZ);
    outT->Branch("LongMuonEndX",       &tLongMuonEndX);
    outT->Branch("LongMuonEndY",       &tLongMuonEndY);
    outT->Branch("LongMuonEndZ",       &tLongMuonEndZ);
    // Selected muons
    outT->Branch("Angle",           &tAngle);
    outT->Branch("FrontSeparation", &tFrontSeparation);
    outT->Branch("BackSeparation",  &tBackSeparation);

    for (std::vector<std::unique_ptr<TridentMuMu> >::const_iterator
	   tridentIt = allTridentIt->second.begin();
	 tridentIt != allTridentIt->second.end();
	 ++tridentIt) {

      //tTridentMuMu = *(*tridentIt).release();
      tEventEnergy = (*tridentIt)->GetEventEnergy();
      tNumIdentifiedMus = (*tridentIt)->NumIdentifiedMus();

      // Get particle objects
      Particle nu = (*tridentIt)->GetNu();
      std::vector<Particle> mus = (*tridentIt)->GetIdentifiedMus();
      Particle muPlus = (*tridentIt)->GetMuPlus();
      Particle muMinus = (*tridentIt)->GetMuMinus();

      Particle longMuon, shortMuon;
      if (tNumIdentifiedMus > 1) {
	// In the instance of multiple identified muons, need to pick two
	// Right now, take the longest track
	std::map<double,Particle> muLength;
	for (std::vector<Particle>::iterator muIt = mus.begin(); muIt != mus.end(); ++muIt)
	  muLength[muIt->Length()] = *muIt;
	std::map<double,Particle>::reverse_iterator muIt = muLength.rbegin();
	longMuon = muIt->second;
	shortMuon = std::next(muIt)->second;
      }
      else if (tNumIdentifiedMus == 1)
	longMuon = mus.at(0);

      tNuID         = nu.ID();
      tNuPDG        = nu.PDG();
      tNuStartX     = nu.Start().X();
      tNuStartY     = nu.Start().Y();
      tNuStartZ     = nu.Start().Z();
      tNuDirectionX = nu.Direction().X();
      tNuDirectionY = nu.Direction().Y();
      tNuDirectionZ = nu.Direction().Z();

      tNumParticles = 0;

      tMuPlusID         = muPlus.ID();
      tMuPlusPDG        = muPlus.PDG();
      tMuPlusEnergy     = muPlus.Energy();
      tMuPlusdEdx       = muPlus.dEdx();
      tMuPlusLength     = muPlus.Length();
      tMuPlusDirectionX = muPlus.Direction().X();
      tMuPlusDirectionY = muPlus.Direction().Y();
      tMuPlusDirectionZ = muPlus.Direction().Z();
      tMuPlusStartX     = muPlus.Start().X();
      tMuPlusStartY     = muPlus.Start().Y();
      tMuPlusStartZ     = muPlus.Start().Z();
      tMuPlusEndX       = muPlus.End().X();
      tMuPlusEndY       = muPlus.End().Y();
      tMuPlusEndZ       = muPlus.End().Z();

      tMuMinusID         = muMinus.ID();
      tMuMinusPDG        = muMinus.PDG();
      tMuMinusEnergy     = muMinus.Energy();
      tMuMinusdEdx       = muMinus.dEdx();
      tMuMinusLength     = muMinus.Length();
      tMuMinusDirectionX = muMinus.Direction().X();
      tMuMinusDirectionY = muMinus.Direction().Y();
      tMuMinusDirectionZ = muMinus.Direction().Z();
      tMuMinusStartX     = muMinus.Start().X();
      tMuMinusStartY     = muMinus.Start().Y();
      tMuMinusStartZ     = muMinus.Start().Z();
      tMuMinusEndX       = muMinus.End().X();
      tMuMinusEndY       = muMinus.End().Y();
      tMuMinusEndZ       = muMinus.End().Z();

      tLongMuonID         = longMuon.ID();
      tLongMuonPDG        = longMuon.PDG();
      tLongMuonEnergy     = longMuon.Energy();
      tLongMuondEdx       = longMuon.dEdx();
      tLongMuonLength     = longMuon.Length();
      tLongMuonDirectionX = longMuon.Direction().X();
      tLongMuonDirectionY = longMuon.Direction().Y();
      tLongMuonDirectionZ = longMuon.Direction().Z();
      tLongMuonStartX     = longMuon.Start().X();
      tLongMuonStartY     = longMuon.Start().Y();
      tLongMuonStartZ     = longMuon.Start().Z();
      tLongMuonEndX       = longMuon.End().X();
      tLongMuonEndY       = longMuon.End().Y();
      tLongMuonEndZ       = longMuon.End().Z();

      tShortMuonID         = shortMuon.ID();
      tShortMuonPDG        = shortMuon.PDG();
      tShortMuonEnergy     = shortMuon.Energy();
      tShortMuondEdx       = shortMuon.dEdx();
      tShortMuonLength     = shortMuon.Length();
      tShortMuonDirectionX = shortMuon.Direction().X();
      tShortMuonDirectionY = shortMuon.Direction().Y();
      tShortMuonDirectionZ = shortMuon.Direction().Z();
      tShortMuonStartX     = shortMuon.Start().X();
      tShortMuonStartY     = shortMuon.Start().Y();
      tShortMuonStartZ     = shortMuon.Start().Z();
      tShortMuonEndX       = shortMuon.End().X();
      tShortMuonEndY       = shortMuon.End().Y();
      tShortMuonEndZ       = shortMuon.End().Z();

      tAngle = longMuon.Direction().
	Angle(shortMuon.Direction()) * 180/TMath::Pi();
      tFrontSeparation = (longMuon.Start()-
			  shortMuon.Start()).Mag();
      tBackSeparation = (longMuon.End()-
			 shortMuon.End()).Mag();

      outT->Fill();

    } // tree

    outF->cd();
    outT->Write();
    delete outT;

  } // signal/bg

  // Write out all the monitor hists
  outF->cd();
  for (unsigned int obj = 0; obj < monitorHists->GetEntries(); ++obj)
    monitorHists->At(obj)->Write();

  outF->Close();
  delete outF;

}
