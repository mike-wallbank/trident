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
#include "TObjArray.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraph.h"

// framework
#include "Particle.h"
#include "dunend/src/MCHit.h"
#include "dunend/src/MCTrack.h"
#include "dunend/src/MCParticle.h"
#include "dunend/src/Event.h"

// Run options
const int kMaxParticles = 1000;
const bool kMakeMonitorHists = true;
const bool kSaveTree = true;

// Selection options
const float kVertexDifference = 100.;

class TridentMuMu {
public:

  TridentMuMu();

  void AddNu(Particle nu);
  void AddMuPlus(Particle mu);
  void AddMuMinus(Particle mu);
  void AddTrack(Particle track);
  double GetEventEnergy();
  int GetEventNumber();
  TVector3 GetEventVertex();
  double GetVertexActivity();
  double GetVertexFraction();
  Particle GetNu();
  Particle GetMuPlus();
  Particle GetMuMinus();
  Particle GetTrack(int track);
  std::vector<Particle> GetTracks();
  std::pair<Particle, Particle> GetCandidateMus();
  double MuonAngle();
  double MuonStartSep();
  double MuonEndSep();
  double MuonParticleSep();
  int NumMuPlus();
  int NumMuMinus();
  int NumTracks();
  void SetEventEnergy(double energy);
  void SetEventNumber(int event);
  void SetEventVertex(TVector3 vertex);
  void SetVertexActivity(double energy);
  void SetVertexFraction(double vertexFraction);

private:

  int fEventNumber;
  double fEventEnergy;
  double fVertexActivity;
  double fVertexFraction;
  TVector3 fEventVertex;
  Particle fNu;
  Particle fMuPlus;
  Particle fMuMinus;
  std::vector<Particle> fTracks;
  std::map<double,Particle> fMuMap;

};

TridentMuMu::TridentMuMu() {
  fEventNumber = -1;
  fEventEnergy = 0;
  fEventVertex = TVector3(-9999.,-9999.,-9999.);
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

void TridentMuMu::AddTrack(Particle mu) {
  fTracks.push_back(mu);
  fMuMap[mu.Length()] = mu;
}

double TridentMuMu::GetEventEnergy() {
  return fEventEnergy;
}

int TridentMuMu::GetEventNumber() {
  return fEventNumber;
}

TVector3 TridentMuMu::GetEventVertex() {
  return fEventVertex;
}

double TridentMuMu::GetVertexActivity() {
  return fVertexActivity;
}

double TridentMuMu::GetVertexFraction() {
  return fVertexFraction;
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

Particle TridentMuMu::GetTrack(int track) {
  return fTracks.at(track);
}

std::pair<Particle, Particle> TridentMuMu::GetCandidateMus() {
  std::pair<Particle, Particle> candidateMus;
  if (this->NumTracks() < 2) {
    std::cout << "Warning: trying to obtain candidate muons from fewer than two tracks" << std::endl;
    return candidateMus;
  }
  candidateMus.first = fMuMap.rbegin()->second;
  candidateMus.second = std::next(fMuMap.rbegin())->second;
  return candidateMus;
}

std::vector<Particle> TridentMuMu::GetTracks() {
  return fTracks;
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

int TridentMuMu::NumTracks() {
  return fTracks.size();
}

void TridentMuMu::SetEventEnergy(double energy) {
  fEventEnergy = energy;
}

void TridentMuMu::SetEventNumber(int event) {
  fEventNumber = event;
}

void TridentMuMu::SetEventVertex(TVector3 vertex) {
  fEventVertex = vertex;
}

void TridentMuMu::SetVertexActivity(double energy) {
  fVertexActivity = energy;
}

void TridentMuMu::SetVertexFraction(double vertexFraction) {
  fVertexFraction = vertexFraction;
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

/// Returns the fraction of energy deposited by the track near the vertex
double TrackVertexFraction(const MCTrack& track, const TVector3& vertex);

/// Write out the event into a TTree for external analysis
void WriteData(TFile* outFile,
	       const std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > >& tridents);

/// Write out the monitor hists
void WriteHists(TFile* outFile,
		const TObjArray* monitorHists);

/// Run main function
void trident_mumu(int n_events = -1,
		  unsigned int nskip = 0,
		  bool run_signal = true,
		  int run_background = -1,
		  std::string inFilePath = "/pnfs/dune/persistent/users/jmalbos/Trident/data/sim/mumu/",
		  std::string outFileName = "TridentMuMuOut.root") {

  std::cout << std::endl
	    << "--------------------------------------------------------------------------------------" << std::endl
	    << "Running trident_mumu:" << std::endl
	    << "  Number of events: " << n_events << std::endl
	    << "  Number of events to skip: " << nskip << std::endl
	    << "  Running signal? " << run_signal << std::endl
	    << "  Running background files (-1: all): " << run_background << std::endl
	    << "  Input file path: " << inFilePath << std::endl
	    << "  Output file name: " << outFileName << std::endl
	    << "--------------------------------------------------------------------------------------" << std::endl
	    << std::endl;

  gStyle->SetOptStat(0);

  // gROOT->ProcessLine(".L Particle.cxx+");
  // gSystem->AddIncludePath("dunend/src/");
  // gSystem->Load("dunend/build/src/libEvent");
  // gSystem->Load("Particle_cxx");

  // Load file and tree
  TFile* file_sig = nullptr;
  TTree* tree_sig = nullptr;
  TChain* tree_bg = nullptr;
  if (run_signal) {
    file_sig = new TFile(Form("%s/mumu.sgn.coh.00.g4.root", inFilePath.c_str()), "READ");
    tree_sig = (TTree*)file_sig->Get("Event");
  }
  if (run_background > -2) {
    tree_bg = new TChain("Event");
    for (unsigned int bg_file = 0; bg_file < 100; ++bg_file)
      if (run_background == -1 or (run_background > -1 and bg_file == (unsigned int)run_background))
	tree_bg->Add(Form("%s/mumu.bkg.%02d.g4.root", inFilePath.c_str(), bg_file));
  }

  // Map of trident objects to fill
  std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > > allTridents;

  // Any histograms to fill
  TObjArray* monitorHists = nullptr;
  if (kMakeMonitorHists)
    monitorHists = DefineDataProducts();

  // Run the module over signal and background events
  if (tree_sig != nullptr) {
    unsigned int sig_events_to_process = n_events == -1 ? tree_sig->GetEntries()
      : TMath::Min((unsigned int)nskip+n_events, (unsigned int)tree_sig->GetEntries());
    allTridents["Signal"] = std::move(ProcessEvents(tree_sig, sig_events_to_process, nskip, "Signal", monitorHists));
  }
  if (tree_bg != nullptr) {
    unsigned int bg_events_to_process = n_events == -1 ? tree_bg->GetEntries()
      : TMath::Min((unsigned int)nskip+n_events, (unsigned int)tree_bg->GetEntries());
    allTridents["Background"] = std::move(ProcessEvents(tree_bg,  bg_events_to_process, nskip, "Background", monitorHists));
  }

  // Open output file
  TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");

  // Save a tree with info if required
  if (kSaveTree)
    WriteData(outFile, allTridents);
  if (kMakeMonitorHists)
    WriteHists(outFile, monitorHists);

  // Close signal file
  if (file_sig != nullptr) {
    if (file_sig->IsOpen())
      file_sig->Close();
    delete file_sig;
  }

  // Close output file
  outFile->Close();
  delete outFile;

  return;

}

TObjArray* DefineDataProducts() {

  TObjArray* objs = new TObjArray();
  TH1D* hTrackLengthSig		= new TH1D("TrackLengthSig",";Track length (mm);",100,0,5000);
  TH1D* hTrackEnergySig		= new TH1D("TrackEnergySig",";Track energy (MeV);",100,0,10000);
  TH1D* hLongTrackLengthSig	= new TH1D("LongTrackLengthSig",";Track length (mm);",100,0,5000);
  TH1D* hLongTrackEnergySig	= new TH1D("LongTrackEnergySig",";Track energy (MeV);",100,0,10000);
  TH1D* hShortTrackLengthSig	= new TH1D("ShortTrackLengthSig",";Track length (mm);",100,0,5000);
  TH1D* hShortTrackEnergySig	= new TH1D("ShortTrackEnergySig",";Track energy (MeV);",100,0,10000);
  TH1D* hTrackLengthBg		= new TH1D("TrackLengthBg",";Track length (mm);",100,0,5000);
  TH1D* hTrackEnergyBg		= new TH1D("TrackEnergyBg",";Track energy (MeV);",100,0,10000);
  TH1D* hLongTrackLengthBg	= new TH1D("LongTrackLengthBg",";Track length (mm);",100,0,5000);
  TH1D* hLongTrackEnergyBg	= new TH1D("LongTrackEnergyBg",";Track energy (MeV);",100,0,10000);
  TH1D* hShortTrackLengthBg	= new TH1D("ShortTrackLengthBg",";Track length (mm);",100,0,5000);
  TH1D* hShortTrackEnergyBg	= new TH1D("ShortTrackEnergyBg",";Track energy (MeV);",100,0,10000);
  TH2D* hVertexActivitySig	= new TH2D("VertexActivitySig",";Distance from vertex (mm);Deposited energy (MeV);",40,0,200,40,0,400);
  TH2D* hVertexActivityBg       = new TH2D("VertexActivityBg",";Distance from vertex (mm);Deposited energy (MeV);",40,0,200,40,0,400);
  TH1D* hVertexEnergySig     = new TH1D("VertexEnergySig",";Distance from vertex (mm);Deposited energy (MeV);",40,0,200);
  TH1D* hVertexEnergyBg      = new TH1D("VertexEnergyBg",";Distance from vertex (mm);Deposited energy (MeV);",40,0,200);
  TProfile* hVertexEnergyFracSig = new TProfile("VertexEnergyFracSig",";Distance from vertex (mm);Fraction of event energy",40,0,200,0,1);
  TProfile* hVertexEnergyFracBg  = new TProfile("VertexEnergyFracBg",";Distance from vertex (mm);Fraction of event energy",40,0,200,0,1);

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
  objs->Add(hVertexEnergySig);
  objs->Add(hVertexEnergyBg);
  objs->Add(hVertexEnergyFracSig);
  objs->Add(hVertexEnergyFracBg);

  return objs;

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

  const char* fSigLabel = (label == "Signal") ? "Sig" : "Bg";

  for (unsigned int event = nskip; event < n_events; ++event) {

    if (event % 100 == 0)
      std::cout << "\r" << label << ": processing event " << event << "/" << n_events
		<< " (" << (event)*100/(double)n_events << "%)" << std::flush;

    tree->GetEntry(event);

    // Save event information
    std::map<int,Particle*> particleMap;
    std::vector<MCHit> hits;
    std::map<int, int> particleAncestry;
    double eventEnergy = 0.;
    TVector3 eventVertex = evt->GetVertexPositionAndTime().Vect();

    // Look through MCParticles
    for (std::vector<MCParticle>::iterator particleIt = evt->MCParticleContainer.begin();
	 particleIt != evt->MCParticleContainer.end();
	 ++particleIt) {
      //particleMap[particleIt->GetMCID()] = Particle(particleIt->GetMCID());
      Particle* particle = new Particle(particleIt->GetMCID());
      particle->SetParticleStart(particleIt->GetInitialPositionAndTime().Vect());
      particle->SetPDG(particleIt->GetPDGCode());
      particle->SetInitialMomentum(particleIt->GetInitialMomentum());
      particle->SetFinalMomentum(particleIt->GetFinalMomentum());
      particle->SetDirection(particleIt->GetInitialMomentum().Unit());
      particle->SetAssociatedTrack(false);
      particleMap[particleIt->GetMCID()] = particle;
      particleAncestry[particleIt->GetMCID()] = particleIt->GetAncestorID();
    }

    // Look through MCTracks
    double trackEnergy;
    std::map<int, double> familyEnergy;
    for (std::vector<MCTrack>::iterator trackIt = evt->MCTrackContainer.begin();
	 trackIt != evt->MCTrackContainer.end();
	 ++trackIt) {
      trackEnergy = TrackDepositedEnergy(*trackIt);
      particleMap[trackIt->GetMCID()]->SetAssociatedTrack(true);
      particleMap[trackIt->GetMCID()]->SetDepositedEnergy(trackEnergy);
      particleMap[trackIt->GetMCID()]->SetdEdx(TrackdEdx(*trackIt, 30.));
      particleMap[trackIt->GetMCID()]->SetLength(TrackLength(*trackIt));
      //particleMap[trackIt->GetMCID()]->SetDirection(TrackDirection(*trackIt));
      particleMap[trackIt->GetMCID()]->SetStart(TrackStart(*trackIt));
      particleMap[trackIt->GetMCID()]->SetEnd(TrackEnd(*trackIt));
      particleMap[trackIt->GetMCID()]->SetVertexFraction(TrackVertexFraction(*trackIt, eventVertex));
      eventEnergy += trackEnergy;
      familyEnergy[particleAncestry[trackIt->GetMCID()]] += trackEnergy;
      std::vector<MCHit> trackHits = trackIt->MCHitContainer;
      for (std::vector<MCHit>::const_iterator trackHitIt = trackHits.begin();
      	   trackHitIt != trackHits.end();
      	   ++trackHitIt)
      	hits.push_back(*trackHitIt);
    }

    // Set family energy
    for (std::map<int, Particle*>::iterator particleIt = particleMap.begin(); particleIt != particleMap.end(); ++particleIt)
      particleIt->second->SetFamilyEnergy(familyEnergy[particleIt->second->ID()]);

    // Save all information about the trident candidate event
    TridentMuMu tridentMuMu;
    tridentMuMu.SetEventNumber(event);
    tridentMuMu.SetEventEnergy(eventEnergy);
    tridentMuMu.SetEventVertex(eventVertex);

    // Look through hits
    double vertexActivity = 0., vertexFraction = 0.;
    std::vector<double> vertexDepEnergy(40, 0.);
    TVector3 pos;
    double hitEnergy, vertexDist;
    for (std::vector<MCHit>::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
      pos = hitIt->GetStartPositionAndTime().Vect();
      hitEnergy = hitIt->GetEnergyDeposit();
      vertexDist = (pos-eventVertex).Mag();
      if (vertexDist < kVertexDifference)
	vertexActivity += hitEnergy;
      if (vertexDist < 200.)
	vertexDepEnergy[vertexDist/5] += hitEnergy;
      if (kMakeMonitorHists) {
	((TH1D*)monitorHists->FindObject(Form("VertexEnergy%s", fSigLabel)))->Fill(vertexDist, hitEnergy);
	((TH2D*)monitorHists->FindObject(Form("VertexActivity%s", fSigLabel)))->Fill(vertexDist, hitEnergy);
      }
    }
    if (kMakeMonitorHists)
      for (unsigned int bin = 0; bin < 40; ++bin)
	((TProfile*)monitorHists->FindObject(Form("VertexEnergyFrac%s", fSigLabel)))->Fill(bin*5., vertexDepEnergy[bin]/eventEnergy);
    vertexFraction = vertexActivity / eventEnergy;
    tridentMuMu.SetVertexActivity(vertexActivity);
    tridentMuMu.SetVertexFraction(vertexFraction);

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

      // Selection info
      //if (particleIt->second->AssociatedTrack() and particleIt->second->Length() > 500 and particleIt->second->DepositedEnergy() > 100)
      if (particleIt->second->AssociatedTrack() and particleIt->second->Length() > 20 and particleIt->second->FamilyEnergy() > 100)
	tridentMuMu.AddTrack(*(particleIt->second));

      particleLength[particleIt->second->Length()] = particleIt->second;

    }

    if (kMakeMonitorHists and particleLength.size() >= 2) {
      Particle* longestParticle = particleLength.rbegin()->second;
      Particle* nextLongestParticle = std::next(particleLength.rbegin())->second;
      ((TH1D*)monitorHists->FindObject(Form("TrackLength%s", fSigLabel)))->Fill(longestParticle->Length());
      ((TH1D*)monitorHists->FindObject(Form("TrackLength%s", fSigLabel)))->Fill(nextLongestParticle->Length());
      ((TH1D*)monitorHists->FindObject(Form("TrackEnergy%s", fSigLabel)))->Fill(longestParticle->DepositedEnergy());
      ((TH1D*)monitorHists->FindObject(Form("TrackEnergy%s", fSigLabel)))->Fill(nextLongestParticle->DepositedEnergy());
      ((TH1D*)monitorHists->FindObject(Form("LongTrackLength%s", fSigLabel)))->Fill(longestParticle->Length());
      ((TH1D*)monitorHists->FindObject(Form("LongTrackEnergy%s", fSigLabel)))->Fill(longestParticle->DepositedEnergy());
      ((TH1D*)monitorHists->FindObject(Form("ShortTrackLength%s", fSigLabel)))->Fill(nextLongestParticle->Length());
      ((TH1D*)monitorHists->FindObject(Form("ShortTrackEnergy%s", fSigLabel)))->Fill(nextLongestParticle->DepositedEnergy());
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
  double trackLength_simple = (firstHit.GetStartPositionAndTime().Vect() - lastHit.GetStopPositionAndTime().Vect()).Mag();

  double trackLength = 0.;
  for (std::vector<MCHit>::const_iterator hitIt = trackHits.begin(); hitIt != trackHits.end(); ++hitIt) {
    if (hitIt != trackHits.begin())
      trackLength += (hitIt->GetStartPositionAndTime().Vect() - std::next(hitIt, -1)->GetStopPositionAndTime().Vect()).Mag();
    trackLength += (hitIt->GetStopPositionAndTime().Vect() - hitIt->GetStartPositionAndTime().Vect()).Mag();
  }

  return trackLength;

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

double TrackVertexFraction(const MCTrack& track, const TVector3& vertex) {

  std::vector<MCHit> trackHits = track.MCHitContainer;

  double vertexEnergy = 0.;
  for (std::vector<MCHit>::const_iterator trackHitIt = trackHits.begin(); trackHitIt != trackHits.end(); ++trackHitIt)
    if ((trackHitIt->GetStartPositionAndTime().Vect()-vertex).Mag() < kVertexDifference)
      vertexEnergy += trackHitIt->GetEnergyDeposit();

  return vertexEnergy/TrackDepositedEnergy(track);

}

void WriteData(TFile* outFile,
	       const std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > >& tridents) {

  for (std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > >::const_iterator
	 allTridentIt = tridents.begin();
       allTridentIt != tridents.end();
       ++allTridentIt) {

    TTree* outTree = new TTree(("TridentMuMu"+allTridentIt->first).c_str(), "");//,"TridentMuMu"+tridentIt->first);

    // Tree variables
    int tEvent;
    TridentMuMu tTridentMuMu;
    double tEventEnergy;
    int tNumTracks;
    int tNumParticles;
    double tVertexActivity;
    double tVertexFraction;
    // Nu
    int tNuID, tNuPDG;
    double tNuStartX, tNuStartY, tNuStartZ;
    double tNuDirectionX, tNuDirectionY, tNuDirectionZ;
    // Particles
    int tID[kMaxParticles], tPDG[kMaxParticles];
    double tEnergy[kMaxParticles], tdEdx[kMaxParticles],
      tLength[kMaxParticles], tVertexFrac[kMaxParticles],
      tDirectionX[kMaxParticles], tDirectionY[kMaxParticles], tDirectionZ[kMaxParticles],
      tStartX[kMaxParticles], tStartY[kMaxParticles], tStartZ[kMaxParticles],
      tEndX[kMaxParticles], tEndY[kMaxParticles], tEndZ[kMaxParticles],
      tInitMomX[kMaxParticles], tInitMomY[kMaxParticles], tInitMomZ[kMaxParticles],
      tFinalMomX[kMaxParticles], tFinalMomY[kMaxParticles], tFinalMomZ[kMaxParticles];
    // Mu plus (true)
    int tMuPlusID, tMuPlusPDG;
    double tMuPlusEnergy, tMuPlusdEdx,
      tMuPlusFamilyEnergy,
      tMuPlusLength, tMuPlusVertexFrac,
      tMuPlusDirectionX, tMuPlusDirectionY, tMuPlusDirectionZ,
      tMuPlusStartX, tMuPlusStartY, tMuPlusStartZ,
      tMuPlusEndX, tMuPlusEndY, tMuPlusEndZ,
      tMuPlusInitMomX, tMuPlusInitMomY, tMuPlusInitMomZ,
      tMuPlusFinalMomX, tMuPlusFinalMomY, tMuPlusFinalMomZ;
    // Mu minus (true)
    int tMuMinusID, tMuMinusPDG;
    double tMuMinusEnergy, tMuMinusdEdx,
      tMuMinusFamilyEnergy,
      tMuMinusLength, tMuMinusVertexFrac,
      tMuMinusDirectionX, tMuMinusDirectionY, tMuMinusDirectionZ,
      tMuMinusStartX, tMuMinusStartY, tMuMinusStartZ,
      tMuMinusEndX, tMuMinusEndY, tMuMinusEndZ,
      tMuMinusInitMomX, tMuMinusInitMomY, tMuMinusInitMomZ,
      tMuMinusFinalMomX, tMuMinusFinalMomY, tMuMinusFinalMomZ;
    // Short muon (selected)
    int tShortMuonID, tShortMuonPDG;
    double tShortMuonEnergy, tShortMuondEdx,
      tShortMuonFamilyEnergy,
      tShortMuonLength, tShortMuonVertexFrac,
      tShortMuonDirectionX, tShortMuonDirectionY, tShortMuonDirectionZ,
      tShortMuonStartX, tShortMuonStartY, tShortMuonStartZ,
      tShortMuonEndX, tShortMuonEndY, tShortMuonEndZ,
      tShortMuonInitMomX, tShortMuonInitMomY, tShortMuonInitMomZ,
      tShortMuonFinalMomX, tShortMuonFinalMomY, tShortMuonFinalMomZ;
    // Long muon (selected)
    int tLongMuonID, tLongMuonPDG;
    double tLongMuonEnergy, tLongMuondEdx,
      tLongMuonFamilyEnergy,
      tLongMuonLength, tLongMuonVertexFrac,
      tLongMuonDirectionX, tLongMuonDirectionY, tLongMuonDirectionZ,
      tLongMuonStartX, tLongMuonStartY, tLongMuonStartZ,
      tLongMuonEndX, tLongMuonEndY, tLongMuonEndZ,
      tLongMuonInitMomX, tLongMuonInitMomY, tLongMuonInitMomZ,
      tLongMuonFinalMomX, tLongMuonFinalMomY, tLongMuonFinalMomZ;
    // Selected muons
    double tAngle, tFrontSeparation, tBackSeparation, tAngle_s;

    // Set branch vars
    // Event
    //outTree->Branch("TridentMuMu",      tTridentMuMu);
    outTree->Branch("Event",          &tEvent);
    outTree->Branch("EventEnergy",    &tEventEnergy);
    outTree->Branch("NumTracks",      &tNumTracks);
    outTree->Branch("NumParticles",   &tNumParticles);
    outTree->Branch("VertexActivity", &tVertexActivity);
    outTree->Branch("VertexFraction", &tVertexFraction);
    // Nu
    outTree->Branch("NuID",         &tNuID);
    outTree->Branch("NuPDG",        &tNuPDG);
    outTree->Branch("NuStartX",     &tNuStartX);
    outTree->Branch("NuStartY",     &tNuStartY);
    outTree->Branch("NuStartZ",     &tNuStartZ);
    outTree->Branch("NuDirectionX", &tNuDirectionX);
    outTree->Branch("NuDirectionY", &tNuDirectionY);
    outTree->Branch("NuDirectionZ", &tNuDirectionZ);
    // All particles
    outTree->Branch("ID",         &tID,         "ID[NumParticles]/I");
    outTree->Branch("PDG",        &tPDG,        "PDG[NumParticles]/I");
    outTree->Branch("Energy",     &tEnergy,     "PDG[NumParticles]/F");
    outTree->Branch("dEdx",       &tdEdx,       "dEdx[NumParticles]/F");
    outTree->Branch("Length",     &tLength,     "Length[NumParticles]/F");
    outTree->Branch("VertexFrac", &tVertexFrac, "VertexFrac[NumParticles]/F");
    outTree->Branch("DirectionX", &tDirectionX, "MuPlusDirectionX[NumParticles]/F");
    outTree->Branch("DirectionY", &tDirectionY, "MuPlusDirectionY[NumParticles]/F");
    outTree->Branch("DirectionZ", &tDirectionZ, "MuPlusDirectionZ[NumParticles]/F");
    outTree->Branch("StartX",     &tStartX,     "MuPlusStartX[NumParticles]/F");
    outTree->Branch("StartY",     &tStartY,     "MuPlusStartY[NumParticles]/F");
    outTree->Branch("StartZ",     &tStartZ,     "MuPlusStartZ[NumParticles]/F");
    outTree->Branch("EndX",       &tEndX,       "MuPlusEndX[NumParticles]/F");
    outTree->Branch("EndY",       &tEndY,       "MuPlusEndY[NumParticles]/F");
    outTree->Branch("EndZ",       &tEndZ,       "MuPlusEndZ[NumParticles]/F");
    outTree->Branch("InitMomX",   &tInitMomX,   "MuPlusInitMomX[NumParticles]/F");
    outTree->Branch("InitMomY",   &tInitMomY,   "MuPlusInitMomY[NumParticles]/F");
    outTree->Branch("InitMomZ",   &tInitMomZ,   "MuPlusInitMomZ[NumParticles]/F");
    outTree->Branch("FinalMomX",  &tFinalMomX,  "MuPlusFinalMomX[NumParticles]/F");
    outTree->Branch("FinalMomY",  &tFinalMomY,  "MuPlusFinalMomY[NumParticles]/F");
    outTree->Branch("FinalMomZ",  &tFinalMomZ,  "MuPlusFinalMomZ[NumParticles]/F");
    // Mu plus (true)
    outTree->Branch("MuPlusID",         &tMuPlusID);
    outTree->Branch("MuPlusPDG",        &tMuPlusPDG);
    outTree->Branch("MuPlusEnergy",     &tMuPlusEnergy);
    outTree->Branch("MuPlusdEdx",       &tMuPlusdEdx);
    outTree->Branch("MuPlusFamilyEnergy", &tMuPlusFamilyEnergy);
    outTree->Branch("MuPlusLength",     &tMuPlusLength);
    outTree->Branch("MuPlusVertexFrac", &tMuPlusVertexFrac);
    outTree->Branch("MuPlusDirectionX", &tMuPlusDirectionX);
    outTree->Branch("MuPlusDirectionY", &tMuPlusDirectionY);
    outTree->Branch("MuPlusDirectionZ", &tMuPlusDirectionZ);
    outTree->Branch("MuPlusStartX",     &tMuPlusStartX);
    outTree->Branch("MuPlusStartY",     &tMuPlusStartY);
    outTree->Branch("MuPlusStartZ",     &tMuPlusStartZ);
    outTree->Branch("MuPlusEndX",       &tMuPlusEndX);
    outTree->Branch("MuPlusEndY",       &tMuPlusEndY);
    outTree->Branch("MuPlusEndZ",       &tMuPlusEndZ);
    outTree->Branch("MuPlusInitMomX",   &tMuPlusInitMomX);
    outTree->Branch("MuPlusInitMomY",   &tMuPlusInitMomY);
    outTree->Branch("MuPlusInitMomZ",   &tMuPlusInitMomZ);
    outTree->Branch("MuPlusFinalMomX",  &tMuPlusFinalMomX);
    outTree->Branch("MuPlusFinalMomY",  &tMuPlusFinalMomY);
    outTree->Branch("MuPlusFinalMomZ",  &tMuPlusFinalMomZ);
    // Mu minus (true)
    outTree->Branch("MuMinusID",         &tMuMinusID);
    outTree->Branch("MuMinusPDG",        &tMuMinusPDG);
    outTree->Branch("MuMinusEnergy",     &tMuMinusEnergy);
    outTree->Branch("MuMinusdEdx",       &tMuMinusdEdx);
    outTree->Branch("MuMinusFamilyEnergy", &tMuMinusFamilyEnergy);
    outTree->Branch("MuMinusLength",     &tMuMinusLength);
    outTree->Branch("MuMinusVertexFrac", &tMuMinusVertexFrac);
    outTree->Branch("MuMinusDirectionX", &tMuMinusDirectionX);
    outTree->Branch("MuMinusDirectionY", &tMuMinusDirectionY);
    outTree->Branch("MuMinusDirectionZ", &tMuMinusDirectionZ);
    outTree->Branch("MuMinusStartX",     &tMuMinusStartX);
    outTree->Branch("MuMinusStartY",     &tMuMinusStartY);
    outTree->Branch("MuMinusStartZ",     &tMuMinusStartZ);
    outTree->Branch("MuMinusEndX",       &tMuMinusEndX);
    outTree->Branch("MuMinusEndY",       &tMuMinusEndY);
    outTree->Branch("MuMinusEndZ",       &tMuMinusEndZ);
    outTree->Branch("MuMinusInitMomX",   &tMuMinusInitMomX);
    outTree->Branch("MuMinusInitMomY",   &tMuMinusInitMomY);
    outTree->Branch("MuMinusInitMomZ",   &tMuMinusInitMomZ);
    outTree->Branch("MuMinusFinalMomX",  &tMuMinusFinalMomX);
    outTree->Branch("MuMinusFinalMomY",  &tMuMinusFinalMomY);
    outTree->Branch("MuMinusFinalMomZ",  &tMuMinusFinalMomZ);
    // Short muon (selected)
    outTree->Branch("ShortMuonID",         &tShortMuonID);
    outTree->Branch("ShortMuonPDG",        &tShortMuonPDG);
    outTree->Branch("ShortMuonEnergy",     &tShortMuonEnergy);
    outTree->Branch("ShortMuondEdx",       &tShortMuondEdx);
    outTree->Branch("ShortMuonFamilyEnergy", &tShortMuonFamilyEnergy);
    outTree->Branch("ShortMuonLength",     &tShortMuonLength);
    outTree->Branch("ShortMuonVertexFrac", &tShortMuonVertexFrac);
    outTree->Branch("ShortMuonDirectionX", &tShortMuonDirectionX);
    outTree->Branch("ShortMuonDirectionY", &tShortMuonDirectionY);
    outTree->Branch("ShortMuonDirectionZ", &tShortMuonDirectionZ);
    outTree->Branch("ShortMuonStartX",     &tShortMuonStartX);
    outTree->Branch("ShortMuonStartY",     &tShortMuonStartY);
    outTree->Branch("ShortMuonStartZ",     &tShortMuonStartZ);
    outTree->Branch("ShortMuonEndX",       &tShortMuonEndX);
    outTree->Branch("ShortMuonEndY",       &tShortMuonEndY);
    outTree->Branch("ShortMuonEndZ",       &tShortMuonEndZ);
    outTree->Branch("ShortMuonInitMomX",   &tShortMuonInitMomX);
    outTree->Branch("ShortMuonInitMomY",   &tShortMuonInitMomY);
    outTree->Branch("ShortMuonInitMomZ",   &tShortMuonInitMomZ);
    outTree->Branch("ShortMuonFinalMomX",  &tShortMuonFinalMomX);
    outTree->Branch("ShortMuonFinalMomY",  &tShortMuonFinalMomY);
    outTree->Branch("ShortMuonFinalMomZ",  &tShortMuonFinalMomZ);
    // Long muon (selected)
    outTree->Branch("LongMuonID",         &tLongMuonID);
    outTree->Branch("LongMuonPDG",        &tLongMuonPDG);
    outTree->Branch("LongMuonEnergy",     &tLongMuonEnergy);
    outTree->Branch("LongMuondEdx",       &tLongMuondEdx);
    outTree->Branch("LongMuonFamilyEnergy", &tLongMuonFamilyEnergy);
    outTree->Branch("LongMuonLength",     &tLongMuonLength);
    outTree->Branch("LongMuonVertexFrac", &tLongMuonVertexFrac);
    outTree->Branch("LongMuonDirectionX", &tLongMuonDirectionX);
    outTree->Branch("LongMuonDirectionY", &tLongMuonDirectionY);
    outTree->Branch("LongMuonDirectionZ", &tLongMuonDirectionZ);
    outTree->Branch("LongMuonStartX",     &tLongMuonStartX);
    outTree->Branch("LongMuonStartY",     &tLongMuonStartY);
    outTree->Branch("LongMuonStartZ",     &tLongMuonStartZ);
    outTree->Branch("LongMuonEndX",       &tLongMuonEndX);
    outTree->Branch("LongMuonEndY",       &tLongMuonEndY);
    outTree->Branch("LongMuonEndZ",       &tLongMuonEndZ);
    outTree->Branch("LongMuonInitMomX",   &tLongMuonInitMomX);
    outTree->Branch("LongMuonInitMomY",   &tLongMuonInitMomY);
    outTree->Branch("LongMuonInitMomZ",   &tLongMuonInitMomZ);
    outTree->Branch("LongMuonFinalMomX",  &tLongMuonFinalMomX);
    outTree->Branch("LongMuonFinalMomY",  &tLongMuonFinalMomY);
    outTree->Branch("LongMuonFinalMomZ",  &tLongMuonFinalMomZ);
    // Selected muons
    outTree->Branch("Angle",           &tAngle);
    outTree->Branch("Angle_s",         &tAngle_s);
    outTree->Branch("FrontSeparation", &tFrontSeparation);
    outTree->Branch("BackSeparation",  &tBackSeparation);

    for (std::vector<std::unique_ptr<TridentMuMu> >::const_iterator
	   tridentIt = allTridentIt->second.begin();
	 tridentIt != allTridentIt->second.end();
	 ++tridentIt) {

      //tTridentMuMu = *(*tridentIt).release();
      tEvent = (*tridentIt)->GetEventNumber();
      tEventEnergy = (*tridentIt)->GetEventEnergy();
      tNumTracks = (*tridentIt)->NumTracks();
      tVertexActivity = (*tridentIt)->GetVertexActivity();
      tVertexFraction = (*tridentIt)->GetVertexFraction();

      // Get particle objects
      Particle nu = (*tridentIt)->GetNu();
      std::vector<Particle> mus = (*tridentIt)->GetTracks();
      Particle muPlus = (*tridentIt)->GetMuPlus();
      Particle muMinus = (*tridentIt)->GetMuMinus();

      Particle longMuon, shortMuon;
      if (tNumTracks > 1) {
	// In the instance of multiple tracks, need to pick two
	// Right now, take the longest track
	std::map<double,Particle> muLength;
	for (std::vector<Particle>::iterator muIt = mus.begin(); muIt != mus.end(); ++muIt)
	  muLength[muIt->Length()] = *muIt;
	std::map<double,Particle>::reverse_iterator muIt = muLength.rbegin();
	longMuon = muIt->second;
	shortMuon = std::next(muIt)->second;
      }
      else if (tNumTracks == 1)
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
      tMuPlusEnergy     = muPlus.DepositedEnergy();
      tMuPlusdEdx       = muPlus.dEdx();
      tMuPlusFamilyEnergy = muPlus.FamilyEnergy();
      tMuPlusLength     = muPlus.Length();
      tMuPlusVertexFrac = muPlus.VertexFraction();
      tMuPlusDirectionX = muPlus.Direction().X();
      tMuPlusDirectionY = muPlus.Direction().Y();
      tMuPlusDirectionZ = muPlus.Direction().Z();
      tMuPlusStartX     = muPlus.Start().X();
      tMuPlusStartY     = muPlus.Start().Y();
      tMuPlusStartZ     = muPlus.Start().Z();
      tMuPlusEndX       = muPlus.End().X();
      tMuPlusEndY       = muPlus.End().Y();
      tMuPlusEndZ       = muPlus.End().Z();
      tMuPlusInitMomX   = muPlus.InitialMomentum().X();
      tMuPlusInitMomY   = muPlus.InitialMomentum().Y();
      tMuPlusInitMomZ   = muPlus.InitialMomentum().Z();
      tMuPlusFinalMomX  = muPlus.FinalMomentum().X();
      tMuPlusFinalMomY  = muPlus.FinalMomentum().Y();
      tMuPlusFinalMomZ  = muPlus.FinalMomentum().Z();

      tMuMinusID         = muMinus.ID();
      tMuMinusPDG        = muMinus.PDG();
      tMuMinusEnergy     = muMinus.DepositedEnergy();
      tMuMinusdEdx       = muMinus.dEdx();
      tMuMinusFamilyEnergy = muMinus.FamilyEnergy();
      tMuMinusLength     = muMinus.Length();
      tMuMinusVertexFrac = muMinus.VertexFraction();
      tMuMinusDirectionX = muMinus.Direction().X();
      tMuMinusDirectionY = muMinus.Direction().Y();
      tMuMinusDirectionZ = muMinus.Direction().Z();
      tMuMinusStartX     = muMinus.Start().X();
      tMuMinusStartY     = muMinus.Start().Y();
      tMuMinusStartZ     = muMinus.Start().Z();
      tMuMinusEndX       = muMinus.End().X();
      tMuMinusEndY       = muMinus.End().Y();
      tMuMinusEndZ       = muMinus.End().Z();
      tMuMinusInitMomX   = muMinus.InitialMomentum().X();
      tMuMinusInitMomY   = muMinus.InitialMomentum().Y();
      tMuMinusInitMomZ   = muMinus.InitialMomentum().Z();
      tMuMinusFinalMomX  = muMinus.FinalMomentum().X();
      tMuMinusFinalMomY  = muMinus.FinalMomentum().Y();
      tMuMinusFinalMomZ  = muMinus.FinalMomentum().Z();

      tLongMuonID         = longMuon.ID();
      tLongMuonPDG        = longMuon.PDG();
      tLongMuonEnergy     = longMuon.DepositedEnergy();
      tLongMuondEdx       = longMuon.dEdx();
      tLongMuonFamilyEnergy = longMuon.FamilyEnergy();
      tLongMuonLength     = longMuon.Length();
      tLongMuonVertexFrac = longMuon.VertexFraction();
      tLongMuonDirectionX = longMuon.Direction().X();
      tLongMuonDirectionY = longMuon.Direction().Y();
      tLongMuonDirectionZ = longMuon.Direction().Z();
      tLongMuonStartX     = longMuon.Start().X();
      tLongMuonStartY     = longMuon.Start().Y();
      tLongMuonStartZ     = longMuon.Start().Z();
      tLongMuonEndX       = longMuon.End().X();
      tLongMuonEndY       = longMuon.End().Y();
      tLongMuonEndZ       = longMuon.End().Z();
      tLongMuonInitMomX   = longMuon.InitialMomentum().X();
      tLongMuonInitMomY   = longMuon.InitialMomentum().Y();
      tLongMuonInitMomZ   = longMuon.InitialMomentum().Z();
      tLongMuonFinalMomX  = longMuon.FinalMomentum().X();
      tLongMuonFinalMomY  = longMuon.FinalMomentum().Y();
      tLongMuonFinalMomZ  = longMuon.FinalMomentum().Z();

      tShortMuonID         = shortMuon.ID();
      tShortMuonPDG        = shortMuon.PDG();
      tShortMuonEnergy     = shortMuon.DepositedEnergy();
      tShortMuondEdx       = shortMuon.dEdx();
      tShortMuonFamilyEnergy = shortMuon.FamilyEnergy();
      tShortMuonLength     = shortMuon.Length();
      tShortMuonVertexFrac = shortMuon.VertexFraction();
      tShortMuonDirectionX = shortMuon.Direction().X();
      tShortMuonDirectionY = shortMuon.Direction().Y();
      tShortMuonDirectionZ = shortMuon.Direction().Z();
      tShortMuonStartX     = shortMuon.Start().X();
      tShortMuonStartY     = shortMuon.Start().Y();
      tShortMuonStartZ     = shortMuon.Start().Z();
      tShortMuonEndX       = shortMuon.End().X();
      tShortMuonEndY       = shortMuon.End().Y();
      tShortMuonEndZ       = shortMuon.End().Z();
      tShortMuonInitMomX   = shortMuon.InitialMomentum().X();
      tShortMuonInitMomY   = shortMuon.InitialMomentum().Y();
      tShortMuonInitMomZ   = shortMuon.InitialMomentum().Z();
      tShortMuonFinalMomX  = shortMuon.FinalMomentum().X();
      tShortMuonFinalMomY  = shortMuon.FinalMomentum().Y();
      tShortMuonFinalMomZ  = shortMuon.FinalMomentum().Z();

      tAngle = longMuon.Direction().
	Angle(shortMuon.Direction()) * 180/TMath::Pi();
      tAngle_s = (longMuon.End()-longMuon.Start()).Unit().
	Angle((shortMuon.End()-shortMuon.Start()).Unit()) * 180/TMath::Pi();
      tFrontSeparation = (longMuon.Start()-
			  shortMuon.Start()).Mag();
      tBackSeparation = (longMuon.End()-
			 shortMuon.End()).Mag();

      outTree->Fill();

    } // tree

    outFile->cd();
    outTree->Write();
    delete outTree;

  } // signal/bg

  return;

}

void WriteHists(TFile* outFile,
		const TObjArray* monitorHists) {

  // Write out all the monitor hists
  if (monitorHists) {
    outFile->cd();
    for (int obj = 0; obj < monitorHists->GetEntries(); ++obj)
      monitorHists->At(obj)->Write();
  }

  return;

}
