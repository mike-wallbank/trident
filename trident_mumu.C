// Load dynamic library containing the definitions of the event classes
#include "Rtypes.h"
R__LOAD_LIBRARY(dunend/build/src/libEvent);
R__LOAD_LIBRARY(Particle_cxx.so);

// #include "dunend/src/Event.h"
// #include "Particle.h"

// C++
#include <iostream>
#include <algorithm>
#include <memory>

// ROOT
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TEfficiency.h"

// framework
#include "Particle.h"
#include "dunend/src/MCHit.h"
#include "dunend/src/MCTrack.h"
#include "dunend/src/MCParticle.h"
#include "dunend/src/Event.h"

const int kMaxParticles = 1000;

class TridentMuMu {
public:

  TridentMuMu();

  void AddNu(Particle nu);
  void AddMuPlus(Particle mu);
  void AddMuMinus(Particle mu);
  void AddIdentifiedMu(Particle mu);
  double GetEventEnergy();
  Particle GetNu();
  Particle GetMuPlus();
  Particle GetMuMinus();
  Particle GetMu(int mu);
  std::vector<Particle> GetMus();
  double MuonAngle();
  double MuonSep();
  double MuonTrackSep();
  int NumMuPlus();
  int NumMuMinus();
  int NumIdentifiedMus();
  void SetEventEnergy(double energy);

private:

  double fEventEnergy;
  Particle fNu;
  Particle fMuPlus;
  Particle fMuMinus;
  std::vector<Particle> fIdentifiedMus;
  std::map<double,Particle> fMuMap;

};

TridentMuMu::TridentMuMu() {
  fEventEnergy = 0;
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

std::vector<Particle> TridentMuMu::GetMus() {
  return fIdentifiedMus;
}

double TridentMuMu::MuonAngle() {
  return (fMuMap.rbegin()->second.Direction().Angle(std::next(fMuMap.rbegin())->second.Direction()))*180/TMath::Pi();
}

double TridentMuMu::MuonSep() {
  return (fMuMap.rbegin()->second.Start()-std::next(fMuMap.rbegin())->second.Start()).Mag();
}

double TridentMuMu::MuonTrackSep() {
  return (fMuMap.rbegin()->second.StartTrack()-std::next(fMuMap.rbegin())->second.StartTrack()).Mag();
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

std::map<int,std::pair<std::string,int> > PDGMap = {{13,{"#mu",0}},{211,{"#pi",1}},{11,{"e",2}},{2212,{"p",3}},{321,{"K",4}}};

/// Process the events in one of the data samples
std::vector<std::unique_ptr<TridentMuMu> > ProcessEvents(TTree* tree,
							 unsigned int n_events,
							 unsigned int nskip,
							 std::string label);

/// Returns the length of the track
double TrackLength(const MCTrack& track);

/// Returns the energy of the track
double TrackEnergy(const MCTrack& track);

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
	       const std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > >& tridents);

/// Run main function
void trident_mumu(int n_events = -1, unsigned int nskip = 0, std::string outFile = "TridentMuMuOut.root") {

  // gROOT->ProcessLine(".L Particle.cxx+");
  //gSystem->Load("dunend/build/src/libEvent.dylib");
  gStyle->SetOptStat(0);

  // Load file and tree
  TFile file_sig("ArSMmumuDec11.root");
  TFile file_bg("ArSMmumu_backg.root");
  TTree* tree_sig = (TTree*)file_sig.Get("Event");
  TTree* tree_bg = (TTree*)file_bg.Get("Event");

  unsigned int sig_events_to_process = n_events == -1 ? tree_sig->GetEntriesFast()
    : TMath::Min((unsigned int)nskip+n_events, (unsigned int)tree_sig->GetEntriesFast());
  unsigned int bg_events_to_process = n_events == -1 ? tree_bg->GetEntriesFast()
    : TMath::Min((unsigned int)nskip+n_events, (unsigned int)tree_bg->GetEntriesFast());

  std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > > allTridents;

  allTridents["Signal"]     = std::move(ProcessEvents(tree_sig, sig_events_to_process, nskip, "Signal"));
  allTridents["Background"] = std::move(ProcessEvents(tree_bg,  bg_events_to_process,  nskip, "Background"));

  WriteData(outFile, allTridents);

  file_sig.Close();
  file_bg.Close();
  
}

std::vector<std::unique_ptr<TridentMuMu> > ProcessEvents(TTree* tree,
							 unsigned int n_events,
							 unsigned int nskip,
							 std::string label) {

  // Sync tree with empty event
  Event* evt = 0;
  tree->SetBranchAddress("Event", &evt);

  std::vector<std::unique_ptr<TridentMuMu> > tridentMuMus;

  for (unsigned int event = nskip; event < n_events; ++event) {

    if (event % 1000 == 0)
      std::cout << "Processing event " << event << std::endl;

    tree->GetEntry(event);

    std::map<int,Particle*> particleMap;
    double eventEnergy = 0.;

    for (std::vector<MCParticle>::iterator particleIt = evt->MCParticleContainer.begin();
	 particleIt != evt->MCParticleContainer.end();
	 ++particleIt) {
      //particleMap[particleIt->GetMCID()] = Particle(particleIt->GetMCID());
      Particle* particle = new Particle(particleIt->GetMCID());
      particle->SetStart(particleIt->GetInitialPositionAndTime().Vect());
      particle->SetPDG(particleIt->GetPDGCode());
      particle->SetAssociatedTrack(false);
      particleMap[particleIt->GetMCID()] = particle;
    }
    for (std::vector<MCTrack>::iterator trackIt = evt->MCTrackContainer.begin();
	 trackIt != evt->MCTrackContainer.end();
	 ++trackIt) {
      particleMap[trackIt->GetMCID()]->SetAssociatedTrack(true);
      particleMap[trackIt->GetMCID()]->SetEnergy(TrackEnergy(*trackIt));
      particleMap[trackIt->GetMCID()]->SetdEdx(TrackdEdx(*trackIt, 100.));
      particleMap[trackIt->GetMCID()]->SetLength(TrackLength(*trackIt));
      particleMap[trackIt->GetMCID()]->SetDirection(TrackDirection(*trackIt));
      particleMap[trackIt->GetMCID()]->SetStartTrack(TrackStart(*trackIt));
      eventEnergy += TrackEnergy(*trackIt);
    }

    // Look for the components of the trident
    TridentMuMu tridentMuMu;
    tridentMuMu.SetEventEnergy(eventEnergy);

    // Look through particles
    int nParticlesLength = 0, nParticlesLength1 = 0;
    for (std::map<int,Particle*>::const_iterator particleIt = particleMap.begin();
	 particleIt != particleMap.end();
	 ++particleIt) {
      if (particleIt->second->Length())
	++nParticlesLength;
      if (particleIt->second->Length() > 1.)
	++nParticlesLength1;
      if (label == "Signal") {
	if (particleIt->second->PDG() == 14)
	  tridentMuMu.AddNu(*(particleIt->second));
	if (particleIt->second->PDG() == 13 and tridentMuMu.NumMuPlus() == 0)
	  tridentMuMu.AddMuPlus(*(particleIt->second));
	if (particleIt->second->PDG() == -13 and tridentMuMu.NumMuMinus() == 0)
	  tridentMuMu.AddMuMinus(*(particleIt->second));
      }
      if (particleIt->second->AssociatedTrack() and particleIt->second->Length() > 50 and particleIt->second->Energy() > 100)
	tridentMuMu.AddIdentifiedMu(*(particleIt->second));
    }

    tridentMuMus.push_back(std::make_unique<TridentMuMu>(tridentMuMu));

  } // events

  return tridentMuMus;

}

double TrackLength(const MCTrack& track) {

  std::vector<MCHit> trackHits = track.MCHitContainer;
  MCHit firstHit = *trackHits.begin();
  MCHit lastHit = *trackHits.rbegin();

  double trackLength = (firstHit.GetPositionAndTime().Vect() - lastHit.GetPositionAndTime().Vect()).Mag();

  return trackLength/10.;

}

double TrackEnergy(const MCTrack& track) {

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
  TVector3 lastHitPos = trackHits.begin()->GetPositionAndTime().Vect();
  int nhit = 0;
  for (std::vector<MCHit>::iterator hitIt = trackHits.begin(); hitIt != trackHits.end() and dx < dEdxlength; ++hitIt) {
    TVector3 hitPos = hitIt->GetPositionAndTime().Vect();
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

  return trackHits.begin()->GetPositionAndTime().Vect();

}

TVector3 TrackEnd(const MCTrack& track) {

  std::vector<MCHit> trackHits = track.MCHitContainer;

  return trackHits.rbegin()->GetPositionAndTime().Vect();

}

TVector3 TrackDirection(const MCTrack& track) {

  TVector3 direction;
  
  std::vector<MCHit> trackHits = track.MCHitContainer;
  
  if (trackHits.size() < 2)
    direction = TVector3(0,0,0);

  else if (trackHits.size() == 2) {
    TVector3 hit1 = trackHits[0].GetPositionAndTime().Vect();
    TVector3 hit2 = trackHits[1].GetPositionAndTime().Vect();
    direction = (hit2-hit1).Unit();
  }

  else if (trackHits.size() > 2) {
    TVector3 hit1 = trackHits[0].GetPositionAndTime().Vect();
    TVector3 hit2 = trackHits[trackHits.size()-1].GetPositionAndTime().Vect();
    direction = (hit2-hit1).Unit();
  }

  return direction;

}

void WriteData(std::string outFile,
	       const std::map<std::string,std::vector<std::unique_ptr<TridentMuMu> > >& tridents) {

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
    double tNuStart[3];
    // Particles
    int tID[kMaxParticles], tPDG[kMaxParticles];
    double tEnergy[kMaxParticles], tdEdx[kMaxParticles],
      tLength[kMaxParticles], tMomentum[kMaxParticles][3],
      tStart[kMaxParticles][3], tEnd[kMaxParticles][3];
    // Mu plus (true)
    int tMuPlusID, tMuPlusPDG;
    double tMuPlusEnergy, tMuPlusdEdx,
      tMuPlusLength, tMuPlusMomentum[3],
      tMuPlusStart[3], tMuPlusEnd[3];
    // Mu minus (true)
    int tMuMinusID, tMuMinusPDG;
    double tMuMinusEnergy, tMuMinusdEdx,
      tMuMinusLength, tMuMinusMomentum[3],
      tMuMinusStart[3], tMuMinusEnd[3];
    // Short muon (selected)
    int tShortMuonID, tShortMuonPDG;
    double tShortMuonEnergy, tShortMuondEdx,
      tShortMuonLength, tShortMuonMomentum[3],
      tShortMuonStart[3], tShortMuonEnd[3];
    // Long muon (selected)
    int tLongMuonID, tLongMuonPDG;
    double tLongMuonEnergy, tLongMuondEdx,
      tLongMuonLength, tLongMuonMomentum[3],
      tLongMuonStart[3], tLongMuonEnd[3];

    // Set branch vars
    // Event
    //outT->Branch("TridentMuMu",        tTridentMuMu);
    outT->Branch("EventEnergy",      &tEventEnergy);
    outT->Branch("NumIdentifiedMus", &tNumIdentifiedMus);
    outT->Branch("NumParticles",     &tNumParticles);
    // Nu
    outT->Branch("NuID",               &tNuID);
    outT->Branch("NuPDG",              &tNuPDG);
    outT->Branch("NuStart",            &tNuStart, "tNuStart[3]/F");
    // All particles
    outT->Branch("ID",           &tID,       "ID[NumParticles]/I");
    outT->Branch("PDG",          &tPDG,      "PDG[NumParticles]/I");
    outT->Branch("Energy",       &tEnergy,   "PDG[NumParticles]/F");
    outT->Branch("dEdx",         &tdEdx,     "dEdx[NumParticles]/F");
    outT->Branch("Length",       &tLength,   "Length[NumParticles]/F");
    outT->Branch("Momentum",     &tMomentum, "MuPlusMomentum[NumParticles][3]/F");
    outT->Branch("Start",        &tStart,    "MuPlusStart[NumParticles][3]/F");
    outT->Branch("End",          &tEnd,      "MuPlusEnd[NumParticles][3]/F");
    // Mu plus (true)
    outT->Branch("MuPlusID",           &tMuPlusID);
    outT->Branch("MuPlusPDG",          &tMuPlusPDG);
    outT->Branch("MuPlusEnergy",       &tMuPlusEnergy);
    outT->Branch("MuPlusdEdx",         &tMuPlusdEdx);
    outT->Branch("MuPlusLength",       &tMuPlusLength);
    outT->Branch("MuPlusMomentum",     &tMuPlusMomentum, "MuPlusMonentum[3]/F");
    outT->Branch("MuPlusStart",        &tMuPlusStart,    "MuPlusStart[3]/F");
    outT->Branch("MuPlusEnd",          &tMuPlusEnd,      "MuPlusEnd[3]/F");
    // Mu minus (true)
    outT->Branch("MuMinusID",          &tMuMinusID);
    outT->Branch("MuMinusPDG",         &tMuMinusPDG);
    outT->Branch("MuMinusEnergy",      &tMuMinusEnergy);
    outT->Branch("MuMinusdEdx",        &tMuMinusdEdx);
    outT->Branch("MuMinusLength",      &tMuMinusLength);
    outT->Branch("MuMinusMomentum",    &tMuMinusMomentum, "MuMinusMonentum[3]/F");
    outT->Branch("MuMinusStart",       &tMuMinusStart,    "MuMinusStart[3]/F");
    outT->Branch("MuMinusEnd",         &tMuMinusEnd,      "MuMinusEnd[3]/F");
    // Short muon (selected)
    outT->Branch("ShortMuonID",        &tShortMuonID);
    outT->Branch("ShortMuonPDG",       &tShortMuonPDG);
    outT->Branch("ShortMuonEnergy",    &tShortMuonEnergy);
    outT->Branch("ShortMuondEdx",      &tShortMuondEdx);
    outT->Branch("ShortMuonLength",    &tShortMuonLength);
    outT->Branch("ShortMuonMomentum",  &tShortMuonMomentum, "ShortMuonMonentum[3]/F");
    outT->Branch("ShortMuonStart",     &tShortMuonStart,    "ShortMuonStart[3]/F");
    outT->Branch("ShortMuonEnd",       &tShortMuonEnd,      "ShortMuonEnd[3]/F");
    // Long muon (selected)
    outT->Branch("LongMuonID",         &tLongMuonID);
    outT->Branch("LongMuonPDG",        &tLongMuonPDG);
    outT->Branch("LongMuonEnergy",     &tLongMuonEnergy);
    outT->Branch("LongMuondEdx",       &tLongMuondEdx);
    outT->Branch("LongMuonLength",     &tLongMuonLength);
    outT->Branch("LongMuonMomentum",   &tLongMuonMomentum, "LongMuonMonentum[3]/F");
    outT->Branch("LongMuonStart",      &tLongMuonStart,    "LongMuonStart[3]/F");
    outT->Branch("LongMuonEnd",        &tLongMuonEnd,      "LongMuonEnd[3]/F");

    for (std::vector<std::unique_ptr<TridentMuMu> >::const_iterator
	   tridentIt = allTridentIt->second.begin();
	 tridentIt != allTridentIt->second.end();
	 ++tridentIt) {
      
      //tTridentMuMu = *(*tridentIt).release();
      tEventEnergy = (*tridentIt)->GetEventEnergy();
      tNumIdentifiedMus = (*tridentIt)->NumIdentifiedMus();

      // Get particle objects
      Particle nu = (*tridentIt)->GetNu();
      std::vector<Particle> mus = (*tridentIt)->GetMus();
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

      tNuID = nu.ID();
      tNuPDG = nu.PDG();
      tNuStart[0] = nu.Start().X();
      tNuStart[1] = nu.Start().Y();
      tNuStart[2] = nu.Start().Z();

      tNumParticles = 0;

      tMuPlusID = muPlus.ID();
      tMuPlusPDG = muPlus.PDG();
      tMuPlusEnergy = muPlus.Energy();
      tMuPlusdEdx = muPlus.dEdx();
      tMuPlusLength = muPlus.Length();
      tMuPlusMomentum[0] = muPlus.Momentum().X();
      tMuPlusMomentum[1] = muPlus.Momentum().Y();
      tMuPlusMomentum[2] = muPlus.Momentum().Z();
      tMuPlusStart[0] = muPlus.Start().X();
      tMuPlusStart[1] = muPlus.Start().Y();
      tMuPlusStart[2] = muPlus.Start().Z();
      tMuPlusEnd[0] = muPlus.End().X();
      tMuPlusEnd[1] = muPlus.End().Y();
      tMuPlusEnd[2] = muPlus.End().Z();

      tMuMinusID = muMinus.ID();
      tMuMinusPDG = muMinus.PDG();
      tMuMinusEnergy = muMinus.Energy();
      tMuMinusdEdx = muMinus.dEdx();
      tMuMinusLength = muMinus.Length();
      tMuMinusMomentum[0] = muMinus.Momentum().X();
      tMuMinusMomentum[1] = muMinus.Momentum().Y();
      tMuMinusMomentum[2] = muMinus.Momentum().Z();
      tMuMinusStart[0] = muMinus.Start().X();
      tMuMinusStart[1] = muMinus.Start().Y();
      tMuMinusStart[2] = muMinus.Start().Z();
      tMuMinusEnd[0] = muMinus.End().X();
      tMuMinusEnd[1] = muMinus.End().Y();
      tMuMinusEnd[2] = muMinus.End().Z();

      tLongMuonID = longMuon.ID();
      tLongMuonPDG = longMuon.PDG();
      tLongMuonEnergy = longMuon.Energy();
      tLongMuondEdx = longMuon.dEdx();
      tLongMuonLength = longMuon.Length();
      tLongMuonMomentum[0] = longMuon.Momentum().X();
      tLongMuonMomentum[1] = longMuon.Momentum().Y();
      tLongMuonMomentum[2] = longMuon.Momentum().Z();
      tLongMuonStart[0] = longMuon.Start().X();
      tLongMuonStart[1] = longMuon.Start().Y();
      tLongMuonStart[2] = longMuon.Start().Z();
      tLongMuonEnd[0] = longMuon.End().X();
      tLongMuonEnd[1] = longMuon.End().Y();
      tLongMuonEnd[2] = longMuon.End().Z();

      tShortMuonID = shortMuon.ID();
      tShortMuonPDG = shortMuon.PDG();
      tShortMuonEnergy = shortMuon.Energy();
      tShortMuondEdx = shortMuon.dEdx();
      tShortMuonLength = shortMuon.Length();
      tShortMuonMomentum[0] = shortMuon.Momentum().X();
      tShortMuonMomentum[1] = shortMuon.Momentum().Y();
      tShortMuonMomentum[2] = shortMuon.Momentum().Z();
      tShortMuonStart[0] = shortMuon.Start().X();
      tShortMuonStart[1] = shortMuon.Start().Y();
      tShortMuonStart[2] = shortMuon.Start().Z();
      tShortMuonEnd[0] = shortMuon.End().X();
      tShortMuonEnd[1] = shortMuon.End().Y();
      tShortMuonEnd[2] = shortMuon.End().Z();

      outT->Fill();

    }

    outF->cd();
    outT->Write();
    delete outT;

  } // signal/bg

  outF->Close();
  delete outF;

}
