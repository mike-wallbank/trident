#include "Particle.h"

// Default constructor
Particle::Particle() {
  fID = -1;
  fPDG = -1;
  fAssociatedTrack = false;
  fEnergy = 0.;
  fStart = TVector3(0,0,0);
  fLength = 0.;
  fDirection = TVector3(0,0,0);
  fMomentum = TVector3(0,0,0);
  fdEdx = 0.;
}

Particle::Particle(int id) {
  fID = id;
}

void Particle::SetPDG(int pdg) {
  fPDG = pdg;
}

void Particle::SetAssociatedTrack(bool track) {
  fAssociatedTrack = track;
}

void Particle::SetEnergy(double energy) {
  fEnergy = energy;
}

void Particle::SetStart(TVector3 start) {
  fStart = start;
}

void Particle::SetEnd(TVector3 end) {
  fEnd = end;
}

void Particle::SetStartTrack(TVector3 start) {
  fStartTrack = start;
}

void Particle::SetLength(double length) {
  fLength = length;
}

void Particle::SetDirection(TVector3 direction) {
  fDirection = direction;
}

void Particle::SetMomentum(TVector3 momentum) {
  fMomentum = momentum;
}

void Particle::SetdEdx(double dEdx) {
  fdEdx = dEdx;
}

int Particle::ID() { return fID; }
int Particle::PDG() { return fPDG; }
bool Particle::AssociatedTrack() { return fAssociatedTrack; }
double Particle::Energy() { return fEnergy; }
TVector3 Particle::Start() { return fStart; }
TVector3 Particle::End() { return fEnd; }
TVector3 Particle::StartTrack() { return fStartTrack; }
double Particle::Length() { return fLength; }
TVector3 Particle::Direction() { return fDirection; }
TVector3 Particle::Momentum() { return fMomentum; }
double Particle::dEdx() { return fdEdx; }
