#include "Particle.h"

// Default constructor
Particle::Particle() {
  fID              = -1;
  fPDG             = -1;
  fAssociatedTrack = false;
  fEnergy          = -9999.;
  fDepositedEnergy = -9999.;
  fFamilyEnergy    = -9999.;
  fStart           = TVector3(-9999.,-9999.,-9999.);
  fEnd             = TVector3(-9999., -9999., -9999.);
  fInitialMomentum = TVector3(-9999., -9999., -9999.);
  fFinalMomentum   = TVector3(-9999., -9999., -9999.);
  fParticleStart   = TVector3(-9999., -9999., -9999.);
  fLength          = -9999.;
  fDirection       = TVector3(-9999.,-9999.,-9999.);
  fdEdx            = -9999.;
  fVertexFraction  = -9999.;
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

void Particle::SetDepositedEnergy(double energy) {
  fDepositedEnergy = energy;
}

void Particle::SetFamilyEnergy(double energy) {
  fFamilyEnergy = energy;
}

void Particle::SetParticleStart(TVector3 start) {
  fParticleStart = start;
}

void Particle::SetStart(TVector3 start) {
  fStart = start;
}

void Particle::SetEnd(TVector3 end) {
  fEnd = end;
}

void Particle::SetInitialMomentum(TVector3 mom) {
  fInitialMomentum = mom;
}

void Particle::SetFinalMomentum(TVector3 mom) {
  fFinalMomentum = mom;
}

void Particle::SetLength(double length) {
  fLength = length;
}

void Particle::SetDirection(TVector3 direction) {
  fDirection = direction;
}

void Particle::SetdEdx(double dEdx) {
  fdEdx = dEdx;
}

void Particle::SetVertexFraction(double vertexFraction) {
  fVertexFraction = vertexFraction;
}

int Particle::ID() { return fID; }
int Particle::PDG() { return fPDG; }
bool Particle::AssociatedTrack() { return fAssociatedTrack; }
double Particle::Energy() { return fEnergy; }
double Particle::DepositedEnergy() { return fDepositedEnergy; }
double Particle::FamilyEnergy() { return fFamilyEnergy; }
TVector3 Particle::Start() { return fStart; }
TVector3 Particle::End() { return fEnd; }
TVector3 Particle::InitialMomentum() { return fInitialMomentum; }
TVector3 Particle::FinalMomentum() { return fFinalMomentum; }
TVector3 Particle::ParticleStart() { return fParticleStart; }
double Particle::Length() { return fLength; }
TVector3 Particle::Direction() { return fDirection; }
double Particle::dEdx() { return fdEdx; }
double Particle::VertexFraction() { return fVertexFraction; }
