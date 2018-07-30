#ifndef PARTICLE_H
#define PARTICLE_H

#include "TMath.h"
#include "TVector3.h"

class Particle {
public:

  Particle();
  Particle(int id);
  ~Particle() { }

  void SetPDG(int pdg);
  void SetAssociatedTrack(bool track);
  void SetEnergy(double energy);
  void SetDepositedEnergy(double energy);
  void SetFamilyEnergy(double energy);
  void SetStart(TVector3 start);
  void SetEnd(TVector3 end);
  void SetInitialMomentum(TVector3 mom);
  void SetFinalMomentum(TVector3 mom);
  void SetParticleStart(TVector3 start);
  void SetLength(double length);
  void SetDirection(TVector3 direction);
  void SetdEdx(double dEdx);
  void SetVertexFraction(double vertexFraction);

  int ID();
  int PDG();
  bool AssociatedTrack();
  double Energy();
  double DepositedEnergy();
  double FamilyEnergy();
  TVector3 Start();
  TVector3 End();
  TVector3 InitialMomentum();
  TVector3 FinalMomentum();
  TVector3 ParticleStart();
  double Length();
  TVector3 Direction();
  double dEdx();
  double VertexFraction();

private:

  int fID;
  int fPDG;
  bool fAssociatedTrack;
  double fEnergy;
  double fDepositedEnergy;
  double fFamilyEnergy;
  TVector3 fStart;
  TVector3 fEnd;
  TVector3 fInitialMomentum;
  TVector3 fFinalMomentum;
  TVector3 fParticleStart;
  double fLength;
  TVector3 fDirection;
  double fdEdx;
  double fVertexFraction;

};

#endif
