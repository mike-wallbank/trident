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
  void SetStart(TVector3 start);
  void SetEnd(TVector3 end);
  void SetStartTrack(TVector3 start);
  void SetLength(double length);
  void SetDirection(TVector3 direction);
  void SetMomentum(TVector3 momentum);
  void SetdEdx(double dEdx);

  int ID();
  int PDG();
  bool AssociatedTrack();
  double Energy();
  TVector3 Start();
  TVector3 End();
  TVector3 StartTrack();
  double Length();
  TVector3 Direction();
  TVector3 Momentum();
  double dEdx();

private:

  int fID;
  int fPDG;
  bool fAssociatedTrack;
  double fEnergy;
  TVector3 fStart;
  TVector3 fEnd;
  TVector3 fStartTrack;
  double fLength;
  TVector3 fDirection;
  TVector3 fMomentum;
  double fdEdx;

};

#endif
