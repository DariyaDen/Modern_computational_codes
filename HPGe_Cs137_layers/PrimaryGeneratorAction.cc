#include "PrimaryGeneratorAction.hh"
#include "DetectorParameters.hh"

#include "G4ParticleGun.hh"
#include "G4Gamma.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  using namespace DetectorParameters;

  fParticleGun = new G4ParticleGun(1);
  fParticleGun->SetParticleDefinition(G4Gamma::GammaDefinition());
  fParticleGun->SetParticleEnergy(sourceEnergy);

  // Monoenergetic Cs-137 gamma line.
  // The source is placed sourceDistance from the +Z crystal face and aimed
  // toward the crystal (-Z). The detector cap is not yet modelled, so this
  // is a distance from the crystal face, not from the physical detector cap.
  fParticleGun->SetParticlePosition(
    G4ThreeVector(0.0, 0.0, H / 2.0 + sourceDistance));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  fParticleGun->GeneratePrimaryVertex(event);
}
