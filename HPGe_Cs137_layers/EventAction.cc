#include "EventAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

EventAction::EventAction()
  : fOutput("edep.csv", std::ios::out)
{
  if (fOutput.is_open()) {
    fOutput << "event_id,edep_keV\n";
  } else {
    G4cerr << "ERROR: cannot open edep.csv for writing." << G4endl;
  }
}

EventAction::~EventAction()
{
  if (fOutput.is_open()) fOutput.close();
}

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
  if (fOutput.is_open()) {
    fOutput << event->GetEventID() << ',' << fEdep/keV << '\n';
  }
}
