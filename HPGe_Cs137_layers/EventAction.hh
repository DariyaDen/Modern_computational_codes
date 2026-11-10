#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <fstream>

class G4Event;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override;
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;
  void AddEdep(G4double edep) { fEdep += edep; }

private:
  G4double fEdep = 0.0;
  std::ofstream fOutput;
};

#endif
