#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class EventAction;
class G4Step;

class SteppingAction : public G4UserSteppingAction
{
public:
  explicit SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}
  ~SteppingAction() override = default;
  void UserSteppingAction(const G4Step*) override;

private:
  EventAction* fEventAction = nullptr;
};

#endif
