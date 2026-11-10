#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::Build() const
{
  auto* generator = new PrimaryGeneratorAction;
  auto* eventAction = new EventAction;
  SetUserAction(generator);
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction(eventAction));
}
