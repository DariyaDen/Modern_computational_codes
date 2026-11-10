#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    auto* preVolume =
        step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

    if (!preVolume) return;

    if (preVolume->GetName() == "HPGeActive") {
        fEventAction->AddEdep(step->GetTotalEnergyDeposit());
    }
}
