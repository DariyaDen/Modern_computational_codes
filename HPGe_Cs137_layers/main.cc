#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "FTFP_BERT.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4SystemOfUnits.hh"

int main(int argc, char** argv)
{
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new FTFP_BERT);
  runManager->SetUserInitialization(new ActionInitialization);

  auto* uiManager = G4UImanager::GetUIpointer();

  if (argc > 1) {
    G4String command = "/control/execute ";
    uiManager->ApplyCommand(command + argv[1]);
  } else {
    auto* ui = new G4UIExecutive(argc, argv);
    auto* visManager = new G4VisExecutive;
    visManager->Initialize();
    uiManager->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete visManager;
    delete ui;
  }

  delete runManager;
  return 0;
}
