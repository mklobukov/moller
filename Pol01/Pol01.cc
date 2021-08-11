//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file polarisation/Pol01/Pol01.cc
/// \brief Main program of the polarisation/Pol01 example
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "FTFP_BERT.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // default output root file name. can be overridden via cmd params
  G4String outputFileName = "data.root";

  // default cutoff energy for moller data collection. can be overridden via cmd params
  G4double cutoffEnergy = 1.;

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
  if (argv[2])
  {
      outputFileName = argv[2];
  }
  if (argv[3])
  {
      cutoffEnergy = atof(argv[3]);
  }

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  // CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine());

#ifdef G4MULTITHREADED
  auto* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(2 * G4Threading::G4GetNumberOfCores());
  G4cout << "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......\n";
  G4cout << "MULTITHREADED MODE \n";
  G4cout << "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......\n";

#else 
  // Construct a serial run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  G4cout << "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......\n";
  G4cout << "SINGLE THREADED MODE \n";
  G4cout << "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......\n";
#endif

  // Import seeds collection:
  std::queue<std::vector<long>> seeds;
  std::ifstream infile("seeds");
  long seedIndex, seedOne, seedTwo;
  while (infile >> seedIndex >> seedOne >> seedTwo)
  {
      seeds.push(std::vector<long> {seedIndex, seedOne, seedTwo});
  }

  // set mandatory initialization classes
  DetectorConstruction* det;
  PrimaryGeneratorAction* prim;
  runManager->SetUserInitialization(det = new DetectorConstruction);
  //runManager->SetUserInitialization(new FTFP_BERT(1)); // physics list
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserAction(prim = new PrimaryGeneratorAction(det));

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // set user action classes
  RunAction* run;
  runManager->SetUserAction(run = new RunAction(det,prim, outputFileName));
  runManager->SetUserAction(new EventAction(run, seeds));
  runManager->SetUserAction(new SteppingAction(det,prim,run,cutoffEnergy));

  // get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)   // Define UI terminal for interactive mode
    {
      ui->SessionStart();
      delete ui;
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  // job termination
  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
