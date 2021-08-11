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
/// \file polarisation/Pol01/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4PolarizationHelper.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
                               PrimaryGeneratorAction* prim, RunAction* ruAct, G4double energyCutoff)
 : G4UserSteppingAction(),
   fDetector(det), fPrimary(prim), fRunAction(ruAct), fLowerEnergyCutoff(1)
{ 
    fLowerEnergyCutoff = energyCutoff;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
#pragma region Toggle simulation setups. Rebuild for these to take effect. 
    bool collectMollerData = false;
    bool collectPhotonPolarization = true;
    bool recordDetectorData = false;
#pragma endregion 

    G4double thetaCOM = -1; // will stay -1 if no moller scatter
    G4double thetaCOM2;

    G4int trackID = aStep->GetTrack()->GetTrackID();
    G4int parentTrackId = aStep->GetTrack()->GetParentID();
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* endPoint = aStep->GetPostStepPoint();

    G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
    G4double primaryKineticEnergy = endPoint->GetKineticEnergy();

    fRunAction->CountProcesses(procName);
    G4Track* aTrack = aStep->GetTrack();
    const G4ParticleDefinition* part =
        aTrack->GetDynamicParticle()->GetDefinition();

    if (collectMollerData == false)
    {
        if (part->GetParticleName() == "e-" &&
            aStep->GetNumberOfSecondariesInCurrentStep() == 1 &&
            aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "e-" &&
            primaryKineticEnergy > 1000. &&
            aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetKineticEnergy() > 1000.)
        {
            // Calculate lab scattering angle for the primary
            G4ThreeVector direction = endPoint->GetMomentumDirection().unit();
            G4ThreeVector beamDirection =
                fPrimary->GetParticleGun()->GetParticleMomentumDirection().unit();
            G4double costheta = direction * beamDirection; // dot product of two unit vectors
            G4double scatterAngle = std::acos(costheta); // returns angle in [0, pi]

            // Calculate lab scattering angle for the secondary
            direction = aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetMomentumDirection().unit();
            G4double costheta2 = direction * beamDirection; // dot product of two unit vectors

            G4double E = endPoint->GetTotalEnergy() * costheta +
                aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetTotalEnergy() * costheta2;

            G4cout << "E1 * cos(theta1) ) + E2 * cos(theta2) = " << E << "\n";

#pragma region Print general info of this step 
            //G4cout << "Posiion of box: " << fDetector->GetBox()->GetObjectTranslation() << "\n";
            //G4cout << "Prepoint volume: " << prePoint->GetTouchableHandle()->GetVolume()->GetName() << "\n";
            //G4cout << "Prepoint: " << prePoint->GetPosition() << "\n";
            //G4cout <<  "Endpoint volume: " << endPoint->GetTouchableHandle()->GetVolume()->GetName() << "\n";
            //G4cout << "Endpoint: " << endPoint->GetPosition() << "\n";
            //G4cout << "Origin: " << aTrack->GetLogicalVolumeAtVertex()->GetName() << "\n";
            //G4cout << aStep->GetTrack()->GetTrackID() << " " << aStep->GetTrack()->GetParentID() << "\n";
            //G4cout << "Energies: " << prePoint->GetKineticEnergy() / (GeV) << " " << endPoint->GetKineticEnergy() / (GeV) << "\n";
            //// G4cout << trackID << "\n";
            //G4cout << part->GetParticleName() << "\n";
#pragma endregion 

            // Calculate COM angle 
            G4double m = endPoint->GetMass(); 
            G4double EL = prePoint->GetTotalEnergy(); 

            G4cout << "mass : " << m << ", incident E: " << EL << "\n";
            G4double cosThetaCOM = (EL - costheta * EL - m) / (-EL + costheta * EL - m);
            G4double cosThetaCOM2 = (EL - costheta2 * EL - m) / (-EL + costheta2 * EL - m);
            thetaCOM = acos(cosThetaCOM);
            thetaCOM2 = acos(cosThetaCOM2);
            fRunAction->SetAngleCOM(thetaCOM);

            G4double E2 = endPoint->GetKineticEnergy();

            // Calculate COM relativistic energy
            G4double energyCOM = 2 * sqrt(m * (EL + m) / 2.0);

            G4cout << "COM ANGLE: " << thetaCOM << "\n";
            G4cout << "Scatter ANGLE: " << scatterAngle << "\n";
            G4cout << "Energy: " << E2 << std::endl;
            fRunAction->FillData(
                part,
                endPoint->GetKineticEnergy(),
                procName,
                thetaCOM,
                cosThetaCOM,
                scatterAngle,
                EL,
                energyCOM,
                acos(costheta2),
                thetaCOM2,
                cosThetaCOM2,
                aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetKineticEnergy());

            // collect and store seeds
            //long seedIndex = CLHEP::HepRandom::getTheSeed();
            //long seedOne = CLHEP::HepRandom::getTheSeeds()[0];
            //long seedTwo = CLHEP::HepRandom::getTheSeeds()[1];

            ////std::ofstream seedsFile;
            ////seedsFile.open("seedsMinus", std::ios_base::app);
            ////seedsFile << seedIndex << " " << seedOne << " " << seedTwo << "\n";
            ////seedsFile.close();

            //G4cout << seedIndex << " " << seedOne << " " << seedTwo << std::endl;
        }
    }

    if (collectPhotonPolarization == true) {
        if (part->GetParticleName() == "e-" &&
            aStep->GetNumberOfSecondariesInCurrentStep() == 1 &&
            aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "gamma" &&
            aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetKineticEnergy() > 100.)
        {
            // record stuff for secondary photon
            G4double photonEnergy = aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetTotalEnergy();
            G4double photonPol = aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetPolarization().z();
            const G4ParticleDefinition* photon = aStep->GetSecondaryInCurrentStep()->at(0)->GetParticleDefinition();
            
            fRunAction->FillBrehmData(
                endPoint->GetTotalEnergy(),
                photonEnergy,
                endPoint->GetPolarization().z(),
                photonPol);
        }
    }

    // This step generates a lot of data. Toggle the conditional to enable/disable
    if (recordDetectorData == false) {
        if (prePoint->GetTouchableHandle()->GetVolume() == fDetector->GetBox() &&
            endPoint->GetTouchableHandle()->GetVolume() == fDetector->GetWorld() &&
            endPoint->GetTotalEnergy() > fLowerEnergyCutoff)
        {
            G4ThreeVector position = endPoint->GetPosition();
            G4ThreeVector direction = endPoint->GetMomentumDirection();
            G4double kinEnergy = endPoint->GetKineticEnergy();

            G4ThreeVector beamDirection =
                fPrimary->GetParticleGun()->GetParticleMomentumDirection();
            G4double polZ = endPoint->GetPolarization().z();
            G4double costheta = direction * beamDirection;
            G4double labScatterAngle = acos(costheta);

            G4double xdir =
                direction * G4PolarizationHelper::GetParticleFrameX(beamDirection);
            G4double ydir =
                direction * G4PolarizationHelper::GetParticleFrameY(beamDirection);

            G4double phi = std::atan2(ydir, xdir);

            // Calculate COM angle 
            G4double m = endPoint->GetMass(); // should be 0.511 MeV
            G4double EL = prePoint->GetTotalEnergy(); // should be 11,000 MeV

            G4double cosThetaCOM = (EL - costheta * EL - m) / (-EL + costheta * EL - m);
            thetaCOM = acos(cosThetaCOM);

            if (bool lightweight = true) {
                G4String name = part->GetParticleName();
                if (name.compare("gamma") == 0) name = "p";
                fRunAction->FillDataLight(name, endPoint->GetTotalEnergy());
            }
            else
            {
                fRunAction->FillData(part, kinEnergy, costheta, phi, polZ, procName, labScatterAngle);
            }
        }
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......