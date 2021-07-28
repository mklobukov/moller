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
                               PrimaryGeneratorAction* prim, RunAction* ruAct)
 : G4UserSteppingAction(),
   fDetector(det), fPrimary(prim), fRunAction(ruAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    //long seedIndex = CLHEP::HepRandom::getTheSeed();
    //long seedOne = CLHEP::HepRandom::getTheSeeds()[0];
    //long seedTwo = CLHEP::HepRandom::getTheSeeds()[1];
    ////G4cout << "Current seeds: " << "\n";
    //G4cout << seedIndex << " " << seedOne << " " << seedTwo << std::endl;

    G4double thetaCOM = -1; // will stay -1 if no moller scatter
    G4double thetaCOM2;

    G4int trackID = aStep->GetTrack()->GetTrackID();
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* endPoint = aStep->GetPostStepPoint();

    G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
    G4double primaryKineticEnergy = endPoint->GetKineticEnergy();

    fRunAction->CountProcesses(procName);
    G4Track* aTrack = aStep->GetTrack();
    const G4ParticleDefinition* part =
        aTrack->GetDynamicParticle()->GetDefinition();

    if (part->GetParticleName() == "e-" &&
        aStep->GetNumberOfSecondariesInCurrentStep() == 1 &&
        aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetParticleDefinition()->GetParticleName() == "e-" &&
        primaryKineticEnergy > 1000. &&
        aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetKineticEnergy() > 1000.)
    {
        // G4cout << "Primary's Kinetic Energies: " << prePoint->GetKineticEnergy() / (GeV) << " " << endPoint->GetKineticEnergy() / (GeV) << "\n";
        G4cout << "Primary's Total Energies: " << prePoint->GetTotalEnergy() / (GeV) << " " << endPoint->GetTotalEnergy() / (GeV) << "\n";
        // G4cout << "Processes: " << prePoint->GetProcessDefinedStep()->GetProcessName() << " " << endPoint->GetProcessDefinedStep()->GetProcessName() << "\n";
        G4cout << "Primary particle ID: " << aStep->GetTrack()->GetTrackID() << "\n";

        G4cout << "Secondary's energy: " << aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetKineticEnergy() << "\n";
        G4cout << "Current: " << aStep->GetTrack()->GetTrackID() << " Parent: " << aStep->GetTrack()->GetParentID() << "\n";
    }

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

#pragma region Print general info
        //G4cout << "Posiion of box: " << fDetector->GetBox()->GetObjectTranslation() << "\n";
        //G4cout << "Prepoint volume: " << prePoint->GetTouchableHandle()->GetVolume()->GetName() << "\n";
        //G4cout << "Prepoint: " << prePoint->GetPosition() << "\n";
        //G4cout <<  "Endpoint volume: " << endPoint->GetTouchableHandle()->GetVolume()->GetName() << "\n";
        //G4cout << "Endpoint: " << endPoint->GetPosition() << "\n";
        //G4cout << "Origin: " << aTrack->GetLogicalVolumeAtVertex()->GetName() << "\n";
        //G4cout << aStep->GetTrack()->GetTrackID() << " " << aStep->GetTrack()->GetParentID() << "\n";
        //G4cout << "Energies: " << prePoint->GetKineticEnergy() / (GeV) << " " << endPoint->GetKineticEnergy() / (GeV) << "\n";
        //// G4cout << trackID << "\n";
        //G4cout << "Welcome to the event of interest!" << "\n";
        //G4cout << part->GetParticleName() << "\n";
#pragma endregion

        // Calculate COM angle 
        G4double m = endPoint->GetMass(); // should be 0.511 MeV
        G4double EL = prePoint->GetTotalEnergy(); // should be 11,000 MeV

        G4cout << "mass : " << m << ", incident E: " << EL << "\n";
        G4double cosThetaCOM = (EL - costheta * EL - m) / (-EL + costheta * EL - m);
        G4double cosThetaCOM2 = (EL - costheta2 * EL - m) / (-EL + costheta2 * EL - m);
        thetaCOM = acos(cosThetaCOM);
        thetaCOM2 = acos(cosThetaCOM2);
        fRunAction->SetAngleCOM(thetaCOM);

#pragma region Previous attempt at COM Calculation

//// Equations taken from http://www.np.ph.bham.ac.uk/research_resources/programs/ckin/kinematics.pdf
//
//// Calculate vic (before collision):
////      Momentum of particle 0:
//G4double p0 = prePoint->GetMomentum().mag();
//G4cout << "p0=" << p0 << "\n"; G4cout << "P0 vector: " << prePoint->GetMomentum() << "\n";
//G4cout << "p0=" << "\n"; G4cout << "P0 vector: " << prePoint->GetMomentum() << "\n";
////      KE of particle zero
//G4double E0 = prePoint->GetKineticEnergy(); G4cout << "E0=" << E0 << "\n";
////      KE of particle one
//G4double E1 = aStep->GetSecondaryInCurrentStep()->at(0)->GetVertexKineticEnergy(); G4cout << "E1=" << E1 << "\n";

//G4double vic = p0 / (E0 + E1); G4cout << "vic=" << vic << "\n";

//// Calculate gamma_ic 
//// gamma = 1 / sqrt( 1 - (v/c)^2 )
//G4double c = 299792458.0; // meters per sec
//G4double gamma_ic = 1. / sqrt(1. - pow(vic / c, 2));

//// Calculate vfc /c (equation 37):        
////      mass of the first particle before collision:
//G4double m0 = prePoint->GetMass(); G4cout << "Mass incident e-: " << m0 << "\n";

////      mass of the second particle before collision
//G4double m1 = 0.5109989461; //mass of electron in MeV
////      mass of the first particle after collision:
//G4double m2 = endPoint->GetMass(); G4cout << "m2: " << m2 << "\n"; // MeV
////      mass of the second particle after collision:
//G4double m3 = aStep->GetSecondaryInCurrentStep()->at(0)->GetDynamicParticle()->GetMass(); G4cout << "m3: " << m3 << "\n";

////      Get numerator of equation 37:
//G4double numerator = gamma_ic * (m0 + m1) * vic / ((m2 + m3) * c);
////      Get denominator of equation 37:
//G4double denominator = 1. + (gamma_ic * (m0 + m1) * vic / ((m2 + m3) * c));
////      Calculate vfc / c (eq 37):
//G4double vfc_c = sqrt(numerator / denominator);
//

//// Calculate gamma_fc
//G4double gamma_fc = 1. / sqrt(1. - pow(vfc_c, 2));

//// Calculate COM angle (eq 46):
////      Momentum of particle 2:
//G4double p2 = endPoint->GetMomentum().mag();
////      Energy of particle 2:
G4double E2 = endPoint->GetKineticEnergy();
////       Numerator of Eq 46:
//G4double COMnumerator = p2 * c * sin(scatterAngle);
////      Denominator of Eq 46:
//G4double COMdenominator = gamma_fc * (p2 * c * cos(scatterAngle) - vfc_c * E2);
////      COM angle Eq 46:
//G4double COMangle = atan(COMnumerator / COMdenominator);
#pragma endregion - produced too small angle for 5.5 GeV point


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

        //long seedIndex = CLHEP::HepRandom::getTheSeed();
        //long seedOne = CLHEP::HepRandom::getTheSeeds()[0];
        //long seedTwo = CLHEP::HepRandom::getTheSeeds()[1];

        ////std::ofstream seedsFile;
        ////seedsFile.open("seedsMinus", std::ios_base::app);
        ////seedsFile << seedIndex << " " << seedOne << " " << seedTwo << "\n";
        ////seedsFile.close();

        //G4cout << seedIndex << " " << seedOne << " " << seedTwo << std::endl;
  }

    // This step generates a lot of data. Toggle the conditional to enable/disable
    if (bool recordDetectorData = true) {
        if (prePoint->GetTouchableHandle()->GetVolume() == fDetector->GetBox() &&
            endPoint->GetTouchableHandle()->GetVolume() == fDetector->GetWorld() &&
            endPoint->GetTotalEnergy() > sqrt(100000.)) // > 1 MeV
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

            // G4cout << "mass : " << m << ", prePoint E: " << EL << ", endPointE: " << endPoint->GetTotalEnergy() << "\n";
            G4double cosThetaCOM = (EL - costheta * EL - m) / (-EL + costheta * EL - m);
            thetaCOM = acos(cosThetaCOM);

            fRunAction->FillData(part, kinEnergy, costheta, phi, polZ, procName, labScatterAngle);
        }
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......