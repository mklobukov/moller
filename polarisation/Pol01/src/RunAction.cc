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
/// \file polarisation/Pol01/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iomanip>

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string> 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
: G4UserRunAction(),
  fGamma(0), fElectron(0), fPositron(0),
  fDetector(det), fPrimary(prim), fProcCounter(0), fAnalysisManager(0), 
  fTotalEventCount(0),
  fPhotonStats(), fElectronStats(), fPositronStats(), fOtherParticleStats(),
    fAngleCOM(0)
{
  fGamma = G4Gamma::Gamma();
  fElectron = G4Electron::Electron();
  fPositron = G4Positron::Positron();

  fAngleCOM = -1.;

  BookHisto();
  BookNTuples();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

G4double RunAction::GetAngleCOM()
{
    return fAngleCOM;
}

void RunAction::SetAngleCOM(G4double angle)
{
    fAngleCOM = angle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  // save Rndm status
  //  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  //  CLHEP::HepRandom::showEngineStatus();

  time_t t; // t passed as argument in function time()
  struct tm* tt; // decalring variable for localtime()
  time(&t); //passing argument to time()
  tt = localtime(&t);
  G4cout << "Beginning of run time: " << asctime(tt);

  if (fProcCounter) delete fProcCounter;
  fProcCounter = new ProcessesCount;
  fTotalEventCount = 0;
  fPhotonStats.Clear();
  fElectronStats.Clear();
  fPositronStats.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunAction::FillData(const G4ParticleDefinition* particle,
    G4double kinEnergy,
    G4String processName,
    G4double cmAngle,
    G4double cosCmAngle,
    G4double labScatterAngle,
    G4double totalEnergy,
    G4double centerOfMassEnergy,
    G4double secondaryScatterAngle,
    G4double secondaryCOMAngle,
    G4double secondaryCOMAngleCosine,
    G4double secondaryKE)
{
    fAnalysisManager->FillNtupleSColumn(0, 0, particle->GetParticleName());
    fAnalysisManager->FillNtupleSColumn(0, 1, processName);
    fAnalysisManager->FillNtupleDColumn(0, 2, cmAngle);
    fAnalysisManager->FillNtupleDColumn(0, 3, kinEnergy);
    fAnalysisManager->FillNtupleDColumn(0, 4, cosCmAngle);
    fAnalysisManager->FillNtupleDColumn(0, 5, labScatterAngle);    
    fAnalysisManager->FillNtupleDColumn(0, 6, totalEnergy);
    fAnalysisManager->FillNtupleDColumn(0, 7, centerOfMassEnergy);
    fAnalysisManager->FillNtupleDColumn(0, 8, secondaryScatterAngle);
    fAnalysisManager->FillNtupleDColumn(0, 9, secondaryCOMAngle);
    fAnalysisManager->FillNtupleDColumn(0, 10, secondaryCOMAngleCosine);
    fAnalysisManager->FillNtupleDColumn(0, 11, secondaryKE);
    fAnalysisManager->AddNtupleRow(0);
}


void RunAction::FillData(const G4ParticleDefinition* particle,
    G4double kinEnergy, G4double costheta,
    G4double phi,
    G4double longitudinalPolarization,
    G4String processName = "none")
{
    G4int id = -1;
    ParticleStatistics currentParticleStatistics;

    if (particle == fGamma)
    {
        fPhotonStats.FillData(kinEnergy, costheta, longitudinalPolarization);
        currentParticleStatistics = fPhotonStats;
        if (fAnalysisManager) { id = 0; }
    }
    else if (particle == fElectron)
    {
        fElectronStats.FillData(kinEnergy, costheta, longitudinalPolarization);
        currentParticleStatistics = fElectronStats;
        //if(fAnalysisManager) { id = 5; } 
        if (fAnalysisManager) { id = 5; }
    }
    else if (particle == fPositron)
    {
        fPositronStats.FillData(kinEnergy, costheta, longitudinalPolarization);
        currentParticleStatistics = fPositronStats;
        //if(fAnalysisManager) { id = 9; } 
        if (fAnalysisManager) { id = 10; }
    }
    else
    {
        fOtherParticleStats.FillData(kinEnergy, costheta, longitudinalPolarization);
    }

    if (id > -1) {
        fAnalysisManager->FillH1(id, kinEnergy, 1.0);
        fAnalysisManager->FillH1(id + 1, costheta, 1.0);
        fAnalysisManager->FillH1(id + 2, phi, 1.0);
        fAnalysisManager->FillH1(id + 3, longitudinalPolarization, 1.0);
        // fAnalysisManager->FillH1(id + 4, currentParticleStatistics.GetCurrentNumber(), 1.0);
    }

    // G4cout << particle->GetParticleName() << " " << longitudinalPolarization << " " << kinEnergy << " " << currentParticleStatistics.GetCurrentNumber() << "\n";
    fAnalysisManager->FillNtupleSColumn(1, 0, particle->GetParticleName());
    fAnalysisManager->FillNtupleDColumn(1, 1, longitudinalPolarization);
    fAnalysisManager->FillNtupleDColumn(1, 2, kinEnergy);
    fAnalysisManager->FillNtupleDColumn(1, 3, currentParticleStatistics.GetCurrentNumber());
    fAnalysisManager->FillNtupleSColumn(1, 4, processName);
    fAnalysisManager->AddNtupleRow(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BookHisto()
{
  // Always creating analysis manager
  fAnalysisManager = G4AnalysisManager::Instance();
  // fAnalysisManager->SetActivation(true);
  fAnalysisManager->SetVerboseLevel(1);

  // Open root file:
  srand((unsigned)time(0));

  G4String incidentPolarization = std::to_string((int)fPrimary->GetParticleGun()->GetParticlePolarization().z());
  G4String thickness = std::to_string(fDetector->GetBoxSizeZ());
  
  // Randomized assignment of file name
  // fAnalysisManager->OpenFile(thickness + "mm_" + std::to_string(rand()) + ".root");
  fAnalysisManager->OpenFile("data-minus.root");
  fAnalysisManager->SetFirstHistoId(0);
  fAnalysisManager->SetFirstNtupleId(0);

  // Creating an 1-dimensional histograms in the root directory of the tree
  //const G4String id[] = { "h1", "h2", "h3", "h4", "h5", 
  //                        "h6", "h7", "h8", "h9", "h10", "h11", "h12"};
  //const G4String title[] = 
  //              { "Photon Energy distribution",                      //1
  //                "Gamma Cos(Theta) distribution",                  //2
  //                "Gamma Phi angular distribution",                 //3
  //                "Gamma longitudinal Polarization",                //4
  //                "Electron Energy distribution",                   //5
  //                "Electron Cos(Theta) distribution",               //6
  //                "Electron Phi angular distribution",              //7
  //                "Electron longitudinal Polarization",             //8
  //                "Positron Energy distribution",                   //9
  //                "Positron Cos(Theta) distribution",               //10
  //                "Positron Phi angular distribution",              //11
  //                "Positron longitudinal Polarization"              //12
  //               };  
  
    const G4String id[] = { "h1", "h2", "h3", "h4", "h5", 
                          "h6", "h7", "h8", "h9", "h10", "h11", "h12", "h13", "h14", "h15"};
  
    G4cout << "NUMBER OF HISTOGRAMS: " << std::size(id);
    const G4String title[] = 
                { "Photon Energy Distribution",                      //1
                  "Gamma Cos(Theta) Distribution",                  //2
                  "Gamma Phi angular Distribution",                 //3
                  "Gamma longitudinal Polarization",                //4
                  "Number of photons per event",            //5 
                  "Electron Energy Distribution",                   //5
                  "Electron Cos(Theta) Distribution",               //6
                  "Electron Phi angular Distribution",              //7
                  "Electron longitudinal Polarization",             //8
                    "Number of electrons per event",
                  "Positron Energy Distribution",                   //9
                  "Positron Cos(Theta) Distribution",               //10
                  "Positron Phi angular Distribution",              //11
                  "Positron longitudinal Polarization",             //12
                    "Number of positrons per event"
                 };

    //const G4String axisTitle[] =
    //{
    //    "MeV",
    //    "Degree",
    //    "Degree",
    //    "",
    //    ""
    //};

  G4double vmin, vmax;
  G4String binningScheme;

  G4int nbins = 100;
  for(int i = 0; i < std::size(title); i++) {
      binningScheme = "linear";
      // G4int j = i % 5;
      switch (i % 5) {
      case(0):
          vmin = 0.; vmax = 13.5 * GeV;
          binningScheme = "log";
          break;
      case(1): 
          vmin = -1.; vmax = 1.;
          break;
      case(2): 
          vmin = 0.; vmax = pi;
          break;
      case(3):
          vmin = -0.5; vmax = 1.5;
          binningScheme = "log";
          break;
      case(4):
          vmin = -0.5; vmax = 100.;
          binningScheme = "log";
          break;
      default:
          break;
      }

    G4int ih = fAnalysisManager->CreateH1(id[i], title[i], nbins, vmin, vmax); // , axisTitle[i % 5]);
    fAnalysisManager->SetH1Activation(ih, false);
  }
}

void RunAction::BookNTuples() {
    // Always creating analysis manager
    //fAnalysisManager = G4AnalysisManager::Instance();
    //fAnalysisManager->SetActivation(true);
    //fAnalysisManager->SetVerboseLevel(1);

    // fAnalysisManager->OpenFile("ntuples");
    
    // NTuple #0: processes and center of mass angles and energies
    G4int currentNTuple = 0;
    fAnalysisManager->CreateNtuple("Polarized Ionization >1 GeV", "Processes, CM Angles, Energies, Cosines of CM Angle, Lab Scatter angle");
    fAnalysisManager->CreateNtupleSColumn(currentNTuple, "ParticleName");
    fAnalysisManager->CreateNtupleSColumn(currentNTuple, "ProcessName");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "CMAngle");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "PrimaryKineticEnergy");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "CosineCMAngle");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "LabScatterAngle");    
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "TotalEnergy");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "CenterOfMassEnergy");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "SecondaryScatterAngle");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "SecondaryCOMAngle");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "SecondaryCOMAngleCosine");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "SecondaryKineticEnergy");
    fAnalysisManager->FinishNtuple(currentNTuple);

    // NTuple #1: data for particles emerging from target
    currentNTuple++;
    fAnalysisManager->CreateNtuple("Particles emerging from target", "Number, Energy, Polarization");
    fAnalysisManager->CreateNtupleSColumn(currentNTuple, "ParticleName");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "Count");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "Energy");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "Polarization");
    fAnalysisManager->CreateNtupleSColumn(currentNTuple, "Process");
    fAnalysisManager->FinishNtuple(currentNTuple);

    // ntuple # 2: center of mass angle vs number of particles per event
    currentNTuple++;
    fAnalysisManager->CreateNtuple("COM angle and particles per event", "COM, nparticles");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "COMAngle");
    fAnalysisManager->CreateNtupleDColumn(currentNTuple, "NScattered");
    fAnalysisManager->FinishNtuple(currentNTuple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SaveHistoAndNtuples(G4int nevents)
{
  if(fAnalysisManager) {
    //G4double norm = 1.0/G4double(nevents);
    //for(int i=0; i<15; ++i) {
    //  fAnalysisManager->ScaleH1(i, norm);
    //}

    fAnalysisManager->Write();
    fAnalysisManager->CloseFile();
    delete fAnalysisManager;
    fAnalysisManager = 0;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName)
{
  // is the process already counted ?
  // *AS* change to std::map?!
  size_t nbProc = fProcCounter->size();
  size_t i = 0;
  while ((i<nbProc)&&((*fProcCounter)[i]->GetName()!=procName)) i++;
  if (i == nbProc) fProcCounter->push_back( new ProcessCount(procName));
  
  (*fProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
    time_t t; // t passed as argument in function time()
    struct tm* tt; // decalring variable for localtime()
    time(&t); //passing argument to time()
    tt = localtime(&t);
    G4cout << "End of run time " << asctime(tt);

  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  G4int  prec = G4cout.precision(5);
    
  G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = 
                            fPrimary->GetParticleGun()->GetParticleDefinition();
  G4String Particle = particle->GetParticleName();    
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();

  G4cout << "\tIncident e- beam polarization = " << fPrimary->GetParticleGun()->GetParticlePolarization() << "\n";
  G4cout << "\tTarget polarization = " << fDetector->GetTargetPolarization() << "\n";
  G4cout << "\tThe run consists of " << NbOfEvents << " "<< Particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(fDetector->GetBoxSizeZ(),"Length") << "of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  //frequency of processes
  G4cout << "\n Process calls frequency --->\n";
  for (size_t i=0; i< fProcCounter->size();i++) {
     G4String procName = (*fProcCounter)[i]->GetName();
     G4int    count    = (*fProcCounter)[i]->GetCounter(); 
     G4cout << "\t" << procName << " = " << count<<"\n";
  }
  
  if (fTotalEventCount == 0) return;
  
  G4cout<<" Gamma: \n";
  fPhotonStats.PrintResults(fTotalEventCount);

  G4cout<<" Electron: \n";
  fElectronStats.PrintResults(fTotalEventCount);
  G4double asymmetry = fElectronStats.CalculateAsymmetry(
      fPrimary->GetParticleGun()->GetParticlePolarization().Z,
      fTotalEventCount);
  G4cout << "Asymmetry = " << asymmetry << "\n";

  G4cout<<" Positron: \n";
  fPositronStats.PrintResults(fTotalEventCount);

  G4cout << "Other particles: \n";
  fOtherParticleStats.PrintResults(fTotalEventCount);

  G4cout << "Total number of scattered/produced particles: " <<
      fElectronStats.GetTotalNumber() + fPhotonStats.GetTotalNumber() +
      fPositronStats.GetTotalNumber() + fOtherParticleStats.GetTotalNumber() << "\n";

  //cross check from G4EmCalculator
  //  G4cout << "\n Verification from G4EmCalculator. \n";  
  //  G4EmCalculator emCal;

  //restore default format         
  G4cout.precision(prec);         

  // write out histograms  
  SaveHistoAndNtuples(NbOfEvents);

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

G4double RunAction::ParticleStatistics::CalculateAsymmetry(G4double incidentPolarization, G4int numberOfIncidentParticles) {
    G4cout << "Incident polarization = " << incidentPolarization << "\n";
    G4cout << "# incident particles = " << numberOfIncidentParticles << "\n";

    G4int numScatteredWithSamePol;
    std::vector<G4int> polarizationCounts = GetPolarizations();

    G4cout << "Pols: " << polarizationCounts.at(0) << " " << polarizationCounts.at(1) << " " << polarizationCounts.at(2) <<  "\n";
    switch ((int)incidentPolarization) {
    case(-1): numScatteredWithSamePol = polarizationCounts.at(0); break;
    case(0): numScatteredWithSamePol = polarizationCounts.at(1); break;
    case(1): numScatteredWithSamePol = polarizationCounts.at(2); break;
    default: numScatteredWithSamePol = INFINITY;
    }

    G4cout << " NUM Scattered = " << numScatteredWithSamePol << "\n";

    G4double asymmetry = ((G4double)numberOfIncidentParticles - (G4double)numScatteredWithSamePol) /
        ((G4double)numberOfIncidentParticles + (G4double)numScatteredWithSamePol);
    return asymmetry;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EventFinished()
{
  ++fTotalEventCount;

  // fill every fifth histogram with particle counts for the given event:
  // ids are 4, 9, 14
  fAnalysisManager->FillH1(4, fPhotonStats.GetCurrentNumber(), 1.0);
  fAnalysisManager->FillH1(9, fElectronStats.GetCurrentNumber(), 1.0);
  fAnalysisManager->FillH1(14, fPositronStats.GetCurrentNumber(), 1.0);

  if (fAngleCOM > -1)
  {
      fAnalysisManager->FillNtupleDColumn(2, 0, fAngleCOM);
      fAnalysisManager->FillNtupleDColumn(
          2,
          1,
          fElectronStats.GetCurrentNumber() + fPhotonStats.GetCurrentNumber() +
          fPositronStats.GetCurrentNumber() + fOtherParticleStats.GetCurrentNumber());
      fAnalysisManager->AddNtupleRow(2);

      fAngleCOM = -1; // reset the flag
  }

  // fill the scattering data here...
  fPhotonStats.EventFinished();
  fElectronStats.EventFinished();
  fPositronStats.EventFinished();
}

G4int RunAction::ParticleStatistics::GetCurrentNumber()
{
    return fCurrentNumber;
}

RunAction::ParticleStatistics RunAction::GetElectronStatistics() {
    return fElectronStats;
}

G4int RunAction::ParticleStatistics::GetTotalNumber()
{
    return fTotalNumber;
}

std::vector<G4int> RunAction::ParticleStatistics::GetPolarizations()
{
    static const G4int arr[] = { fNumNegativePol, fNumZeroPol, fNumPositivePol };
    std::vector<G4int> vec(arr, arr + sizeof(arr) / sizeof(arr[0]));
    return vec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::ParticleStatistics::ParticleStatistics()
  : fCurrentNumber(0),
    fTotalNumber(0), fTotalNumber2(0),
    fSumEnergy(0), fSumEnergy2(0),
    fSumPolarization(0), fSumPolarization2(0),
    fSumCosTheta(0), fSumCosTheta2(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::ParticleStatistics::~ParticleStatistics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::EventFinished()
{
  fTotalNumber+=fCurrentNumber;
  fTotalNumber2+=fCurrentNumber*fCurrentNumber;
  fCurrentNumber=0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics:: FillData(G4double kinEnergy, 
                                              G4double costheta,
                                              G4double longitudinalPolarization)
{
    switch ((int)longitudinalPolarization)
    {
    case(0): fNumZeroPol++; break;
    case(1): fNumPositivePol++; break;
    case(-1): fNumNegativePol++; break;
    default: break;
    }

  ++fCurrentNumber;
  fSumEnergy+=kinEnergy;
  fSumEnergy2+=kinEnergy*kinEnergy;
  fSumPolarization+=longitudinalPolarization;
  fSumPolarization2+=longitudinalPolarization*longitudinalPolarization;
  fSumCosTheta+=costheta;
  fSumCosTheta2+=costheta*costheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::PrintResults(G4int totalNumberOfEvents)
{
  G4cout << "Total across all events: " << fTotalNumber << "\n";
  G4cout << "Total events: " << totalNumberOfEvents << "\n";

  G4cout<<"Mean Number per Event :"
        <<G4double(fTotalNumber)/G4double(totalNumberOfEvents)<<"\n";
  if (fTotalNumber==0) fTotalNumber=1;
  G4double energyMean=fSumEnergy/fTotalNumber;
  G4double energyRms=std::sqrt(fSumEnergy2/fTotalNumber-energyMean*energyMean);
  G4cout<<"Mean Energy :"<< G4BestUnit(energyMean,"Energy")
        <<" +- "<<G4BestUnit(energyRms,"Energy")<<"\n";
  G4double polarizationMean=fSumPolarization/fTotalNumber;
  G4double polarizationRms=
    std::sqrt(fSumPolarization2/fTotalNumber-polarizationMean*polarizationMean);
  G4cout<<"Mean Polarization :"<< polarizationMean
        <<" +- "<<polarizationRms<<"\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::Clear()
{
  fCurrentNumber=0;
  fTotalNumber=fTotalNumber2=0;
  fSumEnergy=fSumEnergy2=0;
  fSumPolarization=fSumPolarization2=0;
  fSumCosTheta=fSumCosTheta2=0;
  fNumNegativePol = fNumZeroPol = fNumPositivePol = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
