// -----------------------------------------------------------------------------
//  G4_QPIX | SteppingAction.cpp
//
//  Definition of detector geometry and materials.
//   * Author: Justo Martin-Albo
//   * Creation date: 14 Aug 2019
// -----------------------------------------------------------------------------

#include "SteppingAction.h"

#include "AnalysisData.h"
#include "AnalysisManager.h"
#include "G4VProcess.hh"
#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpBoundaryProcess.hh>
#include "G4RunManager.hh"
#include "G4Scintillation.hh"


SteppingAction::SteppingAction(): G4UserSteppingAction()
{

}


SteppingAction::~SteppingAction()
{
}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
    G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();
    G4Track* track = step->GetTrack();
    AnalysisManager::Instance();

    if (step->GetPostStepPoint()->GetProcessDefinedStep() != 0){
      event.AddProcess(step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
    } else {
      event.AddProcess("User Limit");
    }

    static G4OpBoundaryProcess* boundary = 0;

    if (!boundary) { // the pointer is not defined yet
        // Get the list of processes defined for the optical photon
        // and loop through it to find the optical boundary process.
        G4ProcessVector* pv = pdef->GetProcessManager()->GetProcessList();
        for (size_t i=0; i<pv->size(); i++) {
            if ((*pv)[i]->GetProcessName() == "OpBoundary") {
                boundary = (G4OpBoundaryProcess*) (*pv)[i];
                break;
            }
        }
    }


    if (pdef != G4OpticalPhoton::Definition()) {
        G4SteppingManager *sMg = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
        G4StepStatus stepStatus = sMg->GetfStepStatus();
        if (stepStatus != fAtRestDoItProc) {
            G4ProcessVector *PostStepProc = sMg->GetfPostStepDoItVector();
            size_t MaxSteps = sMg->GetMAXofPostStepLoops();
            for (int stp = 0; stp < MaxSteps; stp++) {
                if ((*PostStepProc)[stp]->GetProcessName() == "Scintillation") {
                    G4Scintillation *ScintProc = (G4Scintillation *) (*PostStepProc)[stp];
                    G4int num_photons = ScintProc->GetNumPhotons();
                    //std::cout << "Scintilation "<< num_photons <<std::endl;

                    if (num_photons > 0) {
                        G4MaterialPropertiesTable *MPT = track->GetMaterial()->GetMaterialPropertiesTable();
                        G4double t1, t2 = 0;
                        G4int singlets, triplets = 0;
                        t1 = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
                        t2 = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT2);
                        singlets = floor(MPT->GetConstProperty(kSCINTILLATIONYIELD1) * num_photons);
                        triplets = ceil(MPT->GetConstProperty(kSCINTILLATIONYIELD2) * num_photons);


                        //std::cout << "Scintilation "<< num_photons <<" Amount of Singlets " <<singlets <<" Triplets " << triplets <<std::endl;
#ifdef With_Opticks
                        if (singlets > 0)
                            U4::CollectGenstep_DsG4Scintillation_r4695(track, step, singlets, 0, t1);
                        if (triplets > 0)
                            U4::CollectGenstep_DsG4Scintillation_r4695(track, step, triplets, 1, t2);
#endif

                    }

                }
            }
        }
    }

    if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        if (boundary->GetStatus() == Detection ){
            G4String detector_name = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
            // Only care about optical photons
                auto run= G4RunManager::GetRunManager();
                G4int eventID=run->GetCurrentEvent()->GetEventID();

                // Get the name of the volume where the track ended
                G4String volName = track->GetStep()->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

                // Check if it died in the silicon detector
                G4ThreeVector pos = track->GetStep()->GetPostStepPoint()->GetPosition();
                G4double time = track->GetStep()->GetPostStepPoint()->GetGlobalTime();
                G4double wavelength=1239.8/step->GetPostStepPoint()->GetTotalEnergy()*CLHEP::eV; //nm
                AnalysisManager * analysisManager = AnalysisManager::Instance();
                analysisManager->AddG4PhotonHits(eventID , pos.x() / CLHEP::cm, pos.y() / CLHEP::cm,pos.z() / CLHEP::cm,time / CLHEP::ns,wavelength);


        }
    }

}
