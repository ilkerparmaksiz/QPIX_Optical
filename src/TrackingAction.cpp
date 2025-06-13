// -----------------------------------------------------------------------------
//  TrackingAction.cpp
//
//
//   * Author: Everybody is an author!
//   * Creation date: 4 August 2020
// -----------------------------------------------------------------------------

#include "TrackingAction.h"

// Q-Pix includes
#include "MCParticle.h"
#include "MCTruthManager.h"

// GEANT4 includes
#include "G4TrackingManager.hh"

// C++ includes
#include <iostream>
#include <fstream>
#include "G4OpticalPhoton.hh"
TrackingAction::TrackingAction()
{}

TrackingAction::~TrackingAction()
{}

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
    // get MC truth manager
    MCTruthManager * mc_truth_manager = MCTruthManager::Instance();
    auto PDefi = track->GetDynamicParticle()->GetParticleDefinition();
    /*if(PDefi==G4OpticalPhoton::Definition()){
        std::cout << "photon"<<std::endl;
    }*/
    // create new MCParticle object
    MCParticle * particle = new MCParticle();
    particle->SetTrackID(track->GetTrackID());
    particle->SetParentTrackID(track->GetParentID());
    particle->SetPDGCode(track->GetDefinition()->GetPDGEncoding());
    particle->SetMass(track->GetDynamicParticle()->GetMass());
    particle->SetCharge(track->GetDynamicParticle()->GetCharge());
    particle->SetGlobalTime(track->GetGlobalTime() / CLHEP::ns);
    particle->SetTotalOccupancy(track->GetDynamicParticle()->GetTotalOccupancy());

    particle->SetInitialPosition(
        TLorentzVector(
            track->GetPosition().x() / CLHEP::cm,
            track->GetPosition().y() / CLHEP::cm,
            track->GetPosition().z() / CLHEP::cm,
            track->GetGlobalTime()   / CLHEP::ns
        )
    );


    particle->SetInitialMomentum(
        TLorentzVector(
            track->GetMomentum().x() / CLHEP::MeV,
            track->GetMomentum().y() / CLHEP::MeV,
            track->GetMomentum().z() / CLHEP::MeV,
            track->GetTotalEnergy()  / CLHEP::MeV
        )
    );

    // add track ID to parent MC particle
    // we might need to deal with cases where some particles aren't tracked (?)
    // we can use a try block for that if need be
    if (track->GetParentID() > 0)
    {
        // get parent MC particle
        MCParticle * parent_particle = mc_truth_manager->GetMCParticle(track->GetParentID());
        parent_particle->AddDaughter(track->GetTrackID());
    }

    // add MC particle to MC truth manager
    mc_truth_manager->AddMCParticle(particle);
    //G4cout << "Added MCParticle" << G4endl;
}

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
    MCTruthManager* mc_truth_manager = MCTruthManager::Instance();
    MCParticle* particle = mc_truth_manager->GetMCParticle(track->GetTrackID());
    auto* PDefi = track->GetDynamicParticle()->GetParticleDefinition();

    // Only care about optical photons
    if (PDefi == G4OpticalPhoton::Definition()) {
        // Get the name of the volume where the track ended
        G4String volName = track->GetStep()->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

        // Check if it died in the silicon detector
        if (volName == "Photon_detector.physical") {
            G4ThreeVector pos = track->GetStep()->GetPostStepPoint()->GetPosition();
            G4double time = track->GetStep()->GetPostStepPoint()->GetGlobalTime();

            // Save in your MCParticle if you want
            particle->SetFinalPosition(
                TLorentzVector(
                    pos.x() / CLHEP::cm,
                    pos.y() / CLHEP::cm,
                    pos.z() / CLHEP::cm,
                    time   / CLHEP::ns
                )
            );

            // Print or store
            std::cout << "[Photon Detection] at Silicon volume:" << std::endl;
            std::cout << " → Position (cm): " << pos / CLHEP::cm << std::endl;
            std::cout << " → Time (ns): " << time / CLHEP::ns << std::endl;
        }
    }

    // Set process name
    if (track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep() != nullptr) {
        particle->SetProcess(track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
    } else {
        particle->SetProcess("User Limit");
    }
}
