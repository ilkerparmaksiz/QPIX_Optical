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
    // get MC truth manager
    MCTruthManager * mc_truth_manager = MCTruthManager::Instance();

    // get MC particle
    MCParticle * particle = mc_truth_manager->GetMCParticle(track->GetTrackID());
    auto PDefi = track->GetDynamicParticle()->GetParticleDefinition();
    if(PDefi==G4OpticalPhoton::Definition()){
        //std::cout << "photon"<<std::endl;
        particle->SetFinalPosition(
                TLorentzVector(
                        track->GetStep()->GetPostStepPoint()->GetPosition().x() / CLHEP::cm,
                        track->GetStep()->GetPostStepPoint()->GetPosition().y()  / CLHEP::cm,
                        track->GetStep()->GetPostStepPoint()->GetPosition().z()  / CLHEP::cm,
                        track->GetStep()->GetPostStepPoint()->GetGlobalTime()   / CLHEP::ns
                )
        );
    }

    // set process
    if (track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep() != 0) {
      particle->SetProcess(track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
    } else {
      particle->SetProcess("User Limit");
    }
}

