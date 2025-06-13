// -----------------------------------------------------------------------------
//  G4_QPIX | EventAction.cpp
//
//
//   * Author: Everybody is an author!
//   * Creation date: 15 Aug 2019
// -----------------------------------------------------------------------------

#include "EventAction.h"

// Q-Pix includes
#include "AnalysisData.h"
#include "AnalysisManager.h"
#include "ConfigManager.h"
#include "MCTruthManager.h"

// GEANT4 includes
#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "../cfg/config.h"
#ifdef With_Opticks
#  include "SEvt.hh"
    #  include "G4CXOpticks.hh"
    #  include "QSim.hh"
namespace {G4Mutex opticks_mt =G4MUTEX_INITIALIZER;}
#endif

EventAction::EventAction():
  G4UserEventAction()
{
}


EventAction::~EventAction()
{
}


void EventAction::BeginOfEventAction(const G4Event*)
{
}


void EventAction::EndOfEventAction(const G4Event* g4event)
{
    event_id_offset_ = ConfigManager::GetEventIDOffset();
    energy_threshold_ = ConfigManager::GetEnergyThreshold();

    // get MC truth manager
    MCTruthManager * mc_truth_manager = MCTruthManager::Instance();

    // get map of particles from MC truth manager
    const MCTruthManager::MCParticleMap& particleMap = mc_truth_manager->GetMCParticleMap();


    double energy_deposited = 0.;

    // add particle to analysis manager
    for (auto const& p : particleMap)
    {
        const MCParticle* particle = p.second;
        energy_deposited += particle->EnergyDeposited();
        // std::cout << "Energy deposited by particle PDG (" << particle->PDGCode() << "): " << particle->EnergyDeposited() << std::endl;
    }

    // don't save event if total energy deposited is below the energy threshold
    if (energy_deposited < energy_threshold_)
    {
        // reset event variables
        event.EventReset();

        // reset event in MC truth manager
        mc_truth_manager->EventReset();

        return;
    }


    // set event number
    event.SetEvent(g4event->GetEventID() + event_id_offset_);

    // add initial generator particles to analysis manager
    for (auto const& particle : mc_truth_manager->GetInitialGeneratorParticles())
    {
        event.AddInitialGeneratorParticle(particle);
    }

    // add initial generator particles to analysis manager
    for (auto const& particle : mc_truth_manager->GetIntermediateGeneratorParticles())
    {
      event.AddIntermediateGeneratorParticle(particle);
    }

    // add final generator particles to analysis manager
    for (auto const& particle : mc_truth_manager->GetFinalGeneratorParticles())
    {
        event.AddFinalGeneratorParticle(particle);
    }

    // add particle to analysis manager
    for (auto const& p : particleMap)
    {
        const MCParticle* particle = p.second;

        event.AddMCParticle(particle);
    }

    // write event to ROOT file and reset event variables
    AnalysisManager * analysisManager = AnalysisManager::Instance();
    analysisManager->EventFill(event);

#ifdef With_Opticks

        G4AutoLock lock(&opticks_mt);
        G4CXOpticks * g4cx=G4CXOpticks::Get();

        G4int eventID=g4event->GetEventID();
        G4int ngenstep=SEvt::GetNumGenstepFromGenstep(0);
        G4int nphotons=SEvt::GetNumPhotonCollected(0);
        G4int hits;



        // Simulate the photons
        if(ngenstep>0){
            //std::cout<<g4cx->desc()<<std::endl;
            //std::cout<<"--- G4Optickx ---" << g4cx->descSimulate() <<std::endl;
            //std::cout<< "Simulating Photons " <<ngenstep <<std::endl;
            g4cx->simulate(eventID,0); // For Simulation
            cudaDeviceSynchronize();

            hits=SEvt::GetNumHit(0);
            //std::cout << "DefaultEventAction Hits " << hits<<std::endl;
            if(hits>0) analysisManager->AddOPhotonHits();
            std::cout<<"Event " <<eventID <<" Simulating with Opticks nphotons "<< nphotons << " nsteps " << ngenstep << " Hits " <<SEvt::GetNumHit(0) << std::endl;
            G4CXOpticks::Get()->reset(eventID);

        }


        //G4cout<<" Opticks End of Event Action" <<G4endl;

#endif


    event.EventReset();

    // reset event in MC truth manager
    mc_truth_manager->EventReset();
}

