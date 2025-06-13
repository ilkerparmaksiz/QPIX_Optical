// -----------------------------------------------------------------------------
//  G4Basic | RunAction.h
//
//
//   * Author: Everybody is an author!
//   * Creation date: 15 Aug 2019
// -----------------------------------------------------------------------------

#include "RunAction.h"

// Q-Pix includes
#include "AnalysisManager.h"
#include "ConfigManager.h"
#include "MARLEYManager.h"
#include "MCTruthManager.h"
#include "DetectorConstruction.h"
#include "GENIEManager.h"

// GEANT4 includes
#include "G4Box.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"

// C++ includes
#include <filesystem>


RunAction::RunAction()
  : G4UserRunAction()
{
  ConfigManager::Instance();
}


RunAction::~RunAction()
{
}


void RunAction::BeginOfRunAction(const G4Run* g4run)
{
    ConfigManager::Instance();
    //ConfigManager::Print();

    inputFile_ = ConfigManager::GetInputFile();
    outputFile_ = ConfigManager::GetOutputFile();
    marleyJson_ = ConfigManager::GetMarleyJson();
    generator_ = ConfigManager::GetGenerator();
    genieFormat_ = ConfigManager::GetGenieFormat();
    multirun_ = ConfigManager::GetMultirun();
    particleType_ = ConfigManager::GetParticleType();


     G4StrUtil::to_lower(particleType_);
     G4StrUtil::to_lower(generator_);
     G4StrUtil::to_lower(genieFormat_);
    //ConfigManager::Print();

    if (generator_ == "marley") {

        // get MARLEY manager
        MARLEYManager * marleyManager = MARLEYManager::Instance();
    
    
        // configure and create MARLEY generator
        marleyManager->Initialize();
    
        if (particleType_ == "supernova")
        {
        }

    } else if (generator_ == "genie")
    {

        // get ROOT manager
        GENIEManager *genieManager=GENIEManager::Instance();
        if (genieManager->Initialize())
        {
            G4int NumberEventsInTheFile=(G4int)genieManager->GetNEntries();
            G4cout << "nEntries = " << NumberEventsInTheFile << G4endl;
        } else
        {
            G4Exception("[RunAction]","[BeginOfRunAction]",G4ExceptionSeverity::FatalException ,"RootManager is not properly initialized. Check to see if the following file exist /Inputs/ReadFrom_Root_Path in macros/ROOTRead.macro   ") ;

        }
    }

    G4String root_output_path = outputFile_;
    if (multirun_)
    {

        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << g4run->GetRunID();
        std::string runStr_(ss.str());
        G4String parentPath_ = outputFile_.substr(0, outputFile_.find_last_of('/'));
        G4String baseName_ = outputFile_.substr(outputFile_.find_last_of('/')+1, outputFile_.length());
        G4String stem_ = baseName_.substr(0, baseName_.find_last_of('.'));
        G4String extension_ = baseName_.substr(baseName_.find_last_of('.'), baseName_.length());
        
        root_output_path = parentPath_;
        root_output_path += "/";
        root_output_path += stem_;
        root_output_path += "_";
        root_output_path += runStr_;
        root_output_path += extension_;

    }

    // get run number
    AnalysisManager * analysis_manager = AnalysisManager::Instance();

    analysis_manager->Book(root_output_path);
    event.SetRun(g4run->GetRunID());

    // reset event variables
    event.EventReset();

    // get MC truth manager
    MCTruthManager * mc_truth_manager = MCTruthManager::Instance();

    // reset event in MC truth manager
    mc_truth_manager->EventReset();
}


void RunAction::EndOfRunAction(const G4Run*)
{
    // instantiate classes to allow for messengers
    AnalysisManager * analysis_manager = AnalysisManager::Instance();
    GENIEManager::Instance();


    inputFile_ = ConfigManager::GetInputFile();
    outputFile_ = ConfigManager::GetOutputFile();
    marleyJson_ = ConfigManager::GetMarleyJson();
    generator_ = ConfigManager::GetGenerator();
    genieFormat_ = ConfigManager::GetGenieFormat();
    multirun_ = ConfigManager::GetMultirun();
    particleType_ = ConfigManager::GetParticleType();

    // save detector dimensions as metadata
    if (G4Threading::IsMasterThread()){
      analysis_manager->FillMetadata();
    }
   

    // save run to ROOT file
    analysis_manager->Save();
}

