// -----------------------------------------------------------------------------
//  AnalysisManager.cpp
//
//  Class definition of the analysis manager
//   * Author: Everybody is an author!
//   * Creation date: 4 August 2020
// -----------------------------------------------------------------------------

// Class includes
#include "AnalysisManager.h"

// Qpix includes
#include "AnalysisData.h"
#include "ConfigManager.h"
#include "../cfg/config.h"
// Geant4 includes
#include "G4AutoLock.hh"

// ROOT includes
#include "TFile.h"
#include "TTree.h"


#ifdef With_Opticks
#include "SEvt.hh"
#include "G4CXOpticks.hh"
#include "OpticksPhoton.hh"
#include "OpticksGenstep.h"
#include "QSim.hh"
#endif

namespace {
  G4Mutex bookMutex = G4MUTEX_INITIALIZER;
  G4Mutex saveMutex = G4MUTEX_INITIALIZER;
  G4Mutex fillMutex = G4MUTEX_INITIALIZER;
  G4Mutex instMutex = G4MUTEX_INITIALIZER;
  G4Mutex metaMutex = G4MUTEX_INITIALIZER;
}

AnalysisManager * AnalysisManager::instance_ = 0;

//-----------------------------------------------------------------------------
AnalysisManager::AnalysisManager()
  : tfile_(0), metadata_(0), event_tree_(0), G4_Optical_tree_(0), Opticks_Optical_tree_(0),  G4event_id(0)
{
#ifdef G4ANALYSIS_USE
#endif

}

//-----------------------------------------------------------------------------
AnalysisManager::~AnalysisManager()
{
#ifdef G4ANALYSIS_USE
#endif
}

//-----------------------------------------------------------------------------
AnalysisManager * AnalysisManager::Instance()
{
  G4AutoLock instLock(&instMutex);
  if (instance_ == 0) instance_ = new AnalysisManager();
  return instance_;
}

//-----------------------------------------------------------------------------
void AnalysisManager::Book(const std::string& file_path)
{

  if (!G4Threading::IsMasterThread()) return; // only run Book() for the master thread
  Reset();
  // Check if tfile_ is a null pointer, if so, create ROOT output file
  if (tfile_ == 0) {
    //G4cout << "Testing file_path.data() for null pointer. file_path.data() = " << file_path.data() << G4endl;
    tfile_ = new TFile(file_path.data(), "recreate", "qpix");
  }

  // check if metadata_ is a null pointer, if so, create metadata tree
  if (metadata_ == 0)
  {
    metadata_ = new TTree("metadata", "metadata");

    metadata_->Branch("detector_length_x", &detector_length_x_, "detector_length_x/D");
    metadata_->Branch("detector_length_y", &detector_length_y_, "detector_length_y/D");
    metadata_->Branch("detector_length_z", &detector_length_z_, "detector_length_z/D");
    metadata_->Branch("use_HD_detector_configuration", &useHDDetectorConfiguration_);
  }

  // check if event_tree_ is a null pointer, if so, create event tree
  if (event_tree_ == 0) 
  {
    event_tree_ = new TTree("event_tree", "event tree");

    event_tree_->Branch("run",   &event.run_,   "run/I");
    event_tree_->Branch("event", &event.event_, "event/I");

    event_tree_->Branch("generator_initial_number_particles",  &event.generator_initial_number_particles_, "generator_initial_number_particles/I");
    event_tree_->Branch("generator_initial_particle_x",        &event.generator_initial_particle_x_);
    event_tree_->Branch("generator_initial_particle_y",        &event.generator_initial_particle_y_);
    event_tree_->Branch("generator_initial_particle_z",        &event.generator_initial_particle_z_);
    event_tree_->Branch("generator_initial_particle_t",        &event.generator_initial_particle_t_);
    event_tree_->Branch("generator_initial_particle_px",       &event.generator_initial_particle_px_);
    event_tree_->Branch("generator_initial_particle_py",       &event.generator_initial_particle_py_);
    event_tree_->Branch("generator_initial_particle_pz",       &event.generator_initial_particle_pz_);
    event_tree_->Branch("generator_initial_particle_energy",   &event.generator_initial_particle_energy_);
    event_tree_->Branch("generator_initial_particle_pdg_code", &event.generator_initial_particle_pdg_code_);
    event_tree_->Branch("generator_initial_particle_mass",     &event.generator_initial_particle_mass_);
    event_tree_->Branch("generator_initial_particle_charge",   &event.generator_initial_particle_charge_);

    event_tree_->Branch("generator_intermediate_number_particles",  &event.generator_intermediate_number_particles_, "generator_intermediate_number_particles/I");
    event_tree_->Branch("generator_intermediate_particle_x",        &event.generator_intermediate_particle_x_);
    event_tree_->Branch("generator_intermediate_particle_y",        &event.generator_intermediate_particle_y_);
    event_tree_->Branch("generator_intermediate_particle_z",        &event.generator_intermediate_particle_z_);
    event_tree_->Branch("generator_intermediate_particle_t",        &event.generator_intermediate_particle_t_);
    event_tree_->Branch("generator_intermediate_particle_px",       &event.generator_intermediate_particle_px_);
    event_tree_->Branch("generator_intermediate_particle_py",       &event.generator_intermediate_particle_py_);
    event_tree_->Branch("generator_intermediate_particle_pz",       &event.generator_intermediate_particle_pz_);
    event_tree_->Branch("generator_intermediate_particle_energy",&event.generator_intermediate_particle_energy_);
    event_tree_->Branch("generator_intermediate_particle_pdg_code", &event.generator_intermediate_particle_pdg_code_);
    event_tree_->Branch("generator_intermediate_particle_mass",     &event.generator_intermediate_particle_mass_);
    event_tree_->Branch("generator_intermediate_particle_charge",   &event.generator_intermediate_particle_charge_);

    event_tree_->Branch("generator_final_number_particles",    &event.generator_final_number_particles_,   "generator_final_number_particles/I");
    event_tree_->Branch("generator_final_particle_x",          &event.generator_final_particle_x_);
    event_tree_->Branch("generator_final_particle_y",          &event.generator_final_particle_y_);
    event_tree_->Branch("generator_final_particle_z",          &event.generator_final_particle_z_);
    event_tree_->Branch("generator_final_particle_t",          &event.generator_final_particle_t_);
    event_tree_->Branch("generator_final_particle_px",         &event.generator_final_particle_px_);
    event_tree_->Branch("generator_final_particle_py",         &event.generator_final_particle_py_);
    event_tree_->Branch("generator_final_particle_pz",         &event.generator_final_particle_pz_);
    event_tree_->Branch("generator_final_particle_energy",     &event.generator_final_particle_energy_);
    event_tree_->Branch("generator_final_particle_pdg_code",   &event.generator_final_particle_pdg_code_);
    event_tree_->Branch("generator_final_particle_mass",       &event.generator_final_particle_mass_);
    event_tree_->Branch("generator_final_particle_charge",     &event.generator_final_particle_charge_);

    event_tree_->Branch("number_particles", &event.number_particles_, "number_particles/I");
    event_tree_->Branch("number_hits",      &event.number_hits_,      "number_hits/L");

    event_tree_->Branch("energy_deposit",   &event.energy_deposit_,   "energy_deposit/D");

    event_tree_->Branch("particle_track_id",        &event.particle_track_id_);
    event_tree_->Branch("particle_parent_track_id", &event.particle_parent_track_id_);
    event_tree_->Branch("particle_pdg_code",        &event.particle_pdg_code_);
    event_tree_->Branch("particle_mass",            &event.particle_mass_);
    event_tree_->Branch("particle_charge",          &event.particle_charge_);
    event_tree_->Branch("particle_process_key",     &event.particle_process_key_);
    event_tree_->Branch("particle_total_occupancy", &event.particle_total_occupancy_);
    event_tree_->Branch("particle_initial_x",       &event.particle_initial_x_);
    event_tree_->Branch("particle_initial_y",       &event.particle_initial_y_);
    event_tree_->Branch("particle_initial_z",       &event.particle_initial_z_);
    event_tree_->Branch("particle_initial_t",       &event.particle_initial_t_);
    event_tree_->Branch("particle_final_x",       &event.particle_final_x_);
    event_tree_->Branch("particle_final_y",       &event.particle_final_y_);
    event_tree_->Branch("particle_final_z",       &event.particle_final_z_);
    event_tree_->Branch("particle_final_t",       &event.particle_final_t_);
    event_tree_->Branch("particle_initial_px",      &event.particle_initial_px_);
    event_tree_->Branch("particle_initial_py",      &event.particle_initial_py_);
    event_tree_->Branch("particle_initial_pz",      &event.particle_initial_pz_);
    event_tree_->Branch("particle_initial_energy",  &event.particle_initial_energy_);

    event_tree_->Branch("particle_number_daughters",  &event.particle_number_daughters_);
    event_tree_->Branch("particle_daughter_track_id", &event.particle_daughter_track_ids_);

    event_tree_->Branch("hit_track_id",       &event.hit_track_id_);
    event_tree_->Branch("hit_pdg_code",       &event.hit_pdg_code_);
    event_tree_->Branch("hit_start_x",        &event.hit_start_x_);
    event_tree_->Branch("hit_start_y",        &event.hit_start_y_);
    event_tree_->Branch("hit_start_z",        &event.hit_start_z_);
    event_tree_->Branch("hit_start_t",        &event.hit_start_t_);
    event_tree_->Branch("hit_end_x",          &event.hit_end_x_);
    event_tree_->Branch("hit_end_y",          &event.hit_end_y_);
    event_tree_->Branch("hit_end_z",          &event.hit_end_z_);
    event_tree_->Branch("hit_end_t",          &event.hit_end_t_);
    event_tree_->Branch("hit_energy_deposit", &event.hit_energy_deposit_);
    event_tree_->Branch("hit_length",         &event.hit_length_);

      // Optical Tree
      if(G4_Optical_tree_ == 0){
          G4_Optical_tree_ = new TTree("G4_Photons", "Geant4 Photon Hit Information");
          //G4_Optical_tree_->Branch("run",   &event.run_,   "run/I");
          G4_Optical_tree_->Branch("event_id", &G4event_id);
          G4_Optical_tree_->Branch("photon_hit_x",    &G4_photon_hit_x);
          G4_Optical_tree_->Branch("photon_hit_y",    &G4_photon_hit_y);
          G4_Optical_tree_->Branch("photon_hit_z",    &G4_photon_hit_z);
          G4_Optical_tree_->Branch("photon_hit_t",    &G4_photon_hit_t);
          G4_Optical_tree_->Branch("photon_hit_wavelength",    &G4_photon_hit_wavelength);
      }
#ifdef With_Opticks
      // Opticks Tree
      if(Opticks_Optical_tree_ == 0){
          Opticks_Optical_tree_ = new TTree("Opticks_Photons", "Opticks Photon Hit Information");
          //Opticks_Optical_tree_->Branch("run",   &event.run_,   "run/I");
          Opticks_Optical_tree_->Branch("event_id", &Opticks_event_id);
          Opticks_Optical_tree_->Branch("photon_hit_x",    &Opticks_photon_hit_x);
          Opticks_Optical_tree_->Branch("photon_hit_y",    &Opticks_photon_hit_y);
          Opticks_Optical_tree_->Branch("photon_hit_z",    &Opticks_photon_hit_z);
          Opticks_Optical_tree_->Branch("photon_hit_t",    &Opticks_photon_hit_t);
          Opticks_Optical_tree_->Branch("photon_hit_wavelength",    &Opticks_photon_hit_wavelength);
      }
#endif
  }
}

void AnalysisManager::Reset() {
    G4event_id.clear();
    G4event_id.shrink_to_fit();
    G4_photon_hit_x.clear();
    G4_photon_hit_x.shrink_to_fit();
    G4_photon_hit_y.clear();
    G4_photon_hit_y.shrink_to_fit();
    G4_photon_hit_z.clear();
    G4_photon_hit_z.shrink_to_fit();
    G4_photon_hit_t.clear();
    G4_photon_hit_t.shrink_to_fit();
    G4_photon_hit_wavelength.clear();
    G4_photon_hit_wavelength.shrink_to_fit();

    Opticks_event_id.clear();
    Opticks_event_id.shrink_to_fit();
    Opticks_photon_hit_x.clear();
    Opticks_photon_hit_x.shrink_to_fit();
    Opticks_photon_hit_y.clear();
    Opticks_photon_hit_y.shrink_to_fit();

    Opticks_photon_hit_z.clear();
    Opticks_photon_hit_z.shrink_to_fit();
    Opticks_photon_hit_t.clear();
    Opticks_photon_hit_t.shrink_to_fit();
    Opticks_photon_hit_wavelength.clear();
    Opticks_photon_hit_wavelength.shrink_to_fit();

}

//-----------------------------------------------------------------------------
void AnalysisManager::Save()
{
  if (!G4Threading::IsMasterThread()) {return;}  // only run Save() with the master thread

  //G4AutoLock saveLock(&saveMutex);

  // write TTree objects to file and close file
  tfile_->cd();
  metadata_->Write();
  event_tree_->Write();
  G4_Optical_tree_->Write();
#ifdef With_Opticks
  Opticks_Optical_tree_->Write();
#endif
  tfile_->Close();

}

//-----------------------------------------------------------------------------
void AnalysisManager::EventFill(const AnalysisData& rhs)
{
  G4AutoLock fillLock(&fillMutex);
  // fill TTree objects per event
  event = rhs;
  event_tree_->Fill();
  G4_Optical_tree_->Fill();
#ifdef With_Opticks

  Opticks_Optical_tree_->Fill();
#endif
  Reset();
}

//-----------------------------------------------------------------------------
void AnalysisManager::FillMetadata()
{
  G4AutoLock metaLock(&metaMutex);
  detector_length_x_ = ConfigManager::GetDetectorWidth();
  detector_length_y_ = ConfigManager::GetDetectorHeight();
  detector_length_z_ = ConfigManager::GetDetectorLength();
  useHDDetectorConfiguration_ = ConfigManager::GetUseHDDetectorConfiguration();
  metadata_->Fill();
}
void AnalysisManager::AddG4PhotonHits(G4int eventID, double x ,double y, double z,double t, double wavelength ) {
    //std::cout << "Saving Photon Info" << std::endl;
    G4AutoLock booklock(&bookMutex);

    G4_photon_hit_x.push_back(x);
    G4_photon_hit_y.push_back(y);
    G4_photon_hit_z.push_back(z);
    G4_photon_hit_t.push_back(t);
    G4_photon_hit_wavelength.push_back(wavelength);
    G4event_id.push_back(eventID);


}

void AnalysisManager::AddOPhotonHits(G4int eventID,double x ,double y, double z,double t, double wavelength ) {
    G4AutoLock saveLock(&saveMutex);
#ifdef With_Opticks
    SEvt* sev             = SEvt::Get_EGPU();
    std::vector<sphoton> hits;
    sphoton::Get(hits, sev->getHit());

    for (auto & hit : hits){
        Opticks_event_id.push_back(eventID);
        Opticks_photon_hit_x.push_back(hit.pos.x);
        Opticks_photon_hit_y.push_back(hit.pos.y);
        Opticks_photon_hit_z.push_back(hit.pos.z);
        Opticks_photon_hit_t.push_back(hit.time);
        Opticks_photon_hit_wavelength.push_back(hit.wavelength);
    }
    hits.clear();
    hits.shrink_to_fit();
#endif
}

