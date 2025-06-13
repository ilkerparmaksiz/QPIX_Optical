// -----------------------------------------------------------------------------
//  AnalysisManager.h
//
//  Class definition of the analysis manager
//   * Author: Everybody is an author!
//   * Creation date: 4 August 2020
// -----------------------------------------------------------------------------

#ifndef AnalysisManager_h
#define AnalysisManager_h 1

// Qpix includes
#include "AnalysisData.h"

// GEANT4 includes
#include "globals.hh"

// C++ includes
#include <map>
#include <set>

class TFile;
class TTree;
class AnalysisData;

class AnalysisManager {

  public:
    ~AnalysisManager();

    void Book(const std::string&);
    void Save();
    void EventFill(const AnalysisData&);
    void FillMetadata();

    void AddG4PhotonHits(G4int eventID , double x ,double y, double z,double t, double wavelength );
    void AddOPhotonHits();

    static AnalysisManager* Instance();

    AnalysisData event;

  private:

    AnalysisManager();

    static AnalysisManager * instance_;

    // ROOT objects
    TFile * tfile_;
    TTree * metadata_;
    TTree * event_tree_;
    TTree * G4_Optical_tree_;
    TTree * Opticks_Optical_tree_;

    std::vector< int > G4event_id;
    std::vector< double > G4_photon_hit_x;
    std::vector< double > G4_photon_hit_y;
    std::vector< double > G4_photon_hit_z;
    std::vector< double > G4_photon_hit_t;
    std::vector< double > G4_photon_hit_wavelength;

    std::vector< int >  Opticks_event_id;
    std::vector< double > Opticks_photon_hit_x;
    std::vector< double > Opticks_photon_hit_y;
    std::vector< double > Opticks_photon_hit_z;
    std::vector< double > Opticks_photon_hit_t;
    std::vector< double > Opticks_photon_hit_wavelength;

    void AddInitialGeneratorParticle(GeneratorParticle const *);
    void AddFinalGeneratorParticle(GeneratorParticle const *);
    void Reset();

  private:

    std::set< std::string > process_names_;

    // variables that will go into the metadata tree
    double detector_length_x_;
    double detector_length_y_;
    double detector_length_z_;
    bool useHDDetectorConfiguration_;



};

#endif
