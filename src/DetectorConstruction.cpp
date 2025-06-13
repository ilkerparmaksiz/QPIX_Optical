// -----------------------------------------------------------------------------
//  G4_QPIX | DetectorConstruction.cpp
//
//  Definition of detector geometry and materials.
//   * Author: Justo Martin-Albo
//   * Creation date: 14 Aug 2019
// -----------------------------------------------------------------------------

#include "ConfigManager.h"
#include "DetectorConstruction.h"
#include "TrackingSD.h"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "OpticalMaterialProperties.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4OpticalSurface.hh"
#include "../cfg/config.h"
#include "G4LogicalSkinSurface.hh"
#ifdef With_Opticks
#include "G4CXOpticks.hh"
#endif
DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction()
{
}

DetectorConstruction::~DetectorConstruction()
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get Detector Geometry first to dicatate world_size
  G4ThreeVector Photon_detector_dim;

  if(ConfigManager::GetUseHDDetectorConfiguration())
  {
    // DETECTOR HD CONFIGURATION //////////////////////////////////////////////
    // resemble an APA size
    ConfigManager::SetDetectorWidth(2.3*CLHEP::m);   // detector_x
    ConfigManager::SetDetectorHeight(6.0*CLHEP::m);  // detector_y
    ConfigManager::SetDetectorLength(3.6*CLHEP::m);  // detector_z

  } else {

    // DETECTOR VD CONFIGURATION //////////////////////////////////////////////
    ConfigManager::SetDetectorHeight(13.0*CLHEP::m); // detector_y
    ConfigManager::SetDetectorLength(6.5*CLHEP::m);  // detector_z
    ConfigManager::SetDetectorWidth(20.0*CLHEP::m);  // detector_x



  }
  Photon_detector_dim = {ConfigManager::GetDetectorWidth(),ConfigManager::GetDetectorHeight(),15*CLHEP::mm};

  // WORLD /////////////////////////////////////////////////

  G4double world_size = std::max({ConfigManager::GetDetectorHeight(),ConfigManager::GetDetectorLength(),ConfigManager::GetDetectorWidth()})*CLHEP::m;
  G4Material* world_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

  G4Material* Si_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  G4Box* world_solid_vol =
    new G4Box("world.solid", world_size/2., world_size/2., world_size/2.);

  G4LogicalVolume* world_logic_vol =
    new G4LogicalVolume(world_solid_vol, world_mat, "world.logical");
  world_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VPhysicalVolume* world_phys_vol =
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
                      world_logic_vol, "world.physical", 0, false, 0, true);
                      

  std::cout << " Detector configuration is: " << ConfigManager::GetUseHDDetectorConfiguration() << std::endl;

  // DETECTOR
  G4Material* detector_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
  detector_mat->SetMaterialPropertiesTable(OpticalMaterialProperties::OpticalLAr());
  G4Box* detector_solid_vol =
    new G4Box("detector.solid", ConfigManager::GetDetectorWidth()/2., ConfigManager::GetDetectorHeight()/2., ConfigManager::GetDetectorLength()/2.);

  G4LogicalVolume* detector_logic_vol =
            new G4LogicalVolume(detector_solid_vol, detector_mat, "detector.logical");

    G4Box* PhotonDetect_solid_vol = new G4Box("Photon_detector.solid", Photon_detector_dim.x()/2., Photon_detector_dim.y()/2., Photon_detector_dim.z()/2.);

  G4LogicalVolume* Photon_detector_logic_vol =
            new G4LogicalVolume(PhotonDetect_solid_vol, Si_mat, "Photon_detector.logical");
 detector_logic_vol->SetVisAttributes(G4Colour(1,1,1,0.3));
 Photon_detector_logic_vol->SetVisAttributes(G4Color(0,1,0));

  G4ThreeVector offset(ConfigManager::GetDetectorWidth()/2., ConfigManager::GetDetectorHeight()/2., ConfigManager::GetDetectorLength()/2.);

  new G4PVPlacement(0, offset,
                    detector_logic_vol, "detector.physical", world_logic_vol, false, 0, true);

  new G4PVPlacement(0, G4ThreeVector(0,0.,325*CLHEP::cm/2),
                    Photon_detector_logic_vol, "Photon_detector.physical", detector_logic_vol, false, 0, true);
  //////////////////////////////////////////////////////////
G4OpticalSurface * ops= new G4OpticalSurface("SiliconeDetector",unified,polished,dielectric_metal);
ops->SetMaterialPropertiesTable(OpticalMaterialProperties::PerfectDetector());
G4LogicalSkinSurface * lss = new G4LogicalSkinSurface("SiliconeDetector",detector_logic_vol,ops);


#ifdef With_Opticks
    //std::cout <<"Setting our detector geometry with opticks" <<std::endl;
    G4CXOpticks::SetGeometry(world_phys_vol);
    //std::cout << SEventConfig::Desc() <<std::endl;
#endif


  return world_phys_vol;
}

void DetectorConstruction::ConstructSDandField()
{
  // SENSITIVE DETECTOR ////////////////////////////////////

  TrackingSD* tracking_sd = new TrackingSD("/G4QPIX/TRACKING", "TrackingHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(tracking_sd);

  G4LogicalVolume* detector_logic_vol =
    G4LogicalVolumeStore::GetInstance()->GetVolume("detector.logical");

  SetSensitiveDetector(detector_logic_vol, tracking_sd);

  //////////////////////////////////////////////////////////
}
