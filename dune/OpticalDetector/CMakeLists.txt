#add_subdirectory(BackTrackingTests)

art_make( BASENAME_ONLY ALLOW_UNDERSCORES
                LIBRARY_NAME dune_OpticalDetector
                LIB_LIBRARIES 
                                lardataobj_RecoBase
                                larcorealg_Geometry
                                larcore_Geometry_Geometry_service
                                larana_OpticalDetector_OpHitFinder
                                nusimdata::SimulationBase
                                fhiclcpp::fhiclcpp
                                messagefacility::MF_MessageLogger
                                ART_FRAMEWORK_CORE
                                ART_FRAMEWORK_PRINCIPAL
                                ART_FRAMEWORK_SERVICES_REGISTRY
                                ART_ROOT_IO_TFILE_SUPPORT
                                ROOT::Core
                                ART_ROOT_IO_TFILESERVICE_SERVICE
                                ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
                                art::Persistency_Common canvas
                                art::Utilities canvas
                                cetlib::cetlib cetlib_except 
                                ROOT_BASIC_LIB_LIST
                                
                                Boost::filesystem
  
                MODULE_LIBRARIES
                                dune_DuneObj
                                lardataalg_DetectorInfo
                                lardataobj_RecoBase
                                larana_OpticalDetector
                                larana_OpticalDetector_OpHitFinder
                                larcorealg_Geometry
                                larcore_Geometry_Geometry_service
                                dune_OpticalDetector
                                larsim_MCCheater_ParticleInventoryService_service
                                larsim_MCCheater_PhotonBackTrackerService_service 
                                nusimdata::SimulationBase
                                nurandom::RandomUtils_NuRandomService_service
                                fhiclcpp::fhiclcpp
                                messagefacility::MF_MessageLogger
                                ART_FRAMEWORK_CORE
                                ART_FRAMEWORK_PRINCIPAL
                                ART_FRAMEWORK_SERVICES_REGISTRY
                                ART_ROOT_IO_TFILE_SUPPORT
                                ROOT::Core
                                ART_ROOT_IO_TFILESERVICE_SERVICE
                                ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
                                ART_PERSISTENCY_PROVENANCE
                                art::Persistency_Common canvas
                                art::Utilities canvas
                                cetlib::cetlib cetlib_except
                                CLHEP
                                ROOT_BASIC_LIB_LIST
                                
                                Boost::filesystem

                SERVICE_LIBRARIES
                                larcorealg_Geometry
                                larcore_Geometry_Geometry_service
                                dune_OpticalDetector
                                fhiclcpp::fhiclcpp
                                messagefacility::MF_MessageLogger
                                ART_FRAMEWORK_CORE
                                ART_FRAMEWORK_PRINCIPAL
                                ART_FRAMEWORK_SERVICES_REGISTRY
                                ART_ROOT_IO_TFILE_SUPPORT
                                ROOT::Core
                                ART_ROOT_IO_TFILESERVICE_SERVICE
                                art::Persistency_Common canvas
                                art::Utilities canvas
                                cetlib::cetlib cetlib_except
                                CLHEP
                                ROOT_BASIC_LIB_LIST
                                
                                Boost::filesystem

)

install_headers()
install_fhicl()
install_source()
install_fw( LIST sensl_1.dat
  SUBDIRNAME dune/OpticalDetector )
