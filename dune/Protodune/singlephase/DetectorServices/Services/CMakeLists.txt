cet_build_plugin(DetectorPropertiesServiceProtoDUNEsp art::service 
	          ProtoDUNEspDataProviders
              lardataalg_DetectorInfo
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              art::Framework_Principal
              ifdh::ifdh
	      ART_ROOT_IO_ROOTDB
	      ART_FRAMEWORK_SERVICES_REGISTRY
	      art::Utilities
	      canvas::canvas
	      SQLITE3
	      messagefacility::MF_MessageLogger
	      cetlib::cetlib cetlib_except
	      ROOT_BASIC_LIB_LIST
              
)

install_headers()
install_fhicl()
install_source()
