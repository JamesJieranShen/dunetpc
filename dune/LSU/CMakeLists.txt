
art_make(
  LIBRARY_NAME DUNELSU
  LIB_LIBRARIES
  ${ROOT_BASIC_LIB_LIST}
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${ART_PERSISTENCY_COMMON}
  canvas
  art_Persistency_Provenance
  ${ART_UTILITIES}
  ${MF_MESSAGELOGGER}
  cetlib cetlib_except
  #larpandoraobj
  larpandora_LArPandoraInterface
  lardataobj_RecoBase
  lardata_Utilities
  larcorealg_Geometry
  nusimdata_SimulationBase
  #larsim_MCCheater_PhotonBackTrackerService_service           
  larsim_MCCheater_BackTrackerService_service           
  larsim_MCCheater_ParticleInventoryService_service
  larreco_RecoAlg
  dune_DuneObj
  ProtoDUNEAnaUtils
  MODULE_LIBRARIES
  DUNELSU
  # For CalorimetryHack
  larreco_Calorimetry
  lardataobj_AnalysisBase
  lardata_ArtDataHelper
  )
install_headers()
install_fhicl()
install_source()
