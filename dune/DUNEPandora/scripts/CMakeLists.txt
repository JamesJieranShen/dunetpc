
# The xml files here are used to initialize pandora
# They need to be found via FW_SEARCH_PATH

# install xml files
file(GLOB xml_file *.xml )
install(FILES ${xml_file} DESTINATION ${dunetpc_SCRIPTS_DIR})

