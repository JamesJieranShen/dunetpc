# dunetpc/dune/ArtSupport/Tool/CMakeLists.txt

# David Adams
# August 2017

# Add plugin for each tool.

cet_build_plugin(TestTool   art::tool
                art::Utilities canvas
                cetlib::cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )
