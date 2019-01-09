////////////////////////////////////////////////////////////////////////
// Name: TFileMetadataDUNE.h
//
// A struct datatype to hold the metadata information as it is extracted
// from various places.
//
// Created: 1-Nov-2017,  T. Junk
//   Based on the MicroBooNE version by S. Gollapinni
//
////////////////////////////////////////////////////////////////////////
#ifndef TFILEMETADATADUNE_H
#define TFILEMETADATADUNE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include <iostream>
#include <fstream>

using namespace std;

namespace util{

  class TFileMetadataDUNE
  {
  public:
    TFileMetadataDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~TFileMetadataDUNE();
      
    void reconfigure(fhicl::ParameterSet const& p);
    
    struct metadata {
      std::tuple<std::string, std::string, std::string> fapplication;
      //no crc information yet
      //std::vector<std::string> fcrcs;
      std::string fdata_tier;
      time_t fend_time;
      unsigned int fevent_count=0;
      std::string ffile_format;
      std::string ffile_type;
      art::EventNumber_t ffirst_event=0;
      std::string fgroup; 
      art::EventNumber_t flast_event=0;
      std::set<std::string> fParents;
      std::vector<std::tuple<art::RunNumber_t,art::SubRunNumber_t,std::string>> fruns;
      time_t fstart_time=0; 
      std::string fMCGenerators;			  
      std::string fMCOscillationP;      
      std::string fMCTriggerListVersion;		  
      std::string fMCBeamEnergy;			  
      std::string fMCBeamFluxID;			  
      std::string fMCName;	            
      std::string fMCDetectorType;      
      std::string fMCNeutrinoFlavors;	
      std::string fMCMassHierarchy;     
      std::string fMCMiscellaneous;     
      std::string fMCGeometryVersion;	
      std::string fMCOverlay;		
      std::string fDataRunMode;			  
      std::string fDataDetectorType;		
      std::string fDataName;			
      std::string fStageName;           
    };
    
    metadata md;
    std::set<art::SubRunID> fSubRunNumbers;
    std::map<std::string,std::string> mdmap;    

  private:
   
    // Callbacks.
    void postBeginJob();
    void postOpenFile(std::string const& fn);
    void postEvent(art::Event const& ev, art::ScheduleContext);
    void postBeginSubRun(art::SubRun const& subrun);
    void postEndJob();

    // Private member functions.

    // Data members.

    // Fcl parameters.
    bool fGenerateTFileMetadata;  
    std::string frunType;                     
    std::string fJSONFileName;

  }; // class TFileMetadataDUNE

} //namespace utils
	
DECLARE_ART_SERVICE(util::TFileMetadataDUNE, LEGACY)
	
#endif
