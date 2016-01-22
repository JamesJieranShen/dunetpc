// -*- mode: c++; c-basic-offset: 2; -*-
////////////////////////////////////////////////////////////////////////
// Class:       PTBToOffline
// Module Type: producer
// File:        PTBToOffline_module.cc
//
// Repackage PTB data into External Trigger data products
//
// Karl Warburton k.warburton@sheffield.ac.uk
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <map>

//lbne-artdaq includes
#include "lbne-raw-data/Overlays/PennMilliSlice.hh"
#include "lbne-raw-data/Overlays/PennMicroSlice.hh"
#include "artdaq-core/Data/Fragments.hh"

//larsoft includes
#include "RawData/raw.h"
#include "RawData/OpDetWaveform.h"
#include "Geometry/Geometry.h"

//daqinput35t includes

#include "utilities/UnpackFragment.h"
#include "PennToOffline.h"

namespace DAQToOffline {
  class PTBToOffline;
}

class DAQToOffline::PTBToOffline : public art::EDProducer {
public:
  explicit PTBToOffline(fhicl::ParameterSet const & pset);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PTBToOffline(PTBToOffline const &) = delete;
  PTBToOffline(PTBToOffline &&) = delete;
  PTBToOffline & operator = (PTBToOffline const &) = delete;
  PTBToOffline & operator = (PTBToOffline &&) = delete;
  void produce(art::Event & evt) override;
  void reconfigure(fhicl::ParameterSet const& pset);
  void printParameterSet();

private:

  std::string fFragType;
  std::string fRawDataLabel;
  std::string fOutputDataLabel;
  double      fNOvAClockFrequency; //MHz

  std::map<int,int> PTBChannelMap;
  std::string       fPTBChannelMapFile;

  std::pair <std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::CounterWordSize> >,
	     std::pair<lbne::PennMicroSlice::Payload_Header::short_nova_timestamp_t, std::bitset<TypeSizes::TriggerWordSize> > > PrevTimeStampWords;
  //long        first_FirstSample;
  //double      first_TimeStamp;
  //long        first_InternalSample;
  //double      first_InternalTimeStamp;

  //bool fDebug;
  //raw::Compress_t        fCompression;      ///< compression type to use
  //unsigned int           fZeroThreshold;    ///< Zero suppression threshold


};


DAQToOffline::PTBToOffline::PTBToOffline(fhicl::ParameterSet const & pset)
{

  this->reconfigure(pset);

  produces< std::vector<raw::ExternalTrigger> > (fOutputDataLabel);  

  //first_FirstSample = -1;
  //first_TimeStamp = -1;
}

void DAQToOffline::PTBToOffline::reconfigure(fhicl::ParameterSet const& pset){

  fFragType           = pset.get<std::string>("FragType");
  fRawDataLabel       = pset.get<std::string>("RawDataLabel");
  fOutputDataLabel    = pset.get<std::string>("OutputDataLabel");
  fNOvAClockFrequency = pset.get<double>("NOvAClockFrequency"); // in MHz
  fPTBChannelMapFile  = pset.get<std::string>("PTBChannelMapFile","ptb_channel_map_dune35t.txt");

  printParameterSet();

  BuildPTBChannelMap(fPTBChannelMapFile, PTBChannelMap);
}

void DAQToOffline::PTBToOffline::printParameterSet(){

  mf::LogDebug("PTBToOffline") << "===================================="   << "\n"
			       << "Parameter Set"                          << "\n"
			       << "===================================="   << "\n"
			       << "fFragType:        " << fFragType        << "\n"
			       << "fRawDataLabel:    " << fRawDataLabel    << "\n"
			       << "fOutputDataLabel: " << fOutputDataLabel << "\n"
			       << "===================================="   << "\n";
}


void DAQToOffline::PTBToOffline::produce(art::Event & evt)
{

  art::Handle<artdaq::Fragments> rawFragments;
  evt.getByLabel(fRawDataLabel, fFragType, rawFragments);

  // Check if there is PTB data in this event
  // Don't crash code if not present, just don't save anything
  try { rawFragments->size(); }
  catch(std::exception e) {
    mf::LogWarning("PTBToOffline") << "WARNING: Raw PTB data not found in event " << evt.event();
    std::vector<raw::ExternalTrigger> Triggers;
    evt.put(std::make_unique<std::vector<raw::ExternalTrigger>>(std::move(Triggers)), fOutputDataLabel);
    return;
  }

  // Check that the data is valid
  if(!rawFragments.isValid()){
    mf::LogError("PTBToOffline") << "Run: " << evt.run()
				 << ", SubRun: " << evt.subRun()
				 << ", Event: " << evt.event()
				 << " is NOT VALID";
    throw cet::exception("raw NOT VALID");
    return;
  }

  auto triggers = PennFragmentToExternalTrigger(*rawFragments, PTBChannelMap, PrevTimeStampWords);

  evt.put(std::make_unique<decltype(triggers)>(std::move(triggers)), fOutputDataLabel);


}

DEFINE_ART_MODULE(DAQToOffline::PTBToOffline)


