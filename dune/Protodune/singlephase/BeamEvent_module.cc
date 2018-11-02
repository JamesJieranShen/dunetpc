////////////////////////////////////////////////////////////////////////
// Class:       BeamEvent
// Plugin Type: producer (art v2_08_03)
// File:        BeamEvent_module.cc
//
// Generated at Thu Nov  2 22:57:41 2017 by Jonathan Paley using cetskelgen
// from cetlib version v3_01_01.
//
// Edited by Jake Calcutt
////////////////////////////////////////////////////////////////////////


/*///////////////////////////////////////////////////////////////////////
 To-Do: 
       * Now that we're setting the epsilon of the beamfolder,
         we can get rid of my hacky 'FetchWithRetries' functions

       * Need to save all momenta combinations to the beamevt

       * Need to save the current of the magnet to the beamevt
         so we can do any offline fixes to the momentum reconstruction
         !!! Requires change to the ProtoDUNEBeamEvent class !!!

       * Make print statements prettier

       * Fix the tracking portion with real values for monitor positions

       * Make TOF matching more robust. Perhaps search 'downward'
       
*/////////////////////////////////////////////////////////////////////////





#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "IFBeam_service.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "dune/DuneObj/ProtoDUNEBeamSpill.h"
#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "dune/DuneObj/ProtoDUNETimeStamp.h"
#include <bitset>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <limits>

#include "TTree.h"
#include "TH2F.h"
#include "TVectorD.h"
#include "TPolyMarker.h"

namespace proto {
  class BeamEvent;
}

enum tofChan{
  k1A2A,
  k1B2A,
  k1A2B,
  k1B2B
};

typedef std::numeric_limits< double > dbl;

class proto::BeamEvent : public art::EDProducer {
public:
  explicit BeamEvent(fhicl::ParameterSet const & p);
  //  virtual ~BeamEvent();

  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamEvent(BeamEvent const &) = delete;
  BeamEvent(BeamEvent &&) = delete;
  BeamEvent & operator = (BeamEvent const &) = delete;
  BeamEvent & operator = (BeamEvent &&) = delete;

  // Required functions.
  void reconfigure(fhicl::ParameterSet const & p);
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;


  //std::bitset<sizeof(double)*CHAR_BIT> toBinary(const long num);  
  uint64_t joinHighLow(double,double);

  TVector3 ConvertProfCoordinates(double x, double y, double z, double zOffset);
  void BeamMonitorBasisVectors(); 
  bool rotated = false;
  void RotateMonitorVector(TVector3 &vec); 

  uint64_t GetRawDecoderInfo(art::Event &);
  void TimeIn(art::Event &, uint64_t);
  void GetSpillInfo(art::Event &);
  void MatchBeamToTPC();
  void MatchS11ToGen();
  void SetBeamEvent(); 
  void SetCKovInfo(); 

  void MakeTrack(size_t);
  void MomentumSpec(size_t);
  double MomentumCosTheta(double, double, double);

  double GetPosition(std::string, int);
  TVector3 ProjectToTPC(TVector3, TVector3);
  double GetPairedPosition(std::string, size_t);
 
  void  InitXBPFInfo(beam::ProtoDUNEBeamSpill *);
  void  parseGeneralXBPF(std::string, uint64_t, size_t);
  void  parseXBPF(uint64_t);

  void  parseXTOF(uint64_t);
  void  parseXCET(uint64_t);


  void getS11Info(uint64_t);

  //template<class T> 
  //T FetchWithRetries(long long, std::string, int);
  //
  //template<class T> 
  //T FetchWithRetriesDown(long long, std::string, int);

  std::vector<double> FetchWithRetries(long long, std::string, int, long long & );
  std::vector<double> FetchWithRetriesDown(long long, std::string, int, long long &);

  std::vector<double> FetchAndReport(long long, std::string);
   
private:
  
  TTree * fOutTree;
  TH1F  * fFullMomentum;
  TH1F  * fCutMomentum;

  TTree * fGenTrigTree;
  TTree * fXTOF1ATree;
  TTree * fXTOF1BTree;
  TTree * fXTOF2ATree;
  TTree * fXTOF2BTree;


  double fGenTrigFrac;
  double fGenTrigSec;
  double fGenTrigCoarse; 
  double fXTOF1AFrac;
  double fXTOF1ACoarse; 
  double fXTOF1BFrac;
  double fXTOF1BCoarse; 
  double fXTOF2AFrac;
  double fXTOF2ACoarse; 
  double fXTOF2BFrac;
  double fXTOF2BCoarse; 

  double fXTOF1ASec;
  double fXTOF1BSec;
  double fXTOF2ASec;
  double fXTOF2BSec;
  std::vector<double>* diff2A;
  std::vector<double>* diff2B;
 

  long long int eventTime;
  double SpillStart;
  ULong_t SpillStart_alt;
  bool SpillStartValid;
  double PrevStart;
  double SpillEnd;
  double SpillOffset;
  double ActiveTriggerTime;
  long long RDTSTime;
  double RDTSTimeSec;
  double PrevRDTSTimeSec;  
  double RDTSTimeNano; 

  long long valid_fetch_time; 
  long long spill_valid_fetch_time;

  double s11Sec, s11Nano;

  int RDTSTrigger;

  double acqTime;
  double acqStampMBPL;
  int HLTWord;
  long long int HLTTS;
  int BeamOn;
  int BITrigger;
  int Upstream;
  int C1;
  int C2;
  int BP1;
  int BP2;
  int BP3;
  int BP4;


  int eventNum;
  int runNum;
  int subRunNum;
  double CKov1Pressure;
  double CKov2Pressure;
  double CKov1Efficiency;
  double CKov2Efficiency;
  
  TVector3 fBMBasisX;
  TVector3 fBMBasisY;
  TVector3 fBMBasisZ;

  // Declare member data here.
  double  fTimeWindow;
  std::string fBundleName;
  std::string fOutputLabel;
  std::string fURLStr;
  double fBFEpsilon;
  int fIFBeamDebug;
  uint64_t fFixedTime;
  //std::vector< uint64_t > fMultipleTimes;

  std::string firstUpstreamName;
  std::string secondUpstreamName;
  std::string firstDownstreamName;
  std::string secondDownstreamName;

  std::string firstBPROF1;
  std::string secondBPROF1;
  std::string BPROF2;
  std::string BPROF3;
  double      fBeamBend;

  std::vector< std::string > fDevices;
  std::map<std::string, std::string > fDeviceTypes;
  std::map< std::string, double > fFiberDimension;

  int fNRetries;

  // Names of the CERN Beam Devices
  // Names have the form of a prefix + device
  
  // Prefixes for different devices classes:
  // Beam Positions (BPF)
  // Time of flights (TOF)
  // and Cerenkov counters (CET)
  std::string fXBPFPrefix;
  std::string fXTOFPrefix;
  std::string fXCETPrefix;

  // Device Names for each device

  // Time of Flights
  std::string fTOF1;
  std::string fTOF1A, fTOF1B;
  std::string fTOF2;
  std::string fTOF2A, fTOF2B;

  // Cerenkovs
  std::string fCKov1;
  std::string fCKov2;

  double fRotateMonitorXZ;
  double fRotateMonitorYZ;

  double fFirstTrackingProfZ;
  double fSecondTrackingProfZ;
  double fNP04FrontZ;  
  double fBeamX, fBeamY, fBeamZ;

  bool   fForceNewFetch;
  bool   fMatchTime;
  bool   fForceRead;

  bool   fSaveOutTree;
  bool   fDebugTOFs;
  bool   fDebugMomentum;

  bool   gotGeneralTrigger;
  bool   gotTOFs;
  bool   gotCurrent;

  double fTimingCalibration;
  double fCalibrationTolerance;
  double fOffsetTAI;
  int    fFetchOffset;

  double fS11DiffUpper; 
  double fS11DiffLower; 
  double fRDTSToS11Upper; 
  double fRDTSToS11Lower; 

  int    fOffsetCTBtoRDTS;
  int    fToleranceCTBtoRDTS;

  beam::ProtoDUNEBeamEvent * beamevt;
  beam::ProtoDUNEBeamEvent prev_beamevt;

  beam::ProtoDUNEBeamSpill * beamspill;
  beam::ProtoDUNEBeamSpill prev_beamspill;

  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifb;

  art::Handle< std::vector<raw::RDTimeStamp> > RDTimeStampHandle;

  uint64_t validTimeStamp;

  double L1=1.980, L2=1.69472, L3=2.11666;
  double magnetLen, magnetField;
  std::vector<double> current;

  //Hardware Parameters for magnetic field stuff
  double mag_P1 = 5.82044830e-3;
  double mag_P2 = 0.;
  double mag_P3 = -4.68880000e-6;
  double mag_P4 = 324.573967;

};

// Constructor
proto::BeamEvent::BeamEvent(fhicl::ParameterSet const & p)
{
  // Declare products this module will provide
  produces<std::vector<beam::ProtoDUNEBeamEvent>>();  

  // Configure/Reconfigure
  this->reconfigure(p);
}
// END Constructor
////////////////////////

std::vector<double> proto::BeamEvent::FetchAndReport(long long time, std::string name){

  //Note! Sometimes this won't retrieve the data from the database. We'll need
  //      catch the exception outside of this and handle it according to whichever
  //      device we're trying to retrieve from.

  std::vector<double> theResult;
  LOG_INFO("BeamEvent") << "Trying to grab from folder: " << name << "\n";
  LOG_INFO("BeamEvent") << "At Time: " << time << "\n";    

  theResult = bfp->GetNamedVector(time, name);

  LOG_INFO("BeamEvent") << "Successfully fetched " << time << "\n";

  return theResult;
}

////////////////////////
// Fetch Method
//template <class T> 
//T proto::BeamEvent::FetchWithRetries(long long time, std::string name, int nRetry){
std::vector<double> proto::BeamEvent::FetchWithRetries(long long time, std::string name, int nRetry, long long & returned_fetch_time){
  std::vector<double> theResult;
  
  LOG_INFO("BeamEvent") << "\n";
  long long newTime;
  //Search at and above time given with nRetries
  //Will later search below, just in case the event time is actually greater than
  for(newTime = time; newTime < time + nRetry; ++newTime){
    LOG_INFO("BeamEvent") << "Trying to grab from folder: " << name << "\n";
    LOG_INFO("BeamEvent") << "At Time: " << newTime << "\n";    
      try{
        theResult = bfp->GetNamedVector(newTime, name);
        LOG_INFO("BeamEvent") << "Successfully fetched " << newTime << "\n";
        returned_fetch_time = newTime;
        return theResult;
      }
    catch(std::exception e){
      LOG_INFO("BeamEvent") << "Could not fetch with time " << newTime << "\n";      
    }
  }

  //Try a final time. Let it crash if it doesn't work
  LOG_INFO("BeamEvent") << "Trying a final time to grab from folder: " << name << "\n";
  LOG_INFO("BeamEvent") << "At time: " << newTime << "\n";
  theResult = bfp->GetNamedVector(newTime, name);
    LOG_INFO("BeamEvent") << "Successfully fetched" << "\n";
  returned_fetch_time = newTime;
  LOG_INFO("BeamEvent") << "\n";
  return theResult; 
}
// END FetchWithRetries
////////////////////////

////////////////////////
// Fetch Method
//template <class T> 
//T proto::BeamEvent::FetchWithRetriesDown(long long time, std::string name, int nRetry){
std::vector<double> proto::BeamEvent::FetchWithRetriesDown(long long time, std::string name, int nRetry, long long & returned_fetch_time){
  std::vector<double> theResult;
  
  LOG_INFO("BeamEvent") << "\n";
  long long newTime;

  //Now search below
  for(newTime = time; newTime > time - nRetry - 1; --newTime){
    LOG_INFO("BeamEvent") << "Trying to grab from folder: " << name << "\n";
    LOG_INFO("BeamEvent") << "At Time: " << newTime << "\n";    
    try{
      theResult = bfp->GetNamedVector(newTime, name);
      LOG_INFO("BeamEvent") << "Successfully fetched " << newTime << "\n";
      returned_fetch_time = newTime;
      return theResult;
    }
    catch(std::exception e){
//      LOG_INFO("BeamEvent") << "Could not fetch with time " << newTime << "\n";      
    }
  }
  
  //Try a final time. Let it crash if it doesn't work
  LOG_INFO("BeamEvent") << "Trying a final time to grab from folder: " << name << "\n";
  LOG_INFO("BeamEvent") << "At time: " << newTime << "\n";
  theResult = bfp->GetNamedVector(newTime, name);
  returned_fetch_time = newTime;
  LOG_INFO("BeamEvent") << "Successfully fetched" << "\n";
  LOG_INFO("BeamEvent") << "\n";
  return theResult; 
}
// END FetchWithRetries
////////////////////////


//Gets the Timing and CTB raw-decoded info.
//Finds the triggers, and looks for a valid trigger word
//(i.e. coming from beam)
//
//Returns the timestamp of the high level trigger.
uint64_t proto::BeamEvent::GetRawDecoderInfo(art::Event & e){
  LOG_INFO("BeamEvent") << "\n";
  LOG_INFO("BeamEvent") << "Getting Raw Decoder Info" << "\n";
  e.getByLabel("timingrawdecoder","daq",RDTimeStampHandle);
  LOG_INFO("BeamEvent") << "RDTS valid? " << RDTimeStampHandle.isValid() << "\n";
  for (auto const & RDTS : *RDTimeStampHandle){
    LOG_INFO("BeamEvent") << "High: " << RDTS.GetTimeStamp_High() << "\n";
    LOG_INFO("BeamEvent") << "Low: " << RDTS.GetTimeStamp_Low() << "\n"; 

    uint64_t high = RDTS.GetTimeStamp_High();
    uint64_t low  = RDTS.GetTimeStamp_Low();

    high = high << 32; 
    uint64_t joined = (high | low);

    LOG_INFO("BeamEvent") << "Raw Decoder Timestamp: " << joined << "\n";

    RDTSTime = joined;
    RDTSTrigger = RDTS.GetFlags();
    LOG_INFO("BeamEvent") << "Trigger: " << RDTSTrigger << "\n"; 

    //Separates seconds portion of the ticks 
    //From the nanoseconds
    long long RDTSTickSec = (RDTSTime * 2) / (int)(TMath::Power(10,8));
    RDTSTickSec = RDTSTickSec * (int)(TMath::Power(10,8)) / 2;
    long long RDTSTickNano = RDTSTime - RDTSTickSec;
  
    //Units are 20 nanoseconds ticks
    RDTSTimeSec  = 20.e-9 * RDTSTickSec;
    RDTSTimeNano = 20.    * RDTSTickNano;

  }

  auto const CTBHandle = e.getValidHandle< std::vector< raw::ctb::pdspctb > >("ctbrawdecoder:daq");
  LOG_INFO("BeamEvent") << "CTB valid? " << CTBHandle.isValid() << "\n";
  if(CTBHandle.isValid()){
    auto const & CTB = (*CTBHandle)[0];

    bool noHLT = true;

    LOG_INFO("BeamEvent") << "NTriggers: " << CTB.GetNTriggers() << "\n";
    for (size_t nTrig = 0; nTrig < CTB.GetNTriggers(); ++nTrig){

      raw::ctb::Trigger ctbTrig = CTB.GetTrigger(nTrig);      
      uint32_t  theType  = ctbTrig.word_type;
      ULong64_t theWord  = ctbTrig.trigger_word;
      ULong64_t theTS    = ctbTrig.timestamp;
      //LOG_INFO("BeamEvent").precision(21);
        
     
      if (theType == 2 ){

        long long deltaCTBtoRDTS = RDTSTime - (long long)theTS;        

        LOG_INFO("BeamEvent") << "Type 2. deltaT: " << deltaCTBtoRDTS << "\n";

        if( deltaCTBtoRDTS <= (fOffsetCTBtoRDTS + fToleranceCTBtoRDTS) 
        &&  deltaCTBtoRDTS >= (fOffsetCTBtoRDTS - fToleranceCTBtoRDTS) ){
        
          LOG_INFO("BeamEvent") << "Found the High Level Trigger" << "\n";
       
          HLTWord = theWord;
          HLTTS = theTS;
          LOG_INFO("BeamEvent") << HLTTS << "\n";

          //The High Level Trigger consists of 8 bits
          //HLT7 -> HLT0
          std::bitset<8> theHLT(theWord);
          LOG_INFO("BeamEvent") << "High Level Trigger: " << theHLT << "\n";

          
          //HLT5 corresponds to excluding Low Level Triggers
          //from the Beamline
          //So return 0, we'll skip this event
          if (theHLT[5]) {         
            noHLT = false;
            LOG_INFO("BeamEvent") << "HLT 5 activated. This is a Beam-Excluded event." << "\n";
            break;
          }
          else if (theHLT[0] && (theHLT.count() == 1)) {
            noHLT = false;
            LOG_INFO("BeamEvent") << "Only HLT 0 activated. This is just a random trigger. No beamline info was activated." << "\n";
            break;
          }
          else{
            noHLT = false;
            LOG_INFO("BeamEvent") << "Found valid beam event." << "\n";
            break;
          }
        }
      }       
    }

    if(noHLT){
      //This happens sometimes
      //Just skip the event
      LOG_INFO("BeamEvent") << "No High Level Trigger Found!" << "\n";
      return 0;
    }
    else{

      //Now check the channel statuses        
      LOG_INFO("BeamEvent") << "ChStatuses: " << CTB.GetNChStatuses() << "\n";
      for(size_t iStat = 0; iStat < CTB.GetNChStatuses(); ++ iStat){
        raw::ctb::ChStatus theStatus = CTB.GetChStatuse(iStat);

        uint32_t the_beam_hi    = theStatus.beam_hi; 
        uint32_t the_beam_lo    = theStatus.beam_lo; 
        long long the_timestamp = theStatus.timestamp; 
        
        LOG_INFO("BeamEvent") << "Timestamp : " << the_timestamp << "\n";
        int delta = HLTTS - the_timestamp;
        if( delta < 2 && delta >= 0 ){
          LOG_INFO("BeamEvent") << "Found Channel status matching  HLT timestamp" << "\n";
          LOG_INFO("BeamEvent") << "beam_hi   : " << std::bitset<5>(the_beam_hi) << "\n"; 
          LOG_INFO("BeamEvent") << "beam_lo   : " << std::bitset<4>(the_beam_lo) << "\n"; 
    
          std::bitset<4> beam_lo(the_beam_lo);
          std::bitset<5> beam_hi(the_beam_hi);
    
          BeamOn     =  beam_lo[0];
          BITrigger  =  beam_lo[1];
          Upstream   =  beam_lo[2];
          C1         =  beam_lo[3];
          C2         =  beam_hi[0];
          BP1        =  beam_hi[1];
          BP2        =  beam_hi[2];
          BP3        =  beam_hi[3];
          BP4        =  beam_hi[4];
    
    
          LOG_INFO("BeamEvent") << "%%%%Decoding the beam channels%%%" << "\n";
          LOG_INFO("BeamEvent") << "Beam On:    " << BeamOn    << "\n"
                    << "BI Trigger: " << BITrigger << "\n"
                    << "Upstream:   " << Upstream  << "\n"
                    << "C1:         " << C1        << "\n"
                    << "C2:         " << C2        << "\n"
                    << "BP1:        " << BP1       << "\n"
                    << "BP2:        " << BP2       << "\n"
                    << "BP3:        " << BP3       << "\n"
                    << "BP4:        " << BP4       << "\n";
          LOG_INFO("BeamEvent") << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << "\n" << std::endl;        
    

          //This means the beamline wasn't triggered
          if(!BITrigger) return 0;
    
          return HLTTS;

        }
      }
      return 0;
    }
  }
  LOG_INFO("BeamEvent") << "Error! Invalid CTB Handle!" << "\n";
  return 0;

}

void proto::BeamEvent::GetSpillInfo(art::Event & e){

    auto const PDTStampHandle = e.getValidHandle< std::vector< dune::ProtoDUNETimeStamp > > ("timingrawdecoder:daq");  
    std::vector< dune::ProtoDUNETimeStamp > PDTStampVec(*PDTStampHandle); 
    auto PDTStamp = PDTStampVec[0];            
          
    UInt_t    ver        = PDTStamp.getVersion();   
    double RunStart   = 2.e-8*PDTStamp.getLastRunStart();            
    SpillStart = 2.e-8*PDTStamp.getLastSpillStart();          
    SpillEnd   = 2.e-8*PDTStamp.getLastSpillEnd();            
    SpillStart_alt = PDTStamp.getLastSpillStart();


    if( 0ul == ~(SpillStart_alt) ) SpillStartValid = false;
    else SpillStartValid = true;


    LOG_INFO("BeamEvent") << PDTStamp.getLastSpillStart() << " " << SpillStart_alt << " " << SpillStartValid << "\n";
          
    std::cout.precision(dbl::max_digits10);
    LOG_INFO("BeamEvent") << "Version:     " << ver        << "\n";      
    LOG_INFO("BeamEvent") << "Run:         " << RunStart   << "\n";      
          
    LOG_INFO("BeamEvent") << "Spill Start: " << std::setw(20) << SpillStart <<  "\n";      
    std::cout << "Spill Start: " << SpillStart <<  std::endl;      
    LOG_INFO("BeamEvent") << "Spill End:   " << SpillEnd   <<  "\n";      
          
          
          
}

void proto::BeamEvent::TimeIn(art::Event & e, uint64_t time){

    /////Now look at the acqStamp coming out of IFBeam
    try{
      std::vector<double> acqStamp = FetchAndReport(time, "dip/acc/NORTH/NP04/POW/MBPL022699:acqStamp[]"); 

      acqStampMBPL   = 1.e-9 * joinHighLow(acqStamp[0],   acqStamp[1]); 
      std::cout.precision(dbl::max_digits10);
      std::cout << "MBPL: " << acqStampMBPL << std::endl;

      //Assign the calibration offset
      SpillOffset = SpillStart - acqStampMBPL;

    }
    catch( std::exception e){
      LOG_WARNING("BeamEvent") << "Could not get Spill time to time in\n";
    }
    
}

void proto::BeamEvent::MatchS11ToGen(){
  
  LOG_INFO("BeamEvent") << "Matching S11 To Gen" << "\n";
  for(size_t iT = 0; iT < beamspill->GetNT0(); ++iT){
    double GenTrigSec  = beamspill->GetT0(iT).first;
    double GenTrigNano = beamspill->GetT0(iT).second;

    double diff = GenTrigSec - s11Sec;
    diff += 1.e-9*(GenTrigNano - s11Nano);
    //LOG_INFO("BeamEvent") << iT << " " << diff  << "\n";

    if( fS11DiffLower < diff && diff < fS11DiffUpper ){
      LOG_INFO("BeamEvent") << "Found matching S11 and GenTrig!" << "\n";
      LOG_INFO("BeamEvent") << "diff: " << diff << "\n";

      beamspill->SetActiveTrigger( iT ); 
      beamevt->SetActiveTrigger( iT );
      beamevt->SetT0( beamspill->GetT0( iT ) );
//      LOG_INFO("BeamEvent") << "Found at: " << iT << "\n";
//      LOG_INFO("BeamEvent") << "Previous Active Trigger: " << beamspill->GetActiveTrigger() << "\n";
//      LOG_INFO("BeamEvent") << "Set event T0: " << beamevt->GetT0Sec() << " " << beamevt->GetT0Nano() << "\n";
      return;
    }
  }
  LOG_INFO("BeamEvent") << "Could not find matching time " << "\n";
  beamspill->SetUnmatched();
}


void proto::BeamEvent::MatchBeamToTPC(){

  LOG_INFO("BeamEvent") << "Matching in time between Beamline and TPC!!!" << "\n"; 
  for(size_t iT = 0; iT < beamspill->GetNT0(); ++iT){
    
    double GenTrigSec  = beamspill->GetT0(iT).first;
    double GenTrigNano = beamspill->GetT0(iT).second;

    //Separates seconds portion of the ticks 
    //From the nanoseconds
    long long RDTSTickSec = (RDTSTime * 2) / (int)(TMath::Power(10,8));
    RDTSTickSec = RDTSTickSec * (int)(TMath::Power(10,8)) / 2;
    long long RDTSTickNano = RDTSTime - RDTSTickSec;

    //Units are 20 nanoseconds ticks
    double RDTSTimeSec  = 20.e-9 * RDTSTickSec;
    double RDTSTimeNano = 20.    * RDTSTickNano;

 //   LOG_INFO("BeamEvent") << "RDTS: " << RDTSTimeSec << " " << RDTSTimeNano << "\n";

    double diffSec = RDTSTimeSec - GenTrigSec - SpillOffset;
    LOG_INFO("BeamEvent") << "RDTSTimeSec - GenTrigSec " << RDTSTimeSec - GenTrigSec << "\n";
    double diffNano = 1.e-09*(RDTSTimeNano - GenTrigNano);
    LOG_INFO("BeamEvent") << "diff: " << diffSec << " " << diffNano << "\n";

    double diff = diffSec + diffNano; 
//    LOG_INFO("BeamEvent") << diff << "\n";
  
    if( ( fTimingCalibration - fCalibrationTolerance < diff ) && (fTimingCalibration + fCalibrationTolerance > diff) ){
      LOG_INFO("BeamEvent") << "FOUND MATCHING TIME!!!" << "\n";
      LOG_INFO("BeamEvent") << "diff: " << diff << "\n";


      beamspill->SetActiveTrigger( iT ); 
      beamevt->SetActiveTrigger( iT );
      beamevt->SetT0( beamspill->GetT0( iT ) );
      LOG_INFO("BeamEvent") << "Set event T0: " << beamevt->GetT0Sec() << " " << beamevt->GetT0Nano() << "\n";
      return;
    }
  }

  LOG_INFO("BeamEvent") << "Could not find matching time " << "\n";
  beamspill->SetUnmatched();
}

void proto::BeamEvent::SetCKovInfo(){

  beam::CKov theCKov = beamspill->GetCKov0();
  theCKov.trigger = C1;
  beamevt->SetCKov0( theCKov );

  theCKov = beamspill->GetCKov1();
  theCKov.trigger = C2;
  beamevt->SetCKov1( theCKov );

}

void proto::BeamEvent::SetBeamEvent(){

  if( !beamspill->CheckIsMatched() ){
    LOG_INFO("BeamEvent") << "art Event is unmatched to Beam Spill " << "\n";
    return;
  }

  LOG_INFO("BeamEvent") << "Setting beam event for matched event" << "\n";

  size_t activeTrigger = beamspill->GetActiveTrigger();
  beamevt->SetActiveTrigger( activeTrigger );
  
  LOG_INFO("BeamEvent") << "SetActiveTrigger " << beamevt->GetActiveTrigger() << "\n";

  beamevt->SetT0( beamspill->GetT0( activeTrigger ) );

  LOG_INFO("BeamEvent") << "Set T0  " << beamevt->GetT0Sec() << " " << beamevt->GetT0Nano() << "\n" << "\n"; 


  LOG_INFO("BeamEvent") << "Setting FBM statuses" << "\n";

  for( size_t i = 0; i < fDevices.size(); ++i){
    std::string theName = fDevices[i];
    beamevt->SetFBMTrigger( theName, beamspill->GetFBM(theName, activeTrigger) );    
    LOG_INFO("BeamEvent") << "beamevt monitor " << theName << " has "
              << beamevt->GetActiveFibers( theName ).size() << " active Fibers" << "\n";
  }

  LOG_INFO("BeamEvent") << "Setting TOF info for beamevt " << "\n";
  beamevt->SetTOF0Trigger( beamspill->GetTOF0( activeTrigger ) );
  beamevt->SetTOF1Trigger( beamspill->GetTOF1( activeTrigger ) );
  beamevt->DecodeTOF();
  beamevt->SetTOFChan(  beamspill->GetTOFChan( activeTrigger ) );
  LOG_INFO("BeamEvent") << "beamevt has TOF " << beamevt->GetTOF() 
            << " and TOFChan "    << beamevt->GetTOFChan() << "\n";


  LOG_INFO("BeamEvent") << "Finished adding info to beamevt " << "\n";


}

////////////////////////
// Producer Method (reads in the event and derives values)
void proto::BeamEvent::produce(art::Event & e){
  //Reset 
  acqTime = 0;
  acqStampMBPL = 0;
  HLTWord = 0;
  HLTTS = 0;
  BeamOn      = -1;
  BITrigger   = -1;
  RDTSTrigger = -1;
  Upstream    = -1;
  C1          = -1;
  C2          = -1;
  BP1         = -1;
  BP2         = -1;
  BP3         = -1;
  BP4         = -1; 
  SpillStart  = -1;
  SpillEnd    = -1;
  SpillOffset = -1;

  s11Nano     = -1.;
  s11Sec      = -1.;

  ActiveTriggerTime = -1;
  RDTSTime   = 0;

//  bool usedEventTime = false;


  eventNum = e.event();
  runNum = e.run();
  subRunNum = e.subRun();

  if( e.time().timeHigh() == 0 ) eventTime = e.time().timeLow();
  else eventTime = e.time().timeHigh();


   
  // Create a new beam event (note the "new" here)  
  beamevt = new beam::ProtoDUNEBeamEvent();
  beamspill = new beam::ProtoDUNEBeamSpill();

  // Get the coordinate system conversions
  if(!rotated) BeamMonitorBasisVectors();

  validTimeStamp = GetRawDecoderInfo(e);
  



  //Get Spill Information
  //This stores Spill Start, Spill End, 
  //And Prev Spill Start
  GetSpillInfo(e);

  
  //Check if we have a valid beam trigger
  //If not, just place an empty beamevt
  //and move on
  //
  //Also check if we've gotten good spill info
  //
  //Or if we're forcing to read out the Beamline Info
  if( ( (RDTSTrigger == 12) ) || fForceRead ){

    //Start getting beam event info
    uint64_t fetch_time = uint64_t( RDTSTime * 2e-8 ) + fFetchOffset;
    uint64_t fetch_time_down = uint64_t( RDTSTime * 2e-8 );
    LOG_INFO("BeamEvent") << "RDTSTime: " <<  uint64_t( RDTSTime * 2e-8 ) << "\n";

    //Check if we are still using the same spill information.
    //
    //If it's a new spill: Get new info from the database 
    //Each 'parse' command fetches new data from the database
    //
    //If it's the same spill then just pass the old BeamEvent
    //Object. Its 'active trigger' info will be updated below
    //
    //Also: Check if the SpillStart is invalid. If the RDTSTime is  
    //Greater than 5, it's in a new spill, so fetch again
    //
    //Can be overridden with a flag from the fcl
    if( ( ( PrevStart != SpillStart) || ( !SpillStartValid && ( abs(RDTSTimeSec - PrevRDTSTimeSec) > 5 ) ) )
    || fForceNewFetch){
      LOG_INFO("BeamEvent") << "New spill or forced new fetch. Getting new beamspill info" << "\n";

      // Parse the Time of Flight Counter data for the list
      // of times that we are using
      parseXTOF(fetch_time);
      

      // Parse the Beam postion counter information for the list
      // of time that we are using
      InitXBPFInfo(beamspill);

      if( gotGeneralTrigger ){ 
        parseXBPF(fetch_time);
      }

      parseXCET(fetch_time);

      try{
        current = FetchAndReport(fetch_time_down, "dip/acc/NORTH/NP04/POW/MBPL022699:current");
        gotCurrent = true;
        LOG_INFO("BeamEvent") << "Current: " << current[0] << "\n";
      }
      catch( std::exception e){
        LOG_WARNING("BeamEvent") << "Could not get magnet current\n";
        gotCurrent = false;
      }
   
      //Set PrevStart to SpillStart here
      //so that we don't skip in the case 
      //the first event in the spill did not
      //a good beamline trigger
      PrevStart = SpillStart;
      PrevRDTSTimeSec = RDTSTimeSec;

    }
    else{
      LOG_INFO("BeamEvent") << "Same spill. Reusing beamspill info" << "\n";
      LOG_INFO("BeamEvent") << prev_beamspill.GetNT0() << "\n";
      *beamspill = prev_beamspill;
      LOG_INFO("BeamEvent") << beamspill->GetNT0() << "\n";
    }

    LOG_INFO("BeamEvent") << "NGoodParticles: " << beamspill->GetNT0()            << "\n";
    LOG_INFO("BeamEvent") << "NTOF0: "          << beamspill->GetNTOF0Triggers()  << "\n";
    LOG_INFO("BeamEvent") << "NTOF1: "          << beamspill->GetNTOF1Triggers()  << "\n";
    LOG_INFO("BeamEvent") << "acqTime: "        << acqTime                      << "\n";
    LOG_INFO("BeamEvent") << "NXBPF: "          << beamspill->GetNFBMTriggers(fDevices[0]) << "\n";

    if( fMatchTime ){
   
      //Now do the matching in Time:
      //Get the conversion from the ProtoDUNE Timing system
      //To the one in the SPS.
      //

      if(SpillStartValid){
        TimeIn(e, fetch_time_down);
        LOG_INFO("BeamEvent") << "SpillOffset " << SpillOffset << "\n";
        
        //If not successfully timed in, 
        //Oh well. It won't cause a crash
        //But it won't be matched. 
        MatchBeamToTPC();
      }
      else{
        try{
          getS11Info(fetch_time); 
        }
        catch( std::exception e ){
          LOG_WARNING("BeamEvent") << "Could not get S11 Info\n";
        }

        //Again, it won't crash, but it won't match          
        MatchS11ToGen();
      }


      if( beamspill->CheckIsMatched() ){
        std::pair<double,double> theTime = beamspill->GetT0(beamspill->GetActiveTrigger());
        ActiveTriggerTime = theTime.first + theTime.second*1.e-9;
        LOG_INFO("BeamEvent") << "Trigger: " << beamspill->GetActiveTrigger() << " " << ActiveTriggerTime << "\n";       
      
        //Pass the information to the beamevent
        SetBeamEvent();

        MakeTrack( beamspill->GetActiveTrigger() );
        LOG_INFO("BeamEvent") << "Added " << beamevt->GetNBeamTracks() << " tracks to the beam spill" << "\n";

        //Momentum
        if( gotCurrent ){
          MomentumSpec( beamspill->GetActiveTrigger() ); 
        }
        LOG_INFO("BeamEvent") << "Got NRecoBeamMomenta: " << beamevt->GetNRecoBeamMomenta() << "\n";
      }
    }


    //Pass beamspill to the next event;
    //Erase the Track and Reco Momentum info
    prev_beamspill = *beamspill;
    prev_beamspill.ClearBeamTracks();
    prev_beamspill.ClearRecoBeamMomenta();
    prev_beamspill.SetUnmatched();

  }
  //Start of a new spill, but the first event was 
  //Not a beam trigger. In this case, it would not
  //have been filled with info in the block above
  //So let's make it empty so we aren't putting 
  //old spill info in the new event

  //Or. If this was not a beam trigger, and the RDTSTime
  //Comes from a new spill while the SpillStart was invalid, pass it. 
  else if( ( ( PrevStart != SpillStart ) || ( !SpillStartValid && ( abs(RDTSTimeSec - PrevRDTSTimeSec) > 5 ) ) ) 
    && RDTSTrigger != 12 ){
    prev_beamspill = *beamspill;   
  }
  
  beamevt->SetBITrigger(BITrigger);
  beamevt->SetTimingTrigger(RDTSTrigger);
  SetCKovInfo();
  beamevt->SetSpillStart(SpillStart);
  beamevt->SetSpillOffset(SpillOffset);
  beamevt->SetCTBTimestamp( HLTTS );
  beamevt->SetRDTimestamp( RDTSTime );

  std::unique_ptr<std::vector<beam::ProtoDUNEBeamEvent> > beamData(new std::vector<beam::ProtoDUNEBeamEvent>);
  beamData->push_back(beam::ProtoDUNEBeamEvent(*beamevt));
  e.put(std::move(beamData));
  delete beamevt;
  delete beamspill;
 
  // Write out the to tree
  if( fSaveOutTree )fOutTree->Fill();
 
}
// END BeamEvent::produce
////////////////////////

void proto::BeamEvent::InitXBPFInfo(beam::ProtoDUNEBeamSpill * beamspill){
  // Places a dummy trigger vector for each device

  // Make a vector the names of each of the devices being readout
  std::vector<std::string> monitors;
  size_t nDev = 0; 

  for(size_t id = 0; id < fDevices.size(); ++id){
    std::string name = fDevices[id];
    // Put the current device name on the list
    monitors.push_back(name);
    nDev++;
  }

  beamspill->InitFBMs(monitors);
}
// END BeamEvent::InitXBPFInfo
////////////////////////

void proto::BeamEvent::getS11Info(uint64_t time){
  LOG_INFO("BeamEvent") << "Getting S11 Info " << "\n";
//  uint64_t time = uint64_t( RDTSTime * 2e-8 + 5 );


  std::vector<double> coarseS11         = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/S11:coarse[]");
  std::vector<double> fracS11           = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/S11:frac[]"); 
  std::vector<double> secondsS11        = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/S11:seconds[]"); 
  std::vector<double> timestampCountS11 = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/S11:timestampCount"); 
  
  int s11Count = (int)timestampCountS11[0];

  LOG_INFO("BeamEvent") << "Found " << s11Count << " returned signals" << "\n";

  //Separates seconds portion of the ticks 
  //From the nanoseconds
  long long RDTSTickSec = (RDTSTime * 2) / (int)(TMath::Power(10,8));
  RDTSTickSec = RDTSTickSec * (int)(TMath::Power(10,8)) / 2;
  long long RDTSTickNano = RDTSTime - RDTSTickSec;

  //Units are 20 nanoseconds ticks
  double RDTSTimeSec  = 20.e-9 * RDTSTickSec;
  double RDTSTimeNano = 20.    * RDTSTickNano;


  for(int i = 0; i < s11Count; ++i){
    double nano = 8.*coarseS11[i] + fracS11[i]/512.;

    double diffSec = secondsS11[2*i + 1] - RDTSTimeSec - fOffsetTAI;
    double diffNano = nano - RDTSTimeNano;

    //LOG_INFO("BeamEvent") << i << " " << secondsS11[2*i + 1] << " " << RDTSTimeSec << "\n";
    //LOG_INFO("BeamEvent") << "\t" << nano << " " << RDTSTimeNano << "\n";
    LOG_INFO("BeamEvent") << i << " diffSec " << diffSec << "\n"; 
    LOG_INFO("BeamEvent") << i << " diffNano " << diffNano << "\n"; 

    double diff = diffSec + 1.e-9*diffNano;
    if( fRDTSToS11Lower < diff && diff < fRDTSToS11Upper ){
      LOG_INFO("BeamEvent") << "FOUND Match between S11 and RDTS" << "\n";  
      s11Nano = nano;
      s11Sec  = secondsS11[2*i + 1] - fOffsetTAI;
      return;
    }
  }

  LOG_WARNING("BeamEvent") << "Could not match RDTS to S11\n";

}


void proto::BeamEvent::parseXTOF(uint64_t time){

  std::vector<double> coarseGeneralTrigger;         
  std::vector<double> fracGeneralTrigger;           
  std::vector<double> acqStampGeneralTrigger;       
  std::vector<double> secondsGeneralTrigger;        
  std::vector<double> timestampCountGeneralTrigger; 
  

  try{
    LOG_INFO("BeamEvent") << "Getting General trigger info " << "\n";

    coarseGeneralTrigger         = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:coarse[]");
    fracGeneralTrigger           = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:frac[]"); 
    acqStampGeneralTrigger       = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:acqStamp[]"); 
    secondsGeneralTrigger        = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:seconds[]"); 
    timestampCountGeneralTrigger = FetchAndReport(time, "dip/acc/NORTH/NP04/BI/XBTF/GeneralTrigger:timestampCount"); 

    LOG_INFO("BeamEvent") << "timestampCounts: " << timestampCountGeneralTrigger[0] << "\n";
    gotGeneralTrigger = true;

    uint64_t low = (uint64_t)acqStampGeneralTrigger[1];
    uint64_t high = (uint64_t)acqStampGeneralTrigger[0];
    high = high << 32;
    acqTime = ( high | low ) / 1000000000.; 

  }
  catch(std::exception e){
    LOG_WARNING("BeamEvent") << "Could get GeneralTrigger information!!" << "\n";
    gotGeneralTrigger = false;
    return;
  }
  
  std::vector<double> coarseTOF1A;  
  std::vector<double> fracTOF1A;    
  std::vector<double> secondsTOF1A; 

  std::vector<double> coarseTOF1B;  
  std::vector<double> fracTOF1B;    
  std::vector<double> secondsTOF1B; 

  std::vector<double> coarseTOF2A;  
  std::vector<double> fracTOF2A;    
  std::vector<double> secondsTOF2A; 

  std::vector<double> coarseTOF2B;  
  std::vector<double> fracTOF2B;    
  std::vector<double> secondsTOF2B; 

  try{
    LOG_INFO("BeamEvent") << "Getting TOF1A info: " << fTOF1 << "\n";
    coarseTOF1A  = FetchAndReport(time, fXTOFPrefix + fTOF1A + ":coarse[]");
    fracTOF1A    = FetchAndReport(time, fXTOFPrefix + fTOF1A + ":frac[]");
    secondsTOF1A = FetchAndReport(time, fXTOFPrefix + fTOF1A + ":seconds[]");
    LOG_INFO("BeamEvent") << "Size of coarse,frac: " << coarseTOF1A.size() << " " << fracTOF1A.size() << "\n"; 

    LOG_INFO("BeamEvent") << "Getting TOF1B info: " << fTOF1 << "\n";
    coarseTOF1B  = FetchAndReport(time, fXTOFPrefix + fTOF1B + ":coarse[]");
    fracTOF1B    = FetchAndReport(time, fXTOFPrefix + fTOF1B + ":frac[]");
    secondsTOF1B = FetchAndReport(time, fXTOFPrefix + fTOF1B + ":seconds[]");
    LOG_INFO("BeamEvent") << "Size of coarse,frac: " << coarseTOF1B.size() << " " << fracTOF1B.size() << "\n"; 

    LOG_INFO("BeamEvent") << "Getting TOF2A info: " << fTOF2 << "\n";
    coarseTOF2A  = FetchAndReport(time, fXTOFPrefix + fTOF2A + ":coarse[]");
    fracTOF2A    = FetchAndReport(time, fXTOFPrefix + fTOF2A + ":frac[]");
    secondsTOF2A = FetchAndReport(time, fXTOFPrefix + fTOF2A + ":seconds[]");
    LOG_INFO("BeamEvent") << "Size of coarse,frac: " << coarseTOF2A.size() << " " << fracTOF2A.size() << "\n"; 

    LOG_INFO("BeamEvent") << "Getting TOF2B info: " << fTOF2 << "\n";
    coarseTOF2B  = FetchAndReport(time, fXTOFPrefix + fTOF2B + ":coarse[]");
    fracTOF2B    = FetchAndReport(time, fXTOFPrefix + fTOF2B + ":frac[]");
    secondsTOF2B = FetchAndReport(time, fXTOFPrefix + fTOF2B + ":seconds[]");
    LOG_INFO("BeamEvent") << "Size of coarse,frac: " << coarseTOF2B.size() << " " << fracTOF2B.size() << "\n"; 

    gotTOFs = true;
  }
  catch(std::exception e){
    LOG_WARNING("BeamEvent") << "Could get TOF information!!" << "\n";
    gotTOFs = false;
  }

  std::vector<std::pair<double,double>> unorderedGenTrigTime;
  std::vector<std::pair<double,double>> unorderedTOF1ATime;
  std::vector<std::pair<double,double>> unorderedTOF1BTime;
  std::vector<std::pair<double,double>> unorderedTOF2ATime;
  std::vector<std::pair<double,double>> unorderedTOF2BTime;


  for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
    //LOG_INFO("BeamEvent") << i << " " << secondsGeneralTrigger[2*i + 1] << " "  << 8.*coarseGeneralTrigger[i] + fracGeneralTrigger[i]/512. << "\n";
//    LOG_INFO("BeamEvent") << "\t" << std::setw(15) << secondsGeneralTrigger[2*i + 1] + 1.e-9*(8.*coarseGeneralTrigger[i] + fracGeneralTrigger[i]/512.) << "\n";

    //2*i + 1 because the format is weird
    fGenTrigSec    = secondsGeneralTrigger[2*i + 1];
    fGenTrigCoarse = coarseGeneralTrigger[i];
    fGenTrigFrac   = fracGeneralTrigger[i];

    unorderedGenTrigTime.push_back( std::make_pair(fGenTrigSec, (fGenTrigCoarse*8. + fGenTrigFrac/512.)) );

    if (fGenTrigCoarse == 0.0 && fGenTrigFrac == 0.0 && fGenTrigSec == 0.0) break;
    if(fDebugTOFs)fGenTrigTree->Fill();
  }

  if( gotTOFs ){
    for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
//      LOG_INFO("BeamEvent") << "TOF1A " << i << " " << secondsTOF1A[2*i+1] << " "  << 8.*coarseTOF1A[i] +  fracTOF1A[i]/512. << "\n";
      fXTOF1ACoarse = coarseTOF1A[i];
      fXTOF1AFrac   = fracTOF1A[i];
      fXTOF1ASec    = secondsTOF1A[2*i + 1];

      if(fXTOF1ASec < secondsTOF1A[1]) break; 

      if (fXTOF1ACoarse == 0.0 && fXTOF1AFrac == 0.0 && fXTOF1ASec == 0.0) break;
      unorderedTOF1ATime.push_back(std::make_pair(fXTOF1ASec, (fXTOF1ACoarse*8. + fXTOF1AFrac/512.)) );

      if(fDebugTOFs)fXTOF1ATree->Fill();
    }
      
    for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
//      LOG_INFO("BeamEvent") << "TOF1B " << i << " " << secondsTOF1B[2*i+1] << " "  << 8.*coarseTOF1B[i] + fracTOF1B[i]/512. << "\n";
      fXTOF1BCoarse = coarseTOF1B[i];
      fXTOF1BFrac   = fracTOF1B[i];
      fXTOF1BSec    = secondsTOF1B[2*i + 1];

      if(fXTOF1BSec < secondsTOF1B[1]) break; 

      if (fXTOF1BCoarse == 0.0 && fXTOF1BFrac == 0.0 && fXTOF1BSec == 0.0) break;
      unorderedTOF1BTime.push_back(std::make_pair(fXTOF1BSec, (fXTOF1BCoarse*8. + fXTOF1BFrac/512.)) );
      if(fDebugTOFs)fXTOF1BTree->Fill();
    }

    for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
//      LOG_INFO("BeamEvent") << "TOF2A " << i << " " << secondsTOF2A[2*i+1] << " "  << 8.*coarseTOF2A[i] +  fracTOF2A[i]/512. << "\n";
      fXTOF2ACoarse = coarseTOF2A[i];
      fXTOF2AFrac   = fracTOF2A[i];
      fXTOF2ASec    = secondsTOF2A[2*i + 1];

      if(fXTOF2ASec < secondsTOF2A[1]) break; 

      if (fXTOF2ACoarse == 0.0 && fXTOF2AFrac == 0.0 && fXTOF2ASec == 0.0) break;
      unorderedTOF2ATime.push_back(std::make_pair(fXTOF2ASec, (fXTOF2ACoarse*8. + fXTOF2AFrac/512.)) );

      if(diff2A->size()) diff2A->clear();
      for(size_t j = 0; j < timestampCountGeneralTrigger[0]; ++j){
        diff2A->push_back( 1.e9*(unorderedTOF2ATime.back().first - unorderedGenTrigTime[j].first) + (unorderedTOF2ATime.back().second - unorderedGenTrigTime[j].second) );
      }
      
      if(fDebugTOFs)fXTOF2ATree->Fill();

    }  

    for(size_t i = 0; i < timestampCountGeneralTrigger[0]; ++i){
//      LOG_INFO("BeamEvent") << "TOF2B " << i << " " << secondsTOF2B[2*i+1] << " "  << 8.*coarseTOF2B[i] +  fracTOF2B[i]/512. << "\n";
      fXTOF2BCoarse = coarseTOF2B[i];
      fXTOF2BFrac   = fracTOF2B[i];
      fXTOF2BSec    = secondsTOF2B[2*i + 1];

      if(fXTOF2BSec < secondsTOF2B[1]) break; 

      if (fXTOF2BCoarse == 0.0 && fXTOF2BFrac == 0.0 && fXTOF2BSec == 0.0) break;
      unorderedTOF2BTime.push_back(std::make_pair(fXTOF2BSec, (fXTOF2BCoarse*8. + fXTOF2BFrac/512.) ));
      if(diff2B->size()) diff2B->clear();
      for(size_t j = 0; j < timestampCountGeneralTrigger[0]; ++j){
        diff2B->push_back( 1.e9*(unorderedTOF2BTime.back().first - unorderedGenTrigTime[j].first) + (unorderedTOF2BTime.back().second - unorderedGenTrigTime[j].second) );
      }
      if(fDebugTOFs)fXTOF2BTree->Fill();
    }
  }


  for(size_t iT = 0; iT < unorderedGenTrigTime.size(); ++iT){
    
    bool found_TOF = false;

    double the_gen_sec = unorderedGenTrigTime[iT].first;
    double the_gen_ns = unorderedGenTrigTime[iT].second;

    double the_TOF1_sec = -1.;
    double the_TOF2_sec = -1.;
    double the_TOF1_ns = -1.;
    double the_TOF2_ns = -1.;

    //1A2A = 0; 1B2A = 1, 1A2B = 2,  1B2B = 3
    //Add 1 for 1B, add 2 for 2B
    int channel = 0;

    if( gotTOFs ){
      for(size_t iT1A = 0; iT1A < unorderedTOF1ATime.size(); ++iT1A){
        double TOF1A_sec = unorderedTOF1ATime[iT1A].first;
        double TOF1A_ns = unorderedTOF1ATime[iT1A].second;

        double delta_1A = 1.e9*(the_gen_sec - TOF1A_sec) + the_gen_ns - TOF1A_ns;
        if(delta_1A < 0.){
         // LOG_INFO("BeamEvent") << "Passed TOF1A" << "\n";
          //TOF1A_passed = true;
          break;
        }

        if( delta_1A > 500. ) continue;     

        //If here, then 0 < delta_1A < 500ns
        //Check the TOF2 times. 
        
        for(size_t iT2A = 0; iT2A < unorderedTOF2ATime.size(); ++iT2A){
          double TOF2A_sec = unorderedTOF2ATime[iT2A].first;
          double TOF2A_ns = unorderedTOF2ATime[iT2A].second;

          double delta = 1.e9*(TOF2A_sec - TOF1A_sec) + TOF2A_ns - TOF1A_ns;
          
          //Overtaken 1A
          if(delta < 0.){
            continue; 
          }

          if( delta > 500. ){
            break;
          }
          else{

            //If here, then TOF1 is within 500 ns below TOF2

//            LOG_INFO("BeamEvent") << "Found matching TOF2A and TOF1A" << "\n";
            found_TOF = true;

            the_TOF2_sec = TOF2A_sec;
            the_TOF2_ns  = TOF2A_ns;

            the_TOF1_sec = TOF1A_sec;
            the_TOF1_ns  = TOF1A_ns;

            channel = k1A2A;

            break;
          }       
        }

        if(!found_TOF){
          for(size_t iT2B = 0; iT2B < unorderedTOF2BTime.size(); ++iT2B){
            double TOF2B_sec = unorderedTOF2BTime[iT2B].first;
            double TOF2B_ns = unorderedTOF2BTime[iT2B].second;

            double delta = 1.e9*(TOF2B_sec - TOF1A_sec) + TOF2B_ns - TOF1A_ns;
            
            //Overtaken 1A
            if(delta < 0.){
              continue;
            }

            if( delta > 500. ){
              break;
            }
            else{
              //If here, then TOF1 is within 500 ns below TOF2

//              LOG_INFO("BeamEvent") << "Found matching TOF2B and TOF1A" << "\n";
              found_TOF = true;

              the_TOF2_sec = TOF2B_sec;
              the_TOF2_ns  = TOF2B_ns;
  
              the_TOF1_sec = TOF1A_sec;
              the_TOF1_ns  = TOF1A_ns;

              channel = k1A2B;

              break;
            }       
          }
        }

        if(found_TOF) break;
      }

      //Now check 1B with 2A and 2B
      if(!found_TOF){
        
        for(size_t iT1B = 0; iT1B < unorderedTOF1BTime.size(); ++iT1B){
          double TOF1B_sec = unorderedTOF1BTime[iT1B].first;
          double TOF1B_ns = unorderedTOF1BTime[iT1B].second;

          double delta_1B = 1.e9*(the_gen_sec - TOF1B_sec) + the_gen_ns - TOF1B_ns;
          if(delta_1B < 0.){
           // LOG_INFO("BeamEvent") << "Passed TOF1B" << "\n";
            //TOF1B_passed = true;
            break;
          }

          if( delta_1B > 500. ) continue;     

          //If here, then 0 < delta_1B < 500ns
          //Check the TOF2 times. 
          
          for(size_t iT2A = 0; iT2A < unorderedTOF2ATime.size(); ++iT2A){
            double TOF2A_sec = unorderedTOF2ATime[iT2A].first;
            double TOF2A_ns = unorderedTOF2ATime[iT2A].second;

            double delta = 1.e9*(TOF2A_sec - TOF1B_sec) + TOF2A_ns - TOF1B_ns;
            
            //Overtaken 1B
            if(delta < 0.){
              continue;
            }

            if( delta > 500. ){
              break;
            }
            else{

              //If here, then TOF1 is within 500 ns below TOF2

//              LOG_INFO("BeamEvent") << "Found matching TOF2A and TOF1B" << "\n";
              found_TOF = true;

              the_TOF2_sec = TOF2A_sec;
              the_TOF2_ns  = TOF2A_ns;

              the_TOF1_sec = TOF1B_sec;
              the_TOF1_ns  = TOF1B_ns;

              channel = k1B2A;

              break;
            }       
          }

          if(!found_TOF){
            for(size_t iT2B = 0; iT2B < unorderedTOF2BTime.size(); ++iT2B){
              double TOF2B_sec = unorderedTOF2BTime[iT2B].first;
              double TOF2B_ns = unorderedTOF2BTime[iT2B].second;

              double delta = 1.e9*(TOF2B_sec - TOF1B_sec) + TOF2B_ns - TOF1B_ns;
              
              //Overtaken 1B
              if(delta < 0.){
                continue;
              }

              if( delta > 500. ){
                break;
              }
              else{
                //If here, then TOF1 is within 500 ns below TOF2

//                LOG_INFO("BeamEvent") << "Found matching TOF2B and TOF1B" << "\n";
                found_TOF = true;

                the_TOF2_sec = TOF2B_sec;
                the_TOF2_ns  = TOF2B_ns;
  
                the_TOF1_sec = TOF1B_sec;
                the_TOF1_ns  = TOF1B_ns;

                channel = k1B2B;

                break;
              }       
            }
          }
          
          if(found_TOF) break;
        }    
      }
    }

    if(found_TOF){
      //Convert from TAI to UTC at this point

//      LOG_INFO("BeamEvent") << "Adding matched tof" << "\n";

      beamspill->AddT0(std::make_pair(the_gen_sec - fOffsetTAI, the_gen_ns));
      beamspill->AddTOF0Trigger(std::make_pair(the_TOF1_sec - fOffsetTAI, the_TOF1_ns));
      beamspill->AddTOF1Trigger(std::make_pair(the_TOF2_sec - fOffsetTAI, the_TOF2_ns));
      beamspill->AddTOFChan(channel);        
    }
    else{
      //Add dummy

//      LOG_INFO("BeamEvent") << "Adding unmatched tof" << "\n";

      beamspill->AddT0(std::make_pair(the_gen_sec - fOffsetTAI, the_gen_ns));
      beamspill->AddTOF0Trigger(std::make_pair(0., 0.));
      beamspill->AddTOF1Trigger(std::make_pair(0., 0.));
      beamspill->AddTOFChan(-1);        
    }

  }

}
// END BeamEvent::parseXTOF
////////////////////////

void proto::BeamEvent::parseXCET(uint64_t time){
  if(fCKov1 != ""){  

    try{
      std::vector< double > pressureCKov1   = FetchAndReport(time, fXCETPrefix + fCKov1 + ":pressure");
      CKov1Pressure = pressureCKov1[0];
    }
    catch( std::exception e){
      LOG_WARNING("BeamEvent") << "Could not get Cerenkov 1 Pressure\n";
      CKov1Pressure = 0.; 
    }

  }
  
  if(fCKov2 != ""){
    try{
      std::vector< double > pressureCKov2   = FetchAndReport(time, fXCETPrefix + fCKov2 + ":pressure");
      CKov2Pressure = pressureCKov2[0];
    }
    catch( std::exception e){
      LOG_WARNING("BeamEvent") << "Could not get Cerenkov 2 Pressure\n";
      CKov2Pressure = 0.; 
    }  
  }

  beam::CKov CKov1Status, CKov2Status;

  CKov1Status.pressure  = CKov1Pressure;
  CKov1Status.trigger   = C1;
  beamspill->SetCKov0( CKov1Status );


  CKov2Status.pressure  = CKov2Pressure;
  CKov2Status.trigger   = C2;
  beamspill->SetCKov1( CKov2Status );

}
// END BeamEvent::parseXCET
////////////////////////



////////////////////////
// 
void proto::BeamEvent::parseGeneralXBPF(std::string name, uint64_t time, size_t ID){

  // Retrieve the number of counts in the BPF
  std::vector<double> counts;

  try{
    counts = FetchAndReport(time, fXBPFPrefix + name + ":countsRecords[]");
  }
  catch( std::exception e){
    LOG_WARNING("BeamEvent") << "Could not fetch " << fXBPFPrefix + name + ":countsRecords[]\n";
    return;
  }
  
  LOG_INFO("BeamEvent") << "Counts: " << counts[1] << "\n";

  std::vector<double> data;
  try{
    data = FetchAndReport(time, fXBPFPrefix + name + ":eventsData[]");
  }
  catch( std::exception e){
    LOG_WARNING("BeamEvent") << "Could not fetch " << fXBPFPrefix + name + ":eventsData[]\n";
    return;
  }

  // If the number of counts is larger than the number of general triggers
  // make note
  if(counts[1] != beamspill->GetNT0()){
    LOG_WARNING("BeamEvent") << "WARNING MISMATCH " << counts[1] << " " << beamspill->GetNT0() << "\n";
  }
  
  
  beam::FBM fbm;
  fbm.ID = ID;
    
  //Use this just in case any are out of sync?
  //Shouldn't be, but just to be safe...
  //Helps cut down on time
  std::vector<size_t> leftOvers;      
  for(size_t lo = 0; lo < beamspill->GetNT0(); ++lo){
    leftOvers.push_back(lo);
  }
 
  for(size_t i = 0; i < counts[1]; ++i){      
      std::cout << "Count: " << i << std::endl;
    
    for(int j = 0; j < 10; ++j){
      double theData = data[20*i + (2*j + 1)];
      std::cout << std::setw(15) << theData ;
      if(j < 4){
	fbm.timeData[j] = theData;           
      }else{
	fbm.fiberData[j - 4] = theData;
      } 
    } 

    // Check the time data for corruption
    if(fbm.timeData[1] < .0000001){
//      std::cout << "Skipping bad time" << std::endl;
      continue;
    } 

//    std::cout << "Replacing  " << i << std::endl;
//    beamspill->ReplaceFBMTrigger(name, fbm, i);
//    leftOvers.erase( std::find(leftOvers.begin(), leftOvers.end(), i ) );

    // Timestamp is in units of sec + 8ns ticks
//    fbm.timeStamp = fbm.timeData[3] + fbm.timeData[2]*8.e-9;  
//    std::cout << fbm.timeData[3] << " " << fbm.timeData[2]*8.e-9 << " " << fbm.timeData[1] << " " << 8.e-9*fbm.timeData[0] << std::endl;
    
    //Go through the valid Good Particles, and emplace the FBM 
 //   std::cout << "Checking " << beamspill->GetNT0() << " triggers " << leftOvers.size() << std::endl;

    for(std::vector<size_t>::iterator ip = leftOvers.begin(); ip != leftOvers.end(); ++ip){

      // Compute the time delta between the timeStamp and the T0, see if it's less than 500ns away

      double delta = (beamspill->GetT0Sec(*ip)  - (fbm.timeData[3] - fOffsetTAI) );
      delta += 1.e-9*( beamspill->GetT0Nano(*ip) - 8.*fbm.timeData[2] ); 

 //     std::cout << beamspill->GetT0Sec(*ip) << " " << beamspill->GetT0Nano(*ip)*1.e-9 << " " << delta << std::endl;

      if( delta < -1.e-7 && delta > -10.e-7 ){

	 std::cout << "Replacing  " << *ip << std::endl;
	beamspill->ReplaceFBMTrigger(name, fbm, *ip);
	leftOvers.erase(ip);
	break;
      } 
    } 


 //   std::cout << std::endl;
  } 

  if( leftOvers.size() ){
    LOG_WARNING("BeamEvent") << "Warning! Could not match to Good Particles: " << "\n";
    for( size_t ip = 0; ip < leftOvers.size(); ++ip){
      LOG_WARNING("BeamEvent") << leftOvers[ip] << " ";
    }
    LOG_WARNING("BeamEvent") << "\n";
  }

  for(size_t i = 0; i < beamspill->GetNFBMTriggers(name); ++i){
    beamspill->DecodeFibers(name,i);
  } 

}

////////////////////////
// 
void proto::BeamEvent::parseXBPF(uint64_t time){
  for(size_t d = 0; d < fDevices.size(); ++d){
    std::string name = fDevices[d];
    parseGeneralXBPF(name, time, d);
  }  
}
// END BeamEvent::parseXBFP
////////////////////////


void proto::BeamEvent::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  
  if( fDebugMomentum ){
    fFullMomentum = tfs->make<TH1F>("FullMomentum", "", 150, 0, 15.);
    fCutMomentum = tfs->make<TH1F>("CutMomentum", "", 150, 0, 15.);
  }

  if( fSaveOutTree ){
    fOutTree = tfs->make<TTree>("tree", "lines"); 
    fOutTree->Branch("Time", &eventTime);
    fOutTree->Branch("RDTS", &RDTSTime);
    fOutTree->Branch("Event", &eventNum);
    fOutTree->Branch("SpillStart", &SpillStart);
    fOutTree->Branch("SpillEnd", &SpillEnd);
    fOutTree->Branch("SpillOffset", &SpillOffset);
    fOutTree->Branch("ActiveTriggerTime", &ActiveTriggerTime);
    fOutTree->Branch("Run",   &runNum);
    fOutTree->Branch("Subrun", &subRunNum);
    fOutTree->Branch("Pressure1", &CKov1Pressure);
    fOutTree->Branch("Eff1", &CKov1Efficiency);
    fOutTree->Branch("Pressure2", &CKov2Pressure);
    fOutTree->Branch("Eff2", &CKov2Efficiency);
    fOutTree->Branch("acqTime",        &acqTime);
    fOutTree->Branch("acqStampMBPL",   &acqStampMBPL);
    fOutTree->Branch("HLTWord",     &HLTWord);
    fOutTree->Branch("HLTTS",       &HLTTS);
    fOutTree->Branch("BeamOn",      &BeamOn);
    fOutTree->Branch("BITrigger",   &BITrigger);
    fOutTree->Branch("Upstream",    &Upstream);
    fOutTree->Branch("C1",          &C1);
    fOutTree->Branch("C2",          &C2);
    fOutTree->Branch("BP1",         &BP1);
    fOutTree->Branch("BP2",         &BP2);
    fOutTree->Branch("BP3",         &BP3);
    fOutTree->Branch("BP4",         &BP4);
    fOutTree->Branch("s11Sec",      &s11Sec);
    fOutTree->Branch("s11Nano",      &s11Nano);
  }    

  diff2A = new std::vector<double>();
  diff2B = new std::vector<double>();
  if( fDebugTOFs ){
    fGenTrigTree = tfs->make<TTree>("GenTrig","");
    fGenTrigTree->Branch("coarse", &fGenTrigCoarse);
    fGenTrigTree->Branch("frac", &fGenTrigFrac);
    fGenTrigTree->Branch("sec", &fGenTrigSec);
  
    fXTOF1ATree = tfs->make<TTree>("TOF1A","");
    fXTOF1ATree->Branch("coarse", &fXTOF1ACoarse);
    fXTOF1ATree->Branch("frac", &fXTOF1AFrac);
    fXTOF1ATree->Branch("sec", &fXTOF1ASec);
  
    fXTOF1BTree = tfs->make<TTree>("TOF1B","");
    fXTOF1BTree->Branch("coarse", &fXTOF1BCoarse);
    fXTOF1BTree->Branch("frac", &fXTOF1BFrac);
    fXTOF1BTree->Branch("sec", &fXTOF1BSec);
  
    fXTOF2ATree = tfs->make<TTree>("TOF2A","");
    fXTOF2ATree->Branch("coarse", &fXTOF2ACoarse);
    fXTOF2ATree->Branch("frac", &fXTOF2AFrac);
    fXTOF2ATree->Branch("sec", &fXTOF2ASec);
    fXTOF2ATree->Branch("diff2A", &diff2A);
  
    fXTOF2BTree = tfs->make<TTree>("TOF2B","");
    fXTOF2BTree->Branch("coarse", &fXTOF2BCoarse);
    fXTOF2BTree->Branch("frac", &fXTOF2BFrac);
    fXTOF2BTree->Branch("sec", &fXTOF2BSec);
    fXTOF2BTree->Branch("diff2B", &diff2B);
  }

  // Read in and cache the beam bundle folder for a specific time

  ifbeam_ns::BeamFolder::_debug = fIFBeamDebug;

  bfp = ifb->getBeamFolder(fBundleName,fURLStr,fTimeWindow);
  LOG_INFO("BeamEvent") << "%%%%%%%%%% Got beam folder %%%%%%%%%%\n"; 

  bfp->set_epsilon( fBFEpsilon );
  LOG_INFO("BeamEvent") << "%%%%%%%%%% Set beam epislon " << fBFEpsilon << " %%%%%%%%%%\n";


 
}

void proto::BeamEvent::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::endJob()
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void proto::BeamEvent::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fBundleName  = p.get<std::string>("BundleName");
  fOutputLabel = p.get<std::string>("OutputLabel");
  fURLStr      = p.get<std::string>("URLStr");
  fBFEpsilon   = p.get<double>("BFEpsilon");
  fIFBeamDebug = p.get<int>("IFBeamDebug");
  fNRetries    = p.get<int>("NRetries");
  fTimeWindow  = p.get<double>("TimeWindow");
  fFixedTime   = p.get<uint64_t>("FixedTime");
  //fMultipleTimes = p.get< std::vector<uint64_t> >("MultipleTimes");

  fDevices     = p.get< std::vector< std::string > >("Devices");    

  //For Tracking/////
  firstUpstreamName    = p.get< std::string >("FirstUpstream");
  secondUpstreamName   = p.get< std::string >("SecondUpstream");
  firstDownstreamName  = p.get< std::string >("FirstDownstream");
  secondDownstreamName = p.get< std::string >("SecondDownstream");
  ///////////////////

  //For Momentum Spectrometry////
  firstBPROF1  = p.get< std::string >("FirstBPROF1");
  secondBPROF1 = p.get< std::string >("SecondBPROF1");
  BPROF2       = p.get< std::string >("BPROF2");
  BPROF3       = p.get< std::string >("BPROF3");
  fBeamBend    = p.get< double >("BeamBend");
/*  L1           = 2.004;//(m)
  L2           = 1.718*cos(fBeamBend);//(m)
  L3           = 2.728*cos(fBeamBend);//(m)
*/
  magnetLen    = 1.;//(m)
  magnetField  = 1000.;//()
  ///////////////////////////////

  std::vector< std::pair<std::string, std::string> >  tempTypes = p.get<std::vector< std::pair<std::string, std::string> >>("DeviceTypes");
  fDeviceTypes  = std::map<std::string, std::string>(tempTypes.begin(), tempTypes.end() );

  //Deminsion of Fibers
  std::vector< std::pair<std::string, double> > tempFiberDims = p.get< std::vector<std::pair<std::string,double> > >("Dimension");
  fFiberDimension = std::map<std::string, double>(tempFiberDims.begin(), tempFiberDims.end());

  
  //XTOF devices 
  fTOF1 = p.get< std::string >("TOF1");
  fTOF2 = p.get< std::string >("TOF2");

  fTOF1A = fTOF1 + "A";
  fTOF1B = fTOF1 + "B";
  fTOF2A = fTOF2 + "A";
  fTOF2B = fTOF2 + "B";

  fCKov1 = p.get< std::string >("CKov1");
  fCKov2 = p.get< std::string >("CKov2");


  fXBPFPrefix      = p.get<std::string>("XBPFPrefix");
  fXTOFPrefix      = p.get<std::string>("XTOFPrefix");
  fXCETPrefix      = p.get<std::string>("XCETPrefix");


  //New parameters to match Leigh's
  fRotateMonitorXZ = p.get<double>("RotateMonitorXZ");
  fRotateMonitorYZ = p.get<double>("RotateMonitorYZ");

  fFirstTrackingProfZ  = p.get<double>("FirstTrackingProfZ");
  fSecondTrackingProfZ = p.get<double>("SecondTrackingProfZ");
  fNP04FrontZ          = p.get<double>("NP04FrontZ");  

  fBeamX               = p.get<double>("BeamX");
  fBeamY               = p.get<double>("BeamY");
  fBeamZ               = p.get<double>("BeamZ");

  fForceNewFetch       = p.get<bool>("ForceNewFetch");
  fMatchTime           = p.get<bool>("MatchTime");
  fForceRead           = p.get<bool>("ForceRead");


  fTimingCalibration      = p.get<double>("TimingCalibration");
  fCalibrationTolerance   = p.get<double>("CalibrationTolerance");
  fOffsetTAI              = p.get<double>("OffsetTAI");

  fFetchOffset            = p.get<int>("FetchOffset");

  fS11DiffUpper           = p.get<double>("S11DiffUpper");
  fS11DiffLower           = p.get<double>("S11DiffLower");

  fRDTSToS11Upper         = p.get<double>("RDTSToS11Upper");
  fRDTSToS11Lower         = p.get<double>("RDTSToS11Lower");

  fOffsetCTBtoRDTS        = p.get<int>("OffsetCTBtoRDTS");
  fToleranceCTBtoRDTS     = p.get<int>("ToleranceCTBtoRDTS");

  fSaveOutTree            = p.get<bool>("SaveOutTree");
  fDebugMomentum          = p.get<bool>("DebugMomentum");
  fDebugTOFs              = p.get<bool>("DebugTOFs");

}

uint64_t proto::BeamEvent::joinHighLow(double high, double low){

  uint64_t low64 = (uint64_t)low;

  uint64_t high64 = (uint64_t)high;

  high64 = high64 << 32;
  uint64_t joined = high64 | low64;

  return joined; 
}

TVector3 proto::BeamEvent::ConvertProfCoordinates(double x, double y, double z, double zOffset){
  double off = fNP04FrontZ - zOffset;

  TVector3 old(x,y,z);

  double newX = x*fBMBasisX.X() + y*fBMBasisY.X() + /*(z-zOffset)*fBMBasisZ.X()*/ + off*fabs(fBMBasisZ.X());
  double newY = x*fBMBasisX.Y() + y*fBMBasisY.Y() + /*(z-zOffset)*fBMBasisZ.Y()*/ + off*fabs(fBMBasisZ.Y());
  double newZ = x*fBMBasisX.Z() + y*fBMBasisY.Z() + /*(z-zOffset)              */ - off*fabs(fBMBasisZ.Z());

  newX += fBeamX*10.;
  newY += fBeamY*10.;
  newZ += fBeamZ*10.;

  TVector3 result(newX/10., newY/10., newZ/10.);
  return result;
}

void proto::BeamEvent::BeamMonitorBasisVectors(){
  fBMBasisX = TVector3(1.,0.,0.);
  fBMBasisY = TVector3(0.,1.,0.);
  fBMBasisZ = TVector3(0.,0.,1.);
  RotateMonitorVector(fBMBasisX);
  RotateMonitorVector(fBMBasisY);
  RotateMonitorVector(fBMBasisZ);

  rotated = true;
}

void proto::BeamEvent::RotateMonitorVector(TVector3 &vec){
  vec.RotateY(fRotateMonitorXZ * TMath::Pi()/180.);
  vec.RotateX(fRotateMonitorYZ * TMath::Pi()/180.);
}

void proto::BeamEvent::MakeTrack(size_t theTrigger){
  
  LOG_DEBUG("BeamEvent") << "Making Track for time: " << beamspill->GetT0Sec(theTrigger) << " " << beamspill->GetT0Nano(theTrigger) << "\n";

  //Get the active fibers from the upstream tracking XBPF
  std::vector<short> firstUpstreamFibers  = beamspill->GetActiveFibers(firstUpstreamName, theTrigger);
  std::vector<short> secondUpstreamFibers = beamspill->GetActiveFibers(secondUpstreamName, theTrigger);

  LOG_DEBUG("BeamEvent") << firstUpstreamName << " has " << firstUpstreamFibers.size() << " active fibers"<< "\n";
  for(size_t i = 0; i < firstUpstreamFibers.size(); ++i){
    LOG_DEBUG("BeamEvent") << firstUpstreamFibers[i] << " ";
  }
  LOG_DEBUG("BeamEvent") << "\n";

  LOG_DEBUG("BeamEvent") << secondUpstreamName << " has " << secondUpstreamFibers.size() << " active fibers" << "\n";
  for(size_t i = 0; i < secondUpstreamFibers.size(); ++i){
    LOG_DEBUG("BeamEvent") << secondUpstreamFibers[i] << " ";
  }
  LOG_DEBUG("BeamEvent") << "\n";
  //////////////////////////////////////////////

  //Get the active fibers from the downstream tracking XBPF
  std::vector<short> firstDownstreamFibers = beamspill->GetActiveFibers(firstDownstreamName, theTrigger);
  std::vector<short> secondDownstreamFibers = beamspill->GetActiveFibers(secondDownstreamName, theTrigger);

  LOG_DEBUG("BeamEvent") << firstDownstreamName << " has " << firstDownstreamFibers.size() << " active fibers" << "\n";
  for(size_t i = 0; i < firstDownstreamFibers.size(); ++i){
    LOG_DEBUG("BeamEvent") << firstDownstreamFibers[i] << " ";
  }
  LOG_DEBUG("BeamEvent") << "\n";

  LOG_DEBUG("BeamEvent") << secondDownstreamName << " has " << secondDownstreamFibers.size() << " active fibers" << "\n";
  for(size_t i = 0; i < secondDownstreamFibers.size(); ++i){
    LOG_DEBUG("BeamEvent") << secondDownstreamFibers[i] << " ";
  }
  LOG_DEBUG("BeamEvent") << "\n";
  //////////////////////////////////////////////

  if( (firstUpstreamFibers.size() < 1) || (secondUpstreamFibers.size() < 1) || (firstDownstreamFibers.size() < 1) || (secondDownstreamFibers.size() < 1) ){
    LOG_DEBUG("BeamEvent") << "Warning, at least one empty Beam Profiler. Not making track" << "\n";
    return;
  }
/*  else if( (firstUpstreamFibers.size() > 5) || (secondUpstreamFibers.size() > 5) || (firstDownstreamFibers.size() > 5) || (secondDownstreamFibers.size() > 5) ){
    LOG_DEBUG("BeamEvent") << "Warning, too many (>5) active fibers in at least one Beam Profiler. Not making track" << "\n";
    return;
  }*/

  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  std::vector< std::pair<size_t, size_t> > upstreamPairedFibers;
  std::vector< std::pair<size_t, size_t> > downstreamPairedFibers;
  std::string firstUpstreamType    = fDeviceTypes[firstUpstreamName]; 
  std::string secondUpstreamType   = fDeviceTypes[secondUpstreamName]; 
  std::string firstDownstreamType  = fDeviceTypes[firstDownstreamName]; 
  std::string secondDownstreamType = fDeviceTypes[secondDownstreamName]; 

  std::vector< TVector3 > upstreamPositions;
  std::vector< TVector3 > downstreamPositions; 
  
  //Pair the upstream fibers together
  LOG_DEBUG("BeamEvent") << "Upstream" << "\n";
  for(size_t iF1 = 0; iF1 < firstUpstreamFibers.size(); ++iF1){
    
    size_t firstFiber = firstUpstreamFibers[iF1];

    for(size_t iF2 = 0; iF2 < secondUpstreamFibers.size(); ++iF2){
      size_t secondFiber = secondUpstreamFibers[iF2];

      LOG_DEBUG("BeamEvent") << "Paired: " << firstFiber << " " << secondFiber << "\n"; 
      upstreamPairedFibers.push_back(std::make_pair(firstFiber, secondFiber));

      if (iF2 < secondUpstreamFibers.size() - 1){
        if (secondUpstreamFibers[iF2] == (secondUpstreamFibers[iF2 + 1] - 1)) ++iF2;
      }
    }

    if (iF1 < firstUpstreamFibers.size() - 1){
      if (firstUpstreamFibers[iF1] == (firstUpstreamFibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  for(size_t iF = 0; iF < upstreamPairedFibers.size(); ++iF){
    
    std::pair<size_t,size_t> thePair = upstreamPairedFibers.at(iF);
 
    if(firstUpstreamType == "horiz" && secondUpstreamType == "vert"){
      double xPos = GetPosition(firstUpstreamName, thePair.first);
      double yPos = GetPosition(secondUpstreamName, thePair.second);
      
      LOG_DEBUG("BeamEvent") << "normal " << xPos << " " << yPos <<  "\n";
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
      LOG_DEBUG("BeamEvent") << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << "\n";
      upstreamPositions.push_back( posInDet );
    }
    else if(firstUpstreamType == "vert" && secondUpstreamType == "horiz"){
      double yPos = GetPosition(firstUpstreamName, thePair.first);
      double xPos = GetPosition(secondUpstreamName, thePair.second);
      LOG_DEBUG("BeamEvent") << "normal " << xPos << " " << yPos <<  "\n";

      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fFirstTrackingProfZ);
      LOG_DEBUG("BeamEvent") << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << "\n";
      upstreamPositions.push_back( posInDet );
    }

  }
  
  LOG_DEBUG("BeamEvent") << "Downstream" << "\n";
  //Pair the downstream fibers together
  for(size_t iF1 = 0; iF1 < firstDownstreamFibers.size(); ++iF1){
    
    size_t firstFiber = firstDownstreamFibers[iF1];

    for(size_t iF2 = 0; iF2 < secondDownstreamFibers.size(); ++iF2){
      size_t secondFiber = secondDownstreamFibers[iF2];

      LOG_DEBUG("BeamEvent") << "Paired: " << firstFiber << " " << secondFiber << "\n"; 
      downstreamPairedFibers.push_back(std::make_pair(firstFiber, secondFiber));

      if (iF2 < secondDownstreamFibers.size() - 1){
        if (secondDownstreamFibers[iF2] == (secondDownstreamFibers[iF2 + 1] - 1)) ++iF2;
      }
    }

    if (iF1 < firstDownstreamFibers.size() - 1){
      if (firstDownstreamFibers[iF1] == (firstDownstreamFibers[iF1 + 1] - 1)) ++iF1;
    }
  }

  for(size_t iF = 0; iF < downstreamPairedFibers.size(); ++iF){
    
    std::pair<size_t,size_t> thePair = downstreamPairedFibers.at(iF);
 
    if(firstDownstreamType == "horiz" && secondDownstreamType == "vert"){
      double xPos = GetPosition(firstDownstreamName, thePair.first);
      double yPos = GetPosition(secondDownstreamName, thePair.second);

      LOG_DEBUG("BeamEvent") << "normal " << xPos << " " << yPos <<  "\n";
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
      LOG_DEBUG("BeamEvent") << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << "\n";
      downstreamPositions.push_back( posInDet );
    }
    else if(firstDownstreamType == "vert" && secondDownstreamType == "horiz"){
      double yPos = GetPosition(firstDownstreamName, thePair.first);
      double xPos = GetPosition(secondDownstreamName, thePair.second);

      LOG_DEBUG("BeamEvent") << "normal " << xPos << " " << yPos <<  "\n";
      TVector3 posInDet = ConvertProfCoordinates(xPos,yPos,0.,fSecondTrackingProfZ);
      LOG_DEBUG("BeamEvent") << posInDet.X() << " " << posInDet.Y() << " " << posInDet.Z() << "\n";
      downstreamPositions.push_back( posInDet );
    }

  }
 
  //Just for creating tracks
  std::vector< std::vector<double> > dummy;
  std::vector<double> mom(3, util::kBogusD);
  /// 

  for(size_t iU = 0; iU < upstreamPositions.size(); ++iU){
    for(size_t iD = 0; iD < downstreamPositions.size(); ++iD){
      std::vector<TVector3> thePoints;
      thePoints.push_back(upstreamPositions.at(iU));
      thePoints.push_back(downstreamPositions.at(iD));

      //Now project the last point to the TPC face
      thePoints.push_back( ProjectToTPC(thePoints[0],thePoints[1]) );    

      
      std::vector<TVector3> theMomenta;
      //Just push back the unit vector for each point 
      //For now.
      //Eventually, use momentum from curvature?
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );
      theMomenta.push_back( ( downstreamPositions.at(iD) - upstreamPositions.at(iU) ).Unit() );

      recob::Track * tempTrack = new recob::Track(thePoints, theMomenta, dummy, mom, 1);      
      beamevt->AddBeamTrack( *tempTrack );
    }
  }

}

void proto::BeamEvent::MomentumSpec(size_t theTrigger){
  
  LOG_INFO("BeamEvent") << "Doing momentum spectrometry for trigger " << beamspill->GetT0Sec(theTrigger) << " " << beamspill->GetT0Nano(theTrigger) << "\n";

  double LB = mag_P1*fabs(current[0]);
  double deltaI = fabs(current[0]) - mag_P4;
  if(deltaI>0) LB+= mag_P3*deltaI*deltaI;

  //Get the active fibers from the upstream tracking XBPF
  std::string firstBPROF1Type    = fDeviceTypes[firstBPROF1]; 
  std::string secondBPROF1Type   = fDeviceTypes[secondBPROF1]; 
  std::vector<short> BPROF1Fibers;
  std::string BPROF1Name;

  if (firstBPROF1Type == "horiz" && secondBPROF1Type == "vert"){
    BPROF1Fibers = beamspill->GetActiveFibers(firstBPROF1, theTrigger);
    BPROF1Name = firstBPROF1;

    LOG_INFO("BeamEvent") << firstBPROF1 << " has " << BPROF1Fibers.size() << " active fibers" << "\n";
    for(size_t i = 0; i < BPROF1Fibers.size(); ++i){
      LOG_INFO("BeamEvent") << BPROF1Fibers[i] << " ";
    }
    LOG_INFO("BeamEvent") << "\n";

  }
  else if(secondBPROF1Type == "horiz" && firstBPROF1Type == "vert"){
    BPROF1Fibers = beamspill->GetActiveFibers(secondBPROF1, theTrigger);
    BPROF1Name = secondBPROF1;

    LOG_INFO("BeamEvent") << secondBPROF1 << " has " << BPROF1Fibers.size() << " active fibers" << "\n";
    for(size_t i = 0; i < BPROF1Fibers.size(); ++i){
      LOG_INFO("BeamEvent") << BPROF1Fibers[i] << " ";
    }
    LOG_INFO("BeamEvent") << "\n";
  }
  else{
    LOG_INFO("BeamEvent") << "Error: type is not correct" << "\n";
    return;
  }


  //////////////////////////////////////////////

  if( (BPROF1Fibers.size() < 1) ){
    LOG_INFO("BeamEvent") << "Warning, at least one empty Beam Profiler. Not checking momentum" << "\n";
    return;
  }
  //We have the active Fibers, now go through them.
  //Skip the second of any adjacents 
  std::vector< short > strippedFibers; 
  LOG_INFO("BeamEvent") << "BPROF1" << "\n";
  for(size_t iF1 = 0; iF1 < BPROF1Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF1Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[iF1] == (BPROF1Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
  
  double X1, X2, X3;

  //X-direction convention is reversed for the spectrometer   
  X1 = -1.*GetPosition(BPROF1Name, strippedFibers[0])/1.E3;

  //BPROF2////
  //
  std::vector<short>  BPROF2Fibers = beamspill->GetActiveFibers(BPROF2, theTrigger);
  if( (BPROF2Fibers.size() < 1) ){
    LOG_INFO("BeamEvent") << "Warning, at least one empty Beam Profiler. Not checking momentum" << "\n";
    return;
  }
  LOG_INFO("BeamEvent") << BPROF2 << " has " << BPROF2Fibers.size() << " active fibers at time " << beamspill->GetFiberTime(BPROF2,theTrigger) << "\n";
  for(size_t i = 0; i < BPROF2Fibers.size(); ++i){
    LOG_INFO("BeamEvent") << BPROF2Fibers[i] << " ";
  }
  LOG_INFO("BeamEvent") << "\n";

  strippedFibers.clear(); 
  LOG_INFO("BeamEvent") << "BPROF2" << "\n";
  for(size_t iF1 = 0; iF1 < BPROF2Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF2Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF2Fibers.size() - 1){
      if (BPROF2Fibers[iF1] == (BPROF2Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  X2 = -1.*GetPosition(BPROF2, strippedFibers[0])/1.E3; 
  //X2 = X2*cos(TMath::Pi() - fBeamBend/1000.);
  ////////////

  //BPROF3////
  //
  std::vector<short>  BPROF3Fibers = beamspill->GetActiveFibers(BPROF3, theTrigger);
  if( (BPROF3Fibers.size() < 1) ){
    LOG_INFO("BeamEvent") << "Warning, at least one empty Beam Profiler. Not checking momentum" << "\n";
    return;
  }
  LOG_INFO("BeamEvent") << BPROF3 << " has " << BPROF3Fibers.size() << " active fibers at time " << beamspill->GetFiberTime(BPROF3,theTrigger) << "\n";
  for(size_t i = 0; i < BPROF3Fibers.size(); ++i){
    LOG_INFO("BeamEvent") << BPROF3Fibers[i] << " ";
  }
  LOG_INFO("BeamEvent") << "\n";

  strippedFibers.clear(); 
  LOG_INFO("BeamEvent") << "BPROF3" << "\n";
  for(size_t iF1 = 0; iF1 < BPROF3Fibers.size(); ++iF1){
    
    size_t Fiber = BPROF3Fibers[iF1];
    strippedFibers.push_back(Fiber);

    if (iF1 < BPROF3Fibers.size() - 1){
      if (BPROF3Fibers[iF1] == (BPROF3Fibers[iF1 + 1] - 1)) ++iF1;
    }
  }
 
  X3 = -1.*GetPosition(BPROF3, strippedFibers[0])/1.E3; 
  ////////////
  
  double cosTheta = MomentumCosTheta(X1,X2,X3);
  double momentum = 299792458*LB/(1.E9 * acos(cosTheta));
  beamevt->AddRecoBeamMomentum(momentum);
  LOG_INFO("BeamEvent") << "Momentum: " << momentum << "\n"; 

  if( (BPROF1Fibers.size() == 1) && (BPROF2Fibers.size() == 1) && (BPROF3Fibers.size() == 1) ){
    LOG_INFO("BeamEvent") << "Filling Cut Momentum Spectrum" << "\n";
    if( fDebugMomentum ) fCutMomentum->Fill(momentum);
  }


  LOG_INFO("BeamEvent") << "Getting all trio-wise hits" << "\n";
  LOG_INFO("BeamEvent") << "N1,N2,N3 " << BPROF1Fibers.size()
            << " "         << BPROF2Fibers.size() 
            << " "         << BPROF3Fibers.size() << "\n";
  for(size_t i1 = 0; i1 < BPROF1Fibers.size(); ++i1){
    double x1,x2,x3;

    x1 = -1.*GetPosition(BPROF1Name, BPROF1Fibers[i1])/1.E3;
    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Add .5 mm
        x1 += .0005;
      }
    }

    for(size_t i2 = 0; i2 < BPROF2Fibers.size(); ++i2){
      x2 = -1.*GetPosition(BPROF2, BPROF2Fibers[i2])/1.E3;
      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Add .5 mm
          x2 += .0005;
        }
      }

      for(size_t i3 = 0; i3 < BPROF3Fibers.size(); ++i3){
        LOG_INFO("BeamEvent") << "\t" << i1 << " " << i2 << " " << i3 << "\n";
        x3 = -1.*GetPosition(BPROF3, BPROF3Fibers[i3])/1.E3;
        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Add .5 mm
            x3 += .0005;
          }
        }

        double cosTheta_full = MomentumCosTheta(x1,x2,x3);        
        double momentum_full = 299792458*LB/(1.E9 * acos(cosTheta_full));
        if( fDebugMomentum ) fFullMomentum->Fill(momentum_full);

        //NEED TO MOVE BEAMEVENT HERE EVENTUALLY

        if (i3 < BPROF3Fibers.size() - 1){
          if (BPROF3Fibers[i3] == (BPROF3Fibers[i3 + 1] - 1)){
            //Skip the next
            ++i3;
          }
        }
        
      }

      if (i2 < BPROF2Fibers.size() - 1){
        if (BPROF2Fibers[i2] == (BPROF2Fibers[i2 + 1] - 1)){
          //Skip the next
          ++i2;
        }
      }
    }

    if (i1 < BPROF1Fibers.size() - 1){
      if (BPROF1Fibers[i1] == (BPROF1Fibers[i1 + 1] - 1)){
        //Skip the next
        ++i1;
      }
    }
  }
}

double proto::BeamEvent::MomentumCosTheta(double X1, double X2, double X3){
/*
  double a = cos(fBeamBend)*(X3*L2 - X2*L3) / (L3 - L2);

  double cosTheta = (a + X1)*( (L3 - L2)*tan(fBeamBend) + (X2 - X3)*cos(fBeamBend) )+ (L3 - L2)*L1 ;

  double denomTerm1, denomTerm2, denom; 
  denomTerm1 = sqrt( L1*L1 + (a + X1)*(a + X1) );
  denomTerm2 = sqrt( 
               (L3 - L2)*(L3 - L2) 
             + TMath::Power( ( (L3 - L2)*tan(fBeamBend) 
                             + (X2 - X3)*cos(fBeamBend) ), 2 ) );
                            
  denom = denomTerm1*denomTerm2;
  cosTheta = cosTheta / denom;

  //
 
  double a = (X2*L3 - X3*L2) / (L3 - L2);
   
  double cosTheta = (a - X1)*(X3 - X2)*cos(fBeamBend) + (L3 - L2 )*( L1 + (a - X1)*tan(fBeamBend) );
    
  double denomTerm1, denomTerm2, denom; 
  denomTerm1 = sqrt( L1*L1 + (a - X1)*(a - X1) );
  denomTerm2 = sqrt( 
               (L3 - L2)*(L3 - L2) 
             + TMath::Power( ( (L3 - L2)*tan(fBeamBend) 
                             + (X3 - X2)*cos(fBeamBend) ), 2 ) );
  
  denom = denomTerm1*denomTerm2;
  cosTheta = cosTheta / denom;
 */

    //double a = L2*tan(fBeamBend) + X2*cos(fBeamBend) - ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) )*( L2 - X2*sin(fBeamBend) );
    double a =  ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) )*( L2 - X2*sin(fBeamBend) );
    a = a / (L3 - X3*sin(fBeamBend) - L2 + X2*sin(fBeamBend) );
    a = L2*tan(fBeamBend) + X2*cos(fBeamBend) - a;
  
    double numTerm = (a - X1)*(L3 - L2)*tan(fBeamBend) + (a - X1)*(X3 - X2)*cos(fBeamBend) + L1*( (L3 - L2) - (X3 - X2)*sin(fBeamBend) );
  
    double denomTerm1, denomTerm2, denom;
    denomTerm1 = sqrt( L1*L1 + (a - X1)*(a - X1) );
    denomTerm2 = sqrt( TMath::Power( ( (L3 - L2)*tan(fBeamBend) + (X3 - X2)*cos(fBeamBend) ),2)
                     + TMath::Power( ( (L3 - L2)                - (X3 - X2)*sin(fBeamBend) ),2) );
    denom = denomTerm1 * denomTerm2; 

    double cosTheta = numTerm/denom;
  
  return cosTheta;
}

TVector3 proto::BeamEvent::ProjectToTPC(TVector3 firstPoint, TVector3 secondPoint){
  TVector3 dR = (secondPoint - firstPoint);
  
  double deltaZ = -1.*secondPoint.Z();
  double deltaX = deltaZ * (dR.X() / dR.Z());
  double deltaY = deltaZ * (dR.Y() / dR.Z());

  TVector3 lastPoint = secondPoint + TVector3(deltaX, deltaY, deltaZ);
  return lastPoint;
}

double proto::BeamEvent::GetPosition(std::string deviceName, int fiberIdx){
  //NEEDS WORK
  if(fiberIdx > 192){ LOG_WARNING("BeamEvent") << "Please select fiber in range [0,191]" << "\n"; return -1.;}
  double size = fFiberDimension[deviceName];
  //double size = 1.;
  
  //Define 0th fiber as farthest positive. Last fiber is farthest negative. Center is between 96 and 97 
  double pos = size*(96 - fiberIdx) - size/2.;
  return pos;
}


DEFINE_ART_MODULE(proto::BeamEvent)
