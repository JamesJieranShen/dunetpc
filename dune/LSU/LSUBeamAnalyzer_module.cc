////////////////////////////////////////////////////////////////////////
// Class:       LSUBeamAnalyzer
// Module Type: analyzer
// File:        LSUBeamAnalyzer_module.cc
//
// Generated at Mon Feb  1 13:51:14 2016 by Andrew Olivier using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/FileBlock.h"

//LArSoft includes
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/TriggerData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "dunetpc/dune/LSU/TrajectoryInterpExtrapAlg.h"
#include "dunetpc/dune/LSU/MCBeamOrCosmicAlg.h"
#include "dunetpc/dune/LSU/PlaneIntersectionFinder.h"
#include "dunetpc/dune/LSU/MotherDaughterWalkerAlg.h"
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEDataUtils.h"
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "dunetpc/dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "larsim/MCCheater/BackTrackerService.h"

//ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"

//c++ includes
#include <string>
#include <iomanip>
#include <algorithm> //for std::sort
#include <cmath>

#define DEFAULTNEG -99999999
#define DEFAULTPOS 99999999
#define MAXTOFS 10
#define MAXTRACKS 1000
#define MAXDAUGHTER 25
#define MAXMCPARTS 10000
#define MAXBEAMTRACKS 300
#define MAXBEAMMOMS 300
#define MAXIDES 150000
#define MAXPFBEAMPRIM 10
#define MAXZINT 95
#define MAXLINT 20

#define MCHARGEDPION 139.57018 // MeV/c^2
#define MPROTON 938.2720813 // MeV/c^2
#define KINLOSTBEFORETPC 0.0 //MeV; from LArIAT pion total cross section group
#define KINLOSTBEFORETPCPROTON 0.0 //MeV; from LArIAT pion total cross section group

#define setHistTitles(hist,xtitle,ytitle) hist->GetXaxis()->SetTitle(xtitle); hist->GetYaxis()->SetTitle(ytitle);

namespace lana {
  class LSUBeamAnalyzer;
}

class lana::LSUBeamAnalyzer : public art::EDAnalyzer {
public:
  explicit LSUBeamAnalyzer(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LSUBeamAnalyzer(LSUBeamAnalyzer const &) = delete;
  LSUBeamAnalyzer(LSUBeamAnalyzer &&) = delete;
  LSUBeamAnalyzer & operator = (LSUBeamAnalyzer const &) = delete;
  LSUBeamAnalyzer & operator = (LSUBeamAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run const & r) override;
  void beginSubRun(art::SubRun const & sr) override;
  void endJob() override;
  void endRun(art::Run const & r) override;
  void endSubRun(art::SubRun const & sr) override;
  void reconfigure(fhicl::ParameterSet const & p);
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  // Declare member data here.

  //parameters read from fcl file
  art::InputTag fBeamEventTag;
  art::InputTag fRawTriggerTag;

  /////////////////////////////
  //Tree variables
  TTree* tree;

  bool isMC;
  UInt_t runNumber;
  UInt_t subRunNumber;
  UInt_t eventNumber;
  std::string infilename;
  Float_t xWC; // WC track projected to TPC front face, x coord in cm
  Float_t yWC; // WC track projected to TPC front face, y coord in cm
  Float_t thetaWC; // WC track theta
  Float_t phiWC; // WC track phi
  Float_t pzWC; // WC track momentum (z component) MeV/c
  Float_t pWC; // WC track momentum MeV/c
  Float_t eWC; // WC track energy assuming charged pion mass MeV
  Float_t kinWC; // WC track kinetic energy assuming charged pion mass MeV
  Float_t kinWCInTPC; // WC track kinetic energy assuming charged pion mass, and subtracting constant E MeV
  Float_t eWCProton; // WC track energy assuming proton mass MeV
  Float_t kinWCProton; // WC track kinetic energy assuming proton mass MeV
  Float_t kinWCInTPCProton; // WC track kinetic energy assuming proton mass, and subtracting constant E MeV
  Float_t yKinkWC; // WC track y-kink
  UInt_t nHitsWC; // WC track n hits
  Float_t xWC4Hit; // WC4 hit position x coord in cm
  Float_t yWC4Hit; // WC4 hit position y coord in cm
  Float_t zWC4Hit; // WC4 hit position z coord in cm

  UInt_t nBeamTracks;
  Float_t beamTrackXFrontTPC[MAXBEAMTRACKS];
  Float_t beamTrackYFrontTPC[MAXBEAMTRACKS];
  Float_t beamTrackTheta[MAXBEAMTRACKS];
  Float_t beamTrackPhi[MAXBEAMTRACKS];
  Float_t beamTrackMom[MAXBEAMTRACKS];
  Float_t beamTrackEPion[MAXBEAMTRACKS];
  Float_t beamTrackEProton[MAXBEAMTRACKS];
  Float_t beamTrackKinPion[MAXBEAMTRACKS];
  Float_t beamTrackKinProton[MAXBEAMTRACKS];
  Float_t beamTrackTruePDG[MAXBEAMTRACKS];

  UInt_t nBeamMom;
  Float_t beamMom[MAXBEAMMOMS];

  UInt_t nBeamEvents;
  Int_t BITrigger;
  bool BITriggerMatched;
  UInt_t BIActiveTrigger;
  Float_t TOF;
  Int_t TOFChan;
  UInt_t nTOFs;
  float TOFs[MAXTOFS];
  Int_t TOFChans[MAXTOFS];
  UInt_t TOFusTrigs[MAXTOFS];
  UInt_t TOFdsTrigs[MAXTOFS];

  Int_t CKov0Status;
  Int_t CKov1Status;
  Float_t CKov0Time;
  Float_t CKov1Time;
  Float_t CKov0Pressure;
  Float_t CKov1Pressure;

  bool triggerIsBeam;
  UInt_t triggerBits;

  protoana::ProtoDUNEDataUtils fDataUtil;
 
  void ResetTreeVars();
};

lana::LSUBeamAnalyzer::LSUBeamAnalyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  infilename("NoInFilenameFound"),
  fDataUtil(p.get<fhicl::ParameterSet>("DataUtil"))
 // More initializers here.
{
  this->reconfigure(p);
}

void lana::LSUBeamAnalyzer::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  ResetTreeVars();

  art::ServiceHandle<geo::Geometry> geom;

  runNumber = e.run();
  subRunNumber = e.subRun();
  eventNumber = e.id().event();
  isMC = !(e.isRealData());

  // Triggers
  art::Handle< std::vector<raw::Trigger> > triggerHandle;
  if(e.isRealData())
  {
    e.getByLabel(fRawTriggerTag,triggerHandle);
  }

  // Beamline info
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if(e.isRealData())
  {
    auto beamHand = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
    if(beamHand.isValid())
    {
      art::fill_ptr_vector(beamVec, beamHand);
    }
  }

  const bool PRINTBEAMEVENT=true;
  nBeamEvents = beamVec.size();
  nBeamTracks = 0;
  nBeamMom = 0;
  for(size_t iBeamEvent=0; iBeamEvent < beamVec.size(); iBeamEvent++)
  {
    beam::ProtoDUNEBeamEvent beamEvent = *(beamVec.at(iBeamEvent));
    //art::Ptr<beam::ProtoDUNEBeamEvent> beamEvent = beamVec.at(iBeamEvent);
    if(PRINTBEAMEVENT)
    {
      std::cout << "PiAbsSel BeamEvent " << iBeamEvent << ": \n";
      std::cout << "  CTB Timestamp: " << beamEvent.GetCTBTimestamp() << "\n";
      std::cout << "  BI Trigger: " << beamEvent.GetBITrigger() << "\n";
      std::cout << "  Active Trigger: " << beamEvent.GetActiveTrigger() << "\n";
      std::cout << "  Is Trigger Matched: " << beamEvent.CheckIsMatched() << "\n";
      std::cout << "  TOF: " << beamEvent.GetTOF() << "\n";
      std::cout << "  CKov0Status: " << beamEvent.GetCKov0Status() << "\n";
      std::cout << "  CKov1Status: " << beamEvent.GetCKov1Status() << "\n";
      std::cout << "  CKov0Time: " << beamEvent.GetCKov0Time() << "\n";
      std::cout << "  CKov1Time: " << beamEvent.GetCKov1Time() << "\n";
      std::cout << "  CKov0Pressure: " << beamEvent.GetCKov0Pressure() << "\n";
      std::cout << "  CKov1Pressure: " << beamEvent.GetCKov1Pressure() << "\n";
      std::cout << "  Beam Momenta:\n";
    }

    triggerBits = beamEvent.GetTimingTrigger();
    BITrigger = beamEvent.GetBITrigger();
    BITriggerMatched = beamEvent.CheckIsMatched();
    BIActiveTrigger = beamEvent.GetActiveTrigger();
    TOF = beamEvent.GetTOF();
    TOFChan = beamEvent.GetTOFChan();
    CKov0Status = beamEvent.GetCKov0Status();
    CKov1Status = beamEvent.GetCKov1Status();
    CKov0Time = beamEvent.GetCKov0Time();
    CKov1Time = beamEvent.GetCKov1Time();
    CKov0Pressure = beamEvent.GetCKov0Pressure();
    CKov1Pressure = beamEvent.GetCKov1Pressure();
    const bool sameNTracksAsMom = beamEvent.GetNBeamTracks() == beamEvent.GetNRecoBeamMomenta();
    for(size_t iMom=0; iMom < beamEvent.GetNRecoBeamMomenta(); iMom++)
    {
      beamMom[nBeamMom] = beamEvent.GetRecoBeamMomentum(iMom);
      nBeamMom++;
      if(PRINTBEAMEVENT) std::cout << "    " << beamEvent.GetRecoBeamMomentum(iMom) << "\n";
    }
    for(size_t iTrack=0; iTrack < beamEvent.GetNBeamTracks(); iTrack++)
    {
      const recob::Track & track =  beamEvent.GetBeamTrack(iTrack);
      if(nBeamTracks >= MAXBEAMTRACKS)
      {
        throw cet::exception("TooManyBeamTracks","Too many beam tracks in this event, ran out of room in array");
      }
      beamTrackXFrontTPC[nBeamTracks] = track.End().X();
      beamTrackYFrontTPC[nBeamTracks] = track.End().Y();
      beamTrackTheta[nBeamTracks] = track.EndDirection().Theta();
      beamTrackPhi[nBeamTracks] = track.EndDirection().Phi();
      if(sameNTracksAsMom) 
      {
        beamTrackMom[nBeamTracks] = beamEvent.GetRecoBeamMomentum(iTrack);
        beamTrackEPion[nBeamTracks] = sqrt(pow(beamTrackMom[nBeamTracks],2)+pow(MCHARGEDPION,2));
        beamTrackKinPion[nBeamTracks] = beamTrackEPion[nBeamTracks] - MCHARGEDPION;
        beamTrackEProton[nBeamTracks] = sqrt(pow(beamTrackMom[nBeamTracks],2)+pow(MPROTON,2));
        beamTrackKinProton[nBeamTracks] = beamTrackEProton[nBeamTracks] - MPROTON;
      }
      nBeamTracks++;

      xWC = track.End().X();
      yWC = track.End().Y();
      thetaWC = track.EndDirection().Theta();
      phiWC = track.EndDirection().Phi();
    
      if(PRINTBEAMEVENT)
      {
        std::cout << "  Beam Track: "<< iTrack <<"\n";
        std::cout << "    N Points:  " << track.NumberTrajectoryPoints() << "\n";
        std::cout << "    Start Pos: " << track.Vertex().X()
                                  << "  " << track.Vertex().Y()
                                  << "  " << track.Vertex().Z() << "\n";
        std::cout << "    End Pos:   " << track.End().X()
                                  << "  " << track.End().Y()
                                  << "  " << track.End().Z() << "\n";
        std::cout << "    Start Theta: " << track.VertexDirection().Theta()*180/CLHEP::pi << " deg\n";
        std::cout << "    Start Phi:   " << track.VertexDirection().Phi()*180/CLHEP::pi << " deg\n";
        std::cout << "    End Theta:   " << track.EndDirection().Theta()*180/CLHEP::pi << " deg\n";
        std::cout << "    End Theta:   " << track.EndDirection().Phi()*180/CLHEP::pi << " deg\n";
      }
    }

    const auto tofs = beamEvent.GetTOFs(); // double
    const auto tofChans = beamEvent.GetTOFChans(); // int
    const auto usTOFTrigs = beamEvent.GetUpstreamTriggers(); // size_t
    const auto dsTOFTrigs = beamEvent.GetDownstreamTriggers(); // size_t
    nTOFs = tofs.size();
    for(size_t iTOF=0; iTOF < nTOFs; iTOF++)
    {
        std::cout << "iTOF: " << iTOF;
        std::cout << " TOF: " << tofs.at(iTOF);
        TOFs[iTOF] = tofs.at(iTOF);
        if(tofChans.size() > iTOF)
        {
          std::cout << " chan: " << tofChans.at(iTOF);
          TOFChans[iTOF] = tofChans.at(iTOF);
        }
        if(usTOFTrigs.size() > iTOF)
        {
          std::cout << " usTrig: " << usTOFTrigs.at(iTOF);
          TOFusTrigs[iTOF] = usTOFTrigs.at(iTOF);
        }
        if(dsTOFTrigs.size() > iTOF)
        {
          std::cout << " dsTrig: " << dsTOFTrigs.at(iTOF);
          TOFdsTrigs[iTOF] = dsTOFTrigs.at(iTOF);
        }
        std::cout << std::endl;
    }
  }

  // Get Trigger info
  if (triggerHandle.isValid() && triggerHandle->size() > 0)
  {
    triggerBits = triggerHandle->at(0).TriggerBits();
  }

  if(e.isRealData())
  {
    // For data we can see if this event comes from a beam trigger
    triggerIsBeam = fDataUtil.IsBeamTrigger(e);
    if(triggerIsBeam){
      std::cout << "This data event has a beam trigger" << std::endl;
    }
  }

  tree->Fill();

} // analyze function

void lana::LSUBeamAnalyzer::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  // Book tree
  tree = tfs->make<TTree>("tree","tree");

  tree->Branch("isMC",&isMC,"isMC/O");
  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("subRunNumber",&subRunNumber,"subRunNumber/i");
  tree->Branch("eventNumber",&eventNumber,"eventNumber/i");
  tree->Branch("infilename",&infilename);

  tree->Branch("xWC",&xWC,"xWC/F");
  tree->Branch("yWC",&yWC,"yWC/F");
  tree->Branch("thetaWC",&thetaWC,"thetaWC/F");
  tree->Branch("phiWC",&phiWC,"phiWC/F");
  //tree->Branch("pzWC",&pzWC,"pzWC/F");
  tree->Branch("pWC",&pWC,"pWC/F");
  tree->Branch("eWC",&eWC,"eWC/F");
  tree->Branch("kinWC",&kinWC,"kinWC/F");
  //tree->Branch("kinWCInTPC",&kinWCInTPC,"kinWCInTPC/F");
  //tree->Branch("eWCProton",&eWCProton,"eWCProton/F");
  //tree->Branch("kinWCProton",&kinWCProton,"kinWCProton/F");
  //tree->Branch("yKinkWC",&yKinkWC,"yKinkWC/F");
  //tree->Branch("nHitsWC",&nHitsWC,"nHitsWC/i");
  //tree->Branch("xWC4Hit",&xWC4Hit,"xWC4Hit/F");
  //tree->Branch("yWC4Hit",&yWC4Hit,"yWC4Hit/F");
  //tree->Branch("zWC4Hit",&zWC4Hit,"zWC4Hit/F");

  tree->Branch("nBeamEvents",&nBeamEvents,"nBeamEvents/i");
  tree->Branch("BITrigger",&BITrigger,"BITrigger/I");
  tree->Branch("BITriggerMatched",&BITriggerMatched,"BITriggerMatched/O");
  //tree->Branch("BIActiveTrigger",&BIActiveTrigger,"BIActiveTrigger/i");
  tree->Branch("TOF",&TOF,"TOF/F");
  tree->Branch("TOFChan",&TOFChan,"TOFChan/I");
  tree->Branch("nTOFs",&nTOFs,"nTOFs/i");
  tree->Branch("TOFs",&TOFs,"TOFs[nTOFs]/F");
  tree->Branch("TOFChans",&TOFChans,"TOFChans[nTOFs]/I");
  tree->Branch("TOFusTrigs",&TOFusTrigs,"TOFusTrigs[nTOFs]/i");
  tree->Branch("TOFdsTrigs",&TOFdsTrigs,"TOFdsTrigs[nTOFs]/i");

  tree->Branch("CKov0Status",&CKov0Status,"CKov0Status/I");
  tree->Branch("CKov1Status",&CKov1Status,"CKov1Status/I");
  tree->Branch("CKov0Time",&CKov0Time,"CKov0Time/F");
  tree->Branch("CKov1Time",&CKov1Time,"CKov1Time/F");
  tree->Branch("CKov0Pressure",&CKov0Pressure,"CKov0Pressure/F");
  tree->Branch("CKov1Pressure",&CKov1Pressure,"CKov1Pressure/F");

  tree->Branch("nBeamMom",&nBeamMom,"nBeamMom/i");
  tree->Branch("beamMom",&beamMom,"beamMom[nBeamMom]/F");

  tree->Branch("nBeamTracks",&nBeamTracks,"nBeamTracks/i");
  tree->Branch("beamTrackXFrontTPC",&beamTrackXFrontTPC,"beamTrackXFrontTPC[nBeamTracks]/F");
  tree->Branch("beamTrackYFrontTPC",&beamTrackYFrontTPC,"beamTrackYFrontTPC[nBeamTracks]/F");
  tree->Branch("beamTrackTheta",&beamTrackTheta,"beamTrackTheta[nBeamTracks]/F");
  tree->Branch("beamTrackPhi",&beamTrackPhi,"beamTrackPhi[nBeamTracks]/F");
  tree->Branch("beamTrackMom",&beamTrackMom,"beamTrackMom[nBeamTracks]/F");
  tree->Branch("beamTrackEPion",&beamTrackEPion,"beamTrackEPion[nBeamTracks]/F");
  tree->Branch("beamTrackEProton",&beamTrackEProton,"beamTrackEProton[nBeamTracks]/F");
  tree->Branch("beamTrackKinPion",&beamTrackKinPion,"beamTrackKinPion[nBeamTracks]/F");
  tree->Branch("beamTrackKinProton",&beamTrackKinProton,"beamTrackKinProton[nBeamTracks]/F");
  tree->Branch("beamTrackTruePDG",&beamTrackTruePDG,"beamTrackTruePDG[nBeamTracks]/F");


  tree->Branch("triggerIsBeam",&triggerIsBeam,"triggerIsBeam/O");
  tree->Branch("triggerBits",&triggerBits,"triggerBits/i");

}

void lana::LSUBeamAnalyzer::beginRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void lana::LSUBeamAnalyzer::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void lana::LSUBeamAnalyzer::endJob()
{
  // Implementation of optional member function here.
}

void lana::LSUBeamAnalyzer::endRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void lana::LSUBeamAnalyzer::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void lana::LSUBeamAnalyzer::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fBeamEventTag = p.get<art::InputTag>("BeamEventTag");
  fRawTriggerTag = p.get<art::InputTag>("RawTriggerTag");
}

void lana::LSUBeamAnalyzer::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
  infilename = "INVALID";
}

void lana::LSUBeamAnalyzer::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lana::LSUBeamAnalyzer::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
  infilename = fb.fileName();
}

void lana::LSUBeamAnalyzer::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lana::LSUBeamAnalyzer::ResetTreeVars() 
{
  isMC = false;
  runNumber = 0;
  subRunNumber = 0;
  eventNumber = 0;

  xWC = DEFAULTNEG;
  yWC = DEFAULTNEG;
  thetaWC = DEFAULTNEG;
  phiWC = DEFAULTNEG;
  pzWC = DEFAULTNEG;
  pWC = DEFAULTNEG;
  eWC = DEFAULTNEG;
  kinWC = DEFAULTNEG;
  kinWCInTPC = DEFAULTNEG;
  eWCProton = DEFAULTNEG;
  kinWCProton = DEFAULTNEG;
  kinWCInTPCProton = DEFAULTNEG;
  yKinkWC = DEFAULTNEG;
  nHitsWC = 0;
  xWC4Hit = DEFAULTNEG;
  yWC4Hit = DEFAULTNEG;
  zWC4Hit = DEFAULTNEG;

  nBeamTracks = 0;
  for(size_t iTrack=0; iTrack < MAXBEAMTRACKS; iTrack++)
  {
    beamTrackXFrontTPC[iTrack] = DEFAULTNEG;
    beamTrackYFrontTPC[iTrack] = DEFAULTNEG;
    beamTrackTheta[iTrack] = DEFAULTNEG;
    beamTrackPhi[iTrack] = DEFAULTNEG;
    beamTrackMom[iTrack] = DEFAULTNEG;
    beamTrackEPion[iTrack] = DEFAULTNEG;
    beamTrackEProton[iTrack] = DEFAULTNEG;
    beamTrackKinPion[iTrack] = DEFAULTNEG;
    beamTrackKinProton[iTrack] = DEFAULTNEG;
    beamTrackTruePDG[iTrack] = DEFAULTNEG;
  }

  nBeamMom = 0;
  for(size_t iTrack=0; iTrack < MAXBEAMMOMS; iTrack++)
  {
    beamMom[iTrack] = DEFAULTNEG;
  }
 
  nBeamEvents = 0;
  BITrigger = DEFAULTNEG;
  BITriggerMatched = false;
  BIActiveTrigger = 0;
  TOF = DEFAULTNEG;
  TOFChan = DEFAULTNEG;
  nTOFs = 0;
  for(size_t iTOF; iTOF < MAXTOFS; iTOF++)
  {
    TOFs[iTOF] = DEFAULTNEG;
    TOFChans[iTOF] = DEFAULTNEG;
    TOFusTrigs[iTOF] = 0;
    TOFdsTrigs[iTOF] = 0;
  }

  triggerIsBeam = false;
  triggerBits = 0;

  CKov0Status = DEFAULTNEG;
  CKov1Status = DEFAULTNEG;
  CKov0Time = DEFAULTNEG;
  CKov1Time = DEFAULTNEG;
  CKov0Pressure = DEFAULTNEG;
  CKov1Pressure = DEFAULTNEG;

} // ResetTreeVars

DEFINE_ART_MODULE(lana::LSUBeamAnalyzer)
