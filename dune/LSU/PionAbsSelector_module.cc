////////////////////////////////////////////////////////////////////////
// Class:       PionAbsSelector
// Module Type: analyzer
// File:        PionAbsSelector_module.cc
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
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Cluster.h"
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
#include "dunetpc/dune/LSU/GenVectorHelper.h"
#include "dunetpc/dune/LSU/BackTrackerAlg.h"
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEDataUtils.h"
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEBeamlineUtils.h"
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "dunetpc/dune/Protodune/Analysis/ProtoDUNETruthUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "dunetpc/dune/DuneObj/ProtoDUNEBeamEvent.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

//ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/Point3D.h" // for 3D genvector points
#include "Math/Vector3D.h" // for 3D genvectors
#include "Math/Vector4Dfwd.h" // for 4D genvectors
//#include "Math/GenVector/DisplacementVector3D.h"
//#include "Math/GenVector/Polar3D.h"
#include "Math/GenVector/VectorUtil.h"

//c++ includes
#include <string>
#include <iomanip>
#include <algorithm> //for std::sort
#include <cmath>

const auto DEFAULTNEG = -99999999;
const auto DEFAULTPOS = 99999999;
const auto MAXTOFS = 10;
const auto MAXTRACKS = 1000;
const auto MAXDAUGHTER = 25;
const auto MAXMCPARTS = 10000;
const auto MAXTRUETRAJPOINTS = 1000;
const auto MAXBEAMTRACKS = 300;
const auto MAXBEAMMOMS = 300;
const auto MAXIDES = 20000;
const auto MAXPFSECTRKS = 25;
//const auto MAXPFSECSHWRS = 25;
const auto MAXZINT = 95;
const auto MAXLINT = 20;
const auto NTOFCHANS = 4;

const auto MMUON = 105.6583745; // MeV/c^2
const auto MCHARGEDPION = 139.57018; // MeV/c^2
const auto MPROTON = 938.2720813; // MeV/c^2
const auto KINLOSTBEFORETPC = 0.0; //MeV; from LArIAT pion total cross section group
const auto KINLOSTBEFORETPCPROTON = 0.0; //MeV; from LArIAT pion total cross section group

const auto ZSTARTOFTPC = -0.49375; //cm

namespace lana {
  class PionAbsSelector;
  struct PiAbsSecondary;
  //struct TreeSecondary;
}

struct lana::PiAbsSecondary //add more data members here as needed
{
  //data members
  const art::Ptr<recob::Track> fTrackPtr;
  const std::string fPidName;

  //constructor
  PiAbsSecondary(const art::Ptr<recob::Track>& track, const std::string pid): fTrackPtr(track), fPidName(pid)
  {}
};

template <class T, class U, class V>
inline void setHistTitles(T hist, U xtitle, V ytitle)
{
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetYaxis()->SetTitle(ytitle);
}

template <class T>
float findAverage(const T first, const T last)
{
  T tmpItr = first;
  size_t N = 0;
  double sum = 0.;
  for (;tmpItr<last;tmpItr++)
  {
    sum += *tmpItr;
    N++;
  }
  if (N==0) 
  {
    return DEFAULTNEG;
  }
  return sum / N;
}

template <class T>
float findQuantile(const float quantile, const T first, const T last)
{
  std::vector<double> tempVector(first,last);
  const size_t N = tempVector.size();
  if (N == 0)
  {
    return DEFAULTNEG;
  }
  if (N == 1)
  {
    return tempVector.at(0);
  }
  std::sort(tempVector.begin(),tempVector.end());
  const double iFloat = N*quantile;
  const size_t iFloor = std::floor(iFloat);
  const double iRemainder = iFloat - iFloor;
  if (iFloat == 0.)
  {
    return tempVector.at(0);
  }
  if (iFloat >= N)
  {
    return tempVector.at(N-1);
  }
  if (iRemainder == 0.) // exactly inbetween values
  {
    return 0.5*(tempVector.at(iFloor) + tempVector.at(iFloor-1));
  }
  return tempVector.at(iFloor);
}

template <class T, class C>
float findAverageResRangeFunc(const T resRangesBegin, const T resRangesEnd, const T first, const T last, C func)
{
  std::vector<double> tmpVec;
  T resRangeItr = resRangesBegin;
  T tmpItr = first;
  for(;resRangeItr < resRangesEnd && tmpItr < last; resRangeItr++, tmpItr++)
  {
    if(func(*resRangeItr))
        tmpVec.push_back(*tmpItr);
  }
  return findAverage(tmpVec.begin(),tmpVec.end());
}

template <class T, class C>
float findQuantileResRangeFunc(const float quantile, const T resRangesBegin, const T resRangesEnd, const T first, const T last, C func)
{
  std::vector<double> tmpVec;
  T resRangeItr = resRangesBegin;
  T tmpItr = first;
  for(;resRangeItr < resRangesEnd && tmpItr < last; resRangeItr++, tmpItr++)
  {
    if(func(*resRangeItr))
        tmpVec.push_back(*tmpItr);
  }
  return findQuantile(quantile,tmpVec.begin(),tmpVec.end());
}

class lana::PionAbsSelector : public art::EDAnalyzer {
public:
  explicit PionAbsSelector(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PionAbsSelector(PionAbsSelector const &) = delete;
  PionAbsSelector(PionAbsSelector &&) = delete;
  PionAbsSelector & operator = (PionAbsSelector const &) = delete;
  PionAbsSelector & operator = (PionAbsSelector &&) = delete;

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
  bool fDoPMAAnalysis;
  bool fDoPandoraAnalysis;
  art::InputTag fTruePartLabel; //The name of the module that produced simb::MCParticle objects
  art::InputTag fBeamTruthTag; //The name of the module that produced beam simb::MCTruth objects
  art::InputTag fCosmicTruthTag; //The name of the module that produced cosmic simb::MCTruth objects
  art::InputTag fSimChanLabel; //The name of the module that produced sim::SimChannel objects
  art::InputTag fTrackLabel; //The name of the module that produced recob::Track objects
  art::InputTag fCaloLabel; //The name of the module that produced anab::Calorimetry objects
  art::InputTag fBeamEventTag;
  art::InputTag fBeamEventTagOld;
  art::InputTag fRawTriggerTag;
  //art::InputTag fVertexLabel; //The name of the module that produced recob::Vertex objects
  art::InputTag fLikelihoodPIDTag;
  art::InputTag fPIDATag;
  art::InputTag fTrackHitTag;
  art::InputTag fAllHitTag;
  art::InputTag fPFParticleTag;
  art::InputTag fPFTrackTag;
  art::InputTag fPFShowerTag;
  art::InputTag fPFCaloTag;
  art::InputTag fPandoraHitTag;

  unsigned int fCaloPlane; //The plane we are using for calorimetry

  double fFiducialPrimaryXMin;
  double fFiducialPrimaryXMax;
  double fFiducialPrimaryYMin;
  double fFiducialPrimaryYMax;
  double fFiducialPrimaryZMin;
  double fFiducialPrimaryZMax;

  double fTrackMatchZMax;
  double fTrackMatchDeltaXMinData;
  double fTrackMatchDeltaXMaxData;
  double fTrackMatchDeltaXMinMC;
  double fTrackMatchDeltaXMaxMC;
  double fTrackMatchDeltaYMinData;
  double fTrackMatchDeltaYMaxData;
  double fTrackMatchDeltaYMinMC;
  double fTrackMatchDeltaYMaxMC;
  double fTrackMatchAngleMaxDeg;
  double fFlangeCenterX;
  double fFlangeCenterY;
  double fFlangeCenterZ;
  double fFlangeRadiusCut;

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
  Float_t xWC1Hit; // Fiber Tracker hit position x coord in cm
  Float_t yWC1Hit; // Fiber Tracker hit position y coord in cm
  Float_t zWC1Hit; // Fiber Tracker hit position z coord in cm
  Float_t xWC2Hit; // Fiber Tracker hit position x coord in cm
  Float_t yWC2Hit; // Fiber Tracker hit position y coord in cm
  Float_t zWC2Hit; // Fiber Tracker hit position z coord in cm

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
  float TOFsByChan[NTOFCHANS];
  UInt_t TOFusTrigsByChan[NTOFCHANS];
  UInt_t TOFdsTrigsByChan[NTOFCHANS];

  Int_t CKov0Status;
  Int_t CKov1Status;
  Float_t CKov0Time;
  Float_t CKov1Time;
  Float_t CKov0Pressure;
  Float_t CKov1Pressure;

  bool BIElectron0p3GeV;
  bool BIMuon0p3GeV;
  bool BIPion0p3GeV;
  bool BIKaon0p3GeV;
  bool BIProton0p3GeV;
  bool BIDeuteron0p3GeV;
  bool BIElectron0p5GeV;
  bool BIMuon0p5GeV;
  bool BIPion0p5GeV;
  bool BIKaon0p5GeV;
  bool BIProton0p5GeV;
  bool BIDeuteron0p5GeV;
  bool BIElectron1GeV;
  bool BIMuon1GeV;
  bool BIPion1GeV;
  bool BIKaon1GeV;
  bool BIProton1GeV;
  bool BIDeuteron1GeV;
  bool BIElectron2GeV;
  bool BIMuon2GeV;
  bool BIPion2GeV;
  bool BIKaon2GeV;
  bool BIProton2GeV;
  bool BIDeuteron2GeV;
  bool BIElectron3GeV;
  bool BIMuon3GeV;
  bool BIPion3GeV;
  bool BIKaon3GeV;
  bool BIProton3GeV;
  bool BIDeuteron3GeV;
  bool BIElectron6GeV;
  bool BIMuon6GeV;
  bool BIPion6GeV;
  bool BIKaon6GeV;
  bool BIProton6GeV;
  bool BIDeuteron6GeV;
  bool BIElectron7GeV;
  bool BIMuon7GeV;
  bool BIPion7GeV;
  bool BIKaon7GeV;
  bool BIProton7GeV;
  bool BIDeuteron7GeV;

  UInt_t nBeamTracksOld;
  UInt_t nBeamMomOld;
  Float_t beamMomOld[MAXBEAMMOMS];
  Float_t beamTrackXFrontTPCOld[MAXBEAMTRACKS];
  Float_t beamTrackYFrontTPCOld[MAXBEAMTRACKS];
  Float_t beamTrackThetaOld[MAXBEAMTRACKS];
  Float_t beamTrackPhiOld[MAXBEAMTRACKS];

  UInt_t nBeamEventsOld;
  Int_t BITriggerOld;
  bool BITriggerMatchedOld;
  UInt_t BIActiveTriggerOld;
  Float_t TOFOld;
  Int_t CKov0StatusOld;
  Int_t CKov1StatusOld;

  bool BIElectron0p3GeVOld;
  bool BIMuon0p3GeVOld;
  bool BIPion0p3GeVOld;
  bool BIKaon0p3GeVOld;
  bool BIProton0p3GeVOld;
  bool BIDeuteron0p3GeVOld;
  bool BIElectron0p5GeVOld;
  bool BIMuon0p5GeVOld;
  bool BIPion0p5GeVOld;
  bool BIKaon0p5GeVOld;
  bool BIProton0p5GeVOld;
  bool BIDeuteron0p5GeVOld;
  bool BIElectron1GeVOld;
  bool BIMuon1GeVOld;
  bool BIPion1GeVOld;
  bool BIKaon1GeVOld;
  bool BIProton1GeVOld;
  bool BIDeuteron1GeVOld;
  bool BIElectron2GeVOld;
  bool BIMuon2GeVOld;
  bool BIPion2GeVOld;
  bool BIKaon2GeVOld;
  bool BIProton2GeVOld;
  bool BIDeuteron2GeVOld;
  bool BIElectron3GeVOld;
  bool BIMuon3GeVOld;
  bool BIPion3GeVOld;
  bool BIKaon3GeVOld;
  bool BIProton3GeVOld;
  bool BIDeuteron3GeVOld;
  bool BIElectron6GeVOld;
  bool BIMuon6GeVOld;
  bool BIPion6GeVOld;
  bool BIKaon6GeVOld;
  bool BIProton6GeVOld;
  bool BIDeuteron6GeVOld;
  bool BIElectron7GeVOld;
  bool BIMuon7GeVOld;
  bool BIPion7GeVOld;
  bool BIKaon7GeVOld;
  bool BIProton7GeVOld;
  bool BIDeuteron7GeVOld;

  bool triggerIsBeam;
  UInt_t triggerBits;
  Int_t nGoodFEMBs[6];

  UInt_t nPrimaryParticleCandidates; // Number of MCParticles passing primary particle cuts

  Int_t trueCategory; // int showing which category this event is in
  Int_t trueEndProcess; // int showing which end process
  Int_t truePrimaryPDG; // int showing PDG code of primary
  Int_t truePrimaryTrackID; // int showing Track ID code of primary particle
  bool trueSignalT;
  bool trueSignalNT;
  UInt_t trueNDaughters; // number of daughter particles
  UInt_t nSecTracks; // number of daughters output info about 
  UInt_t trueNSecondaryChPions; // number of secondary pi+ and pi-
  UInt_t trueNSecondaryPiZeros; // number of secondary pi0
  UInt_t trueNSecondaryProtons; // number of secondary protons
  UInt_t trueNSecondaryOppChPions; // number of secondary pi+/- with opposite charge of primary (double chex)
  UInt_t trueNSecondaryMuons; // number of secondary muons
  Float_t trueStartX;  // primary MCParticle start position in cm
  Float_t trueStartY;
  Float_t trueStartZ;
  Float_t trueStartT; // ns
  Float_t trueEndX; // primary MCParticle end position in cm
  Float_t trueEndY;
  Float_t trueEndZ;
  Float_t trueEndT; // ns
  Float_t trueStartTheta; // primary MCParticle initial theta radians
  Float_t trueStartPhi; // primary MCParticle initial phi radians
  Float_t trueSecondToEndTheta; // primary MCParticle 2nd to last point theta radians (end is 0)
  Float_t trueSecondToEndPhi; // primary MCParticle 2nd to last point phi radians (end is 0)
  Float_t trueStartMom; // primary MCParticle initial momentum MeV/c
  Float_t trueStartE; // primary MCParticle initial energy MeV
  Float_t trueStartKin; // primary MCParticle initial kinetic energy MeV
  Float_t trueEndMom; // primary MCParticle end momentum MeV/c
  Float_t trueEndE; // primary MCParticle end energy MeV
  Float_t trueEndKin; // primary MCParticle end kinetic energy MeV
  Float_t trueEndMaxAngleToSec; // primary MCParticle end dir max angle to secondary start (only for +/-(pi/mu/p))
  Float_t trueSecondToEndMom;// primary MCParticle second to last point momentum MeV/c
  Float_t trueSecondToEndE;// primary MCParticle second to last point energy MeV
  Float_t trueSecondToEndKin;// primary MCParticle second to last point kinetic energy MeV
  Float_t trueAngleStartDirSecondToEndDir; // Angle between start dir and end dir in radians
  Int_t trueSecondPDG[MAXDAUGHTER]; // secondary Pdg Code
  Float_t trueSecondKin[MAXDAUGHTER]; // secondary Kinetic Energy
  Float_t trueSecondAngleToPrimEnd[MAXDAUGHTER]; // secondary angle to prim end
  Float_t trueXFrontTPC; // the starting trajectory projected to the TPC, x coord
  Float_t trueYFrontTPC; // the starting trajectory projected to the TPC, y coord
  Float_t trueTrajPointFrontTPCDist; // distance from front of TPC to closest traj point on true MCPart cm
  Float_t trueKinFrontTPC; // kinetic energy of primaryParticle at traj point closest to front of TPC MeV

  UInt_t trueTrajNPoints;
  Float_t trueTrajPointX[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointY[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointZ[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointT[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointPx[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointPy[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointPz[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointE[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointMom[MAXTRUETRAJPOINTS];
  Float_t trueTrajPointKin[MAXTRUETRAJPOINTS];
  UChar_t trueTrajPointProc[MAXTRUETRAJPOINTS]; // 0 unknown, 1 hadElastic 2-7 inelastic depend on hadron
  Float_t trueTrajPointKinkAngle[MAXTRUETRAJPOINTS]; // angle between this point mom and last point mom

  UInt_t nMCParts;
  bool mcPartIsBeam[MAXMCPARTS];
  bool mcPartIsPrimary[MAXMCPARTS];
  Int_t mcPartTrackID[MAXMCPARTS];
  Int_t mcPartPDG[MAXMCPARTS];
  Float_t mcPartStartX[MAXMCPARTS];
  Float_t mcPartStartY[MAXMCPARTS];
  Float_t mcPartStartZ[MAXMCPARTS];
  Float_t mcPartStartT[MAXMCPARTS];
  Float_t mcPartEndX[MAXMCPARTS];
  Float_t mcPartEndY[MAXMCPARTS];
  Float_t mcPartEndZ[MAXMCPARTS];
  Float_t mcPartEndT[MAXMCPARTS];
  Float_t mcPartStartTheta[MAXMCPARTS];
  Float_t mcPartStartPhi[MAXMCPARTS];
  Float_t mcPartXFrontTPC[MAXMCPARTS];
  Float_t mcPartYFrontTPC[MAXMCPARTS];
  Float_t mcPartStartMom[MAXMCPARTS];
  Float_t mcPartStartE[MAXMCPARTS];
  Float_t mcPartStartKin[MAXMCPARTS];
  Float_t mcPartEndMom[MAXMCPARTS];
  Float_t mcPartEndE[MAXMCPARTS];
  Float_t mcPartEndKin[MAXMCPARTS];
  Float_t mcPartDeltaAngle[MAXMCPARTS];

  UInt_t nIDEs;
  Float_t simIDENumElectrons[MAXIDES]; // number of electrons at wire for this time tick
  Float_t simIDEEnergy[MAXIDES]; // number energy deposited for this time tick in MeV
  Float_t simIDEX[MAXIDES]; // true X position cm
  Float_t simIDEY[MAXIDES]; // true Y position cm
  Float_t simIDEZ[MAXIDES]; // true Z position cm
  Float_t simIDETrajDistance[MAXIDES]; // simIDe distance from true trajectory cm
  Float_t simIDETrajKin[MAXIDES]; // simIDE KE interpolated from trajectory MeV
  Float_t simIDEWireZ[MAXIDES]; // wire z position cm
  UInt_t simIDEChannel[MAXIDES]; // raw::ChannelID_t for this IDE
  Float_t simIDETDC[MAXIDES]; // TDC tick number for this IDE
  Float_t simIDEPartKin[MAXIDES]; // primary particle KE at this point MeV
  Float_t simIDEEnergySum; // sum of all simIDEEnergy

  std::vector<UInt_t> zWireChannel;
  std::vector<float> zWireWireZ;

  std::vector<UInt_t> zWireNHits;
  std::vector<float> zWireResRange;
  std::vector<float> zWiredEdx;
  std::vector<float> zWiredQdx;
  std::vector<float> zWirePitch;
  std::vector<float> zWireX;
  std::vector<float> zWireY;
  std::vector<float> zWireZ;
  std::vector<float> zWirePartKin;
  std::vector<float> zWirePartKinProton;
  Float_t zWireEnergySum;
  
  std::vector<float> zWireTrueNumElectrons;
  std::vector<float> zWireTrueEnergy;
  std::vector<float> zWireTrueX;
  std::vector<float> zWireTrueY;
  std::vector<float> zWireTrueZ;
  std::vector<float> zWireTruedZ;
  std::vector<float> zWireTruedR;
  std::vector<float> zWireTrueTDC;
  std::vector<float> zWireTrueTrajDistance;
  std::vector<float> zWireTrueTrajPitch;
  std::vector<float> zWireTrueTrajKin;
  std::vector<float> zWireTruePartKin;
  Float_t zWireTrueEnergySum;

  UInt_t nTracks;
  UInt_t nTracksInFirstZ[MAXZINT]; // the number of tracks with a space point in (0,i) cm where i is the index
  UInt_t nTracksLengthLt[MAXLINT]; // the number of tracks with length less than i cm where i is the index
  Float_t trackStartX[MAXTRACKS];
  Float_t trackStartY[MAXTRACKS];
  Float_t trackStartZ[MAXTRACKS];
  Float_t trackStartTheta[MAXTRACKS];
  Float_t trackStartPhi[MAXTRACKS];
  Float_t trackEndX[MAXTRACKS];
  Float_t trackEndY[MAXTRACKS];
  Float_t trackEndZ[MAXTRACKS];
  Float_t trackLength[MAXTRACKS];
  Float_t trackXFrontTPC[MAXTRACKS];
  Float_t trackYFrontTPC[MAXTRACKS];
  Float_t trackCaloKin[MAXTRACKS];
  Float_t trackLLHPion[MAXTRACKS];
  Float_t trackLLHProton[MAXTRACKS];
  Float_t trackLLHMuon[MAXTRACKS];
  Float_t trackLLHKaon[MAXTRACKS];
  Float_t trackPIDA[MAXTRACKS];
  Float_t trackStartDistToPrimTrkEnd[MAXTRACKS];
  Float_t trackEndDistToPrimTrkEnd[MAXTRACKS];
  Float_t trackClosestDistToPrimTrkEnd[MAXTRACKS];
  bool trackStartClosestToPrimTrkEnd[MAXTRACKS];
  Int_t trackTrueID[MAXTRACKS];
  Int_t trackTrueMotherID[MAXTRACKS];
  Int_t trackTruePdg[MAXTRACKS];
  bool trackTrueIsBeam[MAXTRACKS];
  Float_t trackTrueKin[MAXTRACKS];
  Float_t trackTrueEndKin[MAXTRACKS];
  Float_t trackTrueTrajLen[MAXTRACKS];
  Float_t trackTrueChargePurity[MAXTRACKS];
  Float_t trackTrueChargeEfficiencyU[MAXTRACKS];
  Float_t trackTrueChargeEfficiencyV[MAXTRACKS];
  Float_t trackTrueChargeEfficiencyZ[MAXTRACKS];

  Float_t trackTrueStartX[MAXTRACKS];
  Float_t trackTrueStartY[MAXTRACKS];
  Float_t trackTrueStartZ[MAXTRACKS];
  Float_t trackTrueStartT[MAXTRACKS];
  Float_t trackTrueEndX[MAXTRACKS];
  Float_t trackTrueEndY[MAXTRACKS];
  Float_t trackTrueEndZ[MAXTRACKS];
  Float_t trackTrueEndT[MAXTRACKS];
  Float_t trackTrueXFrontTPC[MAXTRACKS];
  Float_t trackTrueYFrontTPC[MAXTRACKS];

  Int_t iBestMatch;
  Float_t trackMatchDeltaX[MAXTRACKS];
  Float_t trackMatchDeltaY[MAXTRACKS];
  Float_t trackMatchDeltaR[MAXTRACKS];
  Float_t trackMatchDeltaAngle[MAXTRACKS];
  Float_t trackMatchLowestZ[MAXTRACKS];
  UInt_t nMatchedTracks; // number of tracks in this event that pass matching criteria

  bool primTrkIsMatchPrimary; // primary TPC track is matched to the primary MCParticle
  bool primTrkIsMatchPrimaryDaughter; // primary TPC track is matched to a daughter of the primary MCPart
  bool primTrkIsMatchBeam; // primary TPC track is matched to a beam MCParticle
  bool primTrkIsMatchAPrimary; // primary TPC track is matched to a primary (like a halo muon)

  Float_t primTrkStartMomTrking; // primary TPC track initial momentum from tracking
  Float_t primTrkStartTheta;  // primary TPC track initial theta
  Float_t primTrkStartPhi; // primary TPC track initial phi
  Float_t primTrkLength; // primary TPC track length cm
  Float_t primTrkStartX; // primary TPC track start x cm
  Float_t primTrkStartY; // primary TPC track start y cm
  Float_t primTrkStartZ; // primary TPC track start z cm
  Float_t primTrkEndX; // primary TPC track End x cm
  Float_t primTrkEndY; // primary TPC track End y cm
  Float_t primTrkEndZ; // primary TPC track End z cm
  Float_t primTrkXFrontTPC;
  Float_t primTrkYFrontTPC;
  Float_t primTrkCaloRange; // Total track length from calo
  bool primTrkEndInFid; // ends in primary fiducial

  Float_t primTrkCaloKin; // primary TPC track calo KE MeV
  Float_t primTrkEndKin; // kinWCInTPC - primTrkCaloKin MeV
  Float_t primTrkEndKinFid; // same as primTrkEndKin, except if end out of fid = DEFAULTNEG
  Float_t primTrkKinInteract; // The kinetic energy at the last hit on the track. If out of fid, = DEFAULTNEG
  Float_t primTrkKinInteractProton; // The kinetic energy at the last hit on the track. If out of fid, = DEFAULTNEG
  Float_t primTrkNHits; // The number of hits in the calo object in the collection plane

  bool primTrkResRangesFlipped; // true if resRange larger at back than front
  Float_t primTrkdEdxMedianLast3Hits; // median dE/dx of last 3 hits
  Float_t primTrkdEdxAverageLast3Hits; // avg dE/dx of last 3 hits
  Float_t primTrkdEdxMedianLast5Hits;
  Float_t primTrkdEdxAverageLast5Hits;
  Float_t primTrkdEdxMedianLast7Hits;
  Float_t primTrkdEdxAverageLast7Hits;
  Float_t primTrkdEdxMedianRRL1; // median dE/dx of hits with residual range < 1 cm
  Float_t primTrkdEdxAverageRRL1;
  Float_t primTrkdEdxMedianRRL2;
  Float_t primTrkdEdxAverageRRL2;
  Float_t primTrkdEdxMedianRRL3;
  Float_t primTrkdEdxAverageRRL3;
  Float_t primTrkdEdxMedianRRL5;
  Float_t primTrkdEdxAverageRRL5;
  Float_t primTrkdEdxMedianRRL7;
  Float_t primTrkdEdxAverageRRL7;
  Float_t primTrkdEdxMedianRRL3G1; // median dE/dx of hits with 1cm < residual range < 3 cm
  Float_t primTrkdEdxAverageRRL3G1;
  Float_t primTrkdEdxMedianRRL5G1;
  Float_t primTrkdEdxAverageRRL5G1;
  Float_t primTrkdEdxMedianRRL7G1;
  Float_t primTrkdEdxAverageRRL7G1;

  Float_t primTrkLLHPion; // log likelihood is pion
  Float_t primTrkLLHProton;
  Float_t primTrkLLHMuon;
  Float_t primTrkLLHKaon;
  Float_t primTrkPIDA;

  std::vector<float> primTrkdEdxs; // dE/dx values for primary track MeV/cm
  std::vector<float> primTrkResRanges; // Residual range values for primary track cm
  std::vector<float> primTrkRangeSoFars; // Range in TPC so far for primary track cm
  std::vector<float> primTrkPitches; // pitch values for every prim track hit cm
  std::vector<float> primTrkIBackwards; // Counting down from size -1 to 0
  std::vector<float> primTrkXs; // hit x positions in cm
  std::vector<float> primTrkYs; // hit y positions in cm
  std::vector<float> primTrkZs; // hit z positions in cm
  std::vector<float> primTrkKins; // Calc Kin from kinWCInTPC and integrating dE/dx
  std::vector<float> primTrkKinsProton; // Calc Kin from kinWCInTPCProton and integrating dE/dx
  std::vector<bool> primTrkInFids; // Is the calo point in the primary fiducial region
  std::vector<float> primTrkKinsTrue; // interpolated true Kin at calo location
  std::vector<float> primTrkDistToTrueTraj; // distance from calo location to true trajectory
  std::vector<float> primTrkDistToTrueTrajPoint; // distance from calo location to true trajectory point

  Int_t nSecTrk; // number of secondary tracks with start dist to prim end < 2.5 cm
  Int_t iSecTrkID[MAXTRACKS];  //Track numbers of secondaries
  bool SecTrkPID[MAXTRACKS];   //Secondary okay for min LLR calc
  Int_t iSecMin;  // secondary with min LLR
  Float_t nSecLLRProtonToPionMax;
  Float_t nSecLLRProtonToPionMin;

  std::vector<float> secMinTrkdEdxs; // dE/dx values for iSecMin track MeV/cm
  std::vector<float> secMinTrkResRanges; // Residual range values for iSecMin track cm
  std::vector<float> secMinTrkPitches; // pitch values for every secMin track hit cm
  std::vector<float> secMinTrkIBackwards; // Counting down from size -1 to 0
  std::vector<float> secMinTrkXs; // hit x positions in cm
  std::vector<float> secMinTrkYs; // hit y positions in cm
  std::vector<float> secMinTrkZs; // hit z positions in cm
  std::vector<bool> secMinTrkInFids; // Is the calo point in the iSecMin fiducial region

  Int_t nSecTrkLLRG0; // number of secondary tracks with pi-p LLR > 0
  Int_t nSecTrkLLRG100;
  Int_t nSecTrkLLRG200;
  Int_t nSecTrkLLRG300;
  Int_t nSecTrkLLRG400;
  Int_t nSecTrkLLRG500;
  Int_t nSecTrkLLRG600;
  Int_t nSecTrkLLRG700;
  Int_t nSecTrkLLRG800;
  Int_t nSecTrkLLRG900;
  Int_t nSecTrkPIDAL8; // number of secondary tracks with PIDA > 8
  Int_t nSecTrkPIDAL10;
  Int_t nSecTrkPIDAL14;
  Int_t nSecTrkPIDAL16;
  Int_t nSecTrkPIDAL18;

  // PFParticle/Pandora stuff
  UInt_t PFNBeamSlices;
  UInt_t PFBeamPrimNDaughters;
  UInt_t PFBeamPrimNDaughterTracks;
  UInt_t PFBeamPrimNDaughterShowers;
  Int_t PFBeamPrimPDG;
  bool PFBeamPrimIsTracklike;
  bool PFBeamPrimIsShowerlike;
  float PFBeamPrimBeamCosmicScore;
  float PFBeamPrimXFrontTPC;
  float PFBeamPrimYFrontTPC;
  float PFBeamPrimXFrontTPCTrackEnd;
  float PFBeamPrimYFrontTPCTrackEnd;
  float PFBeamPrimStartX;
  float PFBeamPrimStartY;
  float PFBeamPrimStartZ;
  float PFBeamPrimEndX;
  float PFBeamPrimEndY;
  float PFBeamPrimEndZ;
  float PFBeamPrimStartTheta;
  float PFBeamPrimStartPhi;
  float PFBeamPrimEndTheta;
  float PFBeamPrimEndPhi;
  float PFBeamPrimAngleToBeamTrk;
  float PFBeamPrimEndAngleToBeamTrk;
  float PFBeamPrimTrkLen;
  float PFBeamPrimTrkMaxKink;
  float PFBeamPrimTrkStartEndDirAngle;
  float PFBeamPrimMaxAngleToSecTrks;
  float PFBeamPrimShwrLen;
  float PFBeamPrimShwrOpenAngle;
  float PFBeamPrimCSDAMomMu;
  float PFBeamPrimCSDAMomProton;
  float PFBeamPrimCSDAKinMu;
  float PFBeamPrimCSDAKinProton;
  float PFBeamPrimdEdxAverageLast3Hits;
  float PFBeamPrimdEdxAverageLast5Hits;
  float PFBeamPrimdEdxAverageLast7Hits;

  std::vector<float> PFBeamPrimResRanges;
  std::vector<float> PFBeamPrimdEdxs;
  std::vector<float> PFBeamPrimdQdxs;
  std::vector<float> PFBeamPrimPitches;
  std::vector<float> PFBeamPrimXs;
  std::vector<float> PFBeamPrimYs;
  std::vector<float> PFBeamPrimZs;
  std::vector<Int_t> PFBeamPrimZWires;
  std::vector<float> PFBeamPrimKins;
  std::vector<float> PFBeamPrimKinsProton;
  float PFBeamPrimKinInteract;
  float PFBeamPrimKinInteractProton;

  Int_t PFBeamPrimTrueTrackID;
  bool PFBeamPrimTrueTrackIDIsNeg;
  Int_t PFBeamPrimTrueMotherTrackID;
  Int_t PFBeamPrimTruePDG;
  float PFBeamPrimTrueKin;
  float PFBeamPrimTrueEndKin;
  UInt_t PFBeamPrimTrueTrajLen;
  bool PFBeamPrimTrueIsBeam;
  float PFBeamPrimTrueStartX;
  float PFBeamPrimTrueStartY;
  float PFBeamPrimTrueStartZ;
  float PFBeamPrimTrueStartT;
  float PFBeamPrimTrueEndX;
  float PFBeamPrimTrueEndY;
  float PFBeamPrimTrueEndZ;
  float PFBeamPrimTrueEndT;
  float PFBeamPrimTrueMaxKink;

  float PFBeamSecTrkLen[MAXPFSECTRKS];
  float PFBeamSecTrkdEdxAverageLast3Hits[MAXPFSECTRKS];
  float PFBeamSecTrkdEdxAverageLast5Hits[MAXPFSECTRKS];
  float PFBeamSecTrkdEdxAverageLast7Hits[MAXPFSECTRKS];
  float PFBeamSecTrkDistToPrimEnd[MAXPFSECTRKS];
  float PFBeamSecTrkAngleToPrimEnd[MAXPFSECTRKS];

  /////////////////////////////
  //Histograms

  TH2F* deltaXYTPCBeamlineHist;
  TH2F* deltaXYTPCBeamlineOnlyInFlangeHist;
  TH2F* deltaXYTPCBeamlineOnlyInFirst25cmHist;
  TH2F* deltaXYTPCBeamlineOnlyInFlangeInFirst25cmHist;

  TH1F* deltaAngleTPCBeamlineHist;
  TH1F* deltaAngleTPCBeamlineOnlyInFlangeHist;
  TH1F* deltaAngleTPCBeamlineOnlyInFirst25cmHist;
  TH1F* deltaAngleTPCBeamlineOnlyInFlangeInFirst25cmHist;

  std::map<UInt_t,std::pair<UInt_t,float>> myChannelToWireMap;

  art::ServiceHandle<cheat::BackTrackerService> bt;

  protoana::ProtoDUNEDataUtils fDataUtil;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtilsOld;
 
  //Internal functions
  const art::Ptr<recob::Track> MatchRecoToTruthOrWCTrack(const std::vector<art::Ptr<recob::Track>>& tracks, bool isData); 
  const art::Ptr<simb::MCParticle> ProcessMCParticles(const std::vector<art::Ptr<simb::MCParticle>>& truePartVec, 
        const pdana::MCBeamOrCosmicAlg * const beamOrCosmic, 
        const art::Event& e); // returns primary particle candidate
  void ProcessAllTracks(const std::vector<art::Ptr<recob::Track>>& trackVec,
        const std::vector<art::Ptr<simb::MCParticle>>& truePartVec,
        const art::FindManyP<anab::Calorimetry>& tracksCaloVec, 
        const art::FindManyP<recob::Hit>& fmHitsForTracks,
        const pdana::MCBeamOrCosmicAlg * const beamOrCosmic,
        const art::Event& e);
  void ProcessPrimaryTrack(const std::vector<art::Ptr<recob::Track>>& trackVec,
        const art::Ptr<simb::MCParticle>& primaryParticle,
        const art::FindManyP<anab::Calorimetry>& tracksCaloVec, 
        const pdana::MotherDaughterWalkerAlg& motherDaughterWalker,
        const art::Event& e); // returns primary track candidate
  void ProcessPFParticles(const art::Event& e,
        //const art::Ptr<simb::MCParticle>& primaryParticle,
        //const pdana::MotherDaughterWalkerAlg& motherDaughterWalker,
        const pdana::MCBeamOrCosmicAlg * const beamOrCosmic
        );

  void ResetTreeVars();

  template <class T>
  bool InPrimaryFiducial(const T & pos); // use a type with X(), Y(), Z() methods

};

template <class T>
bool lana::PionAbsSelector::InPrimaryFiducial(const T & pos)
{
  bool isInX = pos.X() > fFiducialPrimaryXMin && pos.X() < fFiducialPrimaryXMax;
  bool isInY = pos.Y() > fFiducialPrimaryYMin && pos.Y() < fFiducialPrimaryYMax;
  bool isInZ = pos.Z() > fFiducialPrimaryZMin && pos.Z() < fFiducialPrimaryZMax;
  bool isInAll = isInX && isInY && isInZ;
  return isInAll;
}

lana::PionAbsSelector::PionAbsSelector(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  infilename("NoInFilenameFound"),
  fDataUtil(p.get<fhicl::ParameterSet>("DataUtil")),
  fBeamlineUtils(p.get<fhicl::ParameterSet>("BeamlineUtils")),
  fBeamlineUtilsOld(p.get<fhicl::ParameterSet>("BeamlineUtilsOld"))
 // More initializers here.
{
  this->reconfigure(p);
}

void lana::PionAbsSelector::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  ResetTreeVars();

  runNumber = e.run();
  subRunNumber = e.subRun();
  eventNumber = e.id().event();
  isMC = !(e.isRealData());

  //Get needed data products
  std::vector<art::Ptr<simb::MCParticle>> truePartVec;
  if(!e.isRealData())
  {
    auto truePartHand = e.getValidHandle<std::vector<simb::MCParticle>>(fTruePartLabel);
    if(truePartHand.isValid())
    {
      art::fill_ptr_vector(truePartVec, truePartHand);
    }
  }

  pdana::MCBeamOrCosmicAlg* beamOrCosmic = NULL;
  if(!e.isRealData())
  {
    beamOrCosmic = new pdana::MCBeamOrCosmicAlg(e,fTruePartLabel,fBeamTruthTag,fCosmicTruthTag);
  }
  pdana::MotherDaughterWalkerAlg motherDaughterWalker(e,fTruePartLabel);

  auto trackHand = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
  std::vector<art::Ptr<recob::Track>> trackVec;
  if(trackHand.isValid())
  {
    art::fill_ptr_vector(trackVec, trackHand);
  }

  // One element for each track. Each one of those is a vector of calos for each wire plane
  const auto tracksCaloVec = art::FindManyP<anab::Calorimetry>(trackVec, e, fCaloLabel);

  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if(e.isRealData())
  {
    auto beamHand = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
    const auto prov = beamHand.provenance();
    const auto& provTag = prov->inputTag();
    std::cout << "BeamEvent provenance inputTag: '" << provTag.encode() << "'\n";
    //prov->write(std::cout);
    //std::cout << "\nBeamEvent provenance done\n";
    if(beamHand.isValid())
    {
      art::fill_ptr_vector(beamVec, beamHand);
    }
  }

  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVecOld;
  if(e.isRealData())
  {
    auto beamHandOld = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTagOld);
    const auto prov = beamHandOld.provenance();
    const auto& provTag = prov->inputTag();
    std::cout << "BeamEvent Old provenance inputTag: '" << provTag.encode() << "'\n";
    //prov->write(std::cout);
    //std::cout << "\nBeamEvent Old provenance done\n";
    if(beamHandOld.isValid())
    {
      art::fill_ptr_vector(beamVecOld, beamHandOld);
    }
  }

  art::FindManyP<anab::ParticleID>  fmpida(trackHand, e, fPIDATag);
  art::FindManyP<anab::ParticleID>  fmlhpid(trackHand, e, fLikelihoodPIDTag);

  // For track MC truth matching
  art::FindManyP<recob::Hit> fmHitsForTracks(trackHand, e, fTrackHitTag);

  std::vector<art::Ptr<recob::Hit>> allHitsVec;
  if (!e.isRealData())
  {
    auto allHitsHand = e.getValidHandle<std::vector<recob::Hit>>(fAllHitTag);
    if(allHitsHand.isValid())
    {
      art::fill_ptr_vector(allHitsVec, allHitsHand);
    }
  }

  // Triggers
  art::Handle< std::vector<raw::Trigger> > triggerHandle;
  if(e.isRealData())
  {
    e.getByLabel(fRawTriggerTag,triggerHandle);
  }

  // Beamline info
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
      std::cout << "PiAbsSel BeamEvent " << iBeamEvent << " for "<< e.id() << ": \n";
      std::cout << "  CTB Timestamp: " << beamEvent.GetCTBTimestamp() << "\n";
      std::cout << "  BI Trigger: " << beamEvent.GetBITrigger() << "\n";
      std::cout << "  Active Trigger: " << beamEvent.GetActiveTrigger() << "\n";
      std::cout << "  Is Trigger Matched: " << beamEvent.CheckIsMatched() << "\n";
      std::cout << "  TOF: " << beamEvent.GetTOF() << "\n";
      std::cout << "  CKov0Status: " << beamEvent.GetCKov0Status() << "\n";
      std::cout << "  CKov1Status: " << beamEvent.GetCKov1Status() << "\n";
      //std::cout << "  CKov0Time: " << beamEvent.GetCKov0Time() << "\n"; // don't think these are filled
      //std::cout << "  CKov1Time: " << beamEvent.GetCKov1Time() << "\n"; // don't think these are filled
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
    // Backwards so that the first one is the one we record
    for(int jTOF=nTOFs-1; jTOF >= 0; jTOF--)
    {
      const size_t iTOF=jTOF;
      if(tofChans.size() > iTOF)
      {
        const auto & chan = tofChans.at(iTOF);
        if (chan >= 0 && chan < NTOFCHANS)
        {
          TOFsByChan[chan] = tofs.at(iTOF);
          if(usTOFTrigs.size() > iTOF)
          {
            TOFusTrigsByChan[chan] = usTOFTrigs.at(iTOF);
          }
          if(dsTOFTrigs.size() > iTOF)
          {
            TOFdsTrigsByChan[chan] = dsTOFTrigs.at(iTOF);
          }
        }
      }
    } // for jTOF backwards

    const bool sameNTracksAsMom = beamEvent.GetNBeamTracks() == beamEvent.GetNRecoBeamMomenta();
    for(size_t iMom=0; iMom < beamEvent.GetNRecoBeamMomenta(); iMom++)
    {
      beamMom[nBeamMom] = beamEvent.GetRecoBeamMomentum(iMom);
      nBeamMom++;
      if(PRINTBEAMEVENT) std::cout << "    " << beamEvent.GetRecoBeamMomentum(iMom) << "\n";
      if(iMom == 0)
      {
        pWC = beamEvent.GetRecoBeamMomentum(iMom)*1000.;
      }
    }
    for(size_t iTrack=0; iTrack < beamEvent.GetNBeamTracks(); iTrack++)
    {
      const recob::Track & track =  beamEvent.GetBeamTrack(iTrack);
      if(nBeamTracks >= MAXBEAMTRACKS)
      {
        mf::LogError("PionAbsSelector") << "Event " << e.id() << " not stored because too many beam tracks";
        return;
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

      if(iTrack == 0)
      {
        xWC = track.End().X();
        yWC = track.End().Y();
        xWC1Hit = track.LocationAtPoint(0).X();
        yWC1Hit = track.LocationAtPoint(0).Y();
        zWC1Hit = track.LocationAtPoint(0).Z();
        xWC2Hit = track.LocationAtPoint(1).X();
        yWC2Hit = track.LocationAtPoint(1).Y();
        zWC2Hit = track.LocationAtPoint(1).Z();
        thetaWC = track.EndDirection().Theta();
        phiWC = track.EndDirection().Phi();
        track.NPoints();
      }

      if(PRINTBEAMEVENT)
      {
        std::cout << "  Beam Track: "<< iTrack <<"\n";
        std::cout << "    N Points:  " << track.NumberTrajectoryPoints() << "\n";
        std::cout << "    Start Pos: " << track.Vertex().X()
                                  << "  " << track.Vertex().Y()
                                  << "  " << track.Vertex().Z() << "\n";
        std::cout << "    Mid Pos:   " << track.LocationAtPoint(1).X()
                                  << "  " << track.LocationAtPoint(1).Y()
                                  << "  " << track.LocationAtPoint(1).Z() << "\n";
        std::cout << "    End Pos:   " << track.End().X()
                                  << "  " << track.End().Y()
                                  << "  " << track.End().Z() << "\n";
        std::cout << "    Start Theta: " << track.VertexDirection().Theta()*180/CLHEP::pi << " deg\n";
        std::cout << "    Start Phi:   " << track.VertexDirection().Phi()*180/CLHEP::pi << " deg\n";
        std::cout << "    End Theta:   " << track.EndDirection().Theta()*180/CLHEP::pi << " deg\n";
        std::cout << "    End Phi:     " << track.EndDirection().Phi()*180/CLHEP::pi << " deg\n";
      }
    }

    // CERN official beam PID
    //const auto& beamPIDs0p3GeV = fBeamlineUtils.GetPIDCandidates(beamEvent,0.3);
    //const auto& beamPIDs0p5GeV = fBeamlineUtils.GetPIDCandidates(beamEvent,0.5);
    const auto& beamPIDs1GeV = fBeamlineUtils.GetPIDCandidates(beamEvent,1);
    const auto& beamPIDs2GeV = fBeamlineUtils.GetPIDCandidates(beamEvent,2);
    const auto& beamPIDs3GeV = fBeamlineUtils.GetPIDCandidates(beamEvent,3);
    const auto& beamPIDs6GeV = fBeamlineUtils.GetPIDCandidates(beamEvent,6);
    const auto& beamPIDs7GeV = fBeamlineUtils.GetPIDCandidates(beamEvent,7);

    //BIElectron0p3GeV = beamPIDs0p3GeV.electron;
    //BIMuon0p3GeV = beamPIDs0p3GeV.muon;
    //BIPion0p3GeV = beamPIDs0p3GeV.pion;
    //BIKaon0p3GeV = beamPIDs0p3GeV.kaon;
    //BIProton0p3GeV = beamPIDs0p3GeV.proton;
    //BIDeuteron0p3GeV = beamPIDs0p3GeV.deuteron;
    //BIElectron0p5GeV = beamPIDs0p5GeV.electron;
    //BIMuon0p5GeV = beamPIDs0p5GeV.muon;
    //BIPion0p5GeV = beamPIDs0p5GeV.pion;
    //BIKaon0p5GeV = beamPIDs0p5GeV.kaon;
    //BIProton0p5GeV = beamPIDs0p5GeV.proton;
    //BIDeuteron0p5GeV = beamPIDs0p5GeV.deuteron;
    BIElectron1GeV = beamPIDs1GeV.electron;
    BIMuon1GeV = beamPIDs1GeV.muon;
    BIPion1GeV = beamPIDs1GeV.pion;
    BIKaon1GeV = beamPIDs1GeV.kaon;
    BIProton1GeV = beamPIDs1GeV.proton;
    BIDeuteron1GeV = beamPIDs1GeV.deuteron;
    BIElectron2GeV = beamPIDs2GeV.electron;
    BIMuon2GeV = beamPIDs2GeV.muon;
    BIPion2GeV = beamPIDs2GeV.pion;
    BIKaon2GeV = beamPIDs2GeV.kaon;
    BIProton2GeV = beamPIDs2GeV.proton;
    BIDeuteron2GeV = beamPIDs2GeV.deuteron;
    BIElectron3GeV = beamPIDs3GeV.electron;
    BIMuon3GeV = beamPIDs3GeV.muon;
    BIPion3GeV = beamPIDs3GeV.pion;
    BIKaon3GeV = beamPIDs3GeV.kaon;
    BIProton3GeV = beamPIDs3GeV.proton;
    BIDeuteron3GeV = beamPIDs3GeV.deuteron;
    BIElectron6GeV = beamPIDs6GeV.electron;
    BIMuon6GeV = beamPIDs6GeV.muon;
    BIPion6GeV = beamPIDs6GeV.pion;
    BIKaon6GeV = beamPIDs6GeV.kaon;
    BIProton6GeV = beamPIDs6GeV.proton;
    BIDeuteron6GeV = beamPIDs6GeV.deuteron;
    BIElectron7GeV = beamPIDs7GeV.electron;
    BIMuon7GeV = beamPIDs7GeV.muon;
    BIPion7GeV = beamPIDs7GeV.pion;
    BIKaon7GeV = beamPIDs7GeV.kaon;
    BIProton7GeV = beamPIDs7GeV.proton;
    BIDeuteron7GeV = beamPIDs7GeV.deuteron;

    if(PRINTBEAMEVENT)
    {
      std::cout << "CERN BIG PID: e  mu pi k  p  d  \n"
                << "1 GeV:        " << beamPIDs1GeV.electron
                            << "  " << beamPIDs1GeV.muon
                            << "  " << beamPIDs1GeV.pion
                            << "  " << beamPIDs1GeV.kaon
                            << "  " << beamPIDs1GeV.proton
                            << "  " << beamPIDs1GeV.deuteron
                                    << "\n"
                << "2 GeV:        " << beamPIDs2GeV.electron
                            << "  " << beamPIDs2GeV.muon
                            << "  " << beamPIDs2GeV.pion
                            << "  " << beamPIDs2GeV.kaon
                            << "  " << beamPIDs2GeV.proton
                            << "  " << beamPIDs2GeV.deuteron
                                    << "\n"
                << "3 GeV:        " << beamPIDs3GeV.electron
                            << "  " << beamPIDs3GeV.muon
                            << "  " << beamPIDs3GeV.pion
                            << "  " << beamPIDs3GeV.kaon
                            << "  " << beamPIDs3GeV.proton
                            << "  " << beamPIDs3GeV.deuteron
                                    << "\n"
                << "6 GeV:        " << beamPIDs6GeV.electron
                            << "  " << beamPIDs6GeV.muon
                            << "  " << beamPIDs6GeV.pion
                            << "  " << beamPIDs6GeV.kaon
                            << "  " << beamPIDs6GeV.proton
                            << "  " << beamPIDs6GeV.deuteron
                                    << "\n"
                << "7 GeV:        " << beamPIDs7GeV.electron
                            << "  " << beamPIDs7GeV.muon
                            << "  " << beamPIDs7GeV.pion
                            << "  " << beamPIDs7GeV.kaon
                            << "  " << beamPIDs7GeV.proton
                            << "  " << beamPIDs7GeV.deuteron
                                    << "\n";
    }
  }


  nBeamEventsOld = beamVecOld.size();
  nBeamTracksOld = 0;
  nBeamMomOld = 0;
  for(size_t iBeamEvent=0; iBeamEvent < beamVecOld.size(); iBeamEvent++)
  {
    beam::ProtoDUNEBeamEvent beamEvent = *(beamVecOld.at(iBeamEvent));
    if(PRINTBEAMEVENT)
    {
      std::cout << "PiAbsSel BeamEvent (OLD)" << iBeamEvent << " for " << e.id() << ": \n";
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

    BITriggerOld = beamEvent.GetBITrigger();
    BITriggerMatchedOld = beamEvent.CheckIsMatched();
    BIActiveTriggerOld = beamEvent.GetActiveTrigger();
    TOFOld = beamEvent.GetTOF();
    CKov0StatusOld = beamEvent.GetCKov0Status();
    CKov1StatusOld = beamEvent.GetCKov1Status();

    for(size_t iMom=0; iMom < beamEvent.GetNRecoBeamMomenta(); iMom++)
    {
      beamMomOld[nBeamMomOld] = beamEvent.GetRecoBeamMomentum(iMom);
      nBeamMomOld++;
      if(PRINTBEAMEVENT) std::cout << "    " << beamEvent.GetRecoBeamMomentum(iMom) << "\n";
    }
    for(size_t iTrack=0; iTrack < beamEvent.GetNBeamTracks(); iTrack++)
    {
      const recob::Track & track =  beamEvent.GetBeamTrack(iTrack);
      if(nBeamTracks >= MAXBEAMTRACKS)
      {
        mf::LogError("PionAbsSelector") << "Event " << e.id() << " not stored because too many old beam tracks";
        return;
      }
      beamTrackXFrontTPCOld[nBeamTracksOld] = track.End().X();
      beamTrackYFrontTPCOld[nBeamTracksOld] = track.End().Y();
      beamTrackThetaOld[nBeamTracksOld] = track.EndDirection().Theta();
      beamTrackPhiOld[nBeamTracksOld] = track.EndDirection().Phi();
      nBeamTracksOld++;

      if(PRINTBEAMEVENT)
      {
        std::cout << "  Old Beam Track: "<< iTrack <<"\n";
        std::cout << "    N Points:  " << track.NumberTrajectoryPoints() << "\n";
        std::cout << "    Start Pos: " << track.Vertex().X()
                                  << "  " << track.Vertex().Y()
                                  << "  " << track.Vertex().Z() << "\n";
        std::cout << "    Mid Pos:   " << track.LocationAtPoint(1).X()
                                  << "  " << track.LocationAtPoint(1).Y()
                                  << "  " << track.LocationAtPoint(1).Z() << "\n";
        std::cout << "    End Pos:   " << track.End().X()
                                  << "  " << track.End().Y()
                                  << "  " << track.End().Z() << "\n";
        std::cout << "    Start Theta: " << track.VertexDirection().Theta()*180/CLHEP::pi << " deg\n";
        std::cout << "    Start Phi:   " << track.VertexDirection().Phi()*180/CLHEP::pi << " deg\n";
        std::cout << "    End Theta:   " << track.EndDirection().Theta()*180/CLHEP::pi << " deg\n";
        std::cout << "    End Phi:     " << track.EndDirection().Phi()*180/CLHEP::pi << " deg\n";
      }
    }

    // CERN official beam PID
    //const auto& beamPIDs0p3GeV = fBeamlineUtilsOld.GetPIDCandidates(beamEvent,0.3);
    //const auto& beamPIDs0p5GeV = fBeamlineUtilsOld.GetPIDCandidates(beamEvent,0.5);
    const auto& beamPIDs1GeV = fBeamlineUtilsOld.GetPIDCandidates(beamEvent,1);
    const auto& beamPIDs2GeV = fBeamlineUtilsOld.GetPIDCandidates(beamEvent,2);
    const auto& beamPIDs3GeV = fBeamlineUtilsOld.GetPIDCandidates(beamEvent,3);
    const auto& beamPIDs6GeV = fBeamlineUtilsOld.GetPIDCandidates(beamEvent,6);
    const auto& beamPIDs7GeV = fBeamlineUtilsOld.GetPIDCandidates(beamEvent,7);

    //BIElectron0p3GeVOld = beamPIDs0p3GeV.electron;
    //BIMuon0p3GeVOld = beamPIDs0p3GeV.muon;
    //BIPion0p3GeVOld = beamPIDs0p3GeV.pion;
    //BIKaon0p3GeVOld = beamPIDs0p3GeV.kaon;
    //BIProton0p3GeVOld = beamPIDs0p3GeV.proton;
    //BIDeuteron0p3GeVOld = beamPIDs0p3GeV.deuteron;
    //BIElectron0p5GeVOld = beamPIDs0p5GeV.electron;
    //BIMuon0p5GeVOld = beamPIDs0p5GeV.muon;
    //BIPion0p5GeVOld = beamPIDs0p5GeV.pion;
    //BIKaon0p5GeVOld = beamPIDs0p5GeV.kaon;
    //BIProton0p5GeVOld = beamPIDs0p5GeV.proton;
    //BIDeuteron0p5GeVOld = beamPIDs0p5GeV.deuteron;
    BIElectron1GeVOld = beamPIDs1GeV.electron;
    BIMuon1GeVOld = beamPIDs1GeV.muon;
    BIPion1GeVOld = beamPIDs1GeV.pion;
    BIKaon1GeVOld = beamPIDs1GeV.kaon;
    BIProton1GeVOld = beamPIDs1GeV.proton;
    BIDeuteron1GeVOld = beamPIDs1GeV.deuteron;
    BIElectron2GeVOld = beamPIDs2GeV.electron;
    BIMuon2GeVOld = beamPIDs2GeV.muon;
    BIPion2GeVOld = beamPIDs2GeV.pion;
    BIKaon2GeVOld = beamPIDs2GeV.kaon;
    BIProton2GeVOld = beamPIDs2GeV.proton;
    BIDeuteron2GeVOld = beamPIDs2GeV.deuteron;
    BIElectron3GeVOld = beamPIDs3GeV.electron;
    BIMuon3GeVOld = beamPIDs3GeV.muon;
    BIPion3GeVOld = beamPIDs3GeV.pion;
    BIKaon3GeVOld = beamPIDs3GeV.kaon;
    BIProton3GeVOld = beamPIDs3GeV.proton;
    BIDeuteron3GeVOld = beamPIDs3GeV.deuteron;
    BIElectron6GeVOld = beamPIDs6GeV.electron;
    BIMuon6GeVOld = beamPIDs6GeV.muon;
    BIPion6GeVOld = beamPIDs6GeV.pion;
    BIKaon6GeVOld = beamPIDs6GeV.kaon;
    BIProton6GeVOld = beamPIDs6GeV.proton;
    BIDeuteron6GeVOld = beamPIDs6GeV.deuteron;
    BIElectron7GeVOld = beamPIDs7GeV.electron;
    BIMuon7GeVOld = beamPIDs7GeV.muon;
    BIPion7GeVOld = beamPIDs7GeV.pion;
    BIKaon7GeVOld = beamPIDs7GeV.kaon;
    BIProton7GeVOld = beamPIDs7GeV.proton;
    BIDeuteron7GeVOld = beamPIDs7GeV.deuteron;

    if(PRINTBEAMEVENT)
    {
      std::cout << "My   old PID: e  mu pi k  p  d  \n"
                << "1 GeV:        " << beamPIDs1GeV.electron
                            << "  " << beamPIDs1GeV.muon
                            << "  " << beamPIDs1GeV.pion
                            << "  " << beamPIDs1GeV.kaon
                            << "  " << beamPIDs1GeV.proton
                            << "  " << beamPIDs1GeV.deuteron
                                    << "\n"
                << "2 GeV:        " << beamPIDs2GeV.electron
                            << "  " << beamPIDs2GeV.muon
                            << "  " << beamPIDs2GeV.pion
                            << "  " << beamPIDs2GeV.kaon
                            << "  " << beamPIDs2GeV.proton
                            << "  " << beamPIDs2GeV.deuteron
                                    << "\n"
                << "3 GeV:        " << beamPIDs3GeV.electron
                            << "  " << beamPIDs3GeV.muon
                            << "  " << beamPIDs3GeV.pion
                            << "  " << beamPIDs3GeV.kaon
                            << "  " << beamPIDs3GeV.proton
                            << "  " << beamPIDs3GeV.deuteron
                                    << "\n"
                << "6 GeV:        " << beamPIDs6GeV.electron
                            << "  " << beamPIDs6GeV.muon
                            << "  " << beamPIDs6GeV.pion
                            << "  " << beamPIDs6GeV.kaon
                            << "  " << beamPIDs6GeV.proton
                            << "  " << beamPIDs6GeV.deuteron
                                    << "\n"
                << "7 GeV:        " << beamPIDs7GeV.electron
                            << "  " << beamPIDs7GeV.muon
                            << "  " << beamPIDs7GeV.pion
                            << "  " << beamPIDs7GeV.kaon
                            << "  " << beamPIDs7GeV.proton
                            << "  " << beamPIDs7GeV.deuteron
                                    << "\n";
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

    for(size_t k=0; k < 6; k++)
    {
      nGoodFEMBs[k] = fDataUtil.GetNActiveFembsForAPA(e, k);
    }
  }

  const art::Ptr<simb::MCParticle> primaryParticle = ProcessMCParticles(truePartVec,beamOrCosmic,e);

  // Either got pWC from beam or from primaryParticle, so now is the time to do this
  eWC = sqrt(pWC*pWC+MCHARGEDPION*MCHARGEDPION); // assume charged pion in MeV
  kinWC = eWC - MCHARGEDPION; // assume charged pion in MeV
  kinWCInTPC = kinWC - KINLOSTBEFORETPC;

  // for proton
  eWCProton = sqrt(pWC*pWC+MPROTON*MPROTON);
  kinWCProton = eWCProton - MPROTON;
  kinWCInTPCProton = kinWCProton - KINLOSTBEFORETPCPROTON;

  if(fDoPMAAnalysis)
  {
    ProcessAllTracks(trackVec, truePartVec,tracksCaloVec,fmHitsForTracks,beamOrCosmic,e);

    ProcessPrimaryTrack(trackVec,primaryParticle,tracksCaloVec,motherDaughterWalker,e);

    if (isMC && nTracksInFirstZ[2] >= 1 && nTracksInFirstZ[14] < 4 && nTracksLengthLt[5] < 3
              && iBestMatch >= 0 && nMatchedTracks == 1
              && primTrkEndX > 5.4 && primTrkEndX < 42.9
              && primTrkEndY > -15. && primTrkEndX < 15.
              && primTrkEndY > 5. && primTrkEndX < 85.
        )
    {
      mf::LogInfo("PiAbsRecoAccuracy") << std::fixed << std::setprecision(1) 
              << "N_pi+/-: " << trueNSecondaryChPions
              << " N_pi0: " << trueNSecondaryPiZeros
              << " N_p: " << trueNSecondaryProtons
              << " N Secondary Tracks: " << nSecTrk
              << std::setprecision(2)
              << " Min Proton/Pion LHR: " << nSecLLRProtonToPionMin;
    }
  } // if fDoPMAAnalysis

  if (fDoPandoraAnalysis) ProcessPFParticles(e,beamOrCosmic);
  
  tree->Fill();

  if (fDoPMAAnalysis)
  {
    const float flangeRadiusCutSquared = pow(fFlangeRadiusCut,2);
    for (size_t iTrack=0; iTrack < nTracks; iTrack++)
    {
      for (size_t iBeamTrack=0; iBeamTrack < nBeamTracks; iBeamTrack++)
      {
        const float dx = trackXFrontTPC[iTrack] - beamTrackXFrontTPC[iBeamTrack];
        const float dy = trackYFrontTPC[iTrack] - beamTrackYFrontTPC[iBeamTrack];
        TVector3 trackDir;
        TVector3 beamTrackDir;
        trackDir.SetMagThetaPhi(1.,trackStartTheta[iTrack],trackStartPhi[iTrack]);
        beamTrackDir.SetMagThetaPhi(1.,beamTrackTheta[iTrack],beamTrackPhi[iTrack]);
        const float dAngle = trackDir.Angle(beamTrackDir)*180./CLHEP::pi;
        deltaXYTPCBeamlineHist->Fill(dx,dy);
        deltaAngleTPCBeamlineHist->Fill(dAngle);
        const float r2FlangeBeamTrack = pow(beamTrackXFrontTPC[iBeamTrack]-fFlangeCenterX,2)
                                      +pow(beamTrackYFrontTPC[iBeamTrack]-fFlangeCenterY,2);
        const float r2FlangeTrack = pow(trackXFrontTPC[iTrack]-fFlangeCenterX,2)
                                     +pow(trackYFrontTPC[iTrack]-fFlangeCenterY,2);
        const bool beamTrackInFlange = r2FlangeBeamTrack < flangeRadiusCutSquared;
        const bool trackInFlange = r2FlangeTrack < flangeRadiusCutSquared;
        const bool trackInFirst25cm = std::min(trackStartZ[iTrack],trackEndZ[iTrack]) < 25.;
        //std::cout << "iTrack: " << iTrack << " iMCPart: " << iMCPart
        //          << " trackInFirst25cm " << trackInFirst25cm
        //          << " trackInFlange " << trackInFlange
        //          << " beamTrackInFlange " << beamTrackInFlange
        //          << " r2FlangeBeamTrack " << r2FlangeBeamTrack
        //          << " r2FlangeTrack " << r2FlangeTrack
        //          << " trackMinZ " << std::min(trackStartZ[iTrack],trackEndZ[iTrack])
        //          << std::endl;
        if (beamTrackInFlange && trackInFlange)
        {
          deltaXYTPCBeamlineOnlyInFlangeHist->Fill(dx,dy);
          deltaAngleTPCBeamlineOnlyInFlangeHist->Fill(dAngle);
          if (trackInFirst25cm)
          {
            deltaXYTPCBeamlineOnlyInFlangeInFirst25cmHist->Fill(dx,dy);
            deltaAngleTPCBeamlineOnlyInFlangeInFirst25cmHist->Fill(dAngle);
          }
        }
        if (trackInFirst25cm)
        {
          deltaXYTPCBeamlineOnlyInFirst25cmHist->Fill(dx,dy);
          deltaAngleTPCBeamlineOnlyInFirst25cmHist->Fill(dAngle);
        }
      } // for iMCPart
    } // for iTrack
  } // if fDoPMAAnalysis

  if(beamOrCosmic) delete beamOrCosmic;

} // analyze function

void lana::PionAbsSelector::beginJob()
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
  tree->Branch("kinWCInTPC",&kinWCInTPC,"kinWCInTPC/F");
  tree->Branch("eWCProton",&eWCProton,"eWCProton/F");
  tree->Branch("kinWCProton",&kinWCProton,"kinWCProton/F");
  tree->Branch("kinWCInTPCProton",&kinWCInTPCProton,"kinWCInTPCProton/F");
  tree->Branch("xWC1Hit",&xWC1Hit,"xWC1Hit/F");
  tree->Branch("yWC1Hit",&yWC1Hit,"yWC1Hit/F");
  tree->Branch("zWC1Hit",&zWC1Hit,"zWC1Hit/F");
  tree->Branch("xWC2Hit",&xWC2Hit,"xWC2Hit/F");
  tree->Branch("yWC2Hit",&yWC2Hit,"yWC2Hit/F");
  tree->Branch("zWC2Hit",&zWC2Hit,"zWC2Hit/F");

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
  tree->Branch("TOFsByChan",&TOFsByChan,"TOFsByChan[4]/F");
  tree->Branch("TOFusTrigsByChan",&TOFusTrigsByChan,"TOFusTrigsByChan[4]/i");
  tree->Branch("TOFdsTrigsByChan",&TOFdsTrigsByChan,"TOFdsTrigsByChan[4]/i");

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

  tree->Branch("BIElectron0p3GeV",&BIElectron0p3GeV,"BIElectron0p3GeV/O");
  tree->Branch("BIMuon0p3GeV",&BIMuon0p3GeV,"BIMuon0p3GeV/O");
  tree->Branch("BIPion0p3GeV",&BIPion0p3GeV,"BIPion0p3GeV/O");
  tree->Branch("BIKaon0p3GeV",&BIKaon0p3GeV,"BIKaon0p3GeV/O");
  tree->Branch("BIProton0p3GeV",&BIProton0p3GeV,"BIProton0p3GeV/O");
  tree->Branch("BIDeuteron0p3GeV",&BIDeuteron0p3GeV,"BIDeuteron0p3GeV/O");
  tree->Branch("BIElectron0p5GeV",&BIElectron0p5GeV,"BIElectron0p5GeV/O");
  tree->Branch("BIMuon0p5GeV",&BIMuon0p5GeV,"BIMuon0p5GeV/O");
  tree->Branch("BIPion0p5GeV",&BIPion0p5GeV,"BIPion0p5GeV/O");
  tree->Branch("BIKaon0p5GeV",&BIKaon0p5GeV,"BIKaon0p5GeV/O");
  tree->Branch("BIProton0p5GeV",&BIProton0p5GeV,"BIProton0p5GeV/O");
  tree->Branch("BIDeuteron0p5GeV",&BIDeuteron0p5GeV,"BIDeuteron0p5GeV/O");
  tree->Branch("BIElectron1GeV",&BIElectron1GeV,"BIElectron1GeV/O");
  tree->Branch("BIMuon1GeV",&BIMuon1GeV,"BIMuon1GeV/O");
  tree->Branch("BIPion1GeV",&BIPion1GeV,"BIPion1GeV/O");
  tree->Branch("BIKaon1GeV",&BIKaon1GeV,"BIKaon1GeV/O");
  tree->Branch("BIProton1GeV",&BIProton1GeV,"BIProton1GeV/O");
  tree->Branch("BIDeuteron1GeV",&BIDeuteron1GeV,"BIDeuteron1GeV/O");
  tree->Branch("BIElectron2GeV",&BIElectron2GeV,"BIElectron2GeV/O");
  tree->Branch("BIMuon2GeV",&BIMuon2GeV,"BIMuon2GeV/O");
  tree->Branch("BIPion2GeV",&BIPion2GeV,"BIPion2GeV/O");
  tree->Branch("BIKaon2GeV",&BIKaon2GeV,"BIKaon2GeV/O");
  tree->Branch("BIProton2GeV",&BIProton2GeV,"BIProton2GeV/O");
  tree->Branch("BIDeuteron2GeV",&BIDeuteron2GeV,"BIDeuteron2GeV/O");
  tree->Branch("BIElectron3GeV",&BIElectron3GeV,"BIElectron3GeV/O");
  tree->Branch("BIMuon3GeV",&BIMuon3GeV,"BIMuon3GeV/O");
  tree->Branch("BIPion3GeV",&BIPion3GeV,"BIPion3GeV/O");
  tree->Branch("BIKaon3GeV",&BIKaon3GeV,"BIKaon3GeV/O");
  tree->Branch("BIProton3GeV",&BIProton3GeV,"BIProton3GeV/O");
  tree->Branch("BIDeuteron3GeV",&BIDeuteron3GeV,"BIDeuteron3GeV/O");
  tree->Branch("BIElectron6GeV",&BIElectron6GeV,"BIElectron6GeV/O");
  tree->Branch("BIMuon6GeV",&BIMuon6GeV,"BIMuon6GeV/O");
  tree->Branch("BIPion6GeV",&BIPion6GeV,"BIPion6GeV/O");
  tree->Branch("BIKaon6GeV",&BIKaon6GeV,"BIKaon6GeV/O");
  tree->Branch("BIProton6GeV",&BIProton6GeV,"BIProton6GeV/O");
  tree->Branch("BIDeuteron6GeV",&BIDeuteron6GeV,"BIDeuteron6GeV/O");
  tree->Branch("BIElectron7GeV",&BIElectron7GeV,"BIElectron7GeV/O");
  tree->Branch("BIMuon7GeV",&BIMuon7GeV,"BIMuon7GeV/O");
  tree->Branch("BIPion7GeV",&BIPion7GeV,"BIPion7GeV/O");
  tree->Branch("BIKaon7GeV",&BIKaon7GeV,"BIKaon7GeV/O");
  tree->Branch("BIProton7GeV",&BIProton7GeV,"BIProton7GeV/O");
  tree->Branch("BIDeuteron7GeV",&BIDeuteron7GeV,"BIDeuteron7GeV/O");

  tree->Branch("nBeamEventsOld",&nBeamEventsOld,"nBeamEventsOld/i");
  tree->Branch("BITriggerOld",&BITriggerOld,"BITriggerOld/I");
  tree->Branch("BITriggerMatchedOld",&BITriggerMatchedOld,"BITriggerMatchedOld/O");
  tree->Branch("BIActiveTriggerOld",&BIActiveTriggerOld,"BIActiveTriggerOld/i");
  tree->Branch("TOFOld",&TOFOld,"TOFOld/F");

  tree->Branch("CKov0StatusOld",&CKov0StatusOld,"CKov0StatusOld/I");
  tree->Branch("CKov1StatusOld",&CKov1StatusOld,"CKov1StatusOld/I");

  tree->Branch("nBeamMomOld",&nBeamMomOld,"nBeamMomOld/i");
  tree->Branch("beamMomOld",&beamMomOld,"beamMomOld[nBeamMomOld]/F");

  tree->Branch("nBeamTracksOld",&nBeamTracksOld,"nBeamTracksOld/i");
  tree->Branch("beamTrackXFrontTPCOld",&beamTrackXFrontTPCOld,"beamTrackXFrontTPCOld[nBeamTracksOld]/F");
  tree->Branch("beamTrackYFrontTPCOld",&beamTrackYFrontTPCOld,"beamTrackYFrontTPCOld[nBeamTracksOld]/F");
  tree->Branch("beamTrackThetaOld",&beamTrackThetaOld,"beamTrackThetaOld[nBeamTracksOld]/F");
  tree->Branch("beamTrackPhiOld",&beamTrackPhiOld,"beamTrackPhiOld[nBeamTracksOld]/F");

  tree->Branch("BIElectron0p3GeVOld",&BIElectron0p3GeVOld,"BIElectron0p3GeVOld/O");
  tree->Branch("BIMuon0p3GeVOld",&BIMuon0p3GeVOld,"BIMuon0p3GeVOld/O");
  tree->Branch("BIPion0p3GeVOld",&BIPion0p3GeVOld,"BIPion0p3GeVOld/O");
  tree->Branch("BIKaon0p3GeVOld",&BIKaon0p3GeVOld,"BIKaon0p3GeVOld/O");
  tree->Branch("BIProton0p3GeVOld",&BIProton0p3GeVOld,"BIProton0p3GeVOld/O");
  tree->Branch("BIDeuteron0p3GeVOld",&BIDeuteron0p3GeVOld,"BIDeuteron0p3GeVOld/O");
  tree->Branch("BIElectron0p5GeVOld",&BIElectron0p5GeVOld,"BIElectron0p5GeVOld/O");
  tree->Branch("BIMuon0p5GeVOld",&BIMuon0p5GeVOld,"BIMuon0p5GeVOld/O");
  tree->Branch("BIPion0p5GeVOld",&BIPion0p5GeVOld,"BIPion0p5GeVOld/O");
  tree->Branch("BIKaon0p5GeVOld",&BIKaon0p5GeVOld,"BIKaon0p5GeVOld/O");
  tree->Branch("BIProton0p5GeVOld",&BIProton0p5GeVOld,"BIProton0p5GeVOld/O");
  tree->Branch("BIDeuteron0p5GeVOld",&BIDeuteron0p5GeVOld,"BIDeuteron0p5GeVOld/O");
  tree->Branch("BIElectron1GeVOld",&BIElectron1GeVOld,"BIElectron1GeVOld/O");
  tree->Branch("BIMuon1GeVOld",&BIMuon1GeVOld,"BIMuon1GeVOld/O");
  tree->Branch("BIPion1GeVOld",&BIPion1GeVOld,"BIPion1GeVOld/O");
  tree->Branch("BIKaon1GeVOld",&BIKaon1GeVOld,"BIKaon1GeVOld/O");
  tree->Branch("BIProton1GeVOld",&BIProton1GeVOld,"BIProton1GeVOld/O");
  tree->Branch("BIDeuteron1GeVOld",&BIDeuteron1GeVOld,"BIDeuteron1GeVOld/O");
  tree->Branch("BIElectron2GeVOld",&BIElectron2GeVOld,"BIElectron2GeVOld/O");
  tree->Branch("BIMuon2GeVOld",&BIMuon2GeVOld,"BIMuon2GeVOld/O");
  tree->Branch("BIPion2GeVOld",&BIPion2GeVOld,"BIPion2GeVOld/O");
  tree->Branch("BIKaon2GeVOld",&BIKaon2GeVOld,"BIKaon2GeVOld/O");
  tree->Branch("BIProton2GeVOld",&BIProton2GeVOld,"BIProton2GeVOld/O");
  tree->Branch("BIDeuteron2GeVOld",&BIDeuteron2GeVOld,"BIDeuteron2GeVOld/O");
  tree->Branch("BIElectron3GeVOld",&BIElectron3GeVOld,"BIElectron3GeVOld/O");
  tree->Branch("BIMuon3GeVOld",&BIMuon3GeVOld,"BIMuon3GeVOld/O");
  tree->Branch("BIPion3GeVOld",&BIPion3GeVOld,"BIPion3GeVOld/O");
  tree->Branch("BIKaon3GeVOld",&BIKaon3GeVOld,"BIKaon3GeVOld/O");
  tree->Branch("BIProton3GeVOld",&BIProton3GeVOld,"BIProton3GeVOld/O");
  tree->Branch("BIDeuteron3GeVOld",&BIDeuteron3GeVOld,"BIDeuteron3GeVOld/O");
  tree->Branch("BIElectron6GeVOld",&BIElectron6GeVOld,"BIElectron6GeVOld/O");
  tree->Branch("BIMuon6GeVOld",&BIMuon6GeVOld,"BIMuon6GeVOld/O");
  tree->Branch("BIPion6GeVOld",&BIPion6GeVOld,"BIPion6GeVOld/O");
  tree->Branch("BIKaon6GeVOld",&BIKaon6GeVOld,"BIKaon6GeVOld/O");
  tree->Branch("BIProton6GeVOld",&BIProton6GeVOld,"BIProton6GeVOld/O");
  tree->Branch("BIDeuteron6GeVOld",&BIDeuteron6GeVOld,"BIDeuteron6GeVOld/O");
  tree->Branch("BIElectron7GeVOld",&BIElectron7GeVOld,"BIElectron7GeVOld/O");
  tree->Branch("BIMuon7GeVOld",&BIMuon7GeVOld,"BIMuon7GeVOld/O");
  tree->Branch("BIPion7GeVOld",&BIPion7GeVOld,"BIPion7GeVOld/O");
  tree->Branch("BIKaon7GeVOld",&BIKaon7GeVOld,"BIKaon7GeVOld/O");
  tree->Branch("BIProton7GeVOld",&BIProton7GeVOld,"BIProton7GeVOld/O");
  tree->Branch("BIDeuteron7GeVOld",&BIDeuteron7GeVOld,"BIDeuteron7GeVOld/O");

  tree->Branch("triggerIsBeam",&triggerIsBeam,"triggerIsBeam/O");
  tree->Branch("triggerBits",&triggerBits,"triggerBits/i");
  tree->Branch("nGoodFEMBs",&nGoodFEMBs,"nGoodFEMBs[6]/I");

  tree->Branch("nPrimaryParticleCandidates",&nPrimaryParticleCandidates,"nPrimaryParticleCandidates/i");

  tree->Branch("trueCategory",&trueCategory,"trueCategory/I");
  tree->Branch("trueEndProcess",&trueEndProcess,"trueEndProcess/I");
  tree->Branch("truePrimaryPDG",&truePrimaryPDG,"truePrimaryPDG/I");
  tree->Branch("truePrimaryTrackID",&truePrimaryTrackID,"truePrimaryTrackID/I");
  tree->Branch("trueSignalT",&trueSignalT,"trueSignalT/O");
  tree->Branch("trueSignalNT",&trueSignalNT,"trueSignalNT/O");
  tree->Branch("trueNDaughters",&trueNDaughters,"trueNDaughters/i");
  tree->Branch("nSecTracks",&nSecTracks,"nSecTracks/i");
  tree->Branch("trueNSecondaryChPions",&trueNSecondaryChPions,"trueNSecondaryChPions/i");
  tree->Branch("trueNSecondaryPiZeros",&trueNSecondaryPiZeros,"trueNSecondaryPiZeros/i");
  tree->Branch("trueNSecondaryProtons",&trueNSecondaryProtons,"trueNSecondaryProtons/i");
  tree->Branch("trueNSecondaryOppChPions",&trueNSecondaryOppChPions,"trueNSecondaryOppChPions/i");
  tree->Branch("trueNSecondaryMuons",&trueNSecondaryMuons,"trueNSecondaryMuons/i");
  tree->Branch("trueStartX",&trueStartX,"trueStartX/F");
  tree->Branch("trueStartY",&trueStartY,"trueStartY/F");
  tree->Branch("trueStartZ",&trueStartZ,"trueStartZ/F");
  tree->Branch("trueStartT",&trueStartT,"trueStartT/F");
  tree->Branch("trueEndX",&trueEndX,"trueEndX/F");
  tree->Branch("trueEndY",&trueEndY,"trueEndY/F");
  tree->Branch("trueEndZ",&trueEndZ,"trueEndZ/F");
  tree->Branch("trueEndT",&trueEndT,"trueEndT/F");
  tree->Branch("trueStartTheta",&trueStartTheta,"trueStartTheta/F");
  tree->Branch("trueStartPhi",&trueStartPhi,"trueStartPhi/F");
  tree->Branch("trueSecondToEndTheta",&trueSecondToEndTheta,"trueSecondToEndTheta/F");
  tree->Branch("trueSecondToEndPhi",&trueSecondToEndPhi,"trueSecondToEndPhi/F");
  tree->Branch("trueStartMom",&trueStartMom,"trueStartMom/F");
  tree->Branch("trueStartE",&trueStartE,"trueStartE/F");
  tree->Branch("trueStartKin",&trueStartKin,"trueStartKin/F");
  tree->Branch("trueEndMom",&trueEndMom,"trueEndMom/F");
  tree->Branch("trueEndE",&trueEndE,"trueEndE/F");
  tree->Branch("trueEndKin",&trueEndKin,"trueEndKin/F");
  tree->Branch("trueEndMaxAngleToSec",&trueEndMaxAngleToSec,"trueEndMaxAngleToSec/F");
  tree->Branch("trueSecondToEndMom",&trueSecondToEndMom,"trueSecondToEndMom/F");
  tree->Branch("trueSecondToEndE",&trueSecondToEndE,"trueSecondToEndE/F");
  tree->Branch("trueSecondToEndKin",&trueSecondToEndKin,"trueSecondToEndKin/F");
  tree->Branch("trueAngleStartDirSecondToEndDir",&trueAngleStartDirSecondToEndDir,"trueAngleStartDirSecondToEndDir/F");
  tree->Branch("trueSecondPDG",&trueSecondPDG,"trueSecondPDG[nSecTracks]/I");
  tree->Branch("trueSecondKin",&trueSecondKin,"trueSecondKin[nSecTracks]/F");
  tree->Branch("trueSecondAngleToPrimEnd",&trueSecondAngleToPrimEnd,"trueSecondAngleToPrimEnd[nSecTracks]/F");
  tree->Branch("trueXFrontTPC",&trueXFrontTPC,"trueXFrontTPC/F");
  tree->Branch("trueYFrontTPC",&trueYFrontTPC,"trueYFrontTPC/F");
  tree->Branch("trueTrajPointFrontTPCDist",&trueTrajPointFrontTPCDist,"trueTrajPointFrontTPCDist/F");
  tree->Branch("trueKinFrontTPC",&trueKinFrontTPC,"trueKinFrontTPC/F");

  tree->Branch("trueTrajNPoints",&trueTrajNPoints,"trueTrajNPoints/i");
  tree->Branch("trueTrajPointX",&trueTrajPointX,"trueTrajPointX[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointY",&trueTrajPointY,"trueTrajPointY[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointZ",&trueTrajPointZ,"trueTrajPointZ[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointT",&trueTrajPointT,"trueTrajPointT[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointPx",&trueTrajPointPx,"trueTrajPointPx[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointPy",&trueTrajPointPy,"trueTrajPointPy[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointPz",&trueTrajPointPz,"trueTrajPointPz[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointE",&trueTrajPointE,"trueTrajPointE[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointMom",&trueTrajPointMom,"trueTrajPointMom[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointKin",&trueTrajPointKin,"trueTrajPointKin[trueTrajNPoints]/F");
  tree->Branch("trueTrajPointProc",&trueTrajPointProc,"trueTrajPointProc[trueTrajNPoints]/b");
  tree->Branch("trueTrajPointKinkAngle",&trueTrajPointKinkAngle,"trueTrajPointKinkAngle[trueTrajNPoints]/F");

  tree->Branch("nMCParts",&nMCParts,"nMCParts/i");
  tree->Branch("mcPartIsBeam",&mcPartIsBeam,"mcPartIsBeam[nMCParts]/O");
  tree->Branch("mcPartIsPrimary",&mcPartIsPrimary,"mcPartIsPrimary[nMCParts]/O");
  tree->Branch("mcPartTrackID",&mcPartTrackID,"mcPartTrackID[nMCParts]/I");
  tree->Branch("mcPartPDG",&mcPartPDG,"mcPartPDG[nMCParts]/I");
  tree->Branch("mcPartStartX",&mcPartStartX,"mcPartStartX[nMCParts]/F");
  tree->Branch("mcPartStartY",&mcPartStartY,"mcPartStartY[nMCParts]/F");
  tree->Branch("mcPartStartZ",&mcPartStartZ,"mcPartStartZ[nMCParts]/F");
  tree->Branch("mcPartStartT",&mcPartStartT,"mcPartStartT[nMCParts]/F");
  tree->Branch("mcPartEndX",&mcPartEndX,"mcPartEndX[nMCParts]/F");
  tree->Branch("mcPartEndY",&mcPartEndY,"mcPartEndY[nMCParts]/F");
  tree->Branch("mcPartEndZ",&mcPartEndZ,"mcPartEndZ[nMCParts]/F");
  tree->Branch("mcPartEndT",&mcPartEndT,"mcPartEndT[nMCParts]/F");
  tree->Branch("mcPartStartTheta",&mcPartStartTheta,"mcPartStartTheta[nMCParts]/F");
  tree->Branch("mcPartStartPhi",&mcPartStartPhi,"mcPartStartPhi[nMCParts]/F");
  tree->Branch("mcPartXFrontTPC",&mcPartXFrontTPC,"mcPartXFrontTPC[nMCParts]/F");
  tree->Branch("mcPartYFrontTPC",&mcPartYFrontTPC,"mcPartYFrontTPC[nMCParts]/F");
  tree->Branch("mcPartStartMom",&mcPartStartMom,"mcPartStartMom[nMCParts]/F");
  tree->Branch("mcPartStartE",&mcPartStartE,"mcPartStartE[nMCParts]/F");
  tree->Branch("mcPartStartKin",&mcPartStartKin,"mcPartStartKin[nMCParts]/F");
  tree->Branch("mcPartEndMom",&mcPartEndMom,"mcPartEndMom[nMCParts]/F");
  tree->Branch("mcPartEndE",&mcPartEndE,"mcPartEndE[nMCParts]/F");
  tree->Branch("mcPartEndKin",&mcPartEndKin,"mcPartEndKin[nMCParts]/F");
  tree->Branch("mcPartDeltaAngle",&mcPartDeltaAngle,"mcPartDeltaAngle[nMCParts]/F");

  tree->Branch("nIDEs",&nIDEs,"nIDEs/i");
  tree->Branch("simIDENumElectrons",&simIDENumElectrons,"simIDENumElectrons[nIDEs]/F");
  tree->Branch("simIDEEnergy",&simIDEEnergy,"simIDEEnergy[nIDEs]/F");
  tree->Branch("simIDEX",&simIDEX,"simIDEX[nIDEs]/F");
  tree->Branch("simIDEY",&simIDEY,"simIDEY[nIDEs]/F");
  tree->Branch("simIDEZ",&simIDEZ,"simIDEZ[nIDEs]/F");
  tree->Branch("simIDETrajDistance",&simIDETrajDistance,"simIDETrajDistance[nIDEs]/F");
  tree->Branch("simIDETrajKin",&simIDETrajKin,"simIDETrajKin[nIDEs]/F");
  tree->Branch("simIDEWireZ",&simIDEWireZ,"simIDEWireZ[nIDEs]/F");
  tree->Branch("simIDEChannel",&simIDEChannel,"simIDEChannel[nIDEs]/i");
  tree->Branch("simIDETDC",&simIDETDC,"simIDETDC[nIDEs]/F");
  tree->Branch("simIDEPartKin",&simIDEPartKin,"simIDEPartKin[nIDEs]/F");
  tree->Branch("simIDEEnergySum",&simIDEEnergySum,"simIDEEnergySum/F");
  
  tree->Branch("nTracks",&nTracks,"nTracks/i");
  tree->Branch("nTracksInFirstZ",&nTracksInFirstZ,("nTracksInFirstZ["+std::to_string(MAXZINT)+"]/i").c_str());
  tree->Branch("nTracksLengthLt",&nTracksLengthLt,("nTracksLengthLt["+std::to_string(MAXLINT)+"]/i").c_str());
  tree->Branch("trackStartX",&trackStartX,"trackStartX[nTracks]/F");
  tree->Branch("trackStartY",&trackStartY,"trackStartY[nTracks]/F");
  tree->Branch("trackStartZ",&trackStartZ,"trackStartZ[nTracks]/F");
  tree->Branch("trackStartTheta",&trackStartTheta,"trackStartTheta[nTracks]/F");
  tree->Branch("trackStartPhi",&trackStartPhi,"trackStartPhi[nTracks]/F");
  tree->Branch("trackEndX",&trackEndX,"trackEndX[nTracks]/F");
  tree->Branch("trackEndY",&trackEndY,"trackEndY[nTracks]/F");
  tree->Branch("trackEndZ",&trackEndZ,"trackEndZ[nTracks]/F");
  tree->Branch("trackLength",&trackLength,"trackLength[nTracks]/F");
  tree->Branch("trackXFrontTPC",&trackXFrontTPC,"trackXFrontTPC[nTracks]/F");
  tree->Branch("trackYFrontTPC",&trackYFrontTPC,"trackYFrontTPC[nTracks]/F");
  tree->Branch("trackCaloKin",&trackCaloKin,"trackCaloKin[nTracks]/F");
  tree->Branch("trackLLHPion",&trackLLHPion,"trackLLHPion[nTracks]/F");
  tree->Branch("trackLLHProton",&trackLLHProton,"trackLLHProton[nTracks]/F");
  tree->Branch("trackLLHMuon",&trackLLHMuon,"trackLLHMuon[nTracks]/F");
  tree->Branch("trackLLHKaon",&trackLLHKaon,"trackLLHKaon[nTracks]/F");
  tree->Branch("trackPIDA",&trackPIDA,"trackPIDA[nTracks]/F");
  tree->Branch("trackStartDistToPrimTrkEnd",&trackStartDistToPrimTrkEnd,"trackStartDistToPrimTrkEnd[nTracks]/F");
  tree->Branch("trackEndDistToPrimTrkEnd",&trackEndDistToPrimTrkEnd,"trackEndDistToPrimTrkEnd[nTracks]/F");
  tree->Branch("trackClosestDistToPrimTrkEnd",&trackClosestDistToPrimTrkEnd,"trackClosestDistToPrimTrkEnd[nTracks]/F");
  tree->Branch("trackStartClosestToPrimTrkEnd",&trackStartClosestToPrimTrkEnd,"trackStartClosestToPrimTrkEnd[nTracks]/O");
  tree->Branch("trackTrueID",&trackTrueID,"trackTrueID[nTracks]/I");
  tree->Branch("trackTrueMotherID",&trackTrueMotherID,"trackTrueMotherID[nTracks]/I");
  tree->Branch("trackTruePdg",&trackTruePdg,"trackTruePdg[nTracks]/I");
  tree->Branch("trackTrueIsBeam",&trackTrueIsBeam,"trackTrueIsBeam[nTracks]/O");
  tree->Branch("trackTrueKin",&trackTrueKin,"trackTrueKin[nTracks]/F");
  tree->Branch("trackTrueEndKin",&trackTrueEndKin,"trackTrueEndKin[nTracks]/F");
  tree->Branch("trackTrueTrajLen",&trackTrueTrajLen,"trackTrueTrajLen[nTracks]/F");
  tree->Branch("trackTrueChargePurity",&trackTrueChargePurity,"trackTrueChargePurity[nTracks]/F");
  tree->Branch("trackTrueChargeEfficiencyU",&trackTrueChargeEfficiencyU,"trackTrueChargeEfficiencyU[nTracks]/F");
  tree->Branch("trackTrueChargeEfficiencyV",&trackTrueChargeEfficiencyV,"trackTrueChargeEfficiencyV[nTracks]/F");
  tree->Branch("trackTrueChargeEfficiencyZ",&trackTrueChargeEfficiencyZ,"trackTrueChargeEfficiencyZ[nTracks]/F");
  tree->Branch("trackTrueStartX",&trackTrueStartX,"trackTrueStartX[nTracks]/F");
  tree->Branch("trackTrueStartY",&trackTrueStartY,"trackTrueStartY[nTracks]/F");
  tree->Branch("trackTrueStartZ",&trackTrueStartZ,"trackTrueStartZ[nTracks]/F");
  tree->Branch("trackTrueStartT",&trackTrueStartT,"trackTrueStartT[nTracks]/F");
  tree->Branch("trackTrueEndX",&trackTrueEndX,"trackTrueEndX[nTracks]/F");
  tree->Branch("trackTrueEndY",&trackTrueEndY,"trackTrueEndY[nTracks]/F");
  tree->Branch("trackTrueEndZ",&trackTrueEndZ,"trackTrueEndZ[nTracks]/F");
  tree->Branch("trackTrueEndT",&trackTrueEndT,"trackTrueEndT[nTracks]/F");
  tree->Branch("trackTrueXFrontTPC",&trackTrueXFrontTPC,"trackTrueXFrontTPC[nTracks]/F");
  tree->Branch("trackTrueYFrontTPC",&trackTrueYFrontTPC,"trackTrueYFrontTPC[nTracks]/F");

  tree->Branch("iBestMatch",&iBestMatch,"iBestMatch/I");
  tree->Branch("trackMatchDeltaX",&trackMatchDeltaX,"trackMatchDeltaX[nTracks]/F");
  tree->Branch("trackMatchDeltaY",&trackMatchDeltaY,"trackMatchDeltaY[nTracks]/F");
  tree->Branch("trackMatchDeltaR",&trackMatchDeltaR,"trackMatchDeltaR[nTracks]/F");
  tree->Branch("trackMatchDeltaAngle",&trackMatchDeltaAngle,"trackMatchDeltaAngle[nTracks]/F");
  tree->Branch("trackMatchLowestZ",&trackMatchLowestZ,"trackMatchLowestZ[nTracks]/F");
  tree->Branch("nMatchedTracks",&nMatchedTracks,"nMatchedTracks/i");

  tree->Branch("primTrkIsMatchPrimary",&primTrkIsMatchPrimary,"primTrkIsMatchPrimary/O");
  tree->Branch("primTrkIsMatchPrimaryDaughter",&primTrkIsMatchPrimaryDaughter,"primTrkIsMatchPrimaryDaughter/O");
  tree->Branch("primTrkIsMatchBeam",&primTrkIsMatchBeam,"primTrkIsMatchBeam/O");
  tree->Branch("primTrkIsMatchAPrimary",&primTrkIsMatchAPrimary,"primTrkIsMatchAPrimary/O");

  tree->Branch("primTrkStartMomTrking",&primTrkStartMomTrking,"primTrkStartMomTrking/F");
  tree->Branch("primTrkStartTheta",&primTrkStartTheta,"primTrkStartTheta/F");
  tree->Branch("primTrkStartPhi",&primTrkStartPhi,"primTrkStartPhi/F");
  tree->Branch("primTrkLength",&primTrkLength,"primTrkLength/F");
  tree->Branch("primTrkStartX",&primTrkStartX,"primTrkStartX/F");
  tree->Branch("primTrkStartY",&primTrkStartY,"primTrkStartY/F");
  tree->Branch("primTrkStartZ",&primTrkStartZ,"primTrkStartZ/F");
  tree->Branch("primTrkEndX",&primTrkEndX,"primTrkEndX/F");
  tree->Branch("primTrkEndY",&primTrkEndY,"primTrkEndY/F");
  tree->Branch("primTrkEndZ",&primTrkEndZ,"primTrkEndZ/F");
  tree->Branch("primTrkXFrontTPC",&primTrkXFrontTPC,"primTrkXFrontTPC/F");
  tree->Branch("primTrkYFrontTPC",&primTrkYFrontTPC,"primTrkYFrontTPC/F");
  tree->Branch("primTrkCaloRange",&primTrkCaloRange,"primTrkCaloRange/F");
  tree->Branch("primTrkEndInFid",&primTrkEndInFid,"primTrkEndInFid/O");

  tree->Branch("primTrkCaloKin",&primTrkCaloKin,"primTrkCaloKin/F");
  tree->Branch("primTrkEndKin",&primTrkEndKin,"primTrkEndKin/F");
  tree->Branch("primTrkEndKinFid",&primTrkEndKinFid,"primTrkEndKinFid/F");
  tree->Branch("primTrkKinInteract",&primTrkKinInteract,"primTrkKinInteract/F");
  tree->Branch("primTrkKinInteractProton",&primTrkKinInteractProton,"primTrkKinInteractProton/F");
  tree->Branch("primTrkNHits",&primTrkNHits,"primTrkNHits/F");
  tree->Branch("primTrkLLHPion",&primTrkLLHPion,"primTrkLLHPion/F");
  tree->Branch("primTrkLLHProton",&primTrkLLHProton,"primTrkLLHProton/F");
  tree->Branch("primTrkLLHMuon",&primTrkLLHMuon,"primTrkLLHMuon/F");
  tree->Branch("primTrkLLHKaon",&primTrkLLHKaon,"primTrkLLHKaon/F");
  tree->Branch("primTrkPIDA",&primTrkPIDA,"primTrkPIDA/F");

  tree->Branch("primTrkResRangesFlipped",&primTrkResRangesFlipped,"primTrkResRangesFlipped/O");
  tree->Branch("primTrkdEdxMedianLast3Hits",&primTrkdEdxMedianLast3Hits,"primTrkdEdxMedianLast3Hits/F");
  tree->Branch("primTrkdEdxAverageLast3Hits",&primTrkdEdxAverageLast3Hits,"primTrkdEdxAverageLast3Hits/F");
  tree->Branch("primTrkdEdxMedianLast5Hits",&primTrkdEdxMedianLast5Hits,"primTrkdEdxMedianLast5Hits/F");
  tree->Branch("primTrkdEdxAverageLast5Hits",&primTrkdEdxAverageLast5Hits,"primTrkdEdxAverageLast5Hits/F");
  tree->Branch("primTrkdEdxMedianLast7Hits",&primTrkdEdxMedianLast7Hits,"primTrkdEdxMedianLast7Hits/F");
  tree->Branch("primTrkdEdxAverageLast7Hits",&primTrkdEdxAverageLast7Hits,"primTrkdEdxAverageLast7Hits/F");
  tree->Branch("primTrkdEdxMedianRRL1",&primTrkdEdxMedianRRL1,"primTrkdEdxMedianRRL1/F");
  tree->Branch("primTrkdEdxAverageRRL1",&primTrkdEdxAverageRRL1,"primTrkdEdxAverageRRL1/F");
  tree->Branch("primTrkdEdxMedianRRL2",&primTrkdEdxMedianRRL2,"primTrkdEdxMedianRRL2/F");
  tree->Branch("primTrkdEdxAverageRRL2",&primTrkdEdxAverageRRL2,"primTrkdEdxAverageRRL2/F");
  tree->Branch("primTrkdEdxMedianRRL3",&primTrkdEdxMedianRRL3,"primTrkdEdxMedianRRL3/F");
  tree->Branch("primTrkdEdxAverageRRL3",&primTrkdEdxAverageRRL3,"primTrkdEdxAverageRRL3/F");
  tree->Branch("primTrkdEdxMedianRRL5",&primTrkdEdxMedianRRL5,"primTrkdEdxMedianRRL5/F");
  tree->Branch("primTrkdEdxAverageRRL5",&primTrkdEdxAverageRRL5,"primTrkdEdxAverageRRL5/F");
  tree->Branch("primTrkdEdxMedianRRL7",&primTrkdEdxMedianRRL7,"primTrkdEdxMedianRRL7/F");
  tree->Branch("primTrkdEdxAverageRRL7",&primTrkdEdxAverageRRL7,"primTrkdEdxAverageRRL7/F");
  tree->Branch("primTrkdEdxMedianRRL5G1",&primTrkdEdxMedianRRL5G1,"primTrkdEdxMedianRRL5G1/F");
  tree->Branch("primTrkdEdxAverageRRL5G1",&primTrkdEdxAverageRRL5G1,"primTrkdEdxAverageRRL5G1/F");
  tree->Branch("primTrkdEdxMedianRRL7G1",&primTrkdEdxMedianRRL7G1,"primTrkdEdxMedianRRL7G1/F");
  tree->Branch("primTrkdEdxAverageRRL7G1",&primTrkdEdxAverageRRL7G1,"primTrkdEdxAverageRRL7G1/F");

  tree->Branch("primTrkdEdxs",&primTrkdEdxs);
  tree->Branch("primTrkResRanges",&primTrkResRanges);
  tree->Branch("primTrkRangeSoFars",&primTrkRangeSoFars);
  tree->Branch("primTrkPitches",&primTrkPitches);
  tree->Branch("primTrkIBackwards",&primTrkIBackwards);
  tree->Branch("primTrkXs",&primTrkXs);
  tree->Branch("primTrkYs",&primTrkYs);
  tree->Branch("primTrkZs",&primTrkZs);
  tree->Branch("primTrkKins",&primTrkKins);
  tree->Branch("primTrkKinsProton",&primTrkKinsProton);
  tree->Branch("primTrkInFids",&primTrkInFids);
  tree->Branch("primTrkKinsTrue",&primTrkKinsTrue);
  tree->Branch("primTrkDistToTrueTraj",&primTrkDistToTrueTraj);
  tree->Branch("primTrkDistToTrueTrajPoint",&primTrkDistToTrueTrajPoint);

  tree->Branch("nSecTrk",&nSecTrk,"nSecTrk/I");
  tree->Branch("iSecTrkID",&iSecTrkID,"iSecTrkID[nTracks]/I");
  tree->Branch("SecTrkPID",&SecTrkPID,"SecTrkPID[nTracks]/O");
  tree->Branch("iSecMin",&iSecMin,"iSecMin/I");
  tree->Branch("nSecLLRProtonToPionMax",&nSecLLRProtonToPionMax,"nSecLLRProtonToPionMax/F");
  tree->Branch("nSecLLRProtonToPionMin",&nSecLLRProtonToPionMin,"nSecLLRProtonToPionMin/F");

  tree->Branch("secMinTrkdEdxs",&secMinTrkdEdxs);
  tree->Branch("secMinTrkResRanges",&secMinTrkResRanges);
  tree->Branch("secMinTrkPitches",&secMinTrkPitches);
  tree->Branch("secMinTrkIBackwards",&secMinTrkIBackwards);
  tree->Branch("secMinTrkXs",&secMinTrkXs);
  tree->Branch("secMinTrkYs",&secMinTrkYs);
  tree->Branch("secMinTrkZs",&secMinTrkZs);
  tree->Branch("secMinTrkInFids",&secMinTrkInFids);

  tree->Branch("nSecTrkLLRG0",&nSecTrkLLRG0,"nSecTrkLLRG0/I");
  tree->Branch("nSecTrkLLRG100",&nSecTrkLLRG100,"nSecTrkLLRG100/I");
  tree->Branch("nSecTrkLLRG200",&nSecTrkLLRG200,"nSecTrkLLRG200/I");
  tree->Branch("nSecTrkLLRG300",&nSecTrkLLRG300,"nSecTrkLLRG300/I");
  tree->Branch("nSecTrkLLRG400",&nSecTrkLLRG400,"nSecTrkLLRG400/I");
  tree->Branch("nSecTrkLLRG500",&nSecTrkLLRG500,"nSecTrkLLRG500/I");
  tree->Branch("nSecTrkLLRG600",&nSecTrkLLRG600,"nSecTrkLLRG600/I");
  tree->Branch("nSecTrkLLRG700",&nSecTrkLLRG700,"nSecTrkLLRG700/I");
  tree->Branch("nSecTrkLLRG800",&nSecTrkLLRG800,"nSecTrkLLRG800/I");
  tree->Branch("nSecTrkLLRG900",&nSecTrkLLRG900,"nSecTrkLLRG900/I");
  tree->Branch("nSecTrkPIDAL8",&nSecTrkPIDAL8,"nSecTrkPIDAL8/I");
  tree->Branch("nSecTrkPIDAL10",&nSecTrkPIDAL10,"nSecTrkPIDAL10/I");
  tree->Branch("nSecTrkPIDAL14",&nSecTrkPIDAL14,"nSecTrkPIDAL14/I");
  tree->Branch("nSecTrkPIDAL16",&nSecTrkPIDAL16,"nSecTrkPIDAL16/I");
  tree->Branch("nSecTrkPIDAL18",&nSecTrkPIDAL18,"nSecTrkPIDAL18/I");

  tree->Branch("PFNBeamSlices",&PFNBeamSlices,"PFNBeamSlices/i");
  tree->Branch("PFBeamPrimNDaughters",&PFBeamPrimNDaughters,"PFBeamPrimNDaughters/i");
  tree->Branch("PFBeamPrimNDaughterTracks",&PFBeamPrimNDaughterTracks,"PFBeamPrimNDaughterTracks/i");
  tree->Branch("PFBeamPrimNDaughterShowers",&PFBeamPrimNDaughterShowers,"PFBeamPrimNDaughterShowers/i");
  tree->Branch("PFBeamPrimPDG",&PFBeamPrimPDG,"PFBeamPrimPDG/I");
  tree->Branch("PFBeamPrimIsTracklike",&PFBeamPrimIsTracklike,"PFBeamPrimIsTracklike/O");
  tree->Branch("PFBeamPrimIsShowerlike",&PFBeamPrimIsShowerlike,"PFBeamPrimIsShowerlike/O");
  tree->Branch("PFBeamPrimBeamCosmicScore",&PFBeamPrimBeamCosmicScore,"PFBeamPrimBeamCosmicScore/F");
  tree->Branch("PFBeamPrimXFrontTPC",&PFBeamPrimXFrontTPC,"PFBeamPrimXFrontTPC/F");
  tree->Branch("PFBeamPrimYFrontTPC",&PFBeamPrimYFrontTPC,"PFBeamPrimYFrontTPC/F");
  tree->Branch("PFBeamPrimXFrontTPCTrackEnd",&PFBeamPrimXFrontTPCTrackEnd,"PFBeamPrimXFrontTPCTrackEnd/F");
  tree->Branch("PFBeamPrimYFrontTPCTrackEnd",&PFBeamPrimYFrontTPCTrackEnd,"PFBeamPrimYFrontTPCTrackEnd/F");
  tree->Branch("PFBeamPrimStartX",&PFBeamPrimStartX,"PFBeamPrimStartX/F");
  tree->Branch("PFBeamPrimStartY",&PFBeamPrimStartY,"PFBeamPrimStartY/F");
  tree->Branch("PFBeamPrimStartZ",&PFBeamPrimStartZ,"PFBeamPrimStartZ/F");
  tree->Branch("PFBeamPrimEndX",&PFBeamPrimEndX,"PFBeamPrimEndX/F");
  tree->Branch("PFBeamPrimEndY",&PFBeamPrimEndY,"PFBeamPrimEndY/F");
  tree->Branch("PFBeamPrimEndZ",&PFBeamPrimEndZ,"PFBeamPrimEndZ/F");
  tree->Branch("PFBeamPrimStartTheta",&PFBeamPrimStartTheta,"PFBeamPrimStartTheta/F");
  tree->Branch("PFBeamPrimStartPhi",&PFBeamPrimStartPhi,"PFBeamPrimStartPhi/F");
  tree->Branch("PFBeamPrimEndTheta",&PFBeamPrimEndTheta,"PFBeamPrimEndTheta/F");
  tree->Branch("PFBeamPrimEndPhi",&PFBeamPrimEndPhi,"PFBeamPrimEndPhi/F");
  tree->Branch("PFBeamPrimAngleToBeamTrk",&PFBeamPrimAngleToBeamTrk,"PFBeamPrimAngleToBeamTrk/F");
  tree->Branch("PFBeamPrimEndAngleToBeamTrk",&PFBeamPrimEndAngleToBeamTrk,"PFBeamPrimEndAngleToBeamTrk/F");
  tree->Branch("PFBeamPrimTrkLen",&PFBeamPrimTrkLen,"PFBeamPrimTrkLen/F");
  tree->Branch("PFBeamPrimTrkMaxKink",&PFBeamPrimTrkMaxKink,"PFBeamPrimTrkMaxKink/F");
  tree->Branch("PFBeamPrimTrkStartEndDirAngle",&PFBeamPrimTrkStartEndDirAngle,"PFBeamPrimTrkStartEndDirAngle/F");
  tree->Branch("PFBeamPrimMaxAngleToSecTrks",&PFBeamPrimMaxAngleToSecTrks,"PFBeamPrimMaxAngleToSecTrks/F");
  tree->Branch("PFBeamPrimShwrLen",&PFBeamPrimShwrLen,"PFBeamPrimShwrLen/F");
  tree->Branch("PFBeamPrimShwrOpenAngle",&PFBeamPrimShwrOpenAngle,"PFBeamPrimShwrOpenAngle/F");
  tree->Branch("PFBeamPrimCSDAMomMu",&PFBeamPrimCSDAMomMu,"PFBeamPrimCSDAMomMu/F");
  tree->Branch("PFBeamPrimCSDAMomProton",&PFBeamPrimCSDAMomProton,"PFBeamPrimCSDAMomProton/F");
  tree->Branch("PFBeamPrimCSDAKinMu",&PFBeamPrimCSDAKinMu,"PFBeamPrimCSDAKinMu/F");
  tree->Branch("PFBeamPrimCSDAKinProton",&PFBeamPrimCSDAKinProton,"PFBeamPrimCSDAKinProton/F");
  tree->Branch("PFBeamPrimdEdxAverageLast3Hits",&PFBeamPrimdEdxAverageLast3Hits,"PFBeamPrimdEdxAverageLast3Hits/F");
  tree->Branch("PFBeamPrimdEdxAverageLast5Hits",&PFBeamPrimdEdxAverageLast5Hits,"PFBeamPrimdEdxAverageLast5Hits/F");
  tree->Branch("PFBeamPrimdEdxAverageLast7Hits",&PFBeamPrimdEdxAverageLast7Hits,"PFBeamPrimdEdxAverageLast7Hits/F");

  tree->Branch("PFBeamPrimResRanges",&PFBeamPrimResRanges);
  tree->Branch("PFBeamPrimdEdxs",&PFBeamPrimdEdxs);
  tree->Branch("PFBeamPrimdQdxs",&PFBeamPrimdQdxs);
  tree->Branch("PFBeamPrimPitches",&PFBeamPrimPitches);
  tree->Branch("PFBeamPrimXs",&PFBeamPrimXs);
  tree->Branch("PFBeamPrimYs",&PFBeamPrimYs);
  tree->Branch("PFBeamPrimZs",&PFBeamPrimZs);
  tree->Branch("PFBeamPrimZWires",&PFBeamPrimZWires);
  tree->Branch("PFBeamPrimKins",&PFBeamPrimKins);
  tree->Branch("PFBeamPrimKinsProton",&PFBeamPrimKinsProton);
  tree->Branch("PFBeamPrimKinInteract",&PFBeamPrimKinInteract,"PFBeamPrimKinInteract/F");
  tree->Branch("PFBeamPrimKinInteractProton",&PFBeamPrimKinInteractProton,"PFBeamPrimKinInteractProton/F");

  tree->Branch("PFBeamPrimTrueTrackID",&PFBeamPrimTrueTrackID,"PFBeamPrimTrueTrackID/I");
  tree->Branch("PFBeamPrimTrueTrackIDIsNeg",&PFBeamPrimTrueTrackIDIsNeg,"PFBeamPrimTrueTrackIDIsNeg/O");
  tree->Branch("PFBeamPrimTrueMotherTrackID",&PFBeamPrimTrueMotherTrackID,"PFBeamPrimTrueMotherTrackID/I");
  tree->Branch("PFBeamPrimTruePDG",&PFBeamPrimTruePDG,"PFBeamPrimTruePDG/I");
  tree->Branch("PFBeamPrimTrueKin",&PFBeamPrimTrueKin,"PFBeamPrimTrueKin/F");
  tree->Branch("PFBeamPrimTrueEndKin",&PFBeamPrimTrueEndKin,"PFBeamPrimTrueEndKin/F");
  tree->Branch("PFBeamPrimTrueTrajLen",&PFBeamPrimTrueTrajLen,"PFBeamPrimTrueTrajLen/i");
  tree->Branch("PFBeamPrimTrueIsBeam",&PFBeamPrimTrueIsBeam,"PFBeamPrimTrueIsBeam/O");
  tree->Branch("PFBeamPrimTrueStartX",&PFBeamPrimTrueStartX,"PFBeamPrimTrueStartX/F");
  tree->Branch("PFBeamPrimTrueStartY",&PFBeamPrimTrueStartY,"PFBeamPrimTrueStartY/F");
  tree->Branch("PFBeamPrimTrueStartZ",&PFBeamPrimTrueStartZ,"PFBeamPrimTrueStartZ/F");
  tree->Branch("PFBeamPrimTrueStartT",&PFBeamPrimTrueStartT,"PFBeamPrimTrueStartT/F");
  tree->Branch("PFBeamPrimTrueEndX",&PFBeamPrimTrueEndX,"PFBeamPrimTrueEndX/F");
  tree->Branch("PFBeamPrimTrueEndY",&PFBeamPrimTrueEndY,"PFBeamPrimTrueEndY/F");
  tree->Branch("PFBeamPrimTrueEndZ",&PFBeamPrimTrueEndZ,"PFBeamPrimTrueEndZ/F");
  tree->Branch("PFBeamPrimTrueEndT",&PFBeamPrimTrueEndT,"PFBeamPrimTrueEndT/F");
  tree->Branch("PFBeamPrimTrueMaxKink",&PFBeamPrimTrueMaxKink,"PFBeamPrimTrueMaxKink/F");


  tree->Branch("PFBeamSecTrkLen",&PFBeamSecTrkLen,"PFBeamSecTrkLen[PFBeamPrimNDaughterTracks]/F");
  tree->Branch("PFBeamSecTrkdEdxAverageLast3Hits",&PFBeamSecTrkdEdxAverageLast3Hits,"PFBeamSecTrkdEdxAverageLast3Hits[PFBeamPrimNDaughterTracks]/F");
  tree->Branch("PFBeamSecTrkdEdxAverageLast5Hits",&PFBeamSecTrkdEdxAverageLast5Hits,"PFBeamSecTrkdEdxAverageLast5Hits[PFBeamPrimNDaughterTracks]/F");
  tree->Branch("PFBeamSecTrkdEdxAverageLast7Hits",&PFBeamSecTrkdEdxAverageLast7Hits,"PFBeamSecTrkdEdxAverageLast7Hits[PFBeamPrimNDaughterTracks]/F");
  tree->Branch("PFBeamSecTrkDistToPrimEnd",&PFBeamSecTrkDistToPrimEnd,"PFBeamSecTrkDistToPrimEnd[PFBeamPrimNDaughterTracks]/F");
  tree->Branch("PFBeamSecTrkAngleToPrimEnd",&PFBeamSecTrkAngleToPrimEnd,"PFBeamSecTrkAngleToPrimEnd[PFBeamPrimNDaughterTracks]/F");

  ////////////////////////////////////////
  // Book histograms

  deltaXYTPCBeamlineHist = tfs->make<TH2F>("deltaXYTPCBeamline","",500,-500,500,500,-500,500);
  setHistTitles(deltaXYTPCBeamlineHist,"#Delta x TPC Track - Beamline Track","#Delta y TPC Track - Beamline Track");
  deltaXYTPCBeamlineOnlyInFlangeHist = tfs->make<TH2F>("deltaXYTPCBeamlineOnlyInFlange","",500,-500,500,500,-500,500);
  setHistTitles(deltaXYTPCBeamlineOnlyInFlangeHist,"#Delta x TPC Track - Beamline Track","#Delta y TPC Track - Beamline Track");
  deltaXYTPCBeamlineOnlyInFirst25cmHist = tfs->make<TH2F>("deltaXYTPCBeamlineOnlyInFirst25cm","",500,-500,500,500,-500,500);
  setHistTitles(deltaXYTPCBeamlineOnlyInFirst25cmHist,"#Delta x TPC Track - Beamline Track","#Delta y TPC Track - Beamline Track");
  deltaXYTPCBeamlineOnlyInFlangeInFirst25cmHist = tfs->make<TH2F>("deltaXYTPCBeamlineOnlyInFlangeInFirst25cm","",500,-500,500,500,-500,500);
  setHistTitles(deltaXYTPCBeamlineOnlyInFlangeInFirst25cmHist,"#Delta x TPC Track - Beamline Track","#Delta y TPC Track - Beamline Track");

  deltaAngleTPCBeamlineHist = tfs->make<TH1F>("deltaAngleTPCBeamline","",500,0,180);
  setHistTitles(deltaAngleTPCBeamlineHist,"#Delta #alpha TPC Track - Beamline Track","TPC-Beamline Track Pairs / Bin");
  deltaAngleTPCBeamlineOnlyInFlangeHist = tfs->make<TH1F>("deltaAngleTPCBeamlineOnlyInFlange","",500,0,180);
  setHistTitles(deltaAngleTPCBeamlineOnlyInFlangeHist,"#Delta #alpha TPC Track - Beamline Track","TPC-Beamline Track Pairs / Bin");
  deltaAngleTPCBeamlineOnlyInFirst25cmHist = tfs->make<TH1F>("deltaAngleTPCBeamlineOnlyInFirst25cm","",500,0,180);
  setHistTitles(deltaAngleTPCBeamlineOnlyInFirst25cmHist,"#Delta #alpha TPC Track - Beamline Track","TPC-Beamline Track Pairs / Bin");
  deltaAngleTPCBeamlineOnlyInFlangeInFirst25cmHist = tfs->make<TH1F>("deltaAngleTPCBeamlineOnlyInFlangeInFirst25cm","",500,0,180);
  setHistTitles(deltaAngleTPCBeamlineOnlyInFlangeInFirst25cmHist,"#Delta #alpha TPC Track - Beamline Track","TPC-Beamline Track Pairs / Bin");
  
  art::ServiceHandle<geo::Geometry> geom;
  size_t iZWire = 0;
  for(const auto& tpcid: geom->IterateTPCIDs())
  {
    if (tpcid.TPC != 1 && tpcid.TPC != 5 && tpcid.TPC != 9 ) continue;
    //std::cout << "TPC: "<< tpcid << " is valid: "<<tpcid.isValid << std::endl;
    for(const auto& planeid: geom->IteratePlaneIDs(tpcid))
    {
      const auto& view = geom->View(planeid);
      if(view != geo::kZ) continue;
      //const auto& sigType = geom->SignalType(planeid);
      //const auto& plane = geom->Plane(planeid);
      //std::cout << "  Plane: "<< planeid << " is valid: "<< planeid.isValid <<" sigtype: " << sigType << std::endl;
      for(const auto & wireid: geom->IterateWireIDs(planeid))
      {
          const auto & wire = geom->Wire(wireid);
          const auto & channel = geom->PlaneWireToChannel(wireid);
          const auto & center = wire.GetCenter(); // geo::Point
          std::cout << "    Wire: "<<wireid<<" Channel: "<<channel 
                  << " Z: " << center.Z() 
                  << std::endl; 
          myChannelToWireMap.emplace(channel,std::make_pair(iZWire,center.Z()));
          iZWire++;
          zWireChannel.push_back(channel);
          zWireWireZ.push_back(center.Z());
      } // for wireid
    } // for planeid
  }// for tpcid

  // The size is the number of wires
  zWireTrueNumElectrons.resize(iZWire);
  zWireTrueEnergy.resize(iZWire);
  zWireTrueX.resize(iZWire);
  zWireTrueY.resize(iZWire);
  zWireTrueZ.resize(iZWire);
  zWireTruedZ.resize(iZWire);
  zWireTruedR.resize(iZWire);
  zWireTrueTrajDistance.resize(iZWire);
  zWireTrueTrajPitch.resize(iZWire);
  zWireTrueTrajKin.resize(iZWire);
  zWireTrueTDC.resize(iZWire);
  zWireTruePartKin.resize(iZWire);

  zWireNHits.resize(iZWire);
  zWireResRange.resize(iZWire);
  zWiredEdx.resize(iZWire);
  zWiredQdx.resize(iZWire);
  zWirePitch.resize(iZWire);
  zWireX.resize(iZWire);
  zWireY.resize(iZWire);
  zWireZ.resize(iZWire);
  zWirePartKin.resize(iZWire);
  zWirePartKinProton.resize(iZWire);

  tree->Branch("zWireChannel",&zWireChannel);
  tree->Branch("zWireWireZ",&zWireWireZ);

  tree->Branch("zWireNHits",&zWireNHits);
  tree->Branch("zWireResRange",&zWireResRange);
  tree->Branch("zWiredEdx",&zWiredEdx);
  tree->Branch("zWiredQdx",&zWiredQdx);
  tree->Branch("zWirePitch",&zWirePitch);
  tree->Branch("zWireX",&zWireX);
  tree->Branch("zWireY",&zWireY);
  tree->Branch("zWireZ",&zWireZ);
  tree->Branch("zWirePartKin",&zWirePartKin);
  tree->Branch("zWirePartKinProton",&zWirePartKinProton);
  tree->Branch("zWireEnergySum",&zWireEnergySum,"zWireEnergySum/F");

  tree->Branch("zWireTrueNumElectrons",&zWireTrueNumElectrons);
  tree->Branch("zWireTrueEnergy",&zWireTrueEnergy);
  tree->Branch("zWireTrueX",&zWireTrueX);
  tree->Branch("zWireTrueY",&zWireTrueY);
  tree->Branch("zWireTrueZ",&zWireTrueZ);
  tree->Branch("zWireTruedZ",&zWireTruedZ);
  tree->Branch("zWireTruedR",&zWireTruedR);
  tree->Branch("zWireTrueTrajDistance",&zWireTrueTrajDistance);
  tree->Branch("zWireTrueTrajPitch",&zWireTrueTrajPitch);
  tree->Branch("zWireTrueTrajKin",&zWireTrueTrajKin);
  tree->Branch("zWireTrueTDC",&zWireTrueTDC);
  tree->Branch("zWireTruePartKin",&zWireTruePartKin);
  tree->Branch("zWireTrueEnergySum",&zWireTrueEnergySum,"zWireTrueEnergySum/F");
}

void lana::PionAbsSelector::beginRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void lana::PionAbsSelector::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void lana::PionAbsSelector::endJob()
{
  // Implementation of optional member function here.
}

void lana::PionAbsSelector::endRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void lana::PionAbsSelector::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void lana::PionAbsSelector::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fDoPMAAnalysis = p.get<bool>("DoPMAAnalysis");
  fDoPandoraAnalysis = p.get<bool>("DoPandoraAnalysis");
  fTruePartLabel = p.get<art::InputTag>("TruePartLabel");
  fBeamTruthTag = p.get<art::InputTag>("BeamTruthTag");
  fCosmicTruthTag = p.get<art::InputTag>("CosmicTruthTag");
  fSimChanLabel = p.get<art::InputTag>("SimChanLabel");
  fTrackLabel = p.get<art::InputTag>("TrackLabel");
  fCaloLabel = p.get<art::InputTag>("CaloLabel");
  fBeamEventTag = p.get<art::InputTag>("BeamEventTag");
  fBeamEventTagOld = p.get<art::InputTag>("BeamEventTagOld");
  fRawTriggerTag = p.get<art::InputTag>("RawTriggerTag");
  //fVertexLabel = p.get<art::InputTag>("VertexLabel");
  fLikelihoodPIDTag = p.get<art::InputTag>("LikelihoodPIDTag");
  fPIDATag = p.get<art::InputTag>("PIDATag");
  fTrackHitTag = p.get<art::InputTag>("TrackHitTag");
  fAllHitTag = p.get<art::InputTag>("AllHitTag");
  fPFParticleTag = p.get<art::InputTag>("PFParticleTag");
  fPFTrackTag = p.get<art::InputTag>("PFTrackTag");
  fPFShowerTag = p.get<art::InputTag>("PFShowerTag");
  fPFCaloTag = p.get<art::InputTag>("PFCaloTag");
  fPandoraHitTag = p.get<art::InputTag>("PandoraHitTag");
  
  fCaloPlane = p.get<unsigned int>("CaloPlane");

  fFiducialPrimaryXMin = p.get<double>("FiducialPrimaryXMin");
  fFiducialPrimaryXMax = p.get<double>("FiducialPrimaryXMax");
  fFiducialPrimaryYMin = p.get<double>("FiducialPrimaryYMin");
  fFiducialPrimaryYMax = p.get<double>("FiducialPrimaryYMax");
  fFiducialPrimaryZMin = p.get<double>("FiducialPrimaryZMin");
  fFiducialPrimaryZMax = p.get<double>("FiducialPrimaryZMax");

  fTrackMatchZMax = p.get<double>("TrackMatchZMax");
  fTrackMatchDeltaXMinData = p.get<double>("TrackMatchDeltaXMinData");
  fTrackMatchDeltaXMaxData = p.get<double>("TrackMatchDeltaXMaxData");
  fTrackMatchDeltaXMinMC = p.get<double>("TrackMatchDeltaXMinMC");
  fTrackMatchDeltaXMaxMC = p.get<double>("TrackMatchDeltaXMaxMC");
  fTrackMatchDeltaYMinData = p.get<double>("TrackMatchDeltaYMinData");
  fTrackMatchDeltaYMaxData = p.get<double>("TrackMatchDeltaYMaxData");
  fTrackMatchDeltaYMinMC = p.get<double>("TrackMatchDeltaYMinMC");
  fTrackMatchDeltaYMaxMC = p.get<double>("TrackMatchDeltaYMaxMC");
  fTrackMatchAngleMaxDeg = p.get<double>("TrackMatchAngleMaxDeg");
  fFlangeCenterX = p.get<double>("FlangeCenterX");
  fFlangeCenterY = p.get<double>("FlangeCenterY");
  fFlangeCenterZ = p.get<double>("FlangeCenterZ");
  fFlangeRadiusCut = p.get<double>("FlangeRadiusCut");
}

void lana::PionAbsSelector::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
  infilename = "INVALID";
}

void lana::PionAbsSelector::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lana::PionAbsSelector::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
  infilename = fb.fileName();
}

void lana::PionAbsSelector::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

const art::Ptr<recob::Track> lana::PionAbsSelector::MatchRecoToTruthOrWCTrack(const std::vector<art::Ptr<recob::Track>>& tracks, bool isData)
{
  if (tracks.size() == 0)
  {
    return art::Ptr<recob::Track>();
  }
  // Get true/wctrack vars
  float xTrue;
  float yTrue;
  float phiTrue;
  float thetaTrue;
  if(isData)
  {
    if(nBeamTracks == 0) return art::Ptr<recob::Track>();
    xTrue = xWC;
    yTrue = yWC;
    phiTrue = phiWC;
    thetaTrue = thetaWC;
  }
  else
  {
    if(trueXFrontTPC < -90000000.) return art::Ptr<recob::Track>();
    xTrue = trueXFrontTPC;
    yTrue = trueYFrontTPC;
    phiTrue = trueStartPhi;
    thetaTrue = trueStartTheta;
  }
  TVector3 dirTrue(1,0,0);
  dirTrue.SetPhi(phiTrue);
  dirTrue.SetTheta(thetaTrue);
  // Now get differences w/ tracks
  //size_t iBestMatchAngle = DEFAULTNEG;
  //size_t iBestMatchR = DEFAULTNEG;
  float deltaRBest = DEFAULTPOS;
  //float deltaAngleBest = DEFAULTPOS;
  assert(nTracks == tracks.size());
  for(size_t iTrack=0; iTrack < nTracks; iTrack++)
  {
    const auto& track = tracks[iTrack];
    float zTrackStart = 1e10;
    if (track->NumberTrajectoryPoints() > 0)
    {
      zTrackStart = track->Vertex().Z();
      float zTrackEnd = track->End().Z();
      if (zTrackStart > zTrackEnd) 
      {
        zTrackStart = zTrackEnd;
      }
    }
    const TVector3 trackIntersectionPoint = lsu::trackZPlane(ZSTARTOFTPC,*track);
    float xTrack = trackXFrontTPC[iTrack];
    float yTrack = trackYFrontTPC[iTrack];
    const float deltaX = xTrack - xTrue;
    const float deltaY = yTrack - yTrue;
    const float deltaR = sqrt(deltaX*deltaX + deltaY*deltaY);
    const auto& trackStartDir = track->DirectionAtPoint(0);
    const float deltaAngle = ROOT::Math::VectorUtil::Angle(trackStartDir,dirTrue);
    if(iTrack <= MAXTRACKS)
    {
      trackMatchDeltaX[iTrack] = deltaX;
      trackMatchDeltaY[iTrack] = deltaY;
      trackMatchDeltaR[iTrack] = deltaR;
      trackMatchDeltaAngle[iTrack] = deltaAngle;
      trackMatchLowestZ[iTrack] = zTrackStart;
    }
    else
    {
      mf::LogError("TooManyTracks") << "N TPC Tracks > MAXTRACKS, not all matches being considered";
    }

    if (zTrackStart > fTrackMatchZMax) continue; // Only allow these tracks for matching
    if (sqrt(pow(xTrue-fFlangeCenterX,2)+pow(yTrue-fFlangeCenterY,2)) > fFlangeRadiusCut) continue;
    if (sqrt(pow(xTrack-fFlangeCenterX,2)+pow(yTrack-fFlangeCenterY,2)) > fFlangeRadiusCut) continue;

    if (isData && (deltaY < fTrackMatchDeltaYMinData || deltaY > fTrackMatchDeltaYMaxData)) continue; // Only allow these tracks for matching
    if (!isData && (deltaY < fTrackMatchDeltaYMinMC || deltaY > fTrackMatchDeltaYMaxMC)) continue; // Only allow these tracks for matching
    if (isData && (deltaX < fTrackMatchDeltaXMinData || deltaX > fTrackMatchDeltaXMaxData)) continue; // Only allow these tracks for matching
    if (!isData && (deltaX < fTrackMatchDeltaXMinMC || deltaX > fTrackMatchDeltaXMaxMC)) continue; // Only allow these tracks for matching
    if (deltaAngle > fTrackMatchAngleMaxDeg) continue; // Only allow these tracks for matching

    nMatchedTracks++;
    if(deltaR < deltaRBest)
    {
      iBestMatch = iTrack;
      //iBestMatchR = iTrack;
      deltaRBest = deltaR;
      //deltaAngleBest = deltaAngle;
    }
    //if(deltaAngle < deltaAngleBest)
    //{
    //  //iBestMatchAngle = iTrack;
    //  deltaRBest = deltaR;
    //  deltaAngleBest = deltaAngle;
    //}
  } // for iTrack
//  std::cout << "At end of TrackMatch: " << tracks.size() << " " << iBestMatch << std::endl;

  if(iBestMatch<0)
  {
    return art::Ptr<recob::Track>();
  }
  else
  {
    return tracks[iBestMatch];
  }

} // MatchRecoToTruthOrWCTrack

void lana::PionAbsSelector::ResetTreeVars() 
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
  xWC1Hit = DEFAULTNEG;
  yWC1Hit = DEFAULTNEG;
  zWC1Hit = DEFAULTNEG;
  xWC2Hit = DEFAULTNEG;
  yWC2Hit = DEFAULTNEG;
  zWC2Hit = DEFAULTNEG;

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
  for(size_t iTOF=0; iTOF < MAXTOFS; iTOF++)
  {
    TOFs[iTOF] = DEFAULTNEG;
    TOFChans[iTOF] = DEFAULTNEG;
    TOFusTrigs[iTOF] = 0;
    TOFdsTrigs[iTOF] = 0;
  }
  for(size_t iTOF=0; iTOF < NTOFCHANS; iTOF++)
  {
    TOFsByChan[iTOF] = DEFAULTNEG;
    TOFusTrigsByChan[iTOF] = 0;
    TOFdsTrigsByChan[iTOF] = 0;
  }

  BIElectron0p3GeV = false;
  BIMuon0p3GeV = false;
  BIPion0p3GeV = false;
  BIKaon0p3GeV = false;
  BIProton0p3GeV = false;
  BIDeuteron0p3GeV = false;
  BIElectron0p5GeV = false;
  BIMuon0p5GeV = false;
  BIPion0p5GeV = false;
  BIKaon0p5GeV = false;
  BIProton0p5GeV = false;
  BIDeuteron0p5GeV = false;
  BIElectron1GeV = false;
  BIMuon1GeV = false;
  BIPion1GeV = false;
  BIKaon1GeV = false;
  BIProton1GeV = false;
  BIDeuteron1GeV = false;
  BIElectron2GeV = false;
  BIMuon2GeV = false;
  BIPion2GeV = false;
  BIKaon2GeV = false;
  BIProton2GeV = false;
  BIDeuteron2GeV = false;
  BIElectron3GeV = false;
  BIMuon3GeV = false;
  BIPion3GeV = false;
  BIKaon3GeV = false;
  BIProton3GeV = false;
  BIDeuteron3GeV = false;
  BIElectron6GeV = false;
  BIMuon6GeV = false;
  BIPion6GeV = false;
  BIKaon6GeV = false;
  BIProton6GeV = false;
  BIDeuteron6GeV = false;
  BIElectron7GeV = false;
  BIMuon7GeV = false;
  BIPion7GeV = false;
  BIKaon7GeV = false;
  BIProton7GeV = false;
  BIDeuteron7GeV = false;

  nBeamTracksOld = 0;
  for(size_t iTrack=0; iTrack < MAXBEAMTRACKS; iTrack++)
  {
    beamTrackXFrontTPCOld[iTrack] = DEFAULTNEG;
    beamTrackYFrontTPCOld[iTrack] = DEFAULTNEG;
    beamTrackThetaOld[iTrack] = DEFAULTNEG;
    beamTrackPhiOld[iTrack] = DEFAULTNEG;
  }

  nBeamMomOld = 0;
  for(size_t iTrack=0; iTrack < MAXBEAMMOMS; iTrack++)
  {
    beamMomOld[iTrack] = DEFAULTNEG;
  }

  nBeamEventsOld = 0;
  BITriggerOld = DEFAULTNEG;
  BITriggerMatchedOld = false;
  BIActiveTriggerOld = 0;
  TOFOld = DEFAULTNEG;
  CKov0StatusOld = DEFAULTNEG;
  CKov1StatusOld = DEFAULTNEG;

  BIElectron0p3GeVOld = false;
  BIMuon0p3GeVOld = false;
  BIPion0p3GeVOld = false;
  BIKaon0p3GeVOld = false;
  BIProton0p3GeVOld = false;
  BIDeuteron0p3GeVOld = false;
  BIElectron0p5GeVOld = false;
  BIMuon0p5GeVOld = false;
  BIPion0p5GeVOld = false;
  BIKaon0p5GeVOld = false;
  BIProton0p5GeVOld = false;
  BIDeuteron0p5GeVOld = false;
  BIElectron1GeVOld = false;
  BIMuon1GeVOld = false;
  BIPion1GeVOld = false;
  BIKaon1GeVOld = false;
  BIProton1GeVOld = false;
  BIDeuteron1GeVOld = false;
  BIElectron2GeVOld = false;
  BIMuon2GeVOld = false;
  BIPion2GeVOld = false;
  BIKaon2GeVOld = false;
  BIProton2GeVOld = false;
  BIDeuteron2GeVOld = false;
  BIElectron3GeVOld = false;
  BIMuon3GeVOld = false;
  BIPion3GeVOld = false;
  BIKaon3GeVOld = false;
  BIProton3GeVOld = false;
  BIDeuteron3GeVOld = false;
  BIElectron6GeVOld = false;
  BIMuon6GeVOld = false;
  BIPion6GeVOld = false;
  BIKaon6GeVOld = false;
  BIProton6GeVOld = false;
  BIDeuteron6GeVOld = false;
  BIElectron7GeVOld = false;
  BIMuon7GeVOld = false;
  BIPion7GeVOld = false;
  BIKaon7GeVOld = false;
  BIProton7GeVOld = false;
  BIDeuteron7GeVOld = false;

  triggerIsBeam = false;
  triggerBits = 0;
  for(size_t k=0; k < 6; k++)
  {
    nGoodFEMBs[k] = DEFAULTNEG;
  }

  CKov0Status = DEFAULTNEG;
  CKov1Status = DEFAULTNEG;
  CKov0Time = DEFAULTNEG;
  CKov1Time = DEFAULTNEG;
  CKov0Pressure = DEFAULTNEG;
  CKov1Pressure = DEFAULTNEG;

  nPrimaryParticleCandidates = 0;

  trueCategory = DEFAULTNEG;
  trueEndProcess = DEFAULTNEG;
  truePrimaryPDG = DEFAULTNEG;
  truePrimaryTrackID = DEFAULTNEG;
  trueSignalT = false;
  trueSignalNT = false;
  trueNDaughters = 0;
  nSecTracks = 0;
  trueNSecondaryChPions = 0;
  trueNSecondaryPiZeros = 0;
  trueNSecondaryProtons = 0;
  trueNSecondaryOppChPions = 0;
  trueNSecondaryMuons = 0;
  trueStartX = DEFAULTNEG;
  trueStartY = DEFAULTNEG;
  trueStartZ = DEFAULTNEG;
  trueStartT = DEFAULTNEG;
  trueEndX = DEFAULTNEG;
  trueEndY = DEFAULTNEG;
  trueEndZ = DEFAULTNEG;
  trueEndT = DEFAULTNEG;
  trueStartTheta = DEFAULTNEG;
  trueStartPhi = DEFAULTNEG;
  trueSecondToEndTheta = DEFAULTNEG;
  trueSecondToEndPhi = DEFAULTNEG;
  trueStartMom = DEFAULTNEG;
  trueStartE = DEFAULTNEG;
  trueStartKin = DEFAULTNEG;
  trueEndMom = DEFAULTNEG;
  trueEndE = DEFAULTNEG;
  trueEndKin = DEFAULTNEG;
  trueEndMaxAngleToSec = DEFAULTNEG;
  trueSecondToEndMom = DEFAULTNEG;
  trueSecondToEndE = DEFAULTNEG;
  trueSecondToEndKin = DEFAULTNEG;
  trueAngleStartDirSecondToEndDir = DEFAULTNEG;
  for(size_t iSec=0; iSec < MAXDAUGHTER; iSec++)
  {
    trueSecondPDG[iSec] = DEFAULTNEG;
    trueSecondKin[iSec] = DEFAULTNEG;
    trueSecondAngleToPrimEnd[iSec] = DEFAULTNEG;
  }
  trueXFrontTPC = DEFAULTNEG;
  trueYFrontTPC = DEFAULTNEG;
  trueTrajPointFrontTPCDist = DEFAULTNEG;
  trueKinFrontTPC = DEFAULTNEG;

  trueTrajNPoints = 0;
  for(size_t i=0; i < MAXTRUETRAJPOINTS; i++)
  {
    trueTrajPointX[i] = DEFAULTNEG;
    trueTrajPointY[i] = DEFAULTNEG;
    trueTrajPointZ[i] = DEFAULTNEG;
    trueTrajPointT[i] = DEFAULTNEG;
    trueTrajPointPx[i] = DEFAULTNEG;
    trueTrajPointPy[i] = DEFAULTNEG;
    trueTrajPointPz[i] = DEFAULTNEG;
    trueTrajPointE[i] = DEFAULTNEG;
    trueTrajPointMom[i] = DEFAULTNEG;
    trueTrajPointKin[i] = DEFAULTNEG;
    trueTrajPointProc[i] = 0;
    trueTrajPointKinkAngle[i] = DEFAULTNEG;
  }

  nMCParts = 0;
  for(size_t i=0; i < MAXMCPARTS; i++)
  {
    mcPartIsBeam[i] = false;
    mcPartTrackID[i] = DEFAULTNEG;
    mcPartIsPrimary[i] = false;
    mcPartPDG[i] = DEFAULTNEG;
    mcPartStartX[i] = DEFAULTNEG;
    mcPartStartY[i] = DEFAULTNEG;
    mcPartStartZ[i] = DEFAULTNEG;
    mcPartStartT[i] = DEFAULTNEG;
    mcPartEndX[i] = DEFAULTNEG;
    mcPartEndY[i] = DEFAULTNEG;
    mcPartEndZ[i] = DEFAULTNEG;
    mcPartEndT[i] = DEFAULTNEG;
    mcPartStartTheta[i] = DEFAULTNEG;
    mcPartStartPhi[i] = DEFAULTNEG;
    mcPartXFrontTPC[i] = DEFAULTNEG;
    mcPartYFrontTPC[i] = DEFAULTNEG;
    mcPartStartMom[i] = DEFAULTNEG;
    mcPartStartE[i] = DEFAULTNEG;
    mcPartStartKin[i] = DEFAULTNEG;
    mcPartEndMom[i] = DEFAULTNEG;
    mcPartEndE[i] = DEFAULTNEG;
    mcPartEndKin[i] = DEFAULTNEG;
    mcPartDeltaAngle[i] = DEFAULTNEG;
  }

  nIDEs = 0;
  for(size_t iIDE=0; iIDE < MAXIDES; iIDE++)
  {
    simIDENumElectrons[iIDE] = DEFAULTNEG;
    simIDEEnergy[iIDE] = DEFAULTNEG;
    simIDEX[iIDE] = DEFAULTNEG;
    simIDEY[iIDE] = DEFAULTNEG;
    simIDEZ[iIDE] = DEFAULTNEG;
    simIDETrajDistance[iIDE] = DEFAULTNEG;
    simIDETrajKin[iIDE] = DEFAULTNEG;
    simIDEWireZ[iIDE] = DEFAULTNEG;
    simIDEChannel[iIDE] = 0;
    simIDETDC[iIDE] = DEFAULTNEG;
    simIDEPartKin[iIDE] = DEFAULTNEG;
  }
  simIDEEnergySum = DEFAULTNEG;

  for(size_t iWire=0; iWire < zWireTrueNumElectrons.size(); iWire++)
  {
    zWireNHits.at(iWire) = 0;
    zWireResRange.at(iWire) = DEFAULTNEG;
    zWiredEdx.at(iWire) = DEFAULTNEG;
    zWiredQdx.at(iWire) = DEFAULTNEG;
    zWirePitch.at(iWire) = DEFAULTNEG;
    zWireX.at(iWire) = DEFAULTNEG;
    zWireY.at(iWire) = DEFAULTNEG;
    zWireZ.at(iWire) = DEFAULTNEG;
    zWirePartKin.at(iWire) = DEFAULTNEG;
    zWirePartKinProton.at(iWire) = DEFAULTNEG;

    zWireTrueNumElectrons.at(iWire) = 0.;
    zWireTrueEnergy.at(iWire) = DEFAULTNEG;
    zWireTrueX.at(iWire) = DEFAULTNEG;
    zWireTrueY.at(iWire) = DEFAULTNEG;
    zWireTrueZ.at(iWire) = DEFAULTNEG;
    zWireTruedZ.at(iWire) = DEFAULTNEG;
    zWireTruedR.at(iWire) = DEFAULTNEG;
    zWireTrueTrajDistance.at(iWire) = DEFAULTNEG;
    zWireTrueTrajPitch.at(iWire) = DEFAULTNEG;
    zWireTrueTrajKin.at(iWire) = DEFAULTNEG;
    zWireTrueTDC.at(iWire) = DEFAULTNEG;
    zWireTruePartKin.at(iWire) = DEFAULTNEG;
  }
  zWireEnergySum = DEFAULTNEG;
  zWireTrueEnergySum = DEFAULTNEG;

  nTracks = 0;
  for(size_t z=0; z < MAXZINT; z++)
  {
    nTracksInFirstZ[z] = 0;
  }
  for(size_t l=0; l < MAXLINT; l++)
  {
    nTracksLengthLt[l] = 0;
  }

  iBestMatch = DEFAULTNEG;
  nMatchedTracks = 0;
  for(size_t iTrack=0; iTrack < MAXTRACKS; iTrack++)
  {
    trackMatchDeltaX[iTrack] = DEFAULTPOS;
    trackMatchDeltaY[iTrack] = DEFAULTPOS;
    trackMatchDeltaR[iTrack] = DEFAULTPOS;
    trackMatchDeltaAngle[iTrack] = DEFAULTPOS;
    trackMatchLowestZ[iTrack] = DEFAULTPOS;

    trackStartX[iTrack] = DEFAULTNEG;
    trackStartY[iTrack] = DEFAULTNEG;
    trackStartZ[iTrack] = DEFAULTNEG;
    trackStartTheta[iTrack] = DEFAULTNEG;
    trackStartPhi[iTrack] = DEFAULTNEG;
    trackEndX[iTrack] = DEFAULTNEG;
    trackEndY[iTrack] = DEFAULTNEG;
    trackEndZ[iTrack] = DEFAULTNEG;
    trackLength[iTrack] = DEFAULTNEG;
    trackXFrontTPC[iTrack] = DEFAULTNEG;
    trackYFrontTPC[iTrack] = DEFAULTNEG;
    trackCaloKin[iTrack] = DEFAULTNEG;
    trackLLHPion[iTrack] = DEFAULTNEG;
    trackLLHProton[iTrack] = DEFAULTNEG;
    trackLLHMuon[iTrack] = DEFAULTNEG;
    trackLLHKaon[iTrack] = DEFAULTNEG;
    trackPIDA[iTrack] = DEFAULTNEG;
    trackStartDistToPrimTrkEnd[iTrack] = DEFAULTPOS;
    trackEndDistToPrimTrkEnd[iTrack] = DEFAULTPOS;
    trackClosestDistToPrimTrkEnd[iTrack] = DEFAULTPOS;
    trackStartClosestToPrimTrkEnd[iTrack] = false;
    trackTrueID[iTrack] = DEFAULTNEG;
    trackTrueMotherID[iTrack] = DEFAULTNEG;
    trackTruePdg[iTrack] = DEFAULTNEG;
    trackTrueIsBeam[iTrack] = false;
    trackTrueKin[iTrack] = DEFAULTNEG;
    trackTrueEndKin[iTrack] = DEFAULTNEG;
    trackTrueTrajLen[iTrack] = DEFAULTNEG;
    trackTrueChargePurity[iTrack] = DEFAULTNEG;
    trackTrueChargeEfficiencyU[iTrack] = DEFAULTNEG;
    trackTrueChargeEfficiencyV[iTrack] = DEFAULTNEG;
    trackTrueChargeEfficiencyZ[iTrack] = DEFAULTNEG;
    trackTrueStartX[iTrack] = DEFAULTNEG;
    trackTrueStartY[iTrack] = DEFAULTNEG;
    trackTrueStartZ[iTrack] = DEFAULTNEG;
    trackTrueStartT[iTrack] = DEFAULTNEG;
    trackTrueEndX[iTrack] = DEFAULTNEG;
    trackTrueEndY[iTrack] = DEFAULTNEG;
    trackTrueEndZ[iTrack] = DEFAULTNEG;
    trackTrueEndT[iTrack] = DEFAULTNEG;
    trackTrueXFrontTPC[iTrack] = DEFAULTNEG;
    trackTrueYFrontTPC[iTrack] = DEFAULTNEG;
    iSecTrkID[iTrack] = DEFAULTNEG;
    SecTrkPID[iTrack] = false;
  }

  primTrkIsMatchPrimary = false;
  primTrkIsMatchPrimaryDaughter = false;
  primTrkIsMatchBeam = false;
  primTrkIsMatchAPrimary = false;

  primTrkStartMomTrking = DEFAULTNEG;
  primTrkStartTheta = DEFAULTNEG;
  primTrkStartPhi = DEFAULTNEG;
  primTrkLength = DEFAULTNEG;
  primTrkStartX = DEFAULTNEG;
  primTrkStartY = DEFAULTNEG;
  primTrkStartZ = DEFAULTNEG;
  primTrkEndX = DEFAULTNEG;
  primTrkEndY = DEFAULTNEG;
  primTrkEndZ = DEFAULTNEG;
  primTrkXFrontTPC = DEFAULTNEG;
  primTrkYFrontTPC = DEFAULTNEG;
  primTrkCaloRange = DEFAULTNEG;
  primTrkEndInFid = false;

  primTrkCaloKin = DEFAULTNEG;
  primTrkEndKin = DEFAULTNEG;
  primTrkEndKinFid = DEFAULTNEG;
  primTrkKinInteract = DEFAULTNEG;
  primTrkKinInteractProton = DEFAULTNEG;
  primTrkNHits = DEFAULTNEG;
  primTrkLLHPion = DEFAULTNEG;
  primTrkLLHProton = DEFAULTNEG;
  primTrkLLHMuon = DEFAULTNEG;
  primTrkLLHKaon = DEFAULTNEG;
  primTrkPIDA = DEFAULTNEG;

  primTrkResRangesFlipped = false;
  primTrkdEdxMedianLast3Hits = DEFAULTNEG;
  primTrkdEdxAverageLast3Hits = DEFAULTNEG;
  primTrkdEdxMedianLast5Hits = DEFAULTNEG;
  primTrkdEdxAverageLast5Hits = DEFAULTNEG;
  primTrkdEdxMedianLast7Hits = DEFAULTNEG;
  primTrkdEdxAverageLast7Hits = DEFAULTNEG;
  primTrkdEdxMedianRRL1 = DEFAULTNEG;
  primTrkdEdxAverageRRL1 = DEFAULTNEG;
  primTrkdEdxMedianRRL2 = DEFAULTNEG;
  primTrkdEdxAverageRRL2 = DEFAULTNEG;
  primTrkdEdxMedianRRL3 = DEFAULTNEG;
  primTrkdEdxAverageRRL3 = DEFAULTNEG;
  primTrkdEdxMedianRRL5 = DEFAULTNEG;
  primTrkdEdxAverageRRL5 = DEFAULTNEG;
  primTrkdEdxMedianRRL7 = DEFAULTNEG;
  primTrkdEdxAverageRRL7 = DEFAULTNEG;
  primTrkdEdxMedianRRL3G1 = DEFAULTNEG;
  primTrkdEdxAverageRRL3G1 = DEFAULTNEG;
  primTrkdEdxMedianRRL5G1 = DEFAULTNEG;
  primTrkdEdxAverageRRL5G1 = DEFAULTNEG;
  primTrkdEdxMedianRRL7G1 = DEFAULTNEG;
  primTrkdEdxAverageRRL7G1 = DEFAULTNEG;

  primTrkdEdxs.clear();
  primTrkResRanges.clear();
  primTrkRangeSoFars.clear();
  primTrkPitches.clear();
  primTrkIBackwards.clear();
  primTrkXs.clear();
  primTrkYs.clear();
  primTrkZs.clear();
  primTrkKins.clear();
  primTrkKinsProton.clear();
  primTrkInFids.clear();
  primTrkKinsTrue.clear();
  primTrkDistToTrueTraj.clear();
  primTrkDistToTrueTrajPoint.clear();

  nSecTrk = 0;
  iSecMin = DEFAULTNEG;
  nSecLLRProtonToPionMax = DEFAULTNEG;
  nSecLLRProtonToPionMin = DEFAULTPOS;
  secMinTrkdEdxs.clear();
  secMinTrkResRanges.clear();
  secMinTrkPitches.clear();
  secMinTrkIBackwards.clear();
  secMinTrkXs.clear();
  secMinTrkYs.clear();
  secMinTrkZs.clear();
  secMinTrkInFids.clear();

  nSecTrkLLRG0 = 0;
  nSecTrkLLRG100 = 0;
  nSecTrkLLRG200 = 0;
  nSecTrkLLRG300 = 0;
  nSecTrkLLRG400 = 0;
  nSecTrkLLRG500 = 0;
  nSecTrkLLRG600 = 0;
  nSecTrkLLRG700 = 0;
  nSecTrkLLRG800 = 0;
  nSecTrkLLRG900 = 0;
  nSecTrkPIDAL8 = 0; 
  nSecTrkPIDAL10 = 0;
  nSecTrkPIDAL14 = 0;
  nSecTrkPIDAL16 = 0;
  nSecTrkPIDAL18 = 0;

  PFNBeamSlices = 0;
  PFBeamPrimNDaughters = 0;
  PFBeamPrimNDaughterTracks = 0;
  PFBeamPrimNDaughterShowers = 0;
  PFBeamPrimPDG = DEFAULTNEG;
  PFBeamPrimIsTracklike = false;
  PFBeamPrimIsShowerlike = false;
  PFBeamPrimBeamCosmicScore = DEFAULTNEG;
  PFBeamPrimXFrontTPC = DEFAULTNEG;
  PFBeamPrimYFrontTPC = DEFAULTNEG;
  PFBeamPrimXFrontTPCTrackEnd = DEFAULTNEG;
  PFBeamPrimYFrontTPCTrackEnd = DEFAULTNEG;
  PFBeamPrimStartX = DEFAULTNEG;
  PFBeamPrimStartY = DEFAULTNEG;
  PFBeamPrimStartZ = DEFAULTNEG;
  PFBeamPrimEndX = DEFAULTNEG;
  PFBeamPrimEndY = DEFAULTNEG;
  PFBeamPrimEndZ = DEFAULTNEG;
  PFBeamPrimStartTheta = DEFAULTNEG;
  PFBeamPrimStartPhi = DEFAULTNEG;
  PFBeamPrimEndTheta = DEFAULTNEG;
  PFBeamPrimEndPhi = DEFAULTNEG;
  PFBeamPrimAngleToBeamTrk = DEFAULTNEG;
  PFBeamPrimEndAngleToBeamTrk = DEFAULTNEG;
  PFBeamPrimTrkLen = DEFAULTNEG;
  PFBeamPrimTrkMaxKink = DEFAULTNEG;
  PFBeamPrimTrkStartEndDirAngle = DEFAULTNEG;
  PFBeamPrimMaxAngleToSecTrks = DEFAULTNEG;
  PFBeamPrimShwrLen = DEFAULTNEG;
  PFBeamPrimShwrOpenAngle = DEFAULTNEG;
  PFBeamPrimCSDAMomMu = DEFAULTNEG;
  PFBeamPrimCSDAMomProton = DEFAULTNEG;
  PFBeamPrimCSDAKinMu = DEFAULTNEG;
  PFBeamPrimCSDAKinProton = DEFAULTNEG;
  PFBeamPrimdEdxAverageLast3Hits = DEFAULTNEG;
  PFBeamPrimdEdxAverageLast5Hits = DEFAULTNEG;
  PFBeamPrimdEdxAverageLast7Hits = DEFAULTNEG;

  PFBeamPrimResRanges.clear();
  PFBeamPrimdEdxs.clear();
  PFBeamPrimdQdxs.clear();
  PFBeamPrimPitches.clear();
  PFBeamPrimXs.clear();
  PFBeamPrimYs.clear();
  PFBeamPrimZs.clear();
  PFBeamPrimZWires.clear();
  PFBeamPrimKins.clear();
  PFBeamPrimKinsProton.clear();
  PFBeamPrimKinInteract = DEFAULTNEG;
  PFBeamPrimKinInteractProton = DEFAULTNEG;

  PFBeamPrimTrueTrackID = DEFAULTNEG;
  PFBeamPrimTrueTrackIDIsNeg = false;
  PFBeamPrimTrueMotherTrackID = DEFAULTNEG;
  PFBeamPrimTruePDG = DEFAULTNEG;
  PFBeamPrimTrueKin = DEFAULTNEG;
  PFBeamPrimTrueEndKin = DEFAULTNEG;
  PFBeamPrimTrueTrajLen = 0;
  PFBeamPrimTrueIsBeam = false;
  PFBeamPrimTrueStartX = DEFAULTNEG;
  PFBeamPrimTrueStartY = DEFAULTNEG;
  PFBeamPrimTrueStartZ = DEFAULTNEG;
  PFBeamPrimTrueStartT = DEFAULTNEG;
  PFBeamPrimTrueEndX = DEFAULTNEG;
  PFBeamPrimTrueEndY = DEFAULTNEG;
  PFBeamPrimTrueEndZ = DEFAULTNEG;
  PFBeamPrimTrueEndT = DEFAULTNEG;
  PFBeamPrimTrueMaxKink = DEFAULTNEG;

  for(size_t iSec=0; iSec < MAXPFSECTRKS; iSec++)
  {
    PFBeamSecTrkLen[iSec] = DEFAULTNEG;
    PFBeamSecTrkdEdxAverageLast3Hits[iSec] = DEFAULTNEG;
    PFBeamSecTrkdEdxAverageLast5Hits[iSec] = DEFAULTNEG;
    PFBeamSecTrkdEdxAverageLast7Hits[iSec] = DEFAULTNEG;
    PFBeamSecTrkDistToPrimEnd[iSec] = DEFAULTNEG;
    PFBeamSecTrkAngleToPrimEnd[iSec] = DEFAULTNEG;
  }

} // ResetTreeVars

const art::Ptr<simb::MCParticle> lana::PionAbsSelector::ProcessMCParticles(const std::vector<art::Ptr<simb::MCParticle>>& truePartVec, 
        const pdana::MCBeamOrCosmicAlg * const beamOrCosmic, 
        const art::Event& e) // returns primary particle candidate
{
  art::ServiceHandle<geo::Geometry> geom;

  //Get MCParticle Variables and find MCC11 primaryParticleCandidates
  std::vector<art::Ptr<simb::MCParticle> > primaryParticleCandidates;
  std::vector<size_t> primaryParticleCandidateIs;
  for(const auto& truth:(truePartVec))
  {
    if (truth->PdgCode() != 2112 && truth->PdgCode() < 1000000000)
    {
      bool isBeam = false;
      if (beamOrCosmic)
      {
        isBeam= beamOrCosmic->isBeam(truth);
      }
      mf::LogInfo("MCParticle") << std::fixed << std::setprecision(1) 
          << "TrackId: "<<truth->TrackId()
          <<" MotherId: "<<truth->Mother()
          <<" PDG "<< truth->PdgCode()
          <<" KE "<<1000*(truth->E()-truth->Mass()) 
          <<" len: " << truth->Trajectory().TotalLength() 
          <<" start: " << truth->Vx() << ", " << truth->Vy() << ", " << truth->Vz()
          <<" end: " << truth->EndX() << ", " << truth->EndY() << ", " << truth->EndZ()
          <<" Process: " << truth->Process()
          <<" EndProcess: " << truth->EndProcess();

      if (nMCParts >= MAXMCPARTS)
      {
        mf::LogError("MCParticle") << "Too many MCParticles in event to record.";
        continue;
      }

      mcPartIsBeam[nMCParts] = isBeam;
      mcPartIsPrimary[nMCParts] = truth->Process() == "primary";
      mcPartTrackID[nMCParts] = truth->TrackId();
      mcPartPDG[nMCParts] = truth->PdgCode();
      mcPartStartX[nMCParts] = truth->Vx();
      mcPartStartY[nMCParts] = truth->Vy();
      mcPartStartZ[nMCParts] = truth->Vz();
      mcPartStartT[nMCParts] = truth->T();
      mcPartEndX[nMCParts] = truth->EndX();
      mcPartEndY[nMCParts] = truth->EndY();
      mcPartEndZ[nMCParts] = truth->EndZ();
      mcPartEndT[nMCParts] = truth->EndT();
      mcPartStartTheta[nMCParts] = truth->Momentum().Theta();
      mcPartStartPhi[nMCParts] = truth->Momentum().Phi();
      mcPartStartMom[nMCParts] = 1000.*truth->Momentum().Vect().Mag();
      mcPartStartE[nMCParts] = 1000.*truth->E();
      mcPartStartKin[nMCParts] = 1000.*(truth->E()-truth->Mass());
      mcPartEndMom[nMCParts] = 1000.*truth->EndMomentum().Vect().Mag();
      mcPartEndE[nMCParts] = 1000.*truth->EndE();
      mcPartEndKin[nMCParts] = 1000.*(truth->EndE()-truth->Mass());
      mcPartDeltaAngle[nMCParts] = truth->Momentum().Vect().Angle(truth->EndMomentum().Vect());

      const TVector3 particleFrontTPCPoint = lsu::mcPartStartZPlane(ZSTARTOFTPC,*truth);
      mcPartXFrontTPC[nMCParts] = particleFrontTPCPoint.X();
      mcPartYFrontTPC[nMCParts] = particleFrontTPCPoint.Y();

      // for MCC11 primary particle should have start T = 0
      if(mcPartIsBeam[nMCParts] && mcPartIsPrimary[nMCParts]
            && fabs(mcPartStartT[nMCParts]) < 1e-6
        )
      {
        primaryParticleCandidates.push_back(truth);
        primaryParticleCandidateIs.push_back(nMCParts);
      }

      nMCParts++;
    } // if not neturon or nucleus
  } // for true (mcPart)

  // Debugging printing of primaryParticleCandidates
  for(size_t iTruth=0; iTruth < primaryParticleCandidateIs.size(); iTruth++)
  {
    const auto & iMCPart = primaryParticleCandidateIs[iTruth];
    const auto & mcPart = primaryParticleCandidates[iTruth];
    const auto &  trueSecondToEndMomVec4 = lsu::toVector4(std::prev(std::prev(mcPart->Trajectory().end()))->second);
    const auto & trueSecondToEndMom = trueSecondToEndMomVec4.R()*1000.; // in MeV/c
    std::cout << "PrimaryParticle Candidate for: "<< eventNumber <<"\n"
              << "  iMCPart:     " << iMCPart << "\n"
              << "  TrackID:     " << mcPartTrackID[iMCPart] << "\n"
              //<< "  IsPrimary:   " << mcPartIsPrimary[iMCPart] << "\n"
              << "  PDG:         " << mcPartPDG[iMCPart] << "\n"
              << "  Process:     " << mcPart->Process() << "\n"
              << "  EndProcess:  " << mcPart->EndProcess() << "\n"
              << "  StartX:      " << mcPartStartX[iMCPart] << "\n"
              << "  StartY:      " << mcPartStartY[iMCPart] << "\n"
              << "  StartZ:      " << mcPartStartZ[iMCPart] << "\n"
              << "  StartT:      " << mcPartStartT[iMCPart] << "\n"
              << "  EndX:        " << mcPartEndX[iMCPart] << "\n"
              << "  EndY:        " << mcPartEndY[iMCPart] << "\n"
              << "  EndZ:        " << mcPartEndZ[iMCPart] << "\n"
              << "  StartMom:    " << mcPartStartMom[iMCPart] << "\n"
              << "  2ndToEndMom: " << trueSecondToEndMom << "\n"
              << "  StartTheta:  " << mcPartStartTheta[iMCPart]*180/CLHEP::pi << " deg\n"
              << "  StartPhi:    " << mcPartStartPhi[iMCPart]*180/CLHEP::pi << " deg\n"
              << "  XFrontTPC:   " << mcPartXFrontTPC[iMCPart] << "\n"
              << "  YFrontTPC:   " << mcPartYFrontTPC[iMCPart] << "\n"
              ;
    if(nBeamTracks >= MAXBEAMTRACKS)
    {
      throw cet::exception("TooManyBeamTracks","Too many beam tracks in this event, ran out of room in array");
    }
    beamTrackXFrontTPC[nBeamTracks] = mcPartXFrontTPC[iMCPart];
    beamTrackYFrontTPC[nBeamTracks] = mcPartYFrontTPC[iMCPart];
    beamTrackTheta[nBeamTracks] = mcPartStartTheta[iMCPart];
    beamTrackPhi[nBeamTracks] = mcPartStartPhi[iMCPart];
    beamTrackMom[nBeamTracks] = mcPartStartMom[iMCPart]*1e-3; // beam in GeV/c
    beamTrackEPion[nBeamTracks] = sqrt(pow(beamTrackMom[nBeamTracks],2)+pow(MCHARGEDPION,2));
    beamTrackKinPion[nBeamTracks] = beamTrackEPion[nBeamTracks] - MCHARGEDPION;
    beamTrackEProton[nBeamTracks] = sqrt(pow(beamTrackMom[nBeamTracks],2)+pow(MPROTON,2));
    beamTrackKinProton[nBeamTracks] = beamTrackEProton[nBeamTracks] - MPROTON;
    beamTrackTruePDG[nBeamTracks] = mcPartPDG[iMCPart];
    nBeamTracks++;
  }

  // Select primaryParticle
  art::Ptr<simb::MCParticle> primaryParticle;
  nPrimaryParticleCandidates = primaryParticleCandidates.size();
  if(nPrimaryParticleCandidates > 0)
  {
    primaryParticle = primaryParticleCandidates.at(0);
  }

  recob::Track::Point_t trueStartPos;
  recob::Track::Point_t trueEndPos;
  recob::Track::Point_t trueSecondToEndPos;
  ROOT::Math::XYZTVector trueStartMomVec4;
  ROOT::Math::XYZTVector trueEndMomVec4;
  ROOT::Math::XYZTVector trueSecondToEndMomVec4;

  if(primaryParticle)
  {
    trueStartPos = lsu::toPoint(primaryParticle->Trajectory().begin()->first.Vect());
    trueEndPos = lsu::toPoint(std::prev(primaryParticle->Trajectory().end())->first.Vect());
    trueSecondToEndPos = lsu::toPoint(std::prev(std::prev(primaryParticle->Trajectory().end()))->first.Vect());

    trueStartMomVec4 = lsu::toVector4(primaryParticle->Trajectory().begin()->second);
    trueEndMomVec4 = lsu::toVector4(std::prev(primaryParticle->Trajectory().end())->second);
    trueSecondToEndMomVec4 = lsu::toVector4(std::prev(std::prev(primaryParticle->Trajectory().end()))->second);

    trueStartX = trueStartPos.X();
    trueStartY = trueStartPos.Y();
    trueStartZ = trueStartPos.Z();
    trueStartT = primaryParticle->T();
    xWC2Hit = trueStartPos.X();
    yWC2Hit = trueStartPos.Y();
    zWC2Hit = trueStartPos.Z();
    trueEndX = trueEndPos.X();
    trueEndY = trueEndPos.Y();
    trueEndZ = trueEndPos.Z();
    trueEndT = primaryParticle->EndT();
    trueStartTheta = trueStartMomVec4.Vect().Theta();
    trueStartPhi = trueStartMomVec4.Vect().Phi();
    trueSecondToEndTheta = trueSecondToEndMomVec4.Vect().Theta();
    trueSecondToEndPhi = trueSecondToEndMomVec4.Vect().Phi();
    trueAngleStartDirSecondToEndDir = ROOT::Math::VectorUtil::Angle(trueStartMomVec4,trueSecondToEndMomVec4);
    thetaWC = trueStartMomVec4.Vect().Theta();
    phiWC = trueStartMomVec4.Vect().Phi();
    trueStartMom = trueStartMomVec4.R()*1000.; //in MeV/c
    trueStartE = trueStartMomVec4.E()*1000.;// in MeV
    trueStartKin = trueStartE - primaryParticle->Mass()*1000.; // in MeV
    trueEndMom = trueEndMomVec4.R()*1000.; // in MeV/c
    trueEndE = trueEndMomVec4.E()*1000.; // in MeV
    trueEndKin = trueEndE - primaryParticle->Mass()*1000.; // in MeV
    trueSecondToEndMom = trueSecondToEndMomVec4.R()*1000.; // in MeV/c
    trueSecondToEndE = trueSecondToEndMomVec4.E()*1000.; // in MeV
    trueSecondToEndKin = trueSecondToEndE - primaryParticle->Mass()*1000.; // in MeV
    std::string processStr = primaryParticle->Process();
    std::string endProcessStr = primaryParticle->EndProcess();
    trueNDaughters = primaryParticle->NumberDaughters();
    truePrimaryPDG = primaryParticle->PdgCode();
    truePrimaryTrackID = primaryParticle->TrackId();
    nSecTracks = 0;
    for(size_t iDaughter=0; iDaughter < trueNDaughters; iDaughter++)
    {
      int daughterTrackID = primaryParticle->Daughter(iDaughter);
      for(const auto& truth:(truePartVec))
      {
        if(truth->TrackId() == daughterTrackID)
        {
          const int& pdgid = truth->PdgCode();
          const unsigned& abspdgid = abs(pdgid);
          const TLorentzVector& secMom4 = TLorentzVector(truth->Trajectory().begin()->second);
          const auto& secMom3 = lsu::toVector(secMom4);
          const float& primSecAngle = ROOT::Math::VectorUtil::Angle(lsu::toVector(trueEndMomVec4),secMom3);
          const float& secKin = (secMom4.E()-truth->Mass())*1000.;
          if(nSecTracks < MAXDAUGHTER && secKin > 3.0) 
          {
            trueSecondPDG[nSecTracks] = pdgid;
            trueSecondKin[nSecTracks] = secKin;
            trueSecondAngleToPrimEnd[nSecTracks] = primSecAngle;
            nSecTracks++;
          }
          if(abspdgid == 211 || abspdgid == 2212 || abspdgid == 13)
          {
            trueEndMaxAngleToSec = std::max(trueEndMaxAngleToSec,primSecAngle);
          }
          if(abspdgid == 211)
          {
            trueNSecondaryChPions++;
          }
          else if(abspdgid == 2212)
          {
            if(truth->Trajectory().TotalLength() > 1.)
            {
                trueNSecondaryProtons++;
            }
          }
          else if(pdgid == 111)
          {
            trueNSecondaryPiZeros++;
          }
          else if(pdgid == (-truePrimaryPDG))
          {
            trueNSecondaryOppChPions++;
          }
          else if(pdgid == 13)
          {
            trueNSecondaryMuons++;
          }
          break;
        }
      }
    }

	if(trueNSecondaryChPions == 0) trueSignalT = true;
	if(trueNSecondaryChPions > 0) trueSignalNT = true;
	if(endProcessStr == "primary") // created via particle gun or something
	   {trueEndProcess = 0;}
	if(endProcessStr == "pi-Inelastic")
	   {trueEndProcess = 1;}
	if(endProcessStr == "neutronInelastic")
	   {trueEndProcess = 2;}
	if(endProcessStr == "hadElastic")
	   {trueEndProcess = 3;}
	if(endProcessStr == "nCapture")
	   {trueEndProcess = 4;}
	if(endProcessStr == "CHIPSNuclearCaptureAtRest")
	   {trueEndProcess = 5;}
	if(endProcessStr == "Decay")
	   {trueEndProcess = 6;}
	if(endProcessStr == "KaonZeroLInelastic")
	   {trueEndProcess = 7;}
	if(endProcessStr == "CoulombScat")
	   {trueEndProcess = 8;}
	if(endProcessStr == "muMinusCaptureAtRest")
	   {trueEndProcess = 9;}
	if(endProcessStr == "protonInelastic")
	   {trueEndProcess = 10;}
	if(endProcessStr == "kaon+Inelastic")
	   {trueEndProcess = 11;}
	if(endProcessStr == "hBertiniCaptureAtRest")
	   {trueEndProcess = 12;}
	if(endProcessStr == "pi+Inelastic")
	   {trueEndProcess = 13;}
	if(endProcessStr == "LArVoxelReadoutScoringProcess") // just ionizes
	   {trueEndProcess = 14;}
	if(endProcessStr == "CoupledTransportation") // exits world
	   {trueEndProcess = 15;}
	if(endProcessStr == "annihil") // positron
	   {trueEndProcess = 16;}
	if(endProcessStr == "FastScintillation") // just ionizes?
	   {trueEndProcess = 17;}
    mf::LogInfo("PrimaryParticle") <<"process: "<<processStr <<", endProcess: "<< endProcessStr
                        <<", trueEndProcess: "<<trueEndProcess
                        <<", nDaughters: "<< trueNDaughters
                        <<", trueNSecondaryChPions: "<< trueNSecondaryChPions
                        <<", trueNSecondaryPiZeros: "<< trueNSecondaryPiZeros
                        <<", trueNSecondaryProtons: "<< trueNSecondaryProtons
                        << std::setprecision(4) 
                        << ", End: "<<trueEndX<<", "<<trueEndY<<", "<<trueEndZ;

    bool interactedBeforeTPC = false;
    bool interactedOutsideTPC = false;
    if (trueEndZ < 0.) interactedBeforeTPC = true;
    if (trueEndZ > 695.
        || trueEndX < -360. || trueEndX > 360.
        || trueEndY < 0. || trueEndY > 608.
        ) interactedOutsideTPC = true;

    // Set trueCategory
    if (abs(truePrimaryPDG) != 211)
    {
      if (abs(truePrimaryPDG) == 11)
      {
        trueCategory = 11; // electron
      }
      else if (abs(truePrimaryPDG) == 2212)
      {
        trueCategory = 12; // proton
      }
      else if (abs(truePrimaryPDG) == 13)
      {
        trueCategory = 13; // muon
      }
      else if (abs(truePrimaryPDG) == 321)
      {
        trueCategory = 14; // kaon
      }
      else
      {
        trueCategory = 15; // other non-pion primary PDG
      }
    }
    else if (trueEndProcess == 15)
    {
      trueCategory = 8; // through-going
    }
    else if (interactedBeforeTPC)
    {
      trueCategory = 7; // interacted before TPC
    }
    else if (interactedOutsideTPC)
    {
      trueCategory = 6; // interacted outside TPC
    }
    else if (trueEndProcess == 13 || trueEndProcess == 1)
    {
      if (trueNSecondaryOppChPions > 0)
      {
        trueCategory = 4; // pion double charge exchange
      }
      else if(trueNSecondaryChPions > 0)
      {
        trueCategory = 1; // pion inelastic scattering
      }
      else
      {
        if (trueNSecondaryPiZeros > 0)
        {
          trueCategory = 3; // pion charge exchange
        }
        else
        {
          trueCategory = 2; // pion absorption
        }
      }
    }
    else if (trueEndProcess == 6)
    {
      trueCategory = 10; // decay-in-flight
    }
    else if (trueEndProcess == 14 || trueEndProcess == 17) // stopping
    {
      if (trueNSecondaryMuons == 1)
      {
        trueCategory = 9; // Pion Decay At Rest
      }
      else
      {
        trueCategory = 16; // Other Stopping
      }
    }
    else
    {
      trueCategory = 0; // unknown
    }

    //std::cout << "trueEndProcess: " << trueEndProcess << " trueCategory: " << trueCategory << std::endl;

    const TVector3 particleFrontTPCPoint = lsu::mcPartStartZPlane(ZSTARTOFTPC,*primaryParticle);
    trueXFrontTPC = particleFrontTPCPoint.X();
    trueYFrontTPC = particleFrontTPCPoint.Y();
    xWC = trueXFrontTPC;
    yWC = trueYFrontTPC;
    pzWC = trueStartMomVec4.Z()*1000.; //in MeV/c
    pWC = trueStartMom; //in MeV/c

    mctrue::TrajectoryInterpExtrapAlg trajInterpAlg;
    size_t iClosestTrajPointToFrontTPC;
    TLorentzVector momAtFrontTPC;
    double trueTrajPointFrontTPCDistDbl;
    //const auto& closestTrajPointToFrontTPC = trajInterpAlg.pointClosestToPlane(primaryParticle->Trajectory(),ZSTARTOFTPC,momAtFrontTPC,iClosestTrajPointToFrontTPC,trueTrajPointFrontTPCDistDbl);
    trajInterpAlg.pointClosestToPlane(primaryParticle->Trajectory(),ZSTARTOFTPC,momAtFrontTPC,iClosestTrajPointToFrontTPC,trueTrajPointFrontTPCDistDbl);
    trueTrajPointFrontTPCDist = trueTrajPointFrontTPCDistDbl;
    trueKinFrontTPC = 1000.*(momAtFrontTPC.E()-momAtFrontTPC.M()); // MeV
    //std::cout << "closest to front of TPC: iClosest: " << iClosestTrajPointToFrontTPC << " distance: " << trueTrajPointFrontTPCDist << " trueKinFrontTPC: " << trueKinFrontTPC << " z: " << closestTrajPointToFrontTPC.Z() << std::endl;
    //for(size_t iTP=0; iTP < primaryParticle->NumberTrajectoryPoints(); iTP++)
    //{
    //  if(iClosestTrajPointToFrontTPC == iTP) std::cout << "This is TPC front point: " << std::endl;
    //  std::cout << "primaryParticle TrajPoint  iTP: "<<iTP<<" z: " << primaryParticle->Vz(iTP) <<"        KE: "<< 1000.*(primaryParticle->E(iTP) - primaryParticle->Mass()) << std::endl;
    //}

    // Store all trajectory info
    const auto& trajectory = primaryParticle->Trajectory();
    const auto& trajProcMap = trajectory.TrajectoryProcesses(); // I think tells you if elastic interactions happened
    trueTrajNPoints = 0;
    for (size_t iTP=0; iTP < trajectory.size(); iTP++)
    {
      if (trajectory.Z(iTP) > 800. || iTP >= MAXTRUETRAJPOINTS) break;
      trueTrajPointX[iTP] = trajectory.X(iTP);
      trueTrajPointY[iTP] = trajectory.Y(iTP);
      trueTrajPointZ[iTP] = trajectory.Z(iTP);
      trueTrajPointT[iTP] = trajectory.T(iTP);
      const auto& trajMom = trajectory.Momentum(iTP);
      trueTrajPointPx[iTP] = trajMom.Px();
      trueTrajPointPy[iTP] = trajMom.Py();
      trueTrajPointPz[iTP] = trajMom.Pz();
      trueTrajPointE[iTP] = trajMom.E();
      trueTrajPointMom[iTP] = trajMom.Vect().Mag();
      trueTrajPointKin[iTP] = trajMom.E()-trajMom.M();
      unsigned char procCode = 0;
      for (const auto& iPtAndCode: trajProcMap)
      {
        if(iPtAndCode.first == iTP)
        {
          procCode = iPtAndCode.second;
          break;
        } // if iPt == iTP
      } // for iPtAndCode
      trueTrajPointProc[iTP] = procCode;
      if(iTP > 0 && iTP < trajectory.size()-1)
      {
        trueTrajPointKinkAngle[iTP] = trajectory.Momentum(iTP-1).Angle(trajMom.Vect());
      }
      //std::cout << "Trajectory Point " << iTP << " Position: " << trajectory.X(iTP) << "  " << trajectory.Y(iTP) <<"  "<< trajectory.Z(iTP) 
      //          << "  KE: " << trajMom.E() - trajMom.M() << " GeV Theta: "<< trajMom.Theta() << " Phi: " << trajMom.Phi() << " Kink: " << trueTrajPointKinkAngle[iTP] 
      //          << " process: "<<(unsigned) procCode  << " \""<< trajectory.KeyToProcess(procCode) << "\"\n";
      trueTrajNPoints++;
    } // for iTP

    // Get SimChannel Info
    protoana::ProtoDUNETruthUtils pdTruthUtils;
    const auto& ideinfos = pdTruthUtils.GetTotalChannelEnergyFromParticle(*primaryParticle,e,true);
    simIDEEnergySum = 0.;
    for (const auto& ideinfo: ideinfos) {
      const auto& [channel, tdc, ide] = ideinfo; // unsigned int, unsigned short, sim::IDE
      //std::cout << "Channel: "<< channel << " tdc: "<< tdc << " energy: " << ide.energy << " MeV  x: "<<ide.x << "  z: "<<ide.z<<"\n";
      if (nIDEs >= MAXIDES - 1)
      {
          throw cet::exception("TooManyIDEs","Too many IDEs in this event, ran out of room in array");
      }

      simIDENumElectrons[nIDEs] = ide.numElectrons;
      simIDEEnergy[nIDEs] = ide.energy; // MeV
      simIDEX[nIDEs] = ide.x; // cm
      simIDEY[nIDEs] = ide.y;
      simIDEZ[nIDEs] = ide.z;
      simIDEChannel[nIDEs] = channel;
      simIDETDC[nIDEs] = tdc;
      simIDETrajDistance[nIDEs] = 
      simIDEPartKin[nIDEs] = trueKinFrontTPC - simIDEEnergySum;
      simIDEEnergySum += ide.energy;

      TLorentzVector interpolatedMom;
      double distToTraj;
      trajInterpAlg.pointOfClosestApproach(primaryParticle->Trajectory(),TVector3(ide.x,ide.y,ide.z),distToTraj,interpolatedMom);
      const float interpPitch = geom->WirePitch(geo::kZ)/cos(interpolatedMom.Theta());
      simIDETrajDistance[nIDEs] = distToTraj;
      simIDETrajKin[nIDEs] = (interpolatedMom.Vect().Mag() - interpolatedMom.M())*1000.;
      if(myChannelToWireMap.count(channel) > 0)
      {
        simIDEWireZ[nIDEs] = myChannelToWireMap[channel].second;
        const auto& iZWire = myChannelToWireMap[channel].first;
        zWireTrueNumElectrons.at(iZWire) = ide.numElectrons;
        zWireTrueEnergy.at(iZWire) = ide.energy;
        zWireTrueX.at(iZWire) = ide.x;
        zWireTrueY.at(iZWire) = ide.y;
        zWireTrueZ.at(iZWire) = ide.z;
        zWireTrueTDC.at(iZWire) = tdc;
        zWireTrueTrajDistance.at(iZWire) = simIDETrajDistance[nIDEs];
        zWireTrueTrajKin.at(iZWire) = simIDETrajKin[nIDEs];
        zWireTrueTrajPitch.at(iZWire) = interpPitch;
      }
      //std::cout << "primaryParticle IDE z: " << ide.z <<"  energy: "<<ide.energy<<"  Esum: "<<simIDEEnergySum<<"  trueKinFront: "<<trueKinFrontTPC<<"  part kin: " <<simIDEPartKin[nIDEs]<<std::endl;
      nIDEs++;
    } // for ide

    // now z-wire kin energy
    zWireTrueEnergySum = 0.;
    long lastWire = -1;
    for (long iZWire=0; iZWire < (long) zWireTrueEnergy.size(); iZWire++)
    {
      if(zWireTrueEnergy.at(iZWire) >= 0.) lastWire = iZWire;
    }
    for (long iZWire=0; iZWire <= lastWire; iZWire++)
    {
      zWireTruePartKin.at(iZWire) = trueKinFrontTPC - zWireTrueEnergySum;
      if(zWireTrueEnergy.at(iZWire) >= 0.)
      {
        zWireTrueEnergySum += zWireTrueEnergy.at(iZWire);
      }
      if(iZWire < ((long) zWireTrueZ.size() -1) && zWireTrueZ.at(iZWire) > -1000. && zWireTrueZ.at(iZWire+1) > -1000.)
      {
        const recob::Track::Point_t thisPoint(zWireTrueX.at(iZWire),zWireTrueY.at(iZWire),zWireTrueZ.at(iZWire));
        const recob::Track::Point_t nextPoint(zWireTrueX.at(iZWire+1),zWireTrueY.at(iZWire+1),zWireTrueZ.at(iZWire+1));
        const auto& dZ = nextPoint.Z()-thisPoint.Z();
        zWireTruedZ.at(iZWire) = dZ;
        auto&& dR = (nextPoint-thisPoint).R();
        if(dZ < 0.) dR *= -1.;
        zWireTruedR.at(iZWire) = dR;
      }
    } // for iZWire

  } // if primaryParticle

  return primaryParticle;

} // ProcessMCParticles

void lana::PionAbsSelector::ProcessAllTracks(const std::vector<art::Ptr<recob::Track>>& trackVec, 
        const std::vector<art::Ptr<simb::MCParticle>>& truePartVec,
        const art::FindManyP<anab::Calorimetry>& tracksCaloVec, 
        const art::FindManyP<recob::Hit>& fmHitsForTracks,
        const pdana::MCBeamOrCosmicAlg * const beamOrCosmic,
        const art::Event& e)
{
  nTracks = trackVec.size();
  if (nTracks > MAXTRACKS)
  {
    mf::LogError("Track") << "Too many tracks in event to record: "
        << nTracks << " > " << MAXTRACKS << ", not saving event to tree.";
    return;
  }
  for(size_t iTrack=0; iTrack < trackVec.size(); iTrack++)
  {
    const art::Ptr<recob::Track> track = trackVec[iTrack];
    const size_t nPoints = track->NumberTrajectoryPoints();

    if (nPoints > 0)
    {
      trackStartX[iTrack] = track->Vertex().X();
      trackStartY[iTrack] = track->Vertex().Y();
      trackStartZ[iTrack] = track->Vertex().Z();
      trackStartTheta[iTrack] = track->VertexDirection().Theta();
      trackStartPhi[iTrack] = track->VertexDirection().Phi();
      if (nPoints > 1)
      {
        trackEndX[iTrack] = track->End().X();
        trackEndY[iTrack] = track->End().Y();
        trackEndZ[iTrack] = track->End().Z();
      }
      else
      {
        trackEndX[iTrack] = track->Vertex().X();
        trackEndY[iTrack] = track->Vertex().Y();
        trackEndZ[iTrack] = track->Vertex().Z();
      }
      trackLength[iTrack] = track->Length();
    }

    const TVector3 trackIntersectionPoint = lsu::trackZPlane(ZSTARTOFTPC,*track);
    trackXFrontTPC[iTrack] = trackIntersectionPoint.X();
    trackYFrontTPC[iTrack] = trackIntersectionPoint.Y();

    double minz=1000.;
    for(size_t iPoint=0; iPoint < nPoints; iPoint++)
    {
      const double pointZ = track->LocationAtPoint(iPoint).Z();
      if(pointZ < minz)
      {
        minz = pointZ;
      }
    } // for iPoint
    for(size_t z=0; z < MAXZINT; z++)
    {
      if(minz < z)
      {
        // the number of tracks with a space point in [0,i] cm where i is the index
        nTracksInFirstZ[z]++;
      }
    } // for z
    const double trkLength = track->Length();
    for(size_t l=0; l < MAXLINT; l++)
    {
      if(trkLength < l)
      {
        // the number of tracks with length less than i cm where i is the index
        nTracksLengthLt[l]++;
      }
    } // for l
    // calo info
    const auto calos = tracksCaloVec.at(iTrack);
    for(const auto& calo:calos)
    {
      if(calo->PlaneID().Plane == fCaloPlane)
      {
          trackCaloKin[iTrack] = calo->KineticEnergy();
          //for(size_t cRangeIt = 0; cRangeIt < calo->ResidualRange().size() && cRangeIt < calo->dEdx().size(); cRangeIt++)
          //{
          //  trackResRanges[iTrk].push_back(primTrkCalo->ResidualRange().at(cRangeIt));
          //  trackdEdxs[iTrk].push_back(primTrkCalo->dEdx().at(cRangeIt));
          //}
        } // if plane == fCaloPlane
    } //for calo in caloVec
    //const auto lhpids = fmlhpid.at(iTrack);
    //for (const auto& lhpid: lhpids)
    //{
    //  if(lhpid->PlaneID().Plane == fCaloPlane)
    //  {
    //    trackLLHPion[iTrack] = lhpid->Chi2Pion();
    //    trackLLHProton[iTrack] = lhpid->Chi2Proton();
    //    trackLLHMuon[iTrack] = lhpid->Chi2Muon();
    //    trackLLHKaon[iTrack] = lhpid->Chi2Kaon();
    //  } // if lhpid plane == fCaloPlane
    //} // for lhpid
    //const auto pidas = fmpida.at(iTrack);
    //for (const auto& pida: pidas)
    //{
    //  if(pida->PlaneID().Plane == fCaloPlane)
    //  {
    //    trackPIDA[iTrack] = pida->PIDA();
    //  } // if pida plane == fCaloPlane
    //} // for pida
    //// Match track to MCParticle
    mf::LogInfo("Track") << std::fixed << std::setprecision(1)
        << "Track: "<< iTrack 
        << " len: " << trkLength 
        <<"  start: " << trackStartX[iTrack] << ", " << trackStartY[iTrack] << ", " << trackStartZ[iTrack]
        <<"  end: " << trackEndX[iTrack] << ", " << trackEndY[iTrack] << ", " << trackEndZ[iTrack]
        ;//<< std::endl;
    if (isMC && fmHitsForTracks.isValid())
    {
	    std::vector< art::Ptr<recob::Hit> > trackHits = fmHitsForTracks.at(iTrack);
	    std::vector< art::Ptr<recob::Hit> > trackHitsU;
	    std::vector< art::Ptr<recob::Hit> > trackHitsV;
	    std::vector< art::Ptr<recob::Hit> > trackHitsZ;
        for (const auto & trackHit : trackHits)
        {
          if (trackHit->View() == geo::kU) trackHitsU.push_back(trackHit);
          else if (trackHit->View() == geo::kV) trackHitsV.push_back(trackHit);
          else if (trackHit->View() == geo::kZ) trackHitsZ.push_back(trackHit);
        }
	    int TrackID = DEFAULTNEG;
        std::set<int> setOfTrackIDs = bt->GetSetOfTrackIds(trackHits);
        std::set<int> setOfPosTrackIDs;
        for (const auto & trackID : setOfTrackIDs)
        {
          //std::cout << "Track: "<<iTrack<<" has TrackID: "<< trackID << std::endl;
          setOfPosTrackIDs.insert(abs(trackID));
        }
        std::set<int> thisTrackIDSet;
        for (const auto & trackID : setOfPosTrackIDs)
        {
          thisTrackIDSet.clear();
          thisTrackIDSet.insert(trackID);
          thisTrackIDSet.insert(-trackID);
          float purity = bt->HitChargeCollectionPurity(thisTrackIDSet,trackHits);
          if (purity > trackTrueChargePurity[iTrack])
          {
              TrackID = trackID;
              trackTrueChargePurity[iTrack] = purity;
          }
        }
        thisTrackIDSet.clear();
        if (TrackID > DEFAULTNEG)
        {
          trackTrueID[iTrack] = TrackID;
          thisTrackIDSet.insert(TrackID);
          thisTrackIDSet.insert(-TrackID);
          // In LArIAT, view V is collection, U is induction.
          // Was taking 90% of the time spent in analyze!!! (valgrind callgrind)
          //trackTrueChargeEfficiencyU[iTrack] = bt->HitChargeCollectionEfficiency(
          //                        thisTrackIDSet,trackHitsU,allHitsVec,geo::kU);
          //trackTrueChargeEfficiencyV[iTrack] = bt->HitChargeCollectionEfficiency(
          //                        thisTrackIDSet,trackHitsV,allHitsVec,geo::kV);
          //trackTrueChargeEfficiencyZ[iTrack] = bt->HitChargeCollectionEfficiency(
          //                        thisTrackIDSet,trackHitsZ,allHitsVec,geo::kZ);

          // Find the MCParticle corresponding to this trackID
          art::Ptr<simb::MCParticle> particle;
          for(auto mcpart: truePartVec)
          {
            if (mcpart->TrackId() == TrackID)
            {
                particle = mcpart;
                break;
            }
          }
          if (particle.isNull())
          {
            std::string message = "Couldn't find MCParticle for Track: ";
            message.append(std::to_string(iTrack));
            message.append(" TrackID: ");
            message.append(std::to_string(TrackID));
            throw cet::exception("MCParticleNotFound",message);
          } // if can't find MCParticle for highest Track ID

          // This is where you can do some analysis of the true particle and compare it to the reco
          trackTrueMotherID[iTrack] = particle->Mother();
          trackTruePdg[iTrack] = particle->PdgCode();
          if (beamOrCosmic)
          {
            trackTrueIsBeam[iTrack] = beamOrCosmic->isBeam(particle);
          }
          trackTrueKin[iTrack] = 1000*(particle->E()-particle->Mass());
          trackTrueEndKin[iTrack] = 1000*(particle->EndE()-particle->Mass());
          trackTrueTrajLen[iTrack] = particle->Trajectory().TotalLength();

          trackTrueStartX[iTrack] = particle->Vx();
          trackTrueStartY[iTrack] = particle->Vy();
          trackTrueStartZ[iTrack] = particle->Vz();
          trackTrueStartT[iTrack] = particle->T();

          trackTrueEndX[iTrack] = particle->EndX();
          trackTrueEndY[iTrack] = particle->EndY();
          trackTrueEndZ[iTrack] = particle->EndZ();
          trackTrueEndT[iTrack] = particle->EndT();

          const TVector3 particleFrontTPCPoint = lsu::mcPartStartZPlane(ZSTARTOFTPC,*particle);
          trackTrueXFrontTPC[iTrack] = particleFrontTPCPoint.X();
          trackTrueYFrontTPC[iTrack] = particleFrontTPCPoint.Y();

          //std::cout << "Truth matching info: " << std::endl;
          //std::cout << "  iTrack:                 " << iTrack << std::endl;
          //std::cout << "  ID:                     " << trackTrueID[iTrack] << std::endl;
          //std::cout << "  MotherID:               " << trackTrueMotherID[iTrack] << std::endl;
          //std::cout << "  PDG:                    " << trackTruePdg[iTrack] << std::endl;
          //std::cout << "  IsBeam:                 " << trackTrueIsBeam[iTrack] << std::endl;
          //std::cout << "  KE:                     " << trackTrueKin[iTrack] << std::endl;
          //std::cout << "  EndKE:                  " << trackTrueEndKin[iTrack] << std::endl;
          //std::cout << "  TrajLen                 " << trackTrueTrajLen[iTrack] << std::endl;
          //std::cout << "  ChargePurity:           " << trackTrueChargePurity[iTrack] << std::endl;
          //std::cout << "  ChargeEfficiencyU:      " << trackTrueChargeEfficiencyU[iTrack] << std::endl;
          //std::cout << "  ChargeEfficiencyV:      " << trackTrueChargeEfficiencyV[iTrack] << std::endl;
          //std::cout << "  ChargeEfficiencyZ:      " << trackTrueChargeEfficiencyZ[iTrack] << std::endl;
          //std::cout << "  HitPurity:              " << bt->HitCollectionPurity(thisTrackIDSet,trackHits) << std::endl;
          //std::cout << "  HitEfficiencyU:         " << bt->HitCollectionEfficiency(thisTrackIDSet,trackHitsU,allHitsVec,geo::kU) << std::endl;
          //std::cout << "  HitEfficiencyV:         " << bt->HitCollectionEfficiency(thisTrackIDSet,trackHitsV,allHitsVec,geo::kV) << std::endl;
          //std::cout << "  HitEfficiencyZ:         " << bt->HitCollectionEfficiency(thisTrackIDSet,trackHitsZ,allHitsVec,geo::kZ) << std::endl;
          //std::cout <<"   Process:                " << particle->Process() << std::endl;
          //std::cout <<"   EndProcess:             " << particle->EndProcess() << std::endl;
          //std::cout <<"   Momentum:               " << particle->Momentum().Vect().Mag() << std::endl;
          //std::cout <<"   EndMomentum:            " << particle->EndMomentum().Vect().Mag() << std::endl;

        } // if TrackID >= 0
        else
        {
            mf::LogError("lana::PiAbsSelector::ProcessAllTracks") << "Run: "<<e.run()<<" Event " << e.id().event() << " Couldn't find TrackID. File: "<<infilename<<" Track hits: "<<trackHits.size()<<" Z plane hits: "<<trackHitsZ.size()<<"\n";
            //std::string message = "Couldn't find TrackID for Track ";
            //message.append(std::to_string(iTrack));
            //throw cet::exception("TrackIDNotFound",message);
        }
    } // if isMC && fmHItsForTracks.isValid

  } // for track

} // ProcessAllTracks()

void lana::PionAbsSelector::ProcessPrimaryTrack(const std::vector<art::Ptr<recob::Track>>& trackVec,
        const art::Ptr<simb::MCParticle>& primaryParticle,
        const art::FindManyP<anab::Calorimetry>& tracksCaloVec, 
        const pdana::MotherDaughterWalkerAlg& motherDaughterWalker,
        const art::Event& e)
{
  //Match the primary track to the WCTrack or MCParticle
  const art::Ptr<recob::Track> primaryTrack = MatchRecoToTruthOrWCTrack(trackVec,e.isRealData()); // also sets deltaX, deltaY, etc.
  //Get Primary Track Variables
  if(iBestMatch < 0)
  {
    mf::LogWarning("PionAbsorption") << "Event " << e.id().event() << " thrown out because a primary track could not be identified.\n";
  }
  else if (primaryTrack.isNull())
  {
    mf::LogError("PionAbsorption") << "Event " << e.id().event() << " thrown out because a primary track pointer is unexpectedly null.\n";
  }
  else if (primaryTrack->NumberTrajectoryPoints() == 0)
  {
    mf::LogWarning("PionAbsorption") << "Event " << e.id().event() << " thrown out because a primary track had no trajectory points.\n";
  }
  else // good primaryTrack
  {
    recob::Track::Point_t primTrkStart;
    recob::Track::Point_t primTrkEnd;
    recob::Track::Vector_t primTrkStartDir;
    recob::Track::Vector_t primTrkEndDir;
    primTrkStart = primaryTrack->LocationAtPoint(0);
    primTrkEnd = primaryTrack->LocationAtPoint(primaryTrack->NumberTrajectoryPoints()-1);
    primTrkStartDir = primaryTrack->DirectionAtPoint(0);
    primTrkEndDir = primaryTrack->DirectionAtPoint(primaryTrack->NumberTrajectoryPoints()-1);
    primTrkStartMomTrking = primaryTrack->MomentumAtPoint(0)*1000.; // MeV/c
    primTrkStartTheta = primTrkStartDir.Theta();
    primTrkStartPhi = primTrkStartDir.Phi();
    primTrkLength = primaryTrack->Length();
    primTrkStartX = primTrkStart.X();
    primTrkStartY = primTrkStart.Y();
    primTrkStartZ = primTrkStart.Z();
    primTrkEndX = primTrkEnd.X();
    primTrkEndY = primTrkEnd.Y();
    primTrkEndZ = primTrkEnd.Z();
    primTrkXFrontTPC = trackTrueXFrontTPC[iBestMatch];
    primTrkYFrontTPC = trackTrueYFrontTPC[iBestMatch];
    primTrkEndInFid = InPrimaryFiducial(primTrkEnd);

    primTrkCaloKin = trackCaloKin[iBestMatch];
    primTrkEndKin = kinWCInTPC - primTrkCaloKin;
    if(primTrkEndKin < 0.) primTrkEndKin = 0;
    if(primTrkEndInFid) primTrkEndKinFid = primTrkEndKin;
    primTrkLLHPion = trackLLHPion[iBestMatch];
    primTrkLLHProton = trackLLHProton[iBestMatch];
    primTrkLLHMuon = trackLLHMuon[iBestMatch];
    primTrkLLHKaon = trackLLHKaon[iBestMatch];
    primTrkPIDA = trackPIDA[iBestMatch];

    mf::LogInfo("PrimaryTrack") << "Primary Track start point: (" << std::fixed << std::setprecision(1) << primTrkStart.X() << "," << primTrkStart.Y() << "," << primTrkStart.Z() << ")";
    mf::LogInfo("PrimaryTrack") << "Primary Track end point: (" << std::fixed << std::setprecision(1) << primTrkEnd.X() << "," << primTrkEnd.Y() << "," << primTrkEnd.Z() << ")";

    // Primary Track Truth Match Info
    if(isMC)
    {
      primTrkIsMatchPrimary = trackTrueID[iBestMatch] == truePrimaryTrackID;
      if (!primTrkIsMatchPrimary)
      {
        const auto thisMatchedParticle = motherDaughterWalker.getParticle(trackTrueID[iBestMatch]);
        if(thisMatchedParticle.isNull())
        {
            throw cet::exception("CantFindMCParticle","Couldn't find MCParticle for best match TrackID");
        }
        const auto greatestGrandmother = motherDaughterWalker.getGreatestGrandmother(*thisMatchedParticle);
        
        if (greatestGrandmother.isNonnull())
        {
          primTrkIsMatchPrimaryDaughter = greatestGrandmother->TrackId() == truePrimaryTrackID;
        }
      }
      primTrkIsMatchBeam = trackTrueIsBeam[iBestMatch];
      primTrkIsMatchAPrimary = trackTrueMotherID[iBestMatch] == 0;
    }

    // Primary Track Calorimetry
    const auto primTrkCalos = tracksCaloVec.at(iBestMatch);
    for(const auto& primTrkCalo:primTrkCalos)
    {
      if(primTrkCalo->PlaneID().Plane == fCaloPlane)
      {
          primTrkCaloRange = primTrkCalo->Range();
          double intE = 0.;
          size_t IBackwards = primTrkCalo->dEdx().size()-1;
          mctrue::TrajectoryInterpExtrapAlg trajInterpAlg;
          for(size_t cRangeIt = 0; cRangeIt < primTrkCalo->ResidualRange().size() && cRangeIt < primTrkCalo->dEdx().size(); cRangeIt++)
          {
            primTrkResRanges.push_back(primTrkCalo->ResidualRange().at(cRangeIt));
            primTrkRangeSoFars.push_back(primTrkLength-primTrkCalo->ResidualRange().at(cRangeIt));
            primTrkdEdxs.push_back(primTrkCalo->dEdx().at(cRangeIt));
            primTrkPitches.push_back(primTrkCalo->TrkPitchVec().at(cRangeIt));
            primTrkIBackwards.push_back(IBackwards);
            IBackwards--;
            const auto thisPoint = primTrkCalo->XYZ().at(cRangeIt); // PositionVector3D
            primTrkXs.push_back(thisPoint.X());
            primTrkYs.push_back(thisPoint.Y());
            primTrkZs.push_back(thisPoint.Z());

            bool thisInFid = InPrimaryFiducial(thisPoint);
            primTrkInFids.push_back(thisInFid);
            double thisKin = kinWCInTPC-intE;
            primTrkKins.push_back(thisKin);
            double thisKinProton = kinWCInTPCProton-intE;
            primTrkKinsProton.push_back(thisKinProton);
            if(thisInFid)
            {
              primTrkKinInteract = thisKin;
              primTrkKinInteractProton = thisKinProton;
            }
            else
            {
              primTrkKinInteract = DEFAULTNEG;
              primTrkKinInteractProton = DEFAULTNEG;
            }
            intE += primTrkCalo->dEdx().at(cRangeIt)*primTrkCalo->TrkPitchVec().at(cRangeIt);

            // true trajectory matching
            if(primaryParticle)
            {
              double distanceToTraj;
              TLorentzVector trueMomVec;
              size_t iClosestTrajPoint;
              double distanceToClosestTrajPoint;
              const TVector3 truePosVec = trajInterpAlg.pointOfClosestApproach(primaryParticle->Trajectory(),thisPoint,distanceToTraj,trueMomVec,iClosestTrajPoint,distanceToClosestTrajPoint);
              double trueKin = trueMomVec.E() - trueMomVec.M();
              trueKin *= 1000.;
              primTrkKinsTrue.push_back(trueKin);
              primTrkDistToTrueTraj.push_back(distanceToTraj);
              primTrkDistToTrueTrajPoint.push_back(distanceToClosestTrajPoint);
            } // if primaryParticle
          } // for cRangeIt
          primTrkResRangesFlipped = (primTrkResRanges.front() - primTrkResRanges.back()) < 0.;
          // look into the end of the track
          const size_t dEdxSize = primTrkdEdxs.size();
          const std::vector<float>::const_iterator dEdxBegin = primTrkdEdxs.begin();
          const std::vector<float>::const_iterator dEdxEnd = primTrkdEdxs.end();
          if (dEdxSize >= 3)
          {
            primTrkdEdxMedianLast3Hits = findQuantile(0.5,dEdxEnd - 3, dEdxEnd);
            primTrkdEdxAverageLast3Hits = findAverage(dEdxEnd - 3, dEdxEnd);
          }
          if (dEdxSize >= 5)
          {
            primTrkdEdxMedianLast5Hits = findQuantile(0.5,dEdxEnd - 5, dEdxEnd);
            primTrkdEdxAverageLast5Hits = findAverage(dEdxEnd - 5, dEdxEnd);
          }
          if (dEdxSize >= 7)
          {
            primTrkdEdxMedianLast7Hits = findQuantile(0.5,dEdxEnd - 7, dEdxEnd);
            primTrkdEdxAverageLast7Hits = findAverage(dEdxEnd - 7, dEdxEnd);
          }
          const std::vector<float>::const_iterator resRangesBegin = primTrkResRanges.begin();
          const std::vector<float>::const_iterator resRangesEnd = primTrkResRanges.end();
          if(!primTrkResRangesFlipped)
          {
            primTrkdEdxMedianRRL1 = findQuantileResRangeFunc(0.5,resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 1;});
            primTrkdEdxAverageRRL1 = findAverageResRangeFunc(resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 1;});
            primTrkdEdxMedianRRL2 = findQuantileResRangeFunc(0.5,resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 2;});
            primTrkdEdxAverageRRL2 = findAverageResRangeFunc(resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 2;});
            primTrkdEdxMedianRRL3 = findQuantileResRangeFunc(0.5,resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 3;});
            primTrkdEdxAverageRRL3 = findAverageResRangeFunc(resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 3;});
            primTrkdEdxMedianRRL5 = findQuantileResRangeFunc(0.5,resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 5;});
            primTrkdEdxAverageRRL5 = findAverageResRangeFunc(resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 5;});
            primTrkdEdxMedianRRL7 = findQuantileResRangeFunc(0.5,resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 7;});
            primTrkdEdxAverageRRL7 = findAverageResRangeFunc(resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x < 7;});
            primTrkdEdxMedianRRL3G1 = findQuantileResRangeFunc(0.5,resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x>1. && x < 3;});
            primTrkdEdxAverageRRL3G1 = findAverageResRangeFunc(resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x>1. && x < 3;});
            primTrkdEdxMedianRRL5G1 = findQuantileResRangeFunc(0.5,resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x>1. && x < 5;});
            primTrkdEdxAverageRRL5G1 = findAverageResRangeFunc(resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x>1. && x < 5;});
            primTrkdEdxMedianRRL7G1 = findQuantileResRangeFunc(0.5,resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x>1. && x < 7;});
            primTrkdEdxAverageRRL7G1 = findAverageResRangeFunc(resRangesBegin,resRangesEnd,dEdxBegin,dEdxEnd,[](auto x){return x>1. && x < 7;});
          }
        } // if plane == fCaloPlane
    } //for calo in caloVec
    primTrkNHits = primTrkdEdxs.size();


    // Find secondary tracks
    for(size_t iTrack=0; iTrack < nTracks; iTrack++)
    { 
      const auto track = trackVec[iTrack];
      if(track != primaryTrack)
      {
        //mf::LogInfo("SecondaryTrack") << "Got a secondary track to analyze that is not the primary in event " << e.id().event() << ".";
        const auto secTrkStart = track->LocationAtPoint(0);
        const auto secTrkEnd = track->LocationAtPoint(track->NumberTrajectoryPoints()-1);
        mf::LogInfo("SecondaryTrack") << std::fixed << std::setprecision(1) 
                << "Start: (" << secTrkStart.X() << "," << secTrkStart.Y() << "," << secTrkStart.Z() << ") "
                << "End: (" << secTrkEnd.X() << "," << secTrkEnd.Y() << "," << secTrkEnd.Z() << ")";
        const float secTrkStartDistToPrimEnd = (secTrkStart-primTrkEnd).R();
        const float secTrkEndDistToPrimEnd = (secTrkEnd-primTrkEnd).R();
        trackStartDistToPrimTrkEnd[iTrack] = secTrkStartDistToPrimEnd;
        trackEndDistToPrimTrkEnd[iTrack] = secTrkEndDistToPrimEnd;
        trackClosestDistToPrimTrkEnd[iTrack] = secTrkEndDistToPrimEnd;
        trackStartClosestToPrimTrkEnd[iTrack] = false;
	if(secTrkStartDistToPrimEnd <= secTrkEndDistToPrimEnd)
	{
          trackClosestDistToPrimTrkEnd[iTrack] = secTrkStartDistToPrimEnd;
          trackStartClosestToPrimTrkEnd[iTrack] = true;
	}
	  
        const auto secTrkCalos = tracksCaloVec.at(iTrack);
        for(const auto& secTrkCalo:secTrkCalos)
        {
          if(secTrkCalo->PlaneID().Plane == fCaloPlane)
	  {
            size_t secTrkdEdxs = secTrkCalo->dEdx().size();
	    if(secTrkdEdxs > 4) SecTrkPID[iTrack] = true;
	  }
	}

        const double secTrkLLR = trackLLHPion[iTrack] - trackLLHProton[iTrack];
        const double secTrkLLRPro = trackLLHProton[iTrack] - trackLLHPion[iTrack];
        if(trackClosestDistToPrimTrkEnd[iTrack] < 2.5)
        {
          iSecTrkID[nSecTrk] = iTrack;
          nSecTrk++;
          if(SecTrkPID[iTrack]) 
          {
	    if(secTrkLLRPro > nSecLLRProtonToPionMax) nSecLLRProtonToPionMax = secTrkLLRPro;
	    if(secTrkLLRPro < nSecLLRProtonToPionMin) 
            {
	      nSecLLRProtonToPionMin = secTrkLLRPro;
	      iSecMin = iTrack;
            }
	  }
          if(secTrkLLR > 0.) nSecTrkLLRG0++;
          if(secTrkLLR > 100.) nSecTrkLLRG100++;
          if(secTrkLLR > 200.) nSecTrkLLRG200++;
          if(secTrkLLR > 300.) nSecTrkLLRG300++;
          if(secTrkLLR > 400.) nSecTrkLLRG400++;
          if(secTrkLLR > 500.) nSecTrkLLRG500++;
          if(secTrkLLR > 600.) nSecTrkLLRG600++;
          if(secTrkLLR > 700.) nSecTrkLLRG700++;
          if(secTrkLLR > 800.) nSecTrkLLRG800++;
          if(secTrkLLR > 900.) nSecTrkLLRG900++;
          if(trackPIDA[iTrack] < 8 ) nSecTrkPIDAL8++;
          if(trackPIDA[iTrack] < 10 ) nSecTrkPIDAL10++;
          if(trackPIDA[iTrack] < 14 ) nSecTrkPIDAL14++;
          if(trackPIDA[iTrack] < 16 ) nSecTrkPIDAL16++;
          if(trackPIDA[iTrack] < 18 ) nSecTrkPIDAL18++;
        }
      } //if this track is not the primary
    } // for track in trackVec

//    // SecMin Calorimetry
//    int itc = tracksCaloVec.size();
//    if(iSecMin >= 0 && iSecMin < itc)
//    {
//      const auto secMinTrkCalos = tracksCaloVec.at(iSecMin);
//      for(const auto& secMinTrkCalo:secMinTrkCalos)
//      {
//        if(secMinTrkCalo->PlaneID().Plane == fCaloPlane)
//        {
//          size_t IBackwards = secMinTrkCalo->dEdx().size()-1;
//          for(size_t cRangeIt = 0; cRangeIt < secMinTrkCalo->ResidualRange().size() && cRangeIt < secMinTrkCalo->dEdx().size(); cRangeIt++)
//          {
//            secMinTrkResRanges.push_back(secMinTrkCalo->ResidualRange().at(cRangeIt));
//            secMinTrkdEdxs.push_back(secMinTrkCalo->dEdx().at(cRangeIt));
//            secMinTrkPitches.push_back(secMinTrkCalo->TrkPitchVec().at(cRangeIt));
//            secMinTrkIBackwards.push_back(IBackwards);
//            IBackwards--;
//            const auto thisPoint = secMinTrkCalo->XYZ().at(cRangeIt); // PositionVector3D
//            secMinTrkXs.push_back(thisPoint.X());
//            secMinTrkYs.push_back(thisPoint.Y());
//            secMinTrkZs.push_back(thisPoint.Z());
//            bool thisInFid = InPrimaryFiducial(thisPoint);
//            secMinTrkInFids.push_back(thisInFid);
//          } // for cRangeIt
//        } // if plane == fCaloPlane
//      } // for calo in caloVec
//    } // have good secondary

  } // good primaryTrack
  
} // ProcessPrimaryTrack

void lana::PionAbsSelector::ProcessPFParticles(const art::Event& e,
        //const art::Ptr<simb::MCParticle>& primaryParticle,
        //const pdana::MotherDaughterWalkerAlg& motherDaughterWalker,
        const pdana::MCBeamOrCosmicAlg * const beamOrCosmic
        )
{
  /**
  LArPandoraHelper larPandoraHelper;
  lar_pandora::PFParticleVector allPFParticles;
  lar_pandora::PFParticlesToMetadata pfPartsToMetadata;
  larPandoraHelper.CollectPFParticleMetadata(e,fPFParticleTag.encode(),allPFParticles,pfPartsToMetadata);
  lar_pandora::PFParticleMap idToPFParticleMap; // map int to Ptr<PFParticle>
  larPandoraHelper.BuildPFParticleMap(allPFParticles,pfPartsToTracks);
  lar_pandora::TrackVector allPFTracks; // vector<Ptr<Track>>
  lar_pandora::PFParticlesToTracks pfPartsToTracks; // map of Ptr<PFParticle> to Ptr<Track>
  larPandoraHelper.CollectTracks(e,fPFParticleTag.encode(),allPFTracks,pfPartsToTracks);
  lar_pandora::ShowerVector allPFShowers; // vector<Ptr<Shower>>
  lar_pandora::PFParticlesToShowers pfPartsToShowers; // map of Ptr<PFParticle> to Ptr<Shower>
  larPandoraHelper.CollectShowers(e,fPFParticleTag.encode(),allPFShowers,pfPartsToShowers);
  **/

  protoana::ProtoDUNEPFParticleUtils pfPartUtils;
  std::cout << "Number of primary PFParticles: " << pfPartUtils.GetNumberPrimaryPFParticle(e,fPFParticleTag.encode()) << std::endl;
  const std::vector<const recob::PFParticle*> pfFromBeamSlice = pfPartUtils.GetPFParticlesFromBeamSlice(e,fPFParticleTag.encode());

  // All this just to get the calos and truth
  auto allPFTrackHand = e.getValidHandle<std::vector<recob::Track>>(fPFTrackTag);
  std::vector<art::Ptr<recob::Track>> allPFTrackVec;
  if(allPFTrackHand.isValid())
  {
    art::fill_ptr_vector(allPFTrackVec, allPFTrackHand);
  }
  art::FindManyP<recob::Hit> fmHitsForPFTracks(allPFTrackHand, e, fPFTrackTag);
  //art::FindManyP<recob::Hit, recob::TrackHitMeta> fmHitsForPFTracksMeta(allPFTrackHand, e,fPFTrackTag);
  art::FindManyP<anab::Calorimetry>  fmPFCalo(allPFTrackHand, e, fPFCaloTag);
  //auto allPFShowerHand = e.getValidHandle<std::vector<recob::Shower>>(fPFShowerTag);
  //std::vector<art::Ptr<recob::Shower>> allPFShowerVec;
  //if(allPFShowerHand.isValid())
  //{
  //  art::fill_ptr_vector(allPFShowerVec, allPFShowerHand);
  //}
  //art::FindManyP<recob::Hit> fmHitsForPFShowers(allPFShowerHand, e, fPFShowerTag);
  auto allPandoraHitHand = e.getValidHandle<std::vector<recob::Hit>>(fPandoraHitTag);
  std::vector<art::Ptr<recob::Hit>> allPandoraHits;
  art::fill_ptr_vector(allPandoraHits, allPandoraHitHand);

  trkf::TrackMomentumCalculator trackMomCalc;

  PFNBeamSlices = pfFromBeamSlice.size();
  std::cout << "Number of primary beam PFParticles: " << PFNBeamSlices << std::endl;
  for(size_t iPF=0; iPF < PFNBeamSlices; iPF++)
  {
    const recob::PFParticle* pfBeamPart = pfFromBeamSlice.at(iPF);
    PFBeamPrimPDG = pfBeamPart->PdgCode();
    PFBeamPrimNDaughters = pfBeamPart->NumDaughters();
    const bool isPFParticleTracklike = pfPartUtils.IsPFParticleTracklike(*pfBeamPart);
    const bool isPFParticleShowerlike = pfPartUtils.IsPFParticleShowerlike(*pfBeamPart);
    PFBeamPrimIsTracklike = isPFParticleTracklike;
    PFBeamPrimIsShowerlike = isPFParticleShowerlike;
    PFBeamPrimBeamCosmicScore = pfPartUtils.GetBeamCosmicScore(*pfBeamPart,e,fPFParticleTag.encode());


    std::cout << "Beam PFParticle: "<< iPF 
                << " is primary: " << pfBeamPart->IsPrimary() 
                << " PDG: " << pfBeamPart->PdgCode()
                << " is track-like: "<<isPFParticleTracklike
                << " is shower-like: "<<isPFParticleShowerlike
                << std::endl;

    bool primStartValid = false;
    bool primEndValid = false;
    recob::Track::Point_t pfBeamPrimStart;
    recob::Track::Point_t pfBeamPrimEnd;
    recob::Track::Vector_t pfBeamPrimStartDir;
    recob::Track::Vector_t pfBeamPrimEndDir;
    if(isPFParticleTracklike)
    {
      const recob::Track* pfTrack = pfPartUtils.GetPFParticleTrack(*pfBeamPart, e, fPFParticleTag.encode(),fPFTrackTag.encode());
      
      if(pfTrack)
      {
        pfBeamPrimStart = pfTrack->Start();
        pfBeamPrimStartDir = pfTrack->StartDirection();
        primStartValid = true;

        if(pfTrack->NPoints() > 1)
        {
          primEndValid = true;
          pfBeamPrimEnd = pfTrack->End();
          pfBeamPrimEndDir = pfTrack->EndDirection();
          PFBeamPrimTrkLen = pfTrack->Length();
          PFBeamPrimCSDAMomMu = trackMomCalc.GetTrackMomentum(PFBeamPrimTrkLen,13) * 1000; // MeV/c
          PFBeamPrimCSDAMomProton = trackMomCalc.GetTrackMomentum(PFBeamPrimTrkLen,2212) * 1000; // MeV/c
          PFBeamPrimCSDAKinMu = sqrt(pow(PFBeamPrimCSDAMomMu,2)+pow(MMUON,2))-MMUON; // MeV
          PFBeamPrimCSDAKinProton = sqrt(pow(PFBeamPrimCSDAMomProton,2)+pow(MPROTON,2))-MPROTON; // MeV
        }

        std::cout << "PFBeamPrimTrk ID: " << pfTrack->ID() << std::endl;
        std::cout << "PFBeamPrimTrk Start Pos: "
                    << pfBeamPrimStart.X() << ", " << pfBeamPrimStart.Y() << ", " << pfBeamPrimStart.Z()
                    << " End Pos: "
                    << pfBeamPrimEnd.X() << ", " << pfBeamPrimEnd.Y() << ", " << pfBeamPrimEnd.Z()
                    << std::endl;
        std::cout << "PFBeamPrimTrk Start Dir: "
                    << pfBeamPrimStartDir.X() << ", " << pfBeamPrimStartDir.Y() << ", " << pfBeamPrimStartDir.Z()
                    << " End Dir: "
                    << pfBeamPrimEndDir.X() << ", " << pfBeamPrimEndDir.Y() << ", " << pfBeamPrimEndDir.Z()
                    << std::endl;
        std::cout << "PFBeamPrimTrk Start Theta & Phi (deg): " << pfBeamPrimStartDir.Theta()*180./CLHEP::pi << ", " << pfBeamPrimStartDir.Phi()*180./CLHEP::pi << std::endl;
        std::cout << "PFBeamPrimTrk End Theta & Phi (deg): " << pfBeamPrimEndDir.Theta()*180./CLHEP::pi << ", " << pfBeamPrimEndDir.Phi()*180./CLHEP::pi << std::endl;

        if(primEndValid && pfBeamPrimStart.Z() > pfBeamPrimEnd.Z())
        {
            const auto tmpPoint = pfBeamPrimStart;
            pfBeamPrimStart = pfBeamPrimEnd;
            pfBeamPrimEnd = tmpPoint;
            const auto tmpVect = pfBeamPrimStartDir;
            pfBeamPrimStartDir = pfBeamPrimEndDir;
            pfBeamPrimEndDir = tmpVect;
            pfBeamPrimStartDir *= -1.;
            pfBeamPrimEndDir *= -1.;
        }
        std::cout << "After possibly flipping:" << std::endl;
        std::cout << "PFBeamPrimTrk Start Pos: "
                    << pfBeamPrimStart.X() << ", " << pfBeamPrimStart.Y() << ", " << pfBeamPrimStart.Z()
                    << " End Pos: "
                    << pfBeamPrimEnd.X() << ", " << pfBeamPrimEnd.Y() << ", " << pfBeamPrimEnd.Z()
                    << std::endl;
        std::cout << "PFBeamPrimTrk Start Dir: "
                    << pfBeamPrimStartDir.X() << ", " << pfBeamPrimStartDir.Y() << ", " << pfBeamPrimStartDir.Z()
                    << " End Dir: "
                    << pfBeamPrimEndDir.X() << ", " << pfBeamPrimEndDir.Y() << ", " << pfBeamPrimEndDir.Z()
                    << std::endl;
        std::cout << "PFBeamPrimTrk Start Theta & Phi (deg): " << pfBeamPrimStartDir.Theta()*180./CLHEP::pi << ", " << pfBeamPrimStartDir.Phi()*180./CLHEP::pi << std::endl;
        std::cout << "PFBeamPrimTrk End Theta & Phi (deg): " << pfBeamPrimEndDir.Theta()*180./CLHEP::pi << ", " << pfBeamPrimEndDir.Phi()*180./CLHEP::pi << std::endl;

        const auto& pfTrackTraj = pfTrack->Trajectory();
        for(size_t iTP=0; iTP < pfTrackTraj.NumberTrajectoryPoints()-1; iTP++)
        {
          if(pfTrackTraj.HasValidPoint(iTP) && pfTrackTraj.HasValidPoint(iTP+1))
          {
            float thisKink = ROOT::Math::VectorUtil::Angle(pfTrackTraj.DirectionAtPoint(iTP),pfTrackTraj.DirectionAtPoint(iTP+1));
            //const auto& trajFlags = pfTrackTraj.FlagsAtPoint(iTP); // seems empty for pandoraTrack
            PFBeamPrimTrkMaxKink = std::max(PFBeamPrimTrkMaxKink,thisKink);
          }
        } // for iTP

        art::ServiceHandle<geo::Geometry> geom;
        const auto& pfTrackHits = fmHitsForPFTracks.at(pfTrack->ID());
        if(isMC) // Now MCTruth Matching
        {
          //std::map<size_t,size_t> nHitsPerPlane;
          //for(size_t iHit=0; iHit < pfTrackHits.size(); iHit++)
          //{
          //  const auto& iPlane = pfTrackHits.at(iHit)->WireID().Plane;
          //  if(nHitsPerPlane.count(iPlane) > 0)
          //  {
          //      nHitsPerPlane.at(iPlane)++;
          //  }
          //  else
          //  {
          //      nHitsPerPlane.emplace(iPlane,1);
          //  }
          //}
          //for(auto& nHitsPer: nHitsPerPlane)
          //{
          //  std::cout << "For plane: " << nHitsPer.first << " there are " << nHitsPer.second << std::endl;
          //}
          lsu::BackTrackerAlg btAlg;
          const auto pfTrackMCPart = btAlg.getMCParticle(pfTrackHits,e,PFBeamPrimTrueTrackIDIsNeg);
          const auto& pfTrackMCTraj = pfTrackMCPart.Trajectory();
          PFBeamPrimTrueTrackID = pfTrackMCPart.TrackId();
          PFBeamPrimTrueMotherTrackID = pfTrackMCPart.Mother();
          PFBeamPrimTruePDG = pfTrackMCPart.PdgCode();
          PFBeamPrimTrueKin = 1000*(pfTrackMCPart.E()-pfTrackMCPart.Mass());
          if(pfTrackMCTraj.size()>2)
          {
            PFBeamPrimTrueEndKin = 1000*(std::prev(std::prev(pfTrackMCTraj.end()))->first.E()-pfTrackMCPart.Mass());
          }
          PFBeamPrimTrueTrajLen = pfTrackMCTraj.TotalLength();
          if (beamOrCosmic)
          {
            PFBeamPrimTrueIsBeam = beamOrCosmic->isBeam(pfTrackMCPart);
          }

          PFBeamPrimTrueStartX = pfTrackMCPart.Vx();
          PFBeamPrimTrueStartY = pfTrackMCPart.Vy();
          PFBeamPrimTrueStartZ = pfTrackMCPart.Vz();
          PFBeamPrimTrueStartT = pfTrackMCPart.T();

          PFBeamPrimTrueEndX = pfTrackMCPart.EndX();
          PFBeamPrimTrueEndY = pfTrackMCPart.EndY();
          PFBeamPrimTrueEndZ = pfTrackMCPart.EndZ();
          PFBeamPrimTrueEndT = pfTrackMCPart.EndT();

          for(size_t iTP=0; iTP < pfTrackMCTraj.size()-1; iTP++)
          {
            float thisKink = ROOT::Math::VectorUtil::Angle(
                                    lsu::toVector(pfTrackMCTraj.Momentum(iTP)),
                                    lsu::toVector(pfTrackMCTraj.Momentum(iTP+1))
                    );
            PFBeamPrimTrueMaxKink = std::max(PFBeamPrimTrueMaxKink,thisKink);
          } // for iTP

          std::cout << "PFPrimaryTrack MCPart: TrackId: "<< PFBeamPrimTrueTrackID 
                    <<" PDG: " << pfTrackMCPart.PdgCode() 
                    <<"  Mom: "<< pfTrackMCPart.Momentum().Vect().Mag() 
                    << " isBeam: "<<PFBeamPrimTrueIsBeam
                    <<" process: "<<pfTrackMCPart.Process() 
                    << " end process: "<<pfTrackMCPart.EndProcess()
                    << std::endl;

        } // if isMC

        // For truth calo info
        std::vector<art::Ptr<sim::SimChannel>> simChanVec;
        if(isMC)
        {
          auto simChanHand = e.getValidHandle<std::vector<sim::SimChannel> >("largeant");
          art::fill_ptr_vector(simChanVec, simChanHand);
        } // isMC

        const auto& shouldBeThisTrackToo = allPFTrackVec.at(pfTrack->ID());
        if(pfTrack->ID() != shouldBeThisTrackToo->ID() 
            || pfTrack->NPoints() != shouldBeThisTrackToo->NPoints())
        {
          throw cet::exception("PionAbsSelector","track->ID() isn't the track index in the list for primary track!");
        }
        const auto& pfTrackCalos = fmPFCalo.at(pfTrack->ID());
        for(const auto& pfTrackCalo:pfTrackCalos)
        {
          if(pfTrackCalo->PlaneID().Plane == fCaloPlane)
          {
            const auto& dEdxSize = pfTrackCalo->dEdx().size();
            //const auto& dEdxBegin = pfTrackCalo->dEdx().begin();
            const auto& dEdxEnd = pfTrackCalo->dEdx().end();
            if (dEdxSize >= 3)
            {
              PFBeamPrimdEdxAverageLast3Hits = findAverage(dEdxEnd - 3, dEdxEnd);
            }
            if (dEdxSize >= 5)
            {
              PFBeamPrimdEdxAverageLast5Hits = findAverage(dEdxEnd - 5, dEdxEnd);
            }
            if (dEdxSize >= 7)
            {
              PFBeamPrimdEdxAverageLast7Hits = findAverage(dEdxEnd - 7, dEdxEnd);
            }

            std::cout << "N Valid Traj points:: " << pfTrackTraj.CountValidPoints() << std::endl;
            std::cout << "pfTrackHits.size(): " << pfTrackHits.size() << std::endl;
            std::cout << "  first Hit WireID: " << pfTrackHits.front()->WireID() << std::endl;
            std::cout << "  last Hit WireID:  " << pfTrackHits.back()->WireID() << std::endl;
            std::cout  << "Calo dEdxSize: " << dEdxSize 
                        << " XYZ.size(): " << pfTrackCalo->XYZ().size()
                        << " TpIndices.size(): " << pfTrackCalo->TpIndices().size() << std::endl;
            for(size_t cRangeIt = 0; cRangeIt < pfTrackCalo->ResidualRange().size() 
                                && cRangeIt < dEdxSize; cRangeIt++)
            {
              PFBeamPrimResRanges.push_back(pfTrackCalo->ResidualRange().at(cRangeIt));
              PFBeamPrimdEdxs.push_back(pfTrackCalo->dEdx().at(cRangeIt));
              PFBeamPrimdQdxs.push_back(pfTrackCalo->dQdx().at(cRangeIt));
              PFBeamPrimPitches.push_back(pfTrackCalo->TrkPitchVec().at(cRangeIt));
              const auto& thisPoint = pfTrackCalo->XYZ().at(cRangeIt); // PositionVector3D

              PFBeamPrimXs.push_back(thisPoint.X());
              PFBeamPrimYs.push_back(thisPoint.Y());
              PFBeamPrimZs.push_back(thisPoint.Z());
              int iZWire = DEFAULTNEG;

              if(cRangeIt < pfTrackCalo->TpIndices().size()) 
              {
                const size_t& hitIndex = pfTrackCalo->TpIndices().at(cRangeIt);
                const auto& thisHit = allPandoraHits.at(hitIndex);
                const auto& thisHitChan = thisHit->Channel();
                //const auto& thisHitWireID = thisHit->WireID();
                //const auto& thisHitWire = geom->Wire(thisHitWireID);
                //const auto& locationAtPoint = pfTrackTraj.LocationAtPoint(hitIndex);
                //std::cout << "calo i: " << cRangeIt 
                //            << " hitIndex: "<< hitIndex
                //            << std::endl;
                //std::cout << "    calo XYZ: "
                //    << thisPoint.X() << ", " << thisPoint.Y() << ", " << thisPoint.Z()
                //    << std::endl;
                //std::cout << "    traj XYZ: "
                //    << locationAtPoint.X() << ", " << locationAtPoint.Y() << ", " << locationAtPoint.Z()
                //    << std::endl;
                //std::cout << "    Hit channel: " << thisHit->Channel() 
                //    << " wire: " << thisHitWireID
                //    << " wire Z: " << thisHitWire.GetCenter().Z()
                //    << " peak time: " << thisHit->PeakTime()
                //    << " time RMS: " << thisHit->RMS()
                //    << std::endl;
                if(myChannelToWireMap.count(thisHitChan) > 0)
                {
                  iZWire = myChannelToWireMap[thisHitChan].first;
                  zWireNHits.at(iZWire)++;
                  zWireResRange.at(iZWire) = pfTrackCalo->ResidualRange().at(cRangeIt);
                  zWiredEdx.at(iZWire) = pfTrackCalo->dEdx().at(cRangeIt);
                  zWiredQdx.at(iZWire) = pfTrackCalo->dQdx().at(cRangeIt);
                  zWirePitch.at(iZWire) = pfTrackCalo->TrkPitchVec().at(cRangeIt);
                  zWireX.at(iZWire) = thisPoint.X();
                  zWireY.at(iZWire) = thisPoint.Y();
                  zWireZ.at(iZWire) = thisPoint.Z();
                }
                if(isMC && false)
                {
                  for(const auto& simChan: simChanVec)
                  {
                    if(simChan->Channel() == thisHitChan)
                    {
                      for(const auto& TDCIDEs: simChan->TDCIDEMap())
                      {
                        unsigned short TDC = TDCIDEs.first;
                        const std::vector<sim::IDE>& IDEs = TDCIDEs.second;
                        for(const auto& IDE: IDEs)
                        {
                          if(PFBeamPrimTrueTrackID == abs(IDE.trackID))
                          {
                            std::cout << "    TDC: " << int(TDC) 
                                    << " TrackID: " << IDE.trackID
                                    << " Energy: " << IDE.energy
                                    << " XYZ: " << IDE.x << ", " << IDE.y << ", " << IDE.z
                                    << std::endl;
                          } // if PFBeamPrimTrueTrackID == IDE.trackID
                        } // for IDE
                      } // for TDCIDE
                    } // if simChan->Channel == hitChannel
                  } // for simChan
                } // isMC
              } // if cRangeIt < TpIndices.size()
              PFBeamPrimZWires.push_back(iZWire);
            } // for cRangeIt
          } // if Plane is fCaloPlane
        } // for pfTrackCalo
        //std::cout << "Final Calos in Tree:\n";
        if(PFBeamPrimZs.size()>1)
        {
          if(PFBeamPrimZs.front() > PFBeamPrimZs.back())
          {
            std::cout << "Calorimetry was backwards so reversing!\n";
            std::reverse(PFBeamPrimXs.begin(),PFBeamPrimXs.end());
            std::reverse(PFBeamPrimYs.begin(),PFBeamPrimYs.end());
            std::reverse(PFBeamPrimZs.begin(),PFBeamPrimZs.end());
            std::reverse(PFBeamPrimResRanges.begin(),PFBeamPrimResRanges.end());
            std::reverse(PFBeamPrimdEdxs.begin(),PFBeamPrimdEdxs.end());
            std::reverse(PFBeamPrimdQdxs.begin(),PFBeamPrimdQdxs.end());
            std::reverse(PFBeamPrimPitches.begin(),PFBeamPrimPitches.end());
            std::reverse(PFBeamPrimZWires.begin(),PFBeamPrimZWires.end());
          }

          double intE = 0.;
          for(size_t iCaloHit=0; iCaloHit < PFBeamPrimdEdxs.size();iCaloHit++)
          {
              const double thisKin = kinWCInTPC-intE;
              PFBeamPrimKins.push_back(thisKin);
              const double thisKinProton = kinWCInTPCProton-intE;
              PFBeamPrimKinsProton.push_back(thisKinProton);
              PFBeamPrimKinInteract = thisKin;
              PFBeamPrimKinInteractProton = thisKinProton;
              intE += PFBeamPrimdEdxs.at(iCaloHit)*PFBeamPrimPitches.at(iCaloHit);
              //std::cout << "Calo i: " << iCaloHit
              //          << "    iZWire: " << PFBeamPrimZWires.at(iCaloHit) 
              //          << "    z: " << PFBeamPrimZs.at(iCaloHit) 
              //          << "    res range: " << PFBeamPrimResRanges.at(iCaloHit)
              //          << "    dEdx: " << PFBeamPrimdEdxs.at(iCaloHit)
              //          << "    pitch: " << PFBeamPrimPitches.at(iCaloHit)
              //          << "    Kin: " << PFBeamPrimKins.at(iCaloHit)
              //          << std::endl;
          }
        }
        // now z-wire kin energy
        zWireEnergySum = 0.;
        int lastZWire = -1;
        for (long iZWire=0; iZWire < (long) zWiredEdx.size(); iZWire++)
        {
          if(zWiredEdx.at(iZWire) >= 0.) lastZWire = iZWire;
        }
        for (long iZWire=0; iZWire <= lastZWire; iZWire++)
        {
          zWirePartKin.at(iZWire) = kinWCInTPC - zWireEnergySum;
          zWirePartKinProton.at(iZWire) = kinWCInTPCProton - zWireEnergySum;
          if(zWiredEdx.at(iZWire) >= 0.)
          {
            zWireEnergySum += zWiredEdx.at(iZWire) * zWirePitch.at(iZWire);
          }
        } // for iZWire
      } // if pfTrack
    } // if isPFParticleTracklike
    if(isPFParticleShowerlike)
    {
      const recob::Shower* pfShower = pfPartUtils.GetPFParticleShower(*pfBeamPart, e, fPFParticleTag.encode(),fPFShowerTag.encode());
      if(pfShower)
      {
        pfBeamPrimStart = lsu::toPoint(pfShower->ShowerStart());
        pfBeamPrimStartDir = lsu::toVector(pfShower->Direction());
        primStartValid = true;

        if(pfShower->has_open_angle())
        {
          PFBeamPrimShwrOpenAngle = pfShower->OpenAngle();
        }
        if(pfShower->has_length())
        {
          PFBeamPrimShwrLen = pfShower->Length();
        }
        std::cout << "PFBeamPrimShwr Start Pos: "
                    << pfBeamPrimStart.X() << ", " << pfBeamPrimStart.Y() << ", " << pfBeamPrimStart.Z()
                    << std::endl;
        std::cout << "PFBeamPrimShwr Start Dir: "
                    << pfBeamPrimStartDir.X() << ", " << pfBeamPrimStartDir.Y() << ", " << pfBeamPrimStartDir.Z()
                    << std::endl;
        std::cout << "    Shower: "
                  << " has open angle: " << pfShower->has_open_angle()
                  << " has length: " << pfShower->has_length()
                  << " open angle: " << pfShower->OpenAngle()
                  << " length: " << pfShower->Length()
                  << std::endl;
        //if(isMC) // Now MCTruth Matching
        //{
        //  const auto& pfShowerHits = fmHitsForPFShowers.at(pfShower->ID());
        //  lsu::BackTrackerAlg btAlg;
        //  bool trackIDIsNeg = false;
        //  const auto pfShowerMCPart = btAlg.getMCParticle(pfShowerHits,e,trackIDIsNeg);
        //  bool isTrueBeam = false;
        //  if (beamOrCosmic)
        //  {
        //    isTrueBeam = beamOrCosmic->isBeam(pfShowerMCPart);
        //  }
        //  std::cout << "PFPrimaryShower MCPart: PDG: " << pfShowerMCPart.PdgCode() <<"  Mom: "<< pfShowerMCPart.Momentum().Vect().Mag() << " isBeam: "<<isTrueBeam<<" process: "<<pfShowerMCPart.Process() << " end process: "<<pfShowerMCPart.EndProcess()<< std::endl;

        //} // if isMC
      } // if pfShower
    } // if isPFParticleShowerlike

    const ROOT::Math::Polar3DVector dirVecBeamTrk(1.,thetaWC,phiWC);
    if (primStartValid)
    {
      PFBeamPrimStartX = pfBeamPrimStart.X();
      PFBeamPrimStartY = pfBeamPrimStart.Y();
      PFBeamPrimStartZ = pfBeamPrimStart.Z();
      PFBeamPrimStartTheta = pfBeamPrimStartDir.Theta();
      PFBeamPrimStartPhi = pfBeamPrimStartDir.Phi();
      const TVector3 pfFrontTPCPoint = lsu::lineZPlane(ZSTARTOFTPC,pfBeamPrimStart,pfBeamPrimStartDir);
      PFBeamPrimXFrontTPC = pfFrontTPCPoint.X();
      PFBeamPrimYFrontTPC = pfFrontTPCPoint.Y();
      PFBeamPrimAngleToBeamTrk = ROOT::Math::VectorUtil::Angle(pfBeamPrimStartDir,dirVecBeamTrk);
    }
    if (primEndValid)
    {
      PFBeamPrimEndX = pfBeamPrimEnd.X();
      PFBeamPrimEndY = pfBeamPrimEnd.Y();
      PFBeamPrimEndZ = pfBeamPrimEnd.Z();
      PFBeamPrimEndTheta = pfBeamPrimEndDir.Theta();
      PFBeamPrimEndPhi = pfBeamPrimEndDir.Phi();
      const TVector3 pfFrontTPCPointTrackEnd = lsu::lineZPlane(ZSTARTOFTPC,pfBeamPrimEnd,pfBeamPrimEndDir);
      PFBeamPrimXFrontTPCTrackEnd = pfFrontTPCPointTrackEnd.X();
      PFBeamPrimYFrontTPCTrackEnd = pfFrontTPCPointTrackEnd.Y();
      PFBeamPrimTrkStartEndDirAngle = ROOT::Math::VectorUtil::Angle(pfBeamPrimStartDir,pfBeamPrimEndDir);
      PFBeamPrimEndAngleToBeamTrk = ROOT::Math::VectorUtil::Angle(pfBeamPrimEndDir,dirVecBeamTrk);
    }

    const auto& pfBeamDaughterTracks = pfPartUtils.GetPFParticleDaughterTracks(*pfBeamPart,e,fPFParticleTag.encode(),fPFTrackTag.encode());
    const auto& pfBeamDaughterShowers = pfPartUtils.GetPFParticleDaughterShowers(*pfBeamPart,e,fPFParticleTag.encode(),fPFShowerTag.encode());
    PFBeamPrimNDaughterTracks = pfBeamDaughterTracks.size();
    PFBeamPrimNDaughterShowers = pfBeamDaughterShowers.size();
    
    std::cout << "    N daughters: " << pfBeamPart->NumDaughters()
                << " N daughter tracks: " << pfBeamDaughterTracks.size()
                << " N daughter showers: " << pfBeamDaughterShowers.size()
                << std::endl;
    std::cout << "    Beam Prim Start:           "
                << "(" << pfBeamPrimStart.X()
                << "," << pfBeamPrimStart.Y()
                << "," << pfBeamPrimStart.Z()
                << ")"
                << std::endl;
    std::cout << "    Vertex:           "
                << "(" << pfBeamPrimEnd.X()
                << "," << pfBeamPrimEnd.Y()
                << "," << pfBeamPrimEnd.Z()
                << ")"
                << std::endl;
    for(size_t iSec=0; iSec < pfBeamDaughterTracks.size(); iSec++)
    {
        const auto& pfBeamDaughterTrack = pfBeamDaughterTracks.at(iSec);
        std::cout << "Daughter "<<iSec<<" track len: " << pfBeamDaughterTrack->Length() << std::endl;
        PFBeamSecTrkLen[iSec] = pfBeamDaughterTrack->Length();
        if(primEndValid) // there was an end point (and is primary track)
        {
          PFBeamSecTrkDistToPrimEnd[iSec] = (pfBeamDaughterTrack->Start()-pfBeamPrimEnd).R();

          PFBeamSecTrkAngleToPrimEnd[iSec] = ROOT::Math::VectorUtil::Angle(pfBeamPrimEndDir,pfBeamDaughterTrack->StartDirection());
          PFBeamPrimMaxAngleToSecTrks = std::max(PFBeamSecTrkAngleToPrimEnd[iSec],PFBeamPrimMaxAngleToSecTrks);
        }

        const auto& shouldBeThisTrackToo = allPFTrackVec.at(pfBeamDaughterTrack->ID());
        if(pfBeamDaughterTrack->ID() != shouldBeThisTrackToo->ID() 
            || pfBeamDaughterTrack->NPoints() != shouldBeThisTrackToo->NPoints())
        {
          throw cet::exception("PionAbsSelector","track->ID() isn't the track index in the list for daughter track!");
        }
        const auto& pfBeamDaughterTrackCalos = fmPFCalo.at(pfBeamDaughterTrack->ID());

        for(const auto& pfTrackCalo:pfBeamDaughterTrackCalos)
        {
          if(pfTrackCalo->PlaneID().Plane == fCaloPlane)
          {
            const auto& dEdxSize = pfTrackCalo->dEdx().size();
            //const auto& dEdxBegin = pfTrackCalo->dEdx().begin();
            const auto& dEdxEnd = pfTrackCalo->dEdx().end();
            if (dEdxSize >= 3)
            {
              PFBeamSecTrkdEdxAverageLast3Hits[iSec] = findAverage(dEdxEnd - 3, dEdxEnd);
            }
            if (dEdxSize >= 5)
            {
              PFBeamSecTrkdEdxAverageLast5Hits[iSec] = findAverage(dEdxEnd - 5, dEdxEnd);
            }
            if (dEdxSize >= 7)
            {
              PFBeamSecTrkdEdxAverageLast7Hits[iSec] = findAverage(dEdxEnd - 7, dEdxEnd);
            }
          } // if Plane is fCaloPlane
        } // for pfTrackCalo
    } // for pfBeamDaughterTrack
    
    break; // only look at first PFSlice
  } // for iPF
} // ProcessPFParticles

DEFINE_ART_MODULE(lana::PionAbsSelector)
