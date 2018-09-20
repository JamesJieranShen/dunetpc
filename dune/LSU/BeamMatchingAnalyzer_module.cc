////////////////////////////////////////////////////////////////////////
// Class:       BeamMatchingAnalyzer
// Module Type: analyzer
// File:        BeamMatchingAnalyzer_module.cc
//
// Created from BeamMatchingAnalyzer_module.cc on 2018-09-18 by Justin Hugon
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
#define MAXIDES 150000
#define MAXZINT 95
#define MAXLINT 20

#define MCHARGEDPION 139.57018 // MeV/c^2
#define MPROTON 938.2720813 // MeV/c^2
#define KINLOSTBEFORETPC 0.0 //MeV; from LArIAT pion total cross section group
#define KINLOSTBEFORETPCPROTON 0.0 //MeV; from LArIAT pion total cross section group

#define setHistTitles(hist,xtitle,ytitle) hist->GetXaxis()->SetTitle(xtitle); hist->GetYaxis()->SetTitle(ytitle);

namespace lana {
  class BeamMatchingAnalyzer;
}

class lana::BeamMatchingAnalyzer : public art::EDAnalyzer {
public:
  explicit BeamMatchingAnalyzer(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamMatchingAnalyzer(BeamMatchingAnalyzer const &) = delete;
  BeamMatchingAnalyzer(BeamMatchingAnalyzer &&) = delete;
  BeamMatchingAnalyzer & operator = (BeamMatchingAnalyzer const &) = delete;
  BeamMatchingAnalyzer & operator = (BeamMatchingAnalyzer &&) = delete;

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
  art::InputTag fTruePartLabel; //The name of the module that produced simb::MCParticle objects
  art::InputTag fBeamTruthTag; //The name of the module that produced beam simb::MCTruth objects
  art::InputTag fCosmicTruthTag; //The name of the module that produced cosmic simb::MCTruth objects
  art::InputTag fTrackLabel; //The name of the module that produced recob::Track objects
  art::InputTag fTrackHitTag;
  art::InputTag fAllHitTag;

  /////////////////////////////
  //Tree variables
  TTree* tree;

  bool isMC;
  UInt_t runNumber;
  UInt_t subRunNumber;
  UInt_t eventNumber;

  UInt_t nMCParts;
  bool mcPartIsBeam;
  bool mcPartIsPrimary;
  Int_t mcPartTrackID;
  Int_t mcPartPDG;
  Int_t mcPartGreatestGrandmotherTrackID;
  Int_t mcPartGreatestGrandmotherPDG;
  Float_t mcPartStartX;
  Float_t mcPartStartY;
  Float_t mcPartStartZ;
  Float_t mcPartStartT;
  Float_t mcPartEndX;
  Float_t mcPartEndY;
  Float_t mcPartEndZ;
  Float_t mcPartEndT;
  Float_t mcPartStartTheta;
  Float_t mcPartStartPhi;
  Float_t mcPartXFrontTPC;
  Float_t mcPartYFrontTPC;
  Float_t mcPartStartMom;
  Float_t mcPartStartE;
  Float_t mcPartStartKin;
  Float_t mcPartEndMom;
  Float_t mcPartEndE;
  Float_t mcPartEndKin;
  Float_t mcPartDeltaAngle;

  UInt_t nTracks;
  Float_t trackStartX;
  Float_t trackStartY;
  Float_t trackStartZ;
  Float_t trackStartTheta;
  Float_t trackStartPhi;
  Float_t trackEndX;
  Float_t trackEndY;
  Float_t trackEndZ;
  Float_t trackLength;
  Float_t trackXFrontTPC;
  Float_t trackYFrontTPC;
  Int_t trackTrueID;
  Int_t trackTrueMotherID;
  Int_t trackTrueGreatestGrandmotherTrackID;
  Int_t trackTrueGreatestGrandmotherPDG;
  Int_t trackTruePdg;
  bool trackTrueIsBeam;
  Float_t trackTrueKin;
  Float_t trackTrueEndKin;
  Float_t trackTrueTrajLen;
  Float_t trackTrueChargePurity;
  Float_t trackTrueChargeEfficiencyU;
  Float_t trackTrueChargeEfficiencyV;
  Float_t trackTrueChargeEfficiencyZ;

  Float_t trackTrueStartX;
  Float_t trackTrueStartY;
  Float_t trackTrueStartZ;
  Float_t trackTrueStartT;
  Float_t trackTrueEndX;
  Float_t trackTrueEndY;
  Float_t trackTrueEndZ;
  Float_t trackTrueEndT;
  Float_t trackTrueXFrontTPC;
  Float_t trackTrueYFrontTPC;
  Float_t trackMCPartAngle;

  art::ServiceHandle<cheat::BackTrackerService> bt;
 
  void ResetTreeVars();
  void ResetTreeTrackVars();

};


lana::BeamMatchingAnalyzer::BeamMatchingAnalyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  this->reconfigure(p);
}

void lana::BeamMatchingAnalyzer::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  art::ServiceHandle<geo::Geometry> geom;

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

  pdana::MCBeamOrCosmicAlg* beamOrCosmic;
  if(!e.isRealData())
  {
    beamOrCosmic = new pdana::MCBeamOrCosmicAlg(e,fTruePartLabel,fBeamTruthTag,fCosmicTruthTag);
  }

  auto trackHand = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
  std::vector<art::Ptr<recob::Track>> trackVec;
  if(trackHand.isValid())
  {
    art::fill_ptr_vector(trackVec, trackHand);
  }

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

  pdana::MotherDaughterWalkerAlg motherDaughterWalker(e,fTruePartLabel);

  //Get MCParticle Variables
  art::Ptr<simb::MCParticle> primaryParticle;
  for(const auto& truth:(truePartVec))
  {
    if (truth->PdgCode() == 2112 && truth->PdgCode() >= 1000000000) continue;
    ResetTreeVars();

    runNumber = e.run();
    subRunNumber = e.subRun();
    eventNumber = e.id().event();
    isMC = !(e.isRealData());

    bool isBeam = false;
    if (beamOrCosmic)
    {
      isBeam= beamOrCosmic->isBeam(truth);
    }
    if (!isBeam) continue;
    if (nMCParts >= MAXMCPARTS)
    {
      mf::LogError("MCParticle") << "Too many MCParticles in event to record.";
      continue;
    }

    mcPartIsBeam = isBeam;
    mcPartIsPrimary = truth->Process() == "primary";
    mcPartTrackID = truth->TrackId();
    mcPartPDG = truth->PdgCode();
    const auto greatestGrandmotherTrue = motherDaughterWalker.getGreatestGrandmother(*truth);
    if (greatestGrandmotherTrue.isNonnull())
    {
      mcPartGreatestGrandmotherTrackID = greatestGrandmotherTrue->TrackId();
      mcPartGreatestGrandmotherPDG = greatestGrandmotherTrue->PdgCode();
    }
    mcPartStartX = truth->Vx();
    mcPartStartY = truth->Vy();
    mcPartStartZ = truth->Vz();
    if (mcPartStartZ > 20) continue;
    mcPartStartT = truth->T();
    mcPartEndX = truth->EndX();
    mcPartEndY = truth->EndY();
    mcPartEndZ = truth->EndZ();
    mcPartEndT = truth->EndT();
    mcPartStartTheta = truth->Momentum().Theta();
    mcPartStartPhi = truth->Momentum().Phi();
    mcPartStartMom = 1000.*truth->Momentum().Vect().Mag();
    mcPartStartE = 1000.*truth->E();
    mcPartStartKin = 1000.*(truth->E()-truth->Mass());
    mcPartEndMom = 1000.*truth->EndMomentum().Vect().Mag();
    mcPartEndE = 1000.*truth->EndE();
    mcPartEndKin = 1000.*(truth->EndE()-truth->Mass());
    mcPartDeltaAngle = truth->Momentum().Vect().Angle(truth->EndMomentum().Vect());
    TVector3 mcPartDir;
    mcPartDir.SetMagThetaPhi(1.,mcPartStartTheta,mcPartStartPhi);

    const TVector3 particleFrontTPCPoint = lsu::mcPartStartZPlane(0,*truth);
    mcPartXFrontTPC = particleFrontTPCPoint.X();
    mcPartYFrontTPC = particleFrontTPCPoint.Y();
    if (mcPartXFrontTPC > 100 || mcPartXFrontTPC < -200 || mcPartYFrontTPC < 300 || mcPartYFrontTPC > 500) continue;
    nMCParts++;

    //Get Variables using all tracks
    nTracks = trackVec.size();
    if (nTracks > MAXTRACKS)
    {
      mf::LogError("Track") << "Too many tracks in event to record: "
          << nTracks << " > " << MAXTRACKS << ", not saving event to tree.";
      return;
    }
    for(size_t iTrack=0; iTrack < trackVec.size(); iTrack++)
    {
      ResetTreeTrackVars();
      const art::Ptr<recob::Track> track = trackVec.at(iTrack);
      const size_t nPoints = track->NumberTrajectoryPoints();

      const TVector3 trackFrontTPCPoint = lsu::trackZPlane(0,*track);
      trackXFrontTPC = trackFrontTPCPoint.X();
      trackYFrontTPC = trackFrontTPCPoint.Y();
      if (trackXFrontTPC > 100 || trackXFrontTPC < -200 || trackYFrontTPC < 300 || trackYFrontTPC > 500) continue;
      if (nPoints > 0)
      {
        trackStartX = track->Vertex().X();
        trackStartY = track->Vertex().Y();
        trackStartZ = track->Vertex().Z();
        trackStartTheta = track->VertexDirection().Theta();
        trackStartPhi = track->VertexDirection().Phi();
        if (nPoints > 1)
        {
          trackEndX = track->End().X();
          trackEndY = track->End().Y();
          trackEndZ = track->End().Z();
        }
        else
        {
          trackEndX = track->Vertex().X();
          trackEndY = track->Vertex().Y();
          trackEndZ = track->Vertex().Z();
        }
        trackLength = track->Length();
        
        TVector3 trackDir;
        trackDir.SetMagThetaPhi(1.,trackStartTheta,trackStartPhi);
        trackMCPartAngle = trackDir.Angle(mcPartDir);
      }

      if (isMC && fmHitsForTracks.isValid())
      {
          int TrackID = DEFAULTNEG;
          std::set<int> thisTrackIDSet;
          std::vector< art::Ptr<recob::Hit> > trackHits = fmHitsForTracks.at(iTrack);
          std::set<int> setOfTrackIDs = bt->GetSetOfTrackIds(trackHits);
          std::vector< art::Ptr<recob::Hit> > trackHitsU;
          std::vector< art::Ptr<recob::Hit> > trackHitsV;
          std::vector< art::Ptr<recob::Hit> > trackHitsZ;
          for (const auto & trackHit : trackHits)
          {
              if (trackHit->View() == geo::kU) trackHitsU.push_back(trackHit);
              else if (trackHit->View() == geo::kV) trackHitsV.push_back(trackHit);
              else if (trackHit->View() == geo::kZ) trackHitsZ.push_back(trackHit);
          }
          for (const auto & trackID : setOfTrackIDs)
          {
              //std::cout << "Trying trackID: "<< trackID << std::endl;
              thisTrackIDSet.clear();
              thisTrackIDSet.insert(trackID);
              float purity = bt->HitChargeCollectionPurity(thisTrackIDSet,trackHits);
              if (purity > trackTrueChargePurity)
              {
                  TrackID = trackID;
                  trackTrueChargePurity = purity;
              }
          }
          thisTrackIDSet.clear();
          if (TrackID >= 0)
          {
            trackTrueID = TrackID;
            thisTrackIDSet.insert(TrackID);
            // In LArIAT, view V is collection, U is induction.
            trackTrueChargeEfficiencyU = bt->HitChargeCollectionEfficiency(
                                    thisTrackIDSet,trackHitsU,allHitsVec,geo::kU);
            trackTrueChargeEfficiencyV = bt->HitChargeCollectionEfficiency(
                                    thisTrackIDSet,trackHitsV,allHitsVec,geo::kV);
            trackTrueChargeEfficiencyZ = bt->HitChargeCollectionEfficiency(
                                    thisTrackIDSet,trackHitsZ,allHitsVec,geo::kZ);

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
              std::cout<<"Warning: Couldn't find MCParticle for TrackID: "<<TrackID<<" trying to find: " <<abs(TrackID)<<"\n";
              for(auto mcpart: truePartVec)
              {
                  if (mcpart->TrackId() == abs(TrackID))
                  {
                      particle = mcpart;
                      break;
                  }
              }
              if (particle.isNull())
              {
                std::string message = "Couldn't find MCParticle for TrackID: ";
                message.append(std::to_string(TrackID));
                throw cet::exception("MCParticleNotFound",message);
              }
            } // if can't find MCParticle for highest Track ID

            // This is where you can do some analysis of the true particle and compare it to the reco
            trackTrueMotherID = particle->Mother();
            trackTruePdg = particle->PdgCode();
            if (beamOrCosmic)
            {
              trackTrueIsBeam = beamOrCosmic->isBeam(particle);
            }
            trackTrueKin = 1000*(particle->E()-particle->Mass());
            trackTrueEndKin = 1000*(particle->EndE()-particle->Mass());
            trackTrueTrajLen = particle->Trajectory().TotalLength();

            trackTrueStartX = particle->Vx();
            trackTrueStartY = particle->Vy();
            trackTrueStartZ = particle->Vz();
            trackTrueStartT = particle->T();

            trackTrueEndX = particle->EndX();
            trackTrueEndY = particle->EndY();
            trackTrueEndZ = particle->EndZ();
            trackTrueEndT = particle->EndT();

            const TVector3 particleFrontTPCPoint = lsu::mcPartStartZPlane(0,*particle);
            trackTrueXFrontTPC = particleFrontTPCPoint.X();
            trackTrueYFrontTPC = particleFrontTPCPoint.Y();

            const auto greatestGrandmother = motherDaughterWalker.getGreatestGrandmother(*particle);
            if (greatestGrandmother.isNonnull())
            {
              trackTrueGreatestGrandmotherTrackID = greatestGrandmother->TrackId();
              trackTrueGreatestGrandmotherPDG = greatestGrandmother->PdgCode();
            }

            //std::cout << "Truth matching info: " << std::endl;
            //std::cout << "  iTrack:                 " << iTrack << std::endl;
            //std::cout << "  ID:                     " << trackTrueID << std::endl;
            //std::cout << "  MotherID:               " << trackTrueMotherID << std::endl;
            //std::cout << "  PDG:                    " << trackTruePdg << std::endl;
            //std::cout << "  IsBeam:                 " << trackTrueIsBeam << std::endl;
            //std::cout << "  KE:                     " << trackTrueKin << std::endl;
            //std::cout << "  EndKE:                  " << trackTrueEndKin << std::endl;
            //std::cout << "  TrajLen                 " << trackTrueTrajLen << std::endl;
            //std::cout << "  ChargePurity:           " << trackTrueChargePurity << std::endl;
            //std::cout << "  ChargeEfficiencyU:      " << trackTrueChargeEfficiencyU << std::endl;
            //std::cout << "  ChargeEfficiencyV:      " << trackTrueChargeEfficiencyV << std::endl;
            //std::cout << "  ChargeEfficiencyZ:      " << trackTrueChargeEfficiencyZ << std::endl;
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
              std::cout<<"Error: Found negative TrackID: "<<TrackID<<" for Track "<<iTrack<<"\n";
          }
      } // if isMC && fmHItsForTracks.isValid

      tree->Fill();

    } // for track

  } // for true (mcPart)

  if(beamOrCosmic != NULL) delete beamOrCosmic;

} // analyze function

void lana::BeamMatchingAnalyzer::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  // Book tree
  tree = tfs->make<TTree>("tree","tree");

  tree->Branch("isMC",&isMC,"isMC/O");
  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("subRunNumber",&subRunNumber,"subRunNumber/i");
  tree->Branch("eventNumber",&eventNumber,"eventNumber/i");

  //tree->Branch("nMCParts",&nMCParts,"nMCParts/i");
  tree->Branch("mcPartIsBeam",&mcPartIsBeam,"mcPartIsBeam/O");
  tree->Branch("mcPartIsPrimary",&mcPartIsPrimary,"mcPartIsPrimary/O");
  tree->Branch("mcPartTrackID",&mcPartTrackID,"mcPartTrackID/I");
  tree->Branch("mcPartPDG",&mcPartPDG,"mcPartPDG/I");
  tree->Branch("mcPartGreatestGrandmotherTrackID",&mcPartGreatestGrandmotherTrackID,"mcPartGreatestGrandmotherTrackID/I");
  tree->Branch("mcPartGreatestGrandmotherPDG",&mcPartGreatestGrandmotherPDG,"mcPartGreatestGrandmotherPDG/I");
  tree->Branch("mcPartStartX",&mcPartStartX,"mcPartStartX/F");
  tree->Branch("mcPartStartY",&mcPartStartY,"mcPartStartY/F");
  tree->Branch("mcPartStartZ",&mcPartStartZ,"mcPartStartZ/F");
  tree->Branch("mcPartStartT",&mcPartStartT,"mcPartStartT/F");
  tree->Branch("mcPartEndX",&mcPartEndX,"mcPartEndX/F");
  tree->Branch("mcPartEndY",&mcPartEndY,"mcPartEndY/F");
  tree->Branch("mcPartEndZ",&mcPartEndZ,"mcPartEndZ/F");
  tree->Branch("mcPartEndT",&mcPartEndT,"mcPartEndT/F");
  tree->Branch("mcPartStartTheta",&mcPartStartTheta,"mcPartStartTheta/F");
  tree->Branch("mcPartStartPhi",&mcPartStartPhi,"mcPartStartPhi/F");
  tree->Branch("mcPartXFrontTPC",&mcPartXFrontTPC,"mcPartXFrontTPC/F");
  tree->Branch("mcPartYFrontTPC",&mcPartYFrontTPC,"mcPartYFrontTPC/F");
  tree->Branch("mcPartStartMom",&mcPartStartMom,"mcPartStartMom/F");
  tree->Branch("mcPartStartE",&mcPartStartE,"mcPartStartE/F");
  tree->Branch("mcPartStartKin",&mcPartStartKin,"mcPartStartKin/F");
  tree->Branch("mcPartEndMom",&mcPartEndMom,"mcPartEndMom/F");
  tree->Branch("mcPartEndE",&mcPartEndE,"mcPartEndE/F");
  tree->Branch("mcPartEndKin",&mcPartEndKin,"mcPartEndKin/F");
  tree->Branch("mcPartDeltaAngle",&mcPartDeltaAngle,"mcPartDeltaAngle/F");

  //tree->Branch("nTracks",&nTracks,"nTracks/i");
  tree->Branch("trackStartX",&trackStartX,"trackStartX/F");
  tree->Branch("trackStartY",&trackStartY,"trackStartY/F");
  tree->Branch("trackStartZ",&trackStartZ,"trackStartZ/F");
  tree->Branch("trackStartTheta",&trackStartTheta,"trackStartTheta/F");
  tree->Branch("trackStartPhi",&trackStartPhi,"trackStartPhi/F");
  tree->Branch("trackEndX",&trackEndX,"trackEndX/F");
  tree->Branch("trackEndY",&trackEndY,"trackEndY/F");
  tree->Branch("trackEndZ",&trackEndZ,"trackEndZ/F");
  tree->Branch("trackLength",&trackLength,"trackLength/F");
  tree->Branch("trackXFrontTPC",&trackXFrontTPC,"trackXFrontTPC/F");
  tree->Branch("trackYFrontTPC",&trackYFrontTPC,"trackYFrontTPC/F");
  tree->Branch("trackTrueID",&trackTrueID,"trackTrueID/I");
  tree->Branch("trackTrueMotherID",&trackTrueMotherID,"trackTrueMotherID/I");
  tree->Branch("trackTrueGreatestGrandmotherTrackID",&trackTrueGreatestGrandmotherTrackID,"trackTrueGreatestGrandmotherTrackID/I");
  tree->Branch("trackTrueGreatestGrandmotherPDG",&trackTrueGreatestGrandmotherPDG,"trackTrueGreatestGrandmotherPDG/I");
  tree->Branch("trackTruePdg",&trackTruePdg,"trackTruePdg/I");
  tree->Branch("trackTrueIsBeam",&trackTrueIsBeam,"trackTrueIsBeam/O");
  tree->Branch("trackTrueKin",&trackTrueKin,"trackTrueKin/F");
  tree->Branch("trackTrueEndKin",&trackTrueEndKin,"trackTrueEndKin/F");
  tree->Branch("trackTrueTrajLen",&trackTrueTrajLen,"trackTrueTrajLen/F");
  tree->Branch("trackTrueChargePurity",&trackTrueChargePurity,"trackTrueChargePurity/F");
  tree->Branch("trackTrueChargeEfficiencyU",&trackTrueChargeEfficiencyU,"trackTrueChargeEfficiencyU/F");
  tree->Branch("trackTrueChargeEfficiencyV",&trackTrueChargeEfficiencyV,"trackTrueChargeEfficiencyV/F");
  tree->Branch("trackTrueChargeEfficiencyZ",&trackTrueChargeEfficiencyZ,"trackTrueChargeEfficiencyZ/F");
  tree->Branch("trackTrueStartX",&trackTrueStartX,"trackTrueStartX/F");
  tree->Branch("trackTrueStartY",&trackTrueStartY,"trackTrueStartY/F");
  tree->Branch("trackTrueStartZ",&trackTrueStartZ,"trackTrueStartZ/F");
  tree->Branch("trackTrueStartT",&trackTrueStartT,"trackTrueStartT/F");
  tree->Branch("trackTrueEndX",&trackTrueEndX,"trackTrueEndX/F");
  tree->Branch("trackTrueEndY",&trackTrueEndY,"trackTrueEndY/F");
  tree->Branch("trackTrueEndZ",&trackTrueEndZ,"trackTrueEndZ/F");
  tree->Branch("trackTrueEndT",&trackTrueEndT,"trackTrueEndT/F");
  tree->Branch("trackTrueXFrontTPC",&trackTrueXFrontTPC,"trackTrueXFrontTPC/F");
  tree->Branch("trackTrueYFrontTPC",&trackTrueYFrontTPC,"trackTrueYFrontTPC/F");

  tree->Branch("trackMCPartAngle",&trackMCPartAngle,"trackMCPartAngle/F");

}

void lana::BeamMatchingAnalyzer::beginRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::endJob()
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::endRun(art::Run const & r)
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::endSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fTruePartLabel = p.get<art::InputTag>("TruePartLabel");
  fBeamTruthTag = p.get<art::InputTag>("BeamTruthTag");
  fCosmicTruthTag = p.get<art::InputTag>("CosmicTruthTag");
  fTrackLabel = p.get<art::InputTag>("TrackLabel");
  fTrackHitTag = p.get<art::InputTag>("TrackHitTag");
  fAllHitTag = p.get<art::InputTag>("AllHitTag");
}

void lana::BeamMatchingAnalyzer::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lana::BeamMatchingAnalyzer::ResetTreeVars() 
{
  isMC = false;
  runNumber = 0;
  subRunNumber = 0;
  eventNumber = 0;

  nMCParts = 0;
  mcPartIsBeam = false;
  mcPartTrackID = DEFAULTNEG;
  mcPartIsPrimary = false;
  mcPartPDG = DEFAULTNEG;
  mcPartGreatestGrandmotherTrackID = DEFAULTNEG;
  mcPartGreatestGrandmotherPDG = DEFAULTNEG;
  mcPartStartX = DEFAULTNEG;
  mcPartStartY = DEFAULTNEG;
  mcPartStartZ = DEFAULTNEG;
  mcPartStartT = DEFAULTNEG;
  mcPartEndX = DEFAULTNEG;
  mcPartEndY = DEFAULTNEG;
  mcPartEndZ = DEFAULTNEG;
  mcPartEndT = DEFAULTNEG;
  mcPartStartTheta = DEFAULTNEG;
  mcPartStartPhi = DEFAULTNEG;
  mcPartXFrontTPC = DEFAULTNEG;
  mcPartYFrontTPC = DEFAULTNEG;
  mcPartStartMom = DEFAULTNEG;
  mcPartStartE = DEFAULTNEG;
  mcPartStartKin = DEFAULTNEG;
  mcPartEndMom = DEFAULTNEG;
  mcPartEndE = DEFAULTNEG;
  mcPartEndKin = DEFAULTNEG;
  mcPartDeltaAngle = DEFAULTNEG;

  nTracks = 0;

  trackStartX = DEFAULTNEG;
  trackStartY = DEFAULTNEG;
  trackStartZ = DEFAULTNEG;
  trackStartTheta = DEFAULTNEG;
  trackStartPhi = DEFAULTNEG;
  trackEndX = DEFAULTNEG;
  trackEndY = DEFAULTNEG;
  trackEndZ = DEFAULTNEG;
  trackLength = DEFAULTNEG;
  trackXFrontTPC = DEFAULTNEG;
  trackYFrontTPC = DEFAULTNEG;
  trackTrueID = DEFAULTNEG;
  trackTrueMotherID = DEFAULTNEG;
  trackTrueGreatestGrandmotherTrackID = DEFAULTNEG;
  trackTrueGreatestGrandmotherPDG = DEFAULTNEG;
  trackTruePdg = DEFAULTNEG;
  trackTrueIsBeam = false;
  trackTrueKin = DEFAULTNEG;
  trackTrueEndKin = DEFAULTNEG;
  trackTrueTrajLen = DEFAULTNEG;
  trackTrueChargePurity = DEFAULTNEG;
  trackTrueChargeEfficiencyU = DEFAULTNEG;
  trackTrueChargeEfficiencyV = DEFAULTNEG;
  trackTrueChargeEfficiencyZ = DEFAULTNEG;
  trackTrueStartX = DEFAULTNEG;
  trackTrueStartY = DEFAULTNEG;
  trackTrueStartZ = DEFAULTNEG;
  trackTrueStartT = DEFAULTNEG;
  trackTrueEndX = DEFAULTNEG;
  trackTrueEndY = DEFAULTNEG;
  trackTrueEndZ = DEFAULTNEG;
  trackTrueEndT = DEFAULTNEG;
  trackTrueXFrontTPC = DEFAULTNEG;
  trackTrueYFrontTPC = DEFAULTNEG;
  trackMCPartAngle = DEFAULTNEG;
}

void lana::BeamMatchingAnalyzer::ResetTreeTrackVars() 
{
  nTracks = 0;

  trackStartX = DEFAULTNEG;
  trackStartY = DEFAULTNEG;
  trackStartZ = DEFAULTNEG;
  trackStartTheta = DEFAULTNEG;
  trackStartPhi = DEFAULTNEG;
  trackEndX = DEFAULTNEG;
  trackEndY = DEFAULTNEG;
  trackEndZ = DEFAULTNEG;
  trackLength = DEFAULTNEG;
  trackXFrontTPC = DEFAULTNEG;
  trackYFrontTPC = DEFAULTNEG;
  trackTrueID = DEFAULTNEG;
  trackTrueMotherID = DEFAULTNEG;
  trackTrueGreatestGrandmotherTrackID = DEFAULTNEG;
  trackTrueGreatestGrandmotherPDG = DEFAULTNEG;
  trackTruePdg = DEFAULTNEG;
  trackTrueIsBeam = false;
  trackTrueKin = DEFAULTNEG;
  trackTrueEndKin = DEFAULTNEG;
  trackTrueTrajLen = DEFAULTNEG;
  trackTrueChargePurity = DEFAULTNEG;
  trackTrueChargeEfficiencyU = DEFAULTNEG;
  trackTrueChargeEfficiencyV = DEFAULTNEG;
  trackTrueChargeEfficiencyZ = DEFAULTNEG;
  trackTrueStartX = DEFAULTNEG;
  trackTrueStartY = DEFAULTNEG;
  trackTrueStartZ = DEFAULTNEG;
  trackTrueStartT = DEFAULTNEG;
  trackTrueEndX = DEFAULTNEG;
  trackTrueEndY = DEFAULTNEG;
  trackTrueEndZ = DEFAULTNEG;
  trackTrueEndT = DEFAULTNEG;
  trackTrueXFrontTPC = DEFAULTNEG;
  trackTrueYFrontTPC = DEFAULTNEG;
  trackMCPartAngle = DEFAULTNEG;
}

DEFINE_ART_MODULE(lana::BeamMatchingAnalyzer)
