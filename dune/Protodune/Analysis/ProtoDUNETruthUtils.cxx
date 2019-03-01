#include "dune/Protodune/Analysis/ProtoDUNETruthUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "art/Framework/Principal/Event.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

#include <algorithm>

#include "TVector3.h"

protoana::ProtoDUNETruthUtils::ProtoDUNETruthUtils(){

}

protoana::ProtoDUNETruthUtils::~ProtoDUNETruthUtils(){

}

// Function to find the best matched true particle to a reconstructed particle. In case of problems, returns a null pointer
const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetMCParticleFromRecoTrack(const recob::Track &track, art::Event const & evt, std::string trackModule) const{

  const simb::MCParticle* mcParticle = 0x0;

  // We must have MC for this module to make sense
  if(evt.isRealData()) return mcParticle;

  // Get the reconstructed tracks
  auto allRecoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  // We need the association between the tracks and the hits
  const art::FindManyP<recob::Hit> findTrackHits(allRecoTracks, evt, trackModule);

  unsigned int trackIndex = track.ID();

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<int, double> trkIDE;
  for (auto const & h : findTrackHits.at(trackIndex))
  {
    for (auto const & ide : bt_serv->HitToTrackIDEs(h)) // loop over std::vector<sim::TrackIDE>
    {
        trkIDE[ide.trackID] += ide.energy; // sum energy contribution by each track ID
    }
  }

  int best_id = 0;
  double tot_e = 0, max_e = 0;
  for (auto const & contrib : trkIDE)
  {
    tot_e += contrib.second;     // sum total energy in these hits
    if (contrib.second > max_e)  // find track ID corresponding to max energy
    {
        max_e = contrib.second;
        best_id = contrib.first;
    }
  }

  if ((max_e > 0) && (tot_e > 0)) // ok, found something reasonable
  {
    if (best_id < 0)            // NOTE: negative ID means this is EM activity
    {                           // caused by track with the same but positive ID
//        best_id = -best_id;     // --> we'll find mother MCParticle of these hits
      return mcParticle;
    }
    mcParticle = pi_serv->TrackIdToParticle_P(best_id); // MCParticle corresponding to track ID
  }

  return mcParticle;
}

const simb::MCParticle* protoana::ProtoDUNETruthUtils::MatchPduneMCtoG4( const simb::MCParticle & pDunePart, const art::Event & evt )
{  // Function that will match the protoDUNE MC particle to the Geant 4 MC particle, and return the matched particle (or a null pointer).

  // Find the energy of the procided MC Particle
  double pDuneEnergy = pDunePart.E();
  
  // Get list of the g4 particles. plist should be a std::map< int, simb::MCParticle* >
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  const sim::ParticleList & plist = pi_serv->ParticleList();
  
  // Check if plist is empty
  if ( !plist.size() ) {
    std::cerr << "\n\n#####################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "sim::ParticleList from cheat::ParticleInventoryService is empty\n"
              << "A null pointer will be returned\n"
              << "#####################################\n\n";
    return nullptr;
  }
  
  // Go through the list of G4 particles
  for ( auto partIt = plist.begin() ; partIt != plist.end() ; partIt++ ) {
    
    const simb::MCParticle* pPart = partIt->second;
    if (!pPart) {
      std::cerr << "\n\n#####################################\n"
                << "\nEvent " << evt.id().event() << "\n"
                << "GEANT particle #" << partIt->first << " returned a null pointer\n"
                << "This is not necessarily bad. It just means at least one\n"
                << "of the G4 particles returned a null pointer. It may well\n"
                << "have still matched a PD particle and a G4 particle.\n"
                << "#####################################\n\n";
      continue;
    }
    
    // If the initial energy of the g4 particle is very close to the energy of the protoDUNE particle, call it a day and have a cuppa.
    if ( (pDunePart.PdgCode() == pPart->PdgCode()) && fabs(pPart->E() - pDuneEnergy) < 0.00001 ) {
      return pPart;
    }
    
  }  // G4 particle list loop end.
  
  std::cout << "No G4 particle was matched for Event " << evt.id().event() << ". Null pointer returned\n";
  return nullptr;
  
}  // End MatchPduneMCtoG4

const simb::MCParticle* protoana::ProtoDUNETruthUtils::GetGeantGoodParticle(const simb::MCTruth &genTruth, const art::Event &evt) const{

  // Get the good particle from the MCTruth
  simb::MCParticle goodPart;
  bool found = false;
  for(int t = 0; t < genTruth.NParticles(); ++t){
    simb::MCParticle part = genTruth.GetParticle(t);
    if(part.Process() == "primary"){
      goodPart = part;
      found = true;
      break;
    }
  }

  if(!found){
    std::cerr << "No good particle found, returning null pointer" << std::endl;
    return nullptr;
  }

  // Now loop over geant particles to find the one that matches
  // Get list of the g4 particles. plist should be a std::map< int, simb::MCParticle* >
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;
  const sim::ParticleList & plist = pi_serv->ParticleList();

  for(auto const part : plist){
    if((goodPart.PdgCode() == part.second->PdgCode()) && fabs(part.second->E() - goodPart.E()) < 1e-5){
      return part.second;
    }
  } 

  // If we get here something has gone wrong
  std::cerr << "No G4 version of the good particle was found, returning null pointer" << std::endl;
  return nullptr;

}

// Converting times in LArSoft can be a bit of a minefield. These functions convert true times in ns
// to pandora times in ns
const float protoana::ProtoDUNETruthUtils::ConvertTrueTimeToPandoraTimeNano(const simb::MCParticle &part) const{
  return ConvertTrueTimeToPandoraTimeNano(part.T());
}

const float protoana::ProtoDUNETruthUtils::ConvertTrueTimeToPandoraTimeNano(const float trueTime) const{
  return 1000. * ConvertTrueTimeToPandoraTimeMicro(trueTime);
}

// Microsecond versions
const float protoana::ProtoDUNETruthUtils::ConvertTrueTimeToPandoraTimeMicro(const simb::MCParticle &part) const{
  return ConvertTrueTimeToPandoraTimeMicro(part.T());
}

const float protoana::ProtoDUNETruthUtils::ConvertTrueTimeToPandoraTimeMicro(const float trueTime) const{

  // Use the clocks service to account for the offset between the Geant4 time and the electronics clock
  auto const* detclock = lar::providerFrom<detinfo::DetectorClocksService>();

  return detclock->G4ToElecTime(trueTime);
}

std::vector<std::tuple<raw::ChannelID_t,unsigned short, sim::IDE> > protoana::ProtoDUNETruthUtils::GetIDEsFromParticle(const simb::MCParticle & part, const art::Event & evt, const bool includeNegativeTrackID) const{

  art::ServiceHandle<geo::Geometry> geom;
  std::vector<std::tuple<raw::ChannelID_t,unsigned short, sim::IDE> > result;
  const auto& partTrackID = part.TrackId();

  std::vector<art::Ptr<sim::SimChannel>> simChanVec;
  auto simChanHand = evt.getValidHandle<std::vector<sim::SimChannel> >("largeant");
  art::fill_ptr_vector(simChanVec, simChanHand);
  for(const auto& channel: simChanVec)
  {
    const auto& channelNumber = channel->Channel();
    auto signalType = geom->SignalType(channelNumber);
    if(signalType != geo::kCollection) continue;
    for(const auto& TDCIDEs: channel->TDCIDEMap())
    {
      unsigned short TDC = TDCIDEs.first;
      const std::vector<sim::IDE>& IDEs = TDCIDEs.second;
      for(const auto& IDE: IDEs)
      {
        if(IDE.trackID == partTrackID)
        {
          result.push_back(std::make_tuple(channelNumber,TDC,IDE));
        }
        else if (includeNegativeTrackID && (-IDE.trackID == partTrackID))
        {
          result.push_back(std::make_tuple(channelNumber,TDC,IDE));
        }
      }
    } // for TDCIDE
  } // for channel
  return result;
}

std::vector<std::tuple<raw::ChannelID_t,unsigned short, sim::IDE> > protoana::ProtoDUNETruthUtils::GetIDEsFromParticleSortZ(const simb::MCParticle & part, const art::Event & evt, const bool includeNegativeTrackID) const{
  auto result = GetIDEsFromParticle(part,evt,includeNegativeTrackID);
  std::sort(result.begin(),result.end(),[](const auto& a, const auto& b){return std::get<2>(a).z < std::get<2>(b).z;});
  return result;
}

std::vector<std::tuple<raw::ChannelID_t,float, sim::IDE> > protoana::ProtoDUNETruthUtils::GetTotalChannelEnergyFromParticle(const simb::MCParticle & part, const art::Event & evt, const bool includeNegativeTrackID) const{

  const auto& ideinfos = GetIDEsFromParticle(part,evt,includeNegativeTrackID);
  std::map<raw::ChannelID_t,std::vector<std::pair<unsigned short, sim::IDE> > > ideChannelMap;
  for (const auto& ideinfo: ideinfos) 
  {
    const auto& [channel, tdc, ide] = ideinfo;
    if (ideChannelMap.count(channel) == 0)
    {
      ideChannelMap.emplace(channel,std::vector<std::pair<unsigned short, sim::IDE> >());
    }
    ideChannelMap.at(channel).push_back(std::make_pair(tdc,ide));
  }
  std::vector<std::tuple<raw::ChannelID_t,float, sim::IDE> >  result;
  for(const auto& ideChannel: ideChannelMap)
  {
    const auto& channel = ideChannel.first;
    float xAvg = 0;
    float yAvg = 0;
    float zAvg = 0;
    float tdcAvg = 0;
    float energy = 0;
    float numElectrons = 0;
    for(const auto& tdcide: ideChannel.second)
    {
      tdcAvg += tdcide.first;
      const auto& tdc = tdcide.second;
      xAvg += tdc.x;
      yAvg += tdc.y;
      zAvg += tdc.z;
      energy += tdc.energy;
      numElectrons += tdc.numElectrons;
    }
    const size_t nIDEs = ideChannel.second.size();
    xAvg /= nIDEs;
    yAvg /= nIDEs;
    zAvg /= nIDEs;
    tdcAvg /= nIDEs;
    sim::IDE avgIDE(-1,numElectrons,energy,xAvg,yAvg,zAvg);
    result.push_back(std::make_tuple(channel,tdcAvg,avgIDE));
  }
  std::sort(result.begin(),result.end(),[](const auto& a, const auto& b){return std::get<2>(a).z < std::get<2>(b).z;});
  return result;
}

// Get process key.
int protoana::ProtoDUNETruthUtils::GetProcessKey(std::string process){

  if(process.compare("primary") == 0)                    return 0;
  else if(process.compare("hadElastic") == 0)            return 1;
  else if(process.compare("pi-Inelastic") == 0)          return 2;
  else if(process.compare("pi+Inelastic") == 0)          return 3;
  else if(process.compare("kaon-Inelastic") == 0)        return 4;
  else if(process.compare("kaon+Inelastic") == 0)        return 5;
  else if(process.compare("protonInelastic") == 0)       return 6;
  else if(process.compare("neutronInelastic") == 0)      return 7;
  else if(process.compare("kaon0SInelastic") == 0)       return 8;
  else if(process.compare("kaon0LInelastic") == 0)       return 9;
  else if(process.compare("lambdaInelastic") == 0)       return 10;
  else if(process.compare("omega-Inelastic") == 0)       return 11;
  else if(process.compare("sigma+Inelastic") == 0)       return 12;
  else if(process.compare("sigma-Inelastic") == 0)       return 13;
  else if(process.compare("sigma0Inelastic") == 0)       return 14;
  else if(process.compare("xi-Inelastic") == 0)          return 15;
  else if(process.compare("xi0Inelastic") == 0)          return 16;
  else if(process.compare("anti_protonInelastic") == 0)  return 20;
  else if(process.compare("anti_neutronInelastic") == 0) return 21;
  else if(process.compare("anti_lambdaInelastic") == 0)  return 22;
  else if(process.compare("anti_omega-Inelastic") == 0)  return 23;
  else if(process.compare("anti_sigma+Inelastic") == 0)  return 24;
  else if(process.compare("anti_sigma-Inelastic") == 0)  return 25;
  else if(process.compare("anti_xi-Inelastic") == 0)     return 26;
  else if(process.compare("anti_xi0Inelastic") == 0)     return 27;

  else if(process.compare("Decay") == 0)                 return 30;
  else if(process.compare("FastScintillation") == 0)     return 31;
  else if(process.compare("nKiller") == 0)               return 32; // Remove unwanted neutrons: neutron kinetic energy threshold (default 0) or time limit for neutron track
  else if(process.compare("nCapture") == 0)              return 33; // Neutron capture

  else if(process.compare("compt") == 0)                 return 40; // Compton Scattering
  else if(process.compare("rayleigh") == 0)              return 41; // Rayleigh Scattering
  else if(process.compare("phot") == 0)                  return 42; // Photoelectric Effect
  else if(process.compare("conv") == 0)                  return 43; // Pair production
  else if(process.compare("CoupledTransportation") == 0) return 44; //
  
  else return -1;
}

// Get estimated particle energy deposit. The G4 trackId must be provided 
double protoana::ProtoDUNETruthUtils::GetDepEnergyMC(const art::Event &evt, geo::GeometryCore const * fGeom, int trackid, int whichview) const {

  double edep = 0.0;

  art::Handle< std::vector<sim::SimChannel> > simchannelHandle;
  if(evt.getByLabel("largeant", simchannelHandle)){
    // Loop over sim channels
    for(auto const& simchannel : (*simchannelHandle)){
      // Only keep channels in the selected view
      if(fGeom->View(simchannel.Channel()) != whichview) continue;
      // Get all time slices
      auto const& alltimeslices = simchannel.TDCIDEMap();
      // Loop over time slices
      for(auto const& tslice : alltimeslices){
	auto const& simide = tslice.second;
	// Loop over energy deposits
	for(auto const& eDep : simide){
	  if(eDep.trackID != trackid) continue;
	  edep += eDep.energy;
	}
      }
    }
  }
  
  return edep;
  
}

// Get first trajectory point in TPC active volume
int protoana::ProtoDUNETruthUtils::GetFirstTrajectoryPointInTPCActiveVolume(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh){

  int firstpoint = -999;
  for(unsigned int i = 0; i < mcpart.NumberTrajectoryPoints(); ++i) {
    if(mcpart.Vx(i) >= tpcactiveXlow && mcpart.Vx(i) <= tpcactiveXhigh && mcpart.Vy(i) >= tpcactiveYlow && mcpart.Vy(i) <= tpcactiveYhigh && mcpart.Vz(i) >= tpcactiveZlow && mcpart.Vz(i) <= tpcactiveZhigh){

      firstpoint = i;
      break;
    }
  }

  return firstpoint;
}

// Get MC Particle length in TPC active volume
double protoana::ProtoDUNETruthUtils::GetMCParticleLengthInTPCActiveVolume(const simb::MCParticle& mcpart, double tpcactiveXlow, double tpcactiveXhigh, double tpcactiveYlow, double tpcactiveYhigh, double tpcactiveZlow, double tpcactiveZhigh){

  double length = 0.0;
  int firstpoint = GetFirstTrajectoryPointInTPCActiveVolume(mcpart, tpcactiveXlow, tpcactiveXhigh, tpcactiveYlow, tpcactiveYhigh, tpcactiveZlow, tpcactiveZhigh);

  if(firstpoint < 0) return length;

  TVector3 pos =  mcpart.Position(firstpoint).Vect();
  for(unsigned int i = firstpoint+1; i < mcpart.NumberTrajectoryPoints(); ++i) {
    if(mcpart.Vx(i) >= tpcactiveXlow && mcpart.Vx(i) <= tpcactiveXhigh && mcpart.Vy(i) >= tpcactiveYlow && mcpart.Vy(i) <= tpcactiveYhigh && mcpart.Vz(i) >= tpcactiveZlow && mcpart.Vz(i) <= tpcactiveZhigh){
      
      pos -= mcpart.Position(i).Vect();
      length += pos.Mag();
      pos = mcpart.Position(i).Vect();
    }
  }

  return length;
}
