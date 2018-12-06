#ifndef BACKTRACKERALG_H
#define BACKTRACKERALG_H

#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"

namespace lsu
{
  class BackTrackerAlg
  {
   public:
    //Give a collection of hits, get an MCParticle
    const simb::MCParticle getMCParticle(const std::vector<art::Ptr<recob::Hit> >& hits, art::Event const & e);
    const simb::MCParticle getMCParticle(const std::vector<art::Ptr<recob::Hit> >& hits, art::Event const & e, bool& isNeg);
    //Give a track, get an MCParticle
    const simb::MCParticle getMCParticle(art::Ptr<recob::Track> Track,const art::Event & e);
    const simb::MCParticle getMCParticle(art::Ptr<recob::Track> Track,const art::Event & e, bool& isNeg);
   private:
    void getHCPTrackID(std::set<int> setOfTrackIDs,int &highestChargePurityTrackID,double &highestChargePurity,art::ServiceHandle<cheat::BackTrackerService> bt,const std::vector<art::Ptr<recob::Hit> > hits);
    void finishGetting(int highestChargePurityTrackID,std::vector<art::Ptr<simb::MCParticle> > mcparticleList, art::Ptr<simb::MCParticle> &mcparticle);
  };
}
#endif
