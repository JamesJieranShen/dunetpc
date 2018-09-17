#include "dunetpc/dune/LSU/MCBeamOrCosmicAlg.h"

pdana::MCBeamOrCosmicAlg::MCBeamOrCosmicAlg(art::Event const & event, const art::InputTag mcPartTag, const art::InputTag beamTruthTag, const art::InputTag cosmicTruthTag):
  beamTruthHandle(event.getValidHandle<std::vector<simb::MCTruth> >(beamTruthTag)),
  cosmicTruthHandle(event.getValidHandle<std::vector<simb::MCTruth> >(cosmicTruthTag)),
  beamMatchedParticles(beamTruthHandle, event, mcPartTag),
  cosmicMatchedParticles(cosmicTruthHandle, event, mcPartTag)
{
  //beamTruthHandle = event.getValidHandle<std::vector<simb::MCTruth> >(beamTruthTag);
  //cosmicTruthHandle = event.getValidHandle<std::vector<simb::MCTruth> >(cosmicTruthTag);
  //beamMatchedParticles = art::FindManyP<simb::MCParticle>(beamTruthHandle, e, mcPartTag);
  //cosmicMatchedParticles = art::FindManyP<simb::MCParticle>(cosmicTruthHandle, e, mcPartTag);
  if(beamTruthHandle->size() == 0)
  {
    throw cet::exception("NoTruth","There are 0 elements in the beam truth for MCBeamOrCosmicAlg");
  }
  if(beamMatchedParticles.size() > 1)
  {
    throw cet::exception("TooManyTruths","There are more than one element in the beam truth, MCBeamOrCosmicAlg assumes 1");
  }
  if(cosmicTruthHandle->size() == 0)
  {
    throw cet::exception("NoTruth","There are 0 elements in the cosmic truth for MCBeamOrCosmicAlg");
  }
  if(cosmicMatchedParticles.size() > 1)
  {
    throw cet::exception("TooManyTruths","There are more than one element in the cosmic truth, MCBeamOrCosmicAlg assumes 1");
  }
}

bool pdana::MCBeamOrCosmicAlg::isBeam(simb::MCParticle const & particle)
{
  const int trackId = particle.TrackId();
  return isBeam(trackId);
}

bool pdana::MCBeamOrCosmicAlg::isBeam(art::Ptr<simb::MCParticle> const & particle)
{
  const int trackId = particle->TrackId();
  return isBeam(trackId);
}

bool pdana::MCBeamOrCosmicAlg::isBeam(const int trackId)
{
  std::vector<art::Ptr<simb::MCParticle>> const & beamParts = beamMatchedParticles.at(0);
  std::vector<art::Ptr<simb::MCParticle>> const & cosmicParts = cosmicMatchedParticles.at(0);
  for(const auto beamPart: beamParts)
  {
    if (beamPart->TrackId() == trackId)
    {
      return true;
    }
  }
  for(const auto cosmicPart: cosmicParts)
  {
    if (cosmicPart->TrackId() == trackId)
    {
      return false;
    }
  }
  throw cet::exception("ParticleNotInTruth","MCBeamOrCosmicAlg::isBeam can't find the particle in either truth");
}
