#include "dunetpc/dune/LSU/MotherDaughterWalkerAlg.h"

pdana::MotherDaughterWalkerAlg::MotherDaughterWalkerAlg(art::Event const & event, const art::InputTag mcPartTag)
{
  if(!event.isRealData())
  {
    auto allParticlesHandle = event.getValidHandle<std::vector<simb::MCParticle> >(mcPartTag);
    art::fill_ptr_vector(allParticlePtrs, allParticlesHandle);
  }
}

const art::Ptr<simb::MCParticle> pdana::MotherDaughterWalkerAlg::getGreatestGrandmother(simb::MCParticle const & particle) const
{
  art::Ptr<simb::MCParticle> result;
  if (particle.Process() == "primary") return result;
  result = getMother(particle);
  if (result->Process() == "primary") return result;
  return getGreatestGrandmother(*result);
}

const art::Ptr<simb::MCParticle> pdana::MotherDaughterWalkerAlg::getMother(simb::MCParticle const & particle) const
{
  art::Ptr<simb::MCParticle> result;
  if (particle.Process() == "primary") return result;
  const int trackId = particle.Mother();
  result = getParticle(trackId);
  return result;
}

const std::vector<art::Ptr<simb::MCParticle> > pdana::MotherDaughterWalkerAlg::getDaughters(simb::MCParticle const & particle) const
{
  std::vector<art::Ptr<simb::MCParticle>> result;
  const size_t nDaughters = particle.NumberDaughters();
  for(size_t iDaughter=0; iDaughter < nDaughters; iDaughter++)
  {
    const int daughterID = particle.Daughter(iDaughter);
    result.push_back(getParticle(daughterID));
  }
  return result;
}

const art::Ptr<simb::MCParticle> pdana::MotherDaughterWalkerAlg::getParticle(const int TrackID) const
{
  art::Ptr<simb::MCParticle> result;
  for( auto const & mcPart: allParticlePtrs)
  {
    const int partTrackID = mcPart->TrackId();
    if (TrackID == partTrackID)
    {
      result = mcPart;
      break;
    }
  }
  return result;
}

