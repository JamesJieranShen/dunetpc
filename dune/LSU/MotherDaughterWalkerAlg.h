#ifndef MOTHERDAUGHTERWALKERALG_H
#define MOTHERDAUGHTERWALKERALG_H

///////////////////////////////////////////////////////
// Looks through the list of MCParticles to find the
// mothers and daughters of a given MCParticle
//
// Justin Hugon jhugon@fnal.gov  2018-09-18
///////////////////////////////////////////////////////

#include <vector>

#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/Ptr.h" 

namespace pdana {
  class MotherDaughterWalkerAlg;
}

/** 
 *
 * Find the mother and daughter MCParticles of a given MCParticle
 *
 */
class pdana::MotherDaughterWalkerAlg
{
  public:

    /** \brief Constructor
     *
     * Use in analyze method of an EDAnalyzer (or produce for a producer, etc.)
     * mcPartTag must be the input tag (usually just label string) for an
     * std::vector<simb::MCParticle>, usually largeant.
     *
     */
    MotherDaughterWalkerAlg(art::Event const & event, const art::InputTag mcPartTag);

    /** \brief Find the greatest mother or grandmother of the particle
     *
     * Looks through all of the MCParticles for the mother, it's other, and so on
     * until finding a particle that is primary. If the input particle is primary,
     * will return an emtpy pointer
     *
     */
    const art::Ptr<simb::MCParticle> getGreatestGrandmother(simb::MCParticle const & particle) const;

//    /** \brief Find all daughters and granddaughters of a particle
//     *
//     * Looks through all of the MCParticles for the daughters, their daughters, and so on.
//     * returns a vector of all of them
//     */
//    const std::vector<art::Ptr<simb::MCParticle> > getDaughtersGranddaughters(simb::MCParticle const & particle) const;


    /** \brief Find the given particles mother
     *
     * Looks through all of the MCParticles for the mother
     * If the particle is primary, will return an empty pointer
     * (the returned pointer's isNull() method will return true)
     *
     */
    const art::Ptr<simb::MCParticle> getMother(simb::MCParticle const & particle) const;

    /** \brief Find the given particles daughters
     *
     * Looks through all of the MCParticles for daughters.
     * Some pointers may be null if we couldn't find the particle for a TrackID
     *
     */
    const std::vector<art::Ptr<simb::MCParticle> > getDaughters(simb::MCParticle const & particle) const;

    /** \brief Find particle with given TrackID
     *
     * Looks through all of the MCParticles for the one with TrackID.
     * Will return an empty pointer if it can't find the particle
     *
     */
    const art::Ptr<simb::MCParticle> getParticle(const int TrackID) const;

  private:

    std::vector<art::Ptr<simb::MCParticle> > allParticlePtrs;

};

#endif
