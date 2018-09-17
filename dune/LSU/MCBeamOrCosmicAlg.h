#ifndef MCBEAMORCOSMICALG_H
#define MCBEAMORCOSMICALG_H

///////////////////////////////////////////////////////
// Uses associations to find out if a simb::MCParticle
// came from the beam set of generator particles or 
// the cosmic set
//
// Justin Hugon jhugon@fnal.gov  2018-09-13
///////////////////////////////////////////////////////

#include <vector>

#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/FindManyP.h"

namespace pdana {
  class MCBeamOrCosmicAlg;
}

/** 
 *
 * Get whether MCParticles come from the beam or cosmic generator
 *
 */
class pdana::MCBeamOrCosmicAlg
{
  public:

    /** \brief Constructor
     *
     * Use in analyze method of an EDAnalyzer (or produce for a producer, etc.)
     * mcPartTag must be the input tag (usually just label string) for an
     * art::Assns<simb::MCTruth,simb::MCParticle,void>, usually largeant.
     * beamTruthTag and cosmicTruthTag must be input tags of the 
     * std::vector<simb::MCTruth> consisting of the generator beam and 
     * cosmic particles respectively.
     *
     */
    MCBeamOrCosmicAlg(art::Event const & event, const art::InputTag mcPartTag, 
                    const art::InputTag beamTruthTag, const art::InputTag cosmicTruthTag);

    /** \brief Find if MCParticle is from a beam or cosmic
     *
     * Returns true if this MCParticle comes from the beam 
     * truth, and false if it came from a cosmic truth.
     *
     */
    bool isBeam(simb::MCParticle const & particle);

    bool isBeam(art::Ptr<simb::MCParticle> const & particle);


    /** \brief Find if trackId is from a beam or cosmic
     *
     * Uses trackId (given by MCParticle::TrackId(),
     * sim::IDE::trackID, or various methods of
     * cheat::BackTrackerService. Returns true if
     * the trackId corresponds to a beam particle,
     * and false if the trackId corresponds to a
     * cosmic particle.
     *
     */
    bool isBeam(const int trackId);

  private:

    art::ValidHandle<std::vector<simb::MCTruth> > beamTruthHandle;
    art::ValidHandle<std::vector<simb::MCTruth> > cosmicTruthHandle;
    art::FindManyP<simb::MCParticle> beamMatchedParticles;
    art::FindManyP<simb::MCParticle> cosmicMatchedParticles;

};

#endif
