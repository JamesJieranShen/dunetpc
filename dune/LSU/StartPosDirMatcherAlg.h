#ifndef STARTPOSDIRMATCHERALG_H
#define STARTPOSDIRMATCHERALG_H
/*!
 * Title:   Start Position & Direction Matcher Algorithim Class
 * Author:  Justin Hugon (jhugon@fnal.gov)
 *
 * Description: Algorithm for matching MCParticles to reconstructed tracks 
 *              based on the start position and direction of the track and 
 *              MCParticle.
 * Input:       simb::MCParticle, vector<art::Ptr<recob::Track>>, and 
 *              maximum angle for match (double)
 * Output:      Matched track (art::Ptr<recob::Track>) and angle in 
 *              radians between start direction of track and MCParticle 
 *              (double)
*/

//Framework includes
#include "canvas/Persistency/Common/Ptr.h"

//LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"

//c++ includes
#include <vector>
#include <tuple>

namespace mctrue {
  class StartPosDirMatcherAlg;
}

/// Matches Reco Tracks to MCParticles
/**
 * This is a class to make matching reco objects to MCParticles
 * easier. The starting position and direction of the tracks
 * and MCParticle are used for matching
 *
 */
class mctrue::StartPosDirMatcherAlg
{
  public:
    /// Find best match track given MCParticle and set of tracks
    /**
     * Find the best matching track to the given mcParticle.
     * The best match is found by finding the recoTrack with
     * one end closest to the start of the mcParticle trajectory
     * and with angle between the trajectory direction and
     * track direction at that point less than maxAngleDeg.
     * maxAngleDeg is the minimum angle between the two lines
     * i.e. always <= 90 deg.
     */
    inline const art::Ptr<recob::Track> getBestMatch(
                simb::MCParticle const& mcParticle, 
                std::vector<art::Ptr<recob::Track>> const& recoTracks,
                double const& maxAngleDeg)
    {
      double d, a;
      return getBestMatch(mcParticle,recoTracks,maxAngleDeg,d,a);
    };

    /// Find best match track given MCParticle and set of tracks, with diagnostic info
    /**
     * Find the best matching track to the given mcParticle.
     * The best match is found by finding the recoTrack with
     * one end closest to the start of the mcParticle trajectory
     * and with angle between the trajectory direction and
     * track direction at that point less than maxAngleDeg.
     * maxAngleDeg is the minimum angle between the two lines
     * i.e. always <= 90 deg.
     *
     * distance is set to the distance between the track endpoint
     * and the mcParticle start point.
     *
     * angleDeg is set to the angle between the track direction at
     * the matched endpoint and the mcParticle start direction.
     * It is always <= 90 deg.
     */
    const art::Ptr<recob::Track> getBestMatch(
                simb::MCParticle const& mcParticle, 
                std::vector<art::Ptr<recob::Track>> const& recoTracks,
                double const& maxAngleDeg, double& distance, double& angleDeg);

    // Messing with new c++17 structured binding
    // should be able to do:
    // const auto& [track, distance, angleDeg] = getBestMatch(particle,tracks,maxAngleDeg);
    inline const std::tuple<art::Ptr<recob::Track>, double, double> getBestMatchTup(
                simb::MCParticle const& mcParticle, 
                std::vector<art::Ptr<recob::Track>> const& recoTracks,
                double const& maxAngleDeg)
    {
        double distance;
        double angleDeg;
        const auto track = getBestMatch(mcParticle,recoTracks,maxAngleDeg,distance,angleDeg);
        return std::make_tuple(track,distance,angleDeg);
    };
};

#endif
