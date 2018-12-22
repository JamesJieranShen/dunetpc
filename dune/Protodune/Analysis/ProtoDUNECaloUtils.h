#ifndef PROTODUNE_CALO_UTILS_H
#define PROTODUNE_CALO_UTILS_H

///////////////////////////////////////////////////////////////////
// ProtoDUNECaloUtils
//  - Class to help analysers access useful calorimetry information
// 
// Justin Hugon jhugon@fnal.gov
///////////////////////////////////////////////////////////////////

#include <tuple>
#include <vector>
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "art/Framework/Principal/Event.h"

namespace protoana {

  class ProtoDUNECaloUtils {

  public:

    /// Get the collection plane calo and associated trajectory points and hits for a track
    const std::tuple<art::Ptr<anab::Calorimetry>,std::vector<recob::Track::TrajectoryPoint_t>,std::vector<art::Ptr<recob::Hit>>> TrackCalosWithTrajAndHit(const recob::Track& track, const art::Event& evt, const art::InputTag& trackTag, const art::InputTag& caloTag) const;

    /// Same as above, but get the associated sim::SimChannel, too
    const std::tuple<art::Ptr<anab::Calorimetry>,std::vector<recob::Track::TrajectoryPoint_t>,std::vector<art::Ptr<recob::Hit>>,std::vector<art::Ptr<sim::SimChannel>>> TrackCalosWithTrajHitAndChan(const recob::Track& track, const art::Event& evt, const art::InputTag& trackTag, const art::InputTag& caloTag, const art::InputTag& simTag) const;

  };

}

#endif

