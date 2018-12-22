#include "dune/Protodune/Analysis/ProtoDUNECaloUtils.h"

#include "canvas/Persistency/Common/FindManyP.h"

//#include "larsim/MCCheater/BackTrackerService.h"

const std::tuple<art::Ptr<anab::Calorimetry>,std::vector<recob::Track::TrajectoryPoint_t>,std::vector<art::Ptr<recob::Hit>>> protoana::ProtoDUNECaloUtils::TrackCalosWithTrajAndHit(const recob::Track& track, const art::Event& evt, const art::InputTag& trackTag, const art::InputTag& caloTag) const
{
  const auto& allTracks = evt.getValidHandle<std::vector<recob::Track>>(trackTag);
  const auto& fmCalo = art::FindManyP<anab::Calorimetry>(allTracks, evt, caloTag);
  const auto& fmHit = art::FindManyP<recob::Hit>(allTracks, evt, trackTag);
  const auto& caloVec = fmCalo.at(track.ID());
  const auto& hitVec = fmHit.at(track.ID());

  art::Ptr<anab::Calorimetry> resultCalo;
  std::vector<recob::Track::TrajectoryPoint_t> resultTPs;
  std::vector<art::Ptr<recob::Hit>> resultHits;
  for(const auto& calo: caloVec)
  {
    if(calo->PlaneID().Plane == 2)
    {
      resultCalo = calo;
      for(const auto& tpIndex: calo->TpIndices())
      {
        resultTPs.push_back(track.TrajectoryPoint(tpIndex));
        resultHits.push_back(hitVec.at(tpIndex));
      } // for tpIndex
    } // if Plane is 2
  } // for calo
  return std::make_tuple(resultCalo,resultTPs,resultHits);
}

const std::tuple<art::Ptr<anab::Calorimetry>,std::vector<recob::Track::TrajectoryPoint_t>,std::vector<art::Ptr<recob::Hit>>,std::vector<art::Ptr<sim::SimChannel>>> protoana::ProtoDUNECaloUtils::TrackCalosWithTrajHitAndChan(const recob::Track& track, const art::Event& evt, const art::InputTag& trackTag, const art::InputTag& caloTag, const art::InputTag& simTag) const
{
  const auto& allTracks = evt.getValidHandle<std::vector<recob::Track>>(trackTag);
  const auto& fmCalo = art::FindManyP<anab::Calorimetry>(allTracks, evt, caloTag);
  const auto& fmHit = art::FindManyP<recob::Hit>(allTracks, evt, trackTag);
  const auto& caloVec = fmCalo.at(track.ID());
  const auto& hitVec = fmHit.at(track.ID());

  std::vector<art::Ptr<sim::SimChannel>> allSimChannels;
  const auto& simChanHand = evt.getValidHandle<std::vector<sim::SimChannel>>(simTag);
  art::fill_ptr_vector(allSimChannels, simChanHand);

  art::Ptr<anab::Calorimetry> resultCalo;
  std::vector<recob::Track::TrajectoryPoint_t> resultTPs;
  std::vector<art::Ptr<recob::Hit>> resultHits;
  std::vector<art::Ptr<sim::SimChannel>> resultSimChannels;
  for(const auto& calo: caloVec)
  {
    if(calo->PlaneID().Plane == 2)
    {
      resultCalo = calo;
      for(const auto& tpIndex: calo->TpIndices())
      {
        resultTPs.push_back(track.TrajectoryPoint(tpIndex));
        const auto& hit = hitVec.at(tpIndex);
        const auto& hitChan = hit->Channel();
        resultHits.push_back(hit);
        if (!evt.isRealData())
        {
          for(const auto& simChan: allSimChannels)
          {
            if(simChan->Channel() == hitChan)
            {
              resultSimChannels.push_back(simChan);
            } // if simChan->Channel == hitChan
          } // for simChan
        } // if !evt.isRealData()
      } // for tpIndex
    } // if Plane is 2
  } // for calo
  return std::make_tuple(resultCalo,resultTPs,resultHits,resultSimChannels);
}
