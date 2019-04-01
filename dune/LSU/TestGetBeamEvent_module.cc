#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft includes
#include "dunetpc/dune/DuneObj/ProtoDUNEBeamEvent.h"

namespace lana {
  class TestGetBeamEvent;
  struct PiAbsSecondary;
  //struct TreeSecondary;
}

class lana::TestGetBeamEvent : public art::EDAnalyzer {
public:
  explicit TestGetBeamEvent(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TestGetBeamEvent(TestGetBeamEvent const &) = delete;
  TestGetBeamEvent(TestGetBeamEvent &&) = delete;
  TestGetBeamEvent & operator = (TestGetBeamEvent const &) = delete;
  TestGetBeamEvent & operator = (TestGetBeamEvent &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p);

private:

  // Declare member data here.

  //parameters read from fcl file
  art::InputTag fBeamEventTag;
  art::InputTag fBeamEventTagOld;
  bool fIsFirstEvent;
};

lana::TestGetBeamEvent::TestGetBeamEvent(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), fIsFirstEvent(true)
 // More initializers here.
{
  this->reconfigure(p);
}

void lana::TestGetBeamEvent::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  std::cout << "TestGetBeamEvent: " << e.id() << std::endl;

  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
  if(e.isRealData())
  {
    auto beamHand = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTag);
    if (fIsFirstEvent)
    {
      const auto prov = beamHand.provenance();
      const auto& provTag = prov->inputTag();
      std::cout << "TestGetBeamEvent: BeamEvent provenance inputTag: '" << provTag.encode() << "'\n";
      //prov->write(std::cout);
      //std::cout << "\nBeamEvent provenance done\n";
    }
    //if(beamHand.isValid())
    //{
    //  art::fill_ptr_vector(beamVec, beamHand);
    //}
  }

  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVecOld;
  if(e.isRealData())
  {
    auto beamHandOld = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamEventTagOld);
    if (fIsFirstEvent)
    {
      const auto prov = beamHandOld.provenance();
      const auto& provTag = prov->inputTag();
      std::cout << "TestGetBeamEvent: BeamEvent Old provenance inputTag: '" << provTag.encode() << "'\n";
      //prov->write(std::cout);
      //std::cout << "\nBeamEvent Old provenance done\n";
    }
    //if(beamHandOld.isValid())
    //{
    //  art::fill_ptr_vector(beamVecOld, beamHandOld);
    //}
  }

  if(fIsFirstEvent) fIsFirstEvent = false;
} // analyze function

void lana::TestGetBeamEvent::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fBeamEventTag = p.get<art::InputTag>("BeamEventTag");
  fBeamEventTagOld = p.get<art::InputTag>("BeamEventTagOld");
}

DEFINE_ART_MODULE(lana::TestGetBeamEvent)
