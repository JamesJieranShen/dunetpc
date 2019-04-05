#include "dune/Protodune/Analysis/ProtoDUNETrackUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TFile.h"
#include "TH1F.h"

protoana::ProtoDUNETrackUtils::ProtoDUNETrackUtils(){

}

protoana::ProtoDUNETrackUtils::~ProtoDUNETrackUtils(){

}

std::vector<anab::CosmicTag> protoana::ProtoDUNETrackUtils::GetRecoTrackCosmicTag(const recob::Track &track, art::Event const &evt, std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  unsigned int trackIndex = track.ID();

  // Convert to std::vector<anab::CosmicTag> from std::vector<art::Ptr<anab::CosmicTag>>
  std::vector<anab::CosmicTag> trackTags;

  try{
    const art::FindManyP<anab::CosmicTag> findCosmicTags(recoTracks,evt,trackModule);
    for(unsigned int t = 0; t < findCosmicTags.at(trackIndex).size(); ++t){
      trackTags.push_back((*(findCosmicTags.at(trackIndex)[t])));
    }
  }
  catch(...){
//    std::cerr << "Product not found - returning empty vector" << std::endl;
  }

  return trackTags;    
}

std::vector<anab::T0> protoana::ProtoDUNETrackUtils::GetRecoTrackT0(const recob::Track &track, art::Event const &evt, std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);

  unsigned int trackIndex = track.ID();

  // Convert to std::vector<anab::T0> from std::vector<art::Ptr<anab::T0>>
  std::vector<anab::T0> trackT0s;
  
  try{
    const art::FindManyP<anab::T0> findTrackT0s(recoTracks,evt,trackModule);
    for(unsigned int t = 0; t < findTrackT0s.at(trackIndex).size(); ++t){
      trackT0s.push_back((*(findTrackT0s.at(trackIndex)[t])));
    }
  }
  catch(...){
//    std::cerr << "Product not found - returning empty vector" << std::endl;
  }
  
  return trackT0s;

}

// Get the Calorimetry(s) from a given reco track
std::vector<anab::Calorimetry> protoana::ProtoDUNETrackUtils::GetRecoTrackCalorimetry(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  std::vector<anab::Calorimetry> caloInfo;
  
  try{
    const art::FindManyP<anab::Calorimetry> findCalorimetry(recoTracks,evt,caloModule);
    std::vector<art::Ptr<anab::Calorimetry>> theseCalos = findCalorimetry.at(track.ID());

    for( auto calo : theseCalos){
      caloInfo.push_back(*calo);
    }
  }
  catch(...){
    std::cerr << "No calorimetry object found... returning empty vector" << std::endl;
  }

  return caloInfo;
}

// Get the hits from a given reco track
const std::vector<const recob::Hit*> protoana::ProtoDUNETrackUtils::GetRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,trackModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(track.ID());

  std::vector<const recob::Hit*> trackHits;

  for(const art::Ptr<recob::Hit> hit : inputHits){

    trackHits.push_back(hit.get());

  }

  return trackHits;  

}

const std::vector<const recob::Hit*> protoana::ProtoDUNETrackUtils::GetRecoTrackHitsFromPlane(const recob::Track &track, art::Event const &evt, const std::string trackModule, unsigned int planeID ) const{

  std::vector<const recob::Hit*> trackHits;
  if( planeID > 2 ){
    std::cout << "Please input plane 0, 1, or 2" << std::endl;
    return trackHits;
  }

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  art::FindManyP<recob::Hit> findHits(recoTracks,evt,trackModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(track.ID());

  for(const art::Ptr<recob::Hit> hit : inputHits){
    unsigned int thePlane = hit.get()->WireID().asPlaneID().Plane;
    if( thePlane != planeID ) continue;
       
    trackHits.push_back(hit.get());

  }

  return trackHits;  

}

std::vector< float >  protoana::ProtoDUNETrackUtils::CalibrateCalorimetry(  const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, const fhicl::ParameterSet &ps ) {


  unsigned int planeID = ps.get< unsigned int >( "PlaneID" );
   
  double betap  = ps.get< double >( "betap"  );
  double Rho    = ps.get< double >( "Rho"    );
  double Efield = ps.get< double >( "Efield" );
  double Wion   = ps.get< double >( "Wion"   );
  double alpha  = ps.get< double >( "alpha"  );
  double norm_factor = ps.get< double >( "norm_factor" );
  double calib_factor = ps.get< double >( "calib_factor" );
  std::string X_correction_name = ps.get< std::string >( "X_correction" );
  TFile X_correction_file = TFile( X_correction_name.c_str(), "OPEN" );
  TH1F * X_correction_hist = (TH1F*)X_correction_file.Get( "dqdx_X_correction_hist" );


  std::vector< float > calibrated_dEdx;

  //Get the Calorimetry vector from the track
  std::vector< anab::Calorimetry > caloVector = GetRecoTrackCalorimetry( track, evt, trackModule, caloModule ); 
  
  size_t calo_position;
  bool found_plane = false;
  for( size_t i = 0; i < caloVector.size(); ++i ){
     unsigned int thePlane = caloVector.at(i).PlaneID().Plane;
     if( thePlane == planeID ){
       calo_position = i;
       found_plane = true;
       break;
     }
  }

  if( !found_plane ){
    std::cout << "Could not find the correct plane in the calorimetry vector" << std::endl;
    return calibrated_dEdx;
  }

  std::vector< float > dQdX = caloVector.at( calo_position).dQdx();
  auto theXYZPoints = caloVector.at( calo_position).XYZ();

  //Get the hits from the track from a specific plane
  const std::vector< const recob::Hit* > hits = GetRecoTrackHitsFromPlane( track, evt, trackModule, planeID ); 
  if( hits.size() == 0 ){
    std::cout << "Got empty hits vector" << std::endl;
    return calibrated_dEdx;
  }

  //Do Ajib's correction 
  for( size_t i = 0; i < dQdX.size(); ++i ){ 
    float hit_x = theXYZPoints[i].X();
    int X_bin = X_correction_hist->FindBin( hit_x );
    float X_correction = X_correction_hist->GetBinContent(X_bin);

    float corrected_dq_dx = dQdX[i] * X_correction * norm_factor;
    float scaled_corrected_dq_dx = corrected_dq_dx / calib_factor;
    float cal_de_dx = calc_dEdX( scaled_corrected_dq_dx,  betap,  Rho,  Efield,  Wion,  alpha );
 
    calibrated_dEdx.push_back( cal_de_dx );
  }


  return calibrated_dEdx;
}

float protoana::ProtoDUNETrackUtils::calc_dEdX(double dqdx, double betap, double Rho, double Efield, double Wion, double alpha){
  return (exp(dqdx*(betap/(Rho*Efield)*Wion))-alpha)/(betap/(Rho*Efield));  
}


// Get the hits from a given reco track
unsigned int protoana::ProtoDUNETrackUtils::GetNumberRecoTrackHits(const recob::Track &track, art::Event const &evt, const std::string trackModule) const{

  return GetRecoTrackHits(track,evt,trackModule).size();

}

// Get the PID from a given track
std::vector<anab::ParticleID> protoana::ProtoDUNETrackUtils::GetRecoTrackPID(const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string pidModule) const{

  auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
  std::vector<anab::ParticleID> pidvec;

  try{
    const art::FindManyP<anab::ParticleID> findPID(recoTracks,evt,pidModule);
    std::vector<art::Ptr<anab::ParticleID>> thePID = findPID.at(track.ID());
  
    for( auto pid : thePID){
      pidvec.push_back(*pid);
    }
  }
  catch(...){
    std::cerr << "No track PID object found... returning empty vector" << std::endl;
  }

  return pidvec;

}

protoana::BrokenTrack protoana::ProtoDUNETrackUtils::IsBrokenTrack( const recob::Track &track, art::Event const &evt, const std::string trackModule, const std::string caloModule, const fhicl::ParameterSet & BrokenTrackPars, const fhicl::ParameterSet & CalorimetryPars ){
  
  BrokenTrack theBrokenTrack;
  theBrokenTrack.Valid = false;
  
  double fBrokenTrackZ_low  = BrokenTrackPars.get< double >( "BrokenTrackZ_low" );
  double fBrokenTrackZ_high = BrokenTrackPars.get< double >( "BrokenTrackZ_high" );

  double fStitchTrackZ_low  = BrokenTrackPars.get< double >( "StitchTrackZ_low" );
  double fStitchTrackZ_high = BrokenTrackPars.get< double >( "StitchTrackZ_high" );
  
  double fStitchXTol = BrokenTrackPars.get< double >( "StitchXTol" );
  double fStitchYTol = BrokenTrackPars.get< double >( "StitchYTol" );

  ////Check the end of the track
  double endZ = track.Trajectory().End().Z();
  double endX = track.Trajectory().End().X();
  double endY = track.Trajectory().End().Y();
  if( fBrokenTrackZ_low < endZ && endZ < fBrokenTrackZ_high ){

    const auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(trackModule);
    for( auto const & tr : *recoTracks ){          

      //Skip the track in question 
      if( tr.ID() == track.ID() ) continue;

      double stitchStartZ = tr.Trajectory().Start().Z();
      if( fStitchTrackZ_low < stitchStartZ && stitchStartZ < fStitchTrackZ_high ){
        double deltaX = fabs(endX - tr.Trajectory().Start().X());
        double deltaY = fabs(endY - tr.Trajectory().Start().Y());

        std::cout << "Possible stitching track: " << stitchStartZ << " " << deltaX << " " << deltaY << std::endl;

        if( deltaX < fStitchXTol && deltaY < fStitchYTol ){


          //Get the cosine of the angle between them
          auto stitchDir = tr.StartDirection();
          double stitch_cos_theta =  stitchDir.X()*track.EndDirection().X() + stitchDir.Y()*track.EndDirection().Y() + stitchDir.Z()*track.EndDirection().Z() ;
          std::cout << "Cos_theta " << stitch_cos_theta << std::endl;


          unsigned int planeID = CalorimetryPars.get< unsigned int >( "PlaneID");
          //Get the calorimetries, calibrate, and combine
          std::vector< anab::Calorimetry > broken_calos = GetRecoTrackCalorimetry(track, evt, trackModule, caloModule);

          bool found_broken_calo = false; 
          size_t calo_position= 0;
          for( size_t i = 0; i < broken_calos.size(); ++i ){
            unsigned int thePlane = broken_calos.at(i).PlaneID().Plane;
            if( thePlane == planeID ){
              found_broken_calo = true;
              calo_position = i;
              break;
            }
          }
          
          //If no calorimetry object is found for this track, return the default BrokenTrack.
          if( !found_broken_calo ) return theBrokenTrack; 

          auto broken_range = broken_calos.at( calo_position ).ResidualRange();
          auto broken_dQdx  = broken_calos.at( calo_position ).dQdx();
          std::vector< float > broken_cal_dEdx = CalibrateCalorimetry(  track, evt, trackModule, caloModule, CalorimetryPars );

          calo_position = 0;
          std::vector< anab::Calorimetry > stitch_calos = GetRecoTrackCalorimetry(tr, evt, trackModule, caloModule);

          bool found_stitch_calo = false;
          for( size_t i = 0; i < stitch_calos.size(); ++i ){
            unsigned int thePlane = stitch_calos.at(i).PlaneID().Plane;
            if( thePlane == planeID ){
              found_stitch_calo = true;
              calo_position = i;
              break;
            }
          }

          //If no calorimetry object is found for this track, try the next 
          if( !found_stitch_calo ) continue; 

          auto stitch_range = stitch_calos.at( calo_position ).ResidualRange();
          auto stitch_dQdx  = stitch_calos.at( calo_position ).dQdx();
          std::vector< float > stitch_cal_dEdx = CalibrateCalorimetry(  tr, evt, trackModule, caloModule, CalorimetryPars );

          //piece them together in order       
          std::vector< float > combined_range, combined_dQdx, combined_dEdx;
          

          combined_range = stitch_range;
          if( stitch_range[0] > stitch_range.back() ){
            std::cout << "Adding range: " << stitch_range[0] << std::endl;
            for( size_t i = 0 ; i < broken_range.size(); ++i ){
              combined_range.push_back( broken_range[i] + stitch_range[0] );
            }
          }
          else{
            std::cout << "Adding range: " << stitch_range[0] << std::endl;
            for( size_t i = 0 ; i < broken_range.size(); ++i ){
              combined_range.push_back( broken_range[i] + stitch_range.back() );
            }
          }

          for( size_t i = 0; i < combined_range.size(); ++i ){
            std::cout << combined_range[i] << std::endl;
          }

          combined_dQdx = stitch_dQdx;
          combined_dQdx.insert( combined_dQdx.end(), broken_dQdx.begin(), broken_dQdx.end() );
          combined_dEdx = stitch_cal_dEdx;
          combined_dEdx.insert( combined_dEdx.end(), broken_cal_dEdx.begin(), broken_cal_dEdx.end() );
/*          float total_stitch_range;
          if( stitch_range[0] > stitch_range.back() ){
            total_stitch_range = stitch_range[0]; 
          }
          else{
            total_stitch_range = stitch_range.back();
          }

          std::cout << stitch_range[0] << " " << stitch_range.back() << " " << total_stitch_range << std::endl;

          if( broken_range[0] > broken_range.back() ){
            for( size_t i = 0; i < broken_range.size(); ++i ){
              combined_range.push_back( broken_range.at(i) + total_stitch_range );
              combined_dQdx.push_back(  broken_dQdx.at(i) );
              combined_dEdx.push_back(  broken_cal_dEdx.at(i) );
            }
          }
          else{
            for( size_t i = broken_range.size() - 1; i >= 0; --i ){
              combined_range.push_back( broken_range.at(i) + total_stitch_range );
              combined_dQdx.push_back(  broken_dQdx.at(i) );
              combined_dEdx.push_back(  broken_cal_dEdx.at(i) );
            }
          }

          if( stitch_range[0] > stitch_range.back() ){
            combined_range.insert( combined_range.end(), stitch_range.begin(), stitch_range.end() );
            combined_dQdx.insert( combined_dQdx.end(), stitch_dQdx.begin(), stitch_dQdx.end() );
            combined_dEdx.insert( combined_dEdx.end(), stitch_cal_dEdx.begin(), stitch_cal_dEdx.end() );
          }
          else{
            combined_range.insert( combined_range.end(), stitch_range.rbegin(), stitch_range.rend() );
            combined_dQdx.insert( combined_dQdx.end(), stitch_dQdx.rbegin(), stitch_dQdx.rend() );
            combined_dEdx.insert( combined_dEdx.end(), stitch_cal_dEdx.rbegin(), stitch_cal_dEdx.rend() );
          }
*/          

          theBrokenTrack.firstTrack = &track;
          theBrokenTrack.secondTrack = &tr;
          theBrokenTrack.CosTheta = stitch_cos_theta; 
          theBrokenTrack.Combined_ResidualRange = combined_range;
          theBrokenTrack.Combined_dQdx = combined_dQdx;
          theBrokenTrack.Combined_dEdx = combined_dEdx;
          theBrokenTrack.Valid = true;

          return theBrokenTrack;
        }
      }
    }
  }
  return theBrokenTrack;
}
