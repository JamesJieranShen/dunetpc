
// Saves single electron pointing resolution to a TTree.
// Author: James Shen

#ifndef ElectronInfo_H
#define ElectronInfo_H 1

// Framework includes

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Provenance/EventID.h"

// LArSoft includes

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ROOT includes

#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TTree.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TEfficiency.h"
#include "TF2.h"
#include "Math/Functor.h"
#include "TPolyLine3D.h"
#include "Math/Vector3D.h"
#include "TFile.h"

// C++ includes

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

namespace {

	class ElectronInfo : public art::EDAnalyzer {

		public:

			// Functions that art/LArSoft need
			ElectronInfo(fhicl::ParameterSet const&);
			void analyze(art::Event const&) override;
			void endJob() override;

		private:

			// Parameters we'll read from the fcl-file
			std::string fSimulationLabel;
			std::string fTrackModuleLabel;
			std::string fHitsModuleLabel;
			std::string fOpFlashLabel;

			// Tree and tree variables
			TTree * tr;

			double elec_true_x;
			double elec_true_y;
			double elec_true_z;
			double elec_true_en;

			double elec_x_ndf;
			double elec_y_ndf;
			double elec_z_ndf;
			double cosAngle_ndf;
			double elec_x_df;
			double elec_y_df;
			double elec_z_df;
			double cosAngle_df;
			double elec_reco_en_noDC;
			double elec_reco_en_DC;
			double max_flash_time;
			double trk_charge_U;
			double trk_charge_V;
			double trk_charge_Z;
        
            double trk_charge_noDC;
            double trk_charge_DC;

			// Counters
			int has_tracks;
			int opflashcounter;
			int longest_is_prim_largest_en_frac;
			int longest_has_prim_first_pt;
			double maxmindist;

			// String for finding primary process
			std::string prim;

	};

}

#endif 

namespace {

	DEFINE_ART_MODULE(ElectronInfo)

}

namespace {

	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	// Constructor
	ElectronInfo::ElectronInfo(fhicl::ParameterSet const& parameterSet)
		: 
			EDAnalyzer(parameterSet)
	{

		// Read the fcl-file
		fSimulationLabel = parameterSet.get< std::string >("SimulationLabel");
		fTrackModuleLabel = parameterSet.get< std::string >("TrackModuleLabel");
		fHitsModuleLabel = parameterSet.get< std::string >("HitsModuleLabel");
		fOpFlashLabel = parameterSet.get< std::string >("OpFlashLabel");

		// Create ServiceHandle for TFileService output
		art::ServiceHandle<art::TFileService> tfs;

		// Make tree and branches
		tr = tfs->make<TTree>("tr", "tr");
		
		tr->Branch("elec_true_x",&elec_true_x,"elec_true_x/D");
		tr->Branch("elec_true_y",&elec_true_y,"elec_true_y/D");
		tr->Branch("elec_true_z",&elec_true_z,"elec_true_z/D");
		
		tr->Branch("elec_x_ndf",&elec_x_ndf,"elec_x_ndf/D");
		tr->Branch("elec_y_ndf",&elec_y_ndf,"elec_y_ndf/D");
		tr->Branch("elec_z_ndf",&elec_z_ndf,"elec_z_ndf/D");
		tr->Branch("cosAngle_ndf", &cosAngle_ndf, "cosAngle_ndf/D");
		
		tr->Branch("elec_x_df",&elec_x_df,"elec_x_df/D");
		tr->Branch("elec_y_df",&elec_y_df,"elec_y_df/D");
		tr->Branch("elec_z_df",&elec_z_df,"elec_z_df/D");
		tr->Branch("cosAngle_df", &cosAngle_df, "cosAngle_df/D");
		
		tr->Branch("elec_reco_en_noDC",&elec_reco_en_noDC,"elec_reco_en_noDC/D");
		tr->Branch("elec_reco_en_DC",&elec_reco_en_DC,"elec_reco_en_DC/D");
		tr->Branch("elec_true_en",&elec_true_en,"elec_true_en/D");
		tr->Branch("max_flash_time",&max_flash_time,"max_flash_time/D");
		tr->Branch("trk_charge_U",&trk_charge_U,"trk_charge_U/D");
		tr->Branch("trk_charge_V",&trk_charge_V,"trk_charge_V/D");
		tr->Branch("trk_charge_Z",&trk_charge_Z,"trk_charge_Z/D");
                
        tr->Branch("trk_charge_noDC", &trk_charge_noDC, "trk_charge_noDC/D");
        tr->Branch("trk_charge_DC", &trk_charge_DC, "trk_charge_DC/D");
		//tr->Branch("",&,"/");

		// Counters
		has_tracks = 0;
		opflashcounter = 0;
		longest_is_prim_largest_en_frac = 0;
		longest_has_prim_first_pt = 0;
		maxmindist = 0.0;

		prim = "primary";

	}

	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	void ElectronInfo::analyze(art::Event const& event)
	{

		// Initialize tree variables to nonsense values to make events
		// with no reco'ed track/flashes obvious and possible to cut and
		// initialize trk_charge to 0
		elec_true_x = -10.0;
		elec_true_y = -10.0;
		elec_true_z = -10.0;
		
		elec_x_ndf = -10.0;
		elec_y_ndf = -10.0;
		elec_z_ndf = -10.0;
		cosAngle_ndf= 2;


		elec_x_df = -10.0;
		elec_y_df = -10.0;
		elec_z_df = -10.0;
		cosAngle_df = 2;
		
		elec_reco_en_noDC = -10.0;
		elec_reco_en_DC = -10.0;
		elec_true_en = -10.0;
		max_flash_time = -3000.0;
		trk_charge_U = 0.0;
		trk_charge_V = 0.0;
		trk_charge_Z = 0.0;


		//-----------------------------------------------------------------------------------------------------
		// True neutrino direction and energy and electron energy
		//-----------------------------------------------------------------------------------------------------
	
		TVector3 elec_true_dir;
		if(fSimulationLabel != "largeant"){
			std::cout << "Simulation Label is not largeant but is " << 
			fSimulationLabel << ". Error!" <<std::endl;
		}
		if(fSimulationLabel == "largeant") {
		
			auto particleHandle
				= event.getValidHandle<std::vector<simb::MCParticle>> (fSimulationLabel);
				// electron
			for (auto const& particle : (*particleHandle)){
				if ( particle.Process() == "primary" && particle.PdgCode() == 11 ) {
	
					elec_true_x = particle.Px()/particle.P();
					elec_true_y = particle.Py()/particle.P();
					elec_true_z = particle.Pz()/particle.P();
					elec_true_en = 1000*(particle.E() - particle.Mass());
	
				}
	
			}
			elec_true_dir.SetXYZ(elec_true_x, elec_true_y, elec_true_z);
	
		}//  end truth writing


		//-----------------------------------------------------------------------------------------------------
		// Finding longest track
		//-----------------------------------------------------------------------------------------------------

		auto tracklistHandle
			= event.getValidHandle<std::vector<recob::Track>> (fTrackModuleLabel);

		int numtracks = tracklistHandle->size();
		if( numtracks != 0 ) ++has_tracks;
		//std::cout << "numtracks = " << numtracks << std::endl;

		int longest_track = -1;    
		double longest_length = 0;

		for (int i=0; i<numtracks; ++i) {

			if ( tracklistHandle->at(i).Length() > longest_length ) {

				longest_track = i;
				longest_length = tracklistHandle->at(i).Length();

			}

		}

		//-----------------------------------------------------------------------------------------------------
		// How often is longest track the primary electron? 
		//-----------------------------------------------------------------------------------------------------
		/*
		// Create BackTrackerService and ParticleInventoryService objects
		art::ServiceHandle<cheat::BackTrackerService> bt;
		art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

		// Get hits from longest track (bc there's no way to directly get truth track info from track)
		art::FindManyP<recob::Hit> fmth(tracklistHandle, event, fTrackModuleLabel);

		if(fmth.isValid() && numtracks!=0) {

		double prim_energy;
		double max_en = 0.0;
		int max_en_trk = -1;
		double min_dist = 9999.0;
		int min_dist_trk = -1;

		TVector3 start_pos_truth;

		if(fSimulationLabel == "largeant") {

		auto particleHandle
		= event.getValidHandle<std::vector<simb::MCParticle>> (fSimulationLabel);

		// Find primary electron and save its true starting position

		for (auto const& particle : (*particleHandle)) {
		if(particle.PdgCode() == 11 && prim.compare(particle.Process()) == 0) {
		start_pos_truth.SetXYZ(particle.Vx(),particle.Vy(),particle.Vz());
		}
		}

		}

		// Loop through tracks to find which one has the most energy from the primary electron

		for(int i=0; i<numtracks; ++i) {

		prim_energy = 0.0;

		// Loop through hits in track to add up primary electron's energy in track

		for(unsigned int h=0; h<fmth.at(i).size(); ++h) {

		// Loop through trackIDEs of hit to get energy of primary electron in hit

		for(auto const & trackide : bt->HitToTrackIDEs(fmth.at(i)[h])){

		// Use particle inventory service to get MCParticle from Geant track ID (different from longest_track)
		const simb::MCParticle* part = pi_serv->TrackIdToParticle_P(trackide.trackID);

		if(part->PdgCode() == 11 && prim.compare(part->Process()) == 0) {

		prim_energy += trackide.energy;

		}

		} // End of loop through TrackIDEs

		// Check if this hit is closest to true starting position
		// (Checking alternate definition of "true" primary track - closest to true
		// starting position)

		std::vector<double> hit = bt->HitToXYZ(fmth.at(i)[h]);
		double dist = std::sqrt(std::pow((hit[0] - start_pos_truth[0]),2)
		+ std::pow((hit[1] - start_pos_truth[1]),2)
		+ std::pow((hit[2] - start_pos_truth[2]),2));

		if(dist < min_dist) {
		min_dist = dist;
		min_dist_trk = i;
		}

	} // End of loop through hits

	if(prim_energy > max_en) {
		max_en_trk = i;
		max_en = prim_energy;
	}

	} // End of loop through tracks

	if(min_dist > maxmindist) maxmindist = min_dist;
	if(max_en_trk==longest_track) ++longest_is_prim_largest_en_frac;
	if(min_dist_trk==longest_track) ++longest_has_prim_first_pt;

	} // End of if(fmth.isValid() && numtracks!=0)
	*/
		//-----------------------------------------------------------------------------------------------------
		// Daughter flipping on longest track
		//-----------------------------------------------------------------------------------------------------
	


	TVector3 dir_reco_longest_dflip;
	TVector3 dir_reco_longest_flipped;
	TVector3 longest_start;
	TVector3 longest_end;
	double avg_daughter_cos_start = 0.0;
	double avg_daughter_cos_end = 0.0;
	int num_daughter_tracks = 0;
	TVector3 daughter_start;
	daughter_start.SetXYZ(0,0,0);

	// Get Handle of tracks if not done already
	//tracklistHandle
	//  = event.getValidHandle<std::vector<recob::Track>> (fTrackModuleLabel);

	if( numtracks != 0 ) {

		// Save reconstructed direction of electron
		dir_reco_longest_dflip.SetXYZ(tracklistHandle->at(longest_track).StartDirection().X(),
				tracklistHandle->at(longest_track).StartDirection().Y(),
				tracklistHandle->at(longest_track).StartDirection().Z());
		// save to tree:
		elec_x_ndf = dir_reco_longest_dflip.X();
		elec_y_ndf = dir_reco_longest_dflip.Y();
		elec_z_ndf = dir_reco_longest_dflip.Z();
		cosAngle_ndf = TMath::Cos(dir_reco_longest_dflip.Angle(elec_true_dir));
		// Save negative of end direction of electron -- this will be the direction
		// if we decide to flip it
		dir_reco_longest_flipped.SetXYZ(-tracklistHandle->at(longest_track).EndDirection().X(),
				-tracklistHandle->at(longest_track).EndDirection().Y(),
				-tracklistHandle->at(longest_track).EndDirection().Z());
		// Save start and end points of electron track for use in daughter flipping
		longest_start.SetXYZ(tracklistHandle->at(longest_track).Start().X(),
				tracklistHandle->at(longest_track).Start().Y(),
				tracklistHandle->at(longest_track).Start().Z());
		longest_end.SetXYZ(tracklistHandle->at(longest_track).End().X(),
				tracklistHandle->at(longest_track).End().Y(),
				tracklistHandle->at(longest_track).End().Z());

		// Loop through all tracks except for the longest one (daughter tracks). Find avg cos of
		// angles between electron direction and vector from the the vertex of 
		// the track to the daughter track, for each end.

		for ( int i=0; i<numtracks; ++i ) {

			if( i != longest_track ) {

				++num_daughter_tracks;

				daughter_start.SetXYZ(tracklistHandle->at(i).Start().X(),
						tracklistHandle->at(i).Start().Y(),
						tracklistHandle->at(i).Start().Z());

				avg_daughter_cos_start += TMath::Cos(dir_reco_longest_dflip.Angle(daughter_start-longest_start));
				avg_daughter_cos_end += TMath::Cos(dir_reco_longest_flipped.Angle(daughter_start-longest_end));

			}

		} // End loop through daughter tracks

		avg_daughter_cos_start /= num_daughter_tracks;
		avg_daughter_cos_end /= num_daughter_tracks;

		// If avg cos from what reco called the "end" side of the track is
		// is larger than what reco called the "start" side of the track,
		// flip the track.

		if ( avg_daughter_cos_start < avg_daughter_cos_end ) {

			dir_reco_longest_dflip = dir_reco_longest_flipped;

		}

		// Save post-flipping reco'ed electron direction in tree
		elec_x_df = dir_reco_longest_dflip.X();
		elec_y_df = dir_reco_longest_dflip.Y();
		elec_z_df = dir_reco_longest_dflip.Z();

		cosAngle_df = TMath::Cos(dir_reco_longest_dflip.Angle(elec_true_dir));
	}// end daughter flipping
	//std::cout<<"True Dir: " << elec_true_dir.X() << "\t" << elec_true_dir.Y() << "\t" <<  elec_true_dir.Z() << std::endl;
	//std::cout<<"ndf Dir: " << elec_x_ndf << "\t" << elec_y_ndf << "\t" <<  elec_z_ndf << std::endl;
	//std::cout<<"df Dir: " << elec_x_df << "\t" << elec_y_df << "\t" <<  elec_z_df << std::endl;
	//std::cout << "df cosAngle: " << cosAngle_df << " ndf cosAngle: "<< cosAngle_ndf<< std::endl;


	//-----------------------------------------------------------------------------------------------------
	//Total event charge -- for energy reco for clean events
	//-----------------------------------------------------------------------------------------------------

	auto hitListHandle = event.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);

	double tot_charge_U = 0.0;
	double tot_charge_V = 0.0;
	double tot_charge_Z = 0.0;
	double totalevtCharge = 0.0;

	// Loop over hits and add up charge for each wire plane

	for(size_t iHit = 0; iHit < hitListHandle->size(); ++iHit){

		art::Ptr<recob::Hit> hitPtr(hitListHandle, iHit);

		if(hitPtr->View() == geo::kU) {
			tot_charge_U += hitPtr->Integral();
		}
		if(hitPtr->View() == geo::kV) {
			tot_charge_V += hitPtr->Integral();
		}
		if(hitPtr->View() == geo::kZ) {
			tot_charge_Z += hitPtr->Integral();
		}

	}

	// Save charge from wire plane with the most charge as
	// the total event charge

	if(tot_charge_Z>=tot_charge_U) {
		totalevtCharge = tot_charge_Z;
	}
	else {
		totalevtCharge = tot_charge_U;
	}
	if(tot_charge_V>totalevtCharge) {
		totalevtCharge = tot_charge_V;
	}

	//-----------------------------------------------------------------------------------------------------
	// Primary track charge -- for energy reco for events w bkgs/noise
	//-----------------------------------------------------------------------------------------------------

	double trk_charge = 0.0;

	art::FindManyP<recob::Hit> fmth(tracklistHandle, event, fTrackModuleLabel);

	if ( fmth.isValid() && numtracks != 0) {

		//std::cout << "fmth is valid. Size: " << fmth.size() << "\n";

		// Get the hits of the longest track
		std::vector< art::Ptr<recob::Hit> > vhit = fmth.at(longest_track);

		//std::cout << "vhit size: " << vhit.size() << "\n";

		// Loop through hits and add up charge for each wire plane

		for (size_t i = 0; i < vhit.size(); ++i){

			if (vhit[i]->WireID().Plane == geo::kU) trk_charge_U += vhit[i]->Integral();
			if (vhit[i]->WireID().Plane == geo::kV) trk_charge_V += vhit[i]->Integral();
			if (vhit[i]->WireID().Plane == geo::kZ) trk_charge_Z += vhit[i]->Integral();

		}

	}

	//std::cout << "trk_charge_U = " << trk_charge_U << std::endl;
	//std::cout << "trk_charge_V = " << trk_charge_V << std::endl;
	//std::cout << "trk_charge_Z = " << trk_charge_Z << std::endl;

	// Save charge from wire plane with the most charge as
	// the track charge

	if(trk_charge_Z > trk_charge_U) trk_charge = trk_charge_Z;
	else trk_charge = trk_charge_U;
	if(trk_charge_V > trk_charge) trk_charge = trk_charge_V;

	//std::cout << "trk_charge = " << trk_charge << std::endl;

	//-----------------------------------------------------------------------------------------------------
	// Energy reco - for clean events and no cut on opflash time
	//-----------------------------------------------------------------------------------------------------
	/*    
		  double noDCb = 31.1126;
		  double noDCm = 90.3919;
		  elec_reco_en_noDC = (totalevtCharge-noDCb)/noDCm;
		  elec_reco_en_DC = elec_reco_en_noDC;

		  if(numtracks != 0) {

	//setting up detectorClocks
	auto const* detectorClocks = lar::providerFrom< detinfo::DetectorClocksService >();

	//handle for the opFlash
	auto opflashHandle = event.getValidHandle<std::vector<recob::OpFlash>>(fOpFlashLabel);
	double photonFlashTime = 0.0;
	std::vector<double> totalPEs; std::vector<float> flashtimes;

	//loop thru flashes and make vectors of PEs and flash times
	for (std::size_t iFlash = 0; iFlash < opflashHandle->size(); ++iFlash) {
	art::Ptr<recob::OpFlash> opFlashPtr(opflashHandle, iFlash);
	totalPEs.emplace_back(opFlashPtr->TotalPE());
	//std::cout << "opFlashPtr->TotalPE() = " << opFlashPtr->TotalPE() << std::endl;
	flashtimes.emplace_back(opFlashPtr->Time()); //units us
	//std::cout << "opFlashPtr->Time() = " << opFlashPtr->Time() << std::endl;
	}

	//find max opflash and save time of that flash
	if(opflashHandle->size() > 0){
	int index = std::distance(totalPEs.begin(), std::max_element(totalPEs.begin(), totalPEs.end()));
	photonFlashTime = flashtimes[index];
	max_flash_time = photonFlashTime;
	}

	//declare variables
	std::vector<double> hitPeakTimes;
	double driftRecobTime = 0.0;

	//get track information
	//tracklistHandle
	//= event.getValidHandle<std::vector<recob::Track>> (fTrackModuleLabel);

	//find hits associated w tracks
	//commented out bc created above
	//art::FindMany<recob::Hit> fmth(tracklistHandle, event, fTrackModuleLabel);
	//art::FindManyP<recob::Hit> fmth(tracklistHandle, event, fTrackModuleLabel);

	//hits in longest track
	//std::vector<const recob::Hit*> hits = fmth.at(longest_track);
	std::vector< art::Ptr<recob::Hit> > hits = fmth.at(longest_track);

	//now we can set the time if there is information available
	if(hits.size() > 0 && opflashHandle->size() > 0){ //check that there are hits associated with track

	driftRecobTime = detectorClocks->TPCTick2Time(hits[0]->PeakTime()) - photonFlashTime;

	}

	auto const* detectorProperties =
	lar::providerFrom< detinfo::DetectorPropertiesService >();
	float tau = detectorProperties->ElectronLifetime(); //microseconds

	//drift-corrected charge
	double totalevtChargeDC = totalevtCharge*exp(driftRecobTime/tau);

	double DCb = 57.8044;
	double DCm = 127.382;

	if(driftRecobTime != 0) elec_reco_en_DC = (totalevtChargeDC-DCb)/DCm;

	}
	*/
	//-----------------------------------------------------------------------------------------------------
	// Energy reco - uses track charge and has cut on opflash time
	//-----------------------------------------------------------------------------------------------------

	double noDCb = 31.1126;
	double noDCm = 90.3919;
	elec_reco_en_noDC = (trk_charge-noDCb)/noDCm;
	elec_reco_en_DC = elec_reco_en_noDC;

	if(numtracks != 0) {

		// Set up detectorClocks
		// Modified for lar v9
		//auto const detectorClockData= lar::providerFrom< detinfo::DetectorClocksService >();
		auto const detectorClockData = 
			art::ServiceHandle<detinfo::DetectorClocksService const>() -> DataFor(event);
		auto const * detectorClocks = &detectorClockData;

		// Handle for the opFlash
		auto opflashHandle = event.getValidHandle<std::vector<recob::OpFlash>>(fOpFlashLabel);
		double photonFlashTime = 0.0;
		std::vector<double> totalPEs; std::vector<float> flashtimes;

		//std::cout << "opflashHandle->size() = " << opflashHandle->size() << std::endl;
		if(opflashHandle->size() > 0) ++opflashcounter;

		// Loop thru flashes and make vectors of PEs and flash times
		for (std::size_t iFlash = 0; iFlash < opflashHandle->size(); ++iFlash) {
			art::Ptr<recob::OpFlash> opFlashPtr(opflashHandle, iFlash);
			totalPEs.emplace_back(opFlashPtr->TotalPE());
			//std::cout << "opFlashPtr->TotalPE() = " << opFlashPtr->TotalPE() << std::endl;
			flashtimes.emplace_back(opFlashPtr->Time()); //units us
			//std::cout << "opFlashPtr->Time() = " << opFlashPtr->Time() << std::endl;
		}

		// Declare variables
		std::vector<double> hitPeakTimes;
		double driftRecobTime = 0.0;

		// Get track information
		//tracklistHandle
		//= event.getValidHandle<std::vector<recob::Track>> (fTrackModuleLabel);

		// Find hits associated w tracks
		//commented out bc created above
		//art::FindMany<recob::Hit> fmth(tracklistHandle, event, fTrackModuleLabel);
		//art::FindManyP<recob::Hit> fmth(tracklistHandle, event, fTrackModuleLabel);

		// Hits in longest track
		//std::vector<const recob::Hit*> hits = fmth.at(longest_track);
		std::vector< art::Ptr<recob::Hit> > hits = fmth.at(longest_track);

		// Now we can set the time if there is information available
		if(hits.size() > 0 && opflashHandle->size() > 0){ //check that there are hits associated with track

			double hitTime = detectorClocks->TPCTick2Time(hits[0]->PeakTime());

			// Find max opflash within physical constraints and save time of that flash
			while(1) {
				int index = std::distance(totalPEs.begin(), std::max_element(totalPEs.begin(), totalPEs.end()));
				//std::cout << "hitTime = " << hitTime << std::endl;
				//std::cout << "max flashtime = " << flashtimes[index] << std::endl;

				// If the max PE opflash makes sense physically, save as photonFlashTime
				if( flashtimes[index] <= hitTime && hitTime - flashtimes[index] < 2400.0) {
					photonFlashTime = flashtimes[index];
					max_flash_time = photonFlashTime;
					//std::cout << "photonFlashTime = " << photonFlashTime << std::endl;
					break;
				}
				// Else break loop if that was the last opflash
				else if( totalPEs.size()==1 ) break;
				// Else remove max PE opflash in order to find second highest PE opflash
				else{
					totalPEs.erase(totalPEs.begin() + index);
					flashtimes.erase(flashtimes.begin() + index);
				}
			}

			// Set driftRecobTime here
			driftRecobTime = hitTime - photonFlashTime;

		}

		//std::cout << "driftRecobTime = " << driftRecobTime << std::endl;

		// modified for lar v9
		//auto const* detectorProperties =
		//		lar::providerFrom< detinfo::DetectorPropertiesService >();

		auto const detectorPropertiesData = art::ServiceHandle<detinfo::DetectorPropertiesService const>() -> DataFor(event);
		//auto const * detectorProperties = &detectorPropertiesData;

		double tau = detectorPropertiesData.ElectronLifetime(); //microseconds

		// Drift-corrected charge
		//double totalevtChargeDC = totalevtCharge*exp(driftRecobTime/tau);
        trk_charge_noDC = trk_charge;
		trk_charge_DC = trk_charge * exp(driftRecobTime/tau);

		double DCb = 57.8044;
		double DCm = 127.382;

		// Drift-corrected energy
		//if(driftRecobTime != 0) elec_reco_en = (totalevtChargeDC-DCb)/DCm;
		if(driftRecobTime != 0) elec_reco_en_DC = (trk_charge_DC-DCb)/DCm;

		//std::cout << "elec_reco_en_noDC = " << elec_reco_en_noDC << std::endl;
		//std::cout << "elec_reco_en_DC = " << elec_reco_en_DC << std::endl;

	}


	//-----------------------------------------------------------------------------------------------------
	// Fill tree - last thing in analyze function, at the end of each event analysis
	//-----------------------------------------------------------------------------------------------------

	tr->Fill();

	}

	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------------------

	void ElectronInfo::endJob(){

		std::cout << "has_tracks = " << has_tracks << std::endl;
		std::cout << "opflashcounter = " << opflashcounter << std::endl;
		std::cout << "longest_is_prim_largest_en_frac = " << longest_is_prim_largest_en_frac << std::endl;
		std::cout << "longest_has_prim_first_pt = " << longest_has_prim_first_pt << std::endl;
		std::cout << "maxmindist = " << maxmindist << std::endl;

	}

}

