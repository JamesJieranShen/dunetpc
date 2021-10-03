/**
 * Retrieve and calculate relevant infromation from the art reco file. Save all information in a tree.
 * @author James Shen <jieran.shen@duke.edu>
 * Based on the original PointResTree module by AJ Roeth <ajroeth110@gmail.com>
 * 
 * Specifically, this module does the following:
 *      - Retrieve relevant MC information for the neutrino and the electron.
 *      - Resolve track directional ambiguity via daughter flipping.
 *      - Reconstruct electron energy using charge information from the planes, drift corrected by photon flashes.
 * 
 * Structure of the tree:
 *      - TVector3 truth_e_dir: truth direction of electron;
 *      - TVector3 truth_nu_dir: truth direction of neutrino;
 *      - double truth_nu_en: truth energy of neutrino;
 *      - double truth_e_en: truth energy of electron;
 *      - int nu_pdg: PDG of neutrino;
 *      - int Nhits: number of hits in the event;
 *      - int NTrks: number of tracks for reco;
 *      - double charge_U: charge from U plane;
 *      - double charge_V: charge from V plane;
 *      - double charge_Z: charge from Z plane;
 *      - double charge_corrected: corrected charge used for energy reco;
 *      - double drift_time: drift time of electron;
 *      - TVector3 reco_e_dir: reco direction of electron;
 *      - double reco_e_en: reco energy of electron;
 *      - double reco_distance distance between reco start position and truth start position
 *      - int primary_trk_id: index of primary track;
 *      - The following branches are parallel arrays of track info.
 *          - vector<double> trk_length: length of tracks;
 *          - vector<double> trk_start_x, trk_start_y, trk_start_z: unit vectors of start positions;
 *          - vector<double> trk_end_x, trk_end_y, trk_end_z: unit vectors of end positions;
 *          - vector<double> trk_start_dir_x, trk_start_dir_y, trk_start_dir_z: unit vectors of start directions;
 *          - vector<double> trk_end_dir_x, trk_end_dir_y, trk_end_dir_z: unit vectors of end directions;
 * 
 * */

#ifndef PointResTree_H
#define PointResTree_H 1

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
#include "dune/AnaUtils/DUNEAnaHitUtils.h"

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
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"
#include "TPolyLine3D.h"
#include "Math/Vector3D.h"
#include "TFile.h"
#include "TPrincipal.h"

// C++ includes

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

namespace dune {
    struct genFinder{
    private:
        typedef std::pair<int, std::string> track_id_to_string;
        std::vector<track_id_to_string> track_id_map;
        std::set<std::string> generator_names;
        bool isSorted = false;

    public:
        void sort_now(){
            std::sort(this->track_id_map.begin(),
                      this->track_id_map.end(),
                      [](const auto &a, const auto &b)
                      { return (a.first < b.first); });
            isSorted = true;
        }

        void add(const int &track_id, const std::string &gname){
            this->track_id_map.push_back(std::make_pair(track_id, gname));
            generator_names.emplace(gname);
            isSorted = false;
        }

        bool has_gen(std::string gname){
            return static_cast<bool>(generator_names.count(gname));
        };

        std::string get_gen(int tid){
            if (!isSorted)
                this->sort_now();
            return std::lower_bound(track_id_map.begin(),
                                    track_id_map.end(),
                                    tid,
                                    [](const auto &a, const auto &b)
                                    { return (a.first < b); })
                ->second;
        };
    };

    class PointResTree : public art::EDAnalyzer {

		public:

			PointResTree(fhicl::ParameterSet const&);
            void reconfigure(fhicl::ParameterSet const&);
            void beginJob() override;
			void analyze(art::Event const&) override;
			void endJob() override;
        
        private:
            using XYZVector=ROOT::Math::XYZVector;

            //NOT IN TREE:
            genFinder* gf;
            std::string primary_trk_label;
            Int_t primary_trk_idx; // index of the primary track IN THE TRACK LABEL
            XYZVector com_position;

            // art module labels
            std::string fSimulationLabel;
			std::vector<std::string> fTrackModuleLabel;
			std::string fHitsModuleLabel;
			std::string fOpFlashLabel;
            std::string fHitFdModuleLabel;
			std::string fHitToSpacePointLabel;

            // Tree & branches
            TH1D * hit_hist;
			TH1D * good_hit_hist;
            TH1D * hit_distance;
            TH1D * trk_length_hist;
            TTree * tr;
            XYZVector truth_nu_dir, truth_e_dir;
            Double_t truth_nu_en, truth_e_en;
            Double_t truth_en_deposition, truth_charge_deposition;
            Int_t nu_pdg;
            XYZVector truth_e_position;

            Int_t NParticles, NHits, NTrks;
            Double_t charge_U, charge_V, charge_Z;
            std::vector<Double_t> charge_Z_dist;
            Double_t charge_corrected;
            Double_t cut_charge_loss;
            Double_t drift_time;
            
            XYZVector reco_e_dir;
            Double_t reco_e_en;
            XYZVector reco_e_position;


            Double_t reco_distance;
            Int_t primary_trk_id;

            std::vector<Double_t> trk_length;
            std::vector<Double_t> trk_charge;
            std::vector<Double_t> trk_start_x, trk_start_y, trk_start_z;
            std::vector<Double_t> trk_end_x, trk_end_y, trk_end_z;
            std::vector<Double_t> trk_start_dir_x, trk_start_dir_y, trk_start_dir_z;
            std::vector<Double_t> trk_end_dir_x, trk_end_dir_y, trk_end_dir_z;


            // helper functions
            void reset_variables();
            void writeMCTruths_marley(art::Event const&);
            void writeMCTruths_largeant(art::Event const&);
            void calculate_electron_energy(art::Event const&);
            void calculate_electron_direction();
            void get_primary_track();
            Bool_t distance_cut(art::Event const&, recob::Hit const&, std::vector<art::Ptr<recob::SpacePoint>>);

            /* Function to write hit as a distance of truth distance. */
            void write_multi_distance(art::Event const&, recob::Hit const&, std::vector<art::Ptr<recob::SpacePoint>>);
            
	}; // class PointResTree

}// namespace dune

#endif 

namespace dune {
	DEFINE_ART_MODULE(PointResTree)
}

dune::PointResTree::PointResTree(fhicl::ParameterSet const& parameterSet)
        :EDAnalyzer(parameterSet){
    this->reconfigure(parameterSet);
}

void dune::PointResTree::reconfigure(fhicl::ParameterSet const& parameterSet){
    fSimulationLabel = parameterSet.get< std::string >("SimulationLabel");
	fTrackModuleLabel = parameterSet.get< std::vector<std::string> >("TrackModuleLabel");
	fHitsModuleLabel = parameterSet.get< std::string >("HitsModuleLabel");
	fOpFlashLabel = parameterSet.get< std::string >("OpFlashLabel");
	fHitToSpacePointLabel = parameterSet.get<std::string>("HitToSpacePointLabel");
}

void dune::PointResTree::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    hit_hist = tfs->make<TH1D>("hit_hist", "Hit Amplitudes", 100, 0, 50);
	good_hit_hist = tfs->make<TH1D>("good_hit_hist", "Hit Amplitude with associated IDE", 100, 0, 50);
    hit_distance = tfs->make<TH1D>("hit_distance", "Distance from hit to start position of electron", 100, 0, 400);
    trk_length_hist = tfs->make<TH1D>("trk_length_hist", "Length of Tracks", 100, 0, 50);

    tr = tfs->make<TTree>("tr", "tr");

    tr->Branch("truth_e_dir", &truth_e_dir);
    tr->Branch("truth_nu_dir", &truth_nu_dir);
    tr->Branch("truth_nu_en", &truth_nu_en);
    tr->Branch("truth_e_en", &truth_e_en);

    tr->Branch("truth_en_deposition", &truth_en_deposition);
    tr->Branch("truth_charge_deposition", &truth_charge_deposition);
    tr->Branch("truth_e_position", &truth_e_position);

    tr->Branch("nu_pdg", &nu_pdg);

	tr->Branch("NParticles", &NParticles);
    tr->Branch("NHits", &NHits);
    tr->Branch("NTrks", &NTrks);
    tr->Branch("charge_U", &charge_U);
    tr->Branch("charge_V", &charge_V);
    tr->Branch("charge_Z", &charge_Z);
    tr->Branch("charge_Z_dist", &charge_Z_dist);
    tr->Branch("charge_corrected", &charge_corrected);
    tr->Branch("cut_charge_loss", &cut_charge_loss);
    tr->Branch("drift_time", &drift_time);

    tr->Branch("reco_e_dir", &reco_e_dir);
    tr->Branch("reco_e_en", &reco_e_en);
    tr->Branch("reco_e_position", &reco_e_position);

    tr->Branch("reco_distance", &reco_distance);
    tr->Branch("primary_trk_id", &primary_trk_id);

    tr->Branch("trk_length", &trk_length);
    tr->Branch("trk_charge", &trk_charge);
    tr->Branch("trk_start_x", &trk_start_x);
    tr->Branch("trk_start_y", &trk_start_y);
    tr->Branch("trk_start_z", &trk_start_z);

    tr->Branch("trk_end_x", &trk_end_x);
    tr->Branch("trk_end_y", &trk_end_y);
    tr->Branch("trk_end_z", &trk_end_z);

    tr->Branch("trk_start_dir_x", &trk_start_dir_x);
    tr->Branch("trk_start_dir_y", &trk_start_dir_y);
    tr->Branch("trk_start_dir_z", &trk_start_dir_z);

    tr->Branch("trk_end_dir_x", &trk_end_dir_x);
    tr->Branch("trk_end_dir_y", &trk_end_dir_y);
    tr->Branch("trk_end_dir_z", &trk_end_dir_z);

}// beginJob()


void dune::PointResTree::endJob(){
    std::cout<<"Ending Job" << std::endl;
} //endJob()


void dune::PointResTree::analyze(art::Event const& event){
    reset_variables();
    gf = new genFinder();
    // ... Create a map of track IDs to generator labels
    //Get a list of generator names.
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mcHandles;
    event.getManyByType(mcHandles);
    std::vector<std::pair<int, std::string>> track_id_to_label;

    for (auto const &mcHandle : mcHandles)
    {
        const std::string &sModuleLabel = mcHandle.provenance()->moduleLabel();
        art::FindManyP<simb::MCParticle> findMCParts(mcHandle, event, fSimulationLabel);
        std::vector<art::Ptr<simb::MCParticle>> mcParts = findMCParts.at(0);
        for (const art::Ptr<simb::MCParticle> ptr : mcParts)
        {
            int track_id = ptr->TrackId();
            gf->add(track_id, sModuleLabel);
        }
    }

    // Get MC Truths
    if(fSimulationLabel == "marley")
        writeMCTruths_marley(event);
    if(fSimulationLabel == "largeant")
        writeMCTruths_largeant(event);

    // get track info
    // determine primary track: IF pandora, choose longest track
    // IF PMA, use nearest neighbor
    primary_trk_label="NA";
    Double_t max_trk_length = 0;
    for(std::string track_label:fTrackModuleLabel){
        auto trk_list = event.getValidHandle<std::vector<recob::Track>>(track_label);
        NTrks += trk_list->size();
        art::FindManyP<recob::Hit> hitsFromTracks(trk_list, event, track_label);

        for(int i = 0; i < (int)trk_list->size(); i++){// using index loop to get track idx
            recob::Track const& trk = trk_list->at(i);
            trk_length.push_back(trk.Length());
            trk_length_hist->Fill(trk.Length());
            if(track_label=="pandoraTrack" && trk.Length() > max_trk_length){
                max_trk_length = trk.Length();
                primary_trk_id = trk_length.size() - 1; // most recent entry
                primary_trk_label = track_label;
                primary_trk_idx = i;
            }
            

            // double distance = sqrt(pow(trk.Start().X() - truth_e_position.X(), 2) +
            //                         pow(trk.Start().Y() - truth_e_position.Y(), 2) + 
            //                         pow(trk.Start().Z() - truth_e_position.Z(), 2));
            // std::cout<<distance<<std::endl;
            trk_start_x.push_back(trk.Start().X());
            trk_start_y.push_back(trk.Start().Y());
            trk_start_z.push_back(trk.Start().Z());
            trk_end_x.push_back(trk.End().X());
            trk_end_y.push_back(trk.End().Y());
            trk_end_z.push_back(trk.End().Z());

            trk_start_dir_x.push_back(trk.StartDirection().X());
            trk_start_dir_y.push_back(trk.StartDirection().Y());
            trk_start_dir_z.push_back(trk.StartDirection().Z());
            trk_end_dir_x.push_back(trk.EndDirection().X());
            trk_end_dir_y.push_back(trk.EndDirection().Y());
            trk_end_dir_z.push_back(trk.EndDirection().Z());

            // loop through all charges associated with this track to accumulate charge
            auto hits = hitsFromTracks.at(i);
            float q = 0;
            for(auto hit: hits){
                if(hit->View() == geo::kZ){
                    q += hit->Integral();
                }
            }
            trk_charge.push_back(q);

        } //end loop through trks
    }// loop through track labels

    if(NTrks!=0 && primary_trk_id == -1){ // tracks exist yet primary track not initialized, meaning all tracks are pmtrack
        get_primary_track();

    }
    //std::cout << std::distance(trk_charge.begin(),std::max_element(trk_charge.begin(), trk_charge.end()))<<std::endl;
    //std::cout << primary_trk_id << std::endl;
    calculate_electron_direction(); 
    auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    // auto tracklistHandle = event.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);
    // art::FindManyP<recob::Hit> hitsFromTracks(tracklistHandle, event, fTrackModuleLabel);
    // art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    
    // // look at gaushits FIXME:
    // auto gaushit_handle = event.getValidHandle<std::vector<recob::Hit>>("gaushit");
    // for(int i = 0; i < (int)gaushit_handle->size(); i++){
    //     auto const hit = gaushit_handle->at(i);
    //     if(hit.View() == geo::kU || hit.View() == geo::kV)
    //         std::cout<<(hit.WireID().TPC%2) << std::endl;
    // }

    // get hit info
    auto hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
    std::vector<art::Ptr<recob::Hit>> hit_list;
    art::fill_ptr_vector(hit_list, hit_handle);
    auto const detectorPropertiesData = 
            art::ServiceHandle<detinfo::DetectorPropertiesService const>() -> DataFor(event);
	art::ServiceHandle<cheat::BackTrackerService> bt;
    float NHits_Z = 0;
    float NspHits_Z = 0;
    NHits = hit_list.size();
    double uncut_charge = 0;

    for(int i = 0; i<NHits; i++){
        auto const hit = hit_list.at(i);
        //if(hit->PeakAmplitude()<3.0) continue;
        
        
        auto const trk_IDEs = bt->HitToTrackIDEs(clockData, hit);
        //if(hit->View() == geo::kZ){
            hit_hist->Fill(hit->PeakAmplitude());
            if(trk_IDEs.size() != 0) {
				good_hit_hist->Fill(hit->PeakAmplitude());
			}

        //}
        if(hit->PeakAmplitude()>4.8){
            // for induction planes, just sum all hits (not used)
            if(hit->View() == geo::kU) charge_U += hit->Integral();
            if(hit->View() == geo::kV) charge_V += hit->Integral();
            
            if(hit->View() == geo::kZ) uncut_charge += hit->Integral();
            // for collection plane, apply distance cut
			if(hit->View() == geo::kZ){
                auto spacepoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(hit, event, 
				        fHitsModuleLabel, fHitToSpacePointLabel);
                NHits_Z++;
                if(spacepoints.size()) NspHits_Z++;
                write_multi_distance(event, *hit, spacepoints);
                if(distance_cut(event, *hit, spacepoints)){
                    // for(auto trkide : trk_IDEs){
                    //     auto particle = pi_serv->TrackIdToParticle(trkide.trackID);
                    //     std::cout<<"pdg:\t" << particle.PdgCode() << "\tProcess:\t" << particle.Process() << std::endl;
                    // }
                    continue;
                }
                else charge_Z+=hit->Integral();
            }
        }
        
    } //end loop through hits
    cut_charge_loss = charge_Z/uncut_charge;
    // std::cout<<"Spacepoint reco hits / total hits: " << NspHits_Z/NHits_Z << std::endl;
    // std::cout<<NHits_Z/NHits<< std::endl;
    //hit_distance->Fill(NspHits_Z/NHits_Z);
	reco_distance = (truth_e_position - reco_e_position).R();

    
    
    calculate_electron_energy(event);
    
    delete gf;
    tr->Fill();
} //analyze

//=================================================================================================
// Private Helper methods

/**
 * reset all tree variables to impossible values
 * */
void dune::PointResTree::reset_variables(){
    truth_nu_dir.SetXYZ(-2, -2, -2);
    truth_e_dir.SetXYZ(-2, -2, -2);
    reco_e_dir.SetXYZ(-2, -2, -2);
    
    reco_e_en = truth_nu_en = truth_e_en = -10;
    nu_pdg = 0;

    NParticles = NHits = NTrks = 0;
    charge_U = charge_V = charge_Z = charge_corrected = 0;
    charge_Z_dist.resize(20, 0);
    std::fill(charge_Z_dist.begin(), charge_Z_dist.end(), 0);
    cut_charge_loss = 0;
    drift_time = -10;
    reco_distance = -10;
    primary_trk_id = -1;
    truth_en_deposition = 0;
    truth_charge_deposition = 0;

    trk_length.clear();
    trk_charge.clear();
    trk_start_x.clear();
    trk_start_y.clear();
    trk_start_z.clear();
    trk_end_x.clear();
    trk_end_y.clear();
    trk_end_z.clear();
    trk_start_dir_x.clear();
    trk_start_dir_y.clear();
    trk_start_dir_z.clear();
    trk_end_dir_x.clear();
    trk_end_dir_y.clear();
    trk_end_dir_z.clear();


} //reset_variables


/**
 * Set Neutrino and Electron truths to tree variables, assuming marley generator
 * Wrties truth_nu_dir, truth_e_dir, truth_nu_en, truthe_e_en, nu_pdg
 * @param event art event object
 * */
void dune::PointResTree::writeMCTruths_marley(art::Event const& event){
    Bool_t found_nu, found_e = kFALSE;
    auto particleHandle_nu = event.getValidHandle<std::vector<simb::MCTruth>>(fSimulationLabel);
    for (auto const& particle : (*particleHandle_nu)){
        if(particle.NeutrinoSet()){ // entry is a neutrino
            found_nu = kTRUE;
            auto const neutrino = particle.GetNeutrino().Nu();
            truth_nu_en = neutrino.E()*1000;
            nu_pdg = neutrino.PdgCode();

            truth_nu_dir.SetXYZ(neutrino.Px(),
                                neutrino.Py(),
                                neutrino.Pz());
            truth_nu_dir/=neutrino.P();
        }
    }

    auto particleHandle_elec = event.getValidHandle<std::vector<simb::MCParticle>> ("largeant");
    for(auto const& particle : (*particleHandle_elec)){
        if(particle.Process()=="primary" && particle.PdgCode()==11){
            found_e = kTRUE;
            truth_e_en = 1000*(particle.E() - particle.Mass());
            
            truth_e_dir.SetXYZ( particle.Px(),
                                particle.Py(),
                                particle.Pz());
            truth_e_dir/=particle.P();
        }
    }
    if (!found_nu)
        std::cerr<<"WARNING: Neutrino MC info not found" << std::endl;
    if (!found_e)
        std::cerr<<"WARNING: Electron MC info not found" << std::endl;
    
}


/**
 * Set Neutrino and Electron truths to tree variables, assuming geant4/NuE generator
 * Wrties truth_nu_dir, truth_e_dir, truth_nu_en, truthe_e_en, nu_pdg
 * @param event art event object
 * */
void dune::PointResTree::writeMCTruths_largeant(art::Event const& event){
    Bool_t found_nu, found_e = kFALSE;
    auto particleHandle
				= event.getValidHandle<std::vector<simb::MCParticle>> (fSimulationLabel);
	NParticles = 0;


    for (auto const& particle : (*particleHandle)) {
        
		if(gf->get_gen(particle.TrackId())!="generator") continue; // do not track radiologicals
		// neutrino
		if ( particle.Process() == "primary" && 
                (particle.PdgCode() == 12 
                    || particle.PdgCode() == 14
                    || particle.PdgCode() == 16
                    || particle.PdgCode() == -12 
                    || particle.PdgCode() == -14
                    || particle.PdgCode() == -16) ) {
            found_nu = kTRUE;
            truth_nu_en = 1000*(particle.E() - particle.Mass());

            nu_pdg = particle.PdgCode();

            truth_nu_dir.SetXYZ(particle.Px(),
                                particle.Py(),
                                particle.Pz());
            truth_nu_dir/=particle.P();
        }
		// electron
        if ( particle.Process() == "primary" && particle.PdgCode() == 11 ) {
            found_e = kTRUE;
            truth_e_en = 1000*(particle.E() - particle.Mass());
            truth_e_dir.SetXYZ( particle.Px(),
                                particle.Py(),
                                particle.Pz());
            truth_e_dir/=particle.P();
            truth_e_position = particle.Trajectory().Position(0).Vect();
			NParticles++;

        }
    }//endfor

    auto channel_list = event.getValidHandle<std::vector<sim::SimChannel>>("elecDrift");
    for(auto const& channel:(*channel_list)){
        auto const TDCIDEs = channel.TDCIDEMap();
        for(auto const& TDCIDE: TDCIDEs){
            auto simIDEs = TDCIDE.second;
            for(auto& IDE:simIDEs){
                truth_en_deposition+=IDE.energy;
                truth_charge_deposition+=IDE.numElectrons;
            }
        }
    }
	if(NParticles!=1)
		std::cerr<<"WARNING: More than one primary electron is found" << std::endl;
    if (!found_nu)
        std::cerr<<"WARNING: Neutrino MC info not found" << std::endl;
    if (!found_e)
        std::cerr<<"WARNING: Electron MC info not found" << std::endl;
}


/**
 * Given accumulated charges on the U, V, and Z planes, reconstruct the energy of the electron.
 * Use Optical Detector information for drift correction.
 * sets reco_e_en, charge_corrected, drift_time.
 * */
void dune::PointResTree::calculate_electron_energy(art::Event const& event){
    // use collection plane charge for reconstruction
    Double_t charge_raw = charge_Z;
    // drift correction
    if(NTrks == 0 || primary_trk_id == -1){ // no tracks found, do not do drift correction
        charge_corrected = charge_raw;
        return;
    }

    auto const detectorClockData = 
            art::ServiceHandle<detinfo::DetectorClocksService const>() -> DataFor(event);
    auto opflashHandle = event.getValidHandle<std::vector<recob::OpFlash>>(fOpFlashLabel);
    std::vector<double> totalPEs; std::vector<float> flashtimes;
    //std::cout<<opflashHandle->size()<<std::endl;
    // if(opflashHandle->size() == 0){ // do not do drift correction
    //     charge_corrected = charge_raw;
    //     return; 
    // }

    for(unsigned int i = 0; i < opflashHandle->size(); i++){
        auto flash = opflashHandle->at(i);
        totalPEs.push_back(flash.TotalPE());
        flashtimes.push_back(flash.Time()); //units us
    }

    std::vector<double> hitPeakTimes;
    drift_time = 0.0;
    // Create FindManyP object to find hits associated with Tracks
    auto tracklistHandle = event.getValidHandle<std::vector<recob::Track>>(primary_trk_label);
    art::FindManyP<recob::Hit> hitsFromTracks(tracklistHandle, event, primary_trk_label);
    // hits in primary track:
    auto hits_primary = hitsFromTracks.at(primary_trk_idx); // pointers to hits
    // for(auto hit : hits_primary){
    //     std::cout<<detectorClockData.TPCTick2Time(hit->PeakTime()) << std::endl;
    // }
    Double_t hitTime = detectorClockData.TPCTick2Time(hits_primary[0]->PeakTime());
    Double_t current_drift_time = 0;
    while(totalPEs.size() != 0){
        int max_flash_idx = std::distance(totalPEs.begin(), 
                std::max_element(totalPEs.begin(), totalPEs.end()));
        current_drift_time = hitTime - flashtimes[max_flash_idx];
		//std::cout<<current_drift_time<<std::endl;
        if(current_drift_time>0 && current_drift_time < 2400.0){ // drift time found
            drift_time = current_drift_time;
            break;
        }
        // if the flashtime found is not physical
        totalPEs.erase(totalPEs.begin() + max_flash_idx);
        flashtimes.erase(flashtimes.begin() + max_flash_idx);
    }
    // Double_t APA_distance = abs(truth_e_position.X());
    auto const detectorPropertiesData = 
            art::ServiceHandle<detinfo::DetectorPropertiesService const>() ->DataFor(event);
    Double_t electron_lifetime = detectorPropertiesData.ElectronLifetime(); //us, probably 10400
    // Double_t drift_velocity = detectorPropertiesData.DriftVelocity();
    // drift_time = APA_distance/drift_velocity; //using truth!!
    charge_corrected = charge_raw * ROOT::Math::exp(drift_time/electron_lifetime);

    //TODO: Calculate reco energy
}


/**
 * Reads the track direction information, calculate the direction of the primary track.
 * Sets reco_e_dir and reco_e_position. Assumes NTrk, track directions, and primary_trk_id is set.
 * */
void dune::PointResTree::calculate_electron_direction(){
    // std::cerr << "Daughter Flipping" << std::endl;
    // std::cerr << "NTrk: " << NTrks << std::endl;
    // std::cerr << "Vector size: " << trk_start_x.size() << std::endl;
    // std::cerr << "Primary Index:" << primary_trk_id << std::endl;
    if (primary_trk_id < 0) return; // if no trk is in reco, do nothing
    // double dist_thresh = 400;
    using XYZVector=ROOT::Math::XYZVector;
        //reset primary position after index update
    XYZVector primary_start(trk_start_x.at(primary_trk_id),
                            trk_start_y.at(primary_trk_id),
                            trk_start_z.at(primary_trk_id));
    
    XYZVector primary_end(      trk_end_x.at(primary_trk_id),
                                trk_end_y.at(primary_trk_id),
                                trk_end_z.at(primary_trk_id));
    XYZVector primary_forward(  trk_start_dir_x.at(primary_trk_id),
                                trk_start_dir_y.at(primary_trk_id),
                                trk_start_dir_z.at(primary_trk_id));
    XYZVector primary_backward( -trk_end_dir_x.at(primary_trk_id),
                                -trk_end_dir_y.at(primary_trk_id),
                                -trk_end_dir_z.at(primary_trk_id));

    
	
	Double_t sum_cos_forward = 0.0;
    Double_t sum_cos_backward = 0.0;
    XYZVector current_trk;
    XYZVector difference_start, difference_end; 
    // loop through all tracks to find avg cos between primary track and all daughter tracks.
	for (int i = 0; i < NTrks; i++)
    {
		if(i==primary_trk_id) continue; // don't include primary track
        current_trk.SetXYZ(trk_start_x.at(i), trk_start_y.at(i), trk_start_z.at(i));
        //if(trk_charge[i]<100) continue;
        // if((current_trk-com_position).R() > dist_thresh) {
        //     //std::cout<<(current_trk-com_position).R()<<std::endl;
        //     continue; // if track too far away, don't include
        // }
        difference_start = current_trk-primary_start;
        difference_end = current_trk - primary_end;
        sum_cos_forward += ROOT::Math::VectorUtil::CosTheta(difference_start, primary_forward);
        sum_cos_backward += ROOT::Math::VectorUtil::CosTheta(difference_end, primary_backward);

        // if(difference_start.mag2()==0){
        //     sum_cos_forward += 1;
        // }else{
        //     difference_start /= difference_start.R(); //noramlize to mag of 1
        //     sum_cos_forward+= primary_forward.Dot(difference_start);
        // }

		// difference_end = current_trk - primary_end;
        // if(difference_end.mag2()==0){
        //     sum_cos_backward += 1;
        // }else{
        //     difference_end /= difference_end.R();
        //     sum_cos_backward+= primary_backward.Dot(difference_end);
        // }
    }
    //std::cout<<sum_cos_forward << std::endl;
    //std::cout<<sum_cos_backward<<std::endl;
    if(sum_cos_backward > sum_cos_forward){
        reco_e_dir = primary_backward;
        reco_e_position = primary_end;
		//std::cout<<"Flipped"<<std::endl;
    }else{
        reco_e_dir = primary_forward;
        reco_e_position = primary_start;
    }

}


/**
 * Assuming all tracks have no length-wise significance, determine the primary track.
 * */
void dune::PointResTree::get_primary_track(){
    primary_trk_label="pmtracktc";

    Double_t com_thresh=45000; //cm
    Double_t com_x = 0;
    Double_t com_y = 0;
    Double_t com_z = 0;
    Double_t length_sum = 0;
    for (int i = 0; i < NTrks; i++){
        com_x += trk_start_x[i]*trk_length[i];
        com_y += trk_start_y[i]*trk_length[i];
        com_z += trk_start_z[i]*trk_length[i];
        length_sum += trk_length[i];
    }
    com_position.SetXYZ(com_x/length_sum, 
                            com_y/length_sum, 
                            com_z/length_sum);

    double max_trk_length = 0;
    for(int i = 0; i < NTrks; i++){
        double distance = sqrt(pow(trk_start_x[i] - com_position.X(), 2) +
                                pow(trk_start_y[i] - com_position.Y(), 2) + 
                                pow(trk_start_z[i] - com_position.Z(), 2));
        if(trk_length[i] > max_trk_length &&
            distance < com_thresh){
            max_trk_length = trk_length[i];
            primary_trk_id = i;
        }
    } 
    primary_trk_idx = primary_trk_id; // since there will be only one product
    
}


/**
 * Determine if a hit should be included in the energy sum.
 * Apply distance cut to hits with associated spacepoints, 2D cut (WireID and time) for hits without spacepoints.
 * @return True if should be cut, false if should be included.
 * */
Bool_t dune::PointResTree::distance_cut(art::Event const& event, recob::Hit const& hit, std::vector<art::Ptr<recob::SpacePoint>> spacepoints){
    if(primary_trk_id<0) return kFALSE; // don't cut anything if there are no tracks
    double dist_th = 14.0 * 5; // 14 is radiation length
    Double_t distance = 0;
    //art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    //auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    //art::ServiceHandle<cheat::BackTrackerService> bt;

    if(spacepoints.size()){
        // average of all spacepoint recos
        for(auto spacepoint:spacepoints){
            auto pos = XYZVector(spacepoint->XYZ()[0],
                                spacepoint->XYZ()[1],
                                spacepoint->XYZ()[2]);
            distance += (pos - reco_e_position).R();
        }
        distance /= spacepoints.size();
    }else{ // no spacepoint found, do 2D reconstruction (X, Z)
        geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
        auto const detProperties = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(event);
        auto wireID = hit.WireID();
        auto plane = geom.Plane(wireID);
        auto pos = plane.Wire(wireID).GetCenter<geo::Point_t>();
        plane.DriftPoint(pos, -abs(detProperties.ConvertTicksToX(hit.PeakTime(), wireID))); // drift away from plane
        distance = std::pow(pos.X() - reco_e_position.X(), 2)+std::pow(pos.Z() - reco_e_position.Z(), 2);
        distance = sqrt(distance);
    }
    hit_distance->Fill(distance);

    return(distance > dist_th);

}


void dune::PointResTree::write_multi_distance(art::Event const& event, recob::Hit const& hit, std::vector<art::Ptr<recob::SpacePoint>> spacepoints){
    if(primary_trk_id<0) return; // don't cut anything if there are no tracks
    Double_t distance = 0;
    //art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    //auto const& clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    //art::ServiceHandle<cheat::BackTrackerService> bt;

    if(spacepoints.size()){
        // average of all spacepoint recos
        for(auto spacepoint:spacepoints){
            auto pos = XYZVector(spacepoint->XYZ()[0],
                                spacepoint->XYZ()[1],
                                spacepoint->XYZ()[2]);
            distance += (pos - truth_e_position).R();
        }
        distance /= spacepoints.size();
    }else{ // no spacepoint found, do 2D reconstruction (X, Z)
        geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
        auto const detProperties = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(event);
        auto wireID = hit.WireID();
        auto plane = geom.Plane(wireID);
        auto pos = plane.Wire(wireID).GetCenter<geo::Point_t>();
        plane.DriftPoint(pos, -abs(detProperties.ConvertTicksToX(hit.PeakTime(), wireID))); // drift away from plane
        distance = std::pow(pos.X() - truth_e_position.X(), 2)+std::pow(pos.Z() - truth_e_position.Z(), 2);
        distance = sqrt(distance);
    }

    // write charge to respective distance cut bins
    for(int i = 0; i < 20; i++){
        double cut = 14.0*(i+1); // each bin signifies a cut applied at 14*bin
        if(distance<cut)
            charge_Z_dist.at(i)+=hit.Integral();
    }

}