#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"

#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Base.h"
#include "../NanoCORE/tqdm.h"

#include <iostream>
#include <iomanip>

#define SUM(vec) std::accumulate((vec).begin(), (vec).end(), 0);
#define SUM_GT(vec,num) std::accumulate((vec).begin(), (vec).end(), 0, [](float x,float y){return ((y > (num)) ? x+y : x); });
#define COUNT_GT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x > (num); });
#define COUNT_LT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x < (num); });

#define H1(name,nbins,low,high) TH1F *h_##name = new TH1F(#name,#name,nbins,low,high);

// #define DEBUG

struct debugger { template<typename T> debugger& operator , (const T& v) { cerr<<v<<" "; return *this; } } dbg;
#ifdef DEBUG
    #define debug(args...) do {cerr << #args << ": "; dbg,args; cerr << endl;} while(0)
#else
    #define debug(args...)
#endif

using namespace std;
using namespace tas;

// Add argument for genEventSumw (obtained using the function defined in doAll_wwz)
//
// Scaling factor (for MC) is xsec*lumi/GenEventSumw ---------> Need to add cross section information (maybe as a function in doAll_wwz)

int ScanChain(TChain *ch, TString process, TString year, double xsec, double GenEventSumw, bool isMC) {

    TFile* f1 = new TFile("output_"+process+".root", "RECREATE");

    // define kinematic variables for tree
    short sel_mask = 0;  // selection bitmask to keep track of cuts
    float pt_l_Z1 = 0.;  // pT of leading Z candidate lepton
    float pt_l_Z2 = 0.;  // pT of subleading Z candidate lepton
    float pt_l_W1 = 0.;  // pT of leading W candidate lepton
    float pt_l_W2 = 0.;  // pT of subleading W candidate lepton
    double evt_weight = 0.;  // per-event weight

    // Define the output tree
    TTree tree_out("tree","");
    tree_out.Branch("pt_l_Z1",&pt_l_Z1);
    tree_out.Branch("pt_l_Z2",&pt_l_Z2);
    tree_out.Branch("pt_l_W1",&pt_l_W1);
    tree_out.Branch("pt_l_W2",&pt_l_W2);
    tree_out.Branch("evt_weight",&evt_weight);

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();
    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);
    tqdm bar;

    // Luminosities
    double lumi = 0.;
    if ( year == "2016" ) lumi = 35.9;
    if ( year == "2017" ) lumi = 41.3;
    if ( year == "2018" ) lumi = 59.7;

    double factor;

    if ( isMC ) factor = xsec*lumi/GenEventSumw;
    if ( !isMC ) factor = 1.0;

    // set configuration parameters
    //gconf.year = 2017;

    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile *file = TFile::Open( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        TString filename(currentFile->GetTitle());

        tree->SetCacheSize(128*1024*1024);
        tree->SetCacheLearnEntries(100);

        nt.Init(tree);

        // Event loop
        for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

            nt.GetEntry(event);
            tree->LoadTree(event);

            nEventsTotal++;
            bar.progress(nEventsTotal, nEventsChain);

	    // Cuts, filling hists, branches, etc
	    
////////////// Pre-selection////////////////////////////////////////////////////////////////////////
            // Muons
            std::vector<int> pre_muon_idx;

            if (nt.nMuon() < 1) continue;

            for (int i = 0; i < nt.nMuon(); i++){
                if (!nt.Muon_looseId().at(i)) continue;
		if (!(std::abs(nt.Muon_eta().at(i)) < 2.4)) continue;
		if (!(nt.Muon_pt().at(i) > 10.)) continue;
                if (!(nt.Muon_pfIsoId().at(i) >= 1)) continue;
		pre_muon_idx.push_back(i);  // store the index of muons passing pre-selection cuts
	    } 
	  
            // Electrons
            std::vector<int> pre_ele_idx;

            if (nt.nElectron() < 1) continue;
	 
	    for (int j = 0; i < nt.nElectron(); j++){
		if (!nt.Electron_mvaFall17V2noIso_WPL().at(j)) continue;
		if (!(std::abs(nt.Electron_eta().at(j)) < 2.5)) continue;
		if (!(nt.Electron_pt().at(j) > 10.)) continue;
		if (!(nt.Electron_pfRelIso03_all().at(j) > 0.4)) continue;
		pre_ele_idx.push_back(j); // store electron index if it passes pre-selection
	    }

	    
            // Cut on the number of pre-selected muons and electrons. We want exactly 4
            if (pre_muon_idx.size() + pre_ele_idx.size() == 4){
		sel_mask |= 1UL << 0;  // set the first bit of the selmask to be true
	    }

	    
//////////////Lepton ID categorization//////////////////////////////////////////////////////////////
        
            std::vector<int> mu_Z_cand;
	    std::vector<int> mu_W_cand;
            if ( pre_muon_idx.size() > 0 ){
	    // Loop over muons looking for W and Z candidates
	    	for (int mu = 0; mu < pre_muon_idx.size(); mu++){
		    //Calculate IP_3D and sigma_IP_3D
		    float mu_sip_3D = nt.Muon_sip3d().at(pre_muon_idx[mu]);
		    float mu_ip_3D = nt.Muon_ip3d().at(pre_muon_idx[mu]);
		    bool mu_frac_sel = (std::abs(mu_ip_3D/mu_sip_3D) < 4);
		    if ( !(nt.Muon_mediumId().at(pre_muon_idx[mu]) && mu_frac_sel) ) continue;
		    if ( !(nt.Muon_pfRelIso04_all().at(pre_muon_idx[mu]) < 0.25) ) continue;
			//store index as Z candidate
			mu_Z_cand.push_back(pre_muon_idx[mu]);
		    if ( !(nt.Muon_pfRelIso04_all().at(pre_muon_idx[mu]) < 0.15) ) continue;
			//store index as W candidate
			mu_W_cand.push_back(pre_muon_idx[mu]);
	        }
	    }

	    std::vector<int> el_Z_cand;
	    std::vector<int> el_W_cand;
	    if ( pre_ele_idx.size() > 0 ){
		for (int el = 0; el < pre_ele_idx.size(); el++){
		    float el_sip_3D = nt.Electron_sip3d().at(pre_ele_idx[el]);
		    float el_ip_3D = nt.Electron_ip3d().at(pre_ele_idx[el]);
		    bool el_frac_sel = (std::abs(el_ip_3D/el_sip_3D) < 4);
		    if (!el_frac_sel) continue;
		    if ( nt.Electron_mvaFall17V2noIso_WPL().at(pre_ele_idx[el] ) el_Z_cand.push_back(pre_ele_idx[el]);
		    if ( nt.Electron_mvaFall17V2Iso_WP90().at(pre_ele_idx[el] ) el_W_cand.push_back(pre_ele_idx[el]);
		}
	    }

///////////////W and Z Candidate Selection///////////////////////////////////////////////////////
                       

            float m_Z = 91.2;
            std::vector<pair<float,int>> electron_cand;
	    std::vector<pair<float,int>> muon_cand;
            if ( el_Z_cand.size() > 1 ){  
	    	std::vector<pair<int,int>> electron_Z_cands = getCandPairs_el(el_Z_cand);
		electron_cand = get_el_resonance(electron_Z_cands);
            }
            if ( mu_Z_cand.size() > 1 ){
	    	std::vector<pair<int,int>> muon_Z_cands = getCandPairs_mu(mu_Z_cand);
		muon_cand = get_mu_resonance(muon_Z_cands);
            }

	    //Select the Z candidate closest to the Z mass
	    float m_Z_cand = 0.0;
	    std::vector<pair<int,int>> Z_cands_mu;
	    std::vector<pair<int,int>> Z_cands_el;
	    
	    if ( std::abs(muon_cand[0].first - m_Z) < std::abs(electron_cand[0].first - m_Z) ){
		 if ( std::abs(muon_cand[0].first - m_Z) > 10 ) continue;
	         Z_cands_mu.emplace_back(muon_Z_cands[muon_cand[0].second]);
	    } 
  
	    if ( std::abs(muon_cand[0].first - m_Z) > std::abs(electron_cand[0].first - m-Z) ){
		 if ( std::abs(electron_cand[0].first - m_Z) > 10 ) continue;
		 Z_cands_el.emplace_back(electron_Z_cands[electron_cand[0].second]);
	    }

            bool Z_mm = false;
	    bool Z_ee = false;

	    if ( Z_cands_mu.size() > 1 && Z_cands_el.size == 0 ){ 
		 Z_mm = true;
		 pt_l_Z1 = nt.Muon_pt().at(Z_cands_mu[0]);
		 pt_l_Z2 = nt.Muon_pt().at(Z_cands_mu[1]);
            }
	    if ( Z_cands_el.size() > 1 && Z_cands_mu.size == 0 ){ 
		 Z_ee = true;
		 pt_l_Z1 = nt.Electron_pt().at(Z_cands_el[0]);
		 pt_l_Z2 = nt.Electron_pt().at(Z_cands_el[1]);
	    }
	    if ( !(Z_mm) && !(Z_ee) ) continue;
	  
            // Now look for the W candidates

	    // Flags for same-flavor, opposite flavor
	    /*bool same_flav = false;
	    bool opp_flav = false;
            
	    std::vector<pair<int,int>> W_cands;
            int W_id_1;
	    int W_id_2;
	    std::vector<int> W_leptons_el;
	    std::vector<int> W_leptons_mu;

	    if ( Z_mm ){
	       W_cands = get_W_cands(false, Z_cands_mu, el_W_cand, mu_W_cand);
               if ( Muon_pt().at(Z_cands_mu[0].first) < 25 || Muon_pt().at(Z_cands_mu[0].second) < 15 ) continue; 
	       W_id_1 = W_cands[0].second;
	       W_id_2 = W_cands[1].second;
	       if ( W_id_1 + W_id_2 == 0 ){
	            same_flav = true;
		    if ( std::abs(W_id_1) == 13 ){
			 W_leptons_mu.push_back(W_cands[0].first);
			 W_leptons_mu.push_back(W_cands[1].first);
		    }
		    if ( std::abs(W_id_1) == 11 ){
			 W_leptons_el.push_back(W_cands[0].first);
			 W_leptons_el.push_back(W_cands[1].first);
		    }   
	       }
	       if ( std::abs(W_id_1 + W_id_2) == 2 ){
		    opp_flav = true;
		    if ( std::abs(W_id_1) > std::abs(W_id_2) ){
			 W_leptons_mu.push_back(W_cands[0].first);
			 W_leptons_el.push_back(W_cands[1].first);
		    }
		    else{
			 W_leptons_el.push_back(W_cands[0].first);
			 W_leptons_mu.push_back(W_cands[1].first);
		    }
	       }
	       else continue;
	    }
	
	    if ( Z_ee ){
	       W_cands = get_W_cands(true, Z_cands_el, el_W_cand, mu_W_cand);
	       if ( Electron_pt().at(Z_cands_el[0].first) < 25 || Electron_pt().at(Z_cands_el[0].second) < 15 ) continue;
	       W_id_1 = W_cands[0].second;
               W_id_2 = W_cands[1].second;
               if ( W_id_1 + W_id_2 == 0 ){
                    same_flav = true;
                    if ( std::abs(W_id_1) == 13 ){
                         W_leptons_mu.push_back(W_cands[0].first);
                         W_leptons_mu.push_back(W_cands[1].first);
                    }   
                    if ( std::abs(W_id_1) == 11 ){
                         W_leptons_el.push_back(W_cands[0].first);
                         W_leptons_el.push_back(W_cands[1].first);
                    }   
               }
               if ( std::abs(W_id_1 + W_id_2) == 2 ){
                    opp_flav = true;
                    if ( std::abs(W_id_1) > std::abs(W_id_2) ){
                         W_leptons_mu.push_back(W_cands[0].first);
                         W_leptons_el.push_back(W_cands[1].first);
                    }
                    else{
                         W_leptons_el.push_back(W_cands[0].first);
                         W_leptons_mu.push_back(W_cands[1].first);
                    }
               }
               else continue;
	    }

	    
	    // W-lepton pt requirements
	    if ( W_leptons_el.size() == 2 ){
	         pt_l_W1 = nt.Electron_pt().at(W_leptons_el[0]);
		 pt_l_W2 = nt.Electron_pt().at(W_leptons_el[1]);
		 if ( pt_l_W1 < 25 || pt_l_W2 < 15 ) continue;	 
	    }

	    if ( W_leptons_mu.size() == 2 ){
		 pt_l_W1 = nt.Muon_pt().at(W_leptons_mu[0]);
		 pt_l_W2 = nt.Muon_pt().at(W_leptons_mu[1]);
		 if ( pt_l_W1 < 25 || pt_l_W2 < 15 ) continue;
	    }
            
	    if ( W_leptons_el.size() == 1 && W_leptons_mu.size() == 1 ){
                 pt_l_W1 = std::max(nt.Muon_pt().at(W_leptons_mu[0]),nt.Electron_pt().at(W_leptons_el[0]));
		 pt_l_W2 = std::min(nt.Muon_pt().at(W_leptons_mu[0]),nt.Electron_pt().at(W_leptons_el[0]));
		 if ( pt_l_W1 < 25 || pt_l_W2 < 15 ) continue;
	    }
            */
            
///////////////b-tagged-jet veto////////////////////////////////////////////////////////////////////////////

	    	    

            float weight = genWeight();
	    evt_weight = weight*factor;

	    tree_out.Fill();

        } // Event loop

        delete file;


    } // File loop
    bar.finish();

    f1->Write();
    f1->Close();
    return 0;
}
