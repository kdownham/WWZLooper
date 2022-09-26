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
#include "../NanoCORE/WWZSelections.h"

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

    std::cout << "Process = " << process << " , year = " << year << " , cross section (fb) = " << xsec << " , genEventSumw = " << GenEventSumw << std::endl;

    TFile* f1 = new TFile("output/output_"+process+".root", "RECREATE");

    // define kinematic variables for tree
    short sel_mask = 0;  // selection bitmask to keep track of cuts
    float pt_l_Z1 = 0.;  // pT of leading Z candidate lepton
    float pt_l_Z2 = 0.;  // pT of subleading Z candidate lepton
    float pt_l_W1 = 0.;  // pT of leading W candidate lepton
    float pt_l_W2 = 0.;  // pT of subleading W candidate lepton
    double evt_weight = 0.;  // per-event weight
    float m_ll = 0.; // invariant mass of Z candidate leptons 
    bool opp_flav = false;
    bool same_flav = true;

    double sumOfEventWeights = 0.;
    double sumSQEventWeights = 0.;
    double sumW_el_EventWeights = 0.;
    double sumW_mu_EventWeights = 0.;
    double sumZ_mm_EventWeights = 0.;
    double sumZ_ee_EventWeights = 0.;
    double sumZEventWeights = 0.;
    double sum_sf_EventWeights = 0.;
    double sum_of_EventWeights = 0.;
    double sum_pre_sf_EventWeights = 0.;
    double sum_pre_of_EventWeights = 0.;
    double sumW_mm_EventWeights = 0.;
    double sumW_ee_EventWeights = 0.;
    double sumW_em_EventWeights = 0.;
    double nEvents_4mu = 0.;
    double nEvents_4el = 0.;
    double nEvents_3mu = 0.;
    double nEvents_2mu = 0.;
    double nEvents_1mu = 0.;
    double nEvents_Zee_WeWe = 0.;
    double nEvents_Zee_WeWm = 0.;
    double nEvents_Zee_WmWm = 0.;
    double nEvents_Zmm_WeWe = 0.;
    double nEvents_Zmm_WeWm = 0.;
    double nEvents_Zmm_WmWm = 0.;
    double nEvents_1pmu_1pel = 0.;
    double nEvents_2pmu_0el = 0.;
    double nEvents_0mu_2pel = 0.;

    bool debug = false;

    // Define the output tree
    TTree tree_out("evt","");
    tree_out.Branch("pt_l_Z1",&pt_l_Z1);
    tree_out.Branch("pt_l_Z2",&pt_l_Z2);
    tree_out.Branch("pt_l_W1",&pt_l_W1);
    tree_out.Branch("pt_l_W2",&pt_l_W2);
    tree_out.Branch("evt_weight",&evt_weight);
    tree_out.Branch("m_ll",&m_ll);
    tree_out.Branch("opp_flav",&opp_flav);
    tree_out.Branch("same_flav",&same_flav);

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
    if ( year == "2018" ) lumi = 136.9;			// 59.7 is actual value 

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

	    float weight = genWeight();
	    evt_weight = weight*factor;

	    // Cuts, filling hists, branches, etc
	    
////////////// Pre-selection////////////////////////////////////////////////////////////////////////
            // Muons
            std::vector<int> pre_muon_idx;

	    if (debug) std::cout << "debug 1" << endl;

            if (nt.nMuon() < 1) continue;

            for (int i = 0; i < nt.nMuon(); i++){
                if (!nt.Muon_looseId().at(i)) continue;
		if (!(std::abs(nt.Muon_eta().at(i)) < 2.4)) continue;
		if (!(nt.Muon_pt().at(i) > 10.)) continue;
                if (!(nt.Muon_pfIsoId().at(i) >= 1)) continue;
		pre_muon_idx.push_back(i);  // store the index of muons passing pre-selection cuts
	    } 
	  
	    if (debug) std::cout << "debug 2" << endl;

            // Electrons
            std::vector<int> pre_ele_idx;

            if (nt.nElectron() < 1) continue;
	 
	    for (int j = 0; j < nt.nElectron(); j++){
		if (!nt.Electron_mvaFall17V2noIso_WPL().at(j)) continue;
		if (!(std::abs(nt.Electron_eta().at(j)) < 2.5)) continue;
		if (!(nt.Electron_pt().at(j) > 10.)) continue;
		if (!(nt.Electron_pfRelIso03_all().at(j) < 0.4)) continue;
		pre_ele_idx.push_back(j); // store electron index if it passes pre-selection
	    }

	    if (debug) std::cout << "debug 3" << endl;	    

            // Cut on the number of pre-selected muons and electrons. We want exactly 4
            if (!(pre_muon_idx.size() + pre_ele_idx.size() == 4)){
		//sel_mask |= 1UL << 0;  // set the first bit of the selmask to be true
		continue;
	    }

	    if (pre_muon_idx.size() == 4) nEvents_4mu += evt_weight;
	    if (pre_muon_idx.size() == 3) nEvents_3mu += evt_weight;
	    if (pre_muon_idx.size() == 2) nEvents_2mu += evt_weight;
	    if (pre_muon_idx.size() == 1) nEvents_1mu += evt_weight;
	    if (pre_muon_idx.size() == 0) nEvents_4el += evt_weight;

	    if (debug) std::cout << "debug 4" << endl;
	    
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

	    if (debug) std::cout << "debug 5" << endl;

	    std::vector<int> el_Z_cand;
	    std::vector<int> el_W_cand;
	    if ( pre_ele_idx.size() > 0 ){
		for (int el = 0; el < pre_ele_idx.size(); el++){
		    float el_sip_3D = nt.Electron_sip3d().at(pre_ele_idx[el]);
		    float el_ip_3D = nt.Electron_ip3d().at(pre_ele_idx[el]);
		    bool el_frac_sel = (std::abs(el_ip_3D/el_sip_3D) < 4);
		    if (!el_frac_sel) continue;
		    if ( nt.Electron_mvaFall17V2noIso_WPL().at(pre_ele_idx[el]) ) el_Z_cand.push_back(pre_ele_idx[el]);
		    if ( nt.Electron_mvaFall17V2Iso_WP90().at(pre_ele_idx[el]) ) el_W_cand.push_back(pre_ele_idx[el]);
		}
	    }

	    if ( el_Z_cand.size() > 1 ){
		 if ( mu_W_cand.size() == 0) nEvents_Zee_WeWe += evt_weight;
		    // 4 el candidate event
		 if ( mu_W_cand.size() == 1 ) nEvents_Zee_WeWm += evt_weight;
		    // 3 el 1 mu candidate event
		 if ( mu_W_cand.size() == 2 ) nEvents_Zee_WmWm += evt_weight;
		    // 2 el 2 mu candidate event
	    }

	    if ( mu_Z_cand.size() > 1 ){
		 if ( el_W_cand.size() == 0 ) nEvents_Zmm_WmWm += evt_weight;
		    // 4 mu candidate event
	 	 if ( el_W_cand.size() == 1 ) nEvents_Zmm_WeWm += evt_weight;
		    // 3 mu candidate event
		 if ( el_W_cand.size() == 2 ) nEvents_Zmm_WeWe += evt_weight;
		    // 2 mu candidate event
	    }
	    
///////////////W and Z Candidate Selection///////////////////////////////////////////////////////
                       
	    if (debug) std::cout << "debug 6" << endl;

            float m_Z = 91.2;
            std::vector<pair<float,int>> electron_cand;
	    std::vector<pair<float,int>> muon_cand;
	    std::vector<pair<int,int>> electron_Z_cands;
            std::vector<pair<int,int>> muon_Z_cands;
            if ( el_Z_cand.size() > 1 ){  
	    	electron_Z_cands = getCandPairs_el(el_Z_cand);  // Stores indices of electron cands
		electron_cand = get_el_resonance(electron_Z_cands); // stores mass and index of electron candidates in electron_Z_cands
            }
            if ( mu_Z_cand.size() > 1 ){
	    	muon_Z_cands = getCandPairs_mu(mu_Z_cand);
		muon_cand = get_mu_resonance(muon_Z_cands);
            }

	    if (debug) std::cout << "debug 7" << endl;

	    //Select the Z candidate closest to the Z mass
	    float m_Z_cand = 0.0;
	    std::vector<pair<int,int>> Z_cands_mu;
	    std::vector<pair<int,int>> Z_cands_el;

	    
	    // Only works if both electron cand and muon cand are nonzero size
	    if ( muon_cand.size() > 0 && electron_cand.size() > 0 ){
	        nEvents_1pmu_1pel += evt_weight;
	    	if ( std::abs(muon_cand[0].first - m_Z) < std::abs(electron_cand[0].first - m_Z) ){
			 if ( std::abs(muon_cand[0].first - m_Z) > 10 ) continue;
	         	 Z_cands_mu.emplace_back(muon_Z_cands[muon_cand[0].second]);
	    	} 
  
	    	if ( std::abs(muon_cand[0].first - m_Z) > std::abs(electron_cand[0].first - m_Z) ){
		 	if ( std::abs(electron_cand[0].first - m_Z) > 10 ) continue;
		 	Z_cands_el.emplace_back(electron_Z_cands[electron_cand[0].second]);
	    	}
	    }

	    // Now do the same if we only have muons or electrons in the final state
	    if ( muon_cand.size() > 0 && electron_cand.size() == 0 ){
		nEvents_2pmu_0el += evt_weight;
		if ( std::abs(muon_cand[0].first - m_Z) > 10 ) continue;
		Z_cands_mu.emplace_back(muon_Z_cands[muon_cand[0].second]);
	    }
	   

	    if ( muon_cand.size() == 0 && electron_cand.size() > 0 ){
		nEvents_0mu_2pel += evt_weight;
                if ( std::abs(electron_cand[0].first - m_Z) > 10 ) continue;
                Z_cands_el.emplace_back(electron_Z_cands[electron_cand[0].second]);
            }


	    if (debug) std::cout << "debug 8" << endl;

            bool Z_mm = false;
	    bool Z_ee = false;

	    double mll = 0.;

	    

	    if ( Z_cands_mu.size() > 0 && Z_cands_el.size() == 0 ){ 
		 sumZ_mm_EventWeights += evt_weight;
		 Z_mm = true;
		 pt_l_Z1 = nt.Muon_pt().at(Z_cands_mu[0].first);
		 pt_l_Z2 = nt.Muon_pt().at(Z_cands_mu[0].second);
		 mll = (nt.Muon_p4().at(Z_cands_mu[0].first)+nt.Muon_p4().at(Z_cands_mu[0].second)).M();
            }
	    if ( Z_cands_el.size() > 0 && Z_cands_mu.size() == 0 ){ 
		 sumZ_ee_EventWeights += evt_weight;	 
		 Z_ee = true;
		 pt_l_Z1 = nt.Electron_pt().at(Z_cands_el[0].first);
		 pt_l_Z2 = nt.Electron_pt().at(Z_cands_el[0].second);
		 mll = (nt.Electron_p4().at(Z_cands_el[0].first)+nt.Electron_p4().at(Z_cands_el[0].second)).M();
	    }
	    if ( !(Z_mm) && !(Z_ee) ) continue;
	    
	    m_ll = mll;

	    sumZEventWeights += evt_weight;
	  
	    if (debug) std::cout << "debug 9" << endl;

            // Now look for the W candidates

	    // Flags for same-flavor, opposite flavor
	    bool same_flav = false;
	    bool opp_flav = false;
            
	    std::vector<pair<int,int>> W_cands;
            int W_id_1;
	    int W_id_2;
	    std::vector<int> W_leptons_el;
	    std::vector<int> W_leptons_mu;

	    if ( Z_mm ){
	       //sumZ_mm_EventWeights += evt_weight;
	       W_cands = get_W_cands(false, Z_cands_mu, el_W_cand, mu_W_cand);
               if ( Muon_pt().at(Z_cands_mu[0].first) < 25 || Muon_pt().at(Z_cands_mu[0].second) < 15 ) continue; 
	       if (W_cands.size() < 2) continue;
	       W_id_1 = W_cands[0].second;
	       W_id_2 = W_cands[1].second;
	       if ( W_id_1 + W_id_2 == 0 ){
	            same_flav = true;
		    sum_pre_sf_EventWeights += evt_weight;
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
		    sum_pre_of_EventWeights += evt_weight;
		    if ( std::abs(W_id_1) > std::abs(W_id_2) ){
			 W_leptons_mu.push_back(W_cands[0].first);
			 W_leptons_el.push_back(W_cands[1].first);
		    }
		    else{
			 W_leptons_el.push_back(W_cands[0].first);
			 W_leptons_mu.push_back(W_cands[1].first);
		    }
	       }
	       sumW_mu_EventWeights += evt_weight;
	       if ( !( (W_id_1 + W_id_2 == 0) || (std::abs(W_id_1+W_id_2) == 2) )) continue;
	    }
	
	    //sumW_mu_EventWeights += evt_weight;
	
	    if (debug) std::cout << "debug 10" << endl;

	    if ( Z_ee ){
	       //sumZ_ee_EventWeights += evt_weight;
	       W_cands = get_W_cands(true, Z_cands_el, el_W_cand, mu_W_cand);
	       if ( Electron_pt().at(Z_cands_el[0].first) < 25 || Electron_pt().at(Z_cands_el[0].second) < 15 ) continue;
	       if (debug) std::cout << "debug 10.1" << endl;
	       if (W_cands.size() < 2) continue;
	       W_id_1 = W_cands[0].second;
               W_id_2 = W_cands[1].second;
	       if (debug) std::cout << "debug 10.2" << endl;
               if ( W_id_1 + W_id_2 == 0 ){
                    same_flav = true;
		    sum_pre_sf_EventWeights += evt_weight;
                    if ( std::abs(W_id_1) == 13 ){
                         W_leptons_mu.push_back(W_cands[0].first);
                         W_leptons_mu.push_back(W_cands[1].first);
                    }   
                    if ( std::abs(W_id_1) == 11 ){
                         W_leptons_el.push_back(W_cands[0].first);
                         W_leptons_el.push_back(W_cands[1].first);
                    }   
               }
	       if (debug) std::cout << "debug 10.3" << endl;
               if ( std::abs(W_id_1 + W_id_2) == 2 ){
                    opp_flav = true;
		    sum_pre_of_EventWeights += evt_weight;
                    if ( std::abs(W_id_1) > std::abs(W_id_2) ){
                         W_leptons_mu.push_back(W_cands[0].first);
                         W_leptons_el.push_back(W_cands[1].first);
                    }
                    else{
                         W_leptons_el.push_back(W_cands[0].first);
                         W_leptons_mu.push_back(W_cands[1].first);
                    }
               }
	       sumW_el_EventWeights += evt_weight;
               if ( !( (W_id_1 + W_id_2 == 0) || (std::abs(W_id_1+W_id_2) == 2) )) continue;
	    }


	    //if (opp_flav) sum_pre_of_EventWeights += evt_weight;
            //if (same_flav) sum_pre_sf_EventWeights += evt_weight;

	    //sumW_el_EventWeights += evt_weight;

	    if (debug) std::cout << "debug 11" << endl;	    

	    // W-lepton pt requirements
	    
	    bool opposite_flavor = false;
	    bool same_flavor = false;
	    
	    if ( W_leptons_el.size() == 2 && W_leptons_mu.size() == 0 ){
		 same_flavor = true;
		 sumW_ee_EventWeights += evt_weight;
	         pt_l_W1 = nt.Electron_pt().at(W_leptons_el[0]);
		 pt_l_W2 = nt.Electron_pt().at(W_leptons_el[1]);
		 if ( pt_l_W1 < 25 || pt_l_W2 < 15 ) continue;	 
	    }

	    if ( W_leptons_mu.size() == 2 && W_leptons_el.size() == 0 ){
		 same_flavor = true;
		 sumW_mm_EventWeights += evt_weight;
		 pt_l_W1 = nt.Muon_pt().at(W_leptons_mu[0]);
		 pt_l_W2 = nt.Muon_pt().at(W_leptons_mu[1]);
		 if ( pt_l_W1 < 25 || pt_l_W2 < 15 ) continue;
	    }
            
	    if ( W_leptons_el.size() == 1 && W_leptons_mu.size() == 1 ){
		 opposite_flavor = true;
		 sumW_em_EventWeights += evt_weight;
                 pt_l_W1 = std::max(nt.Muon_pt().at(W_leptons_mu[0]),nt.Electron_pt().at(W_leptons_el[0]));
		 pt_l_W2 = std::min(nt.Muon_pt().at(W_leptons_mu[0]),nt.Electron_pt().at(W_leptons_el[0]));
		 if ( pt_l_W1 < 25 || pt_l_W2 < 15 ) continue;
	    }

	    opp_flav = opposite_flavor;
	    same_flav = same_flavor;
            
            if (debug) std::cout << "debug 12" << endl;
///////////////b-tagged-jet veto////////////////////////////////////////////////////////////////////////////

	    int nBjets = 0;
	    for (int jet_idx = 0; jet_idx < nt.nJet(); jet_idx++){
		 if ( nt.Jet_btagDeepB().at(jet_idx) > 0.1208 && nt.Jet_pt().at(jet_idx) > 20. ){
		      nBjets++;
		 }
	    }	  

	    if (nBjets > 0) continue; 

///////////////signal region definitions//////////////////////////////////////////////////////////////////// 

	    // Construct MT2 variable

	    if ( same_flavor ){
		 //if ( m_ll < 100. ){
		    // Impose MT2 cut
		 sum_sf_EventWeights += evt_weight;

		 //}
		 //else{
		 //   sum_sf_EventWeights += evt_weight;
		 //}	 
	    }

	    if ( opposite_flavor ){
		 sum_of_EventWeights += evt_weight;
	    }

	    
            //float weight = genWeight();
	    //evt_weight = weight*factor;
	    sumOfEventWeights += evt_weight;
	    sumSQEventWeights += std::pow(evt_weight,2.0);

	    tree_out.Fill();

        } // Event loop
        delete file;


    } // File loop
    bar.finish();

    std::cout << "Number of events passing selection = " << sumOfEventWeights << " +/- " << std::sqrt(sumSQEventWeights) << endl;
    //std::cout << "Number of events passing Z selection = " << sumZEventWeights << endl;
    //std::cout << "Number of events in Z to mu mu category = " << sumZ_mm_EventWeights << endl;
    //std::cout << "Number of events in Z to el el category = " << sumZ_ee_EventWeights << endl;
    //std::cout << "N events in Z to mu mu category passing W selection = " << sumW_mu_EventWeights << endl;
    //std::cout << "N events in Z to el el category passing W selection = " << sumW_el_EventWeights << endl;
    //std::cout << "Number of events in pre same-flavor category = " << sum_pre_sf_EventWeights << endl;
    //std::cout << "Number of events in pre opposite-flavor category = " << sum_pre_of_EventWeights << endl;
    std::cout << "Number of events in same-flavor category = " << sum_sf_EventWeights  << endl;
    std::cout << "Number of events in opposite-flavor category = " << sum_of_EventWeights << endl;
    //std::cout << "Number of events in mu-mu category = " << sumW_mm_EventWeights << endl;
    //std::cout << "Number of events in el-el category = " << sumW_ee_EventWeights << endl;
    //std::cout << "Number of events in el-mu category = " << sumW_em_EventWeights << endl;
    //std::cout << "N events 4mu = " << nEvents_4mu << endl;
    //std::cout << "N events 3mu = " << nEvents_3mu << endl;
    //std::cout << "N events 2mu = " << nEvents_2mu << endl;
    //std::cout << "N events 1mu = " << nEvents_1mu << endl;
    //std::cout << "N events 0mu = " << nEvents_4el << endl;
    //std::cout << "N events Zee WeWe = " << nEvents_Zee_WeWe << endl;
    //std::cout << "N events Zee WeWm = " << nEvents_Zee_WeWm << endl;
    //std::cout << "N events Zee WmWm = " << nEvents_Zee_WmWm << endl;
    //std::cout << "N events Zmm WeWe = " << nEvents_Zmm_WeWe << endl;
    //std::cout << "N events Zmm WeWm = " << nEvents_Zmm_WeWm << endl;
    //std::cout << "N events Zmm WmWm = " << nEvents_Zmm_WmWm << endl; 
    //std::cout << "N events 2+ muons, 0 electrons = " << nEvents_2pmu_0el << endl;
    //std::cout << "N events 1+ muons, 1+ electrons = " << nEvents_1pmu_1pel << endl;
    //std::cout << "N events 0 muons, 2+ electrons = " << nEvents_0mu_2pel << endl;

    f1->Write();
    f1->Close();
    return 0;
}
