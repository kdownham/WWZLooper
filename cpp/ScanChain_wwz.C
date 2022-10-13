#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include <math.h>
//#include "Math/GenVector/Boost.h"
#include "Math/Boost.h"
#include "../rooutil/rooutil.h"
#include "../rooutil/calc.h"

// For .. gconf???
#include "../NanoCORE/SSSelections.cc"
#include "../NanoCORE/MetSelections.cc"

#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Base.h"
#include "../NanoCORE/tqdm.h"
#include "../NanoCORE/WWZSelections.h"
#include "../NanoCORE/lester_mt2_bisect.h"

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

    TFile* f1 = new TFile("output/output_"+process+"_"+year+".root", "RECREATE");

    // define kinematic variables for tree
    short sel_mask = 0;  // selection bitmask to keep track of cuts
    double pt_l_Z1;  // pT of leading Z candidate lepton
    double pt_l_Z2;  // pT of subleading Z candidate lepton
    double pt_l_W1;  // pT of leading W candidate lepton
    double pt_l_W2;  // pT of subleading W candidate lepton
    double evt_weight;  // per-event weight
    double m_ll; // invariant mass of Z candidate leptons 
    bool opp_flav;
    bool same_flav;
    double MT2;
    double MET_px;
    double MET_py;
    double MET;

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
    tree_out.Branch("MT2",&MT2);
    tree_out.Branch("opp_flav",&opp_flav);
    tree_out.Branch("same_flav",&same_flav);
    tree_out.Branch("MET_px",&MET_px);
    tree_out.Branch("MET_py",&MET_py);
    tree_out.Branch("MET",&MET);

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
    if ( year == "2018" ) lumi = 59.8;         //136.9;			// 59.7 is actual value 

    double factor;

    if ( isMC ) factor = xsec*lumi/GenEventSumw;
    if ( !isMC ) factor = 1.0;

    // set configuration parameters
    //gconf.year = 2018;

    //float bWPloose = gconf.WP_DeepFlav_loose;
    //std::cout << "Deep Flavor loose working point = " << bWPloose << endl;

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
	
	    bool isEvent = false;

	    //if ( nt.run() == 1 && nt.luminosityBlock() == 1801 && nt.event() == 1800002 ){
	//	 isEvent = true;
	  //  }

	    bool triggerReqs = false;

	    if ( year == "2016" ) triggerReqs = (nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL() || nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ());
	    if ( year == "2017" ) triggerReqs = (nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() || nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() | nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ());
	    if ( year == "2018" ) triggerReqs = (nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() || nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() || nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ());

	    if (!triggerReqs) continue;

	    if ( isEvent ) std::cout << "Passed trigger requirement!" << endl;

            nEventsTotal++;
            bar.progress(nEventsTotal, nEventsChain);

	    float weight = genWeight();
	    evt_weight = weight*factor;

	    // Cuts, filling hists, branches, etc
	    
////////////// Pre-selection////////////////////////////////////////////////////////////////////////
            // Muons
            std::vector<int> lep_ids;
	    std::vector<int> lep_idxs;
	    std::vector<LorentzVector> lep_p4s;
	    std::tuple<std::vector<int>,std::vector<int>,std::vector<LorentzVector>> lepton_p4s;

	    if (debug) std::cout << "debug 1" << endl;


            for (int i = 0; i < nt.nMuon(); i++){
                if (!nt.Muon_looseId().at(i)) continue;
		if (isEvent) std::cout << "Passed muon loose ID" << endl;
		if (!(std::abs(nt.Muon_eta().at(i)) < 2.4)) continue;
		if (isEvent) std::cout << "Passed muon eta" << endl;
		if (!(nt.Muon_pt().at(i) > 10.)) continue;
		if (isEvent) std::cout << "Passed muon pt" << endl;
                if (!(nt.Muon_pfIsoId().at(i) >= 1)) continue;
		if (isEvent) std::cout << "Passed muon pfIso" << endl;
		if (isEvent) std::cout << "Muon pT = " << Muon_pt()[i] << endl;
		if (isEvent) std::cout << "Muon eta = " << Muon_eta()[i] << endl;
		if (isEvent) std::cout << "Muon phi = " << Muon_phi()[i] << endl;
		if (isEvent) std::cout << "Muon pdgId = " << Muon_pdgId()[i] << endl;
		lep_idxs.push_back(i);
		lep_ids.push_back(Muon_pdgId()[i]);
		lep_p4s.push_back(Muon_p4()[i]);	
	    } 
	  
	    if (debug) std::cout << "debug 2" << endl;

            // Electrons
	 
	    for (int j = 0; j < nt.nElectron(); j++){
		if (!nt.Electron_mvaFall17V2noIso_WPL().at(j)) continue;
		if (isEvent) std::cout << "Passed electron noIso WPL" << endl;
		if (!(std::abs(nt.Electron_eta().at(j)) < 2.5)) continue;
		if (isEvent) std::cout << "Passed electron eta" << endl;
		if (!(nt.Electron_pt().at(j) > 10.)) continue;
		if (isEvent) std::cout << "Passed electron pt" << endl;
		if (!(nt.Electron_pfRelIso03_all().at(j) < 0.4)) continue;
		if (isEvent) std::cout << "Passed electron pfRelIso" << endl;
		if (isEvent) std::cout << "Electron pT = " << Electron_pt()[j] << endl;
                if (isEvent) std::cout << "Electron eta = " << Electron_eta()[j] << endl;
                if (isEvent) std::cout << "Electron phi = " << Electron_phi()[j] << endl;
                if (isEvent) std::cout << "Electron pdgId = " << Electron_pdgId()[j] << endl;
		lep_idxs.push_back(j);
                lep_ids.push_back(Electron_pdgId()[j]);
                lep_p4s.push_back(Electron_p4()[j]);
	    }

	    std::get<0>(lepton_p4s) = lep_idxs;
	    std::get<1>(lepton_p4s) = lep_ids;
	    std::get<2>(lepton_p4s) = lep_p4s;

	    if (debug) std::cout << "debug 3" << endl;	    

	    //int size = std::tuple_size<decltype(lepton_p4s)>::value;
	    int size = lep_idxs.size();
	    if (isEvent) std::cout << "Number of common veto ID leptons = " << size << endl;

            // Cut on the number of pre-selected muons and electrons. We want exactly 4
            if (!( size == 4) ){
		//sel_mask |= 1UL << 0;  // set the first bit of the selmask to be true
		continue;
	    }

	    if ( isEvent ) std::cout << "Passed 4-lepton requirement!" << endl;

	    if (debug) std::cout << "debug 4" << endl;

	    auto [ has_sfos, idx1, idx2, m_leps ] = RooUtil::Calc::pickZcandidateIdxs(std::get<1>(lepton_p4s),std::get<2>(lepton_p4s));

	    if (isEvent){
		std::cout << "has sfos = " << has_sfos << endl;
		std::cout << "Index 1 = " << idx1 << endl;
		std::cout << "Index 2 = " << idx2 << endl;
		std::cout << "mass = " << m_leps << endl;
	    }

	    if (!has_sfos) continue;

	    if (isEvent) std::cout << "Has sf os pair" << endl;

            // Define vectors containing indices for lepton Z and W candidates	   
            // Maybe do this separately for electrons and Muons so I can keep track of the indices......... 
            // idxs correspond to the index of the lepton_p4 vector that makes up the Z candidate
            bool Z_ee = false;
	    bool Z_mm = false;
	    std::vector<int> el_Z_cand;
	    std::vector<int> mu_Z_cand;
	    std::vector<int> el_W_cand;
	    std::vector<int> mu_W_cand;
	    if ( std::abs(std::get<1>(lepton_p4s)[idx1]) == 11 ){
		el_Z_cand.push_back(std::get<0>(lepton_p4s)[idx1]);
		el_Z_cand.push_back(std::get<0>(lepton_p4s)[idx2]);	
		Z_ee = true;
	    }
	    if ( std::abs(std::get<1>(lepton_p4s)[idx1]) == 13 ){
		mu_Z_cand.push_back(std::get<0>(lepton_p4s)[idx1]);
                mu_Z_cand.push_back(std::get<0>(lepton_p4s)[idx2]);
		Z_mm = true;
	    }

            // Indices in respective bins (muon or electron) 
	    if (isEvent) std::cout << "Z lepton indices = " << std::get<0>(lepton_p4s)[idx1] <<  " , " << std::get<0>(lepton_p4s)[idx2] << endl;
	    if (isEvent) std::cout << "Z lepton ids = " << std::get<1>(lepton_p4s)[idx1] << " , " << std::get<1>(lepton_p4s)[idx2] << endl;

	    std::vector<pair<int,int>> Z_cands;
	    Z_cands.emplace_back(std::get<0>(lepton_p4s)[idx1],std::get<1>(lepton_p4s)[idx1]);
	    Z_cands.emplace_back(std::get<0>(lepton_p4s)[idx2],std::get<1>(lepton_p4s)[idx2]);

	    // Now we have our Z candidates. Time to find the W candidates.
	    if ( Z_mm ){
	         // Loop over the stored indices
	         if (isEvent) std::cout << "Z_mm = true" << endl;
	         for ( int i = 0; i < lep_idxs.size(); i++ ){
		       int index_mm = std::get<0>(lepton_p4s)[i];
		       int flav_mm = std::get<1>(lepton_p4s)[i];
		       // If you have the same index and flavor as a Z candidate, skip to the next lepton
		       if ( index_mm == Z_cands[0].first && flav_mm == Z_cands[0].second ) continue;
		       if ( index_mm == Z_cands[1].first && flav_mm == Z_cands[1].second ) continue;
		       if ( std::abs(flav_mm) == 11 ) el_W_cand.push_back(std::get<0>(lepton_p4s)[i]);
		       if ( std::abs(flav_mm) == 13 ) mu_W_cand.push_back(std::get<0>(lepton_p4s)[i]);
		 }
	    }
	    if ( Z_ee ){
		 if (isEvent) std::cout << "Z_ee = true" << endl;
		 for ( int j = 0; j < lep_idxs.size(); j++ ){
		       int index_ee = std::get<0>(lepton_p4s)[j];
                       int flav_ee = std::get<1>(lepton_p4s)[j];
		       if ( index_ee == Z_cands[0].first && flav_ee == Z_cands[0].second ) continue;
                       if ( index_ee == Z_cands[1].first && flav_ee == Z_cands[1].second ) continue;
                       if ( std::abs(flav_ee) == 11 ) el_W_cand.push_back(std::get<0>(lepton_p4s)[j]);
                       if ( std::abs(flav_ee) == 13 ) mu_W_cand.push_back(std::get<0>(lepton_p4s)[j]);      
		 }
	    } 
	    
	    
//////////////Lepton ID categorization//////////////////////////////////////////////////////////////
           

	    if ( mu_Z_cand.size() > 0 ){
		 if (isEvent) std::cout << "Number of mu Z candidates = " << mu_Z_cand.size() << endl;
		 int nmuZ = 0;
		 for ( int mu_Z = 0; mu_Z < mu_Z_cand.size(); mu_Z++ ){
		       //if ( !Muon_mediumId()[mu_Z_cand[mu_Z]] ) continue;
		       //if (isEvent) std::cout << "Passed muon medium id" << endl;
		       //if (isEvent) std::cout << "Muon sip3d = " << Muon_sip3d()[mu_Z_cand[mu_Z]] << endl;
		       if ( !(std::abs(Muon_ip3d()[mu_Z_cand[mu_Z]]/Muon_sip3d()[mu_Z_cand[mu_Z]]) < 4) ) continue;
		       if (isEvent) std::cout << "Passed muon sip3d cut" << endl;
		       //if ( !(Muon_pfRelIso04_all()[mu_Z_cand[mu_Z]] < 0.25 )) continue;
		       //if (isEvent) std::cout << "Passed muon pfRelIso cut" << endl;
		       //if ( !(Muon_pfIsoId()[mu_Z_cand[mu_Z]] >= 2) ) continue;
		       if ( !(Muon_pfIsoId()[mu_Z_cand[mu_Z]] >= 1) ) continue;
		       nmuZ++;
		 }
		 if (isEvent) std::cout << "nmuZ = " << nmuZ << " , number of Z cand mu's = " << mu_Z_cand.size() << endl;
		 if ( nmuZ < mu_Z_cand.size() ) continue;
	    }

	    if (isEvent) std::cout << "Passes mu Z cand cuts" << endl;

	    if ( mu_W_cand.size() > 0 ){
		 if (isEvent) std::cout << "Mu cand W size > 0" << endl;
		 int nmuW = 0;
		 for ( int mu_W = 0; mu_W < mu_W_cand.size(); mu_W++ ){
		       //if ( !Muon_mediumId()[mu_W_cand[mu_W]] ) continue;
		       //if (isEvent) std::cout << "Passed muon medium ID" << endl;
                       if ( !(std::abs(Muon_ip3d()[mu_W_cand[mu_W]]/Muon_sip3d()[mu_W_cand[mu_W]]) < 4) ) continue;
		       if (isEvent) std::cout << "Passed muon sip3d cut" << endl;
                       //if ( !(Muon_pfRelIso04_all()[mu_W_cand[mu_W]] < 0.15 )) continue;
		       //if (isEvent) std::cout << "Passed muon pfRelIso cut" << endl;
                       //if ( !(Muon_pfIsoId()[mu_W_cand[mu_W]] >= 3) ) continue;
		       if ( !(Muon_pfIsoId()[mu_W_cand[mu_W]] >= 1) ) continue;
		       nmuW++;
		 }
		 if (isEvent) std::cout << "nmuW = " << nmuW << " , number of W cand mu's = " << mu_W_cand.size() << endl;
		 if ( nmuW < mu_W_cand.size() ) continue;
	    }
 
	    if (isEvent) std::cout << "Passes mu W cand cuts" << endl;	

	    if (debug) std::cout << "debug 5" << endl;

	    if ( el_Z_cand.size() > 0 ){
		 int nelZ = 0;
		 for ( int el_Z = 0; el_Z < el_Z_cand.size(); el_Z++ ){
		       if (!(std::abs(Electron_ip3d()[el_Z_cand[el_Z]]/Electron_sip3d()[el_Z_cand[el_Z]]) < 4) ) continue;
		       //if (!Electron_mvaFall17V2noIso_WPL()[el_Z_cand[el_Z]]) continue;
		       //if (!(Electron_pfRelIso03_all()[el_Z_cand[el_Z]] < 0.2)) continue;
		       nelZ++;
		 }
		 if ( nelZ < el_Z_cand.size() ) continue;
	    }

	    if (isEvent) std::cout << "Passes el Z cand cuts" << endl;

	    if ( el_W_cand.size() > 0 ){
		 if (isEvent) std::cout << "El cand W size > 0" << endl;
		 int nelW = 0;
		 for ( int el_W = 0; el_W < el_W_cand.size(); el_W++ ){
                       if (!(std::abs(Electron_ip3d()[el_W_cand[el_W]]/Electron_sip3d()[el_W_cand[el_W]]) < 4) ) continue;
		       if (isEvent) std::cout << "Electron passed sip cut!" << endl;
		       //if (!(Electron_pfRelIso03_all()[el_W_cand[el_W]] < 0.2)) continue;
                       //if (!Electron_mvaFall17V2Iso_WP90()[el_W_cand[el_W]]) continue;
                       nelW++;
                 }
		 if (isEvent) std::cout << "nelW = " << nelW << " , number of W cand el's = " << el_W_cand.size() << endl;
                 if ( nelW < el_W_cand.size() ) continue;
	    }


	    if (isEvent) std::cout << "Passes el W cand cuts" << endl;

	    
///////////////W and Z Candidate Selection///////////////////////////////////////////////////////
                       

	    bool Wee = ( el_W_cand.size() == 2 && mu_W_cand.size() == 0 );
	    bool Wmm = ( mu_W_cand.size() == 2 && el_W_cand.size() == 0 );
	    bool Wem = ( mu_W_cand.size() == 1 && el_W_cand.size() == 1 );

	    if (isEvent) std::cout << "Number of W cand el's = " << el_W_cand.size() << " , Number of W cand mu's = " << mu_W_cand.size() << endl;


	    float mll;

	    //double px_l_Z1, py_l_Z1, px_l_Z2, py_l_Z2;
	    LorentzVector p4_l_Z1, p4_l_Z2;

	    if ( Z_mm ){ 
		 pt_l_Z1 = nt.Muon_pt().at(mu_Z_cand[0]);
		 p4_l_Z1 = nt.Muon_p4().at(mu_Z_cand[0]);
		 pt_l_Z2 = nt.Muon_pt().at(mu_Z_cand[1]);
		 p4_l_Z2 = nt.Muon_p4().at(mu_Z_cand[1]);
		 mll = (nt.Muon_p4().at(mu_Z_cand[0])+nt.Muon_p4().at(mu_Z_cand[1])).M();
		 if (isEvent) std::cout << "m_ll = " << mll << " , m_leps = " << m_leps << endl;
            }
	    if ( Z_ee ){ 
		 pt_l_Z1 = nt.Electron_pt().at(el_Z_cand[0]);
                 p4_l_Z1 = nt.Electron_p4().at(el_Z_cand[0]);
		 pt_l_Z2 = nt.Electron_pt().at(el_Z_cand[1]);
                 p4_l_Z2 = nt.Electron_p4().at(el_Z_cand[1]);
		 mll = (nt.Electron_p4().at(el_Z_cand[0])+nt.Electron_p4().at(el_Z_cand[1])).M();
		 if (isEvent) std::cout << "m_ll = " << mll << " , m_leps = " << m_leps << endl;
	    }
	    if ( !(Z_mm) && !(Z_ee) ) continue;


	    if ( isEvent ) std::cout << "Passed Z candidate requirement!" << endl;
	    
	    m_ll = m_leps;

	    if ( std::abs(m_ll - 91.2) > 10. ) continue;

	    if ( isEvent ) std::cout << "Passed Z candidate within 10 GeV of Z mass requirement!" << endl;

	    sumZEventWeights += evt_weight;
	  
	    if (debug) std::cout << "debug 9" << endl;

            // Now look for the W candidates
            if ( !(Wee || Wmm || Wem) ) continue;
	    if (isEvent) std::cout << "W's go to el or mu!" << endl;

	    LorentzVector p4_l_W1;
            LorentzVector p4_l_W2;
	    
	    if (Wee){
		if (isEvent) std::cout << "WWee" << endl;
		if ( Electron_pt()[el_W_cand[0]] < 25. || Electron_pt()[el_W_cand[1]] < 15. ) continue;
		if (isEvent) std::cout << "W goes to electrons!" << endl;
		p4_l_W1 = Electron_p4()[el_W_cand[0]];
                p4_l_W2 = Electron_p4()[el_W_cand[1]];
                pt_l_W1 = Electron_pt()[el_W_cand[0]];
                pt_l_W2 = Electron_pt()[el_W_cand[1]];
	    }
            
	    if (Wmm){
		if (isEvent) std::cout << "WWmumu" << endl;
		if ( Muon_pt()[mu_W_cand[0]] < 25. || Muon_pt()[mu_W_cand[1]] < 15. ) continue;
		if (isEvent) std::cout << "W goes to electrons!" << endl;
		p4_l_W1 = Muon_p4()[mu_W_cand[0]];
		p4_l_W2 = Muon_p4()[mu_W_cand[1]];
		pt_l_W1 = Muon_pt()[mu_W_cand[0]];
		pt_l_W2 = Muon_pt()[mu_W_cand[1]];
	    }	

	    if (Wem){
		if (isEvent) std::cout << "WWemu" << endl;
		if (isEvent) std::cout << "pT of electron W cand = " << Electron_pt()[el_W_cand[0]] << endl;
		if (isEvent) std::cout << "pT of muon W cand = " << Muon_pt()[mu_W_cand[0]] << endl;
		double pt_l_W1 = std::max(Electron_pt()[el_W_cand[0]],Muon_pt()[mu_W_cand[0]]);
		double pt_l_W2 = std::min(Electron_pt()[el_W_cand[0]],Muon_pt()[mu_W_cand[0]]);
		if ( pt_l_W1 < 25. || pt_l_W2 < 15. ) continue;
		if (isEvent) std::cout << "Opposite flavor channel!" << endl;
		if ( Electron_pt()[el_W_cand[0]] > Muon_pt()[mu_W_cand[0]] ){
		     p4_l_W1 = Electron_p4()[el_W_cand[0]];
		     p4_l_W2 = Muon_p4()[mu_W_cand[0]];   
		}
		if ( Electron_pt()[el_W_cand[0]] < Muon_pt()[mu_W_cand[0]] ){
		     p4_l_W2 = Electron_p4()[el_W_cand[0]];
                     p4_l_W1 = Muon_p4()[mu_W_cand[0]];
		}
	    }
	    //sumW_mu_EventWeights += evt_weight;
	    
	    double m_Wcands = (p4_l_W1+p4_l_W2).M();
	    if (isEvent){
		std::cout << "Px of W1 lepton = " << p4_l_W1.Px() << endl;
		std::cout << "Py of W2 lepton = " << p4_l_W2.Py() << endl;
	    }
	
	    if (debug) std::cout << "debug 10" << endl;

            bool opp_flav = Wem;
	    bool same_flav = (Wee || Wmm);

	    if (opp_flav) sum_pre_of_EventWeights += evt_weight;
            if (same_flav) sum_pre_sf_EventWeights += evt_weight;


            
            if (debug) std::cout << "debug 12" << endl;
///////////////b-tagged-jet veto////////////////////////////////////////////////////////////////////////////
	    
	    int nBjets = 0;
	    for (int jet_idx = 0; jet_idx < nt.nJet(); jet_idx++){
		 /*if (isEvent){
		     std::cout << "Jet pT = " << nt.Jet_pt()[jet_idx] << endl;
 		     std::cout << "Jet btag ID = " << nt.Jet_btagDeepFlavB()[jet_idx] << endl;
		     std::cout << "Jet eta = " << nt.Jet_eta()[jet_idx] << endl;
		     std::cout << "Jet phi = " << nt.Jet_phi()[jet_idx] << endl;
		 }*/
		 // Overlap removal
		 if ( ROOT::Math::VectorUtil::DeltaR(nt.Jet_p4()[jet_idx],p4_l_Z1) < 0.4 ) continue;
		 if ( ROOT::Math::VectorUtil::DeltaR(nt.Jet_p4()[jet_idx],p4_l_Z2) < 0.4 ) continue;
		 if ( ROOT::Math::VectorUtil::DeltaR(nt.Jet_p4()[jet_idx],p4_l_W1) < 0.4 ) continue;
		 if ( ROOT::Math::VectorUtil::DeltaR(nt.Jet_p4()[jet_idx],p4_l_W2) < 0.4 ) continue;
		
		 /*if (isEvent){
		     std::cout << "Jet pT (post overlap removal) = " << nt.Jet_pt()[jet_idx] << endl;
		     std::cout << "Jet btag ID (post overlap removal) = " << nt.Jet_btagDeepFlavB()[jet_idx] << endl;
		     std::cout << "Jet eta (post overlap removal) = " << nt.Jet_eta()[jet_idx] << endl;
                     std::cout << "Jet phi (post overlap removal) = " << nt.Jet_phi()[jet_idx] << endl;
		 }*/	
	
		 if ( nt.Jet_btagDeepB().at(jet_idx) > 0.1208 && nt.Jet_p4()[jet_idx].Pt() > 20. && std::abs(nt.Jet_p4()[jet_idx].Eta()) < 2.4 ){
		      nBjets++;
		 }
	    }	  
	    if (isEvent) std::cout << "Number of b jets = "<< nBjets << endl;
	    if (nBjets != 0) continue;
	    if (isEvent) std::cout << "b-tagged jet requirement passed!" << endl;
	    

///////////////signal region definitions//////////////////////////////////////////////////////////////////// 
	    //TLorentzVector p4_lepW1(p4_l_W1.E(),p4_l_W1.Px(),p4_l_W1.Py(),p4_l_W1.Pz());
	    //TLorentzVector p4_lepW2(p4_l_W2.E(),p4_l_W2.Px(),p4_l_W2.Py(),p4_l_W2.Pz());

	    // Construct MT2 variable
	    double mVisA = 0.; // Mass of visible object on side A
	    double mVisB = 0.; // Mass of visible object on side B
	    
	    double chiA = 0.;
	    double chiB = 0.;

	    double desiredPrecisionOnMt2 = 0.;

	    LorentzVector MET_p4(nt.MET_pt(),0.,nt.MET_phi(),0.);

	    LorentzVector rest_WW = p4_l_W1 + p4_l_W2 + MET_p4;
	    auto beta_from_miss_reverse(-rest_WW.BoostToCM());
	    //TVector3 beta_from_miss(-beta_from_miss_reverse.X(),-beta_from_miss_reverse.Y(),-beta_from_miss_reverse.Z());
	    auto beta_from_miss = -beta_from_miss_reverse;

	    ROOT::Math::Boost xf(beta_from_miss.X(),beta_from_miss.Y(),beta_from_miss.Z());

	    LorentzVector boosted_p4_l_W1 = xf*p4_l_W1;
	    LorentzVector boosted_p4_l_W2 = xf*p4_l_W2;
	    MET_p4 = xf*MET_p4;

	    if (isEvent){

	        std::cout << "Inputs into MT2: " << endl;
		std::cout << "boosted_p4_l_W1.Px = " << boosted_p4_l_W1.Px() << endl;
		std::cout << "boosted_p4_l_W1.Py = " << boosted_p4_l_W1.Py() << endl;
		std::cout << "boosted_p4_l_W2.Px = " << boosted_p4_l_W2.Px() << endl;
                std::cout << "boosted_p4_l_W2.Py = " << boosted_p4_l_W2.Py() << endl;
		std::cout << "MET Px = " << MET_p4.Px() << endl;
		std::cout << "MET Py = " << MET_p4.Py() << endl;

	    }

	    asymm_mt2_lester_bisect::disableCopyrightMessage();
	    MT2 = asymm_mt2_lester_bisect::get_mT2( mVisA, boosted_p4_l_W1.Px(), boosted_p4_l_W1.Py(), mVisB, boosted_p4_l_W2.Px(), boosted_p4_l_W2.Py(), MET_p4.Px(), MET_p4.Py(), chiA, chiB, desiredPrecisionOnMt2);

	    // Construct pT_4L
	    double lep_px_sum = p4_l_Z1.Px() + p4_l_Z2.Px() + p4_l_W1.Px() + p4_l_W2.Px();
	    double lep_py_sum = p4_l_Z1.Py() + p4_l_Z2.Py() + p4_l_W1.Py() + p4_l_W2.Py();

	    double pt_4l = std::sqrt(std::pow(lep_px_sum,2.0)+std::pow(lep_py_sum,2.0));

	    bool sf_pass = false;
	    bool of_pass = false;
	    /*
	    if ( same_flav ){
		// Bin A
		if ( nt.MET_pt() > 120. ) sf_pass = true;
		if ( 70. < nt.MET_pt() < 120. ){
		     // Bin B
		     if ( pt_4l > 70. ) sf_pass = true;
		     // Bin C
		     if ( 40. < pt_4l < 70. ) sf_pass = true;
		     if ( pt_4l <= 40. ) continue;
		}
		if ( nt.MET_pt() <= 70. ) continue;
		      sf_pass = true;
	    }
	    */

	    //if (isEvent)

	    if ( same_flav ){
		 continue;
	    }

	    if (isEvent){
		std::cout << "MT2 = " << MT2 << endl;
		std::cout << "m_Wcands = " << m_Wcands << endl;
	    }
	    if ( opp_flav ){
		 if ( m_Wcands < 100. && MT2 < 25.) continue;
		 
		 if ( m_Wcands < 100. && MT2 > 25.) of_pass = true;
		      
                 if ( m_Wcands >= 100. ){
		      of_pass = true;
                 }
	    }

	
            //if ( sf_pass ) sum_sf_EventWeights += evt_weight;
	    if ( of_pass ) sum_of_EventWeights += evt_weight;

	    if (isEvent) std::cout << "Passed full selection!" << endl;
            //float weight = genWeight();
	    //evt_weight = weight*factor;
	    sumOfEventWeights += evt_weight;
	    sumSQEventWeights += std::pow(evt_weight,2.0);

	    std::cout << nt.run() << ":" << nt.luminosityBlock() << ":" << nt.event() << endl;

	    tree_out.Fill();

        } // Event loop
        delete file;


    } // File loop
    bar.finish();

    std::cout << "Number of events passing selection = " << sumOfEventWeights << " +/- " << std::sqrt(sumSQEventWeights) << endl;
    std::cout << "Number of events (scaled to full Run 2) passing selection = " << (136.9/59.8)*sumOfEventWeights << endl;
    //std::cout << "Number of events passing Z selection = " << sumZEventWeights << endl;
    //std::cout << "Number of events in Z to mu mu category = " << sumZ_mm_EventWeights << endl;
    //std::cout << "Number of events in Z to el el category = " << sumZ_ee_EventWeights << endl;
    //std::cout << "N events in Z to mu mu category passing W selection = " << sumW_mu_EventWeights << endl;
    //std::cout << "N events in Z to el el category passing W selection = " << sumW_el_EventWeights << endl;
    //std::cout << "Number of events in pre same-flavor category = " << sum_pre_sf_EventWeights << endl;
    //std::cout << "Number of events in pre opposite-flavor category = " << sum_pre_of_EventWeights << endl;
    //std::cout << "Number of events in same-flavor category = " << sum_sf_EventWeights  << endl;
    //std::cout << "Number of events in opposite-flavor category = " << sum_of_EventWeights << endl;
    /*
    std::cout << "Number of events in mu-mu category = " << sumW_mm_EventWeights << endl;
    std::cout << "Number of events in el-el category = " << sumW_ee_EventWeights << endl;
    std::cout << "Number of events in el-mu category = " << sumW_em_EventWeights << endl;
    std::cout << "N events 4mu = " << nEvents_4mu << endl;
    std::cout << "N events 3mu = " << nEvents_3mu << endl;
    std::cout << "N events 2mu = " << nEvents_2mu << endl;
    std::cout << "N events 1mu = " << nEvents_1mu << endl;
    std::cout << "N events 0mu = " << nEvents_4el << endl;
    std::cout << "N events Zee WeWe = " << nEvents_Zee_WeWe << endl;
    std::cout << "N events Zee WeWm = " << nEvents_Zee_WeWm << endl;
    std::cout << "N events Zee WmWm = " << nEvents_Zee_WmWm << endl;
    std::cout << "N events Zmm WeWe = " << nEvents_Zmm_WeWe << endl;
    std::cout << "N events Zmm WeWm = " << nEvents_Zmm_WeWm << endl;
    std::cout << "N events Zmm WmWm = " << nEvents_Zmm_WmWm << endl; 
    std::cout << "N events 2+ muons, 0 electrons = " << nEvents_2pmu_0el << endl;
    std::cout << "N events 1+ muons, 1+ electrons = " << nEvents_1pmu_1pel << endl;
    std::cout << "N events 0 muons, 2+ electrons = " << nEvents_0mu_2pel << endl;
    */
    f1->Write();
    f1->Close();
    return 0;
}
