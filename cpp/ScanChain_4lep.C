#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"
#include "TLorentzVector.h"
#include <math.h>
//#include "Math/GenVector/Boost.h"
#include "Math/Boost.h"


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
    double m_Wcands;

    double sum_of_EventWeights;
    double sum_sf_EventWeights;
    double sumOfEventWeights;
    double sumSQEventWeights;
    double sumWCut1;
    double sumWCut2;
    double sumWCut3;
    double sumWCut4;
    double sumWCut4_1;
    double sumWCut4_2;
    double sumWCut4_3;
    double sumWCut4_4;
    double sumWCut4_5;
    double sumWCut4_6;
    double sumWCut4_7;
    double sumWCut4_8;
    double sumWCut4_9;
    double sumWCut4_10;
    double sumWCut4_11;
    double sumWCut4_12;
    double sumWCut5;
    double sumWCut6;
    double sumWCut7;

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
    tree_out.Branch("m_Wcands",&m_Wcands);

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

	    bool triggerReqs = false;

	    if ( year == "2016" ) triggerReqs = (nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL() || nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ());
	    if ( year == "2017" ) triggerReqs = (nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() || nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() | nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ());
	    if ( year == "2018" ) triggerReqs = (nt.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() || nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL() || nt.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() || nt.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ());

	    if (!triggerReqs) continue;

            nEventsTotal++;
            bar.progress(nEventsTotal, nEventsChain);

	    float weight = genWeight();
	    evt_weight = weight*factor;

	    // Cuts, filling hists, branches, etc
	    
////////////// Pre-selection////////////////////////////////////////////////////////////////////////
            // Muons
            if (debug) std::cout << "1" << endl;
            std::vector<pair<int,int>> cand_leps;
	    for ( int i = 0; i < nt.nMuon(); i++ ){
		  if (!nt.Muon_looseId().at(i)) continue;
		  if (!(std::abs(nt.Muon_eta().at(i)) < 2.4)) continue;
		  if (!(nt.Muon_pt().at(i) > 10.)) continue;
		  if (!(nt.Muon_pfIsoId().at(i) >= 1)) continue;
		  cand_leps.emplace_back(i,nt.Muon_pdgId().at(i));
	    }

	    // Electrons
	    if (debug) std::cout << "2" << endl;
	    for ( int j = 0; j < nt.nElectron(); j++ ){
		  if (!nt.Electron_mvaFall17V2noIso_WPL().at(j)) continue;
		  if (!(std::abs(nt.Electron_eta().at(j)) < 2.5 )) continue;
		  if (!(nt.Electron_pt().at(j) > 10.)) continue;
		  if (!(nt.Electron_pfRelIso03_all().at(j) < 0.4)) continue;
		  cand_leps.emplace_back(j,nt.Electron_pdgId().at(j));
            }

	    // Require exactly 4 candidate leptons
	    if (debug) std::cout << "3" << endl;
	    sumWCut1 += evt_weight;
	    if ( cand_leps.size() != 4 ) continue;
	    sumWCut2 += evt_weight;

///////////////Get W and Z candidates///////////////////////////////////////////////////////////////////

	    std::vector<pair<int,int>> Z_cands;
	    //int Z_cand_pdgId;
	    std::vector<pair<int,int>> W_cands;
	    double Z_cand_mass = 99999999999999.;

	    for ( int lep = 0; lep < cand_leps.size(); lep++ ){
		  LorentzVector p4_lep1;
		  int lep1_id = cand_leps[lep].second;
		  if (debug) std::cout << "loop 1" << endl;
		  if ( std::abs(lep1_id) == 11 ){
		       p4_lep1 = nt.Electron_p4().at(cand_leps[lep].first);
		  }
		  if (debug) std::cout << "loop 2" << endl;
		  if ( std::abs(lep1_id) == 13 ){
		       p4_lep1 = nt.Muon_p4().at(cand_leps[lep].first);
		  }
		  if (debug) std::cout << "loop 3" << endl;
		  for ( int lep2 = lep+1; lep2 < cand_leps.size(); lep2++ ){	
		       if (debug) std::cout << "loop 4" << endl;
		       LorentzVector p4_lep2;
		       int lep2_id = cand_leps[lep2].second;
		       if ( std::abs(lep2_id) == 11 ){ 
			    if ( std::abs(lep1_id) == 13 ) continue;
			    p4_lep2 = nt.Electron_p4().at(cand_leps[lep2].first);
		       }
		       if (debug) std::cout << "loop 5" << endl;
		       if ( std::abs(lep2_id) == 13 ){
			    if ( std::abs(lep1_id) == 11 ) continue;
			    p4_lep2 = nt.Muon_p4().at(cand_leps[lep2].first);
		       }

		       if ( lep1_id + lep2_id != 0 ) continue;
			
		       double m_Zcands = (p4_lep1+p4_lep2).M();
		       if ( std::abs(m_Zcands - 91.2) > 10. ) continue;

		       if ( std::abs(m_Zcands - 91.2) < std::abs(Z_cand_mass - 91.2) ){
			    Z_cands.clear();
			    Z_cands.emplace_back(lep,lep1_id);
			    Z_cands.emplace_back(lep2,lep2_id);
			    Z_cand_mass = m_Zcands;
			    if (debug) std::cout << "Z cand mass = " << Z_cand_mass << endl;
		       }

		  }
	    }
	    if (debug) std::cout << "4" << endl;

	    if ( Z_cands.size() != 2 ) continue;           // veto event if there is no Z candidate
	    sumWCut3 += evt_weight;	    

	    // Get W candidates
	    for ( int w = 0; w < cand_leps.size(); w++ ){
		  int w_idx = cand_leps[w].first;
		  int w_id = cand_leps[w].second;
		  if ( w_idx == Z_cands[0].first && w_id == Z_cands[0].second ) continue;
		  if ( w_idx == Z_cands[1].first && w_id == Z_cands[1].second ) continue;
		  W_cands.emplace_back(w_idx,w_id);	  
	    }

	    if (debug) std::cout << "5" << endl;

	    if ( W_cands.size() != 2 ) continue;	     // veto event if there are no W candidates
	    sumWCut4 += evt_weight;
///////////////Cuts for W and Z candidates//////////////////////////////////////////////////////////////////

	    double ip_3d_1;
	    double ip_3d_2;
	    double sip_3d_1;
	    double sip_3d_2;
	
	    int nWmu = 0;
	    int nWel = 0;

	    bool opp_flav = false;
	    bool same_flav = false;

	    LorentzVector p4_l_Z1;
	    LorentzVector p4_l_Z2;
	    LorentzVector p4_l_W1;
	    LorentzVector p4_l_W2;

	    // Z cands are electrons
	    if ( std::abs(Z_cands[0].second) == 11 ){
	 	 ip_3d_1 = nt.Electron_ip3d()[Z_cands[0].first];
		 ip_3d_2 = nt.Electron_ip3d()[Z_cands[1].first];
		 sip_3d_1 = nt.Electron_sip3d()[Z_cands[0].first];                    
                 sip_3d_2 = nt.Electron_sip3d()[Z_cands[1].first];
		 if ( !(std::abs(ip_3d_1/sip_3d_1) < 4) || !(std::abs(ip_3d_2/sip_3d_2) < 4) ) continue;
		 //sumWCut4_1 += evt_weight;	 
		 if ( !nt.Electron_mvaFall17V2noIso_WPL()[Z_cands[0].first] || !nt.Electron_mvaFall17V2noIso_WPL()[Z_cands[1].first] ) continue;
		 //sumWCut4_2 += evt_weight;
		 if ( !(nt.Electron_pt()[Z_cands[0].first] > 25.) || !(nt.Electron_pt()[Z_cands[1].first] > 15.) ) continue;     
		 pt_l_Z1 = nt.Electron_pt()[Z_cands[0].first];
	         pt_l_Z2 = nt.Electron_pt()[Z_cands[1].first];
		 p4_l_Z1 = nt.Electron_p4()[Z_cands[0].first];
		 p4_l_Z2 = nt.Electron_p4()[Z_cands[1].first];
	    } 
	    //sumWCut4_3 += evt_weight;
	    // Z cands are muons
	    if ( std::abs(Z_cands[0].second) == 13 ){
		 ip_3d_1 = nt.Muon_ip3d()[Z_cands[0].first];                    
                 ip_3d_2 = nt.Muon_ip3d()[Z_cands[1].first];
                 sip_3d_1 = nt.Muon_sip3d()[Z_cands[0].first];
                 sip_3d_2 = nt.Muon_sip3d()[Z_cands[1].first];
                 if ( !(std::abs(ip_3d_1/sip_3d_1) < 4) || !(std::abs(ip_3d_2/sip_3d_2) < 4) ) continue;
		 //sumWCut4_4 += evt_weight;
		 if ( !nt.Muon_mediumId()[Z_cands[0].first] || !nt.Muon_mediumId()[Z_cands[1].first] ) continue;
		 //sumWCut4_5 += evt_weight;
		 if ( !(nt.Muon_pfRelIso04_all()[Z_cands[0].first] < 0.25) || !(nt.Muon_pfRelIso04_all()[Z_cands[1].first] < 0.25) ) continue;
		 //sumWCut4_6 += evt_weight;
		 if ( !(nt.Muon_pt()[Z_cands[0].first] > 25.) || !(nt.Muon_pt()[Z_cands[1].first] > 15.) ) continue;
		 pt_l_Z1 = nt.Muon_pt()[Z_cands[0].first];
                 pt_l_Z2 = nt.Muon_pt()[Z_cands[1].first];
                 p4_l_Z1 = nt.Muon_p4()[Z_cands[0].first];
                 p4_l_Z2 = nt.Muon_p4()[Z_cands[1].first];
	    }
	    //sumWCut4_7 += evt_weight;
            // Leading W is electron
	    if ( std::abs(W_cands[0].second) == 11 ){
		 if ( !nt.Electron_mvaFall17V2Iso_WP90()[W_cands[0].first] ) continue;
		 if ( !(std::abs(nt.Electron_ip3d()[W_cands[0].first]/nt.Electron_sip3d()[W_cands[0].first]) < 4) ) continue; 
		 if ( !(nt.Electron_pt()[W_cands[0].first] > 25.) ) continue;
		 nWel++;
		 pt_l_W1 = nt.Electron_pt()[W_cands[0].first];
		 p4_l_W1 = nt.Electron_p4()[W_cands[0].first];
	    }
	    //sumWCut4_8 += evt_weight;
	    // Subleading W is electron
	    if ( std::abs(W_cands[1].second) == 11 ){
                 if ( !nt.Electron_mvaFall17V2Iso_WP90()[W_cands[1].first] ) continue;
                 if ( !(std::abs(nt.Electron_ip3d()[W_cands[1].first]/nt.Electron_sip3d()[W_cands[1].first]) < 4) ) continue;
                 if ( !(nt.Electron_pt()[W_cands[1].first] > 15.) ) continue;
                 nWel++;
		 pt_l_W2 = nt.Electron_pt()[W_cands[1].first];
                 p4_l_W2 = nt.Electron_p4()[W_cands[1].first];
            }
	    //sumWCut4_9 += evt_weight;
	    // Leading W is muon
	    if ( std::abs(W_cands[0].second) == 13 ){
                 if ( !nt.Muon_mediumId()[W_cands[0].first] ) continue;
		 if ( !(nt.Muon_pfRelIso04_all()[W_cands[0].first] < 0.15) ) continue;
                 if ( !(std::abs(nt.Muon_ip3d()[W_cands[0].first]/nt.Muon_sip3d()[W_cands[0].first]) < 4) ) continue;
                 if ( !(nt.Muon_pt()[W_cands[0].first] > 25.) ) continue;
                 nWmu++;
		 pt_l_W1 = nt.Muon_pt()[W_cands[0].first];
                 p4_l_W1 = nt.Muon_p4()[W_cands[0].first];
            }
	    //sumWCut4_10 += evt_weight;
	    // Subleading W is muon
	    if ( std::abs(W_cands[1].second) == 13 ){
                 if ( !nt.Muon_mediumId()[W_cands[1].first] ) continue;
                 if ( !(nt.Muon_pfRelIso04_all()[W_cands[1].first] < 0.15) ) continue;
                 if ( !(std::abs(nt.Muon_ip3d()[W_cands[1].first]/nt.Muon_sip3d()[W_cands[1].first]) < 4) ) continue;
                 if ( !(nt.Muon_pt()[W_cands[1].first] > 15.) ) continue;
                 nWmu++;
		 pt_l_W2 = nt.Muon_pt()[W_cands[1].first];
                 p4_l_W2 = nt.Muon_p4()[W_cands[1].first];
            }

	    sumWCut5 += evt_weight;

            if ( nWmu == 2 || nWel == 2 ) same_flav = true;
	    if ( nWmu == 1 && nWel == 1 ) opp_flav = true;

///////////////Low mass resonance veto//////////////////////////////////////////////////////////////////////

	    double m_Z1Z2 = (p4_l_Z1+p4_l_Z2).M();
	    double m_W1W2 = (p4_l_W1+p4_l_W2).M();
	    double m_Z1W2 = (p4_l_Z1+p4_l_W2).M();
            double m_W1Z2 = (p4_l_W1+p4_l_Z2).M(); 

	    if ( m_Z1Z2 < 12. || m_W1W2 < 12. || m_Z1W2 < 12. || m_W1Z2 < 12. ) continue;
             
	    sumWCut6 += evt_weight;          
///////////////b-tagged-jet veto////////////////////////////////////////////////////////////////////////////
	    
	    int nBjets = 0;
	    for (int jet_idx = 0; jet_idx < nt.nJet(); jet_idx++){
		 if ( nt.Jet_btagDeepB().at(jet_idx) > 0.1208 && nt.Jet_pt().at(jet_idx) > 20. ){
		      nBjets++;
		 }
	    }	  

	    if (nBjets > 0) continue;
	    
	    sumWCut7 += evt_weight;
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

	    asymm_mt2_lester_bisect::disableCopyrightMessage();
	    MT2 = asymm_mt2_lester_bisect::get_mT2( mVisA, boosted_p4_l_W1.Px(), boosted_p4_l_W1.Py(), mVisB, boosted_p4_l_W2.Px(), boosted_p4_l_W2.Py(), MET_p4.Px(), MET_p4.Py(), chiA, chiB, desiredPrecisionOnMt2);

	    // Construct pT_4L
	    double lep_px_sum = p4_l_Z1.Px() + p4_l_Z2.Px() + p4_l_W1.Px() + p4_l_W2.Px();
	    double lep_py_sum = p4_l_Z1.Py() + p4_l_Z2.Py() + p4_l_W1.Py() + p4_l_W2.Py();

	    double pt_4l = std::sqrt(std::pow(lep_px_sum,2.0)+std::pow(lep_py_sum,2.0));

	    bool sf_pass = false;
	    bool of_pass = false;

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

	    if ( opp_flav ){
		 if ( m_Wcands < 100. && MT2 < 25.) continue;
		 
		 if ( m_Wcands < 100. && MT2 > 25.) of_pass = true;
		      
                 if ( m_Wcands >= 100. ){
		      of_pass = true;
                 }
	    }

	
            if ( sf_pass ) sum_sf_EventWeights += evt_weight;
	    if ( of_pass ) sum_of_EventWeights += evt_weight;

	    
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
    std::cout << "Number of events in same-flavor category = " << sum_sf_EventWeights  << endl;
    std::cout << "Number of events in opposite-flavor category = " << sum_of_EventWeights << endl;
    std::cout << "Number of events passing cut 1 = " << sumWCut1 << endl;
    std::cout << "Number of events passing cut 2 = " << sumWCut2 << endl;
    std::cout << "Number of events passing cut 3 = " << sumWCut3 << endl;
    std::cout << "Number of events passing cut 4 = " << sumWCut4 << endl;
    /*
    std::cout << "Number of events passing cut 4.1 = " << sumWCut4_1 << endl;
    std::cout << "Number of events passing cut 4.2 = " << sumWCut4_2 << endl;
    std::cout << "Number of events passing cut 4.3 = " << sumWCut4_3 << endl;
    std::cout << "Number of events passing cut 4.4 = " << sumWCut4_4 << endl;
    std::cout << "Number of events passing cut 4.5 = " << sumWCut4_5 << endl;
    std::cout << "Number of events passing cut 4.6 = " << sumWCut4_6 << endl;
    std::cout << "Number of events passing cut 4.7 = " << sumWCut4_7 << endl;
    std::cout << "Number of events passing cut 4.8 = " << sumWCut4_8 << endl;
    std::cout << "Number of events passing cut 4.9 = " << sumWCut4_9 << endl;
    std::cout << "Number of events passing cut 4.10 = " << sumWCut4_10 << endl;
    */
    std::cout << "Number of events passing cut 5 = " << sumWCut5 << endl;
    std::cout << "Number of events passing cut 6 = " << sumWCut6 << endl;
    std::cout << "Number of events passing cut 7 = " << sumWCut7 << endl;
    f1->Write();
    f1->Close();
    return 0;
}
