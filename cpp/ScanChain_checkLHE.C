#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"
#include "TH2.h"

#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Base.h"
#include "../NanoCORE/tqdm.h"

// For.. gconf?
#include "../NanoCORE/SSSelections.cc"
#include "../NanoCORE/MetSelections.cc"

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

int ScanChain(TChain *ch) {

    TFile* f1 = new TFile("output_testLHE.root", "RECREATE");
    //float pt = 0;
    //TTree tree_out("tree","");
    //tree_out.Branch("pt", &pt);
    H1(nLeps,6,0,6);
    H1(nCandMu,6,0,6);
    H1(nCandEl,6,0,6);
    H1(nCandLeps,6,0,6);
    H1(fourLeps,2,0,2);
    H1(mu_pt,50,0,500);
    H1(el_pt,50,0,500);
    H1(mu_eta,50,-4.,4.);
    H1(el_eta,50,-4.,4.);
    H1(m_Z1,50,0,300);
    H1(m_Z2,50,0,300);
    //H2(mZ1,mZ2,50,0,300,50,0,300);
    H1(m_ee,50,0,300);
    H1(m_mm,50,0,300);
    TH2* h_Z1Z2 = new TH2D("h2","h2",50,0.,300.,50,0.,300.);

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();
    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);
    tqdm bar;

    //gconf.year = 2018;
    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile *file = TFile::Open( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        TString filename(currentFile->GetTitle());

        tree->SetCacheSize(128*1024*1024);
        tree->SetCacheLearnEntries(100);

        nt.Init(tree);

        for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

            nt.GetEntry(event);
            tree->LoadTree(event);


            nEventsTotal++;
            bar.progress(nEventsTotal, nEventsChain);

	    int nLHE_el = 0;
	    std::vector<int> lhe_el;
	    int nLHE_mu = 0;
	    std::vector<int> lhe_mu;

	    int nCand_mu = 0;
            int nCand_el = 0;

	    bool four_leps = false;

            for ( int i = 0; i < nt.nLHEPart(); i++ ){
		 if ( !(nt.LHEPart_status().at(i) == 1) ) continue;
		 if ( std::abs(nt.LHEPart_pdgId().at(i)) == 11 ){
		      lhe_el.push_back(i);
		      h_el_pt->Fill(nt.LHEPart_pt().at(i));
		      h_el_eta->Fill(nt.LHEPart_eta().at(i));
		      if ( nt.LHEPart_pt().at(i) > 10 && std::abs(nt.LHEPart_eta().at(i)) < 2.5 ){
			   nCand_el++;	
		      }
		 }
		 if ( std::abs(nt.LHEPart_pdgId().at(i)) == 13 ){
		      lhe_mu.push_back(i);
		      h_mu_pt->Fill(nt.LHEPart_pt().at(i));
		      h_mu_eta->Fill(nt.LHEPart_eta().at(i));
		      if ( nt.LHEPart_pt().at(i) > 10 && std::abs(nt.LHEPart_eta().at(i)) < 2.4 ){
                           nCand_mu++;
                      }
		 }
	    }

	    nLHE_el = lhe_el.size();
	    nLHE_mu = lhe_mu.size();

	    if (nCand_mu+nCand_el >=4){
		four_leps = true;
	    }

	    float m_Z1 = 0;
	    float m_Z2 = 0;
	    float m_ee = 0;
	    float m_mm = 0;

	    if ( nLHE_el > 3 ){
		 m_Z1 = (nt.LHEPart_p4().at(lhe_el[0])+nt.LHEPart_p4().at(lhe_el[1])).M();
		 m_Z2 = (nt.LHEPart_p4().at(lhe_el[2])+nt.LHEPart_p4().at(lhe_el[3])).M(); 
	    }

	    if ( nLHE_mu > 3 ){
                 m_Z1 = (nt.LHEPart_p4().at(lhe_mu[0])+nt.LHEPart_p4().at(lhe_mu[1])).M();
                 m_Z2 = (nt.LHEPart_p4().at(lhe_mu[2])+nt.LHEPart_p4().at(lhe_mu[3])).M();
            }

	    if ( nLHE_mu > 1 && nLHE_el > 1 ){
		 m_ee = (nt.LHEPart_p4().at(lhe_el[0])+nt.LHEPart_p4().at(lhe_el[1])).M();
		 m_mm = (nt.LHEPart_p4().at(lhe_mu[0])+nt.LHEPart_p4().at(lhe_mu[1])).M();	 
	    }

	    h_nLeps->Fill(nLHE_el+nLHE_mu);
	    h_nCandMu->Fill(nCand_mu);
	    h_nCandEl->Fill(nCand_el);
	    h_nCandLeps->Fill(nCand_mu+nCand_el);
	    h_fourLeps->Fill(four_leps);
	    if (m_Z1 > 0.0) h_m_Z1->Fill(m_Z1);
	    if (m_Z2 > 0.0) h_m_Z2->Fill(m_Z2);
	    if (m_ee > 0.0) h_m_ee->Fill(m_ee);
	    if (m_mm > 0.0) h_m_mm->Fill(m_mm);

	    if (m_Z1 > 0.0 && m_Z2 > 0.0) h_Z1Z2->Fill(m_Z1,m_Z2);


        } // Event loop

        delete file;


    } // File loop
    bar.finish();

    f1->Write();
    f1->Close();
    return 0;
}
