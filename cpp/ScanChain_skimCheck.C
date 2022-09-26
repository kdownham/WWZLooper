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

    TFile* f1 = new TFile("skimCheck_output/output_skimTest.root", "RECREATE");
    //float pt = 0;
    //TTree tree_out("tree","");
    //tree_out.Branch("pt", &pt);
    //H1(isSkim,2,0,1);

    int nEventsTotal = 0;
    int nEventsPassing = 0;
    int nEventsFourMu = 0;
    int nEventsTwoMuTwoEl = 0;
    int nEventsFourEl = 0;
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
	    //bool isSkimEvent = false;
            int nCand_leps = 0;
            int nCand_mu = 0;
	    int nCand_ele = 0;

	    //Loop over muons to find muon candidates
	    for ( int i = 0; i < nt.nMuon(); i++ ){
		 if ( !(nt.Muon_pt().at(i) > 10.) ) continue;
		 if ( !(std::abs(nt.Muon_eta().at(i)) < 2.4) ) continue;
		 if ( !nt.Muon_looseId().at(i) ) continue;
		 if ( !(nt.Muon_pfIsoId().at(i) >= 1) ) continue;
		 nCand_mu++;
	    }    

            // Do the same with electrons
            for ( int k = 0; k < nt.nElectron(); k++ ){
		 if ( !(nt.Electron_pt().at(k) > 10.) ) continue;
		 if ( !(std::abs(nt.Electron_eta().at(k)) < 2.5) ) continue;
		 if ( !nt.Electron_mvaFall17V2noIso_WPL().at(k) ) continue;
		 if ( !(nt.Electron_pfRelIso03_all().at(k) < 0.4) ) continue;
		 nCand_ele++;
	    }

	    nCand_leps = nCand_mu + nCand_ele;
	    if ( nCand_leps < 4 ) continue;
	    nEventsPassing++;

	    if ( nCand_mu > 3 ) nEventsFourMu++;
	    if ( nCand_ele > 3 ) nEventsFourEl++;
	    if ( nCand_mu > 1 && nCand_ele > 1 ) nEventsTwoMuTwoEl++;

	    //h_isSkim->Fill(isSkimEvent);

        } // Event loop

        delete file;


    } // File loop

    bar.finish();

    std::cout << nEventsPassing << " events passed skim requirements, out of " << nEventsTotal << " total events" << endl;
    std::cout << "Number of events with at least 4 candidate muons = " << nEventsFourMu << endl;
    std::cout << "Number of events with at least 4 candidate electrons = " << nEventsFourEl << endl;
    std::cout << "Number of events with at least 2 candidate muons and 2 candidate electrons = " << nEventsTwoMuTwoEl << endl;

    f1->Write();
    f1->Close();
    return 0;
}
