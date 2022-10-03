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

    TFile* f1 = new TFile("output_checkGenPart.root", "RECREATE");
    //float pt = 0;
    //TTree tree_out("tree","");
    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();
    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);
    tqdm bar;

    int nEvents_ZHWW = 0;
    int nEvents_WHWW = 0;
    int nEvents_WHZZ = 0;
    int nEvents_ZHZZ = 0;

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

	    int V_abspdgid = abs(nt.GenPart_pdgId().at(2));
	    int HXX_abspdgid = 0;
	    int nWs = 0;
            int nZs = 0;

	    for (unsigned int igen = 0; igen < nt.GenPart_pdgId().size(); ++igen){
		 // Get the particles with higgs as mother that is not higgs
		 if ( std::abs(nt.GenPart_pdgId().at(igen)) == 25 ) continue;
		 if ( std::abs(nt.GenPart_pdgId().at(igen)) == 24 && std::abs(nt.GenPart_pdgId().at(nt.GenPart_genPartIdxMother().at(igen))) == 25 ){
			nWs++;
	         }
		 if ( std::abs(nt.GenPart_pdgId().at(igen)) == 23 && std::abs(nt.GenPart_pdgId().at(nt.GenPart_genPartIdxMother().at(igen))) == 25 ){
                        nZs++;
                 }
	    }

	    if ( nWs == 2 && V_abspdgid == 23 ) nEvents_ZHWW++;
	    if ( nWs == 2 && V_abspdgid == 24 ) nEvents_WHWW++;
	    if ( nZs == 2 && V_abspdgid == 23 ) nEvents_ZHZZ++;
	    if ( nZs == 2 && V_abspdgid == 24 ) nEvents_WHZZ++;


        } // Event loop

        delete file;


    } // File loop
    bar.finish();

    std::cout << "Number of ZHWW events = " << nEvents_ZHWW << endl;
    std::cout << "Number of WHWW events = " << nEvents_WHWW << endl;
    std::cout << "Number of ZHZZ events = " << nEvents_ZHZZ << endl;
    std::cout << "Number of WHZZ events = " << nEvents_WHZZ << endl;

    f1->Write();
    f1->Close();
    return 0;
}
