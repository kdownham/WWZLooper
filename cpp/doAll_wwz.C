#include "Riostream.h"

{

    gROOT->ProcessLine(".L ../NanoCORE/NANO_CORE.so");
    gROOT->ProcessLine(".L ScanChain_wwz.C+");
    TChain *ch = new TChain("Events");
    TString baseDir = "/ceph/cms/store/user/skimOutput/wwz_4lep/";

    ch->Add(baseDir+"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_wwz_4lep/*");

    ScanChain(ch,"WWZ","2018",get_xsec("WWZ"),getSumOfGenEventSumw(ch),true);

}
// Called in ScanChain function ----- See how Z prime does this
double getSumOfGenEventSumw(TChain *chaux){
  double genEventSumw, sumOfGenEventSumw=0.0;
  chaux->SetBranchAddress("genEventSumw",&genEventSumw);
  for (unsigned int run = 0; run < chaux->GetEntriesFast(); run++){
      chaux->GetEntry(run);
      sumOfGenEventSumw += genEventSumw;
  }

  return sumOfGenEventSumw;
}

double get_xsec(TString sample){
     double xsec;
     std::ifstream f("data/xsecs.txt");
     std::string line;
     while (std::getline(f, line)){
	TString col1;
	double col2;
	std::istringstream ss(line);
	ss >> col1 >> col2;
	if(sample == col1){
	   xsec = 1000.*col2;  // Convert from pb to fb
	}
     }
      
     return xsec; 
}
