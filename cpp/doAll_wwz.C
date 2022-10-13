R__LOAD_LIBRARY(../NanoCORE/NANO_CORE.so)
#include "ScanChain_wwz.C+"
#include "Riostream.h"

double getSumOfGenEventSumw(TChain *chaux);
double get_xsec(TString sample);

void doAll_wwz(){

    //gROOT->ProcessLine(".L ../NanoCORE/NANO_CORE.so");
    //gROOT->ProcessLine(".L ScanChain_wwz.C+");
    TChain *ch = new TChain("Events");
    TChain *ch_aux = new TChain("Runs");
    TChain *ch2 = new TChain("Events");
    TChain *ch2_aux = new TChain("Runs");
    TChain *ch3 = new TChain("Events");
    TChain *ch3_aux = new TChain("Runs");
    TChain *ch4 = new TChain("Events");
    TChain *ch4_aux = new TChain("Runs");
    TString baseDir = "/ceph/cms/store/user/kdownham/skimOutput/WWZ_4L/";

    TString WWZ_dir_18 = "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_WWZ_4L";
    TString WWZ_ext1_dir_18 = "WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2_NANOAODSIM_WWZ_4L";
    TString GGZH_dir_18 = "GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_4L";
    TString HZJ_dir_18 = "HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_4L";
    TString WWZJets_dir_18 = "WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_4L"; 

    ch->Add(baseDir+WWZ_dir_18+"/*");
    ch->Add(baseDir+WWZ_ext1_dir_18+"/*");
    ch_aux->Add(baseDir+WWZ_dir_18+"/*");
    ch_aux->Add(baseDir+WWZ_ext1_dir_18+"/*");
    //ScanChain(ch,"WWZ","2018",get_xsec("WWZ"),getSumOfGenEventSumw(ch_aux),true);
    
    ch2->Add(baseDir+GGZH_dir_18+"/*");
    ch2_aux->Add(baseDir+GGZH_dir_18+"/*");
    //ScanChain(ch2,"ggZH","2018",get_xsec("ggZH"),getSumOfGenEventSumw(ch2_aux),true);

    ch3->Add(baseDir+HZJ_dir_18+"/*");
    ch3_aux->Add(baseDir+HZJ_dir_18+"/*");
    //ScanChain(ch3,"HZJ","2018",get_xsec("HZJ"),getSumOfGenEventSumw(ch3_aux),true);

    ch4->Add(baseDir+WWZJets_dir_18+"/output_7.root");
    ch4_aux->Add(baseDir+WWZJets_dir_18+"/output_7.root");
    ScanChain(ch4,"WWZJets","2018",get_xsec("WWZJets"),getSumOfGenEventSumw(ch4_aux),true);
    
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
