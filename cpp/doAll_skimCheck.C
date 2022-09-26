{

    gROOT->ProcessLine(".L ../NanoCORE/NANO_CORE.so");
    gROOT->ProcessLine(".L ScanChain_skimCheck.C+");
    TChain *ch = new TChain("Events");
    TString baseDir = "/ceph/cms/store/user/kdownham/skimOutput/WWZ_4L/";

    ch->Add(baseDir+"GluGluZH_HToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_4L/*");

    ScanChain(ch);

}

