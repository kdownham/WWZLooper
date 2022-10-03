{

    gROOT->ProcessLine(".L ../NanoCORE/NANO_CORE.so");
    gROOT->ProcessLine(".L ScanChain_skimCheck.C+");
    TChain *ch = new TChain("Events");
    TChain *ch2 = new TChain("Events");
    TChain *ch3 = new TChain("Events");
    TChain *ch4 = new TChain("Events");
    TString baseDir = "/ceph/cms/store/user/kdownham/skimOutput/WWZ_4L/";

    ch->Add(baseDir+"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1_NANOAODSIM_WWZ_4L/*");
    ch->Add(baseDir+"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v2_NANOAODSIM_WWZ_4L/*");

    ScanChain(ch);

    ch2->Add(baseDir+"GluGluZH_HToWWTo2L2Nu_M125_13TeV_powheg_pythia8_TuneCP5_PSweights_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_4L/*");

    ScanChain(ch2);

    ch3->Add(baseDir+"HZJ_HToWWTo2L2Nu_ZTo2L_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_4L/*");

    ScanChain(ch3);

    ch4->Add(baseDir+"VHToNonbb_M125_TuneCP5_13TeV-amcatnloFXFX_madspin_pythia8_RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2_NANOAODSIM_WWZ_4L/*");

    ScanChain(ch4);

}

