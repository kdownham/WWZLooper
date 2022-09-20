{

    gROOT->ProcessLine(".L ../NanoCORE/NANO_CORE.so");
    gROOT->ProcessLine(".L ScanChain_checkLHE.C+");
    TChain *ch = new TChain("Events");
    TString baseDir = "~/Triboson/VVVNanoLooper/mc/";

    ch->Add(baseDir+"RunIISummer20UL18NanoAODv9/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/7075899E-49EC-3B4F-BA70-877BC8E8C8CF.root");

    ScanChain(ch);

}

