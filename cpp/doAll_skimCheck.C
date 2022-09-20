{

    gROOT->ProcessLine(".L ../NanoCORE/NANO_CORE.so");
    gROOT->ProcessLine(".L ScanChain_skimCheck.C+");
    TChain *ch = new TChain("Events");
    TString baseDir = "~/Triboson/NanoSkimmer/Run2018/";

    ch->Add(baseDir+"7075899E-49EC-3B4F-BA70-877BC8E8C8CF_Skim.root");

    ScanChain(ch);

}

