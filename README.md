
## This is a WIP of course

### Environment & setup
```bash
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_9_2_8
cmsenv
cd -
cd NanoCORE
make -j8 # takes about 90 seconds
cd -
```

### Python dependencies
```bash
pip install uproot --user
pip install backports.lzma --user
```


### Test file
```bash
# from /TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM
xrdcp root://xrootd.unl.edu//store/mc/RunIIFall17NanoAOD/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/508B6DBB-9742-E811-BA65-A4BF0112BCB0.root .
```

