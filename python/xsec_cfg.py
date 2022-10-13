# run like
# cmsRun xsec_cfg.py datasetName=/WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM maxEvents=100000

from sys import version_info as sys_version_info

def check_module_exists(module_name):
   #print("Sys version: ",sys_version_info)
   impexcpt = None
   if sys_version_info < (3, 0):
      # python 2
      from pkgutil import find_loader as find_module_loader
      impexcpt = ImportError
   elif sys_version_info <= (3, 3):
      # python 3.0 to 3.3
      from importlib import find_loader as find_module_loader
   elif sys_version_info >= (3, 4):
      # python 3.4 and above
      from importlib import util
      find_module_loader = util.find_spec
   res = None
   if impexcpt is not None:
      try:
         res = find_module_loader(module_name)
      except impexcpt:
         res = None
   else:
      res = find_module_loader(module_name)
   return res is not None


#if check_module_exists("IvyFramework.IvyDataTools.cmseostools"):
   #from IvyFramework.IvyDataTools.cmseostools import listFiles
   #from IvyFramework.IvyDataTools.cmseostools import findParent
if check_module_exists("cmseostools"):
   from cmseostools import listFiles
   from cmseostools import findParent
else:
   from cmseostools import listFiles
   from cmseostools import findParent


def get_max_files(DAS_name, max_files, dbase) :
   result = []
   file_names = listFiles(DAS_name, dbase)
   nfiles = len(file_names)
   if max_files>0:
      nfiles=min(len(file_names), max_files)
   for i in range(nfiles):
      if file_names[i]:
         result.append(file_names[i])
   return result


import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

options.register('datasetName',
		'',
		VarParsing.multiplicity.singleton,
		VarParsing.varType.string,
		"DAS-style name of a primary dataset, e.g. /ZZTo4L_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_v2_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM")
options.register('localdir', False, mytype=VarParsing.varType.bool)
options.register('useFileList', False, mytype=VarParsing.varType.bool)
options.register('maxfiles', 50, mytype=VarParsing.varType.int)
options.register('disableDuplicateCheck', False, mytype=VarParsing.varType.bool)

options.parseArguments()

# Always operate over MINIAOD
dsetname = str(options.datasetName)
if "NANOAOD" in dsetname:
   dsetname = findParent(dsetname)

print("Running over data set {}".format(dsetname))

input_is_already_flist = options.useFileList

dbase="dbs"
if input_is_already_flist:
   dbase="csflist"
elif options.localdir:
   dbase="local"

process = cms.Process('XSec')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents),
)

process.load('FWCore.MessageService.MessageLogger_cfi')
# How often do you want to see the "Begin processing the .."
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source ("PoolSource",
   fileNames = cms.untracked.vstring(*get_max_files(dsetname, options.maxfiles, dbase)),
)
if options.disableDuplicateCheck:
   process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.xsec = cms.EDAnalyzer("GenXSecAnalyzer")

process.ana = cms.Path(process.xsec)
process.schedule = cms.Schedule(process.ana)
