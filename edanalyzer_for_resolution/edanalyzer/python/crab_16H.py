from WMCore.Configuration import Configuration

store_dir = 'SingleMuon_16H'
sample_name = 'SingleMuon_16H'

MIN_DSET = '/SingleMuon/Run2016H-17Jul2018-v1/MINIAOD'
#RAW_DSET = '/SingleMuon/Run2016B-v2/RAW'

config = Configuration()

config.section_('General')
config.General.requestName = 'jmeTriggerNTuple_'+sample_name
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName = 'for_making_res_plots.py'
config.JobType.inputFiles = []
config.JobType.pyCfgParams = ['output='+sample_name+'.root']
config.JobType.maxJobRuntimeMin = 1000
config.JobType.maxMemoryMB = 2000
#config.JobType.numCores = 4

config.section_('Data')
config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.splitting = 'FileBased'
config.Data.inputDataset = MIN_DSET
#config.Data.secondaryInputDataset = RAW_DSET
config.Data.outLFNDirBase = '/store/user/dvats/for_res/'+store_dir+'/'+sample_name
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.lumiMask = '/eos/home-d/dvats/jet_eff_plots_area/CMSSW_10_2_22/src/JMETriggerAnalysis/NTuplizers/test/crab/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'

config.section_('Site')
config.Site.storageSite = 'T2_IN_TIFR'
if config.Data.ignoreLocality:
   config.Site.whitelist = ['T2_CH_CERN', 'T2_DE_*']

config.section_('User')
