import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
#            '/store/data/Run2018A/SingleMuon/MINIAOD/17Sep2018-v2/00000/0034E0CE-45EA-2249-A559-62746A83F4F1.root'
#'/store/data/Run2017F/SingleMuon/MINIAOD/31Mar2018-v1/00000/00245C8B-8C37-E811-B942-002590E2F9D4.root'
'/store/data/Run2016H/SingleMuon/MINIAOD/17Jul2018-v1/00000/02C005C4-998B-E811-A5CA-001E67DDD348.root'
                ),
#                                    skipEvents=cms.untracked.uint32(50)
                            )

process.demo = cms.EDAnalyzer('edanalyzer',

                            JetCollectionTag = cms.InputTag("slimmedJets"),
                            TriggerObjectCollectionTag = cms.InputTag("slimmedPatTrigger")
                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hltJetMetNtuple.root")
)

process.p = cms.Path(process.demo)
