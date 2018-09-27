import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v16', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#,SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

#'/store/data/Run2017B/Charmonium/MINIAOD/PromptReco-v1/000/297/050/00000/183B4680-3356-E711-B33A-02163E014487.root',
'/store/data/Run2016H/Charmonium/MINIAOD/03Feb2017_ver3-v1/100000/285996B8-64EC-E611-BF01-008CFAFBEFE8.root',
 )
)

#process.load("myAnalyzers.JPsiKsPAT.miniAODmuonsRootupler_cfi")
#process.rootuple.dimuons = cms.InputTag('miniaodPATMuonsWithTrigger') 
process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_Dimuon20_Jpsi_v*',
                                                                        'HLT_Dimuon16_Jpsi_v*',
                                                                        'HLT_DoubleMu4_3_Jpsi_Displaced_v*',
                                                                        'HLT_DoubleMu4_JpsiTrk_Displaced_v*',
                                                                       # 'HLT_DoubleMu4_Jpsi_Displaced_v*'                                   
                                                                       ),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )
process.load("BUAA-CMS-group.JPsiUpsilonPhi.slimmedMuonsTriggerMatcher2016_cfi")

process.rootuple = cms.EDAnalyzer('jpsipipi',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          isMC = cms.bool(False),
                          )
process.rootuple.dimuons = cms.InputTag('slimmedMuonsWithTrigger')

process.TFileService = cms.Service("TFileService",
       #fileName = cms.string('Rootuple_Jpsi_2017-MiniAOD.root'),
       fileName = cms.string('Rootuple_JpsiPiPi.root'),
)

process.mySequence = cms.Sequence(
                                   process.triggerSelection *
                                   process.slimmedMuonsWithTriggerSequence *
                                   process.rootuple
                                   )
process.p = cms.Path(process.triggerSelection*process.slimmedMuonsWithTriggerSequence*process.rootuple)