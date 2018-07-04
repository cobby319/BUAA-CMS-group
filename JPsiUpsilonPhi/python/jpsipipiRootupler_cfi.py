import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('jpsipipi',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          isMC = cms.bool(False),
                          )
