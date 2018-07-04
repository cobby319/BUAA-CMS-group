import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('jpsiupsgamma',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          isMC = cms.bool(False),
                          )
