##########################
#                        #
#  2017-01-16 MC         #
#                        #
##########################
dataset = {
'/BcToJPsiBcPt8Y2p5_MuNoCut_13TeV-bcvegpy2-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
}

nfiles = {
}

filesPerJob = {
}

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand

    def submit(config):
        res = crabCommand('submit', config = config)

    from CRABClient.UserUtilities import config
    config = config()

    name = '_data'
    config.General.workArea = 'crab_MC'
    config.General.transferLogs = False
    config.General.transferOutputs = True
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'test_MC.py'
    config.JobType.pyCfgParams = ['DataProcessing=mc2016']
    config.JobType.inputFiles   = ['data','rcdata.2016.v3']
  #  config.JobType.inputFiles   = ['data','rcdata.2016.v3']
  #  config.JobType.pyCfgParams = ['DataProcessing=legacy_rerecodata']
#    config.JobType.pyCfgParams = ['DataProcessing=rerecodata']
#    config.JobType.pyCfgParams = ['DataProcessing=promptdata']
    config.Data.inputDBS = 'global'
    #config.Data.splitting = 'Automatic'
    config.Data.splitting = 'EventAwareLumiBased'
    config.Data.unitsPerJob = 20000
    config.Data.publication = False
    config.Data.ignoreLocality = False
    config.Data.outLFNDirBase = '/store/user/hanwen/BcToJPsiBcPt8Y2p5'
 #   config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
#    config.Data.lumiMask = '/user/wenxing/FINAL_code_170407/CMSSW_8_0_26_patch1/src/UserCode/IIHETree/test/resubmit_DoubleEG.json'
#    config.Data.lumiMask = '/user/wenxing/FINAL_code_170407/CMSSW_8_0_26_patch1/src/UserCode/IIHETree/test/resubmit_SingleEG.json'
    config.Site.storageSite = 'T2_BE_IIHE'

    for sample in dataset:
        config.General.requestName = sample
        config.Data.inputDataset = dataset[sample]
#        config.Data.outputDatasetTag = sample
        submit(config)
