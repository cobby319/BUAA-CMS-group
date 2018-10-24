##########################
#                        #
#  2017-01-16 MC         #
#                        #
##########################
dataset = {
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

    name = '_mc'
    config.General.workArea = 'crab_projects_TTsys_20180426'
    config.General.transferLogs = False
    config.General.transferOutputs = True
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'IIHE.py'
#    config.JobType.psetName = 'IIHE_PL.py'
    config.JobType.pyCfgParams = ['DataProcessing=mc2016']
    config.JobType.inputFiles   = ['data','rcdata.2016.v3']
    config.Data.inputDBS = 'global'
    #config.Data.splitting = 'FileBased'
    #config.Data.unitsPerJob = 1
    config.Data.splitting = 'EventAwareLumiBased'
    config.Data.unitsPerJob = 20000
    config.Data.publication = False
    config.Data.ignoreLocality = False
    config.Data.outLFNDirBase = '/store/user/wenxing/TTbar_sys_ext_sample_2016'
    config.Site.storageSite = 'T2_BE_IIHE'

    for sample in dataset:
        config.General.requestName = sample
        config.Data.inputDataset = dataset[sample]
#        config.Data.outputDatasetTag = sample
        submit(config)

