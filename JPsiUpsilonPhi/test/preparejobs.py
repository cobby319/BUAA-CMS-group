import sys
import re
import os

def printPath(path):  
    dirname = os.getcwd()
    if not os.path.exists("OUTPUTS/"+type+"/JOBS/catalog_"+type+".txt"): 
        os.mknod("OUTPUTS/"+type+"/JOBS/catalog_"+type+".txt")
    file=open("OUTPUTS/"+type+"/JOBS/catalog_"+type+".txt",'a')
    for dirs in os.listdir(path):
        for fs in os.listdir(path+'/'+dirs):
            if (os.path.isfile(path+'/'+dirs+'/'+fs)):
                file.write(path[21:]+'/'+dirs+'/'+fs)
                file.write("\n")
    file.close()


def prepare_job_script(theCatalog,jobID):
    global cfgDirectory  
    
    if os.path.exists(cfgDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'_cfg.py'): 
        os.remove(cfgDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'_cfg.py')
    os.mknod(cfgDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'_cfg.py') 

    if os.path.exists(jobsDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'.sh'): 
        os.remove(jobsDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'.sh')
    os.mknod(jobsDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'.sh') 

    configFile = open(cfgDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'_cfg.py','w')
    configLines  = ''
    configLines +="import FWCore.ParameterSet.Config as cms\n"
    configLines +="process = cms.Process('Rootuple')\n"
    configLines +="process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')\n"
    configLines +="process.load('Configuration.StandardSequences.Services_cff')\n"
    configLines +="process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')\n"
    configLines +="process.load('FWCore.MessageService.MessageLogger_cfi')\n"
    configLines +="process.load('Configuration.EventContent.EventContent_cff')\n"
    configLines +="process.load('Configuration.StandardSequences.GeometryRecoDB_cff')\n"
    configLines +="process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')\n"
    configLines +="process.load('Configuration.StandardSequences.EndOfProcess_cff')\n"
    configLines +="process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')\n"
    configLines +="from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag\n"
    if('Run2016H' in theCatalog):
        configLines +="process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v16', '')\n"
    else:
        configLines +="process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')\n"
    configLines +="process.MessageLogger.cerr.FwkReport.reportEvery = 200\n"
    configLines +="process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )\n"
    configLines +="process.options.allowUnscheduled = cms.untracked.bool(True)\n"
    configLines +="process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1),SkipEvent = cms.untracked.vstring('ProductNotFound'))\n"
    configLines +="process.source = cms.Source('PoolSource',\n"
    configLines +="    fileNames = cms.untracked.vstring(\n"
    configLines +=("'"+theCatalog+"'"+"\n")
    configLines +=",)\n"
    configLines +=")\n"
    configLines +="process.rootuple = cms.EDAnalyzer('jpsipipi',\n"
    configLines +='                          dimuons = cms.InputTag("slimmedMuons"),\n'
    configLines +='                          Trak = cms.InputTag("packedPFCandidates"),\n'
    configLines +='                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),\n'
    configLines +="                          isMC = cms.bool(False),\n"
    configLines +="                          )\n"
    configLines +='process.TFileService = cms.Service("TFileService",\n'
    configLines +="       fileName = cms.string('Rootuple_"+theCatalog[11:19]+"_"+str(jobID)+".root'),\n"
    configLines +=")\n"
    configLines +="process.p = cms.Path(process.rootuple)\n"
    configFile.write(configLines)
    configFile.close()



    scriptFile = open(jobsDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'.sh','w')
    scriptLines = ''
    scriptLines += 'source $VO_CMS_SW_DIR/cmsset_default.sh\n'
    scriptLines += 'export SCRAM_ARCH=slc6_amd64_gcc530\n'
#    scriptLines += 'export BUILD_ARCH=slc6_amd64_gcc530\n'
#    scriptLines += 'export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n'
#    scriptLines += 'export XRD_NETWORKSTACK=IPv4\n'
    scriptLines += ('export INITDIR='+os.getcwd()+'\n')
    scriptLines += ('cd $INITDIR\n')
    scriptLines += 'eval `scramv1 runtime -sh`\n'
    scriptLines += 'cd -\n'
#    scriptLines += 'ulimit -c 0;\n'
    scriptLines += 'if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi;\n'
#    scriptLines += 'cp '+os.getcwd()+'/makentupleproduce .;\n'
    scriptLines += 'hostname ;\n'
    iteFileInJob=0
#    for aFile in listFiles:
    scriptLines += ("date;\n")
#        scriptLines += ("dccp "+aFile+" inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n")
    scriptLines += ("cmsRun /user/hanwen/bphytest/CMSSW_8_0_26_patch1/src/BUAA-CMS-group/JPsiUpsilonPhi/test/"+cfgDirectory+"/runOnBatch_ntupleproduce_"+str(jobID)+"_cfg.py;\n")
#        scriptLines += ("rm inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n\n")
#        iteFileInJob = iteFileInJob+1
#    scriptLines += ('$ROOTSYS/bin/hadd output_'+name+"_"+str(jobID)+".root theOutput_"+name+"_"+str(jobID)+"_*.root;\n\n")
    scriptLines += ("cp Rootuple_JpsiPiPi_2016MiniAOD_"+str(jobID)+".root $INITDIR/"+cfgDirectory+"\n")
    scriptFile.write(scriptLines)
    scriptFile.close()

    jobsFiles = open('OUTPUTS/'+type+'/sendJobs_'+type+'.cmd',"a")
    jobsFiles.write("qsub -l walltime=20:00:00 -j oe "+os.getcwd()+'/'+jobsDirectory+'/runOnBatch_'+'ntupleproduce'+'_'+str(jobID)+'.sh\n')
    jobsFiles.close()

def main():

    jobId = 0
    global jobsDirectory
    global cfgDirectory
    if not os.path.exists('OUTPUTS/'+type+'/JOBS'): 
        os.makedirs('OUTPUTS/'+type+'/JOBS')
    if not os.path.exists('OUTPUTS/'+type+'/MERGED'):
        os.makedirs('OUTPUTS/'+type+'/MERGED')
    if not os.path.exists('OUTPUTS/'+type+'/CONFIGS'):
        os.makedirs('OUTPUTS/'+type+'/CONFIGS')
    if type == "mc_ele":
        printPath('/pnfs/iihe/cms/store/user/hanwen/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-001-Ele/CRAB_PrivateMC/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-TagProbePruner-Ele/180408_045629')
    if type == "mc_mu":
        printPath('/pnfs/iihe/cms/store/user/hanwen/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-001/CRAB_PrivateMC/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-TagProbePruner-Mu/180408_034333')    
    if type == "single_ele":
        printPath('/pnfs/iihe/cms/store/user/hanwen/SingleElectron001/CRAB_PrivateMC/crab_SingleElectron-TagProbePruner-Ele/180501_034315')
        printPath('/pnfs/iihe/cms/store/user/hanwen/SingleElectron002/CRAB_PrivateMC/crab_SingleElectron-TagProbePruner-Ele/180501_034717')
        printPath('/pnfs/iihe/cms/store/user/hanwen/SingleElectron003/CRAB_PrivateMC/crab_SingleElectron-TagProbePruner-Ele/180501_044954')
    if type == "single_mu":
        printPath('/pnfs/iihe/cms/store/user/hanwen/SingleMuon-001/CRAB_PrivateMC/crab_SingleMuon-TagProbePruner-Mu/180408_042749')
        printPath('/pnfs/iihe/cms/store/user/hanwen/SingleMuon-002/CRAB_PrivateMC/crab_SingleMuon-TagProbePruner-Mu/180408_042911') 
        printPath('/pnfs/iihe/cms/store/user/hanwen/SingleMuon-003/CRAB_PrivateMC/crab_SingleMuon-TagProbePruner-Mu/180408_042947')
    if type == "DoubleMuon_Run2016":
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2')
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016C/DoubleMuon/MINIAOD/03Feb2017-v1')
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1')
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016E/DoubleMuon/MINIAOD/03Feb2017-v1')
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016F/DoubleMuon/MINIAOD/03Feb2017-v1')
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1')
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016H/DoubleMuon/MINIAOD/03Feb2017_ver3-v1')
    if type == "Charm_Run2016":
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016E/Charmonium/MINIAOD/23Sep2016-v1')
        printPath('/pnfs/iihe/cms/ph/sc4/store/data/Run2016H/Charmonium/MINIAOD/03Feb2017_ver3-v1')
    jobsDirectory ='OUTPUTS/'+type+'/JOBS' 
    cfgDirectory = 'OUTPUTS/'+type +'/CONFIGS' 
    datasetFile = open("OUTPUTS/"+type+"/JOBS/catalog_"+type+".txt",'r')
    datasets = datasetFile.readlines()
    for line in datasets:
        jobId = jobId + 1
        prepare_job_script(line.strip(),jobId)
def harvestJobs():
    os.system("hadd "+cfgDirectory+"/Merged_"+type+".root "++cfgDirectory+"/*.root")
if __name__ == '__main__':
    type = sys.argv[1]
    harvest = sys.argv[2]
   # for arg in sys.argv:  
    #    print arg
    if not harvest == "harvest":
        main()
    if harvest == "harvest":
        harvestJobs()



