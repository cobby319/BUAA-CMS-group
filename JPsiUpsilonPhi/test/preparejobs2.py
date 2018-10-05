import sys
import re
import os

def printPath(path):  
    if not os.path.exists("OUTPUTS/"+type+"/JOBS/catalog_"+type+".txt"): 
        os.mknod("OUTPUTS/"+type+"/JOBS/catalog_"+type+".txt")
    catalogfile=open("OUTPUTS/"+type+"/JOBS/catalog_"+type+".txt",'a')
    for file in os.listdir(path):
        if os.path.isfile(path+'/'+file):
            if '.root' in file:
                catalogfile.write('dcap://maite.iihe.ac.be'+path+'/'+file)
                catalogfile.write("\n")
        else :
            printPath(path+'/'+file)
    catalogfile.close()


def prepare_job_script(theCatalog,jobID):
    global cfgDirectory  
    

    if os.path.exists(jobsDirectory+'/runOnBatch_'+'scantuple'+'_'+str(jobID)+'.sh'): 
        os.remove(jobsDirectory+'/runOnBatch_'+'scantuple'+'_'+str(jobID)+'.sh')
    os.mknod(jobsDirectory+'/runOnBatch_'+'scantuple'+'_'+str(jobID)+'.sh') 



    scriptFile = open(jobsDirectory+'/runOnBatch_'+'scantuple'+'_'+str(jobID)+'.sh','w')
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
    scriptLines += 'cp '+os.getcwd()+'/scantuple .;\n'
    scriptLines += 'hostname ;\n'
    iteFileInJob=0
#    for aFile in listFiles:
    scriptLines += ("date;\n")
    scriptLines += ("./scantuple input="+theCatalog+" output=histos_"+str(jobID)+".root;\n")
 #   scriptLines += ("cmsRun /user/hanwen/bphytest/CMSSW_8_0_26_patch1/src/BUAA-CMS-group/JPsiUpsilonPhi/test/"+cfgDirectory+"/runOnBatch_scantuple_"+str(jobID)+"_cfg.py;\n")
#        scriptLines += ("rm inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n\n")
#        iteFileInJob = iteFileInJob+1
#    scriptLines += ('$ROOTSYS/bin/hadd output_'+name+"_"+str(jobID)+".root theOutput_"+name+"_"+str(jobID)+"_*.root;\n\n")
    scriptLines += ("cp histos_"+str(jobID)+".root $INITDIR/"+cfgDirectory+"\n")
    scriptFile.write(scriptLines)
    scriptFile.close()

    jobsFiles = open('OUTPUTS/'+type+'/sendJobs_'+type+'.cmd',"a")
    jobsFiles.write("qsub -l walltime=20:00:00 -j oe "+os.getcwd()+'/'+jobsDirectory+'/runOnBatch_'+'scantuple'+'_'+str(jobID)+'.sh\n')
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
    if type == "testold":
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium1005')
    if type == "test":
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium') 
    jobsDirectory ='OUTPUTS/'+type+'/JOBS' 
    cfgDirectory = 'OUTPUTS/'+type +'/CONFIGS' 
    datasetFile = open("OUTPUTS/"+type+"/JOBS/catalog_"+type+".txt",'r')
    datasets = datasetFile.readlines()
    for line in datasets:
        jobId = jobId + 1
        prepare_job_script(line.strip(),jobId)
def harvestJobs():
    global cfgDirectory 
    cfgDirectory = 'OUTPUTS/'+type +'/CONFIGS'
    os.system("hadd "+cfgDirectory+"/Merged.root "+cfgDirectory+"/*.root")
    os.system("mv "+cfgDirectory+"/Merged.root .")
def dotheFit(path):
    os.system("python bphyfit.py "+type+" "+path)
def drawDistributions():
    os.system("python plottery/distributions.py")
if __name__ == '__main__':
    type = sys.argv[1]
    harvest = sys.argv[2]
   # for arg in sys.argv:  
    #    print arg
    if not harvest == "harvest":
        main()
    if harvest == "harvest":
        harvestJobs()
    if harvest == "fit":
        dotheFit("Merged.root")
    if harvest == "plot":
        drawDistributions()


