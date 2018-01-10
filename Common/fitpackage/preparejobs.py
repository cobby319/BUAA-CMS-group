import sys
import re
import os
    
def prepare_job_script(i,j):
    global outputDirectory  
    os.mknod('JOBS/'+'mkhis'+'_'+'pt_'+str(i)+'eta_'+str(j)+'.sh') 
    scriptFile = open('JOBS/'+'mkhis'+'_'+'pt_'+str(i)+'eta_'+str(j)+'.sh','w')
    scriptLines = ''
    scriptLines += 'source $VO_CMS_SW_DIR/cmsset_default.sh\n'
    scriptLines += 'export SCRAM_ARCH=slc6_amd64_gcc530\n'
#    scriptLines += 'export BUILD_ARCH=slc6_amd64_gcc530\n'
#    scriptLines += 'export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n'
#    scriptLines += 'export XRD_NETWORKSTACK=IPv4\n'
    scriptLines += ('export INITDIR='+os.getcwd()+'\n')
    scriptLines += ('cd $INITDIR\n')
    scriptLines += 'eval `scramv1 runtime -sh`\n'
#    scriptLines += 'ulimit -c 0;\n'
    #scriptLines += 'cp '+os.getcwd()+'/makehistos .;\n'
    scriptLines += 'hostname ;\n'
    iteFileInJob=0
#    for aFile in listFiles:
    scriptLines += ("date;\n")
#        scriptLines += ("dccp "+aFile+" inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n")
    scriptLines += ("./makehistos pt="+str(i)+" eta="+str(j)+";\n")
    scriptLines += ("mv "+"histos_muons_data_pt"+str(i)+"_eta"+str(j)+".root Outputs;\n")

#        scriptLines += ("rm inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n\n")
#        iteFileInJob = iteFileInJob+1
#    scriptLines += ('$ROOTSYS/bin/hadd output_'+name+"_"+str(jobID)+".root theOutput_"+name+"_"+str(jobID)+"_*.root;\n\n")
    scriptFile.write(scriptLines)
    scriptFile.close()

    jobsFiles = open("sendJobs_"+outputDirectory+".cmd","a")
    jobsFiles.write("qsub -l walltime=20:00:00 -j oe "+os.getcwd()+'/JOBS/'+'mkhis'+'_'+'pt_'+str(i)+'eta_'+str(j)+'.sh\n')
    jobsFiles.close()
def main():


    global outputDirectory
    outputDirectory = "Outputs"
    for i in range(7):
        for j in range(7):
            prepare_job_script(i,j)
if __name__ == '__main__':
    main()


