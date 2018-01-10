import os
if not os.path.exists('resendjobs.cmd'):
    os.mknod('resendjobs.cmd')
else:
    os.remove('resendjobs.cmd')
    os.mknod('resendjobs.cmd')
file=open('resendjobs.cmd','r+')
def printPath(path):
    for i in range(7):
        for j in range(7):
            if not(os.path.exists("/user/hanwen/work/CMSSW_8_0_26_patch1/src/shears_old/fit/test/Outputs/histos_muons_mc_pt"+str(i)+"_eta"+str(j)+".root")):
                    file.write("qsub -l walltime=20:00:00 -j oe /storage_mnt/storage/user/hanwen/work/CMSSW_8_0_26_patch1/src/shears_old/fit/test/JOBS/mkhis_pt_"+str(i)+"eta_"+str(j)+".sh")
                    file.write("\n")
    file.close()
if __name__ == '__main__':
    printPath("/user/hanwen/work/CMSSW_8_0_26_patch1/src/shears_old/fit/test/Outputs")
