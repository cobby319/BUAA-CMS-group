import sys
import re
import os

def printPath(path):  
    dirname = os.getcwd()
    if not os.path.exists("catalog_"+type+".txt"): 
        os.mknod("catalog_"+type+".txt")
    file=open("catalog_"+type+".txt",'a')
    for dirs in os.listdir(path):
        for fs in os.listdir(path+'/'+dirs):
            if (os.path.isfile(path+'/'+dirs+'/'+fs)):
                file.write('dcap://maite.iihe.ac.be'+path+'/'+dirs+'/'+fs)
                file.write("\n")
    file.close()



def main():
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
    if type == "Run2016":
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium/crab_Charmonium_Run2016B-03Feb2017_ver2-v2/180714_004045')
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium/crab_Charmonium_Run2016C-03Feb2017-v1/180714_004026')
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium/crab_Charmonium_Run2016D-03Feb2017-v1/180714_004054')
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium/crab_Charmonium_Run2016E-03Feb2017-v1/180714_004007')
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium/crab_Charmonium_Run2016F-03Feb2017-v1/180714_004017')
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium/crab_Charmonium_Run2016G-03Feb2017-v1/180714_004036')
        printPath('/pnfs/iihe/cms/store/user/hanwen/Charmonium/crab_Charmonium_Run2016H-03Feb2017_ver3-v1/180713_140833')


    

if __name__ == '__main__':
    type = sys.argv[1]
    main()
    



