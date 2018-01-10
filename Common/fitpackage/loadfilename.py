#!/usr/bin/python  
# -*- coding:utf8 -*-  
  
import os  
allFileNum = 0 
dirname = os.getcwd()
if not os.path.exists('catalog.txt'): 
    os.mknod('catalog.txt')
else:
    os.remove('catalog.txt')
    os.mknod('catalog.txt')
file=open('catalog.txt','r+')
def printPath(path):  
    for dirs in os.listdir(path):
        for fs in os.listdir(path+'/'+dirs):
            if (os.path.isfile(path+'/'+dirs+'/'+fs)):
                file.write('dcap://maite.iihe.ac.be'+path+'/'+dirs+'/'+fs)
                file.write("\n")
    file.close()
if __name__ == '__main__':  
    printPath('/pnfs/iihe/cms/store/user/hanwen/DYJetsToLL_M-50/CRAB_PrivateMC/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-TagProbePruner-Mu/171205_234824')


    
    
    
