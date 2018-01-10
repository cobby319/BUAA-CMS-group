#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"

TChain *chain = new TChain("tpTree/fitter_tree");
//TFile *myFile = new TFile("histos_muons_data.root","RECREATE");
TString baseAcceptCut = "tag_pt>25&&abs(tag_eta)<2.1&&mass>60&&mass<120";
TString theCutToTest = "Tight2012==1";
TString shortNameCut = "isTight";

int nbPtBins = 7;
float ptBins[8] = {20, 25,30,40, 50,60,120, 200};
int nbEtaBins =7 ;
float etaBins[8] = {0, 0.2,0.3,0.9,1.2,1.6,2.1,2.4};
