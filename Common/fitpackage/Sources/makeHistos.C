#include "../Includes/makeHistos.h"

void doTheHistos(TString nameHisto, TString testCut, TString baseCut,int i, int j){
    std::cout << "fill the histo:" << nameHisto << std::endl;
    TString fileName = Form("histos_muons_mc_pt%i_eta%i.root",i,j);
    TFile *myFile = new TFile(fileName,"RECREATE");
    
    TH1F *passing = new TH1F("passing","",40,60,120);
    chain->Draw("mass>>passing",baseCut+"&&"+testCut);
    myFile->cd();
    passing->Write(nameHisto+"_pass");

    
    TH1F *failling = new TH1F("failling","",40,60,120);
    chain->Draw("mass>>failling",baseCut+"&&(!("+testCut+"))");
    myFile->cd();
    failling->Write(nameHisto+"_fail");

    
    TH1F *all = new TH1F("all","",40,60,120);
    all->Add(passing,failling);
    myFile->cd();
    all->Write(nameHisto+"_all");
    delete all;    
    delete failling;
    delete passing;
}


void makeHistos(TString nameFile,int i, int j){



    chain->Add(nameFile);


    


         TString theCut = Form("pt>%f&&pt<%f&&abs(eta)>%f&&abs(eta)<%f",ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1]);
         doTheHistos(shortNameCut+Form("_%i_%i",i,j),theCutToTest,baseAcceptCut+"&&"+theCut,i,j);



    
    
    
    
}
