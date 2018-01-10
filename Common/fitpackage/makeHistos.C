#include "../Includes/makeHistos.h"

void doTheHistos(TString nameHisto, TString testCut, TString baseCut){
    std::cout << "fill the histo:" << nameHisto << std::endl;
    
    
    TH1F *passing = new TH1F("passing","",40,60,120);
    chain->Draw("mass>>passing",baseCut+"&&"+testCut);
    myFile->cd();
    passing->Write(nameHisto+"_pass");
    delete passing;
    
    TH1F *failling = new TH1F("failling","",40,60,120);
    chain->Draw("mass>>failling",baseCut+"&&(!("+testCut+"))");
    myFile->cd();
    failling->Write(nameHisto+"_fail");
    delete failling;
    
    TH1F *all = new TH1F("all","",40,60,120);
    chain->Draw("mass",baseCut);
    myFile->cd();
    all->Write(nameHisto+"_all");
    delete all;    
}


void makeHistos(TString nameFile){



    chain->Add(nameFile);


    
    for (int i = 0 ; i < nbPtBins; i++){
        for (int j = 0 ; j < nbEtaBins ; j++){
            TString theCut = Form("pt>%f&&pt<%f&&abs(eta)>%f&&abs(eta)<%f",ptBins[i], ptBins[i+1], etaBins[j], etaBins[j+1]);
            doTheHistos(shortNameCut+Form("_%i_%i",i,j),theCutToTest,baseAcceptCut+"&&"+theCut);

        }
    }
    
    
    
    
}
