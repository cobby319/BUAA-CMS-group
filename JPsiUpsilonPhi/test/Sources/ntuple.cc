#define ntuple_cxx
#include "../Includes/ntuple.h"
#include "../Includes/SmartSelectionMonitor.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ntuple::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L ntuple.C
//      Root > ntuple t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   TFile *outFile = new TFile(outputFile_,"RECREATE");
   SmartSelectionMonitor mon;
   TH1F *h =(TH1F*) mon.addHistogram(new TH1F("eventflow",";;Events",6,0,6));
   h->GetXaxis()->SetBinLabel(1,"skimmed");
   h->GetXaxis()->SetBinLabel(2,"muonIDsoft");
   h->GetXaxis()->SetBinLabel(3,"J/PsiLxy");
   h->GetXaxis()->SetBinLabel(4,"J+PfitChi2");
   h->GetXaxis()->SetBinLabel(5,"J+2PfitChi2");
   h->GetXaxis()->SetBinLabel(6,"cosine<P,r>");

   mon.addHistogram(new TH1F("Num_J/Psi",";N_{J/#psi};Events",5,0,5));
   mon.addHistogram(new TH1F("M_J/Psi",";M_{J/#psi};Events",50,2.8,3.3));
   mon.addHistogram(new TH1F("pT_J/Psi",";p_{T,J/Psi};Events",120,0,120)); 


 
   mon.addHistogram(new TH1F("M_J/PsiPicut4.2",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
   mon.addHistogram(new TH1F("M_J/PsiPicut4.25",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
   mon.addHistogram(new TH1F("M_J/PsiPicut4.3",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
   mon.addHistogram(new TH1F("M_J/PsiPicut4.4",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
   mon.addHistogram(new TH1F("M_J/PsiPicut4.7",";m_{J/psi,Pi};Events",40,3.5,4.3)); 
   mon.addHistogram(new TH1F("M_J/PsiPicut5.0",";m_{J/psi,Pi};Events",40,3.5,4.3)); 




   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      double weight = 1.;
      mon.fillHisto("Num_J/Psi","tot",nJ,weight);
      
      // if (Cut(ientry) < 0) continue;
   }
   outFile->Write();
   outFile->Close();
}
