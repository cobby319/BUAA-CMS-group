#define ntuple_cxx
#include "../Includes/ntuple.h"
#include "../Includes/SmartSelectionMonitor.h"
#include "../Includes/SmartSelectionMonitor_hzz.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

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
   
   SmartSelectionMonitor_hzz mon;
   mon.declareHistos_jpsipipi();
   /*TH1F *h =(TH1F*) mon.addHistogram(new TH1F("eventflow",";;Events",6,0,6));
   h->GetXaxis()->SetBinLabel(1,"skimmed");
   h->GetXaxis()->SetBinLabel(2,"muonIDsoft");
   h->GetXaxis()->SetBinLabel(3,"J/PsiLxy");
   h->GetXaxis()->SetBinLabel(4,"J+PfitChi2");
   h->GetXaxis()->SetBinLabel(5,"J+2PfitChi2");
   h->GetXaxis()->SetBinLabel(6,"cosine<P,r>");*/





   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry % 10000 ==0) cout << jentry << " of " << nentries << endl;
      double weight = 1.;
      mon.fillHisto("Num_J/Psi","tot",nJ,weight);;

      auto smallestchi2 = std::min_element(Pi_vertexchisq2->begin(), Pi_vertexchisq2->end());
      int piN =std::distance(Pi_vertexchisq2->begin(), smallestchi2);
      cout <<  "selected pi is " << piN << endl;
      auto largestlxy = std::max_element(J_lxy->begin(), J_lxy->end());
      int jpsiN =std::distance(J_lxy->begin(), largestlxy);
      cout <<  "selected jpsi is " << jpsiN << endl;
      TLorentzVector jpsi, pion1,pion2;
      jpsi.SetXYZM(J_px->at(jpsiN),J_py->at(jpsiN),J_pz->at(jpsiN),J_mass->at(jpsiN)); 
      pion1.SetPtEtaPhiE(Pi_pt1->at(piN),Pi_eta1->at(piN),Pi_phi1->at(piN),Pi_e1->at(piN));
      pion2.SetPtEtaPhiE(Pi_pt2->at(piN),Pi_eta2->at(piN),Pi_phi2->at(piN),Pi_e2->at(piN));
      mon.fillHisto("M_J/Psi","tot",jpsi.M(),weight);;
      mon.fillHisto("pT_J/Psi","tot",jpsi.Pt(),weight);
      float jpipi_mass = (jpsi+pion1+pion2).M();
      if (jpipi_mass>4.1 &&jpipi_mass<5)  cout <<  "jpp mass is  " << jpipi_mass << endl;
      if (jpipi_mass>4.1 &&jpipi_mass<4.2)        mon.fillHisto("M_J/PsiPicut4.2","tot",(jpsi+pion1).M(),weight);
      else if (jpipi_mass>4.2 &&jpipi_mass<4.25)  mon.fillHisto("M_J/PsiPicut4.25","tot",(jpsi+pion1).M(),weight);
      else if (jpipi_mass>4.25 &&jpipi_mass<4.3)  mon.fillHisto("M_J/PsiPicut4.3","tot",(jpsi+pion1).M(),weight);
      else if (jpipi_mass>4.3 &&jpipi_mass<4.4)   mon.fillHisto("M_J/PsiPicut4.4","tot",(jpsi+pion1).M(),weight);
      else if (jpipi_mass>4.4 &&jpipi_mass<4.7)   mon.fillHisto("M_J/PsiPicut4.7","tot",(jpsi+pion1).M(),weight);
      else if (jpipi_mass>4.7 &&jpipi_mass<5.0)   mon.fillHisto("M_J/PsiPicut5.0","tot",(jpsi+pion1).M(),weight);
      // if (Cut(ientry) < 0) continue;
   }
   TFile* outFile=TFile::Open(outputFile_,"recreate");
    mon.Write();
    outFile->Close();
}
