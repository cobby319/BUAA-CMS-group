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
      if(jentry % 10000 ==0) {
         cout <<  "selected piN-1 is " << piN -1 << "and selected chi2 at piN-1 is "<< Pi_vertexchisq2->at(piN-1)<<endl;
         cout <<  "piN  is " << piN << "and selected chi2 at piN is "<< Pi_vertexchisq2->at(piN)<< endl;
         }
      auto minchi2J = std::min_element(J_vertexchi2->begin(), J_vertexchi2->end());
      int jpsiN =std::distance(J_vertexchi2->begin(), minchi2J) ;
      //cout <<  "selected jpsi is " << jpsiN << endl;
      TLorentzVector jpsi, pion1,pion2;
      float Pion_mass = 0.13957061;
      jpsi.SetXYZM(J_px->at(jpsiN),J_py->at(jpsiN),J_pz->at(jpsiN),J_mass->at(jpsiN)); 
      pion1.SetPtEtaPhiM(Pi_pt1->at(piN),Pi_eta1->at(piN),Pi_phi1->at(piN),Pion_mass);
      pion2.SetPtEtaPhiM(Pi_pt2->at(piN),Pi_eta2->at(piN),Pi_phi2->at(piN),Pion_mass);
      mon.fillHisto("M_JPsi","tot",jpsi.M(),weight);
      mon.fillHisto("pT_JPsi","tot",jpsi.Pt(),weight);
      mon.fillHisto("Deta JpsiPi1","tot",jpsi.Eta()-pion1.Eta(),weight); 
      mon.fillHisto("Deta JpsiPi2","tot",jpsi.Eta()-pion2.Eta(),weight); 
      mon.fillHisto("Dphi JpsiPi1","tot",jpsi.Phi()-pion1.Phi(),weight); 
      mon.fillHisto("Dphi JpsiPi2","tot",jpsi.Phi()-pion2.Phi(),weight); 
      mon.fillHisto("lxy_jpsipipi","tot",JPiPi_lxy->at(piN),weight); 
      float jpipi_mass = (jpsi+pion1+pion2).M();
      float px,py,pz;
      float rx,ry,rz;
      px = (jpsi+pion1+pion2).Px();
      py = (jpsi+pion1+pion2).Py();
      pz = (jpsi+pion1+pion2).Pz();
      rx = JPiPi_x->at(piN);
      ry = JPiPi_y->at(piN);
      rz = JPiPi_z->at(piN);
      float cosine = (px*rx+py*ry+pz*rz)/((px*px+py*py+pz*pz)*(rx*rx+ry*ry+rz*rz));
      mon.fillHisto("cosine of P&r","total",cosine,weight);
      float pipi_mass =(pion1+pion2).M();
      if (pipi_mass>0.81 && pipi_mass<0.97) continue;
      if (pipi_mass>1.01 && pipi_mass<1.03) continue;
      if (pipi_mass<0.35) continue;
      mon.fillHisto("M_JPsiPiPi3-8","total",jpipi_mass,weight);
      if (jpipi_mass>4.0 &&jpipi_mass<5.0)   mon.fillHisto("M_JPsiPiPi4-5","total",jpipi_mass,weight);
      mon.fillHisto("M_JpsiPi1&M_JpsiPi2","total",(jpsi+pion1).M(),(jpsi+pion2).M(),weight);

     // if (jpipi_mass>4.1 &&jpipi_mass<5)  cout <<  "jpp mass is  " << jpipi_mass << endl;
      if (jpipi_mass>4.1 &&jpipi_mass<4.2)         mon.fillHisto("M_JPsiPicut4.1-4.2","total",(jpsi+pion1).M(),weight);
      if (jpipi_mass>4.2 &&jpipi_mass<4.25)         mon.fillHisto("M_JPsiPicut4.2-4.25","total",(jpsi+pion1).M(),weight);
      if (jpipi_mass>4.25 &&jpipi_mass<4.3)          mon.fillHisto("M_JPsiPicut4.25-4.3","total",(jpsi+pion1).M(),weight);
      if (jpipi_mass>4.3 &&jpipi_mass<4.4)          mon.fillHisto("M_JPsiPicut4.3-4.4","total",(jpsi+pion1).M(),weight);
      if (jpipi_mass>4.4 &&jpipi_mass<4.7)          mon.fillHisto("M_JPsiPicut4.4-4.7","total",(jpsi+pion1).M(),weight);
      if (jpipi_mass>4.7 &&jpipi_mass<5.0)        mon.fillHisto("M_JPsiPicut4.7-5.0","total",(jpsi+pion1).M(),weight);
      //if (jpipi_mass>4.24 &&jpipi_mass<4.28)        mon.fillHisto("M_JPsiPicut4.24-4.28","total",(jpsi+pion1).M(),weight);

      // if (Cut(ientry) < 0) continue;
   }
   TFile* outFile=TFile::Open(outputFile_,"recreate");
   TDirectoryFile *dir = new TDirectoryFile("histos","histos","",outFile);
   gDirectory->cd("histos");
   mon.Write();
   outFile->Close();
}
