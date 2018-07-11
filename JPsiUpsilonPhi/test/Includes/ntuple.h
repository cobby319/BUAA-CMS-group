//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 11 13:09:57 2018 by ROOT version 5.34/36
// from TTree ntuple/ J/psi ntuple
// found on file: C:/Users/allen/Documents/testnew4.root
//////////////////////////////////////////////////////////

#ifndef ntuple_h
#define ntuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class ntuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TString         outputFile_;
   UInt_t          nJ;
   UInt_t          nPiPair;
   vector<float>   *J_mass;
   vector<float>   *J_px;
   vector<float>   *J_py;
   vector<float>   *J_pz;
   vector<float>   *J_px1;
   vector<float>   *J_py1;
   vector<float>   *J_pz1;
   vector<int>     *J_charge1;
   vector<float>   *J_px2;
   vector<float>   *J_py2;
   vector<float>   *J_pz2;
   vector<int>     *J_charge2;
   vector<float>   *mumC2;
   vector<int>     *mumNHits;
   vector<int>     *mumNPHits;
   vector<float>   *mupC2;
   vector<int>     *mupNHits;
   vector<int>     *mupNPHits;
   vector<float>   *mumdxy;
   vector<float>   *mupdxy;
   vector<float>   *mumdz;
   vector<float>   *mupdz;
   vector<float>   *muon_dca;
   vector<bool>    *mu1soft;
   vector<bool>    *mu2soft;
   vector<bool>    *mu1tight;
   vector<bool>    *mu2tight;
   vector<bool>    *mu1PF;
   vector<bool>    *mu2PF;
   vector<bool>    *mu1loose;
   vector<bool>    *mu2loose;
   vector<float>   *J_lxy;
   vector<float>   *J_lxyErr;
   vector<float>   *Pi_dJP;
   vector<float>   *JPi_lxy;
   vector<float>   *JPi_lxyErr;
   vector<float>   *JPiPi_lxy;
   vector<float>   *JPiPi_lxyErr;
   vector<int>     *Pi_nhits1;
   vector<int>     *Pi_npixelhits1;
   vector<int>     *Pi_nhits2;
   vector<int>     *Pi_npixelhits2;
   vector<float>   *Pi_eta1;
   vector<float>   *Pi_eta2;
   vector<float>   *Pi_phi1;
   vector<float>   *Pi_phi2;
   vector<float>   *Pi_pt1;
   vector<float>   *Pi_pt2;
   vector<float>   *Pi_e1;
   vector<float>   *Pi_e2;
   vector<int>     *Pi_charge1;
   vector<int>     *Pi_charge2;
   vector<float>   *Pi_dxy1;
   vector<float>   *Pi_dxy2;
   vector<float>   *Pi_dxyerr1;
   vector<float>   *Pi_dxyerr2;
   vector<float>   *Pi_vertexchisq1;
   vector<float>   *Pi_vertexchisq2;

   // List of branches
   TBranch        *b_nJ;   //!
   TBranch        *b_nPiPair;   //!
   TBranch        *b_J_mass;   //!
   TBranch        *b_J_px;   //!
   TBranch        *b_J_py;   //!
   TBranch        *b_J_pz;   //!
   TBranch        *b_J_px1;   //!
   TBranch        *b_J_py1;   //!
   TBranch        *b_J_pz1;   //!
   TBranch        *b_J_charge1;   //!
   TBranch        *b_J_px2;   //!
   TBranch        *b_J_py2;   //!
   TBranch        *b_J_pz2;   //!
   TBranch        *b_J_charge2;   //!
   TBranch        *b_mumC2;   //!
   TBranch        *b_mumNHits;   //!
   TBranch        *b_mumNPHits;   //!
   TBranch        *b_mupC2;   //!
   TBranch        *b_mupNHits;   //!
   TBranch        *b_mupNPHits;   //!
   TBranch        *b_mumdxy;   //!
   TBranch        *b_mupdxy;   //!
   TBranch        *b_mumdz;   //!
   TBranch        *b_mupdz;   //!
   TBranch        *b_muon_dca;   //!
   TBranch        *b_mu1soft;   //!
   TBranch        *b_mu2soft;   //!
   TBranch        *b_mu1tight;   //!
   TBranch        *b_mu2tight;   //!
   TBranch        *b_mu1PF;   //!
   TBranch        *b_mu2PF;   //!
   TBranch        *b_mu1loose;   //!
   TBranch        *b_mu2loose;   //!
   TBranch        *b_J_lxy;   //!
   TBranch        *b_J_lxyErr;   //!
   TBranch        *b_Pi_dJP;   //!
   TBranch        *b_JPi_lxy;   //!
   TBranch        *b_JPi_lxyErr;   //!
   TBranch        *b_JPiPi_lxy;   //!
   TBranch        *b_JPiPi_lxyErr;   //!
   TBranch        *b_Pi_nhits1;   //!
   TBranch        *b_Pi_npixelhits1;   //!
   TBranch        *b_Pi_nhits2;   //!
   TBranch        *b_Pi_npixelhits2;   //!
   TBranch        *b_Pi_eta1;   //!
   TBranch        *b_Pi_eta2;   //!
   TBranch        *b_Pi_phi1;   //!
   TBranch        *b_Pi_phi2;   //!
   TBranch        *b_Pi_pt1;   //!
   TBranch        *b_Pi_pt2;   //!
   TBranch        *b_Pi_e1;   //!
   TBranch        *b_Pi_e2;   //!
   TBranch        *b_Pi_charge1;   //!
   TBranch        *b_Pi_charge2;   //!
   TBranch        *b_Pi_dxy1;   //!
   TBranch        *b_Pi_dxy2;   //!
   TBranch        *b_Pi_dxyerr1;   //!
   TBranch        *b_Pi_dxyerr2;   //!
   TBranch        *b_Pi_vertexchisq1;   //!
   TBranch        *b_Pi_vertexchisq2;   //!

   ntuple(TString fileName, TString outputFile);
   virtual ~ntuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntuple_cxx
ntuple::ntuple(TString fileName, TString outputFile) : fChain(0) 
{
   outputFile_ = outputFile;
   TChain * chain = new TChain("rootuple/ntuple","");
  chain->Add(fileName);
  TTree *tree = chain;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
#ifdef SINGLE_TREE
TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("rootuple/ntuple",tree);
#else // SINGLE_TREE
// The following code should be used if you want this class to access a chain
//       // of trees.
TChain * chain = new TChain("rootuple/ntuple","");
      chain->Add("DoubleMuon_jpsipipi.root/rootuple/ntuple");
      tree = chain;
#endif // SINGLE_TREE

}
   Init(tree);
}
ntuple::~ntuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ntuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   J_mass = 0;
   J_px = 0;
   J_py = 0;
   J_pz = 0;
   J_px1 = 0;
   J_py1 = 0;
   J_pz1 = 0;
   J_charge1 = 0;
   J_px2 = 0;
   J_py2 = 0;
   J_pz2 = 0;
   J_charge2 = 0;
   mumC2 = 0;
   mumNHits = 0;
   mumNPHits = 0;
   mupC2 = 0;
   mupNHits = 0;
   mupNPHits = 0;
   mumdxy = 0;
   mupdxy = 0;
   mumdz = 0;
   mupdz = 0;
   muon_dca = 0;
   mu1soft = 0;
   mu2soft = 0;
   mu1tight = 0;
   mu2tight = 0;
   mu1PF = 0;
   mu2PF = 0;
   mu1loose = 0;
   mu2loose = 0;
   J_lxy = 0;
   J_lxyErr = 0;
   Pi_dJP = 0;
   JPi_lxy = 0;
   JPi_lxyErr = 0;
   JPiPi_lxy = 0;
   JPiPi_lxyErr = 0;
   Pi_nhits1 = 0;
   Pi_npixelhits1 = 0;
   Pi_nhits2 = 0;
   Pi_npixelhits2 = 0;
   Pi_eta1 = 0;
   Pi_eta2 = 0;
   Pi_phi1 = 0;
   Pi_phi2 = 0;
   Pi_pt1 = 0;
   Pi_pt2 = 0;
   Pi_e1 = 0;
   Pi_e2 = 0;
   Pi_charge1 = 0;
   Pi_charge2 = 0;
   Pi_dxy1 = 0;
   Pi_dxy2 = 0;
   Pi_dxyerr1 = 0;
   Pi_dxyerr2 = 0;
   Pi_vertexchisq1 = 0;
   Pi_vertexchisq2 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nJ", &nJ, &b_nJ);
   fChain->SetBranchAddress("nPiPair", &nPiPair, &b_nPiPair);
   fChain->SetBranchAddress("J_mass", &J_mass, &b_J_mass);
   fChain->SetBranchAddress("J_px", &J_px, &b_J_px);
   fChain->SetBranchAddress("J_py", &J_py, &b_J_py);
   fChain->SetBranchAddress("J_pz", &J_pz, &b_J_pz);
   fChain->SetBranchAddress("J_px1", &J_px1, &b_J_px1);
   fChain->SetBranchAddress("J_py1", &J_py1, &b_J_py1);
   fChain->SetBranchAddress("J_pz1", &J_pz1, &b_J_pz1);
   fChain->SetBranchAddress("J_charge1", &J_charge1, &b_J_charge1);
   fChain->SetBranchAddress("J_px2", &J_px2, &b_J_px2);
   fChain->SetBranchAddress("J_py2", &J_py2, &b_J_py2);
   fChain->SetBranchAddress("J_pz2", &J_pz2, &b_J_pz2);
   fChain->SetBranchAddress("J_charge2", &J_charge2, &b_J_charge2);
   fChain->SetBranchAddress("mumC2", &mumC2, &b_mumC2);
   fChain->SetBranchAddress("mumNHits", &mumNHits, &b_mumNHits);
   fChain->SetBranchAddress("mumNPHits", &mumNPHits, &b_mumNPHits);
   fChain->SetBranchAddress("mupC2", &mupC2, &b_mupC2);
   fChain->SetBranchAddress("mupNHits", &mupNHits, &b_mupNHits);
   fChain->SetBranchAddress("mupNPHits", &mupNPHits, &b_mupNPHits);
   fChain->SetBranchAddress("mumdxy", &mumdxy, &b_mumdxy);
   fChain->SetBranchAddress("mupdxy", &mupdxy, &b_mupdxy);
   fChain->SetBranchAddress("mumdz", &mumdz, &b_mumdz);
   fChain->SetBranchAddress("mupdz", &mupdz, &b_mupdz);
   fChain->SetBranchAddress("muon_dca", &muon_dca, &b_muon_dca);
   fChain->SetBranchAddress("mu1soft", &mu1soft, &b_mu1soft);
   fChain->SetBranchAddress("mu2soft", &mu2soft, &b_mu2soft);
   fChain->SetBranchAddress("mu1tight", &mu1tight, &b_mu1tight);
   fChain->SetBranchAddress("mu2tight", &mu2tight, &b_mu2tight);
   fChain->SetBranchAddress("mu1PF", &mu1PF, &b_mu1PF);
   fChain->SetBranchAddress("mu2PF", &mu2PF, &b_mu2PF);
   fChain->SetBranchAddress("mu1loose", &mu1loose, &b_mu1loose);
   fChain->SetBranchAddress("mu2loose", &mu2loose, &b_mu2loose);
   fChain->SetBranchAddress("J_lxy", &J_lxy, &b_J_lxy);
   fChain->SetBranchAddress("J_lxyErr", &J_lxyErr, &b_J_lxyErr);
   fChain->SetBranchAddress("Pi_dJP", &Pi_dJP, &b_Pi_dJP);
   fChain->SetBranchAddress("JPi_lxy", &JPi_lxy, &b_JPi_lxy);
   fChain->SetBranchAddress("JPi_lxyErr", &JPi_lxyErr, &b_JPi_lxyErr);
   fChain->SetBranchAddress("JPiPi_lxy", &JPiPi_lxy, &b_JPiPi_lxy);
   fChain->SetBranchAddress("JPiPi_lxyErr", &JPiPi_lxyErr, &b_JPiPi_lxyErr);
   fChain->SetBranchAddress("Pi_nhits1", &Pi_nhits1, &b_Pi_nhits1);
   fChain->SetBranchAddress("Pi_npixelhits1", &Pi_npixelhits1, &b_Pi_npixelhits1);
   fChain->SetBranchAddress("Pi_nhits2", &Pi_nhits2, &b_Pi_nhits2);
   fChain->SetBranchAddress("Pi_npixelhits2", &Pi_npixelhits2, &b_Pi_npixelhits2);
   fChain->SetBranchAddress("Pi_eta1", &Pi_eta1, &b_Pi_eta1);
   fChain->SetBranchAddress("Pi_eta2", &Pi_eta2, &b_Pi_eta2);
   fChain->SetBranchAddress("Pi_phi1", &Pi_phi1, &b_Pi_phi1);
   fChain->SetBranchAddress("Pi_phi2", &Pi_phi2, &b_Pi_phi2);
   fChain->SetBranchAddress("Pi_pt1", &Pi_pt1, &b_Pi_pt1);
   fChain->SetBranchAddress("Pi_pt2", &Pi_pt2, &b_Pi_pt2);
   fChain->SetBranchAddress("Pi_e1", &Pi_e1, &b_Pi_e1);
   fChain->SetBranchAddress("Pi_e2", &Pi_e2, &b_Pi_e2);
   fChain->SetBranchAddress("Pi_charge1", &Pi_charge1, &b_Pi_charge1);
   fChain->SetBranchAddress("Pi_charge2", &Pi_charge2, &b_Pi_charge2);
   fChain->SetBranchAddress("Pi_dxy1", &Pi_dxy1, &b_Pi_dxy1);
   fChain->SetBranchAddress("Pi_dxy2", &Pi_dxy2, &b_Pi_dxy2);
   fChain->SetBranchAddress("Pi_dxyerr1", &Pi_dxyerr1, &b_Pi_dxyerr1);
   fChain->SetBranchAddress("Pi_dxyerr2", &Pi_dxyerr2, &b_Pi_dxyerr2);
   fChain->SetBranchAddress("Pi_vertexchisq1", &Pi_vertexchisq1, &b_Pi_vertexchisq1);
   fChain->SetBranchAddress("Pi_vertexchisq2", &Pi_vertexchisq2, &b_Pi_vertexchisq2);
   Notify();
}

Bool_t ntuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntuple_cxx
