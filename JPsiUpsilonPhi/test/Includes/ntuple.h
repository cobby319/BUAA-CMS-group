//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 11 13:09:57 2018 by ROOT version 5.34/36
// from TTree ntuple/ J/psi ntuple
// found on file: C:/Users/allen/Documents/testnew4.root
//////////////////////////////////////////////////////////

#ifndef ntuple_h
#define ntuple_h

#include <iostream>
#include <fstream>
#include <thread>
#include <TROOT.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include "../Common/ObjectSelection.h"
// Header file for the classes stored in the TTree if any.
#include "vector"
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class ntuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   int maxEvents_;
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
   vector<float>   *J_vertexchi2;
   vector<float>   *Pi_dJP;
   vector<float>   *JPi_lxy;
   vector<float>   *JPi_lxyErr;
   vector<float>   *JPiPi_lxy;
   vector<float>   *JPiPi_lxyErr;
   vector<float>   *JPiPi_x;
   vector<float>   *JPiPi_y;
   vector<float>   *JPiPi_z;
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
   vector<float>       *Pi1_hcalFraction;
   vector<float>       *Pi2_hcalFraction;
   vector<float>       *Pi1_vertexNdof;
   vector<float>       *Pi2_vertexNdof;
   vector<float>       *Pi1_vertexNchi2;
   vector<float>       *Pi2_vertexNchi2;
   vector<float>       *Pi1_lambda;
   vector<float>       *Pi2_lambda;
   vector<float>       *Pi1_lambdaError;
   vector<float>       *Pi2_lambdaError;
   vector<float>       *Pi1_qoverp;
   vector<float>       *Pi2_qoverp;
   vector<float>       *Pi1_qoverpError;
   vector<float>       *Pi2_qoverpError;
   vector<float>       *Pi1_validTkFraction;
   vector<float>       *Pi2_validTkFraction;
   vector<int>         *Pi1_numberOfMothers;
   vector<int>         *Pi2_numberOfMothers;
   vector<int>         *Pi1_numberOfSourceCandidatePtrs;
   vector<int>         *Pi2_numberOfSourceCandidatePtrs;
   vector<int>         *Pi1_pdgId;
   vector<int>         *Pi2_pdgId;
   vector<int>         *Pi1_numberOfValidHitsOnTrack;
   vector<int>         *Pi2_numberOfValidHitsOnTrack;
   vector<int>         *Pi1_innerDetId;
   vector<int>         *Pi2_innerDetId;
   vector<bool>        *Pi1_innerOk;
   vector<bool>        *Pi2_innerOk;
   vector<bool>        *Pi1_isCaloMuon;
   vector<bool>        *Pi2_isCaloMuon;
   vector<bool>        *Pi1_isConvertedPhoton;
   vector<bool>        *Pi2_isConvertedPhoton;
   vector<bool>        *Pi1_isElectron;
   vector<bool>        *Pi2_isElectron;
   vector<bool>        *Pi1_isMuon ;
   vector<bool>        *Pi2_isMuon;
   vector<bool>        *Pi1_isPhoton;
   vector<bool>        *Pi2_isPhoton ;
   vector<bool>        *Pi1_isGlobalMuon;
   vector<bool>        *Pi2_isGlobalMuon;
   vector<bool>        *Pi1_isJet;
   vector<bool>        *Pi2_isJet;
   vector<bool>        *Pi1_isLonglived;
   vector<bool>        *Pi2_isLonglived;
   vector<bool>        *Pi1_massConstraint;
   vector<bool>        *Pi2_massConstraint;
   vector<float>       *J_Prob;
   vector<float>       *JPi_Prob;
   vector<float>       *JPiPi_Prob;

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
   TBranch        *b_J_vertexchi2; //!
   TBranch        *b_Pi_dJP;   //!
   TBranch        *b_JPi_lxy;   //!
   TBranch        *b_JPi_lxyErr;   //!
   TBranch        *b_JPiPi_lxy;   //!
   TBranch        *b_JPiPi_lxyErr;   //!
   TBranch        *b_JPiPi_x;
   TBranch        *b_JPiPi_y;
   TBranch        *b_JPiPi_z;
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
   TBranch        *b_Pi1_hcalFraction;
   TBranch        *b_Pi2_hcalFraction;
   TBranch        *b_Pi1_vertexNdof;
   TBranch        *b_Pi2_vertexNdof;
   TBranch        *b_Pi1_vertexNchi2;
   TBranch        *b_Pi2_vertexNchi2;
   TBranch        *b_Pi1_lambda;
   TBranch        *b_Pi2_lambda;
   TBranch        *b_Pi1_lambdaError;
   TBranch        *b_Pi2_lambdaError;
   TBranch        *b_Pi1_qoverp;
   TBranch        *b_Pi2_qoverp;
   TBranch        *b_Pi1_qoverpError;
   TBranch        *b_Pi2_qoverpError;
   TBranch        *b_Pi1_validTkFraction;
   TBranch        *b_Pi2_validTkFraction;
   TBranch        *b_Pi1_numberOfMothers;
   TBranch        *b_Pi2_numberOfMothers;
   TBranch        *b_Pi1_numberOfSourceCandidatePtrs;
   TBranch        *b_Pi2_numberOfSourceCandidatePtrs;
   TBranch        *b_Pi1_pdgId;
   TBranch        *b_Pi2_pdgId;
   TBranch        *b_Pi1_numberOfValidHitsOnTrack;
   TBranch        *b_Pi2_numberOfValidHitsOnTrack;
   TBranch        *b_Pi1_innerDetId;
   TBranch        *b_Pi2_innerDetId;
   TBranch        *b_Pi1_innerOk;
   TBranch        *b_Pi2_innerOk;
   TBranch        *b_Pi1_isCaloMuon;
   TBranch        *b_Pi2_isCaloMuon;
   TBranch        *b_Pi1_isConvertedPhoton;
   TBranch        *b_Pi2_isConvertedPhoton;
   TBranch        *b_Pi1_isElectron;
   TBranch        *b_Pi2_isElectron;
   TBranch        *b_Pi1_isMuon ;
   TBranch        *b_Pi2_isMuon;
   TBranch        *b_Pi1_isPhoton;
   TBranch        *b_Pi2_isPhoton ;
   TBranch        *b_Pi1_isGlobalMuon;
   TBranch        *b_Pi2_isGlobalMuon;
   TBranch        *b_Pi1_isJet;
   TBranch        *b_Pi2_isJet;
   TBranch        *b_Pi1_isLonglived;
   TBranch        *b_Pi2_isLonglived;
   TBranch        *b_Pi1_massConstraint;
   TBranch        *b_Pi2_massConstraint;
   TBranch        *b_J_Prob;
   TBranch        *b_JPi_Prob;
   TBranch        *b_JPiPi_Prob;

   ntuple(TString fileName, TString outputFile, int skipFile, int maxFiles,int CatalogOrFile);
   virtual ~ntuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     FillTheTChain(TChain *theChain, TString theInputCatalog, int skipFiles, int maxFiles);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntuple_cxx
ntuple::ntuple(TString fileName, TString outputFile, int skipFile, int maxFiles,int CatalogOrFile) : fChain(0) 
{
   
  outputFile_ = outputFile;

  
   TChain * chain = new TChain("rootuple/ntuple","");
   if      (CatalogOrFile)  FillTheTChain(chain, fileName, skipFile, maxFiles);
   else if (!CatalogOrFile)   chain->Add(fileName);
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
      chain->Add("test.root/rootuple/ntuple");
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

void ntuple::FillTheTChain(TChain *theChain, TString theInputCatalog, int skipFiles, int maxFiles){
  cout << "catalog name=" << theInputCatalog << endl;

  std::ifstream f(theInputCatalog);
  if(!f.good()){
    std::cerr << "Failed to open file "<< theInputCatalog << "!\n";
    return;
  }


  int iline = 0;
  int nfiles = 0;
  std::string firstFile_ = "";
  while(f.good()){
    ++iline;
    std::string l;
    std::string::size_type p;

    std::getline(f, l);

    //trim white spaces:
    p = l.find_first_not_of(" \t");
    if(p!=std::string::npos) l.erase(0, p);
    p = l.find_last_not_of(" \t\n\r");
    if(p!=std::string::npos) l.erase(p + 1);
    else l.clear();

    //skip empty lines and comment lines:
    if (!l.size() || l[0] == '#' || l[0] == '*') continue;

    //extract first column (file name):
    p = l.find_first_of(" \t");
    if(p!=std::string::npos) l.erase(p);

    //sanity check:
    const char ext[6] = ".root";

    if(l.size() < sizeof(ext) || l.substr(l.size() - sizeof(ext) + 1) != ext){
      std::cerr << "Line " << iline << " of catalog file " << theInputCatalog << " was skipped.\n";
      continue;
    }



    if(skipFiles <= 0){
      ++nfiles;
      if((maxFiles > 0) &&  (nfiles > maxFiles)) break;
      std::cout << "Add file " << l.c_str() << " to the list of input files.\n";
      theChain->Add(l.c_str());
      if(firstFile_.size()==0) firstFile_ = l;
    } else{
    --skipFiles;
  }
}

return ;

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
   J_vertexchi2 =0;
   Pi_dJP = 0;
   JPi_lxy = 0;
   JPi_lxyErr = 0;
   JPiPi_lxy = 0;
   JPiPi_lxyErr = 0;
   JPiPi_x = 0;
   JPiPi_y = 0;
   JPiPi_z = 0;
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
   Pi1_numberOfSourceCandidatePtrs = 0;
   Pi2_numberOfSourceCandidatePtrs = 0;
   Pi1_hcalFraction = 0;
   Pi2_hcalFraction = 0;
   Pi1_vertexNdof = 0;
   Pi2_vertexNdof = 0;
   Pi1_vertexNchi2 = 0;
   Pi2_vertexNchi2 = 0;
   Pi1_lambda = 0;
   Pi2_lambda = 0;
   Pi1_lambdaError = 0;
   Pi2_lambdaError = 0;
   Pi1_qoverp = 0;
   Pi2_qoverp = 0;
   Pi1_qoverpError = 0;
   Pi2_qoverpError = 0;
   Pi1_validTkFraction = 0;
   Pi2_validTkFraction = 0;
   Pi1_numberOfMothers = 0;
   Pi2_numberOfMothers = 0;
   Pi1_pdgId = 0;
   Pi2_pdgId = 0;
   Pi1_numberOfValidHitsOnTrack = 0;
   Pi2_numberOfValidHitsOnTrack = 0;
   Pi1_innerDetId = 0;
   Pi2_innerDetId = 0;
   Pi1_innerOk = 0;
   Pi2_innerOk = 0;
   Pi1_isCaloMuon = 0;
   Pi2_isCaloMuon = 0;
   Pi1_isConvertedPhoton = 0;
   Pi2_isConvertedPhoton = 0;
   Pi1_isElectron = 0;
   Pi2_isElectron = 0;
   Pi1_isMuon  = 0;
   Pi2_isMuon = 0;
   Pi1_isPhoton = 0;
   Pi2_isPhoton  = 0;
   Pi1_isGlobalMuon = 0;
   Pi2_isGlobalMuon = 0;
   Pi1_isJet = 0;
   Pi2_isJet = 0;
   Pi1_isLonglived = 0;
   Pi2_isLonglived = 0;
   Pi1_massConstraint = 0;
   Pi2_massConstraint = 0;
   J_Prob = 0;
   JPi_Prob = 0;
   JPiPi_Prob = 0;
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
   fChain->SetBranchAddress("J_vertexchi2",&J_vertexchi2,&b_J_vertexchi2);
   fChain->SetBranchAddress("Pi_dJP", &Pi_dJP, &b_Pi_dJP);
   fChain->SetBranchAddress("JPi_lxy", &JPi_lxy, &b_JPi_lxy);
   fChain->SetBranchAddress("JPi_lxyErr", &JPi_lxyErr, &b_JPi_lxyErr);
   fChain->SetBranchAddress("JPiPi_lxy", &JPiPi_lxy, &b_JPiPi_lxy);
   fChain->SetBranchAddress("JPiPi_lxyErr", &JPiPi_lxyErr, &b_JPiPi_lxyErr);
   fChain->SetBranchAddress("JPiPi_x",&JPiPi_x,&b_JPiPi_x);
   fChain->SetBranchAddress("JPiPi_y",&JPiPi_y,&b_JPiPi_y);
   fChain->SetBranchAddress("JPiPi_z",&JPiPi_z,&b_JPiPi_z);
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
   fChain->SetBranchAddress("Pi1_numberOfSourceCandidatePtrs".&Pi1_numberOfSourceCandidatePtrs,&b_Pi1_numberOfSourceCandidatePtrs);
   fChain->SetBranchAddress("Pi2_numberOfSourceCandidatePtrs",&Pi2_numberOfSourceCandidatePtrs,&b_Pi2_numberOfSourceCandidatePtrs);
   fChain->SetBranchAddress("Pi1_hcalFraction",&Pi1_hcalFraction,&b_Pi1_hcalFraction);
   fChain->SetBranchAddress("Pi2_hcalFraction",&Pi2_hcalFraction,&b_Pi2_hcalFraction);
   fChain->SetBranchAddress("Pi1_vertexNdof",&Pi1_vertexNdof,&b_Pi1_vertexNdof);
   fChain->SetBranchAddress("Pi2_vertexNdof",&Pi2_vertexNdof,&b_Pi2_vertexNdof);
   fChain->SetBranchAddress("Pi1_vertexNchi2",&Pi1_vertexNchi2,&b_Pi1_vertexNchi2);
   fChain->SetBranchAddress("Pi2_vertexNchi2",&Pi2_vertexNchi2,&b_Pi2_vertexNchi2);
   fChain->SetBranchAddress("Pi1_lambda",&Pi1_lambda,&b_Pi1_lambda);
   fChain->SetBranchAddress("Pi2_lambda",&Pi2_lambda,&b_Pi2_lambda);
   fChain->SetBranchAddress("Pi1_lambdaError",&Pi1_lambdaError,&b_Pi1_lambdaError);
   fChain->SetBranchAddress("Pi2_lambdaError",&Pi2_lambdaError,&b_Pi2_lambdaError);
   fChain->SetBranchAddress("Pi1_qoverp",&Pi1_qoverp,&b_Pi1_qoverp);
   fChain->SetBranchAddress("Pi2_qoverp",&Pi2_qoverp,&b_Pi2_qoverp);
   fChain->SetBranchAddress("Pi1_qoverpError",&Pi1_qoverpError,&b_Pi1_qoverpError);
   fChain->SetBranchAddress("Pi2_qoverpError",&Pi2_qoverpError,&b_Pi2_qoverpError);
   fChain->SetBranchAddress("Pi1_validTkFraction",&Pi1_validTkFraction,&b_Pi1_validTkFraction);
   fChain->SetBranchAddress("Pi2_validTkFraction",&Pi2_validTkFraction,&b_Pi2_validTkFraction);
   fChain->SetBranchAddress("Pi1_numberOfMothers",&Pi1_numberOfMothers,&b_Pi1_numberOfMothers);
   fChain->SetBranchAddress("Pi2_numberOfMothers",&Pi2_numberOfMothers,&b_Pi2_numberOfMothers);
   fChain->SetBranchAddress("Pi1_pdgId",&Pi1_pdgId,&b_Pi1_pdgId);
   fChain->SetBranchAddress("Pi2_pdgId",&Pi2_pdgId,&b_Pi2_pdgId);
   fChain->SetBranchAddress("Pi1_numberOfValidHitsOnTrack",&Pi1_numberOfValidHitsOnTrack,&b_Pi1_numberOfValidHitsOnTrack);
   fChain->SetBranchAddress("Pi2_numberOfValidHitsOnTrack",&Pi2_numberOfValidHitsOnTrack,&b_Pi2_numberOfValidHitsOnTrack);
   fChain->SetBranchAddress("Pi1_innerDetId",&Pi1_innerDetId,&b_Pi1_innerDetId);
   fChain->SetBranchAddress("Pi2_innerDetId",&Pi2_innerDetId,&b_Pi2_innerDetId);
   fChain->SetBranchAddress("Pi1_innerOk",&Pi1_innerOk,&b_Pi1_innerOk);
   fChain->SetBranchAddress("Pi2_innerOk",&Pi2_innerOk,&b_Pi2_innerOk);
   fChain->SetBranchAddress("Pi1_isCaloMuon",&Pi1_isCaloMuon,&b_Pi1_isCaloMuon);
   fChain->SetBranchAddress("Pi2_isCaloMuon",&Pi2_isCaloMuon,&b_Pi2_isCaloMuon);
   fChain->SetBranchAddress("Pi1_isConvertedPhoton",&Pi1_isConvertedPhoton,&b_Pi1_isConvertedPhoton);
   fChain->SetBranchAddress("Pi2_isConvertedPhoton",&Pi2_isConvertedPhoton,&b_Pi2_isConvertedPhoton);
   fChain->SetBranchAddress("Pi1_isElectron",&Pi1_isElectron,&b_Pi1_isElectron);
   fChain->SetBranchAddress("Pi2_isElectron",&Pi2_isElectron,&b_Pi2_isElectron);
   fChain->SetBranchAddress("Pi1_isMuon",&Pi1_isMuon,&b_Pi1_isMuon);
   fChain->SetBranchAddress("Pi2_isMuon",&Pi2_isMuon,&b_Pi2_isMuon);
   fChain->SetBranchAddress("Pi1_isPhoton",&Pi1_isPhoton,&b_Pi1_isPhoton);
   fChain->SetBranchAddress("Pi2_isPhoton",&Pi2_isPhoton,&b_Pi2_isPhoton);
   fChain->SetBranchAddress("Pi1_isGlobalMuon",&Pi1_isGlobalMuon,&b_Pi1_isGlobalMuon);
   fChain->SetBranchAddress("Pi2_isGlobalMuon",&Pi2_isGlobalMuon,&b_Pi2_isGlobalMuon);
   fChain->SetBranchAddress("Pi1_isJet",&Pi1_isJet,&b_Pi1_isJet);
   fChain->SetBranchAddress("Pi2_isJet",&Pi2_isJet,&b_Pi2_isJet);
   fChain->SetBranchAddress("Pi1_isLonglived",&Pi1_isLonglived,&b_Pi1_isLonglived);
   fChain->SetBranchAddress("Pi2_isLonglived",&Pi2_isLonglived,&b_Pi2_isLonglived);
   fChain->SetBranchAddress("Pi1_massConstraint",&Pi1_massConstraint,&b_Pi1_massConstraint);
   fChain->SetBranchAddress("Pi2_massConstraint",&Pi2_massConstraint,&b_Pi2_massConstraint);
   fChain->SetBranchAddress("J_Prob",&J_Prob,&b_J_Prob);
   fChain->SetBranchAddress("JPi_Prob",&JPi_Prob,&b_JPi_Prob);
   fChain->SetBranchAddress("JPiPi_Prob",&JPiPi_Prob,&b_JPiPi_Prob);
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
