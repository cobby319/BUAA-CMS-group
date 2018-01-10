//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 27 08:16:15 2017 by ROOT version 6.06/01
// from TTree fitter_tree/fitter_tree
// found on file: output_tptree_2.root
//////////////////////////////////////////////////////////

#ifndef fitter_tree_h
#define fitter_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class fitter_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         mass;
   Float_t         tag_pt;
   Float_t         pt;
   Float_t         tag_eta;
   Float_t         eta;
   Float_t         tag_phi;
   Float_t         phi;
   Float_t         tag_abseta;
   Float_t         abseta;
   Float_t         tag_charge;
   Float_t         charge;
   Float_t         tkIso;
   Int_t           Tight2012;
   Int_t           tag_HltMatch;
   Int_t           HltMatch;

   // List of branches
   TBranch        *b_mass;   //!
   TBranch        *b_tag_pt;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_tag_eta;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_tag_phi;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_tag_abseta;   //!
   TBranch        *b_abseta;   //!
   TBranch        *b_tag_charge;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_tkIso;   //!
   TBranch        *b_Tight2012;   //!
   TBranch        *b_tag_HltMatch;   //!
   TBranch        *b_HltMatch;   //!

   fitter_tree(TTree *tree=0);
   virtual ~fitter_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fitter_tree_cxx
fitter_tree::fitter_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_tptree_2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_tptree_2.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("output_tptree_2.root:/tpTree");
      dir->GetObject("fitter_tree",tree);

   }
   Init(tree);
}

fitter_tree::~fitter_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fitter_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fitter_tree::LoadTree(Long64_t entry)
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

void fitter_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("tag_pt", &tag_pt, &b_tag_pt);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("tag_eta", &tag_eta, &b_tag_eta);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("tag_phi", &tag_phi, &b_tag_phi);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("tag_abseta", &tag_abseta, &b_tag_abseta);
   fChain->SetBranchAddress("abseta", &abseta, &b_abseta);
   fChain->SetBranchAddress("tag_charge", &tag_charge, &b_tag_charge);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("tkIso", &tkIso, &b_tkIso);
   fChain->SetBranchAddress("Tight2012", &Tight2012, &b_Tight2012);
   fChain->SetBranchAddress("tag_HltMatch", &tag_HltMatch, &b_tag_HltMatch);
   fChain->SetBranchAddress("HltMatch", &HltMatch, &b_HltMatch);
   Notify();
}

Bool_t fitter_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fitter_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fitter_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fitter_tree_cxx
