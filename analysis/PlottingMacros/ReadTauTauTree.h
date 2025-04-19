//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 13 17:23:41 2023 by ROOT version 6.26/06
// from TTree output_tree/
// found on file: RootFile/mc_GammaUPC_EleMu_ArashNew_13thSep_moreEvents_merged.root
//////////////////////////////////////////////////////////

#ifndef ReadTauTauTree_h
#define ReadTauTauTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ReadTauTauTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         Tau_Pt1;
   Float_t         Tau_Pt;
   Float_t         Tau_Eta;
   Float_t         Tau_Phi;
   Float_t         Tau_E;
   Float_t         Gen_Ele_Pt;
   Float_t         Gen_Ele_Eta;
   Float_t         Gen_Ele_Phi;
   Float_t         Gen_Ele_E;
   Float_t         Gen_Mu_Pt;
   Float_t         Gen_Mu_Eta;
   Float_t         Gen_Mu_Phi;
   Float_t         Gen_Mu_E;
   Float_t         Gen_vSum_M;
   Float_t         Gen_vSum_Pt;
   Float_t         Gen_Scalar_Pt;
   Float_t         Gen_vSum_Eta;
   Float_t         Gen_vSum_Phi;
   Float_t         Gen_vSum_Rapidity;
   Float_t         Gen_EleMu_dphi;
   Float_t         Gen_EleMu_acop;
   Int_t           EleMuEvents;
   Float_t         Gen_Mu_Pt2;
   Float_t         Gen_Mu_Eta2;
   Float_t         Gen_Mu_Phi2;
   Float_t         Tau_Pt2;
   Float_t         Gen_vSum_M2;
   Float_t         Gen_vSum_Pt2;
   Float_t         Gen_Ele_Pt2;
   Float_t         Gen_Ele_Eta2;
   Float_t         Gen_Ele_Phi2;
   Int_t           run;
   Int_t           ls;
   Int_t           evtnb;
   Float_t         Ele_Pt;
   Float_t         Ele_Eta;
   Float_t         Ele_Phi;
   Float_t         Ele_E;
   Int_t           Ele_Charge;
   Float_t         Mu_Pt;
   Float_t         Mu_Eta;
   Float_t         Mu_Phi;
   Float_t         Mu_E;
   Int_t           Mu_Charge;
   Float_t         vSum_M;
   Float_t         vSum_Pt;
   Float_t         Scalar_Pt;
   Float_t         vSum_Eta;
   Float_t         vSum_Phi;
   Float_t         vSum_Rapidity;
   Float_t         EleMu_dphi;
   Float_t         EleMu_acop;
   Int_t           ok_neuexcl;
   Int_t           ok_chexcl;
   Int_t           ok_chexcl_tracks;
   Int_t           ok_chexcl_goodtracks;
   Int_t           ok_trigger;
   Float_t         dEta_Gen_Reco;
   Float_t         dPt_Gen_Reco;
   Float_t         dPhi_Gen_Reco;
   Float_t         DeltaR_ele;
   Float_t         DeltaR_mu;
   Int_t           RecoEleMu;
   Int_t           Check;

   // List of branches
   TBranch        *b_Tau_Pt1;   //!
   TBranch        *b_Tau_Pt;   //!
   TBranch        *b_Tau_Eta;   //!
   TBranch        *b_Tau_Phi;   //!
   TBranch        *b_Tau_E;   //!
   TBranch        *b_Gen_Ele_Pt;   //!
   TBranch        *b_Gen_Ele_Eta;   //!
   TBranch        *b_Gen_Ele_Phi;   //!
   TBranch        *b_Gen_Ele_E;   //!
   TBranch        *b_Gen_Mu_Pt;   //!
   TBranch        *b_Gen_Mu_Eta;   //!
   TBranch        *b_Gen_Mu_Phi;   //!
   TBranch        *b_Gen_Mu_E;   //!
   TBranch        *b_Gen_vSum_M;   //!
   TBranch        *b_Gen_vSum_Pt;   //!
   TBranch        *b_Gen_Scalar_Pt;   //!
   TBranch        *b_Gen_vSum_Eta;   //!
   TBranch        *b_Gen_vSum_Phi;   //!
   TBranch        *b_Gen_vSum_Rapidity;   //!
   TBranch        *b_Gen_EleMu_dphi;   //!
   TBranch        *b_Gen_EleMu_acop;   //!
   TBranch        *b_EleMuEvents;   //!
   TBranch        *b_Gen_Mu_Pt2;   //!
   TBranch        *b_Gen_Mu_Eta2;   //!
   TBranch        *b_Gen_Mu_Phi2;   //!
   TBranch        *b_Tau_Pt2;   //!
   TBranch        *b_Gen_vSum_M2;   //!
   TBranch        *b_Gen_vSum_Pt2;   //!
   TBranch        *b_Gen_Ele_Pt2;   //!
   TBranch        *b_Gen_Ele_Eta2;   //!
   TBranch        *b_Gen_Ele_Phi2;   //!
   TBranch        *b_run;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_evtnb;   //!
   TBranch        *b_Ele_Pt;   //!
   TBranch        *b_Ele_Eta;   //!
   TBranch        *b_Ele_Phi;   //!
   TBranch        *b_Ele_E;   //!
   TBranch        *b_Ele_Charge;   //!
   TBranch        *b_Mu_Pt;   //!
   TBranch        *b_Mu_Eta;   //!
   TBranch        *b_Mu_Phi;   //!
   TBranch        *b_Mu_E;   //!
   TBranch        *b_Mu_Charge;   //!
   TBranch        *b_vSum_M;   //!
   TBranch        *b_vSum_Pt;   //!
   TBranch        *b_Scalar_Pt;   //!
   TBranch        *b_vSum_Eta;   //!
   TBranch        *b_vSum_Phi;   //!
   TBranch        *b_vSum_Rapidity;   //!
   TBranch        *b_EleMu_dphi;   //!
   TBranch        *b_EleMu_acop;   //!
   TBranch        *b_ok_neuexcl;   //!
   TBranch        *b_ok_chexcl;   //!
   TBranch        *b_ok_chexcl_tracks;   //!
   TBranch        *b_ok_chexcl_goodtracks;   //!
   TBranch        *b_ok_trigger;   //!
   TBranch        *b_dEta_Gen_Reco;   //!
   TBranch        *b_dPt_Gen_Reco;   //!
   TBranch        *b_dPhi_Gen_Reco;   //!
   TBranch        *b_DeltaR_ele;   //!
   TBranch        *b_DeltaR_mu;   //!
   TBranch        *b_RecoEleMu;   //!
   TBranch        *b_Check;   //!

   ReadTauTauTree(TTree *tree=0);
   virtual ~ReadTauTauTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ReadTauTauTree_cxx
ReadTauTauTree::ReadTauTauTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootFile/mc_GammaUPC_EleMu_ArashNew_13thSep_moreEvents_merged.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootFile/mc_GammaUPC_EleMu_ArashNew_13thSep_moreEvents_merged.root");
      }
      f->GetObject("output_tree",tree);

   }
   Init(tree);
}

ReadTauTauTree::~ReadTauTauTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ReadTauTauTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ReadTauTauTree::LoadTree(Long64_t entry)
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

void ReadTauTauTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("Tau_Pt1", &Tau_Pt1, &b_Tau_Pt1);
   fChain->SetBranchAddress("Tau_Pt", &Tau_Pt, &b_Tau_Pt);
   fChain->SetBranchAddress("Tau_Eta", &Tau_Eta, &b_Tau_Eta);
   fChain->SetBranchAddress("Tau_Phi", &Tau_Phi, &b_Tau_Phi);
   fChain->SetBranchAddress("Tau_E", &Tau_E, &b_Tau_E);
   fChain->SetBranchAddress("Gen_Ele_Pt", &Gen_Ele_Pt, &b_Gen_Ele_Pt);
   fChain->SetBranchAddress("Gen_Ele_Eta", &Gen_Ele_Eta, &b_Gen_Ele_Eta);
   fChain->SetBranchAddress("Gen_Ele_Phi", &Gen_Ele_Phi, &b_Gen_Ele_Phi);
   fChain->SetBranchAddress("Gen_Ele_E", &Gen_Ele_E, &b_Gen_Ele_E);
   fChain->SetBranchAddress("Gen_Mu_Pt", &Gen_Mu_Pt, &b_Gen_Mu_Pt);
   fChain->SetBranchAddress("Gen_Mu_Eta", &Gen_Mu_Eta, &b_Gen_Mu_Eta);
   fChain->SetBranchAddress("Gen_Mu_Phi", &Gen_Mu_Phi, &b_Gen_Mu_Phi);
   fChain->SetBranchAddress("Gen_Mu_E", &Gen_Mu_E, &b_Gen_Mu_E);
   fChain->SetBranchAddress("Gen_vSum_M", &Gen_vSum_M, &b_Gen_vSum_M);
   fChain->SetBranchAddress("Gen_vSum_Pt", &Gen_vSum_Pt, &b_Gen_vSum_Pt);
   fChain->SetBranchAddress("Gen_Scalar_Pt", &Gen_Scalar_Pt, &b_Gen_Scalar_Pt);
   fChain->SetBranchAddress("Gen_vSum_Eta", &Gen_vSum_Eta, &b_Gen_vSum_Eta);
   fChain->SetBranchAddress("Gen_vSum_Phi", &Gen_vSum_Phi, &b_Gen_vSum_Phi);
   fChain->SetBranchAddress("Gen_vSum_Rapidity", &Gen_vSum_Rapidity, &b_Gen_vSum_Rapidity);
   fChain->SetBranchAddress("Gen_EleMu_dphi", &Gen_EleMu_dphi, &b_Gen_EleMu_dphi);
   fChain->SetBranchAddress("Gen_EleMu_acop", &Gen_EleMu_acop, &b_Gen_EleMu_acop);
   fChain->SetBranchAddress("EleMuEvents", &EleMuEvents, &b_EleMuEvents);
   fChain->SetBranchAddress("Gen_Mu_Pt2", &Gen_Mu_Pt2, &b_Gen_Mu_Pt2);
   fChain->SetBranchAddress("Gen_Mu_Eta2", &Gen_Mu_Eta2, &b_Gen_Mu_Eta2);
   fChain->SetBranchAddress("Gen_Mu_Phi2", &Gen_Mu_Phi2, &b_Gen_Mu_Phi2);
   fChain->SetBranchAddress("Tau_Pt2", &Tau_Pt2, &b_Tau_Pt2);
   fChain->SetBranchAddress("Gen_vSum_M2", &Gen_vSum_M2, &b_Gen_vSum_M2);
   fChain->SetBranchAddress("Gen_vSum_Pt2", &Gen_vSum_Pt2, &b_Gen_vSum_Pt2);
   fChain->SetBranchAddress("Gen_Ele_Pt2", &Gen_Ele_Pt2, &b_Gen_Ele_Pt2);
   fChain->SetBranchAddress("Gen_Ele_Eta2", &Gen_Ele_Eta2, &b_Gen_Ele_Eta2);
   fChain->SetBranchAddress("Gen_Ele_Phi2", &Gen_Ele_Phi2, &b_Gen_Ele_Phi2);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("evtnb", &evtnb, &b_evtnb);
   fChain->SetBranchAddress("Ele_Pt", &Ele_Pt, &b_Ele_Pt);
   fChain->SetBranchAddress("Ele_Eta", &Ele_Eta, &b_Ele_Eta);
   fChain->SetBranchAddress("Ele_Phi", &Ele_Phi, &b_Ele_Phi);
   fChain->SetBranchAddress("Ele_E", &Ele_E, &b_Ele_E);
   fChain->SetBranchAddress("Ele_Charge", &Ele_Charge, &b_Ele_Charge);
   fChain->SetBranchAddress("Mu_Pt", &Mu_Pt, &b_Mu_Pt);
   fChain->SetBranchAddress("Mu_Eta", &Mu_Eta, &b_Mu_Eta);
   fChain->SetBranchAddress("Mu_Phi", &Mu_Phi, &b_Mu_Phi);
   fChain->SetBranchAddress("Mu_E", &Mu_E, &b_Mu_E);
   fChain->SetBranchAddress("Mu_Charge", &Mu_Charge, &b_Mu_Charge);
   fChain->SetBranchAddress("vSum_M", &vSum_M, &b_vSum_M);
   fChain->SetBranchAddress("vSum_Pt", &vSum_Pt, &b_vSum_Pt);
   fChain->SetBranchAddress("Scalar_Pt", &Scalar_Pt, &b_Scalar_Pt);
   fChain->SetBranchAddress("vSum_Eta", &vSum_Eta, &b_vSum_Eta);
   fChain->SetBranchAddress("vSum_Phi", &vSum_Phi, &b_vSum_Phi);
   fChain->SetBranchAddress("vSum_Rapidity", &vSum_Rapidity, &b_vSum_Rapidity);
   fChain->SetBranchAddress("EleMu_dphi", &EleMu_dphi, &b_EleMu_dphi);
   fChain->SetBranchAddress("EleMu_acop", &EleMu_acop, &b_EleMu_acop);
   fChain->SetBranchAddress("ok_neuexcl", &ok_neuexcl, &b_ok_neuexcl);
   fChain->SetBranchAddress("ok_chexcl", &ok_chexcl, &b_ok_chexcl);
   fChain->SetBranchAddress("ok_chexcl_tracks", &ok_chexcl_tracks, &b_ok_chexcl_tracks);
   fChain->SetBranchAddress("ok_chexcl_goodtracks", &ok_chexcl_goodtracks, &b_ok_chexcl_goodtracks);
   fChain->SetBranchAddress("ok_trigger", &ok_trigger, &b_ok_trigger);
   fChain->SetBranchAddress("dEta_Gen_Reco", &dEta_Gen_Reco, &b_dEta_Gen_Reco);
   fChain->SetBranchAddress("dPt_Gen_Reco", &dPt_Gen_Reco, &b_dPt_Gen_Reco);
   fChain->SetBranchAddress("dPhi_Gen_Reco", &dPhi_Gen_Reco, &b_dPhi_Gen_Reco);
   fChain->SetBranchAddress("DeltaR_ele", &DeltaR_ele, &b_DeltaR_ele);
   fChain->SetBranchAddress("DeltaR_mu", &DeltaR_mu, &b_DeltaR_mu);
   fChain->SetBranchAddress("RecoEleMu", &RecoEleMu, &b_RecoEleMu);
   fChain->SetBranchAddress("Check", &Check, &b_Check);
   Notify();
}

Bool_t ReadTauTauTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ReadTauTauTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ReadTauTauTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ReadTauTauTree_cxx
