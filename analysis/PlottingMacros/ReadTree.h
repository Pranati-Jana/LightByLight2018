//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 24 12:55:55 2025 by ROOT version 6.28/06
// from TTree output_tree/
// found on file: RootFile_From_selectEleMuEvents/NewRunningHFThreshold90Percent/data_obs_EleMu.root
//////////////////////////////////////////////////////////

#ifndef ReadTree_h
#define ReadTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ReadTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           ls;
   Int_t           evtnb;
   Float_t         Ele_Pt;
   Float_t         Ele_Eta;
   Float_t         Ele_Phi;
   Float_t         Ele_E;
   Int_t           Ele_charge;
   Float_t         Mu_Pt;
   Float_t         Mu_Eta;
   Float_t         Mu_Phi;
   Float_t         Mu_E;
   Int_t           Mu_charge;
   Float_t         vSum_M;
   Float_t         vSum_Pt;
   Float_t         Scalar_Pt;
   Float_t         vSum_Eta;
   Float_t         vSum_Phi;
   Float_t         vSum_Rapidity;
   Float_t         EleMu_dphi;
   Float_t         EleMu_acop;
   Int_t           ok_neuexcl;
   Int_t           ok_neuexcl_HFonly;
   Int_t           ok_chexcl;
   Int_t           ok_chexcl_extratracks;
   Int_t           ok_chexcl_goodtracks;
   Int_t           ok_trigger;
   Int_t           ok_singleMuon;
   Int_t           ok_singleEG3;
   Int_t           ok_singleEG5;
   Float_t         deltaRtrack_electron;
   Float_t         deltaRtrack_muon;
   Float_t         zdc_energy_pos;
   Float_t         zdc_energy_neg;
   Int_t           nTracks;
   Float_t         track_Pt;
   Float_t         track_Eta;
   Float_t         track_Phi;
   Int_t           track_Charge;
   Int_t           nExtrk;
   Float_t         vSum_Pz;
   Float_t         tower_eta;
   Float_t         tower_phi;
   Float_t         tower_em;
   Float_t         tower_had;
   Float_t         tower_deltaRmuon_HE;
   Float_t         tower_deltaRmuon_HB;
   Float_t         tracks_deltaRElectron;
   Float_t         tower_HF_Plus;
   Float_t         tower_HF_Minus;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_evtnb;   //!
   TBranch        *b_Ele_Pt;   //!
   TBranch        *b_Ele_Eta;   //!
   TBranch        *b_Ele_Phi;   //!
   TBranch        *b_Ele_E;   //!
   TBranch        *b_Ele_charge;   //!
   TBranch        *b_Mu_Pt;   //!
   TBranch        *b_Mu_Eta;   //!
   TBranch        *b_Mu_Phi;   //!
   TBranch        *b_Mu_E;   //!
   TBranch        *b_Mu_charge;   //!
   TBranch        *b_vSum_M;   //!
   TBranch        *b_vSum_Pt;   //!
   TBranch        *b_Scalar_Pt;   //!
   TBranch        *b_vSum_Eta;   //!
   TBranch        *b_vSum_Phi;   //!
   TBranch        *b_vSum_Rapidity;   //!
   TBranch        *b_EleMu_dphi;   //!
   TBranch        *b_EleMu_acop;   //!
   TBranch        *b_ok_neuexcl;   //!
   TBranch        *b_ok_neuexcl_HFonly;   //!
   TBranch        *b_ok_chexcl;   //!
   TBranch        *b_ok_chexcl_extratracks;   //!
   TBranch        *b_ok_chexcl_goodtracks;   //!
   TBranch        *b_ok_trigger;   //!
   TBranch        *b_ok_singleMuon;   //!
   TBranch        *b_ok_singleEG3;   //!
   TBranch        *b_ok_singleEG5;   //!
   TBranch        *b_deltaRtrack_electron;   //!
   TBranch        *b_deltaRtrack_muon;   //!
   TBranch        *b_zdc_energy_pos;   //!
   TBranch        *b_zdc_energy_neg;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_track_Pt;   //!
   TBranch        *b_track_Eta;   //!
   TBranch        *b_track_Phi;   //!
   TBranch        *b_track_Charge;   //!
   TBranch        *b_nExtrk;   //!
   TBranch        *b_vSum_Pz;   //!
   TBranch        *b_tower_eta;   //!
   TBranch        *b_tower_phi;   //!
   TBranch        *b_tower_em;   //!
   TBranch        *b_tower_had;   //!
   TBranch        *b_tower_deltaRmuon_HE;   //!
   TBranch        *b_tower_deltaRmuon_HB;   //!
   TBranch        *b_tracks_deltaRElectron;   //!
   TBranch        *b_tower_HF_Plus;   //!
   TBranch        *b_tower_HF_Minus;   //!

   ReadTree(TTree *tree=0);
   virtual ~ReadTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ReadTree_cxx
ReadTree::ReadTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootFile_From_selectEleMuEvents/NewRunningHFThreshold90Percent/data_obs_EleMu.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootFile_From_selectEleMuEvents/NewRunningHFThreshold90Percent/data_obs_EleMu.root");
      }
      f->GetObject("output_tree",tree);

   }
   Init(tree);
}

ReadTree::~ReadTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ReadTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ReadTree::LoadTree(Long64_t entry)
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

void ReadTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("evtnb", &evtnb, &b_evtnb);
   fChain->SetBranchAddress("Ele_Pt", &Ele_Pt, &b_Ele_Pt);
   fChain->SetBranchAddress("Ele_Eta", &Ele_Eta, &b_Ele_Eta);
   fChain->SetBranchAddress("Ele_Phi", &Ele_Phi, &b_Ele_Phi);
   fChain->SetBranchAddress("Ele_E", &Ele_E, &b_Ele_E);
   fChain->SetBranchAddress("Ele_charge", &Ele_charge, &b_Ele_charge);
   fChain->SetBranchAddress("Mu_Pt", &Mu_Pt, &b_Mu_Pt);
   fChain->SetBranchAddress("Mu_Eta", &Mu_Eta, &b_Mu_Eta);
   fChain->SetBranchAddress("Mu_Phi", &Mu_Phi, &b_Mu_Phi);
   fChain->SetBranchAddress("Mu_E", &Mu_E, &b_Mu_E);
   fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
   fChain->SetBranchAddress("vSum_M", &vSum_M, &b_vSum_M);
   fChain->SetBranchAddress("vSum_Pt", &vSum_Pt, &b_vSum_Pt);
   fChain->SetBranchAddress("Scalar_Pt", &Scalar_Pt, &b_Scalar_Pt);
   fChain->SetBranchAddress("vSum_Eta", &vSum_Eta, &b_vSum_Eta);
   fChain->SetBranchAddress("vSum_Phi", &vSum_Phi, &b_vSum_Phi);
   fChain->SetBranchAddress("vSum_Rapidity", &vSum_Rapidity, &b_vSum_Rapidity);
   fChain->SetBranchAddress("EleMu_dphi", &EleMu_dphi, &b_EleMu_dphi);
   fChain->SetBranchAddress("EleMu_acop", &EleMu_acop, &b_EleMu_acop);
   fChain->SetBranchAddress("ok_neuexcl", &ok_neuexcl, &b_ok_neuexcl);
   fChain->SetBranchAddress("ok_neuexcl_HFonly", &ok_neuexcl_HFonly, &b_ok_neuexcl_HFonly);
   fChain->SetBranchAddress("ok_chexcl", &ok_chexcl, &b_ok_chexcl);
   fChain->SetBranchAddress("ok_chexcl_extratracks", &ok_chexcl_extratracks, &b_ok_chexcl_extratracks);
   fChain->SetBranchAddress("ok_chexcl_goodtracks", &ok_chexcl_goodtracks, &b_ok_chexcl_goodtracks);
   fChain->SetBranchAddress("ok_trigger", &ok_trigger, &b_ok_trigger);
   fChain->SetBranchAddress("ok_singleMuon", &ok_singleMuon, &b_ok_singleMuon);
   fChain->SetBranchAddress("ok_singleEG3", &ok_singleEG3, &b_ok_singleEG3);
   fChain->SetBranchAddress("ok_singleEG5", &ok_singleEG5, &b_ok_singleEG5);
   fChain->SetBranchAddress("deltaRtrack_electron", &deltaRtrack_electron, &b_deltaRtrack_electron);
   fChain->SetBranchAddress("deltaRtrack_muon", &deltaRtrack_muon, &b_deltaRtrack_muon);
   fChain->SetBranchAddress("zdc_energy_pos", &zdc_energy_pos, &b_zdc_energy_pos);
   fChain->SetBranchAddress("zdc_energy_neg", &zdc_energy_neg, &b_zdc_energy_neg);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("track_Pt", &track_Pt, &b_track_Pt);
   fChain->SetBranchAddress("track_Eta", &track_Eta, &b_track_Eta);
   fChain->SetBranchAddress("track_Phi", &track_Phi, &b_track_Phi);
   fChain->SetBranchAddress("track_Charge", &track_Charge, &b_track_Charge);
   fChain->SetBranchAddress("nExtrk", &nExtrk, &b_nExtrk);
   fChain->SetBranchAddress("vSum_Pz", &vSum_Pz, &b_vSum_Pz);
   fChain->SetBranchAddress("tower_eta", &tower_eta, &b_tower_eta);
   fChain->SetBranchAddress("tower_phi", &tower_phi, &b_tower_phi);
   fChain->SetBranchAddress("tower_em", &tower_em, &b_tower_em);
   fChain->SetBranchAddress("tower_had", &tower_had, &b_tower_had);
   fChain->SetBranchAddress("tower_deltaRmuon_HE", &tower_deltaRmuon_HE, &b_tower_deltaRmuon_HE);
   fChain->SetBranchAddress("tower_deltaRmuon_HB", &tower_deltaRmuon_HB, &b_tower_deltaRmuon_HB);
   fChain->SetBranchAddress("tracks_deltaRElectron", &tracks_deltaRElectron, &b_tracks_deltaRElectron);
   fChain->SetBranchAddress("tower_HF_Plus", &tower_HF_Plus, &b_tower_HF_Plus);
   fChain->SetBranchAddress("tower_HF_Minus", &tower_HF_Minus, &b_tower_HF_Minus);
   Notify();
}

Bool_t ReadTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ReadTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ReadTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ReadTree_cxx
