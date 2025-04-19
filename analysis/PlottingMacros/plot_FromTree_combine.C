#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include "ReadTree.C"
#include "CMS_lumi.C"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TCut.h"
#include "TChain.h"
#include "THStack.h"
#include "TString.h"
#include "tnp_muon_UPC_PbPb2018.h"

int xlo = 1;
int xhi = 2;
int nbin = 8;


void make_canvas(TCanvas *&);
void make_hist(TH1D *&, Color_t , int , int);
void make_canvas_ratio(TCanvas *&);
void make_hist_ratio(TH1D *&, Color_t, int );
//TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D* , TH1D*, TH1D*, TH1D*, TH1D*, double , double , double , double , double , double , const char *, bool); 

TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, TH1D*, double , double, double , double , double, double, const char *, bool); 

const int nSample = 24;
//const char *Sample[nSample]={"data_obs","GAMMA_UPC","MUMU_FSR","gamgamToTauTau0", "gamgamToTauTau0.01", "gamgamToTauTau0.02","gamgamToTauTau0.03", "gamgamToTauTau0.04", "gamgamToTauTau0.05", "gamgamToTauTau0.06", "gamgamToTauTau0.07", "gamgamToTauTau0.08", "gamgamToTauTau0.09", "gamgamToTauTau0.1", "gamgamToTauTau-0.01", "gamgamToTauTau-0.02", "gamgamToTauTau-0.03", "gamgamToTauTau-0.04", "gamgamToTauTau-0.05", "gamgamToTauTau-0.06", "gamgamToTauTau-0.07", "gamgamToTauTau-0.08", "gamgamToTauTau-0.09", "gamgamToTauTau-0.1"};
const char *Sample[nSample]={"data_obs","GAMMA_UPC","MUMU_FSR","gamgamToTauTau0", "gamgamToTauTau0.01", "gamgamToTauTau0.02","gamgamToTauTau0.03", "gamgamToTauTau0.04", "gamgamToTauTau0.05", "gamgamToTauTau0.06", "gamgamToTauTau0.07", "gamgamToTauTau0.08", "gamgamToTauTau0.09", "gamgamToTauTau0.1", "gamgamToTauTau-0.01", "gamgamToTauTau-0.02", "gamgamToTauTau-0.03", "gamgamToTauTau-0.04", "gamgamToTauTau-0.05", "gamgamToTauTau-0.06", "gamgamToTauTau-0.07", "gamgamToTauTau-0.08", "gamgamToTauTau-0.09", "gamgamToTauTau-0.1"};
//const char *Sample[nSample]={"Data"};

const char *dir = "figures_eta2p2";

//double luminosity = 425.91;/
//double luminosity = 1647.228136;
double luminosity = 1700.0;
//double crossSectionGammaUPC = 1060.8959;//GammaUPC
double crossSectionGammaUPC =  1;
double generatedEvtGammaUPC = 5290000;//GammaUPC
double ScaleFactor_new = 0.943* //0.952//electron reco+ID //ele pt 2.5GeV
                       0.96*
                       0.9*  //HF
                       0.99;  // neutral exclusivity
                      //charge exclusivity

double Other_gUPC = (luminosity * crossSectionGammaUPC * ScaleFactor_new)/generatedEvtGammaUPC;//gUPC
//double Other_gUPC = 1;

//double crossSectionFSR = 139.7;//MC_FSR_Arash
double crossSectionFSR = 1;
double generatedEvtFSR = 2471810;//FSR

double Other_mumuFSR = (luminosity * crossSectionFSR*ScaleFactor_new)/generatedEvtFSR;//MC_Arash_FSR

//double crossSectionUPCGen = 1000;
double crossSectionUPCGen = 1;
double generatedEvtUPCGen = 10000000;
double Other_UPCGen = (luminosity * crossSectionUPCGen*ScaleFactor_new)/generatedEvtUPCGen;//


void BinLogX(TH1* h);

 const int nPtbins=10;
 //double Ptbin[nPtbins]={1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5};
 //double Ptbin[nPtbins]={0.0,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,20.5};
 double Ptbin[nPtbins]={1.5,3.5,5.5,7.5,9.5,11.5,13.5,15.5,17.5,20.5};
 const int nPtbin= sizeof(Ptbin)/sizeof(double) - 1;
void plot_FromTree_combine(){
  cout << "other_UPCGen:" << Other_UPCGen << endl;
  //const double wt[nSample] = {1}; 
  double ET;
  double PT;
//  int idx = -1;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  string mc_name = "SC";
  TFile *outf= new TFile("RootFile_For_combine/ElectronMuon_6.4GeVHFcut_ChFF_NormalizedTo1mbCrossSection_2GeVbinning.root","recreate");  
  TChain *tree[nSample];

  //// histograms
  TH1D* hmu_pt[nSample], *hmu_pt_1sigmaUp_sys[nSample], *hmu_pt_1sigmaDown_sys[nSample], *hmu_pt_1sigmaUp_stat[nSample], *hmu_pt_1sigmaDown_stat[nSample], *hele_pt[nSample], *hele_pt_1sigmaUp_sys[nSample], *hele_pt_1sigmaDown_sys[nSample], *hele_pt_1sigmaUp_stat[nSample], *hele_pt_1sigmaDown_stat[nSample];
  TH1D* hmu_pt_1sigmaUp_sys_stat[nSample], *hmu_pt_1sigmaDown_sys_stat[nSample], *hele_pt_1sigmaUp_sys_stat[nSample], *hele_pt_1sigmaDown_sys_stat[nSample];
  for (int i = 0; i < nSample; i++){
//  for (int i = 0; i < 4; i++){
    cout << "nSample:" << i  << endl;
   tree[i] = new TChain("output_tree");
//   tree[i]->Add(Form("/home/pranati/Documents/GammaGammaToTauTau/Workspace/plottingMacros/RootFile_From_selectEleMuEvents/11thDecember_From_selectEleMuEvents_Tree/%s_EleMu.root", Sample[i]));   
   //tree[i]->Add(Form("/Users/pranatijana/Documents/GammaGamma_to_Tautau/Workspace/plottingMacros/RootFile_For_combine/a_Tau_points/%s_EleMu.root",Sample[i]));  
 //tree[i]->Add(Form("/Users/pranatijana/Documents/GammaGamma_to_Tautau/Workspace/plottingMacros/RootFile_For_combine/Test/%s_EleMu.root",Sample[i]));  
   //tree[i]->Add(Form("/Users/pranatijana/Documents/GammaGamma_to_Tautau/Workspace/plottingMacros/RootFile_From_selectEleMuEvents/ExtraTrack_ExclusivityAcopCut_NotApplied/21aTau_Points_UPCGen_ChargedExcElectronFromLByLMuonDeltaR0.001_CaloTowerInfo/%s_EleMu.root",Sample[i]));
   //tree[i]->Add(Form("/Users/pranatijana/Documents/GammaGamma_to_Tautau/Workspace/plottingMacros/RootFile_From_selectEleMuEvents/loose_ElectronMatchingToTower_deltaRMatchingTowerToMuon/aTau_DeltaR1HE3/%s_EleMu.root",Sample[i]));
   //tree[i]->Add(Form("RootFile_From_selectEleMuEvents/MuonEndcapPt1.5GeV_HE1To3DeltaRMuonHE1/%s_EleMu.root", Sample[i]));
   tree[i]->Add(Form("RootFile_From_selectEleMuEvents/a_Tau_Points_6.4GeVHFcut/%s_EleMu.root", Sample[i]));
   hmu_pt[i] = new TH1D(Form("%s", Sample[i]),"",nPtbin, Ptbin);

   if(i==2){
   hmu_pt_1sigmaUp_sys[i] = new TH1D(Form("%s_mumu_fsr_bkg_sys_SFUp", Sample[i]),"",nPtbin, Ptbin);
   hmu_pt_1sigmaDown_sys[i] = new TH1D(Form("%s_mumu_fsr_bkg_sys_SFDown", Sample[i]),"",nPtbin, Ptbin);
   hmu_pt_1sigmaUp_stat[i] = new TH1D(Form("%s_mumu_fsr_bkg_stat_SFUp", Sample[i]),"",nPtbin, Ptbin);
   hmu_pt_1sigmaDown_stat[i] = new TH1D(Form("%s_mumu_fsr_bkg_stat_SFDown", Sample[i]),"",nPtbin,Ptbin);
   hmu_pt_1sigmaUp_sys_stat[i] = new TH1D(Form("%s_mumu_fsr_bkg_sys_stat_SFUp", Sample[i]),"",nPtbin,Ptbin);
   hmu_pt_1sigmaDown_sys_stat[i] = new TH1D(Form("%s_mumu_fsr_bkg_sys_stat_SFDown", Sample[i]),"",nPtbin,Ptbin);
   }
  else{
   hmu_pt_1sigmaUp_sys[i] = new TH1D(Form("%s_sys_SFUp", Sample[i]),"",nPtbin, Ptbin);
   hmu_pt_1sigmaDown_sys[i] = new TH1D(Form("%s_sys_SFDown", Sample[i]),"",nPtbin, Ptbin);
   hmu_pt_1sigmaUp_stat[i] = new TH1D(Form("%s_stat_SFUp", Sample[i]),"",nPtbin, Ptbin);
   hmu_pt_1sigmaDown_stat[i] = new TH1D(Form("%s_stat_SFDown", Sample[i]),"",nPtbin,Ptbin);
   hmu_pt_1sigmaUp_sys_stat[i] = new TH1D(Form("%s_sys_stat_SFUp", Sample[i]),"",nPtbin,Ptbin);
   hmu_pt_1sigmaDown_sys_stat[i] = new TH1D(Form("%s_sys_stat_SFDown", Sample[i]),"",nPtbin,Ptbin);
  }
  
  hele_pt[i] = new TH1D(Form("%s", Sample[i]),"",nPtbin, Ptbin);
   if(i==2){
   hele_pt_1sigmaUp_sys[i]= new TH1D(Form("%s_mumu_fsr_bkg_sys_SFUp", Sample[i]),"",nPtbin, Ptbin);
   hele_pt_1sigmaDown_sys[i]= new TH1D(Form("%s_mumu_fsr_bkg_sys_SFDown", Sample[i]),"",nPtbin, Ptbin);
   hele_pt_1sigmaUp_stat[i]= new TH1D(Form("%s_mumu_fsr_bkg_stat_SFUp", Sample[i]),"",nPtbin, Ptbin);
   hele_pt_1sigmaDown_stat[i]= new TH1D(Form("%s_mumu_fsr_bkg_stat_SFDown", Sample[i]),"",nPtbin, Ptbin);
   hele_pt_1sigmaUp_sys_stat[i] = new TH1D(Form("%s_mumu_fsr_bkg_sys_stat_SFUp", Sample[i]),"",nPtbin,Ptbin);
   hele_pt_1sigmaDown_sys_stat[i] = new TH1D(Form("%s_mumu_fsr_bkg_sys_stat_SFDown", Sample[i]),"",nPtbin,Ptbin);

   }
   else{
   hele_pt_1sigmaUp_sys[i]= new TH1D(Form("%s_sys_SFUp", Sample[i]),"",nPtbin, Ptbin);
   hele_pt_1sigmaDown_sys[i]= new TH1D(Form("%s_sys_SFDown", Sample[i]),"",nPtbin, Ptbin);
   hele_pt_1sigmaUp_stat[i]= new TH1D(Form("%s_stat_SFUp", Sample[i]),"",nPtbin, Ptbin);
   hele_pt_1sigmaDown_stat[i]= new TH1D(Form("%s_stat_SFDown", Sample[i]),"",nPtbin, Ptbin);
   hele_pt_1sigmaUp_sys_stat[i] = new TH1D(Form("%s_sys_stat_SFUp", Sample[i]),"",nPtbin,Ptbin);
   hele_pt_1sigmaDown_sys_stat[i] = new TH1D(Form("%s_sys_stat_SFDown", Sample[i]),"",nPtbin,Ptbin);
   }

   //  hmu_pt[i] = new TH1D(Form("%s", Sample[i]),"",23,2.5,25.5);
  //  hmu_pt_1sigmaUp_sys[i] = new TH1D(Form("%s_sys_SFUp", Sample[i]),"",23,2.5,25.5);
  //  hmu_pt_1sigmaDown_sys[i] = new TH1D(Form("%s_sys_SFDown", Sample[i]),"",23,2.5,25.5);
  //  hmu_pt_1sigmaUp_stat[i] = new TH1D(Form("%s_stat_SFUp", Sample[i]),"",23,2.5,25.5);
  //  hmu_pt_1sigmaDown_stat[i] = new TH1D(Form("%s_stat_SFDown", Sample[i]),"",23,2.5,25.5);

  //  hele_pt[i] = new TH1D(Form("%s", Sample[i]),"",23,2.5,25.5);
  //  if(i==2){
  //  hele_pt_1sigmaUp_sys[i]= new TH1D(Form("%s_mumu_fsr_bkg_sys_SFUp", Sample[i]),"",23,2.5,25.5);
  //  hele_pt_1sigmaDown_sys[i]= new TH1D(Form("%s_mumu_fsr_bkg_sys_SFDown", Sample[i]),"",23,2.5,25.5);
  //  hele_pt_1sigmaUp_stat[i]= new TH1D(Form("%s_mumu_fsr_bkg_stat_SFUp", Sample[i]),"",23,2.5,25.5);
  //  hele_pt_1sigmaDown_stat[i]= new TH1D(Form("%s_mumu_fsr_bkg_stat_SFDown", Sample[i]),"",23,2.5,25.5);
  //  }
  //  else{
  //   hele_pt_1sigmaUp_sys[i]= new TH1D(Form("%s_sys_SFUp", Sample[i]),"",23,2.5,25.5);
  //  hele_pt_1sigmaDown_sys[i]= new TH1D(Form("%s_sys_SFDown", Sample[i]),"",23,2.5,25.5);
  //  hele_pt_1sigmaUp_stat[i]= new TH1D(Form("%s_stat_SFUp", Sample[i]),"",23,2.5,25.5);
  //  hele_pt_1sigmaDown_stat[i]= new TH1D(Form("%s_stat_SFDown", Sample[i]),"",23,2.5,25.5);
  //  }
  

   cout << "nSample:" << i  << endl;
   cout << "file " << Sample[i] << ":" << tree[i]->GetEntries()  << endl;
   ReadTree  treeR(tree[i]);
   treeR.fChain->SetBranchStatus("*",1);
   if (treeR.fChain == 0) return;
   Long64_t nentries = treeR.fChain->GetEntriesFast();
   cout << tree[i]->GetName() << "    " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
  
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<30;jentry++) { 
      Long64_t ientry_evt = treeR.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = treeR.fChain->GetEntry(jentry);   nbytes += nb;


       if(treeR.Ele_Pt < 2.5) continue;
     if(treeR.ok_neuexcl!=1) continue;
    // if(i==0){
     if(treeR.ok_neuexcl_HFonly!=1) continue;
     if(treeR.tower_HF_Minus > 6.4) continue;
     if(treeR.tower_HF_Plus > 6.4) continue;
     //}
     if(treeR.ok_chexcl_goodtracks!=1) continue;
     //if(treeR.ok_chexcl_extratracks!=1) continue;
     //if(treeR.nTracks >4) continue;
     if(treeR.EleMu_acop <= 0.01) continue;
     if(treeR.EleMu_acop > 0.15) continue;
      //SF
      ET = treeR.Mu_Eta;
      PT = treeR.Mu_Pt;

     

      //std::vector<double> CrossSection = {1,1060.8959,139.7,466.087,478.525,491.875,506.371,521.758,538.059,555.433,573.662,593.409,614.035,635.410,454.508,443.656,433.879,424.608,416.405,408.804,402.353,396.321,391.171,386.746};
    // std::vector<double> CrossSection = {1,1060.8959,139.7,466.087,478.525,491.875,506.371,521.758,538.059,555.433,573.662,593.409,614.035,635.410,454.508,443.656,433.879,424.608,416.405,408.804,402.353,396.321,391.171,386.746};
     //std::vector<double> CrossSection = {1,1060.8959,139.7,564.89,579.9,596.1,613.8,632.66,652.74,673.61,696.4,720.29,745.5,772.4,550.56,537.5,525.6,514.7,504.8,496.0,488.37,481.4,475.4,470.49};
      std::vector<double> CrossSection = {1,1000,139.7,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000};
     // std::vector<double> CrossSection = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      double norm_gUPC = tnp_weight_softid_upc_pbpb( PT, ET, 0) * tnp_weight_trigger_upc_pbpb( PT, ET, 0) * Other_gUPC;
      double norm_mumuFSR = tnp_weight_softid_upc_pbpb( PT, ET, 0) * tnp_weight_trigger_upc_pbpb( PT, ET, 0) * Other_mumuFSR;
      double norm_UPCGen = tnp_weight_softid_upc_pbpb( PT, ET, 0) * tnp_weight_trigger_upc_pbpb( PT, ET, 0) * Other_UPCGen;
     std::vector<double> wt = {1, norm_gUPC, norm_mumuFSR, norm_UPCGen, norm_UPCGen, norm_UPCGen,norm_UPCGen*1.19,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen,norm_UPCGen, norm_UPCGen,norm_UPCGen};
    //std::vector<double> wt = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

      double norm_gUPC_Sys_Up = tnp_weight_softid_upc_pbpb( PT, ET, -1) * tnp_weight_trigger_upc_pbpb( PT, ET, -1) * Other_gUPC;
      double norm_mumuFSR_Sys_Up = tnp_weight_softid_upc_pbpb( PT, ET, -1) * tnp_weight_trigger_upc_pbpb( PT, ET, -1) * Other_mumuFSR;
      double norm_UPCGen_Sys_Up = tnp_weight_softid_upc_pbpb( PT, ET, -1) * tnp_weight_trigger_upc_pbpb( PT, ET, -1) * Other_UPCGen;
 
      std::vector<double> wtSysUp = {1, norm_gUPC_Sys_Up, norm_mumuFSR_Sys_Up, norm_UPCGen_Sys_Up, norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up,norm_UPCGen_Sys_Up, norm_UPCGen_Sys_Up};
     //std::vector<double> wtSysUp = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      

      double norm_gUPC_Sys_Down = tnp_weight_softid_upc_pbpb( PT, ET, -2) * tnp_weight_trigger_upc_pbpb( PT, ET, -2) * Other_gUPC;
      double norm_mumuFSR_Sys_Down = tnp_weight_softid_upc_pbpb( PT, ET, -2) * tnp_weight_trigger_upc_pbpb( PT, ET, -2) * Other_mumuFSR;
      double norm_UPCGen_Sys_Down = tnp_weight_softid_upc_pbpb( PT, ET, -2) * tnp_weight_trigger_upc_pbpb( PT, ET, -2)*Other_UPCGen;
      std::vector<double> wtSysDown = {1, norm_gUPC_Sys_Down, norm_mumuFSR_Sys_Down, norm_UPCGen_Sys_Down, norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down,norm_UPCGen_Sys_Down, norm_UPCGen_Sys_Down};
      // std::vector<double> wtSysDown = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

      double norm_gUPC_Stat_Up = tnp_weight_softid_upc_pbpb( PT, ET, +1) * tnp_weight_trigger_upc_pbpb( PT, ET, +1) * Other_gUPC;
      double norm_mumuFSR_Stat_Up = tnp_weight_softid_upc_pbpb( PT, ET, +1) * tnp_weight_trigger_upc_pbpb( PT, ET, +1) * Other_mumuFSR;
      double norm_UPCGen_Stat_Up = tnp_weight_softid_upc_pbpb( PT, ET, +1) * tnp_weight_trigger_upc_pbpb( PT, ET, +1) *Other_UPCGen;
      std::vector<double> wtStatUp = {1, norm_gUPC_Stat_Up, norm_mumuFSR_Stat_Up, norm_UPCGen_Stat_Up, norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up,norm_UPCGen_Stat_Up, norm_UPCGen_Stat_Up};
      //std::vector<double> wtStatUp = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

      double norm_gUPC_Stat_Down = tnp_weight_softid_upc_pbpb( PT, ET, +2) * tnp_weight_trigger_upc_pbpb( PT, ET, +2) * Other_gUPC;
      double norm_mumuFSR_Stat_Down = tnp_weight_softid_upc_pbpb( PT, ET, +2) * tnp_weight_trigger_upc_pbpb( PT, ET, +2) * Other_mumuFSR;
      double norm_UPCGen_Stat_Down = tnp_weight_softid_upc_pbpb( PT, ET, +2) * tnp_weight_trigger_upc_pbpb( PT, ET, +2) * Other_UPCGen;
      std::vector<double> wtStatDown = {1, norm_gUPC_Stat_Down, norm_mumuFSR_Stat_Down, norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down,norm_UPCGen_Stat_Down};
     // std::vector<double> wtStatDown = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

     ///Stat and sys combined////
      /*
      double norm_gUPC_Sys_Stat_Up = sqrt(norm_gUPC_Sys_Up*norm_gUPC_Sys_Up +  norm_gUPC_Stat_Up*norm_gUPC_Stat_Up );
      double norm_gUPC_Sys_Stat_Down = sqrt(norm_gUPC_Sys_Down*norm_gUPC_Sys_Down +  norm_gUPC_Stat_Down*norm_gUPC_Stat_Down);
     
      double norm_UPCGen_Sys_Stat_Up = sqrt(norm_UPCGen_Sys_Up*norm_UPCGen_Sys_Up + norm_UPCGen_Stat_Up*norm_UPCGen_Stat_Up);
      double norm_UPCGen_Sys_Stat_Down = sqrt(norm_UPCGen_Sys_Down*norm_UPCGen_Sys_Down + norm_UPCGen_Stat_Down*norm_UPCGen_Stat_Down);
      
      double norm_mumuFSR_Sys_Stat_Up = sqrt(norm_mumuFSR_Sys_Up*norm_mumuFSR_Sys_Up + norm_mumuFSR_Stat_Up*norm_mumuFSR_Stat_Up);
      double norm_mumuFSR_Sys_Stat_Down = sqrt(norm_mumuFSR_Sys_Down*norm_mumuFSR_Sys_Down + norm_mumuFSR_Stat_Down*norm_mumuFSR_Stat_Down);
      
       std::vector<double> wtSysStatUp = {1, norm_gUPC_Sys_Stat_Up, norm_mumuFSR_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up*1.19, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up};   
       std::vector<double> wtSysStatDown = {1, norm_gUPC_Sys_Stat_Down, norm_mumuFSR_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down*1.19, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down};   
      */

        float muon_weight_sys_trigger = tnp_weight_softid_upc_pbpb( PT, ET, 0)*max(abs(tnp_weight_trigger_upc_pbpb(PT, ET,-1)-tnp_weight_trigger_upc_pbpb( PT, ET, 0)),abs(tnp_weight_trigger_upc_pbpb(PT, ET,-2)-tnp_weight_trigger_upc_pbpb( PT, ET, 0)));  
        float muon_weight_sys_soft = tnp_weight_trigger_upc_pbpb( PT, ET, 0)*max(abs(tnp_weight_softid_upc_pbpb(PT, ET,-1)-tnp_weight_softid_upc_pbpb( PT, ET, 0)),abs(tnp_weight_softid_upc_pbpb(PT, ET,-2)-tnp_weight_softid_upc_pbpb( PT, ET, 0)));
       
        float muon_weight_stat_trigger = tnp_weight_softid_upc_pbpb( PT, ET, 0)*max(abs(tnp_weight_trigger_upc_pbpb(PT, ET,+1)-tnp_weight_trigger_upc_pbpb( PT, ET, 0)),abs(tnp_weight_trigger_upc_pbpb(PT, ET,+2)-tnp_weight_trigger_upc_pbpb( PT, ET, 0)));  
        float muon_weight_stat_soft = tnp_weight_trigger_upc_pbpb( PT, ET, 0)*max(abs(tnp_weight_softid_upc_pbpb(PT, ET,+1)-tnp_weight_softid_upc_pbpb( PT, ET, 0)),abs(tnp_weight_softid_upc_pbpb(PT, ET,+2)-tnp_weight_softid_upc_pbpb( PT, ET, 0)));
        // #### Adding up all Muon SF uncertainties ####
        float  muon_weight_error = sqrt(pow(muon_weight_stat_trigger,2)+pow(muon_weight_stat_soft,2)+pow(muon_weight_sys_trigger,2)+pow(muon_weight_sys_soft,2));
      
        float muon_weight = tnp_weight_softid_upc_pbpb( PT, ET, 0) * tnp_weight_trigger_upc_pbpb( PT, ET, 0);

        float norm_gUPC_Sys_Stat_Up = Other_gUPC * (muon_weight+muon_weight_error);
        float norm_gUPC_Sys_Stat_Down = Other_gUPC * (muon_weight-muon_weight_error);
        
        float norm_UPCGen_Sys_Stat_Up = Other_UPCGen * (muon_weight+muon_weight_error);
        float norm_UPCGen_Sys_Stat_Down = Other_UPCGen * (muon_weight-muon_weight_error);

        float norm_mumuFSR_Sys_Stat_Up = Other_mumuFSR * (muon_weight+muon_weight_error);
        float norm_mumuFSR_Sys_Stat_Down = Other_mumuFSR * (muon_weight-muon_weight_error);
      
       std::vector<double> wtSysStatUp = {1, norm_gUPC_Sys_Stat_Up, norm_mumuFSR_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up, norm_UPCGen_Sys_Stat_Up};   
       std::vector<double> wtSysStatDown = {1, norm_gUPC_Sys_Stat_Down, norm_mumuFSR_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down, norm_UPCGen_Sys_Stat_Down};   
     
      cout << "norm_UPCGen nominal:" << norm_UPCGen << endl;
     
      cout << "norm_UPCGen_Sys_Stat_Up:" << norm_UPCGen_Sys_Stat_Up << endl;
      cout << "norm_UPCGen_Sys_Stat_Down:" << norm_UPCGen_Sys_Stat_Down << endl;
      // cout << "Muon Pt:" << treeR.Mu_Pt << endl;
      // cout << "Muon Eta:" << treeR.Mu_Eta << endl;
      // cout << "softID:" << tnp_weight_softid_upc_pbpb( PT, ET, 0) << endl;
      // cout << "TriggerID:" << tnp_weight_trigger_upc_pbpb( PT, ET, 0) << endl;
      // cout << "softID stat up:" << tnp_weight_softid_upc_pbpb( PT, ET, 1) << endl;
      // cout << "Other_gUPC:" << Other_gUPC << endl;
      // cout << "wt:" << norm_gUPC << endl;
      // cout << "wtUP:" << norm_gUPC_Sys_Up << endl;
      // cout << "wtDown:" << norm_gUPC_Sys_Down << endl;
      // cout << "wt mumuFSR:" << norm_mumuFSR << endl;
      // cout << "wt mumuFSRUp:" << norm_mumuFSR_Sys_Up  << endl;
      // cout << "wt mumuFSRDown:" << norm_mumuFSR_Sys_Down << endl;
     // if(treeR.Ele_Pt < 2.5) continue;
    
    
     //float binwidth1 = hmu_pt[i]->GetBinWidth(hmu_pt[i]->GetXaxis()->FindFixBin(treeR.Mu_Pt));
     float binwidth1 = 1;
      //cout << "binwidth:" << binwidth1 << endl;
      hmu_pt[i]->Fill(treeR.Mu_Pt,wt[i]*CrossSection[i]/binwidth1);
      //float binwidth2 = hele_pt[i]->GetBinWidth(hele_pt[i]->GetXaxis()->FindFixBin(treeR.Ele_Pt));
      float binwidth2 =1;
      hele_pt[i]->Fill(treeR.Ele_Pt,wt[i]*CrossSection[i]/binwidth2);
     // if(i > 0){
     // hmu_pt[i]->Fill(treeR.Mu_Pt,wt[i]);
   
      hmu_pt_1sigmaUp_sys[i]->Fill(treeR.Mu_Pt,wtSysUp[i]*CrossSection[i]/binwidth1);
      hmu_pt_1sigmaDown_sys[i]->Fill(treeR.Mu_Pt,wtSysDown[i]*CrossSection[i]/binwidth1);
      hmu_pt_1sigmaUp_stat[i]->Fill(treeR.Mu_Pt,wtStatUp[i]*CrossSection[i]/binwidth1);
      hmu_pt_1sigmaDown_stat[i]->Fill(treeR.Mu_Pt,wtStatDown[i]*CrossSection[i]/binwidth1);
      hmu_pt_1sigmaUp_sys_stat[i]->Fill(treeR.Mu_Pt,wtSysStatUp[i]*CrossSection[i]/binwidth1);
      hmu_pt_1sigmaDown_sys_stat[i]->Fill(treeR.Mu_Pt,wtSysStatDown[i]*CrossSection[i]/binwidth1);
//      cout << "Test" << endl;
     // hele_pt[i]->Fill(treeR.Ele_Pt,wt[i]);
      hele_pt_1sigmaUp_sys[i]->Fill(treeR.Ele_Pt,wtSysUp[i]*CrossSection[i]/binwidth2);
      hele_pt_1sigmaDown_sys[i]->Fill(treeR.Ele_Pt,wtSysDown[i]*CrossSection[i]/binwidth2);
      hele_pt_1sigmaUp_stat[i]->Fill(treeR.Ele_Pt,wtStatUp[i]*CrossSection[i]/binwidth2);
      hele_pt_1sigmaDown_stat[i]->Fill(treeR.Ele_Pt,wtStatDown[i]*CrossSection[i]/binwidth2);
      hele_pt_1sigmaUp_sys_stat[i]->Fill(treeR.Ele_Pt,wtSysStatUp[i]*CrossSection[i]/binwidth2);
      hele_pt_1sigmaDown_sys_stat[i]->Fill(treeR.Ele_Pt,wtSysStatDown[i]*CrossSection[i]/binwidth2);
    
   
   }	   

  } 

   
   TCanvas*c2 = new TCanvas();
//   hmu_PtEta[2]->Draw("colz");
   TCanvas*c3 = new TCanvas();
  // hmu_pt[2]->Draw("HIST");
   cout << "test1:" << endl; 
  TCanvas* cc1 = new TCanvas("Muon_pt_Rebin","Muon pT",254,411,639,592); 
  PlotHistsAndRatio(cc1, hmu_pt[0], hmu_pt[3], hmu_pt[2], 2.5,22.5,0,35, -1, 1,"#tau_{#mu} p_{T} in [GeV/c]", 0);
  


  TCanvas* cc2 = new TCanvas("Electron_pt_Rebin","Electron pT",254,411,639,592);
  PlotHistsAndRatio(cc2, hele_pt[0], hele_pt[3], hele_pt[2], 2.0,14.0,0,28, -1, 1,"#tau_{e} p_{T} in [GeV/c]", 0);




  cout << "test2:" << endl; 
  for (int i = 0; i < 24; i++){
  //cout << " Events electron pT: " << Sample[i] <<  " :" << hele_pt[i]->Integral(1,5) << endl;
   cout << " Events muon pT Sys_Stat_Up: " << Sample[i] <<  " :" << hmu_pt_1sigmaUp_sys_stat[i]->Integral(1,20) << endl;
  cout << " Events muon pT: " << Sample[i] <<  " :" << hmu_pt[i]->Integral(1,6) << endl;
 
  cout << " Events muon pT Sys_Stat_Down: " << Sample[i] <<  " :" << hmu_pt_1sigmaDown_sys_stat[i]->Integral(1,20) << endl;
  //cout << " Events muon pT: " << Sample[i] <<  " :" << hmu_pt[i]->Integral(1,9) << endl;
  // cout << " Events muon pT UP_sys: " << Sample[i] <<  " :" << hmu_pt_1sigmaUp_sys[i]->Integral(1,40) << endl;
  // cout << " Events muon pT Down_sys: " << Sample[i] <<  " :" << hmu_pt_1sigmaDown_sys[i]->Integral(1,40) << endl;
  // cout << " Events muon pT Up_stat: " << Sample[i] <<  " :" << hmu_pt_1sigmaUp_stat[i]->Integral(1,40) << endl;
  // cout << " Events muon pT Down_stat: " << Sample[i] <<  " :" << hmu_pt_1sigmaDown_stat[i]->Integral(1,40) << endl;

 }
 TString dirName = "muon_pt";
 TString dirName2 = "electron_pt";
 outf->mkdir(dirName);
 outf->mkdir(dirName2);
 //outf->cd();
 outf->cd(dirName);
 for( int i=0; i < 24; i++){
 hmu_pt[i]->Write();

 }
 for( int i=1; i < 24; i++){
 hmu_pt_1sigmaUp_sys[i]->Write();
 hmu_pt_1sigmaDown_sys[i]->Write();
 hmu_pt_1sigmaDown_stat[i]->Write();
 hmu_pt_1sigmaUp_stat[i]->Write();
 hmu_pt_1sigmaUp_sys_stat[i]->Write();
 hmu_pt_1sigmaDown_sys_stat[i]->Write();


 }
// outf->Write();
outf->cd(); 
outf->cd(dirName2);
for( int i=0; i < 24; i++){
 hele_pt[i]->Write();
 }
 for( int i=1; i < 24; i++){
 hele_pt_1sigmaUp_sys[i]->Write();
 hele_pt_1sigmaDown_sys[i]->Write();
 hele_pt_1sigmaDown_stat[i]->Write();
 hele_pt_1sigmaUp_stat[i]->Write();
 hele_pt_1sigmaUp_sys_stat[i]->Write();
 hele_pt_1sigmaDown_sys_stat[i]->Write();
 }

 outf->Close();

}

TCanvas* PlotHistsAndRatio(TCanvas* c1, TH1D* hdata, TH1D* hgammaUPC, TH1D* hmumuFSR, double hxmin, double hxmax, double hymin, double hymax, double rmax, double rmin,  const char *ytitle, bool iflogy){
  
  float T = 0.08;
  float B = 0.14; 
  float L = 0.14;
  float R = 0.04;
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd(); 

  THStack *hs = new THStack("hs"," ");
  if(iflogy)gPad->SetLogy();


  hs->SetMaximum(hymax);
  hs->SetMinimum(hymin);
  
  make_hist(hmumuFSR, kYellow, 21, 4);
  hs->Add(hmumuFSR);
 
  make_hist(hgammaUPC,kAzure+1, 21, 4);
  //hs->SetLineColor(kRed);
  hs->Add(hgammaUPC);

 // make_hist(hmumuFSR, kYellow, 21);
 // hs->Add(hmumuFSR);
  //hs->GetXaxis()->SetTitle(ytitle);
  //hs->GetYaxis()->SetTitle("Events");
  hs->Draw("hist");
  
  hs->GetXaxis()->SetTitle(ytitle);
  hs->GetYaxis()->SetTitle("Events");
//  hs->SetLineWidth(2);
   
  make_hist(hdata, kBlack, 21, 2);
  hdata->SetLineColor(kBlack);
  hdata->Draw("e1same");

 // make_hist(hmumuFSR, kBlue, 21);
 // hmumuFSR->SetLineColor(kBlue);
 // hmumuFSR->Draw("histsame");
 
  ///
  TLegend *leg2=new TLegend(0.55,0.60,0.90,0.91);
  //TLegend *leg2=new TLegend(0.55,0.60,0.90,0.91);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(20); 
  leg2->AddEntry(hdata,"Data","pl");
  leg2->AddEntry(hgammaUPC,"Signal gUPC","f");
  leg2->AddEntry(hmumuFSR,"#gamma #gamma #rightarrow #mu^{+}#mu^{-} #gamma","f");
  leg2->Draw();
  

       TLatex *mark = new TLatex();
       mark->SetNDC(true);
       double startY = 0.92;
       mark->SetTextFont(42);
       mark->SetTextSize(0.035);
       mark->DrawLatex(0.15,startY,"#bf{CMS} #it{Work in progresss}");
       mark->DrawLatex(0.65,startY,"#scale[0.8]{1642.79 #mub^{-1} (5.02 TeV)}");
       mark->Draw();
  

  //Ratio plot
  /*
  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.06, 1, 0.29);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad2->SetGridx();
  pad2->Draw();
  pad2->cd();
  TString strRatio=hgammaUPC->GetName();
  strRatio+="_over_";
  strRatio+=hgammaUPC->GetName();
  TH1D* hratio=(TH1D*)hdata->Clone(strRatio);
  hratio->Sumw2();  
  
  hratio->Divide(hgammaUPC);
  hratio->SetTitle("");
  hratio->GetXaxis()->SetTitle(ytitle);
  hratio->GetYaxis()->SetTitle("Data/gUPC");
  hratio->SetLabelSize(0.1, "XYZ");
  hratio->SetLabelFont(42, "XYZ");
  hratio->SetLabelOffset(0.007, "XYZ");
  hratio->SetTitleSize(0.11, "XYZ");
  hratio->GetXaxis()->SetTitleOffset(1.0);
  hratio->GetYaxis()->SetTitleOffset(0.3);
//  pad2->cd();
  //pad2->SetGridy(1);
 // hratio->SetLineWidth(2);
  hratio->Draw("p");
  *///to remove ratio

  c1->Update();
  TString cName=c1->GetName();
  cName+=".png";
 // c1->SaveAs("Plots_From_selectEleMuEvents/11thDecember/"+cName);
  //c1->SaveAs("Test/"+cName);
  TString c2Name=c1->GetName();
  c2Name+=".pdf";
  c1->SaveAs("Plots_From_selectEleMuEvents/tmp/"+c2Name);
  return c1;
}
  

void make_canvas(TCanvas *& canvas){

  int W = 600;
  int H = 670;
  float T = 0.08;
  float B = 0.2; 
  float L = 0.18;
  float R = 0.04;

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L );
  canvas->SetRightMargin( R );
  canvas->SetTopMargin( T );
  canvas->SetBottomMargin( B );
  canvas->SetTickx(0);
  canvas->SetTicky(0);
}

void make_canvas_ratio(TCanvas *& can){

  int W = 700;
  int H = 600;
  float T = 0.08;
  float B = 0.14; 
  float L = 0.14;
  float R = 0.04;

  can->SetFillColor(0);
  can->SetBorderMode(0);
  can->SetFrameFillStyle(0);
  can->SetFrameBorderMode(0);
  can->SetLeftMargin( L );
  can->SetRightMargin( R );
  can->SetTopMargin( T );
  can->SetBottomMargin( B );
  can->SetTickx(0);
  can->SetTicky(0);
}


void make_hist(TH1D *& hist, Color_t kcolor, int kstyle, int lwidth){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.1, "XYZ");
 
  hist->SetMarkerColor(kcolor);
  hist->SetFillColor(kcolor);
  hist->SetLineColor(kBlack);
  hist->SetLineColor(kBlack);
  //hist->SetLineWidth(2);
  hist->SetMarkerStyle(kstyle);
  hist->SetBinErrorOption(TH1::kPoisson); 
  hist->SetLineWidth(1);
}

void make_hist_ratio(TH1D *& hist, Color_t kcolor, int kstyle){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.08, "XYZ");
  hist->SetLabelSize(0.1, "XYZ");

  // For the axis titles:
  hist->SetTitleFont(42, "XYZ");
  hist->SetTitleSize(0.11, "XYZ");
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleOffset(0.4);
  hist->SetMarkerColor(kcolor);
  hist->SetFillColor(kcolor);
  hist->SetMarkerStyle(kstyle);
   
}

