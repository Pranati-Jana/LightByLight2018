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
//#include "tnp_electron.h"

int xlo = 1;
int xhi = 2;
int nbin = 8;


void make_canvas(TCanvas *&);
void make_hist(TH1D *&, Color_t , int , int);
void make_canvas_ratio(TCanvas *&);
void make_hist_ratio(TH1D *&, Color_t, int );
//TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D* , TH1D*, TH1D*, TH1D*, TH1D*, double , double , double , double , double , double , const char *, bool); 
// const int nPtbins=6;
//  double Ptbin[nPtbins]={2.5,3.5,4.5,5.5,9.5,15.5};
//  const int nPtbin= sizeof(Ptbin)/sizeof(double) - 1;
TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, TH1D*, double , double, double , double , double, double, const char *,const char *, bool); 

const int nSample = 3;
//const char *Sample[nSample]={"Data_OneFourthLumi_DeltaR1HE3","gamgamToTauTau0_LByLCuts_DeltaR1HE3","MUMU_FSR_LByLCuts_DeltaR1HE3"};
//const char *Sample[nSample]={"data_obs","GAMMA_UPC","MUMU_FSR"};
const char *Sample[nSample]={"data_obs","gamgamToTauTau0","MUMU_FSR"};


const char *dir = "figures_eta2p2";

double luminosity = 1700.00; //ub
//double luminosity = 425.91;//
double crossSectionGammaUPC = 1060.8959;//GammaUPC
double generatedEvtGammaUPC = 5290000;//GammaUPC

double ScaleFactorSC= 0.943*0.961237*0.9*0.99;
double Other_gUPC = (luminosity * crossSectionGammaUPC * ScaleFactorSC)/generatedEvtGammaUPC;//gUPC

double crossSectionFSR = 139.7;//MC_FSR_Arash
double generatedEvtFSR = 2471810;//FSR

double Other_mumuFSR = (luminosity * crossSectionFSR*ScaleFactorSC)/generatedEvtFSR;//MC_Arash_FSR
//double Other_UPCGen = (luminosity * 484 *ScaleFactorSC)/10000000;
double Other_UPCGen = (luminosity * 564.89 *ScaleFactorSC)/10000000;
//double Other_UPCGen =1;
 const int nPtbins=18;

 //double Ptbin[nPtbins]={1.5,2.5,3.5,4.5,5.5,22.5};
 double Ptbin[nPtbins]={0.0,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,20.5};
 const int nPtbin= sizeof(Ptbin)/sizeof(double) - 1;
void BinLogX(TH1* h);


void plot_FromTree(){
  
  //const double wt[nSample] = {1};
 cout << "Other gUPC:" << Other_gUPC << endl;
 cout << "Other mumuFSR:" << Other_mumuFSR << endl;
  double ET;
  double PT;
  double ETe;
  double PTe;

  int idx = 0;
    
	
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  string mc_name = "SC";
  TFile *outf= new TFile("EleMu_histos_gUPC_signal_0-1.5-2.5-3.5-4.5-5.5-6.5-7.5-8.5-9.5-10.5-11.5-12.5-13.5-14.5-15.5-16.5-20.5GeV.root","recreate");  
//  TFile* f[nSample];
  TChain *tree[nSample];

  //// histograms
  TH1D* hmu_pt[nSample], *hele_pt[nSample], *helemu_pt[nSample], *helemuScalar_pt[nSample], *helemu_mass[nSample], *helemu_aco[nSample], *hele_eta[nSample], *hele_phi[nSample], *hmu_eta[nSample], *hmu_phi[nSample], *helemu_rapidity[nSample], *helemu_phi[nSample], *helemu_pz[nSample], *hZDC_XnXn[nSample], *hDeltaR_TrkElectron[nSample], *hDeltaR_TrkMuon[nSample], *hnTracks[nSample];
  TH2D* hmu_PtEta[nSample], *hZDC[nSample], *hDeltaR_HE[nSample];
  TH1D* h_LeadingHE_Energy[nSample], *h_DeltaRMuon_LeadingHE[nSample];
  
  for (int i = 0; i < nSample; i++){
    cout << "nSample:" << i  << endl;
   tree[i] = new TChain("output_tree");
   //tree[i]->Add(Form("/home/pranati/Documents/GammaGammaToTauTau/Workspace/plottingMacros/RootFile_From_selectEleMuEvents/11thDecember_From_selectEleMuEvents_Tree/%s_EleMu.root", Sample[i]));   
   //tree[i]->Add(Form("RootFile_From_selectEleMuEvents/AllSelections_FromSelectEleMuScript/%s_EleMu.root", Sample[i])); 
   //tree[i]->Add(Form("RootFile_From_selectEleMuEvents/NewCutFlowHistadded/%s_EleMu.root", Sample[i]));  
  //tree[i]->Add(Form("/Users/pranatijana/Documents/GammaGamma_to_Tautau/Workspace/plottingMacros/RootFile_From_selectEleMuEvents/ExtraTrack_ExclusivityAcopCut_NotApplied/loose_ElectronMatchingToTower/%s_EleMu.root",Sample[i]));
   //tree[i]->Add(Form("/Users/pranatijana/Documents/GammaGamma_to_Tautau/Workspace/plottingMacros/RootFile_From_selectEleMuEvents/ExtraTrack_ExclusivityAcopCut_NotApplied/loose_ElectronMatchingToTower_deltaRMatchingTowerToMuon/%s_EleMu.root",Sample[i]));
  //tree[i]->Add(Form("/Users/pranatijana/Documents/GammaGamma_to_Tautau/Workspace/plottingMacros/RootFile_From_selectEleMuEvents/Leading_HE/%s_EleMu.root",Sample[i]));
  // tree[i]->Add(Form("RootFile_From_selectEleMuEvents/ExtraTrack_ExclusivityAcopCut_NotApplied/21aTau_Points_UPCGen_ChargedExcElectronFromLByLMuonDeltaR0.001_CaloTowerInfo/%s_EleMu.root", Sample[i])); 
  //tree[i]->Add(Form("RootFile_From_selectEleMuEvents/loose_ElectronMatchingToTower_deltaRMatchingTowerToMuon/%s_EleMu.root", Sample[i]));
 //tree[i]->Add(Form("RootFile_From_selectEleMuEvents/MuonEndcapPt1.5GeV_HE1To3DeltaRMuonHE1/%s_EleMu.root", Sample[i]));
//tree[i]->Add(Form("RootFile_From_selectEleMuEvents/NewRunningHFThreshold90Percent/%s_EleMu.root", Sample[i]));
tree[i]->Add(Form("RootFile_From_selectEleMuEvents/a_Tau_Points_6.4GeVHFcut/%s_EleMu.root", Sample[i]));
  // hmu_pt[i] = new TH1D(Form("hmu_pt%s", Sample[i]),"",19,1.5,20.5);
   hmu_pt[i] = new TH1D(Form("hmu_pt"),"",nPtbin, Ptbin);
   hmu_eta[i] = new TH1D(Form("hmu_eta"),"",13,-2.4,2.4);//12
   hmu_phi[i] = new TH1D(Form("hmu_phi"),"",13,-3.5,3.5);//8
   hele_pt[i] = new TH1D(Form("hele_pt"),"",16,2.5,18.5);
   //hele_pt[i] = new TH1D(Form("hele_pt%s", Sample[i]),"",nPtbin,Ptbin);
   hele_eta[i] = new TH1D(Form("hele_eta"),"",11,-2.2,2.2);//11
   hele_phi[i] = new TH1D(Form("hele_phi"),"",13,-3.5,3.5);//8
   helemu_pt[i] = new TH1D(Form("helemu_pt"),"",14,0.0,14.0);
   helemuScalar_pt[i] = new TH1D(Form("helemuScalar_pt"),"",11,0.0,22.0);
   helemu_mass[i] = new TH1D(Form("helemu_mass"),"",11,2,24.0);
   helemu_aco[i] = new TH1D(Form("helemu_aco"),"",20,0.0,0.2);
   helemu_rapidity[i] = new TH1D(Form("helemu_rapidity"),"",11,-2.5,2.5);
   helemu_phi[i] = new TH1D(Form("helemu_phi"),"",10,2.6,3.2);
   helemu_pz[i] = new TH1D(Form("helemu_pz"),"",22,-30,30);
   hmu_PtEta[i] = new TH2D(Form("hmu_PtEta"), "",14,2.5,16.5,180,-2.5,2.5);
   hZDC[i]      = new TH2D(Form("hZDC%s", Sample[i]), "",87,0,130000,43,0,60000);
   hZDC_XnXn[i] = new TH1D(Form("hZDC_XnXn%s", Sample[i]),"",3,0,3);
   hDeltaR_TrkElectron[i] = new TH1D(Form("hDeltaR_TrkElectron%s", Sample[i]),"",30,0,0.05);
   hDeltaR_TrkMuon[i] = new TH1D(Form("hDeltaR_TrkMuon%s", Sample[i]),"",10,0,0.01);
   hnTracks[i] = new TH1D(Form("hnTracks%s", Sample[i]),"",10,0,10);
   h_LeadingHE_Energy[i] = new TH1D(Form("h_LeadingHE_Energy%s", Sample[i]),"",10,0,10);
   h_DeltaRMuon_LeadingHE[i]= new TH1D(Form("h_DeltaRMuon_LeadingHE%s", Sample[i]),"",60,0,6);
   hDeltaR_HE[i] = new TH2D(Form("hDeltaR_HE%s", Sample[i]), "",60,0,6,100,0,10);
////////////////////////
// hmu_pt[i] = new TH1D(Form("hmu_pt%s", Sample[i]),"",nPtbin, Ptbin);
//    hmu_eta[i] = new TH1D(Form("hmu_eta%s", Sample[i]),"",13,-2.4,2.4);//12
//    hmu_phi[i] = new TH1D(Form("hmu_phi%s", Sample[i]),"",13,-3.5,3.5);//8
//    hele_pt[i] = new TH1D(Form("hele_pt%s", Sample[i]),"",16,2.5,18.5);
//    //hele_pt[i] = new TH1D(Form("hele_pt%s", Sample[i]),"",nPtbin,Ptbin);
//    hele_eta[i] = new TH1D(Form("hele_eta%s", Sample[i]),"",11,-2.2,2.2);//11
//    hele_phi[i] = new TH1D(Form("hele_phi%s", Sample[i]),"",13,-3.5,3.5);//8
//    helemu_pt[i] = new TH1D(Form("helemu_pt%s", Sample[i]),"",14,0.0,14.0);
//    helemuScalar_pt[i] = new TH1D(Form("helemuScalar_pt%s", Sample[i]),"",11,0.0,22.0);
//    helemu_mass[i] = new TH1D(Form("helemu_mass%s", Sample[i]),"",11,2,24.0);
//    helemu_aco[i] = new TH1D(Form("helemu_aco%s", Sample[i]),"",20,0.0,0.2);
//    helemu_rapidity[i] = new TH1D(Form("helemu_rapidity%s", Sample[i]),"",11,-2.5,2.5);
//    helemu_phi[i] = new TH1D(Form("helemu_phi%s", Sample[i]),"",10,2.6,3.2);
//    helemu_pz[i] = new TH1D(Form("helemu_pz%s", Sample[i]),"",22,-30,30);
//    hmu_PtEta[i] = new TH2D(Form("hmu_PtEta%s", Sample[i]), "",14,2.5,16.5,180,-2.5,2.5);
//    hZDC[i]      = new TH2D(Form("hZDC%s", Sample[i]), "",87,0,130000,43,0,60000);
//    hZDC_XnXn[i] = new TH1D(Form("hZDC_XnXn%s", Sample[i]),"",3,0,3);
//    hDeltaR_TrkElectron[i] = new TH1D(Form("hDeltaR_TrkElectron%s", Sample[i]),"",30,0,0.05);
//    hDeltaR_TrkMuon[i] = new TH1D(Form("hDeltaR_TrkMuon%s", Sample[i]),"",10,0,0.01);
//    hnTracks[i] = new TH1D(Form("hnTracks%s", Sample[i]),"",10,0,10);
//    h_LeadingHE_Energy[i] = new TH1D(Form("h_LeadingHE_Energy%s", Sample[i]),"",10,0,10);
//    h_DeltaRMuon_LeadingHE[i]= new TH1D(Form("h_DeltaRMuon_LeadingHE%s", Sample[i]),"",60,0,6);
//    hDeltaR_HE[i] = new TH2D(Form("hDeltaR_HE%s", Sample[i]), "",60,0,6,100,0,10);
   //////////

   cout << "nSample:" << i  << endl;
   cout << "file " << Sample[i] << ":" << tree[i]->GetEntries()  << endl;
   ReadTree  treeR(tree[i]);
   treeR.fChain->SetBranchStatus("*",1);
   if (treeR.fChain == 0) return;
   Long64_t nentries = treeR.fChain->GetEntriesFast();
   cout << tree[i]->GetName() << "    " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
  
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<62;jentry++) {
      Long64_t ientry_evt = treeR.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = treeR.fChain->GetEntry(jentry);   nbytes += nb;

      //SF
      ET = treeR.Mu_Eta;
      PT = treeR.Mu_Pt;
      ETe = treeR.Ele_Eta;
      PTe = treeR.Ele_Pt;

      //double norm_gUPC = tnp_weight_softid_upc_pbpb( PT, ET, idx) * tnp_weight_trigger_upc_pbpb( PT, ET, idx) * Other_gUPC * tnp_electron( PTe, ETe, idx);
     
     
      //std::vector<double> wt = {1, 1, 1};
      //if(treeR.Ele_Pt < 2.5) continue;
     if(treeR.ok_neuexcl!=1) continue;
     if(treeR.ok_neuexcl_HFonly!=1) continue;
     if(treeR.ok_chexcl_goodtracks!=1) continue;
    //if(treeR.ok_chexcl_extratracks!=1) continue;
    // if(treeR.nTracks <=2) continue;
     if(treeR.EleMu_acop < 0.01) continue; 
     if(treeR.EleMu_acop > 0.15) continue;
     //if(treeR.tower_had > 3.0) continue;
      if(i==0){
        if(treeR.tower_HF_Minus > 6.4) continue;
        if(treeR.tower_HF_Plus > 6.4) continue;
     //  if( treeR.zdc_energy_pos > 4200 || treeR.zdc_energy_neg > 5000) continue;
      } 
      double norm_gUPC = tnp_weight_softid_upc_pbpb( PT, ET, idx) * tnp_weight_trigger_upc_pbpb( PT, ET, idx) * Other_gUPC;
     // double norm_gUPC =  Other_gUPC;
      //double norm_mumuFSR = tnp_weight_softid_upc_pbpb( PT, ET, idx) * tnp_weight_trigger_upc_pbpb( PT, ET, idx) * Other_mumuFSR * tnp_electron( PTe, ETe, idx);
      double norm_mumuFSR = tnp_weight_softid_upc_pbpb( PT, ET, idx) * tnp_weight_trigger_upc_pbpb( PT, ET, idx) * Other_mumuFSR;
      double norm_UPCGen = tnp_weight_softid_upc_pbpb( PT, ET, idx) * tnp_weight_trigger_upc_pbpb( PT, ET, idx) * Other_UPCGen;
    //  if(PT > 1.5 && PT < 2.5){
    //   cout << "Trigger SF:" << tnp_weight_trigger_upc_pbpb( PT, ET, idx) << endl;
    //   cout << "SoftID SF:" << tnp_weight_softid_upc_pbpb( PT, ET, idx) << endl;

    //  }

     std::vector<double> wt = {1, norm_UPCGen, norm_mumuFSR};
     //std::vector<double> wt = {1, 1, 1};
      // if((treeR.Mu_Eta) > 1.2 && treeR.Mu_Pt > 2.5 ){
      // cout << "Muon Pt:" << treeR.Mu_Pt << endl;
      // cout << "Muon Eta:" << treeR.Mu_Eta << endl;
      //cout << "Electron Pt:" << treeR.Ele_Pt << endl;
      // cout << "softID:" << tnp_weight_softid_upc_pbpb( PT, ET, idx) << endl;
      // cout << "TriggerID:" << tnp_weight_trigger_upc_pbpb( PT, ET, idx) << endl;
      // cout << "Other_gUPC:" << Other_gUPC << endl;
      // cout << "Electron SF:" << tnp_electron( PTe, ETe, idx) << endl;
      // cout << "wt:" << norm_gUPC << endl;
      // cout << "wt mumuFSR:" << norm_mumuFSR << endl;
        //    }
      
      hmu_pt[i]->Fill(treeR.Mu_Pt,wt[i]);
      hmu_eta[i]->Fill(treeR.Mu_Eta,wt[i]);
      hmu_phi[i]->Fill(treeR.Mu_Phi,wt[i]);
      float binwidth2 = hele_pt[i]->GetBinWidth(hele_pt[i]->GetXaxis()->FindFixBin(treeR.Ele_Pt));
      hele_pt[i]->Fill(treeR.Ele_Pt,wt[i]);
      hele_eta[i]->Fill(treeR.Ele_Eta,wt[i]);
      hele_phi[i]->Fill(treeR.Ele_Phi,wt[i]);
      helemu_pt[i]->Fill(treeR.vSum_Pt,wt[i]);
      helemuScalar_pt[i]->Fill(treeR.Scalar_Pt,wt[i]);
      helemu_mass[i]->Fill(treeR.vSum_M,wt[i]);
      helemu_aco[i]->Fill(treeR.EleMu_acop,wt[i]);
      helemu_rapidity[i]->Fill(treeR.vSum_Rapidity,wt[i]);
      helemu_phi[i]->Fill(treeR.EleMu_dphi,wt[i]);
      helemu_pz[i]->Fill(treeR.vSum_Pz,wt[i]);//Change to the selectEleMu code
      hmu_PtEta[i]->Fill(treeR.Mu_Pt,treeR.Mu_Eta,wt[i]);
      hDeltaR_TrkElectron[i]->Fill(treeR.deltaRtrack_electron,wt[i]);
      hDeltaR_TrkMuon[i]->Fill(treeR.deltaRtrack_muon,wt[i]);
      hnTracks[i]->Fill(treeR.nTracks,wt[i]);
      h_LeadingHE_Energy[i]->Fill(treeR.tower_had,wt[i]);
      h_DeltaRMuon_LeadingHE[i]->Fill(treeR.tower_deltaRmuon_HE,wt[i]);
      hDeltaR_HE[i]->Fill(treeR.tower_deltaRmuon_HE,treeR.tower_had,wt[i]);
     
      if(i==0){
        hZDC[i]->Fill(treeR.zdc_energy_pos,treeR.zdc_energy_neg);
        float zdc_plus = treeR.zdc_energy_pos;
        float zdc_minus = treeR.zdc_energy_neg;
        int nNeutronsPlus = 0;
        if (zdc_plus <  1600) nNeutronsPlus = 0;
        else
         nNeutronsPlus = 1;
        int nNeutronsMinus = 0;
        if (zdc_minus <  1600) nNeutronsMinus = 0;
        else nNeutronsMinus = 1;

        if (nNeutronsPlus == 0 && nNeutronsMinus == 0) hZDC_XnXn[0]->Fill(0);  // 0n0n
        if ((nNeutronsPlus == 0 && nNeutronsMinus == 1 ) || (nNeutronsPlus == 1 && nNeutronsMinus == 0)) hZDC_XnXn[0]->Fill(1);  // 0nXn,Xn0n
        if (nNeutronsPlus == 1 && nNeutronsMinus == 1) hZDC_XnXn[0]->Fill(2);
        // hZDC_XnXn[0]->Fill(2);
       // treeR.zdc_energy_pos > 1600 || treeR.zdc_energy_neg > 1600;
      }
     
      

   }	   

  } 

   
   TCanvas*c2 = new TCanvas();
   hmu_PtEta[1]->Draw("colz");
  // c2->SaveAs("Plots_From_selectEleMuEvents/LeadingHE/Muon_Pt_Eta.pdf");

   TLegend *leg=new TLegend(0.55,0.60,0.85,0.91);
  //TLegend *leg2=new TLegend(0.55,0.60,0.90,0.91);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSize(15); 
  leg->AddEntry("Data");
  
  //leg2->Draw();

   /////Data
   TCanvas*c6 = new TCanvas();
   hDeltaR_HE[0]->GetYaxis()->SetTitleOffset(1.4);
   hDeltaR_HE[0]->GetXaxis()->SetTitle("#DeltaR(#tau_{#mu},Leading HE)");
   hDeltaR_HE[0]->GetYaxis()->SetTitle("Leading HE");
   hDeltaR_HE[0]->SetTitle("Data");
   hDeltaR_HE[0]->Draw("colz");
   //leg->Draw();
   //c6->SaveAs("Plots_From_selectEleMuEvents/LByL_ChargedExclusivity_DeltaRMuon1HE3/Data_DeltaR_HE.pdf");
   //c6->SaveAs("Plots_From_selectEleMuEvents/LeadingHE/Data_DeltaR_HE_2Track.pdf");
   ///MC
   TCanvas*c7 = new TCanvas();
   hDeltaR_HE[1]->GetYaxis()->SetTitleOffset(1.4);
   hDeltaR_HE[1]->GetXaxis()->SetTitle("#DeltaR(#tau_{#mu},Leading HE)");
   hDeltaR_HE[1]->GetYaxis()->SetTitle("Leading HE");
   hDeltaR_HE[1]->SetTitle("Signal MC, #gamma #gamma #rightarrow #tau^{+}#tau^{-}");
   hDeltaR_HE[1]->Draw("colz");
  // c7->SaveAs("Plots_From_selectEleMuEvents/LByL_ChargedExclusivity_DeltaRMuon1HE3/MC_DeltaR_HE.pdf");
   //c7->SaveAs("Plots_From_selectEleMuEvents/LeadingHE/MC_DeltaR_HE_2Track.pdf");
   //////
  
   
   TCanvas*c5 = new TCanvas();
   hZDC[0]->GetYaxis()->SetTitleOffset(1.4);
   hZDC[0]->GetXaxis()->SetTitle("ZDC^{+}");
   hZDC[0]->GetYaxis()->SetTitle("ZDC^{-}");
   hZDC[0]->Draw("colz");
   //c5->SaveAs("Plots_From_selectEleMuEvents/tmp/Muon_Pt_Eta.pdf");
   /////////
   //c3->SaveAs("Plots_From_selectEleMuEvents/OnlyRemovingDeltaEtaForElectronAndECALMatching/ZDC.pdf");
   TCanvas*c4 = new TCanvas();
    hZDC_XnXn[0]->SetMarkerColor(kBlue);
    hZDC_XnXn[0]->SetMarkerStyle(20);
    hZDC_XnXn[0]->GetYaxis()->SetRangeUser(0,30);
    hZDC_XnXn[0]->GetXaxis()->SetRangeUser(0,3);
    hZDC_XnXn[0]->GetXaxis()->SetLabelSize(0);
    hZDC_XnXn[0]->GetYaxis()->SetTitle("Events");
   hZDC_XnXn[0]->Draw("P");
   //c4->SaveAs("Plots_From_selectEleMuEvents/OnlyRemovingDeltaEtaForElectronAndECALMatching/ZDC_XnXn.pdf");
 //outf->cd();
 TString dirName = "data";
 TString dirName1 = "signal";
 TString dirName2 = "mumuFSR";
 outf->mkdir(dirName);
 outf->mkdir(dirName1);
 outf->mkdir(dirName2);

 outf->cd(dirName);
 hmu_pt[0]->Write();
 hele_pt[0]->Write();
 helemu_pt[0]->Write();
 helemuScalar_pt[0]->Write();
 helemu_mass[0]->Write();
 helemu_aco[0]->Write();
 hmu_eta[0]->Write();
 hmu_phi[0]->Write();
 hele_eta[0]->Write();
 hele_phi[0]->Write();
 helemu_phi[0]->Write();
 helemu_pz[0]->Write();
 helemu_rapidity[0]->Write();
 
outf->cd(); 
outf->cd(dirName1);
 hmu_pt[1]->Write();
 hele_pt[1]->Write();
 helemu_pt[1]->Write();
 helemuScalar_pt[1]->Write();
 helemu_mass[1]->Write();
 helemu_aco[1]->Write();
 hmu_eta[1]->Write();
 hmu_phi[1]->Write();
 hele_eta[1]->Write();
 hele_phi[1]->Write();
 helemu_phi[1]->Write();
 helemu_pz[1]->Write();
 helemu_rapidity[1]->Write();
 
 outf->cd(); 
 outf->cd(dirName2);
 hmu_pt[2]->Write();
 hele_pt[2]->Write();
 helemu_pt[2]->Write();
 helemuScalar_pt[2]->Write();
 helemu_mass[2]->Write();
 helemu_aco[2]->Write();
 hmu_eta[2]->Write();
 hmu_phi[2]->Write();
 hele_eta[2]->Write();
 hele_phi[2]->Write();
 helemu_phi[2]->Write();
 helemu_pz[2]->Write();
 helemu_rapidity[2]->Write();


//  hmu_pt[0]->Write();
//  hmu_pt[1]->Write();
//  hmu_pt[2]->Write();
 outf->Write();
 
   cout << "test1:" << endl; 
   
  TCanvas* cc1 = new TCanvas("Muon_pt_Full_Lumi_rebin_6.4GeVHFThreshold","Muon pT",254,411,639,592); 
  PlotHistsAndRatio(cc1, hmu_pt[0], hmu_pt[1], hmu_pt[2], 1.5,20.5,0,50, -1, 1,"#tau_{#mu}  p_{T}  [GeV]","Events / bin", 0);//35
  TCanvas* cc2 = new TCanvas("Electron_pt_Full_Lumi_6.4GeVHFThreshold","Electron pT_rebin",254,411,639,592);
  PlotHistsAndRatio(cc2, hele_pt[0], hele_pt[1], hele_pt[2], 2.5,18.5,0,55, -1, 1," #tau_{e} p_{T} [GeV]", "Events / 1GeV", 0);//28
  TCanvas* cc3 = new TCanvas("TauTau_VectorSum_pt_Full_Lumi_6.4GeVHFThreshold","TauTau Vector pT",254,411,639,592);
  PlotHistsAndRatio(cc3, helemu_pt[0], helemu_pt[1], helemu_pt[2], 0,12.0,0,55,-1,1,"Visible #tau_{#mu}#tau_{e}  p_{T} [GeV]","Events / 1GeV", 0);
  TCanvas* cc4 = new TCanvas("TauTau_ScalarSum_pt_Full_Lumi_6.4GeVHFThreshold","TauTau Scalar pT",254,411,639,592);
  PlotHistsAndRatio(cc4, helemuScalar_pt[0], helemuScalar_pt[1], helemuScalar_pt[2], 0,20.0,0,50,-1,1,"Scalar sum p_{T} [GeV]","Events / 2GeV", 0);
  TCanvas* cc5 = new TCanvas("TauTau_mass_Full_Lumi_6.4GeVHFThreshold","TauTau mass",254,411,639,592);
  PlotHistsAndRatio(cc5, helemu_mass[0], helemu_mass[1], helemu_mass[2], 0,24.0,0,40,-1,1,"Visible mass (#tau_{e},#tau_{#mu}) [GeV]","Events / 2GeV", 0);
  
  TCanvas* cc6 = new TCanvas("TauTau_acoplanarity_Full_Lumi_6.4GeVHFThreshold","TauTau acoplanarity",254,411,639,592);
  PlotHistsAndRatio(cc6, helemu_aco[0], helemu_aco[1], helemu_aco[2], 0,0.40,0,30,-1,1,"Acoplanarity_{#phi}(#tau_{e},#tau_{#mu})","Events / 0.01", 0);//15
  TCanvas* cc7 = new TCanvas("Muon_eta_Full_Lumi_6.4GeVHFThreshold","Muon Eta",254,411,639,592);
  PlotHistsAndRatio(cc7, hmu_eta[0], hmu_eta[1], hmu_eta[2], -2.5,2.5,0,25,-1,1,"#tau_{#mu} #eta","Events / 0.37", 0);//18
  TCanvas* cc8 = new TCanvas("Muon_phi_Full_Lumi_6.4GeVHFThreshold","Muon phi",254,411,639,592);
  PlotHistsAndRatio(cc8, hmu_phi[0], hmu_phi[1], hmu_phi[2], -3.5,3.5,0,25,-1,1,"#tau_{#mu} #phi","Events / 0.54", 0);//20
  TCanvas* cc9 = new TCanvas("Electron_eta_Full_Lumi_6.4GeVHFThreshold","Electron Eta",254,411,639,592);
  PlotHistsAndRatio(cc9, hele_eta[0], hele_eta[1], hele_eta[2], -2.5,-2.5,0,30,-1,1,"#tau_{e} #eta","Events / 0.4", 0);//18
  TCanvas* cc10 = new TCanvas("Electron_phi_Full_Lumi_6.4GeVHFThreshold","Electron phi",254,411,639,592);
  PlotHistsAndRatio(cc10, hele_phi[0], hele_phi[1], hele_phi[2], -3.5,-3.5,0,30,-1,1,"#tau_{e} #phi","Events / 0.54", 0);//
  TCanvas* cc11 = new TCanvas("TauTau_phi_Full_Lumi_6.4GeVHFThreshold","TauTau phi",254,411,639,592);
  PlotHistsAndRatio(cc11, helemu_phi[0], helemu_phi[1], helemu_phi[2], 2.6,3.2,0,50,-1,1,"#Delta#phi(e,#mu)","Events / 0.06", 0);
  TCanvas* cc12 = new TCanvas("TauTau_rapidity_Full_Lumi_6.4GeVHFThreshold","TauTau rapidity",254,411,639,592);
  PlotHistsAndRatio(cc12, helemu_rapidity[0], helemu_rapidity[1], helemu_rapidity[2], -2.5,2.5,0,30,-1,1,"y^{#tau_{e}#tau_{#mu}}","Events / 0.454", 0);
  TCanvas* cc13 = new TCanvas("TauTau_pz_Full_Lumi_6.4GeVHFThreshold","TauTau pz",254,411,639,592);
  PlotHistsAndRatio(cc13, helemu_pz[0], helemu_pz[1], helemu_pz[2], -30.0,30.0,0,20,-1,1,"#tau_{e}#tau_{#mu} Pz ","Events / 2.72", 0);
  TCanvas* cc14 = new TCanvas("DeltaR_Electron_Full_Lumi_6.4GeVHFThreshold","DeltaR_Electron",254,411,639,592);
  PlotHistsAndRatio(cc14, hDeltaR_TrkElectron[0], hDeltaR_TrkElectron[1], hDeltaR_TrkElectron[2], 0.0,0.05,0,30,-1,1,"#DeltaR(#tau_{e^{#pm}},trk^{#pm}) ","Events / 0.1", 0);
  TCanvas* cc15 = new TCanvas("DeltaR_Muon_Full_Lumi_6.4GeVHFThreshold","DeltaR_Muon",254,411,639,592);
  PlotHistsAndRatio(cc15, hDeltaR_TrkMuon[0], hDeltaR_TrkMuon[1], hDeltaR_TrkMuon[2], 0.0,0.05,0,30,-1,1,"#DeltaR(#tau_{#mu^{#pm}},trk^{#pm}) ","Events / 0.1", 0);
  TCanvas* cc16 = new TCanvas("nTracks_Full_Lumi_6.4GeVHFThreshold","nTracks",254,411,639,592);
  PlotHistsAndRatio(cc16, hnTracks[0], hnTracks[1], hnTracks[2], 0.0,10,0,200,-1,1,"nTracks ","Events / 1", 0);
  TCanvas* cc17 = new TCanvas("LeadingHE_Full_Lumi_6.4GeVHFThreshold","LeadingHE",254,411,639,592);
  PlotHistsAndRatio(cc17, h_LeadingHE_Energy[0], h_LeadingHE_Energy[1], h_LeadingHE_Energy[2], 0.0,10,0,50,-1,1,"Leading HE [GeV] ","Events / 0.1", 0);
  TCanvas* cc18 = new TCanvas("DeltaR_MuonLeadingHE_Full_Lumi_6.4GeVHFThreshold","nDeltaR_MuonLeadingHE",254,411,639,592);
  PlotHistsAndRatio(cc18, h_DeltaRMuon_LeadingHE[0], h_DeltaRMuon_LeadingHE[1], h_DeltaRMuon_LeadingHE[2], 0.0,6,0,30,-1,1,"#DeltaR(#tau_{#mu^{#pm}},Leading HE)","Events / 0.1", 0);



  //cout << "test2:" << endl; 
  for (int i = 0; i < nSample; i++){
  cout << " Events: " << Sample[i] <<  " :" << helemu_aco[i]->Integral(1,40) << endl;
  cout << " Events from pt: " << Sample[i] <<  " :" << hmu_pt[i]->Integral(1,19) << endl;
 // cout << "Bin content: " << Sample[i] << " : " <<  hmu_pt[i]->GetBinContent(4) << endl;
 }
 
 outf->Close();

}

TCanvas* PlotHistsAndRatio(TCanvas* c1, TH1D* hdata, TH1D* hgammaUPC, TH1D* hmumuFSR, double hxmin, double hxmax, double hymin, double hymax, double rmax, double rmin,  const char *xtitle,  const char *ytitle,bool iflogy){
  
  float T = 0.08;
  float B = 0.14; 
  float L = 0.50;
  float R = 0.04;
 // canvas->SetLeftMargin( R );
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.019);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.05);
  pad1->Draw();
  pad1->cd(); 

  THStack *hs = new THStack("hs"," ");
  if(iflogy)gPad->SetLogy();


  hs->SetMaximum(hymax);
  hs->SetMinimum(hymin);
  //hs->GetXaxis()->SetLabelSize(0);
  make_hist(hmumuFSR, kCyan+2, 0, 2);//kRed-7
  hs->Add(hmumuFSR);
  make_hist(hgammaUPC,kYellow+1, 0, 2);//kCyan+2
  //hs->SetLineColor(kRed);
  //hs->Draw("e1samex0");
  hs->Add(hgammaUPC);
 
  TH1D *hs_uncertainty = (TH1D*)hgammaUPC->Clone();
  for(int iBin=1; iBin<=hgammaUPC->GetNbinsX(); iBin++){
    double error = hgammaUPC->GetBinError(iBin);
    hs_uncertainty->SetBinError(iBin, error);
  } 
  cout << "Test:" << endl;
  hs_uncertainty->SetFillColorAlpha(kBlue-4, 0.5);
  hs_uncertainty->SetFillStyle(3144); 
  hs_uncertainty->SetLineWidth(2);

  hs->Draw("hist");
  hs->Draw("e1same");
  hs_uncertainty->Draw("e2same");
 
  hs->GetXaxis()->SetTitle(xtitle);
  hs->GetYaxis()->SetTitle(ytitle);
  hs->GetXaxis()->SetTitleSize(0.06);
   hs->GetYaxis()->SetTitleSize(0.06);
   hs->GetYaxis()->SetLabelSize(0.06);
   hs->GetXaxis()->SetLabelSize(0.06);
   hs->GetXaxis()->SetLabelSize(0);
//  hs->SetLineWidth(2);


  make_hist(hdata, kBlack, 20, 2);
  hdata->SetLineColor(kBlack);
  hdata->Draw("e1samex0");

 // make_hist(hmumuFSR, kBlue, 21);
 // hmumuFSR->SetLineColor(kBlue);
 // hmumuFSR->Draw("histsame");
 
  ///
  TLegend *leg2=new TLegend(0.55,0.60,0.85,0.91);
  //TLegend *leg2=new TLegend(0.55,0.60,0.90,0.91);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(18); 
  leg2->AddEntry(hdata,"Data","EP");
  leg2->AddEntry(hgammaUPC,"#gamma #gamma #rightarrow #tau^{+}#tau^{-}, ChFF","f");
  leg2->AddEntry(hmumuFSR,"#gamma #gamma #rightarrow #mu^{+}#mu^{-} #gamma","f");
  leg2->AddEntry(hs_uncertainty,"total prefit","f");
  leg2->Draw();
    TLatex *mark = new TLatex();
       mark->SetNDC(true);
       double startY = 0.92;
       mark->SetTextFont(42);
       mark->SetTextSize(0.05);
       mark->DrawLatex(0.15,startY,"#bf{CMS} #it{Preliminary}");
       mark->DrawLatex(0.58,startY,"#scale[0.8]{PbPb, 1.7nb^{-1}  (#sqrt{s_{NN}} = 5.02 TeV)}");
       mark->Draw();

  //Ratio plot
  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.06, 1, 0.29);
  pad2->SetTopMargin(0.1);
  pad2->SetBottomMargin(0.4);
  pad2->SetLeftMargin(0.14);
  pad2->SetRightMargin(0.05);
  pad2->SetGridx();
  pad2->Draw();
  pad2->cd();
  auto line = new TF1("line", "1", -100, 100);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  TH1D* hRatio = (TH1D*) hgammaUPC->Clone();
  hRatio->Add(hmumuFSR);

  TString strRatio=hgammaUPC->GetName();
  strRatio+="_over_";
  strRatio+=hgammaUPC->GetName();
  TH1D* hratio=(TH1D*)hdata->Clone(strRatio);
  hratio->Sumw2();  
  hratio->Divide(hRatio);
  TH1D *hratio_uncertainty = (TH1D*)hgammaUPC->Clone();
  for(int iBin=1; iBin<=hgammaUPC->GetNbinsX(); iBin++){
    double error_bin =  hgammaUPC->GetBinError(iBin);
    double error_ratio = error_bin/hgammaUPC->GetBinContent(iBin);
    hratio_uncertainty->SetBinContent(iBin, 1);
    hratio_uncertainty->SetBinError(iBin, error_ratio);
  } 
  cout << "Test:" << endl;
  hratio_uncertainty->SetFillColorAlpha(kBlue-4, 0.5);
  hratio_uncertainty->SetFillStyle(3144); 
  hratio_uncertainty->SetLineWidth(2);
  

  hratio->SetTitle("");
  hratio->GetXaxis()->SetTitle(xtitle);
  hratio->GetYaxis()->SetTitle("Data/MC");
  hratio->SetLabelSize(0.1, "XYZ");
  hratio->SetLabelFont(42, "XYZ");
  hratio->SetLabelOffset(0.01, "XYZ");
  hratio->SetTitleSize(0.18, "XYZ");
  hratio->GetXaxis()->SetTitleOffset(1.0);
  hratio->GetYaxis()->SetTitleOffset(0.3);
  hratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hratio->GetXaxis()->SetTitleSize(0.18);
   hratio->GetYaxis()->SetTitleSize(0.18);
   hratio->GetYaxis()->SetLabelSize(0.18);
   hratio->GetXaxis()->SetLabelSize(0.18);
   hratio->GetYaxis()->SetNdivisions(505);
  // cout << "hratio bin content: " << hratio->GetBinContent(4) << endl;
//  pad2->cd();
  //pad2->SetGridy(1);
 // hratio->SetLineWidth(2);
  hratio->Draw("pex0");
  hratio_uncertainty->Draw("e2same");
  line->DrawCopy("same");
  
  c1->Update();
  TString cName=c1->GetName();
  cName+=".png";
  //c1->SaveAs("Plots_From_selectEleMuEvents/tmp/"+cName);
  //c1->SaveAs("Test/"+cName);
  TString c2Name=c1->GetName();
  c2Name+=".pdf";
  //c1->SaveAs("Plots_From_selectEleMuEvents/OnlyRemovingDeltaEtaForElectronAndECALMatching/"+c2Name);
  //c1->SaveAs("Plots_From_selectEleMuEvents/RemovingDeltaEtaForElectronAndECALMatching_DeltaRMuonHCAL0.5_OneFourthLumi/"+c2Name);
  c1->SaveAs("Plots_From_selectEleMuEvents/NewRunningHFThreshold90Percent/"+c2Name);
  return c1;
}
  

void make_canvas(TCanvas *& canvas){

  int W = 650;
  int H = 670;
  float T = 0.08;
  float B = 0.2; 
  float L = 0.20;
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
  hist->SetLabelSize(3, "XYZ");
  hist->GetYaxis()->SetTitleSize(18);
  hist->GetXaxis()->SetTitleSize(18);
  hist->GetYaxis()->SetLabelSize(18);
  hist->GetXaxis()->SetLabelSize(18);
 // hist->SetMarkerColor(kcolor);
  hist->SetFillColor(kcolor);
  hist->SetLineColor(kcolor);
  
  //hist->SetLineWidth(2);
  hist->SetMarkerStyle(kstyle);
  hist->SetBinErrorOption(TH1::kPoisson); 
  hist->SetLineWidth(2);

  //hist->Draw("psamex0");
}

void make_hist_ratio(TH1D *& hist, Color_t kcolor, int kstyle){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.08, "XYZ");
  hist->SetLabelSize(0.1, "XYZ");

  // For the axis titles:
  hist->SetTitleFont(42, "XYZ");
  hist->SetTitleSize(0.18, "XYZ");
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleOffset(0.4);
  hist->SetMarkerColor(kcolor);
  hist->SetFillColor(kcolor);
  hist->SetMarkerStyle(kstyle);
   
}

