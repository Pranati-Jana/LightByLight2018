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
#include "ReadeeTree.C"
#include "CMS_lumi.C"
#include "TStyle.h"
#include "TText.h"
#include "TCut.h"
#include "TChain.h"
#include "THStack.h"
#include "TString.h"



int xlo = 1;
int xhi = 2;
int nbin = 8;


void make_canvas(TCanvas *&);
void make_hist(TH1D *&, Color_t , int , int);
void make_canvas_ratio(TCanvas *&);
void make_hist_ratio(TH1D *&, Color_t, int );

TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, double , double, double , double , double, double,const char *,  const char *); 

const int nSample = 2;

const char *Sample[nSample]={"CutFlow_BW_Data_", "CutFlow_BW_mc_photos_sc_dielectron"};
const char *dir = "figures_eta2p2";


 const int nPtbins=7;
 double Ptbin[nPtbins]={2.5,3.5,4.5,5.5,7.5,9.5,22.5};
 const int nPtbin= sizeof(Ptbin)/sizeof(double) - 1;
void BinLogX(TH1* h);


void plot_ggee(){
  

  double ET;
  double PT;
  double ETe;
  double PTe;

  int idx = 0;
    
	
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  string mc_name = "SC";
  TFile *outf= new TFile("ee.root","recreate");  
//  TFile* f[nSample];
  TChain *tree[nSample];

  //// histograms
  TH1D* hmu_pt1[nSample],*hmu_pt2[nSample],*hmumu_pt[nSample], *hmumuScalar_pt[nSample], *hmumu_mass[nSample], *hmumu_aco[nSample],  *hmu_phi1[nSample], *hmu_phi2[nSample], *hmu_eta1[nSample], *hmu_eta2[nSample],  *hmumu_rapidity[nSample], *hmumu_phi[nSample], *hmumu_pz[nSample], *hZDC_XnXn[nSample], *hDeltaR_TrkElectron[nSample], *hDeltaR_TrkMuon[nSample], *hnTracks[nSample];
  TH2D* hmu_PtEta[nSample], *hZDC[nSample], *hDeltaR_HE[nSample];
   TH1D* h_HE_Energy[nSample],*h_DeltaRHEMuon1[nSample], *h_DeltaRHEMuon2[nSample], *hmumu_aco_woNeutralExclusivity[nSample];
   
  for (int i = 0; i < nSample; i++){
    cout << "nSample:" << i  << endl;
   tree[i] = new TChain("output_tree");
 
  tree[i]->Add(Form("/Users/pranatijana/Documents/GammaGamma_to_Tautau/Workspace/plottingMacros/RootFile_From_selectEleMuEvents/NEE_SF_ggToee/%s_HF+9.1HF-8.8_Eta3.root",Sample[i])); // NEE HF only
 
   hmu_pt1[i] = new TH1D(Form("hmu_pt1%s", Sample[i]),"",23,2.5,25.5);
   hmu_eta1[i] = new TH1D(Form("hmu_eta1%s", Sample[i]),"",13,-2.4,2.4);
   hmu_phi1[i] = new TH1D(Form("hmu_phi1%s", Sample[i]),"",13,-3.5,3.5);

   hmu_pt2[i] = new TH1D(Form("hmu_pt2%s", Sample[i]),"",23,2.5,25.5);
   hmu_eta2[i] = new TH1D(Form("hmu_eta2%s", Sample[i]),"",13,-2.4,2.4);
   hmu_phi2[i] = new TH1D(Form("hmu_phi2%s", Sample[i]),"",13,-3.5,3.5);

   hmumu_pt[i] = new TH1D(Form("hmumu_pt%s", Sample[i]),"",14,0.0,14.0);
   hmumuScalar_pt[i] = new TH1D(Form("hmumuScalar_pt%s", Sample[i]),"",22,0.0,22.0);
   hmumu_mass[i] = new TH1D(Form("hmumu_mass%s", Sample[i]),"",52,8,60);
   hmumu_aco[i] = new TH1D(Form("hmumu_aco%s", Sample[i]),"",50,0.0,0.01);
   hmumu_rapidity[i] = new TH1D(Form("hmumu_rapidity%s", Sample[i]),"",11,-2.5,2.5);
   hZDC[i]      = new TH2D(Form("hZDC%s", Sample[i]), "",87,0,130000,43,0,60000);
   h_HE_Energy[i] = new TH1D(Form("h_HE_Energy%s", Sample[i]),"",10,0,10);
   
   hmumu_aco_woNeutralExclusivity[i] = new TH1D(Form("hmumu_aco_woNeutralExclusivity%s", Sample[i]),"",50,0.0,0.01);

   cout << "nSample:" << i  << endl;
   cout << "file " << Sample[i] << ":" << tree[i]->GetEntries()  << endl;
   ReadeeTree  treeR(tree[i]);
   treeR.fChain->SetBranchStatus("*",1);
   if (treeR.fChain == 0) return;
   Long64_t nentries = treeR.fChain->GetEntriesFast();
   cout << tree[i]->GetName() << "    " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
  
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry_evt = treeR.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = treeR.fChain->GetEntry(jentry);   nbytes += nb;

      std::vector<double> wt = {1, 1};
      if(abs(treeR.eleEta_1)  > 2.2) continue;
       if(abs(treeR.eleEta_2) > 2.2) continue;
       if(treeR.elePt_1 < 2.0) continue;
       if(treeR.elePt_2 < 2.0) continue;
     if(treeR.ok_chexcl_extrk!=1) continue;
     if(treeR.vSum_M< 5) continue;
     if(treeR.vSum_Pt>1) continue;
     if(treeR.ele_acop>0.01) continue;
     if(treeR.ok_HFNeeExcl!=1 ) continue;
      
    // // if((treeR.ok_zdcexcl_3n_pos!=1) && (treeR.ok_zdcexcl_3n_neg!=1) ) continue;
      if(i==0){
       if(treeR.zdc_energy_pos > 7000 &&  treeR.zdc_energy_neg > 7000 ) continue;
             
       // if((treeR.ok_zdcexcl_3n_pos!=1) && (treeR.ok_zdcexcl_3n_neg!=1) ) continue;
      }
      
     hmumu_aco_woNeutralExclusivity[i]->Fill(treeR.ele_acop,wt[i]);
   
     if(treeR.ok_neuexcl!=1) continue; 
     
          
      hmu_pt1[i]->Fill(treeR.elePt_1,wt[i]);
      hmu_eta1[i]->Fill(treeR.eleEta_1,wt[i]);
      hmu_phi1[i]->Fill(treeR.elePhi_1,wt[i]);
     
      hmu_pt2[i]->Fill(treeR.elePt_2,wt[i]);
      hmu_eta2[i]->Fill(treeR.eleEta_2,wt[i]);
      hmu_phi2[i]->Fill(treeR.elePhi_2,wt[i]);

      hmumu_pt[i]->Fill(treeR.vSum_Pt,wt[i]);
    
      hmumu_mass[i]->Fill(treeR.vSum_M,wt[i]);
      hmumu_aco[i]->Fill(treeR.ele_acop,wt[i]);
    
      h_HE_Energy[i]->Fill(treeR.tower_HFp,wt[i]);    
      
   }	   

  } 

 for (int i = 0; i < 2; i++){    
    hmu_pt1[i]->Add(hmu_pt2[i]);
    hmu_eta1[i]->Add(hmu_eta2[i]);
    hmu_phi1[i]->Add(hmu_phi2[i]);
   
  }


 outf->cd();
 outf->Write();
   cout << "test1:" << endl; 

    float Norm = 1;
    hmu_pt1[1]->Scale(Norm);
    hmumu_pt[1]->Scale(Norm);

    hmumu_mass[1]->Scale(Norm);
    hmu_eta1[1]->Scale(Norm);
    hmu_phi1[1]->Scale(Norm);
    h_HE_Energy[1]->Scale(Norm);
 /*
 TCanvas* cc1 = new TCanvas("Muon_pt_LeadingHE","Muon pT",254,411,550,592); 
 cc1->SetLeftMargin( 0.20 );
 cc1->SetBottomMargin( 0.15 );
 //PlotHistsAndRatio(cc1, hmu_pt1[0], hmu_pt1[1],  2.5,25.5,0.01,70000, -1, 1,"  p^{#mu}_{T}  [GeV]");//35
  
  TCanvas* cc3 = new TCanvas("MuMu_VectorSum_pt_LeadingHE","MuMu Vector pT",254,411,550,592);
  cc3->SetLeftMargin( 0.20 );
  cc3->SetBottomMargin( 0.15 );
  
   TCanvas* cc5 = new TCanvas("MuMu_mass","MuMu mass",254,411,550,592);
   cc5->SetLeftMargin( 0.20 );
   cc5->SetBottomMargin( 0.15 );
 // PlotHistsAndRatio(cc5, hmumu_mass[0], hmumu_mass[1], 0.1,60.0,0.5,15000,-1,1, "","M_{#mu#mu} [GeV]");
  
  TCanvas* cc6 = new TCanvas("MuMu_acoplanarity_LeadingHE","MuMu acoplanarity",254,411,550,592);
  cc6->SetLeftMargin( 0.20 );
  cc6->SetBottomMargin( 0.15 );
  //PlotHistsAndRatio(cc6, hmumu_aco[0], hmumu_aco[1], 0,0.25,0.01,20000,-1,1,"Acoplanarity_{#phi}(#mu, #mu)");//15
   TCanvas* cc7 = new TCanvas("Muon_eta_LeadingHE","Muon Eta",254,411,550,592);
   cc7->SetLeftMargin( 0.20 );
   cc7->SetBottomMargin( 0.15 );
 // PlotHistsAndRatio(cc7, hmu_eta1[0], hmu_eta1[1], -2.5,2.5,0.01,25000,-1,1," #eta_{#mu}");//18
   TCanvas* cc8 = new TCanvas("Muon_phi_LeadingHE","Muon phi",254,411,550,592);
   cc8->SetLeftMargin( 0.20 );
   cc8->SetBottomMargin( 0.15 );
 // PlotHistsAndRatio(cc8, hmu_phi1[0], hmu_phi1[1],  -3.5,3.5,0.01,25000,-1,1,"#phi_{#mu} ");//20

  TCanvas* cc9 = new TCanvas("HE_Energy_LeadingHE_All","HE Energy",254,411,550,592); 
 cc9->SetLeftMargin( 0.20 );
 cc9->SetBottomMargin( 0.15 );
 //PlotHistsAndRatio(cc9, h_HE_Energy[0], h_HE_Energy[1],  0.0,10,0.01,40000, -1, 1,"  Leading HE  [GeV]");//35

 TCanvas* cc10 = new TCanvas("DeltaRMuonHE_LeadingHE_all","DeltaRMuonHE",254,411,550,592); 
 cc10->SetLeftMargin( 0.20 );
 cc10->SetBottomMargin( 0.15 );
 //PlotHistsAndRatio(cc10, h_DeltaRHEMuon1[0], h_DeltaRHEMuon1[1],  0,6,0.01,20000, -1, 1,"  #DeltaR(#mu, Leading HE)");//35 
*/
// For NEE SF

  TCanvas* cc11 = new TCanvas("ee_acoplanarity_neuexclNoHF_Data","ee acoplanarity",254,411,550,592);
  cc11->SetLeftMargin( 0.20 );
  cc11->SetBottomMargin( 0.15 );
  PlotHistsAndRatio(cc11, hmumu_aco[0], hmumu_aco_woNeutralExclusivity[0], 0,0.25,10,20000,-1,1,"", "Acoplanarity_{#phi}(#mu, #mu)");
 
  TCanvas* cc12 = new TCanvas("ee_acoplanarity_neuexclNoHF_MC","ee acoplanarity",254,411,550,592);
  cc12->SetLeftMargin( 0.20 );
  cc12->SetBottomMargin( 0.15 );
 //PlotHistsAndRatio(cc12, hmumu_aco[1], hmumu_aco_woNeutralExclusivity[1], 0,0.25,30,50000,-1,1,"", "Acoplanarity_{#phi}(e, e^{-})");
  
  TH1D* hratio_mc=(TH1D*)hmumu_aco[1]->Clone();
  hratio_mc->Divide(hmumu_aco_woNeutralExclusivity[1]);

  TH1D* hratio_data=(TH1D*)hmumu_aco[0]->Clone();
  hratio_data->Divide(hmumu_aco_woNeutralExclusivity[0]);
 
  TH1D* hratio_data_mc = (TH1D*)hratio_data->Clone();
  hratio_data_mc->Divide(hratio_mc);

  TCanvas* cc13 = new TCanvas("ee_acoplanarity_neuexclNoHF_DataMC_","ee aco",254,411,550,592);
  cc13->SetLeftMargin( 0.20 );
  cc13->SetBottomMargin( 0.15 );
 // PlotHistsAndRatio(cc13, hratio_data, hratio_mc, 0,0.25,0.5,1.8,-1,1,"", "Acoplanarity_{#phi}(e,e)");
  
 
  for (int i = 0; i < nSample; i++){
  cout << " Events: " << Sample[i] <<  " :" << hmumu_aco[i]->Integral(1,40) << endl;
  cout << " Events from pt: " << Sample[i] <<  " :" << hmu_pt1[i]->Integral(1,40) << endl;
  cout << " Events from Mass: " << Sample[i] <<  " :" << hmumu_mass[i]->Integral(1,40) << endl;
 // cout << "Bin content: " << Sample[i] << " : " <<  hmu_pt1[i]->GetBinContent(4) << endl;
  cout << " Events MC w: " << Sample[i] <<  " :" << hmumu_aco[0]->Integral(1,50) << endl;
  cout << " Events MC wo: " << Sample[i] <<  " :" << hmumu_aco_woNeutralExclusivity[0]->Integral(1,50) << endl;
 }

 outf->Close();

}

TCanvas* PlotHistsAndRatio(TCanvas* c1, TH1D* hdata, TH1D* hgammaUPC, double hxmin, double hxmax, double hymin, double hymax, double rmax, double rmin, const char *xtitle, const char *ytitle){
  
  // float T = 0.08;
  // float B = 0.14; 
  // float L = 0.34;
  // float R = 0.14;
  //  float T = 0.08;
  float B = 0.14; 
  float L = 0.50;
  float R = 0.04;
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.019);
  pad1->SetLeftMargin(0.14);
  pad1->SetRightMargin(0.05);
  pad1->Draw();
  pad1->cd(); 
  gPad->SetLogy();

  THStack *hs = new THStack("hs"," ");
 

  hs->SetMaximum(hymax);
  hs->SetMinimum(hymin);

  make_hist(hgammaUPC, kRed-7, 20, 2);
  hs->Add(hgammaUPC);
  hs->Draw("e1");
  hs->Draw("e1same");
  hs->GetXaxis()->SetTitle(ytitle);
  hs->GetYaxis()->SetTitle("Events / 0.0002");
  // hs->GetYaxis()->SetTitle("Efficiency");
  hs->GetXaxis()->SetTitleSize(0.06);
   hs->GetYaxis()->SetTitleSize(0.05);
   hs->GetYaxis()->SetLabelSize(0.05);
   hs->GetXaxis()->SetLabelSize(0.05);
   hs->GetXaxis()->SetNdivisions(509);
   hs->GetXaxis()->SetLabelSize(0);

   
  make_hist(hdata, kBlack, 20, 2);
  hdata->SetLineColor(kBlack);
  hdata->Draw("e1samex0");
 
  ///
  //TLegend *leg2=new TLegend(0.55,0.10,0.85,0.21);
  TLegend *leg2=new TLegend(0.20,0.60,0.90,0.91);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(15); 
  leg2->AddEntry(hdata,"Data with neutral exclusivity (No HF) ","EP");
  leg2->AddEntry(hgammaUPC,"Data without neutral exclusivity (No HF) ","EP");
   // leg2->AddEntry(hdata,"Data efficiency   ","EP");
     //leg2->AddEntry(hgammaUPC,"MC efficiency ","EP");

  leg2->Draw();
  

       TLatex *mark = new TLatex();
       mark->SetNDC(true);
       double startY = 0.92;
       mark->SetTextFont(42);
       mark->SetTextSize(0.04);
       mark->DrawLatex(0.15,startY,"#bf{CMS} #it{Private work}");
       mark->DrawLatex(0.50,startY,"#scale[0.8]{PbPb, (#sqrt{s_{NN}} = 5.02 TeV)}");
       mark->Draw();
  
  //c1->Update();
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
  //hRatio->Add(hmumuFSR);

  TString strRatio=hgammaUPC->GetName();
  strRatio+="_over_";
  strRatio+=hgammaUPC->GetName();
  TH1D* hratio=(TH1D*)hdata->Clone(strRatio);
  hratio->Sumw2();  
  
  hratio->Divide(hRatio);
  hratio->SetTitle("");
  hratio->GetXaxis()->SetTitle(xtitle);
 // hratio->GetYaxis()->SetTitle("Data/MC");
  hratio->GetYaxis()->SetTitle("Efficiency");
  hratio->GetXaxis()->SetTitle("Acoplanarity_{#phi}(e, e^{-})");
  hratio->SetLabelSize(0.1, "XYZ");
  hratio->SetLabelFont(42, "XYZ");
  hratio->SetLabelOffset(0.01, "XYZ");
  hratio->SetTitleSize(0.18, "XYZ");
  hratio->GetXaxis()->SetTitleOffset(1.0);
  hratio->GetYaxis()->SetTitleOffset(0.3);
  hratio->GetYaxis()->SetRangeUser(0.5,1.2);
  hratio->GetXaxis()->SetTitleSize(0.18);
   hratio->GetYaxis()->SetTitleSize(0.18);
   hratio->GetYaxis()->SetLabelSize(0.18);
   hratio->GetXaxis()->SetLabelSize(0.18);
   hratio->GetYaxis()->SetNdivisions(505);
   hratio->GetXaxis()->SetNdivisions(505);
   cout << "hratio bin content: " << hratio->GetBinContent(4) << endl;
//  pad2->cd();
  //pad2->SetGridy(1);
 // hratio->SetLineWidth(2);
 hratio->Fit("pol0");
  hratio->Draw("pex0");
  line->DrawCopy("same");
  c1->Update();
  TString cName=c1->GetName();
  cName+=".png";
 
  TString c2Name=c1->GetName();
  c2Name+=".pdf";

  //c1->SaveAs("Plots_From_selectEleMuEvents/ForTriggerHFveto_SF_Study_ggmumu/"+c2Name);
  c1->SaveAs("Plots_From_selectEleMuEvents/NEE_Dielectron/"+cName);
  return c1;
}
  

void make_canvas(TCanvas *& canvas){

  int W = 550;
  int H = 670;
  float T = 0.08;
  float B = 0.2; 
  float L = 0.50;
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
  hist->SetLabelOffset(0.01, "XYZ");
  hist->SetLabelSize(3, "XYZ");
  hist->GetYaxis()->SetTitleSize(18);
  hist->GetXaxis()->SetTitleSize(18);
  hist->GetYaxis()->SetLabelSize(15);
  hist->GetXaxis()->SetLabelSize(15);
  hist->SetMarkerColor(kcolor);
  //hist->SetFillColor(kcolor);
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

