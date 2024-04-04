#include <iostream>
#include <fstream>
#include <sstream>
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double syst_uncertainty = 0.068;
double syst_uncertainty_mc = 0.068;

void make_hist(TH1D *&, string, string, int, float, Color_t, int, int) ;
void make_hist_syst(TH1D *&, string, string, int, float, Color_t, int, int) ;
TH1D* get_stat_uncertainty_hist(TH1D *input_hist, double uncertainty = syst_uncertainty);
TH1D* get_stat_uncertainty_hist_pt(TH1D *input_hist);
TH1D* get_stat_uncertainty_hist_mass(TH1D *input_hist);
TH1D* get_stat_uncertainty_hist_rap(TH1D *input_hist);
TH1D* get_stat_uncertainty_hist_cos(TH1D *input_hist);
TH1D* get_stat_uncertainty_hist_ratio(TH1D *input_hist);
const char *dir  = "fig_AN/fig_AN_Paper";
void printBinCont(TH1D *hist);
void prepare_canvas(TCanvas *canvas);
TH1D* get_ratio(TH1D *hist_num, TH1D* hist_den);


int mc_gen_color = kCyan+2;
int data_unfold_color = kBlack;
int mc_reco_color = kBlack;
int data_reco_color = kOrange+1;
int gUPC_color = kBlue+2;
int SL_color = kRed-7;//kBlue+2;

void draw_QED_gUPC(int method =1, string mc_name = "SC")
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  //gStyle->SetPadLeftMargin(0.01);
  gStyle->SetPadRightMargin(0.005);
  char MethID1[100];
  if(method==1)
  {
  
  }
  if(method==2)
  {
 
  }
  if(method==3)
  {
  
  }
  
  TFile *f1 = new TFile(Form("unfolding_histograms_Bayes_QEDSC%s_gUPC.root", mc_name.c_str()),"r");
  TFile *f2 = new TFile(Form("unfolding_histograms_Bayes_QEDSLSC_gUPC.root"),"r");
  TFile *f3 = new TFile(Form("Comparison/gammaUPC_PY8FSR_predictions_QED.root"),"r");
   
  TH1D *hMuMu_pt_SL                = (TH1D*)f2->Get("hGenPt_xSec");
  TH1D *hMuMu_rapidity_SL          = (TH1D*)f2->Get("hGenRap_xSec");
  TH1D *hMuMu_invmass_SL           = (TH1D*)f2->Get("hGenInvmass_xSec");
  TH1D *hMuMu_costhetastar_SL      = (TH1D*)f2->Get("hGenCosthetastar_xSec");

  TH1D *hMuMu_pt_gUPC          = (TH1D*)f3->Get("hgammaUPC_QED_PtPair");
  TH1D *hMuMu_rapidity_gUPC    = (TH1D*)f3->Get("hgammaUPC_QED_Rap");
  TH1D *hMuMu_invmass_gUPC     = (TH1D*)f3->Get("hgammaUPC_QED_Invmass");
  TH1D *hMuMu_costhetastar_gUPC     = (TH1D*)f3->Get("hgammaUPC_QED_CosThetaStar");

  hMuMu_pt_gUPC->Scale(0.000001);
  hMuMu_rapidity_gUPC->Scale(0.000001);
  hMuMu_invmass_gUPC->Scale(0.000001);
  hMuMu_costhetastar_gUPC->Scale(0.000001); 

  hMuMu_pt_SL->Scale(0.001);
  hMuMu_rapidity_SL->Scale(0.001);
  hMuMu_invmass_SL->Scale(0.001);
  hMuMu_costhetastar_SL->Scale(0.001);


  TFile *outf= new TFile(Form("Test.root"),"recreate"); 
  if(!f1){
    cout<<"File not found: "<<Form("unfolding_histograms_%s.root", MethID1)<<endl;
    exit(0);
  }
  cout<<"   "<<f1->GetName()<<endl;
  cout<<"   "<<f2->GetName()<<endl;
  cout<<"   "<<f3->GetName()<<endl;
  
  TH1D *hMuMu_pt_gen          = (TH1D*)f1->Get("hGenPt_xSec");
  TH1D *hMuMu_pt_recomc       = (TH1D*)f1->Get("hRecoMCPt_xSec");
  TH1D *hMuMu_pt_recodata     = (TH1D*)f1->Get("hRecoDataPt_xSec");
  TH1D *hMuMu_pt_unfodata     = (TH1D*)f1->Get("hUnfoDataPt_xSec");
  TH1D *hMuMu_pt_unfomc       = (TH1D*)f1->Get("hUnfoMCPt_xSec");

  hMuMu_pt_gen->Scale(0.001);
  hMuMu_pt_recomc->Scale(0.001);
  hMuMu_pt_recodata->Scale(0.001);
  hMuMu_pt_unfodata->Scale(0.001);
  double renorm1 = 271.5/hMuMu_pt_unfodata->Integral("width");
  hMuMu_pt_unfodata->Scale(renorm1);
  hMuMu_pt_unfomc->Scale(0.001);
 
 
  TH1D *hMuMu_rapidity_gen          = (TH1D*)f1->Get("hGenRap_xSec");
  TH1D *hMuMu_rapidity_recomc       = (TH1D*)f1->Get("hRecoMCRap_xSec");
  TH1D *hMuMu_rapidity_recodata     = (TH1D*)f1->Get("hRecoDataRap_xSec");
  TH1D *hMuMu_rapidity_unfodata     = (TH1D*)f1->Get("hUnfoDataRap_Invert_xSec");
  TH1D *hMuMu_rapidity_unfomc       = (TH1D*)f1->Get("hUnfoMCRap_Invert_xSec");
  
  hMuMu_rapidity_gen->Scale(0.001);
  hMuMu_rapidity_recomc->Scale(0.001);
  hMuMu_rapidity_recodata->Scale(0.001);
  hMuMu_rapidity_unfodata->Scale(0.001);
  double renorm2 = 271.5/hMuMu_rapidity_unfodata->Integral("width");
  hMuMu_rapidity_unfodata->Scale(renorm2);
  cout << "rap cross section:" << hMuMu_rapidity_unfodata->Integral("width") << endl;
  hMuMu_rapidity_unfomc->Scale(0.001);
  
  TH1D *hMuMu_invmass_gen          = (TH1D*)f1->Get("hGenInvmass_xSec");
  TH1D *hMuMu_invmass_recomc       = (TH1D*)f1->Get("hRecoMCInvmass_xSec");
  TH1D *hMuMu_invmass_recodata     = (TH1D*)f1->Get("hRecoDataInvmass_xSec");
  TH1D *hMuMu_invmass_unfodata     = (TH1D*)f1->Get("hUnfoDataInvmass_Invert_xSec");
  TH1D *hMuMu_invmass_unfomc       = (TH1D*)f1->Get("hUnfoMCInvmass_Invert_xSec");

  hMuMu_invmass_gen->Scale(0.001);
  hMuMu_invmass_recomc->Scale(0.001);
  hMuMu_invmass_recodata->Scale(0.001);
  hMuMu_invmass_unfodata->Scale(0.001);
  double renorm3 = 271.5/hMuMu_invmass_unfodata->Integral("width");
  hMuMu_invmass_unfodata->Scale(renorm3);
  hMuMu_invmass_unfomc->Scale(0.001);
 
  TH1D *hMuMu_costhetastar_gen          = (TH1D*)f1->Get("hGenCosthetastar_xSec");
  TH1D *hMuMu_costhetastar_recomc       = (TH1D*)f1->Get("hRecoMCCosthetastar_xSec");
  TH1D *hMuMu_costhetastar_recodata     = (TH1D*)f1->Get("hRecoDataCosthetastar_xSec");
  TH1D *hMuMu_costhetastar_unfodata     = (TH1D*)f1->Get("hUnfoDataCosthetastar_Invert_xSec");
  TH1D *hMuMu_costhetastar_unfomc       = (TH1D*)f1->Get("hUnfoMCCosthetastar_Invert_xSec");

  hMuMu_costhetastar_gen->Scale(0.001);
  hMuMu_costhetastar_recomc->Scale(0.001); 
  hMuMu_costhetastar_recodata->Scale(0.001);
  hMuMu_costhetastar_unfodata->Scale(0.001);
  double renorm4 = 271.5/hMuMu_costhetastar_unfodata->Integral("width");
  hMuMu_costhetastar_unfodata->Scale(renorm4);
  hMuMu_costhetastar_unfomc->Scale(0.001);
  cout<<"ok"<<endl;
  

  TLegend *leg1=new TLegend(0.30,0.72,0.90,0.9);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(43);
  leg1->SetTextSize(16);
  leg1->AddEntry(hMuMu_pt_unfomc,"Unfolded MC (Bayes, 1 iteration)","pl");
  leg1->AddEntry(hMuMu_pt_gen, "Gen-level MC (Superchic+Photos)","pl");
  leg1->AddEntry(hMuMu_pt_recomc,"Reconstructed MC (Superchic+Photos)","pl");
  

 
  TLegend *leg3=new TLegend(0.20,0.6,0.91,0.9);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextFont(43);
  leg3->SetTextSize(18);
  leg3->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  leg3->AddEntry(hMuMu_rapidity_unfodata,"Data","EP");
  leg3->AddEntry(hMuMu_rapidity_gUPC,  "gamma-UPC@NLO+PY8","pl");
  leg3->AddEntry(hMuMu_rapidity_gen,  "SUPERCHIC+PHOTOS++","pl");
  leg3->AddEntry(hMuMu_rapidity_SL,  "STARLIGHT","pl");
  

  TLegend *leg6=new TLegend(0.30,0.72,0.90,0.9);
  leg6->SetFillColor(0);
  leg6->SetBorderSize(0);
  leg6->SetFillStyle(0);
  leg6->SetTextFont(43);
  leg6->SetTextSize(16);
  leg6->AddEntry(hMuMu_rapidity_unfomc,"Unfolded MC (Matrix inversion)","pl");
  leg6->AddEntry(hMuMu_rapidity_gen, "Gen-level MC (Superchic+Photos)","pl");
  leg6->AddEntry(hMuMu_rapidity_recomc,"Reconstructed MC (Superchic+Photos)","pl");
  
  TLegend *leg4=new TLegend(0.20,0.6,0.91,0.9);
  leg4->SetFillColor(0);
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetTextFont(43);
  leg4->SetTextSize(18);
  leg4->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  leg4->AddEntry(hMuMu_invmass_unfodata,"Data","EP");
  leg4->AddEntry(hMuMu_invmass_gUPC,  "gamma-UPC@NLO+PY8","pl");
  leg4->AddEntry(hMuMu_invmass_gen,  "SUPERCHIC+PHOTOS++","pl");
  leg4->AddEntry(hMuMu_invmass_SL,  "STARLIGHT","pl");
  
   
  TLegend *leg5=new TLegend(0.30,0.72,0.90,0.9);
  leg5->SetFillColor(0);
  leg5->SetBorderSize(0);
  leg5->SetFillStyle(0);
  leg5->SetTextFont(43);
  leg5->SetTextSize(18);
  leg5->AddEntry(hMuMu_invmass_unfomc,"Unfolded MC (Matrix inversion)","pl");
  leg5->AddEntry(hMuMu_invmass_gen, "Gen-level MC (Superchic+Photos)","pl");
  leg5->AddEntry(hMuMu_invmass_recomc,"Reconstructed MC (Superchic+Photos)","pl");
  
   TLatex *mark = new TLatex();
       mark->SetNDC(true);
       double startY = 0.90;
       mark->SetTextFont(42);
       mark->SetTextSize(0.035);
       mark->DrawLatex(0.15,startY,"#bf{CMS} #it{Preliminary}");
       mark->DrawLatex(0.55,startY,"#scale[0.8]{1647.23 #mub^{-1} (PbPb @ 5.02 TeV)}");
       mark->Draw();

  auto line = new TF1("line", "1", -100, 100);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);


  int W = 600;
  int H = 500;
  
  // pT

 TCanvas* c2 = new TCanvas("c2","Delectron pt",50,50,W,H);
  prepare_canvas(c2);
  TPad *pad7 = new TPad("pad7", "pad7", 0, 0.3, 1, 1.0);
  pad7->SetBottomMargin(0.012);
  pad7->SetRightMargin(0.08); 
  pad7->SetTopMargin(0.08);
  pad7->Draw();
  pad7->cd();
  gPad->SetLogy();
  hMuMu_pt_gen->SetMaximum(10000);
  hMuMu_pt_gen->SetMinimum(0.02);
  
 string x_title = "p_{T}^{ee} (GeV) ";
 string  y_title = "d#sigma^{ee}/dp_{T}^{ee} (#mub / GeV)";
 
  make_hist(hMuMu_pt_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  make_hist(hMuMu_pt_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hMuMu_pt_gUPC,x_title,y_title,2,1.0,gUPC_color,1,2);
  make_hist(hMuMu_pt_SL,x_title,y_title,2,1.0,SL_color,1,2);
  auto hMuMu_pt_unfodata_syst = get_stat_uncertainty_hist_pt(hMuMu_pt_unfodata);
  auto hMuMu_pt_gen_syst = get_stat_uncertainty_hist_pt(hMuMu_pt_gen);
  auto hMuMu_pt_gUPC_syst = get_stat_uncertainty_hist_pt(hMuMu_pt_gUPC);
  auto hMuMu_pt_SL_syst = get_stat_uncertainty_hist_pt(hMuMu_pt_SL );
  hMuMu_pt_gen->GetXaxis()->SetLabelSize(0);
  hMuMu_pt_gen->Draw("histx0");
  //hMuMu_rpt_gen_syst->Draw("e2same");
 // hMuMu_pt_gUPC_syst->Draw("e2same");
  hMuMu_pt_gUPC->Draw("histsamex0");
  hMuMu_pt_gUPC->Draw("psame");
  //hMuMu_pt_SL_syst->Draw("e2same");
  hMuMu_pt_SL->Draw("histsamex0");
  hMuMu_pt_SL->Draw("psame");
  hMuMu_pt_unfodata_syst->Draw("e2same");
  hMuMu_pt_unfodata->Draw("psameex0");
  /**************************Legend************************************************/
  TH1D *hMuMu_pt_data_ratio = get_ratio(hMuMu_pt_unfodata, hMuMu_pt_gen);
  TH1D *hMuMu_pt_data_gUPC_ratio = get_ratio(hMuMu_pt_unfodata, hMuMu_pt_gUPC);
  TH1D *hMuMu_pt_data_SL_ratio = get_ratio(hMuMu_pt_unfodata, hMuMu_pt_SL);

   make_hist(hMuMu_pt_data_gUPC_ratio,x_title,y_title,20,1.0,gUPC_color,1,2);
   make_hist(hMuMu_pt_data_SL_ratio,x_title,y_title,22,1.0,SL_color,1,2);
   
   hMuMu_pt_data_ratio->SetLineColor(mc_gen_color);
   hMuMu_pt_data_ratio->SetMarkerColor(mc_gen_color);

  auto hMuMu_pt_data_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_pt_unfodata_syst);
  auto hMuMu_pt_data_gUPC_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_pt_gUPC);
  auto hMuMu_pt_data_SL_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_pt_SL);

  make_hist(hMuMu_pt_data_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  make_hist(hMuMu_pt_data_gUPC_ratio_syst,x_title,y_title,0,1.0,gUPC_color,2,2);
  make_hist(hMuMu_pt_data_SL_ratio_syst,x_title,y_title,0,1.0,SL_color,2,2);

  TLegend *leg2=new TLegend(0.30,0.6,0.91,0.9);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(16);
  leg2->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  //leg2->AddEntry("","#gamma#gamma#rigtharrow e^{+}e^{-} (PbPb 5.02 TeV)");
 // leg2->SetEntrySeparation(1.5); 
  leg2->AddEntry(hMuMu_pt_unfodata,"Data","EP");
  leg2->AddEntry(hMuMu_pt_data_gUPC_ratio,  "gamma-UPC/MG5 + FSR (PY8)","pl");
  leg2->AddEntry(hMuMu_pt_data_ratio,  "SUPERCHIC 3.03 + FSR (PHOTOS++)","pl");
  leg2->AddEntry(hMuMu_pt_data_SL_ratio,  "STARLIGHT 3.13","pl");
  leg2->Draw();

  c2->cd();

  TLatex *mark8 = new TLatex();
  mark8->SetTextSize(0.037); 
  mark8->SetTextFont(42);
  
   mark8->DrawLatex(0.16, 0.96, "#bf{CMS} #it{Preliminary}");
   mark8->DrawLatex(0.55, 0.96, "#scale[0.8]{PbPb, 1.65 nb^{-1} (#sqrt{s_{NN}} = 5.02 TeV)}");
  TPad *pad8 = new TPad("pad8", "pad8", 0, 0.06, 1, 0.28);
  pad8->SetTopMargin(0.04);
  pad8->SetBottomMargin(0.45);
  pad8->Draw();
  pad8->cd();
  pad8->SetRightMargin(0.08);
  
  

  hMuMu_pt_data_ratio->GetXaxis()->SetTitle("p_{T}^{ee} (GeV)");
  hMuMu_pt_data_ratio->GetYaxis()->SetTitle("Data/MC");
  
  hMuMu_pt_data_ratio->SetMaximum(1.45);
  hMuMu_pt_data_ratio->SetMinimum(0.65);
  hMuMu_pt_data_ratio->SetMarkerStyle(21);
  hMuMu_pt_data_ratio->GetXaxis()->SetLabelOffset(0.03);
  hMuMu_pt_data_ratio->Draw("pex0"); 
  hMuMu_pt_data_gUPC_ratio->Draw("pex0same"); 
  hMuMu_pt_data_SL_ratio->Draw("pex0same"); 
  hMuMu_pt_data_ratio_syst->Draw("e2same");

  line->DrawCopy("same");

  c2->Print(Form("./%s/data_pt_QED_gUPCPythia8_Comparison_30thMarch.pdf",dir));
  c2->Print(Form("./%s/data_pt_QED_gUPCPythia8_Comparison_ForAN_Paper_30thMarch.C",dir));
  //Data Rapidity

  TCanvas* c4 = new TCanvas("c4","Dielectron rapidity",50,50,W,H);
  prepare_canvas(c4);
  TPad *pad9 = new TPad("pad9", "pad9", 0, 0.3, 1, 1.0);
  pad9->SetBottomMargin(0.02);
  pad9->SetRightMargin(0.08); 
  pad9->SetTopMargin(0.08);
  pad9->Draw();
  pad9->cd();
 // gPad->SetLogy();
  hMuMu_rapidity_gen->SetMaximum(185);//500
  hMuMu_rapidity_gen->SetMinimum(0);
  
  x_title = "y^{ee} ";
  y_title = "d#sigma^{ee}/dy^{ee} (#mub)";
 
  make_hist(hMuMu_rapidity_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  make_hist(hMuMu_rapidity_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hMuMu_rapidity_gUPC,x_title,y_title,2,1.0,gUPC_color,1,2);
  make_hist(hMuMu_rapidity_SL,x_title,y_title,2,1.0,SL_color,1,2);
  auto hMuMu_rapidity_unfodata_syst = get_stat_uncertainty_hist_rap(hMuMu_rapidity_unfodata);
  auto hMuMu_rapidity_gen_syst = get_stat_uncertainty_hist_rap(hMuMu_rapidity_gen);
  auto hMuMu_rapidity_gUPC_syst = get_stat_uncertainty_hist_rap(hMuMu_rapidity_gUPC );
  auto hMuMu_rapidity_SL_syst = get_stat_uncertainty_hist_rap(hMuMu_rapidity_SL );
  hMuMu_rapidity_gen->GetXaxis()->SetLabelSize(0);
  hMuMu_rapidity_gen->Draw("histx0");
  hMuMu_rapidity_gen->Draw("psamex0");
  //hMuMu_rapidity_gen_syst->Draw("e2same");
 // hMuMu_rapidity_gUPC_syst->Draw("e2same");
  hMuMu_rapidity_gUPC->Draw("histsamex0");
  hMuMu_rapidity_gUPC->Draw("psamex0");
  //hMuMu_rapidity_SL_syst->Draw("e2same");
  hMuMu_rapidity_SL->Draw("histsamex0");
  hMuMu_rapidity_SL->Draw("psamex0");
  hMuMu_rapidity_unfodata_syst->Draw("e2same");
  hMuMu_rapidity_unfodata->Draw("psameex0");
  leg2->Draw();
  
  c4->cd();
  TLatex *mark9 = new TLatex();
  mark9->SetTextSize(0.037); 
  mark9->SetTextFont(42);
  
   mark9->DrawLatex(0.16, 0.96, "#bf{CMS} #it{Preliminary}");
   mark9->DrawLatex(0.55, 0.96, "#scale[0.8]{PbPb, 1.65 nb^{-1} (#sqrt{s_{NN}} = 5.02 TeV)}");
  TPad *pad10 = new TPad("pad10", "pad10", 0, 0.06, 1, 0.28);
  pad10->SetTopMargin(0.04);
  pad10->SetBottomMargin(0.45);
  pad10->Draw();
  pad10->cd();
  pad10->SetRightMargin(0.08);
  TH1D *hMuMu_rapidity_data_ratio = get_ratio(hMuMu_rapidity_unfodata, hMuMu_rapidity_gen);
  TH1D *hMuMu_rapidity_data_gUPC_ratio = get_ratio(hMuMu_rapidity_unfodata, hMuMu_rapidity_gUPC);
  TH1D *hMuMu_rapidity_data_SL_ratio = get_ratio(hMuMu_rapidity_unfodata, hMuMu_rapidity_SL);

   make_hist(hMuMu_rapidity_data_gUPC_ratio,x_title,y_title,20,1.0,gUPC_color,1,2);
   make_hist(hMuMu_rapidity_data_SL_ratio,x_title,y_title,22,1.0,SL_color,1,2);

   hMuMu_rapidity_data_ratio->SetLineColor(mc_gen_color);
   hMuMu_rapidity_data_ratio->SetMarkerColor(mc_gen_color);

  auto hMuMu_rapidity_data_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_rapidity_unfodata_syst);
  auto hMuMu_rapidity_data_gUPC_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_rapidity_gUPC);
  auto hMuMu_rapidity_data_SL_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_rapidity_SL);

  make_hist(hMuMu_rapidity_data_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  make_hist(hMuMu_rapidity_data_gUPC_ratio_syst,x_title,y_title,0,1.0,gUPC_color,2,2);
  make_hist(hMuMu_rapidity_data_SL_ratio_syst,x_title,y_title,0,1.0,SL_color,2,2);

  hMuMu_rapidity_data_ratio->GetXaxis()->SetTitle("y^{ee}");
  hMuMu_rapidity_data_ratio->GetYaxis()->SetTitle("Data/MC");
  
  hMuMu_rapidity_data_ratio->SetMaximum(1.3);
  hMuMu_rapidity_data_ratio->SetMinimum(0.6);
  hMuMu_rapidity_data_ratio->SetMarkerStyle(21);
  hMuMu_rapidity_data_ratio->GetXaxis()->SetLabelOffset(0.03);
  hMuMu_rapidity_data_ratio->Draw("pex0");
  
   hMuMu_rapidity_data_gUPC_ratio->Draw("pex0same"); 
  hMuMu_rapidity_data_SL_ratio->Draw("pex0same"); 
 
  hMuMu_rapidity_data_ratio_syst->Draw("e2same");
  
  line->DrawCopy("same");
 
  c4->Print(Form("./%s/data_rapidity_QED_gUPCPythia8_Comparison_30thMarch.pdf",dir));
  c4->Print(Form("./%s/data_rapidity_QED_gUPCPythia8_Comparison_ForAN_Paper_30thMarch.C",dir));
  ///Data Inv mass
  TCanvas* c6 = new TCanvas("c6","Dielectron mass",50,50,W,H);
  prepare_canvas(c6);
  TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
  pad11->SetBottomMargin(0.01);
  pad11->SetRightMargin(0.08); 
  pad11->SetTopMargin(0.08);
  pad11->Draw();
  pad11->cd();
  gPad->SetLogy();
  // hMuMu_invmass_unfodata->SetMaximum(200000);
  // hMuMu_invmass_unfodata->SetMinimum(1E-1);

  x_title = "m_{ee} (GeV)";
  y_title =   "d#sigma^{ee}/dm^{ee} (#mub / GeV)";
 
  make_hist(hMuMu_invmass_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  make_hist(hMuMu_invmass_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hMuMu_invmass_gUPC,x_title,y_title,2,1.0,gUPC_color,1,2);
  make_hist(hMuMu_invmass_SL,x_title,y_title,2,1.0,SL_color,1,2);
  auto hMuMu_invmass_unfodata_syst = get_stat_uncertainty_hist_mass(hMuMu_invmass_unfodata);
  auto hMuMu_invmass_gen_syst = get_stat_uncertainty_hist_mass(hMuMu_invmass_gen);
  auto hMuMu_invmass_gUPC_syst = get_stat_uncertainty_hist_mass(hMuMu_invmass_gUPC );
  auto hMuMu_invmass_SL_syst = get_stat_uncertainty_hist_mass(hMuMu_invmass_SL );

  //hMuMu_invmass_gen->Draw("histsamex0");
  hMuMu_invmass_gen->Draw("histx0");
  hMuMu_invmass_gen->Draw("px0same");
 // hMuMu_invmass_gen_syst->Draw("e2same");
 // hMuMu_invmass_gUPC_syst->Draw("e2same");
  hMuMu_invmass_gUPC->Draw("histsamex0");
  hMuMu_invmass_gUPC->Draw("psamex0");
 // hMuMu_invmass_SL_syst->Draw("e2same");
  hMuMu_invmass_SL->Draw("histsamex0");
  hMuMu_invmass_SL->Draw("psamex0");
  hMuMu_invmass_unfodata_syst->Draw("e2same");
  hMuMu_invmass_unfodata->Draw("psameex0");
  leg2->Draw();
  
  c6->cd();
  TLatex *mark10 = new TLatex();
  mark10->SetTextSize(0.037); 
  mark10->SetTextFont(42);
  
   mark10->DrawLatex(0.16, 0.96, "#bf{CMS} #it{Preliminary}");
   mark10->DrawLatex(0.55, 0.96, "#scale[0.8]{PbPb, 1.65 nb^{-1} (#sqrt{s_{NN}} = 5.02 TeV)}");
  TPad *pad12 = new TPad("pad12", "pad12", 0, 0.06, 1, 0.28);
  pad12->SetTopMargin(0.052);
  pad12->SetBottomMargin(0.45);
  pad12->Draw();
  pad12->cd();
  pad12->SetRightMargin(0.08);
  TH1D *hMuMu_invmass_data_ratio = get_ratio(hMuMu_invmass_unfodata, hMuMu_invmass_gen);
  TH1D *hMuMu_invmass_data_gUPC_ratio = get_ratio(hMuMu_invmass_unfodata, hMuMu_invmass_gUPC);
  TH1D *hMuMu_invmass_data_SL_ratio = get_ratio(hMuMu_invmass_unfodata, hMuMu_invmass_SL);

   make_hist(hMuMu_invmass_data_gUPC_ratio,x_title,y_title,20,1.0,gUPC_color,1,2);
   make_hist(hMuMu_invmass_data_SL_ratio,x_title,y_title,22,1.0,SL_color,1,2);

   hMuMu_invmass_data_ratio->SetLineColor(mc_gen_color);
   hMuMu_invmass_data_ratio->SetMarkerColor(mc_gen_color);

  auto hMuMu_invmass_data_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_invmass_unfodata_syst);
  auto hMuMu_invmass_data_gUPC_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_invmass_gUPC);
  auto hMuMu_invmass_data_SL_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_invmass_SL);
  
  make_hist(hMuMu_invmass_data_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  make_hist(hMuMu_invmass_data_gUPC_ratio_syst,x_title,y_title,0,1.0,gUPC_color,2,2);
  make_hist(hMuMu_invmass_data_SL_ratio_syst,x_title,y_title,0,1.0,SL_color,2,2);

  hMuMu_invmass_data_ratio->GetXaxis()->SetTitle("m^{ee} (GeV)");
  hMuMu_invmass_data_ratio->GetYaxis()->SetTitle("Data/MC");
  hMuMu_invmass_data_ratio->SetMarkerStyle(21);
  hMuMu_invmass_data_ratio->SetMaximum(2.0);
  hMuMu_invmass_data_ratio->GetXaxis()->SetLabelOffset(0.03);
  hMuMu_invmass_data_ratio->Draw("pex0");
  hMuMu_invmass_data_gUPC_ratio->Draw("pex0same"); 
  hMuMu_invmass_data_SL_ratio->Draw("pex0same"); 
  
  hMuMu_invmass_data_ratio_syst->Draw("e2same");

 
 // hMuMu_invmass_data_SL_ratio_syst->Draw("e2same");
  //hMuMu_invmass_data_gUPC_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
 
  c6->Print(Form("./%s/data_invmass_QED_gUPCPythia8_30thMarch.pdf",dir));
  c6->Print(Form("./%s/data_invmass_QED_gUPCPythia8_ForAN_Paper_30thMarch.C",dir));

  //costhetastar
   TCanvas* c8 = new TCanvas("c8","Dielectron costhetastar",50,50,W,H);
  prepare_canvas(c8);
  TPad *pad13 = new TPad("pad13", "pad13", 0, 0.3, 1, 1.0);
  pad13->SetBottomMargin(0.0);
  pad13->SetRightMargin(0.08); 
  pad13->SetTopMargin(0.08);
  pad13->Draw();
  pad13->cd();
  //gPad->SetLogy();
  hMuMu_costhetastar_gen->SetMaximum(610);
  hMuMu_costhetastar_gen->SetMinimum(5);
  
  x_title = "|cos#theta*|^{ee} ";
  y_title = "d#sigma^{ee}/d|cos#theta*|^{ee} (#mub)";
 
  make_hist(hMuMu_costhetastar_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  make_hist(hMuMu_costhetastar_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hMuMu_costhetastar_gUPC,x_title,y_title,2,1.0,gUPC_color,1,2);
  make_hist(hMuMu_costhetastar_SL,x_title,y_title,2,1.0,SL_color,1,2);
  auto hMuMu_costhetastar_unfodata_syst = get_stat_uncertainty_hist_cos(hMuMu_costhetastar_unfodata);
  auto hMuMu_costhetastar_gen_syst = get_stat_uncertainty_hist_cos(hMuMu_costhetastar_gen);
  auto hMuMu_costhetastar_gUPC_syst = get_stat_uncertainty_hist_cos(hMuMu_costhetastar_gUPC );
  auto hMuMu_costhetastar_SL_syst = get_stat_uncertainty_hist_cos(hMuMu_costhetastar_SL );

  hMuMu_costhetastar_gen->Draw("histx0");
  hMuMu_costhetastar_gen->Draw("psamex0");
  //hMuMu_costhetastar_gen_syst->Draw("e2same");
 // hMuMu_costhetastar_gUPC_syst->Draw("e2same");
  hMuMu_costhetastar_gUPC->Draw("histsamex0");
  hMuMu_costhetastar_gUPC->Draw("psamex0");
  //hMuMu_costhetastar_SL_syst->Draw("e2same");
  hMuMu_costhetastar_SL->Draw("histsamex0");
  hMuMu_costhetastar_SL->Draw("psamex0");
  hMuMu_costhetastar_unfodata_syst->Draw("e2same");
  hMuMu_costhetastar_unfodata->Draw("psameex0");
  leg2->Draw();
  
  c8->cd();
  TLatex *mark11 = new TLatex();
  mark11->SetTextSize(0.037); 
  mark11->SetTextFont(42);
  
   mark11->DrawLatex(0.16, 0.96, "#bf{CMS} #it{Preliminary}");
   mark11->DrawLatex(0.55, 0.96, "#scale[0.8]{PbPb, 1.65 nb^{-1} (#sqrt{s_{NN}} = 5.02 TeV)}");
  TPad *pad14 = new TPad("pad14", "pad14", 0, 0.06, 1, 0.28);
  pad14->SetTopMargin(0.04);
  pad14->SetBottomMargin(0.47);
  pad14->Draw();
  pad14->cd();
  pad14->SetRightMargin(0.08);
  TH1D *hMuMu_costhetastar_data_ratio = get_ratio(hMuMu_costhetastar_unfodata, hMuMu_costhetastar_gen);
  TH1D *hMuMu_costhetastar_data_gUPC_ratio = get_ratio(hMuMu_costhetastar_unfodata, hMuMu_costhetastar_gUPC);
  TH1D *hMuMu_costhetastar_data_SL_ratio = get_ratio(hMuMu_costhetastar_unfodata, hMuMu_costhetastar_SL);

   make_hist(hMuMu_costhetastar_data_gUPC_ratio,x_title,y_title,20,1.0,gUPC_color,1,2);
   make_hist(hMuMu_costhetastar_data_SL_ratio,x_title,y_title,22,1.0,SL_color,1,2);

   hMuMu_costhetastar_data_ratio->SetLineColor(mc_gen_color);
   hMuMu_costhetastar_data_ratio->SetMarkerColor(mc_gen_color);

  auto hMuMu_costhetastar_data_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_costhetastar_unfodata_syst);
  auto hMuMu_costhetastar_data_gUPC_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_costhetastar_gUPC);
  auto hMuMu_costhetastar_data_SL_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_costhetastar_SL);

  make_hist(hMuMu_costhetastar_data_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  make_hist(hMuMu_costhetastar_data_gUPC_ratio_syst,x_title,y_title,0,1.0,gUPC_color,2,2);
  make_hist(hMuMu_costhetastar_data_SL_ratio_syst,x_title,y_title,0,1.0,SL_color,2,2);

  hMuMu_costhetastar_data_ratio->GetXaxis()->SetTitle("|cos#theta*|^{ee}");
  hMuMu_costhetastar_data_ratio->GetYaxis()->SetTitle("Data/MC");
  hMuMu_costhetastar_data_ratio->SetMarkerStyle(21);
  hMuMu_costhetastar_data_ratio->SetMaximum(1.3);
  hMuMu_costhetastar_data_ratio->SetMinimum(0.7);
  hMuMu_costhetastar_data_ratio->GetXaxis()->SetLabelOffset(0.03);
  hMuMu_costhetastar_data_ratio->Draw("pex0");
  
   hMuMu_costhetastar_data_gUPC_ratio->Draw("pex0same"); 
  hMuMu_costhetastar_data_SL_ratio->Draw("pex0same"); 
  //hMuMu_costhetastar_data_SL_ratio->Draw("histsame"); 
 
  hMuMu_costhetastar_data_ratio_syst->Draw("e2same");
  
  line->DrawCopy("same");
 
  c8->Print(Form("./%s/data_costhetastar_QED_gUPCPythia8_Comparison_30thMarch.pdf",dir));
  c8->Print(Form("./%s/data_costhetastar_QED_gUPCPythia8_Comparison_ForAN_Paper_30thMarch.C",dir));
 

  printBinCont(hMuMu_pt_unfodata);
  printBinCont(hMuMu_rapidity_unfodata);
  printBinCont(hMuMu_invmass_unfodata);
  printBinCont(hMuMu_costhetastar_unfodata);
  outf->cd();
// hist_invmass_data_ratio->Write();
// hist_rapidity_data_ratio->Write();
// hist_pt_data_ratio->Write();
//hist_costhetastar_data_ratio->Write();
//outf->Write();
outf->Close();
  
}

void printBinCont(TH1D *hist){
  
  double bincont = 0;
  
  for (int i=1; i<=hist->GetNbinsX(); i++) {
    
    cout << "******bin:" << i << " " << hist->GetBinContent(i) << "+/-" <<  hist->GetBinError(i) << endl;
    bincont +=  hist->GetBinContent(i);
    
  }
  cout << "Total bincontent:" << bincont << endl;
}


void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}

void make_hist(TH1D *& hist, string xtitle, string ytitle, int kstyle, float ksize, Color_t kcolor, int lstyle, int lwidth){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(40, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.05, "XYZ");
  
  // For the axis titles:
  hist->SetTitleFont(42, "XYZ");
  hist->SetTitleSize(0.06, "XYZ");
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleOffset(1.0);
  hist->SetMarkerColor(kcolor);
  hist->SetLineColor(kcolor);
  hist->SetMarkerStyle(kstyle);
  hist->SetMarkerSize(ksize);
  hist->GetYaxis()->SetTitle(ytitle.c_str());
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->SetLineStyle(lstyle);
  hist->SetLineWidth(lwidth);
}

void make_hist_syst(TH1D *& hist, string xtitle, string ytitle, int kstyle, float ksize, Color_t kcolor, int lstyle, int lwidth){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(40, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.05, "XYZ");
  
  // For the axis titles:
  hist->SetTitleFont(42, "XYZ");
  hist->SetTitleSize(0.06, "XYZ");
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetTitleOffset(1.0);
  hist->SetFillColor(kcolor);
  hist->SetFillStyle(3005);
  hist->SetMarkerStyle(kstyle);
  hist->SetMarkerSize(ksize);
  hist->GetYaxis()->SetTitle(ytitle.c_str());
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->SetLineStyle(lstyle);
  hist->SetLineWidth(lwidth);
}

TH1D* get_stat_uncertainty_hist(TH1D *input_hist, double uncertainty = syst_uncertainty){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  
  for(int iBin=1; iBin<=output_hist->GetNbinsX(); iBin++){
    double value = output_hist->GetBinContent(iBin);
    output_hist->SetBinError(iBin, uncertainty * value);
  }
  
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.5);

  return output_hist;
}

//Ratio uncertainty plot
TH1D* get_stat_uncertainty_hist_ratio(TH1D *input_hist){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  output_hist->Divide(output_hist);
  for(int i=1; i<=input_hist->GetNbinsX(); i++){
      if(input_hist->GetBinContent(i) == 0){
        output_hist->SetBinError(i,0);
      }
      // else{
      //   output_hist->SetBinContent(i,1);
      // }
  }
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.5);
  output_hist->SetFillStyle(3005); 
  output_hist->SetLineWidth(2);
  return output_hist;
}

// pt uncertainty 
TH1D* get_stat_uncertainty_hist_pt(TH1D *input_hist){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  
  for(int iBin=1; iBin<=output_hist->GetNbinsX(); iBin++){
    double value = output_hist->GetBinContent(iBin);
    double xbincentre = output_hist->GetBinCenter(iBin);
    if(xbincentre < 0.15){
    output_hist->SetBinError(iBin, sqrt(0.01*0.01 + 0.068*0.068) * value);
    }
    else if(xbincentre >= 0.15 &&  xbincentre < 0.7){
      output_hist->SetBinError(iBin, sqrt(0.10*0.10 + 0.068*0.068) * value);
    }
    else if(xbincentre >=0.7 && xbincentre < 1){
      output_hist->SetBinError(iBin, sqrt(0.15*0.15 + 0.068*0.068) * value);
    }
  }
  
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.5);
    output_hist->SetFillStyle(3005); 
  output_hist->SetLineWidth(2);
  return output_hist;
}
// mass uncertainty
TH1D* get_stat_uncertainty_hist_mass(TH1D *input_hist){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  
  for(int iBin=1; iBin<=output_hist->GetNbinsX(); iBin++){
    double value = output_hist->GetBinContent(iBin);
    double xbincentre = output_hist->GetBinCenter(iBin);
    if(xbincentre > 0.0 &&  xbincentre < 30){
    output_hist->SetBinError(iBin, sqrt(0.03*0.03 + 0.068*0.068) * value);
    }
    else if(xbincentre >= 30 &&  xbincentre < 50){
      output_hist->SetBinError(iBin, sqrt(0.10*0.10 + 0.068*0.068) * value);
    }
    else if(xbincentre >=50 && xbincentre < 100){
      output_hist->SetBinError(iBin, sqrt(0.2*0.2 + 0.068*0.068) * value);
    }
    }
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.5);
  output_hist->SetFillStyle(3005); 
  output_hist->SetLineWidth(2);
  return output_hist;
}
//rapidity uncertainty
TH1D* get_stat_uncertainty_hist_rap(TH1D *input_hist){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  
  for(int iBin=1; iBin<=output_hist->GetNbinsX(); iBin++){
    double value = output_hist->GetBinContent(iBin);
    double xbincentre = output_hist->GetBinCenter(iBin);
    if(abs(xbincentre) < 1){
    output_hist->SetBinError(iBin, sqrt(0.0 + 0.068*0.068) * value);
    }
    else if( abs(xbincentre) >= 1 &&  abs(xbincentre) < 2.2 ){
      output_hist->SetBinError(iBin, sqrt(0.05*0.05 + 0.068*0.068) * value);
    }
    else if(abs(xbincentre) > 2.2){
      output_hist->SetBinError(iBin, 0.0 * value);
    }
    
    }
    output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.5);
      output_hist->SetFillStyle(3005); 
  output_hist->SetLineWidth(2);
  return output_hist;
  }
  
  TH1D* get_stat_uncertainty_hist_cos(TH1D *input_hist){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  
  for(int iBin=1; iBin<=output_hist->GetNbinsX(); iBin++){
    double value = output_hist->GetBinContent(iBin);
    double xbincentre = output_hist->GetBinCenter(iBin);
    if(abs(xbincentre) < 1){
    output_hist->SetBinError(iBin, sqrt(0.025*0.025 + 0.068*0.068) * value);
    }
  }
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.5);
  output_hist->SetFillStyle(3005); 
  output_hist->SetLineWidth(2);
  return output_hist;
  
}
  

void prepare_canvas(TCanvas *canvas){
  
  float T = 0.08;
  float B = 0.14;
  float L = 0.14;
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
  canvas->cd();
}

TH1D* get_ratio(TH1D *hist_num, TH1D* hist_den){
  auto hist_ratio = (TH1D*)hist_num->Clone();
  hist_ratio->Divide(hist_den);

  hist_ratio->SetMinimum(0.0);
  hist_ratio->SetMaximum(2.0);

  hist_ratio->GetXaxis()->SetTitleFont(43);
  hist_ratio->GetXaxis()->SetTitleSize(22);
  hist_ratio->GetXaxis()->SetTitleOffset(1.0);
  hist_ratio->GetXaxis()->SetLabelFont(43);
  hist_ratio->GetXaxis()->SetLabelSize(18);
  hist_ratio->GetXaxis()->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("UnfoData/Gen");
  hist_ratio->GetYaxis()->SetTitleFont(43);
  hist_ratio->GetYaxis()->SetTitleSize(18);
  hist_ratio->GetYaxis()->SetTitleOffset(2);
  hist_ratio->GetYaxis()->SetLabelFont(43);
  hist_ratio->GetYaxis()->SetLabelSize(18);
  hist_ratio->GetYaxis()->SetNdivisions(5);
  return hist_ratio;
}


