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
#include <TLatex.h>
double syst_uncertainty = 0.24;
double syst_uncertainty_mc = 0.24;

void make_hist(TH1D *&, string, string, int, float, Color_t, int, int) ;
TH1D* get_stat_uncertainty_hist(TH1D *input_hist, double uncertainty = syst_uncertainty);
TH1D* get_stat_uncertainty_hist_ratio(TH1D *input_hist);
void printBinCont(TH1D *hist);
void prepare_canvas(TCanvas *canvas);
TH1D* get_ratio(TH1D *hist_num, TH1D* hist_den);

int mc_gen_color = kCyan+2;
int data_unfold_color = kBlack;
int mc_reco_color = kBlack;
int data_reco_color = kOrange+1;
int data_unfold_invert_color = kRed;
//int gUPC_color = kRed-7;//kBlue+2;
int gUPC_color = kBlue+2;
void draw_LbL_gUPC(int method =1)
{
  
  //  gStyle->SetErrorX(0); 
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
  
  string mc_name = "SC";

   //gUPC
  TFile *f2= new TFile(Form("Comparison/gammaUPC_predictions_LbL.root"),"r");
  TH1D *hist_pt_gUPC          = (TH1D*)f2->Get("hgammaUPC_LbL_PtPair");
  TH1D *hist_rapidity_gUPC    = (TH1D*)f2->Get("hgammaUPC_LbL_Rap");
  TH1D *hist_invmass_gUPC     = (TH1D*)f2->Get("hgammaUPC_LbL_Invmass");
  hist_pt_gUPC->Scale(0.001);
  hist_rapidity_gUPC->Scale(0.001);
  hist_invmass_gUPC->Scale(0.001);

///
 //TFile *f1= new TFile(Form("unfolding_histograms_Bayes_%s_Final.root",mc_name.c_str()),"r");
  TFile *f1= new TFile(Form("unfolding_histograms_Bayes_%s.root",mc_name.c_str()),"r");
  cout<<"   "<<f1->GetName()<<endl;
  TFile *outf= new TFile(Form("gUPC_Comparison.root"),"recreate"); 
  

  TH1D *hist_pt_gen          = (TH1D*)f1->Get("hGenPt_xSec");
  TH1D *hist_pt_recomc       = (TH1D*)f1->Get("hRecoMCPt_xSec");
  TH1D *hist_pt_recodata     = (TH1D*)f1->Get("hRecoDataPt_xSec");
  TH1D *hist_pt_unfodata     = (TH1D*)f1->Get("hUnfoDataPt_xSec");
  TH1D *hist_pt_unfomc       = (TH1D*)f1->Get("hUnfoMCPt_xSec");

  
  TH1D *hist_rapidity_gen          = (TH1D*)f1->Get("hGenRap_xSec");
  TH1D *hist_rapidity_recomc       = (TH1D*)f1->Get("hRecoMCRap_xSec");
  TH1D *hist_rapidity_recodata     = (TH1D*)f1->Get("hRecoDataRap_xSec");
  TH1D *hist_rapidity_unfodata     = (TH1D*)f1->Get("hUnfoDataRap_xSec");
  TH1D *hist_rapidity_unfomc       = (TH1D*)f1->Get("hUnfoMCRap_xSec");
    
  double renorm1 = 107/hist_rapidity_unfodata->Integral("width"); 
  hist_rapidity_unfodata->Scale(renorm1);
  
  TH1D *hist_invmass_gen          = (TH1D*)f1->Get("hGenInvmass_xSec");
  TH1D *hist_invmass_recomc       = (TH1D*)f1->Get("hRecoMCInvmass_xSec");
  TH1D *hist_invmass_recodata     = (TH1D*)f1->Get("hRecoDataInvmass_xSec");
  TH1D *hist_invmass_unfodata     = (TH1D*)f1->Get("hUnfoDataInvmass_xSec");
  TH1D *hist_invmass_unfomc       = (TH1D*)f1->Get("hUnfoMCInvmass_xSec");
  double renorm2 = 107/hist_invmass_unfodata->Integral("width"); 
  hist_invmass_unfodata->Scale(renorm2);

  //Costhetastar
  TH1D *hist_costhetastar_gen          = (TH1D*)f1->Get("hGenCosthetastar_xSec");
  TH1D *hist_costhetastar_recomc       = (TH1D*)f1->Get("hRecoMCCosthetastar_xSec");
  TH1D *hist_costhetastar_recodata     = (TH1D*)f1->Get("hRecoDataCosthetastar_xSec");
  TH1D *hist_costhetastar_unfodata     = (TH1D*)f1->Get("hUnfoDataCosthetastar_xSec");
  TH1D *hist_costhetastar_unfomc       = (TH1D*)f1->Get("hUnfoMCCosthetastar_xSec");

  TH1D *hist_invmass_unfodata_Invert       = (TH1D*)f1->Get("hUnfoDataInvmass_Invert_xSec");
  TH1D *hist_rapidity_unfodata_Invert     = (TH1D*)f1->Get("hUnfoDataRap_Invert_xSec");
  TH1D *hist_pt_unfodata_Invert     = (TH1D*)f1->Get("hUnfoDataPt_Invert_xSec");
  TH1D *hist_costhetastar_unfodata_Invert     = (TH1D*)f1->Get("hUnfoDataCosthetastar_Invert_xSec");

  TH1D *hist_invmass_unfomc_Invert       = (TH1D*)f1->Get("hUnfoMCInvmass_Invert_xSec");
  TH1D *hist_rapidity_unfomc_Invert     = (TH1D*)f1->Get("hUnfoMCRap_Invert_xSec");
  TH1D *hist_pt_unfomc_Invert     = (TH1D*)f1->Get("hUnfoMCPt_Invert_xSec");
  TH1D *hist_costhetastar_unfomc_Invert     = (TH1D*)f1->Get("hUnfoMCCosthetastar_Invert_xSec");

  cout<<"ok"<<endl;
  
  TLegend *leg1=new TLegend(0.30,0.72,0.90,0.9);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(43);
  leg1->SetTextSize(20);
  leg1->AddEntry(hist_pt_unfomc,"Unfolded MC (Bayes,iteration 4)","pl");
  leg1->AddEntry(hist_pt_gen, "Gen-level MC (Superchic)","pl");
  leg1->AddEntry(hist_pt_recomc, "Reconstructed MC (Superchic)", "pl");
 
  
  TLegend *leg2=new TLegend(0.30,0.68,0.90,0.9);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(20);
  
  leg2->AddEntry(hist_pt_unfodata,"Unfolded data (Bayes,1 iteration)","pl");
  leg2->AddEntry(hist_pt_gen, "SUPERCHIC", "pl");
  leg2->AddEntry(hist_pt_gUPC, "gamma-UPC@NLO", "pl");
  
  

  

 TLegend *leg4=new TLegend(0.40,0.72,0.90,0.9);
  leg4->SetFillColor(0);
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetTextFont(43);
  leg4->SetTextSize(20);
  leg4->AddEntry(hist_rapidity_unfomc,"Unfolded MC (Matrix Inversion)","pl");
  leg4->AddEntry(hist_rapidity_gen, "Gen-level MC (Superchic)","pl");
  leg4->AddEntry(hist_rapidity_recomc, "Reconstructed MC (Superchic)", "pl");
 

  int W = 600;
  int H = 670;//500

  auto line = new TF1("line", "1", -100, 100);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);

  TLatex *mark = new TLatex();
       mark->SetNDC(true);
       double startY = 0.90;
       mark->SetTextFont(42);
       mark->SetTextSize(0.035);
       mark->DrawLatex(0.16,startY,"#bf{CMS} #it{Preliminary}");
       mark->DrawLatex(0.55,startY,"#scale[0.8]{1.64 nb^{-1} (PbPb @ 5.02 TeV)}");
       //mark->Draw();
  


  
  //Pt  
  TCanvas* c2 = new TCanvas("c2","Diphoton pT",50,50,W,H);
  prepare_canvas(c2);
  TPad *pad3 = new TPad("pad3", "pad3", 0, 0.3, 1, 1.0);
  pad3->SetBottomMargin(0);
  pad3->SetRightMargin(0.08); 
  pad3->Draw();
  pad3->cd();
  pad3->SetLogy();  
   string x_title = "p_{T}^{#gamma#gamma} (GeV)";
  string y_title = "d#sigma_{#gamma#gamma}/dp_{T}^{#gamma#gamma} (nb/0.2GeV)";

  hist_pt_unfodata->SetMaximum(2*pow(10.,3.));
  hist_pt_unfodata->SetMinimum(5*pow(10.,-1.));
  make_hist(hist_pt_gen,x_title,y_title,20,1.0,mc_gen_color,1,2);
  make_hist(hist_pt_unfodata,x_title,y_title,20,1.0,data_unfold_color,2,2);
  make_hist(hist_pt_unfodata_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2);
  make_hist(hist_pt_gUPC,x_title,y_title,20,1.0,kGreen,2,2);
  auto hist_pt_unfodata_syst = get_stat_uncertainty_hist(hist_pt_unfodata);
  auto hist_pt_gen_syst = get_stat_uncertainty_hist(hist_pt_gen, syst_uncertainty_mc);
  auto hist_pt_gUPC_syst = get_stat_uncertainty_hist(hist_pt_gUPC, syst_uncertainty_mc);
  hist_pt_unfodata->SetBinErrorOption(TH1::kPoisson);
  hist_pt_unfodata_syst->Draw("e2");
  hist_pt_gen->Draw("histsamex0");
  hist_pt_gen_syst->Draw("e2same");
  hist_pt_gUPC_syst->Draw("e2same");
  hist_pt_unfodata->Draw("psameex0");
  hist_pt_unfodata_Invert->Draw("psameex0");
  hist_pt_gUPC->Draw("histsamex0");
  hist_pt_gUPC->Draw("psameex0");
  leg2->Draw();
  c2->cd();
  TLatex *mark6 = new TLatex();
  mark6->SetTextSize(0.035); 
  mark6->SetTextFont(42);
  
   mark6->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark6->DrawLatex(0.55, 0.97, "#scale[0.8]{PbPb, 1.65 nb^{-1} (#sqrt{S_{NN}} = 5.02 TeV)}");
  TPad *pad4 = new TPad("pad4", "pad4", 0, 0.06, 1, 0.28);
  pad4->SetTopMargin(0);
  pad4->SetBottomMargin(0.3);
  pad4->SetRightMargin(0.08); 
 // pad4->SetGridx();
  pad4->Draw();
  pad4->cd();
  TH1D *hist_pt_data_ratio = get_ratio(hist_pt_unfodata, hist_pt_gen);
  TH1D *hist_pt_data_invert_ratio = get_ratio(hist_pt_unfodata_Invert, hist_pt_gen);
  TH1D *hist_pt_data_gUPC_ratio = get_ratio(hist_pt_unfodata_Invert, hist_pt_gUPC);
  auto hist_pt_data_ratio_syst = get_stat_uncertainty_hist(hist_pt_data_ratio);
  auto hist_pt_data_gUPC_ratio_syst = get_stat_uncertainty_hist(hist_pt_data_gUPC_ratio);

  hist_pt_data_invert_ratio->GetXaxis()->SetTitle("p_{T}^{#gamma#gamma} (GeV)");
  hist_pt_data_invert_ratio->GetYaxis()->SetRangeUser(0.0,3.0);
  hist_pt_data_invert_ratio->Draw("pex0same");
  hist_pt_data_ratio->Draw("pex0same");
  hist_pt_data_gUPC_ratio->Draw("pex0same");
  hist_pt_data_ratio_syst->Draw("e2same");
  hist_pt_data_gUPC_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
  c2->Print("fig_Bayes_Inverse/data_pt_gUPC_Comparison.pdf");

  // Rapidity  
  //Rapidity
  TCanvas* c4 = new TCanvas("c4","Diphoton rapidity",50,50,W,H);
  prepare_canvas(c4);
  TPad *pad9 = new TPad("pad9", "pad9", 0, 0.3, 1, 1.0);
  pad9->SetBottomMargin(0.017);
  pad9->SetRightMargin(0.08); 
  pad9->SetTopMargin(0.08);
  pad9->Draw();
  pad9->cd();
  //gPad->SetLogy();
  hist_rapidity_unfodata->SetMaximum(60);
  hist_rapidity_unfodata->SetMinimum(0.0);

 
  x_title = "y^{#gamma#gamma}";
  y_title = "d#sigma^{#gamma#gamma}/dy^{#gamma#gamma} (nb)";
 
  make_hist(hist_rapidity_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  make_hist(hist_rapidity_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hist_rapidity_gUPC,x_title,y_title,2,1.0,gUPC_color,1,2);
 
  auto hist_rapidity_unfodata_syst = get_stat_uncertainty_hist(hist_rapidity_unfodata,syst_uncertainty);
  auto hist_rapidity_gen_syst = get_stat_uncertainty_hist(hist_rapidity_gen);
  auto hist_rapidity_gUPC_syst = get_stat_uncertainty_hist(hist_rapidity_gUPC );
  hist_rapidity_unfodata_syst->GetXaxis()->SetLabelSize(0);
  hist_rapidity_unfodata_syst->Draw("e2");
  hist_rapidity_unfodata->Draw("psameex0");
  hist_rapidity_gen->Draw("histsame");
  hist_rapidity_gen->Draw("psameex0");

  hist_rapidity_unfodata->Draw("pe1sameex0");
  hist_rapidity_gUPC->Draw("histsame");
  hist_rapidity_gUPC->Draw("psamex0");
  
  TH1D *hist_rapidity_data_ratio = get_ratio(hist_rapidity_unfodata, hist_rapidity_gen);
  TH1D *hist_rapidity_data_gUPC_ratio = get_ratio(hist_rapidity_unfodata, hist_rapidity_gUPC);
  // make_hist(hMuMu_rapidity_data_ratio,x_title,y_title,20,1.0,mc_gen_color,2,2);
   make_hist(hist_rapidity_data_gUPC_ratio,x_title,y_title,22,1.0,gUPC_color,1,2);
   hist_rapidity_data_ratio->SetLineColor(mc_gen_color);
   hist_rapidity_data_ratio->SetMarkerColor(mc_gen_color);

  auto hist_rapidity_data_ratio_syst = get_stat_uncertainty_hist_ratio(hist_rapidity_unfodata_syst);
  auto hist_rapidity_data_gUPC_ratio_syst = get_stat_uncertainty_hist_ratio(hist_rapidity_gUPC);

  make_hist(hist_rapidity_data_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  make_hist(hist_rapidity_data_gUPC_ratio_syst,x_title,y_title,0,1.0,gUPC_color,2,2);

  TLegend *leg3=new TLegend(0.40,0.68,0.90,0.9);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextFont(43);
  leg3->SetTextSize(20);
  
  leg3->SetHeader("#gamma#gamma#rightarrow #gamma#gamma");
  leg3->AddEntry(hist_rapidity_unfodata,"Data ","EP");
  leg3->AddEntry(hist_rapidity_data_gUPC_ratio, "gamma-UPC@NLO", "pl");
  leg3->AddEntry(hist_rapidity_data_ratio, "SUPERCHIC 3.03", "pl");
  leg3->Draw();
  c4->cd();
  TLatex *mark9 = new TLatex();
  mark9->SetTextSize(0.035); 
  mark9->SetTextFont(42);
  
   mark9->DrawLatex(0.16, 0.96, "#bf{CMS} #it{Preliminary}");
   mark9->DrawLatex(0.52, 0.96, "#scale[0.8]{PbPb, 1.65 nb^{-1} (#sqrt{s_{NN}} = 5.02 TeV)}");
  TPad *pad10 = new TPad("pad10", "pad10", 0, 0.06, 1, 0.28);
  pad10->SetTopMargin(0.047);
  pad10->SetBottomMargin(0.45);
 // pad10->SetGridx();
  pad10->Draw();
  pad10->cd();
  pad10->SetRightMargin(0.08);
  
  
  hist_rapidity_data_ratio->GetXaxis()->SetTitle("y^{#gamma#gamma}");
  hist_rapidity_data_ratio->GetYaxis()->SetTitle("Data/MC");
  
  hist_rapidity_data_ratio->SetMaximum(3.0);
  hist_rapidity_data_ratio->GetXaxis()->SetLabelOffset(0.03);
  hist_rapidity_data_ratio->Draw("pex0");
  hist_rapidity_data_gUPC_ratio->Draw("pex0same"); 
  hist_rapidity_data_ratio_syst->Draw("e2same");
 
  line->DrawCopy("same");
 
  c4->Print("fig_AN/fig_AN_Paper/data_rapidity_LByL_gUPC_30thMarch.pdf");
  c4->Print("fig_AN/fig_AN_Paper/data_rapidity_LByL_gUPC_ForAN_Paper_30thMarch.C");
  
  //invmass  
  //////////
 TCanvas* c6 = new TCanvas("c6","Diphoton mass",50,50,W,H);
  prepare_canvas(c6);
  TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
  pad11->SetBottomMargin(0.017);
  pad11->SetRightMargin(0.08); 
  pad11->SetTopMargin(0.08);
  pad11->Draw();
  pad11->cd();
  gPad->SetLogy();
   //gStyle->SetErrorX(0);
 hist_invmass_unfodata->SetMaximum(100);
 hist_invmass_unfodata->SetMinimum(0.3);
 hist_invmass_unfodata->SetBinError(4,1.84);
//  hist_invmass_unfodata->SetBinError(3,1.5641);
//  hist_invmass_unfodata->SetBinError(2,3.51786);
//  hist_invmass_unfodata->SetBinError(1,9.57971);
  x_title = "m^{#gamma#gamma} (GeV)";
  y_title = "d#sigma^{#gamma#gamma}/dm^{#gamma#gamma} (nb / GeV)";
 
  make_hist(hist_invmass_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  make_hist(hist_invmass_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hist_invmass_gUPC,x_title,y_title,2,1.0,gUPC_color,1,2);
 
  auto hist_invmass_unfodata_syst = get_stat_uncertainty_hist(hist_invmass_unfodata,syst_uncertainty);
  auto hist_invmass_gen_syst = get_stat_uncertainty_hist(hist_invmass_gen);
  auto hist_invmass_gUPC_syst = get_stat_uncertainty_hist(hist_invmass_gUPC );
  hist_invmass_unfodata_syst->GetXaxis()->SetRangeUser(5.0,21.0);
  hist_invmass_unfodata_syst->GetXaxis()->SetLabelSize(0);
  hist_invmass_unfodata->SetBinErrorOption(TH1::kPoisson);
 //hist_invmass_gen->SetMinimum(0.2);
    
  int ibin = 4;
  double err_low = hist_invmass_unfodata->GetBinErrorLow(ibin);
  double err_up = hist_invmass_unfodata->GetBinErrorUp(ibin);
  cout << "err_low" << err_low << endl;
  cout << "err_up" << err_up << endl;

  hist_invmass_unfodata_syst->Draw("e2");
  hist_invmass_gen->Draw("psame");
  hist_invmass_gen->Draw("histsame");
  hist_invmass_gUPC->Draw("histsame");
  hist_invmass_gUPC->Draw("psamex0");
  hist_invmass_unfodata->Draw("P0X0E0E1SAMEX0");

  leg3->Draw();
  c6->cd();
  TLatex *mark10 = new TLatex();
  mark10->SetTextSize(0.035); 
  mark10->SetTextFont(42);
  
   mark10->DrawLatex(0.16, 0.96, "#bf{CMS} #it{Preliminary}");
   mark10->DrawLatex(0.52, 0.96, "#scale[0.8]{PbPb, 1.65 nb^{-1} (#sqrt{s_{NN}} = 5.02 TeV)}");
  TPad *pad12 = new TPad("pad12", "pad12", 0, 0.06, 1, 0.28);
  pad12->SetTopMargin(0.047);
  pad12->SetBottomMargin(0.45);
 // pad12->SetGridx();
  pad12->Draw();
  pad12->cd();
  pad12->SetRightMargin(0.08);
  TH1D *hist_invmass_data_ratio = get_ratio(hist_invmass_unfodata, hist_invmass_gen);
  TH1D *hist_invmass_data_gUPC_ratio = get_ratio(hist_invmass_unfodata, hist_invmass_gUPC);

  // make_hist(hMuMu_invmass_data_ratio,x_title,y_title,20,1.0,mc_gen_color,2,2);
   make_hist(hist_invmass_data_gUPC_ratio,x_title,y_title,22,1.0,gUPC_color,1,2);
   hist_invmass_data_ratio->SetLineColor(mc_gen_color);
   hist_invmass_data_ratio->SetMarkerColor(mc_gen_color);
  
  auto hist_invmass_data_ratio_syst = get_stat_uncertainty_hist_ratio(hist_invmass_unfodata_syst);
  auto hist_invmass_data_gUPC_ratio_syst = get_stat_uncertainty_hist_ratio(hist_invmass_gUPC);
  
  make_hist(hist_invmass_data_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  make_hist(hist_invmass_data_gUPC_ratio_syst,x_title,y_title,0,1.0,gUPC_color,2,2);
  
  hist_invmass_data_ratio->GetXaxis()->SetTitle("m^{#gamma#gamma} (GeV)");
  hist_invmass_data_ratio->GetYaxis()->SetTitle("Data/MC");
  hist_invmass_data_ratio->GetXaxis()->SetRangeUser(5.0,21.0); 
  hist_invmass_data_ratio->SetMaximum(3.0);
  hist_invmass_data_ratio->SetMinimum(0.1);
  hist_invmass_data_ratio->GetXaxis()->SetLabelOffset(0.03);
  hist_invmass_data_ratio->Draw("pex0");
  hist_invmass_data_gUPC_ratio->Draw("pex0same"); 
  hist_invmass_data_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
  c6->Print("fig_AN/fig_AN_Paper/data_invmass_LByL_gUPC_Comparison_30thMarch.pdf");
  // c6->Print("fig_AN/fig_AN_Paper/data_invmass_LByL_gUPC_Comparison_ForAN_Paper_30thMarch.C");
  //  c6->Print("fig_AN/fig_AN_Paper/data_invmass_LByL_gUPC_Comparison_ForAN_Paper_30thMarch.root");


////Costhetastar
  TCanvas* c8 = new TCanvas("c8","Diphoton costhetastar",50,50,W,H);
  prepare_canvas(c8);
  TPad *pad15 = new TPad("pad15", "pad15", 0, 0.3, 1, 1.0);
  pad15->SetBottomMargin(0.0);
  pad15->SetRightMargin(0.08); 
  pad15->Draw();
  pad15->cd();
  hist_costhetastar_unfodata->SetMaximum(100);//50
  hist_costhetastar_unfodata->SetMinimum(0);
  make_hist(hist_costhetastar_gen,x_title,y_title,20,1.0,mc_gen_color,1,2);
  make_hist(hist_costhetastar_unfodata,x_title,y_title,20,1.0,data_unfold_color,2,2);
  auto hist_costhetastar_unfodata_syst = get_stat_uncertainty_hist(hist_costhetastar_unfodata);
  auto hist_costhetastar_gen_syst = get_stat_uncertainty_hist(hist_costhetastar_gen, syst_uncertainty_mc);
  hist_costhetastar_unfodata_syst->Draw("e2");
  make_hist(hist_costhetastar_unfodata_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2); 
  hist_costhetastar_unfodata->Draw("psameex0");
   hist_costhetastar_gen->Draw("histsamex0");
    hist_costhetastar_gen_syst->Draw("e2same");
    hist_costhetastar_recomc->Draw("ehistsamex0");
// hist_invmass_unfodata_Invert->Draw("psameex0");
   hist_costhetastar_unfodata->Draw("psameex0");
 // hist_invmass_unfodata_Invert->Draw("psameex0");
  make_hist(hist_costhetastar_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  hist_costhetastar_recodata->Draw("ehistsamex0");
  leg3->Draw();
  
  c8->cd();
  TLatex *mark12 = new TLatex();
  mark12->SetTextSize(0.035); 
  mark12->SetTextFont(42);
  
   mark12->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark12->DrawLatex(0.55, 0.97, "#scale[0.8]{1.65 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad16 = new TPad("pad16", "pad16", 0, 0.06, 1, 0.28);
  pad16->SetTopMargin(0);
  pad16->SetBottomMargin(0.32);
 // pad12->SetGridx();
  pad16->Draw();
  pad16->cd();
  pad16->SetRightMargin(0.08);
  TH1D *hist_costhetastar_data_ratio = get_ratio(hist_costhetastar_unfodata, hist_costhetastar_gen);
  TH1D *hist_costhetastar_data_invert_ratio = get_ratio(hist_costhetastar_unfodata_Invert, hist_costhetastar_gen);
  auto hist_costhetastar_data_ratio_syst = get_stat_uncertainty_hist(hist_costhetastar_data_ratio);
  hist_costhetastar_data_ratio->GetXaxis()->SetTitle("|cos#theta^{*}|_{#gamma#gamma}");
  hist_costhetastar_data_ratio->Draw("pex0");
  hist_costhetastar_data_ratio->SetMaximum(3.0);
//  hist_invmass_data_invert_ratio->Draw("pex0same");
  hist_costhetastar_data_ratio->Draw("pex0same");
  hist_costhetastar_data_ratio_syst->Draw("e2same");
   line->DrawCopy("same");
  // c8->Print("fig_Bayes_Inverse/data_costhetastar_ZDCAND3n_SCSL.pdf");
  // c8->Print("fig_Bayes_Inverse/data_costhetastar_ZDCAND3n_SCSL.png");
  //c8->Print("fig_Bayes_Inverse/data_costhetastar_ZDCAND3n.C");
  //////////////
  //Comparison with GammaUPC
  



  
  ///////
  printBinCont(hist_pt_unfodata);
  printBinCont(hist_rapidity_unfodata);
  // printBinCont(hist_rapidity_unfodata_syst);
  // printBinCont(hist_rapidity_data_ratio_syst);
  printBinCont(hist_invmass_unfodata);

  // printBinCont(hist_invmass_unfodata_syst);
  // printBinCont(hist_invmass_data_ratio_syst);
  printBinCont(hist_costhetastar_unfodata);
  printBinCont(hist_invmass_gen);

  outf->cd();
  hist_pt_gen->Write();
  hist_pt_gen_syst->Write();
  hist_pt_unfodata->Write();
  hist_pt_unfodata_syst->Write();
  hist_pt_data_ratio->Write();
  hist_pt_data_ratio_syst->Write();
  hist_rapidity_gen->Write();
  hist_rapidity_gen_syst->Write();
  hist_rapidity_unfodata->Write();
  hist_rapidity_unfodata_syst->Write();
  hist_rapidity_data_ratio->Write();
  hist_rapidity_data_ratio_syst->Write();
  hist_invmass_gen->Write();
  hist_invmass_gen_syst->Write();
  hist_invmass_unfodata->Write();
  hist_invmass_unfodata_syst->Write();
  hist_invmass_data_ratio->Write();
  hist_invmass_data_ratio_syst->Write();
  hist_costhetastar_gen->Write();
  hist_costhetastar_gen_syst->Write();
  hist_costhetastar_unfodata->Write();
  hist_costhetastar_unfodata_syst->Write();
  hist_costhetastar_data_ratio->Write();
  hist_costhetastar_data_ratio_syst->Write();

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
  tex->SetTextColor(data_reco_color);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}

void make_hist(TH1D *& hist, string xtitle, string ytitle, int kstyle, float ksize, Color_t kcolor, int lstyle, int lwidth){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.05, "XYZ");
  
  // For the axis titles:
  hist->SetTitleFont(42, "XYZ");
  hist->SetTitleSize(0.06, "XYZ");
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->SetMarkerColor(kcolor);
  hist->SetLineColor(kcolor);
  hist->SetMarkerStyle(kstyle);
  hist->SetMarkerSize(ksize);
  hist->GetYaxis()->SetTitle(ytitle.c_str());
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->SetLineStyle(lstyle);
  hist->SetLineWidth(lwidth);
//  hist->SetMaximum(300);  
}

TH1D* get_stat_uncertainty_hist(TH1D *input_hist, double uncertainty = syst_uncertainty){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  
  for(int iBin=1; iBin<=output_hist->GetNbinsX(); iBin++){
    double value = output_hist->GetBinContent(iBin);
    output_hist->SetBinError(iBin, uncertainty * value);
  } 
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.5);
   output_hist->SetFillStyle(3005); 
  output_hist->SetLineWidth(2);
  return output_hist;
}

TH1D* get_stat_uncertainty_hist_ratio(TH1D *input_hist){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  output_hist->Divide(output_hist);
   for(int i=1; i<=input_hist->GetNbinsX(); i++){
       if(input_hist->GetBinContent(i) <= 0){
         output_hist->SetBinError(i,0);
       }
      //  else{
      //    output_hist->SetBinContent(i,1);
      //  }
  }
  
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.5);
  output_hist->SetFillStyle(3005); 
  output_hist->SetLineWidth(2);
  return output_hist;
}
void prepare_canvas(TCanvas *canvas){
  float T = 0.08;
  float B = 0.18;
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

TH1D* get_ratio(TH1D *hist_num, TH1D* hist_den){
  auto hist_ratio = (TH1D*)hist_num->Clone();
  hist_ratio->Divide(hist_den);

  hist_ratio->SetMinimum(0.0);
  hist_ratio->SetMaximum(4.0);

  hist_ratio->GetXaxis()->SetTitleFont(43);
  hist_ratio->GetXaxis()->SetTitleSize(20);
  hist_ratio->GetXaxis()->SetTitleOffset(1.0);
  hist_ratio->GetXaxis()->SetLabelFont(43);
  hist_ratio->GetXaxis()->SetLabelSize(20);
  hist_ratio->GetXaxis()->SetTitle("");
  hist_ratio->GetYaxis()->SetTitle("UnfoData/Gen");
  hist_ratio->GetYaxis()->SetTitleFont(43);
  hist_ratio->GetYaxis()->SetTitleSize(22);
  hist_ratio->GetYaxis()->SetTitleOffset(2);
  hist_ratio->GetYaxis()->SetLabelFont(43);
  hist_ratio->GetYaxis()->SetLabelSize(20);
  hist_ratio->GetYaxis()->SetNdivisions(5);
  

  return hist_ratio;
}

