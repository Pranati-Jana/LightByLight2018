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
TH1D* get_stat_uncertainty_hist(TH1D *input_hist, double uncertainty = syst_uncertainty);
TH1D* get_stat_uncertainty_hist_pt(TH1D *input_hist);
TH1D* get_stat_uncertainty_hist_mass(TH1D *input_hist);
TH1D* get_stat_uncertainty_hist_rap(TH1D *input_hist);
TH1D* get_stat_uncertainty_hist_ratio(TH1D *input_hist);
TH1D* get_stat_uncertainty_hist_cos(TH1D *input_hist);

const char *dir  = "fig_AN";
void printBinCont(TH1D *hist);
void prepare_canvas(TCanvas *canvas);
TH1D* get_ratio(TH1D *hist_num, TH1D* hist_den);


int mc_gen_color = kBlue;
int data_unfold_color = kRed;
int mc_reco_color = kBlack;
int data_reco_color = kOrange+1;

void draw_unfold_QED_SC(int method =1, string mc_name = "SC")
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
  
  TFile *f1 = new TFile(Form("unfolding_histograms_Bayes_QEDSC%s.root", mc_name.c_str()),"r");
  TFile *f2 = new TFile(Form("unfolding_histograms_Bayes_QEDSCSC_gUPC.root"),"r");
  TH1D  *h1 = (TH1D*)f2->Get("hUnfoDataPt_xSec");
  h1->Scale(0.001);
  TFile *outf= new TFile(Form("After_Unfolding_QEDSC_costheta.root"),"recreate"); 
  

  if(!f1){
    cout<<"File not found: "<<Form("unfolding_histograms_%s.root", MethID1)<<endl;
    exit(0);
  }
  
  //TFile *f1= new TFile(Form("unfolding_QEDSL_histograms_%s.root",MethID1),"r");
  cout<<"   "<<f1->GetName()<<endl;
  
  TH1D *hMuMu_pt_gen          = (TH1D*)f1->Get("hGenPt_xSec");
  TH1D *hMuMu_pt_recomc       = (TH1D*)f1->Get("hRecoMCPt_xSec");
  TH1D *hMuMu_pt_recodata     = (TH1D*)f1->Get("hRecoDataPt_xSec");
  TH1D *hMuMu_pt_unfodata     = (TH1D*)f1->Get("hUnfoDataPt_xSec");
  TH1D *hMuMu_pt_unfomc       = (TH1D*)f1->Get("hUnfoMCPt_xSec");

  hMuMu_pt_gen->Scale(0.001);
  hMuMu_pt_recomc->Scale(0.001);
  hMuMu_pt_recodata->Scale(0.001);
  hMuMu_pt_unfodata->Scale(0.001);
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
  hMuMu_costhetastar_unfomc->Scale(0.001);

  
  cout<<"ok"<<endl;
  
  TLegend *leg1=new TLegend(0.25,0.72,0.90,0.9);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(43);
  leg1->SetTextSize(16);
  leg1->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  leg1->AddEntry(hMuMu_pt_unfomc,"Unfolded MC (Bayes, 3 iteration)","pl");
  leg1->AddEntry(hMuMu_pt_gen, "Gen-level MC (Superchic 3.02+Photos)","pl");
  leg1->AddEntry(hMuMu_pt_recomc,"Reconstructed MC (Superchic 3.02+Photos)","pl");
  
  
  TLegend *leg2=new TLegend(0.25,0.68,0.90,0.9);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(10);
  leg2->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  leg2->AddEntry(hMuMu_pt_recodata,"Reconstructed data","pl");
  leg2->AddEntry(hMuMu_pt_unfodata,"Unfolded data (Bayes, 3 iteration)","pl");
  //leg2->AddEntry(hMuMu_pt_unfodata,"Unfolded data (Matrix inversion)","pl");
  leg2->AddEntry(hMuMu_pt_recomc, "Reconstructed MC (Superchic 3.02+Photos)","pl");
  leg2->AddEntry(hMuMu_pt_gen,  "Gen-level MC (Superchic 3.02+Photos) ","pl");
  //leg2->AddEntry(hMuMu_pt_gen,  "Gen MC (SC 3.02+FSR(PHOTOS) reweighted to gamma-UPC+FSR(PY8) ) ","pl");
  leg2->SetEntrySeparation(0.5); 

  TLegend *leg3=new TLegend(0.25,0.68,0.90,0.9);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextFont(43);
  leg3->SetTextSize(16);
  leg3->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  leg3->AddEntry(hMuMu_rapidity_recodata,"Reconstructed data","pl");
  leg3->AddEntry(hMuMu_rapidity_unfodata,"Unfolded data (Matrix inversion)","pl");
  leg3->AddEntry(hMuMu_rapidity_recomc, "Reconstructed MC (Superchic 3.02+Photos)","pl");
  leg3->AddEntry(hMuMu_rapidity_gen,  "Gen-level MC (Superchic 3.02+Photos)","pl");

  TLegend *leg6=new TLegend(0.25,0.72,0.90,0.9);
  leg6->SetFillColor(0);
  leg6->SetBorderSize(0);
  leg6->SetFillStyle(0);
  leg6->SetTextFont(43);
  leg6->SetTextSize(16);
  leg6->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  leg6->AddEntry(hMuMu_rapidity_unfomc,"Unfolded MC (Matrix inversion)","pl");
  leg6->AddEntry(hMuMu_rapidity_gen, "Gen-level MC (Superchic 3.02+Photos)","pl");
  leg6->AddEntry(hMuMu_rapidity_recomc,"Reconstructed MC (Superchic 3.02+Photos)","pl");
  
  TLegend *leg4=new TLegend(0.25,0.68,0.90,0.9);
  leg4->SetFillColor(0);
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetTextFont(43);
  leg4->SetTextSize(16);
  leg4->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  leg4->AddEntry(hMuMu_invmass_recodata,"Reconstructed data","pl");
  //leg4->AddEntry(hMuMu_invmass_unfodata,"Unfolded data (Bayes, 1 iteration)","pl");
  leg4->AddEntry(hMuMu_invmass_unfodata,"Unfolded data (Matrix inversion)","pl");
  leg4->AddEntry(hMuMu_invmass_recomc, "Reconstructed MC (Superchic 3.02+Photos)","pl");
  leg4->AddEntry(hMuMu_invmass_gen,  "Gen-level MC (Superchic 3.02+Photos)","pl");

  TLegend *leg5=new TLegend(0.25,0.72,0.90,0.9);
  leg5->SetFillColor(0);
  leg5->SetBorderSize(0);
  leg5->SetFillStyle(0);
  leg5->SetTextFont(43);
  leg5->SetTextSize(16);
  leg5->SetHeader("#gamma#gamma#rightarrow e^{+}e^{-}");
  leg5->AddEntry(hMuMu_invmass_unfomc,"Unfolded MC (Matrix inversion)","pl");
  leg5->AddEntry(hMuMu_invmass_gen, "Gen-level MC (Superchic 3.02+Photos)","pl");
  leg5->AddEntry(hMuMu_invmass_recomc,"Reconstructed MC (Superchic 3.02+Photos)","pl");
  
  
  

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
  TCanvas* c1 = new TCanvas("c1","dielectron pT",50,50,W,H);
  prepare_canvas(c1);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->SetRightMargin(0.08);
  pad1->SetTopMargin(0.07);
  pad1->Draw();
  pad1->cd();
 // c1->cd();
  gPad->SetLogy();
  
  string x_title = "p_{T}^{ee} (GeV)";
  string y_title = "#frac{d#sigma_{ee}}{dp_{T}^{ee}} (#mub / 0.142 GeV)";
  hMuMu_pt_unfomc->SetMaximum(10000);
 // hMuMu_pt_unfomc->SetMaximum(2*pow(10.,4.));
  //hMuMu_pt_unfomc->SetMinimum(5*pow(10.,1.));
   hMuMu_pt_unfomc->SetMinimum(0.02);
  //hMuMu_pt_unfomc->SetMinimum(10);
  
  cout<<"pt bin width:" <<hMuMu_pt_unfodata->GetXaxis()->GetBinWidth(1)<<endl;

  make_hist(hMuMu_pt_unfomc,x_title,y_title,20,1.0,data_unfold_color,1,1);
  hMuMu_pt_unfomc->Draw("p");
  make_hist(hMuMu_pt_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);

  auto hMuMu_pt_gen_syst = get_stat_uncertainty_hist_pt(hMuMu_pt_unfomc);
  make_hist(hMuMu_pt_gen_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);

  auto hMuMu_pt_recomc_syst = get_stat_uncertainty_hist(hMuMu_pt_recomc);
  make_hist(hMuMu_pt_recomc_syst,x_title,y_title,0,1.0,mc_reco_color,2,2);

  hMuMu_pt_unfomc->Draw("p");
  hMuMu_pt_gen->Draw("histsame");
  hMuMu_pt_gen->Draw("psame");
  hMuMu_pt_gen_syst->Draw("samee2");
  
  make_hist(hMuMu_pt_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hMuMu_pt_recomc->Draw("ehistsamex0");
  hMuMu_pt_recomc_syst->Draw("samee2");
 
  leg1->Draw();

  c1->cd();
  TLatex *mark5 = new TLatex();
  mark5->SetTextSize(0.035); 
  mark5->SetTextFont(42);
   mark5->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark5->DrawLatex(0.62, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.06, 1, 0.28);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.4);
  pad2->SetRightMargin(0.08);
  //pad2->SetGridx();
  pad2->Draw();
  pad2->cd();
  TH1D *hist_pt_closure_ratio = get_ratio(hMuMu_pt_unfomc, hMuMu_pt_gen);
  
  auto hist_pt_closure_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_pt_gen_syst);
  make_hist(hist_pt_closure_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  //make_hist(hist_pt_closure_ratio,x_title,y_title,20,1.0,mc_gen_color,2,2);
  hist_pt_closure_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_pt_closure_ratio->GetXaxis()->SetTitle("p_{T}^{ee} (GeV)");
  hist_pt_closure_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hist_pt_closure_ratio->Draw("pex0");
  //hist_pt_closure_ratio_syst->Draw("e2same");
  hist_pt_closure_ratio->Draw("pex0same");
  hist_pt_closure_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
  mark->Draw();

  c1->Print(Form("./%s/closure_pt_SC_ZDCAND3n_7bin_newUncertainty.pdf",dir));
  
  
  TCanvas* c2 = new TCanvas("c2","dielectron pT",50,50,W,H);
  prepare_canvas(c2);
  TPad *pad3 = new TPad("pad3", "pad3", 0, 0.3, 1, 1.0);
  pad3->SetBottomMargin(0);
  pad3->SetRightMargin(0.08);
    pad3->SetTopMargin(0.07);
  pad3->Draw();
  pad3->cd();
  
  gPad->SetLogy();
  hMuMu_pt_unfodata->SetMaximum(10000);
//  hMuMu_pt_unfodata->SetMaximum(2*pow(10.,4.));
  hMuMu_pt_unfodata->SetMinimum(0.02);
  make_hist(hMuMu_pt_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hMuMu_pt_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  auto hMuMu_pt_unfodata_syst = get_stat_uncertainty_hist(hMuMu_pt_unfodata,0.068);
  auto hMuMu_pt_recodata_syst = get_stat_uncertainty_hist(hMuMu_pt_recodata);
  make_hist(hMuMu_pt_recodata_syst,x_title,y_title,0,1.0,data_reco_color,2,2);

  hMuMu_pt_unfodata_syst->Draw("e2");
  cout<<"Unfolded data pt integral: "<<hMuMu_pt_unfodata->Integral()<<" mub"<<endl;
  hMuMu_pt_unfodata->Draw("psamex0");
  hMuMu_pt_gen->Draw("histsamex0");
  hMuMu_pt_gen->Draw("psame");
  hMuMu_pt_recomc->Draw("ehistsamex0");
  hMuMu_pt_recodata->Draw("ehistsamex0");
  hMuMu_pt_recodata_syst->Draw("samee2");
  hMuMu_pt_recomc_syst->Draw("samee2");
  //CMS_lumi( c2, 1, 10 );
  leg2->Draw();
  c2->cd();
  TLatex *mark6 = new TLatex();
  mark6->SetTextSize(0.035); 
  mark6->SetTextFont(42);
  
   mark6->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark6->DrawLatex(0.62, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad4 = new TPad("pad4", "pad4", 0, 0.06, 1, 0.28);
  pad4->SetTopMargin(0);
  pad4->SetBottomMargin(0.4);
  pad4->SetRightMargin(0.08);
 // pad4->SetGridx();
  pad4->Draw();
  pad4->cd();

  TH1D *hist_pt_data_ratio = get_ratio(hMuMu_pt_unfodata, hMuMu_pt_gen);
  auto hist_pt_data_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_pt_unfodata_syst);
   make_hist(hist_pt_data_ratio_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);
  hist_pt_data_ratio->GetXaxis()->SetTitle("p_{T}^{ee} (GeV)");
  hist_pt_data_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_pt_data_ratio->GetYaxis()->SetRangeUser(0.0,5);
  hist_pt_data_ratio->Draw("pex0same");
  hist_pt_data_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
  
  c2->Print(Form("./%s/data_pt_SC_ZDCAND3n_Bayes_7bin.pdf",dir));

  // rapidity
  TCanvas* c3 = new TCanvas("c3","dielectron rapidity",50,50,W,H);
  prepare_canvas(c3);
  //c3->cd();
  TPad *pad5 = new TPad("pad5", "pad5", 0, 0.3, 1, 1.0);
  pad5->SetBottomMargin(0);
  pad5->SetRightMargin(0.08);
  pad5->SetTopMargin(0.07);
  pad5->Draw();
  pad5->cd();

  gPad->SetLogy();  
  x_title = "y_{ee}";
  y_title = "#frac{d#sigma_{ee}}{dy_{ee}} (#mub / 0.4)";
  
  hMuMu_rapidity_unfomc->SetMaximum(3000);
  hMuMu_rapidity_unfomc->SetMinimum(0.08);
  //hMuMu_rapidity_unfomc->GetXaxis()->SetRangeUser(0,1.0);
  make_hist(hMuMu_rapidity_unfomc,x_title,y_title,20,1.0,data_unfold_color,1,1);
  hMuMu_rapidity_unfomc->Draw("p");
  make_hist(hMuMu_rapidity_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  hMuMu_rapidity_gen->Draw("histsame");
  hMuMu_rapidity_gen->Draw("psame");
  
  auto hMuMu_rapidity_gen_syst = get_stat_uncertainty_hist_rap(hMuMu_rapidity_unfomc);
  make_hist(hMuMu_rapidity_gen_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);

  auto hMuMu_rapidity_recomc_syst = get_stat_uncertainty_hist(hMuMu_rapidity_recomc);
  make_hist(hMuMu_rapidity_recomc_syst,x_title,y_title,0,1.0,mc_reco_color,2,2);

  hMuMu_rapidity_gen_syst->Draw("samee2");
  

  make_hist(hMuMu_rapidity_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hMuMu_rapidity_recomc->Draw("ehistsamex0");
  hMuMu_rapidity_recomc_syst->Draw("samee2");
 
  //CMS_lumi( c3, 1, 10 );
  leg6->Draw();

  c3->cd();
  TLatex *mark7 = new TLatex();
  mark7->SetTextSize(0.035); 
  mark7->SetTextFont(42);
  
   mark7->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark7->DrawLatex(0.62, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad6 = new TPad("pad6", "pad6", 0, 0.06, 1, 0.28);
  pad6->SetTopMargin(0);
  pad6->SetBottomMargin(0.42);
  pad6->SetRightMargin(0.08);
//  pad6->SetGridx();
  pad6->Draw();
  pad6->cd();
  TH1D *hist_rapidity_closure_ratio = get_ratio(hMuMu_rapidity_unfomc, hMuMu_rapidity_gen);
  auto hist_rapidity_closure_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_rapidity_gen_syst);
  make_hist(hist_rapidity_closure_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  hist_rapidity_closure_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_rapidity_closure_ratio->GetXaxis()->SetTitle("y_{ee}");
  hist_rapidity_closure_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hist_rapidity_closure_ratio->Draw("pex0");
  hist_rapidity_closure_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
  //c3->Print(Form("./%s/closure_rapidity_SC_ZDCAND3n_11bin_17thMarch.pdf",dir));
  
  
  TCanvas* c4 = new TCanvas("c4","dielectron rap",50,50,W,H);
  prepare_canvas(c4);
  TPad *pad7 = new TPad("pad7", "pad7", 0, 0.3, 1, 1.0);
  //pad7->SetTopMargin(0.0);
  pad7->SetBottomMargin(0.0);
  pad7->SetTopMargin(0.07);
  pad7->SetRightMargin(0.08);
  pad7->Draw();
  pad7->cd();  

  gPad->SetLogy();
  hMuMu_rapidity_unfodata->SetMaximum(3000);
  hMuMu_rapidity_unfodata->SetMinimum(0.08);
  make_hist(hMuMu_rapidity_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hMuMu_rapidity_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  auto hMuMu_rapidity_unfodata_syst = get_stat_uncertainty_hist_rap(hMuMu_rapidity_unfodata);
  auto hMuMu_rapidity_recodata_syst = get_stat_uncertainty_hist(hMuMu_rapidity_recodata);
  make_hist(hMuMu_rapidity_recodata_syst,x_title,y_title,0,1.0,data_reco_color,2,2);

  hMuMu_rapidity_unfodata_syst->Draw("e2");
  cout<<"Unfolded data rapidity integral: "<<hMuMu_rapidity_unfodata->Integral()<<" mub"<<endl;
  hMuMu_rapidity_unfodata->Draw("psameex0");
  hMuMu_rapidity_gen->Draw("histsamex0");
  hMuMu_rapidity_gen->Draw("psame");
  //hMuMu_rapidity_gen_syst->Draw("samee2");
  hMuMu_rapidity_recomc->Draw("ehistsamex0");

  hMuMu_rapidity_recodata->Draw("ehistsamex0");
  hMuMu_rapidity_recodata_syst->Draw("samee2");
  hMuMu_rapidity_recomc_syst->Draw("samee2");
  //CMS_lumi( c4, 1, 10 );
  leg3->Draw();

  c4->cd();
  TLatex *mark8 = new TLatex();
  mark8->SetTextSize(0.035); 
  mark8->SetTextFont(42);
 
   mark8->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark8->DrawLatex(0.62, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad8 = new TPad("pad8", "pad8", 0, 0.06, 1, 0.28);
  pad8->SetTopMargin(0);
  pad8->SetBottomMargin(0.42);
  pad8->SetRightMargin(0.08);
  //pad8->SetGridx();
  pad8->Draw();
  pad8->cd();
  TH1D *hist_rapidity_data_ratio = get_ratio(hMuMu_rapidity_unfodata, hMuMu_rapidity_gen);
  auto hist_rapidity_data_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_rapidity_unfodata_syst);
  make_hist(hist_rapidity_data_ratio_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);
  hist_rapidity_data_ratio->GetXaxis()->SetTitle("y_{ee}");
  hist_rapidity_data_ratio_syst->GetXaxis()->SetTitle("y_{ee}");
  hist_rapidity_data_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_rapidity_data_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
// hist_rapidity_data_ratio_syst->SetMaximum(2.0);
  hist_rapidity_data_ratio->Draw("pex0");
  hist_rapidity_data_ratio->Draw("pex0same");
  hist_rapidity_data_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
 // c4->Print(Form("./%s/data_rapidity_SC_ZDCAND3n_11bin_17thMarch.pdf",dir));
  
  //invmass
  TCanvas* c5 = new TCanvas("c5","dielectron invmass",50,50,W,H);
  prepare_canvas(c5);
  TPad *pad9 = new TPad("pad9", "pad9", 0, 0.3, 1, 1.0);
  pad9->SetBottomMargin(0);
  pad9->SetRightMargin(0.08);
  pad9->SetTopMargin(0.07);
  pad9->Draw();
  pad9->cd();  

  gPad->SetLogy();
  
  x_title = "m_{ee} (GeV)";
  y_title = "#frac{d#sigma_{ee}}{dm_{ee}} (#mub / bin width)";
  
  hMuMu_invmass_unfomc->SetMaximum(100);
  hMuMu_invmass_unfomc->SetMinimum(0.00002);
  make_hist(hMuMu_invmass_unfomc,x_title,y_title,20,1.0,data_unfold_color,1,1);
  hMuMu_invmass_unfomc->Draw("p");
  make_hist(hMuMu_invmass_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  auto hMuMu_invmass_gen_syst = get_stat_uncertainty_hist_mass(hMuMu_invmass_unfomc);
  auto hMuMu_invmass_recomc_syst = get_stat_uncertainty_hist(hMuMu_invmass_recomc);
  make_hist(hMuMu_invmass_recomc_syst,x_title,y_title,0,1.0,mc_reco_color,2,2);
  hMuMu_invmass_gen_syst->Draw("samee2");
  hMuMu_invmass_recomc_syst->Draw("samee2");
  hMuMu_invmass_gen->Draw("histsame");
  hMuMu_invmass_gen->Draw("psame");
  make_hist(hMuMu_invmass_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hMuMu_invmass_recomc->Draw("ehistsamex0");
  
  //CMS_lumi( c5, 1, 10 );
  leg5->Draw();

  c5->cd();
  TLatex *mark9 = new TLatex();
  mark9->SetTextSize(0.035); 
  mark9->SetTextFont(42);

   mark9->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark9->DrawLatex(0.62, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad10 = new TPad("pad10", "pad10", 0, 0.06, 1, 0.28);
  pad10->SetTopMargin(0);
  pad10->SetBottomMargin(0.4);
  pad10->SetRightMargin(0.08);
 // pad10->SetGridx();
  pad10->Draw();
  pad10->cd();
  TH1D *hist_invmass_closure_ratio = get_ratio(hMuMu_invmass_unfomc, hMuMu_invmass_gen);
  auto hist_invmass_closure_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_invmass_gen_syst);
   make_hist(hist_invmass_closure_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  hist_invmass_closure_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_invmass_closure_ratio->GetXaxis()->SetTitle("m_{ee} (GeV)");
  hist_invmass_closure_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hist_invmass_closure_ratio->Draw("pex0");
  hist_invmass_closure_ratio_syst->Draw("e2same");
  line->DrawCopy("same");

 // c5->Print(Form("./%s/closure_invmass_SC_ZDCAND3n_21bin_17thMarch.pdf",dir));
  
  TCanvas* c6 = new TCanvas("c6","dielectron invmass",50,50,W,H);
  prepare_canvas(c6);
  TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
  pad11->SetBottomMargin(0.01);
  pad11->SetRightMargin(0.08);
   pad11->SetTopMargin(0.07);
  pad11->Draw();
  pad11->cd();
  gPad->SetLogy();
  hMuMu_invmass_unfodata->SetMaximum(100);
  hMuMu_invmass_unfodata->SetMinimum(.00002);
  
  cout<<"m_inv bin width:" <<hMuMu_invmass_unfodata->GetXaxis()->GetBinWidth(1)<<endl;
  
  make_hist(hMuMu_invmass_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  auto hMuMu_invmass_unfodata_syst = get_stat_uncertainty_hist_mass(hMuMu_invmass_unfodata);

  make_hist(hMuMu_invmass_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
 
  auto hMuMu_invmass_recodata_syst = get_stat_uncertainty_hist(hMuMu_invmass_recodata);
  make_hist(hMuMu_invmass_recodata_syst,x_title,y_title,0,1.0,data_reco_color,2,2);

  hMuMu_invmass_unfodata_syst->Draw("e2");
  
  cout<<"Unfolded data mass integral: "<<hMuMu_invmass_unfodata->Integral()<<" mub"<<endl;
  
  hMuMu_invmass_unfodata->Draw("psameex0");
  hMuMu_invmass_gen->Draw("histsamex0");
  hMuMu_invmass_gen->Draw("psame");
  hMuMu_invmass_recodata_syst->Draw("e2same");
  hMuMu_invmass_recomc->Draw("ehistsamex0");
  hMuMu_invmass_recomc_syst->Draw("samee2");
  hMuMu_invmass_recodata->Draw("ehistsamex0");
  
  leg4->Draw();

  c6->cd();
  TLatex *mark10 = new TLatex();
  mark10->SetTextSize(0.035); 
  mark10->SetTextFont(42);
  
   mark10->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark10->DrawLatex(0.62, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad12 = new TPad("pad12", "pad12", 0, 0.06, 1, 0.28);
  pad12->SetTopMargin(0.05);
  pad12->SetBottomMargin(0.45);
  pad12->Draw();
  pad12->cd();
  pad12->SetRightMargin(0.08);
  TH1D *hist_invmass_data_ratio = get_ratio(hMuMu_invmass_unfodata, hMuMu_invmass_gen);
  auto hist_invmass_data_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_invmass_unfodata_syst);
  make_hist(hist_invmass_data_ratio_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);
  hist_invmass_data_ratio->GetXaxis()->SetTitle("m_{ee} (GeV)");
  hist_invmass_data_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_invmass_data_ratio->Draw("pex0");
  hist_invmass_data_ratio->GetYaxis()->SetRangeUser(0.0,1.5);
  hist_invmass_data_ratio->Draw("pex0same");
  hist_invmass_data_ratio_syst->Draw("e2same");
   line->DrawCopy("same");

 // c6->Print(Form("./%s/data_invmass_SC_ZDCAND3n_21bin_17thMarch.pdf",dir));


  //costhetastar
  TCanvas* c7 = new TCanvas("c7","dielectron Costhetastar",50,50,W,H);
  prepare_canvas(c7);
  //c3->cd();
  TPad *pad13 = new TPad("pad13", "pad13", 0, 0.3, 1, 1.0);
  pad13->SetBottomMargin(0);
  pad13->SetRightMargin(0.08);
  pad13->SetTopMargin(0.07);
  pad13->Draw();
  pad13->cd();

  gPad->SetLogy();  
  x_title = "|cos#theta*|_{ee}";
  y_title = "#frac{d#sigma_{ee}}{d|cos#theta*|_{ee}} (#mub / 0.1)";
  
  hMuMu_costhetastar_unfomc->SetMaximum(10000);
  hMuMu_costhetastar_unfomc->SetMinimum(0.08);
  
  make_hist(hMuMu_costhetastar_unfomc,x_title,y_title,20,1.0,data_unfold_color,1,1);
  hMuMu_costhetastar_unfomc->Draw("p");
  make_hist(hMuMu_costhetastar_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  hMuMu_costhetastar_gen->Draw("histsame");
  hMuMu_costhetastar_gen->Draw("psame");
  
  auto hMuMu_costhetastar_gen_syst = get_stat_uncertainty_hist_cos(hMuMu_costhetastar_unfomc);
  make_hist(hMuMu_costhetastar_gen_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);

  auto hMuMu_costhetastar_recomc_syst = get_stat_uncertainty_hist(hMuMu_costhetastar_recomc);
  make_hist(hMuMu_costhetastar_recomc_syst,x_title,y_title,0,1.0,mc_reco_color,2,2);

  hMuMu_costhetastar_gen_syst->Draw("samee2");
  

  make_hist(hMuMu_costhetastar_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hMuMu_costhetastar_recomc->Draw("ehistsamex0");
  hMuMu_costhetastar_recomc_syst->Draw("samee2");
 

  leg6->Draw();

  c7->cd();
  TLatex *mark11 = new TLatex();
  mark11->SetTextSize(0.035); 
  mark11->SetTextFont(42);
  
   mark11->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark11->DrawLatex(0.62, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad14 = new TPad("pad14", "pad14", 0, 0.06, 1, 0.28);
  pad14->SetTopMargin(0);
  pad14->SetBottomMargin(0.42);
  pad14->SetRightMargin(0.08);

  pad14->Draw();
  pad14->cd();
  TH1D *hist_costhetastar_closure_ratio = get_ratio(hMuMu_costhetastar_unfomc, hMuMu_costhetastar_gen);
  auto hist_costhetastar_closure_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_rapidity_gen_syst);
  make_hist(hist_costhetastar_closure_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  hist_costhetastar_closure_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_costhetastar_closure_ratio->GetXaxis()->SetTitle("|cos#theta*|_{ee}");
  hist_costhetastar_closure_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hist_costhetastar_closure_ratio->Draw("pex0");
  hist_costhetastar_closure_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
  //c7->Print(Form("./%s/closure_costhetastar_SC_ZDCAND3n_10bin_17thMarch.pdf",dir));
  
  
  TCanvas* c8 = new TCanvas("c8","dielectron costhetastar data",50,50,W,H);
  prepare_canvas(c8);
  TPad *pad15 = new TPad("pad15", "pad15", 0, 0.3, 1, 1.0);

  pad15->SetBottomMargin(0.0);
  pad15->SetTopMargin(0.07);
  pad15->SetRightMargin(0.08);
  pad15->Draw();
  pad15->cd();  

  gPad->SetLogy();
  hMuMu_costhetastar_unfodata->SetMaximum(10000);
  hMuMu_costhetastar_unfodata->SetMinimum(0.08);
  make_hist(hMuMu_costhetastar_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  make_hist(hMuMu_costhetastar_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  auto hMuMu_costhetastar_unfodata_syst = get_stat_uncertainty_hist_cos(hMuMu_costhetastar_unfodata);
  auto hMuMu_costhetastar_recodata_syst = get_stat_uncertainty_hist(hMuMu_costhetastar_recodata);
  make_hist(hMuMu_costhetastar_recodata_syst,x_title,y_title,0,1.0,data_reco_color,2,2);

  hMuMu_costhetastar_unfodata_syst->Draw("e2");
  cout<<"Unfolded data costhetastar integral: "<<hMuMu_costhetastar_unfodata->Integral()<<" mub"<<endl;
  hMuMu_costhetastar_unfodata->Draw("psameex0");
  hMuMu_costhetastar_gen->Draw("histsamex0");
  hMuMu_costhetastar_gen->Draw("psame");
  //hMuMu_rapidity_gen_syst->Draw("samee2");
  hMuMu_costhetastar_recomc->Draw("ehistsamex0");

  hMuMu_costhetastar_recodata->Draw("ehistsamex0");
  hMuMu_costhetastar_recodata_syst->Draw("samee2");
  hMuMu_costhetastar_recomc_syst->Draw("samee2");
  //CMS_lumi( c4, 1, 10 );
  leg3->Draw();

  c8->cd();
  TLatex *mark12 = new TLatex();
  mark12->SetTextSize(0.035); 
  mark12->SetTextFont(42);
 
   mark12->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark12->DrawLatex(0.62, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad16 = new TPad("pad16", "pad16", 0, 0.06, 1, 0.28);
  pad16->SetTopMargin(0);
  pad16->SetBottomMargin(0.42);
  pad16->SetRightMargin(0.08);

  pad16->Draw();
  pad16->cd();
  TH1D *hist_costhetastar_data_ratio = get_ratio(hMuMu_costhetastar_unfodata, hMuMu_costhetastar_gen);
  auto hist_costhetastar_data_ratio_syst = get_stat_uncertainty_hist_ratio(hMuMu_costhetastar_unfodata_syst);
  make_hist(hist_costhetastar_data_ratio_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);
  hist_costhetastar_data_ratio->GetXaxis()->SetTitle("|cos#theta*|_{ee}");
  hist_costhetastar_data_ratio_syst->GetXaxis()->SetTitle("|cos#theta*|_{ee}");
  hist_costhetastar_data_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_costhetastar_data_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hist_costhetastar_data_ratio->Draw("pex0");
  hist_costhetastar_data_ratio->Draw("pex0same");
  hist_costhetastar_data_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
 // c8->Print(Form("./%s/data_costhetastar_SC_ZDCAND3n_10bin_17thMarch.pdf",dir));
  /////
  
  

  printBinCont(hMuMu_pt_unfodata);
  printBinCont(hMuMu_rapidity_unfodata);
  printBinCont(hMuMu_invmass_unfodata);
  printBinCont(hMuMu_costhetastar_unfodata);
  outf->cd();
// hist_invmass_data_ratio->Write();
// hist_rapidity_data_ratio->Write();
  hist_pt_data_ratio->Write();
//hist_costhetastar_data_ratio->Write();
//outf->Write();
outf->Close();
TCanvas* c = new TCanvas("c","ratio",50,50,W,H);
TH1D *hist_SCgUPC_ratio = get_ratio(h1,hMuMu_pt_unfodata);
hist_SCgUPC_ratio->GetYaxis()->SetRangeUser(0.6,1.1);
hist_SCgUPC_ratio->GetYaxis()->SetTitle("Unfolded gUPC 3 itr/Unfolded Superchic 3 itr");
gPad->SetGridy();
hist_SCgUPC_ratio->Draw("pex0");
line->DrawCopy("same");
//c->SaveAs("itr_gUPC3.png");

  
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

TH1D* get_stat_uncertainty_hist(TH1D *input_hist, double uncertainty = syst_uncertainty){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  
  for(int iBin=1; iBin<=output_hist->GetNbinsX(); iBin++){
    double value = output_hist->GetBinContent(iBin);
    output_hist->SetBinError(iBin, uncertainty * value);
  }
  
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.8);
  output_hist->SetFillStyle(3013); 
  output_hist->SetLineWidth(2);
  return output_hist;
}
//ratio uncertainty
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
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.8);
  output_hist->SetFillStyle(3013); 
  output_hist->SetLineWidth(2);
  return output_hist;
}

//pt uncertainty for main plot
TH1D* get_stat_uncertainty_hist_pt(TH1D *input_hist){
  
  TH1D *output_hist = (TH1D*)input_hist->Clone();
  
  for(int iBin=1; iBin<=output_hist->GetNbinsX(); iBin++){
    double value = output_hist->GetBinContent(iBin);
    double xbincentre = output_hist->GetBinCenter(iBin);
    if(xbincentre < 0.3){
    output_hist->SetBinError(iBin, sqrt( 0.068*0.068) * value);
    //cout << "value:"<<  value << endl;
    }
    else if(xbincentre >= 0.3 &&  xbincentre <= 0.7){
      output_hist->SetBinError(iBin, sqrt(0.0*0.0 + 0.068*0.068) * value);
      //cout << "value:"<<  value << endl;
    }
    else if(xbincentre >0.7 && xbincentre < 1){
      output_hist->SetBinError(iBin, sqrt(0.0*0.0 + 0.068*0.068) * value);
     // cout << "value:"<<  value << endl;
    }
  }
  
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.8);
  output_hist->SetFillStyle(3013); 
  output_hist->SetLineWidth(2);
  return output_hist;
}

// mass uncertainty for main plot
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
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.8);
  output_hist->SetFillStyle(3013); 
  output_hist->SetLineWidth(2);
  return output_hist;
}
//rapidity uncertainty for main plot
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
    output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.8);
        output_hist->SetFillStyle(3013); 
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
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.8);
  output_hist->SetFillStyle(3013); 
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
  hist_ratio->GetXaxis()->SetTitleSize(20);
  hist_ratio->GetXaxis()->SetTitleOffset(1.0);
  hist_ratio->GetXaxis()->SetLabelFont(43);
  hist_ratio->GetXaxis()->SetLabelSize(18);
  hist_ratio->GetXaxis()->SetTitle("");
  //hist_ratio->GetYaxis()->SetTitle("UnfoData/Gen");
  hist_ratio->GetYaxis()->SetTitleFont(43);
  hist_ratio->GetYaxis()->SetTitleSize(15);
  hist_ratio->GetYaxis()->SetTitleOffset(2);
  hist_ratio->GetYaxis()->SetLabelFont(43);
  hist_ratio->GetYaxis()->SetLabelSize(15);
  hist_ratio->GetYaxis()->SetNdivisions(5);


  return hist_ratio;
}


