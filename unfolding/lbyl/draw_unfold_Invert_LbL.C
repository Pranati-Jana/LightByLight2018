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
double syst_uncertainty = 0.23;
double syst_uncertainty_mc = 0.23;

void make_hist(TH1D *&, string, string, int, float, Color_t, int, int) ;

TH1D* get_stat_uncertainty_hist(TH1D *input_hist, double uncertainty = syst_uncertainty);
TH1D* get_stat_uncertainty_hist_ratio(TH1D *input_hist);
void printBinCont(TH1D *hist);
void prepare_canvas(TCanvas *canvas);
TH1D* get_ratio(TH1D *hist_num, TH1D* hist_den);

int mc_gen_color = kBlue;
int data_unfold_color = kRed;
int mc_reco_color = kBlack;
int data_reco_color = kOrange+1;
int data_unfold_invert_color = kRed;
void draw_unfold_Invert_LbL(int method =1)
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

  TFile *f1= new TFile(Form("unfolding_histograms_Bayes_%s.root",mc_name.c_str()),"r");
  cout<<"   "<<f1->GetName()<<endl;
  TFile *outf= new TFile(Form("After_Unfolding_LByL.root"),"recreate"); 
  
  //TFile *f2= new TFile(Form("unfolding_histograms_invmass_%s_%s.root",MethID1,mc_name.c_str()),"r");
  //cout<<"   "<<f2->GetName()<<endl;
  
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
    
  
  TH1D *hist_invmass_gen          = (TH1D*)f1->Get("hGenInvmass_xSec");
  TH1D *hist_invmass_recomc       = (TH1D*)f1->Get("hRecoMCInvmass_xSec");
  TH1D *hist_invmass_recodata     = (TH1D*)f1->Get("hRecoDataInvmass_xSec");
  TH1D *hist_invmass_unfodata     = (TH1D*)f1->Get("hUnfoDataInvmass_xSec");
  TH1D *hist_invmass_unfomc       = (TH1D*)f1->Get("hUnfoMCInvmass_xSec");
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
  leg1->SetTextSize(15);
  leg1->SetHeader("#gamma#gamma#rightarrow #gamma#gamma");
  leg1->AddEntry(hist_pt_unfomc,"Unfolded MC (Bayes,iteration 4)","pl");
  leg1->AddEntry(hist_pt_gen, "Gen-level MC (Superchic)","pl");
  leg1->AddEntry(hist_pt_recomc, "Reconstructed MC (Superchic 3.02)", "pl");
  //leg1->AddEntry(hist_pt_unfomc,"Unfolded MC(Bayes,)","pl");
  
  TLegend *leg2=new TLegend(0.30,0.68,0.90,0.9);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(15);
  leg2->SetHeader("#gamma#gamma#rightarrow #gamma#gamma");
  leg2->AddEntry(hist_pt_recodata,"Reconstructed data","pl");
  //leg2->AddEntry(hist_pt_unfodata,"Unfolded data (RooBayes)","pl");
  //leg2->AddEntry(hist_rapidity_unfodata_Invert,"Unfolded data (Matrix Inversion)","pl");
  leg2->AddEntry(hist_pt_unfodata,"Unfolded data (Bayes,1 iteration)","pl");
  leg2->AddEntry(hist_pt_recomc, "Reconstructed MC (Superchic 3.02)", "pl");
  leg2->AddEntry(hist_pt_gen, "Gen-level MC (Superchic 3.02)", "pl");
  
  
  

  TLegend *leg3=new TLegend(0.30,0.68,0.90,0.9);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextFont(43);
  leg3->SetTextSize(15);
  leg3->SetHeader("#gamma#gamma#rightarrow #gamma#gamma");
  leg3->AddEntry(hist_rapidity_recodata,"Reconstructed data","pl");
  //leg3->AddEntry(hist_pt_unfodata,"Unfolded data (Bayes)","pl");
  //leg2->AddEntry(hist_rapidity_unfodata_Invert,"Unfolded data (Matrix Inversion)","pl");
  leg3->AddEntry(hist_rapidity_unfodata,"Unfolded data (Matrix Inversion)","pl");
  leg3->AddEntry(hist_rapidity_recomc, "Reconstructed MC (Superchic 3.02)", "pl");
  leg3->AddEntry(hist_rapidity_gen, "Gen-level MC (Superchic 3.02)", "pl");

 TLegend *leg4=new TLegend(0.30,0.72,0.92,0.95);
  leg4->SetFillColor(0);
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetTextFont(43);
  leg4->SetTextSize(15);
  leg4->SetHeader("#gamma#gamma#rightarrow #gamma#gamma");
  leg4->AddEntry(hist_rapidity_unfomc,"Unfolded MC (Matrix Inversion)","pl");
  leg4->AddEntry(hist_rapidity_gen, "Gen-level MC (Superchic 3.02)","pl");
  leg4->AddEntry(hist_rapidity_recomc, "Reconstructed MC (Superchic 3.02)", "pl");
 

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
       mark->DrawLatex(0.55,startY,"#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
       //mark->Draw();
  
  // pT
  TCanvas* c1 = new TCanvas("c1","DiPhoton pT",50,50,W,H);
  prepare_canvas(c1);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->SetRightMargin(0.08); 
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  string x_title = "p_{T}^{#gamma#gamma} (GeV)";
  string y_title = "#frac{d#sigma_{#gamma#gamma}}{dp_{T}^{#gamma#gamma}} (nb/0.2GeV)";
  
  cout<<"pt bin width:" <<hist_pt_unfodata->GetXaxis()->GetBinWidth(1)<<endl;
  
  make_hist(hist_pt_unfomc,x_title,y_title,20,1.0,data_unfold_color,2,2);
  make_hist(hist_pt_unfomc_Invert,x_title,y_title,2,1.0,data_unfold_invert_color,2,2);
  make_hist(hist_pt_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  auto hist_pt_gen_syst = get_stat_uncertainty_hist(hist_pt_unfomc, syst_uncertainty_mc);
  auto hist_pt_recomc_syst = get_stat_uncertainty_hist(hist_pt_recomc, syst_uncertainty_mc);
  make_hist(hist_pt_recomc_syst,x_title,y_title,0,1.0,mc_reco_color,1,2);
  
  hist_pt_gen_syst->SetMaximum(2*pow(10.,3.));
  hist_pt_gen_syst->SetMinimum(0.25);//5*pow(10.,-1.
   // hist_pt_gen_syst->GetYaxis()->SetRangeUser(0,700);    
  hist_pt_gen_syst->Draw("e2");
  hist_pt_gen->Draw("histsamex0");
  hist_pt_unfomc->Draw("psame");
 // hist_pt_unfomc_Invert->Draw("psame");
  make_hist(hist_pt_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hist_pt_recomc->Draw("ehistsamex0");
  hist_pt_recomc_syst->Draw("e2same");
  leg1->Draw();
  c1->cd();
  TLatex *mark5 = new TLatex();
  mark5->SetTextSize(0.035); 
  mark5->SetTextFont(42);
   mark5->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark5->DrawLatex(0.55, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
//  mark->Draw();
  //c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.06, 1, 0.28);
  pad2->SetTopMargin(0.08);
  pad2->SetBottomMargin(0.3);
  pad2->SetRightMargin(0.08); 
  //pad2->SetGridx();
  pad2->Draw();
  pad2->cd();
//  TH1D *hist_pt_closure_ratio = get_ratio(hist_pt_unfomc, hist_pt_gen);
  TH1D *hist_pt_closure_invert_ratio = get_ratio(hist_pt_unfomc_Invert, hist_pt_gen);
   TH1D *hist_pt_closure_ratio = get_ratio(hist_pt_unfomc, hist_pt_gen);
  auto hist_pt_closure_ratio_syst = get_stat_uncertainty_hist_ratio(hist_pt_gen_syst);
  make_hist(hist_pt_closure_ratio_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);
  hist_pt_closure_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_pt_closure_ratio->GetXaxis()->SetTitle("p_{T}^{#gamma#gamma} (GeV)");
  //hist_pt_closure_ratio->Draw("pex0");
//  hist_pt_closure_ratio_syst->Draw("e2same");
//  hist_pt_closure_invert_ratio->Draw("pex0");
   hist_pt_closure_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
   hist_pt_closure_ratio->Draw("pex0same");
   hist_pt_closure_ratio_syst->Draw("e2same");

  line->DrawCopy("same");
  
  //c1->Print("fig_Bayes_Inverse/closure_pt_ZDCAND3n_SCSL_17thMarch.png");
  c1->Print("fig_Bayes_Inverse/closure_pt_ZDCAND3n_SCSL_newSF.pdf");
  
  
  TCanvas* c2 = new TCanvas("c2","Diphoton pT",50,50,W,H);
  prepare_canvas(c2);
  TPad *pad3 = new TPad("pad3", "pad3", 0, 0.3, 1, 1.0);
  pad3->SetBottomMargin(0);
  pad3->SetRightMargin(0.08); 
  pad3->Draw();
  pad3->cd();
  pad3->SetLogy();  
  hist_pt_unfodata->SetMaximum(2*pow(10.,3.));
  
  
  make_hist(hist_pt_unfodata,x_title,y_title,20,1.0,data_unfold_color,2,2);
  make_hist(hist_pt_unfodata_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2);

  auto hist_pt_unfodata_syst = get_stat_uncertainty_hist(hist_pt_unfodata);
  make_hist(hist_pt_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  auto hist_pt_recodata_syst = get_stat_uncertainty_hist(hist_pt_recodata);
  make_hist(hist_pt_recodata_syst,x_title,y_title,0,1.0,data_reco_color,2,2);
hist_pt_unfodata_syst->SetMinimum(0.25);
  hist_pt_unfodata_syst->Draw("e2");
  hist_pt_gen->Draw("histsamex0");
 //hist_pt_gen_syst->Draw("e2same");
  //hist_pt_unfodata_syst->GetYaxis()->SetRangeUser(0,700); 
  hist_pt_unfodata->Draw("psameex0");
  hist_pt_unfodata_Invert->Draw("psameex0");
  hist_pt_recodata_syst->Draw("e2same");

  hist_pt_recomc->Draw("ehistsamex0");
   hist_pt_recomc_syst->Draw("e2same");
  
  hist_pt_recodata->Draw("ehistsamex0");
  
  //CMS_lumi( c2, 1, 10 );
  leg2->Draw();
  c2->cd();
  TLatex *mark6 = new TLatex();
  mark6->SetTextSize(0.035); 
  mark6->SetTextFont(42);
  
   mark6->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark6->DrawLatex(0.55, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad4 = new TPad("pad4", "pad4", 0, 0.06, 1, 0.28);
  pad4->SetTopMargin(0.08);
  pad4->SetBottomMargin(0.3);
  pad4->SetRightMargin(0.08); 
 // pad4->SetGridx();
  pad4->Draw();
  pad4->cd();
  TH1D *hist_pt_data_ratio = get_ratio(hist_pt_unfodata, hist_pt_gen);
  TH1D *hist_pt_data_invert_ratio = get_ratio(hist_pt_unfodata_Invert, hist_pt_gen);
  auto hist_pt_data_ratio_syst = get_stat_uncertainty_hist_ratio(hist_pt_unfodata_syst);
  make_hist(hist_pt_data_ratio_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);
  hist_pt_data_invert_ratio->GetXaxis()->SetTitle("p_{T}^{#gamma#gamma} (GeV)");
  hist_pt_data_invert_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hist_pt_data_invert_ratio->Draw("pex0same");
  hist_pt_data_ratio->Draw("pex0same");
  hist_pt_data_ratio_syst->Draw("e2same");
  line->DrawCopy("same");
  c2->Print("fig_Bayes_Inverse/data_pt_ZDCAND3n_SCSL_newSF.pdf");
 // c2->Print("fig_Bayes_Inverse/data_pt_ZDCAND3n_SCSL.png");
  //c2->Print("fig_Bayes_Inverse/data_pt_ZDCAND3n.C");

  
  // Rapidity
  TCanvas* c3 = new TCanvas("c3","Dimuon rapidity",50,50,W,H);
  c3->SetLeftMargin(0.15);
  prepare_canvas(c3);
  TPad *pad5 = new TPad("pad5", "pad5", 0, 0.3, 1, 1.0);
  pad5->SetBottomMargin(0);
  pad5->SetRightMargin(0.08); 
  pad5->Draw();
  pad5->cd();
  
  x_title = "y_{#gamma#gamma}";
  y_title = "#frac{d#sigma_{#gamma#gamma}}{dy_{#gamma#gamma}} (nb / 0.147)";
  
  //hist_rapidity_unfomc->GetXaxis()->SetRangeUser(0,1.0);
  make_hist(hist_rapidity_unfomc,x_title,y_title,20,1.0,data_unfold_color,2,2);
  make_hist(hist_rapidity_unfomc_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2);
  
  make_hist(hist_rapidity_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);

  auto hist_rapidity_gen_syst = get_stat_uncertainty_hist(hist_rapidity_unfomc, syst_uncertainty_mc);
  make_hist(hist_rapidity_gen_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);

  auto hist_rapidity_recomc_syst = get_stat_uncertainty_hist(hist_rapidity_recomc, syst_uncertainty_mc);
  make_hist(hist_rapidity_recomc_syst,x_title,y_title,0,1.0,mc_reco_color,1,2);
  
  hist_rapidity_gen_syst->SetMaximum(60);//100
  hist_rapidity_gen_syst->SetMinimum(0);
  
  hist_rapidity_gen_syst->Draw("e2");
  hist_rapidity_gen->Draw("histsame");
  hist_rapidity_gen->Draw("psame");
  hist_rapidity_unfomc->Draw("psame");
  hist_rapidity_unfomc_Invert->Draw("psame");
  make_hist(hist_rapidity_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hist_rapidity_recomc->Draw("ehistsamex0");
  hist_rapidity_recomc_syst->Draw("e2same");
  leg4->Draw();
  c3->cd();
  TLatex *mark7 = new TLatex();
  mark7->SetTextSize(0.035); 
  mark7->SetTextFont(42);
  
   mark7->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark7->DrawLatex(0.55, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad6 = new TPad("pad6", "pad6", 0, 0.06, 1, 0.28);
  pad6->SetTopMargin(0.08);
  pad6->SetBottomMargin(0.3);
  pad6->SetRightMargin(0.08); 
//  pad6->SetGridx();
  pad6->Draw();
  pad6->cd();
  TH1D *hist_rapidity_closure_ratio = get_ratio(hist_rapidity_unfomc, hist_rapidity_gen);
  TH1D *hist_rapidity_closure_invert_ratio = get_ratio(hist_rapidity_unfomc_Invert, hist_rapidity_gen);
  auto hist_rapidity_closure_ratio_syst = get_stat_uncertainty_hist_ratio(hist_rapidity_gen_syst);
  make_hist(hist_rapidity_closure_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  hist_rapidity_closure_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_rapidity_closure_ratio->GetXaxis()->SetTitle("y_{#gamma#gamma}");
  hist_rapidity_closure_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hist_rapidity_closure_ratio->Draw("pex0");
  hist_rapidity_closure_ratio_syst->Draw("e2same");
  hist_rapidity_closure_invert_ratio->Draw("pex0same");
  line->DrawCopy("same");
  c3->Print("fig_Bayes_Inverse/closure_rapidity_ZDCAND3n_SCSL_newSF.pdf");
  //c3->Print("fig_Bayes_Inverse/closure_rapidity_ZDCAND3n_SCSL.png");


  
  TCanvas* c4 = new TCanvas("c4","Diphoton rap data",10,50,W,H);
  prepare_canvas(c4);
  TPad *pad7 = new TPad("pad7", "pad7", 0, 0.3, 1, 1.0);
  
  pad7->SetBottomMargin(0.02);
  pad7->SetRightMargin(0.08); 
  pad7->Draw();
  pad7->cd();
  TLatex *mark1 = new TLatex();
       mark1->SetNDC(true);
       double startY1 = 0.99;
       mark1->SetTextFont(42);
       mark1->SetTextSize(0.035);
       mark1->DrawLatex(0.25,startY1,"#bf{CMS} #it{Preliminary}");
       mark1->DrawLatex(0.55,startY1,"#scale[0.8]{1647.18 #mub^{-1} ( PbPb @ 5.02 TeV)}");
       mark1->Draw();
  hist_rapidity_unfodata->SetMaximum(60);//100
  hist_rapidity_unfodata->SetMinimum(0);
  hist_rapidity_unfodata->GetXaxis()->SetRangeUser(-2.2,2.2);
   make_hist(hist_rapidity_unfodata_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2);
  make_hist(hist_rapidity_unfodata,x_title,y_title,20,1.0,data_unfold_color,2,2);
    make_hist(hist_rapidity_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  auto hist_rapidity_unfodata_syst = get_stat_uncertainty_hist(hist_rapidity_unfodata);
  auto hist_rapidity_recodata_syst = get_stat_uncertainty_hist(hist_rapidity_recodata);
  make_hist(hist_rapidity_recodata_syst,x_title,y_title,0,1.0,data_reco_color,2,2);
   hist_rapidity_unfodata_syst->GetXaxis()->SetLabelSize(0);
  hist_rapidity_unfodata_syst->Draw("e2");
  make_hist(hist_rapidity_unfodata_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2);

  hist_rapidity_gen->Draw("histsamex0");
  hist_rapidity_gen->Draw("psame");
  //hist_rapidity_gen_syst->Draw("e2same");
  //hist_rapidity_unfodata_Invert->Draw("psameex0");
  hist_rapidity_unfodata->Draw("psameex0");
  hist_rapidity_recomc->Draw("ehistsamex0");

  hist_rapidity_recodata->Draw("ehistsamex0");
  hist_rapidity_unfodata_syst->Draw("e2same");
  hist_rapidity_recodata_syst->Draw("e2same");
  hist_rapidity_recomc_syst->Draw("e2same");
  mark1->Draw();
  leg3->Draw();
  c4->cd();
  TLatex *mark8 = new TLatex();
  mark8->SetTextSize(0.035); 
  mark8->SetTextFont(42);
 
   mark8->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark8->DrawLatex(0.55, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad8 = new TPad("pad8", "pad8", 0, 0.06, 1, 0.28);
  pad8->SetTopMargin(0.08);
  pad8->SetBottomMargin(0.31);
  pad8->SetRightMargin(0.08); 
  //pad8->SetGridx();
  pad8->Draw();
  pad8->cd();
  TH1D *hist_rapidity_data_ratio = get_ratio(hist_rapidity_unfodata, hist_rapidity_gen);
  TH1D *hist_rapidity_data_invert_ratio = get_ratio(hist_rapidity_unfodata_Invert, hist_rapidity_gen);
  auto hist_rapidity_data_ratio_syst = get_stat_uncertainty_hist_ratio(hist_rapidity_unfodata_syst);
  make_hist(hist_rapidity_data_ratio_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);
  hist_rapidity_data_ratio->GetXaxis()->SetTitle("y_{#gamma#gamma}");
  hist_rapidity_data_ratio_syst->GetXaxis()->SetTitle("y_{#gamma#gamma}");
  hist_rapidity_data_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
 
  hist_rapidity_data_ratio->Draw("pex0");
 // hist_rapidity_data_invert_ratio->Draw("pex0same");
  hist_rapidity_data_ratio->SetMaximum(3.0);
  hist_rapidity_data_ratio->Draw("pex0same");
  hist_rapidity_data_ratio_syst->Draw("e2same");
  mark1->Draw("same");
  line->DrawCopy("same");
  mark1->Draw();
  c4->Print("fig_Bayes_Inverse/data_rapidity_ZDCAND3n_SCSL_newSF.pdf");
  // c4->Print("fig_Bayes_Inverse/data_rapidity_ZDCAND3n_SCSL.png");
  // c4->Print("fig_Bayes_Inverse/data_rapidity_ZDCAND3n.C");
  
  //invmass
  TCanvas* c5 = new TCanvas("c5","Diphoton invmass gen",50,50,W,H);
  prepare_canvas(c5);
  TPad *pad9 = new TPad("pad9", "pad9", 0, 0.3, 1, 1.0);
  pad9->SetBottomMargin(0);
  pad9->SetRightMargin(0.08); 
  pad9->Draw();
  pad9->cd();
  gPad->SetLogy();
 // c5->cd();
  x_title = "m_{#gamma#gamma} (GeV)";
  y_title = "#frac{d#sigma_{#gamma#gamma}}{dm_{#gamma#gamma}} (nb/4GeV)";
  cout<<"m_inv bin width:" <<hist_invmass_unfodata->GetXaxis()->GetBinWidth(1)<<endl;
  
  make_hist(hist_invmass_unfomc,x_title,y_title,20,1.0,data_unfold_color,2,2);
  make_hist(hist_invmass_unfomc_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2);
  make_hist(hist_invmass_gen,x_title,y_title,2,1.0,mc_gen_color,1,2);
  auto hist_invmass_recomc_syst = get_stat_uncertainty_hist(hist_invmass_recomc, syst_uncertainty_mc);
  make_hist(hist_invmass_recomc_syst,x_title,y_title,0,1.0,mc_reco_color,1,2);
  
  hist_invmass_gen->SetMaximum(100);//50
  hist_invmass_gen->SetMinimum(0.015);
  hist_invmass_gen->Draw("histsame");
  hist_invmass_gen->Draw("psame");
  auto hist_invmass_gen_syst = get_stat_uncertainty_hist(hist_invmass_unfomc, syst_uncertainty_mc);
  hist_invmass_gen_syst->Draw("e2same");
  make_hist(hist_invmass_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hist_invmass_recomc->Draw("ehistsamex0");
  hist_invmass_unfomc->Draw("psamex0");
  hist_invmass_unfomc_Invert->Draw("psamex0");
  hist_invmass_recomc_syst->Draw("e2same");

  leg4->Draw();
  c5->cd();
  TLatex *mark9 = new TLatex();
  mark9->SetTextSize(0.035); 
  mark9->SetTextFont(42);

   mark9->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark9->DrawLatex(0.55, 0.97, "#scale[0.8]{1.647 nb^{-1} ( PbPb @ 5.02 TeV)}");
  TPad *pad10 = new TPad("pad10", "pad10", 0, 0.06, 1, 0.28);
  pad10->SetTopMargin(0.08);
  pad10->SetBottomMargin(0.3);
  pad10->SetRightMargin(0.08); 
 // pad10->SetGridx();
  pad10->Draw();
  pad10->cd();
  TH1D *hist_invmass_closure_ratio = get_ratio(hist_invmass_unfomc, hist_invmass_gen);
  TH1D *hist_invmass_closure_invert_ratio = get_ratio(hist_invmass_unfomc_Invert, hist_invmass_gen);
  auto hist_invmass_closure_ratio_syst = get_stat_uncertainty_hist_ratio(hist_invmass_gen_syst);
   make_hist(hist_invmass_closure_ratio_syst,x_title,y_title,0,1.0,mc_gen_color,2,2);
  hist_invmass_closure_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_invmass_closure_ratio->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
  hist_invmass_closure_ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  hist_invmass_closure_ratio->Draw("pex0");
  hist_invmass_closure_ratio_syst->Draw("e2same");
  hist_invmass_closure_invert_ratio->Draw("pex0same");
  line->DrawCopy("same");
  c5->Print("fig_Bayes_Inverse/closure_invmass_ZDCAND3n_SCSL_newSF.pdf");
 // c5->Print("fig_Bayes_Inverse/closure_invmass_ZDCAND3n_SCSL.png");
  
  TCanvas* c6 = new TCanvas("c6","Diphoton mass",50,50,W,H);
  prepare_canvas(c6);
  TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
  pad11->SetBottomMargin(0.0);
  pad11->SetRightMargin(0.08); 
  pad11->Draw();
  pad11->cd();
  gPad->SetLogy();
   TLatex *mark2 = new TLatex();
       mark2->SetNDC(true);
       double startY2 = 0.90;

       mark2->SetTextFont(42);
       mark2->SetTextSize(0.035);
       mark2->DrawLatex(0.15,startY2,"#bf{CMS} #it{Preliminary}");
       mark2->DrawLatex(0.55,startY2,"#scale[0.8]{1647.18 #mub^{-1} (PbPb @ 5.02 TeV)}");
       mark2->Draw();
  hist_invmass_unfodata->SetMaximum(100);//50
  hist_invmass_unfodata->SetMinimum(0.015);
  make_hist(hist_invmass_unfodata,x_title,y_title,20,1.0,data_unfold_color,2,2);
  auto hist_invmass_unfodata_syst = get_stat_uncertainty_hist(hist_invmass_unfodata);
  hist_invmass_unfodata_syst->Draw("e2");
    make_hist(hist_invmass_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  make_hist(hist_invmass_unfodata_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2); 
   auto hist_invmass_recodata_syst = get_stat_uncertainty_hist(hist_invmass_recodata);
  make_hist(hist_invmass_recodata_syst,x_title,y_title,0,1.0,data_reco_color,2,2);

  hist_invmass_unfodata->Draw("psameex0");
   hist_invmass_gen->Draw("histsamex0");
   hist_invmass_gen->Draw("psame");
   // hist_invmass_gen_syst->Draw("e2same");
    hist_invmass_recomc->Draw("ehistsamex0");
    hist_invmass_recodata_syst->Draw("e2same");
// hist_invmass_unfodata_Invert->Draw("psameex0");
   hist_invmass_unfodata->Draw("psameex0");
 // hist_invmass_unfodata_Invert->Draw("psameex0");

  hist_invmass_recodata->Draw("ehistsamex0");
  hist_invmass_recomc_syst->Draw("e2same");
  leg3->Draw();
  
  c6->cd();
  TLatex *mark10 = new TLatex();
  mark10->SetTextSize(0.035); 
  mark10->SetTextFont(42);
  
   mark10->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark10->DrawLatex(0.55, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
  TPad *pad12 = new TPad("pad12", "pad12", 0, 0.06, 1, 0.28);
  pad12->SetTopMargin(0.08);
  pad12->SetBottomMargin(0.3);
 // pad12->SetGridx();
  pad12->Draw();
  pad12->cd();
  pad12->SetRightMargin(0.08);
  TH1D *hist_invmass_data_ratio = get_ratio(hist_invmass_unfodata, hist_invmass_gen);
  TH1D *hist_invmass_data_invert_ratio = get_ratio(hist_invmass_unfodata_Invert, hist_invmass_gen);
  auto hist_invmass_data_ratio_syst = get_stat_uncertainty_hist_ratio(hist_invmass_unfodata_syst);
  make_hist(hist_invmass_data_ratio_syst,x_title,y_title,0,1.0,data_unfold_color,2,2);
  hist_invmass_data_ratio->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
  hist_invmass_data_ratio->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
  hist_invmass_data_ratio->Draw("pex0");
  hist_invmass_data_ratio->SetMaximum(2.5);
//  hist_invmass_data_invert_ratio->Draw("pex0same");
  hist_invmass_data_ratio->Draw("pex0same");
  hist_invmass_data_ratio_syst->Draw("e2same");
   line->DrawCopy("same");
  c6->Print("fig_Bayes_Inverse/data_invmass_ZDCAND3n_SCSL_newSF.pdf");

  // c6->Print("fig_Bayes_Inverse/data_invmass_ZDCAND3n_SCSL.png");
  // c6->Print("fig_Bayes_Inverse/data_invmass_ZDCAND3n.C");
  
////Costhetastar
 TCanvas* c7 = new TCanvas("c7","Diphoton costhetastar gen",50,50,W,H);
  prepare_canvas(c7);
  TPad *pad13 = new TPad("pad13", "pad13", 0, 0.3, 1, 1.0);
  pad13->SetBottomMargin(0);
  pad13->SetRightMargin(0.08); 
  pad13->Draw();
  pad13->cd();
 // c5->cd();
  x_title = "m_{#gamma#gamma} (GeV)";
  y_title = "d#sigma_{#gamma#gamma}/d|cos#theta*|_{#gamma#gamma} (nb/binwidth)";
  cout<<"costhetastar bin width:" <<hist_costhetastar_unfodata->GetXaxis()->GetBinWidth(1)<<endl;
  
  make_hist(hist_costhetastar_unfomc,x_title,y_title,20,1.0,data_unfold_color,2,2);
  make_hist(hist_costhetastar_unfomc_Invert,x_title,y_title,20,1.0,data_unfold_invert_color,2,2);
  make_hist(hist_costhetastar_gen,x_title,y_title,20,1.0,mc_gen_color,1,2);
  
  hist_costhetastar_gen->SetMaximum(200);//50
  hist_costhetastar_gen->SetMinimum(0);
  hist_costhetastar_gen->Draw("histsame");
  
  auto hist_costhetastar_gen_syst = get_stat_uncertainty_hist(hist_costhetastar_gen, syst_uncertainty_mc);
  hist_costhetastar_gen_syst->Draw("e2same");
  
  make_hist(hist_costhetastar_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hist_costhetastar_recomc->Draw("ehistsamex0");
  hist_costhetastar_unfomc->Draw("psamex0");
  hist_costhetastar_unfomc_Invert->Draw("psamex0");
  
  leg4->Draw();
  c7->cd();
  TLatex *mark11 = new TLatex();
  mark11->SetTextSize(0.035); 
  mark11->SetTextFont(42);

   mark11->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
   mark11->DrawLatex(0.55, 0.97, "#scale[0.8]{1.647 nb^{-1} ( PbPb @ 5.02 TeV)}");
  TPad *pad14 = new TPad("pad11", "pad10", 0, 0.06, 1, 0.28);
  pad14->SetTopMargin(0);
  pad14->SetBottomMargin(0.3);
  pad14->SetRightMargin(0.08); 

  pad14->Draw();
  pad14->cd();
  TH1D *hist_costhetastar_closure_ratio = get_ratio(hist_costhetastar_unfomc, hist_costhetastar_gen);
  TH1D *hist_costhetastar_closure_invert_ratio = get_ratio(hist_costhetastar_unfomc_Invert, hist_costhetastar_gen);
  auto hist_costhetastar_closure_ratio_syst = get_stat_uncertainty_hist(hist_costhetastar_closure_ratio, syst_uncertainty_mc);
  hist_costhetastar_closure_ratio->GetYaxis()->SetTitle("UnfoMC/Gen");
  hist_costhetastar_closure_ratio->GetXaxis()->SetTitle("|cos#theta^{*}|");
  hist_costhetastar_closure_ratio->GetYaxis()->SetRangeUser(0.0,3.0);
 
  hist_costhetastar_closure_ratio->Draw("pex0");
  hist_costhetastar_closure_ratio_syst->Draw("e2same");
  hist_costhetastar_closure_invert_ratio->Draw("pex0same");
  line->DrawCopy("same");
  c7->Print("fig_Bayes_Inverse/closure_costhetastar_ZDCAND3n_newSF.pdf");
 // c7->Print("fig_Bayes_Inverse/closure_costhetastar_ZDCAND3n_SCSL.png");
 ////
  TCanvas* c8 = new TCanvas("c8","Diphoton costhetastar",50,50,W,H);
  prepare_canvas(c8);
  TPad *pad15 = new TPad("pad15", "pad15", 0, 0.3, 1, 1.0);
  pad15->SetBottomMargin(0.0);
  pad15->SetRightMargin(0.08); 
  pad15->Draw();
  pad15->cd();
  hist_costhetastar_unfodata->SetMaximum(250);//50
  hist_costhetastar_unfodata->SetMinimum(0);
  make_hist(hist_costhetastar_unfodata,x_title,y_title,20,1.0,data_unfold_color,2,2);
  auto hist_costhetastar_unfodata_syst = get_stat_uncertainty_hist(hist_costhetastar_unfodata);
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
   mark12->DrawLatex(0.55, 0.97, "#scale[0.8]{1.647 nb^{-1} (PbPb @ 5.02 TeV)}");
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
  c8->Print("fig_Bayes_Inverse/data_costhetatar_ZDCAND3n_SCSL_newSF.pdf");
  //c8->Print("fig_Bayes_Inverse/data_costhetastar_ZDCAND3n_SCSL.png");
 // c8->Print("fig_Bayes_Inverse/data_costhetastar_ZDCAND3n.C");
  //////////////
  //Comparison with GammaUPC
  



  
  ///////
  printBinCont(hist_pt_unfodata);
  printBinCont(hist_rapidity_unfodata);
  printBinCont(hist_invmass_unfodata);
  printBinCont(hist_costhetastar_unfodata);

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
  hist->GetXaxis()->SetTitleOffset(0.9);
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
  
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.8);
  output_hist->SetFillStyle(3013); 
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
  //     else{
  //       output_hist->SetBinContent(i,1);
      // }
  }
  
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.8);
  output_hist->SetFillStyle(3013); 
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
  hist_ratio->GetYaxis()->SetTitleSize(20);
  hist_ratio->GetYaxis()->SetTitleOffset(2);
  hist_ratio->GetYaxis()->SetLabelFont(43);
  hist_ratio->GetYaxis()->SetLabelSize(20);
  hist_ratio->GetYaxis()->SetNdivisions(5);
  

  return hist_ratio;
}

