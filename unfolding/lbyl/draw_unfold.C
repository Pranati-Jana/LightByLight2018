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

double syst_uncertainty = 0.26;
double syst_uncertainty_mc = 0.10;

void make_hist(TH1D *&, string, string, int, float, Color_t, int, int) ;
TH1D* get_stat_uncertainty_hist(TH1D *input_hist, double uncertainty = syst_uncertainty);
void printBinCont(TH1D *hist);
void prepare_canvas(TCanvas *canvas);

int mc_gen_color = kBlue;
int data_unfold_color = kRed;
int mc_reco_color = kBlack;
int data_reco_color = kOrange+1;

void draw_unfold(int method =1)
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
    
//    sprintf(MethID1,"Bayes_unfo");
  }
  if(method==2)
  {
  //  sprintf(MethID1,"Svd_unfo ");
  }
  if(method==3)
  {
    //sprintf(MethID1,"BinByBin_unfo");
  }
  
  string mc_name = "SC";

  TFile *f1= new TFile(Form("unfolding_histograms_Bayes_%s.root",mc_name.c_str()),"r");
  cout<<"   "<<f1->GetName()<<endl;
  
  
  //TFile *f2= new TFile(Form("unfolding_histograms_invmass_%s_%s.root",MethID1,mc_name.c_str()),"r");
  //cout<<"   "<<f2->GetName()<<endl;
  
  TH1D *hMuMu_pt_gen          = (TH1D*)f1->Get("hGenPt_xSec");
  TH1D *hMuMu_pt_recomc       = (TH1D*)f1->Get("hRecoMCPt_xSec");
  TH1D *hMuMu_pt_recodata     = (TH1D*)f1->Get("hRecoDataPt_xSec");
  TH1D *hMuMu_pt_unfodata     = (TH1D*)f1->Get("hUnfoDataPt_xSec");
  TH1D *hMuMu_pt_unfomc       = (TH1D*)f1->Get("hUnfoMCPt_xSec");
  
  TH1D *hMuMu_rapidity_gen          = (TH1D*)f1->Get("hGenRap_xSec");
  TH1D *hMuMu_rapidity_recomc       = (TH1D*)f1->Get("hRecoMCRap_xSec");
  TH1D *hMuMu_rapidity_recodata     = (TH1D*)f1->Get("hRecoDataRap_xSec");
  TH1D *hMuMu_rapidity_unfodata     = (TH1D*)f1->Get("hUnfoDataRap_xSec");
  TH1D *hMuMu_rapidity_unfomc       = (TH1D*)f1->Get("hUnfoMCRap_xSec");
  
  
  TH1D *hMuMu_invmass_gen          = (TH1D*)f1->Get("hGenInvmass_xSec");
  TH1D *hMuMu_invmass_recomc       = (TH1D*)f1->Get("hRecoMCInvmass_xSec");
  TH1D *hMuMu_invmass_recodata     = (TH1D*)f1->Get("hRecoDataInvmass_xSec");
  TH1D *hMuMu_invmass_unfodata     = (TH1D*)f1->Get("hUnfoDataInvmass_xSec");
  TH1D *hMuMu_invmass_unfomc       = (TH1D*)f1->Get("hUnfoMCInvmass_xSec");
  
  cout<<"ok"<<endl;
  
  TLegend *leg1=new TLegend(0.30,0.72,0.90,0.9);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(43);
  leg1->SetTextSize(24);
  leg1->AddEntry(hMuMu_pt_unfomc,"Unfolded MC","pl");
  leg1->AddEntry(hMuMu_pt_gen, "Gen-level MC (Superchic)","pl");
  leg1->AddEntry(hMuMu_pt_recomc, "Reconstructed MC (Superchic)", "pl");
  
  
  TLegend *leg2=new TLegend(0.30,0.68,0.90,0.9);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(24);
  
  leg2->AddEntry(hMuMu_pt_recodata,"Reconstructed data","pl");
  leg2->AddEntry(hMuMu_pt_unfodata,"Unfolded data","pl");
  leg2->AddEntry(hMuMu_pt_recomc, "Reconstructed MC (Superchic)", "pl");
  leg2->AddEntry(hMuMu_pt_gen, "Gen-level MC (Superchic)", "pl");
  
  int W = 600;
  int H = 670;//500
  
  // pT
  TCanvas* c1 = new TCanvas("c1","Dimuon pT",50,50,W,H);
  prepare_canvas(c1);
//  c1->cd();
//  gPad->SetLogy();//Diphototn pT
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  string x_title = "p_{T}^{#gamma#gamma} (GeV)";
  string y_title = "d#sigma_{#gamma#gamma}/dp_{T}^{#gamma#gamma} (nb)";
  
  cout<<"pt bin width:" <<hMuMu_pt_unfodata->GetXaxis()->GetBinWidth(1)<<endl;
  
  make_hist(hMuMu_pt_unfomc,x_title,y_title,20,1.0,data_unfold_color,1,1);
  
  make_hist(hMuMu_pt_gen,x_title,y_title,20,1.0,mc_gen_color,1,2);
  auto hMuMu_pt_gen_syst = get_stat_uncertainty_hist(hMuMu_pt_gen, syst_uncertainty_mc);
  
  hMuMu_pt_gen_syst->SetMaximum(2*pow(10.,3.));
  hMuMu_pt_gen_syst->SetMinimum(5*pow(10.,-2.));
  
  hMuMu_pt_gen_syst->Draw("e2");
  hMuMu_pt_gen->Draw("histsamex0");
  hMuMu_pt_unfomc->Draw("psame");
  make_hist(hMuMu_pt_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hMuMu_pt_recomc->Draw("ehistsamex0");
  leg1->Draw();

  c1->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.06, 1, 0.28);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad2->SetGridx();
  pad2->Draw();
  pad2->cd();
  TString strRatio=hMuMu_pt_gen->GetName();
  strRatio+="_over_";
  strRatio+=hMuMu_pt_gen->GetName();
  TH1D* hratio=(TH1D*)hMuMu_pt_unfomc->Clone(strRatio);
  hratio->Sumw2();

  hratio->Divide(hMuMu_pt_gen);
  hratio->SetTitle("");
  //hratio->GetXaxis()->SetTitle(ytitle);
  hratio->GetYaxis()->SetTitle("Unfo/Gen");
  hratio->SetLabelSize(0.1, "XYZ");
  hratio->SetLabelFont(42, "XYZ");
  hratio->SetLabelOffset(0.007, "XYZ");
  hratio->SetTitleSize(0.11, "XYZ");
  hratio->GetXaxis()->SetTitleOffset(1.0);
  hratio->GetYaxis()->SetTitleOffset(0.3);
  hratio->Draw("p");
  // c1->Update();
  //CMS_lumi( c1, 1, 10 );
  //leg1->Draw();
  c1->Print("7thNov/closure_pt.png");
  c1->Print("7thNov/closure_pt.pdf");
  
  
  TCanvas* c2 = new TCanvas("c2","Dimuon pT",50,50,W,H);
  prepare_canvas(c2);
 // c2->cd();
  //gPad->SetLogy();
  //hMuMu_pt_unfodata->SetMaximum(60);
  //hMuMu_pt_unfodata->SetMinimum(0);
  //gPad->SetLogy();
  TPad *pad3 = new TPad("pad3", "pad3", 0, 0.3, 1, 1.0);
  pad3->SetBottomMargin(0);
  pad3->Draw();
  pad3->cd();
  pad3->SetLogy();  
  hMuMu_pt_unfodata->SetMaximum(2*pow(10.,3.));
  hMuMu_pt_unfodata->SetMinimum(5*pow(10.,-2.));
  
  make_hist(hMuMu_pt_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  auto hMuMu_pt_unfodata_syst = get_stat_uncertainty_hist(hMuMu_pt_unfodata);
  hMuMu_pt_unfodata_syst->Draw("e2");
//  hMuMu_pt_unfodata_syst->GetYaxis()->SetRangeUser(0,500); 
  hMuMu_pt_unfodata->Draw("psameex0");
  hMuMu_pt_gen->Draw("histsamex0");
  hMuMu_pt_gen_syst->Draw("e2same");
  hMuMu_pt_recomc->Draw("ehistsamex0");

  make_hist(hMuMu_pt_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  hMuMu_pt_recodata->Draw("ehistsamex0");
  
  //CMS_lumi( c2, 1, 10 );
  leg2->Draw();
  c2->cd();
  TPad *pad4 = new TPad("pad4", "pad4", 0, 0.06, 1, 0.28);
  pad4->SetTopMargin(0);
  pad4->SetBottomMargin(0.3);
  pad4->SetGridx();
  pad4->Draw();
  pad4->cd();
  TString strRatio2=hMuMu_pt_gen->GetName();
  strRatio2+="_over_";
  strRatio2+=hMuMu_pt_gen->GetName();
  TH1D* hratio2=(TH1D*)hMuMu_pt_unfodata->Clone(strRatio2);
  hratio2->Sumw2();

  hratio2->Divide(hMuMu_pt_gen);
  hratio2->SetTitle("");
  hratio2->GetYaxis()->SetTitle("UnfoData/Gen");
  hratio2->SetLabelSize(0.1, "XYZ");
  hratio2->SetLabelFont(42, "XYZ");
  hratio2->SetLabelOffset(0.007, "XYZ");
  hratio2->SetTitleSize(0.11, "XYZ");
  hratio2->GetXaxis()->SetTitleOffset(1.0);
  hratio2->GetYaxis()->SetTitleOffset(0.3);
  hratio2->GetYaxis()->SetRangeUser(0.0, 2.0);
  hratio2->Draw("p");
  c2->Print("7thNov/data_pt.pdf");
  c2->Print("7thNov/data_pt.png");

  
  // Rapidity
  TCanvas* c3 = new TCanvas("c3","Dimuon rapidity",50,50,W,H);
  prepare_canvas(c3);
  TPad *pad5 = new TPad("pad5", "pad5", 0, 0.3, 1, 1.0);
  pad5->SetBottomMargin(0);
  pad5->Draw();
  pad5->cd();

 // c3->cd();
  
  x_title = "|y_{#gamma#gamma}|";
  y_title = "d#sigma_{#gamma#gamma}/dy_{#gamma#gamma} (nb)";
  
  
  //hMuMu_rapidity_unfomc->GetXaxis()->SetRangeUser(0,1.0);
  make_hist(hMuMu_rapidity_unfomc,x_title,y_title,20,1.0,data_unfold_color,1,1);
  
  make_hist(hMuMu_rapidity_gen,x_title,y_title,20,1.0,mc_gen_color,1,2);
  auto hMuMu_rapidity_gen_syst = get_stat_uncertainty_hist(hMuMu_rapidity_gen, syst_uncertainty_mc);
  
  hMuMu_rapidity_gen_syst->SetMaximum(60);//100
  hMuMu_rapidity_gen_syst->SetMinimum(0);
  
  hMuMu_rapidity_gen_syst->Draw("e2");
  hMuMu_rapidity_gen->Draw("histsame");
  hMuMu_rapidity_unfomc->Draw("psame");
  make_hist(hMuMu_rapidity_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hMuMu_rapidity_recomc->Draw("ehistsamex0");
  
  //CMS_lumi( c3, 1, 10 );
  leg1->Draw();
  c3->cd();
  TPad *pad6 = new TPad("pad6", "pad6", 0, 0.06, 1, 0.28);
  pad6->SetTopMargin(0);
  pad6->SetBottomMargin(0.3);
  pad6->SetGridx();
  pad6->Draw();
  pad6->cd();
  TString strRatio3=hMuMu_rapidity_gen->GetName();
  strRatio3+="_over_";
  strRatio3+=hMuMu_rapidity_gen->GetName();
  TH1D* hratio3=(TH1D*)hMuMu_rapidity_unfomc->Clone(strRatio3);
  hratio3->Sumw2();

  hratio3->Divide(hMuMu_rapidity_gen);
  hratio3->SetTitle("");
  hratio3->GetYaxis()->SetTitle("Unfo/Gen");
  hratio3->SetLabelSize(0.1, "XYZ");
  hratio3->SetLabelFont(42, "XYZ");
  hratio3->SetLabelOffset(0.007, "XYZ");
  hratio3->SetTitleSize(0.11, "XYZ");
  hratio3->GetXaxis()->SetTitleOffset(1.0);
  hratio3->GetYaxis()->SetTitleOffset(0.3);
  //hratio3->GetYaxis()->SetRangeUser(0.0, 2.0);
  hratio3->Draw("p");

  c3->Print("7thNov/closure_rapidity.pdf");
  c3->Print("7thNov/closure_rapidity.png");

  
  TCanvas* c4 = new TCanvas("c4","Dimuon rap",50,50,W,H);
  prepare_canvas(c4);
  TPad *pad7 = new TPad("pad7", "pad7", 0, 0.3, 1, 1.0);
  pad7->SetBottomMargin(0);
  pad7->Draw();
  pad7->cd();

//  c4->cd();
  
  hMuMu_rapidity_unfodata->SetMaximum(60);//100
  hMuMu_rapidity_unfodata->SetMinimum(0);
  make_hist(hMuMu_rapidity_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  auto hMuMu_rapidity_unfodata_syst = get_stat_uncertainty_hist(hMuMu_rapidity_unfodata);
  hMuMu_rapidity_unfodata_syst->Draw("e2");
  
  hMuMu_rapidity_unfodata->Draw("psameex0");
  hMuMu_rapidity_gen->Draw("histsamex0");
  hMuMu_rapidity_gen_syst->Draw("e2same");
  hMuMu_rapidity_recomc->Draw("ehistsamex0");
  make_hist(hMuMu_rapidity_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  hMuMu_rapidity_recodata->Draw("ehistsamex0");
  
  //CMS_lumi( c4, 1, 10 );
  leg2->Draw();
  c4->cd();
  TPad *pad8 = new TPad("pad8", "pad8", 0, 0.06, 1, 0.28);
  pad8->SetTopMargin(0);
  pad8->SetBottomMargin(0.3);
  pad8->SetGridx();
  pad8->Draw();
  pad8->cd();
  TString strRatio4=hMuMu_rapidity_gen->GetName();
  strRatio4+="_over_";
  strRatio4+=hMuMu_rapidity_gen->GetName();
  TH1D* hratio4=(TH1D*)hMuMu_rapidity_unfodata->Clone(strRatio4);
  hratio4->Sumw2();

  hratio4->Divide(hMuMu_rapidity_gen);
  hratio4->SetTitle("");
  hratio4->GetYaxis()->SetTitle("UnfoData/Gen");
  hratio4->SetLabelSize(0.1, "XYZ");
  hratio4->SetLabelFont(42, "XYZ");
  hratio4->SetLabelOffset(0.007, "XYZ");
  hratio4->SetTitleSize(0.11, "XYZ");
  hratio4->GetXaxis()->SetTitleOffset(1.0);
  hratio4->GetYaxis()->SetTitleOffset(0.3);
  hratio4->GetYaxis()->SetRangeUser(0.0, 2.0);
  hratio4->Draw("p");
  c4->Print("7thNov/data_rapidity.pdf");
  c4->Print("7thNov/data_rapidity.png");
  
  //invmass
  TCanvas* c5 = new TCanvas("c5","Dimuon invmass",50,50,W,H);
  prepare_canvas(c5);
  TPad *pad9 = new TPad("pad9", "pad9", 0, 0.3, 1, 1.0);
  pad9->SetBottomMargin(0);
  pad9->Draw();
  pad9->cd();
 // c5->cd();
  
  x_title = "m_{#gamma#gamma} (GeV)";
  y_title = "d#sigma_{#gamma#gamma}/dm_{#gamma#gamma} (nb)";
  
  cout<<"m_inv bin width:" <<hMuMu_invmass_unfodata->GetXaxis()->GetBinWidth(1)<<endl;
  
  make_hist(hMuMu_invmass_unfomc,x_title,y_title,20,1.0,data_unfold_color,1,1);
  make_hist(hMuMu_invmass_gen,x_title,y_title,20,1.0,mc_gen_color,1,2);
  
  hMuMu_invmass_gen->SetMaximum(25);//50
  hMuMu_invmass_gen->SetMinimum(0);
  hMuMu_invmass_gen->Draw("histsame");
  auto hMuMu_invmass_gen_syst = get_stat_uncertainty_hist(hMuMu_invmass_gen, syst_uncertainty_mc);
  hMuMu_invmass_gen_syst->Draw("e2same");
  
  make_hist(hMuMu_invmass_recomc,x_title,y_title,20,1.0,mc_reco_color,1,2);
  hMuMu_invmass_recomc->Draw("ehistsamex0");
  hMuMu_invmass_unfomc->Draw("psamex0");
  
  //CMS_lumi( c5, 1, 10 );
  leg1->Draw();
  c5->cd();
  TPad *pad10 = new TPad("pad10", "pad10", 0, 0.06, 1, 0.28);
  pad10->SetTopMargin(0);
  pad10->SetBottomMargin(0.3);
  pad10->SetGridx();
  pad10->Draw();
  pad10->cd();
  TString strRatio5=hMuMu_invmass_gen->GetName();
  strRatio5+="_over_";
  strRatio5+=hMuMu_invmass_gen->GetName();
  TH1D* hratio5=(TH1D*)hMuMu_invmass_unfomc->Clone(strRatio5);
  hratio5->Sumw2();

  hratio5->Divide(hMuMu_invmass_gen);
  hratio5->SetTitle("");
  hratio5->GetYaxis()->SetTitle("Unfo/Gen");
  hratio5->SetLabelSize(0.1, "XYZ");
  hratio5->SetLabelFont(42, "XYZ");
  hratio5->SetLabelOffset(0.007, "XYZ");
  hratio5->SetTitleSize(0.11, "XYZ");
  hratio5->GetXaxis()->SetTitleOffset(1.0);
  hratio5->GetYaxis()->SetTitleOffset(0.3);
  hratio5->Draw("p");
  c5->Print("7thNov/closure_invmass.pdf");
  c5->Print("7thNov/closure_invmass.png");
  
  TCanvas* c6 = new TCanvas("c6","Dimuon rap",50,50,W,H);
  prepare_canvas(c6);
  TPad *pad11 = new TPad("pad11", "pad11", 0, 0.3, 1, 1.0);
  pad11->SetBottomMargin(0);
  pad11->Draw();
  pad11->cd();
  //c6->cd();

  
  hMuMu_invmass_unfodata->SetMaximum(20);//50
  hMuMu_invmass_unfodata->SetMinimum(0);
  make_hist(hMuMu_invmass_unfodata,x_title,y_title,20,1.0,data_unfold_color,1,2);
  auto hMuMu_invmass_unfodata_syst = get_stat_uncertainty_hist(hMuMu_invmass_unfodata);
  hMuMu_invmass_unfodata_syst->Draw("e2");
  
  hMuMu_invmass_unfodata->Draw("psameex0");
  hMuMu_invmass_gen->Draw("histsamex0");
  hMuMu_invmass_gen_syst->Draw("e2same");
  hMuMu_invmass_recomc->Draw("ehistsamex0");
  make_hist(hMuMu_invmass_recodata,x_title,y_title,22,1.0,data_reco_color,1,2);
  hMuMu_invmass_recodata->Draw("ehistsamex0");
  
  //CMS_lumi( c6, 1, 10 );
  leg2->Draw();
  
  c6->cd();
  TPad *pad12 = new TPad("pad12", "pad12", 0, 0.06, 1, 0.28);
  pad12->SetTopMargin(0);
  pad12->SetBottomMargin(0.3);
  pad12->SetGridx();
  pad12->Draw();
  pad12->cd();
  TString strRatio6=hMuMu_invmass_gen->GetName();
  strRatio6+="_over_";
  strRatio6+=hMuMu_invmass_gen->GetName();
  TH1D* hratio6=(TH1D*)hMuMu_invmass_unfodata->Clone(strRatio6);
  hratio6->Sumw2();

  hratio6->Divide(hMuMu_invmass_gen);
  hratio6->SetTitle("");
  hratio6->GetYaxis()->SetTitle("UnfoData/Gen");
  hratio6->SetLabelSize(0.1, "XYZ");
  hratio6->SetLabelFont(42, "XYZ");
  hratio6->SetLabelOffset(0.007, "XYZ");
  hratio6->SetTitleSize(0.11, "XYZ");
  hratio6->GetXaxis()->SetTitleOffset(1.0);
  hratio6->GetYaxis()->SetTitleOffset(0.3);
  hratio6->GetYaxis()->SetRangeUser(0.0, 1.6);
  hratio6->Draw("p");
//  leg2->Draw();
  c6->Print("7thNov/data_invmass.pdf");
  c6->Print("7thNov/data_invmass.png");
  
  printBinCont(hMuMu_pt_unfodata);
  printBinCont(hMuMu_rapidity_unfodata);
  printBinCont(hMuMu_invmass_unfodata);
  
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
  
  output_hist->SetFillColorAlpha(input_hist->GetLineColor(), 0.2);
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
}
