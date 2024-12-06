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
#include "CMS_lumi.C"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TCut.h"
#include "TChain.h"
#include "THStack.h"

int xlo = 1;
int xhi = 4;
int nbin = 8;




const char *dir = "figures";

#define PI 3.141592653589

void make_hist(TH1D *&, Color_t , int );

const double luminosity       = 1635.123139823; // μb^-1
const double nEventsGeneratedSL = 63957000; // older number with less files 63398400; //starlight
const double xsecGeneratedSL    = 7920; // μb, starlight

const double nEventsGeneratedSC = 67536400; // older number with less files 67262800 superchic 
const double xsecGeneratedSC    = 8827.220; // μb
const double purity = 0.96; 
double scaleFactorsSC = 0.84 *  // NEE
                        0.96 *      // CHE
                        pow(0.98, 2)* // electron reco+ID
                        0.83 *      // HF veto
                        0.99;       // L1 EG trigger

double scaleFactorsSL = 0.84 *  // NEE
                        0.96 *      // CHE
                        pow(0.99, 2)* // electron reco+ID
                        0.83 *      // HF veto
                        0.99;       // L1 EG trigger

double lumiNormSC = xsecGeneratedSC*luminosity*purity*scaleFactorsSC/nEventsGeneratedSC;
double lumiNormSL = xsecGeneratedSL*luminosity*purity*scaleFactorsSL/nEventsGeneratedSL;


const double xsecGeneratedLbL    = 2.59; // μb
const double nEventsGeneratedLbL = 466000; 
double norm_exgg = 0.83*0.99*xsecGeneratedLbL*luminosity/nEventsGeneratedLbL; 


const double xsecGeneratedCEP    = 0.0058; // μb
const double nEventsGeneratedCEP = 668000; 
double norm_cep = xsecGeneratedCEP*luminosity/nEventsGeneratedCEP; 


double minAco = 0.0;
double maxAco = 0.1;

void BinLogX(TH1* h);

void plot_data_mc_log(){
 
  //gROOT->LoadMacro("CMS_lumi.C");
  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
   gStyle->SetOptStat(0);



  TFile *outf= new TFile("histos.root","recreate");  


  TChain *data = new TChain("output_tree");
  data->Add("/eos/cms/store/group/phys_heavyions/osuranyi/lbyl_2018/diphoton/diphoton_data11.root");

  TChain *lbyl = new TChain("output_tree");
  //lbyl->Add("lbyl_ntuple_all.root");
  lbyl->Add("/eos/cms/store/group/phys_heavyions/osuranyi/lbyl_2018/diphoton/diphoton_lbl.root");

 
  TChain *qed = new TChain("output_tree");
  qed->Add("/eos/cms/store/group/phys_heavyions/osuranyi/lbyl_2018/diphoton/diphoton_qed.root");
lumiNormSC = lumiNormSL;

  TChain *cep = new TChain("output_tree");
  cep->Add("/eos/cms/store/group/phys_heavyions/osuranyi/lbyl_2018/diphoton/diphoton_cep.root");


  cout << "norm Exgg " << norm_exgg << endl;

 // const TCut excl = " ngenTrk==0 && nEle==0 ";


 // const TCut calocut = "leadingEmEnergy_EB< 0.55 && leadingEmEnergy_EE < 3.16 && leadingHadEnergy_HB < 2.0 && leadingHadEnergy_HE < 3.0 && leadingHadEnergy_HF_Plus < 4.85 && leadingHadEnergy_HF_Minus < 4.12 ";

   const TCut no_addittional_tower  = " ok_neuexcl==1 ";
   const TCut no_trk  = "ok_chexcl_goodtracks == 1"; //&& ok_chexcl_muons == 1"; //" ok_chexcl==1";
   const TCut apply_cut =  "vSum_M > 5 && vSum_Pt<2 && pho_acop<0.01" && no_trk && no_addittional_tower;
   //const TCut apply_cut =  "vSum_M > 5 && vSum_Pt<1" && no_trk && no_addittional_tower;
  const TCut zdc_cut = "ok_zdcexcl_4n_pos == 1 && ok_zdcexcl_4n_neg == 1";

  const TCut aco_cut   = "vSum_M > 5 && vSum_Pt<2 && phoSwissCross_1 < 0.95 && phoSwissCross_2 < 0.95" && no_trk && no_addittional_tower && zdc_cut;
  const TCut aco_cut_mc   = "vSum_M > 5 && vSum_Pt<2 && phoSwissCross_1 < 0.95 && phoSwissCross_2 < 0.95" && no_trk && no_addittional_tower;

  //// data histogram

  TH1D* hpho_pt_data   = new TH1D("hpho_pt_data","",8,2,10);
  TH1D* hpho_eta_data  = new TH1D("hpho_eta_data","",12,-2.4,2.4);
  TH1D* hpho_phi_data  = new TH1D("hpho_phi_data","",12,-PI,PI);

  TH1D* hpho_pt2_data   = new TH1D("hpho_pt2_data","",8,2,10);
  TH1D* hpho_eta2_data  = new TH1D("hpho_eta2_data","",12,-2.4,2.4);
  TH1D* hpho_phi2_data  = new TH1D("hpho_phi2_data","",12,-PI,PI);
    
  TH1D* hdiphoton_pt_data         = new TH1D("hdiphoton_pt_data","",10,0,2);
  TH1D* hdiphoton_rapidity_data   = new TH1D("hdiphoton_rapidity_data","",12,-2.4,2.4);
  TH1D* hinvmass_data             = new TH1D("hinvmass_data","",30,0,30);
  TH1D* hacoplanarity_data   = new TH1D("hacoplanarity_data","",20,minAco,maxAco);

  //BinLogX(hacoplanarity_data);
  
  TH1D* hntrk = new TH1D("hntrk","",8,2,10);
  cout << " data project ongoing" << endl;
  data->Project(hpho_pt_data->GetName(), "phoEt_1",apply_cut && "ok_zdcexcl==1" );
  data->Project(hpho_eta_data->GetName(), "phoEta_1",apply_cut  && "ok_zdcexcl==1" );
  data->Project(hpho_phi_data->GetName(), "phoPhi_1",apply_cut  && "ok_zdcexcl==1" );
  data->Project(hpho_pt2_data->GetName(), "phoEt_2",apply_cut && "ok_zdcexcl==1" );
  data->Project(hpho_eta2_data->GetName(), "phoEta_2",apply_cut  && "ok_zdcexcl==1" );
  data->Project(hpho_phi2_data->GetName(), "phoPhi_2",apply_cut  && "ok_zdcexcl==1" );

  hpho_pt_data->Add(hpho_pt2_data);
  hpho_eta_data->Add(hpho_eta2_data);
  hpho_phi_data->Add(hpho_phi2_data);
  
  data->Project(hdiphoton_pt_data->GetName(), "vSum_Pt",apply_cut  && "ok_zdcexcl==1" );
  data->Project(hdiphoton_rapidity_data->GetName(), "vSum_Rapidity",apply_cut  && "ok_zdcexcl==1" );
  data->Project(hinvmass_data->GetName(), "vSum_M",apply_cut  /*&& "ok_zdcexcl==1"*/ );
  data->Project(hacoplanarity_data->GetName(),  "pho_acop", aco_cut /*&& "ok_zdcexcl==1"*/ );
  cout << " data project finished" << endl;
 
  // get exclusive gg file


  TH1D* hpho_pt_lbyl   = new TH1D("hpho_pt_lbyl","",8,2,10);
  TH1D* hpho_eta_lbyl  = new TH1D("hpho_eta_lbyl","",12,-2.4,2.4);
  TH1D* hpho_phi_lbyl  = new TH1D("hpho_phi_lbyl","",12,-PI,PI);

  TH1D* hpho_pt2_lbyl   = new TH1D("hpho_pt2_lbyl","",8,2,10);
  TH1D* hpho_eta2_lbyl  = new TH1D("hpho_eta2_lbyl","",12,-2.4,2.4);
  TH1D* hpho_phi2_lbyl  = new TH1D("hpho_phi2_lbyl","",12,-PI,PI);
  
 
  TH1D* hdiphoton_pt_lbyl         = new TH1D("hdiphoton_pt_lbyl","",10,0,2);
  TH1D* hdiphoton_rapidity_lbyl   = new TH1D("hdiphoton_rapidity_lbyl","",12,-2.4,2.4);
  TH1D* hinvmass_lbyl    = new TH1D("hinvmass_lbyl","",30,0,30);
  TH1D* hacoplanarity_lbyl   = new TH1D("hacoplanarity_lbyl","",20,minAco,maxAco);

  //BinLogX(hacoplanarity_lbyl);


  cout << " signal project ongoing" << endl;
  lbyl->Project(hpho_pt_lbyl->GetName(), "phoEt_1",apply_cut );
  lbyl->Project(hpho_eta_lbyl->GetName(), "phoEta_1",apply_cut );
  lbyl->Project(hpho_phi_lbyl->GetName(), "phoPhi_1",apply_cut );
  lbyl->Project(hpho_pt2_lbyl->GetName(), "phoEt_2",apply_cut );
  lbyl->Project(hpho_eta2_lbyl->GetName(), "phoEta_2",apply_cut );
  lbyl->Project(hpho_phi2_lbyl->GetName(), "phoPhi_2",apply_cut );

  
  hpho_pt_lbyl->Add(hpho_pt2_lbyl);
  hpho_eta_lbyl->Add(hpho_eta2_lbyl);
  hpho_phi_lbyl->Add(hpho_phi2_lbyl);
  

  lbyl->Project(hdiphoton_pt_lbyl->GetName(), "vSum_Pt",apply_cut );
  lbyl->Project(hdiphoton_rapidity_lbyl->GetName(), "vSum_Rapidity",apply_cut );
  lbyl->Project(hinvmass_lbyl->GetName(), "vSum_M",apply_cut );
  lbyl->Project(hacoplanarity_lbyl->GetName(),  "pho_acop", aco_cut_mc );

  
  hpho_pt_lbyl->Scale(norm_exgg);   hpho_eta_lbyl->Scale(norm_exgg);   hpho_phi_lbyl->Scale(norm_exgg);       
  hacoplanarity_lbyl->Scale(norm_exgg);
  hdiphoton_pt_lbyl->Scale(norm_exgg);   hdiphoton_rapidity_lbyl->Scale(norm_exgg);     hinvmass_lbyl->Scale(norm_exgg);      
  cout << " signal project finished" << endl;



  // get qed file

  TH1D* hpho_pt_qed   = new TH1D("hpho_pt_qed","",8,2,10);
  TH1D* hpho_eta_qed  = new TH1D("hpho_eta_qed","",12,-2.4,2.4);
  TH1D* hpho_phi_qed  = new TH1D("hpho_phi_qed","",12,-PI,PI);

  TH1D* hpho_pt2_qed   = new TH1D("hpho_pt2_qed","",8,2,10);
  TH1D* hpho_eta2_qed  = new TH1D("hpho_eta2_qed","",12,-2.4,2.4);
  TH1D* hpho_phi2_qed  = new TH1D("hpho_phi2_qed","",12,-PI,PI);
  
  TH1D* hdiphoton_pt_qed         = new TH1D("hdiphoton_pt_qed","",10,0,2);
  TH1D* hdiphoton_rapidity_qed   = new TH1D("hdiphoton_rapidity_qed","",12,-2.4,2.4);
  TH1D* hinvmass_qed    = new TH1D("hinvmass_qed","",30,0,30);
  TH1D* hacoplanarity_qed   = new TH1D("hacoplanarity_qed","",20,minAco,maxAco);

  //BinLogX(hacoplanarity_qed);


  cout << "qed project ongoing" << endl;

  qed->Project(hpho_pt_qed->GetName(), "phoEt_1",apply_cut );
  qed->Project(hpho_eta_qed->GetName(), "phoEta_1",apply_cut );
  qed->Project(hpho_phi_qed->GetName(), "phoPhi_1",apply_cut );
  qed->Project(hpho_pt2_qed->GetName(), "phoEt_2",apply_cut );
  qed->Project(hpho_eta2_qed->GetName(), "phoEta_2",apply_cut );
  qed->Project(hpho_phi2_qed->GetName(), "phoPhi_2",apply_cut );

  hpho_pt_qed->Add(hpho_pt2_qed);
  hpho_eta_qed->Add(hpho_eta2_qed);
  hpho_phi_qed->Add(hpho_phi2_qed);


  qed->Project(hdiphoton_pt_qed->GetName(), "vSum_Pt",apply_cut );
  qed->Project(hdiphoton_rapidity_qed->GetName(), "vSum_Rapidity",apply_cut );
  qed->Project(hinvmass_qed->GetName(), "vSum_M",apply_cut );
  qed->Project(hacoplanarity_qed->GetName(), "pho_acop", aco_cut_mc);
  
  double old_norm = lumiNormSC;
  lumiNormSC = old_norm*2.7; //(hacoplanarity_data->Integral(5,20) - hacoplanarity_lbyl->Integral(5,20)) / hacoplanarity_qed->Integral(5,20);
  std::cout << lumiNormSC/old_norm << std::endl;

  hpho_pt_qed->Scale(lumiNormSC);   hpho_eta_qed->Scale(lumiNormSC);   hpho_phi_qed->Scale(lumiNormSC);     
  hacoplanarity_qed ->Scale(lumiNormSC);   
  hdiphoton_pt_qed->Scale(lumiNormSC);   hdiphoton_rapidity_qed->Scale(lumiNormSC);     hinvmass_qed->Scale(lumiNormSC);    
  cout << "qed project finished" << endl;


 // get cep file


  TH1D* hpho_pt_cep   = new TH1D("hpho_pt_cep","",8,2,10);
  TH1D* hpho_eta_cep  = new TH1D("hpho_eta_cep","",12,-2.4,2.4);
  TH1D* hpho_phi_cep  = new TH1D("hpho_phi_cep","",12,-PI,PI);

  TH1D* hpho_pt2_cep   = new TH1D("hpho_pt2_cep","",8,2,10);
  TH1D* hpho_eta2_cep  = new TH1D("hpho_eta2_cep","",12,-2.4,2.4);
  TH1D* hpho_phi2_cep  = new TH1D("hpho_phi2_cep","",12,-PI,PI);
  
  
  TH1D* hdiphoton_pt_cep         = new TH1D("hdiphoton_pt_cep","",10,0,2);
  TH1D* hdiphoton_rapidity_cep   = new TH1D("hdiphoton_rapidity_cep","",12,-2.4,2.4);
  TH1D* hinvmass_cep    = new TH1D("hinvmass_cep","",30,0,30);
  TH1D* hacoplanarity_cep   = new TH1D("hacoplanarity_cep","",20,minAco,maxAco);

  //BinLogX(hacoplanarity_cep);


  cout << "cep project ongoing" << endl;

  cep->Project(hpho_pt_cep->GetName(), "phoEt_1",apply_cut );
  cep->Project(hpho_eta_cep->GetName(), "phoEta_1",apply_cut );
  cep->Project(hpho_phi_cep->GetName(), "phoPhi_1",apply_cut );
  cep->Project(hpho_pt2_cep->GetName(), "phoEt_2",apply_cut );
  cep->Project(hpho_eta2_cep->GetName(), "phoEta_2",apply_cut );
  cep->Project(hpho_phi2_cep->GetName(), "phoPhi_2",apply_cut );

  hpho_pt_cep->Add(hpho_pt2_cep);
  hpho_eta_cep->Add(hpho_eta2_cep);
  hpho_phi_cep->Add(hpho_phi2_cep);
  

  cep->Project(hdiphoton_pt_cep->GetName(), "vSum_Pt",apply_cut );
  cep->Project(hdiphoton_rapidity_cep->GetName(), "vSum_Rapidity",apply_cut );
  cep->Project(hinvmass_cep->GetName(), "vSum_M",apply_cut );
  cep->Project(hacoplanarity_cep->GetName(),  "pho_acop", aco_cut_mc );
 

  // Normalizing CEP xsec to tail
  norm_cep = 0;(hacoplanarity_data->Integral(1,5) - hacoplanarity_lbyl->Integral(1,5) - hacoplanarity_qed->Integral(1,5)) / hacoplanarity_cep->Integral(1,5);
  //std::cout << hacoplanarity_data->Integral(1,1) << "\t" << hacoplanarity_lbyl->Integral(1,1) << "\t" << hacoplanarity_qed->Integral(1,1) << std::endl;

  hpho_pt_cep->Scale(norm_cep);   hpho_eta_cep->Scale(norm_cep);   hpho_phi_cep->Scale(norm_cep);     
  hacoplanarity_cep ->Scale(norm_cep);  
  hdiphoton_pt_cep->Scale(norm_cep);   hdiphoton_rapidity_cep->Scale(norm_cep);     hinvmass_cep->Scale(norm_cep);      
  cout << "cep project finished" << endl;
  /*hpho_pt_cep->Scale(0.23);   hpho_eta_cep->Scale(0.23);   hpho_phi_cep->Scale(0.23);     
  hdpt_cep->Scale(0.23);      hdeta_cep->Scale(0.23);      hdphi_cep ->Scale(0.23);   hacoplanarity_cep ->Scale(0.23);  
  hdiphoton_pt_cep->Scale(0.23);   hdiphoton_rapidity_cep->Scale(0.23);     hinvmass_cep->Scale(0.23);   */  
  
  
    //int ibin = 20;
    hinvmass_data->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    hdiphoton_pt_data->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    hdiphoton_rapidity_data->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    hpho_pt_data->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    hpho_eta_data->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    hpho_phi_data->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    hacoplanarity_data->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    //hinvmass_data->SetBinError(1,0);
    //hinvmass_data->SetBinError(2,0);

    double err_low = hinvmass_data->GetBinErrorLow(3);
    double err_up = hinvmass_data->GetBinErrorUp(3);

   cout << " low error: " << err_low << " up error:" << err_up << endl;
  
  TLegend *leg2=new TLegend(0.18,0.65,0.60,0.91);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(24);
  leg2->AddEntry(hinvmass_data,"Data ","pl");
  leg2->AddEntry(hinvmass_lbyl,"LbyL #gamma #gamma #rightarrow #gamma #gamma (MC)","f");
  leg2->AddEntry(hinvmass_cep,"CEP (gg #rightarrow #gamma #gamma) + other bkg","f");
  leg2->AddEntry(hinvmass_qed,"QED #gamma #gamma #rightarrow e^{+}e^{-} (MC)","f");

  
  TH1D* hinvmass_empty   = new TH1D("hinvmass_empty","",1,0,30);
  TH1D* hpho_pt_empty   = new TH1D("hpho_pt_empty","",1,0,30);
    
  int W = 700;
  int H = 600;
  
  float T = 0.08;
  float B = 0.14; 
  float L = 0.14;
  float R = 0.04;


  //// Diphoton invmass......................
  THStack *hs = new THStack("hs","Diphoton invariant mass");
  TCanvas* c1 = new TCanvas("c1","Dimuon pT",50,50,W,H);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin( L );
  c1->SetRightMargin( R );
  c1->SetTopMargin( T );
  c1->SetBottomMargin( B );
  c1->SetTickx(0);
  c1->SetTicky(0);
  c1->cd();

  hinvmass_empty->SetMaximum(60);
  hinvmass_empty->SetMinimum(0);
 
  make_hist(hinvmass_empty, kBlack, 21);
  hinvmass_empty->SetLineColor(kBlack);

  hinvmass_empty->GetXaxis()->SetTitle("Diphoton Invariant Mass (GeV)");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  hinvmass_empty->GetYaxis()->SetTitle("# events");
  hinvmass_empty->Draw("p");


  make_hist(hinvmass_qed, kYellow, 21);
  //hinvmass_qed->Draw("psames");
  hs->Add(hinvmass_qed);

  make_hist(hinvmass_cep, 51, 21);
  hinvmass_cep->SetFillStyle(3001);
  //hinvmass_cep->Draw("psames");
  hs->Add(hinvmass_cep);

  make_hist(hinvmass_lbyl, 95, 21);
  //hinvmass_lbyl->Draw("p");
  hs->Add(hinvmass_lbyl);
  hs->Draw("histsame");   hs->GetXaxis()->SetTitle("Diphoton Invariant Mass (GeV)");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  //hs->GetXaxis()->SetRangeUser(5,20);
  hs->GetYaxis()->SetTitle("# events");

  make_hist(hinvmass_data, kBlack, 21);
  hinvmass_data->SetLineColor(kBlack);
  hinvmass_data->Draw("psames");
  
  leg2->Draw();
  CMS_lumi( c1, 104, 33,lumi_PbPb2015 );
  c1->Print(Form("./%s/diphoton_invmass.gif",dir));
  c1->Print(Form("./%s/diphoton_invmass.pdf",dir));
  c1->Print(Form("./%s/diphoton_invmass.C",dir));

  cout << " Data invamss integral:" << hinvmass_data->Integral() << endl;
  cout << " LbyL MC invamss integral:" << hinvmass_lbyl->Integral() << endl;

  //// Diphoton pT....................
  THStack *hs_pt = new THStack("hs_pt","Diphoton pt");
  TCanvas* c2 = new TCanvas("c2","Dimuon pT",50,50,W,H);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillStyle(0);
  c2->SetFrameBorderMode(0);
  c2->SetLeftMargin( L );
  c2->SetRightMargin( R );
  c2->SetTopMargin( T );
  c2->SetBottomMargin( B );
  c2->SetTickx(0);
  c2->SetTicky(0);
  c2->cd();
  //gPad->SetLogy();
  hs_pt->SetMaximum(50);
  hs_pt->SetMinimum(0);
 
  make_hist(hdiphoton_pt_qed, kYellow, 21);
  hs_pt->Add(hdiphoton_pt_qed);

  make_hist(hdiphoton_pt_cep, 51, 21);
  hdiphoton_pt_cep->SetFillStyle(3001);
  hs_pt->Add(hdiphoton_pt_cep);

  make_hist(hdiphoton_pt_lbyl, 95, 21);
  //hdiphoton_pt_lbyl->Draw("p");
  hs_pt->Add(hdiphoton_pt_lbyl);
 
  hs_pt->Draw("hist");   hs_pt->GetXaxis()->SetTitle("Diphoton p_{T} (GeV)");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  hs_pt->GetXaxis()->SetRangeUser(0,1);
  hs_pt->GetYaxis()->SetTitle("# events");

  make_hist(hdiphoton_pt_data, kBlack, 21);
  hdiphoton_pt_data->SetLineColor(kBlack);
  hdiphoton_pt_data->Draw("psame");
  
  leg2->Draw();
  CMS_lumi( c2, 104, 33,lumi_PbPb2015 );
  c2->Print(Form("./%s/diphoton_pt.gif",dir));
  c2->Print(Form("./%s/diphoton_pt.pdf",dir));
  c2->Print(Form("./%s/diphoton_pt.C",dir));


  //// Diphoton rapidity....................
  THStack *hs_rapidity = new THStack("hs_rapidity","Diphoton rapidity");
  TCanvas* c3 = new TCanvas("c3","Dimuon rapidity",50,50,W,H);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetFrameFillStyle(0);
  c3->SetFrameBorderMode(0);
  c3->SetLeftMargin( L );
  c3->SetRightMargin( R );
  c3->SetTopMargin( T );
  c3->SetBottomMargin( B );
  c3->SetTickx(0);
  c3->SetTicky(0);
  c3->cd();
  //gPad->SetLogy();
  hs_rapidity->SetMaximum(30);
  hs_rapidity->SetMinimum(0);
  
  make_hist(hdiphoton_rapidity_qed, kYellow, 21);
  hs_rapidity->Add(hdiphoton_rapidity_qed);

  make_hist(hdiphoton_rapidity_cep, 51, 21);
  hdiphoton_rapidity_cep->SetFillStyle(3001);
  hs_rapidity->Add(hdiphoton_rapidity_cep);

  make_hist(hdiphoton_rapidity_lbyl, 95, 21);
  //hdiphoton_rapidity_lbyl->Draw("p");
  hs_rapidity->Add(hdiphoton_rapidity_lbyl);

  hs_rapidity->Draw("hist");   hs_rapidity->GetXaxis()->SetTitle("Diphoton y");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  //hs_rapidity->GetXaxis()->SetRangeUser(0,2);
  hs_rapidity->GetYaxis()->SetTitle("Events / (0.8)");
 
  make_hist(hdiphoton_rapidity_data, kBlack, 21);
  hdiphoton_rapidity_data->SetLineColor(kBlack);
  hdiphoton_rapidity_data->Draw("psame");
  
  leg2->Draw();
  CMS_lumi( c3, 104, 33,lumi_PbPb2015 );
  c3->Print(Form("./%s/diphoton_rap.gif",dir));
  c3->Print(Form("./%s/diphoton_rap.pdf",dir));
  c3->Print(Form("./%s/diphoton_rap.C",dir));


  /// diphoton acoplanarity ...........................................

  TLegend *leg1=new TLegend(0.23,0.65,0.95,0.91);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(43);
  leg1->SetTextSize(24);
  leg1->AddEntry(hinvmass_data,"Data ","pl");
  leg1->AddEntry(hinvmass_lbyl,"LbyL #gamma #gamma #rightarrow #gamma #gamma (MC)","f");
  leg1->AddEntry(hinvmass_cep,"CEP (gg #rightarrow #gamma #gamma) + other bkg","f");
  leg1->AddEntry(hinvmass_qed,"QED #gamma #gamma #rightarrow e^{+}e^{-} (MC)","f");


  THStack *hs_aco = new THStack("hs_aco","");
  TCanvas* c4 = new TCanvas("c4","Dimuon rapidity",50,50,W,H);
  c4->SetFillColor(0);
  c4->SetBorderMode(0);
  c4->SetFrameFillStyle(0);
  c4->SetFrameBorderMode(0);
  c4->SetLeftMargin( L );
  c4->SetRightMargin( R );
  c4->SetTopMargin( T );
  c4->SetBottomMargin( B );
  c4->SetTickx(0);
  c4->SetTicky(0);
  //c4->SetLimits(0.001,0.1);
  c4->SetLogx();
  c4->cd();
  //gPad->SetLogy();
  hs_aco->SetMaximum(50);
  hs_aco->SetMinimum(0);
  
  make_hist(hacoplanarity_qed, kYellow, 21);
  hs_aco->Add(hacoplanarity_qed);

  make_hist(hacoplanarity_cep, 51, 21);
  hacoplanarity_cep->SetFillStyle(3001);
  hs_aco->Add(hacoplanarity_cep);

  make_hist(hacoplanarity_lbyl, 95, 21);
  hacoplanarity_lbyl->Draw("p");
  hs_aco->Add(hacoplanarity_lbyl);
  
  hs_aco->Draw("hist"); hs_aco->GetXaxis()->SetTitle("Diphoton A_{#phi}");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  //hs_aco->GetXaxis()->SetRangeUser(0,0.2);   
  hs_aco->GetXaxis()->SetNdivisions(509);
  //hs_aco->GetXaxis()->SetRangeUser(0,2);
  hs_aco->GetYaxis()->SetTitle("Events / (0.005)");
  
  make_hist(hacoplanarity_data, kBlack, 21);
  hacoplanarity_data->SetLineColor(kBlack);
  hacoplanarity_data->Draw("p same"); //"p same"
  
  leg1->Draw();
  CMS_lumi( c4, 104, 33,lumi_PbPb2015 );
  c4->Print(Form("./%s/acoplanarity_ged.gif",dir));
  c4->Print(Form("./%s/acoplanarity_ged.pdf",dir));
  c4->Print(Form("./%s/acoplanarity_ged.C",dir));
 


  /// single photon pT ...............................................
  THStack *hs_spt  = new THStack("hs_spt","photon pt");
  TCanvas* c5 = new TCanvas("c5","Single pt",50,50,W,H);
  c5->SetFillColor(0);
  c5->SetBorderMode(0);
  c5->SetFrameFillStyle(0);
  c5->SetFrameBorderMode(0);
  c5->SetLeftMargin( L );
  c5->SetRightMargin( R );
  c5->SetTopMargin( T );
  c5->SetBottomMargin( B );
  c5->SetTickx(0);
  c5->SetTicky(0);
  c5->cd();


  hpho_pt_empty->SetMaximum(100);
  hpho_pt_empty->SetMinimum(0);
 
  make_hist(hpho_pt_empty, kBlack, 21);
  hpho_pt_empty->SetLineColor(kBlack);

  hpho_pt_empty->GetXaxis()->SetTitle("Photon E_{T} (GeV)");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  hpho_pt_empty->GetYaxis()->SetTitle("# events");
  hpho_pt_empty->Draw("p");


  //gPad->SetLogy();
  //hs_spt->SetMaximum(25);
  //hs_spt->SetMinimum(0);
  
  make_hist(hpho_pt_qed, kYellow, 21);
  hs_spt->Add(hpho_pt_qed);

  make_hist(hpho_pt_cep, 51, 21);
  hpho_pt_cep->SetFillStyle(3001);
  hs_spt->Add(hpho_pt_cep);

  make_hist(hpho_pt_lbyl, 95, 21);
  //hpho_pt_lbyl->Draw("p");
  hs_spt->Add(hpho_pt_lbyl);

  hs_spt->Draw("histsame");   hs_spt->GetXaxis()->SetTitle("Photon E_{T} (GeV)");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  //hs_spt->GetXaxis()->SetRangeUser(0,2);
  hs_spt->GetYaxis()->SetTitle("# events");
  make_hist(hpho_pt_data, kBlack, 21);
  hpho_pt_data->SetLineColor(kBlack);
  hpho_pt_data->Draw("psame");
  leg2->Draw();
  CMS_lumi( c5, 104, 33,lumi_PbPb2015 );
  c5->Print(Form("./%s/single_pt.gif",dir));
  c5->Print(Form("./%s/single_pt.pdf",dir));
  c5->Print(Form("./%s/single_pt.C",dir));

  /// single photon eta ...............................................
  THStack *hs_eta = new THStack("hs_eta","photon eta");
  TCanvas* c6 = new TCanvas("c6","Single pt",50,50,W,H);
  c6->SetFillColor(0);
  c6->SetBorderMode(0);
  c6->SetFrameFillStyle(0);
  c6->SetFrameBorderMode(0);
  c6->SetLeftMargin( L );
  c6->SetRightMargin( R );
  c6->SetTopMargin( T );
  c6->SetBottomMargin( B );
  c6->SetTickx(0);
  c6->SetTicky(0);
  c6->cd();
  //gPad->SetLogy();
  hs_eta->SetMaximum(50);
  hs_eta->SetMinimum(0);
  
  make_hist(hpho_eta_qed, kYellow, 21);
  hs_eta->Add(hpho_eta_qed);

  make_hist(hpho_eta_cep, 51, 21);
  hpho_eta_cep->SetFillStyle(3001);
  hs_eta->Add(hpho_eta_cep);

  make_hist(hpho_eta_lbyl, 95, 21);
  //hpho_eta_lbyl->Draw("p");
  hs_eta->Add(hpho_eta_lbyl);

  hs_eta->Draw("hist");   hs_eta->GetXaxis()->SetTitle("Photon #eta");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  //hs_eta->GetXaxis()->SetRangeUser(0,2);
  hs_eta->GetYaxis()->SetTitle("# events");
  make_hist(hpho_eta_data, kBlack, 21);
  hpho_eta_data->SetLineColor(kBlack);
  hpho_eta_data->Draw("psame");
  leg2->Draw();
  CMS_lumi( c6, 104, 33,lumi_PbPb2015 );
  c6->Print(Form("./%s/single_eta.gif",dir));
  c6->Print(Form("./%s/single_eta.pdf",dir));
  c6->Print(Form("./%s/single_eta.C",dir));

  /// single photon phi ...............................................
  THStack *hs_phi = new THStack("hs_phi","photon phi");
  TCanvas* c7 = new TCanvas("c7","Single phi",50,50,W,H);
  c7->SetFillColor(0);
  c7->SetBorderMode(0);
  c7->SetFrameFillStyle(0);
  c7->SetFrameBorderMode(0);
  c7->SetLeftMargin( L );
  c7->SetRightMargin( R );
  c7->SetTopMargin( T );
  c7->SetBottomMargin( B );
  c7->SetTickx(0);
  c7->SetTicky(0);
  c7->cd();
  //gPad->SetLogy();
  hs_phi->SetMaximum(50);
  hs_phi->SetMinimum(0);
  
  make_hist(hpho_phi_qed, kYellow, 21);
  //hpho_phi_qed->Draw("p");
  hs_phi->Add(hpho_phi_qed);

  make_hist(hpho_phi_cep, 51, 21);
  hpho_phi_cep->SetFillStyle(3001);
  //hpho_phi_qed->Draw("psames");
  hs_phi->Add(hpho_phi_cep);

  make_hist(hpho_phi_lbyl, 95, 21);
  //hpho_phi_lbyl->Draw("psames");
  hs_phi->Add(hpho_phi_lbyl);

  hs_phi->Draw("hist");   hs_phi->GetXaxis()->SetTitle("Photon #phi");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  //hs_phi->GetXaxis()->SetRangeUser(0,2);
  hs_phi->GetYaxis()->SetTitle("# events");
  make_hist(hpho_phi_data, kBlack, 21);
  hpho_phi_data->SetLineColor(kBlack);
  hpho_phi_data->Draw("psames");
  leg2->Draw();
  CMS_lumi( c7, 104, 33,lumi_PbPb2015 );
  c7->Print(Form("./%s/single_phi.gif",dir));
  c7->Print(Form("./%s/single_phi.pdf",dir));
  c7->Print(Form("./%s/single_phi.C",dir));

  /*double err, err1,err2,err3 ; 
  hinvmass_data->IntegralAndError(xlo,xhi,err);
  cout << "Data Integral:   "  << hinvmass_data->Integral(xlo,xhi) << "  error " << err << endl;
  hinvmass_lbyl->IntegralAndError(xlo,xhi,err1);
  cout << "Light by light Integral:   "  << hinvmass_lbyl->Integral() << "  error " << err1 << endl;
  *hinvmass_cep->IntegralAndError(xlo,xhi,err3);
  cout << "cep Integral:   "  << hinvmass_cep->Integral(xlo,xhi) << "  error " << err3 << endl;
  hinvmass_qed->IntegralAndError(xlo,xhi,err2);
  cout << "qed Integral:   "  << hinvmass_qed->Integral(xlo,xhi) << "  error " << err2 << endl;*/

 /* double err, err1,err2,err3 ; 
  hacoplanarity_data->IntegralAndError(xlo,xhi,err);
  cout << "Data Integral:   "  << hacoplanarity_data->Integral(xlo,xhi) << "  error " << err << endl;
  hacoplanarity_lbyl->IntegralAndError(xlo,xhi,err1);
  cout << "Light by light Integral:   "  << hacoplanarity_lbyl->Integral(xlo,xhi) << "  error " << err1 << endl;
  hacoplanarity_cep->IntegralAndError(xlo,xhi,err3);
  cout << "cep Integral:   "  << hacoplanarity_cep->Integral(xlo,xhi) << "  error " << err3 << endl;
  hacoplanarity_qed->IntegralAndError(xlo,xhi,err2);
  cout << "qed Integral:   "  << hacoplanarity_qed->Integral(xlo,xhi) << "  error " << err2 << endl;

  outf->cd();
  hinvmass_data->Write();
  hinvmass_lbyl->Write();
  hinvmass_qed->Write();
  hinvmass_cep->Write();

  hdiphoton_pt_data->Write();
  hdiphoton_pt_lbyl->Write();
  hdiphoton_pt_qed->Write();
  hdiphoton_pt_cep->Write();

  hdiphoton_rapidity_data->Write();
  hdiphoton_rapidity_lbyl->Write();
  hdiphoton_rapidity_qed->Write();
  hdiphoton_rapidity_cep->Write();

  hacoplanarity_data->Write();
  hacoplanarity_lbyl->Write();
  hacoplanarity_qed->Write();
  hacoplanarity_cep->Write();*/
 

 
}

void make_hist(TH1D *& hist, Color_t kcolor, int kstyle){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.05, "XYZ");
 
  hist->SetMarkerColor(kcolor);
  hist->SetMarkerSize(1.5);
  hist->SetFillColor(kcolor);
  hist->SetMarkerStyle(kstyle);
  
  //hist->SetTitle("");
}

void BinLogX(TH1* h){
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++){
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}
