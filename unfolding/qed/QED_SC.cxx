//this is a script to make response matrix from LbyL MC 
// Background subtracted from data
// Data unfolded
// All numbers (cross-section, luminosity for data to get cross-section is in nano barn
// Created by Ruchi Chudasama
///////////////////////////////////////////////////////////////////////////////////////
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include "TVectorD.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "THStack.h"
#include "TCut.h"
#include "TAxis.h"
#include "TChain.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldInvert.h"
#endif

#include "ReadqedTree.C"

const int nSample = 1;
const char *Sample[nSample]={"QED_SC"};

void make_canvas(TCanvas *&);
void make_canvas_ratio(TCanvas *&);
void make_hist(TH1D *&, Color_t , int );
void make_hist_ratio(TH1D *&,Color_t, int) ;
Double_t RecoID_SF_Electron( double ETT);

//const double luminosity       = 1635.123139823; // μb^-1
const double luminosity       = 1;

double scaleFactorSC_Old = 0.8477 *  // NEE    21.12.2021
                        0.9322 *      // CHE  21.12.2021
                        pow(0.952, 2)* // electron reco+ID 21.12.2021
                        0.8643 *      // HF veto
                       1.0006;       // L1 EG trigger

const double xsecGeneratedSC    = 8827220; // nb Superchic
const double nEventsGeneratedSC = 59260000;  //Superchic


double norm_LbLSC_other = (1.0089 * 0.8716 * 0.9252 * 0.8487 * xsecGeneratedSC*luminosity)/nEventsGeneratedSC;  //31.12.21 (SF = 1.048, took the square here). 
//double norm_LbLSC_other = (1.0089 * 0.8716 * 0.943*0.943* 0.8487 * xsecGeneratedSC*luminosity)/nEventsGeneratedSC
//double norm_LbLSC = scaleFactorSC*xsecGeneratedSC*luminosity/nEventsGeneratedSC;
//double norm_LbLSC = 1;
double W1 = xsecGeneratedSC*luminosity/nEventsGeneratedSC; 
//double W1 = 1;
//double norm_LbLSC = scaleFactorSC*7920000*luminosity/66750000; //For Starlight
double norm_LbLSL_other = (1.0089 * 0.8716 * 0.9252 * 0.8487 * 7920000*luminosity)/66750000;

//double W1 = 7920000*luminosity/66750000; 

void printOutput(TH1D *hist);
double_t GetDeltaPhi(double dphi, double phi);

float getAcoBinSF(float acop);


void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kMagenta);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}


const int nRapbins=3;
double Rapbin[nRapbins]={0.0, 1.3, 2.2};
const int nRapbin= sizeof(Rapbin)/sizeof(double) - 1;

// const int nMassbins=25;
// double Massbin[nMassbins]={4,6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,62,68,78,88,100};
// const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;
//  const int nMassbins=23;
//  double Massbin[nMassbins]={4,6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,48,54,60,68,78,100};
//  const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;
 
 const int nMassbins=22;
 double Massbin[nMassbins]={5,7,9,11,13,15,17,19,21,23,25,27,31,35,39,43,49,55,61,69,79,100};
 const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;
TH1D* subtractBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep, const char* name, int nbins, float xmin, float xmax);
TH1D* subtractInvBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep, const char* name);
TH1D* getXSecHist(TH1D *hist, const char* name, int nbins, float xmin, float xmax, float lumi);
TH1D* getInvXSecHist(TH1D *hist, const char* name, float lumi);
TH1D* getRapXSecHist(TH1D *hist, const char* name, float lumi);

TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, TH1D*, double , double , double , double , double , double , const char *, bool); 
void QED_SC(int method=1){
 
 // double wt[nSample] = {norm_LbLSC};
  int itr = 1; 
  //gROOT->LoadMacro("CMS_lumi.C");
  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  
  char MethID1[100]; 
  if(method==1){
        
 
    } 
  if(method==2){

    }
  if(method==3){
    

    }

  string mc_name = "SC";


  bool absRap = false;

  TFile *outf= new TFile(Form("unfolding_histograms_Bayes_QEDSC%s_300.root",mc_name.c_str()),"recreate");  
  
  TChain *qed[nSample];
  
  TH1D* hdiele_Pt[nSample], *hdiele_Rapidity[nSample], *hdiele_Invmass[nSample], *hAcoplanarity[nSample], *hdiele_Costhetastar[nSample] ; TH1F *hSF[nSample]; 
  
  TH1D* hdiele_GenPt[nSample], *hdiele_GenRapidity[nSample], *hdiele_GenInvmass[nSample], *hdiele_GenCosthetastar[nSample];
  TH2D* hdiele_RecoGenPt[nSample], *hdiele_RecoGenRapidity[nSample], *hRecoGenInvmass[nSample], *hdiele_RecoGenCosthetastar[nSample];


  
  cout << "define hists"<< endl;
  
  int ptbin = 7;
  int rapbin = 11;

  int i =0;
  //for (int i = 0; i < nSample; i++){
   // TFile *f = new TFile(Form("gammaUPC_over_SC_pTee_scalefactor.root"),"r");
    TFile *f = new TFile(Form("ratio_pTee_gammaUPC_PY8FSR_over_SC_PHOTOS.root"),"r");
    TH1D *hgUPCPY  = (TH1D*)f->Get("ratio_pTee_gammaUPC_PY8FSR_over_SC_PHOTOS");
    qed[i] = new TChain("output_tree");
   // qed[i]->Add(Form("/Users/pranatijana/Documents/LightByLight/Workspace/plottingMacros/After_Dielectron/WithGen/mc_qed_sl_forUnfolding_14thFeb.root"));
   //qed[i]->Add(Form("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/ForDielectron/QEDSC__Costheta.root"));//Old
   qed[i]->Add(Form("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/ForDielectron/QEDSC__Costheatstar_GenReco_6thMarch.root"));
   
    
    hdiele_Pt[i]    = new TH1D(Form("hdiele_Pt%s", Sample[i]),"",ptbin,0,1);
    hdiele_Rapidity[i]   = new TH1D(Form("hdiele_Rapidity%s",Sample[i]),"",rapbin,-2.2,2.2);
    //hdiele_Invmass[i]   = new TH1D(Form("hdiele_Invmass%s",Sample[i]),"",50,0,100);
    hdiele_Invmass[i]   = new TH1D(Form("hdiele_Invmass%s",Sample[i]),"",nMassbin, Massbin);
    hdiele_Costhetastar[i]   = new TH1D(Form("hdiele_Costhetastar%s",Sample[i]),"",10,0,1);

    hAcoplanarity[i]   = new TH1D(Form("hAcoplanarity%s",Sample[i]),"",50,0,0.1);
    
    hdiele_GenPt[i]         = new TH1D(Form("hdiele_GenPt%s", Sample[i]),"",ptbin,0,1);
    hdiele_GenRapidity[i]   = new TH1D(Form("hdiele_GenRapidity%s",Sample[i]),"",rapbin,-2.2,2.2);
    //hdiele_GenInvmass[i]          = new TH1D(Form("hdiele_GenInvmass%s",Sample[i]),"",50,0,100);
    hdiele_GenInvmass[i]          = new TH1D(Form("hdiele_GenInvmass%s",Sample[i]),"",nMassbin, Massbin);
    hdiele_GenCosthetastar[i]          = new TH1D(Form("hdiele_GenCosthetastar%s",Sample[i]),"",10,0,1);
    
    hdiele_RecoGenPt[i]         = new TH2D(Form("hdiele_RecoGenPt%s", Sample[i]),"",ptbin,0,1,ptbin,0,1);
    hdiele_RecoGenRapidity[i]   = new TH2D(Form("hdiele_RecoGenRapidity%s",Sample[i]),"",rapbin,-2.2,2.2,rapbin,-2.2,2.2);
    //hRecoGenInvmass[i]          = new TH2D(Form("hRecoGenInvmass%s",Sample[i]),"",50,0,100,50,0,100);
    hRecoGenInvmass[i]          = new TH2D(Form("hRecoGenInvmass%s",Sample[i]),"",nMassbin, Massbin,nMassbin, Massbin);
    hdiele_RecoGenCosthetastar[i]   = new TH2D(Form("hdiele_RecoGenCosthetastar%s",Sample[i]),"",10,0,1,10,0,1);

  // ===================   Cross section histograms   ======================     
 //TH1D* hUnfoldPt_data = new TH1D("hUnfoldPt_data", "Pt data unfolded",ptbin,0,1);
  TH1D* hGenPt_xSec              = new TH1D("hGenPt_xSec", "Pt gen MC",ptbin,0,1);
  TH1D* hGenRap_xSec             = new TH1D("hGenRap_xSec","Rap gen MC",rapbin,-2.2,2.2);
  //TH1D* hGenInvmass_xSec         = new TH1D("hGenInvmass_xSec","Invmass gen MC",50,0,100);
  TH1D* hGenInvmass_xSec         = new TH1D("hGenInvmass_xSec","Invmass gen MC",nMassbin, Massbin);
  TH1D* hGenCosthetastar_xSec         = new TH1D("hGenCosthetastar_xSec","Costhetastar gen MC",10,0,1);

  TH1D* hRecoMCPt_xSec              = new TH1D("hRecoMCPt_xSec", "Pt reco MC",ptbin,0,1);
  TH1D* hRecoMCRap_xSec             = new TH1D("hRecoMCRap_xSec","Rap reco MC",rapbin,-2.2,2.2);
  //TH1D* hRecoMCInvmass_xSec         = new TH1D("hRecoMCInvmass_xSec","Invmass reco MC",50,0,100);
  TH1D* hRecoMCInvmass_xSec         = new TH1D("hRecoMCInvmass_xSec","Invmass reco MC",nMassbin, Massbin);
  TH1D* hRecoMCCosthetastar_xSec         = new TH1D("hRecoMCCosthetastar_xSec","Costhetastar reco MC",10,0,1);

  TH1D* hRecoDataPt_xSec              = new TH1D("hRecoDataPt_xSec", "Pt reco Data",ptbin,0,1);
  TH1D* hRecoDataRap_xSec             = new TH1D("hRecoDataRap_xSec","Rap reco Data",rapbin,-2.2,2.2);
 // TH1D* hRecoDataInvmass_xSec         = new TH1D("hRecoDataInvmass_xSec","Invmass reco Data",50,0,100);
  TH1D* hRecoDataInvmass_xSec         = new TH1D("hRecoDataInvmass_xSec","Invmass reco Data",nMassbin, Massbin);
  TH1D* hRecoDataCosthetastar_xSec         = new TH1D("hRecoDataCosthetastar_xSec","Costhetastar reco Data",10,0,1);

  TH1D* hUnfoMCPt_xSec              = new TH1D("hUnfoMCPt_xSec", "Pt unfo MC",ptbin,0,1);
  TH1D* hUnfoMCRap_xSec             = new TH1D("hUnfoMCRap_xSec","Rap unfo MC",rapbin,-2.2,2.2);
  // TH1D* hUnfoMCInvmass_xSec         = new TH1D("hUnfoMCInvmass_xSec","Invmass unfo MC",50,0,100);
  // TH1D* hUnfoMCInvmass_Invert_xSec         = new TH1D("hUnfoMCInvmass_Invert_xSec","Invmass unfo Invert MC",50,0,100);
  TH1D* hUnfoMCInvmass_xSec         = new TH1D("hUnfoMCInvmass_xSec","Invmass unfo MC",nMassbin, Massbin);
  TH1D* hUnfoMCCosthetastar_xSec         = new TH1D("hUnfoMCCosthetastar_xSec","Costhetastar unfo MC",10,0,1);

  TH1D* hUnfoMCInvmass_Invert_xSec         = new TH1D("hUnfoMCInvmass_Invert_xSec","Invmass unfo Invert MC",nMassbin, Massbin);
  TH1D* hUnfoMCPt_Invert_xSec              = new TH1D("hUnfoMCPt_Invert_xSec", "Pt unfo Invert MC",ptbin,0,1);
  TH1D* hUnfoMCRap_Invert_xSec             = new TH1D("hUnfoMCRap_Invert_xSec","Rap unfo Invert MC",rapbin,-2.2,2.2);
  TH1D* hUnfoMCCosthetastar_Invert_xSec             = new TH1D("hUnfoMCCosthetastar_Invert_xSec","Costhetastar unfo Invert MC",rapbin,-2.2,2.2);

  TH1D* hUnfoDataPt_xSec              = new TH1D("hUnfoDataPt_xSec", "Pt unfo Data",ptbin,0,1);
  TH1D* hUnfoDataRap_xSec             = new TH1D("hUnfoDataRap_xSec","Rap unfo Data",rapbin,-2.2,2.2);
  // TH1D* hUnfoDataInvmass_xSec         = new TH1D("hUnfoDataInvmass_xSec","Invmass unfo Data",50,0,100);
  // TH1D* hUnfoDataInvmass_Invert_xSec         = new TH1D("hUnfoDataInvmass_Invert_xSec","Invmass unfo Data",50,0,100);
  TH1D* hUnfoDataInvmass_xSec         = new TH1D("hUnfoDataInvmass_xSec","Invmass unfo Data",nMassbin, Massbin);
   TH1D* hUnfoDataCosthetastar_xSec         = new TH1D("hUnfoDataCosthetastar_xSec","Costhetastar unfo Data",10,0,1);

  TH1D* hUnfoDataInvmass_Invert_xSec         = new TH1D("hUnfoDataInvmass_Invert_xSec","Invmass unfo Invert Data",nMassbin, Massbin);
  TH1D* hUnfoDataPt_Invert_xSec              = new TH1D("hUnfoDataPt_Invert_xSec", "Pt unfo Invert Data",ptbin,0,1);
  TH1D* hUnfoDataRap_Invert_xSec             = new TH1D("hUnfoDataRap_Invert_xSec","Rap unfo Invert Data",rapbin,-2.2,2.2);
  TH1D* hUnfoDataCosthetastar_Invert_xSec             = new TH1D("hUnfoDataCosthetastar_Invert_xSec","Costhetastar unfo Invert Data",10,0,1);
   
    RooUnfoldResponse responsePt (hdiele_Pt[i], hdiele_GenPt[i]);
   // RooUnfoldResponse responsePt (hdiele_Pt[i], hdiele_GenPt[i],hdiele_RecoGenPt[i]);
    RooUnfoldResponse responseRapidity (hdiele_Rapidity[i], hdiele_GenRapidity[i],hdiele_RecoGenRapidity[i]);
    RooUnfoldResponse responseInvmass (hdiele_Invmass[i], hdiele_GenInvmass[i],hRecoGenInvmass[i]);
    RooUnfoldResponse responseCosthetastar (hdiele_Costhetastar[i], hdiele_GenCosthetastar[i], hdiele_RecoGenCosthetastar[i]);
    cout << "file " << qed[i]->GetEntries()  << endl;
    ReadqedTree  qedR(qed[i]);
    qedR.fChain->SetBranchStatus("*",1);

    if (qedR.fChain == 0) return;
    
    Long64_t nentries = qedR.fChain->GetEntriesFast();
    cout << qed[i]->GetName() << "    " << nentries << endl;
    
    Long64_t nbytes = 0, nb = 0;
    //for (Long64_t jentry=0; jentry<1000;jentry++) {
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry_evt = qedR.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = qedR.fChain->GetEntry(jentry);   nbytes += nb;
     

      // "===================      Fill gen Info here ======================     
     
     if(abs(qedR.gen_Eta1) > 2.2 || abs(qedR.gen_Eta2) > 2.2 || qedR.gen_Pt1 < 2.0 || qedR.gen_Pt2 < 2.0 || qedR.gen_diEle_M < 5 || qedR.gen_diEle_Pt > 1) continue;
      if(abs(qedR.gen_diEle_Rapidity) > 2.2) continue;
      double ele_dphi = GetDeltaPhi(qedR.gen_Phi1, qedR.gen_Phi2);
      double ele_gen_acop = 1 - (GetDeltaPhi(qedR.gen_Phi1, qedR.gen_Phi2)/3.141592653589);
      if(ele_gen_acop > 0.01) continue; 

       // cout << "BinNumber1:" <<hdiele_GenPt[i]->GetXaxis()->FindFixBin(qedR.gen_diEle_Pt) << endl;
      int binNumber1 = hdiele_GenPt[i]->GetXaxis()->FindFixBin(qedR.gen_diEle_Pt);
       // double reweight1 = hgUPCPY->GetBinContent(binNumber1);
       // cout << "reweight:" << reweight1 << endl;
      double reweight1 =1;
      hdiele_GenPt[i]->Fill(qedR.gen_diEle_Pt,W1*reweight1);
      
      if(absRap)hdiele_GenRapidity[i]->Fill(abs(qedR.gen_diEle_Rapidity),W1);
      if(!absRap)hdiele_GenRapidity[i]->Fill(qedR.gen_diEle_Rapidity,W1);
      hdiele_GenInvmass[i]->Fill(qedR.gen_diEle_M,W1);
      hdiele_GenCosthetastar[i]->Fill(abs(qedR.gen_diEle_Costhetastar),W1);

      // ===================   Apply analysis cuts on reco MC ======================           
      bool isPassing =  qedR.ok_neuexcl == 1 && qedR.ok_chexcl_extrk == 1 
                        && qedR.vSum_M > 5 && qedR.vSum_Pt < 1 
                        && abs(qedR.eleEta_1) < 2.2 && abs(qedR.eleEta_2) < 2.2 && qedR.elePt_1 > 2.0 && qedR.elePt_2 > 2.0
                        
                        && qedR.ele_acop < 0.01;

  if(qedR.ok_trigger == 1) { //trigger 
	if(isPassing){

    double  ET1 = qedR.eleEta_1;
    double  ET2 = qedR.eleEta_2;
    //SC
    double lumiNorm_SC = RecoID_SF_Electron(ET1) * RecoID_SF_Electron(ET2) * norm_LbLSC_other;  
    std::vector<double> wt = {lumiNorm_SC };
 
   //For SL
    // double lumiNorm_SL = RecoID_SF_Electron(ET1) * RecoID_SF_Electron(ET2) * norm_LbLSL_other;  
    // std::vector<double> wt = {lumiNorm_SL};
  
 // cout << "BinNumber2:" <<hdiele_GenPt[i]->GetXaxis()->FindFixBin(qedR.gen_diEle_Pt) << endl;
  int binNumber2 = hdiele_GenPt[i]->GetXaxis()->FindFixBin(qedR.gen_diEle_Pt);
  //double reweight2 = hgUPCPY->GetBinContent(binNumber2);
  //cout << "reweight2:" << reweight2 << endl;
  double reweight2=1;

	hdiele_Pt[i]->Fill(qedR.vSum_Pt,wt[i]);//reweighting to gUPC
  
	if(absRap)hdiele_Rapidity[i]->Fill(abs(qedR.vSum_Rapidity),wt[i]);
	if(!absRap)hdiele_Rapidity[i]->Fill(qedR.vSum_Rapidity,wt[i]);
	hdiele_Invmass[i]->Fill(qedR.vSum_M,wt[i]);
  hdiele_Costhetastar[i]->Fill(abs(qedR.costhetastar),wt[i]);
	
  
	hdiele_RecoGenPt[i]->Fill(qedR.vSum_Pt, qedR.gen_diEle_Pt*reweight2, wt[i]);//reweighting to gUPC

	if(absRap)hdiele_RecoGenRapidity[i]->Fill(abs(qedR.vSum_Rapidity),abs(qedR.gen_diEle_Rapidity),wt[i]);
	if(!absRap)hdiele_RecoGenRapidity[i]->Fill(qedR.vSum_Rapidity,qedR.gen_diEle_Rapidity,wt[i]);
	hRecoGenInvmass[i]->Fill(qedR.vSum_M,qedR.gen_diEle_M,wt[i]);
  hdiele_RecoGenCosthetastar[i]->Fill(abs(qedR.costhetastar),abs(qedR.gen_diEle_Costhetastar),wt[i]);
  
  responsePt.Fill(qedR.vSum_Pt,qedR.gen_diEle_Pt*reweight2,wt[i]); //reweighting to gUPC
  
  if(absRap)responseRapidity.Fill(abs(qedR.vSum_Rapidity),abs(qedR.gen_diEle_Rapidity),wt[i]);
	if(!absRap)responseRapidity.Fill(qedR.vSum_Rapidity,qedR.gen_diEle_Rapidity,wt[i]);
	responseInvmass.Fill(qedR.vSum_M,qedR.gen_diEle_M,wt[i]);
  responseCosthetastar.Fill(abs(qedR.costhetastar),abs(qedR.gen_diEle_Costhetastar),wt[i]);
	} //isPassing analysis cuts
	else {
  //cout << "BinNumber3:" <<hdiele_GenPt[i]->GetXaxis()->FindFixBin(qedR.gen_diEle_Pt) << endl;
  int binNumber3 = hdiele_GenPt[i]->GetXaxis()->FindFixBin(qedR.gen_diEle_Pt);
  //double reweight3 = hgUPCPY->GetBinContent(binNumber3);
  double reweight3 =1;
	responsePt.Miss(qedR.gen_diEle_Pt,W1*reweight3);
  
	if(absRap)responseRapidity.Miss(abs(qedR.gen_diEle_Rapidity),W1);
	if(!absRap)responseRapidity.Miss(qedR.gen_diEle_Rapidity,W1);
	responseInvmass.Miss(qedR.gen_diEle_M,W1);
  responseCosthetastar.Miss(abs(qedR.gen_diEle_Costhetastar),W1);
	}
      } //trigger
	else {
  int binNumber4 = hdiele_GenPt[i]->GetXaxis()->FindFixBin(qedR.gen_diEle_Pt);
  //double reweight4 = hgUPCPY->GetBinContent(binNumber4);
  double reweight4 =1;
	responsePt.Miss(qedR.gen_diEle_Pt, W1*reweight4);

	if(absRap)responseRapidity.Miss(abs(qedR.gen_diEle_Rapidity),W1);
	if(!absRap)responseRapidity.Miss(qedR.gen_diEle_Rapidity,W1);
	responseInvmass.Miss(qedR.gen_diEle_M,W1);
  responseCosthetastar.Miss(abs(qedR.gen_diEle_Costhetastar),W1);
	}
      
      
    } //entry
 // } // for 4 files

// TCanvas*C =  new TCanvas();
 //hdiele_RecoGenRapidity[0]->Draw("Colz");
 //C->Print("fig_Bayes_Inverse/Invmass.pdf");
cout << "gen entries pass:" << hdiele_GenPt[0]->GetEntries() << endl; 
cout << "reco entries pass:" << hdiele_Pt[0]->GetEntries() << endl; 
cout << "gen N bin :" << hdiele_GenPt[0]->GetNbinsX() << endl; 
cout << "reco N bin:" << hdiele_Pt[0]->GetNbinsX() << endl;
/*******************************Condition*number*********************************/
cout << "Inverse1" << endl;
TDecompSVD *svd= new TDecompSVD (responsePt.Mresponse());
auto singular_values = svd->GetSig();
//cout << "Test1" << endl;
svd->Print();
double singular_value_min;
      for (int i = 0; i < 7 ; i++)
      {
        //cout << "Test2" << endl;
        cout << singular_values[i] << endl;
      }

//cout << "condition number: "  << singular_values[0]/singular_value_min << endl;
cout << "Inverse2" << endl;
/*********************************************************************************************/
  cout<< "======================================Response matrix========================="<<endl;
  char MethID2[100]; 
  //if(method==1)
   // {

      RooUnfoldBayes unfold(&responsePt, hdiele_Pt[0], 3);
      RooUnfoldBayes unfold2(&responseRapidity, hdiele_Rapidity[0], itr);
      RooUnfoldBayes unfold3(&responseInvmass, hdiele_Invmass[0], 1);
      RooUnfoldBayes unfold4(&responseCosthetastar, hdiele_Costhetastar[0], 1);

      RooUnfoldInvert unfold_Invert(&responsePt, hdiele_Pt[0]);
      RooUnfoldInvert unfold2_Invert(&responseRapidity, hdiele_Rapidity[0]);
      RooUnfoldInvert unfold3_Invert(&responseInvmass, hdiele_Invmass[0]);
      RooUnfoldInvert unfold4_Invert(&responseCosthetastar, hdiele_Costhetastar[0]);


/********************Matrix-plotting************/
        TCanvas*c4 =  new TCanvas("c4", "" ,800,600);
        c4->SetRightMargin(0.15);
        TMatrixD  responseMatrix1 = responsePt.Mresponse();
        TH2D *hResponseMatrix = new TH2D("hResponseMatrix", "Response Matrix ;Reco dielectron p_{T} bin index;Gen dielectron p_{T} bin index", responseMatrix1.GetNrows(), 0, responseMatrix1.GetNrows(), responseMatrix1.GetNcols(), 0, responseMatrix1.GetNcols());
        for (int i = 0; i < responseMatrix1.GetNrows(); ++i) {
         // cout << "rows:" <<  responseMatrix1.GetNrows() << endl;
             for (int j = 0; j < responseMatrix1.GetNcols(); ++j) {
                  hResponseMatrix->SetBinContent(i + 1, j + 1, responseMatrix1[i][j]);
          //        cout << "cols:" <<  responseMatrix1.GetNcols() << endl;
                  }
        }
        hResponseMatrix->GetXaxis()->SetNdivisions(5);
        hResponseMatrix->GetYaxis()->SetNdivisions(5);
        hResponseMatrix->Draw("textcolz");

        c4->Update();
        
        c4->SaveAs("fig_AN/Pt_SC_ele.pdf");
        //
        TCanvas*c5 =  new TCanvas("c5","Response Matrix2",800,600);
        c5->SetRightMargin(0.15);
        TMatrixD  responseMatrix2 = responseRapidity.Mresponse();
        TH2D *hResponseMatrix2 = new TH2D("hResponseMatrix2", "Response Matrix2;Reco dielectron rapidity bin index;Gen dielectron rapidity bin index", responseMatrix2.GetNrows(), 0, responseMatrix2.GetNrows(), responseMatrix2.GetNcols(), 0, responseMatrix2.GetNcols());
        for (int i = 0; i < responseMatrix2.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix2.GetNcols(); ++j) {
                  hResponseMatrix2->SetBinContent(i + 1, j + 1, responseMatrix2[i][j]);
                  }
        }
        hResponseMatrix2->GetXaxis()->SetNdivisions(4);
        hResponseMatrix2->GetYaxis()->SetNdivisions(4);
        hResponseMatrix2->Draw("textcolz");
        TLatex *mark6 = new TLatex();
        mark6->SetTextSize(0.035); 
        mark6->SetTextFont(42);
        // mark6->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
        // mark6->DrawLatex(0.65, 0.97, "#scale[0.8]{1.647 nb^{-1} (5.02 TeV)}");

        c5->Update();
        c5->SaveAs("fig_AN/Rapidity_SC_ele.pdf");
        //
        TCanvas*c6 =  new TCanvas("c6","Response Matrix3", 800,600);
        c6->SetRightMargin(0.15);
        TMatrixD  responseMatrix3 = responseInvmass.Mresponse();
        TH2D *hResponseMatrix3 = new TH2D("hResponseMatrix3", "Response Matrix3;Reco dielectron inv. mass bin index;Gen dielectron inv. mass bin index", responseMatrix3.GetNrows(), 0, responseMatrix3.GetNrows(), responseMatrix3.GetNcols(), 0, responseMatrix3.GetNcols());
        for (int i = 0; i < responseMatrix3.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix3.GetNcols(); ++j) {
                  hResponseMatrix3->SetBinContent(i + 1, j + 1, responseMatrix3[i][j]);
                  }
        }
        TLatex *mark8 = new TLatex();
      mark8->SetTextSize(0.035); 
      mark8->SetTextFont(42);
 
       mark8->DrawLatex(0.16, 0.97, "#bf{CMS} #it{Preliminary}");
       mark8->DrawLatex(0.65, 0.97, "#scale[0.8]{1.647 nb^{-1} (5.02 TeV)}");
       hResponseMatrix3->Draw("textcolz");
       c6->Update();
       c6->SaveAs("fig_AN/InvMass_SC_ele.pdf");
//
 TCanvas*c7 =  new TCanvas("c7","Response Matrix4", 800,600);
        c7->SetRightMargin(0.15);
        TMatrixD  responseMatrix4 = responseCosthetastar.Mresponse();
        TH2D *hResponseMatrix4 = new TH2D("hResponseMatrix4", "Response Matrix4;Reco dielectron Costhetatstar bin index;Gen dielectron Costhetastar bin index", responseMatrix4.GetNrows(), 0, responseMatrix4.GetNrows(), responseMatrix4.GetNcols(), 0, responseMatrix4.GetNcols());
        for (int i = 0; i < responseMatrix4.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix4.GetNcols(); ++j) {
                  hResponseMatrix4->SetBinContent(i + 1, j + 1, responseMatrix4[i][j]);
                  }
        }
      
       hResponseMatrix4->Draw("textcolz");
       c7->Update();
       c7->SaveAs("fig_AN/Costhetastar_SC_ele.pdf");

        
/***********************************************/
  
   
 // unfold.IncludeSystematics(2);
  TH1D* hUnfoldPt= (TH1D*) unfold.Hunfold();
  unfold.PrintTable (cout,hdiele_GenPt[0]);
  for(int j=1; j <=hdiele_GenPt[0]->GetNbinsX(); j ++){
  cout << "GetBinentries gen:" << hdiele_GenPt[0]->GetBinContent(j) << endl;
  }
  for(int j=1; j <=hdiele_Pt[0]->GetNbinsX(); j ++){
  cout << "GetBinentries reco:" << hdiele_Pt[0]->GetBinContent(j) << endl;
  }
  TH1D* hfold= (TH1D*) responsePt.ApplyToTruth(hdiele_GenPt[0],"");
  cout << "hfold entries :" << hfold->GetEntries() << endl; 
  
  //unfold2.IncludeSystematics(2);  
  TH1D* hUnfoldRapidity= (TH1D*) unfold2.Hunfold();
  unfold2.PrintTable (cout,hdiele_GenRapidity[0]);
  
  //unfold3.IncludeSystematics(2); 
  TH1D* hUnfoldInvmass= (TH1D*) unfold3.Hunfold();
  unfold3.PrintTable (cout,hdiele_GenInvmass[0]);

  TH1D* hUnfoldCosthetastar= (TH1D*) unfold4.Hunfold();
  unfold4.PrintTable (cout,hdiele_GenCosthetastar[0]);

  cout <<"Bayes MC finished" << endl;

  //unfold_Invert.IncludeSystematics(2);
  TH1D* hUnfoldPt_Invert= (TH1D*) unfold_Invert.Hunfold();
  unfold_Invert.PrintTable (cout,hdiele_GenPt[0]);
  //unfold2_Invert.IncludeSystematics(2);
  TH1D* hUnfoldRapidity_Invert= (TH1D*) unfold2_Invert.Hunfold();
  unfold2_Invert.PrintTable (cout,hdiele_GenRapidity[0]);
 // unfold3_Invert.IncludeSystematics(2);
  TH1D* hUnfoldInvmass_Invert= (TH1D*) unfold3_Invert.Hunfold();
   unfold3_Invert.PrintTable (cout,hdiele_GenInvmass[0]);

  TH1D* hUnfoldCosthetastar_Invert= (TH1D*) unfold4_Invert.Hunfold();
  unfold4_Invert.PrintTable (cout,hdiele_GenCosthetastar[0]);
  /***********************************************Chi2************************************************/
   


  /***************************  Start reading data and subtract background contributions ********************************/

  //TFile *f1=new TFile("dieleton_histos.root","r"); 
  //TFile *f1=new TFile("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/ForDielectron/output_new_ZDCAND4n4n.root ","r"); 
//TFile *f1=new TFile("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/ForDielectron/output_ZDCAND3n.root","r");
//TFile *f1=new TFile("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/ForDielectron/output_ZDCAND3n_costhetastar.root","r");
//TFile *f1=new TFile("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/ForDielectron/Data_sumpt_6bin.root","r");
//TFile *f1=new TFile("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/ForDielectron/Data_sumpt_7bin.root","r");
 TFile *f1=new TFile("/Users/pranatijana/Documents/LightByLight/Workspace/plottingMacros/Data_rap_11bin.root","r");
 
  TH1D *hdiele_Pt_data        = (TH1D*)f1->Get("hSumPt_unfoldData");
  TH1D *hdiele_Rapidity_data  = (TH1D*)f1->Get("hRap_unfoldData"); 
  TH1D *hdiele_Invmass_data   = (TH1D*)f1->Get("hInvmass_unfoldData"); 
  TH1D *hdiele_Costhetastar_data   = (TH1D*)f1->Get("hCosthetastar_unfoldData"); 
      
      RooUnfoldBayes unfold_dataPt(&responsePt, hdiele_Pt_data, 100); // unfold data
      RooUnfoldBayes unfold_dataRapidity(&responseRapidity, hdiele_Rapidity_data, 1);
      RooUnfoldBayes unfold_dataInvmass(&responseInvmass, hdiele_Invmass_data, 1);
      RooUnfoldBayes unfold_dataCosthetastar(&responseCosthetastar, hdiele_Costhetastar_data, 1);
  
      RooUnfoldInvert unfold_dataPt_Invert(&responsePt, hdiele_Pt_data);  
      RooUnfoldInvert unfold_dataRapidity_Invert(&responseRapidity, hdiele_Rapidity_data);
      RooUnfoldInvert unfold_dataInvmass_Invert(&responseInvmass, hdiele_Invmass_data);  
      RooUnfoldInvert unfold_dataCosthetastar_Invert(&responseCosthetastar, hdiele_Costhetastar_data); 
  
  TH1D* hUnfoldPt_data= (TH1D*) unfold_dataPt.Hunfold();
  unfold_dataPt.PrintTable (cout, hdiele_GenPt[0]);
 

   TH1D* hUnfoldRapidity_data= (TH1D*) unfold_dataRapidity_Invert.Hunfold();
   //unfold_dataRapidity_Invert.PrintTable (cout, hdiele_GenRapidity[0]);

   TH1D* hUnfoldInvmass_data= (TH1D*) unfold_dataInvmass.Hunfold();
   unfold_dataInvmass.PrintTable (cout, hdiele_GenInvmass[0]);

   TH1D* hUnfoldCosthetastar_data= (TH1D*) unfold_dataCosthetastar.Hunfold();
   unfold_dataCosthetastar.PrintTable (cout, hdiele_GenCosthetastar[0]);

 
  TH1D* hUnfoldPt_Invert_data= (TH1D*) unfold_dataPt_Invert.Hunfold();
  //unfold_dataPt_Invert.PrintTable (cout, hdiele_GenPt[0]);

  TH1D* hUnfoldRapidity_Invert_data= (TH1D*) unfold_dataRapidity_Invert.Hunfold();
  unfold_dataRapidity_Invert.PrintTable (cout, hdiele_GenRapidity[0]);
  
   TH1D* hUnfoldInvmass_Invert_data= (TH1D*) unfold_dataInvmass_Invert.Hunfold();
   //unfold_dataInvmass_Invert.PrintTable (cout, hdiele_GenInvmass[0]);

   TH1D* hUnfoldCosthetastar_Invert_data= (TH1D*) unfold_dataCosthetastar_Invert.Hunfold();
    unfold_dataCosthetastar_Invert.PrintTable (cout, hdiele_GenCosthetastar[0]);
  /**********************Covarinace matrix*******************************/

  TCanvas*c8 =  new TCanvas("c8","covarianceMatrixPt", 800,600);      
        c8->SetRightMargin(0.15);
        c8->SetTopMargin(0.15);
        TMatrixD  covarianceMatrix_Pt = unfold_dataPt.Eunfold();
        TH2D *hCov_Pt = new TH2D("hCov_Pt", "covarianceMatrix_Pt ;p_{T}^{ee} - Bin index;p_{T}^{ee} -  Bin index", covarianceMatrix_Pt.GetNrows(), 0, covarianceMatrix_Pt.GetNrows(), covarianceMatrix_Pt.GetNcols(), 0, covarianceMatrix_Pt.GetNcols());
        for (int i = 0; i < covarianceMatrix_Pt.GetNrows(); ++i) {
             for (int j = 0; j < covarianceMatrix_Pt.GetNcols(); ++j) {
                  hCov_Pt->SetBinContent(i + 1, j + 1,covarianceMatrix_Pt[i][j]);
                  }
        }
       hCov_Pt->SetTitle("Covariance Matrix - Pt");
       hCov_Pt->Draw("textcolz");
       c8->Update();
       c8->SaveAs("fig_AN/Cov_Pt_ele_SC.pdf");

  TCanvas*c9 =  new TCanvas("c9","covarianceMatrixRapidity", 800,600);
        c9->SetRightMargin(0.15);
        TMatrixD  covarianceMatrix_Rapidity = unfold_dataRapidity_Invert.Eunfold();
        TH2D *hCov_Rapidity = new TH2D("hCov_Rapidity", "covarianceMatrix_Rapidity;y_{ee} - Bin index ; y_{ee} - Bin index", covarianceMatrix_Rapidity.GetNrows(), 0, covarianceMatrix_Rapidity.GetNrows(), covarianceMatrix_Rapidity.GetNcols(), 0, covarianceMatrix_Rapidity.GetNcols());
        for (int i = 0; i < covarianceMatrix_Rapidity.GetNrows(); ++i) {
             for (int j = 0; j < covarianceMatrix_Rapidity.GetNcols(); ++j) {
                  hCov_Rapidity->SetBinContent(i + 1, j + 1,covarianceMatrix_Rapidity[i][j]);
                  }
        }
       hCov_Rapidity->SetTitle("Covariance Matrix - Rapidity");
       hCov_Rapidity->Draw("textcolz");
       c9->Update();
       c9->SaveAs("fig_AN/Cov_Rapidity_ele_SC.pdf");
  
  TCanvas*c10 =  new TCanvas("c10","covarianceMatrixInvmass", 800,600);
        c10->SetRightMargin(0.15);
        TMatrixD  covarianceMatrix_Invmass = unfold_dataInvmass_Invert.Eunfold();
        TH2D *hCov_Invmass= new TH2D("hCov_Invmass", "covarianceMatrix_Invmass;m_{ee}- Bin index;m_{ee}- Bin index", covarianceMatrix_Invmass.GetNrows(), 0, covarianceMatrix_Invmass.GetNrows(), covarianceMatrix_Invmass.GetNcols(), 0, covarianceMatrix_Invmass.GetNcols());
        for (int i = 0; i < covarianceMatrix_Invmass.GetNrows(); ++i) {
             for (int j = 0; j < covarianceMatrix_Invmass.GetNcols(); ++j) {
                  hCov_Invmass->SetBinContent(i + 1, j + 1,covarianceMatrix_Invmass[i][j]);
                  }
        }
       hCov_Invmass->SetTitle("Covariance Matrix - Invarinat mass");
       hCov_Invmass->Draw("textcolz");
       c10->Update();
       c10->SaveAs("fig_AN/Cov_Invmass_ele_SC.pdf");

   TCanvas*c11 =  new TCanvas("c11","covarianceMatrixCosthetastar", 800,600);
        c11->SetRightMargin(0.15);
        TMatrixD  covarianceMatrix_Costhetastar = unfold_dataCosthetastar_Invert.Eunfold();
        TH2D *hCov_Costhetastar= new TH2D("hCov_Costhetastar", "covarianceMatrix_Costhetastar; |cos#theta^{*}|_{ee} - Bin index; |cos#theta^{*}|_{ee} - Bin index", covarianceMatrix_Costhetastar.GetNrows(), 0, covarianceMatrix_Costhetastar.GetNrows(), covarianceMatrix_Costhetastar.GetNcols(), 0, covarianceMatrix_Costhetastar.GetNcols());
        for (int i = 0; i < covarianceMatrix_Costhetastar.GetNrows(); ++i) {
             for (int j = 0; j < covarianceMatrix_Costhetastar.GetNcols(); ++j) {
                  hCov_Costhetastar->SetBinContent(i + 1, j + 1,covarianceMatrix_Costhetastar[i][j]);
                
        }
        }
      hCov_Costhetastar->SetTitle("Covariance Matrix - Costhetastar");
    
       hCov_Costhetastar->Draw("textcolz");
       c11->Update();
       c11->SaveAs("fig_AN/Cov_Costhetastar_ele_SC.pdf");


   

  /************************************************************************/
  cout << "Invert2:" << endl;
  hGenPt_xSec  = getXSecHist(hdiele_GenPt[0], "hGenPt_xSec", ptbin, 0, 1 ,1 );
  hGenRap_xSec = getXSecHist(hdiele_GenRapidity[0], "hGenRap_xSec",rapbin,-2.2,2.2,1 );
  // hGenInvmass_xSec = getXSecHist(hdiele_GenInvmass[0], "hGenInvmass_xSec" ,50,0,100,1 );
  hGenInvmass_xSec = getInvXSecHist(hdiele_GenInvmass[0], "hGenInvmass_xSec" ,1 );
  hGenCosthetastar_xSec = getXSecHist(hdiele_GenCosthetastar[0], "hGenCosthetastar_xSec" ,10,0,1 ,1);

  hRecoMCPt_xSec  = getXSecHist(hdiele_Pt[0], "hRecoMCPt_xSec", ptbin, 0, 1 ,1 );
  hRecoMCRap_xSec = getXSecHist(hdiele_Rapidity[0], "hRecoMCRap_xSec", rapbin,-2.2,2.2,1 );
  //hRecoMCInvmass_xSec = getXSecHist(hdiele_Invmass[0], "hRecoMCInvmass_xSec",50,0,100,1 );
  hRecoMCInvmass_xSec = getInvXSecHist(hdiele_Invmass[0], "hRecoMCInvmass_xSec",1 );
  hRecoMCCosthetastar_xSec = getXSecHist(hdiele_Costhetastar[0], "hRecoMCCosthetastar_xSec",10,0,1,1 );

  hRecoDataPt_xSec  = getXSecHist(hdiele_Pt_data, "hRecoDataPt_xSec", ptbin, 0, 1, 1.64718);
  hRecoDataRap_xSec = getXSecHist(hdiele_Rapidity_data, "hRecoDataRap_xSec",rapbin,-2.2,2.2, 1.64718);
  //hRecoDataInvmass_xSec = getXSecHist(hdiele_Invmass_data, "hRecoDataInvmass_xSec", 50,0,100,1.647);
  hRecoDataInvmass_xSec = getInvXSecHist(hdiele_Invmass_data, "hRecoDataInvmass_xSec", 1.64718);
  hRecoDataCosthetastar_xSec = getXSecHist(hdiele_Costhetastar_data, "hRecoDataCosthetastar_xSec",10,0,1, 1.64718);

  hUnfoMCPt_xSec  = getXSecHist(hUnfoldPt,       "hUnfoMCPt_xSec", ptbin, 0, 1,1 );
  hUnfoMCRap_xSec = getXSecHist(hUnfoldRapidity, "hUnfoMCRap_xSec",rapbin,-2.2,2.2, 1 );
  hUnfoMCInvmass_xSec = getInvXSecHist(hUnfoldInvmass, "hUnfoMCInvmass_xSec", 1 );
  hUnfoMCCosthetastar_xSec  = getXSecHist(hUnfoldCosthetastar,       "hUnfoMCCosthetastar_xSec", 10, 0, 1,1 );  

  hUnfoMCInvmass_Invert_xSec = getInvXSecHist(hUnfoldInvmass_Invert, "hUnfoMCInvmass_Invert_xSec", 1 );
  hUnfoMCPt_Invert_xSec  = getXSecHist(hUnfoldPt_Invert,       "hUnfoMCPt_Invert_xSec", ptbin, 0, 1,1 );
  hUnfoMCRap_Invert_xSec = getXSecHist(hUnfoldRapidity_Invert, "hUnfoMCRap_Invert_xSec", rapbin,-2.2,2.2,1); 
  hUnfoMCCosthetastar_Invert_xSec = getXSecHist(hUnfoldCosthetastar_Invert, "hUnfoMCCosthetastar_Invert_xSec", 10,0,1,1); 
 

  hUnfoDataPt_xSec  = getXSecHist(hUnfoldPt_data, "hUnfoDataPt_xSec", ptbin, 0, 1, 1.64718);
   //hUnfoldPt_data->Scale(1/(1.647*0.05)); 
  hUnfoDataRap_xSec = getXSecHist(hUnfoldRapidity_data, "hUnfoDataRap_xSec",rapbin,-2.2,2.2, 1.64718);
  hUnfoDataInvmass_xSec = getInvXSecHist(hUnfoldInvmass_data, "hUnfoDataInvmass_xSec", 1.64718);
  hUnfoDataCosthetastar_xSec  = getXSecHist(hUnfoldCosthetastar_data, "hUnfoDataCosthetastar_xSec", 10, 0, 1, 1.64718);

  hUnfoDataInvmass_Invert_xSec = getInvXSecHist(hUnfoldInvmass_Invert_data, "hUnfoDataInvmass_Invert_xSec", 1.64718);
  hUnfoDataPt_Invert_xSec  = getXSecHist(hUnfoldPt_Invert_data, "hUnfoDataPt_Invert_xSec", ptbin, 0, 1, 1.64718);
  hUnfoDataRap_Invert_xSec = getXSecHist(hUnfoldRapidity_Invert_data, "hUnfoDataRap_Invert_xSec",rapbin,-2.2,2.2,1.64718);
  hUnfoDataCosthetastar_Invert_xSec  = getXSecHist(hUnfoldCosthetastar_Invert_data, "hUnfoDataCosthetastar_Invert_xSec", 10, 0, 1, 1.64718);
 
  //hUnfoDataRap_xSec->Draw("p");

 
 //hUnfoDataInvmass_Invert_xSec->Draw();

 

/********************/
  cout << " invariant mass reco " << Sample[0] <<  " :" << hdiele_Invmass[0]->Integral() << endl;  
  cout << " invariant mass Gen " << Sample[0] <<  " :" << hdiele_GenInvmass[0]->Integral() << endl;  
  cout << " Rapidity reco " << Sample[0] <<  " :" << hdiele_Rapidity[0]->Integral() << endl;  
  cout << " Pt reco " << Sample[0] <<  " :" << hdiele_Pt[0]->Integral() << endl;  
  cout << " invariant mass data" << Sample[0] <<  " :" << hdiele_Invmass_data->Integral() << endl;
  cout << " Rapidity data " << Sample[0] <<  " :" << hdiele_Rapidity_data->Integral() << endl;
  
  cout << "Cross section Pt Data:" <<  hRecoDataPt_xSec->Integral() << endl;
  cout << "Cross section Rapidity Data:" <<  hRecoDataRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass Data:" << hRecoDataInvmass_xSec->Integral()<< endl;
  cout << "Cross section Pt Reco MC:" <<   hRecoMCPt_xSec->Integral() << endl;
  cout << "Cross section Rapidity Reco MC:" <<  hRecoMCRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass Reco MC:" << hRecoMCInvmass_xSec->Integral()<< endl;
  cout << "Cross section Pt Gen MC:" <<   hGenPt_xSec->Integral() << endl;
  cout << "Cross section Rapidity Gen MC:" <<  hGenRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass Gen  MC:" << hGenInvmass_xSec->Integral()<< endl;
  cout << "Cross section Pt UnfoldedData:" <<   hUnfoDataPt_xSec->Integral() << endl;
  //cout << "Cross section Rapidity UnfoldedData:" <<  hUnfoDataRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass UnfoldedData:" << hUnfoDataInvmass_xSec->Integral()<< endl;
  cout << "Cross section Costhetastar UnfoldedData:" << hUnfoDataCosthetastar_xSec->Integral()<< endl;
  cout << "Cross section Pt UnfoldedData_Invert:" <<   hUnfoDataPt_Invert_xSec->Integral() << endl;
  cout << "Cross section Rapidity UnfoldedData_Invert:" <<  hUnfoDataRap_Invert_xSec->Integral()<< endl;
  cout << "Cross section Invmass UnfoldedData_Invert:" << hUnfoDataInvmass_Invert_xSec->Integral()<< endl;
   cout << "Cross section Costhetastar UnfoldedData_Invert:" <<  hUnfoDataCosthetastar_Invert_xSec->Integral()<< endl;
  
  
  //outf->Write();
//Unfold
double chi2_pt, chi2_mass, chi2_rap, chi2_pt_smear, chi2_mass_smear, chi2_rap_smear;
 for(int j =0; j <= hGenPt_xSec->GetNbinsX()+1; j++ ) {
     if(hGenPt_xSec->GetBinContent(j) > 0){
      chi2_pt+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinError(j) * hGenPt_xSec->GetBinError(j));
     // chi2_pt+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinContent(j));
      //chi2+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinContent(j));
     }
 }
 cout << "Chi2_Pt:" << chi2_pt << endl;
//smear
 for(int j =0; j <= hRecoMCPt_xSec->GetNbinsX()+1; j++ ) {
     if(hRecoMCPt_xSec->GetBinContent(j) > 0){
     // chi2_pt_smear+= ((hRecoMCPt_xSec->GetBinContent(j)-hRecoDataPt_xSec->GetBinContent(j)) * (hRecoMCPt_xSec->GetBinContent(j)-hRecoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinError(j)*hGenPt_xSec->GetBinError(j) + hRecoDataPt_xSec->GetBinError(j)*hRecoDataPt_xSec->GetBinError(j));
       chi2_pt_smear+= ((hRecoMCPt_xSec->GetBinContent(j)-hRecoDataPt_xSec->GetBinContent(j)) * (hRecoMCPt_xSec->GetBinContent(j)-hRecoDataPt_xSec->GetBinContent(j)))/( hRecoMCPt_xSec->GetBinError(j) * hRecoMCPt_xSec->GetBinError(j) );
      //chi2+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinContent(j));
     }
 }
 cout << "Chi2_pt_smear:" << chi2_pt_smear << endl;


//Unfold
 for(int j =0; j <= hGenInvmass_xSec->GetNbinsX(); j++ ) {
     if(hGenInvmass_xSec->GetBinContent(j) > 0){
       chi2_mass+= ((hGenInvmass_xSec->GetBinContent(j)-hUnfoDataInvmass_xSec->GetBinContent(j)) * (hGenInvmass_xSec->GetBinContent(j)-hUnfoDataInvmass_xSec->GetBinContent(j)))/(hGenInvmass_xSec->GetBinError(j)*hUnfoDataInvmass_xSec->GetBinError(j));
      //chi2+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinContent(j));
     }
 }
 cout << "Chi2_Mass:" << chi2_mass << endl;
//smear
 for(int j =0; j <= hRecoMCInvmass_xSec->GetNbinsX(); j++ ) {
     if(hRecoMCInvmass_xSec->GetBinContent(j) > 0){
       chi2_mass_smear+= ((hRecoMCInvmass_xSec->GetBinContent(j)-hRecoDataInvmass_xSec->GetBinContent(j)) * (hRecoMCInvmass_xSec->GetBinContent(j)-hRecoDataInvmass_xSec->GetBinContent(j)))/(hRecoMCInvmass_xSec->GetBinError(j)*hRecoDataInvmass_xSec->GetBinError(j));
      //chi2+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinContent(j));
     }
 }
 cout << "Chi2_Mass_smear:" << chi2_mass_smear << endl;
//Unfold
for(int j =0; j <= hGenRap_xSec->GetNbinsX(); j++ ) {
     if(hGenRap_xSec->GetBinContent(j) > 0){
       chi2_rap+= ((hGenRap_xSec->GetBinContent(j)-hUnfoDataRap_xSec->GetBinContent(j)) * (hGenRap_xSec->GetBinContent(j)-hUnfoDataRap_xSec->GetBinContent(j)))/(hGenRap_xSec->GetBinError(j)*hUnfoDataRap_xSec->GetBinError(j));
      //chi2+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinContent(j));
     }
 }
 cout << "Chi2_Rap:" << chi2_rap << endl;
//Smear
for(int j =0; j <= hRecoMCRap_xSec->GetNbinsX(); j++ ) {
     if(hRecoMCRap_xSec->GetBinContent(j) > 0){
       chi2_rap_smear+= ((hRecoMCRap_xSec->GetBinContent(j)-hRecoDataRap_xSec->GetBinContent(j)) * (hRecoMCRap_xSec->GetBinContent(j)-hRecoDataRap_xSec->GetBinContent(j)))/(hRecoMCRap_xSec->GetBinError(j)*hRecoDataRap_xSec->GetBinError(j));
      //chi2+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinContent(j));
     }
 }
 cout << "Chi2_Rap_smear:" << chi2_rap_smear << endl;
////
 double chi2_pt_t, chi2_pt_smear_t;
 for(int j =0; j <= hdiele_GenPt[0]->GetNbinsX()+1; j++ ) {
     if(hdiele_GenPt[0]->GetBinContent(j) > 0){
      chi2_pt_t+= ((hdiele_GenPt[0]->GetBinContent(j)-hUnfoldPt_data->GetBinContent(j)) * (hdiele_GenPt[0]->GetBinContent(j)-hUnfoldPt_data->GetBinContent(j)))/(hdiele_GenPt[0]->GetBinError(j)*hdiele_GenPt[0]->GetBinError(j));
    // chi2_pt_t+= ((hdiele_GenPt[0]->GetBinContent(j)-hUnfoldPt_data->GetBinContent(j)) * (hdiele_GenPt[0]->GetBinContent(j)-hUnfoldPt_data->GetBinContent(j)))/( hdiele_data->GetBinError(j)*hUnfoldPt_data->GetBinError(j));
     // chi2_pt_t+= ((hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)) * (hGenPt_xSec->GetBinContent(j)-hUnfoDataPt_xSec->GetBinContent(j)))/(hGenPt_xSec->GetBinContent(j));
     }
 }
 cout << "Chi2_Pt_t:" << chi2_pt_t << endl;
//smear
 for(int j =0; j <= hdiele_Pt[0]->GetNbinsX()+1; j++ ) {
     if(hdiele_Pt[0]->GetBinContent(j) > 0){
      chi2_pt_smear_t+= ((hdiele_Pt[0]->GetBinContent(j)-hdiele_Pt_data->GetBinContent(j)) * (hdiele_Pt[0]->GetBinContent(j)-hdiele_Pt_data->GetBinContent(j)))/(hdiele_Pt[0]->GetBinError(j)*hdiele_Pt[0]->GetBinError(j));
    // chi2_pt_smear_t+= ((hdiele_Pt[0]->GetBinContent(j)-hdiele_Pt_data->GetBinContent(j)) * (hdiele_Pt[0]->GetBinContent(j)-hdiele_Pt_data->GetBinContent(j)))/(hdiele_Pt_data->GetBinError(j) * hdiele_Pt_data->GetBinError(j));
      //chi2_pt_smear_t+= ((hdiele_Pt[0]->GetBinContent(j)-hdiele_Pt_data->GetBinContent(j)) * (hdiele_Pt[0]->GetBinContent(j)-hdiele_Pt_data->GetBinContent(j)))/(hdiele_Pt[0]->GetBinContent(j));
     }
 }
 cout << "Chi2_pt_smear_t:" << chi2_pt_smear_t << endl;
////
 float p_Pt = TMath::Prob(chi2_pt,hRecoDataPt_xSec->GetNbinsX());
 cout << "p_Pt:" << p_Pt << endl;

/*************************Chi2 From Covarinace matrix*************************************/

TMatrixD J = unfold_dataPt.Eunfold();
    TMatrixD K = responsePt.Mresponse();
     TMatrixD JJt(TMatrixD::kTransposed, J);
    JJt *= J;
    
    TMatrixD Kt(TMatrixD::kTransposed, K);
   // KtK *= K;
    TMatrixD KJJtKt(K, TMatrixD::kMult, JJt);
     KJJtKt *= Kt;
     TDecompSVD svd1(KJJtKt);
    Double_t rank = svd1.Condition();

    std::cout << "Effective rank of KJJtKtK: " << rank << std::endl;
/**************************************************************************************/
TH1D* diff_hist = (TH1D*)hUnfoldPt_data->Clone();
diff_hist->Add(hdiele_GenPt[0], -1);
// TH1D* diff_hist = (TH1D*)hUnfoDataPt_xSec->Clone();
// diff_hist->Add(hGenPt_xSec, -1);
TMatrixD diff_matrix(1, diff_hist->GetNbinsX());
    for (int i = 0; i < diff_hist->GetNbinsX(); ++i) {
        diff_matrix(0, i) = diff_hist->GetBinContent(i + 1);
    }
TMatrixD diff_matrix_T = diff_matrix.Transpose(diff_matrix);

TMatrixD Cov_Pt_Inverse =  covarianceMatrix_Pt.Invert();
double chi_square = 0;
    for (int i = 0; i < diff_matrix.GetNcols(); ++i) {
        for (int j = 0; j < diff_matrix.GetNcols(); ++j) {
            chi_square += diff_matrix_T(0, i) * Cov_Pt_Inverse(i, j) * diff_matrix(0, j);
        }
    }
cout << "Chi_cov_Pt_Unfold:" << chi_square << endl;

/// Smear Chi2

TH1D* diff_hist_smear = (TH1D*)hdiele_Pt_data->Clone();
diff_hist_smear->Add(hdiele_Pt[0], -1);

// TH1D* diff_hist_smear = (TH1D*)hRecoDataPt_xSec->Clone();
// diff_hist_smear->Add(hRecoMCPt_xSec, -1);

TMatrixD diff_matrix_smear(1, diff_hist_smear->GetNbinsX());

    for (int i = 0; i < diff_hist_smear->GetNbinsX(); ++i) {
        diff_matrix_smear(0, i) = diff_hist_smear->GetBinContent(i + 1);
    }
TMatrixD diff_matrix_smear_T = diff_matrix_smear.Transpose(diff_matrix_smear);

TMatrixD res_Pt = responsePt.Mresponse();
TMatrixD res_Pt_Inverse = res_Pt.Invert();
double chi_square_smear_cov = 0;
    for (int i = 0; i < diff_matrix_smear.GetNcols(); ++i) {
        for (int j = 0; j < diff_matrix_smear.GetNcols(); ++j) {
            chi_square_smear_cov += diff_matrix_smear_T(0, i) * res_Pt_Inverse(i, j) * diff_matrix_smear(0, j);
        }
    }
cout << "Chi_cov_Pt_smear:" << chi_square_smear_cov << endl;
/*****************************************************************************************************************************************/
// TMatrixD y_matrix(1, hdiele_Pt_data->GetNbinsX());
//     for (int i = 0; i < hdiele_Pt_data->GetNbinsX(); ++i) {
//         y_matrix(0, i) = hdiele_Pt_data->GetBinContent(i + 1);
//     }

// TMatrixD lambda_matrix(1, hdiele_GenPt[0]->GetNbinsX());
//     for (int i = 0; i < hdiele_GenPt[0]->GetNbinsX(); ++i) {
//         lambda_matrix(0, i) = hdiele_GenPt[0]->GetBinContent(i + 1);
//     }
// TMatrixD  Matrix_Pt_response = responsePt.Mresponse();
//  TMatrixD K_lambda = Matrix_Pt_response * lambda_matrix;
//  TMatrixD y_minus_K_lambda = y_matrix - K_lambda;
//  TMatrixD V_inverse = covarianceMatrix_Pt.Invert();
   
//  TMatrixD y_minus_K_lambda_T = y_minus_K_lambda.Transpose(y_minus_K_lambda);
// // double chi_square_cov_smear;
// // for (int i = 0; i < y_minus_K_lambda.GetNcols(); ++i) {
// //         for (int j = 0; j < y_minus_K_lambda.GetNcols(); ++j) {
// //             chi_square_cov_smear += y_minus_K_lambda_T(0, i) * V_inverse(i, j) * y_minus_K_lambda(0, j);
// //         }
// //     }
// double chi_square_cov_smear = (y_minus_K_lambda_T * V_inverse * y_minus_K_lambda)(0, 0);

// cout << "Chi_cov_Pt_smear:" << chi_square_cov_smear << endl;


// TH1D* diff_hist_xSec = (TH1D*)hUnfoDataPt_xSec->Clone();
// diff_hist_xSec->Add(hGenPt_xSec, -1);
// TMatrixD diff_matrix_xSec(1, diff_hist_xSec->GetNbinsX());
//     for (int i = 0; i < diff_hist_xSec->GetNbinsX(); ++i) {
//         diff_matrix_xSec(0, i) = diff_hist_xSec->GetBinContent(i + 1);
//     }
//  TMatrixD Cov_Pt_Inverse_xSec =  covarianceMatrix_Pt.Invert();
// double chi_square_xSec = 0;
//     for (int i = 0; i < diff_matrix_xSec.GetNcols(); ++i) {
//         for (int j = 0; j < diff_matrix_xSec.GetNcols(); ++j) {
//             chi_square_xSec += diff_matrix_xSec(0, i) * Cov_Pt_Inverse_xSec(i, j) * diff_matrix_xSec(0, j);
//         }
//     }
// cout << "Chi_cov_Pt_xSec:" << chi_square_xSec << endl;

/*****************************************************************************************/
///
outf->cd();

  hGenPt_xSec->Write();
  hGenRap_xSec->Write();
  hGenInvmass_xSec->Write();
  hGenCosthetastar_xSec->Write();

  hRecoMCPt_xSec->Write();
  hRecoMCRap_xSec->Write();
  hRecoMCInvmass_xSec->Write();
  hRecoMCCosthetastar_xSec->Write();

  hRecoDataPt_xSec->Write();
  hRecoDataRap_xSec->Write();
  hRecoDataInvmass_xSec->Write();
  hRecoDataCosthetastar_xSec->Write();


  hUnfoMCPt_xSec->Write();
  hUnfoMCRap_xSec->Write();
  hUnfoMCInvmass_xSec->Write();
  hUnfoMCCosthetastar_xSec->Write();

  hUnfoMCPt_Invert_xSec->Write();
  hUnfoMCRap_Invert_xSec->Write();
  hUnfoMCInvmass_Invert_xSec->Write();
  hUnfoMCCosthetastar_Invert_xSec->Write();
  


  hUnfoDataPt_xSec->Write();
  hUnfoDataRap_xSec->Write();
  hUnfoDataInvmass_xSec->Write();
  hUnfoDataCosthetastar_xSec->Write();

  hUnfoDataInvmass_Invert_xSec->Write();
  hUnfoDataPt_Invert_xSec->Write();
  hUnfoDataRap_Invert_xSec->Write();
  hUnfoDataCosthetastar_Invert_xSec->Write();
  //hfold->Write();
//outf->Write();
 outf->Close();
  TCanvas*c1 =  new TCanvas();
  hfold->SetLineColor(kRed);
  hfold->Draw("hist");
  c1->SaveAs("hfold.pdf");
   
 
  //hUnfold->Write();
  //hUnfoldRapidity->Write();
  //hUnfoldRapidity_data->Write();
  //hUnfoldInvmass->Write();
  
}

TH1D* subtractBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep, const char* name, int nbins, float xmin, float xmax){

  TH1D *hans = new TH1D(name,"",nbins,xmin,xmax);
  //TH1D* hans = (TH1D*)hdata->Clone();
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hdata->GetNbinsX(); i++) {
      double binCont  = hdata->GetBinContent(i)-hqed->GetBinContent(i)-hcep->GetBinContent(i);
      double binError = sqrt(pow(hdata->GetBinError(i),2) + pow(hqed->GetBinError(i),2)  + pow(hcep->GetBinError(i),2) );


      hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError);
      //cout << "******bin:" << i << " " << hdata->GetBinContent(i) << " " << hqed->GetBinContent(i) << " " << hcep->GetBinContent(i) << "   " <<  hans->GetBinContent(i) << "+/-" <<  binError<< endl;
      //double err=0;
      //for (int ivar=1; ivar<14; ivar++) err += pow(hvari[ivar]->GetBinContent(i)-hvari[0]->GetBinContent(i),2);
      //hans->SetBinError(i,sqrt(err+pow(hans->GetBinError(i),2)+pow(glob_syst,2)));
   }

   return hans;
}

TH1D* subtractInvBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep, const char* name) {
    // Clone the data histogram to create the result histogram
    TH1D *hans = (TH1D*)hdata->Clone(name);

    // Add the histograms with unequal binning
    hdata->Add(hqed, -1);
    hdata->Add(hcep, -1);

    // Set the result histogram values based on the modified data histogram
    // 
    
   for (int i = 1; i <= hdata->GetNbinsX(); i++) {
        double binCont  = hdata->GetBinContent(i);
        double binError = hdata->GetBinError(i);

        // If data bin content is zero, set the hans bin content to 0
        if (binCont == 0) {
            hans->SetBinContent(i, 0);
            hans->SetBinError(i, 0);
        } else {
            hans->SetBinContent(i, binCont);
            hans->SetBinError(i, binError);
        }
    }
    // Reset the original data histogram
    hdata->Add(hqed);
    hdata->Add(hcep);

    return hans;
}


TH1D* getXSecHist(TH1D *hist, const char* name, int nbins, float xmin, float xmax, float luminosity){

  TH1D *hans = new TH1D(name,"",nbins,xmin,xmax);
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hist->GetNbinsX(); i++) {

     double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));
     //double binCont  = hist->GetBinContent(i)/(luminosity);
  
      double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
      //double binError = hist->GetBinError(i)/hist->GetBinWidth(i);
      if(binCont == 0){ hans->SetBinContent(i, 0);
       hans->SetBinError(i, 0);}
	else{
      hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError);
     
  
    }

      cout << "******bin:" << i << " " << hans->GetBinContent(i) << "+/-" <<  binError<< endl;
      //double err=0;
      //for (int ivar=1; ivar<14; ivar++) err += pow(hvari[ivar]->GetBinContent(i)-hvari[0]->GetBinContent(i),2);
      //hans->SetBinError(i,sqrt(err+pow(hans->GetBinError(i),2)+pow(glob_syst,2)));
   }
   return hans;

}

TH1D* getInvXSecHist(TH1D *hist, const char* name, float luminosity){

  TH1D *hans = new TH1D(name,"",nMassbin, Massbin);
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hist->GetNbinsX(); i++) {

  //  double binCont  = hist->GetBinContent(i)/(luminosity);
    double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));

     double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
     // double binError = hist->GetBinError(i)/hist->GetBinWidth(i);
      if(binCont == 0){ hans->SetBinContent(i, 0);
       hans->SetBinError(i, 0);}
	else{
      hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError);
    }

      cout << "******bin:" << i << " " << hans->GetBinContent(i) << "+/-" <<  binError<< endl;
      //double err=0;
      //for (int ivar=1; ivar<14; ivar++) err += pow(hvari[ivar]->GetBinContent(i)-hvari[0]->GetBinContent(i),2);
      //hans->SetBinError(i,sqrt(err+pow(hans->GetBinError(i),2)+pow(glob_syst,2)));
   }
   return hans;

}

TH1D* getRapXSecHist(TH1D *hist, const char* name, float luminosity){

  TH1D *hans = new TH1D(name,"",nRapbin, Rapbin);
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hist->GetNbinsX(); i++) {

    //  double binCont  = hist->GetBinContent(i)/(luminosity);
      double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));

      double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
      //double binError = hist->GetBinError(i)/hist->GetBinWidth(i);
      if(binCont == 0){ hans->SetBinContent(i, 0);
       hans->SetBinError(i, 0);}
	else{
      hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError);
    }

      cout << "******bin:" << i << " " << hans->GetBinContent(i) << "+/-" <<  binError<< endl;
      //double err=0;
      //for (int ivar=1; ivar<14; ivar++) err += pow(hvari[ivar]->GetBinContent(i)-hvari[0]->GetBinContent(i),2);
      //hans->SetBinError(i,sqrt(err+pow(hans->GetBinError(i),2)+pow(glob_syst,2)));
   }
   return hans;

}

void printOutput(TH1D *hist){
 for (int i=1; i<=hist->GetNbinsX(); i++) {
      cout << "******bin:" << i << " " << hist->GetBinContent(i) << "+/-" <<  hist->GetBinError(i) << endl; 
    }
}
void make_canvas(TCanvas *& canvas){

  int W = 600;
  int H = 670;
  float T = 0.08;
  float B = 0.14; 
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


void make_hist(TH1D *& hist, Color_t kcolor, int kstyle){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.05, "XYZ");
 
  hist->SetMarkerColor(kcolor);
  hist->SetFillColor(kcolor);
  hist->SetMarkerStyle(kstyle);
 
}

void make_hist_ratio(TH1D *& hist, Color_t kcolor, int kstyle){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
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


float getAcoBinSF(float acop) {
      if (acop<0.002) return 0.990761;
      else if (acop<0.004) return 0.962675;
      else if (acop<0.006) return 0.861166;
      else if (acop<0.008) return 0.598683;
      else if (acop<0.01) return 0.264046;    
      else return 1;
   
}

Double_t GetDeltaPhi(double dphi, double phi){

     Double_t deltaPhi = dphi - phi;

  if ( deltaPhi > 3.141592653589 )
    deltaPhi = deltaPhi - 2. * 3.141592653589;
  if ( deltaPhi <= -3.141592653589 )
    deltaPhi = deltaPhi + 2. * 3.141592653589;

  if ( TMath::Abs(deltaPhi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
  }

  return TMath::Abs(deltaPhi);
  //return dphi;


}


Double_t RecoID_SF_Electron( double ETT){

    double weight(-1);
    double absETT = abs(ETT);

    if (absETT > 2.4) {
        std::cout << "[WARNING] Electron pseudo-rapidity (" << ETT << ") outside [-2.4, 2.4]" << std::endl;
        return 1;
    }

    if (ETT > -2.2 && ETT < -2) {
            weight = 0.873564;
    }
    else if (ETT >= -2 && ETT < -1.5)
            weight = 0.865897;
    else if (ETT >= -1.5 && ETT < -1)
            weight = 0.956456;
    else if (ETT >= -1 && ETT < -0.5)
            weight = 0.915635;
    else if (ETT >= -0.5 && ETT < 0)
            weight = 0.926652;

    else if (ETT >= 0 && ETT < 0.5)
           weight = 0.960205;
    else if (ETT >= 0.5 && ETT < 1.0)
            weight =  0.948991;
    else if (ETT >= 1.0 && ETT < 1.5)
            weight = 0.978441;
    else if (ETT >= 1.5 && ETT < 2.0)
            weight = 0.975938;
    else if (ETT >= 2.0 && ETT < 2.2)
            weight = 0.859258;

 return weight;
}
