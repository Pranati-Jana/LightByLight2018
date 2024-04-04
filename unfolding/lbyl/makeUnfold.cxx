//this is a script to make response matrix from LbyL MC 
// Background subtracted from data
// Data unfolded
// All numbers (cross-section, luminosity for data to get cross-section is in nano barn
// Edited from the script created by Ruchi Chudasama
///////////////////////////////////////////////////////////////////////////////////////
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
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

#include "ReadLbyLTree.C"

const int nSample = 1;
//const char *Sample[nSample]={"Data","QED_SC","QED_SL"};
const char *Sample[nSample]={"LbL"};

void make_canvas(TCanvas *&);
void make_canvas_ratio(TCanvas *&);
void make_hist(TH1D *&, Color_t , int );
void make_hist_ratio(TH1D *&,Color_t, int) ;
double_t GetDeltaPhi(double dphi, double phi);
void printOutput(TH1D *hist);
double costhetastar_CS(TLorentzVector p1, TLorentzVector pair); 
//const double luminosity       = 1635.123139823; // μb^-1
const double luminosity       = 1;
//const double luminosity       = 1.6722;
//const double luminosity       = 1639.207543; // μb^-1
//const double luminosity       = 1.639207543; // μb^-1
double scaleFactorPhoton_Old = 0.85 * 0.93 * 0.866 * 1.008 * 1.037 * 1.037;

double scaleFactorPhoton_old1 = 0.8477 *  //0.85// NEE    21.12.2021
  0.9322 *      //0.93// CHE  21.12.2021
  pow(0.9771, 2)* //1.048// Photon  reco+ID 21.12.2021
  0.8643 *      //0.866// HF veto
  1.0006;       //1.008// L1 EG trigger
double scaleFactorPhoton  = 0.9758 * 0.9758 *0.946 * 0.946 * 0.8487 * 0.9252* 1.0089 * 0.8716;
//double scaleFactorPhoton  = 0.9758*0.946*1.0089*0.8716*0.9252*0.8487;
const double xsecGeneratedLbLSC    = 2590; // nb Superchic
const double nEventsGeneratedLbLSC = 466000;  //Superchic
double norm_LbLSC = scaleFactorPhoton*xsecGeneratedLbLSC*luminosity/nEventsGeneratedLbLSC;  //31.12.21 (SF = 1.048, took the square here). 
//double norm_LbLSC = 1;
double W1 = xsecGeneratedLbLSC*luminosity/nEventsGeneratedLbLSC; 
//double W1 = 1;
//double W1 = scaleFactorPhoton*xsecGeneratedLbLSC*luminosity/nEventsGeneratedLbLSC;

const double xsecGeneratedLbLMG    = 140.6; // nb Madgraph David
const double nEventsGeneratedLbLMG = 788069; //Madgraph David
double norm_LbLMG = scaleFactorPhoton*xsecGeneratedLbLMG*luminosity/nEventsGeneratedLbLMG;  //31.12.21 (SF = 1.048, took the square here). 



double wt[nSample] = {norm_LbLSC};
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


const int nRapbins=5;
//double Rapbin[nRapbins]={-2.5,-0.75,0.75,2.5};
//double Rapbin[nRapbins]={0.0,0.5,1.0,1.5,2.2};
double Rapbin[nRapbins]={-2.2, -1.1, 0.0, 1.1, 2.2};
const int nRapbin= sizeof(Rapbin)/sizeof(double) - 1;

const int nMassbins=6;
//double Massbin[nMassbins]={5.0,7.0,10.0,30.0};
//double Massbin[nMassbins]={5.0,7.5,10.5,16.0};
//double Massbin[nMassbins]={5.0,7.3,13.0,16.0};
//double Massbin[nMassbins]={5.0,7.0,10.0,16.0};
double Massbin[nMassbins]={5.0,9.0,13.0,17.0, 21.0,25.0};
const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;


TH1D* subtractBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep,TH1D *hqed_sl, const char* name, int nbins, float xmin, float xmax);
TH1D* subtractInvBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep,TH1D *hqed_sl, const char* name);
TH1D* getXSecHist(TH1D *hist, const char* name, int nbins, float xmin, float xmax, float lumi);
TH1D* getInvXSecHist(TH1D *hist, const char* name, float lumi);
TH1D* getRapXSecHist(TH1D *hist, const char* name, float lumi);

TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, TH1D*, double , double , double , double , double , double , const char *, bool); 
void makeUnfold(int method=1){
  int itr = 1; 
  int cosbin = 4;
  int rapbin = 3;
  //gROOT->LoadMacro("CMS_lumi.C");
  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  
  char MethID1[100]; 
  if(method==1){
        
      //sprintf(MethID1,"Bayes_unfo");
    } 
  if(method==2){
    
      //sprintf(MethID1,"Svd_unfo ");
    }
  if(method==3){
    
      //sprintf(MethID1,"BinByBin_unfo");
    }

  string mc_name = "SC";


  bool absRap = false;

  TFile *outf= new TFile(Form("unfolding_histograms_Bayes_%s.root",mc_name.c_str()),"recreate");  
  
  TChain *lbyl[nSample];
  
  TH1D* hdipho_Pt[nSample], *hdipho_Rapidity[nSample], *hdipho_Invmass[nSample], *hAcoplanarity[nSample], *hdipho_Costhetastar[nSample]; TH1F *hSF[nSample]; 
  
  TH1D* hdipho_GenPt[nSample], *hdipho_GenRapidity[nSample], *hdipho_GenInvmass[nSample], *hdipho_GenCosthetastar[nSample];
  TH2D* hdipho_RecoGenPt[nSample], *hdipho_RecoGenRapidity[nSample], *hRecoGenInvmass[nSample], *hdipho_RecoGenCosthetastar[nSample];
    

  //RooUnfoldResponse responseRapidity (hdipho_Rapidity, hdipho_GenRapidity);
  //RooUnfoldResponse responseInvmass (hdipho_Invmass, hdipho_GenInvmass);
  
  cout << "define hists"<< endl;
  
  int i =0;
  //for (int i = 0; i < nSample; i++){
    lbyl[i] = new TChain("output_tree");
//    lbyl[i]->Add(Form("RootFile/lbylSC_diphoton_genInfo.root"));
     cout << "test:" << endl;
   //lbyl[i]->Add(Form("RootFile/14thFeb_LbyL_Unfolding.root"));
   lbyl[i]->Add(Form("RootFile/SC_LbL_Costhetastar.root"));
      
    hdipho_Pt[i]    = new TH1D(Form("hdipho_Pt%s", Sample[i]),"",5,0,1);   
   // hdipho_Rapidity[i]   = new TH1D(Form("hdipho_Rapidity%s",Sample[i]),"",nRapbin, Rapbin);
    hdipho_Rapidity[i]   = new TH1D(Form("hdipho_Rapidity%s",Sample[i]),"",rapbin,-2.2,2.2);
    cout << "test:" << endl; 
    hdipho_Invmass[i]   = new TH1D(Form("hdipho_Invmass%s",Sample[i]),"",nMassbin, Massbin);
   // hdipho_Invmass[i]   = new TH1D(Form("hdipho_Invmass%s",Sample[i]),"",5,5,25);
    hdipho_Costhetastar[i] = new TH1D(Form("hdipho_Costhetastar%s",Sample[i]),"",cosbin,0,1);

    
    hdipho_GenPt[i]         = new TH1D(Form("hdipho_GenPt%s", Sample[i]),"",5,0,1);
   // hdipho_GenRapidity[i]   = new TH1D(Form("hdipho_GenRapidity%s",Sample[i]),"",nRapbin, Rapbin);
    hdipho_GenRapidity[i]   = new TH1D(Form("hdipho_GenRapidity%s",Sample[i]),"",rapbin,-2.2,2.2);
    hdipho_GenInvmass[i]    = new TH1D(Form("hdipho_GenInvmass%s",Sample[i]),"",nMassbin, Massbin);
    hdipho_GenCosthetastar[i] = new TH1D(Form("hdipho_GenCosthetastar%s",Sample[i]),"",cosbin,0,1);
    
    hdipho_RecoGenPt[i]         = new TH2D(Form("hdipho_RecoGenPt%s", Sample[i]),"",5,0,1,5,0,1);
    //hdipho_RecoGenRapidity[i]   = new TH2D(Form("hdipho_RecoGenRapidity%s",Sample[i]),"",nRapbin, Rapbin,nRapbin, Rapbin);
    hdipho_RecoGenRapidity[i]   = new TH2D(Form("hdipho_RecoGenRapidity%s",Sample[i]),"",rapbin,-2.2,2.2,rapbin,-2.2,2.2);
    hRecoGenInvmass[i]          = new TH2D(Form("hRecoGenInvmass%s",Sample[i]),"",nMassbin, Massbin,nMassbin, Massbin);
    hdipho_RecoGenCosthetastar[i]  = new TH2D(Form("hdipho_RecoGenCosthetastar%s",Sample[i]),"",cosbin,0,1,cosbin,0,1);

  // ===================   Cross section histograms   ======================     

  TH1D* hGenPt_xSec              = new TH1D("hGenPt_xSec", "Pt gen MC",5,0,1);
  //TH1D* hGenRap_xSec             = new TH1D("hGenRap_xSec","Rap gen MC",nRapbin, Rapbin);
  TH1D* hGenRap_xSec             = new TH1D("hGenRap_xSec","Rap gen MC",rapbin,-2.2,2.2);
  TH1D* hGenInvmass_xSec         = new TH1D("hGenInvmass_xSec","Invmass gen MC",nMassbin, Massbin);
  TH1D* hGenCosthetastar_xSec  = new TH1D("GenCosthetastar_xSec","Costhetastar gen MC",cosbin,0,1);

  TH1D* hRecoMCPt_xSec              = new TH1D("hRecoMCPt_xSec", "Pt reco MC",5,0,1);
  //TH1D* hRecoMCRap_xSec             = new TH1D("hRecoMCRap_xSec","Rap reco MC",nRapbin, Rapbin);
  TH1D* hRecoMCRap_xSec             = new TH1D("hRecoMCRap_xSec","Rap gen MC",rapbin,-2.2,2.2);
  TH1D* hRecoMCInvmass_xSec         = new TH1D("hRecoMCInvmass_xSec","Invmass reco MC",nMassbin, Massbin);
  TH1D* hRecoMCCosthetastar_xSec    = new TH1D("hRecoMCCosthetastar_xSec","Costhetastar reco MC",cosbin,0,1);

  TH1D* hRecoDataPt_xSec              = new TH1D("hRecoDataPt_xSec", "Pt reco Data",5,0,1);
  //TH1D* hRecoDataRap_xSec             = new TH1D("hRecoDataRap_xSec","Rap reco Data",nRapbin, Rapbin);
  TH1D* hRecoDataRap_xSec             = new TH1D("hRecoDataRap_xSec","Rap gen MC",rapbin,-2.2,2.2);
  TH1D* hRecoDataInvmass_xSec         = new TH1D("hRecoDataInvmass_xSec","Invmass reco Data",nMassbin, Massbin);
  TH1D* hRecoDataCosthetastar_xSec    = new TH1D("hRecoDataCosthetastar_xSec","Costhetastar reco Data",cosbin,0,1);


  TH1D* hUnfoMCPt_xSec              = new TH1D("hUnfoMCPt_xSec", "Pt unfo MC",5,0,1);
  //TH1D* hUnfoMCRap_xSec             = new TH1D("hUnfoMCRap_xSec","Rap unfo MC",nRapbin, Rapbin);
  TH1D* hUnfoMCRap_xSec             = new TH1D("hUnfoMCRap_xSec","Rap gen MC",rapbin,-2.2,2.2);
  TH1D* hUnfoMCInvmass_xSec         = new TH1D("hUnfoMCInvmass_xSec","Invmass unfo MC",nMassbin, Massbin);
  TH1D* hUnfoMCCosthetastar_xSec    = new TH1D("hUnfoMCCosthetastar_xSec","Costhetastar unfo MC",cosbin,0,1);

  TH1D* hUnfoMCPt_Invert_xSec              = new TH1D("hUnfoMCPt_Invert_xSec", "Pt unfo Invert MC",5,0,1);
//  TH1D* hUnfoMCRap_Invert_xSec             = new TH1D("hUnfoMCRap_Invert_xSec","Rap unfo Invert MC",nRapbin, Rapbin);
  TH1D* hUnfoMCRap_Invert_xSec             = new TH1D("hUnfoMCRap_Invert_xSec","Rap unfo Invert MC",rapbin,-2.2, 2.2);
  TH1D* hUnfoMCInvmass_Invert_xSec         = new TH1D("hUnfoMCInvmass_Invert_xSec","Invmass unfo Invert MC",nMassbin, Massbin);
  TH1D* hUnfoMCCosthetastar_Invert_xSec    = new TH1D("hUnfoMCCosthetastar_Invert_xSec","Costhetastar unfo Invert MC",cosbin,0,1);

 TH1D* hUnfoDataPt_xSec              = new TH1D("hUnfoDataPt_xSec", "Pt unfo Data",5,0,1);
  //TH1D* hUnfoDataRap_xSec             = new TH1D("hUnfoDataRap_xSec","Rap unfo Data",nRapbin, Rapbin);
 TH1D* hUnfoDataRap_xSec             = new TH1D("hUnfoDataRap_xSec","Rap gen MC",rapbin,-2.2,2.2);
 TH1D* hUnfoDataInvmass_xSec         = new TH1D("hUnfoDataInvmass_xSec","Invmass unfo Data",nMassbin, Massbin);
 TH1D* hUnfoDataCosthetastar_xSec    = new TH1D("hUnfoDataCosthetastar_xSec","Costhetastar unfo Data",cosbin,0,1);
 
  TH1D* hUnfoDataPt_Invert_xSec              = new TH1D("hUnfoDataPt_Invert_xSec", "Pt unfo Invert Data",5,0,1);
//TH1D* hUnfoDataRap_Invert_xSec             = new TH1D("hUnfoDataRap_Invert_xSec","Rap unfo Invert Data",nRapbin, Rapbin);
  TH1D* hUnfoDataRap_Invert_xSec             = new TH1D("hUnfoDataRap_Invert_xSec","Rap unfo Invert Data",rapbin,-2.2,2.2);
  TH1D* hUnfoDataInvmass_Invert_xSec         = new TH1D("hUnfoDataInvmass_Invert_xSec","Invmass unfo Invert Data",nMassbin, Massbin);
  TH1D* hUnfoDataCosthetastar_Invert_xSec    = new TH1D("hUnfoDataCosthetastar_Invert_xSec","Costhetastar unfo Invert Data",cosbin,0,1);

   //RooUnfoldResponse responsePt (hdipho_Pt[i], hdipho_GenPt[i]);
//   RooUnfoldResponse responsePt (5,0,1,5,0,1);
  // RooUnfoldResponse responseRapidity (hdipho_Rapidity[i], hdipho_GenRapidity[i]);
//    RooUnfoldResponse responseInvmass (hdipho_Invmass[i], hdipho_GenInvmass[i]);
   
    RooUnfoldResponse responsePt (hdipho_Pt[i], hdipho_GenPt[i],hdipho_RecoGenPt[i]);
    RooUnfoldResponse responseRapidity (hdipho_Rapidity[i], hdipho_GenRapidity[i],hdipho_RecoGenRapidity[i]);
    RooUnfoldResponse responseInvmass (hdipho_Invmass[i], hdipho_GenInvmass[i],hRecoGenInvmass[i]);
    RooUnfoldResponse responseCosthetastar (hdipho_Costhetastar[i], hdipho_GenCosthetastar[i],hdipho_RecoGenCosthetastar[i]);
    cout << "file " << lbyl[i]->GetEntries()  << endl;

    ReadLbyLTree  lbylR(lbyl[i]);
    lbylR.fChain->SetBranchStatus("*",1);

    if (lbylR.fChain == 0) return;
    
    Long64_t nentries = lbylR.fChain->GetEntriesFast();
    cout << lbyl[i]->GetName() << "    " << nentries << endl;
    
    Long64_t nbytes = 0, nb = 0;
    //for (Long64_t jentry=0; jentry<10;jentry++) {
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry_evt = lbylR.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = lbylR.fChain->GetEntry(jentry);   nbytes += nb;
     

      // "===================      Fill gen Info here ====================== 
      
      if(abs(lbylR.gen_Eta1) > 2.2 || abs(lbylR.gen_Eta2) > 2.2 || lbylR.gen_Pt1 < 2.0 || lbylR.gen_Pt2 < 2.0 || lbylR.gen_diPho_M < 5 || lbylR.gen_diPho_Pt > 1) continue; 
      //  if(abs(lbylR.gen_diPho_Rapidity) > 2.2) continue;
      double ele_dphi = GetDeltaPhi(lbylR.gen_Phi1, lbylR.gen_Phi2);
      double ele_gen_acop = 1 - (GetDeltaPhi(lbylR.gen_Phi1, lbylR.gen_Phi2)/3.141592653589);
      if(ele_gen_acop > 0.01) continue;
      

  
      hdipho_GenPt[i]->Fill(lbylR.gen_diPho_Pt,W1);
      if(absRap)hdipho_GenRapidity[i]->Fill(abs(lbylR.gen_diPho_Rapidity),W1);
      if(!absRap)hdipho_GenRapidity[i]->Fill(lbylR.gen_diPho_Rapidity,W1);
      hdipho_GenInvmass[i]->Fill(lbylR.gen_diPho_M,W1);
      hdipho_GenCosthetastar[i]->Fill(abs(lbylR.gen_Costhetastar),W1);
 
      // ===================   Apply analysis cuts on reco MC ======================           
  //    bool isPassing = 1;
      bool isPassing = lbylR.phoSwissCross_1 < 0.95 && lbylR.phoSwissCross_2 <  0.95
		        && lbylR.ok_neuexcl == 1 && lbylR.ok_chexcl_goodtracks == 1 && lbylR.ok_chexcl_goodelectrons == 1 
		        && lbylR.vSum_M > 5  && lbylR.vSum_Pt < 1 
 			&& abs(lbylR.phoEta_1) < 2.2 && abs(lbylR.phoEta_2) < 2.2 && lbylR.phoEt_1 > 2.0 && lbylR.phoEt_2 > 2.0		
                        && lbylR.ok_chexcl_muons ==1  
		 	&& lbylR.pho_acop < 0.01 ;


  if(lbylR.ok_trigger == 1) { //trigger 
	if(isPassing){

	hdipho_Pt[i]->Fill(lbylR.vSum_Pt,wt[i]);
	if(absRap)hdipho_Rapidity[i]->Fill(abs(lbylR.vSum_Rapidity),wt[i]);
	if(!absRap)hdipho_Rapidity[i]->Fill(lbylR.vSum_Rapidity,wt[i]);
	hdipho_Invmass[i]->Fill(lbylR.vSum_M,wt[i]);
  hdipho_Costhetastar[i]->Fill(abs(lbylR.Costhetastar),wt[i]);
	
	hdipho_RecoGenPt[i]->Fill(lbylR.vSum_Pt,lbylR.gen_diPho_Pt, wt[i]);
	if(absRap)hdipho_RecoGenRapidity[i]->Fill(abs(lbylR.vSum_Rapidity),abs(lbylR.gen_diPho_Rapidity),wt[i]);
	if(!absRap)hdipho_RecoGenRapidity[i]->Fill(lbylR.vSum_Rapidity,lbylR.gen_diPho_Rapidity,wt[i]);
	hRecoGenInvmass[i]->Fill(lbylR.vSum_M,lbylR.gen_diPho_M,wt[i]);
  hdipho_RecoGenCosthetastar[i]->Fill(abs(lbylR.Costhetastar),abs(lbylR.gen_Costhetastar),wt[i]);
        
  responsePt.Fill(lbylR.vSum_Pt,lbylR.gen_diPho_Pt,wt[i]);
  if(absRap)responseRapidity.Fill(abs(lbylR.vSum_Rapidity),abs(lbylR.gen_diPho_Rapidity),wt[i]);
	if(!absRap)responseRapidity.Fill(lbylR.vSum_Rapidity,lbylR.gen_diPho_Rapidity,wt[i]);
	responseInvmass.Fill(lbylR.vSum_M,lbylR.gen_diPho_M,wt[i]);
  responseCosthetastar.Fill(abs(lbylR.Costhetastar),abs(lbylR.gen_Costhetastar),wt[i]);
	} //isPassing analysis cuts
	else {
	responsePt.Miss(lbylR.gen_diPho_Pt,W1);
	if(absRap)responseRapidity.Miss(abs(lbylR.gen_diPho_Rapidity),W1);
	if(!absRap)responseRapidity.Miss(lbylR.gen_diPho_Rapidity,W1);
	responseInvmass.Miss(lbylR.gen_diPho_M,W1);
  responseCosthetastar.Miss(abs(lbylR.gen_Costhetastar),W1);
	}
      } //trigger
	else {
	responsePt.Miss(lbylR.gen_diPho_Pt, W1);
	if(absRap)responseRapidity.Miss(abs(lbylR.gen_diPho_Rapidity),W1);
	if(!absRap)responseRapidity.Miss(lbylR.gen_diPho_Rapidity,W1);
	responseInvmass.Miss(lbylR.gen_diPho_M,W1);
  responseCosthetastar.Miss(abs(lbylR.gen_Costhetastar),W1);
	}
 
    } //entry
 // } // for 4 files

 TCanvas*C =  new TCanvas();
 hdipho_RecoGenRapidity[0]->Draw("Colz");
 //C->Print("fig_Bayes_Inverse/Invmass.pdf");
cout << "gen entries pass:" << hdipho_GenPt[0]->GetEntries() << endl; 
cout << "reco entries pass:" << hdipho_Pt[0]->GetEntries() << endl; 
/*******************************Condition*number*********************************/
cout << "Inverse1" << endl;
TDecompSVD *svd= new TDecompSVD (responseInvmass.Mresponse());
auto singular_values = svd->GetSig();
//cout << "Test1" << endl;
svd->Print();
double singular_value_min;
      for (int i = 0; i < 4 ; i++)
      {
    //    cout << "Test2" << endl;
        cout << "Condition number:" << singular_values[i] << endl;
       if ( singular_values[i] > pow(10,-15)  ) singular_value_min = singular_values[i]; // the pow(10,-15) > requirement is to suppress singular values from empty bins 
      }

cout << "condition number: "  << singular_values[0]/singular_value_min << endl;
cout << "Inverse2" << endl;
/*********************************************************************************************/
  cout<< "======================================Response matrix========================="<<endl;
  char MethID2[100]; 
  //if(method==1)
   // {

      RooUnfoldBayes unfold(&responsePt, hdipho_Pt[0], itr);
      RooUnfoldInvert unfold2(&responseRapidity, hdipho_Rapidity[0]);
      RooUnfoldInvert unfold3(&responseInvmass, hdipho_Invmass[0]);
      RooUnfoldInvert unfold4(&responseCosthetastar, hdipho_Costhetastar[0]);

      RooUnfoldInvert unfold_Invert(&responsePt, hdipho_Pt[0]);
      RooUnfoldInvert unfold2_Invert(&responseRapidity, hdipho_Rapidity[0]);
      RooUnfoldInvert unfold3_Invert(&responseInvmass, hdipho_Invmass[0]);
      RooUnfoldInvert unfold4_Invert(&responseCosthetastar, hdipho_Costhetastar[0]);


     //   RooUnfoldSvd unfold(&responsePt, hdipho_Pt[0], 4);
       // RooUnfoldSvd unfold2(&responseRapidity, hdipho_Rapidity[0], 4);
      //  RooUnfoldSvd unfold3(&responseInvmass, hdipho_Invmass[0], 4);
/********************Matrix-plotting************/
        TCanvas*c4 =  new TCanvas("c4", "" ,800,600);
        c4->SetRightMargin(0.15);
        TMatrixD  responseMatrix1 = responsePt.Mresponse();
        TH2D *hResponseMatrix = new TH2D("hResponseMatrix", "Response Matrix;Reco diphoton p_{T} bin index;Gen diphoton p_{T} bin index", responseMatrix1.GetNrows(), 0, responseMatrix1.GetNrows(), responseMatrix1.GetNcols(), 0, responseMatrix1.GetNcols());
        for (int i = 0; i < responseMatrix1.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix1.GetNcols(); ++j) {
                  hResponseMatrix->SetBinContent(i + 1, j + 1, responseMatrix1[i][j]);
                  }
        }
        hResponseMatrix->GetXaxis()->SetNdivisions(5);
        hResponseMatrix->GetYaxis()->SetNdivisions(5);
        hResponseMatrix->Draw("textcolz");
        c4->Update();
        c4->SaveAs("fig_Bayes_Inverse/Pt_lbl_ZDCAND3n.pdf");
        //
        TCanvas*c5 =  new TCanvas("c5","Response Matrix2",800,600);
        c5->SetRightMargin(0.15);
        TMatrixD  responseMatrix2 = responseRapidity.Mresponse();
        TH2D *hResponseMatrix2 = new TH2D("hResponseMatrix2", "Response Matrix2;Reco rapidity bin index;Gen rapidity bin index", responseMatrix2.GetNrows(), 0, responseMatrix2.GetNrows(), responseMatrix2.GetNcols(), 0, responseMatrix2.GetNcols());
        for (int i = 0; i < responseMatrix2.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix2.GetNcols(); ++j) {
                  hResponseMatrix2->SetBinContent(i + 1, j + 1, responseMatrix2[i][j]);
                  }
        }
        hResponseMatrix2->GetXaxis()->SetNdivisions(4);
        hResponseMatrix2->GetYaxis()->SetNdivisions(4);
        hResponseMatrix2->Draw("textcolz");
        c5->Update();
        c5->SaveAs("fig_Bayes_Inverse/Rapidity_lbl_ZDCAND3n.pdf");
        //
        TCanvas*c6 =  new TCanvas("c6","Response Matrix3", 800,600);
        c6->SetRightMargin(0.15);
        TMatrixD  responseMatrix3 = responseInvmass.Mresponse();
        TH2D *hResponseMatrix3 = new TH2D("hResponseMatrix3", "Response Matrix3;Reco invariant mass bin index;Gen invariant mass bin index", responseMatrix3.GetNrows(), 0, responseMatrix3.GetNrows(), responseMatrix3.GetNcols(), 0, responseMatrix3.GetNcols());
        for (int i = 0; i < responseMatrix3.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix3.GetNcols(); ++j) {
                  hResponseMatrix3->SetBinContent(i + 1, j + 1, responseMatrix3[i][j]);
                  }
        }
        hResponseMatrix3->Draw("textcolz");
        c6->Update();
        c6->SaveAs("fig_Bayes_Inverse/InvMass_lbl_ZDCAND3n.pdf");
      //
      TCanvas*c7 =  new TCanvas("c7","Response Matrix4", 800,600);
        c7->SetRightMargin(0.15);
        TMatrixD  responseMatrix4 = responseCosthetastar.Mresponse();
        TH2D *hResponseMatrix4 = new TH2D("hResponseMatrix4", "Response Matrix4;Reco Costhetastar bin index;Gen costhetastar bin index", responseMatrix4.GetNrows(), 0, responseMatrix4.GetNrows(), responseMatrix4.GetNcols(), 0, responseMatrix4.GetNcols());
        for (int i = 0; i < responseMatrix4.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix4.GetNcols(); ++j) {
                  hResponseMatrix4->SetBinContent(i + 1, j + 1, responseMatrix4[i][j]);
                  }
        }
        hResponseMatrix4->Draw("textcolz");
        c7->Update();
        c7->SaveAs("fig_Bayes_Inverse/Costhetatstar_lbl_ZDCAND3n.pdf");  
        
/***********************************************/

//  unfold.IncludeSystematics(2);
  TH1D* hUnfoldPt= (TH1D*) unfold.Hunfold();
  unfold.PrintTable (cout,hdipho_GenPt[0]);

  unfold2.IncludeSystematics(2); 
  unfold3.IncludeSystematics(2);
  unfold2_Invert.IncludeSystematics(2);
  unfold3_Invert.IncludeSystematics(2);
  
  TH1D* hUnfoldRapidity= (TH1D*) unfold2.Hunfold();
  unfold2.PrintTable (cout,hdipho_GenRapidity[0]);

  TH1D* hUnfoldInvmass= (TH1D*) unfold3.Hunfold();
  unfold3.PrintTable (cout,hdipho_GenInvmass[0]);

  TH1D* hUnfoldCosthetastar= (TH1D*) unfold4.Hunfold();
  unfold4.PrintTable (cout,hdipho_GenCosthetastar[0]);
  
  TH1D* hUnfoldPt_Invert= (TH1D*) unfold_Invert.Hunfold();
  unfold_Invert.PrintTable (cout,hdipho_GenPt[0]);

  TH1D* hUnfoldRapidity_Invert= (TH1D*) unfold2_Invert.Hunfold();
  unfold2_Invert.PrintTable (cout,hdipho_GenRapidity[0]);

  TH1D* hUnfoldInvmass_Invert= (TH1D*) unfold3_Invert.Hunfold();
  unfold3_Invert.PrintTable (cout,hdipho_GenInvmass[0]);

  TH1D* hUnfoldCosthetastar_Invert= (TH1D*) unfold4_Invert.Hunfold();
  unfold4_Invert.PrintTable (cout,hdipho_GenCosthetastar[0]);

  /***************************  Start reading data and subtract background contributions ********************************/

  //TFile *f1=new TFile("diphoton_histos.root","r"); 
  TFile *f1=new TFile("RootFile/unfolding_histograms/unfoldingHistograms_collisionData.root","r"); 
  TFile *f2=new TFile("RootFile/unfolding_histograms/unfoldingHistograms_qed.root","r"); 
  TFile *f3=new TFile("RootFile/unfolding_histograms/unfoldingHistograms_cep.root","r"); 
  TFile *f4=new TFile("RootFile/unfolding_histograms/unfoldingHistograms_qed_starlight.root","r"); 

cout << "Test1:" << endl;
  TH1D *hdipho_Pt_data        = (TH1D*)f1->Get("unfoldingPhoton_pt");
  //TH1D *hdipho_Rapidity_data  = (TH1D*)f1->Get("unfoldingPhoton_absRap"); 
  TH1D *hdipho_Rapidity_data  = (TH1D*)f1->Get("unfoldingPhoton_rap3"); 
  TH1D *hdipho_Invmass_data   = (TH1D*)f1->Get("unfoldingPhoton_mass"); 
  TH1D *hdipho_Costhetastar_data   = (TH1D*)f1->Get("unfoldingPhoton_costhetastar4");

  TH1D *hdipho_Pt_qed        = (TH1D*)f2->Get("unfoldingPhoton_pt");
  //TH1D *hdipho_Rapidity_qed  = (TH1D*)f2->Get("unfoldingPhoton_absRap"); 
  TH1D *hdipho_Rapidity_qed  = (TH1D*)f2->Get("unfoldingPhoton_rap3"); 
  TH1D *hdipho_Invmass_qed   = (TH1D*)f2->Get("unfoldingPhoton_mass");
  TH1D *hdipho_Costhetastar_qed   = (TH1D*)f2->Get("unfoldingPhoton_costhetastar4");
  // hdipho_Invmass_qed->SetBinContent(4,0);
  // hdipho_Invmass_qed->SetBinContent(5,0);
  
  TH1D *hdipho_Pt_cep        = (TH1D*)f3->Get("unfoldingPhoton_pt");
//  TH1D *hdipho_Rapidity_cep  = (TH1D*)f3->Get("unfoldingPhoton_absRap"); 
  TH1D *hdipho_Rapidity_cep  = (TH1D*)f3->Get("unfoldingPhoton_rap3"); 
  TH1D *hdipho_Invmass_cep   = (TH1D*)f3->Get("unfoldingPhoton_mass"); 
  TH1D *hdipho_Costhetastar_cep   = (TH1D*)f3->Get("unfoldingPhoton_costhetastar4");
  // hdipho_Invmass_cep->SetBinContent(4,0);
  // hdipho_Invmass_cep->SetBinContent(5,0);
  TH1D *hdipho_Pt_qed_sl        = (TH1D*)f4->Get("unfoldingPhoton_pt"); 
 // TH1D *hdipho_Rapidity_qed_sl  = (TH1D*)f4->Get("unfoldingPhoton_absRap");
  TH1D *hdipho_Rapidity_qed_sl  = (TH1D*)f4->Get("unfoldingPhoton_rap3"); 
  TH1D *hdipho_Invmass_qed_sl   = (TH1D*)f4->Get("unfoldingPhoton_mass");
  TH1D *hdipho_Costhetastar_qed_sl   = (TH1D*)f4->Get("unfoldingPhoton_costhetastar4");
  // hdipho_Invmass_qed_sl->SetBinContent(4,0);
  // hdipho_Invmass_qed_sl->SetBinContent(5,0);

  // "===================      Subtract background ======================     
  TH1D* hPt_data_bkgSub       = subtractBkg(hdipho_Pt_data,      hdipho_Pt_qed,      hdipho_Pt_cep, hdipho_Pt_qed_sl , "hPt_data_bkgSub", 5,0,1);
  TH1D* hRapidity_data_bkgSub = subtractBkg(hdipho_Rapidity_data,hdipho_Rapidity_qed,hdipho_Rapidity_cep, hdipho_Rapidity_qed_sl, "hRapidity_data_bkgSub", rapbin,-2.2,2.2);
  TH1D* hInvmass_data_bkgSub  = subtractInvBkg(hdipho_Invmass_data,hdipho_Invmass_qed,  hdipho_Invmass_cep, hdipho_Invmass_qed_sl , "hInvmass_data_bkgSub");
 //  TH1D* hInvmass_data_bkgSub  = subtractBkg(hdipho_Invmass_data,hdipho_Invmass_qed,  hdipho_Invmass_cep, hdipho_Invmass_qed_sl , "hInvmass_data_bkgSub",5,5,25);
  TH1D* hCosthetastar_data_bkgSub = subtractBkg(hdipho_Costhetastar_data,hdipho_Costhetastar_qed,hdipho_Costhetastar_cep, hdipho_Costhetastar_qed_sl, "hCosthetastar_data_bkgSub", cosbin,0,1);
  
  // hdipho_Rapidity_qed->Add(hdipho_Rapidity_cep,1);
  // TH1D* hRapidity_data_bkgSub = (TH1D*)hdipho_Rapidity_data->Clone();
  // hRapidity_data_bkgSub->Add(hdipho_Rapidity_qed,-1);

 //For ()

  // hdipho_Invmass_qed->Add(hdipho_Invmass_cep,1);
  // TH1D* hInvmass_data_bkgSub = (TH1D*)hdipho_Invmass_data->Clone();
  // hInvmass_data_bkgSub->Add(hdipho_Invmass_qed,-1); 
////////////////////////////

/*
  hdipho_Rapidity_qed->Add(hdipho_Rapidity_cep,1);
  TH1D* hRapidity_data_bkgSub = (TH1D*)hdipho_Rapidity_data->Clone();
//  hRapidity_data_bkgSub->Add(hdipho_Rapidity_qed,-1);
  hRapidity_data_bkgSub->Add(hdipho_Rapidity_qed, -1);  

  hdipho_Invmass_cep->Add(hdipho_Invmass_qed,1);
  TH1D* hInvmass_data_bkgSub = (TH1D*)hdipho_Invmass_data->Clone();
  hInvmass_data_bkgSub->Add(hdipho_Invmass_cep,-1);
*/  
//hInvmass_data_bkgSub->Add(hdipho_Invmass_cep, -1);
////////////////////////////  
//  TH1D* hInvmass_data_bkgSub = (TH1D*)hdipho_Invmass_data->Clone();
//  hInvmass_data_bkgSub->Add(hdipho_Invmass_qed, -1);
 // hInvmass_data_bkgSub->Add(hdipho_Invmass_cep, -1);

//cout << "Test2:" << endl;

  // "===================     Unfold data ======================     
  TCanvas *c2 = new TCanvas("c2","Dimuon pt ",1100,500);
  c2->Divide(2,1);
  c2->cd(1);
//      hRecoDataInvmass_xSec = getInvXSecHist(hInvmass_data_bkgSub, "hRecoDataInvmass_xSec", 1.64718);    
      
      RooUnfoldBayes unfold_dataPt(&responsePt, hPt_data_bkgSub, itr); // unfold data
      //RooUnfoldBayes unfold_dataRapidity(&responseRapidity, hRapidity_data_bkgSub, itr);
      RooUnfoldInvert unfold_dataRapidity(&responseRapidity, hRapidity_data_bkgSub);
     // RooUnfoldBayes unfold_dataInvmass(&responseInvmass, hInvmass_data_bkgSub, itr);
      RooUnfoldInvert unfold_dataInvmass(&responseInvmass, hInvmass_data_bkgSub);
      RooUnfoldInvert unfold_dataCosthetastar(&responseCosthetastar, hCosthetastar_data_bkgSub);
  
      RooUnfoldBayes unfold_dataPt_Invert(&responsePt, hPt_data_bkgSub, itr);  
      RooUnfoldInvert unfold_dataRapidity_Invert(&responseRapidity, hRapidity_data_bkgSub);
      RooUnfoldInvert unfold_dataInvmass_Invert(&responseInvmass, hInvmass_data_bkgSub); 
      RooUnfoldInvert unfold_dataCosthetastar_Invert(&responseCosthetastar, hCosthetastar_data_bkgSub);

  //unfold_dataPt.IncludeSystematics(4);
  TH1D* hUnfoldPt_data= (TH1D*) unfold_dataPt.Hunfold();
  unfold_dataPt.PrintTable (cout, hdipho_GenPt[0]);
 // unfold_dataRapidity.IncludeSystematics(4);
   TH1D* hUnfoldRapidity_data= (TH1D*) unfold_dataRapidity.Hunfold();
  unfold_dataRapidity.PrintTable (cout, hdipho_GenRapidity[0]);

  TH1D* hUnfoldInvmass_data= (TH1D*) unfold_dataInvmass.Hunfold();
  unfold_dataInvmass.PrintTable (cout, hdipho_GenInvmass[0]);

  TH1D* hUnfoldCosthetastar_data= (TH1D*) unfold_dataCosthetastar.Hunfold();
  unfold_dataCosthetastar.PrintTable (cout, hdipho_GenCosthetastar[0]);

  cout << "Invert1" << endl;
  
 // unfold_dataPt_Invert.IncludeSystematics(4);
  TH1D* hUnfoldPt_Invert_data= (TH1D*) unfold_dataPt_Invert.Hunfold();
  unfold_dataPt_Invert.PrintTable (cout, hdipho_GenPt[0]);
  
 // unfold_dataRapidity.IncludeSystematics(4); 
  TH1D* hUnfoldRapidity_Invert_data= (TH1D*) unfold_dataRapidity_Invert.Hunfold();
  unfold_dataRapidity_Invert.PrintTable (cout, hdipho_GenRapidity[0]);
 //// unfold_dataRapidity.IncludeSystematics(4);
  TH1D* hUnfoldInvmass_Invert_data= (TH1D*) unfold_dataInvmass_Invert.Hunfold();
  unfold_dataInvmass_Invert.PrintTable (cout, hdipho_GenInvmass[0]);

  TH1D* hUnfoldCosthetastar_Invert_data= (TH1D*) unfold_dataCosthetastar_Invert.Hunfold();
  unfold_dataCosthetastar_Invert.PrintTable (cout, hdipho_GenCosthetastar[0]);

  /*****************Covariance matrix***************************/
    TCanvas*c8 =  new TCanvas("c8","covarianceMatrixPt", 800,600);
        c8->SetRightMargin(0.15);
        TMatrixD  covarianceMatrix_Pt = unfold_dataPt.Eunfold();
        TH2D *hCov_Pt = new TH2D("hCov_Pt", "covarianceMatrix_Pt", covarianceMatrix_Pt.GetNrows(), 0, covarianceMatrix_Pt.GetNrows(), covarianceMatrix_Pt.GetNcols(), 0, covarianceMatrix_Pt.GetNcols());
        for (int i = 0; i < covarianceMatrix_Pt.GetNrows(); ++i) {
             for (int j = 0; j < covarianceMatrix_Pt.GetNcols(); ++j) {
                  hCov_Pt->SetBinContent(i + 1, j + 1,covarianceMatrix_Pt[i][j]);
                  }
        }
      
       hCov_Pt->Draw("textcolz");
       c8->Update();
       c8->SaveAs("fig_Bayes_Inverse/Cov_Pt_lbl.pdf");

  TCanvas*c9 =  new TCanvas("c9","covarianceMatrixRapidity", 800,600);
        c9->SetRightMargin(0.15);
        TMatrixD  covarianceMatrix_Rapidity = unfold_dataRapidity_Invert.Eunfold();
        TH2D *hCov_Rapidity = new TH2D("hCov_Rapidity", "covarianceMatrix_Rapidity; Bin number; Bin number", covarianceMatrix_Rapidity.GetNrows(), 0, covarianceMatrix_Rapidity.GetNrows(), covarianceMatrix_Rapidity.GetNcols(), 0, covarianceMatrix_Rapidity.GetNcols());
        for (int i = 0; i < covarianceMatrix_Rapidity.GetNrows(); ++i) {
             for (int j = 0; j < covarianceMatrix_Rapidity.GetNcols(); ++j) {
                  hCov_Rapidity->SetBinContent(i + 1, j + 1,covarianceMatrix_Rapidity[i][j]);
                  }
        }
      
       hCov_Rapidity->Draw("textcolz");
       c9->Update();
       c9->SaveAs("fig_Bayes_Inverse/Cov_Rapidity_lbl.pdf");
  
  TCanvas*c10 =  new TCanvas("c10","covarianceMatrixInvmass", 800,600);
        c10->SetRightMargin(0.15);
       // c10->SetLeftMargin(0.15);
        TMatrixD  covarianceMatrix_Invmass = unfold_dataInvmass_Invert.Eunfold();
        TH2D *hCov_Invmass= new TH2D("hCov_Invmass", "covarianceMatrix_Invmass; Bin number; Bin number", covarianceMatrix_Invmass.GetNrows(), 0, covarianceMatrix_Invmass.GetNrows(), covarianceMatrix_Invmass.GetNcols(), 0, covarianceMatrix_Invmass.GetNcols());
        for (int i = 0; i < covarianceMatrix_Invmass.GetNrows(); ++i) {
             for (int j = 0; j < covarianceMatrix_Invmass.GetNcols(); ++j) {
                  hCov_Invmass->SetBinContent(i + 1, j + 1,covarianceMatrix_Invmass[i][j]);
                  }
        }
      
       hCov_Invmass->Draw("textcolz");
       c10->Update();
       c10->SaveAs("fig_Bayes_Inverse/Cov_Invmass_lbl.pdf");

   TCanvas*c11 =  new TCanvas("c11","covarianceMatrixCosthetastar", 800,600);
        c11->SetRightMargin(0.15);
        TMatrixD  covarianceMatrix_Costhetastar = unfold_dataCosthetastar_Invert.Eunfold();
        TH2D *hCov_Costhetastar= new TH2D("hCov_Costhetastar", "covarianceMatrix_Costhetastar; Bin number; Bin number", covarianceMatrix_Costhetastar.GetNrows(), 0, covarianceMatrix_Costhetastar.GetNrows(), covarianceMatrix_Costhetastar.GetNcols(), 0, covarianceMatrix_Costhetastar.GetNcols());
        for (int i = 0; i < covarianceMatrix_Costhetastar.GetNrows(); ++i) {
             for (int j = 0; j < covarianceMatrix_Costhetastar.GetNcols(); ++j) {
                  hCov_Costhetastar->SetBinContent(i + 1, j + 1,covarianceMatrix_Costhetastar[i][j]);
                  }
        }
      
       hCov_Costhetastar->Draw("textcolz");
       c11->Update();
       c11->SaveAs("fig_Bayes_Inverse/Cov_Costhetastar_lbl.pdf");

       /****************************************************************/


  


  cout << "Invert2:" << endl;
  // divide histograms by luminosty and binwidth to get the cross-sections. 
  hGenPt_xSec  = getXSecHist(hdipho_GenPt[0], "hGenPt_xSec", 5, 0, 1 ,1 );
  hGenRap_xSec = getXSecHist(hdipho_GenRapidity[0], "hGenRap_xSec",rapbin,-2.2,2.2, 1 );
//  hGenRap_xSec = getRapXSecHist(hdipho_GenRapidity[0], "hGenRap_xSec", 1 );
  hGenInvmass_xSec = getInvXSecHist(hdipho_GenInvmass[0], "hGenInvmass_xSec" ,1 );
  hGenCosthetastar_xSec = getXSecHist(hdipho_GenCosthetastar[0], "hGenCosthetastar_xSec" ,cosbin,0,1,1 );
  

  hRecoMCPt_xSec  = getXSecHist(hdipho_Pt[0], "hRecoMCPt_xSec", 5, 0, 1,1 );
  hRecoMCRap_xSec = getXSecHist(hdipho_Rapidity[0], "hRecoMCRap_xSec", rapbin,-2.2,2.2,1 );
 // hRecoMCRap_xSec = getRapXSecHist(hdipho_Rapidity[0], "hRecoMCRap_xSec", 1 );
  hRecoMCInvmass_xSec = getInvXSecHist(hdipho_Invmass[0], "hRecoMCInvmass_xSec",1 );
  hRecoMCCosthetastar_xSec = getXSecHist(hdipho_Costhetastar[0], "hRecoMCCosthetastar_xSec",cosbin,0,1,1 );
  

/*
  hRecoDataPt_xSec  = getXSecHist(hdipho_Pt_data, "hRecoDataPt_xSec", 10, 0, 2, 1.639);
  //hRecoDataRap_xSec = getXSecHist(hdipho_Rapidity_data, "hRecoDataRap_xSec",2,0,1.6, 1.639);
  hRecoDataRap_xSec = getRapXSecHist(hdipho_Rapidity_data, "hRecoDataRap_xSec", 1.639);
  hRecoDataInvmass_xSec = getInvXSecHist(hdipho_Invmass_data, "hRecoDataInvmass_xSec", 1.639);
*/

  hRecoDataPt_xSec  = getXSecHist(hPt_data_bkgSub, "hRecoDataPt_xSec", 5, 0, 1, 1.64718);
  hRecoDataRap_xSec = getXSecHist(hRapidity_data_bkgSub, "hRecoDataRap_xSec",rapbin,-2.2,2.2, 1.64723);
// hRecoDataRap_xSec = getRapXSecHist(hRapidity_data_bkgSub, "hRecoDataRap_xSec", 1.64718);
  hRecoDataInvmass_xSec = getInvXSecHist(hInvmass_data_bkgSub, "hRecoDataInvmass_xSec", 1.64718);
  hRecoDataCosthetastar_xSec = getXSecHist(hCosthetastar_data_bkgSub, "hRecoDataCosthetastar_xSec",cosbin,0,1,1.64718 );


  hUnfoMCPt_xSec  = getXSecHist(hUnfoldPt,       "hUnfoMCPt_xSec", 5, 0, 1,1 );
  hUnfoMCRap_xSec = getXSecHist(hUnfoldRapidity, "hUnfoMCRap_xSec",rapbin,-2.2,2.2, 1 );
 // hUnfoMCRap_xSec = getRapXSecHist(hUnfoldRapidity, "hUnfoMCRap_xSec", 1 );
  hUnfoMCInvmass_xSec = getInvXSecHist(hUnfoldInvmass, "hUnfoMCInvmass_xSec", 1 );
  hUnfoMCCosthetastar_xSec = getXSecHist(hUnfoldCosthetastar, "hUnfoMCCosthetastar_xSec", cosbin,0,1,1 );

  hUnfoMCPt_Invert_xSec  = getXSecHist(hUnfoldPt_Invert,       "hUnfoMCPt_Invert_xSec", 5, 0, 1,1 );
  hUnfoMCRap_Invert_xSec = getXSecHist(hUnfoldRapidity_Invert, "hUnfoMCRap_Invert_xSec",rapbin,-2.2,2.2, 1 ); 
  hUnfoMCInvmass_Invert_xSec = getInvXSecHist(hUnfoldInvmass_Invert, "hUnfoMCInvmass_Invert_xSec", 1 );
 //hUnfoMCRap_Invert_xSec = getRapXSecHist(hUnfoldRapidity_Invert, "hUnfoMCRap_Invert_xSec", 1 ); 
  hUnfoMCCosthetastar_Invert_xSec = getXSecHist(hUnfoldCosthetastar_Invert, "hUnfoMCCosthetastar_Invert_xSec", cosbin,0,1,1 );
  


  hUnfoDataPt_xSec  = getXSecHist(hUnfoldPt_data, "hUnfoDataPt_xSec", 5, 0, 1, 1.64718);
  hUnfoDataRap_xSec = getXSecHist(hUnfoldRapidity_data, "hUnfoDataRap_xSec",rapbin,-2.2,2.2, 1.64723);
  //hUnfoDataRap_xSec = getRapXSecHist(hUnfoldRapidity_data, "hUnfoDataRap_xSec",  1.64718);
  hUnfoDataInvmass_xSec = getInvXSecHist(hUnfoldInvmass_data, "hUnfoDataInvmass_xSec", 1.64718);
  hUnfoDataCosthetastar_xSec = getXSecHist(hUnfoldCosthetastar_data, "hUnfoDataCosthetastar_xSec",cosbin,0,1, 1.64718);

   hUnfoDataPt_Invert_xSec  = getXSecHist(hUnfoldPt_Invert_data, "hUnfoDataPt_Invert_xSec", 5, 0, 1, 1.64718);
  // hUnfoDataRap_Invert_xSec = getRapXSecHist(hUnfoldRapidity_Invert_data, "hUnfoDataRap_Invert_xSec",  1.64718);
   hUnfoDataRap_Invert_xSec = getXSecHist(hUnfoldRapidity_Invert_data, "hUnfoDataRap_Invert_xSec", rapbin,-2.2,2.2, 1.64718);
  hUnfoDataInvmass_Invert_xSec = getInvXSecHist(hUnfoldInvmass_Invert_data, "hUnfoDataInvmass_Invert_xSec", 1.64718);
  hUnfoDataCosthetastar_Invert_xSec = getXSecHist(hUnfoldCosthetastar_Invert_data, "hUnfoDataCosthetastar_Invert_xSec",cosbin,0,1, 1.64718);
 
/********************/
   cout << " Pt Gen " << Sample[0] <<  " :" << hdipho_GenPt[0]->Integral() << endl; 
   cout << " Rapidity Gen " << Sample[0] <<  " :" << hdipho_GenRapidity[0]->Integral() << endl; 
  cout << " invariant mass Gen " << Sample[0] <<  " :" << hdipho_GenInvmass[0]->Integral() << endl; 
  cout << " Costhetatsar Gen " << Sample[0] <<  " :" << hdipho_GenCosthetastar[0]->Integral() << endl;
  
  cout << " Rapidity reco " << Sample[0] <<  " :" << hdipho_Rapidity[0]->Integral() << endl;  
  cout << " Pt reco " << Sample[0] <<  " :" << hdipho_Pt[0]->Integral() << endl;  
  cout << " invariant mass reco " << Sample[0] <<  " :" << hdipho_Invmass[0]->Integral() << endl;  
  cout << " Costhetastar reco " << Sample[0] <<  " :" << hdipho_Costhetastar[0]->Integral() << endl;

  cout << " invariant mass data" << Sample[0] <<  " :" << hdipho_Invmass_data->Integral() << endl;
  cout << " invariant mass data w/o bkg" << Sample[0] <<  " :" << hInvmass_data_bkgSub->Integral() << endl;
  cout << " Rapidity data " << Sample[0] <<  " :" << hdipho_Rapidity_data->Integral() << endl;
  cout << " Rapidity data w/o bkg " << Sample[0] <<  " :" << hRapidity_data_bkgSub->Integral() << endl;
  cout << " costhetastar data " << Sample[0] <<  " :" << hdipho_Costhetastar_data->Integral() << endl;
  cout << " costhetastar data w/o bkg " << Sample[0] <<  " :" << hCosthetastar_data_bkgSub->Integral() << endl;
  
  cout << " Pt data " << Sample[0] <<  " :" << hdipho_Rapidity_data->Integral() << endl;
  cout << " Pt data w/o bkg " << Sample[0] <<  " :" << hRapidity_data_bkgSub->Integral() << endl;
  cout << " Pt QED_Jeremi " << Sample[0] <<  " :" << hdipho_Rapidity_qed->Integral() << endl;
  cout << " Pt CEP_Jeremi " << Sample[0] <<  " :" << hdipho_Rapidity_cep->Integral() << endl;

  cout << "Cross section Pt Data:" <<  hRecoDataPt_xSec->Integral() << endl;
  cout << "Cross section Rapidity Data:" <<  hRecoDataRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass Data:" << hRecoDataInvmass_xSec->Integral()<< endl;
  cout << "Cross section Costhetastar Data:" << hRecoDataCosthetastar_xSec->Integral()<< endl;

  cout << "Cross section Pt Reco MC:" <<   hRecoMCPt_xSec->Integral() << endl;
  cout << "Cross section Rapidity Reco MC:" <<  hRecoMCRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass Reco MC:" << hRecoMCInvmass_xSec->Integral()<< endl;
  cout << "Cross section Costhetastar Reco MC:" << hRecoMCCosthetastar_xSec->Integral()<< endl;

  cout << "Cross section Pt Gen MC:" <<   hGenPt_xSec->Integral() << endl;
  cout << "Cross section Rapidity Gen MC:" <<  hGenRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass Gen  MC:" << hGenInvmass_xSec->Integral()<< endl;
  cout << "Cross section Costhetastar Gen  MC:" << hGenCosthetastar_xSec->Integral()<< endl;

  cout << "Cross section Pt UnfoldedData:" <<   hUnfoDataPt_xSec->Integral() << endl;
  cout << "Cross section Rapidity UnfoldedData:" <<  hUnfoDataRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass UnfoldedData:" << hUnfoDataInvmass_xSec->Integral()<< endl;
  cout << "Cross section Costhetastar UnfoldedData:" << hUnfoDataCosthetastar_xSec->Integral()<< endl;

  cout << "Cross section Pt UnfoldedData_Invert:" <<   hUnfoDataPt_Invert_xSec->Integral() << endl;
  cout << "Cross section Rapidity UnfoldedData_Invert:" <<  hUnfoDataRap_Invert_xSec->Integral()<< endl;
  cout << "Cross section Invmass UnfoldedData_Invert:" << hUnfoDataInvmass_Invert_xSec->Integral()<< endl;
  cout << "Cross section Costhetastar UnfoldedData_Invert:" << hUnfoDataCosthetastar_Invert_xSec->Integral()<< endl;
  
  outf->cd();
  //outf->Write();
///
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
  
  hUnfoDataPt_Invert_xSec->Write();
  hUnfoDataRap_Invert_xSec->Write();
  hUnfoDataInvmass_Invert_xSec->Write();
  hUnfoDataCosthetastar_Invert_xSec->Write();


  TCanvas*c1 =  new TCanvas();
  hdipho_Rapidity[0]->SetLineColor(kRed);
  hdipho_Rapidity[0]->Draw("p");;
   


  outf->Close();
  //hUnfold->Write();
  //hUnfoldRapidity->Write();
  //hUnfoldRapidity_data->Write();
  //hUnfoldInvmass->Write();
  
}

TH1D* subtractBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep, TH1D *hqed_sl, const char* name, int nbins, float xmin, float xmax){

  TH1D *hans = new TH1D(name,"",nbins,xmin,xmax);
  //TH1D* hans = (TH1D*)hdata->Clone();
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hdata->GetNbinsX(); i++) {
      
      double binCont  = hdata->GetBinContent(i)-hqed->GetBinContent(i)-hcep->GetBinContent(i)-hqed_sl->GetBinContent(i);
     
      double binError = sqrt(pow(hdata->GetBinError(i),2) + pow(hqed->GetBinError(i),2)  + pow(hcep->GetBinError(i),2) + pow(hqed->GetBinError(i),2));

     
      if(binCont <=0){
     
        hans->SetBinContent(i, 0*0);
         hans->SetBinError(i, 0);
      }
      else{
      hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError);
      }
      
      //cout << "******bin:" << i << " " << hdata->GetBinContent(i) << " " << hqed->GetBinContent(i) << " " << hcep->GetBinContent(i) << "   " <<  hans->GetBinContent(i) << "+/-" <<  binError<< endl;
      //double err=0;
      //for (int ivar=1; ivar<14; ivar++) err += pow(hvari[ivar]->GetBinContent(i)-hvari[0]->GetBinContent(i),2);
      //hans->SetBinError(i,sqrt(err+pow(hans->GetBinError(i),2)+pow(glob_syst,2)));
   }

   return hans;
}

TH1D* subtractInvBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep, TH1D *hqed_sl, const char* name) {
    // Clone the data histogram to create the result histogram
    TH1D *hans = (TH1D*)hdata->Clone(name);

    // Add the histograms with unequal binning
    hdata->Add(hqed, -1);
    hdata->Add(hcep, -1);
    hdata->Add(hqed_sl, -1);
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
    hdata->Add(hqed_sl);

    return hans;
}


TH1D* getXSecHist(TH1D *hist, const char* name, int nbins, float xmin, float xmax, float luminosity){

  TH1D *hans = new TH1D(name,"",nbins,xmin,xmax);
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hist->GetNbinsX(); i++) {

     double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));
    //  double binCont  = hist->GetBinContent(i)/(luminosity);

      double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
     // double binError = hist->GetBinError(i)/hist->GetBinWidth(i);
      if(binCont <= 0){ hans->SetBinContent(i, 0);
       hans->SetBinError(i, 0);
       }
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
      if(binCont <=0){
       hans->SetBinContent(i, 0);
       hans->SetBinError(i,0);
    
       }
	else{
     hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError); 
    }

    //  if(i==4){
    //   hans->SetBinError(i,0);
    //  }
      cout << "binerror:" << hans->GetBinError(4) << endl;
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

     //double binCont  = hist->GetBinContent(i)/(luminosity);
     double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));

      double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
     // double binError = hist->GetBinError(i)/hist->GetBinWidth(i);
      if(binCont == 0){ hans->SetBinContent(i, 0);
       hans->SetBinError(i, 0);
       }
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

//Costhetastar
double costhetastar_CS(TLorentzVector p1, TLorentzVector pair)
{
  double costhetastarCS = 2.*(pair.E()*p1.Pz()-pair.Pz()*p1.E())/(pair.M()*pair.Mt());

  return costhetastarCS;
}



