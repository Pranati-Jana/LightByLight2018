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

//const double luminosity       = 1635.123139823; // μb^-1
const double luminosity       = 1;
//const double luminosity       = 1.6722;
//const double luminosity       = 1639.207543; // μb^-1
//const double luminosity       = 1.639207543; // μb^-1
double scaleFactorPhoton_Old = 0.85 * 0.93 * 0.866 * 1.008 * 1.037 * 1.037;

double scaleFactorPhoton = 0.8477 *  //0.85// NEE    21.12.2021
  0.9322 *      //0.93// CHE  21.12.2021
  pow(0.9771, 2)* //1.048// Photon  reco+ID 21.12.2021
  0.8643 *      //0.866// HF veto
  1.0006;       //1.008// L1 EG trigger

//const double xsecGeneratedLbLSC    = 2.59; // μb Superchic
const double xsecGeneratedLbLSC    = 2590; // nb Superchic
const double nEventsGeneratedLbLSC = 466000;  //Superchic
double norm_LbLSC = scaleFactorPhoton*xsecGeneratedLbLSC*luminosity/nEventsGeneratedLbLSC;  //31.12.21 (SF = 1.048, took the square here). 
//double norm_LbLSC = 1;
double W1 = xsecGeneratedLbLSC*luminosity/nEventsGeneratedLbLSC; 
//double W1 = 1;
//double W1 = scaleFactorPhoton*xsecGeneratedLbLSC*luminosity/nEventsGeneratedLbLSC;
//const double xsecGeneratedLbLMG    = 0.1406; // μb Madgraph David
const double xsecGeneratedLbLMG    = 140.6; // nb Madgraph David
const double nEventsGeneratedLbLMG = 788069; //Madgraph David
double norm_LbLMG = scaleFactorPhoton*xsecGeneratedLbLMG*luminosity/nEventsGeneratedLbLMG;  //31.12.21 (SF = 1.048, took the square here). 


void printOutput(TH1D *hist);

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


const int nRapbins=3;
//double Rapbin[nRapbins]={-2.5,-0.75,0.75,2.5};
//double Rapbin[nRapbins]={0.0,0.5,1.0,1.5,2.2};
double Rapbin[nRapbins]={0.0, 0.8, 2.2};
const int nRapbin= sizeof(Rapbin)/sizeof(double) - 1;

const int nMassbins=4;
//double Massbin[nMassbins]={5.0,7.0,10.0,30.0};
//double Massbin[nMassbins]={5.0,7.5,10.5,16.0};
//double Massbin[nMassbins]={5.0,7.3,13.0,16.0};
//double Massbin[nMassbins]={5.0,7.0,10.0,16.0};
double Massbin[nMassbins]={5.0,9.0,13.0,17.0};
const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;


TH1D* subtractBkg(TH1D *hdata, TH1D *hqed, TH1D *hcep, const char* name, int nbins, float xmin, float xmax);

TH1D* getXSecHist(TH1D *hist, const char* name, int nbins, float xmin, float xmax, float lumi);
TH1D* getInvXSecHist(TH1D *hist, const char* name, float lumi);
TH1D* getRapXSecHist(TH1D *hist, const char* name, float lumi);

TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, TH1D*, double , double , double , double , double , double , const char *, bool); 
void makeUnfold(int method=1){
  int itr = 4; 
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


  bool absRap = true;

  TFile *outf= new TFile(Form("unfolding_histograms_Bayes_%s.root",mc_name.c_str()),"recreate");  
  
  TChain *lbyl[nSample];
  
  TH1D* hdipho_Pt[nSample], *hdipho_Rapidity[nSample], *hdipho_Invmass[nSample], *hAcoplanarity[nSample] ; TH1F *hSF[nSample]; 
  
  TH1D* hdipho_GenPt[nSample], *hdipho_GenRapidity[nSample], *hdipho_GenInvmass[nSample];
  TH2D* hdipho_RecoGenPt[nSample], *hdipho_RecoGenRapidity[nSample], *hRecoGenInvmass[nSample];


  //RooUnfoldResponse responseRapidity (hdipho_Rapidity, hdipho_GenRapidity);
  //RooUnfoldResponse responseInvmass (hdipho_Invmass, hdipho_GenInvmass);
  
  cout << "define hists"<< endl;
  
  int i =0;
  //for (int i = 0; i < nSample; i++){
    lbyl[i] = new TChain("output_tree");
    lbyl[i]->Add(Form("RootFile/lbylSC_diphoton_genInfo.root"));


    
    hdipho_Pt[i]    = new TH1D(Form("hdipho_Pt%s", Sample[i]),"",5,0,1);
    hdipho_Rapidity[i]   = new TH1D(Form("hdipho_Rapidity%s",Sample[i]),"",nRapbin, Rapbin);
//    hdipho_Rapidity[i]   = new TH1D(Form("hdipho_Rapidity%s",Sample[i]),"",2,0,2.2);
    hdipho_Invmass[i]   = new TH1D(Form("hdipho_Invmass%s",Sample[i]),"",nMassbin, Massbin);
    hAcoplanarity[i]   = new TH1D(Form("hAcoplanarity%s",Sample[i]),"",50,0,0.1);
    
    hdipho_GenPt[i]         = new TH1D(Form("hdipho_GenPt%s", Sample[i]),"",5,0,1);
    hdipho_GenRapidity[i]   = new TH1D(Form("hdipho_GenRapidity%s",Sample[i]),"",nRapbin, Rapbin);
  //  hdipho_GenRapidity[i]   = new TH1D(Form("hdipho_GenRapidity%s",Sample[i]),"",2,0,2.2);
    hdipho_GenInvmass[i]          = new TH1D(Form("hdipho_GenInvmass%s",Sample[i]),"",nMassbin, Massbin);
    
    hdipho_RecoGenPt[i]         = new TH2D(Form("hdipho_RecoGenPt%s", Sample[i]),"",5,0,1,5,0,1);
    hdipho_RecoGenRapidity[i]   = new TH2D(Form("hdipho_RecoGenRapidity%s",Sample[i]),"",nRapbin, Rapbin,nRapbin, Rapbin);
   // hdipho_RecoGenRapidity[i]   = new TH2D(Form("hdipho_RecoGenRapidity%s",Sample[i]),"",2,0,2.2,2,0,2.2);
    hRecoGenInvmass[i]          = new TH2D(Form("hRecoGenInvmass%s",Sample[i]),"",nMassbin, Massbin,nMassbin, Massbin);


  // ===================   Cross section histograms   ======================     

  TH1D* hGenPt_xSec              = new TH1D("hGenPt_xSec", "Pt gen MC",5,0,1);
  TH1D* hGenRap_xSec             = new TH1D("hGenRap_xSec","Rap gen MC",nRapbin, Rapbin);
  //TH1D* hGenRap_xSec             = new TH1D("hGenRap_xSec","Rap gen MC",2,0,2.2);
  TH1D* hGenInvmass_xSec         = new TH1D("hGenInvmass_xSec","Invmass gen MC",nMassbin, Massbin);

  TH1D* hRecoMCPt_xSec              = new TH1D("hRecoMCPt_xSec", "Pt gen MC",5,0,1);
  TH1D* hRecoMCRap_xSec             = new TH1D("hRecoMCRap_xSec","Rap gen MC",nRapbin, Rapbin);
  //TH1D* hRecoMCRap_xSec             = new TH1D("hRecoMCRap_xSec","Rap gen MC",2,0,2.2);
  TH1D* hRecoMCInvmass_xSec         = new TH1D("hRecoMCInvmass_xSec","Invmass gen MC",nMassbin, Massbin);

  TH1D* hRecoDataPt_xSec              = new TH1D("hRecoDataPt_xSec", "Pt gen MC",5,0,1);
  TH1D* hRecoDataRap_xSec             = new TH1D("hRecoDataRap_xSec","Rap gen MC",nRapbin, Rapbin);
  //TH1D* hRecoDataRap_xSec             = new TH1D("hRecoDataRap_xSec","Rap gen MC",2,0,2.2);
  TH1D* hRecoDataInvmass_xSec         = new TH1D("hRecoDataInvmass_xSec","Invmass gen MC",nMassbin, Massbin);


  TH1D* hUnfoMCPt_xSec              = new TH1D("hUnfoMCPt_xSec", "Pt gen MC",5,0,1);
  TH1D* hUnfoMCRap_xSec             = new TH1D("hUnfoMCRap_xSec","Rap gen MC",nRapbin, Rapbin);
  //TH1D* hUnfoMCRap_xSec             = new TH1D("hUnfoMCRap_xSec","Rap gen MC",2,0,2.2);
  TH1D* hUnfoMCInvmass_xSec         = new TH1D("hUnfoMCInvmass_xSec","Invmass gen MC",nMassbin, Massbin);
  TH1D* hUnfoMCInvmass_Invert_xSec         = new TH1D("hUnfoMCInvmass_Invert_xSec","Invmass gen MC",nMassbin, Massbin);
  TH1D* hUnfoMCPt_Invert_xSec              = new TH1D("hUnfoMCPt_Invert_xSec", "Pt gen MC",5,0,1);
  TH1D* hUnfoMCRap_Invert_xSec             = new TH1D("hUnfoMCRap_Invert_xSec","Rap gen MC",nRapbin, Rapbin);

  TH1D* hUnfoDataPt_xSec              = new TH1D("hUnfoDataPt_xSec", "Pt gen MC",5,0,1);
  TH1D* hUnfoDataRap_xSec             = new TH1D("hUnfoDataRap_xSec","Rap gen MC",nRapbin, Rapbin);
  //TH1D* hUnfoDataRap_xSec             = new TH1D("hUnfoDataRap_xSec","Rap gen MC",2,0,2.2);
  TH1D* hUnfoDataInvmass_xSec         = new TH1D("hUnfoDataInvmass_xSec","Invmass gen MC",nMassbin, Massbin);
  TH1D* hUnfoDataInvmass_Invert_xSec         = new TH1D("hUnfoDataInvmass_Invert_xSec","Invmass gen MC",nMassbin, Massbin);
  TH1D* hUnfoDataPt_Invert_xSec              = new TH1D("hUnfoDataPt_Invert_xSec", "Pt gen MC",5,0,1);
  TH1D* hUnfoDataRap_Invert_xSec             = new TH1D("hUnfoDataRap_Invert_xSec","Rap gen MC",nRapbin, Rapbin);




/*    RooUnfoldResponse responsePt (hdipho_Pt[i], hdipho_GenPt[i]);
    RooUnfoldResponse responseRapidity (hdipho_Rapidity[i], hdipho_GenRapidity[i]);
    RooUnfoldResponse responseInvmass (hdipho_Invmass[i], hdipho_GenInvmass[i]);
 */   
    RooUnfoldResponse responsePt (hdipho_Pt[i], hdipho_GenPt[i],hdipho_RecoGenPt[i]);
    RooUnfoldResponse responseRapidity (hdipho_Rapidity[i], hdipho_GenRapidity[i],hdipho_RecoGenRapidity[i]);
    RooUnfoldResponse responseInvmass (hdipho_Invmass[i], hdipho_GenInvmass[i],hRecoGenInvmass[i]);

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
      if(abs(lbylR.gen_Eta1) > 2.2 || abs(lbylR.gen_Eta2) > 2.2 || lbylR.gen_Pt1 < 2.5 || lbylR.gen_Pt2 < 2.5 || lbylR.gen_diPho_M < 5) continue;
//      if(abs(lbylR.gen_diPho_Rapidity) > 2.2) continue;
      
      hdipho_GenPt[i]->Fill(lbylR.gen_diPho_Pt,W1);
      if(absRap)hdipho_GenRapidity[i]->Fill(abs(lbylR.gen_diPho_Rapidity),W1);
      if(!absRap)hdipho_GenRapidity[i]->Fill(lbylR.gen_diPho_Rapidity,W1);
      hdipho_GenInvmass[i]->Fill(lbylR.gen_diPho_M,W1);

      // ===================   Apply analysis cuts on reco MC ======================           
      bool isPassing = lbylR.phoSwissCross_1 < 0.95 && lbylR.phoSwissCross_2 <  0.95
		        && lbylR.ok_neuexcl == 1 && lbylR.ok_chexcl_goodtracks == 1 && lbylR.ok_chexcl_goodelectrons == 1 
		        && lbylR.vSum_M > 5  && lbylR.vSum_Pt < 1 
 			&& abs(lbylR.phoEta_1) < 2.2 && abs(lbylR.phoEta_2) < 2.2 && lbylR.phoEt_1 > 2.5 && lbylR.phoEt_2 > 2.5		
		 	&& lbylR.pho_acop < 0.01;
			
      if(lbylR.ok_trigger == 1) { //trigger 
	if(isPassing){
	
	hdipho_Pt[i]->Fill(lbylR.vSum_Pt,wt[i]);
	if(absRap)hdipho_Rapidity[i]->Fill(abs(lbylR.vSum_Rapidity),wt[i]);
	if(!absRap)hdipho_Rapidity[i]->Fill(lbylR.vSum_Rapidity,wt[i]);
	hdipho_Invmass[i]->Fill(lbylR.vSum_M,wt[i]);
	
	hdipho_RecoGenPt[i]->Fill(lbylR.vSum_Pt,lbylR.gen_diPho_Pt, wt[i]);
	if(absRap)hdipho_RecoGenRapidity[i]->Fill(abs(lbylR.vSum_Rapidity),abs(lbylR.gen_diPho_Rapidity),wt[i]);
	if(!absRap)hdipho_RecoGenRapidity[i]->Fill(lbylR.vSum_Rapidity,lbylR.gen_diPho_Rapidity,wt[i]);
	hRecoGenInvmass[i]->Fill(lbylR.vSum_M,lbylR.gen_diPho_M,wt[i]);
	
        responsePt.Fill(lbylR.vSum_Pt,lbylR.gen_diPho_Pt,wt[i]);
	if(absRap)responseRapidity.Fill(abs(lbylR.vSum_Rapidity),abs(lbylR.gen_diPho_Rapidity),wt[i]);
	if(!absRap)responseRapidity.Fill(lbylR.vSum_Rapidity,lbylR.gen_diPho_Rapidity,wt[i]);
	responseInvmass.Fill(lbylR.vSum_M,lbylR.gen_diPho_M,wt[i]);

	} //isPassing analysis cuts
	else {
	responsePt.Miss(lbylR.gen_diPho_Pt,W1);
	if(absRap)responseRapidity.Miss(abs(lbylR.gen_diPho_Rapidity),W1);
	if(!absRap)responseRapidity.Miss(lbylR.gen_diPho_Rapidity,W1);
	responseInvmass.Miss(lbylR.gen_diPho_M,W1);
	}
      } //trigger
	else {
	responsePt.Miss(lbylR.gen_diPho_Pt,W1);
	if(absRap)responseRapidity.Miss(abs(lbylR.gen_diPho_Rapidity),W1);
	if(!absRap)responseRapidity.Miss(lbylR.gen_diPho_Rapidity,W1);
	responseInvmass.Miss(lbylR.gen_diPho_M,W1);
	}
      
      
    } //entry
 // } // for 4 files
cout << "gen entries pass:" << hdipho_GenPt[0]->GetEntries() << endl; 
cout << "reco entries pass:" << hdipho_Pt[0]->GetEntries() << endl; 
/*******************************Condition*number*********************************/
cout << "Inverse1" << endl;
TDecompSVD *svd= new TDecompSVD (responseRapidity.Mresponse());
auto singular_values = svd->GetSig();
cout << "Test1" << endl;
svd->Print();
double singular_value_min;
      for (int i = 0; i < 4 ; i++)
      {
        cout << "Test2" << endl;
        cout << singular_values[i] << endl;
//       if ( singular_values[i] > pow(10,-15)  ) singular_value_min = singular_values[i]; // the pow(10,-15) > requirement is to suppress singular values from empty bins 
      }

//cout << "condition number: "  << singular_values[0]/singular_value_min << endl;
cout << "Inverse2" << endl;
/*********************************************************************************************/
  cout<< "======================================Response matrix========================="<<endl;
  char MethID2[100]; 
  //if(method==1)
   // {

      RooUnfoldBayes unfold(&responsePt, hdipho_Pt[0], itr);
      RooUnfoldBayes unfold2(&responseRapidity, hdipho_Rapidity[0], itr);
      RooUnfoldBayes unfold3(&responseInvmass, hdipho_Invmass[0], itr);

      RooUnfoldInvert unfold_Invert(&responsePt, hdipho_Pt[0]);
      RooUnfoldInvert unfold2_Invert(&responseRapidity, hdipho_Rapidity[0]);
      RooUnfoldInvert unfold3_Invert(&responseInvmass, hdipho_Invmass[0]);


     //   RooUnfoldSvd unfold(&responsePt, hdipho_Pt[0], 4);
       // RooUnfoldSvd unfold2(&responseRapidity, hdipho_Rapidity[0], 4);
      //  RooUnfoldSvd unfold3(&responseInvmass, hdipho_Invmass[0], 4);
/********************Matrix-plotting************/
        TCanvas*c4 =  new TCanvas("c4", "" ,800,600);
        c4->SetRightMargin(0.15);
        TMatrixD  responseMatrix1 = responsePt.Eresponse();
        TH2D *hResponseMatrix = new TH2D("hResponseMatrix", "Response Matrix;Reco diphoton p_{T};Gen diphoton p_{T}", responseMatrix1.GetNrows(), 0, responseMatrix1.GetNrows(), responseMatrix1.GetNcols(), 0, responseMatrix1.GetNcols());
        for (int i = 0; i < responseMatrix1.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix1.GetNcols(); ++j) {
                  hResponseMatrix->SetBinContent(i + 1, j + 1, responseMatrix1[i][j]);
                  }
        }
        hResponseMatrix->GetXaxis()->SetNdivisions(5);
        hResponseMatrix->GetYaxis()->SetNdivisions(5);
        hResponseMatrix->Draw("textcolz");
        c4->Update();
        c4->SaveAs("fig_Bayes_Inverse/Pt.png");
        //
        TCanvas*c5 =  new TCanvas("c5","Response Matrix2",800,800);
        c5->SetRightMargin(0.15);
        TMatrixD  responseMatrix2 = responseRapidity.Mresponse();
        TH2D *hResponseMatrix2 = new TH2D("hResponseMatrix2", "Response Matrix2;Reco rapidity;Gen rapidity", responseMatrix2.GetNrows(), 0, responseMatrix2.GetNrows(), responseMatrix2.GetNcols(), 0, responseMatrix2.GetNcols());
        for (int i = 0; i < responseMatrix2.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix2.GetNcols(); ++j) {
                  hResponseMatrix2->SetBinContent(i + 1, j + 1, responseMatrix2[i][j]);
                  }
        }
        hResponseMatrix2->GetXaxis()->SetNdivisions(4);
        hResponseMatrix2->GetYaxis()->SetNdivisions(4);
        hResponseMatrix2->Draw("textcolz");
        c5->Update();
        c5->SaveAs("fig_Bayes_Inverse/Rapidity.png");
        //
        TCanvas*c6 =  new TCanvas("c6","Response Matrix3", 800,600);
        c6->SetRightMargin(0.15);
        TMatrixD  responseMatrix3 = responseInvmass.Mresponse();
        TH2D *hResponseMatrix3 = new TH2D("hResponseMatrix3", "Response Matrix3;Reco invariant mass;Gen invariant mass", responseMatrix3.GetNrows(), 0, responseMatrix3.GetNrows(), responseMatrix3.GetNcols(), 0, responseMatrix3.GetNcols());
        for (int i = 0; i < responseMatrix3.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix3.GetNcols(); ++j) {
                  hResponseMatrix3->SetBinContent(i + 1, j + 1, responseMatrix3[i][j]);
                  }
        }
        hResponseMatrix3->Draw("textcolz");
        c6->Update();
        c6->SaveAs("fig_Bayes_Inverse/InvMass.png");
/***********************************************/

//  unfold.IncludeSystematics(2);
  TH1D* hUnfoldPt= (TH1D*) unfold.Hunfold();
  unfold.PrintTable (cout,hdipho_GenPt[0]);

//  unfold2.IncludeSystematics(2);  
  TH1D* hUnfoldRapidity= (TH1D*) unfold2.Hunfold();
  unfold2.PrintTable (cout,hdipho_GenRapidity[0]);

  TH1D* hUnfoldInvmass= (TH1D*) unfold3.Hunfold();
  unfold3.PrintTable (cout,hdipho_GenInvmass[0]);
  
  TH1D* hUnfoldPt_Invert= (TH1D*) unfold_Invert.Hunfold();
  unfold_Invert.PrintTable (cout,hdipho_GenPt[0]);

  TH1D* hUnfoldRapidity_Invert= (TH1D*) unfold2_Invert.Hunfold();
  unfold2_Invert.PrintTable (cout,hdipho_GenRapidity[0]);

  TH1D* hUnfoldInvmass_Invert= (TH1D*) unfold3_Invert.Hunfold();
  unfold3_Invert.PrintTable (cout,hdipho_GenInvmass[0]);



  /***************************  Start reading data and subtract background contributions ********************************/

  //TFile *f1=new TFile("diphoton_histos.root","r"); 
  TFile *f1=new TFile("RootFile/diphoton_histos.root","r"); 
//  TFile *f1=new TFile("RootFile/diphoton_histos_Invmassbin-5-8-12-28.root","r"); 

  TH1D *hdipho_Pt_data        = (TH1D*)f1->Get("hdipho_ptData");
  //TH1D *hdipho_Rapidity_data  = (TH1D*)f1->Get("hdipho_absRapidityData"); 
  TH1D *hdipho_Rapidity_data  = (TH1D*)f1->Get("hdipho_Rapidity_unfoData"); 
//  TH1D *hdipho_Invmass_data   = (TH1D*)f1->Get("hInvmassData"); 
  TH1D *hdipho_Invmass_data   = (TH1D*)f1->Get("hInvmass_unfoData"); 

  TH1D *hdipho_Pt_qed        = (TH1D*)f1->Get("hdipho_ptQEDSCFSR");
  //TH1D *hdipho_Rapidity_qed  = (TH1D*)f1->Get("hdipho_absRapidityQEDSCFSR"); 
  TH1D *hdipho_Rapidity_qed  = (TH1D*)f1->Get("hdipho_Rapidity_unfoQEDSCFSR"); 
  //TH1D *hdipho_Invmass_qed   = (TH1D*)f1->Get("hInvmassQEDSCFSR"); 
  TH1D *hdipho_Invmass_qed   = (TH1D*)f1->Get("hInvmass_unfoQEDSCFSR"); 

  TH1D *hdipho_Pt_cep        = (TH1D*)f1->Get("hdipho_ptCEPIncoh");
 // TH1D *hdipho_Rapidity_cep  = (TH1D*)f1->Get("hdipho_absRapidityCEPIncoh"); 
  TH1D *hdipho_Rapidity_cep  = (TH1D*)f1->Get("hdipho_Rapidity_unfoCEPIncoh"); 
  //TH1D *hdipho_Invmass_cep   = (TH1D*)f1->Get("hInvmassCEPIncoh"); 
  TH1D *hdipho_Invmass_cep   = (TH1D*)f1->Get("hInvmass_unfoCEPIncoh"); 

  // "===================      Subtract background ======================     
  TH1D* hPt_data_bkgSub       = subtractBkg(hdipho_Pt_data,      hdipho_Pt_qed,      hdipho_Pt_cep, "hPt_data_bkgSub", 5,0,1);
 // TH1D* hRapidity_data_bkgSub = subtractBkg(hdipho_Rapidity_data,hdipho_Rapidity_qed,hdipho_Rapidity_cep,  "hRapidity_data_bkgSub", 2,0,2.2);
  //TH1D* hInvmass_data_bkgSub  = subtractBkg(hdipho_Invmass_data,hdipho_Invmass_qed,  hdipho_Invmass_cep, "hInvmass_data_bkgSub", nMassbin, Massbin);

  hdipho_Rapidity_qed->Add(hdipho_Rapidity_cep,1);
  TH1D* hRapidity_data_bkgSub = (TH1D*)hdipho_Rapidity_data->Clone();
  hRapidity_data_bkgSub->Add(hdipho_Rapidity_qed,-1);
 //For ()

  hdipho_Invmass_qed->Add(hdipho_Invmass_cep,1);
  TH1D* hInvmass_data_bkgSub = (TH1D*)hdipho_Invmass_data->Clone();
  hInvmass_data_bkgSub->Add(hdipho_Invmass_qed,-1); 
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



  // "===================     Unfold data ======================     
  TCanvas *c2 = new TCanvas("c2","Dimuon pt ",1100,500);
  c2->Divide(2,1);
  c2->cd(1);
//      hRecoDataPt_xSec  = getXSecHist(hPt_data_bkgSub, "hRecoDataPt_xSec", 5, 0, 1, 1.64723);
  //    hRecoDataPt_xSec->GetEntries();
      //RooUnfoldBayes unfold_dataPt(&responsePt, hRecoDataPt_xSec, itr); // unfold data
      RooUnfoldBayes unfold_dataPt(&responsePt, hPt_data_bkgSub, itr); // unfold data
      RooUnfoldBayes unfold_dataRapidity(&responseRapidity, hRapidity_data_bkgSub, itr);
      RooUnfoldBayes unfold_dataInvmass(&responseInvmass, hInvmass_data_bkgSub, itr);
  
      RooUnfoldInvert unfold_dataInvmass_Invert(&responseInvmass, hInvmass_data_bkgSub);
      RooUnfoldInvert unfold_dataPt_Invert(&responsePt, hPt_data_bkgSub);  
      RooUnfoldInvert unfold_dataRapidity_Invert(&responseRapidity, hRapidity_data_bkgSub);
  
  TH1D* hUnfoldPt_data= (TH1D*) unfold_dataPt.Hunfold();
  unfold_dataPt.PrintTable (cout, hdipho_GenPt[0]);
 
  TH1D* hUnfoldRapidity_data= (TH1D*) unfold_dataRapidity.Hunfold();
  unfold_dataRapidity.PrintTable (cout, hdipho_GenRapidity[0]);
 
  TH1D* hUnfoldInvmass_data= (TH1D*) unfold_dataInvmass.Hunfold();
  unfold_dataInvmass.PrintTable (cout, hdipho_GenInvmass[0]);

  cout << "Invert1" << endl;
 
  TH1D* hUnfoldInvmass_Invert_data= (TH1D*) unfold_dataInvmass_Invert.Hunfold();
  unfold_dataInvmass_Invert.PrintTable (cout, hdipho_GenInvmass[0]);

  TH1D* hUnfoldPt_Invert_data= (TH1D*) unfold_dataPt_Invert.Hunfold();
  unfold_dataPt_Invert.PrintTable (cout, hdipho_GenPt[0]);
 
  TH1D* hUnfoldRapidity_Invert_data= (TH1D*) unfold_dataRapidity_Invert.Hunfold();
  unfold_dataRapidity_Invert.PrintTable (cout, hdipho_GenRapidity[0]);

  

  cout << "Invert2:" << endl;
  // divide histograms by luminosty and binwidth to get the cross-sections. 
  hGenPt_xSec  = getXSecHist(hdipho_GenPt[0], "hGenPt_xSec", 5, 0, 1 ,1 );
//  hGenRap_xSec = getXSecHist(hdipho_GenRapidity[0], "hGenRap_xSec",2,0,2.2, 1 );
  hGenRap_xSec = getRapXSecHist(hdipho_GenRapidity[0], "hGenRap_xSec", 1 );
  hGenInvmass_xSec = getInvXSecHist(hdipho_GenInvmass[0], "hGenInvmass_xSec" ,1 );

  hRecoMCPt_xSec  = getXSecHist(hdipho_Pt[0], "hRecoMCPt_xSec", 5, 0, 1,1 );
//  hRecoMCRap_xSec = getXSecHist(hdipho_Rapidity[0], "hRecoMCRap_xSec", 2,0,2.2,1 );
  hRecoMCRap_xSec = getRapXSecHist(hdipho_Rapidity[0], "hRecoMCRap_xSec", 1 );
  hRecoMCInvmass_xSec = getInvXSecHist(hdipho_Invmass[0], "hRecoMCInvmass_xSec",1 );

/*
  hRecoDataPt_xSec  = getXSecHist(hdipho_Pt_data, "hRecoDataPt_xSec", 10, 0, 2, 1.639);
  //hRecoDataRap_xSec = getXSecHist(hdipho_Rapidity_data, "hRecoDataRap_xSec",2,0,1.6, 1.639);
  hRecoDataRap_xSec = getRapXSecHist(hdipho_Rapidity_data, "hRecoDataRap_xSec", 1.639);
  hRecoDataInvmass_xSec = getInvXSecHist(hdipho_Invmass_data, "hRecoDataInvmass_xSec", 1.639);
*/

  hRecoDataPt_xSec  = getXSecHist(hPt_data_bkgSub, "hRecoDataPt_xSec", 5, 0, 1, 1.64723);
//  hRecoDataRap_xSec = getXSecHist(hRapidity_data_bkgSub, "hRecoDataRap_xSec",2,0,2.2, 1.64723);
  hRecoDataRap_xSec = getRapXSecHist(hRapidity_data_bkgSub, "hRecoDataRap_xSec", 1.64723);
  hRecoDataInvmass_xSec = getInvXSecHist(hInvmass_data_bkgSub, "hRecoDataInvmass_xSec", 1.64723);


  hUnfoMCPt_xSec  = getXSecHist(hUnfoldPt,       "hUnfoMCPt_xSec", 5, 0, 1,1 );
 // hUnfoMCRap_xSec = getXSecHist(hUnfoldRapidity, "hUnfoMCRap_xSec",2,0,2.2, 1 );
  hUnfoMCRap_xSec = getRapXSecHist(hUnfoldRapidity, "hUnfoMCRap_xSec", 1 );
  hUnfoMCInvmass_xSec = getInvXSecHist(hUnfoldInvmass, "hUnfoMCInvmass_xSec", 1 );
  hUnfoMCInvmass_Invert_xSec = getInvXSecHist(hUnfoldInvmass_Invert, "hUnfoMCInvmass_Invert_xSec", 1 );
  hUnfoMCPt_Invert_xSec  = getXSecHist(hUnfoldPt_Invert,       "hUnfoMCPt_Invert_xSec", 5, 0, 1,1 );
  hUnfoMCRap_Invert_xSec = getRapXSecHist(hUnfoldRapidity_Invert, "hUnfoMCRap_Invert_xSec", 1 ); 
 

  hUnfoDataPt_xSec  = getXSecHist(hUnfoldPt_data, "hUnfoDataPt_xSec", 5, 0, 1, 1.64723);
 // hUnfoDataRap_xSec = getXSecHist(hUnfoldRapidity_data, "hUnfoDataRap_xSec",2,0,2.2, 1.64723);
  hUnfoDataRap_xSec = getRapXSecHist(hUnfoldRapidity_data, "hUnfoDataRap_xSec",  1.64723);
  hUnfoDataInvmass_xSec = getInvXSecHist(hUnfoldInvmass_data, "hUnfoDataInvmass_xSec", 1.64723);
  hUnfoDataInvmass_Invert_xSec = getInvXSecHist(hUnfoldInvmass_Invert_data, "hUnfoDataInvmass_Invert_xSec", 1.64723);
   hUnfoDataPt_Invert_xSec  = getXSecHist(hUnfoldPt_Invert_data, "hUnfoDataPt_Invert_xSec", 5, 0, 1, 1.64723);
   hUnfoDataRap_Invert_xSec = getRapXSecHist(hUnfoldRapidity_Invert_data, "hUnfoDataRap_Invert_xSec",  1.64723);
 
 hUnfoDataRap_xSec->Draw("p");

 new TCanvas();
 hUnfoDataInvmass_Invert_xSec->Draw();
/*  new TCanvas();
 
 hdipho_GenPt[0]->SetLineColor(kGreen);
 hdipho_GenPt[0]->Draw("hist");
 hdipho_Pt_data->SetLineColor(kBlue);
 hdipho_Pt_data->Draw("histsame");

 hUnfoldPt_data->SetLineColor(kRed);
 hUnfoldPt_data->Draw("psame");


  new TCanvas();
 
 hdipho_GenRapidity[0]->SetLineColor(kGreen);
 hdipho_GenRapidity[0]->Draw("hist");
 hdipho_Rapidity_data->SetLineColor(kBlue);
 hdipho_Rapidity_data->Draw("histsame");

 hUnfoldRapidity_data->SetLineColor(kRed);
 hUnfoldRapidity_data->Draw("psame");

  new TCanvas();
 
 hdipho_GenInvmass[0]->SetLineColor(kGreen);
 hdipho_GenInvmass[0]->Draw("hist");
 hdipho_Invmass_data->SetLineColor(kBlue);
 hdipho_Invmass_data->Draw("histsame");

 hUnfoldInvmass_data->SetLineColor(kRed);
 hUnfoldInvmass_data->Draw("psame");*/
 /*******************/
/*  TCanvas*c11 = new TCanvas();
  hPt_data_bkgSub->SetLineColor(kOrange);
  hPt_data_bkgSub->SetMarkerColor(kOrange);
  hPt_data_bkgSub->Draw("p");
  hdipho_Pt[0]->Sumw2();
  hdipho_Pt[0]->SetMarkerColor(kBlack);
  hdipho_Pt[0]->Draw("psame");
  
  c11->SaveAs("Plot_7thNov/DiPho_Pt.png");
  
  TCanvas*c12 = new TCanvas();
  hRapidity_data_bkgSub->SetLineColor(kOrange);
  hRapidity_data_bkgSub->SetMarkerColor(kOrange);
  hRapidity_data_bkgSub->Draw("p");
  hdipho_Rapidity[0]->SetMarkerColor(kBlack);
  hdipho_Rapidity[0]->Draw("psame");
  c12->SaveAs("Plot_7thNov/DiPho_Rap.png");
  TCanvas*c13 = new TCanvas();
  hInvmass_data_bkgSub->SetLineColor(kOrange);
  hInvmass_data_bkgSub->SetMarkerColor(kOrange);
  hInvmass_data_bkgSub->Draw("p");
  hdipho_Invmass[0]->SetMarkerColor(kBlack);
  hdipho_Invmass[0]->Draw("psame");
  c13->SaveAs("Plot_7thNov/DiPho_Mass.png");
*/
/********************/
  cout << " invariant mass reco " << Sample[0] <<  " :" << hdipho_Invmass[0]->Integral() << endl;  
  cout << " invariant mass Gen " << Sample[0] <<  " :" << hdipho_GenInvmass[0]->Integral() << endl;  
  cout << " Rapidity reco " << Sample[0] <<  " :" << hdipho_Rapidity[0]->Integral() << endl;  
  cout << " Pt reco " << Sample[0] <<  " :" << hdipho_Pt[0]->Integral() << endl;  
  cout << " invariant mass data" << Sample[0] <<  " :" << hdipho_Invmass_data->Integral() << endl;
  cout << " invariant mass data w/o bkg" << Sample[0] <<  " :" << hInvmass_data_bkgSub->Integral() << endl;
  cout << " Rapidity data " << Sample[0] <<  " :" << hdipho_Rapidity_data->Integral() << endl;
  cout << " Rapidity data w/o bkg " << Sample[0] <<  " :" << hRapidity_data_bkgSub->Integral() << endl;
  cout << " Pt data " << Sample[0] <<  " :" << hdipho_Pt_data->Integral() << endl;
  cout << " Pt data w/o bkg " << Sample[0] <<  " :" << hPt_data_bkgSub->Integral() << endl;
  
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
  cout << "Cross section Rapidity UnfoldedData:" <<  hUnfoDataRap_xSec->Integral()<< endl;
  cout << "Cross section Invmass UnfoldedData:" << hUnfoDataInvmass_xSec->Integral()<< endl;
  cout << "Cross section Invmass UnfoldedData_Invert:" << hUnfoDataInvmass_Invert_xSec->Integral()<< endl;
  cout << "Cross section Pt UnfoldedData_Invert:" <<   hUnfoDataPt_Invert_xSec->Integral() << endl;
  cout << "Cross section Rapidity UnfoldedData_Invert:" <<  hUnfoDataRap_Invert_xSec->Integral()<< endl;
  
  outf->cd();
  //outf->Write();
  

///
  hGenPt_xSec->Write();
  hGenRap_xSec->Write();
  hGenInvmass_xSec->Write();

  hRecoMCPt_xSec->Write();
  hRecoMCRap_xSec->Write();
  hRecoMCInvmass_xSec->Write();

  hRecoDataPt_xSec->Write();
  hRecoDataRap_xSec->Write();
  hRecoDataInvmass_xSec->Write();
//  hRecoDataInvmass_Invert_xSec->Write();

  hUnfoMCPt_xSec->Write();
  hUnfoMCRap_xSec->Write();
  hUnfoMCInvmass_xSec->Write();
  hUnfoMCPt_Invert_xSec->Write();
  hUnfoMCRap_Invert_xSec->Write();
  hUnfoMCInvmass_Invert_xSec->Write();
  


  hUnfoDataPt_xSec->Write();
  hUnfoDataRap_xSec->Write();
  hUnfoDataInvmass_xSec->Write();
  hUnfoDataInvmass_Invert_xSec->Write();
  hUnfoDataPt_Invert_xSec->Write();
  hUnfoDataRap_Invert_xSec->Write();
  TCanvas*c1 =  new TCanvas();
  hdipho_Rapidity[0]->SetLineColor(kRed);
  hdipho_Rapidity[0]->Draw("p");;
   


  outf->Close();


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




TH1D* getXSecHist(TH1D *hist, const char* name, int nbins, float xmin, float xmax, float luminosity){

  TH1D *hans = new TH1D(name,"",nbins,xmin,xmax);
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hist->GetNbinsX(); i++) {

//     double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));
      double binCont  = hist->GetBinContent(i)/(luminosity);

      double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
      //double binError = hist->GetBinError(i);
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

      double binCont  = hist->GetBinContent(i)/(luminosity);
//      double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));

      double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
      //double binError = hist->GetBinError(i);
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

      double binCont  = hist->GetBinContent(i)/(luminosity);
//      double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));

      double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
      //double binError = hist->GetBinError(i);
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


