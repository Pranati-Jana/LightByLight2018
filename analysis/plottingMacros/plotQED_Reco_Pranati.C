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
#include "ReadQEDTree.C"
//#include "ReadGenEleTree.C"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "THStack.h"
#include "TCut.h"
#include "TAxis.h"
#include "TChain.h"
//#include "../eta2p2/CMS_lumi.C"
#include "MyCanvas.C"

//const double LumiLossHotZDCneg = 0.047433369; // 1034./21799 from QED number LbyL 5Aug 2022 slides
const double LumiLossHotZDCneg = 0;

const double luminosity       = 1647.228136; // μb^-1 from Gabi Feb 13, 2022

const int nSample = 3;
const char *Sample[nSample]={"Data", "QEDSCFSR_Private", "QEDSL"};

int ireg(double et, double eta);
void make_canvas(TCanvas *&);
//void make_canvas_ratio(TCanvas *&);
void make_hist(TH1D *&, Color_t , int );
//void make_hist_ratio(TH1D *&,Color_t, int) ;


const double nEventsGeneratedSC_Old = 59260000;//Superchic + photos
const double xsecGeneratedSC    = 8827.220; // μb





double scaleFactorsSC = 0.8497 * 0.9339 * pow(0.943, 2) * 1.0089 * 0.8696 ;

double lumiNormSC_Private = (xsecGeneratedSC*luminosity*scaleFactorsSC)/nEventsGeneratedSC_Old;
double lumiNormSL = (7920*luminosity*scaleFactorsSC)/66750000;
//double lumiNormSC_Private =1;
//double lumiNormSL =1;

float getAcoBinSF(float acop);
double SFReco(double et, double eta);
double SFReco_uncert(double et, double eta);
double SFTrig(double et, double eta);
double_t GetDeltaPhi(double dphi, double phi);
TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, TH1D*, TH1D*, TH1D*, double , double , double , double , double,double , const char *, bool); 
//TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, TH1D*, double , double , double , double , double, double, const char *, bool); 


TH1D* SFuncert(TTree *tr, const char* name, const char* var, const char* cut, int nbins, double binmin, double binmax);

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kMagenta);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}


void plotQED_Reco_Pranati(){
 

  const double wt[nSample] = {1,lumiNormSC_Private,lumiNormSL};

  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  TFile *outf= new TFile("output_gen.root","recreate");  

  TChain *qed[nSample];

  TH1D* hElePt1[nSample], *hEleEta1[nSample], *hElePhi1[nSample], *hElePt2[nSample], *hEleEta2[nSample], *hElePhi2[nSample]; 
  TH1D* hSumPt[nSample], *hRap[nSample], *hInvmass[nSample], *hAcop[nSample] ; TH1F *hSF[nSample], *hInvmass_err[nSample];
  TH1D* hCosThetaStar[nSample], *hdeltaPhi1[nSample];
 

  for (int i = 0; i < 3; i++){
    qed[i] = new TChain("output_tree");
    qed[i]->Add(Form("After_Dielectron/WithTrigger/%s_qed.root",Sample[i]));

    hElePt1[i]    = new TH1D(Form("hElePt1%s", Sample[i]),"",50,0,50);
    hEleEta1[i]   = new TH1D(Form("hEleEta1%s",Sample[i]),"",22,-2.2,2.2);
    hElePhi1[i]   = new TH1D(Form("hElePhi1%s",Sample[i]),"",16,-4,4);

    hElePt2[i]    = new TH1D(Form("hElePt2%s", Sample[i]),"",50,0,50);
    hEleEta2[i]   = new TH1D(Form("hEleEta2%s",Sample[i]),"",22,-2.2,2.2);
    hElePhi2[i]   = new TH1D(Form("hElePhi2%s",Sample[i]),"",16,-4,4);
 
    hSumPt[i]    = new TH1D(Form("hSumPt%s", Sample[i]),"",20,0,20);
    hRap[i]   = new TH1D(Form("hRap%s",Sample[i]),"",22,-2.2,2.2);
    hInvmass[i]   = new TH1D(Form("hInvmass%s",Sample[i]),"",50,0,100);
    hAcop[i]   = new TH1D(Form("hAcop%s",Sample[i]),"",100,0,0.1);
    hdeltaPhi1[i]  = new TH1D(Form("hdeltaPhi1%s",Sample[i]), "", 21, -1, 4);
     hCosThetaStar[i]   = new TH1D(Form("hCosThetaStar%s",Sample[i]),"",20,0,1); 
    cout << "file " << qed[i]->GetEntries()  << endl;
    ReadQEDTree  qedR(qed[i]);
    qedR.fChain->SetBranchStatus("*",1);
    
    //cout << "QEDSCFSR norm Official:" << lumiNormSC_Official << endl;
   // cout << "QEDSCFSR norm Private:" << lumiNormSC_Private << endl;
    if (qedR.fChain == 0) return;
  
    Long64_t nentries = qedR.fChain->GetEntriesFast();
    cout << qed[i]->GetName() << "    " << nentries << endl;
  
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry_evt = qedR.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = qedR.fChain->GetEntry(jentry);   nbytes += nb;

     // if(qedR.ok_trigger !=1) continue;
      if(qedR.vSum_M < 5) continue; //invmass
      if(qedR.vSum_Pt > 1) continue; //diele pt
      if(abs(qedR.eleEta_1)  > 2.2) continue;
      if(abs(qedR.eleEta_2) > 2.2) continue;
      if(qedR.elePt_1 < 2.0) continue;
      if(qedR.elePt_2 < 2.0) continue;

      if(qedR.ok_neuexcl != 1) continue; 
       if(qedR.ok_chexcl_extrk != 1) continue;       

      if(i==0){
       //  if(qedR.zdc_energy_pos > 10000 || qedR.zdc_energy_neg > 10000 ) continue;
//	 if(qedR.zdc_energy_neg > 10000 ){
////	 cout << "" << qedR.zdc_energy_neg << endl;
          }
	 //cout << "Test" << endl;
      //   }


       hAcop[i]->Fill(qedR.ele_acop,wt[i]);

      if(qedR.ele_acop > 0.01) continue; //acop
	hElePt1[i]->Fill(qedR.elePt_1,wt[i]);
	hEleEta1[i]->Fill(qedR.eleEta_1,wt[i]);
	hElePhi1[i]->Fill(qedR.elePhi_1,wt[i]);
	
	hElePt2[i]->Fill(qedR.elePt_2,wt[i]);
	hEleEta2[i]->Fill(qedR.eleEta_2,wt[i]);
	hElePhi2[i]->Fill(qedR.elePhi_2,wt[i]);
	
	hSumPt[i]->Fill(qedR.vSum_Pt,wt[i]);
	hRap[i]->Fill(qedR.vSum_Rapidity,wt[i]);
	hInvmass[i]->Fill(qedR.vSum_M,wt[i]);
        hCosThetaStar[i]->Fill(abs(qedR.costhetastar),wt[i]);
//	hAcop[i]->Fill(qedR.ele_acop,wt[i]);
	//if( fabs(qedR.eleEta_1) < 1.5 && fabs(qedR.eleEta_2) < 1.5) {
//	if(abs(qedR.eleEta_1) < 1.5) continue;
 //       if(abs(qedR.eleEta_2) < 1.5) continue;
/*	if(qedR.eleCharge_1 == 1){
	hdeltaPhi1[i]->Fill(((qedR.ele_dphi)-(qedR.elePhi_1)));

//	hdeltaPhi1[i]->Fill((qedR.ele_dphi));
	
        }
	else{
	hdeltaPhi1[i]->Fill(((qedR.ele_dphi)-(qedR.elePhi_2)));
//	hdeltaPhi1[i]->Fill((qedR.ele_dphi-(qedR.elePhi_2)));
	//hdeltaPhi1[i]->Fill((qedR.ele_dphi));
	}
	}*/

         if(qedR.eleCharge_1 == 1){
//        hCosdphi2[i]->Fill(((qedR.ele_dphi) - (qedR.elePhi_1)));
     //    cout << "ele charge1:" << qedR.eleCharge_1 << endl;
     //   cout << "dphi:" << qedR.vSum_Phi << endl;
        double dPhi = (qedR.ele_dphi);
        double phi = (qedR.elePhi_1);
        double cosphi1 = GetDeltaPhi(dPhi , phi);

         hdeltaPhi1[i]->Fill(cosphi1);

       // cout << "cosdphi1:" << cosphi1 << endl;
        }
        else{
        // cout << "ele charge2:" << qedR.eleCharge_2 << endl;
        // cout << "dphi:" << qedR.ele_dphi << endl;
         double dPhi = (qedR.ele_dphi);
        double phi = (qedR.elePhi_2);
        double cosphi1 = GetDeltaPhi(dPhi , phi);
         hdeltaPhi1[i]->Fill(cosphi1);
      //  cout << "cosdphi2:" << cosphi1 << endl;
        }


    } //entry
  } // for 4 files

  for (int i = 0; i < 3; i++){    
    hElePt1[i]->Add(hElePt2[i]);
    hEleEta1[i]->Add(hEleEta2[i]);
    hElePhi1[i]->Add(hElePhi2[i]);

  }
  

   int W = 700;
   int H = 600;

   float T = 0.08;
   float B = 0.14;
   float L = 0.14;
   float R = 0.04;
   //
   /////


   MyCanvas mc1("mass_Private","Dielectron invariant mass (GeV)", "Entries", W, H);
   mc1.SetLogy(false);
   mc1.SetYRange(0.1,10000);
   mc1.SetRatioRange(0.1,1.9);
   mc1.SetLegendPosition(0.60,0.68,0.9,0.85);
   mc1.CanvasWithThreeHistogramsRatioPlot(hInvmass[0],hInvmass[1],hInvmass[2],"Data","Superchic+Photos","Starlight","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
  //mc1.CanvasWithThreeHistogramsRatioPlot(hInvmass[0],hInvmass[1],"Data","Superchic+Photos","Starlight","Data/MC",kBlack,kBlue,kFALSE,kTRUE,"EP","hist SAME");
   mc1.PrintCanvas();

   MyCanvas mc2("pt_Private","Dielectron p_{T} (GeV)", "Entries", W, H);
   mc2.SetLogy(false);
   mc2.SetYRange(0.1,6000);
   mc2.SetXRange(0,1);
   mc2.SetRatioRange(0.1,1.9);
   mc2.SetLegendPosition(0.6,0.68,0.9,0.85);
   mc2.CanvasWithThreeHistogramsRatioPlot(hSumPt[0],hSumPt[1],hSumPt[2],"Data","Superchic+photos ","Starlight","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   mc2.PrintCanvas();

   MyCanvas mc3("rap_Private","Dielectron y", "Entries", W, H);
   mc3.SetYRange(0,5000);
   //mc3.SetRatioRange(0.1,1.9);
   mc3.SetLegendPosition(0.16,0.68,0.56,0.85);
   mc3.CanvasWithThreeHistogramsRatioPlot(hRap[0],hRap[1],hRap[2],"Data","Superchic+Photos","Starlight","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   mc3.PrintCanvas();

    MyCanvas mc4("elePt_Private","Electron Pt", "Entries", W, H);
   mc4.SetLogy(false);
   mc4.SetYRange(0.1,50000);
   mc4.SetRatioRange(0.1,1.9);
   mc4.SetLegendPosition(0.6,0.68,0.9,0.85);
   mc4.CanvasWithThreeHistogramsRatioPlot(hElePt1[0],hElePt1[1],hElePt1[2],"Data","Superchic+Photos","Starlight","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   mc4.PrintCanvas();


   MyCanvas mc5("eleEta_Private","Electron Eta", "Entries", W, H);
   mc5.SetYRange(0,5000);
   mc5.SetRatioRange(0.1,1.9);
   mc5.SetLegendPosition(0.16,0.68,0.56,0.85);
   mc5.CanvasWithThreeHistogramsRatioPlot(hEleEta1[0],hEleEta1[1],hEleEta1[2],"Data","Superchic+Photos","Starlight","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   mc5.PrintCanvas();


   MyCanvas mc6("elePhi_Private","Electron Phi", "Entries", W, H);
   mc6.SetYRange(0,5000);
   mc6.SetRatioRange(0.1,1.9);
   mc6.SetLegendPosition(0.16,0.68,0.56,0.85);
   mc6.CanvasWithThreeHistogramsRatioPlot(hElePhi1[0],hElePhi1[1],hElePhi1[2],"Data","Superchic+Photos","Starlight","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   mc6.PrintCanvas();

   MyCanvas mc7("acop_Private","Dielectron acoplanarity", "Entries", W, H);
   mc7.SetLogy(false);
   mc7.SetXRange(0,0.06);
   mc7.SetYRange(0.1,100000);
   mc7.SetRatioRange(0.1,1.9);
   mc7.SetLegendPosition(0.6,0.68,0.9,0.85);
   mc7.CanvasWithThreeHistogramsRatioPlot(hAcop[0],hAcop[1],hAcop[2],"Data","Superchic+Photos","Starlight", "Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   mc7.PrintCanvas();


   MyCanvas mc8("costhetastar_Private","|Cos#theta^{*}|", "Entries", W, H);
   //mc8.SetLogy(true);
   mc8.SetXRange(0,1.0);
   mc8.SetYRange(0,2000);
   mc8.SetRatioRange(0.1,1.9);
   mc8.SetLegendPosition(0.6,0.68,0.9,0.85);
   mc8.CanvasWithThreeHistogramsRatioPlot(hCosThetaStar[0],hCosThetaStar[1],hCosThetaStar[2],"Data","Superchic+Photos","Starlight", "Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   mc8.PrintCanvas();

  for (int i = 0; i < 3; i++){
  cout << " Acop < 0.01 " << Sample[i] <<  " :" << hAcop[i]->Integral(0,10) << endl;  
  cout << "  pt: " << Sample[i] <<  " :" << hSumPt[i]->Integral(0,50) << endl;  
 }

  
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
 // return dphi;


}

