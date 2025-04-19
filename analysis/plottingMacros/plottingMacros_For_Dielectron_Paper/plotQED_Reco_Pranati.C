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
//#include "ReadQEDTree.C"
#include "ReadQEDHFTree.C"
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

//#include "RecoID_SF_Electron.h"

//const double LumiLossHotZDCneg = 0.047433369; // 1034./21799 from QED number LbyL 5Aug 2022 slides
const double LumiLossHotZDCneg = 0;

//const double luminosity       = 1647.2; // μb^-1 from Gabi Feb 13, 2022
const double luminosity       = 1700.0;
//const double luminosity       = 1500.00;
const int nSample = 3;
//const char *Sample[nSample]={"Data", "QEDSCFSR_Private", "QEDSL"};
const char *Sample[nSample]={"CutFlow_BW_Data_", "CutFlow_BW_mc_photos_sc_dielectron", "CutFlow_BW_mc_qed_sl"};

int ireg(double et, double eta);
void make_canvas(TCanvas *&);
//void make_canvas_ratio(TCanvas *&);
void make_hist(TH1D *&, Color_t , int );
//void make_hist_ratio(TH1D *&,Color_t, int) ;


const double nEventsGeneratedSC_Old = 59260000;//Superchic + photos
const double xsecGeneratedSC    = 8827.220; // μb


double scaleFactorsSC =  1.0089 * 0.8716*0.9252 *0.952753;//0.8487  
//double scaleFactorsSC =  1.0089 * 0.8716;

double lumiNormSC_Other = (xsecGeneratedSC*luminosity*scaleFactorsSC)/nEventsGeneratedSC_Old;
double lumiNormSL_Other = (7920*luminosity*scaleFactorsSC)/66750000;
double lumiNormGUPC_Other = (13450*luminosity*scaleFactorsSC)/100000000;//gUPC_1FSR

float getAcoBinSF(float acop);
double SFReco(double et, double eta);
double SFReco_uncert(double et, double eta);
double SFTrig(double et, double eta);
double_t GetDeltaPhi(double dphi, double phi);
Double_t RecoID_SF_Electron( double ETT);


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

  const int nMassbins=22;
 //double Massbin[nMassbins]={4,6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,48,54,60,68,78,100};
 double Massbin[nMassbins]={5,7,9,11,13,15,17,19,21,23,25,27,31,35,39,43,49,55,61,69,79,100};
//  const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;

//David requested for thinner binning 2 times
//  const int nMassbins=43;

//  double Massbin[nMassbins]={5,7,9,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,31,33,35,37,39,41,43,46,49,52,55,58,61,65,69,74,79,90,100};
 
 const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;


const int nPtbins=12;
double Ptbin[nPtbins]={2,4,6,8,10,12,14,18,22,30,40,50};
const int nPtbin= sizeof(Ptbin)/sizeof(double) - 1;
//ZDC bin
const int nZDCbins=6;
double ZDCbin[nZDCbins]={0,2500,5000,7500,10000,30000};
const int nZDCbin= sizeof(ZDCbin)/sizeof(double) - 1;
//
const int nZDCbins1=6;
double ZDCbin1[nZDCbins1]={0,2000,4000,7000,10000,30000};
const int nZDCbin1= sizeof(ZDCbin1)/sizeof(double) - 1;

void plotQED_Reco_Pranati(){
//  const double wt[nSample] = {1,lumiNormSC_Private,lumiNormSL};
cout << "SF_Starlight:" << (7920*1647.2)/66750000 << endl;
cout << "SF_SC:" << 8827.220*1647.2/59260000 << endl;
  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
   gStyle->SetTextFont(42);
  TFile *outf= new TFile("output_dielectron_mass_43bin.root","recreate");
  TFile *outf1= new TFile("Data_mass_43bin.root","recreate"); 


  TChain *qed[nSample];

  TH1D* hElePt1[nSample], *hEleEta1[nSample], *hElePhi1[nSample], *hElePt2[nSample], *hEleEta2[nSample], *hElePhi2[nSample]; 
  TH1D* hSumPt[nSample], *hRap[nSample], *hInvmass[nSample], *hAcop[nSample] ; TH1F *hSF[nSample], *hInvmass_err[nSample];
  TH1D* hCosThetaStar[nSample], *hdeltaPhi1[nSample],  *hInvmass_norm[nSample], *hZDC[nSample], *hZDC_sum[nSample], *hZDC_sum_peak[nSample], *hZDC_sum_max[nSample], *hZDC_sum_peak_max[nSample];
  TH1D* hAco_WOZDC[nSample], *hAco_ZDCAND4n[nSample];
  TH1D* hInvmass_unfold[nSample], *hSumPt_unfold[nSample], *hRap_unfold[nSample], *hCosthetastar_unfold[nSample];
  TH1D* hInvmass_ZDC[nSample];
  for (int i = 0; i < 3; i++){
    
    qed[i] = new TChain("output_tree");
   // qed[i]->Add(Form("After_Dielectron/WithTrigger/%s_qed.root",Sample[i]));
    //qed[i]->Add(Form("After_Dielectron/SameSignElectrons/%s_qed.root",Sample[i]));
    qed[i]->Add(Form("After_Dielectron/WithTrigger/%s_HF+9.1HF-8.8_Eta3.root",Sample[i]));

    hElePt1[i]    = new TH1D(Form("hElePt1%s", Sample[i]),"",nPtbin,Ptbin);
    hEleEta1[i]   = new TH1D(Form("hEleEta1%s",Sample[i]),"",23,-2.2,2.2);
    hElePhi1[i]   = new TH1D(Form("hElePhi1%s",Sample[i]),"",14,-3.5,3.5);

    hElePt2[i]    = new TH1D(Form("hElePt2%s", Sample[i]),"",nPtbin,Ptbin);
    hEleEta2[i]   = new TH1D(Form("hEleEta2%s",Sample[i]),"",23,-2.2,2.2);
    hElePhi2[i]   = new TH1D(Form("hElePhi2%s",Sample[i]),"",14,-3.5,3.5);
 
    hSumPt[i]    = new TH1D(Form("hSumPt%s", Sample[i]),"",20,0,1);
    hRap[i]   = new TH1D(Form("hRap%s",Sample[i]),"",11,-2.2,2.2);
    hInvmass[i]   = new TH1D(Form("hInvmass%s",Sample[i]),"",nMassbins-1, Massbin);
    hAcop[i]   = new TH1D(Form("hAcop%s",Sample[i]),"",14,0.0,0.014);//0.014
    hdeltaPhi1[i]  = new TH1D(Form("hdeltaPhi1%s",Sample[i]), "", 21, -1, 4);
    hCosThetaStar[i]   = new TH1D(Form("hCosThetaStar%s",Sample[i]),"",20,0,1);
    //Unfolding plots 
    hInvmass_unfold[i]  = new TH1D(Form("hInvmass_unfold%s",Sample[i]),"",nMassbin,Massbin);
    hSumPt_unfold[i]  = new TH1D(Form("hSumPt_unfold%s", Sample[i]),"",7,0,1);
    hRap_unfold[i]   = new TH1D(Form("hRap_unfold%s",Sample[i]),"",11,-2.2,2.2);
    hCosthetastar_unfold[i] = new TH1D(Form("hCosthetastar_unfold%s",Sample[i]),"",10,0,1);
    //ZDC binning
     hZDC[i]   = new TH1D(Form("hZDC%s",Sample[i]),"",80,0,200000);
     hZDC_sum[i] = new TH1D(Form("hZDC_sum%s",Sample[i]),"",nZDCbin, ZDCbin);
     hZDC_sum_peak[i] = new TH1D(Form("hZDC_sum_peak%s",Sample[i]),"",nZDCbin1, ZDCbin1);
     hZDC_sum_max[i] = new TH1D(Form("hZDC_sum_max%s",Sample[i]),"",nZDCbin, ZDCbin);
     hZDC_sum_peak_max[i] = new TH1D(Form("hZDC_sum_peak_max%s",Sample[i]),"",nZDCbin1, ZDCbin1);
     hAco_WOZDC[i] = new TH1D(Form("hAco_WOZDC%s",Sample[i]),"",100,0,0.5);
     hAco_ZDCAND4n[i] = new TH1D(Form("hAco_ZDCAND4n%s",Sample[i]),"",100,0,0.1);
     hInvmass_ZDC[i] =  new TH1D(Form("hInvmass_ZDC%s",Sample[i]),"",nMassbin,Massbin);
    
    int zeroNzeroN_Passed = 0, zeroNxN_Passed = 0, zeroNoneN_Passed = 0, oneNoneN_Passed = 0,zeroNtwoN_Passed = 0,twoNtwoN_Passed =0, zeroNthreeN_Passed = 0, threeNthreeN_Passed = 0, zeroNfourN_Passed = 0,fourNfourN_Passed = 0,fourNxN_Passed = 0 ; 
  
    
    cout << "file " << qed[i]->GetEntries()  << endl;
   // ReadQEDTree  qedR(qed[i]);
     ReadQEDHFTree  qedR(qed[i]);
    qedR.fChain->SetBranchStatus("*",1);
    
    if (qedR.fChain == 0) return;
  
    Long64_t nentries = qedR.fChain->GetEntriesFast();
    cout << qed[i]->GetName() << "    " << nentries << endl;
  
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry_evt = qedR.LoadTree(jentry);
      if (ientry_evt <0) break;
      nb = qedR.fChain->GetEntry(jentry);   nbytes += nb;

    //  if(qedR.ok_trigger !=1) continue;
       if(abs(qedR.eleEta_1)  > 2.2) continue;
       if(abs(qedR.eleEta_2) > 2.2) continue;
       if(qedR.elePt_1 < 2.0) continue;
       if(qedR.elePt_2 < 2.0) continue;

       if(qedR.vSum_M < 5) continue; //invmass
   
       if(qedR.ok_chexcl_extrk != 1) continue;  
      if(qedR.ok_neuexcl != 1) continue; 
       if(qedR.ok_HFNeeExcl!=1 ) continue;
          if(qedR.vSum_Pt > 1) continue; //diele pt
          
   
    

      if(i==0){
        if(qedR.zdc_energy_pos > 7000 &&  qedR.zdc_energy_neg > 7000 ) continue;
             
        
      }
    double  ET1 = qedR.eleEta_1;
    double  ET2 = qedR.eleEta_2;
     double lumiNorm_SC = RecoID_SF_Electron(ET1) * RecoID_SF_Electron(ET2) * lumiNormSC_Other;
     double lumiNorm_SL = RecoID_SF_Electron(ET1) * RecoID_SF_Electron(ET2) * lumiNormSL_Other;
    // double lumiNorm_GUPC = RecoID_SF_Electron(ET1) * RecoID_SF_Electron(ET2) * lumiNormGUPC_Other;
    // double lumiNorm_SC =  lumiNormSC_Other;
    // double lumiNorm_SL = lumiNormSL_Other;
    double lumiNorm_GUPC = RecoID_SF_Electron(ET1) * RecoID_SF_Electron(ET2) * lumiNormGUPC_Other;
/*
      cout << "ET1:" << ET1 << "\t" << "\tSF1:" << RecoID_SF_Electron(ET1) << endl;
      cout << "ET2:" << ET2 << "\t" << "\tSF2:" << RecoID_SF_Electron(ET2) << endl;
      cout << "lumiNorm_Other:" << lumiNormSC_Other << endl;
      cout << "lumiNorm_SC:" << lumiNorm_SC << endl;
  */
      std::vector<double> wt = {1, lumiNorm_SC, lumiNorm_SL};
     //std::vector<double> wt = {1, 1, 1};

  hAco_ZDCAND4n[i]->Fill(qedR.ele_acop,wt[i]);

  hAcop[i]->Fill(qedR.ele_acop,wt[i]);

  if(qedR.ele_acop > 0.01) continue; //acop
  
  float binwidth1 = hElePt1[i]->GetBinWidth(hElePt1[i]->GetXaxis()->FindFixBin(qedR.elePt_1));
  //float binwidth1 =1 ;
	hElePt1[i]->Fill(qedR.elePt_1,wt[i]/binwidth1);
	hEleEta1[i]->Fill(qedR.eleEta_1,wt[i]);
	hElePhi1[i]->Fill(qedR.elePhi_1,wt[i]);
	float binwidth2 = hElePt2[i]->GetBinWidth(hElePt2[i]->GetXaxis()->FindFixBin(qedR.elePt_2));
  //float binwidth2 = 1;
	hElePt2[i]->Fill(qedR.elePt_2,wt[i]/binwidth2);
	hEleEta2[i]->Fill(qedR.eleEta_2,wt[i]);
	hElePhi2[i]->Fill(qedR.elePhi_2,wt[i]);
	
	hSumPt[i]->Fill(qedR.vSum_Pt,wt[i]);
	hRap[i]->Fill(qedR.vSum_Rapidity,wt[i]);
  
  float binwidth3 = hInvmass[i]->GetBinWidth(hInvmass[i]->GetXaxis()->FindFixBin(qedR.vSum_M));
  //float binwidth3 = 1;
	hInvmass[i]->Fill(qedR.vSum_M,wt[i]/binwidth3);
  hCosThetaStar[i]->Fill(abs(qedR.costhetastar),wt[i]);
  //unfodling plots
  hInvmass_unfold[i]->Fill(qedR.vSum_M,wt[i]);
  hSumPt_unfold[i]->Fill(qedR.vSum_Pt,wt[i]);
  hRap_unfold[i]->Fill(qedR.vSum_Rapidity,wt[i]);
  hCosthetastar_unfold[i]->Fill(abs(qedR.costhetastar),wt[i]);
  hInvmass_ZDC[i]->Fill(qedR.vSum_M,wt[i]);

         if(qedR.eleCharge_1 == 1){

        double dPhi = (qedR.ele_dphi);
        double phi = (qedR.elePhi_1);
        double cosphi1 = GetDeltaPhi(dPhi , phi);
         hdeltaPhi1[i]->Fill(cosphi1);
        }
        else{ 
        double dPhi = (qedR.ele_dphi);
        double phi = (qedR.elePhi_2);
        double cosphi1 = GetDeltaPhi(dPhi , phi);
         hdeltaPhi1[i]->Fill(cosphi1);
      //  cout << "cosdphi2:" << cosphi1 << endl
        }
     /************************** ZDC************/
     if(i==0){
      hZDC[i]->Fill(qedR.zdc_energy_pos);
      
      if(qedR.zdc_energy_pos + qedR.zdc_energy_neg > 20000){
       hZDC_sum[0]->Fill(20000);
      }
      else{
      hZDC_sum[0]->Fill(qedR.zdc_energy_pos + qedR.zdc_energy_neg);

      } 
      if(qedR.zdc_energy_pos + qedR.zdc_energy_neg > 20000){
       hZDC_sum_peak[0]->Fill(20000);
      }
      else{
       hZDC_sum_peak[0]->Fill(qedR.zdc_energy_pos + qedR.zdc_energy_neg);
      }
      if(qedR.zdc_energy_pos >qedR.zdc_energy_neg){

        if(qedR.zdc_energy_pos>20000){
          hZDC_sum_max[0]->Fill(20000);
          hZDC_sum_peak_max[0]->Fill(20000);
        }
        else{
        hZDC_sum_max[0]->Fill(qedR.zdc_energy_pos);
        hZDC_sum_peak_max[0]->Fill(qedR.zdc_energy_pos);
        }
      }
      else{
        if(qedR.zdc_energy_neg>20000)
        {
        hZDC_sum_max[0]->Fill(20000);
        hZDC_sum_peak_max[0]->Fill(20000);

        }
        else{
        hZDC_sum_max[0]->Fill(qedR.zdc_energy_neg);
        hZDC_sum_peak_max[0]->Fill(qedR.zdc_energy_neg);
        }
      }
     }
///////////////
     
    } //entry
    
  } // for 4 files
//cout << "Starlight: " << hInvmass[2]->GetBinContent(22) << endl;

  for (int i = 0; i < 3; i++){    
    hElePt1[i]->Add(hElePt2[i]);
    hEleEta1[i]->Add(hEleEta2[i]);
    hElePhi1[i]->Add(hElePhi2[i]);
  }
  //Normalization of inv mass and electron pt
  for(int i = 0; i < 3; i++){
   int nBins = hInvmass[i]->GetNbinsX();
   for (int bin = 1; bin <= nBins; ++bin) {
    double binContent = hInvmass[i]->GetBinContent(bin);
    double binWidth = hInvmass[i]->GetBinWidth(bin);
   // hInvmass[i]->SetBinContent(bin, binContent / binWidth);
   }
   int nBinpt = hElePt1[i]->GetNbinsX();
   for (int bin = 1; bin <= nBinpt; ++bin) {
    double binContentpt = hElePt1[i]->GetBinContent(bin);
    double binWidthpt = hElePt1[i]->GetBinWidth(bin);
    //hElePt1[i]->SetBinContent(bin, binContentpt / binWidthpt);
   }

  }
  
  TCanvas* c = new TCanvas();
 hZDC[0]->Draw("hist");
 c->Print("ZDC.png");
TCanvas* c1 = new TCanvas();
 hInvmass_ZDC[0]->Draw("hist");
  hInvmass_ZDC[0]->SetXTitle( "Invariant mass" );
// c1->Print("InvmassZDC.png");


   int W = 800;
   int H = 600;

   float T = 0.08;
   float B = 0.14;
   float L = 0.14;
   float R = 0.04;
   //
   /////
  
  gStyle->SetTextFont(42);
   MyCanvas mc1("Figure_002_g","m^{ee} (GeV)", "#LT Events / GeV #GT ", W, H);
   mc1.SetLogy(false);
   mc1.SetYRange(0.02,10000);//10000
   mc1.SetRatioRange(0.1,1.9);
   mc1.SetLegendPosition(0.33,0.62,0.9,0.83);
   mc1.CanvasWithThreeHistogramsRatioPlot(hInvmass[0],hInvmass[1],hInvmass[2],"Data","SUPERCHIC 3.03 + FSR (PHOTOS++)","STARLIGHT 3.13","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
  //mc1.CanvasWithThreeHistogramsRatioPlot(hInvmass[0],hInvmass[1],"Data","Superchic+Photos","Starlight","Data/MC",kBlack,kBlue,kFALSE,kTRUE,"EP","hist SAME");
  // mc1.PrintCanvas();

   MyCanvas mc2("Figure_002_b","p^{ee}_{T} (GeV)", "Events / 0.05 (GeV^{-1}) ", W, H);
   mc2.SetLogy(false);
   mc2.SetYRange(50,6000);//(50,6000)
   mc2.SetXRange(0,1);
   mc2.SetRatioRange(0.1,1.9);
   mc2.SetLegendPosition(0.36,0.62,0.85,0.82);
   mc2.CanvasWithThreeHistogramsRatioPlot(hSumPt[0],hSumPt[1],hSumPt[2],"Data","SUPERCHIC 3.03 + FSR (PHOTOS++) ","STARLIGHT 3.13","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   //mc2.PrintCanvas();

   MyCanvas mc3("Figure_002_d"," y^{ee}", "Events / 0.4", W, H);
   mc3.SetYRange(0,7000);//(0,7000)
   mc3.SetRatioRange(0.6,1.19);
   mc3.SetLegendPosition(0.20,0.60,0.56,0.84);
   mc3.CanvasWithThreeHistogramsRatioPlot(hRap[0],hRap[1],hRap[2],"Data","SUPERCHIC 3.03 + FSR (PHOTOS++)","STARLIGHT 3.13","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   //mc3.PrintCanvas();

  MyCanvas mc4("Figure_002_a","p^{e}_{T} (GeV)", " #LT Entries / GeV #GT ", W, H);

   mc4.SetLogy(false);
   mc4.SetYRange(0.08,50000);//(0.08,50000)
   mc4.SetRatioRange(0.1,1.9);
   mc4.SetLegendPosition(0.35,0.58,0.85,0.80);
   mc4.CanvasWithThreeHistogramsRatioPlot(hElePt1[0],hElePt1[1],hElePt1[2],"Data","SUPERCHIC 3.03 + FSR (PHOTOS++)","STARLIGHT 3.13","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   //mc4.PrintCanvas();


   MyCanvas mc5("Figure_002_c","#eta^{e}", "Entries / 0.19", W, H);
   mc5.SetYRange(0,5000);//(0,5000)
   mc5.SetRatioRange(0.6,1.35);
   mc5.SetLegendPosition(0.18,0.65,0.56,0.84);
   mc5.CanvasWithThreeHistogramsRatioPlot(hEleEta1[0],hEleEta1[1],hEleEta1[2],"Data","SUPERCHIC 3.03 + FSR (PHOTOS++)","STARLIGHT 3.13","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
  // mc5.PrintCanvas();


   MyCanvas mc6("Figure_002_e","#phi^{e}", "Entries / 0.5", W, H);
   mc6.SetYRange(0,5500);//(0,5500)
   mc6.SetRatioRange(0.85,1.19);
   mc6.SetLegendPosition(0.26,0.62,0.56,0.83);
   mc6.CanvasWithThreeHistogramsRatioPlot(hElePhi1[0],hElePhi1[1],hElePhi1[2],"Data","SUPERCHIC 3.03 + FSR (PHOTOS++)","STARLIGHT 3.13","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
  // mc6.PrintCanvas();

   MyCanvas mc7("Figure_002_f","A_{#phi}^{ee}", "Events / 0.001", W, H);
   mc7.SetLogy(false);
   mc7.SetXRange(0.0,0.014);
   mc7.SetYRange(0.12,100000);//(0.12,100000)
   //mc7.SetYRange(0.0,10);
   mc7.SetRatioRange(0.1,9.8);
   mc7.SetLegendPosition(0.35,0.62,0.9,0.83);
   mc7.CanvasWithThreeHistogramsRatioPlot(hAcop[0],hAcop[1],hAcop[2],"Data","SUPERCHIC 3.03 + FSR (PHOTOS++)","STARLIGHT 3.13", "Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   mc7.PrintCanvas();
   //mc7->SaveAs("ZDCAND2n.C");


   MyCanvas mc8("Figure_002_h","|cos#theta*|^{ee}", "Events / 0.1", W, H);
   //mc8.SetLogy(true);
   mc8.SetXRange(0,1.0);
   mc8.SetYRange(0,2200);//(0,2200)
   mc8.SetRatioRange(0.5,1.19);
   mc8.SetLegendPosition(0.35,0.62,0.9,0.83);
   mc8.CanvasWithThreeHistogramsRatioPlot(hCosThetaStar[0],hCosThetaStar[1],hCosThetaStar[2],"Data","SUPERCHIC 3.03 + FSR (PHOTOS++)","STARLIGHT 3.13", "Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   //mc8.PrintCanvas();
  
  MyCanvas mc9("mass_Private_unfold","Dielectron invariant mass (GeV)", "Events", W, H);
   mc9.SetLogy(false);
   mc9.SetYRange(0.02,10000);
   mc9.SetRatioRange(0.1,1.9);
   mc9.SetLegendPosition(0.55,0.7,0.9,0.85);
   mc9.CanvasWithThreeHistogramsRatioPlot(hInvmass_unfold[0],hInvmass_unfold[1],hInvmass_unfold[2],"Data","SUPERCHIS 3.03 + PHOTOS++","Starlight 3.13 (LO)","Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
  //mc1.CanvasWithThreeHistogramsRatioPlot(hInvmass[0],hInvmass[1],"Data","Superchic+Photos","Starlight","Data/MC",kBlack,kBlue,kFALSE,kTRUE,"EP","hist SAME");
   //mc9.PrintCanvas();

   MyCanvas mc10("acop_UptoNeutralExclusivity","Dielectron acoplanarity", "Event", W, H);
   //mc10.SetLogy(false);
   mc10.SetXRange(0.2,1.0);
   mc10.SetYRange(0.0,10);
   mc10.SetRatioRange(0.1,5);
   mc10.SetLegendPosition(0.6,0.65,0.9,0.85);
   mc10.CanvasWithThreeHistogramsRatioPlot(hAcop[0],hAcop[1],hAcop[2],"Data","SUPERCHIS+PHOTOS++","Starlight 3.13", "Data/MC",kBlack,kBlue,kRed,kFALSE,kTRUE,kTRUE,"EP","hist SAME", "hist SAME");
   //mc10.PrintCanvas();

  for (int i = 0; i < 3; i++){
  cout << " Acop < 0.01 " << Sample[i] <<  " :" << hAcop[i]->Integral(0,10) << endl;  
  cout << "  pt: " << Sample[i] <<  " :" << hSumPt[i]->Integral(0,20) << endl; 
  cout << "Mass:" << Sample[i] <<  " :" << hInvmass[i]->Integral("width") << endl; 
  cout << "Costhetastar:" << Sample[i] <<  " :" << hCosThetaStar[i]->Integral(0,20) << endl;
  //cout << " Acop < 0.01 " << Sample[i] <<  " :" << hAco_WOZDC[i]->Integral(0,10) << endl;
  //cout << " Acop < 0.01 " << Sample[i] <<  " :" << hAco_ZDCAND4n[i]->Integral(0,10) << endl;  
  //cout << "  pt unfold: " << Sample[i] <<  " :" << hSumPt_unfold[i]->Integral(0,40) << endl; 
 // cout << "Mass unfold:" << Sample[i] <<  " :" << hInvmass_unfold[i]->Integral(0,25) << endl; 
//  cout << "Rap unfold:" << Sample[i] <<  " :" << hRap_unfold[i]->Integral(0,23) << endl; 
 //cout << "Costhetastar unfold:" << Sample[i] <<  " :" << hCosthetastar_unfold[i]->Integral(0,10) << endl; 
   cout << "Mass ZDC:" << Sample[i] <<  " :" << hInvmass_ZDC[i]->Integral(0,25) << endl;
   cout << "Single pT:" << Sample[i] <<  " :" << hElePt1[i]->Integral("width") << endl;
   cout << "Single phi:" << Sample[i] <<  " :" << hElePhi1[i]->Integral(0,14) << endl;
   cout << "Single Eta:" << Sample[i] <<  " :" << hEleEta1[i]->Integral(0,23) << endl;
 }
//  cout << "ZDC:" << hZDC[0]-> Integral(0,80) << endl; 
//  cout << "ZDC sum:" << hZDC_sum[0]-> Integral(0,5) << endl;  
//   cout << "ZDC sum:" << hZDC_sum_peak[0]-> Integral(0,5) << endl;
//   cout << "ZDC sum:" << hZDC_sum_peak_max[0]-> Integral(0,5) << endl;
//   cout << "ZDC sum:" << hZDC_sum_max[0]-> Integral(0,5) << endl;

 outf->cd();
 outf->Write();
 outf->Close();

outf1->cd();
hSumPt_unfold[0]->Write();
hInvmass_unfold[0]->Write();
hRap_unfold[0]->Write();
hCosthetastar_unfold[0]->Write();
// hAco_WOZDC[0]->Write();
// hAco_WOZDC[1]->Write();
// hAco_ZDCAND4n[0]->Write();
// hAco_ZDCAND4n[1]->Write();
outf1->Close();

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
 // return dphi
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


