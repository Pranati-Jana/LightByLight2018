#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include "ReadTree.C"
#include "CMS_lumi.C"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TCut.h"
#include "TChain.h"
#include "THStack.h"

int xlo = 1;
int xhi = 2;
int nbin = 8;

const int nMassbins=4;
//double Massbin[nMassbins]={5.0,7.0,10.0,16.0};
double Massbin[nMassbins]={5.0,9.0,13.0,17.0};
//double Massbin[nMassbins]={5.0,7.3,13.0,16.0};//
//double Massbin[nMassbins]={5.0,7.0,10.0,16.0};//2022
const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;

const int nRapbins=3;
//double Rapbin[nRapbins]={-2.2,-0.71,0.71,2.2};
//double Rapbin[nRapbins]={-2.2,-0.30,0.20,2.2};
//double Rapbin[nRapbins]={-2.2,-0.20,0.72,2.2};
//double Rapbin[nRapbins]={0.0,0.5,1.0,1.5,2.2};
double Rapbin[nRapbins]={0.0,0.8,2.2};
const int nRapbin= sizeof(Rapbin)/sizeof(double) - 1;

void make_canvas(TCanvas *&);
void make_hist(TH1D *&, Color_t , int );
TCanvas* PlotStackHists(TCanvas* , TH1D* , TH1D* , TH1D*, TH1D*, double , double , double , double , double , double, const char *, bool); 

//TCanvas* PlotStackHists(TCanvas* , TH1D* , TH1D* , TH1D*, TH1D*, TH1D*, TH1D*, TH1D*, double , double , double , double , double , double , const char *, bool); 

const int nSample = 4;
const char *Sample[nSample]={"Data","LbyL","CEPIncoh","QEDSCFSR"};
//const char *sample[nSample]={"Data"};

const char *dir = "figures_eta2p2";

#define PI 3.141592653589

const double LumiLossHotZDCneg = 0.047433369; // 1034./21799 from QED number LbyL 5Aug 2022 slides
//const double LumiLossHotZDCneg = 0;
const double luminosity       = 1647.228136; // μb^-1 from Gabi Feb 13, 2022
const double nEventsGeneratedSC = 67810000; // older number with less files 67262800 superchic 
const double xsecGeneratedSC    = 8827.220; // μb
const double purity = 0.96; 
//estimate superchic QED normalisation factors
double scaleFactorsSC = 0.8477 *  // 0.8477NEE    31.12.2021
                        0.9322 *  //0.9322    // CHE  31.12.2021
                        pow(0.952, 2)* //0.952 electron reco+ID 05 Oct 2022 from singleEG3 only
                        1.0006 *       // L1 EG trigger 05 Oct 2022 from singleEG5 only
                        0.8643;  // HF veto
double scaleFactorsSC_Old = 0.85 *  // NEE
                        0.93*      // CHE
                        pow(0.976, 2)* // electron reco+ID
                        0.866 *      // HF veto
                        1.008;       // L1 EG trigger
//Official QED+PHOTOS
const double nEventsGeneratedSCFSR_Private = 59260000;//Private
const double nEventsGeneratedSCFSR_Official = 67326400;//Official 
const double xsecGeneratedSC_Official = 8827.220;
//double lumiNormSCFSR = xsecGeneratedSC_Official*luminosity*scaleFactorsSC*(1-LumiLossHotZDCneg)/nEventsGeneratedSCFSR_Official;
double lumiNormSCFSR = xsecGeneratedSC_Official*luminosity*scaleFactorsSC*(1-LumiLossHotZDCneg)/nEventsGeneratedSCFSR_Private*0.96;
double scaleFactorPhoton = 0.8477 *  //0.8477 NEE    21.12.2021
  0.9322 *      // CHE  21.12.2021
  pow(0.9771, 2)* //0.9771 photon reco+ID 21.12.2021
  0.8643 *      //0.8643 HF veto
  1.0006;       //1.0006 L1 EG trigger
double scaleFactorPhoton_Old = 0.85 * //NEE
                           0.93* //CHE
                           pow(1.037,2)* //Photon Reco+ID
                           0.866*   //HF Veto
                           1.008;   //L1 EG
const double xsecGeneratedLbLSC    = 2.59; // μb Superchic
const double nEventsGeneratedLbLSC = 466000;  //Superchic
double norm_LbLSC = scaleFactorPhoton*xsecGeneratedLbLSC*luminosity*(1-LumiLossHotZDCneg)/nEventsGeneratedLbLSC;   
//double norm_LbLSC = 1;   
const double xsecGeneratedCEPIncoh    = 0.001431; // μb
const double nEventsGeneratedCEPIncoh = 500000; 
double norm_cepIncoh = scaleFactorPhoton*xsecGeneratedCEPIncoh*luminosity*(1-LumiLossHotZDCneg)/nEventsGeneratedCEPIncoh; 

double_t GetDeltaPhi(double dphi, double phi);


void BinLogX(TH1* h);


void plot_data_mc_comp_New(bool CEPIncohNorm){
   if(CEPIncohNorm) {norm_cepIncoh = 1;}
  
   const double wt[nSample] = {1, norm_LbLSC, norm_cepIncoh, lumiNormSCFSR};
  // const double wt[nSample] = {1};

  cout << "LumiNorm SC FSR:" << lumiNormSCFSR << endl;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  string mc_name = "SC";


  TFile *outf= new TFile("diphoton_histos.root","recreate");  

  TChain *tree[nSample];

  //// histograms
  TH1D* hpho_pt[nSample], *hpho_eta[nSample], *hpho_phi[nSample], *hpho_pt2[nSample], *hpho_eta2[nSample], *hpho_phi2[nSample]; 
  TH1D* hdipho_pt[nSample], *hdipho_Rapidity[nSample], *hInvmass[nSample], *hAcoplanarity[nSample], *hdipho_cosThetaStar[nSample], *hdipho_absRapidity[nSample], *hInvmass_unfo[nSample], *hdipho_Rapidity_unfo[nSample];
  TH2D* hpho_EtaPhi[nSample], *hpho_EtaPhi2[nSample];
  TH1D* hnMu[nSample];
  TH1D* hnVtx[nSample], *hxVtx[nSample], *hyVtx[nSample], *hzVtx[nSample];  
  TH1D* hDeltaSeedTime[nSample], *hSeedTime_1[nSample], *hSeedTime_2[nSample];
  TH2D* hSeedTime[nSample];
  TH2D* hSigmaIEta_Eta[nSample], *hSigmaIEta_Eta2[nSample];  
  TH1D* hCosPhoton_helicity0[nSample], *hCosPhoton_helicity1[nSample], *hdeltaPhi1[nSample], *hdeltaPhi2[nSample];
  //
  for (int i = 0; i < nSample; i++){
    tree[i] = new TChain("output_tree");
  //  tree[i]->Add(Form("/home/pranati/Documents/LightByLight/Workspace/CutFlowStudy_2018_Again/forPranati/%s_diphoton.root",Sample[i]));
//    tree[i]->Add(Form("/home/pranati/Documents/LightByLight/Workspace/plottingMacros/forPranati/%s_diphoton.root",Sample[i]));
    tree[i]->Add(Form("After_Diphoton/%s_diphoton.root", Sample[i]));
    hpho_pt[i]    = new TH1D(Form("hpho_pt%s", Sample[i]),"",16,2,10);
    hpho_eta[i]   = new TH1D(Form("hpho_eta%s",Sample[i]),"",11,-2.2,2.2);
    hpho_phi[i]   = new TH1D(Form("hpho_phi%s",Sample[i]),"",12,-PI,PI);

    hpho_pt2[i]    = new TH1D(Form("hpho_pt2%s", Sample[i]),"",16,2,10);
    hpho_eta2[i]   = new TH1D(Form("hpho_eta2%s",Sample[i]),"",11,-2.2,2.2);
    hpho_phi2[i]   = new TH1D(Form("hpho_phi2%s",Sample[i]),"",12,-PI,PI);
 
    hdipho_pt[i]    = new TH1D(Form("hdipho_pt%s", Sample[i]),"",5,0,1);
    hdipho_Rapidity[i]   = new TH1D(Form("hdipho_Rapidity%s",Sample[i]),"",11,-2.2,2.2);
    hInvmass[i]   = new TH1D(Form("hInvmass%s",Sample[i]),"",10,0,50);
    hAcoplanarity[i]   = new TH1D(Form("hAcoplanarity%s",Sample[i]),"",20,0,0.1);
    
    hpho_EtaPhi[i]    = new TH2D(Form("hpho_EtaPhi%s", Sample[i]),"",44,-2.2,2.2,48,-PI,PI);
    hpho_EtaPhi2[i]   = new TH2D(Form("hpho_EtaPhi2%s", Sample[i]),"",44,-2.2,2.2,48,-PI,PI);
    //Muon
    hnMu[i]  = new TH1D(Form("hnMu%s", Sample[i]),"",5,0,5);
    hnVtx[i]  = new TH1D(Form("hnVtx%s", Sample[i]),"",5,0,5);
    hxVtx[i]  = new TH1D(Form("hxVtx%s", Sample[i]),"",10,-1000,50);
    hyVtx[i]  = new TH1D(Form("hyVtx%s", Sample[i]),"",50,-1500,50);
    hzVtx[i]  = new TH1D(Form("hzVtx%s", Sample[i]),"",50,-1500,50);
    //SeedTime//Pranati
    hSeedTime[i] = new TH2D(Form("hSeedTime%s", Sample[i]),"",48,-4,4,48,-4,4);
    hDeltaSeedTime[i] = new TH1D(Form("hDeltaSeedTime%s", Sample[i]),"",25,0,5);
    hSeedTime_1[i] = new TH1D(Form("hSeedTime_1%s", Sample[i]),"",24,-4,4);
    hSeedTime_2[i] = new TH1D(Form("hSeedTime_2%s", Sample[i]),"",24,-4,4);
    ////
    hSigmaIEta_Eta[i] =  new TH2D(Form("hSigmaIEta_Eta%s", Sample[i]), "", 40,0,0.08,44,-2.2,2.2);
    hSigmaIEta_Eta2[i] =  new TH2D(Form("hSigmaIEta_Eta2%s", Sample[i]), "", 40,0,0.08,44,-2.2,2.2);
    //
    hdipho_absRapidity[i] = new TH1D(Form("hdipho_absRapidity%s", Sample[i]), "",2,0,2.2);
    hdipho_Rapidity_unfo[i] = new TH1D(Form("hdipho_Rapidity_unfo%s", Sample[i]), "",nRapbin, Rapbin);
    hInvmass_unfo[i] = new TH1D(Form("hInvmass_unfo%s", Sample[i]), "",nMassbin, Massbin);
    //
    cout << "nSample:" << i  << endl;
    hdipho_cosThetaStar[i]   = new TH1D(Form("hdipho_cosThetaStar%s",Sample[i]),"",4,0,1);
    hCosPhoton_helicity0[i]  = new TH1D(Form("hCosPhoton_helicity0%s",Sample[i]), "", 40, -1,1);
    hCosPhoton_helicity1[i]  = new TH1D(Form("hCosPhoton_helicity1%s",Sample[i]), "", 40, -1,1);
    hdeltaPhi1[i]  = new TH1D(Form("hdeltaPhi1%s",Sample[i]), "", 18, 0, 3.15);
    hdeltaPhi2[i]  = new TH1D(Form("hdeltaPhi2%s",Sample[i]), "", 18, 0, 3.15);
    cout << "file " << Sample[i] << ":" << tree[i]->GetEntries()  << endl;
    ReadTree  treeR(tree[i]);
    treeR.fChain->SetBranchStatus("*",1);


    if (treeR.fChain == 0) return;
  
    Long64_t nentries = treeR.fChain->GetEntriesFast();
    cout << tree[i]->GetName() << "    " << nentries << endl;
  
    Long64_t nbytes = 0, nb = 0;
    //for (Long64_t jentry=0; jentry<10;jentry++) {
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry_evt = treeR.LoadTree(jentry);
      if (ientry_evt < 0) break;
      nb = treeR.fChain->GetEntry(jentry);   nbytes += nb;

      if(treeR.phoSwissCross_1 > 0.95) continue;
      if(treeR.phoSwissCross_2 > 0.95) continue;
      if(treeR.ok_neuexcl != 1) continue; //neutral exclusivity
     // if(treeR.ok_chexcl_goodtracks != 1) continue; //charged exclusivity
    if(treeR.ok_chexcl_goodtracks != 1 || treeR.ok_chexcl_goodelectrons !=1 || treeR.ok_chexcl_muons !=1) continue; //charged exclusivity
//    if(treeR.ok_chexcl != 1) continue;
 //  if(treeR.ok_chexcl_tracks != 1 || treeR.ok_chexcl_electrons !=1 ) continue;//Only to check loose exc//
 //     cout <<
     // }
      
     // if(treeR.ok_trigger!=1) continue;
      if(treeR.vSum_M < 5) continue; //invmass
      if(treeR.vSum_Pt > 1) continue; //diphoton pt
      if(abs(treeR.phoEta_1) > 2.2) continue;
      if(abs(treeR.phoEta_2) > 2.2) continue;
      //
      if(treeR.phoEt_1 < 2.5) continue;
      if(treeR.phoEt_2 < 2.5) continue;
      if(abs(treeR.phoSeedTime_1) > 3) continue;
      if(abs(treeR.phoSeedTime_2) > 3) continue;
      //if(treeR.phoSigmaIEta_1 > 0.02) continue;
     // if(treeR.phoSigmaIEta_2 > 0.02) continue;
      //SigmaIEta cut
      
      if ((abs(treeR.phoEta_1) >= 0 && abs(treeR.phoEta_1) < 1.442 && treeR.phoSigmaIEta_1 < 0.02) ||
         (abs(treeR.phoEta_1) > 1.442 && treeR.phoSigmaIEta_1 < 0.06)) {
	       } 
         else {
		   continue;
         }
    if ((abs(treeR.phoEta_2) >= 0 && abs(treeR.phoEta_2) < 1.442 && treeR.phoSigmaIEta_2 < 0.02) ||
       (abs(treeR.phoEta_2) > 1.442 && treeR.phoSigmaIEta_2 < 0.06)) {
               } 
     else {
                   continue;
      }




      if(i==0) {
        // if(treeR.ok_zdcexcl_3n_pos != 1 || treeR.ok_zdcexcl_3n_neg != 1) continue;
        // if(treeR.ok_zdcexcl_1n_pos != 1 || treeR.ok_zdcexcl_1n_neg != 1) continue;//0n0n
        // if(treeR.zdc_energy_pos > 1500 || treeR.zdc_energy_neg > 1500)continue;
        if(treeR.zdc_energy_pos >10000  || treeR.zdc_energy_neg > 10000)continue;
        //if(treeR.zdc_energy_pos >12500  || treeR.zdc_energy_neg > 12500)continue;
         //if(treeR.ok_zdcexcl!= 1) continue; //bad ZDC negative runs removed in ok_zdcexcl variable
       }
 
      //if(i==1)
      //{  if(treeR.ok_trigger!=1) continue; }


      hAcoplanarity[i]->Fill(treeR.pho_acop,wt[i]);

      if(treeR.pho_acop > 0.1) {
	hpho_EtaPhi[i]->Fill(treeR.phoEta_1,treeR.phoPhi_1,wt[i]);
        hpho_EtaPhi2[i]->Fill(treeR.phoEta_2,treeR.phoPhi_2,wt[i]);
        }
      
        if(treeR.pho_acop >  0.01) continue; //acop
       
      	
	hpho_pt[i]->Fill(treeR.phoEt_1,wt[i]);
	hpho_eta[i]->Fill(treeR.phoEta_1,wt[i]);
	hpho_phi[i]->Fill(treeR.phoPhi_1,wt[i]);
	
	hpho_pt2[i]->Fill(treeR.phoEt_2,wt[i]);
	hpho_eta2[i]->Fill(treeR.phoEta_2,wt[i]);
	hpho_phi2[i]->Fill(treeR.phoPhi_2,wt[i]);
	
	hdipho_pt[i]->Fill(treeR.vSum_Pt,wt[i]);
	hdipho_Rapidity[i]->Fill(treeR.vSum_Rapidity,wt[i]);
	hInvmass[i]->Fill(treeR.vSum_M,wt[i]);
  //      cout << "Invmass:" << treeR.vSum_M << endl;
        hdipho_cosThetaStar[i]->Fill(fabs(treeR.costhetastar),wt[i]);
	hCosPhoton_helicity0[i]->Fill(treeR.cos_photon_pair_helicity0, wt[i]);
	hCosPhoton_helicity1[i]->Fill(treeR.cos_photon_pair_helicity1, wt[i]);
        //Muon
        hnMu[i]->Fill(treeR.nMu);	     
  //      cout << "Muon size:" << treeR.nMu << endl;
        hSeedTime_1[i]->Fill(treeR.phoSeedTime_1,wt[i]);
        hSeedTime_2[i]->Fill(treeR.phoSeedTime_2,wt[i]);
        hSeedTime[i]->Fill(treeR.phoSeedTime_1,treeR.phoSeedTime_2,wt[i]);
        hDeltaSeedTime[i]->Fill(fabs(treeR.phoSeedTime_1 - treeR.phoSeedTime_2),wt[i]);
        //
        //cout << "phoSigmaIEta_2:" << treeR.phoSigmaIEta_2 << endl;
        //cout << "phoEta_2:" << treeR.phoEta_2 << endl;
        hSigmaIEta_Eta[i]->Fill(treeR.phoSigmaIEta_1,treeR.phoEta_1,wt[i]);
        hSigmaIEta_Eta2[i]->Fill(treeR.phoSigmaIEta_2,treeR.phoEta_2,wt[i]);
        hdipho_absRapidity[i]->Fill(fabs(treeR.vSum_Rapidity),wt[i]);      
        hdipho_Rapidity_unfo[i]->Fill(fabs(treeR.vSum_Rapidity),wt[i]);
	hInvmass_unfo[i]->Fill(treeR.vSum_M,wt[i]); 
        double dphi = treeR.pho_dphi;
        double phi = treeR.phoPhi_1;	
        double deltaPhi = dphi - phi;
//	cout << "dPhi:" << dphi << endl;
/*
	if(treeR.phoEt_1 > treeR.phoEt_2 ){
	hdeltaPhi1[i]->Fill((treeR.pho_dphi-fabs(treeR.phoPhi_1)), wt[i]);
	cout << "Phi1:" << treeR.phoPhi_1 << endl;
	cout << "deltaPhi1:" <<  (treeR.pho_dphi-(treeR.phoPhi_1))<< endl;
        }
	else{
        hdeltaPhi2[i]->Fill((treeR.pho_dphi-fabs(treeR.phoPhi_2)), wt[i]);
	cout << "Phi2:" << treeR.phoPhi_2 << endl;;
	cout << "deltaPhi2:" <<  (treeR.pho_dphi-(treeR.phoPhi_2))<< endl;
	}
*/
        if(treeR.phoEt_1 > treeR.phoEt_2 ){
	double dPhi = treeR.vSum_Phi;
        double phi = (treeR.phoPhi_1);
        double cosphi1 = GetDeltaPhi(dPhi , phi);

        hdeltaPhi1[i]->Fill(cosphi1, wt[i]);
//	hdeltaPhi1[i]->Fill((treeR.pho_dphi-fabs(treeR.phoPhi_1)), wt[i]);
        }
        else{
	double dPhi = treeR.vSum_Phi;
        double phi = (treeR.phoPhi_2);
        double cosphi1 = GetDeltaPhi(dPhi , phi);

        hdeltaPhi1[i]->Fill(cosphi1,wt[i]);
//	hdeltaPhi1[i]->Fill((treeR.pho_dphi-fabs(treeR.phoPhi_2)), wt[i]);
        }


  } //entry

     if(treeR.pho_acop > 0.09){
//   std::cout << "pho1 et:" << treeR.phoEt_1 << std::endl;

 }
  } // for 4 files
  
  for (int i = 0; i < nSample; i++){    
    hpho_pt[i]->Add(hpho_pt2[i]);
    hpho_eta[i]->Add(hpho_eta2[i]);
    hpho_phi[i]->Add(hpho_phi2[i]);
    hpho_EtaPhi[i]->Add(hpho_EtaPhi2[i]);
    hSigmaIEta_Eta[i]->Add(hSigmaIEta_Eta2[i]);
//    hdeltaPhi1[i]->Add(hdeltaPhi2[i]);
  }
  


 TCanvas* ccc = new TCanvas();
// hdeltaPhi[0]->Draw("P");
 hdeltaPhi1[1]->SetMarkerStyle(20);
 //hdeltaPhi1[1]->GetXaxis->SetTitle("Events");
 hdeltaPhi1[1]->Draw("p");
 ccc->Print("Plot/tmp/deltaPhi.png");
 cout << "DeltaPhi:" << hdeltaPhi1[1]->Integral(0,18) << endl;  

 TCanvas* ccc1 = new TCanvas();
 hnMu[0]->Draw("p");
 ccc1->Print("Plot/tmp/nMu.png");



   if(CEPIncohNorm){
  double old_norm = norm_cepIncoh;
  double norm_cep2 = (hAcoplanarity[0]->Integral(5,20) - hAcoplanarity[1]->Integral(5,20) - hAcoplanarity[3]->Integral(5,20) ) / hAcoplanarity[2]->Integral(5,20);
 
  std::cout << "------------ Scaling CEP Incoherent to tail -----------------------------------------" << std::endl; 
    
  for (int i = 1; i <= hAcoplanarity[0]->GetNbinsX(); i++) std::cout << i << "  " << hAcoplanarity[0]->Integral(i,i)  << std::endl;
  std::cout << "data   " << hAcoplanarity[0]->Integral(5,20)  << std::endl;
  //std::cout << "data   " << h
  std::cout << "LbyL MC   " << hAcoplanarity[1]->Integral(5,20)  << std::endl;
  std::cout << "QED SC FSR  " << hAcoplanarity[3]->Integral(5,20)  << std::endl;
  std::cout << "CEP incoh  " << hAcoplanarity[2]->Integral(5,20)  << std::endl;
  std::cout << "new lumi Norm   " << norm_cep2 << std::endl;
  std::cout << "CEP Incoherent Norm:" << norm_cep2 <<std::endl;
  hpho_pt[2]->Scale(norm_cep2);   hpho_eta[2]->Scale(norm_cep2);   hpho_phi[2]->Scale(norm_cep2);     
  hAcoplanarity[2] ->Scale(norm_cep2);   
  hdipho_pt[2]->Scale(norm_cep2);   hdipho_Rapidity[2]->Scale(norm_cep2);     hInvmass[2]->Scale(norm_cep2);   hdipho_absRapidity[2]->Scale(norm_cep2);    hInvmass_unfo[2]->Scale(norm_cep2);  hdipho_Rapidity_unfo[2]->Scale(norm_cep2); hdipho_cosThetaStar[2]->Scale(norm_cep2); 
  }

 
  
  
  /*int W = 700;
  int H = 600;
  
  float T = 0.08;
  float B = 0.14; 
  float L = 0.14;
  float R = 0.04;*/

  //hInvmass[0]->Draw("p");
  
  //Muon passing
   TCanvas*c = new TCanvas();
   hnMu[0]->Draw("HIST");
//   c->Print("NMuon.png");
/*   TCanvas*c20 = new TCanvas();
   TCanvas*c24 = new TCanvas();
   hSeedTime[0]->Draw("COLZ");
   //c24->Print("SeedTime/SeedTime2D_2.png");
   TCanvas*c25 = new TCanvas();
   hDeltaSeedTime[0]->Draw("HIST");
   //c25->Print("SeedTime/DeltaSeedTime_2.png");
   TCanvas*c26 = new TCanvas();
   hSeedTime_1[0]->Draw("HIST");
   //c26->Print("SeedTime/SeedTime_1.png");
   TCanvas*c27 = new TCanvas();
   hSeedTime_2[0]->Draw("HIST");
   //c27->Print("SeedTime/SeedTime_2.png");
   */
     TCanvas*c20 = new TCanvas();
//     c20->SetOptstat(0);
     hSigmaIEta_Eta[0]->GetXaxis()->SetTitle("SigmaIEtaIEta");
     hSigmaIEta_Eta[0]->GetYaxis()->SetTitle("Photon #eta");
     hSigmaIEta_Eta[0]->Draw("COLZ");
     TCanvas*c21 = new TCanvas();  
     hdipho_cosThetaStar[0]->Draw("P");
     c21->Print("Plot/tmp/costheta.png");
     TCanvas*c22 = new TCanvas();
     hCosPhoton_helicity0[0]->Draw("P");
     c22->Print("Plot/tmp/costhelicity0.png");
     TCanvas*c23 = new TCanvas();
     hCosPhoton_helicity0[0]->Draw("P");
     c23->Print("Plot/tmp/costhelicity1.png");
     //(TCanvas* c1, TH1D* hdata, TH1D* hmc, TH1D* hmc2, double hxmin, double hxmax, double hymin, double hymax, double rymin , double rymax, const char *ytitle, bool iflogy)



  TCanvas* cc1 = new TCanvas("Sum_pt","diphoton pT",254,411,639,592);
  make_canvas(cc1);
  PlotStackHists(cc1, hdipho_pt[0], hdipho_pt[1], hdipho_pt[2], hdipho_pt[3], 0.0,2.0,0,30,0.7,1.8,"Diphoton p_{T}", 0);

  TCanvas* c2 = new TCanvas("Inv_mass","Invmass",254,411,639,592);
  make_canvas(c2);
  PlotStackHists(c2, hInvmass[0], hInvmass[1], hInvmass[2], hInvmass[3], 5,30,0,30,0.6,1.4,"Invariant Mass", 0);

  TCanvas* c3 = new TCanvas("Rapidity","Rapidity",254,411,639,592);
  make_canvas(c3);
  PlotStackHists(c3, hdipho_Rapidity[0], hdipho_Rapidity[1], hdipho_Rapidity[2], hdipho_Rapidity[3], -3,3,0,20,0.7,1.3,"Rapidity", 0);

  TCanvas* c4 = new TCanvas("Photon_pt","Photon pT",254,411,639,592);
  make_canvas(c4);
  PlotStackHists(c4, hpho_pt[0], hpho_pt[1],  hpho_pt[2], hpho_pt[3], 0,50,0.1,30,0.1,1.8,"Photon p_{T}", 0);

  TCanvas* c5 = new TCanvas("Photon_eta","Photon eta",254,411,639,592);
  make_canvas(c5);
  PlotStackHists(c5, hpho_eta[0], hpho_eta[1], hpho_eta[2],  hpho_eta[3], -3,3,0,30,0.7,1.3,"Photon #eta", 0);

  TCanvas* c6 = new TCanvas("Photon_phi","Photon phi",254,411,639,592);
  make_canvas(c6);
  PlotStackHists(c6, hpho_phi[0], hpho_phi[1], hpho_phi[2], hpho_phi[3],  -4,4,0,30,0.7,1.3,"Photon #phi", 0);


  
  TCanvas* c7 = new TCanvas("Acoplanarity","Acoplanarity",2563,306,777,575);
  make_canvas(c7);
  PlotStackHists(c7, hAcoplanarity[0], hAcoplanarity[1], hAcoplanarity[2], hAcoplanarity[3],  0,0.16,0,30,0.,2.5,"A_{#phi}", 0);
  
  TCanvas* c8 = new TCanvas("Abs_rapidity","Abs_rapidity",2563,306,777,575);
  make_canvas(c8);
  PlotStackHists(c8, hdipho_absRapidity[0], hdipho_absRapidity[1], hdipho_absRapidity[2], hdipho_absRapidity[3], 0,1.6,0,30,0.7,1.3,"Rapidity", 0);

   TCanvas* c9 = new TCanvas("Invmass_unfo","Invmass_unfo",2563,306,777,575);
  make_canvas(c9);
  PlotStackHists(c9, hInvmass_unfo[0], hInvmass_unfo[1], hInvmass_unfo[2], hInvmass_unfo[3], 0,120,0,30,0.6,1.4,"Invmass", 0);

  TCanvas* c10 = new TCanvas("Rapidity_unfo","Rapidity_unfo",254,411,639,592);
  make_canvas(c10);
  PlotStackHists(c10, hdipho_Rapidity_unfo[0], hdipho_Rapidity_unfo[1], hdipho_Rapidity_unfo[2], hdipho_Rapidity_unfo[3],  -3,3,0,30,0.7,1.3,"Rapidity", 0);


  TCanvas* c11 = new TCanvas("abs_Costhetastar","abs_Costhetastar",254,411,639,592);
  make_canvas(c11);
  PlotStackHists(c11, hdipho_cosThetaStar[0], hdipho_cosThetaStar[1], hdipho_cosThetaStar[2], hdipho_cosThetaStar[3],  0,1,0,20,0.7,1.3,"|cos#theta*|", 0);

  for (int i = 0; i < nSample; i++){
  cout << " Acop < 0.01 " << Sample[i] <<  " :" << hAcoplanarity[i]->Integral(0,2) << endl;  
  cout << "Rapidity at  Acop < 0.01 " << Sample[i] <<  " :" << hdipho_Rapidity[i]->Integral(0,11) << endl;  
  cout << "pT at  Acop < 0.01 " << Sample[i] <<  " :" << hdipho_pt[i]->Integral(0,10) << endl;  
  cout << " Invmass at Acop < 0.01 " << Sample[i] <<  " :" << hInvmass[i]->Integral(0,30) << endl;  
  cout << " Unfo Invmass at Acop < 0.01 " << Sample[i] <<  " :" << hInvmass_unfo[i]->Integral(0,4) << endl;  
  cout << " Unfo Rapidity at Acop < 0.01 " << Sample[i] <<  " :" << hdipho_Rapidity_unfo[i]->Integral(0,2) << endl;  
 }
 
/*  TLegend *leg2=new TLegend(0.55,0.60,0.90,0.91);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(20); 
  leg2->AddEntry(hAcoplanarity[4],"CEP Incoh (gg #rightarrow #gamma #gamma)","pl");
  leg2->AddEntry(hAcoplanarity[3],"CEP coh (gg #rightarrow #gamma #gamma)","pl");
   
 new TCanvas();
 make_hist(hAcoplanarity[4], kAzure+1, 21);
 hAcoplanarity[4]->GetXaxis()->SetTitle("Acop");//,"dN/dp_{T}^{2} (1/GeV)^{2} 
 hAcoplanarity[4]->GetYaxis()->SetTitle(" Events");
 hAcoplanarity[4]->Draw("p");  
 make_hist(hAcoplanarity[3], kOrange, 21);
 hAcoplanarity[3]->Draw("psame");  
 leg2->Draw();
 
 new TCanvas();
 hpho_EtaPhi[0]->Draw("colz");
 

  
  */

 outf->cd();
 outf->Write();
 outf->Close();
 
}

TCanvas* PlotStackHists(TCanvas* c1, TH1D* hdata, TH1D* hlbyl, TH1D* hcep,TH1D *hfsr, double hxmin, double hxmax, double hymin, double hymax, double rymin , double rymax, const char *ytitle, bool iflogy){
//TCanvas* PlotStackHists(TCanvas* c1, TH1D* hdata, TH1D* hlbyl, TH1D* hqed, TH1D* hcep,TH1D *hfsr, TH1D *hMG5, TH1D *hMG5_2FSR, double hxmin, double hxmax, double hymin, double hymax, double rymin , double rymax, const char *ytitle, bool iflogy){
  
  float T = 0.08;
  float B = 0.14; 
  float L = 0.14;
  float R = 0.04;

  /*c1->Divide(1,2);
  TPad* pad1 = (TPad*)c1->GetPad(1);
  TPad* pad2 = (TPad*)c1->GetPad(2);
  pad1->SetPad(0,0.3,1.0,1.0);
  pad2->SetPad(0,0,  1.0,0.28);
  pad1->SetLeftMargin(0.18);
  pad1->SetTopMargin(0.08);
  pad1->SetRightMargin(0.07);
  pad1->SetBottomMargin(0.01); // All X axis labels and titles are thus cut off
  pad2->SetLeftMargin(0.18);
  pad2->SetTopMargin(0.01);
  pad2->SetRightMargin(0.07);
  pad2->SetBottomMargin(0.45);
  
  pad1->cd();*/

  //hdata->Add(hlbyl,-1); 
  //hdata->Add(hfsr,-1); 
 
  THStack *hs = new THStack("hs"," ");
  if(iflogy)gPad->SetLogy();

  hs->SetMaximum(hymax);
  hs->SetMinimum(hymin);
  //hs->GetXaxis()->SetMinimum(hxmin);
 
  make_hist(hcep, kAzure+1, 21);
  hs->Add(hcep);

  make_hist(hfsr, kYellow, 21);
  hs->Add(hfsr);

 // make_hist(hcep, kAzure+1, 21);
 // hs->Add(hcep);

  make_hist(hlbyl,kOrange+7, 21);
  hs->Add(hlbyl);
 // make_hist(hqed,kRed, 21);
 // hs->Add(hqed);
  
  hs->Draw("hist");   
  hs->GetXaxis()->SetTitle(ytitle);//,"dN/dp_{T}^{2} (1/GeV)^{2} 
  hs->GetYaxis()->SetTitle("# events");

  make_hist(hdata, kBlack, 21);
  hdata->SetLineColor(kBlack);
  hdata->Draw("e1same");
  
   TLegend *leg2=new TLegend(0.55,0.60,0.90,0.91);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(20); 
  leg2->AddEntry(hdata,"Data","pl");
  leg2->AddEntry(hlbyl,"LbyL #gamma #gamma #rightarrow #gamma #gamma (MC)","f");
  leg2->AddEntry(hfsr,"QED SC +Photos","f");
  leg2->AddEntry(hcep,"CEP Incoh (gg #rightarrow #gamma #gamma)","f");
//  leg2->AddEntry(hfsr,"QED SC +Photos","f");
  leg2->Draw();
  
       TLatex *mark = new TLatex();
       mark->SetNDC(true);
       double startY = 0.93;
       mark->SetTextFont(42);
       mark->SetTextSize(0.035);
       mark->DrawLatex(0.18,startY,"#bf{CMS} #it{Preliminary}");
       mark->DrawLatex(0.65,startY,"#scale[0.8]{PbPb, 1647.18 #mub^{-1} (5.02 TeV)}");
       mark->Draw();

  
  

  /*pad2->cd();
  pad2->SetGridy(1);
  hratio->SetLineWidth(2);
  //hratio->SetStats(0);
  hratio->GetXaxis()->SetRangeUser(hxmin,hxmax);
  hratio->GetYaxis()->SetRangeUser(rymin,rymax);
  //hratio->GetXaxis()->SetTitle("Dielectron p_{T}");
  hratio->GetXaxis()->SetTitle(ytitle);
  hratio->GetYaxis()->SetTitle("QED+FSR/QED");
  hratio->GetYaxis()->SetNdivisions(504);
  hratio->GetYaxis()->SetLabelSize(0.11);
  hratio->GetXaxis()->SetLabelSize(0.15);
  hratio->GetXaxis()->SetTitleSize(0.18);
  hratio->GetYaxis()->SetTitleSize(0.15);
  hratio->GetXaxis()->SetTitleOffset(1.1);
  hratio->GetYaxis()->SetTitleOffset(0.31);
  hratio->GetXaxis()->SetTickLength(0.09);
  //hratio->SetMarkerColor(33);
  hratio->SetMarkerColor(kBlack);
  //hratio->SetMarkerStyle(23);
  hratio->SetLineColor(kBlack);
  //hratio->Fit("pol0","W");
 // hratio->Draw("p");*/

  //CMS_lumi( c1, 104, 33,lumi_PbPb2018 );

  c1->Update();
  TString cName=c1->GetName();
  cName+=".png";
  c1->SaveAs("Plot/tmp/"+cName);
  TString c2Name=c1->GetName();
  c2Name+=".pdf";
  c1->SaveAs("Plot/tmp/"+c2Name);
  return c1;
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
  hist->SetLineColor(kBlack);
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


