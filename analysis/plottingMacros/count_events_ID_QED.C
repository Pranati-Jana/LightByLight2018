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
//#include "../../CMS_lumi.C"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "THStack.h"
#include "TCut.h"
#include "TAxis.h"
#include "TChain.h"


const double luminosity       = 1647.228136; // μb^-1




double scaleFactorElectron = 0.8497 * 0.9339 * pow(0.943, 2) * 1.0089 * 0.8696 ;
const double xsecGeneratedLbLSC    = 8827.220; // μb
const double nEventsGeneratedLbLSC = 59260000;
//double normSC_LbL = scaleFactorElectron*xsecGeneratedLbLSC*luminosity/nEventsGeneratedLbLSC;
double normSC_LbL =1;

const double xsecGeneratedCEPIncoh    = 7920; // μb
const double nEventsGeneratedCEPIncoh = 66750000;
//double norm_CEPIncoh = scaleFactorElectron*xsecGeneratedCEPIncoh*luminosity/nEventsGeneratedCEPIncoh;
double norm_CEPIncoh = 1;

void count_events_ID_QED(){
 
  //gROOT->LoadMacro("CMS_lumi.C");
  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
   gStyle->SetOptStat(0);
//
  TChain *data = new TChain("output_tree");
  //data->Add("After_Diphoton/Data_diphoton_latestNtuple_merged.root");
  data->Add("After_Dielectron/WithTrigger/Data_qed.root");
  

  TChain *mc_CEPIncoh = new TChain("output_tree");
  mc_CEPIncoh->Add("After_Dielectron/WithTrigger/QEDSL_qed.root");
  
  TChain *mcsc_lbyl = new TChain("output_tree");
  mcsc_lbyl->Add("After_Dielectron/WithTrigger/QEDSCFSR_Private_qed.root");
//
 

   const TCut no_addittional_tower  = " ok_neuexcl==1";
  const TCut swissCross_photonID = "ok_trigger == 1 && elePt_1 > 2 && elePt_2 > 2 && abs(eleEta_1) < 2.2 && abs(eleEta_2) < 2.2 ";
  
   const TCut no_extrk  = " ok_chexcl_extrk ==1 ";

  TH1D* hmass_data_ID   = new TH1D("hmass_data_ID","",75,0,150);
  TH1D* hmass_data_inv   = new TH1D("hmass_data_inv","",75,0,150);
  TH1D* hmass_data_ch  = new TH1D("hmass_data_ch","",75,0,150);
  TH1D* hmass_data_neu  = new TH1D("hmass_data_neu","",75,0,150);
  TH1D* hmass_data_pt  = new TH1D("hmass_data_pt","",75,0,150);
  TH1D* hmass_data_acop  = new TH1D("hmass_data_acop","",75,0,150);

  data->Project(hmass_data_ID->GetName(), "vSum_M");
  data->Project(hmass_data_inv->GetName(), "vSum_M",  "vSum_M>5");
  data->Project(hmass_data_ch->GetName(), "vSum_M",  no_extrk && "vSum_M>5");
  data->Project(hmass_data_neu->GetName(), "vSum_M", no_extrk && no_addittional_tower && "vSum_M>5" );
  data->Project(hmass_data_pt->GetName(), "vSum_M",  no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1" );
  data->Project(hmass_data_acop->GetName(), "vSum_M", no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1 && ele_acop<0.01" );


  TH1D* hmass_mc_CEPIncoh_ID  = new TH1D("hmass_mc_CEPIncoh_ID","",75,0,150);
  TH1D* hmass_mc_CEPIncoh_inv  = new TH1D("hmass_mc_CEPIncoh_inv","",75,0,150);
  TH1D* hmass_mc_CEPIncoh_ch   = new TH1D("hmass_mc_CEPIncoh_ch","",75,0,150);
  TH1D* hmass_mc_CEPIncoh_neu  = new TH1D("hmass_mc_CEPIncoh_neu","",75,0,150);
  TH1D* hmass_mc_CEPIncoh_pt  = new TH1D("hmass_mc_CEPIncoh_pt","",75,0,150);
  TH1D* hmass_mc_CEPIncoh_acop  = new TH1D("hmass_mc_CEPIncoh_acop","",75,0,150);
  
  mc_CEPIncoh->Project(hmass_mc_CEPIncoh_ID->GetName(),"vSum_M",swissCross_photonID);
  mc_CEPIncoh->Project(hmass_mc_CEPIncoh_inv->GetName(),"vSum_M",swissCross_photonID && "vSum_M>5");
  mc_CEPIncoh->Project(hmass_mc_CEPIncoh_ch->GetName(), "vSum_M",swissCross_photonID && no_extrk && "vSum_M>5");
  mc_CEPIncoh->Project(hmass_mc_CEPIncoh_neu->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5" );
  mc_CEPIncoh->Project(hmass_mc_CEPIncoh_pt->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1" );
  mc_CEPIncoh->Project(hmass_mc_CEPIncoh_acop->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1 && ele_acop<0.01" );

  
  hmass_mc_CEPIncoh_ID->Scale(norm_CEPIncoh);hmass_mc_CEPIncoh_ch->Scale(norm_CEPIncoh);  hmass_mc_CEPIncoh_neu->Scale(norm_CEPIncoh); hmass_mc_CEPIncoh_inv->Scale(norm_CEPIncoh);  
  hmass_mc_CEPIncoh_pt->Scale(norm_CEPIncoh);   hmass_mc_CEPIncoh_acop->Scale(norm_CEPIncoh);

  TH1D* hmass_mcsc_lbyl_ID = new TH1D("hmass_mcsc_lbyl_ID","",75,0,150);  
  TH1D* hmass_mcsc_lbyl_inv  = new TH1D("hmass_mcsc_lbyl_inv","",75,0,150);  
  TH1D* hmass_mcsc_lbyl_ch   = new TH1D("hmass_mcsc_lbyl_ch","",75,0,150);
  TH1D* hmass_mcsc_lbyl_neu  = new TH1D("hmass_mcsc_lbyl_neu","",75,0,150);
  TH1D* hmass_mcsc_lbyl_pt  = new TH1D("hmass_mcsc_lbyl_pt","",75,0,150);
  TH1D* hmass_mcsc_lbyl_acop  = new TH1D("hmass_mcsc_lbyl_acop","",75,0,150);

  mcsc_lbyl->Project(hmass_mcsc_lbyl_ID->GetName(), "vSum_M",swissCross_photonID);
  mcsc_lbyl->Project(hmass_mcsc_lbyl_inv->GetName(), "vSum_M",swissCross_photonID && "vSum_M>5");
  mcsc_lbyl->Project(hmass_mcsc_lbyl_ch->GetName(), "vSum_M",swissCross_photonID && no_extrk && "vSum_M>5");
  mcsc_lbyl->Project(hmass_mcsc_lbyl_neu->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5");
  mcsc_lbyl->Project(hmass_mcsc_lbyl_pt->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1" );
  mcsc_lbyl->Project(hmass_mcsc_lbyl_acop->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1 && ele_acop<0.01" );
  hmass_mcsc_lbyl_ID->Scale(normSC_LbL);hmass_mcsc_lbyl_ch->Scale(normSC_LbL);  hmass_mcsc_lbyl_neu->Scale(normSC_LbL); hmass_mcsc_lbyl_inv->Scale(normSC_LbL);
  hmass_mcsc_lbyl_pt->Scale(normSC_LbL);   hmass_mcsc_lbyl_acop->Scale(normSC_LbL);

  
  cout << "ID cuts   & " <<  hmass_data_ID->Integral() <<"  &  " << hmass_mcsc_lbyl_ID->Integral() <<"  &  " << hmass_mc_CEPIncoh_ID->Integral() << endl;
  cout << "dipho m >5   & " <<  hmass_data_inv->Integral() <<"  &  " << hmass_mcsc_lbyl_inv->Integral() <<"  &  " << hmass_mc_CEPIncoh_inv->Integral() << endl;  
  cout << "chrg excl  & " << hmass_data_ch->Integral() <<"  &  " <<hmass_mcsc_lbyl_ch->Integral() <<"  &  " << hmass_mc_CEPIncoh_ch->Integral() << endl;
  cout << "neu excl   & " << hmass_data_neu->Integral() <<"  &  " <<  hmass_mcsc_lbyl_neu->Integral() <<"  &  " << hmass_mc_CEPIncoh_neu->Integral() << endl;
  cout << "pt < 1: & " << hmass_data_pt->Integral()<<"  &  " << hmass_mcsc_lbyl_pt->Integral() <<"  &  " << hmass_mc_CEPIncoh_pt->Integral() << endl;
  cout << "acop < 0.01:    & "<< hmass_data_acop->Integral() <<"  &  " <<hmass_mcsc_lbyl_acop->Integral() <<"  &  " << hmass_mc_CEPIncoh_acop->Integral() << endl;




}


