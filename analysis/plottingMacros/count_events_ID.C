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

const double nEventsGeneratedSC_FSR = 67326400; // older number with less files 67262800 superchic 
const double xsecGeneratedSC_FSR    = 8827.220; // μb

double scaleFactorsSC_FSR_Old = 0.85 *  // NEE
                        0.93*      // CHE
                        pow(0.976, 2)* // electron reco+ID
                        0.866 *      // HF veto
                        1.008;       // L1 EG trigger

double scaleFactorsSC_FSR_New = 0.8477 * //NEE
                            0.9322 * //CHE
			    pow(0.952,2)*  //Electron reco+ID
			    0.8643 * //HF veto
			    1.0006;  //L1EG

double lumiNormSC_FSR = xsecGeneratedSC_FSR*luminosity*scaleFactorsSC_FSR_New/nEventsGeneratedSC_FSR;

double scaleFactorPhoton_Old = 0.85 * //NEE
                           0.93* //CHE
			   pow(1.037,2)* //Photon Reco+ID
                           0.866*   //HF Veto
                           1.008;   //L1 EG

double scaleFactorPhoton_New = 0.8477* //NEE
                             0.9322*   //CHE
			     pow(0.9771,2)*  //Photon reco+ID
			     0.8643*   //HF veto
			     1.0006;  //L1 EG

const double xsecGeneratedLbLSC    = 2.59; // μb
const double nEventsGeneratedLbLSC = 466000;
double normSC_LbL = scaleFactorPhoton_New*xsecGeneratedLbLSC*luminosity/nEventsGeneratedLbLSC;

const double xsecGeneratedCEPIncoh    = 0.001431; // μb
const double nEventsGeneratedCEPIncoh = 500000;
double norm_CEPIncoh = scaleFactorPhoton_New*xsecGeneratedCEPIncoh*luminosity/nEventsGeneratedCEPIncoh;


void count_events_ID(){
 
  //gROOT->LoadMacro("CMS_lumi.C");
  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
   gStyle->SetOptStat(0);
//
  TChain *data = new TChain("output_tree");
  //data->Add("After_Diphoton/Data_diphoton_latestNtuple_merged.root");
  data->Add("RootFile_CrossCheckWithJeremi/Data_From_Merged_ntuples_0.root");
  
  TChain *mcsc_fsr = new TChain("output_tree");
  mcsc_fsr->Add("After_Diphoton/QEDSCFSR_diphoton.root");

  TChain *mc_CEPIncoh = new TChain("output_tree");
  mc_CEPIncoh->Add("After_Diphoton/CEPIncoh_diphoton.root");
  
  TChain *mcsc_lbyl = new TChain("output_tree");
  mcsc_lbyl->Add("After_Diphoton/LbyL_diphoton.root");
  //mcsc_lbyl->Add("After_Diphoton/mc_diphoton_RuchiNtuple_merged.root");
//
 

   const TCut no_addittional_tower  = " ok_neuexcl==1";
   //const TCut no_extrk  = " ok_chexcl==1 ";
   const TCut swissCross_photonID = "phoSwissCross_1 < 0.95 && phoSwissCross_2 < 0.95 && phoSeedTime_1 < 3 && phoSeedTime_2 < 3";
   const TCut no_extrk  = " ok_chexcl_goodtracks==1 && ok_chexcl_goodelectrons==1 ";

  TH1D* hmass_data_ID   = new TH1D("hmass_data_ID","",75,0,150);
  TH1D* hmass_data_inv   = new TH1D("hmass_data_inv","",75,0,150);
  TH1D* hmass_data_ch  = new TH1D("hmass_data_ch","",75,0,150);
  TH1D* hmass_data_neu  = new TH1D("hmass_data_neu","",75,0,150);
  TH1D* hmass_data_pt  = new TH1D("hmass_data_pt","",75,0,150);
  TH1D* hmass_data_acop  = new TH1D("hmass_data_acop","",75,0,150);

  data->Project(hmass_data_ID->GetName(), "vSum_M", swissCross_photonID);
  data->Project(hmass_data_inv->GetName(), "vSum_M", swissCross_photonID && "vSum_M>5");
  data->Project(hmass_data_ch->GetName(), "vSum_M", swissCross_photonID && no_extrk && "vSum_M>5");
  data->Project(hmass_data_neu->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5" );
  data->Project(hmass_data_pt->GetName(), "vSum_M", swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1" );
  data->Project(hmass_data_acop->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1 && pho_acop<0.01" );

  TH1D* hmass_mcsc_fsr_ID   = new TH1D("hmass_mcsc_fsr_ID","",75,0,150);
  TH1D* hmass_mcsc_fsr_inv   = new TH1D("hmass_mcsc_fsr_inv","",75,0,150);
  TH1D* hmass_mcsc_fsr_ch   = new TH1D("hmass_mcsc_fsr_ch","",75,0,150);
  TH1D* hmass_mcsc_fsr_neu  = new TH1D("hmass_mcsc_fsr_neu","",75,0,150);
  TH1D* hmass_mcsc_fsr_pt  = new TH1D("hmass_mcsc_fsr_pt","",75,0,150);
  TH1D* hmass_mcsc_fsr_acop  = new TH1D("hmass_mcsc_fsr_acop","",75,0,150);
  
  mcsc_fsr->Project(hmass_mcsc_fsr_ID->GetName(),"vSum_M",swissCross_photonID);
  mcsc_fsr->Project(hmass_mcsc_fsr_inv->GetName(),"vSum_M",swissCross_photonID && "vSum_M>5");
  mcsc_fsr->Project(hmass_mcsc_fsr_ch->GetName(), "vSum_M",swissCross_photonID && no_extrk && "vSum_M>5");
  mcsc_fsr->Project(hmass_mcsc_fsr_neu->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5" );
  mcsc_fsr->Project(hmass_mcsc_fsr_pt->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1" );
  mcsc_fsr->Project(hmass_mcsc_fsr_acop->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1 && pho_acop<0.01" );
  hmass_mcsc_fsr_ID->Scale(lumiNormSC_FSR); hmass_mcsc_fsr_ch->Scale(lumiNormSC_FSR);  hmass_mcsc_fsr_neu->Scale(lumiNormSC_FSR); hmass_mcsc_fsr_inv->Scale(lumiNormSC_FSR);  
  hmass_mcsc_fsr_pt->Scale(lumiNormSC_FSR);   hmass_mcsc_fsr_acop->Scale(lumiNormSC_FSR);

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
  mc_CEPIncoh->Project(hmass_mc_CEPIncoh_acop->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1 && pho_acop<0.01" );

  
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
  mcsc_lbyl->Project(hmass_mcsc_lbyl_acop->GetName(), "vSum_M",swissCross_photonID && no_extrk && no_addittional_tower && "vSum_M>5 && vSum_Pt<1 && pho_acop<0.01" );
  hmass_mcsc_lbyl_ID->Scale(normSC_LbL);hmass_mcsc_lbyl_ch->Scale(normSC_LbL);  hmass_mcsc_lbyl_neu->Scale(normSC_LbL); hmass_mcsc_lbyl_inv->Scale(normSC_LbL);
  hmass_mcsc_lbyl_pt->Scale(normSC_LbL);   hmass_mcsc_lbyl_acop->Scale(normSC_LbL);

  
  cout << "ID cuts   & " <<  hmass_data_ID->Integral() <<"  &  " << hmass_mcsc_lbyl_ID->Integral() <<"  &  " << hmass_mcsc_fsr_ID->Integral() << "   &   " << hmass_mc_CEPIncoh_ID->Integral() << endl;
  cout << "dipho m >5   & " <<  hmass_data_inv->Integral() <<"  &  " << hmass_mcsc_lbyl_inv->Integral() <<"  &  " << hmass_mcsc_fsr_inv->Integral() << "   &   " << hmass_mc_CEPIncoh_inv->Integral() << endl;  
  cout << "chrg excl  & " << hmass_data_ch->Integral() <<"  &  " <<hmass_mcsc_lbyl_ch->Integral() <<"  &  " << hmass_mcsc_fsr_ch->Integral() << "   &   "  << hmass_mc_CEPIncoh_ch->Integral() << endl;
  cout << "neu excl   & " << hmass_data_neu->Integral() <<"  &  " <<  hmass_mcsc_lbyl_neu->Integral() <<"  &  " <<  hmass_mcsc_fsr_neu->Integral() << "   &   " << hmass_mc_CEPIncoh_neu->Integral() << endl;
  cout << "pt < 1: & " << hmass_data_pt->Integral()<<"  &  " << hmass_mcsc_lbyl_pt->Integral() <<"  &  " << hmass_mcsc_fsr_pt->Integral() << "  &  " << hmass_mc_CEPIncoh_pt->Integral() << endl;
  cout << "acop < 0.01:    & "<< hmass_data_acop->Integral() <<"  &  " <<hmass_mcsc_lbyl_acop->Integral() <<"  &  " <<hmass_mcsc_fsr_acop->Integral() << "  &  " << hmass_mc_CEPIncoh_acop->Integral() << endl;




}


