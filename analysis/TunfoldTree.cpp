//  SelectQEDEvents
//
//  Created by Ruchi Chudasama on 04/03/2021
//
//  Selects diphoton events in MC and data, as defined in the Analysis Note

#include "Helpers.hpp"
#include "EventProcessor.hpp"
#include "PhysObjectProcessor.hpp"
#include "ConfigManager.hpp"
#include "EventDisplay.hpp"
#include "Logger.hpp"
#include "TLorentzVector.h"
#include "Event.hpp"

int ireg(double et, double eta);
double SF_reco(double et, double eta);
double SF_uncert_reco(double et, double eta);
double SF_trig(double et, double eta);
double SF_uncert_trig(double et, double eta);

const double eleMass = 0.5109989461e-3;
Double_t getDR  ( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
Double_t getDETA(Double_t eta1, Double_t eta2);

float Tau_Pt1;
float Tau_Pt;
float Tau_Eta;
float Tau_Phi;
float Tau_E;
float Tau_Eta2;
float Tau_Phi2;

float  Gen_Ele_Pt;
float  Gen_Ele_Eta;
float  Gen_Ele_Phi;
float  Gen_Ele_E;


float  Gen_Mu_Pt;
float  Gen_Mu_Eta;
float  Gen_Mu_Phi;
float  Gen_Mu_E;

float Gen_vSum_M;
float Gen_vSum_Pt;
float Gen_Scalar_Pt;
float Gen_vSum_Eta;
float Gen_vSum_Phi;
float Gen_vSum_Rapidity;
float Gen_EleMu_dphi;
float Gen_EleMu_acop;
int   EleMuEvents;

float Gen_Mu_Pt2;
float Gen_Mu_Eta2;
float Gen_Mu_Phi2;
float Tau_Pt2;
float Gen_vSum_M2;
float Gen_vSum_Pt2;
float Gen_Ele_Pt2;
float Gen_Ele_Eta2;
float Gen_Ele_Phi2;
int hasMuon;
//
int run;
int ls;
int evtnb;
//int    nPho;
//
float  Ele_Pt;
float  Ele_Eta;
float  Ele_Phi;
float  Ele_E;
int    Ele_Charge;

float  Mu_Pt;
float  Mu_Eta;
float  Mu_Phi;
float  Mu_E;
int    Mu_Charge;

float  vSum_M;
//float vSum_M;
float vSum_Energy;
float vSum_Pt;
float Scalar_Pt;
float vSum_Eta;
float vSum_Phi;
float vSum_Rapidity;
float EleMu_dphi;
float EleMu_acop;
int ok_neuexcl;

int ok_chexcl;
int ok_chexcl_tracks;
int ok_chexcl_goodtracks;
int ok_trigger;
int   nTracks;
float dEta_Gen_Reco;
float dPt_Gen_Reco;
float dPhi_Gen_Reco;
float DeltaR_ele;
float DeltaR_mu;
int RecoEleMu;
int Check;
//int hasMuon;
///initialise gen tree
void InitGenTree(TTree *genTree) {
  genTree->Branch("Tau_Pt1",                 &Tau_Pt1,           "Tau_Pt1/F");
  genTree->Branch("Tau_Pt",                 &Tau_Pt,           "Tau_Pt/F");
  genTree->Branch("Tau_Eta",                &Tau_Eta,          "Tau_Eta/F");
  genTree->Branch("Tau_Eta2",                &Tau_Eta2,          "Tau_Eta2/F");
  genTree->Branch("Tau_Phi2",                &Tau_Phi2,          "Tau_Phi2/F");
  genTree->Branch("Tau_Phi",                &Tau_Phi,          "Tau_Phi/F");
  genTree->Branch("Tau_E",                  &Tau_E,            "Tau_E/F");
  genTree->Branch("Gen_Ele_Pt",              &Gen_Ele_Pt,        "Gen_Ele_Pt/F");
  genTree->Branch("Gen_Ele_Eta",             &Gen_Ele_Eta,       "Gen_Ele_Eta/F");
  genTree->Branch("Gen_Ele_Phi",             &Gen_Ele_Phi,       "Gen_Ele_Phi/F");
  genTree->Branch("Gen_Ele_E",               &Gen_Ele_E,         "Gen_Ele_E/F");


  genTree->Branch("Gen_Mu_Pt",               &Gen_Mu_Pt,       "Gen_Mu_Pt/F");
  genTree->Branch("Gen_Mu_Eta",              &Gen_Mu_Eta,      "Gen_Mu_Eta/F");
  genTree->Branch("Gen_Mu_Phi",              &Gen_Mu_Phi,      "Gen_Mu_Phi/F");
  genTree->Branch("Gen_Mu_E",                &Gen_Mu_E,        "Gen_Mu_E/F");

  genTree->Branch("Gen_vSum_M",              &Gen_vSum_M,      "Gen_vSum_M/F");
  genTree->Branch("Gen_vSum_Pt",             &Gen_vSum_Pt,     "Gen_vSum_Pt/F");
  genTree->Branch("Gen_Scalar_Pt",           &Gen_Scalar_Pt,   "Gen_Scalar_Pt/F");
  genTree->Branch("Gen_vSum_Eta",            &Gen_vSum_Eta,    "Gen_vSum_Eta/F");
  genTree->Branch("Gen_vSum_Phi",            &Gen_vSum_Phi,    "Gen_vSum_Phi/F");
  genTree->Branch("Gen_vSum_Rapidity",       &Gen_vSum_Rapidity,"Gen_vSum_Rapidity/F");
  genTree->Branch("Gen_EleMu_dphi",          &Gen_EleMu_dphi,   "Gen_EleMu_dphi/F");
  genTree->Branch("Gen_EleMu_acop",          &Gen_EleMu_acop,   "Gen_EleMu_acop/F");
  genTree->Branch("EleMuEvents",             &EleMuEvents,      "EleMuEvents/I");
  genTree->Branch("Gen_Mu_Pt2",               &Gen_Mu_Pt2,       "Gen_Mu_Pt2/F");
  genTree->Branch("Gen_Mu_Eta2",               &Gen_Mu_Eta2,       "Gen_Mu_Eta2/F");
  genTree->Branch("Gen_Mu_Phi2",               &Gen_Mu_Phi2,       "Gen_Mu_Phi2/F");
  genTree->Branch("Tau_Pt2",                 &Tau_Pt2,           "Tau_Pt2/F");
  genTree->Branch("Gen_vSum_M2",              &Gen_vSum_M2,      "Gen_vSum_M2/F");
  genTree->Branch("Gen_vSum_Pt2",             &Gen_vSum_Pt2,     "Gen_vSum_Pt2/F");
  genTree->Branch("Gen_Ele_Pt2",               &Gen_Ele_Pt2,     "Gen_Ele_Pt2/F");
  genTree->Branch("Gen_Ele_Eta2",               &Gen_Ele_Eta2,     "Gen_Ele_Eta2/F");
  genTree->Branch("Gen_Ele_Phi2",               &Gen_Ele_Phi2,     "Gen_Ele_Phi2/F");
   genTree->Branch("hasMuon",               &hasMuon,     "hasMuon/I");
}

/// initialise tree
void InitTree(TTree *tr) {
  //
  tr->Branch("Tau_Pt1",                 &Tau_Pt1,           "Tau_Pt1/F");
  tr->Branch("Tau_Pt",                 &Tau_Pt,           "Tau_Pt/F");
  tr->Branch("Tau_Eta",                &Tau_Eta,          "Tau_Eta/F");
  tr->Branch("Tau_Phi",                &Tau_Phi,          "Tau_Phi/F");
  tr->Branch("Tau_E",                  &Tau_E,            "Tau_E/F");
  tr->Branch("Gen_Ele_Pt",              &Gen_Ele_Pt,        "Gen_Ele_Pt/F");
  tr->Branch("Gen_Ele_Eta",             &Gen_Ele_Eta,       "Gen_Ele_Eta/F");
  tr->Branch("Gen_Ele_Phi",             &Gen_Ele_Phi,       "Gen_Ele_Phi/F");
  tr->Branch("Gen_Ele_E",               &Gen_Ele_E,         "Gen_Ele_E/F");


  tr->Branch("Gen_Mu_Pt",               &Gen_Mu_Pt,       "Gen_Mu_Pt/F");
  tr->Branch("Gen_Mu_Eta",              &Gen_Mu_Eta,      "Gen_Mu_Eta/F");
  tr->Branch("Gen_Mu_Phi",              &Gen_Mu_Phi,      "Gen_Mu_Phi/F");
  tr->Branch("Gen_Mu_E",                &Gen_Mu_E,        "Gen_Mu_E/F");

  tr->Branch("Gen_vSum_M",              &Gen_vSum_M,      "Gen_vSum_M/F");
  tr->Branch("Gen_vSum_Pt",             &Gen_vSum_Pt,     "Gen_vSum_Pt/F");
  tr->Branch("Gen_Scalar_Pt",           &Gen_Scalar_Pt,   "Gen_Scalar_Pt/F");
  tr->Branch("Gen_vSum_Eta",            &Gen_vSum_Eta,    "Gen_vSum_Eta/F");
  tr->Branch("Gen_vSum_Phi",            &Gen_vSum_Phi,    "Gen_vSum_Phi/F");
  tr->Branch("Gen_vSum_Rapidity",       &Gen_vSum_Rapidity,"Gen_vSum_Rapidity/F");
  tr->Branch("Gen_EleMu_dphi",          &Gen_EleMu_dphi,   "Gen_EleMu_dphi/F");
  tr->Branch("Gen_EleMu_acop",          &Gen_EleMu_acop,   "Gen_EleMu_acop/F");
  tr->Branch("EleMuEvents",             &EleMuEvents,      "EleMuEvents/I");
  tr->Branch("Gen_Mu_Pt2",               &Gen_Mu_Pt2,       "Gen_Mu_Pt2/F");
  tr->Branch("Gen_Mu_Eta2",               &Gen_Mu_Eta2,       "Gen_Mu_Eta2/F");
  tr->Branch("Gen_Mu_Phi2",               &Gen_Mu_Phi2,       "Gen_Mu_Phi2/F");
  tr->Branch("Tau_Pt2",                 &Tau_Pt2,           "Tau_Pt2/F");
  tr->Branch("Tau_Eta2",                 &Tau_Eta2,           "Tau_Eta2/F");
  tr->Branch("Tau_Phi2",                 &Tau_Phi2,           "Tau_Phi2/F");
 // tr->Branch("Gen_Mu_Pt2",               &Gen_Mu_Pt2,       "Gen_Mu_Pt2/F");
  tr->Branch("Gen_vSum_M2",              &Gen_vSum_M2,      "Gen_vSum_M2/F");
  tr->Branch("Gen_vSum_Pt2",             &Gen_vSum_Pt2,     "Gen_vSum_Pt2/F");
  tr->Branch("Gen_Ele_Pt2",               &Gen_Ele_Pt2,     "Gen_Ele_Pt2/F");
  tr->Branch("Gen_Ele_Eta2",               &Gen_Ele_Eta2,     "Gen_Ele_Eta2/F");
  tr->Branch("Gen_Ele_Phi2",               &Gen_Ele_Phi2,     "Gen_Ele_Phi2/F");
////
  tr->Branch("run",                 &run,           "run/I");
  tr->Branch("ls",                  &ls,            "ls/I");
  tr->Branch("evtnb",               &evtnb,         "evtnb/I");
  //tr->Branch("nPho",                &nPho,          "nPho/I");
  tr->Branch("Ele_Pt",              &Ele_Pt,        "Ele_Pt/F");
  tr->Branch("Ele_Eta",             &Ele_Eta,       "Ele_Eta/F");
  tr->Branch("Ele_Phi",             &Ele_Phi,       "Ele_Phi/F");
  tr->Branch("Ele_E",               &Ele_E,         "Ele_E/F");
  tr->Branch("Ele_Charge",          &Ele_Charge,    "Ele_Charge/I");
 
  tr->Branch("Mu_Pt",               &Mu_Pt,       "Mu_Pt/F");
  tr->Branch("Mu_Eta",              &Mu_Eta,      "Mu_Eta/F");
  tr->Branch("Mu_Phi",              &Mu_Phi,      "Mu_Phi/F");
  tr->Branch("Mu_E",             &Mu_E,     "Mu_E/F");
  tr->Branch("Mu_Charge",        &Mu_Charge,  "Mu_Charge/I");
 
  tr->Branch("vSum_M",              &vSum_M,       "vSum_M/F");
 // tr->Branch("vSum_Energy",         &vSum_diPho_Energy,  "vSum_diPho_Energy/F");
  tr->Branch("vSum_Pt",             &vSum_Pt,      "vSum_Pt/F");
  tr->Branch("Scalar_Pt",           &Scalar_Pt,      "Scalar_Pt/F");
  tr->Branch("vSum_Eta",            &vSum_Eta,     "vSum_Eta/F");
  tr->Branch("vSum_Phi",            &vSum_Phi,     "vSum_Phi/F");
  tr->Branch("vSum_Rapidity",       &vSum_Rapidity,"vSum_Rapidity/F");
  tr->Branch("EleMu_dphi",          &EleMu_dphi,   "EleMu_dphi/F");
  tr->Branch("EleMu_acop",          &EleMu_acop,   "EleMu_acop/F");
  tr->Branch("ok_neuexcl",          &ok_neuexcl,   "ok_neuexcl/I");
  

  tr->Branch("ok_chexcl",           &ok_chexcl,       "ok_chexcl/I");
  tr->Branch("ok_chexcl_tracks",           &ok_chexcl_tracks,       "ok_chexcl_tracks/I");
  tr->Branch("ok_chexcl_goodtracks",           &ok_chexcl_goodtracks,       "ok_chexcl_goodtracks/I");
  tr->Branch("ok_trigger",                 &ok_trigger, "ok_trigger/I");
  tr->Branch("dEta_Gen_Reco",              &dEta_Gen_Reco,  "dEta_Gen_Reco/F");
  tr->Branch("dPt_Gen_Reco",              &dPt_Gen_Reco,  "dPt_Gen_Reco/F");
  tr->Branch("dPhi_Gen_Reco",              &dPhi_Gen_Reco,  "dPhi_Gen_Reco/F");
  tr->Branch("DeltaR_ele",                  &DeltaR_ele, "DeltaR_ele/F");
  tr->Branch("DeltaR_mu",                  &DeltaR_mu, "DeltaR_mu/F");
  tr->Branch("RecoEleMu",                  &RecoEleMu, "RecoEleMu/I");
  tr->Branch("Check",                      &Check,     "Check/I");
  tr->Branch("hasMuon",                      &hasMuon,     "hasMuon/I");
}
void ResetGenVars() {
  Tau_Pt1 = -999;
  Tau_Pt = -999;
  Tau_Eta = -999;
  Tau_Phi = -999;
  Tau_E = -999;
  Gen_Ele_Pt= -999 ;
  Gen_Ele_Eta= -999 ;
  Gen_Ele_Phi= -999 ;
  Gen_Ele_E= -999 ;

 Gen_Mu_Pt= -999 ;
 Gen_Mu_Eta= -999 ;
 Gen_Mu_Phi= -999 ;
 Gen_Mu_E= -999 ;
 Gen_vSum_M = -999;
 Gen_vSum_Pt = -999 ;
 Gen_Scalar_Pt = -999;
 Gen_vSum_Eta = -999;
 Gen_vSum_Phi = -999 ;
 Gen_vSum_Rapidity = -999;
 Gen_EleMu_dphi = -999;
 Gen_EleMu_acop = -999;
 EleMuEvents = 0;
 //dEta_Gen_Reco = -999;
 Gen_Mu_Pt2= -999 ;
 Gen_Mu_Eta2 = -999;
 Gen_Mu_Phi2 = -999;
 Tau_Pt2 = -999;
 Tau_Eta2 = -999;
 Tau_Phi2 = -999;
 Gen_vSum_M2 = -999;
 Gen_vSum_Pt2 = -999;
 Gen_Ele_Pt2 = -999;
 Gen_Ele_Eta2 = -999;
 Gen_Ele_Phi2 = -999;
 hasMuon = 0;


}

void ResetVars() {
  run =0;
  ls=0;
  evtnb =0;
  //nMu = 0;
  Ele_Pt= -999 ;
  Ele_Eta= -999 ;
  Ele_Phi= -999 ;
  Ele_E= -999 ;
  Ele_Charge= 0;
  
  Mu_Pt= -999 ;
  Mu_Eta= -999 ;
  Mu_Phi= -999 ;
  Mu_E= -999 ;
  Mu_Charge= 0;

  vSum_M = 0;
  vSum_Energy = 0 ;
  vSum_Pt = 0 ;
  Scalar_Pt = 0;
  vSum_Eta = 0 ;
  vSum_Phi = 0 ;
  vSum_Rapidity = 0 ;
  EleMu_dphi = 0;
  EleMu_acop = 0;
  ok_neuexcl = 0;

  ok_chexcl = 0;
  ok_chexcl_tracks = 0;
  ok_chexcl_goodtracks = 0;
  ok_trigger = 0;
  dEta_Gen_Reco = -999;
  dPt_Gen_Reco = -999;
  dPhi_Gen_Reco = -999;
  DeltaR_ele = -999;
  DeltaR_mu = -999;
  RecoEleMu = 0;
  Check = 0;
//  hasMuon = 0;
}




/// Checks that number of arguments provided is correct and sets corresponding variables
void ReadInputArguments(int argc, char* argv[],
                        string &configPath, string &inputPath, string &outputPath, string &sampleName)
{
  if(argc != 5){
    Log(0)<<"This app requires 4 parameters:\n";
    Log(0)<<"./TUnfoldTree configPath inputPath outputPath datasetName[Data|QED_SC|QED_SL|LbL|CEP]\n";
    exit(0);
  }
  
  configPath = argv[1];
  inputPath  = argv[2];
  outputPath = argv[3];
  sampleName = argv[4];
}


/// Application starting point
int main(int argc, char* argv[])
{
  string configPath, inputPath, outputPath, sampleName;
  ReadInputArguments(argc, argv, configPath, inputPath, outputPath, sampleName);

  config = ConfigManager(configPath);
  EDataset dataset = datasetForName.at(sampleName);
 
   TH1D *hist       = new TH1D("hist","",9,1,10); 
  TFile *outFile = TFile::Open(outputPath.c_str(), "recreate");
  TTree *tr = new TTree("output_tree","");
  InitTree(tr);

  
  auto events = make_unique<EventProcessor>(inputPath, dataset);
  
  int trigger_passed=0, twoPho=0, twoGoodPho=0, oppCharge=0, neutral_excl=0, charged_excl=0, diphomass_wozdc=0, diphopt_wozdc=0;
  int acop_cut_wozdc=0, diphomass=0, diphopt=0, acop_cut = 0, zdc_excl=0;
  
  // Loop over events
  for(int iEvent=0; iEvent<events->GetNevents(); iEvent++){
    if(iEvent%1000 == 0) Log(0)<<"Processing event "<<iEvent<<"\n";
    //if(iEvent >= config.params("maxEvents")) break;
    
    auto event = events->GetEvent(iEvent);
    ResetGenVars();  
    ResetVars();
    if(sampleName == "QED_SC"){

     ResetGenVars();  

     auto gen = event->GetPhysObjects(EPhysObjType::kGenParticle);
    cout << "gen size:" << gen.size() << endl;
    auto GenTaus_1 = event->GetPhysObjects(EPhysObjType::kGenParticle)[0];
    auto GenTaus_2 = event->GetPhysObjects(EPhysObjType::kGenParticle)[1];
    //cout << "Tau_PID1:" << GenTaus_1->GetPID() << endl;    
    //cout << "Tau_PID2:" << GenTaus_2->GetPID() << endl;    
    auto genDau = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter);
      //   cout << "genDau size:" << genDau.size() << endl;

    auto GenTausD_1 = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter)[0];
    auto GenTausD_2 = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter)[1];
   // auto GenTausD_3 = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter)[2];
   // auto GenTausD_4 = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter)[3];
     
  for(auto Taus : gen){
      Tau_Pt  = Taus->GetPt();
      Tau_Eta = Taus->GetEta();
      Tau_Phi = Taus->GetPhi();
      Tau_E   = Taus->GetEnergy();
      //cout << "TauPt:" << Tau_Pt << endl;
     // cout << "TauEta:" << Tau_Eta << endl;
      int  TauD_PID = Taus->GetPID();
       // cout << "Tau PID:" << Taus->GetPID() << endl;
    // } 

 // bool hasMuon = false;
     
//   for(auto Taus : gen){
 
   int numMuons = 0;
   for (auto TauD : genDau) {
    int TauD_PID = TauD->GetPID();
   // cout << "TauD PID:" << TauD_PID << endl;
   // if((TauD_PID == -13 && TauD_PID == -11) || (TauD_PID == 13 && TauD_PID == 11)) continue;

    if (abs(TauD_PID) == 13 ) {
       numMuons++; // Increment the counter for muons
     // cout << "GenMu_Pt:" << TauD->GetPt() << endl;
     //   cout << "GenMu_Eta:" << TauD->GetEta() << endl;
        Gen_Mu_Pt = TauD->GetPt();
        Gen_Mu_Eta = TauD->GetEta();
        Gen_Mu_Phi = TauD->GetPhi();
        Gen_Mu_E = TauD->GetEnergy();
   
      if(Taus->GetPID() * TauD->GetPID() == 195){
        //cout << "Tau Pt decay to muon:" << Taus->GetPt() << endl;
        Tau_Pt2 = Taus->GetPt();
        Tau_Eta2 = Taus->GetEta();
        Tau_Phi2 = Taus->GetPhi();
      }
   // else{
   //  cout << "Tau Pt decay to muon:" << GenTaus_2.GetPt() << endl;
   // }
          
   }
}

   if (numMuons == 1) {
    // This event has exactly one muon with PID 13
     hasMuon = true;
     cout<< "hasMuon:" << hasMuon << endl;
     cout << "One muon exacty" << endl;
}
/* 
 else{

     hasMuon = false;
}
*/
 }   
//}
 /*   for(auto Taus : gen){
      Tau_Pt  = Taus->GetPt();
      Tau_Eta = Taus->GetEta();
      Tau_Phi = Taus->GetPhi();
      Tau_E   = Taus->GetEnergy();
      cout << "TauPt:" << Tau_Pt << endl;
      cout << "TauEta:" << Tau_Eta << endl;
      // int  TauD_PID = Taus->GetPID();
       cout << "Tau PID:" << Taus->GetPID() << endl;
      
      for(auto TauD : genDau){
       int  TauD_PID = TauD->GetPID();
       cout << "TauD PID:" << TauD_PID << endl;
       if(abs(TauD_PID) == 13){
        hasMuon = true;
        cout << "GenMu_Pt:" << TauD->GetPt() << endl;
        cout << "GenMu_Eta:" << TauD->GetEta() << endl;
        Gen_Mu_Pt = TauD->GetPt();
        Gen_Mu_Eta = TauD->GetEta();
        Gen_Mu_Phi = TauD->GetPhi();
        Gen_Mu_E = TauD->GetEnergy();
    }
   }

    }
*/

/*
    bool hasMuon = false;
    for(auto TauD : genDau){
       int  TauD_PID = TauD->GetPID();
       cout << "TauD PID:" << TauD_PID << endl;
       if(abs(TauD_PID) == 13){
        hasMuon = true;
        cout << "GenMu_Pt:" << TauD->GetPt() << endl;
        cout << "GenMu_Eta:" << TauD->GetEta() << endl;
        Gen_Mu_Pt = TauD->GetPt();
        Gen_Mu_Eta = TauD->GetEta();
        Gen_Mu_Phi = TauD->GetPhi();
        Gen_Mu_E = TauD->GetEnergy();
    }
   }
*/
             

/////////////End of Gen info
    ok_trigger = (event->HasTrigger(kSingleMuOpenNoHF));     
    //ok_trigger = (event->HasTrigger(kSingleMuOpenNoHF));     
    cout << "Test1:" << endl; 
    run = event->GetRunNumber();
    ls = event->GetLumiSection();
    evtnb = event->GetEventNumber();
    cout << "Test2:" << endl;
    /******************Only events with muon****************************/
   // if((event->GetPhysObjects(EPhysObjType::kGoodMuon).size() == 1) && (event->GetPhysObjects(EPhysObjType::kGoodElectron).size() == 1)){
    cout << "Test3:" << endl;
   // auto goodgenTracks = event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack);
    if(event->GetPhysObjects(EPhysObjType::kGoodMuon).size() != 1) continue;    
     cout << "Test4:" << endl;
   // ok_neuexcl = (!event->HasAdditionalTowers());
    auto goodMuon =  event->GetPhysObjects(EPhysObjType::kGoodMuon)[0];
     cout << "Test5:" << endl; 
   //cout << "Muon size:" << goodMuon.size() << endl; 
    // ok_chexcl_goodtracks = (goodgenTracks.size()==2);
    // start filling  information here ........................................
   // for(auto goodMuon :  event->GetPhysObjects(EPhysObjType::kGoodMuon)){ 
   
    Mu_Pt        = goodMuon->GetPt();
    Mu_Eta       = goodMuon->GetEta();
    Mu_Phi       = goodMuon->GetPhi();
    Mu_E         = goodMuon->GetEnergy();
    cout << "Reco Muon Eta" << goodMuon->GetEta() << endl; 
 //  }
  // }//Mu size  
 
  //RecoEleMu = (ok_trigger == 1) && (event->GetPhysObjects(EPhysObjType::kGoodMuon).size() == 1) && (event->GetPhysObjects(EPhysObjType::kGoodElectron).size() == 1);
  RecoEleMu = (ok_trigger == 1) && (event->GetPhysObjects(EPhysObjType::kGoodMuon).size() == 1);
  Check = hasMuon && RecoEleMu;
  if(Check){
   DeltaR_mu = getDR(Tau_Eta2, Tau_Phi2, Mu_Eta, Mu_Phi);
   cout << "DeltaR_Mu:" << DeltaR_mu << endl;
  }

 
   }//sample name
    tr->Fill(); 
  // }
  } //nevents
  outFile->cd();
 // hist->Write();
  outFile->Write();
  outFile->Close();
  
return 0;
} // main loop 


Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
  Double_t theDphi = getDPHI( phi1, phi2);
  Double_t theDeta = eta1 - eta2;
  return TMath::Sqrt ( theDphi*theDphi + theDeta*theDeta);
}

Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;
  
  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 )
    dphi = dphi + 2. * 3.141592653589;
  
  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
  }
  
  return TMath::Abs(dphi);
  //return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2){
  return TMath::Abs(eta1 - eta2);
}

bool outside_HEM(Double_t SCEta, Double_t SCPhi){
  if(SCEta < -1.39 && SCEta > -3 && SCPhi > -1.6 && SCPhi < -0.9) return false;
  
  return true;
}



//_____________________________________________________________________________
//     Angle between direction of P1 in the restframe of the pair (P1+P2)
//    and the direction of the pair (P1+P2) in the labframe.
/*

• Collins‐Soper frame: Simplest anisotropy (purely polar)
• Helicity frame: smallest polar and largest azimuthal anisotropies, indefinite tilt
• Gottfried‐Jackson frame: midway between Collins‐Soper and helicity frames

*/

double cosphotonpair(TLorentzVector p1, TLorentzVector pair, bool helicityFrame )
{

  double cos = 0.;

  // Boost of one photon in the pair direction (i.e. to the rest frame of the pair). The other will be at pi rads from the 1st.
  //TVector3 vpair = pair.Vect();
  TVector3 vpair = pair.BoostVector();
  p1.Boost(-vpair);

  // Angle of the boosted-vector with respect to the direction of the pair (Helicity frame)
  cos = TMath::Cos( p1.Angle( pair.Vect() ));

  if (helicityFrame) return cos;

  // Angle of the boosted-vector with respect to the direction of the z-axis (Gottfrid-Jackson frame)
  TVector3 zaxis = TVector3(0.,0.,1.);
  cos = TMath::Cos( p1.Angle( zaxis ));

  return cos;

}

//_____________________________________________________________________________
//  costheta(angle) of the photon in the Collins-Soper frame
/* 

Eq. (2.54) of https://www-d0.fnal.gov/results/publications_talks/thesis/guo_feng/thesis.pdf (SUNY, 2010)
Eq. page 8 of ATLAS paper: arXiv:2107.09330v1 [hep-ex]
*/

//double LbL_QED_lhe_analysis::costhetastar_CS(TLorentzVector p1, TLorentzVector p2, TLorentzVector pair)
double costhetastar_CS(TLorentzVector p1, TLorentzVector pair)
{

  //double Q2 = pair*pair;
  //double QT = pair.Pt();
  //double pp1 = (p1.E()+p1.Pz())/TMath::Sqrt(2); //P^0 and P^3 represent the energy and z component
  //double pm1 = (p1.E()-p1.Pz())/TMath::Sqrt(2);
  //double pp2 = (p2.E()+p2.Pz())/TMath::Sqrt(2);
  //double pm2 = (p2.E()-p2.Pz())/TMath::Sqrt(2);

  //double costhetastarCS = 2.*(pp1*pm2-pm1*pp2)/TMath::Sqrt(Q2*(Q2+QT*QT));

  //double costhetastarCS = 2.*(pair.E()*p1.Pz()-pair.Pz()*p1.E())/(pair.M()*pair.Mt());

  double costhetastarCS = 2.*(pair.E()*p1.Pz()-pair.Pz()*p1.E())/(pair.M()*pair.Mt());

  return costhetastarCS;

}