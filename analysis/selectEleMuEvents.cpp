//  SelectEleMuEvents_Tree
//
//  Modified Ruchi Chudasama's selectQEDEvents
//
//  Selects eleMu events in MC and data, as defined in the Analysis Note

#include "Helpers.hpp"
#include "EventProcessor.hpp"
#include "PhysObjectProcessor.hpp"
#include "ConfigManager.hpp"
#include "EventDisplay.hpp"
#include "Logger.hpp"
#include "TLorentzVector.h"
#include "Event.hpp"


const double eleMass = 0.5109989461e-3;
Double_t getDR  ( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
Double_t getDETA(Double_t eta1, Double_t eta2);


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
int Ele_charge;

float  Mu_Pt;
float  Mu_Eta;
float  Mu_Phi;
float  Mu_E;
int Mu_charge;

float vSum_M;
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
int ok_chexcl_extratracks;
int ok_chexcl_goodtracks;
int ok_trigger;
int   nTracks;
int ok_singleMuon;
int ok_singleEG3;
int ok_singleEG5;
float deltaRtrack_electron;
float deltaRtrack_muon;
float zdc_energy_pos;
float zdc_energy_neg;
float track_Pt;
float track_Eta;
float track_Phi;
int track_Charge;
int nExtrk;
float vSum_Pz;
float tower_eta;
float tower_phi;
float tower_em;
float tower_had;
float tower_deltaRmuon_HE;
float tower_deltaRmuon_HB;
float tracks_deltaRElectron;
float deltaRtrack_electron_test;
float deltaRtrack_muon_test;
float deltaRtrack_electron_dEta;
float deltaRtrack_electron_dPhi;
float tower_HF_Plus;
float tower_HF_Minus;
/// initialise tree
void InitGenTree(TTree *genTree) {
   genTree->Branch("deltaRtrack_electron_test",                 &deltaRtrack_electron_test, "deltaRtrack_electron_test/F");
   genTree->Branch("deltaRtrack_muon_test",                     &deltaRtrack_muon_test, "deltaRtrack_muon_test/F");
   genTree->Branch("deltaRtrack_electron_dEta",                     &deltaRtrack_electron_dEta, "deltaRtrack_electron_dEta/F");
   genTree->Branch("deltaRtrack_electron_dPhi",                     &deltaRtrack_electron_dPhi, "deltaRtrack_electron_dPhi/F");
}

void InitTree(TTree *tr) {
  //
  tr->Branch("run",                 &run,           "run/I");
  tr->Branch("ls",                  &ls,            "ls/I");
  tr->Branch("evtnb",               &evtnb,         "evtnb/I");
  tr->Branch("Ele_Pt",              &Ele_Pt,        "Ele_Pt/F");
  tr->Branch("Ele_Eta",             &Ele_Eta,       "Ele_Eta/F");
  tr->Branch("Ele_Phi",             &Ele_Phi,       "Ele_Phi/F");
  tr->Branch("Ele_E",               &Ele_E,         "Ele_E/F");
  tr->Branch("Ele_charge",          &Ele_charge,    "Ele_charge/I");
 
  tr->Branch("Mu_Pt",               &Mu_Pt,        "Mu_Pt/F");
  tr->Branch("Mu_Eta",              &Mu_Eta,       "Mu_Eta/F");
  tr->Branch("Mu_Phi",              &Mu_Phi,       "Mu_Phi/F");
  tr->Branch("Mu_E",                &Mu_E,         "Mu_E/F");
  tr->Branch("Mu_charge",           &Mu_charge,    "Mu_charge/I");
  
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
  tr->Branch("ok_chexcl_extratracks",           &ok_chexcl_extratracks,       "ok_chexcl_extratracks/I");
  tr->Branch("ok_chexcl_goodtracks",           &ok_chexcl_goodtracks,       "ok_chexcl_goodtracks/I");
  tr->Branch("ok_trigger",                 &ok_trigger, "ok_trigger/I");
  tr->Branch("ok_singleMuon",                 &ok_trigger, "ok_singleMuon/I");
  tr->Branch("ok_singleEG3",                 &ok_trigger, "ok_singleEG3/I");
  tr->Branch("ok_singleEG5",                 &ok_trigger, "ok_singleEG5/I");
  tr->Branch("deltaRtrack_electron",                 &deltaRtrack_electron, "deltaRtrack_electron/F");
  tr->Branch("deltaRtrack_muon",                 &deltaRtrack_muon, "deltaRtrack_muon/F");
  tr->Branch("zdc_energy_pos",             &zdc_energy_pos,         "zdc_energy_pos/F");
  tr->Branch("zdc_energy_neg",             &zdc_energy_neg,         "zdc_energy_neg/F");
  tr->Branch("nTracks",             &nTracks,         "nTracks/I");
  tr->Branch("track_Pt",             &track_Pt,         "track_Pt/F");
  tr->Branch("track_Eta",             &track_Eta,         "track_Eta/F");
  tr->Branch("track_Phi",             &track_Phi,         "track_Phi/F");
  tr->Branch("track_Charge",             &track_Charge,         "track_Charge/I");
  tr->Branch("nExtrk",              &nExtrk,          "nExtrk/I");
  tr->Branch("vSum_Pz",              &vSum_Pz,       "vSum_Pz/F");
  
  tr->Branch("tower_eta",              &tower_eta,       "tower_eta/F");
  tr->Branch("tower_phi",              &tower_phi,       "tower_phi/F");
  tr->Branch("tower_em",              &tower_em,       "tower_em/F");
  tr->Branch("tower_had",              &tower_had,       "tower_had/F");
  tr->Branch("tower_deltaRmuon_HE",              &tower_deltaRmuon_HE,       "tower_deltaRmuon_HE/F");
  tr->Branch("tower_deltaRmuon_HB",              &tower_deltaRmuon_HB,       "tower_deltaRmuon_HB/F");
  tr->Branch("tracks_deltaRElectron",              &tracks_deltaRElectron,      "tracks_deltaRElectron/F");
  tr->Branch("tower_HF_Plus",              &tower_HF_Plus,       "tower_HF_Plus/F");
  tr->Branch("tower_HF_Minus",              &tower_HF_Minus,       "tower_HF_Minus/F");
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
  Ele_charge = 0;
  
  Mu_Pt= -999 ;
  Mu_Eta= -999 ;
  Mu_Phi= -999 ;
  Mu_E= -999 ;
  Mu_charge = 0;


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
  ok_chexcl_extratracks = 0;
  ok_chexcl_goodtracks = 0;
  ok_trigger = 0;
  ok_singleMuon = 0;
  ok_singleEG3 = 0;
  ok_singleEG5 = 0;
  deltaRtrack_electron = -999;
  deltaRtrack_muon = -999;
  zdc_energy_pos = 0;
  zdc_energy_neg = 0;
  nTracks = -999;
  track_Pt = -999; 
  track_Eta = -999; 
  track_Phi = -999; 
  track_Charge = 0;
  nExtrk = 0;
  vSum_M = -999;
  tower_eta = -999;
  tower_phi = -999;
  tower_em = -999;
  tower_had = -999;
  tower_deltaRmuon_HE = -999;
  tower_deltaRmuon_HB = -999;
  tracks_deltaRElectron = -999;

  deltaRtrack_electron_test = -999;
  deltaRtrack_muon_test = -999;
  deltaRtrack_electron_dEta = -999;
  deltaRtrack_electron_dPhi = -999;
  tower_HF_Plus = -999;
  tower_HF_Minus = -999;
}


/// Checks that number of arguments provided is correct and sets corresponding variables
void ReadInputArguments(int argc, char* argv[],
                        string &configPath, string &inputPath, string &outputPath, string &sampleName)
{
  if(argc != 5){
    Log(0)<<"This app requires 4 parameters:\n";
    Log(0)<<"./selectEleMuEvents configPath inputPath outputPath datasetName[Data|QED_SC|QED_SL|LbL|CEP|GAMMA_UPC|MUMU_FSR]\n";
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
  TTree *genTree = new TTree("gen_tree","");
  TTree *tr = new TTree("output_tree","");
  InitTree(tr);
  InitGenTree(genTree);

  
  auto events = make_unique<EventProcessor>(inputPath, dataset);
  
  int trigger_passed=0, one_goodElectron=0, one_goodElectron_OneGoodMuon=0,oppositelyCharged_ElectronAndMuon=0,  twoGoodGeneralTracks=0, neutral_excl=0, acop=0;
  
  
  // Loop over events
  for(int iEvent=0; iEvent<events->GetNevents(); iEvent++){
  //if(iEvent > 10000) continue; 
   //if(iEvent%10000 == 0) Log(0)<<"Processing event "<<iEvent<<"\n";
    //if(iEvent >= config.params("maxEvents")) break;
   
    auto event = events->GetEvent(iEvent);
    ResetVars();
    //cout <<"After trigger "<<iEvent<<"\n";
    run = event->GetRunNumber();
    ls = event->GetLumiSection();
    evtnb = event->GetEventNumber();
////////
  //  if(event->GetPhysObjects(EPhysObjType::kGoodElectron).size() == 1 & (event->GetPhysObjects(EPhysObjType::kGoodMuon).size() == 1)){
    
 
    if(!(event->HasTrigger(kSingleMuOpenNoHF) || event->HasTrigger(kSingleEG5noHF))) continue;    
    ok_singleMuon = (event->HasTrigger(kSingleMuOpenNoHF));
    ok_singleEG3 = (event->HasTrigger(kSingleEG3noHF));
    ok_singleEG5 = (event->HasTrigger(kSingleEG5noHF));
    
   hist->SetBinContent(1,trigger_passed);  
     //cout << "Event number:" << iEvent << endl;
    trigger_passed++;
    

    if(event->GetPhysObjects(EPhysObjType::kGoodElectron).size() != 1) continue;
    one_goodElectron++;
    hist->SetBinContent(2,one_goodElectron);

    if(event->GetPhysObjects(EPhysObjType::kGoodMuon).size() != 1) continue;
    one_goodElectron_OneGoodMuon++;
    hist->SetBinContent(3,one_goodElectron_OneGoodMuon);
    //

    if(event->GetPhysObjects(EPhysObjType::kGoodMuon)[0]->GetCharge() ==  event->GetPhysObjects(EPhysObjType::kGoodElectron)[0]->GetCharge()) continue;
    oppositelyCharged_ElectronAndMuon++;
    hist->SetBinContent(4,oppositelyCharged_ElectronAndMuon);
 //   cout << "e charge:" << event->GetPhysObjects(EPhysObjType::kGoodElectron)[0]->GetCharge() << "mu charge:" << event->GetPhysObjects(EPhysObjType::kGoodMuon)[0]->GetCharge() << endl;
    //if(event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack).size() !=2) continue;
    twoGoodGeneralTracks++;
    hist->SetBinContent(5,twoGoodGeneralTracks);
    
   // if(event->HasAdditionalTowers()) continue;
    neutral_excl++;
    hist->SetBinContent(6,neutral_excl);




    auto genTracks = event->GetPhysObjects(EPhysObjType::kGeneralTrack);
    auto goodgenTracks = event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack);
    auto goodElectron = event->GetPhysObjects(EPhysObjType::kGoodElectron)[0];
    auto goodMuon = event->GetPhysObjects(EPhysObjType::kGoodMuon)[0];
  //  auto goodgenTracks1 = event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack)[0];
   // auto goodgenTracks2 = event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack)[1];
 

    //if(EleMu_acop < 0.01) continue;
    cout << "Event number:" << iEvent << endl;
   // extra_track_info
    int nextratracks =0;
    for (auto trk : goodgenTracks) {
     track_Pt = trk->GetPt(); 
     track_Eta = trk->GetEta(); 
     track_Phi = trk->GetPhi(); 
     track_Charge = trk->GetCharge(); 
     tracks_deltaRElectron = getDR(track_Eta,track_Phi,goodElectron->GetEta(),goodElectron->GetPhi());
 
     if (getDPHI(trk->GetPhi(), goodElectron->GetPhi())<0.7 && getDETA(trk->GetEta(), goodElectron->GetEta())<0.15) continue;
     if (getDR(trk->GetEta(),trk->GetPhi(), goodMuon->GetEta(),goodMuon->GetPhi()) < 0.1) continue; // For deltaR 
     //tr->Fill(); 
      nextratracks++;
    }
   if(nextratracks>0)Log(0) << "Event :" << iEvent << " extra tracks:" << nextratracks << "\n";    
    nExtrk = nextratracks ;
    ok_chexcl_extratracks = (nExtrk==0);  
    
   // if(ok_chexcl_extratracks!=1) continue;

       nTracks = goodgenTracks.size(); 
     auto trks = event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack);
    for(auto trk: trks){
     double goodgenTracks_Pt = trk->GetPt();
    double goodgenTracks_Eta = trk->GetEta();
    double goodgenTracks_Phi = trk->GetPhi();
    deltaRtrack_electron = getDR(goodgenTracks_Eta,goodgenTracks_Phi,goodElectron->GetEta(),goodElectron->GetPhi());
    deltaRtrack_muon = getDR(goodgenTracks_Eta,goodgenTracks_Phi,goodMuon->GetEta(),goodMuon->GetPhi());
    deltaRtrack_electron_test = getDR(goodgenTracks_Eta,goodgenTracks_Phi,goodElectron->GetEta(),goodElectron->GetPhi());
    deltaRtrack_muon_test = getDR(goodgenTracks_Eta,goodgenTracks_Phi,goodMuon->GetEta(),goodMuon->GetPhi());
    deltaRtrack_electron_dEta = getDETA(trk->GetEta(), goodElectron->GetEta());
    deltaRtrack_electron_dPhi = getDPHI(trk->GetPhi(), goodElectron->GetPhi());
     genTree->Fill();
    }
   ///To Check same sign electron and muon
    ///
    ok_neuexcl = (!event->HasAdditionalTowers());
    ok_chexcl_goodtracks = (goodgenTracks.size()==2);
    if(sampleName == "Data"){
    zdc_energy_pos = event->GetTotalZDCenergyPos(); zdc_energy_neg = event->GetTotalZDCenergyNeg();
   } 



   int ntowers = 0;
   auto towers = event->GetPhysObjects(EPhysObjType::kCaloTower);
   double max_tower_had = -1.0;
   double max_tower_hf_plus = -1.0;
   double max_tower_hf_minus = -1.0;
   for(auto tow: towers){
     //cout << "Tower eta:" << tow->GetEta() << ", Tower energy EM: " << tow->GetEnergyEm() << ", Tower energy Had:" << tow->GetEnergyHad() << endl;
     tower_eta = tow->GetEta();
     tower_phi = tow->GetPhi();
     tower_em = tow->GetEnergyEm();
    double  tower_hf_plus = tow->GetEnergy();
    double  tower_hf_minus = tow->GetEnergy();
     //tower_had = tow->GetEnergyHad();
     if( fabs(tower_eta) > 1.41 && fabs(tower_eta) < 3 && tow->GetEnergyHad()> 0.0 ){
     //cout  <<  "Event no:" << iEvent <<" Had energy HE: "<< tower_had <<  ":deltaR muon and Had: " <<  getDR(tower_eta,tower_phi,goodMuon->GetEta(),goodMuon->GetPhi()) << endl;
     //tower_deltaRmuon_HE = getDR(tower_eta,tower_phi,goodMuon->GetEta(),goodMuon->GetPhi());
     double current_tower_had = tow->GetEnergyHad();
     if(current_tower_had > max_tower_had) { 
      max_tower_had = current_tower_had;
      tower_had = max_tower_had;  
      tower_deltaRmuon_HE = getDR(tower_eta,tower_phi,goodMuon->GetEta(),goodMuon->GetPhi());
     // cout  << " Had energy HE: "<< tower_had <<  ":deltaR muon and Had: " <<  getDR(tower_eta,tower_phi,goodMuon->GetEta(),goodMuon->GetPhi()) << endl;
     }
     }

     if(fabs(tower_eta) > 0 && fabs(tower_eta) < 1.305 && tower_had>2.8){
     tower_deltaRmuon_HB = getDR(tower_eta,tower_phi,goodMuon->GetEta(),goodMuon->GetPhi());
    //cout  << " Had energy HB: "<< tower_had <<  ":deltaR muon and Had: " <<  getDR(tower_eta,tower_phi,goodMuon->GetEta(),goodMuon->GetPhi()) << endl;
     }


    // For HF towers
    if( (tower_eta) > 3.15 && (tower_eta) < 5.2 && tow->GetEnergy()> 0.0 ){
    
     double current_tower_hf_plus = tow->GetEnergy();
     if(current_tower_hf_plus> max_tower_hf_plus) { 
      max_tower_hf_plus = current_tower_hf_plus;
      tower_hf_plus = max_tower_hf_plus;  
      tower_HF_Plus = tower_hf_plus;
      
     }
      
     }
     //cout << "HF Plus:" << tower_HF_Plus << endl;
     // HF minus towers
     if( (tower_eta) > - 5.2 && (tower_eta) < -3.15 && tow->GetEnergy()> 0.0 ){
    
     double current_tower_hf_minus = tow->GetEnergy();
     if(current_tower_hf_minus> max_tower_hf_minus) { 
      max_tower_hf_minus = current_tower_hf_minus;
      tower_hf_minus = max_tower_hf_minus;  
      tower_HF_Minus = tower_hf_minus;
      
     }
     
     }
     //cout << "HF Minus:" << tower_HF_Minus << endl;
     


    ntowers++;
      //tr->Fill(); 
    }



    
     // start filling EleMu information here ........................................
  //  cout << "Ele_Pt:" << goodElectron->GetPt() << endl; 
   
    Ele_Pt       = goodElectron->GetPt();
    Ele_Eta      = goodElectron->GetEta();
    Ele_Phi      = goodElectron->GetPhi();
    Ele_E        = goodElectron->GetEnergy();
    Ele_charge   = goodElectron->GetCharge();
       
    Mu_Pt        = goodMuon->GetPt();
    Mu_Eta       = goodMuon->GetEta();
    Mu_Phi       = goodMuon->GetPhi();
    Mu_E         = goodMuon->GetEnergy();
    Mu_charge    = goodMuon->GetCharge();
    ///
    //cout << "Muon Eta:" << Mu_Eta << endl;
    double eleMass = 0.5109989461e-3;
    double muMass = 105.6583755e-3;
    TLorentzVector ele, mu, EleMu;
    ele.SetPtEtaPhiM(Ele_Pt,Ele_Eta,Ele_Phi,eleMass);
    mu.SetPtEtaPhiM(Mu_Pt,Mu_Eta,Mu_Phi,muMass);
    EleMu = ele + mu;
    vSum_M = EleMu.M();
    vSum_Pt = EleMu.Pt();
    vSum_Phi = EleMu.Phi();    
    Scalar_Pt = (Ele_Pt + Mu_Pt);    
    vSum_Rapidity = EleMu.Rapidity();
    EleMu_dphi = getDPHI(Ele_Phi,Mu_Phi); 
    EleMu_acop = 1 - (EleMu_dphi/3.141592653589);     
    vSum_Pz = EleMu.Pz();
     
 
   /*
    double goodgenTracks1_Pt = goodgenTracks1->GetPt();
    double goodgenTracks1_Eta = goodgenTracks1->GetEta();
    double goodgenTracks1_Phi = goodgenTracks1->GetPhi();
    double goodgenTracks2_Pt = goodgenTracks2->GetPt();
    double goodgenTracks2_Eta = goodgenTracks2->GetEta();
    double goodgenTracks2_Phi = goodgenTracks2->GetPhi();
    


   if(goodgenTracks1->GetCharge() == goodElectron->GetCharge()){
    deltaRtrack_electron = getDR(goodgenTracks1_Eta,goodgenTracks1_Phi,Ele_Eta,Ele_Phi);
  //  cout << "trackcharge:" << goodgenTracks1->GetCharge() << endl;
    }
   else{
    deltaRtrack_muon = getDR(goodgenTracks1_Eta,goodgenTracks1_Phi,Mu_Eta,Mu_Phi);
 }

   if(goodgenTracks2->GetCharge() == goodElectron->GetCharge()){
    deltaRtrack_electron = getDR(goodgenTracks2_Eta,goodgenTracks2_Phi,Ele_Eta,Ele_Phi);
  //   cout << "trackcharge2:" << goodgenTracks2->GetCharge() << endl; 
    }
   else{
    deltaRtrack_muon = getDR(goodgenTracks2_Eta,goodgenTracks2_Phi,Mu_Eta,Mu_Phi);
 //  cout << "trackcharge2:" << goodgenTracks2->GetCharge() << endl; 
  }

*/
 //   if(EleMu_acop > 0.01){
    acop++;
    hist->SetBinContent(7,acop); 
    tr->Fill(); 
  // }
  } //nevents
  Log(0) << "Number of events triggered:" << trigger_passed << "\n" ;
  Log(0) << "One Good electron:" << one_goodElectron << "\n";
  Log(0) << "One Good electron and one good muon:" << one_goodElectron_OneGoodMuon << "\n";
  Log(0) << "Oppositely charged electron and muon:" << oppositelyCharged_ElectronAndMuon << "\n";
  Log(0) << "Two good general tracks:" << twoGoodGeneralTracks << "\n";
  Log(0) << "Neutral exclusivity:" << neutral_excl << "\n";
  
  outFile->cd();
  hist->Write();
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


