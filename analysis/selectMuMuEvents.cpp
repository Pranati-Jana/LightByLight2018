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
float  Mu_Pt1;
float  Mu_Eta1;
float  Mu_Phi1;
float  Mu_E1;
int Mu_charge1;

float  Mu_Pt2;
float  Mu_Eta2;
float  Mu_Phi2;
float  Mu_E2;
int Mu_charge2;

float vSum_M;
//float vSum_M;
float vSum_Energy;
float vSum_Pt;
float Scalar_Pt;
float vSum_Eta;
float vSum_Phi;
float vSum_Rapidity;
float MuMu_dphi;
float MuMu_acop;
int ok_neuexcl;
float E_gamma1;
float E_gamma2;
int ok_chexcl;
int ok_chexcl_extratracks;
int ok_chexcl_goodtracks;
int ok_trigger;
int ok_trigger_noHFveto;
int   nTracks;
int ok_singleMuon;
int ok_singleMuon_noHFveto;
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
float tower_deltaRmuon_HE1;
float tower_deltaRmuon_HE2;
float tower_deltaRmuon_HB;
float tracks_deltaRElectron;
float tower_had_test;
float tower_deltaRmuon_HE1_test;
float tower_deltaRmuon_HE2_test;
/// initialise tree
void InitGenTree(TTree *genTree) {
 genTree->Branch("tower_had_test",              &tower_had_test,       "tower_had_test/F");
   genTree->Branch("tower_deltaRmuon_HE1_test",              &tower_deltaRmuon_HE1_test,       "tower_deltaRmuon_HE1_test/F");
  genTree->Branch("tower_deltaRmuon_HE2_test",              &tower_deltaRmuon_HE2_test,       "tower_deltaRmuon_HE2_test/F");
}
void InitTree(TTree *tr) {
  //
  tr->Branch("run",                 &run,           "run/I");
  tr->Branch("ls",                  &ls,            "ls/I");
  tr->Branch("evtnb",               &evtnb,         "evtnb/I");
  tr->Branch("Mu_Pt1",              &Mu_Pt1,        "Mu_Pt1/F");
  tr->Branch("Mu_Eta1",             &Mu_Eta1,       "Mu_Eta1/F");
  tr->Branch("Mu_Phi1",             &Mu_Phi1,       "Mu_Phi1/F");
  tr->Branch("Mu_E1",               &Mu_E1,         "Mu_E1/F");
  tr->Branch("Mu_charge1",          &Mu_charge1,    "Mu_charge1/I");
 
  tr->Branch("Mu_Pt2",               &Mu_Pt2,        "Mu_Pt2/F");
  tr->Branch("Mu_Eta2",              &Mu_Eta2,       "Mu_Eta2/F");
  tr->Branch("Mu_Phi2",              &Mu_Phi2,       "Mu_Phi2/F");
  tr->Branch("Mu_E2",                &Mu_E2,         "Mu_E2/F");
  tr->Branch("Mu_charge2",           &Mu_charge2,    "Mu_charge2/I");
  
  tr->Branch("vSum_M",              &vSum_M,       "vSum_M/F");
 // tr->Branch("vSum_Energy",         &vSum_diPho_Energy,  "vSum_diPho_Energy/F");
  tr->Branch("vSum_Pt",             &vSum_Pt,      "vSum_Pt/F");
  tr->Branch("Scalar_Pt",           &Scalar_Pt,      "Scalar_Pt/F");
  tr->Branch("vSum_Eta",            &vSum_Eta,     "vSum_Eta/F");
  tr->Branch("vSum_Phi",            &vSum_Phi,     "vSum_Phi/F");
  tr->Branch("vSum_Rapidity",       &vSum_Rapidity,"vSum_Rapidity/F");
  tr->Branch("MuMu_dphi",          &MuMu_dphi,   "MuMu_dphi/F");
  tr->Branch("MuMu_acop",          &MuMu_acop,   "MuMu_acop/F");
  tr->Branch("ok_neuexcl",          &ok_neuexcl,   "ok_neuexcl/I");
  tr->Branch("E_gamma1",          &E_gamma1,   "E_gamma1/I");
  tr->Branch("E_gamma2",          &E_gamma2,   "E_gamma2/I");

  tr->Branch("ok_chexcl",           &ok_chexcl,       "ok_chexcl/I");
  tr->Branch("ok_chexcl_extratracks",           &ok_chexcl_extratracks,       "ok_chexcl_extratracks/I");
  tr->Branch("ok_chexcl_goodtracks",           &ok_chexcl_goodtracks,       "ok_chexcl_goodtracks/I");
  tr->Branch("ok_trigger",                 &ok_trigger, "ok_trigger/I");
  tr->Branch("ok_singleMuon",                 &ok_singleMuon, "ok_singleMuon/I");
  tr->Branch("ok_singleMuon_noHFveto",                 &ok_singleMuon_noHFveto, "ok_singleMuon_noHFveto/I");
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
  tr->Branch("tower_deltaRmuon_HE1",              &tower_deltaRmuon_HE1,       "tower_deltaRmuon_HE1/F");
  tr->Branch("tower_deltaRmuon_HE2",              &tower_deltaRmuon_HE2,       "tower_deltaRmuon_HE2/F");
  tr->Branch("tower_deltaRmuon_HB",              &tower_deltaRmuon_HB,       "tower_deltaRmuon_HB/F");
  tr->Branch("tracks_deltaRElectron",              &tracks_deltaRElectron,      "tracks_deltaRElectron/F");
}


void ResetVars() {
  run =0;
  ls=0;
  evtnb =0;
  //nMu = 0;
  Mu_Pt1= -999 ;
  Mu_Eta1= -999 ;
  Mu_Phi1= -999 ;
  Mu_E1= -999 ;
  Mu_charge1 = 0;
  
  Mu_Pt2= -999 ;
  Mu_Eta2= -999 ;
  Mu_Phi2= -999 ;
  Mu_E2= -999 ;
  Mu_charge2 = 0;


  vSum_M = 0;
  vSum_Energy = 0 ;
  vSum_Pt = 0 ;
  Scalar_Pt = 0;
  vSum_Eta = 0 ;
  vSum_Phi = 0 ;
  vSum_Rapidity = 0 ;
  MuMu_dphi = 0;
  MuMu_acop = 0;
  ok_neuexcl = 0;
  E_gamma1 = -9999;
  E_gamma2 = -9999;

  ok_chexcl = 0;
  ok_chexcl_extratracks = 0;
  ok_chexcl_goodtracks = 0;
  ok_trigger = 0;
  ok_singleMuon = 0;
  ok_singleMuon_noHFveto = 0;
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
  tower_deltaRmuon_HE1 = -999;
  tower_deltaRmuon_HE2 = -999;
  tower_deltaRmuon_HB = -999;
  tracks_deltaRElectron = -999;
  tower_had_test = -999;
  tower_deltaRmuon_HE1_test = -999;
  tower_deltaRmuon_HE2_test = -999;
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
   //if(iEvent > 30) continue; 
   //if(iEvent%10 == 0) Log(0)<<"Processing event "<<iEvent<<"\n";
    //if(iEvent >= config.params("maxEvents")) break;
   
    auto event = events->GetEvent(iEvent);
    ResetVars();
    //cout <<"After trigger "<<iEvent<<"\n";
    run = event->GetRunNumber();
    ls = event->GetLumiSection();
    evtnb = event->GetEventNumber();
////////
  //  if(event->GetPhysObjects(EPhysObjType::kGoodElectron).size() == 1 & (event->GetPhysObjects(EPhysObjType::kGoodMuon).size() == 1)){
    
/*    if(event->GetPhysObjects(EPhysObjType::kGoodElectron).size() != 1) continue;
    cout << "Test3:" << endl;
    if(event->GetPhysObjects(EPhysObjType::kGoodMuon).size() != 1) continue;
    if(event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack).size() !=2) continue;
    if(event->HasAdditionalTowers()) continue;
 */  
    if(!(event->HasTrigger(kSingleMuOpenNoHF))) continue;    
   // if(!(event->HasTrigger(kSingleMuOpenNoHFveto))) continue;    //To calculate HF veto part of singlemuon trigger
    ok_singleMuon_noHFveto = (event->HasTrigger(kSingleMuOpenNoHFveto)); //To calculate HF veto part of singlemuon trigger
    ok_singleMuon = (event->HasTrigger(kSingleMuOpenNoHF));
   
//    cout << "ok_singleMuon_noHFveto:" << ok_singleMuon_noHFveto << endl;   
  //  cout << "ok_singleMuon:" << ok_singleMuon << endl; 
  //  cout << "ok_singleEG3:" << ok_singleEG3 << endl;    
  //  cout << "ok_singleEG5:" << ok_singleEG5 << endl;    
    trigger_passed++;
    

    if(event->GetPhysObjects(EPhysObjType::kGoodMuon).size() != 2) continue;
    one_goodElectron++;
    hist->SetBinContent(2,one_goodElectron);

   // if(event->GetPhysObjects(EPhysObjType::kGoodMuon).size() != 1) continue;
    one_goodElectron_OneGoodMuon++;
    hist->SetBinContent(3,one_goodElectron_OneGoodMuon);
    //
    auto Photons = event->GetPhysObjects(EPhysObjType::kGoodPhoton);
    if(Photons.size()>0) continue;
    auto genTracks = event->GetPhysObjects(EPhysObjType::kGeneralTrack);
    auto goodgenTracks = event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack);
    auto goodMuon1 = event->GetPhysObjects(EPhysObjType::kGoodMuon)[0];
    auto goodMuon2 = event->GetPhysObjects(EPhysObjType::kGoodMuon)[1];
    if(goodMuon1->GetCharge() ==  goodMuon2->GetCharge()) continue;
    oppositelyCharged_ElectronAndMuon++;
    hist->SetBinContent(4,oppositelyCharged_ElectronAndMuon);
 //   cout << "e charge:" << event->GetPhysObjects(EPhysObjType::kGoodElectron)[0]->GetCharge() << "mu charge:" << event->GetPhysObjects(EPhysObjType::kGoodMuon)[0]->GetCharge() << endl;
    if(event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack).size() !=2) continue;
    twoGoodGeneralTracks++;
    hist->SetBinContent(5,twoGoodGeneralTracks);
     
    //if(event->HasAdditionalTowers()) continue;
    
    hist->SetBinContent(6,neutral_excl);

   //Tower energy and eta
   hist->SetBinContent(1,trigger_passed);  
     //cout << "Event number:" << iEvent << endl;

   

    nTracks = goodgenTracks.size(); 
   // extra_track_info
    int nextratracks =0;
    for (auto trk : goodgenTracks) {
     track_Pt = trk->GetPt(); 
     track_Eta = trk->GetEta(); 
     track_Phi = trk->GetPhi(); 
     track_Charge = trk->GetCharge(); 
    // tracks_deltaRElectron = getDR(track_Eta,track_Phi,goodElectron->GetEta(),goodElectron->GetPhi());
 
    //  if (getDPHI(trk->GetPhi(), goodElectron->GetPhi())<0.7 && getDETA(trk->GetEta(), goodElectron->GetEta())<0.15) continue;
    //  if (getDR(trk->GetEta(),trk->GetPhi(), goodMuon->GetEta(),goodMuon->GetPhi()) < 0.001) continue; // For deltaR 
     //tr->Fill(); 
      nextratracks++;
    }
   if(nextratracks>0)Log(0) << "Event :" << iEvent << " extra tracks:" << nextratracks << "\n";    
    nExtrk = nextratracks ;
    ok_chexcl_extratracks = (nExtrk==0);  
   ///To Check same sign electron and muon
    
    ok_neuexcl = (!event->HasAdditionalTowers());
    //if(goodgenTracks.size() !=2) continue;
    ok_chexcl_goodtracks = (goodgenTracks.size()==2);
    if(sampleName == "Data"){
    zdc_energy_pos = event->GetTotalZDCenergyPos(); zdc_energy_neg = event->GetTotalZDCenergyNeg();
   } 
   ///////////
       Mu_Pt1       = goodMuon1->GetPt();
    Mu_Eta1      = goodMuon1->GetEta();
    Mu_Phi1      = goodMuon1->GetPhi();
    Mu_E1        = goodMuon1->GetEnergy();
    Mu_charge1   = goodMuon1->GetCharge();    
    Mu_Pt2        = goodMuon2->GetPt();
    Mu_Eta2       = goodMuon2->GetEta();
    Mu_Phi2       = goodMuon2->GetPhi();
    Mu_E2         = goodMuon2->GetEnergy();
    Mu_charge2    = goodMuon2->GetCharge();
    ///
    //cout << "Muon Eta:" << Mu_Eta1 << ";"<< Mu_Eta2  << endl;
    double eleMass = 0.5109989461e-3;
    double muMass = 105.6583755e-3;
    TLorentzVector mu1, mu2, MuMu;
    mu1.SetPtEtaPhiM(Mu_Pt1,Mu_Eta1,Mu_Phi1,muMass);
    mu2.SetPtEtaPhiM(Mu_Pt2,Mu_Eta2,Mu_Phi2,muMass);
    MuMu = mu1 + mu2;
    vSum_M = MuMu.M();
   // if(vSum_M<8) continue;
   // if(vSum_M>60) continue;
    vSum_Pt = MuMu.Pt();
    vSum_Phi = MuMu.Phi();    
    Scalar_Pt = (Mu_Pt1 + Mu_Pt2);    
    vSum_Rapidity = MuMu.Rapidity();
    MuMu_dphi = getDPHI(Mu_Phi1,Mu_Phi2); 
    MuMu_acop = 1 - (MuMu_dphi/3.141592653589);     
    vSum_Pz = MuMu.Pz();
    E_gamma1 = 0.5*vSum_M*exp(vSum_Rapidity);
    E_gamma2 = 0.5*vSum_M*exp((-1)*vSum_Rapidity);
    neutral_excl++;
    /////////////////
   int ntowers = 0;
   auto towers = event->GetPhysObjects(EPhysObjType::kCaloTower);

   double max_tower_had = -1.0;
   for(auto tow: towers){
      // cout << "tower eta:" <<  tow->GetEta() << endl;
 
   }

 /*  for(auto tow: towers){
    
     tower_eta = tow->GetEta();
     tower_phi = tow->GetPhi();
     tower_em = tow->GetEnergyEm();
     if( fabs(tower_eta) > 1.41 && fabs(tower_eta) < 3 && tow->GetEnergyHad()>0.0 ){
    */
	     /**********Any HE tower > 1.0GeV*/
      /*
      tower_had_test = tow->GetEnergyHad();
      tower_had = tow->GetEnergyHad();
      cout  << " Had energy HE1: "<< tower_had <<  ":deltaR muon and Had: " <<  getDR(tower_eta,tower_phi,goodMuon1->GetEta(),goodMuon1->GetPhi()) << endl;
     cout  << " Had energy HE2: "<< tower_had <<  ":deltaR muon and Had: " <<  getDR(tower_eta,tower_phi,goodMuon2->GetEta(),goodMuon2->GetPhi()) << endl;
     tower_deltaRmuon_HE1 = getDR(tower_eta,tower_phi,goodMuon1->GetEta(),goodMuon1->GetPhi());
     tower_deltaRmuon_HE2 = getDR(tower_eta,tower_phi,goodMuon2->GetEta(),goodMuon2->GetPhi());
     tower_deltaRmuon_HE1_test = getDR(tower_eta,tower_phi,goodMuon1->GetEta(),goodMuon1->GetPhi());
     tower_deltaRmuon_HE2_test = getDR(tower_eta,tower_phi,goodMuon2->GetEta(),goodMuon2->GetPhi());
     genTree->Fill();
     */
    /********************Leading HE towers******************/
  /*   double current_tower_had = tow->GetEnergyHad();
     if(current_tower_had > max_tower_had) { 
      max_tower_had = current_tower_had;
      tower_had = max_tower_had;  
      tower_had_test =  max_tower_had;
      tower_deltaRmuon_HE1 = getDR(tower_eta,tower_phi,goodMuon1->GetEta(),goodMuon1->GetPhi());
      tower_deltaRmuon_HE2 = getDR(tower_eta,tower_phi,goodMuon2->GetEta(),goodMuon2->GetPhi());
  
//        cout  << " Had energy HE1: "<< tower_had <<  ":deltaR muon and Had: " <<  getDR(tower_eta,tower_phi,goodMuon1->GetEta(),goodMuon1->GetPhi()) << endl;
  //    cout  << " Had energy HE2: "<< tower_had <<  ":deltaR muon and Had: " <<  getDR(tower_eta,tower_phi,goodMuon2->GetEta(),goodMuon2->GetPhi()) << endl;
     }
    

     }
   

    //  if(fabs(tower_eta) > 0 && fabs(tower_eta) < 1.305 && tow->GetEnergyHad()>0.0){
    //  //tower_deltaRmuon_HB = getDR(tower_eta,tower_phi,goodMuon1->GetEta(),goodMuon1->GetPhi());
    // //cout  << " Had energy HB: "<< tower_had <<  ":deltaR muon and Had: " <<  getDR(tower_eta,tower_phi,goodMuon1->GetEta(),goodMuon1->GetPhi()) << endl;
    //  }
    // ntowers++;    
    } */
    
     // start filling EleMu information here ........................................
  
     
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


