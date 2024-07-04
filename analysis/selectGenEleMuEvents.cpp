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


//Gen_Tau
//float Tau_Pt;
//
int run;
int ls;
int evtnb;
//int    nPho;

//
float Tau_Pt;
float Tau_Eta;
float Tau_Phi;
float Tau_E;
float TauTau_M;

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

void InitGenTree(TTree *genTree) {
 //    genTree->Branch("Tau_Pt",                 &Tau_Pt,           "Tau_Pt/F");    
}

/// initialise tree
void InitTree(TTree *tr) {
  //GenTau
  tr->Branch("Tau_Pt",                 &Tau_Pt,           "Tau_Pt/F");
  tr->Branch("Tau_Eta",                &Tau_Eta,          "Tau_Eta/F");
  tr->Branch("Tau_Phi",                &Tau_Phi,          "Tau_Phi/F");
  tr->Branch("Tau_E",                  &Tau_E,            "Tau_E/F");
  tr->Branch("TauTau_M",               &TauTau_M,         "TauTau_M/F");
  //
  tr->Branch("run",                 &run,           "run/I");
  tr->Branch("ls",                  &ls,            "ls/I");
  tr->Branch("evtnb",               &evtnb,         "evtnb/I");
  //tr->Branch("nPho",                &nPho,          "nPho/I");
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
  


}


// reset all variables
void ResetGenVars() {
//Tau_Pt = 0;

}


void ResetVars() {
  //GenTau
 // Tau_Pt = -999;
  //
  run =0;
  ls=0;
  evtnb =0;
  //nMu = 0;
  Tau_Pt = -999;
  Tau_Eta = -999;
  Tau_Phi = -999;
  Tau_E = -999;
  TauTau_M = -999;

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
}


/// Checks that number of arguments provided is correct and sets corresponding variables
void ReadInputArguments(int argc, char* argv[],
                        string &configPath, string &inputPath, string &outputPath, string &sampleName)
{
  if(argc != 5){
    Log(0)<<"This app requires 4 parameters:\n";
    Log(0)<<"./selectGenEleMuEvents configPath inputPath outputPath datasetName[Data|QED_SC|QED_SL|LbL|CEP]\n";
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
  // InitGenTree(genTree);

  
  auto events = make_unique<EventProcessor>(inputPath, dataset);
  
  int trigger_passed=0, twoPho=0, twoGoodPho=0, oppCharge=0, neutral_excl=0, charged_excl=0, diphomass_wozdc=0, diphopt_wozdc=0;
  int acop_cut_wozdc=0, diphomass=0, diphopt=0, acop_cut = 0, zdc_excl=0;
  int EleMuCount=0; 
  // Loop over events
  for(int iEvent=0; iEvent<events->GetNevents(); iEvent++){
    if(iEvent%1000 == 0) Log(0)<<"Processing event "<<iEvent<<"\n";
    //if(iEvent >= config.params("maxEvents")) break;
    
    auto event = events->GetEvent(iEvent);
  //  ResetGenVars();  
    ResetVars();
      
    auto gen = event->GetPhysObjects(EPhysObjType::kGenParticle);
    cout << "gen size:" << gen.size() << endl;
   // cout << "gen size:" << gen.size() << endl;
    //auto GenTaus_1 = event->GetPhysObjects(EPhysObjType::kGenParticle)[0];
    //auto GenTaus_2 = event->GetPhysObjects(EPhysObjType::kGenParticle)[1];
    auto GenTaus_1 = event->GetPhysObjects(EPhysObjType::kGenParticle)[0];
    auto GenTaus_2 = event->GetPhysObjects(EPhysObjType::kGenParticle)[1];
    //auto GenTaus_5 = event->GetPhysObjects(EPhysObjType::kGenParticle)[4];
    //auto GenTaus_6 = event->GetPhysObjects(EPhysObjType::kGenParticle)[5];
    //cout << "gen PID1:" << GenTaus_1->GetPID() << endl;
    //cout << "gen Pt1:" << GenTaus_1->GetPt() << endl;
    //cout << "gen PID2:" << GenTaus_2->GetPID() << endl;
    //cout << "gen Pt2:" << GenTaus_2->GetPt() << endl;
    //cout << "gen PID1:" << GenTaus_1->GetPID() << endl;
   // cout << "gen PID2:" << GenTaus_2->GetPID() << endl;
    double taumass = 1.777;
    TLorentzVector tau1, tau2, tautau;
    tau1.SetPtEtaPhiM(GenTaus_1->GetPt(),GenTaus_1->GetEta(),GenTaus_1->GetPhi(),taumass);
    tau2.SetPtEtaPhiM(GenTaus_2->GetPt(),GenTaus_2->GetEta(),GenTaus_2->GetPhi(),taumass);
    tautau = tau1 + tau2 ;
    cout << "ditauMass : " << tautau.M();
    TauTau_M = tautau.M();
    auto genDau = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter);
         cout << "genDau size:" << genDau.size() << endl;
   
    auto GenTausD_1 = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter)[0];
    auto GenTausD_2 = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter)[1];
    auto GenTausD_3 = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter)[2];
    auto GenTausD_4 = event->GetPhysObjects(EPhysObjType::kGenParticleTauDaughter)[3];
   
   // cout << "genDau PID1:" << GenTausD_1->GetPID() << endl;
   // cout << "gendau PID2:" << GenTausD_2->GetPID() << endl;
   // cout << "gendau PID3:" << GenTausD_3->GetPID() << endl;
   // cout << "gendau PID4:" << GenTausD_4->GetPID() << endl;
    for(auto Taus : gen){
      Tau_Pt  = Taus->GetPt();
      Tau_Eta = Taus->GetEta();
      Tau_Phi = Taus->GetPhi();
      Tau_E   = Taus->GetEnergy();
      cout << "TauPt:" << Tau_Pt << endl;
       
    }
    bool hasMuon = false;
    bool hasElectron = false;
//    bool EleMuEvents = false;
    for(auto TauD : genDau){
       int  TauD_PID = TauD->GetPID();
       cout << "TauD PID:" << TauD_PID << endl;
       if(abs(TauD_PID) == 13){
        hasMuon = true;
        Gen_Mu_Pt = TauD->GetPt();
        Gen_Mu_Eta = TauD->GetEta();
        Gen_Mu_Phi = TauD->GetPhi();
        Gen_Mu_E = TauD->GetEnergy();
     }else if (abs(TauD_PID) == 11) {
        hasElectron = true;
        Gen_Ele_Pt = TauD->GetPt();
        Gen_Ele_Eta = TauD->GetEta();
        Gen_Ele_Phi = TauD->GetPhi();
        Gen_Ele_E = TauD->GetEnergy();
        }
     }
      EleMuEvents = hasMuon && hasElectron;
     if ( EleMuEvents ) {
        cout << "EleMu Events:::::::::::::::::::::::::::::::::::::::::" << EleMuEvents << endl;
    //  bool EleMuEvents = true; 
  //     EleMuCount++;
        double eleMass = 0.5109989461e-3;
        double muMass = 105.6583755e-3;
        TLorentzVector ele, mu, EleMu;
        ele.SetPtEtaPhiM(Gen_Ele_Pt,Gen_Ele_Eta,Gen_Ele_Phi,eleMass);
        mu.SetPtEtaPhiM(Gen_Mu_Pt,Gen_Mu_Eta,Gen_Mu_Phi,muMass);
        EleMu = ele + mu;
        Gen_vSum_M = EleMu.M();
        Gen_vSum_Pt = EleMu.Pt();
        Gen_vSum_Phi = EleMu.Phi();    
        Gen_Scalar_Pt = (Gen_Ele_Pt + Gen_Mu_Pt);    
        Gen_vSum_Rapidity = EleMu.Rapidity();
        Gen_EleMu_dphi = getDPHI(Gen_Ele_Phi,Gen_Mu_Phi); 
        Gen_EleMu_acop = 1 - (Gen_EleMu_dphi/3.141592653589); 
        cout << "Ele Phi:" << Gen_Ele_Phi << endl;
        cout << "Mu Phi:" << Gen_Mu_Phi << endl;
        cout << "Ele Pt:" << Gen_Ele_Pt << endl;
        cout << "Mu Pt:" << Gen_Mu_Pt << endl;
        cout << "Ele Eta:" << Gen_Ele_Eta << endl;
        cout << "Mu Eta:" << Gen_Mu_Eta << endl;

 
        cout << " Gen_EleMu_acop : " <<  Gen_EleMu_acop << endl;
        cout << " Gen_EleMu_M : " <<  Gen_vSum_M << endl;
 
    } 
     /*
    //GenTau
    if(event->GetPhysObjects(EPhysObjType::kGenParticleTau).size() == 2); {
    auto GenTaus_1 = event->GetPhysObjects(EPhysObjType::kGenParticleTau)[0];
    auto GenTaus_2 = event->GetPhysObjects(EPhysObjType::kGenParticleTau)[0];
    ResetGenVars(); 
    for(auto Taus : GenTaus){
       Tau_Pt = GenTaus_1->GetPt();
    }
    }
    */
    //cout <<"After trigger "<<iEvent<<"\n";
    //ok_trigger = ((event->HasTrigger(kSingleMuOpenNoHF)) or (event->HasTrigger(kSingleEG5noHF)));     
   // cout << "Test1:" << endl; 
   // run = event->GetRunNumber();
  //  ls = event->GetLumiSection();
  //  evtnb = event->GetEventNumber();
   // cout << "Test2:" << endl;

   // if( (0.01 < Gen_EleMu_acop < 0.4) && Gen_Ele_Pt > 2 && Gen_Mu_Pt > 2.5 && abs(Gen_Ele_Eta < 2.2) && abs(Gen_Mu_Eta < 2.4)){
      
    
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
