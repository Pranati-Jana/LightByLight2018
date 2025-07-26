//  SelectQEDEvents
//
//  Created by Ruchi Chudasama on 24/12/2020
//
//  Selects QED events in MC and data, as defined in the Analysis Note

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
double cosphotonpair(TLorentzVector p1, TLorentzVector pair, bool helicityFrame );
double costhetastar_CS(TLorentzVector p1, TLorentzVector pair); 
float  GetAngle( float vPX1, float vPY1, float vPX2, float vPY2);

float gen_Pt1;
float gen_Eta1;
float gen_Phi1;
float gen_PID1;
float gen_Pt2;
float gen_Eta2;
float gen_Phi2;
float gen_PID2;
float gen_diEle_M;
float gen_diEle_Pt;
float gen_diEle_Rapidity;
float gen_diEle_Costhetastar;
float gen_diEle_PhiModulation;
int run;
int ls;
int evtnb;
int    nEle;
int    eleCharge_1;
float  elePt_1;
float  eleEta_1;
float  elePhi_1;
float  eleSCEt_1;
float  eleSCEta_1;
float  eleSCPhi_1;
float  eleEoverP_1;

int    eleCharge_2;
float  elePt_2;
float  eleEta_2;
float  elePhi_2;
float  eleSCEt_2;
float  eleSCEta_2;
float  eleSCPhi_2;
float  eleEoverP_2;

float vSum_ee_M;
float vSum_ee_Energy;
float vSum_ee_Pt;
float vSum_ee_Eta;
float vSum_ee_Phi;
float vSum_ee_Rapidity;
float ele_PhiModulation;
float ele_dpt;
float ele_deta;
float ele_dphi;
float ele_acop;
int ok_neuexcl;
int ok_zdcexcl;
int ok_chexcl;
int ok_chexcl_extrk;
int nExtrk;
int ok_trigger;
int ok_event;

float zdc_energy_pos;
float zdc_energy_neg;

float costhetastar;
float cos_photon_pair_helicity0;
float cos_photon_pair_helicity1;

/// initialise gen tree
void InitGenTree(TTree *genTree) {
  genTree->Branch("gen_Pt1",  &gen_Pt1, "gen_Pt1/F");
  genTree->Branch("gen_Eta1", &gen_Eta1, "gen_Eta1/F");
  genTree->Branch("gen_Phi1", &gen_Phi1, "gen_Phi1/F");
  genTree->Branch("gen_PID1", &gen_PID1, "gen_PID1/F");

  genTree->Branch("gen_Pt2",  &gen_Pt2, "gen_Pt2/F");
  genTree->Branch("gen_Eta2", &gen_Eta2, "gen_Eta2/F");
  genTree->Branch("gen_Phi2", &gen_Phi2, "gen_Phi2/F");
  genTree->Branch("gen_PID2", &gen_PID2, "gen_PID2/F");

  genTree->Branch("gen_diEle_M",  &gen_diEle_M,  "gen_diEle_M/F");
  genTree->Branch("gen_diEle_Pt", &gen_diEle_Pt, "gen_diEle_Pt/F");
  genTree->Branch("gen_diEle_Rapidity", &gen_diEle_Rapidity, "gen_diEle_Rapidity/F");
  genTree->Branch("gen_diEle_Costhetastar", &gen_diEle_Costhetastar, "gen_diEle_Costhetastar/F");
  genTree->Branch("gen_diEle_PhiModulation", &gen_diEle_PhiModulation, "gen_diEle_PhiModulation/F");
}


/// initialise tree
void InitTree(TTree *tr) {

  tr->Branch("gen_Pt1",  &gen_Pt1, "gen_Pt1/F");
  tr->Branch("gen_Eta1", &gen_Eta1, "gen_Eta1/F");
  tr->Branch("gen_Phi1", &gen_Phi1, "gen_Phi1/F");
  tr->Branch("gen_PID1", &gen_PID1, "gen_PID1/F");

  tr->Branch("gen_Pt2",  &gen_Pt2, "gen_Pt2/F");
  tr->Branch("gen_Eta2", &gen_Eta2, "gen_Eta2/F");
  tr->Branch("gen_Phi2", &gen_Phi2, "gen_Phi2/F");
  tr->Branch("gen_PID2", &gen_PID2, "gen_PID2/F");
  
  tr->Branch("gen_diEle_M",  &gen_diEle_M,  "gen_diEle_M/F");
  tr->Branch("gen_diEle_Pt", &gen_diEle_Pt, "gen_diEle_Pt/F");
  tr->Branch("gen_diEle_Rapidity", &gen_diEle_Rapidity, "gen_diEle_Rapidity/F");
  tr->Branch("gen_diEle_Costhetastar", &gen_diEle_Costhetastar, "gen_diEle_Costhetastar/F");
  tr->Branch("gen_diEle_PhiModulation", &gen_diEle_PhiModulation, "gen_diEle_PhiModulation/F");
  
  tr->Branch("run",                 &run,           "run/I");
  tr->Branch("ls",                  &ls,            "ls/I");
  tr->Branch("evtnb",               &evtnb,         "evtnb/I");
  tr->Branch("nEle",                &nEle,          "nEle/I");
  tr->Branch("eleCharge_1",         &eleCharge_1,   "eleCharge_1/I");
  tr->Branch("elePt_1",             &elePt_1,       "elePt_1/F");
  tr->Branch("eleEta_1",            &eleEta_1,      "eleEta_1/F");
  tr->Branch("elePhi_1",            &elePhi_1,      "elePhi_1/F");
  tr->Branch("eleSCEt_1",           &eleSCEt_1,     "eleSCEt_1/F");
  tr->Branch("eleSCEta_1",          &eleSCEta_1,    "eleSCEta_1/F");
  tr->Branch("eleSCPhi_1",          &eleSCPhi_1,    "eleSCPhi_1/F");
  tr->Branch("eleEoverP_1",         &eleEoverP_1,   "eleEoverP_1/F");
 
  tr->Branch("eleCharge_2",         &eleCharge_2,   "eleCharge_2/I");
  tr->Branch("elePt_2",             &elePt_2,       "elePt_2/F");
  tr->Branch("eleEta_2",            &eleEta_2,      "eleEta_2/F");
  tr->Branch("elePhi_2",            &elePhi_2,      "elePhi_2/F");
  tr->Branch("eleSCEt_2",           &eleSCEt_2,     "eleSCEt_2/F");
  tr->Branch("eleSCEta_2",          &eleSCEta_2,    "eleSCEta_2/F");
  tr->Branch("eleSCPhi_2",          &eleSCPhi_2,    "eleSCPhi_2/F");
  tr->Branch("eleEoverP_2",         &eleEoverP_2,   "eleEoverP_2/F");
  
  tr->Branch("vSum_M",              &vSum_ee_M,       "vSum_ee_M/F");
  tr->Branch("vSum_Energy",         &vSum_ee_Energy,  "vSum_ee_Energy/F");
  tr->Branch("vSum_Pt",             &vSum_ee_Pt,      "vSum_ee_Pt/F");
  tr->Branch("vSum_Eta",            &vSum_ee_Eta,     "vSum_ee_Eta/F");
  tr->Branch("vSum_Phi",            &vSum_ee_Phi,     "vSum_ee_Phi/F");
  tr->Branch("vSum_Rapidity",       &vSum_ee_Rapidity,"vSum_ee_Rapidity/F");
  tr->Branch("ele_PhiModulation",   &ele_PhiModulation,"ele_PhiModulation/F");
  tr->Branch("ele_dpt",             &ele_dpt,         "ele_dpt/F");
  tr->Branch("ele_deta",            &ele_deta,        "ele_deta/F");
  tr->Branch("ele_dphi",            &ele_dphi,        "ele_dphi/F");
  tr->Branch("ele_acop",            &ele_acop,        "ele_acop/F");
  tr->Branch("ok_neuexcl",          &ok_neuexcl,      "ok_neuexcl/I");
  tr->Branch("ok_zdcexcl",          &ok_zdcexcl,      "ok_zdcexcl/I");
  tr->Branch("ok_chexcl",           &ok_chexcl,       "ok_chexcl/I");
  tr->Branch("ok_chexcl_extrk",     &ok_chexcl_extrk, "ok_chexcl_extrk/I");
  tr->Branch("nExtrk",              &nExtrk,          "nExtrk/I");
  tr->Branch("ok_trigger",          &ok_trigger,      "ok_trigger/I");
  tr->Branch("ok_event",            &ok_event,      "ok_event/I");
  tr->Branch("zdc_energy_pos",             &zdc_energy_pos,         "zdc_energy_pos/F");
  tr->Branch("zdc_energy_neg",             &zdc_energy_neg,         "zdc_energy_neg/F");

  tr->Branch("costhetastar",        &costhetastar,        "costhetastar/F");
  tr->Branch("cos_photon_pair_helicity0",     &cos_photon_pair_helicity0,     "cos_photon_pair_helicity0/F");
  tr->Branch("cos_photon_pair_helicity1",     &cos_photon_pair_helicity1,     "cos_photon_pair_helicity1/F");
  //
}

// reset gen variables
void ResetGenVars() {
  gen_Pt1 = 0;
  gen_Eta1 = 0;
  gen_Phi1 = 0;
  gen_PID1 = 0;
  gen_Pt2 = 0;
  gen_Eta2 = 0;
  gen_Phi2 = 0;
  gen_PID2 = 0;
  gen_diEle_M = 0;
  gen_diEle_Pt = 0;
  gen_diEle_Rapidity = 0;
  gen_diEle_Costhetastar = -999;
  gen_diEle_PhiModulation = -999;
}


// reset all variables
void ResetVars() {
  run =0;
  ls=0;
  evtnb =0;
  nEle = 0;
  eleCharge_1= 0 ;
  elePt_1= 0 ;
  eleEta_1= 0 ;
  elePhi_1= 0 ;
  eleSCEt_1= 0 ;
  eleSCEta_1= 0 ;
  eleSCPhi_1= 0 ;
  eleEoverP_1 = 0 ;
  
  eleCharge_2 = 0 ;
  elePt_2 = 0 ;
  eleEta_2 = 0 ;
  elePhi_2 = 0 ;
  eleSCEt_2 = 0 ;
  eleSCEta_2 = 0 ;
  eleSCPhi_2 = 0 ;
  eleEoverP_2 = 0 ;
 
  vSum_ee_M = 0 ;
  vSum_ee_Energy = 0 ;
  vSum_ee_Pt = 0 ;
  vSum_ee_Eta = 0 ;
  vSum_ee_Phi = 0 ;
  vSum_ee_Rapidity = 0 ;
  ele_PhiModulation = 0;
  ele_dpt = 0;
  ele_deta = 0;
  ele_dphi = 0;
  ele_acop = 0;
  ok_neuexcl = 0;
  ok_zdcexcl = 0;
  ok_chexcl = 0;
  ok_chexcl_extrk = 0;
  nExtrk = 0;
  ok_trigger = 1;
  ok_event = 0;
  zdc_energy_pos = 0;
  zdc_energy_neg = 0;
  
  costhetastar = -999;
  cos_photon_pair_helicity0 = -999;
  cos_photon_pair_helicity1 = -999;
}


/// Checks that number of arguments provided is correct and sets corresponding variables
void ReadInputArguments(int argc, char* argv[],
                        string &configPath, string &inputPath, string &outputPath, string &sampleName)
{
  if(argc != 5){
    Log(0)<<"This app requires 4 parameters:\n";
    Log(0)<<"./selectQEDEvents configPath inputPath outputPath datasetName[Data|QED_SC|QED_SL|LbL|CEP]\n";
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
  TH1D *hist_wozdc = new TH1D("hist_wozdc","",9,1,10);
 
  TFile *outFile = TFile::Open(outputPath.c_str(), "recreate");
  TTree *tr = new TTree("output_tree","");
  InitTree(tr);

  
  auto events = make_unique<EventProcessor>(inputPath, dataset);
  
  int event_passed=0, trigger_passed=0, twoGoodEle=0, oppCharge=0, neutral_excl=0, charged_excl=0, dielemass_wozdc=0, dielept_wozdc=0;
  int acop_cut_wozdc=0, dielemass=0, dielept=0, acop_cut = 0, zdc_excl=0;
  
  // Loop over events
  for(int iEvent=0; iEvent<events->GetNevents(); iEvent++){
//  for(int iEvent=0; iEvent<100; iEvent++){
//if(iEvent > 200) continue;
    if(iEvent%10000 == 0) Log(0)<<"Processing event "<<iEvent<<"\n";
    //if(iEvent >= config.params("maxEvents")) break;
    
    auto event = events->GetEvent(iEvent);
    
    ResetGenVars();  
    ResetVars();  
    
   /* if(sampleName == "QED_SL"){   
      // cout << "PID1:" <<  genP1->GetPID() << "  two:" << genP2->GetPID() << "\n";  
      auto genP1 = event->GetPhysObjects(EPhysObjType::kGenParticle)[0]; //3 for SC
      auto genP2 = event->GetPhysObjects(EPhysObjType::kGenParticle)[1]; // 4 for SC

      //cout << "PID1:" <<  genP1->GetPID() << "  two:" << genP2->GetPID() << "\n";     
      */ 
  if(sampleName == "QED_SC" ){
	auto genP1 = event->GetPhysObjects(EPhysObjType::kGenParticle)[3]; //3 for SC // 0 for SL
	auto genP2 = event->GetPhysObjects(EPhysObjType::kGenParticle)[4]; // 4 for SC // 1 for SL
      
      
      ResetGenVars();
      if(abs(genP1->GetPID())== 11 && abs(genP2->GetPID())== 11 ){
	//Log(0) << " PID one:" << genP1->GetPID() << "  two:" << genP2->GetPID() << "\n";
	gen_Pt1  = genP1->GetEt();
	gen_Eta1 = genP1->GetEta();
	gen_Phi1 = genP1->GetPhi();
	gen_PID1 = genP1->GetPID();
	gen_Pt2  = genP2->GetEt();
	gen_Eta2 = genP2->GetEta();
	gen_Phi2 = genP2->GetPhi();
  gen_PID2 = genP2->GetPID();
      //  cout << "Tesrt:" << genP2->GetEt() << endl;	
	TLorentzVector gen_ele1, gen_ele2, gen_diele, gen_ele;
	gen_ele1.SetPtEtaPhiM(gen_Pt1, gen_Eta1, gen_Phi1, 0.000511);
	gen_ele2.SetPtEtaPhiM(gen_Pt2, gen_Eta2, gen_Phi2, 0.000511);
	
	gen_diele = gen_ele1 + gen_ele2;

	gen_diEle_M   = gen_diele.M();
	gen_diEle_Pt  = gen_diele.Pt();
	gen_diEle_Rapidity = gen_diele.Rapidity();
	gen_diEle_Costhetastar = costhetastar_CS(gen_ele1,gen_diele);

  //momentum components to determine the angle between the two vectors

    if(gen_PID1 == 11){
    gen_ele.SetPtEtaPhiM(gen_Pt1, gen_Eta1, gen_Phi1, 0.000511);
    //cout << "neg charge:" << gen_PID1<< ":gen pt1:" << gen_Pt1 << endl;
   }
   else if (gen_PID1 == -11){
    gen_ele.SetPtEtaPhiM(gen_Pt2, gen_Eta2, gen_Phi2, 0.000511);
    //cout << "neg charge:" << gen_PID2 <<  ":gen pt2:" << gen_Pt2 << endl;
   }
        float GenDPx = gen_diele.Pt()*cos(gen_diele.Phi());
        float GenDPy = gen_diele.Pt()*sin(gen_diele.Phi());
        float GenMPx = gen_ele.Pt()*cos(gen_ele.Phi());
        float GenMPy = gen_ele.Pt()*sin(gen_ele.Phi());
///////////////
	//gen_diEle_PhiModulation = GetAngle( GenDPx, GenDPy, GenMPx, GenMPy );  //Gen delta phi between 0 and pi// for SC
  gen_diEle_PhiModulation = GetAngle( GenDPx, GenDPy, GenMPx, GenMPy ); // For SL
  //cout << " Gen phi modulation:" << gen_diEle_PhiModulation << endl;

	
      } // PID
}
        // Check trigger
      event_passed++;
      hist->SetBinContent(1,event_passed);          
      ok_trigger = (event->HasTrigger(kDoubleEG2noHF));
      
      //trigger_passed++;
          
    //  if(!event->HasTrigger(kDoubleEG2noHF)) continue;
      trigger_passed++;
      hist->SetBinContent(2,trigger_passed);
      hist_wozdc->SetBinContent(1,trigger_passed);
      run = event->GetRunNumber();
      ls = event->GetLumiSection();
      evtnb = event->GetEventNumber();
      
      // select two exclusive electrons for QED event selection
      if(event->GetPhysObjects(EPhysObjType::kGoodElectron).size() == 2){
	twoGoodEle++;
	hist->SetBinContent(3,twoGoodEle);     hist_wozdc->SetBinContent(2,twoGoodEle);
	
	auto genTracks = event->GetPhysObjects(EPhysObjType::kGoodGeneralTrack);
	auto electron1 = event->GetPhysObjects(EPhysObjType::kGoodElectron)[0];
	auto electron2 = event->GetPhysObjects(EPhysObjType::kGoodElectron)[1];
	auto caloTower = event->GetPhysObjects(EPhysObjType::kCaloTower);
  auto muons = event->GetPhysObjects(EPhysObjType::kGoodMuon); 

	//if(electron1->GetCharge() == electron2->GetCharge()) continue;
	if(electron1->GetCharge() != electron2->GetCharge()) {
	oppCharge++;
	hist->SetBinContent(4,oppCharge);      hist_wozdc->SetBinContent(3,oppCharge);
	
	// event variables
	ok_neuexcl = (!event->HasAdditionalTowers());
	ok_chexcl  = (genTracks.size()==2);
        
	if(sampleName == "Data"){
        ok_zdcexcl = event->GetTotalZDCenergyPos() < 10000 && event->GetTotalZDCenergyNeg() < 10000;
	zdc_energy_pos = event->GetTotalZDCenergyPos();
        zdc_energy_neg = event->GetTotalZDCenergyNeg();
        }
	// start filling extra track information here ........................................
	int nextratracks =0;
    for (auto trk : genTracks) {
      if (getDPHI(trk->GetPhi(), electron1->GetPhi())<0.7 && getDETA(trk->GetEta(), electron1->GetEta())<0.15) continue;
      if (getDPHI(trk->GetPhi(), electron2->GetPhi())<0.7 && getDETA(trk->GetEta(), electron2->GetEta())<0.15) continue;
      nextratracks++;
    }

        if(nextratracks>0)Log(0) << "Event :" << iEvent << " extra tracks:" << nextratracks << "\n";
	nExtrk = nextratracks ;
	ok_chexcl_extrk = (nExtrk==0 && muons.size()==0);
	
	eleCharge_1  = electron1->GetCharge();
	elePt_1      = electron1->GetPt();
	eleEta_1     = electron1->GetEta();
	elePhi_1     = electron1->GetPhi();
	eleSCEt_1    = electron1->GetEnergySC()*sin(2.*atan(exp(-electron1->GetEtaSC())));
	eleSCEta_1   = electron1->GetEtaSC();
	eleSCPhi_1   = electron1->GetPhiSC();
	eleEoverP_1  = electron1->GetEoverP();
//	cout << "charge:" << eleCharge_1 << ":" << "Pt:" << elePt_1 << endl;
	eleCharge_2  = electron2->GetCharge();
	elePt_2      = electron2->GetPt();
	eleEta_2     = electron2->GetEta();
	elePhi_2     = electron2->GetPhi();
	eleSCEt_2    = electron2->GetEnergySC()*sin(2.*atan(exp(-electron2->GetEtaSC())));
	eleSCEta_2   = electron2->GetEtaSC();
	eleSCPhi_2   = electron2->GetPhiSC();
	eleEoverP_2  = electron2->GetEoverP(); 
	
	TLorentzVector ele1, ele2, diele, ele;
	ele1.SetPtEtaPhiM(elePt_1,eleEta_1,elePhi_1,eleMass);
	ele2.SetPtEtaPhiM(elePt_2,eleEta_2,elePhi_2,eleMass);
	
	diele = ele1 + ele2;
	vSum_ee_M = diele.M();
	vSum_ee_Energy = diele.Energy();
	vSum_ee_Pt  = diele.Pt();
	vSum_ee_Eta = diele.Eta();
	vSum_ee_Phi = diele.Phi();
	vSum_ee_Rapidity = diele.Rapidity();
	
	ele_dpt  = fabs(elePt_1  - elePt_2);
	ele_deta = fabs(eleEta_1 - eleEta_2);
	ele_dphi = getDPHI(elePhi_1,elePhi_2);
	ele_acop = 1- (ele_dphi/3.141592653589);  

   cos_photon_pair_helicity0 = cosphotonpair(ele1, diele, 0); // Boost of one photon in the pair direction (in the rest frame of the pair). The other will be at pi rads from the 1st.
   cos_photon_pair_helicity1 = cosphotonpair(ele1, diele, 1); 
   costhetastar = costhetastar_CS(ele1,diele);
 
   if(eleCharge_1 == -1){
    ele.SetPtEtaPhiM(elePt_1,eleEta_1,elePhi_1,eleMass);
   // cout << "neg charge:" << eleCharge_1 << ":" << "neg Pt:" << elePt_1 << endl;
   }
   else if (eleCharge_1 == 1){
    ele.SetPtEtaPhiM(elePt_2,eleEta_2,elePhi_2,eleMass);
   // cout << "neg charge:" << eleCharge_2 << ":" << "neg Pt:" << elePt_2 << endl;
   }
    //creating the 4vectors for Reconstrutcted MC

                float RecoDPx = diele.Pt()*cos(diele.Phi());
                float RecoDPy = diele.Pt()*sin(diele.Phi());
                float RecoMPx = ele.Pt()*cos(ele.Phi());
                float RecoMPy = ele.Pt()*sin(ele.Phi());

///////////////////////////////////////
		ele_PhiModulation = GetAngle( RecoDPx, RecoDPy, RecoMPx, RecoMPy ); // Reco deltaphi between -pi and pi
   // cout << "Reco phi modulation:" << ele_PhiModulation << endl;

	} //electron charge
   } //if two good electrons
     //samplename
    
    tr->Fill(); 
  } //nevents
  Log(0) << "Number of events triggered:" << trigger_passed << "\n" ;
  outFile->cd();
  hist->Write(); 
  hist_wozdc->Write(); 
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


int ireg(double et, double eta) {
   if (fabs(eta)<1.5) {
      if (et<3) return 1;
      else if (et<4) return 2;
      else if (et<5) return 3;
      else if (et<6) return 4;
      else if (et<7) return 5;
      else if (et<10) return 6;
      else if (et<14) return 7;
      else return 8;
   } else {
      if (et<3) return 9;
      else if (et<4) return 10;
      else if (et<5) return 11;
      else if (et<6) return 12;
      else if (et<7) return 13;
      else if (et<10) return 14;
      else if (et<14) return 15;
      else return 16;
   }
}

double SF_reco(double et, double eta) {
   if (fabs(eta)<1.5) {
      if (et<3) return 3.0043;
      else if (et<4) return 1.1693;
      else if (et<5) return 0.94329;
      else if (et<6) return 0.963789;
      else if (et<7) return 0.966706;
      else if (et<10) return 0.971721;
      else if (et<14) return 0.997961;
      else return 0.974345;
   } else {
      if (et<3) return 3.55878;
      else if (et<4) return 1.23809;
      else if (et<5) return 1.0068;
      else if (et<6) return 1.04464;
      else if (et<7) return 1.09334;
      else if (et<10) return 1.05342;
      else if (et<14) return 1.02881;
      else return 1.10496;
   }
}

double SF_uncert_reco(double et, double eta) {
   if (fabs(eta)<1.5) {
      if (et<3) return 0.438893;
      else if (et<4) return 0.0392326;
      else if (et<5) return 0.011089;
      else if (et<6) return 0.00769222;
      else if (et<7) return 0.00751307;
      else if (et<10) return 0.00577934;
      else if (et<14) return 0.00887761;
      else return 0.0183508;
   } else {
     if (et<3) return 1.46253;
      else if (et<4) return 0.0953307;
      else if (et<5) return 0.041571;
      else if (et<6) return 0.0359209;
      else if (et<7) return 0.0359005;
      else if (et<10) return 0.0295281;
      else if (et<14) return 0.0625413;
      else return 0.198791;
   }
}

double SF_trig(double et, double eta) {
   if (fabs(eta)<1.5) {
      if (et<3) return 1.20499;
      else if (et<4) return 1.08971;
      else if (et<5) return 1.04186;
      else if (et<6) return 1.02313;
      else if (et<7) return 1.02668;
      else if (et<10) return 1.02;
      else if (et<14) return 1.01338;
      else return 1.00665;
   } else {
      if (et<3) return 0.964039;
      else if (et<4) return 0.960585;
      else if (et<5) return 0.961264;
      else if (et<6) return 0.813962;
      else if (et<7) return 0.872495;
      else if (et<10) return 0.857176;
      else if (et<14) return 0.812907;
      else return 1.06977;
   }
}

double SF_uncert_trig(double et, double eta) {
   if (fabs(eta)<1.5) {
      if (et<3) return 0.0650598;
      else if (et<4) return 0.0145146;
      else if (et<5) return 0.0051025;
      else if (et<6) return 0.00298585;
      else if (et<7) return 0.00243375;
      else if (et<10) return 0.00162252;
      else if (et<14) return 0.00201474;
      else return 0.00607588;
   } else {
     if (et<3) return 0.151598;
      else if (et<4) return 0.0556211;
      else if (et<5) return 0.03333;
      else if (et<6) return 0.0318;
      else if (et<7) return 0.0305721;
      else if (et<10) return 0.0275402;
      else if (et<14) return 0.0718492;
      else return 0.12731;
   }
}




double cosphotonpair(TLorentzVector p1, TLorentzVector pair, bool helicityFrame )
{

  double cos = 0.;

  // Boost of one photon in the pair direction (i.e. to the rest frame of the pair). The other will be at pi rads from the 1st.
  //   //TVector3 vpair = pair.Vect();
  TVector3 vpair = pair.BoostVector();
  p1.Boost(-vpair);

  cos = TMath::Cos( p1.Angle( pair.Vect() ));

  if (helicityFrame) return cos;

  TVector3 zaxis = TVector3(0.,0.,1.);
  cos = TMath::Cos( p1.Angle( zaxis ));

  return cos;

}
float GetAngle( float vPX1, float vPY1, float vPX2, float vPY2)
{

   TVector2 DimuonVect;
   TVector2 MuonVect;

   DimuonVect.Set( vPX1, vPY1 );
   MuonVect.Set( vPX2, vPY2 );

   float DimuonMag = sqrt( DimuonVect.X()*DimuonVect.X() + DimuonVect.Y()*DimuonVect.Y()   );
   float MuonMag = sqrt( MuonVect.X()*MuonVect.X() + MuonVect.Y()*MuonVect.Y()   );


   float cosAngle = ( DimuonVect.X()*MuonVect.X() + DimuonVect.Y()*MuonVect.Y() )/( DimuonMag*MuonMag );
   float sinAngle = ( DimuonVect.X()*MuonVect.Y() - DimuonVect.Y()*MuonVect.X() )/( DimuonMag*MuonMag );
   float tanAngle = atan2(sinAngle, cosAngle);

   return tanAngle;

}


double costhetastar_CS(TLorentzVector p1, TLorentzVector pair)
{

double costhetastarCS = 2.*(pair.E()*p1.Pz()-pair.Pz()*p1.E())/(pair.M()*pair.Mt());

  return costhetastarCS;

}


