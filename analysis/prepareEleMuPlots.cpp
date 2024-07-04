//Based on prepareBasicPlots.cpp
#include <iostream>

#include "Helpers.hpp"
#include "EventProcessor.hpp"
#include "PhysObjectProcessor.hpp"
#include "PhysObject.hpp"
#include "ConfigManager.hpp"
#include "EventDisplay.hpp"
#include "Logger.hpp"
#include "tnp_muon_UPC_PbPb2018.h"

  string configPath = "configs/efficiencies_eleNoIsolation_newThresholdsEta2p2.md";
  string outputPath = "Output/basicPlots_test.root";

bool saveCalosFailingNEE = true;
//bool saveTriphotonHists = true;
//bool checkTriggers = false;

//int nThreePhotonEvents = 0;


// Only those datasets will be analyzed
 const vector<EDataset> datasetsToAnalyze = {
   kData,kMCqedSC,kMCgammaUPC,kMCmumuFSR,kMCTauTau,
};

double luminosity = 1642.79;
//SC
double crossSectionSC = 570;//MC_Arash
double generatedEvtSC = 1000000;//SC_Arash

//FSR
double crossSectionFSR = 139.7;//MC_FSR_Arash
double generatedEvtFSR = 2471810;//FSR

//GammaUPC
double crossSectionGammaUPC = 1060.8959;//GammaUPC
//double generatedEvtGammaUPC = 5290000;//GammaUPC
double generatedEvtGammaUPC = 5268607;//GammaUPC//New
//Trigger and reco+ID Scale Factor
double ScaleFactor = 0.952* //electron reco+ID 05 Oct 2022 from singleEG3 only
                      // 1.008* //L1 EG trigger 05 Oct 2022 from singleEG5 only
                      // 0.866* // HF veto
                       0.85;  //Neutral Ex//Muon left                      
//Weight calculation
//double Weight = 1;//Data
//double W  = (luminosity * crossSectionSC/generatedEvtSC;//SC
//double W = (luminosity * crossSectionGammaUPC/generatedEvtGammaUPC;//GammaUPC
//double W = (luminosity * crossSectionFSR)/generatedEvtFSR;//FSR
double W = 1;
double ET;
double PT;
int idx = 0;
//double W = 1;
//double Other = (luminosity * crossSectionSC * ScaleFactor)/generatedEvtSC;
//double Other = (luminosity * crossSectionGammaUPC * ScaleFactor)/generatedEvtGammaUPC;
double Other  = (luminosity * crossSectionFSR * ScaleFactor)/generatedEvtFSR;
//double Weight;
//double Weight = tnp_weight_softid_upc_pbpb( PT, ET, idx);
vector<string> suffixes = {
  "all"
};

vector<tuple<string, int, double, double>> histParams = {
 
   // title                         nBins                min              max
  {"dilepton_acoplanarity"       ,    40          ,     0       ,       0.4  },
  {"muon_pt"                     ,    6          ,     2.5     ,      14.5  },
  {"muon_eta"                    ,    15           ,    -2.5    ,       2.5  },
  {"muon_phi"                    ,    21           ,   -3.5       ,       3.5 },
  {"muon_px"                     ,    28          ,   -14.0       ,      14.0  },
  {"muon_py"                     ,    35          ,   -15.0       ,      20.0  },
  {"muon_pz"                     ,    20          ,   -20.0       ,      20.0  },
  {"e_pt"                        ,    10          ,     2.0       ,      12.0  },
  {"e_eta"                       ,    15          ,    -2.5       ,       2.5  },
  {"e_phi"                       ,    21         ,    -3.5       ,       3.5  },
  {"e_px"                        ,    30          ,   -15.0       ,      15.0  },
  {"e_py"                        ,    30          ,   -15.0       ,      15.0  },
  {"e_pz"                        ,    28          ,   -26.0       ,      30.0  },
  {"dilepton_mass"               ,    5          ,       4       ,      24.0  },
  //{"dilepton_mass"               ,    36          ,       4       ,      22.0  },
  {"dilepton_rapidity"           ,    15          ,    -2.5       ,       2.5  },
  {"dilepton_pt"                 ,    5          ,       0       ,      10.0  },
  //{"dilepton_pt"                 ,    60          ,       0       ,      6.0  },
  {"dilepton_phi"                ,    30          ,     2.6       ,       3.2  },
  {"dilepton_Scalar_Pt"          ,    20          ,       0       ,      20.0  },
  {"dilepton_deltaEleR"          ,    20          ,       0       ,      10.0  },
  {"dilepton_deltaR"             ,    40          ,       0       ,       5.0  },
  {"dilepton_deltaRelectron"     ,    40          ,       0       ,       0.2  },
  {"dilepton_deltaRmuon"         ,    40          ,       0       ,       0.2  },
  {"dilepton_TauTauPz"           ,    30          ,   -30.0       ,      30.0  },
  {"dilepton_TauTauPx"           ,   120          ,   -30.0       ,      30.0  },
  {"dilepton_TauTauPy"           ,   120          ,   -30.0       ,      30.0  },
  {"dilepton_TauTauEta"          ,   10           ,     0.0       ,       5.0  },
  {"dilepton_Pseudorapidity"     ,   10           ,     0.0       ,       5.0  },
 // {"muon_InnerDz"                ,  2200          , -1100.0       ,    1100.0  },
 // {"muon_Dz"                ,  2200          , -1100.0       ,    1100.0  },
  {"track_pt"                    ,   40           ,     0.0       ,       20.0 },
  {"track_eta"                   ,   10           ,    -2.5       ,        2.5 },
  {"track_phi"                   ,   12           ,    -3.0       ,        3.0 },
  {"NTrk"                        ,   25           ,     0.0       ,       25.0 },
};
/*
vector<tuple<string, int, double, double, int, double, double>> histParams2D = {
  // title         nBinsX minX   maxX   nBinsY minY  maxY
     {"MuonPtEta" , 6 , 2.5 , 16.5 , 2 , -2.4 , 2.4 },
  //   {"HoverEmapDen" , 314 , -3.14 , 3.14 , 460 , -2.3 , 2.3 },
  };
*/
//bool endsWith(const std::string &mainStr, const std::string &toMatch)
//{
 // if(mainStr.size() >= toMatch.size() &&
   //  mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(), toMatch) == 0)
   // return true;
  //else
   // return false;
//}

///Muon
void fillMuonHists(Event &event, const map<string, TH1D*> &hists, string datasetName, string suffix="")
{
  for(auto muon : event.GetPhysObjects(EPhysObjType::kGoodMuon)){
     ET = muon->GetEta();
     PT = muon->GetPt();
    double Weight1 = tnp_weight_softid_upc_pbpb( PT, ET, idx);
    double Weight2 = tnp_weight_trigger_upc_pbpb( PT, ET, idx);
//    double W = Other * Weight1 * Weight2;
     
    hists.at("muon_pt_" +suffix+datasetName)->Fill(muon->GetPt(),W);
    hists.at("muon_eta_"+suffix+datasetName)->Fill(muon->GetEta(),W);
    hists.at("muon_phi_"+suffix+datasetName)->Fill(muon->GetPhi(),W);
     cout << "MuonPt:" << muon->GetPt() << endl;    
     cout << "MuonEta:" << muon->GetEta() << endl;    
     cout << "Weight1:" << W << endl; 
    //hists2D.at("MuonPtEta_"+suffix+datasetName)->Fill(muon->GetPt(), muon->GetEta(),W);  

  }
}

//Electron
void fillElectronHists(Event &event, const map<string, TH1D*> &hists, string datasetName, string suffix="")
{
  for(auto electron : event.GetPhysObjects(EPhysObjType::kGoodElectron)){
    double Weight1 = tnp_weight_softid_upc_pbpb( PT, ET, idx);
    double Weight2 = tnp_weight_trigger_upc_pbpb( PT, ET, idx);
 //   double W = Other * Weight1 * Weight2;
    hists.at("e_pt_" +suffix+datasetName)->Fill(electron->GetPt(),W);
    cout<< "ElectronPt:"<< electron->GetPt() <<endl;
    hists.at("e_eta_"+suffix+datasetName)->Fill(electron->GetEta(),W);
    hists.at("e_phi_"+suffix+datasetName)->Fill(electron->GetPhi(),W);
    cout << "Weight2:" << W << endl;
 }
}



//Electron+Muon
void fillEleMuHists(Event &event, const map<string, TH1D*> &hists, string datasetName, string suffix="")
{
   
  TLorentzVector elemu = physObjectProcessor.GetEleMu(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0],
                                                      *event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]);

  TLorentzVector muonpz = physObjectProcessor.Getmu(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0]);

  TLorentzVector elepz = physObjectProcessor.Getele(*event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]);

  double ScalarPz = (physObjectProcessor.Getmu(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0]).Pz() + physObjectProcessor.Getele(*event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]).Pz());
  double TauTauEta = (physObjectProcessor.Getmu(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0]).Eta() - physObjectProcessor.Getele(*event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]).Eta());  

  double Weight1 = tnp_weight_softid_upc_pbpb( PT, ET, idx);
  double Weight2 = tnp_weight_trigger_upc_pbpb( PT, ET, idx);
 // double W = Other * Weight1 * Weight2;
   
  hists.at("dilepton_mass_"     +suffix+datasetName)->Fill(elemu.M(),W);
  hists.at("dilepton_rapidity_" +suffix+datasetName)->Fill(elemu.Rapidity(),W);
  hists.at("dilepton_TauTauEta_" +suffix+datasetName)->Fill(fabs(TauTauEta),W);
  cout << "TauTauEta:" <<fabs(TauTauEta) << endl;
  hists.at("dilepton_pt_"       +suffix+datasetName)->Fill(elemu.Pt(),W);
  hists.at("dilepton_TauTauPz_" +suffix+datasetName)->Fill(elemu.Pz(),W);
  cout << "TauTauPz:" << elemu.Pz() << endl;
  hists.at("dilepton_TauTauPx_" +suffix+datasetName)->Fill(elemu.Px(),W);
  hists.at("dilepton_TauTauPy_" +suffix+datasetName)->Fill(elemu.Py(),W);
  hists.at("e_pz_" +suffix+datasetName)->Fill(elepz.Pz(),W);
  hists.at("e_px_" +suffix+datasetName)->Fill(elepz.Px(),W);
  hists.at("e_py_" +suffix+datasetName)->Fill(elepz.Py(),W);
  hists.at("muon_pz_" +suffix+datasetName)->Fill(muonpz.Pz(),W);
  hists.at("muon_px_" +suffix+datasetName)->Fill(muonpz.Px(),W);
  hists.at("muon_py_" +suffix+datasetName)->Fill(muonpz.Py(),W);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  Delta R betweeen electron and muon
  //for(auto track : event.GetPhysObjects(EPhysObjType::kGoodGeneralTrack)){
    //  for(auto electron : event.GetPhysObjects(EPhysObjType::kGoodElectron)){
      //   if(track->GetCharge() == electron->GetCharge()){
        //     double deltaRelectron = sqrt(pow(electron->GetEta() - track->GetEta(), 2) + pow(electron->GetPhi() - track->GetPhi(), 2));

         // hists.at("dilepton_deltaRelectron_"+suffix+datasetName)->Fill(deltaRelectron,Weight);
         // cout << "trackcharge1:" << track->GetCharge() << endl;
         // cout << "trackEta1:" << track->GetEta() << endl;
         // cout << "trackPhi1:" << track->GetPhi() << endl;
         // cout << "deltaRelectron:" << deltaRelectron << endl;
        // cout << "ElectronCharge:" << electron->GetCharge() << endl;
       //  cout << "ElectronEta1:" << electron->GetEta() << endl;
        // cout << "ElectronPhi1:" << electron->GetPhi() << endl;

    // }
  //  }
 //  }

   // for(auto track : event.GetPhysObjects(EPhysObjType::kGoodGeneralTrack)){
     // for(auto muon : event.GetPhysObjects(EPhysObjType::kGoodMuon)){
       //  if(track->GetCharge() == muon->GetCharge())
      //  {
         //double deltaRmuon = sqrt(pow(muon->GetEta() - track->GetEta(), 2) + pow(muon->GetPhi() - track->GetPhi(), 2));


        // hists.at("dilepton_deltaRmuon_"+suffix+datasetName)->Fill(deltaRmuon,Weight);
        // cout << "trackcharge2:" << track->GetCharge() << endl;
       //  cout << "trackPhi2:" << track->GetPhi() << endl;
        // cout << "trackEta2:" << track->GetEta() << endl;
       //  cout << "deltaRmuon:" << deltaRmuon << endl;
       //  cout << "MuonCharge:" << muon->GetCharge() << endl;
       //  cout << "MuonPhi2:" << muon->GetPhi() << endl;
       //  cout << "MuonEta2:" << muon->GetEta() << endl;
    // }

  // }
// }
//
for(auto track : event.GetPhysObjects(EPhysObjType::kGoodGeneralTrack)){
      for(auto electron : event.GetPhysObjects(EPhysObjType::kGoodElectron)){
         if(track->GetCharge() == electron->GetCharge()){
            double deltaRelectron = sqrt(pow(electron->GetEta() - track->GetEta(), 2) + pow(electron->GetPhi() - track->GetPhi(), 2));
            double Weight1 = tnp_weight_softid_upc_pbpb( PT, ET, idx);       
            double Weight2 = tnp_weight_trigger_upc_pbpb( PT, ET, idx);       
   //         double W = Other * Weight1 * Weight2;
            hists.at("dilepton_deltaRelectron_"+suffix+datasetName)->Fill(deltaRelectron,W);
            cout << "trackcharge1:" << track->GetCharge() << endl;
            cout << "trackEta1:" << track->GetEta() << endl;
            cout << "trackPhi1:" << track->GetPhi() << endl;
            cout << "deltaRelectron:" << deltaRelectron << endl;
            cout << "ElectronCharge:" << electron->GetCharge() << endl;
            cout << "ElectronEta1:" << electron->GetEta() << endl;
            cout << "ElectronPhi1:" << electron->GetPhi() << endl;
           
            }
           }
        }
           


for(auto track : event.GetPhysObjects(EPhysObjType::kGoodGeneralTrack)){
      for(auto muon : event.GetPhysObjects(EPhysObjType::kGoodMuon)){
         if(track->GetCharge() == muon->GetCharge())
        {
double deltaRmuon = sqrt(pow(muon->GetEta() - track->GetEta(), 2) + pow(muon->GetPhi() - track->GetPhi(), 2));

         double Weight1 = tnp_weight_softid_upc_pbpb( PT, ET, idx);
         double Weight2 = tnp_weight_trigger_upc_pbpb( PT, ET, idx);
  //       double W = Other * Weight1 * Weight2;
         hists.at("dilepton_deltaRmuon_"+suffix+datasetName)->Fill(deltaRmuon,W);
         cout << "trackcharge2:" << track->GetCharge() << endl;
         cout << "trackPhi2:" << track->GetPhi() << endl;
         cout << "trackEta2:" << track->GetEta() << endl;
         cout << "deltaRmuon:" << deltaRmuon << endl;
         cout << "MuonCharge:" << muon->GetCharge() << endl;
         cout << "MuonPhi2:" << muon->GetPhi() << endl;
         cout << "MuonEta2:" << muon->GetEta() << endl;
     }

   }
 }
         cout << "Weight3:" << W << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
}

/////////////////Track Hists//////////////////////////////////////////////////////////////

void fillTrackHists(Event &event, const map<string, TH1D*> &hists, EDataset dataset, string suffix="")
{
       string name = datasetName.at(dataset);
        if(suffix != "") suffix += "_";
  //     int nTracks = (int)event.GetPhysObjects(EPhysObjType::kGeneralTrack).size();
//       hists.at("NTrk_"+suffix+name)->Fill(nTracks);
            
          

 // }   
}
////////////////////////////////////////////////////////////////////////////////////////


void fillTauTauHistograms(Event &event, const map<string, TH1D*> &hists, EDataset dataset, vector<string> suffix_list)
{
  string name = datasetName.at(dataset);
  int cutThrough=0;

  double aco = physObjectProcessor.GetAcoplanarity(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0],
                                                   *event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]);


  double delPhi =  physObjectProcessor.GetdeltaPhi(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0],
                                                   *event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]);

  double scalarsumpt =  physObjectProcessor.GetScalarPt(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0],
                                                        *event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]);     
  

  double deltaR = physObjectProcessor.GetDeltaR(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0],
                                                *event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]);
 


  double tautauPz = physObjectProcessor.GetTauTauPz(*event.GetPhysObjects(EPhysObjType::kGoodMuon)[0],
                                                    *event.GetPhysObjects(EPhysObjType::kGoodElectron)[0]);

  int nTracks = (int)event.GetPhysObjects(EPhysObjType::kGoodGeneralTrack).size();
     //  cout << "nTracks = " << nTracks << endl;



  for(string suffix : suffix_list){
   
    if(aco < 0.01 || aco > 0.4) continue;

    if(suffix != "") suffix += "_";

    fillMuonHists(event, hists, name, suffix);
    fillElectronHists(event, hists, name, suffix);
    fillEleMuHists(event, hists, name, suffix);
    double Weight1 = tnp_weight_softid_upc_pbpb( PT, ET, idx);
    double Weight2 = tnp_weight_trigger_upc_pbpb( PT, ET, idx);
 //   double W = Other * Weight1 * Weight2;
     cout << "Weight4:" << W << endl;     
    hists.at("dilepton_acoplanarity_"+suffix+name)->Fill(aco,W);
    hists.at("dilepton_phi_"+suffix+name)->Fill(delPhi,W);
    hists.at("dilepton_Scalar_Pt_"+suffix+name)->Fill(scalarsumpt,W);
    hists.at("dilepton_deltaR_"+suffix+name)->Fill(deltaR,W);
    hists.at("NTrk_"+suffix+name)->Fill(nTracks);
   // hists.at("dilepton_deltaRelectron_"+suffix+name)->Fill(deltaRelectron,W);
   // fillMuonHists(event, hists, name, suffix);
    //fillElectronHists(event, hists, name, suffix);
   // fillEleMuHists(event, hists, name, suffix);
 //   fillTrackHists(event, hists, dataset, suffix);
    //cout << "Weight4:" << Weight << endl;  
  }
}


/// Creates histograms, cut through and event counters for given dataset name, for each
///// histogram specified in `histParams` vector.


void InitializeHistograms(map<string, TH1D*> &hists, string datasetType, string suffix="")
{
  if(suffix != "") suffix = "_" + suffix;

  for(auto &[histName, nBins, min, max] : histParams){
    string title = histName + suffix + "_" + datasetType;

    if(hists.find(title) != hists.end()) continue;
    hists[title] = new TH1D(title.c_str(), title.c_str(), nBins, min, max);
  }
}


int main(int argc, char* argv[])
{
  if(argc != 5){
    cout<<"This app requires 4 parameters.\n";
    cout<<"./prepareBasicPlots configPath inputPath outputPath datasetName[Data|QED_SC|QED_SL|LbL|CEP|GAMMA_UPC|MUMU_FSR|TAUTAU]\n";
    exit(0);
  }
  cout<<"Reading input arguments"<<endl;

  string inputPath = "";
  string sampleName = "";

  //if(argc == 5){
  configPath = argv[1];
  inputPath  = argv[2];
  outputPath = argv[3];
  sampleName = argv[4];
  //}
  cout<<"Config: "<<configPath<<endl;
  cout<<"Input: "<<inputPath<<endl;
  cout<<"Output: "<<outputPath<<endl;
  cout<<"Sample name: "<<sampleName<<endl;


  config = ConfigManager(configPath);

  map<string, TH1D*> hists;
  TFile *outFile = new TFile(outputPath.c_str(), "recreate");


 //if (endsWith(inputPath, "root")){
   // cout << "root file" << endl;
    EDataset dataset = nDatasets;

    if(sampleName == "Data")    dataset = kData;
    if(sampleName == "QED_SC")  dataset = kMCqedSC;
    if(sampleName == "QED_SL")  dataset = kMCqedSL;
    if(sampleName == "LbL")     dataset = kMClbl;
    if(sampleName == "CEP")     dataset = kMCcep;
    if(sampleName == "GAMMA_UPC")     dataset = kMCgammaUPC;
    if(sampleName == "MUMU_FSR")     dataset = kMCmumuFSR;
    if(sampleName == "TAUTAU")     dataset = kMCTauTau;
    auto events = make_unique<EventProcessor>(inputPath, dataset);

    for(string suffix : suffixes){
      InitializeHistograms(hists, sampleName, suffix);
    }


    if(dataset == nDatasets){
      Log(0)<<"ERROR -- unknown dataset name provided: "<<sampleName<<"\n";
      exit(0);
    }
    
    for(int iEvent=0; iEvent<events->GetNevents(); iEvent++){
      if(iEvent%5 == 0)  Log(0)<<"Processing event "<<iEvent<<"\n";
      //if(iEvent%10000 == 0) Log(0)<<"Processing event "<<iEvent<<"\n";
      if(iEvent >= config.params("maxEvents")) break;
      auto event = events->GetEvent(iEvent);
      fillTauTauHistograms(*event, hists, dataset, suffixes);
      //fillDimuonHists(*event, hists, dataset);
    }

  cout << "Writing Histograms" <<endl;
  outFile->cd();
  for(auto &[histName, hist] : hists) hist->Write();


  outFile->Close();
  cout << "Finished" << endl;
}


   
