//  getDiphotonEfficiency
//
//  Created by Jeremi Niedziela on 22/07/2019.
//
// Calculates diphoton efficiency, defined as:
// Eff^{\gamma\gamma} = N^{reco}/N^{gen}
//
// where:
//
// N^{reco}: number of reconstructed diphoton events with E_T^{reco} > 2 GeV, |\eta^{reco}|<2.4,
// passing trigger, identification and exclusivity cuts
//
// N^{gen}: number of generated events with E_T^{gen} > 2 GeV, |\eta^{gen}|<2.4

#include "Helpers.hpp"
#include "EventProcessor.hpp"

// Eta and pt cuts
const double maxEta = 2.4;
const double minEt = 2; // GeV

// You can limit number of events analyzed here:
const int maxEvents = 10000;

int main()
{
  unique_ptr<EventProcessor> eventProcessor(new EventProcessor(kMClbl));
  
  int nGenEvents = 0;
  int nRecEvents = 0;
  int iEvent;
  
  for(iEvent=0; iEvent<eventProcessor->GetNevents(); iEvent++){
    if(iEvent >= maxEvents) break;
    
    Event event = eventProcessor->GetEvent(iEvent);
    
    // Check if gen event is within η and Et limits
    int nGenPhotonsPassing=0;
    
    for(int iGenPhoton=0; iGenPhoton<event.GetNgenPhotons(); iGenPhoton++){
      if((fabs(event.GetGenPhotonEta(iGenPhoton)) < maxEta) &&
         event.GetGenPhotonEt(iGenPhoton) > minEt){
        nGenPhotonsPassing++;
      }
    }
    if(nGenPhotonsPassing == 2) nGenEvents++;
    
    // Check if rec event is within η and Et limits, has LbL trigger and passes other selections
    if(!event.HasLbLTrigger()) continue;
    
    int nRecPhotonSCsPassing = 0;
    for(int iRecPhotonSC=0; iRecPhotonSC<event.GetNphotonSCs(); iRecPhotonSC++){
      if(fabs(event.GetPhotonSCeta(iRecPhotonSC)) < maxEta &&
         event.GetPhotonSCet(iRecPhotonSC) > minEt){
        nRecPhotonSCsPassing++;
      }
    }
    
    if(nRecPhotonSCsPassing != 2) continue;
    
    // Here check also other selections!
    // ...
    
    nRecEvents++;
  }
  
  cout<<"N event analyzed: "<<iEvent<<endl;
  cout<<"N gen events within limits: "<<nGenEvents<<endl;
  cout<<"N rec events passing selection: "<<nRecEvents<<endl;
  
  cout<<"Diphoton efficiency : "<<(double)nRecEvents/nGenEvents<<endl;
  
  return 0;
}
