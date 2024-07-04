//  PhysObjectProcessor.cpp
//
//  Created by Jeremi Niedziela on 29/07/2019.

#include "PhysObjectProcessor.hpp"
#include "ConfigManager.hpp"
#include "Logger.hpp"

PhysObjectProcessor physObjectProcessor = PhysObjectProcessor();

PhysObjectProcessor::PhysObjectProcessor()
{
  
}

PhysObjectProcessor::~PhysObjectProcessor()
{
  
}

double PhysObjectProcessor::GetdeltaPhi(const PhysObject &a, const PhysObject &b)
{
  double phi1 = a.GetPhi();
  double phi2 = b.GetPhi();

    // Make sure that angles are in range [0, 2π)
  while(phi1 < 0)               phi1 += 2*TMath::Pi();
  while(phi1 >= 2*TMath::Pi())  phi1 -= 2*TMath::Pi();
  
  while(phi2 < 0)               phi2 += 2*TMath::Pi();
  while(phi2 >= 2*TMath::Pi())  phi2 -= 2*TMath::Pi();
  double deltaPhi = fabs(phi2-phi1);
  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
  return deltaPhi;

}
//Date:20/08/2022, PxPyPzM Lorentz vector
//TLorentzVector PhysObjectProcessor::GetTauTau(const PhysObject &a, const PhysObject &b)
//{
  //TLorentzVector aVec, bVec;
   //double muMass = 105.6583755e-3;
   //double eleMass = 0.5109989461e-3;

  //aVec.SetPxPyPzE(a.GetPx(), a.GetPy(), a.GetPz(), a.GetEnergy());
  //bVec.SetPxPyPzE(b.GetPx(), b.GetPy(), b.GetPz(), b.GetEnergy());

  //return aVec + bVec;
//}
//
//
//Date:20/08/2022,px,py,pz
double PhysObjectProcessor::GetTauTauPz(const PhysObject &a, const PhysObject &b)
{
  double pz1 = a.GetPz();
  double pz2 = b.GetPz();

  double TauTauPz = (pz1+pz2);
  return TauTauPz;
    
}
   



//
//Date:18/08/2022,For deltaR between eleTrk and muonTrk with generalTrk
//for electron
double PhysObjectProcessor::GetDeltaEleR(const PhysObject &a, const PhysObject &b, const PhysObject &c)
{
  return (sqrt(pow(a.eta-b.eta, 2) + pow(a.phi-b.phi, 2)) or sqrt(pow(a.eta-c.eta, 2) + pow(a.phi-c.phi, 2)));
}

//for muon
double PhysObjectProcessor::GetDeltaMuonR(const PhysObject &a, const PhysObject &b)
{
  return sqrt(pow(a.eta-b.innerEta, 2) + pow(a.phi-b.innerPhi, 2));
}

//
double PhysObjectProcessor::GetDeltaR(const PhysObject &a, const PhysObject &b)
{
  return sqrt(pow(a.eta-b.eta, 2) + pow(a.phi-b.phi, 2));
}

double PhysObjectProcessor::GetDeltaR_SC(const PhysObject &a, const PhysObject &b)
{
  return sqrt(pow(a.etaSC-b.eta, 2) + pow(a.phiSC-b.phi, 2));
}

//Date:22/08/2022, for muon lorentzvector
TLorentzVector PhysObjectProcessor::Getmu(const PhysObject &a)
{
  TLorentzVector aVec;

  double  muMass = 105.6583755e-3;
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, muMass);

  return aVec;
}
//
//Date:22/08/2022, for electron lorentzvectoir
TLorentzVector PhysObjectProcessor::Getele(const PhysObject &a)
{
  TLorentzVector aVec;

  double eleMass = 0.5109989461e-3 ;
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, eleMass);

  return aVec;
}
//
//Date:12/09/2022
TLorentzVector PhysObjectProcessor::Gettrk(const PhysObject &a)
{
  TLorentzVector aVec;

  double trkMass = 0.13957018 ;
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, trkMass);

  return aVec;
}

//Date:20/09/2022, lorentz vector for e and mu track
TLorentzVector PhysObjectProcessor::GetEletrk(const PhysObject &a, const PhysObject &b)
{
  TLorentzVector aVec, bVec;

  double trkMass = 0.5109989461e-3 ;
  double eleMass = 0.5109989461e-3;
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, eleMass);
  bVec.SetPtEtaPhiM(b.pt, b.eta, b.phi, trkMass);

  return aVec + bVec;
 
}

TLorentzVector PhysObjectProcessor::GetMutrk(const PhysObject &a, const PhysObject &b)
{
  TLorentzVector aVec, bVec;

  double trkMass = 105.6583755e-3 ;
  double muMass = 105.6583755e-3;
 
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, muMass);
  bVec.SetPtEtaPhiM(b.pt, b.eta, b.phi, trkMass);

  return aVec + bVec;
}


TLorentzVector PhysObjectProcessor::GetDiphoton(const PhysObject &a, const PhysObject &b)
{
  TLorentzVector aVec, bVec;
  
  aVec.SetPtEtaPhiE(a.GetEt(), a.GetEta(), a.GetPhi(), a.GetEnergy());
  bVec.SetPtEtaPhiE(b.GetEt(), b.GetEta(), b.GetPhi(), b.GetEnergy());

  return aVec + bVec;
}

TLorentzVector PhysObjectProcessor::GetTriphoton(const PhysObject &a, const PhysObject &b, const PhysObject &c)
{
  TLorentzVector aVec, bVec, cVec;
  
  aVec.SetPtEtaPhiE(a.GetEt(), a.GetEta(), a.GetPhi(), a.GetEnergy());
  bVec.SetPtEtaPhiE(b.GetEt(), b.GetEta(), b.GetPhi(), b.GetEnergy());
  cVec.SetPtEtaPhiE(c.GetEt(), c.GetEta(), c.GetPhi(), c.GetEnergy());

  return aVec + bVec + cVec;
}

TLorentzVector PhysObjectProcessor::GetDielectron(const PhysObject &a, const PhysObject &b)
{
  TLorentzVector aVec, bVec;
  
  double eleMass = 0.5109989461e-3;
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, eleMass);
  bVec.SetPtEtaPhiM(b.pt, b.eta, b.phi, eleMass);

  return aVec + bVec;
}

TLorentzVector PhysObjectProcessor::GetDimuon(const PhysObject &a, const PhysObject &b)
{
  TLorentzVector aVec, bVec;
  
  double muMass = 105.6583755e-3;
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, muMass);
  bVec.SetPtEtaPhiM(b.pt, b.eta, b.phi, muMass);

  return aVec + bVec;
}
//Date-16/06/2022
TLorentzVector PhysObjectProcessor::GetEleMu(const PhysObject &a, const PhysObject &b)
{
  TLorentzVector aVec, bVec;

  double muMass = 105.6583755e-3;
  double eleMass = 0.5109989461e-3;
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, muMass);
  bVec.SetPtEtaPhiM(b.pt, b.eta, b.phi, eleMass);

  return aVec + bVec;
}
//Date-12/09/2022
TLorentzVector PhysObjectProcessor::GetEleTrk(const PhysObject &a, const PhysObject &b)
{
  TLorentzVector aVec, bVec;

  double trkMass = 0.13957018;
  double eleMass = 0.5109989461e-3;
  aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, trkMass);
  bVec.SetPtEtaPhiM(b.pt, b.eta, b.phi, eleMass);

  return aVec + bVec;
}

//Date-5/08/202
double  PhysObjectProcessor::GetScalarPt(const PhysObject &a, const PhysObject &b)
{
//  TLorentzVector aVec, bVec;

  //double muMass = 105.6583755e-3;
 // double eleMass = 0.5109989461e-3;
 // aVec.SetPtEtaPhiM(a.pt, a.eta, a.phi, muMass);
  //bVec.SetPtEtaPhiM(b.pt, b.eta, b.phi, eleMass);
 // double ScalarSum = aVec.GetPt() + bVec.GetPt();
  return a.pt + b.pt;
}



double PhysObjectProcessor::GetAcoplanarity(const PhysObject &a, const PhysObject &b)
{
  double phi1 = a.GetPhi();
  double phi2 = b.GetPhi();
  
  // Make sure that angles are in range [0, 2π)
  while(phi1 < 0)               phi1 += 2*TMath::Pi();
  while(phi1 >= 2*TMath::Pi())  phi1 -= 2*TMath::Pi();
  
  while(phi2 < 0)               phi2 += 2*TMath::Pi();
  while(phi2 >= 2*TMath::Pi())  phi2 -= 2*TMath::Pi();
  
  double deltaPhi = fabs(phi2-phi1);
  
  // Make sure that Δφ is in range [0, π]
  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
  
  // Calcualte acoplanarity
  double aco = 1 - deltaPhi/TMath::Pi();
  
  if(aco < 0){
    Log(0)<<"ERROR -- acoplanarity <0, which should never happen!!\n";
  }
  
  return aco;
}

bool PhysObjectProcessor::IsInCrack(const PhysObject &a)
{
  if(fabs(a.GetEtaSC()) > config.params("ecalCrackMin") &&
     fabs(a.GetEtaSC()) < config.params("ecalCrackMax")) return true;
  
  return false;
}

bool PhysObjectProcessor::IsInHEM(const PhysObject &a)
{
  if(a.GetEtaSC() < config.params("ecalHEMmaxEta") &&
     a.GetEtaSC() > -maxEtaEE &&
     a.GetPhiSC() > config.params("ecalHEMminPhi") &&
     a.GetPhiSC() < config.params("ecalHEMmaxPhi")) return true;
  
  return false;
}
