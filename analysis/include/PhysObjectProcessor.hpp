//  PhysObjectProcessor.hpp
//
//  Created by Jeremi Niedziela on 29/07/2019.

#ifndef PhysObjectProcessor_hpp
#define PhysObjectProcessor_hpp

#include "Helpers.hpp"
#include "PhysObject.hpp"

class PhysObjectProcessor;
extern PhysObjectProcessor physObjectProcessor;

/// This class performs operations on instances of PhysObject class
class PhysObjectProcessor {
public:
  /// Default constructor
  PhysObjectProcessor();
  
  /// Default destructor
  ~PhysObjectProcessor();
  
   ///Return TauTauPz
   double GetTauTauPz(const PhysObject &a, const PhysObject &b);
   ///Return TauTau
   ///TLorentzVector GetTauTau(const PhysObject &a, const PhysObject &b);
  ///Returns delaPhi
   double GetdeltaPhi(const PhysObject &a, const PhysObject &b);
  ///Returns DeltaEleR for ele and general trk
   double GetDeltaEleR(const PhysObject &a, const PhysObject &b, const PhysObject &c);
  ///Returns DeltaMuonR for muon and general trk
   double GetDeltaMuonR(const PhysObject &a, const PhysObject &b); 

  /// Returns ΔR=sqrt(Δη^2 + Δφ^2) between two PhysObjects
  double GetDeltaR(const PhysObject &a, const PhysObject &b);
  
  double GetDeltaR_SC(const PhysObject &a, const PhysObject &b);
  ///Returns muon lorentzvector
  TLorentzVector Getmu(const PhysObject &a);
  //Return electron Lorentzvector
  TLorentzVector Getele(const PhysObject &a);
  //REturn track Lorentzvector
  TLorentzVector Gettrk(const PhysObject &a);
  ///Returns a Lorentzvector of eletrk and mutrk
  TLorentzVector GetEletrk(const PhysObject &a, const PhysObject &b);
  TLorentzVector GetMutrk(const PhysObject &a, const PhysObject &b);

  /// Returns a lorentz vector which is a sum of provided objects, assuming they are photons
  TLorentzVector GetDiphoton(const PhysObject &a, const PhysObject &b);
  
  /// Returns a lorentz vector which is a sum of provided objects, assuming they are photons
  TLorentzVector GetTriphoton(const PhysObject &a, const PhysObject &b, const PhysObject &c);
  
  /// Returns a lorentz vector which is a sum of provided objects, assuming they are electrons
  TLorentzVector GetDielectron(const PhysObject &a, const PhysObject &b);
 
  /// Returns a lorentz vector which is a sum of provided objects, assuming they are muons
  TLorentzVector GetDimuon(const PhysObject &a, const PhysObject &b);
  
  ///Returns a lorentz vector which is a sum of provided objects,assuming they are electron+muon
  TLorentzVector GetEleMu(const PhysObject &a, const PhysObject &b);
  ///Return a Lorentzvector which is sum of electron+track
   TLorentzVector GetEleTrk(const PhysObject &a, const PhysObject &b);

  ///Return ScalarPt
  double GetScalarPt(const PhysObject &a, const PhysObject &b);
 
  /// Returns acoplanarity of pair of objects
  double GetAcoplanarity(const PhysObject &a, const PhysObject &b);
  
  /// Checks if physics object falls into HEM region of transition between barrel and endcap
  bool IsInCrack(const PhysObject &a);
  bool IsInHEM(const PhysObject &a);
  
  struct CompareByEnergy{
    bool operator() (const shared_ptr<PhysObject> &o1, const shared_ptr<PhysObject> &o2){
      return (o1->GetEnergy() > o2->GetEnergy());
    }
  };
  
private:
  
  
};

#endif /* PhysObjectProcessor_hpp */
