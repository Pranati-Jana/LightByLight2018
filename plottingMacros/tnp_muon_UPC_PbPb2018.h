#ifndef tnp_muon_UPC_PbPb2018_h
#define tnp_muon_UPC_PbPb2018_h

#include <iostream>

// IN THIS FILE YOU WILL FIND:
// ++++++++++++++
//
// - SoftID: (tnp_weight_softid_upc_pbpb) 
//   * idx = 0: nominal
//   * idx = -1: syst variation,  +1 sigma
//   * idx = -2: syst variation,  -1 sigma
//   * idx = +1: stat variation,  +1 sigma
//   * idx = +2: stat variation,  -1 sigma
//
// - Trigger: (tnp_weight_trigger_upc_pbpb) 
//   * idx = 0:  nominal
//   * idx = -1: syst variation,  +1 sigma
//   * idx = -2: syst variation,  -1 sigma
//   * idx = +1: stat variation,  +1 sigma
//   * idx = +2: stat variation,  -1 sigma

// THE INDIVIDUAL SFs
// ++++++++++++++++++
double tnp_weight_softid_upc_pbpb( const double& PTT, const double& ETT, const int& idx=0);
double tnp_weight_trigger_upc_pbpb(const double& PTT, const double& ETT, const int& idx=0);


///////////////////////////////////////////////////
//                 Soft I D    P b P b           //
///////////////////////////////////////////////////
double tnp_weight_softid_upc_pbpb( double& PTT, const double& ETT, const int& idx)
{
    double weight(-1), unc(0);
    const auto absETT = std::abs(ETT);
    
    if (absETT > 2.4) {
        std::cout << "[WARNING] Muon pseudo-rapidity (" << ETT << ") outside [-2.4, 2.4]" << std::endl;
        return 1;
    }
    
    // nominal
    if (absETT >= 0.0 && absETT < 1.0) {
        if (PTT >= 3.3 && PTT < 6.0)
            weight = 0.982402;
        else if (PTT >= 6.0 && PTT < 30.0)
            weight = 1.0087;
    }
    else if (absETT >= 1.0 && absETT < 1.6) {
        if (PTT >= 1.0 && PTT < 1.5)
            weight = 0.7915;
        else if (PTT >= 1.5 && PTT < 2.0)
            weight = 0.895003;
        else if (PTT >= 2.0 && PTT < 2.5)
            weight = 0.93851;
        else if (PTT >= 2.5 && PTT < 3.5)
            weight = 1.0079;
        else if (PTT >= 3.5 && PTT < 30.0)
            weight = 0.989807;
    }
    else if (absETT >= 1.6 && absETT < 2.4) {
        if (PTT >= 1.0 && PTT < 1.2)
            weight = 0.895651;
        else if (PTT >= 1.2 && PTT < 1.5)
            weight = 0.955272;
        else if (PTT >= 1.5 && PTT < 2.0)
            weight = 0.955527;
        else if (PTT >= 2.0 && PTT < 30.0)
            weight = 0.981941;
    }
    
    // statistics uncertainty
    if (idx == 1 || idx == 2) {
        if (absETT >= 0.0 && absETT < 1.0) {
            if (PTT >= 3.3 && PTT < 6.0)
                unc = 0.0258814;
            else if (PTT >= 6.0 && PTT < 30.0)
                unc = 0.0154683;
        }
        else if (absETT >= 1.0 && absETT < 1.6) {
            if (PTT >= 1.0 && PTT < 1.5)
                unc = 0.0299885;
            else if (PTT >= 1.5 && PTT < 2.0)
                unc = 0.000402623;
            else if (PTT >= 2.0 && PTT < 2.5)
                unc = 0.0436706;
            else if (PTT >= 2.5 && PTT < 3.5)
                unc = 0.0301486;
            else if (PTT >= 3.5 && PTT < 30.0)
                unc = 0.0266832;
        }
        else if (absETT >= 1.6 && absETT < 2.4) {
            if (PTT >= 1.0 && PTT < 1.2)
                unc = 0.0260094;
            else if (PTT >= 1.2 && PTT < 1.5)
                unc = 0.0134018;
            else if (PTT >= 1.5 && PTT < 2.0)
                unc = 0.00722138;
            else if (PTT >= 2.0 && PTT < 30.0)
                unc = 0.012674;
        }
    }
    // systematic uncertainty
    else if (idx == -1 || idx == -2) {
        if (absETT >= 0.0 && absETT < 1.0) {
            if (PTT >= 3.3 && PTT < 6.0)
                unc = 0.00339975;
            else if (PTT >= 6.0 && PTT < 30.0)
                unc = 0.000215163;
        }
        else if (absETT >= 1.0 && absETT < 1.6) {
            if (PTT >= 1.0 && PTT < 1.5)
                unc = 0.00854174;
            else if (PTT >= 1.5 && PTT < 2.0)
                unc = 0.0111685;
            else if (PTT >= 2.0 && PTT < 2.5)
                unc = 0.0111676;
            else if (PTT >= 2.5 && PTT < 3.5)
                unc = 0.0206014;
            else if (PTT >= 3.5 && PTT < 30.0)
                unc = 0.0279548;
        }
        else if (absETT >= 1.6 && absETT < 2.4) {
            if (PTT >= 1.0 && PTT < 1.2)
                unc = 0.0101644;
            else if (PTT >= 1.2 && PTT < 1.5)
                unc = 0.00834835;
            else if (PTT >= 1.5 && PTT < 2.0)
                unc = 0.00323968;
            else if (PTT >= 2.0 && PTT < 30.0)
                unc = 0.00697054;
        }
    }
    else if (idx != 0) {
        std::cout << "[WARNING] Invalid index " << idx << ". Valid values are: -2, -1, 0, 1, 2." << std::endl;
        return 1;
    }
    
    // up variation
    if (idx == 1 || idx == -1)
        weight += unc;
    // down variation
    else if (idx == 2 || idx == -2)
        weight -= unc;
    
    if (weight < 0) {
        std::cout << "[WARNING] Muon (PTT= " << PTT << ", ETT= " << ETT <<") outside of soft ID T&P kinematic range" << std::endl;
        return 1;
    }
    
    return weight;
}


///////////////////////////////////////////////////
//                 Trigger    P b P b            //
//        HLT_HIUPC_SingleMuOpen_NotMBHF2AND     //
///////////////////////////////////////////////////
double tnp_weight_trigger_upc_pbpb(const double& PTT, const double& ETT, const int& idx)
{
    double weight(-1), unc(0);
    const auto absETT = std::abs(ETT);
    
    if (absETT > 2.4) {
        std::cout << "[WARNING] Muon pseudo-rapidity (" << ETT << ") outside [-2.4, 2.4]" << std::endl;
        return 1;
    }
    
    // nominal
    if (absETT >= 0.0 && absETT < 1.1) {
        if (PTT >= 3.3 && PTT < 6.0)
            weight = 1.03745;
        else if (PTT >= 6.0 && PTT < 30.0)
            weight = 1.00452;
    }
    else if (absETT >= 1.1 && absETT < 1.6) {
        //if (PTT >= 2.15 && PTT < 2.7) 
        if (PTT >= 1.5 && PTT < 2.7) 
            weight = 1.06713;
        else if (PTT >= 2.7 && PTT < 3.3)
            weight = 1.0969;
        else if (PTT >= 3.3 && PTT < 4.5)
            weight = 0.953613;
        else if (PTT >= 4.5 && PTT < 30.0)
            weight = 0.993225;
    }
    else if (absETT >= 1.6 && absETT < 2.4) {
        if (PTT >= 1.2 && PTT < 1.5)
            weight = 1.59913;
        else if (PTT >= 1.5 && PTT < 2.0)
            weight = 1.0564;
        else if (PTT >= 2.0 && PTT < 2.5)
            weight = 1.05008;
        else if (PTT >= 2.5 && PTT < 30.0)
            weight = 1.01014;
    }
    
    // statistics uncertainty
    if (idx == 1 || idx == 2) {
        if (absETT >= 0.0 && absETT < 1.1) {
            if (PTT >= 3.3 && PTT < 6.0)
                unc = 0.0252039;
            else if (PTT >= 6.0 && PTT < 30.0)
                unc = 0.022945;
        }
        else if (absETT >= 1.1 && absETT < 1.6) {
            if (PTT >= 2.15 && PTT < 2.7)
                unc = 0.100835;
            else if (PTT >= 2.7 && PTT < 3.3)
                unc = 0.0599291;
            else if (PTT >= 3.3 && PTT < 4.5)
                unc = 0.0427363;
            else if (PTT >= 4.5 && PTT < 30.0)
                unc = 0.028474;
        }
        else if (absETT >= 1.6 && absETT < 2.4) {
            if (PTT >= 1.2 && PTT < 1.5)
                unc = 0.0271744;
            else if (PTT >= 1.5 && PTT < 2.0)
                unc = 0.0205656;
            else if (PTT >= 2.0 && PTT < 2.5)
                unc = 0.0206436;
            else if (PTT >= 2.5 && PTT < 30.0)
                unc = 0.00842663;
        }
    }
    // systematic uncertainty
    else if (idx == -1 || idx == -2) {
        if (absETT >= 0.0 && absETT < 1.1) {
            if (PTT >= 3.3 && PTT < 6.0)
                unc = 0.00302591;
            else if (PTT >= 6.0 && PTT < 30.0)
                unc = 0.0028882;
        }
        else if (absETT >= 1.1 && absETT < 1.6) {
            if (PTT >= 2.15 && PTT < 2.7)
                unc = 0.0401421;
            else if (PTT >= 2.7 && PTT < 3.3)
                unc = 0.0541902;
            else if (PTT >= 3.3 && PTT < 4.5)
                unc = 0.00723339;
            else if (PTT >= 4.5 && PTT < 30.0)
                unc = 0.00580133;
        }
        else if (absETT >= 1.6 && absETT < 2.4) {
            if (PTT >= 1.2 && PTT < 1.5)
                unc = 0.00739062;
            else if (PTT >= 1.5 && PTT < 2.0)
                unc = 0.0152213;
            else if (PTT >= 2.0 && PTT < 2.5)
                unc = 0.0165649;
            else if (PTT >= 2.5 && PTT < 30.0)
                unc = 0.00254733;
        }
    }
    else if (idx != 0) {
        std::cout << "[WARNING] Invalid index " << idx << ". Valid values are: -2, -1, 0, 1, 2." << std::endl;
        return 1;
    }
    
    // up variation
    if (idx == 1 || idx == -1)
        weight += unc;
    // down variation
    else if (idx == 2 || idx == -2)
        weight -= unc;
    
    if (weight < 0) {
        std::cout << "[WARNING] Muon (PTT= " << PTT << ", ETT= " << ETT <<") outside of trigger T&P kinematic range" << std::endl;
        return 1;
    }
    
    return weight;
}
    
    
#endif
