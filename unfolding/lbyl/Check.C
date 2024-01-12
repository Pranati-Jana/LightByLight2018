#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>

#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "THStack.h"
#include "TCut.h"
#include "TAxis.h"
#include "TChain.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "RooUnfoldInvert.h"
#endif


void Check(){

   TFile *f=new TFile("unfolding_histograms_Bayes_SC.root","r");
   TH1D *hdipho_GenInvmass  = (TH1D*)f->Get("hGenInvmass_xSec");
   TH1D *hdipho_Invmass     = (TH1D*)f->Get("hRecoMCInvmass_xSec");


   new TCanvas();
   hdipho_GenInvmass->Draw();
   new TCanvas();
   hdipho_Invmass->Draw();


TH2D *h2D_Invmass = new TH2D("h2D_Invmass", "2D Invmassidity Histogram", 3, 5,17,3,5,17);

for (Int_t i = 1; i <= 3; ++i) {
    for (Int_t j = 1; j <= 3; ++j) {
        Double_t contentGenInvmass = hdipho_GenInvmass->GetBinContent(i);
        Double_t contentRecoInvmass = hdipho_Invmass->GetXaxis()->GetBinCenter(j);
 std::cout << "Bin (" << i << ", " << j << "): " << contentGenInvmass << ", " << contentRecoInvmass << std::endl;
        h2D_Invmass->Fill(contentRecoInvmass, contentGenInvmass);
    }
}
 new TCanvas();
 
  h2D_Invmass->Draw("COLZ");

  RooUnfoldResponse responseInvmass (hdipho_Invmass, hdipho_GenInvmass,h2D_Invmass);
 // RooUnfoldResponse responseInvmass (hdipho_Invmass, hdipho_GenInvmass);
  


  RooUnfoldBayes unfold2(&responseInvmass, hdipho_Invmass, 4);

  TCanvas*c6 =  new TCanvas("c6","Response Matrix3", 800,600);
        c6->SetRightMargin(0.15);
        TMatrixD  responseMatrix3 = responseInvmass.Mresponse();
        TH2D *hResponseMatrix3 = new TH2D("hResponseMatrix3", "Response Matrix3;Reco invariant mass;Gen invariant mass", responseMatrix3.GetNrows(), 0, responseMatrix3.GetNrows(), responseMatrix3.GetNcols(), 0, responseMatrix3.GetNcols());
        for (int i = 0; i < responseMatrix3.GetNrows(); ++i) {
             for (int j = 0; j < responseMatrix3.GetNcols(); ++j) {
                  hResponseMatrix3->SetBinContent(i + 1, j + 1, responseMatrix3[i][j]);
                  }
        }
        hResponseMatrix3->Draw("textcolz");
        c6->Update();
        c6->SaveAs("fig_Bayes_Inverse/InvMass.png");


        TCanvas*c8 =  new TCanvas();
        TMatrixD  responseMatrix5 = responseInvmass.Mresponse(true);
        responseMatrix5.Draw("textcolz");
        c8->SaveAs("fig_Bayes_Inverse/CovMatrixMass.png"); 
  
  TH1D* hUnfoldInvmass= (TH1D*) unfold2.Hunfold();
  unfold2.PrintTable (cout,hdipho_GenInvmass);

  TH1D *hdipho_Invmass_data        = (TH1D*)f->Get("hRecoDataInvmass_xSec");
  RooUnfoldBayes unfold_dataInvmass(&responseInvmass, hdipho_Invmass_data , 4);

  TH1D* hUnfoldInvmass_data= (TH1D*) unfold_dataInvmass.Hunfold();
  unfold_dataInvmass.PrintTable (cout, hdipho_GenInvmass);
}
