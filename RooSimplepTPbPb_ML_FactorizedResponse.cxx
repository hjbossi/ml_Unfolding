//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 279 2011-02-11 18:23:44Z T.J.Adye $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================


// RooSimplepTPbPb_ML_FactorizedResponse.cxx: Script to unfold the ML corrected data with a factorized Response
// Hannah Bossi <hannah.bossi@yale.edu>
// 3/4/2020




#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "TFile.h"
#include "TVectorD.h"

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TPostScript.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TProfile.h"
#inlcude "TRandom3.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldTestHarness2D.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

TH2D* CorrelationHistShape (const TMatrixD& cov,const char* name, const char* title,
		       Int_t na, Int_t nb, Int_t kbin)
{
 
   TH2D* h= new TH2D (name, title, nb, 0, nb, nb, 0, nb);
 
     	  for(int l=0;l<nb;l++){
                 for(int n=0;n<nb;n++){
                int index1=kbin+na*l;
                int index2=kbin+na*n;
   	        Double_t Vv=cov(index1,index1)*cov(index2,index2);
   		 if (Vv>0.0) h->SetBinContent(l+1,n+1,cov(index1,index2)/sqrt(Vv));


  		  }}
  return h;
}

TH2D* CorrelationHistPt (const TMatrixD& cov,const char* name, const char* title,
		       Int_t na, Int_t nb, Int_t kbin)
{
 
   TH2D* h= new TH2D (name, title, na, 0, na, na, 0, na);
 
     	  for(int l=0;l<na;l++){
                 for(int n=0;n<na;n++){

                int index1=l+na*kbin;
                int index2=n+na*kbin;
   	        Double_t Vv=cov(index1,index1)*cov(index2,index2);
   		 if (Vv>0.0)h->SetBinContent(l+1,n+1,cov(index1,index2)/sqrt(Vv));


  		  }}
  return h;
}


   





void Normalize2D(TH2* h)
{
   Int_t nbinsYtmp = h->GetNbinsY();
   const Int_t nbinsY = nbinsYtmp;
   Double_t norm[nbinsY];
   for(Int_t biny=1; biny<=nbinsY; biny++)
     {
       norm[biny-1] = 0;
       for(Int_t binx=1; binx<=h->GetNbinsX(); binx++)
     {
       norm[biny-1] += h->GetBinContent(binx,biny);
     }
     }

   for(Int_t biny=1; biny<=nbinsY; biny++)
     {
       for(Int_t binx=1; binx<=h->GetNbinsX(); binx++)
     {
       if(norm[biny-1]==0)  continue;
       else
         {
  h->SetBinContent(binx,biny,h->GetBinContent(binx,biny)/norm[biny-1]);
  h->SetBinError(binx,biny,h->GetBinError(binx,biny)/norm[biny-1]);
         }
     }
     }
}







TH2D* CorrelationHist (const TMatrixD& cov,const char* name, const char* title,
		       Double_t lo, Double_t hi,Double_t lon,Double_t hin)
{
  Int_t nb= cov.GetNrows();
  Int_t na= cov.GetNcols();
  cout<<nb<<" "<<na<<endl;
  TH2D* h= new TH2D (name, title, nb, 0, nb, na, 0, na);
  h->SetAxisRange (-1.0, 1.0, "Z");
  for(int i=0; i < na; i++)
  for(int j=0; j < nb; j++) {
  Double_t Viijj= cov(i,i)*cov(j,j);
      if (Viijj>0.0) h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
    }
  return h;
}

//==============================================================================
// Example Unfolding
//==============================================================================

void RooSimplepTPbPb_ML_FactorizedResponse(TString cFiles2="filesML.txt")
{
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif
  Int_t difference=1;
  Int_t Ppol=0;
  cout << "==================================== pick up the response matrix for background==========================" << endl;
  ///////////////////parameter setting
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
  
  TRandom3* rand = new TRandom3();

   //***************************************************

  Double_t xbins[13];
  xbins[0]=10;
  xbins[1]=20;
  xbins[2]=30;
  xbins[3]=40;
  xbins[4]=50;
  xbins[5]=60;
  xbins[6]=70;
  xbins[7]=80;
  xbins[8]=100;
  xbins[9]=120;
  xbins[10]=140;
  xbins[11]=190;
  xbins[12]=250;
  //xbins[13]=250;




  
   
  //the raw correlation (data or psuedodata)
  //TH1D *h1raw(0);
  //h1raw=new TH1D("r","raw", 19, 25, 120);

  //detector measure level (reco or hybrid MC)
  TH1D *h1smeared(0);
  h1smeared=new TH1D("smeared","smeared",19,25,120);

  // full range of reco
  TH1D *h1smearedFullRange(0);
  h1smearedFullRange = new TH1D("smearedFull", "smearedFull", 220, -20,200);
    
  //detector measure level no cuts
  TH1D *h1smearednocuts(0);
  h1smearednocuts=new TH1D("smearednocuts","smearednocuts", 19,25,120);
  //true correlations with measured cuts
  TH1D *h1true(0);
  h1true=new TH1D("true","true", 12, xbins);
  //full true correlation
  TH1D *h1fulleff(0);
  h1fulleff=new TH1D("truef","truef", 12, xbins); 
  
  TH2D *hcovariance(0);
  hcovariance=new TH2D("covariance","covariance",10,0.,1.,10,0,1.);

  TH2D *effnum=(TH2D*)h1fulleff->Clone("effnum");
  TH2D *effdenom=(TH2D*)h1fulleff->Clone("effdenom");
 
   effnum->Sumw2();
   effdenom->Sumw2();
   h1smeared->Sumw2();
   h1smearedFullRange->Sumw2(); 
   h1true->Sumw2();
   h1fulleff->Sumw2();

   //branches in the tree that you need in this analysis
   // we need the hybrid Pt to determine what EB scaling factor
   Float_t ptJetMatch, hybridPt, hybridPtData, ptJet;
   Double_t leadingTrackPtData, leadingTrackPt;
   Long64_t pTHardBin; // we are getting the pt hard bin from the tree this time
   Float_t cent;
   Double_t pTRec; // ml corrected data 
   
   //so mcr is correctly normalized to one, not the response.       
   cout<<"cucu"<<endl;
   RooUnfoldResponse response;
   RooUnfoldResponse responsenotrunc;
   response.Setup(h1smeared,h1true);
   responsenotrunc.Setup(h1smearednocuts,h1fulleff);


   TFile*_file0= TFile::Open("Unfolding_NeuralNetwork_Mar26th_FullStats_Det_RecoBinning.root");
   TH1D* h1raw = (TH1D*)_file0->Get("Bayesian_Unfoldediter5");
   h1raw->Sumw2();
   TH1D* effnm       = (TH1D*)_file0->Get("htruth"); // numerator for kinematic efficiecny          
   TH1D* effdm     = (TH1D*)_file0->Get("trueptd");// denominator for kinematic efficiency  
   effnm->Divide(effdm);
   for(Int_t i =0 ; i < h1raw->GetNbinsX()+1; i++){
     h1raw->SetBinContent(i, h1raw->GetBinContent(i)*(1./effnm->GetBinContent(i)));
   }


   // we don't need this instead we take the distribution from the unfolding file from before
   /*
   TFile *input1=TFile::Open("/home/alidock/PredictionTrees/predictionTree_NeuralNetwork_For_LHC15o_R040_032620_FullStats.root");
   TTree *data=(TTree*)input1->Get("NeuralNetwork_For_LHC15o_R040");
   Int_t nEvents=data->GetEntries();
   std::cout << nEvents << std::endl;
   data->SetBranchAddress("Jet_Pt", &hybridPtData); 
   data->SetBranchAddress("Predicted_Jet_Pt", &pTRec);
   data->SetBranchAddress("Jet_TrackPt0", &leadingTrackPtData); 
   for(Int_t i = 0; i < nEvents; i++){
     data->GetEntry(i);
     if(leadingTrackPtData > 100) continue; 
     double EBscale = 1.;
     if(hybridPtData >= -20. && hybridPtData < 10.)       EBscale = 10.;
     else if(hybridPtData >= 10. && hybridPtData < 20.)   EBscale = 10.;
     else if(hybridPtData >= 20. && hybridPtData < 40.)   EBscale = 2.5;
     else if(hybridPtData >= 40. && hybridPtData < 60.)   EBscale = 1.25;
     else if(hybridPtData >= 60. && hybridPtData < 80.)   EBscale = 1.111;
     else if(hybridPtData >= 80. && hybridPtData < 100.)  EBscale = 1.111;
     else if(hybridPtData >= 100. && hybridPtData < 500.) EBscale = 1.0;
     if(pTRec>120 || pTRec<25) continue;
     h1raw->Fill(pTRec, EBscale); 
   }
   */
   // previously derived pT hard bin scaling factors INCLUDING rk. index with PtHardBin Branch
   Double_t scalingFactors[20] = {0.492069, 0.418282, 0.405473, 0.265289, 0.133235, 0.0642945, 0.0229219, 0.008419, 0.00360539, 0.00122459, 0.000494578, 0.000177912, 8.83422e-05, 4.28423e-05, 2.03583e-05, 1.03874e-05, 5.68744e-06, 2.99194e-06, 1.59651e-06, 2.09335e-06};

   TFile *input2=TFile::Open("/home/alidock/PredictionTrees/predictionTree_NeuralNetwork_For_LHC16j5_Embedded_R040_032620_FullStats.root");
   TTree *mc=(TTree*)input2->Get("NeuralNetwork_For_LHC16j5_Embedded_R040"); 
   Int_t nEv=mc->GetEntries(); 
   // get the jet pT predicted by the ml
   mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &ptJet); 
   mc->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &ptJetMatch);
   mc->SetBranchAddress("PtHardBin", &pTHardBin);
   mc->SetBranchAddress("Jet_Pt", &hybridPt);
   mc->SetBranchAddress("Event_Centrality", &cent);
   mc->SetBranchAddress("Jet_TrackPt0", &leadingTrackPt); 
   Int_t countm=0;
   for(int iEntry=0; iEntry< nEv; iEntry++){
     mc->GetEntry(iEntry);
     if (cent > 10) continue;
     if (leadingTrackPt > 100) continue; 
     double scalefactor = scalingFactors[pTHardBin-1];
     double EBscale = 1.;
     // put in if/else statements on the ptJet 
     if(hybridPt >= -20. && hybridPt < 10.)       EBscale = 10.;
     else if(hybridPt >= 10. && hybridPt < 20.)   EBscale = 10.;
     else if(hybridPt >= 20. && hybridPt < 40.)   EBscale = 2.5;
     else if(hybridPt >= 40. && hybridPt < 60.)   EBscale = 1.25;
     else if(hybridPt >= 60. && hybridPt < 80.)   EBscale = 1.111;
     else if(hybridPt >= 80. && hybridPt < 100.)  EBscale = 1.111;
     else if(hybridPt >= 100. && hybridPt < 500.) EBscale = 1.0;

     scalefactor*=EBscale; 
     if(ptJetMatch < 10 || ptJetMatch > 250) continue;
     h1fulleff->Fill(ptJetMatch,scalefactor);  
     h1smearedFullRange->Fill(ptJet, scalefactor);
     h1smearednocuts->Fill(ptJet,scalefactor);  
     responsenotrunc.Fill(ptJet,ptJetMatch,scalefactor);
     if(ptJet>250 || ptJet<10) continue;
     h1smeared->Fill(ptJet,scalefactor);
     //this is the half split to be the response 
     response.Fill(ptJet, ptJetMatch,scalefactor);
     //this is the generator level distribution for the pseudo data or our answer :)
     h1true->Fill(ptJetMatch,scalefactor);
      
   }
 

    
    TH1D *htrueptd=(TH1D*) h1fulleff->Clone("trueptd");
    TH1D *htruept=(TH1D*) h1fulleff->Clone( "truept"); 
 
    //////////efficiencies done////////////////////////////////////
 
    TFile *fout=new TFile (Form("Unfolding_NeuralNetwork_Mar26th_FactorizedResponse_RecoBinning.root"),"RECREATE");
    fout->cd();
    h1raw->SetName("raw");
    h1raw->Write();
    h1smeared->SetName("smeared");
    h1smeared->Write();
    h1smearedFullRange->Write();
    htrueptd->Write();
    h1true->SetName("true");
    h1true->Write();
    response.Write();
    TH1D* htruth = (TH1D*)response.Htruth();
    htruth->SetName("htruth");
    htruth->Write();
    for(int jar=1;jar<10;jar++){
      Int_t iter=jar;
      cout<<"iteration"<<iter<<endl;
      cout<<"==============Unfold h1====================="<<endl;

      RooUnfoldBayes unfold(&response, h1raw, iter, false);    // OR
      TH1D* hunf= (TH1D*) unfold.Hreco(errorTreatment);
      //FOLD BACK
      TH1* hfold = response.ApplyToTruth(hunf, "");

      TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
      htempUnf->SetName(Form("Bayesian_Unfoldediter%d",iter));
      
      TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
      htempFold->SetName(Form("Bayesian_Foldediter%d",iter));        

      htempUnf->Write();
      htempFold->Write();
	  
      ///HERE I GET THE COVARIANCE MATRIX/////
      /*
      if(iter==8){
	TMatrixD covmat= unfold.Ereco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
	for(Int_t k=0;k<h2true->GetNbinsX();k++){
	  TH2D *hCorr= (TH2D*) CorrelationHistShape(covmat, Form("corr%d",k), "Covariance matrix",h2true->GetNbinsX(),h2true->GetNbinsY(),k);
	  TH2D *covshape=(TH2D*)hCorr->Clone("covshape");      
	  covshape->SetName(Form("pearsonmatrix_iter%d_binshape%d",iter,k));
	  covshape->SetDrawOption("colz");
	  covshape->Write();
	}
	
	for(Int_t k=0;k<h2true->GetNbinsY();k++){
	  TH2D *hCorr= (TH2D*) CorrelationHistPt(covmat, Form("corr%d",k), "Covariance matrix",h2true->GetNbinsX(),h2true->GetNbinsY(),k);
	  TH2D *covpt=(TH2D*)hCorr->Clone("covpt");      
	  covpt->SetName(Form("pearsonmatrix_iter%d_binpt%d",iter,k));
	  covpt->SetDrawOption("colz");
	  covpt->Write();
	  }
	  }*/
    }
	  
}
#ifndef __CINT__
int main () { RooSimplepTPbPb_ML_FactorizedResponse(); return 0; }  // Main program when run stand-alone
#endif
