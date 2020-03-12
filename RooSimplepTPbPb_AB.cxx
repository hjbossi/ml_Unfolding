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

/*
  RooSimplepTPbPb_AB.cxx : Unfolds data corrected using the area based method. Data/MC comes directly from jet extractor output.
  Hannah Bossi <hannah.bossi@yale.edu>
  3/11/2020
 */



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

void RooSimplepTPbPb_AB(TString cFiles2="files1.txt")
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

  Double_t xbins[14];
  xbins[0]=5;
  xbins[1]=10;
  xbins[2]=20;
  xbins[3]=30;
  xbins[4]=40;
  xbins[5]=50;
  xbins[6]=60;
  xbins[7]=70;
  xbins[8]=80;
  xbins[9]=100;
  xbins[10]=120;
  xbins[11]=140;
  xbins[12]=190;
  xbins[13]=250;
  
  //the raw correlation (data or psuedodata)
  TH1F *h1raw(0);
  h1raw=new TH1F("r","raw", 17, 35, 120);
  //detector measure level (reco or hybrid MC)
  TH1F *h1smeared(0);
  h1smeared=new TH1F("smeared","smeared",17, 35, 120);

  //detector measure level no cuts
  TH1F *h1smearednocuts(0);
  h1smearednocuts=new TH1F("smearednocuts","smearednocuts", 13, xbins);
  //true correlations with measured cuts
  TH1F *h1true(0);
  h1true=new TH1F("true","true", 13, xbins);
  //full true correlation
  TH1F *h1fulleff(0);
  h1fulleff=new TH1F("truef","truef", 13, xbins); 

  
  TH2F *hcovariance(0);
  hcovariance=new TH2F("covariance","covariance",10,0.,1.,10,0,1.);

  TH2F *effnum=(TH2F*)h1fulleff->Clone("effnum");
  TH2F *effdenom=(TH2F*)h1fulleff->Clone("effdenom");
 
  effnum->Sumw2();
  effdenom->Sumw2();
  h1smeared->Sumw2();
  h1true->Sumw2();
  h1raw->Sumw2();
  h1fulleff->Sumw2();

  //branches in the tree that you need in this analysis
  Float_t ptJet,ptJetMatch, pTRec, ptdet;
  Float_t centData, cent;
  Int_t nEv=0;; 
  //so mcr is correctly normalized to one, not the response.       
  cout<<"cucu"<<endl;

  
  ifstream infile2;
  infile2.open(cFiles2.Data());
  char filename2[300];
  RooUnfoldResponse response;
  RooUnfoldResponse responsenotrunc;
  response.Setup(h1smeared,h1true);
  responsenotrunc.Setup(h1smearednocuts,h1fulleff);

  TFile *input1=TFile::Open("/home/hbossi/alidock/data/Data_031020/AnalysisResults.root");
  TTree *data=(TTree*)input1->Get("JetTree_AliAnalysisTaskJetExtractor_Jets_AKTFullR040_tracks_pT0150_caloClusters_E0300_pt_scheme_Rho_Scaled_Jets");
  Int_t nEvents=data->GetEntries();
  std::cout << nEvents << std::endl;
  data->SetBranchAddress("Jet_Pt", &pTRec);
  data->SetBranchAddress("Event_Centrality", &centData); 
  Int_t numTracksData;
  Float_t jetTrackPtData[400];
  data->SetBranchAddress("Jet_Track_Pt", &jetTrackPtData);
  data->SetBranchAddress("Jet_NumTracks", &numTracksData);
  for(Int_t i = 0; i < nEvents; i++){
    data->GetEntry(i);
    Double_t maxPtData = 0;
    for(Int_t index = 0; index < numTracksData; index++){
      if(jetTrackPtData[index] > maxPtData){
	maxPtData = jetTrackPtData[index];
      }
    }
    //std::cout << "max track pt : " << maxPtData << " JetPt: " << pTRec << std::endl;
    if(maxPtData < 7.0 ) continue;
    if(centData > 10) continue; 
    if(pTRec>120 || pTRec<35) continue;
    h1raw->Fill(pTRec);
  }

  
    while(infile2>>filename2){
    int pthardbin=0;
   
    TFile *input=TFile::Open(filename2);
    TList *list=(TList*) input->Get("AliAnalysisTaskEmcalEmbeddingHelper_histos");
    TList *list2=(TList*) list->FindObject("EventCuts");
    TH1D *hcent=(TH1D*)list2->FindObject("Centrality_selected");
    TProfile *hcross=(TProfile*)list->FindObject("fHistXsection");
    TH1D *htrials=(TH1D*)list->FindObject("fHistTrials");
    TH1D *hpthard=(TH1D*)list->FindObject("fHistPtHard");
    TH1D *hnevent=(TH1D*)list->FindObject("fHistEventCount");
    for(Int_t i=1;i<=htrials->GetNbinsX();i++){
      if(htrials->GetBinContent(i)!=0) pthardbin=i;}
    double pTHardscalefactor=(hcross->Integral(pthardbin,pthardbin)*hcross->GetEntries())/htrials->Integral(pthardbin,pthardbin);
    std::cout << "pT Hard Bin: " << pthardbin << "with scaling factor " << pTHardscalefactor << std::endl;
    TFile *input2=TFile::Open(filename2);
    TTree *mc=(TTree*)input2->Get("JetTree_AliAnalysisTaskJetExtractor_hybridLevelJets_AKTFullR040_tracks_pT0150_caloClustersCombined_E0300_pt_scheme_Rho_Scaled_hybridLevelJets"); 
    Int_t nEv=mc->GetEntries(); 
    Int_t numTracks;
    mc->SetBranchAddress("Jet_Pt", &ptJet);
    mc->SetBranchAddress("Event_Centrality", &cent);
    mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &ptJetMatch);
    //mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &ptdet);
    mc->SetBranchAddress("Jet_NumTracks", &numTracks);
    Float_t jetTrackPt[400];
    mc->SetBranchAddress("Jet_Track_Pt", &jetTrackPt);

    Int_t countm=0;
    for(int iEntry=0; iEntry< nEv; iEntry++){
      mc->GetEntry(iEntry);
      if(cent > 10) continue;
      double scalefactor = pTHardscalefactor;
      double EBscale = 1.;
      // put in if/else statements on the ptJet
      if(ptJet >= -20. && ptJet < 10.)       EBscale = 10.;
      else if(ptJet >= 10. && ptJet < 20.)   EBscale = 10.;
      else if(ptJet >= 20. && ptJet < 40.)   EBscale = 2.5;
      else if(ptJet >= 40. && ptJet < 60.)   EBscale = 1.25;
      else if(ptJet >= 60. && ptJet < 80.)   EBscale = 1.111;
      else if(ptJet >= 80. && ptJet < 100.)  EBscale = 1.111;
      else if(ptJet >= 100. && ptJet < 500.) EBscale = 1.0;

      scalefactor*=EBscale;
      //std::cout << ptJetMatch << std::endl;
      // LTB
      Double_t maxPt = 0;
      for(Int_t index = 0; index < numTracks; index++){
	if(jetTrackPt[index] > maxPt){
	  maxPt = jetTrackPt[index];
	}
      }
      // apply a 7 GeV leading track bias
      if(maxPt < 7.0 ) continue;      
      if(ptJetMatch> 250. || ptJetMatch < 5.) continue;
    
      h1fulleff->Fill(ptJetMatch,scalefactor);  
      h1smearednocuts->Fill(ptJet,scalefactor);  
      responsenotrunc.Fill(ptJet,ptJetMatch,scalefactor);
      
      if(ptJet>120 || ptJet<35) continue;
      h1smeared->Fill(ptJet,scalefactor);
      //this is the half split to be the response 
      response.Fill(ptJet, ptJetMatch,scalefactor);
      //this is the psuedo data!
      //h1raw->Fill(ptJet, scalefactor);
      //this is the generator level distribution for the pseudo data or our answer :)
      h1true->Fill(ptJetMatch,scalefactor);
      
    }}
 

    
    TH1F *htrueptd=(TH1F*) h1fulleff->Clone("trueptd");
    TH1F *htruept=(TH1F*) h1fulleff->Clone( "truept"); 
 
    //////////efficiencies done////////////////////////////////////
    TString filename = "UnfoldingDet_AreaBased_R04_Mar10.root"; 
    TFile *fout=new TFile (filename,"RECREATE");
    fout->cd();
    h1raw->SetName("raw");
    h1raw->Write();
    h1smeared->SetName("smeared");
    h1smeared->Write();
    htrueptd->Write();
    h1true->SetName("true");
    h1true->Write();
    TH1D* htruth = (TH1D*)response.Htruth();
    htruth->SetName("htruth");
    htruth->Write();
    response.Write();


    for(int jar=1;jar<10;jar++){
      Int_t iter=jar;
      cout<<"iteration"<<iter<<endl;
      cout<<"==============Unfold h1====================="<<endl;

      RooUnfoldBayes unfold(&response, h1raw, iter, false);    // OR
      TH1D* hunf= (TH1D*) unfold.Hreco(errorTreatment);
      //FOLD BACK
      TH1* hfold = response.ApplyToTruth (hunf, "");

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
int main () { RooSimplepTPbPb_AB(); return 0; }  // Main program when run stand-alone
#endif
