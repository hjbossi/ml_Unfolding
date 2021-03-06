
// Hannah Bossi <hannah.bossi@cern.ch>
// unfoldingTrueRatio_SVD.C : Plots the ratio of unfolded distributions to k = 4
// 02/13/2020

void unfoldingTrueRatio_SVD(){
  // plotting styles
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //whether or not to apply the efficiency corrections.
  Bool_t correctForEff = kTRUE;
  // draw only the middle four iterations. 
  Bool_t drawMiddle4   = kTRUE; 
  Bool_t plotRecRange  = kFALSE; //whether or not to only show the range we will end up reporting 
  
  // Get the input file and the relevant objects to make the ratio plots
  TFile*_file0       = TFile::Open("../Unfolding_NeuralNetwork_R04_Part_SVD_Mar13th.root");
  TH1F* trueDist     = (TH1F*)_file0->Get("true");
  TH1F* unfold_Iter1 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_1");
  TH1F* unfold_Iter2 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_2");
  TH1F* unfold_Iter3 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_3");
  TH1F* unfold_Iter4 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_4");
  TH1F* unfold_Iter5 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_5");
  TH1F* unfold_Iter6 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_6");
  TH1F* unfold_Iter7 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_7");
  TH1F* unfold_Iter8 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_8");
  TH1F* unfold_Iter9 = (TH1F*)_file0->Get("SVD_Unfolded_KValue_9");
  TH1F* effnum       = (TH1F*)_file0->Get("htruth"); // numerator for kinematic efficiecny
  TH1F* effdenom     = (TH1F*)_file0->Get("trueptd");// denominator for kinematic efficiency
  effnum->Divide(effdenom);
  
  int colors[20] = {kRed+2, kRed-4, kOrange+7, kOrange, kYellow-4, kSpring+10, kSpring, kGreen-3, kGreen+3, kTeal-7, kTeal, kAzure+10, kAzure-4, kBlue+2, kViolet+8, kViolet-1, kMagenta+1, kMagenta-4, kPink+7, kPink-4};
  
  // Define the Canvas
   TCanvas *c = new TCanvas("c", "canvas", 800, 800);

   // Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); 
   pad1->SetGridx();
   pad1->Draw();
   pad1->cd();
   pad1->SetLogy();
   
   trueDist->SetMarkerStyle(20);
   trueDist->SetMarkerColor(kBlue);
   trueDist->SetLineColor(kBlue);
   if(correctForEff){
     for(Int_t i =0 ; i < trueDist->GetNbinsX()+1; i++){
       trueDist->SetBinContent(i, trueDist->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   trueDist->Scale(1./trueDist->Integral(), "width");
   if(plotRecRange)trueDist->GetXaxis()->SetRangeUser(40,120);
   trueDist->Draw();

   unfold_Iter1->SetMarkerStyle(20);
   unfold_Iter1->SetMarkerColor(colors[0]);
   unfold_Iter1->SetLineColor(colors[0]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter1->GetNbinsX()+1; i++){
       unfold_Iter1->SetBinContent(i, unfold_Iter1->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter1->Scale(1./unfold_Iter1->Integral(), "width");
   if(!drawMiddle4) unfold_Iter1->Draw("same");

   unfold_Iter2->SetMarkerStyle(20);
   unfold_Iter2->SetMarkerColor(colors[2]);
   unfold_Iter2->SetLineColor(colors[2]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter2->GetNbinsX()+1; i++){
       unfold_Iter2->SetBinContent(i, unfold_Iter2->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter2->Scale(1./unfold_Iter2->Integral(), "width");
   if(!drawMiddle4) unfold_Iter2->Draw("same");

   unfold_Iter3->SetMarkerStyle(20);
   unfold_Iter3->SetMarkerColor(colors[4]);
   unfold_Iter3->SetLineColor(colors[4]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter3->GetNbinsX()+1; i++){
       unfold_Iter3->SetBinContent(i, unfold_Iter3->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter3->Scale(1./unfold_Iter3->Integral(), "width");
   unfold_Iter3->Draw("same");
   
   unfold_Iter4->SetMarkerStyle(20);
   unfold_Iter4->SetMarkerColor(colors[6]);
   unfold_Iter4->SetLineColor(colors[6]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter4->GetNbinsX()+1; i++){
       unfold_Iter4->SetBinContent(i, unfold_Iter4->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter4->Scale(1./unfold_Iter4->Integral(), "width");
   unfold_Iter4->Draw("same");

   unfold_Iter5->SetMarkerStyle(20);
   unfold_Iter5->SetMarkerColor(colors[8]);
   unfold_Iter5->SetLineColor(colors[8]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter5->GetNbinsX()+1; i++){
       unfold_Iter5->SetBinContent(i, unfold_Iter5->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter5->Scale(1./unfold_Iter5->Integral(), "width");
   unfold_Iter5->Draw("same");

   unfold_Iter6->SetMarkerStyle(20);
   unfold_Iter6->SetMarkerColor(colors[10]);
   unfold_Iter6->SetLineColor(colors[10]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter6->GetNbinsX()+1; i++){
       unfold_Iter6->SetBinContent(i, unfold_Iter6->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter6->Scale(1./unfold_Iter6->Integral(), "width");
   unfold_Iter6->Draw("same");

   unfold_Iter7->SetMarkerStyle(20);
   unfold_Iter7->SetMarkerColor(colors[12]);
   unfold_Iter7->SetLineColor(colors[12]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter7->GetNbinsX()+1; i++){
       unfold_Iter7->SetBinContent(i, unfold_Iter7->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter7->Scale(1./unfold_Iter7->Integral(), "width");
   unfold_Iter7->Draw("same");

   unfold_Iter8->SetMarkerStyle(20);
   unfold_Iter8->SetMarkerColor(colors[14]);
   unfold_Iter8->SetLineColor(colors[14]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter8->GetNbinsX()+1; i++){
       unfold_Iter8->SetBinContent(i, unfold_Iter8->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter8->Scale(1./unfold_Iter8->Integral(), "width");
   if(!drawMiddle4) unfold_Iter8->Draw("same");

   unfold_Iter9->SetMarkerStyle(20);
   unfold_Iter9->SetMarkerColor(colors[16]);
   unfold_Iter9->SetLineColor(colors[16]);
   if(correctForEff){
     for(Int_t i =0 ; i < unfold_Iter9->GetNbinsX()+1; i++){
       unfold_Iter9->SetBinContent(i, unfold_Iter9->GetBinContent(i)*(1./effnum->GetBinContent(i)));
     }
   }
   unfold_Iter9->Scale(1./unfold_Iter9->Integral(), "width");
   if(!drawMiddle4) unfold_Iter9->Draw("same");
   
  // want to plot the ratio between the unfolded distributions and the raw distribution
  // do this for each of the iterationas on the same plot
  // process in general goes as (1) Clone numerator (2) Sumw2 (3) Divide by denominator
  TH1F* unfoldClone_1 = (TH1F*) unfold_Iter1->Clone();
  unfoldClone_1->Sumw2();
  unfoldClone_1->Divide(unfold_Iter5);
  unfoldClone_1->SetMarkerStyle(20);
  unfoldClone_1->SetMarkerColor(colors[0]);

  TH1F* unfoldClone_2 = (TH1F*) unfold_Iter2->Clone();
  unfoldClone_2->Sumw2();
  unfoldClone_2->Divide(unfold_Iter5);
  unfoldClone_2->SetMarkerStyle(20);
  unfoldClone_2->SetMarkerColor(colors[2]);

  TH1F* unfoldClone_3 = (TH1F*) unfold_Iter3->Clone();
  unfoldClone_3->Sumw2();
  unfoldClone_3->Divide(unfold_Iter5);
  unfoldClone_3->SetMarkerStyle(20);
  unfoldClone_3->SetMarkerColor(colors[4]);
  
  TH1F* unfoldClone_4 = (TH1F*) unfold_Iter4->Clone();
  unfoldClone_4->Sumw2();
  unfoldClone_4->Divide(unfold_Iter5);
  unfoldClone_4->SetMarkerStyle(20);
  unfoldClone_4->SetMarkerColor(colors[6]);

  TH1F* unfoldClone_5 = (TH1F*) unfold_Iter5->Clone();
  unfoldClone_5->Sumw2();
  unfoldClone_5->Divide(unfold_Iter5);
  unfoldClone_5->SetMarkerStyle(20);
  unfoldClone_5->SetMarkerColor(colors[8]);
  
  TH1F* unfoldClone_6 = (TH1F*) unfold_Iter6->Clone();
  unfoldClone_6->Sumw2();
  unfoldClone_6->Divide(unfold_Iter5);
  unfoldClone_6->SetMarkerStyle(20);
  unfoldClone_6->SetMarkerColor(colors[10]);

  TH1F* unfoldClone_7 = (TH1F*) unfold_Iter7->Clone();
  unfoldClone_7->Sumw2();
  unfoldClone_7->Divide(unfold_Iter5);
  unfoldClone_7->SetMarkerStyle(20);
  unfoldClone_7->SetMarkerColor(colors[12]);

  TH1F* unfoldClone_8 = (TH1F*) unfold_Iter8->Clone();
  unfoldClone_8->Sumw2();
  unfoldClone_8->Divide(unfold_Iter5);
  unfoldClone_8->SetMarkerStyle(20);
  unfoldClone_8->SetMarkerColor(colors[14]);

  TH1F* unfoldClone_9 = (TH1F*) unfold_Iter9->Clone();
  unfoldClone_9->Sumw2();
  unfoldClone_9->Divide(unfold_Iter5);
  unfoldClone_9->SetMarkerStyle(20);
  unfoldClone_9->SetMarkerColor(colors[16]);

  
  c->cd();
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx();
  pad2->Draw();
  pad2->cd();
  
  unfoldClone_3->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
  unfoldClone_3->GetXaxis()->SetLabelSize(0.08);
  unfoldClone_3->GetXaxis()->SetTitleSize(0.09);
  unfoldClone_3->GetYaxis()->SetTitle("Ratio Unfolded/Unfolded k = 5");
  unfoldClone_3->GetYaxis()->SetLabelSize(0.06);
  unfoldClone_3->GetYaxis()->SetTitleSize(0.08);
  unfoldClone_3->GetYaxis()->SetTitleOffset(0.55);
  unfoldClone_3->GetYaxis()->SetRangeUser(0.8, 1.2);
  if(plotRecRange)unfoldClone_3->GetXaxis()->SetRangeUser(40,120);
  unfoldClone_3->Draw();
  if(!drawMiddle4)unfoldClone_2->Draw("same");
  if(!drawMiddle4)unfoldClone_1->Draw("same");
  unfoldClone_4->Draw("same");
  unfoldClone_5->Draw("same");
  unfoldClone_6->Draw("same");
  unfoldClone_7->Draw("same");
  if(!drawMiddle4)unfoldClone_8->Draw("same");
  if(!drawMiddle4)unfoldClone_9->Draw("same");

  //effnum->Draw();
  TLegend* leg = new TLegend(0.7, 0.5, 0.9, 0.9);
  leg->SetTextSize(0.04);
  leg->AddEntry(trueDist, "True");
  if(!drawMiddle4)leg->AddEntry(unfoldClone_1, "k = 1");
  if(!drawMiddle4)leg->AddEntry(unfoldClone_2, "k = 2");
  leg->AddEntry(unfoldClone_3, "k = 3");
  leg->AddEntry(unfoldClone_4, "k = 4");
  leg->AddEntry(unfoldClone_5, "k = 5");
  leg->AddEntry(unfoldClone_6, "k = 6");
  leg->AddEntry(unfoldClone_7, "k = 7");
  if(!drawMiddle4)leg->AddEntry(unfoldClone_8, "k = 8");
  if(!drawMiddle4)leg->AddEntry(unfoldClone_9, "k = 9");
  pad1->cd();
  leg->Draw("same");


}
