// Hannah Bossi <hannah.bossi@cern.ch>
// refoldedRawRatio.C : Plots the ratio of refolded to raw distributions after unfolding
// 02/13/2020

void refoldedRawRatio(){
  // plotting styles
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  // Get the input file and the relevant objects to make the ratio plots
  TFile*_file0       = TFile::Open("../Unfolding_NeuralNetwork_R04_Det_Mar11th.root");
  TH1F* rawDist      = (TH1F*)_file0->Get("raw");
  TH1F* refold_Iter1 = (TH1F*)_file0->Get("Bayesian_Foldediter1");
  TH1F* refold_Iter2 = (TH1F*)_file0->Get("Bayesian_Foldediter2");
  TH1F* refold_Iter3 = (TH1F*)_file0->Get("Bayesian_Foldediter3");
  TH1F* refold_Iter4 = (TH1F*)_file0->Get("Bayesian_Foldediter4");
  TH1F* refold_Iter5 = (TH1F*)_file0->Get("Bayesian_Foldediter5");
  TH1F* refold_Iter6 = (TH1F*)_file0->Get("Bayesian_Foldediter6");
  TH1F* refold_Iter7 = (TH1F*)_file0->Get("Bayesian_Foldediter7");
  TH1F* refold_Iter8 = (TH1F*)_file0->Get("Bayesian_Foldediter8");
  TH1F* refold_Iter9 = (TH1F*)_file0->Get("Bayesian_Foldediter9");


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
   
   rawDist->SetMarkerStyle(20);
   rawDist->SetMarkerColor(kBlue);
   rawDist->SetLineColor(kBlue);
   rawDist->Scale(1., "width");
   rawDist->Draw("HIST");

   refold_Iter1->SetMarkerStyle(20);
   refold_Iter1->SetMarkerColor(colors[0]);
   refold_Iter1->SetLineColor(colors[0]);
   refold_Iter1->Scale(1., "width");
   refold_Iter1->Draw("same");

   refold_Iter2->SetMarkerStyle(20);
   refold_Iter2->SetMarkerColor(colors[2]);
   refold_Iter2->SetLineColor(colors[2]);
   refold_Iter2->Scale(1., "width");
   refold_Iter2->Draw("same");

   refold_Iter3->SetMarkerStyle(20);
   refold_Iter3->SetMarkerColor(colors[4]);
   refold_Iter3->SetLineColor(colors[4]);
   refold_Iter3->Scale(1., "width");
   refold_Iter3->Draw("same");
   
   refold_Iter4->SetMarkerStyle(20);
   refold_Iter4->SetMarkerColor(colors[6]);
   refold_Iter4->SetLineColor(colors[6]);
   refold_Iter4->Scale(1., "width");
   refold_Iter4->Draw("same");

   refold_Iter5->SetMarkerStyle(20);
   refold_Iter5->SetMarkerColor(colors[8]);
   refold_Iter5->SetLineColor(colors[8]);
   refold_Iter5->Scale(1., "width");
   refold_Iter5->Draw("same");

   refold_Iter6->SetMarkerStyle(20);
   refold_Iter6->SetMarkerColor(colors[10]);
   refold_Iter6->SetLineColor(colors[10]);
   refold_Iter6->Scale(1., "width");
   refold_Iter6->Draw("same");

   refold_Iter7->SetMarkerStyle(20);
   refold_Iter7->SetMarkerColor(colors[12]);
   refold_Iter7->SetLineColor(colors[12]);
   refold_Iter7->Scale(1., "width");
   refold_Iter7->Draw("same");

   refold_Iter8->SetMarkerStyle(20);
   refold_Iter8->SetMarkerColor(colors[14]);
   refold_Iter8->SetLineColor(colors[14]);
   refold_Iter8->Scale(1., "width");
   refold_Iter8->Draw("same");

   refold_Iter9->SetMarkerStyle(20);
   refold_Iter9->SetMarkerColor(colors[16]);
   refold_Iter9->SetLineColor(colors[16]);
   refold_Iter9->Scale(1., "width");
   refold_Iter9->Draw("same");
   
  // want to plot the ratio between the refolded distributions and the raw distribution
  // do this for each of the iterationas on the same plot
  // process in general goes as (1) Clone numerator (2) Sumw2 (3) Divide by denominator
  TH1F* refoldClone_1 = (TH1F*) refold_Iter1->Clone();
  refoldClone_1->Sumw2();
  refoldClone_1->Divide(rawDist);
  refoldClone_1->SetMarkerStyle(20);
  refoldClone_1->SetMarkerColor(colors[0]);

  TH1F* refoldClone_2 = (TH1F*) refold_Iter2->Clone();
  refoldClone_2->Sumw2();
  refoldClone_2->Divide(rawDist);
  refoldClone_2->SetMarkerStyle(20);
  refoldClone_2->SetMarkerColor(colors[2]);

  TH1F* refoldClone_3 = (TH1F*) refold_Iter3->Clone();
  refoldClone_3->Sumw2();
  refoldClone_3->Divide(rawDist);
  refoldClone_3->SetMarkerStyle(20);
  refoldClone_3->SetMarkerColor(colors[4]);
  
  TH1F* refoldClone_4 = (TH1F*) refold_Iter4->Clone();
  refoldClone_4->Sumw2();
  refoldClone_4->Divide(rawDist);
  refoldClone_4->SetMarkerStyle(20);
  refoldClone_4->SetMarkerColor(colors[6]);

  TH1F* refoldClone_5 = (TH1F*) refold_Iter5->Clone();
  refoldClone_5->Sumw2();
  refoldClone_5->Divide(rawDist);
  refoldClone_5->SetMarkerStyle(20);
  refoldClone_5->SetMarkerColor(colors[8]);
  
  TH1F* refoldClone_6 = (TH1F*) refold_Iter6->Clone();
  refoldClone_6->Sumw2();
  refoldClone_6->Divide(rawDist);
  refoldClone_6->SetMarkerStyle(20);
  refoldClone_6->SetMarkerColor(colors[10]);

  TH1F* refoldClone_7 = (TH1F*) refold_Iter7->Clone();
  refoldClone_7->Sumw2();
  refoldClone_7->Divide(rawDist);
  refoldClone_7->SetMarkerStyle(20);
  refoldClone_7->SetMarkerColor(colors[12]);

  TH1F* refoldClone_8 = (TH1F*) refold_Iter8->Clone();
  refoldClone_8->Sumw2();
  refoldClone_8->Divide(rawDist);
  refoldClone_8->SetMarkerStyle(20);
  refoldClone_8->SetMarkerColor(colors[14]);

  TH1F* refoldClone_9 = (TH1F*) refold_Iter9->Clone();
  refoldClone_9->Sumw2();
  refoldClone_9->Divide(rawDist);
  refoldClone_9->SetMarkerStyle(20);
  refoldClone_9->SetMarkerColor(colors[16]);

  c->cd();
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx();
  pad2->Draw();
  pad2->cd();
  
  refoldClone_1->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
  refoldClone_1->GetXaxis()->SetLabelSize(0.08);
  refoldClone_1->GetXaxis()->SetTitleSize(0.09);
  refoldClone_1->GetYaxis()->SetTitle("Ratio Refolded/Raw");
  refoldClone_1->GetYaxis()->SetLabelSize(0.08);
  refoldClone_1->GetYaxis()->SetTitleSize(0.08);
  refoldClone_1->GetYaxis()->SetTitleOffset(0.55);
  refoldClone_1->GetYaxis()->SetRangeUser(0,2);
  refoldClone_1->Draw();
  refoldClone_2->Draw("same");
  refoldClone_3->Draw("same");
  refoldClone_4->Draw("same");
  refoldClone_5->Draw("same");
  refoldClone_6->Draw("same");
  refoldClone_7->Draw("same");
  refoldClone_8->Draw("same");
  refoldClone_9->Draw("same");

  TLegend* leg = new TLegend(0.7, 0.5, 0.9, 0.9);
  leg->SetTextSize(0.04);
  leg->AddEntry(rawDist, "Raw");
  leg->AddEntry(refoldClone_1, "Iteration 1");
  leg->AddEntry(refoldClone_2, "Iteration 2");
  leg->AddEntry(refoldClone_3, "Iteration 3");
  leg->AddEntry(refoldClone_4, "Iteration 4");
  leg->AddEntry(refoldClone_5, "Iteration 5");
  leg->AddEntry(refoldClone_6, "Iteration 6");
  leg->AddEntry(refoldClone_7, "Iteration 7");
  leg->AddEntry(refoldClone_8, "Iteration 8");
  leg->AddEntry(refoldClone_9, "Iteration 9");
  pad1->cd();
  leg->Draw("same");


}
