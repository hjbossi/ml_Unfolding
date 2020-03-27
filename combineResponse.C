/* combineResponse.C: 
 * Reads the response from two root files and multiplies them 
 * writing the new response to a root file.
 * 
 *
 */

void combineResponse(){
  // First get the responses from the files
  TFile *_file0 = TFile::Open("Unfolding_NeuralNetwork_Mar26th_FullStats_Det.root"); 
  TFile *_file1 = TFile::Open("Unfolding_NeuralNetwork_Mar26th_FactorizedResponse.root");

  // get the fluctuations response and the detector response
  RooUnfoldResponse* flucResp = (RooUnfoldResponse*)_file0->Get("smeared_true"); 
  RooUnfoldResponse* detResp  = (RooUnfoldResponse*)_file1->Get("smeared_true"); 
  TH2D* flucHist     = (TH2D*)flucResp->Hresponse();
  TH2D* detHist      = (TH2D*)detResp->Hresponse();
  TH1D* truthHist    = (TH1D*)detResp->Htruth();
  TH1D* measuredHist = (TH1D*)flucResp->Hmeasured();


  // the fluctuation histogram has the same binning that we want in our final response
  // copy this hist and set all elements to zero
  // note this may change depending on intermedeate binning you use
  TH2D* final = (TH2D*)flucHist->Clone();
  final->SetName("final");
  final->Reset();


  // work on the fluctuations hist - Don't need to reweight, done in response automatically :) 
  Int_t numBinsY_Fluctuations = flucHist->GetNbinsY(); 
  Int_t numBinsX_Fluctuations = flucHist->GetNbinsX(); 

  for(Int_t i = 1;i  < numBinsY_Fluctuations+1; i++){
    Double_t numEntriesY = flucHist->ProjectionX(Form("px_%d", i), i, i)->Integral(); 
    for(Int_t j =1;  j < numBinsX_Fluctuations+1; j++){
      Double_t binContent = flucHist->GetBinContent(j,i); 
      flucHist->SetBinContent(j,i, (binContent/numEntriesY));
    }
  }

  // work on the detector hist - Reweighting done in response automatically :) 
  Int_t numBinsY_Det = detHist->GetNbinsY();
  Int_t numBinsX_Det = detHist->GetNbinsX();
  // ------------ General Strategy -----------------------                                                             
  for(Int_t i = 1; i < numBinsY_Det+1; i++){
    Double_t numEntriesY = detHist->ProjectionX(Form("px_%d", i), i, i)->Integral();
    for(Int_t j =1; j < numBinsX_Det+1; j++){
      Double_t binContent = detHist->GetBinContent(j,i);
      detHist->SetBinContent(j,i, (binContent/numEntriesY));
    }
  }
  

  // now do the matrix multiplication, first we need the rows and columns of each matrix
  // this will specify the dimensionality of our multiplication
  Int_t rowsFluc = flucHist->GetNbinsX();
  Int_t colsFluc = flucHist->GetNbinsY(); 
  Int_t rowsDet  = detHist->GetNbinsX(); 
  Int_t colsDet  = detHist->GetNbinsY(); 
  
  std::cout << "Fluctuations Response has " << rowsFluc << " rows  and " << colsFluc << " columns " << std::endl;
  std::cout << "Detector Response has " << rowsDet << " rows  and " << colsDet <<" columns " << std::endl;
  if(rowsDet != colsFluc){
    std::cout << "The dimensionality of your matrices is wrong, try again :) " << std::endl;
  }

  // we are multiplying the detector level response by the fluctuation response (order matters)
  // Multiplying matrix a and b and storing in array mult.
  for(Int_t i = 1; i < rowsFluc+1; ++i){
    for(Int_t j = 1; j < colsDet+1; ++j){
      Double_t element = 0; 
      for(Int_t k = 1; k < colsFluc+1; ++k){
	element += flucHist->GetBinContent(i,k) * detHist->GetBinContent(k,j);
      }
      
      final->SetBinContent(i,j, element);
    }
  }



  // now create the final response by doing projections of final
  RooUnfoldResponse *finalResponse = new RooUnfoldResponse(measuredHist, truthHist, final);

  // now get the full Response for comparison
  TFile *_file2 = TFile::Open("Unfolding_NeuralNetwork_Mar26th_FullStats_Part.root");
  RooUnfoldResponse* fullResp = (RooUnfoldResponse*)_file2->Get("smeared_true");
  TH2D* finalHist = (TH2D*) finalResponse->Hresponse();
  TH2D* fullHist = (TH2D*)fullResp->Hresponse();
  std::cout << fullHist->GetNbinsX() << std::endl;
  //finalHist->Divide(fullHist);
  finalHist->Draw("colz");

  TFile *fout=new TFile (Form("FullResponse_Factorized.root"),"RECREATE");
  fout->cd();
  finalResponse->Write();


} 
