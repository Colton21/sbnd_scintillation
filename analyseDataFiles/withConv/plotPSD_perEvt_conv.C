// Read in the .root event file, take the data and plot sensible histograms regarding the detected photons. 

{
  
  // read in the first (supernova) file
  const char* filename1 = "output_totalLight_conv.root";
  
  TFile * infile1=new TFile(filename1,"READ"); // define the input file
  TH1F * sn_Fp=(TH1F *)infile1->Get("sn_Fp"); 
  TH1F * rn_Fp=(TH1F *)infile1->Get("rn_Fp");

  TH1F * sn_Fp_conv=(TH1F *)infile1->Get("sn_Fp_conv"); 
  TH1F * rn_Fp_conv=(TH1F *)infile1->Get("rn_Fp_conv"); 

  
  
 //plot the histograms
    TCanvas *c1 = new TCanvas("c1","c1");
    sn_Fp->SetTitle("Fraction of prompt photons (sub 153 ns) for 10k events");
    sn_Fp->GetXaxis()->SetTitle("F_{prompt}");
    sn_Fp->GetYaxis()->SetTitle("Events");
    sn_Fp->SetMinimum(0);

    sn_Fp->GetXaxis()->SetTitleSize(0.05);
    sn_Fp->GetXaxis()->SetTitleOffset(0.85);
    sn_Fp->GetYaxis()->SetTitleSize(0.05);
    sn_Fp->GetYaxis()->SetTitleOffset(0.85);

    sn_Fp->SetLineColor(1);
    rn_Fp->SetLineColor(2);

    sn_Fp->Draw();
    rn_Fp->Draw("SAMES");



  TLegend* leg1 = new TLegend(0.2, 0.2, .8, .8);
  leg1->AddEntry("sn_Fp", "e^{-} induced scintillation");
  leg1->AddEntry("rn_Fp", "#alpha induced scintillation");
  leg1->Draw();

  
 //plot the histograms
    TCanvas *c2 = new TCanvas("c2","c2");
    sn_Fp_conv->SetTitle("Fraction of prompt photons (sub 153 ns) for 10k convolved events");
    sn_Fp_conv->GetXaxis()->SetTitle("F_{prompt}");
    sn_Fp_conv->GetYaxis()->SetTitle("Events");
    sn_Fp_conv->SetMinimum(0);

    sn_Fp_conv->GetXaxis()->SetTitleSize(0.05);
    sn_Fp_conv->GetXaxis()->SetTitleOffset(0.85);
    sn_Fp_conv->GetYaxis()->SetTitleSize(0.05);
    sn_Fp_conv->GetYaxis()->SetTitleOffset(0.85);

    sn_Fp_conv->SetLineColor(1);
    rn_Fp_conv->SetLineColor(2);

    sn_Fp_conv->Draw();
    rn_Fp_conv->Draw("SAMES");



  TLegend* leg2 = new TLegend(0.2, 0.2, .8, .8);
  leg2->AddEntry("sn_Fp_conv", "e^{-} scint (#sigma = 10 ns)");
  leg2->AddEntry("rn_Fp_conv", "#alpha scint (#sigma = 10 ns)");
  leg2->Draw();
  
  


  

}

