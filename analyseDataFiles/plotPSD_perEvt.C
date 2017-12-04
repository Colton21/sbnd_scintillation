// Read in the .root event file, take the data and plot sensible histograms regarding the detected photons. 

{
  
  // read in the first (supernova) file
  const char* filename1 = "output_totalLight.root";
  
  TFile * infile1=new TFile(filename1,"READ"); // define the input file
  TH1F * sn_Fp=(TH1F *)infile1->Get("sn_Fp"); //directs towards the graph in the root file
  TH1F * rn_Fp=(TH1F *)infile1->Get("rn_Fp"); //directs towards the graph in the root file
  
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

  TLegend* leg = new TLegend(0.2, 0.2, .8, .8);
  leg->AddEntry("sn_Fp", "e^{-} scint");
  leg->AddEntry("rn_Fp", "#alpha scint");
  leg->Draw();

  

}

