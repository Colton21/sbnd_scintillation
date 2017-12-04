// Read in the .root event file, take the data and plot sensible histograms regarding the detected photons. 

{
  
  // read in the first (supernova) file
  const char* filename1 = "../output_1evt.root";
  TFile * infile1=new TFile(filename1,"READ"); // define the input file
  TH1F * times=(TH1F *)infile1->Get("sn_currentEvt_times"); 

  const char* filename2 = "output_1evt_conv.root";
  TFile * infile2=new TFile(filename2,"READ"); // define the input file
  TH1F * times_conv=(TH1F *)infile1->Get("sn_times_conv"); 



 //plot the histograms
    TCanvas *c1 = new TCanvas("c1","c1");
    times->SetTitle("Photon arrival times for 1 event (no convolution)");
    times->GetXaxis()->SetTitle("Time [#mus]");
    times->GetYaxis()->SetTitle("Photons");
    times->SetMinimum(0);

    times->GetXaxis()->SetTitleSize(0.05);
    times->GetXaxis()->SetTitleOffset(0.85);
    times->GetYaxis()->SetTitleSize(0.05);
    times->GetYaxis()->SetTitleOffset(0.85);

    times->SetLineColor(1);
    
    times->Draw();



 //plot the histograms
    TCanvas *c2 = new TCanvas("c2","c2");
    times_conv->SetTitle("Photon arrival times for 1 event (with convolution)");
    times_conv->GetXaxis()->SetTitle("Time [#mus]");
    times_conv->GetYaxis()->SetTitle("Photons");
    times_conv->SetMinimum(0);

    times_conv->GetXaxis()->SetTitleSize(0.05);
    times_conv->GetXaxis()->SetTitleOffset(0.85);
    times_conv->GetYaxis()->SetTitleSize(0.05);
    times_conv->GetYaxis()->SetTitleOffset(0.85);

    times_conv->SetLineColor(2);

    times_conv->Draw();

  

}

