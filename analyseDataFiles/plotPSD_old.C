// Read in the .root event file, take the data and plot sensible histograms regarding the detected photons. 

{

  const char* filename1 = "10k_sn_randompos2.root";
  TFile * infile1=new TFile(filename1,"READ"); // define the input file
  TTree * data_tree1=(TTree *)infile1->Get("data_tree"); //directs towards to root tree inside the input file

  const char* filename2 = "10k_rn_randompos2.root";
  TFile * infile2=new TFile(filename2,"READ"); // define the input file
  TTree * data_tree2=(TTree *)infile2->Get("data_tree"); //directs towards to root tree inside the input file
  

  /*
Using the TF1::Fit method:
"expo" is an exponential with 2 parameters - f(x) = exp(p0 + p1*x)
   */
  /*
  TF1 *expfn = new TF1("expfn","expo(0)+expo(2)",0,5.5);
  //expfn->SetParameters(10,-25,3,-6.7);
  expfn->SetParameter(1,-160);
  expfn->SetParameter(3,-0.7);
  expfn->SetParLimits(1,-160,-10);
  expfn->SetParLimits(3,-10,0);


  //Set parameter names
  expfn->SetParName(0,"C_{1}");
  expfn->SetParName(1,"-1/#tau_{f}");
  expfn->SetParName(2,"C_{2}");
  expfn->SetParName(3,"-1/#tau_{s}");
  */



 //////////////////////////////////////////////////////////////////////////////////
 ///////// ***** Histogram for z  position of each detected photon ***** //////////   
 //////////////////////////////////////////////////////////////////////////////////

//Make a canvas
  TCanvas *c1 = new TCanvas("c1","c1");
  TH1F *sn_times = new TH1F("sn_times","Photon arrival times 10k events",2500,0,5); //500 = 10ns bins
  TH1F *rn_times = new TH1F("rn_times","Photon arrival times 10k events",2500,0,5); //500 = 10ns bins



  //Labelling
  rn_times->GetXaxis()->SetTitle("Arrival time [#mus]");
  rn_times->GetYaxis()->SetTitle("Number of photons");
  rn_times->SetMinimum(0.1);

  rn_times->GetXaxis()->SetTitleSize(0.05);
  rn_times->GetXaxis()->SetTitleOffset(0.85);
  rn_times->GetYaxis()->SetTitleSize(0.05);
  rn_times->GetYaxis()->SetTitleOffset(0.85);

  c1->SetLogy();

  //set line colors
  sn_times->SetLineColor(1);
  rn_times->SetLineColor(2);

  
  //Draw on the histrogram using the data from the tree
  data_tree1->Draw("data_time >> sn_times");
  data_tree2->Draw("data_time >> rn_times");

  


  double sn_entries = sn_times->GetEntries();
  double rn_entries = rn_times->GetEntries();
 
  double sn_integral = sn_times->Integral();
  double rn_integral = rn_times->Integral();

  cout << "SN integral: " << sn_integral << "\t";
  cout << "SN entries: " << sn_entries << "\n";

  cout << "Rn integral: " << rn_integral << "\t";
  cout << "Rn entries: " << rn_entries << "\n";
  

  /*
  double maxBin = times->GetMaximumBin();
  double maxBinValue = times->GetBinCenter(maxBin);
  double minRange = maxBinValue + 0.001; //add a nanosecond to be safe
  times->Fit("expfn","BR+","",minRange,6);
  gStyle->SetOptFit(111);
  */

  //scale
  rn_times->Scale(1./rn_integral);
  sn_times->Scale(1./sn_integral);
  
  double sn_entries_new = sn_times->GetEntries();
  double rn_entries_new = rn_times->GetEntries();

  double sn_integral_new = sn_times->Integral();
  double rn_integral_new = rn_times->Integral();
  
  cout << "SN integral new: " << sn_integral_new << "\t";
  cout << "SN entries new: " << sn_entries_new << "\n";

  cout << "Rn integral new: " << rn_integral_new << "\t";
  cout << "Rn entries new: " << rn_entries_new << "\n";



  gStyle->SetOptStat(0);
  //Draw the histogram
  rn_times->Draw();
  sn_times->Draw("SAME");
 
  
  
    TLegend* leg = new TLegend(0.2, 0.2, .8, .8);
  leg->AddEntry("sn_times", "Supernova");
  leg->AddEntry("rn_times", "Radon");
//leg->AddEntry(expfn, "exp[C_{1}-1/#tau_{f}*x]+exp[C_{2}-1/#tau_{s}*x]  ");
  leg->Draw();




  ////--------------------- PSD Analysis -----------------------------

  int rn_bins = rn_times->GetSize() - 2; // -2 two as it includes "undeflow" and "overflow"
  int sn_bins = sn_times->GetSize() - 2;

  cout << "rn_bins: " << rn_bins << endl;
  cout << "sn_bins: " << sn_bins << endl;

  int bins;
  if(rn_bins == sn_bins){ bins = rn_bins; }
  else { cerr << "ERROR: Histograms do not have the same number of bins" << endl; }


  double sum_rn = 0;
  double sum_sn = 0;
  int count = 0;
  double temp_sn, temp_rn;
  for(int i = 0; i < bins; i++) {
    temp_sn = sn_times->GetBinContent(i);
    temp_rn = rn_times->GetBinContent(i);

    sum_sn += temp_sn;
    sum_rn += temp_rn;

    sn_times->SetBinContent(i,sum_sn);
    rn_times->SetBinContent(i,sum_rn);
    count++;
  }
  
//Make a canvas
  TCanvas *c2 = new TCanvas("c2","c2");


  //Labelling
  rn_times->GetXaxis()->SetTitle("Arrival time [#mus]");
  rn_times->GetYaxis()->SetTitle("Number of photons");
  rn_times->SetMinimum(0);

  rn_times->GetXaxis()->SetTitleSize(0.05);
  rn_times->GetXaxis()->SetTitleOffset(0.85);
  rn_times->GetYaxis()->SetTitleSize(0.05);
  rn_times->GetYaxis()->SetTitleOffset(0.85);

  //set line colors
  sn_times->SetLineColor(1);
  rn_times->SetLineColor(2);

  gStyle->SetOptStat(0);
  //Draw the histogram
  rn_times->Draw();
  sn_times->Draw("SAME");
 
  
  
    TLegend* leg2 = new TLegend(0.2, 0.2, .8, .8);
  leg2->AddEntry("sn_times", "e^{-} induced scintillation");
  leg2->AddEntry("rn_times", "#alpha induced scintillation");
  leg2->Draw();


  
}

