-// File to read event_file.root, take the data and plot sensible histograms regarding the Ar 39 decays (mainly the positions). 


// Important to note
  // coordinates of dimensions of SBND LAr tank (all in cm):
  // x = [ 0 -> 200 ] IN COLTONS CODE for 1 TPC. IRL active region x for 2 TPCs x_active = [-202.5 -> -2.5 ] & [2.5 -> 200] where cathode [-2.5 -> 2.5]
  // y = [ -200 -> 200 ] middle at 0
  // z = [ 0 -> 500 ]  where 0 = front.

{
  filename = "event_file.root";
  TFile * infile=new TFile(filename,"READ"); // define the input file

  TTree * anatree=(TTree *)infile->Get("event_tree"); //directs towards to root tree inside the input file

  ///////// ***** Histogram for z position of decay event ***** ///////////  

//Make a canvas
  TCanvas *c1 = new TCanvas("c1","c1");
  TH1F *z_pos = new TH1F("z_pos","Ar 39 decays wrt z axis",100,0,500);

  //Labelling
  z_pos->GetXaxis()->SetTitle("Z position [cm]");
  z_pos->GetYaxis()->SetTitle("Number of decay events");
  z_pos->SetMinimum(0);

  //Draw the momentum for every entry in the array where the pdg==11 (electron) and dump it to the h_P histogram
  event_tree->Draw("event_z_pos >> z_pos");

  //Draw the histogram
  z_pos->DrawCopy();

  /////////// ***** For y position of decay event ***** //////////

  TCanvas *c2 = new TCanvas("c2","c2");
  TH1F *y_pos = new TH1F("y_pos", "Ar 39 decays wrt y axis", 80, -200, 200);  

  //Labelling
  y_pos->GetXaxis()->SetTitle("y position [cm]");
  y_pos->GetYaxis()->SetTitle("Number of decay events");
  y_pos->SetMinimum(0);

  //Drawing
  event_tree->Draw("event_y_pos >> y_pos ");
  y_pos->DrawCopy();

  ///////// ***** For x position of decay event ***** //////////

  TCanvas *c3 = new TCanvas("c3,c3");
  TH1F *x_pos = new TH1F("x_pos", "Ar 39 decays wrt x axis", 40, 0, 200);

  //Labelling
  x_pos->GetXaxis()->SetTitle("x position [cm]");
  x_pos->GetYaxis()->SetTitle("Number of decay events");   // # makes look greeky
  x_pos->SetMinimum(0);

  event_tree->Draw("event_x_pos >> x_pos");
  x_pos->DrawCopy(); // because pointers


  ////////// ***** Plot z-y plane ***** //////////

  TCanvas *c4 = new TCanvas("c4,c4");
  TH2F *zyplane_pos = new TH2F("zyplane_pos", "Fixed energy events as viewed from ZY plane", 100, 0, 500, 80, -200, 200);  


  //Labelling
  zyplane_pos->GetXaxis()->SetTitle("z position [cm]");
  zyplane_pos->GetYaxis()->SetTitle("y position [cm]");
  //Colouring in
  zyplane_pos->SetMarkerStyle(1);
  zyplane_pos->SetLineWidth(3);
  zyplane_pos->SetLineColor(4);
  zyplane_pos->SetFillStyle(3003);
  zyplane_pos->SetFillColor(4);

  event_tree->Draw("event_y_pos:event_z_pos >> zyplane_pos"); 
  zyplane_pos->DrawCopy("COLZ"); 


  ///////// ***** Energy of event  ***** //////////

  TCanvas *c5 = new TCanvas("c5,c5");
  TH1F *E = new TH1F("E", "Energy of events",100, 5,6);

  //Labelling
  E->GetXaxis()->SetTitle("Energy of event");
  E->GetYaxis()->SetTitle("Number of events");   
  E->SetMinimum(0);

  event_tree->Draw("event_E >> E");
  E->DrawCopy(); // because pointers


 infile->Close(); // close the input file

}

