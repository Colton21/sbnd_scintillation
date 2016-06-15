#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>  

#include "library_access.h"
#include "utility_functions.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "THStack.h"
#include "TColor.h"
#include "TLegend.h"
#include "TMarker.h"

using namespace std;

//no. SN events
//1000 in DUNE active volume (40kTon) for ~25ms time window
//SBND active volume 112 Ton --> 0.0028 of DUNE --> ~2.8 SN in SBND -->1 tpc so 1.4?
bool sort_function(std::pair<double, int> pair1, std::pair<double, int> pair2)
  { return (pair1.first < pair2.first); }


int main()
{
  gRandom->SetSeed(0);


  const double MassE = 0.510998910; // mass electron - MeV/c^2
  const double Q = 0.565;//Q value of decay - Ar39
  const double t_singlet = 0.000000006; //6ns
  const double t_triplet = 0.0000015; //1.5 us
  const double scint_time_window = 0.00001; //10 us
  //const double scint_time_window = 0.000005; //5 us
  const double Eav = 20.;
  const int scint_yield = 24000;
  const double quantum_efficiency = 0.2;//Expected value
  //const double quantum_efficiency = 0.4;//Test value
  const double activity = 1.; //Ar39 roughly 1 Bq/kg
  const double volume = 112000.; //SBND 112ton LAr
  const double time_window = 0.0012 * 50.; //[s] -- 25ms is the timing window for SN arrival, 1.2ms is measurement window
  const double expected_sn = 2.8;
  const int cut_number = 6;
  
  cout << "Cut: " << cut_number << endl;

  vector<double> timing;
  vector< std::pair< double, int > > total_timing;
  signed long int big_number = 100000000;
  total_timing.reserve(big_number);

  //const int max_events = 10000;
  //const int max_events_sn = 1;
  //int max_events_sn = utility::poisson(expected_sn,gRandom->Uniform(1.),1.);
  const int max_events_sn = 10000;
  const int max_events = activity * volume/2 * time_window;
  //const int max_events = 20;
  cout << "Ar39 events: " << max_events << " in: " << time_window << " s" << endl;


  double temp_energy;

  TF1 *fSpectrum = new TF1("fSpectrum",utility::SpectrumFunction,0,Q,1);//-----Beta decay spectrum
  TF1 *flandau = new TF1("flandau",utility::fsn, 0, 50, 1);//--SN Nu spectrum
  fSpectrum->SetParameter(0, Q);
  flandau->SetParameter(0, Eav);
  TF1 *fScintillation_function = new TF1("Ar Scintillation Timing", utility::Scintillation_function, 0, scint_time_window, 2);
  fScintillation_function->SetParameter(0, t_singlet);
  fScintillation_function->SetParameter(1, t_triplet);

  cout << '\n' << "Poisson fluctuated SN " << max_events_sn << endl;

  LibraryAccess lar_light;
/*
  TCanvas *scint_can = new TCanvas();
  scint_can->cd();
  scint_can->SetLogy();
  fScintillation_function->SetNpx(500);
  fScintillation_function->Draw();
  fScintillation_function->GetXaxis()->SetTitle("Time [s]");
  fScintillation_function->GetYaxis()->SetTitle("Probability");
  scint_can->Print("ar40_scintillation_timing.pdf");

  TCanvas *test1 = new TCanvas("c1","c1",500,600);
  TCanvas *test2 = new TCanvas("c2","c2",500,600);

  test1->cd();
  fSpectrum->GetXaxis()->SetTitle("Energy [MeV]");
  fSpectrum->SetTitle("Ar39 Beta Spectrum");
  fSpectrum->Draw();

  test2->cd();
  flandau->GetXaxis()->SetTitle("Energy [MeV]");
  flandau->SetTitle("Super Nova Neutrino Spectrum");
  flandau->Draw();

  test1->Print("ar39_spectrum.pdf");
  test2->Print("sn_spectrum.pdf");
*/
 //*********** loading the library ****************
  
  //std::string libraryfile = "/home/colton/PhD/Year1/Practice/20160111/opLibArrayPMTswCathd.root";
  std::string libraryfile = "/home/colton/PhD/Year1/Practice/20160418/timing_files/Lib154PMTs8inch_FullFoilsTPB.root";
  //std::string libraryfile = "/home/colton/PhD/Year1/Practice/20160418/timing_files/Lib154PMTs8inch_NoCathodeNoFoils.root";
  //std::string libraryfile = "/home/colton/PhD/Year1/Practice/20160418/timing_files/Lib154PMTs8inch_OnlyCathodeTPB.root";
  //fNVoxels = int(gxSteps*gySteps*gzSteps);
  //fNOpChannels = 308;
  bool reflected = true;
  bool reflT = true;
  lar_light.LoadLibraryFromFile(libraryfile, reflected, reflT);

 //*********** end loading the library *************

 // Reading out positions of PMT from txt file *****
  ifstream myfile;
  myfile.open("posPMTs_setup1.txt");
  if(myfile.is_open()){cout << "File opened successfully" << endl;}
  vector<vector<double>> myfile_data;
  while(!myfile.eof()){
    double num_pmt, x_pmt, y_pmt, z_pmt;
    if(myfile >> num_pmt >> x_pmt >> y_pmt >> z_pmt){
      vector<double> line_data({num_pmt, x_pmt, y_pmt, z_pmt});
      myfile_data.push_back(line_data);
    }
    else{break;}
  }
 // End Reading out positions of PMT from txt file


//Declare Vectors of TMarkers
     vector<TMarker*> candidateMarker_ar;
     vector<TMarker*> candidateMarker_sn;
     vector<TMarker*> candidateMarker_higherE;
     vector<TMarker*> candidateMarker_rd;
     const int pmtMarkerstyle = 24;


 //------Histogram Declarations--------------

  //Histogram for how many pmts see signal from Ar39 decay
  //Each stacked histogram contains both ar39 and sn results for a single cut

  TH1D *h_no_pmt = new TH1D("h_no_pmt", "h_no_pmt" , 42, 0, 42);
  TH1D *h_no_pmt_sn = new TH1D("h_no_pmt_sn", "h_no_pmt_sn" , 42, 0, 42);
  THStack *h_no_pmt_stack = new THStack("h_no_pmt_stack", "Rn and SN Events on PMTs");

  TH1D *h_no_pmt3 = new TH1D("h_no_pmt3", "h_no_pmt3" , 42, 0, 42);
  TH1D *h_no_pmt_sn3 = new TH1D("h_no_pmt_sn3", "h_no_pmt_sn3" , 42, 0, 42);
  THStack *h_no_pmt_stack3 = new THStack("h_no_pmt_stack3", "Rn and SN Events on PMTs");

  TH1D *h_no_pmt10 = new TH1D("h_no_pmt10", "h_no_pmt10" , 42, 0, 42);
  TH1D *h_no_pmt_sn10 = new TH1D("h_no_pmt_sn10", "h_no_pmt_sn10" , 42, 0, 42);
  THStack *h_no_pmt_stack10 = new THStack("h_no_pmt_stack10", "Rn and SN Events on PMTs");

  TH1D *h_no_pmt103 = new TH1D("h_no_pmt103", "h_no_pmt103" , 42, 0, 42);
  TH1D *h_no_pmt_sn103 = new TH1D("h_no_pmt_sn103", "h_no_pmt_sn103" , 42, 0, 42);
  THStack *h_no_pmt_stack103 = new THStack("h_no_pmt_stack10", "Rn and SN Events on PMTs");
//--------------------------------------------------------------------------------

  TH2D *h_no_pmt_distance = new TH2D("h_no_pmt_distance", "h_no_pmt_distance", 50, 0, 50, 20, 0, 200);
  TH2D *h_no_pmt_distance_cut1 = new TH2D("h_no_pmt_distance1", "h_no_pmt_distance1", 50, 0, 50, 20, 0, 200);
  TH2D *h_no_pmt_distance_cut2 = new TH2D("h_no_pmt_distance2", "h_no_pmt_distance2", 50, 0, 50, 20, 0, 200);
  TH2D *h_no_pmt_distance_cut3 = new TH2D("h_no_pmt_distance3", "h_no_pmt_distance3", 50, 0, 50, 20, 0, 200);

  TH2D *h_no_pmt_distance_sn = new TH2D("h_no_pmt_distance_sn", "h_no_pmt_distance_sn", 50, 0, 50, 20, 0, 200);
  TH2D *h_no_pmt_distance_cut1_sn = new TH2D("h_no_pmt_distance1_sn", "h_no_pmt_distance1_sn", 50, 0, 50, 20, 0, 200);
  TH2D *h_no_pmt_distance_cut2_sn = new TH2D("h_no_pmt_distance2_sn", "h_no_pmt_distance2_sn", 50, 0, 50, 20, 0, 200);
  TH2D *h_no_pmt_distance_cut3_sn = new TH2D("h_no_pmt_distance3_sn", "h_no_pmt_distance3_sn", 50, 0, 50, 20, 0, 200);

//--------------------------------------------------------------------------------

//Histograms for the 2d representation of the pmts ------------------
//binning appropriate for dimensions and attempting to get symmetrical
  TH2D *h_pmt_plane0 = new TH2D("h_pmt_plane0", "h_pmt_plane0", 14 ,0, 500, 11, -200, 200);
  TH2D *h_pmt_plane1 = new TH2D("h_pmt_plane1", "h_pmt_plane1", 14 ,0, 500, 11, -200, 200);
  TH2D *h_pmt_plane2 = new TH2D("h_pmt_plane2", "h_pmt_plane2", 14 ,0, 500, 11, -200, 200);

  TH2D *h_pmt_plane02 = new TH2D("h_pmt_plane02", "h_pmt_plane02", 14 ,0, 500, 11, -200, 200);
  TH2D *h_pmt_plane12 = new TH2D("h_pmt_plane12", "h_pmt_plane12", 14 ,0, 500, 11, -200, 200);
  TH2D *h_pmt_plane22 = new TH2D("h_pmt_plane22", "h_pmt_plane22", 14 ,0, 500, 11, -200, 200);

  TH2D *h_pmt_plane_rd0 = new TH2D("h_pmt_plane_rd0", "h_pmt_plane_rd0", 14 ,0, 500, 11, -200, 200);
  TH2D *h_pmt_plane_rd1 = new TH2D("h_pmt_plane_rd1", "h_pmt_plane_rd1", 14 ,0, 500, 11, -200, 200);
  TH2D *h_pmt_plane_rd2 = new TH2D("h_pmt_plane_rd2", "h_pmt_plane_rd2", 14 ,0, 500, 11, -200, 200);

  TH2D *h_pmt_plane_ar0 = new TH2D("h_pmt_plane_ar0", "h_pmt_plane0", 14 ,0, 500, 11, -200, 200);
  TH2D *h_pmt_plane_ar1 = new TH2D("h_pmt_plane_ar1", "h_pmt_plane1", 14 ,0, 500, 11, -200, 200);
  TH2D *h_pmt_plane_ar2 = new TH2D("h_pmt_plane_ar2", "h_pmt_plane2", 14 ,0, 500, 11, -200, 200);

/*  TH2D *h_pmt_plane0 = new TH2D("h_pmt_plane0", "h_pmt_plane0", 12 ,33.3333, 466.667, 9, -166.667, 166.667);
  TH2D *h_pmt_plane1 = new TH2D("h_pmt_plane1", "h_pmt_plane1", 12 ,33.3333, 466.667, 9, -166.667, 166.667);
  TH2D *h_pmt_plane2 = new TH2D("h_pmt_plane2", "h_pmt_plane2", 12 ,33.3333, 466.667, 9, -166.667, 166.667);
*/

//-------------------------------------------------------------------


  //Simple histograms for the change in energy spectrum
  //as cuts are applied.
  //NOTE: The number of bins is important for counting later
  const int num_energy_bins = 100;
  TH1D *h_energy = new TH1D("h_energy", "", num_energy_bins, 0, 0.57);
  TH1D *h_energy_sn = new TH1D("h_energy_sn", "", num_energy_bins, 0, 50);
  TH1D *h_energy1 = new TH1D("h_energy1", "", num_energy_bins, 0, 0.57);
  TH1D *h_energy_sn1 = new TH1D("h_energy_sn1", "", num_energy_bins, 0, 50);
  TH1D *h_energy2 = new TH1D("h_energy2", "", num_energy_bins, 0, 0.57);
  TH1D *h_energy_sn2 = new TH1D("h_energy_sn2", "", num_energy_bins, 0, 50);
  TH1D *h_energy22 = new TH1D("h_energy22", "", num_energy_bins, 0, 0.57);
  TH1D *h_energy_sn22 = new TH1D("h_energy_sn22","",num_energy_bins, 0, 50);

//--------------------------------------------------------------------------

  TH1D *h_pe_time = new TH1D("h_pe_time", "", 30, 0, 30);
  TH1D *h_single = new TH1D("h_single", "", 12000, 0, 1200);


  //---------End histogram code -------------
//-----------------Ar39 Loop-----------------------------
  vector<double> energy_ar_list;
  energy_ar_list.reserve(max_events);
  vector<double> scint_time_list;
  scint_time_list.reserve(max_events);
  vector<int> voxel_list;
  voxel_list.reserve(max_events);

  for(int events = 0; events < max_events; events++)
  {
    energy_ar_list.push_back(fSpectrum->GetRandom());
    //energy_ar_list.push_back(0.565);
    double decay_time = time_window * gRandom->Uniform(1.);
    scint_time_list.push_back(decay_time);
    int rand_voxel = gRandom->Uniform(319999);
    voxel_list.push_back(rand_voxel);//select random voxel -> random position in detector

  }

  int realisticPMT_IDs[60] = {0, 4, 8, 12, 16, 20, 24, 32, 40, 44, 48, 52, 56, 60, 64, 88, 92, 96, 100, 104, 108, 112, 120, 128, 132, 136, 140, 144, 148, 152, 154, 158, 162, 166, 170, 174, 178, 186, 194, 198, 202, 206, 210, 214, 218, 242, 246, 250, 254, 258, 262, 266, 274, 282, 286, 290, 294, 298, 302, 306};

/*
  for(int pmt_loop = 0; pmt_loop < 60; pmt_loop++){
     int num_pmt = realisticPMT_IDs[pmt_loop];
     cout << "------" << num_pmt << endl;
     int events_counter = 0;

     double x_pmt = myfile_data.at(num_pmt).at(1);
     double y_pmt = myfile_data.at(num_pmt).at(2);
     double z_pmt = myfile_data.at(num_pmt).at(3);

     TVector3 optdet (x_pmt/100., y_pmt/100., z_pmt/100.);

     for(int events = 0; events < max_events; events++){

        vector<double> pmt_hits = lar_light.PhotonLibraryAnalyzer(energy_ar_list.at(events), scint_yield, quantum_efficiency, num_pmt, voxel_list.at(events));

	int num_VUV = pmt_hits.at(0);
	int num_VIS = pmt_hits.at(1);
*/
/*
	h_pmt_plane_ar0->Fill(z_pmt, y_pmt, pmt_hits.at(0)+pmt_hits.at(1));/// Total
	h_pmt_plane_ar1->Fill(z_pmt, y_pmt, pmt_hits.at(0));/// VUV
	h_pmt_plane_ar2->Fill(z_pmt, y_pmt, pmt_hits.at(1));/// VIS
*/
/*

	if(num_VUV+num_VIS == 0){continue;}

	double distance_to_pmt = 
	  std::sqrt((pmt_hits.at(2)-x_pmt)*(pmt_hits.at(2)-x_pmt) + 
	  (pmt_hits.at(3)-y_pmt)*(pmt_hits.at(3)-y_pmt) + 
	  (pmt_hits.at(4)-z_pmt)*(pmt_hits.at(4)-z_pmt));
           
		  
	///*************************** ///
	///***** TIMING OPERATIONS *** ///
	///*************************** ///

	vector<double> time_vuv;
	if(num_VUV != 0){ time_vuv = utility::GetVUVTime(distance_to_pmt, num_VUV);}

	if(num_VUV != time_vuv.size())
	{
	  cout << "Param fail" << endl;
	  num_VUV = time_vuv.size();
	}

	for(auto& x: time_vuv)
	{
	  total_timing.push_back(make_pair(x * 0.001 + (scint_time_list.at(events) + fScintillation_function->GetRandom())*1000000., num_pmt));
	}
	time_vuv.clear();
	  
	  
*/
	 /* // this is for foils only on cathode
 	  vector<double> time_vis
	  if(num_VIS != 0){
	    time_vis = utility::GetVisibleTimeOnlyCathode(pmt_hits.at(5), num_VIS);
	    for(auto &y : time_vis)
	    {
	      total_timing.push_back(make_pair(y*0.001+(scint_time_list.at(events) + fScintillation_function->GetRandom())*1000000., num_pmt));
	    }
	    time_vis.clear();
	  }
	  */

	  /// This function expects distances in m
	  /// For full foil configuration only
/*
	  if(num_VIS != 0)
	  {
	    vector<double> time_vis;
	    TVector3 scint_point (pmt_hits.at(2)/100., pmt_hits.at(3)/100., pmt_hits.at(4)/100.);
	    double tmean = utility::TimingParamReflected(scint_point, optdet);

	    //time_vis.reserve(num_VIS);
	    time_vis = utility::GetVisibleTimeFullConfig(pmt_hits.at(5), tmean, distance_to_pmt, num_VIS);

	  /// Loop over number of photons generated
	    for(auto &y : time_vis)
	    {
	      total_timing.push_back(make_pair(y*0.001+(scint_time_list.at(events) + fScintillation_function->GetRandom())*1000000., num_pmt));
	    }
  	    time_vis.clear();
	  }

    }//end looping events

  if(total_timing.size() >= 10000)
  {
    cout << total_timing.size() << endl;
  }

  }//end loop over pmts

  ///For full foils configuration!
  int total_extra_vis = 292416211 * time_window;
  for(int extra_vis = 0; extra_vis < total_extra_vis; extra_vis++)
  {
    double rand_time = gRandom->Uniform(time_window)*1000000.;
    int rand_num = gRandom->Uniform(60);
    int rand_pmt = realisticPMT_IDs[rand_num];
    total_timing.push_back(make_pair(rand_time, rand_pmt));
  }
*/
/*
  ///For cathode only foils!
  int total_extra_vis = 143210695 * time_window;
  for(int extra_vis = 0; extra_vis < total_extra_vis; extra_vis++)
  {
    double rand_time = gRandom->Uniform(time_window)*1000000.;
    int rand_num = gRandom->Uniform(60);
    int rand_pmt = realisticPMT_IDs[rand_num];
    total_timing.push_back(make_pair(rand_time, rand_pmt));
  }
*/
/*
  ///For no foils! - VUV
  int total_extra_vuv = 80003251 * time_window;
  for(int extra_vuv = 0; extra_vuv < total_extra_vuv; extra_vuv++)
  {
    double rand_time = gRandom->Uniform(time_window)*1000000.;
    int rand_num = gRandom->Uniform(60);
    int rand_pmt = realisticPMT_IDs[rand_num];
    total_timing.push_back(make_pair(rand_time, rand_pmt));
  }

  std::sort(total_timing.begin(), total_timing.end(), sort_function);

  energy_ar_list.clear();
  scint_time_list.clear();
  voxel_list.clear();


  int counter = 0; /// Unique PMTs
  int prev_pmt1;
  int prev_pmt2;
  int flag = 0;
  int trigger_flag = 0;
  int new_gamma1;
  bool triggered = false;
  vector<int> pmt_list;


  if(total_timing.size() != 0){
    for(int gamma1 = 0; gamma1 < total_timing.size()-1; gamma1++)
    {
      ///If previous gamma1 saw a trigger, move up 100ns
      if(triggered == true){gamma1 = new_gamma1;}
      triggered = false;

      double t01 = total_timing.at(gamma1).first;
      for(int gamma2 = gamma1+1; gamma2 < total_timing.size(); gamma2++)
      {
        double t02 = total_timing.at(gamma2).first;
        if(total_timing.at(gamma1).second != total_timing.at(gamma2).second)
        {
  	  if(t02 - t01 <= 0.01)
	  {
	    if(counter == 0)
	    {
	      counter++;
	      pmt_list.push_back(gamma1);
	      pmt_list.push_back(gamma2);
	    }
	    if(counter == 1 && total_timing.at(gamma2).second != prev_pmt1 && total_timing.at(gamma2).second != prev_pmt2)
	    {
	      counter++;
	      pmt_list.push_back(gamma2);
	    }/// If another unique PMT
	  }/// if within time

	  else
	  {
	    /// Number of PMT combinations (1,2,3)
	    if(counter >= 2)
	    {
	      //cout << counter << endl;
	      for(int pmt1 = 0; pmt1 < pmt_list.size(); pmt1++)
	      {
	        int evt = pmt_list.at(pmt1); /// number in total_timing
	        int same_pmt_counter = 0;
	        for(int element = evt+1; element < total_timing.size(); element++)
	        {
		  //events within 100 ns
		  if(total_timing.at(element).first - total_timing.at(evt).first <= 0.1)
		  {
		    if(total_timing.at(element).second == total_timing.at(evt).second)/// now see if two events are on the same PMT
		    {
		      same_pmt_counter++;
		    }
		  }/// End loop if events match within 100ns
		  else
		  {
		    ///*************************************
		    ///***** Number of PE on single PMT ****
		    ///*************************************
		    ///***** THIS IS THE CUT - HERE! *******
		    if(same_pmt_counter >= cut_number)
		    {
		      flag++;
		      //cout << same_pmt_counter << endl;
		    }
		    break;
		  }
	        }/// End loop all events again
	      }/// End loop of 3 PMTs which saw signal within 10ns
	    }/// if 3 unique PMTs
	    counter = 0;
	    pmt_list.clear();
	    break;
	  }/// else

  	  prev_pmt1 = total_timing.at(gamma1).second;
	  prev_pmt2 = total_timing.at(gamma2).second;

        }//pm1!=pmt2
      }//gamma2

      /// If all 3 PMTs satisfy the X PE timing
      if(flag >= 3)
      {
        trigger_flag++;
	  /// Next trigger can only happen +100ns
        for(new_gamma1 = gamma1; new_gamma1 < total_timing.size()-1; new_gamma1++)
        {
	  if(total_timing.at(gamma1).first + 0.1 >= time_window * 1000000.){break;}
	  if(total_timing.at(new_gamma1).first >= total_timing.at(gamma1).first + 0.1)
	  {
	    triggered = true;
	    break;
	  } 
        }
      }
      flag = 0;
    }//gamma1
    cout << trigger_flag << endl;
  }// end if total_timing.size() ! = 0

  total_timing.clear();

  TCanvas *can_time = new TCanvas();
  can_time->cd();
  h_pe_time->GetXaxis()->SetTitle("Number of PE within 100ns");
  h_pe_time->Draw();
  can_time->SetLogy();
  can_time->Print("pe_time.pdf");

  TCanvas *can_time2 = new TCanvas();
  can_time2->cd();
  h_single->GetYaxis()->SetTitle("Number of Photoelectrons / 60 PMTs");
  h_single->GetXaxis()->SetTitle("Time [#mu s]");
  h_single->Draw();
  can_time2->Print("pmt_time.pdf");

*/
 //------End Ar39 loop------------------
 //--------------------------------------------------- 

 //-----Supernova loop-------------------
 // Loop for number of times generate random voxel, with decay, considering visibility, energy
     int counter_sn;

     for(int events_sn = 0; events_sn < max_events_sn; events_sn++){

	//if(events_sn % 500 == 0 ){cout << events_sn << endl; }

       double scint_time_list = time_window * gRandom->Uniform(1.);

       vector<double> pmt_hits_sn;
       vector<double> pmt_hits_sn2;
       vector<double> pmt_hits_rd;

	//const double energy_sn = 5.;
	double energy_sn = flandau->GetRandom();
	//We only need the energy_sn2 if we're creating the colour normalised plots with 1 ar39 event and 2 sn events
	//Otherwise this should be commmented out, along with terms using energy_sn2
	const double energy_sn2 = 22.;
	const double energy_rd = 6.9;
        int rand_voxel = gRandom->Uniform(319999);//select random voxel -> random position in detector
	//const double energy_sn = 6.906; //actually the Po216 alpha decay
	//const double energy_sn = 8.955;//Po212 decay
	//cout << "SN LOOP" << endl;
	//cout << "SN Energy: " << energy_sn << endl;
	//cout << "SN Energy2: " << energy_sn2 << endl;

	//Counters for counting number of times a PMT sees 1+ events
	int counter_sn3 = 0;
	int counter_sn10 = 0;
	int counter_sn103 = 0;
	int pmt_counter = 0;

        int realisticPMT_IDs[60] = {0, 4, 8, 12, 16, 20, 24, 32, 40, 44, 48, 52, 56, 60, 64, 88, 92, 96, 100, 104, 108, 112, 120, 128, 132, 136, 140, 144, 148, 152, 154, 158, 162, 166, 170, 174, 178, 186, 194, 198, 202, 206, 210, 214, 218, 242, 246, 250, 254, 258, 262, 266, 274, 282, 286, 290, 294, 298, 302, 306};

	//Loop over all the PMTs
	 for(int loop_pmt = 0; loop_pmt < 60; loop_pmt++){

	  int num_pmt = realisticPMT_IDs[loop_pmt];

          pmt_hits_sn = lar_light.PhotonLibraryAnalyzer(energy_sn, scint_yield, quantum_efficiency, num_pmt, rand_voxel);
          pmt_hits_sn2 = lar_light.PhotonLibraryAnalyzer(energy_sn2, scint_yield, quantum_efficiency, num_pmt, rand_voxel);
          pmt_hits_rd = lar_light.PhotonLibraryAnalyzer(energy_rd, scint_yield, quantum_efficiency, num_pmt, rand_voxel);

          double x_pmt = myfile_data.at(num_pmt).at(1);
          double y_pmt = myfile_data.at(num_pmt).at(2);
          double z_pmt = myfile_data.at(num_pmt).at(3);

          TVector3 optdet (x_pmt/100., y_pmt/100., z_pmt/100.);

	  int num_VUV = pmt_hits_sn.at(0);
	  int num_VIS = pmt_hits_sn.at(1);

	  if(num_VUV+num_VIS == 0){continue;}

	  double distance_to_pmt = 
	    std::sqrt((pmt_hits_sn.at(2)-x_pmt)*(pmt_hits_sn.at(2)-x_pmt) + 
	    (pmt_hits_sn.at(3)-y_pmt)*(pmt_hits_sn.at(3)-y_pmt) + 
	    (pmt_hits_sn.at(4)-z_pmt)*(pmt_hits_sn.at(4)-z_pmt));
           
		  
	///*************************** ///
	///***** TIMING OPERATIONS *** ///
	///*************************** ///

	vector<double> time_vuv;
	if(num_VUV != 0){ time_vuv = utility::GetVUVTime(distance_to_pmt, num_VUV);}

	if(num_VUV != time_vuv.size())
	{
	  cout << "Param fail" << endl;
	  num_VUV = time_vuv.size();
	}

	for(auto& x: time_vuv)
	{
	  total_timing.push_back(make_pair(x * 0.001 + (scint_time_list + fScintillation_function->GetRandom())*1000000., num_pmt));
	}
	time_vuv.clear();
	  
	  

	 /* // this is for foils only on cathode
 	  vector<double> time_vis
	  if(num_VIS != 0){
	    time_vis = utility::GetVisibleTimeOnlyCathode(pmt_hits.at(5), num_VIS);
	    for(auto &y : time_vis)
	    {
	      total_timing.push_back(make_pair(y*0.001+(scint_time_list + fScintillation_function->GetRandom())*1000000., num_pmt));
	    }
	    time_vis.clear();
	  }
	  */

	  /// This function expects distances in m
	  /// For full foil configuration only

	  if(num_VIS != 0)
	  {
	    vector<double> time_vis;
	    TVector3 scint_point (pmt_hits_sn.at(2)/100., pmt_hits_sn.at(3)/100., pmt_hits_sn.at(4)/100.);
	    double tmean = utility::TimingParamReflected(scint_point, optdet);

	    //time_vis.reserve(num_VIS);
	    time_vis = utility::GetVisibleTimeFullConfig(pmt_hits_sn.at(5), tmean, distance_to_pmt, num_VIS);

	  /// Loop over number of photons generated
	    for(auto &y : time_vis)
	    {
	      total_timing.push_back(make_pair(y*0.001+(scint_time_list + fScintillation_function->GetRandom())*1000000., num_pmt));
	    }
  	    time_vis.clear();
	  }

         }//End looping pmts

         std::sort(total_timing.begin(), total_timing.end(), sort_function);

         int counter = 0; /// Unique PMTs
         int prev_pmt1;
         int prev_pmt2;
         int flag = 0;
         int trigger_flag = 0;
         int new_gamma1;
         bool triggered = false;
         vector<int> pmt_list;


         if(total_timing.size() != 0){
           for(int gamma1 = 0; gamma1 < total_timing.size()-1; gamma1++)
           {
             ///If previous gamma1 saw a trigger, move up 100ns
             if(triggered == true){gamma1 = new_gamma1;}
             triggered = false;

             double t01 = total_timing.at(gamma1).first;
             for(int gamma2 = gamma1+1; gamma2 < total_timing.size(); gamma2++)
             {
               double t02 = total_timing.at(gamma2).first;
               if(total_timing.at(gamma1).second != total_timing.at(gamma2).second)
               {
  	         if(t02 - t01 <= 0.01)
	         {
	           if(counter == 0)
	           {
	             counter++;
	             pmt_list.push_back(gamma1);
	             pmt_list.push_back(gamma2);
	           }
	           if(counter == 1 && total_timing.at(gamma2).second != prev_pmt1 && total_timing.at(gamma2).second != prev_pmt2)
	           {
	             counter++;
	             pmt_list.push_back(gamma2);
	           }/// If another unique PMT
	         }/// if within time

	         else
	         {
	           /// Number of PMT combinations (1,2,3)
	           if(counter >= 2)
	           {
	             //cout << counter << endl;
	             for(int pmt1 = 0; pmt1 < pmt_list.size(); pmt1++)
	             {
	               int evt = pmt_list.at(pmt1); /// number in total_timing
	               int same_pmt_counter = 0;
	               for(int element = evt+1; element < total_timing.size(); element++)
	               {
		         //events within 100 ns
		         if(total_timing.at(element).first - total_timing.at(evt).first <= 0.1)
		         {
		           if(total_timing.at(element).second == total_timing.at(evt).second)/// now see if two events are on the same PMT
		           {
		             same_pmt_counter++;
		           }
		         }/// End loop if events match within 100ns
		         else
		         {
		           ///*************************************
		           ///***** Number of PE on single PMT ****
		           ///*************************************
		           ///***** THIS IS THE CUT - HERE! *******
		           if(same_pmt_counter >= cut_number)
		           {
		             flag++;
		           
		           }
		           break;
       		         }
	               }/// End loop all events again
	             }/// End loop of 3 PMTs which saw signal within 10ns
	           }/// if 3 unique PMTs
	           counter = 0;
	           pmt_list.clear();
	           break;
	         }/// else

  	         prev_pmt1 = total_timing.at(gamma1).second;
	         prev_pmt2 = total_timing.at(gamma2).second;

               }//pm1!=pmt2
             }//gamma2

             /// If all 3 PMTs satisfy the X PE timing
             if(flag >= 3)
             {
               trigger_flag++;
	         /// Next trigger can only happen +100ns
               for(new_gamma1 = gamma1; new_gamma1 < total_timing.size()-1; new_gamma1++)
               {
	         if(total_timing.at(gamma1).first + 0.1 >= time_window * 1000000.){break;}
	         if(total_timing.at(new_gamma1).first >= total_timing.at(gamma1).first + 0.1)
	         {
	           triggered = true;
	           break;
	         } 
               }
             }
             flag = 0;
           }//gamma1
           //cout << trigger_flag << endl;
         }// end if total_timing.size() ! = 0


	total_timing.clear();

	if(trigger_flag != 0){counter_sn++;}


       //counter_sn = 0;
       counter_sn3 = 0;
       counter_sn10 = 0;
       counter_sn103 = 0;

       pmt_hits_sn.clear();
       pmt_hits_sn2.clear();
       pmt_hits_rd.clear();
    }/// End loop SN

    cout << "SN Pass: " << counter_sn << endl;

    TCanvas *can_marker_ar = new TCanvas();
    TCanvas *can_marker_sn = new TCanvas();
    TCanvas *can_marker_hE = new TCanvas();
    TCanvas *can_marker_rd = new TCanvas();

    can_marker_ar->cd();

    double z[2] = {0, 500};//{-50, 550};
    double y[2] = {-200, 200};
    TGraph *g0 = new TGraph(2,z,y);
    TGraph *g1 = new TGraph(2,z,y);
    TGraph *g2 = new TGraph(2,z,y);
    TGraph *g3 = new TGraph(2,z,y);

    g0->SetTitle("Ar39 0.5 MeV"); 
    g0->GetYaxis()->SetTitle("y [cm]");
    g0->GetXaxis()->SetTitle("z [cm]");
    g0->Draw("ap");
    for(int i = 0; i < candidateMarker_ar.size(); i++){candidateMarker_ar.at(i)->Draw("same");}
    can_marker_ar->Update();

    can_marker_sn->cd();
    g1->SetTitle("Supernova 5 MeV"); 
    g1->GetYaxis()->SetTitle("y [cm]");
    g1->GetXaxis()->SetTitle("z [cm]");
    g1->Draw("ap");
    for(int i = 0; i < candidateMarker_sn.size(); i++){candidateMarker_sn.at(i)->Draw("same");}
    can_marker_sn->Update();

    can_marker_hE->cd();
    g2->SetTitle("Supernova 22 MeV"); 
    g2->GetYaxis()->SetTitle("y [cm]");
    g2->GetXaxis()->SetTitle("z [cm]");
    g2->Draw("ap");
    for(int i = 0; i < candidateMarker_higherE.size(); i++){candidateMarker_higherE.at(i)->Draw("same");}
    can_marker_hE->Update();

    can_marker_rd->cd();
    g3->SetTitle("Radon 6.9 MeV");
    g3->GetYaxis()->SetTitle("y [cm]");
    g3->GetXaxis()->SetTitle("z [cm]");
    g3->Draw("ap");
    for(int i = 0; i < candidateMarker_rd.size(); i++){candidateMarker_rd.at(i)->Draw("same");}
    can_marker_rd->Update();

    candidateMarker_ar.clear();
    candidateMarker_sn.clear();
    candidateMarker_higherE.clear();
    candidateMarker_rd.clear();
//-----------End Supernova loop--------------------------------------------
//--------------------------------------------------------------------------
//------Modifying colour scale ---------------------------


//Normalising colour scale accross VUV,VIS,TOTAL Ar39 signal

  Int_t ci; //for color index setting
  TColor *colorz;//for color definition with alpha
  ci = TColor::GetColor("#000099");
//---
//Normalising colur scale accross totals: ar39, sn1, sn2, and rd
//uniform across totals
  double temp6 = h_pmt_plane1->GetBinContent(h_pmt_plane1->GetMaximumBin());
  double temp7 = h_pmt_plane_ar1->GetBinContent(h_pmt_plane_ar1->GetMaximumBin());
  double temp8 = h_pmt_plane12->GetBinContent(h_pmt_plane12->GetMaximumBin());
  double temp9 = h_pmt_plane_rd1->GetBinContent(h_pmt_plane_rd1->GetMaximumBin());
  double max2; 
  if(temp6 > temp7 && temp6 > temp8 && temp6 > temp9){max2 = temp6;}
  if(temp7 > temp6 && temp7 > temp8 && temp7 > temp9){max2 = temp7;}
  if(temp8 > temp6 && temp8 > temp7 && temp8 > temp9){max2 = temp8;}
  if(temp9 > temp8 && temp9 > temp7 && temp9 > temp6){max2 = temp9;}

  temp6 = 0; temp7 = 0; temp8 = 0; temp9 = 0;

  temp6 = h_pmt_plane1->GetBinContent(h_pmt_plane1->GetMinimumBin());
  temp7 = h_pmt_plane_ar1->GetBinContent(h_pmt_plane_ar1->GetMinimumBin());
  temp8 = h_pmt_plane12->GetBinContent(h_pmt_plane12->GetMinimumBin());
  temp9 = h_pmt_plane_rd1->GetBinContent(h_pmt_plane_rd1->GetMinimumBin());
  double min2; 
  if(temp6 < temp7 && temp6 < temp8 && temp6 < temp9){min2 = temp6;}
  if(temp7 < temp6 && temp7 < temp8 && temp7 < temp9){min2 = temp7;}
  if(temp8 < temp6 && temp8 < temp7 && temp8 < temp9){min2 = temp8;}
  if(temp9 < temp6 && temp9 < temp7 && temp9 < temp8){min2 = temp9;}

  double min_x = 0.95*min2;
  double max_x = 1.15*max2;

  const Int_t pLevels = 99;
  Double_t __levels[pLevels];

  for(Int_t i = 1; i < pLevels; i++){

    __levels[i] = min_x + (max_x-min_x)/(pLevels-1)*(i);

  }

  __levels[0]=0.0001;
  h_pmt_plane1->SetContour((sizeof(__levels)/sizeof(Double_t)),__levels);
  h_pmt_plane_ar1->SetContour((sizeof(__levels)/sizeof(Double_t)),__levels);
  h_pmt_plane12->SetContour((sizeof(__levels)/sizeof(Double_t)),__levels);
  h_pmt_plane_rd1->SetContour((sizeof(__levels)/sizeof(Double_t)),__levels);

  gStyle->SetPalette(kBird);

  Int_t cj; //for color index setting
  TColor *colorz2;//for color definition with alpha
  cj = TColor::GetColor("#000099");

//--------------------------------------------------------
//****2D Histograms of events "displayed"****************

  TCanvas *can60 = new TCanvas();//supernova event
  TCanvas *can66 = new TCanvas();//optional higher energy sn
  TCanvas *can67 = new TCanvas();//optional higher energy sn -vuv
  TCanvas *can61 = new TCanvas();//sn vuv only
  TCanvas *can62 = new TCanvas();//sn vis only


  h_pmt_plane0->SetLineColor(cj);
  h_pmt_plane0->GetZaxis()->SetRangeUser(min_x,max_x);
  h_pmt_plane02->SetLineColor(cj);
  h_pmt_plane02->GetZaxis()->SetRangeUser(min_x,max_x);
  h_pmt_plane1->SetLineColor(cj);
  h_pmt_plane1->GetZaxis()->SetRangeUser(min_x,max_x);
  h_pmt_plane2->SetLineColor(cj);
  h_pmt_plane2->GetZaxis()->SetRangeUser(min_x,max_x);

  can60->cd();

  h_pmt_plane0->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane0->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane0->SetTitle("PMT Signal Map - Supernova 5MeV - Total");
  h_pmt_plane0->SetStats(kFALSE);
  h_pmt_plane0->Draw("colz");

  can66->cd();

  h_pmt_plane02->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane02->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane02->SetTitle("PMT Signal Map - Supernova 22MeV - Total");
  h_pmt_plane02->SetStats(kFALSE);
  h_pmt_plane02->Draw("colz");

  can67->cd();

  h_pmt_plane12->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane12->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane12->SetTitle("PMT Signal Map - Supernova 22MeV - VUV");
  h_pmt_plane12->SetStats(kFALSE);
  h_pmt_plane12->Draw("colz");

  can61->cd();

  h_pmt_plane1->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane1->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane1->SetTitle("PMT Signal Map - Supernova 5MeV - VUV");
  h_pmt_plane1->SetStats(kFALSE);
  h_pmt_plane1->Draw("colz");

  can62->cd();

  h_pmt_plane2->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane2->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane2->SetTitle("PMT Signal Map - Supernova 5MeV - VIS");
  h_pmt_plane2->SetStats(kFALSE);
  h_pmt_plane2->Draw("colz");

//------------------------------------------------------------
//****2D Histograms of events "displayed"****************

  TCanvas *can63 = new TCanvas();//signal map ar39 total
  TCanvas *can64 = new TCanvas();//ar39 vuv
  TCanvas *can65 = new TCanvas();//ar39 vis

  h_pmt_plane_ar0->SetLineColor(cj);
  h_pmt_plane_ar0->GetZaxis()->SetRangeUser(min_x,max_x);
  h_pmt_plane_ar1->SetLineColor(cj);
  h_pmt_plane_ar1->GetZaxis()->SetRangeUser(min_x,max_x);
  h_pmt_plane_ar2->SetLineColor(cj);
  h_pmt_plane_ar2->GetZaxis()->SetRangeUser(min_x,max_x);

  can63->cd();

  h_pmt_plane_ar0->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane_ar0->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane_ar0->SetTitle("PMT Signal Map - Ar39 0.5 MeV - Total");
  h_pmt_plane_ar0->SetStats(kFALSE);
  h_pmt_plane_ar0->Draw("colz");

  can64->cd();

  h_pmt_plane_ar1->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane_ar1->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane_ar1->SetTitle("PMT Signal Map - Ar39 0.5 MeV - VUV");
  h_pmt_plane_ar1->SetStats(kFALSE);
  h_pmt_plane_ar1->Draw("colz");

  can65->cd();

  h_pmt_plane_ar2->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane_ar2->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane_ar2->SetTitle("PMT Signal Map - Ar39 0.5 MeV - VIS");
  h_pmt_plane_ar2->SetStats(kFALSE);
  h_pmt_plane_ar2->Draw("colz");

//------------------------------------------------------------
//****2D Histograms of events "displayed"****************

  TCanvas *can_rd0 = new TCanvas();//signal map rd total
  TCanvas *can_rd1 = new TCanvas();//ar39 rd
  TCanvas *can_rd2 = new TCanvas();//ar39 rd

  h_pmt_plane_rd0->SetLineColor(cj);
  h_pmt_plane_rd0->GetZaxis()->SetRangeUser(min_x,max_x);
  h_pmt_plane_rd1->SetLineColor(cj);
  h_pmt_plane_rd1->GetZaxis()->SetRangeUser(min_x,max_x);
  h_pmt_plane_rd2->SetLineColor(cj);
  h_pmt_plane_rd2->GetZaxis()->SetRangeUser(min_x,max_x);

  can_rd0->cd();

  h_pmt_plane_rd0->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane_rd0->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane_rd0->SetTitle("PMT Signal Map - Rn 6.9MeV - Total");
  h_pmt_plane_rd0->SetStats(kFALSE);
  h_pmt_plane_rd0->Draw("colz");

  can_rd1->cd();

  h_pmt_plane_rd1->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane_rd1->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane_rd1->SetTitle("PMT Signal Map - Rn 6.9MeV - VUV");
  h_pmt_plane_rd1->SetStats(kFALSE);
  h_pmt_plane_rd1->Draw("colz");

  can_rd2->cd();

  h_pmt_plane_rd2->GetXaxis()->SetTitle("z [cm]");
  h_pmt_plane_rd2->GetYaxis()->SetTitle("y [cm]");
  h_pmt_plane_rd2->SetTitle("PMT Signal Map - Ar39 - VIS");
  h_pmt_plane_rd2->SetStats(kFALSE);
  h_pmt_plane_rd2->Draw("colz");

//****END :: 2D Histograms of events "displayed"****************
//------------------------------------------------------------

//***Energy plotting code for Ar39: no cut, 1 pe, 2pe, 2pe+2pmt
  TCanvas *can_e = new TCanvas();
  can_e->cd();
  h_energy->GetXaxis()->SetTitle("Energy (MeV)");
  h_energy->SetFillColorAlpha(46,1.0);
  h_energy->SetStats(kFALSE);
  h_energy1->SetStats(kFALSE);
  h_energy2->SetStats(kFALSE);
  h_energy22->SetStats(kFALSE);
  //h_energy->SetLineColor(46);
  h_energy->SetFillStyle(3003);
  h_energy->Draw();
  h_energy1->SetFillColor(46);
  h_energy1->Draw("same");
  h_energy2->SetFillColor(49);
  h_energy2->Draw("same");
  h_energy22->SetFillColor(2);
  h_energy22->Draw("same");

  gStyle->SetOptStat(0);
  TLegend *leg = new TLegend(0.63,0.73,0.88,0.88);
  //leg->SetHeader("Ar39");
  leg->AddEntry(h_energy," Truth ", "f");
  leg->AddEntry(h_energy1," > 1 #gamma ", "f");
  leg->AddEntry(h_energy2, " > 3 #gamma", "f");
  leg->AddEntry(h_energy22, " > 5 #gamma & on 3 PMT", "f");
  leg->Draw();
//***END :: Ar39 energy plotting code*****
//--------------------------------------------------------
//***SN energy plotting code - same cuts!****

  TCanvas *can_e_sn = new TCanvas();
  can_e_sn->cd();
  h_energy_sn->GetXaxis()->SetTitle("Energy (MeV)");
  h_energy_sn->SetStats(kFALSE);
  h_energy_sn1->SetStats(kFALSE);
  h_energy_sn2->SetStats(kFALSE);
  h_energy_sn22->SetStats(kFALSE);
  h_energy_sn->Draw();
  h_energy_sn1->SetFillColorAlpha(38,1.0);
  //h_energy_sn1->SetLineColor(38);
  h_energy_sn1->SetFillStyle(3003);
  h_energy_sn1->Draw("same");
  h_energy_sn2->SetFillColor(38);
  h_energy_sn2->Draw("same");
  h_energy_sn22->SetFillColor(9);
  h_energy_sn22->Draw("same");

  TLegend *leg2 = new TLegend(0.63,0.73,0.88,0.88);
  //leg2->SetHeader("SN");
  leg2->AddEntry(h_energy_sn, " Truth ", "f");
  leg2->AddEntry(h_energy_sn1, " > 1 #gamma ", "f");
  leg2->AddEntry(h_energy_sn2, " > 3 #gamma", "f");
  leg2->AddEntry(h_energy_sn22, " > 5 #gamma & on 3 PMT", "f");
  leg2->Draw();

//***END :: Plotting code SN energy ***
//------------------------------------------------
//Counting code to find effectiveness of cuts, also useful for tableing results
//***Cuts are in code before plotting section***
//Ar39-------
/*
  int counter_new = 0;
  for(int i = 0; i<num_energy_bins; i++){if(h_energy->GetBinContent(i) > 0){counter_new += h_energy->GetBinContent(i);}}
  cout << "--Ar39--" << endl;
  cout << max_events << endl;
  cout << counter_new << endl;
  cout << "-----------" << endl;
  counter_new = 0;
  for(int i = 0; i<num_energy_bins; i++){if(h_energy1->GetBinContent(i) > 0){counter_new += h_energy1->GetBinContent(i);}}
  cout << counter_new << endl;
  cout << "-----------" << endl;
  counter_new = 0;
  for(int i = 0; i<num_energy_bins; i++){if(h_energy2->GetBinContent(i) > 0){counter_new += h_energy2->GetBinContent(i);}}
  cout << counter_new << endl;
  cout << "-----------" << endl;
  counter_new = 0;
  for(int i = 0; i<num_energy_bins; i++){if(h_energy22->GetBinContent(i) > 0){counter_new += h_energy22->GetBinContent(i);}}
  cout << counter_new << endl;
  cout << "-----------" << endl;


//SN--------
  int counter_new1 = 0;
  for(int i = 0; i<num_energy_bins; i++){if(h_energy_sn->GetBinContent(i) > 0){counter_new1 += h_energy_sn->GetBinContent(i);}}
  cout << "--SN--" << endl;
  cout << max_events_sn << endl;
  cout << counter_new1 << endl;
  cout << "-----------" << endl;
  counter_new1 = 0;
  for(int i = 0; i<num_energy_bins; i++){if(h_energy_sn1->GetBinContent(i) > 0){counter_new1 += h_energy_sn1->GetBinContent(i);}}
  cout << counter_new1 << endl;
  cout << "-----------" << endl;
  counter_new1 = 0;
  for(int i = 0; i<num_energy_bins; i++){if(h_energy_sn2->GetBinContent(i) > 0){counter_new1 += h_energy_sn2->GetBinContent(i);}}
  cout << counter_new1 << endl;
  cout << "-----------" << endl;
  counter_new1 = 0;
  for(int i = 0; i<num_energy_bins; i++){if(h_energy_sn22->GetBinContent(i) > 0){counter_new1 += h_energy_sn22->GetBinContent(i);}}
  cout << counter_new1 << endl;
  cout << "-----------" << endl;
*/
//-----End counting code-------------------
//----Multiple stacked histograms for the number of PMTs seeing signal before
//----and after certain cuts for both the SN and Ar39


  TCanvas *can50 = new TCanvas();
  TCanvas *can51 = new TCanvas();
  TCanvas *can52 = new TCanvas();


  can50->cd();
  h_no_pmt->GetXaxis()->SetTitle("PMTs with signal");
  //h_no_pmt->Draw();


  can51->cd();
  h_no_pmt_sn->GetXaxis()->SetTitle("PMTs with signal");
  //h_no_pmt_sn->Draw();

//-------------------------------------------
  can52->cd();

  h_no_pmt3->GetXaxis()->SetTitle("PMTs with signal");
  h_no_pmt3->SetLineColor(46);
  h_no_pmt3->SetFillColorAlpha(46,1);
  h_no_pmt3->SetFillStyle(3005);
  h_no_pmt3->SetStats(kFALSE);
  h_no_pmt_stack3->Add(h_no_pmt_sn3);

  h_no_pmt_sn3->SetLineColor(9);
  h_no_pmt_sn3->SetFillColorAlpha(9,1);
  h_no_pmt_sn3->SetFillStyle(3005);
  h_no_pmt_sn3->SetStats(kFALSE);
  h_no_pmt_stack3->Add(h_no_pmt3);


  h_no_pmt_stack3->Draw();
  h_no_pmt_stack3->GetXaxis()->SetTitle("PMTs with signal");
  h_no_pmt_stack3->GetYaxis()->SetTitle("Events");
  //h_no_pmt_stack3->SetTitle("");
//-------------------------------------------
  h_no_pmt103->SetLineColor(46);
  h_no_pmt103->SetFillColorAlpha(2,1);
  h_no_pmt103->SetStats(kFALSE);
  //h_no_pmt103->SetFillStyle(3002);
  h_no_pmt_stack103->Add(h_no_pmt103);

  h_no_pmt_sn103->SetLineColor(9);
  h_no_pmt_sn103->SetFillColorAlpha(9,1);
  //h_no_pmt_sn103->SetFillStyle(3002);
  h_no_pmt_sn103->SetStats(kFALSE);
  h_no_pmt_stack103->Add(h_no_pmt_sn103);


  h_no_pmt_stack103->Draw("same");

//-------------------------------------------
  h_no_pmt10->SetLineColor(46);
  h_no_pmt10->SetFillColorAlpha(2,1);
  h_no_pmt10->SetFillStyle(3002);
  h_no_pmt10->SetStats(kFALSE);
  h_no_pmt_stack10->Add(h_no_pmt10);


  h_no_pmt_sn10->SetLineColor(9);
  h_no_pmt_sn10->SetFillColorAlpha(9,1);
  h_no_pmt_sn10->SetFillStyle(3002);
  h_no_pmt_sn10->SetStats(kFALSE);
  h_no_pmt_stack10->Add(h_no_pmt_sn10);



  h_no_pmt_stack10->Draw("same");
  h_no_pmt_stack10->GetXaxis()->SetTitle("PMTs with signal");
  can52->Update();
//-------------------------------------------
  h_no_pmt->SetLineColor(46);//SetFillColorAlpha(46,0.05);
  h_no_pmt->SetStats(kFALSE);
  //h_no_pmt_stack->Add(h_no_pmt);

  h_no_pmt_sn->SetLineColor(9);//SetFillColorAlpha(9,0.05);
  h_no_pmt_sn->SetStats(kFALSE);
  //h_no_pmt_stack->Add(h_no_pmt_sn);

  can52->SetLogy();

  h_no_pmt_stack->Draw("same");


/*//--Legend detailing cuts
  gStyle->SetOptStat(0);
  TLegend *leg3 = new TLegend(0.62,0.73,0.87,0.88);
  leg3->SetHeader("Ar39");
  //leg3->AddEntry(h_no_pmt," Truth ", "f");
  leg3->AddEntry(h_no_pmt3," > 1 #gamma ", "f");
  leg3->AddEntry(h_no_pmt10, " > 2 #gamma ", "f");
  leg3->AddEntry(h_no_pmt103, " > 5 #gamma & on 3 PMT", "f");
  leg3->Draw();

  TLegend *leg4 = new TLegend(0.62,0.58,0.87,0.72);
  leg4->SetHeader("SN");
  //leg4->AddEntry(h_no_pmt_sn, " Truth ", "f");
  leg4->AddEntry(h_no_pmt_sn3, " > 1 #gamma ", "f");
  leg4->AddEntry(h_no_pmt_sn10, " > 2 #gamma ", "f");
  leg4->AddEntry(h_no_pmt_sn103, " > 5 #gamma & on 3 PMT", "f");
  leg4->Draw();
*/
//--Legend detailing cuts
  gStyle->SetOptStat(0);
  TLegend *leg3 = new TLegend(0.42,0.73,0.67,0.88);
  leg3->SetHeader("Ar39");
  //leg3->AddEntry(h_no_pmt," Truth ", "f");
  leg3->AddEntry(h_no_pmt3," > 1 #gamma ", "f");
  leg3->AddEntry(h_no_pmt10, " > 2 #gamma ", "f");
  leg3->AddEntry(h_no_pmt103, " > 5 #gamma & on 3 PMT", "f");
  leg3->Draw();

  TLegend *leg4 = new TLegend(0.42,0.58,0.67,0.72);
  leg4->SetHeader("SN");
  //leg4->AddEntry(h_no_pmt_sn, " Truth ", "f");
  leg4->AddEntry(h_no_pmt_sn3, " > 1 #gamma ", "f");
  leg4->AddEntry(h_no_pmt_sn10, " > 2 #gamma ", "f");
  leg4->AddEntry(h_no_pmt_sn103, " > 5 #gamma & on 3 PMT", "f");
  leg4->Draw();





  TCanvas *can_pmt_d1 = new TCanvas();
  TCanvas *can_pmt_d2 = new TCanvas();
  //TCanvas *can_pmt_d3 = new TCanvas();
  //TCanvas *can_pmt_d4 = new TCanvas();

  can_pmt_d1->cd();
  h_no_pmt_distance_cut1->SetTitle("Ar39: Distance to Cathode vs No. PMTs w/ Signal - VUV");
  //h_no_pmt_distance->Draw();
  //can_pmt_d2->cd();
  h_no_pmt_distance_cut1->GetXaxis()->SetTitle("Number PMTs w/ Signal");
  h_no_pmt_distance_cut1->GetYaxis()->SetTitle("Distance to Cathode [cm]");
  h_no_pmt_distance_cut1->Draw();
  //can_pmt_d3->cd();
  h_no_pmt_distance_cut2->SetMarkerColor(kBlue);
  h_no_pmt_distance_cut2->SetLineColor(kBlue);
  h_no_pmt_distance_cut2->Draw("same");
  //can_pmt_d4->cd();
  h_no_pmt_distance_cut3->SetMarkerColor(kRed);
  h_no_pmt_distance_cut3->SetLineColor(kRed);
  h_no_pmt_distance_cut3->Draw("same");

  TLegend *leg5 = new TLegend(0.62,0.73,0.87,0.88);
  leg5->SetHeader();
  leg5->AddEntry(h_no_pmt_distance_cut1, " > 1 #gamma ", "l");
  leg5->AddEntry(h_no_pmt_distance_cut2, " > 2 #gamma ", "l");
  leg5->AddEntry(h_no_pmt_distance_cut3, " > 3 #gamma & on 3 PMT", "l");
  leg5->Draw();


  can_pmt_d2->cd();

  h_no_pmt_distance_cut1_sn->SetTitle("SN: Distance to Cathode vs No. PMTs w/ Signal - VUV");
  h_no_pmt_distance_cut1_sn->GetXaxis()->SetTitle("Number PMTs w/ Signal");
  h_no_pmt_distance_cut1_sn->GetYaxis()->SetTitle("Distance to Cathode [cm]");
  h_no_pmt_distance_cut1_sn->Draw();
  h_no_pmt_distance_cut2_sn->SetMarkerColor(kBlue);
  h_no_pmt_distance_cut2_sn->SetLineColor(kBlue);
  h_no_pmt_distance_cut2_sn->Draw("same");
  h_no_pmt_distance_cut3_sn->SetMarkerColor(kRed);
  h_no_pmt_distance_cut3_sn->SetLineColor(kRed);
  h_no_pmt_distance_cut3_sn->Draw("same");

  TLegend *leg6 = new TLegend(0.62,0.73,0.87,0.88);
  leg6->SetHeader();
  leg6->AddEntry(h_no_pmt_distance_cut1_sn, " > 1 #gamma ", "l");
  leg6->AddEntry(h_no_pmt_distance_cut2_sn, " > 2 #gamma ", "l");
  leg6->AddEntry(h_no_pmt_distance_cut3_sn, " > 3 #gamma & on 3 PMT", "l");
  leg6->Draw();

//----------- End plotting codes --------

  //close the pmt location datafile
  myfile.close();

//Save histograms to file

   can_e->Print("ar39_energy.pdf");
   can_e_sn->Print("sn_energy.pdf"); 
   can60->Print("sn_signal_map.pdf");
   can61->Print("sn_signal_map_vuv.pdf");
   can62->Print("sn_signal_map_vis.pdf");
   can63->Print("ar39_signal_map.pdf");
   can64->Print("ar39_signal_map_vuv.pdf");
   can65->Print("ar39_signal_map_vis.pdf");
   can66->Print("sn_signal_map_higherE.pdf");
   can67->Print("sn_signal_map_higherE_vuv.pdf");
   can_rd0->Print("rd_signal_map.pdf");
   can_rd1->Print("rd_signal_map_vuv.pdf");
   can_rd2->Print("rd_signal_map_vis.pdf");
   can52->Print("num_pmts.pdf");
   can_pmt_d1->Print("num_pmts_distance.pdf");
   can_marker_ar->Print("ar_marker.pdf");
   can_marker_sn->Print("sn_marker.pdf");
   can_marker_hE->Print("sn_marker_higherE.pdf");
   can_marker_rd->Print("rd_marker.pdf");
   can_pmt_d2->Print("num_pmts_distance_sn.pdf");

  return 0;

}//end main



