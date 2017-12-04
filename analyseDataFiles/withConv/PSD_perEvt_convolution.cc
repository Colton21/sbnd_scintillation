#include <vector>
#include "unistd.h"
#include <stdio.h>

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

//#include "utility_functions.h"

using namespace std;


 //convolution fn
  pair <vector<double>,vector<double> >
    Convolve(pair <vector<double>,vector<double> >& a, pair <vector<double>,vector<double> >& b)
  {
 
    const unsigned int na = a.first.size();
    const unsigned int nb = b.first.size();
    //if (na < 2 || nb < 2)
    //cout<<"sizes too small";
  
    const double dx = 1.;
    const unsigned int nr = na + nb;
    vector<double> vx;
    vector<double> vy;
    vx.resize(nr);
    vy.resize(nr);
  
    const double xmin = a.first.at(0) + b.first.at(0);
  
    for (int j = 0, n = nr; j < n; ++j) {
      for (int i = max(0, j - int(nb) + 1), m = min(int(na), j + 1); i < m; ++i) {
	const int ji = j - i;
	vy[j] += a.second.at(i) * b.second.at(ji);
      }
      vy[j] *= dx;
      vx[j] = xmin + j*dx;
    }
  
    pair <vector<double>,vector<double> > r = make_pair(vx,vy);
  
    return r;
  }


  //downsample fn
pair <vector<double>,vector<double> >
Downsample(pair <vector<double>,vector<double> >& xy, const double sampling_time, const double scaleY)
{
  const unsigned int n = xy.first.size();
 
  int min_time = floor(xy.first.at(0));
  int max_time = ceil(xy.first.at(n-1));
  if(min_time%int(sampling_time))
    min_time -= 1;
  if(max_time%int(sampling_time))
    max_time += 1;

  const unsigned int m = (max_time - min_time)/sampling_time;

  //cout<<"min_time:"<<min_time<<" max_time:"<<max_time<<"  m:"<<m<<endl;
  vector<double> vx;
  vector<double> vy;
  vx.resize(m);
  vy.resize(m);
  
  for (unsigned int i = 0, j = 0; j < m; j++) {
    double x = min_time + sampling_time * j;
    double y = 0;
    while (xy.first.at(i) < (x + sampling_time)) {
      //cout<<i<<"  "<<j<<"  "<<xy.first.at(i)<<"  "<<x<<"  "<<sampling_time<<"  "<<min_time<<"  "<<sampling_time * j<<" "<<(x + sampling_time)<<endl;
      y += xy.second.at(i);
      if(i < n-1)
	++i;
      else
	break;
    }
    vx[j] = x + sampling_time/2;//moving the time to the center of the bin
    vy[j] = scaleY * y;

  }
 
  pair <vector<double>,vector<double> > res = make_pair(vx,vy);
 
  return res;
  }




int main()
{
  //read in file 1
  const char* filename1 = "../10k_sn_randompos2.root";
  TFile * infile1=new TFile(filename1,"READ"); // define the input file
  TTree * data_tree_sn=(TTree *)infile1->Get("data_tree"); //directs towards to root tree inside the input file
  TTree * event_tree_sn=(TTree *)infile1->Get("event_tree");

  int data_event_sn;
  double data_time_sn;
  data_tree_sn->SetBranchAddress("data_event", &data_event_sn);
  data_tree_sn->SetBranchAddress("data_time", &data_time_sn);

  //read in file 2
  const char* filename2 = "../10k_rn_randompos2.root";
  TFile * infile2=new TFile(filename2,"READ"); // define the input file
  TTree * data_tree_rn=(TTree *)infile2->Get("data_tree"); //directs towards to root tree inside the input file
  TTree * event_tree_rn=(TTree *)infile2->Get("event_tree");

  int data_event_rn;
  double data_time_rn;
  data_tree_rn->SetBranchAddress("data_event", &data_event_rn);
  data_tree_rn->SetBranchAddress("data_time", &data_time_rn);


  /////// Define histograms ///////
  //for F prompt
  TH1F *sn_Fp = new TH1F("sn_Fp","F prompt for SN events",500,0,1); 
  TH1F *rn_Fp = new TH1F("rn_Fp","F prompt for Rn events",500,0,1); 
  TH1F *sn_Fp_conv = new TH1F("sn_Fp_conv","F prompt for SN events (convoluted)",500,0,1); 
  TH1F *rn_Fp_conv = new TH1F("rn_Fp_conv","F prompt for Rn events (convoluted)",500,0,1);

  
  //for photon arrival times per event
  TH1F *sn_currentEvt_times = new TH1F("sn_currentEvt_times","Photon arrival times",2500,0,5); //500 = 10ns bins
  TH1F *rn_currentEvt_times = new TH1F("rn_currentEvt_times","Photon arrival times",2500,0,5); //500 = 10ns bins

  TH1F *sn_times_conv = new TH1F("sn_times_conv","Convolved photon arrival times",2500,0,5000); //500 = 10ns bins
  TH1F *rn_times_conv = new TH1F("rn_times_conv","Convolved photon arrival times",2500,0,5000);
  TH1F *sn_times_conv_currentEvt = new TH1F("sn_times_conv_currentEvt","Convolved photon arrival times",2500,0,5000); //500 = 10ns bins
  TH1F *rn_times_conv_currentEvt = new TH1F("rn_times_conv_currentEvt","Convolved photon arrival times",2500,0,5000);

//SER definition:
  const double sampling_time = 2.;//ns
  TF1 *fser = new TF1("fser", "gaus(0)",-100.,100.);
  double norm = 32.;//adc
  double mean = 0.;
  double sigma = 50.;
  fser->SetParameter(0,norm);
  fser->SetParameter(1,mean);
  fser->SetParameter(2,sigma);
  vector<double> vx;
  vector<double> vy;
  const int npoints = 100;
  
  for(int i=0; i<npoints; i++) {
    double t = i - 0.5*npoints;
    vx.push_back(t);
    vy.push_back(fser->Eval(t));
  }
  pair <vector<double>,vector<double> >  ser = make_pair(vx,vy);


 //////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////// ***** Analysis ***** ///////////////////////////////   
 //////////////////////////////////////////////////////////////////////////////////

  int sn_evts = event_tree_sn->GetEntries();
  int sn_photons = data_tree_sn->GetEntries();
  cout << "SN evts: " << sn_evts << "\t";
  cout << "SN photons: " << sn_photons << endl;

  int rn_evts = event_tree_sn->GetEntries();
  int rn_photons = data_tree_rn->GetEntries();
  cout << "Rn evts: " << rn_evts << "\t";
  cout << "Rn photons: " << rn_photons << endl;


  double sn_I_all_current; double sn_I_all_current_conv;
  double sn_I_153ns_current; double sn_I_153ns_current_conv;
  double sn_Fp_current; double sn_Fp_current_conv;

  double rn_I_all_current; double rn_I_all_current_conv;
  double rn_I_153ns_current; double rn_I_153ns_current_conv;
  double rn_Fp_current; double rn_Fp_current_conv;

  vector<double> current_snevt_times;
  vector<double> current_rnevt_times;




  //////// Main loop ////////
  
  int current_sn_photon = 0;
  int current_rn_photon = 0;

  if(sn_evts != rn_evts){cerr << "ERROR, files do not have the same number of evts.\n Unable to analyse." << endl;}
  else{
    int evts = sn_evts;
    for(int evt = 0; evt < evts; evt++){

      //SN loop
      for(int i = current_sn_photon; i < sn_photons; i++){
	data_tree_sn->GetEntry(i);
	if(data_event_sn == evt){
	  current_snevt_times.push_back(data_time_sn*1000);
	  sn_currentEvt_times->Fill(data_time_sn);
	}
	else {current_sn_photon = i; break;}  
      }
      //Radon loop
      for(int j = current_rn_photon; j < rn_photons; j++){
	data_tree_rn->GetEntry(j);
	if(data_event_rn == evt){
	  current_rnevt_times.push_back(data_time_rn*1000);
	  rn_currentEvt_times->Fill(data_time_rn);
	}
	else {current_rn_photon = j; break;}  
      }
      

      //find the bin with 153 ns center
      //for SN
      int sn_bins = sn_currentEvt_times->GetSize()-2;
      double sn_currentBinTime = 0;
      int sn_bin153ns = 0;
      for(int i = 0; i < sn_bins; i++) {
	sn_currentBinTime = sn_currentEvt_times->GetBinCenter(i);
	if(sn_currentBinTime != 0.153){continue;}
	else{ sn_bin153ns = i; }
      }
     

      //for Radon
      int rn_bins = rn_currentEvt_times->GetSize()-2;
      double rn_currentBinTime = 0;
      int rn_bin153ns = 0;
      for(int j = 0; j < rn_bins; j++) {
	rn_currentBinTime = rn_currentEvt_times->GetBinCenter(j);
	if(rn_currentBinTime != 0.153){continue;}
	else{ rn_bin153ns = j; }
      }

      ////////////////////////
      //NON CONVOLUTED TIMES
      ////////////////////////

      //find intergrals for all light and 0 -> 153 ns
      if(sn_bin153ns == 0 || rn_bin153ns == 0 ){ cerr << "Unable to find bin at 153 ns for event: " << evt << endl; }
      else { 
	sn_I_all_current = sn_currentEvt_times->Integral(); // integrate over 0 -> 5 mu_s	
	sn_I_153ns_current = sn_currentEvt_times->Integral(0,sn_bin153ns);
	
	rn_I_all_current = rn_currentEvt_times->Integral();
	rn_I_153ns_current = rn_currentEvt_times->Integral(0,rn_bin153ns);
      }

      sn_Fp_current = sn_I_153ns_current/sn_I_all_current;
      rn_Fp_current = rn_I_153ns_current/rn_I_all_current;

      sn_Fp->Fill(sn_Fp_current);
      rn_Fp->Fill(rn_Fp_current);

      //reset for next loop
      sn_Fp_current = 0;
      rn_Fp_current = 0;
      //Reset histograms
      sn_currentEvt_times->Reset();
      rn_currentEvt_times->Reset();



      ////////////////////////
      //CONVOLVED TIMES///////
      ////////////////////////

      //convolution analysis here //
      sort(current_snevt_times.begin(), current_snevt_times.end());
      sort(current_rnevt_times.begin(), current_rnevt_times.end());

      vector<double> sn_signaly(current_snevt_times.size(), 1.);         
      pair <vector<double>,vector<double> > sn_signal = make_pair(current_snevt_times,sn_signaly); 
      pair <vector<double>,vector<double> > sn_signal_1ns = Downsample(sn_signal, 1, 1);
      pair <vector<double>,vector<double> > sn_fadc = Convolve(sn_signal_1ns, ser);
      //pair <vector<double>,vector<double> > fadc_2ns = Downsample(fadc, 2, 1);

      vector<double> rn_signaly(current_rnevt_times.size(), 1.);
      pair <vector<double>,vector<double> > rn_signal = make_pair(current_rnevt_times,rn_signaly); 
      pair <vector<double>,vector<double> > rn_signal_1ns = Downsample(rn_signal, 1, 1);
      pair <vector<double>,vector<double> > rn_fadc = Convolve(rn_signal_1ns, ser);

      for(int i=0; i<int(sn_fadc.first.size()); i++) {
	sn_times_conv->Fill(sn_fadc.first.at(i), sn_fadc.second.at(i));
	sn_times_conv_currentEvt->Fill(sn_fadc.first.at(i), sn_fadc.second.at(i));
      }

      for(int i=0; i<int(rn_fadc.first.size()); i++) {
	rn_times_conv->Fill(rn_fadc.first.at(i), rn_fadc.second.at(i));
	rn_times_conv_currentEvt->Fill(rn_fadc.first.at(i), rn_fadc.second.at(i));
      }


      current_snevt_times.clear();
      current_rnevt_times.clear();


      //find the bin with 153 ns center
      //for SN
      int sn_bins_conv = sn_times_conv_currentEvt->GetSize()-2;
      double sn_currentBinTime_conv = 0;
      int sn_bin153ns_conv = 0;
      for(int i = 0; i < sn_bins_conv; i++) {
	sn_currentBinTime_conv = sn_times_conv_currentEvt->GetBinCenter(i);
	if(sn_currentBinTime_conv != 153){continue;}
	else{ sn_bin153ns_conv = i; }
      }
     

      //for Radon
      int rn_bins_conv = rn_times_conv_currentEvt->GetSize()-2;
      double rn_currentBinTime_conv = 0;
      int rn_bin153ns_conv = 0;
      for(int j = 0; j < rn_bins_conv; j++) {
	rn_currentBinTime_conv = rn_times_conv_currentEvt->GetBinCenter(j);
	if(rn_currentBinTime_conv != 153){continue;}
	else{ rn_bin153ns_conv = j; }
      }


      //find intergrals for all light and 0 -> 153 ns
      if(sn_bin153ns_conv == 0 || rn_bin153ns_conv == 0 ){ cerr << "Unable to find bin at 153 ns for event: " << evt << endl; }
      else { 
	sn_I_all_current_conv = sn_times_conv_currentEvt->Integral(); // integrate over 0 -> 5 mu_s	
	sn_I_153ns_current_conv = sn_times_conv_currentEvt->Integral(0,sn_bin153ns_conv);
	
	rn_I_all_current_conv = rn_times_conv_currentEvt->Integral();
	rn_I_153ns_current_conv = rn_times_conv_currentEvt->Integral(0,rn_bin153ns_conv);
      }

      sn_Fp_current_conv = sn_I_153ns_current_conv/sn_I_all_current_conv;
      rn_Fp_current_conv = rn_I_153ns_current_conv/rn_I_all_current_conv;

      sn_Fp_conv->Fill(sn_Fp_current_conv);
      rn_Fp_conv->Fill(rn_Fp_current_conv);

      //reset for next loop
      sn_Fp_current_conv = 0;
      rn_Fp_current_conv = 0;
      
      sn_times_conv_currentEvt->Reset();
      rn_times_conv_currentEvt->Reset();


      //print out when analysis of every 100 events is complete
      if((evt + 1) % 100 == 0) {cout << "Analysis of event " << evt + 1 << " complete." << endl;}
    } // end of evt loop

   


    TFile *output_conv=new TFile("output_conv.root","RECREATE");
  
    output_conv->cd();
    sn_times_conv->Write();
    rn_times_conv->Write();

    sn_Fp->Write();
    rn_Fp->Write();

    sn_Fp_conv->Write();
    rn_Fp_conv->Write();
    output_conv->Close();
  


  }//end of if(sn_evts == rn_evts)



  return 0;

}//end main























