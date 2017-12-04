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


int main()
{
  //read in file 1
  const char* filename1 = "10k_sn_randompos2.root";
  TFile * infile1=new TFile(filename1,"READ"); // define the input file
  TTree * data_tree_sn=(TTree *)infile1->Get("data_tree"); //directs towards to root tree inside the input file
  TTree * event_tree_sn=(TTree *)infile1->Get("event_tree");

  int data_event_sn;
  double data_time_sn;
  data_tree_sn->SetBranchAddress("data_event", &data_event_sn);
  data_tree_sn->SetBranchAddress("data_time", &data_time_sn);

  //read in file 2
  const char* filename2 = "10k_rn_randompos2.root";
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
 
  //for photon arrival times per event
  TH1F *sn_currentEvt_times = new TH1F("sn_currentEvt_times","Photon arrival times",2500,0,5); //500 = 10ns bins
  TH1F *rn_currentEvt_times = new TH1F("rn_currentEvt_times","Photon arrival times",2500,0,5); //500 = 10ns bins



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


  double sn_I_all_current; 
  double sn_I_153ns_current; 
  double sn_Fp_current; 

  double rn_I_all_current; 
  double rn_I_153ns_current; 
  double rn_Fp_current; 






  //////// Main loop ////////
  
  int current_sn_photon = 0;
  int current_rn_photon = 0;

  if(sn_evts != rn_evts){cerr << "ERROR, files do not have the same number of evts.\n Unable to analyse." << endl;}
  else{
    int evts = sn_evts;
    for(int evt = 0; evt < 1; evt++){

      //SN loop
      for(int i = current_sn_photon; i < sn_photons; i++){
	data_tree_sn->GetEntry(i);
	if(data_event_sn == evt){	
	  sn_currentEvt_times->Fill(data_time_sn);
	}
	else {current_sn_photon = i; break;}  
      }
      //Radon loop
      for(int j = current_rn_photon; j < rn_photons; j++){
	data_tree_rn->GetEntry(j);
	if(data_event_rn == evt){
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
      //sn_currentEvt_times->Reset();
      rn_currentEvt_times->Reset();


      //print out when analysis of every 100 events is complete
      if((evt + 1) % 100 == 0) {cout << "Analysis of event " << evt + 1 << " complete." << endl;}
    } // end of evt loop

   


    TFile *output=new TFile("output.root","RECREATE");
  
    output->cd();
    sn_currentEvt_times->Write();
    sn_Fp->Write();
    rn_Fp->Write();

    output->Close();
  


  }//end of if(sn_evts == rn_evts)



  return 0;

}//end main























