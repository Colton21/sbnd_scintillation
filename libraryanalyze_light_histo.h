#ifndef LIBRARYANALYZE_LIGHT_HISTO_H
#define LIBRARYANALYZE_LIGHT_HISTO_H

#include <vector>

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

#include "utility_functions.h"

using namespace std;

///-------------------------------------
///Custom Sorting Function
///-------------------------------------
bool sort_function(std::pair<double, int> pair1, std::pair<double, int> pair2)
{ return (pair1.first < pair2.first); }
///-------------------------------------

///-------------------------------------
//--------What to generate?-------------
///-------------------------------------
  bool gen_ions = false;
  bool supernova = false;
  bool gen_argon = true;
///-------------------------------------
//------Light System Configuration------
///-------------------------------------
//--------------------------------------
  ///Config List:
  ///0 = Full Foils
  ///1 = Cath Foils
  ///2 = VUV only
  const int config = 1;
//--------------------------------------
//--------------------------------------
//--------------------------------------
  std::string libraryfile;
  bool reflected;
  bool reflT;
//--------------------------------------
//--------------------------------------
//TTree branches and data products:
//-------------------------------------
  TFile data_file("test.root", "RECREATE", "Timing PMT File");
  TFile event_file("event_file.root", "RECREATE", "Event File");
  TTree *data_tree = new TTree("data_tree", "data tree");
  TTree *event_tree = new TTree("event_tree", "event tree");
  double data_time;
  int data_pmt;
  int data_event;
  double data_x_pos;
  double data_y_pos;
  double data_z_pos;

  int event_no;
  int event_vox;
  double event_x_pos;
  double event_y_pos;
  double event_z_pos;
  double event_E;
//--------------------------------------
///-------------------------------------
  const double MassE = 0.510998910; // mass electron - MeV/c^2
  const double Q = 0.565;//Q value of decay - Ar39
///-------------------------------------
///-------------------------------------
  const double t_singlet = 0.000000006; //6ns
  const double t_triplet = 0.0000015; //1.5 us
  const double scint_time_window = 0.00001; //10 us
///-------------------------------------
///-------------------------------------
  const double Eav = 20.;//Average energy for SN spectrum
  const double expected_sn = 2.8;//For poisson weighting
///-------------------------------------
///-------------------------------------
  const int scint_yield = 24000;//Scintillation yield of LAr at 500 V/cm
  const double quantum_efficiency = 0.2;//Expected value
  const double activity = 1.; //Ar39 roughly 1 Bq/kg
  const double volume = 112000.; //SBND 112ton LAr
  const double time_window = (0.0012 * 10.);//1.2 [ms] is the readout window
///-------------------------------------
///-------------------------------------
//----Number of Events------------------
///-------------------------------------
  const int max_events_sn = 150;
  //int max_events_sn = utility::poisson(expected_sn,gRandom->Uniform(1.),1.);
  const int max_events = activity * volume/2 * time_window;//Half volume for 1 TPC
  vector<vector<double>> myfile_data;
//--------------------------------------
//--------------------------------------
//---Of type LibraryAccess-------------
  LibraryAccess lar_light;
//--------------------------------------
//--------------------------------------
//--Lists of variables for generating---
  vector<double> energy_ar_list;
  vector<double> scint_time_list;
  vector<int> voxel_list;
//--------------------------------------
//--------------------------------------
//--Predefined number of PMTs for-------
//--SBND configuration------------------
  int realisticPMT_IDs[60] = {0, 4, 8, 12, 16, 20, 24, 32, 40, 44, 48, 52, 56, 60, 64, 88, 92, 96, 100, 104, 108, 112, 120, 128, 132, 136, 140, 144, 148, 152, 154, 158, 162, 166, 170, 174, 178, 186, 194, 198, 202, 206, 210, 214, 218, 242, 246, 250, 254, 258, 262, 266, 274, 282, 286, 290, 294, 298, 302, 306};
//--------------------------------------
//--------------------------------------
//--For timing parameterization---------
  const double signal_t_range = 1000.;


#endif
