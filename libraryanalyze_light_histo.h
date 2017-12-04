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
//--------WHAT to generate?-------------
///-------------------------------------
bool fixed_energy = false; double fixedE = 5.590; //MeV
bool supernova = true;
bool gen_argon = false;
bool gen_radon = false;
///-------------------------------------
//--------WHERE to generate?-------------
///-------------------------------------
bool random_pos = true;
bool fixed_xpos = false; 
bool fixed_pos = false;
double fixedX = 100; // cm (x = 0 at cathode, x = 200 at PMTs)
double fixedY = 0; // cm (y = -200 bottom of TPC, y = 200 top of TPC)
double fixedZ = 250; // cm (z = 0 front of TPC, x = 500 end of TPC)
///-------------------------------------
//--------time cut?-------------
///-------------------------------------
bool cut = false; // NB you can always make time cuts when you're analysing the files - so I tend not to use this
double time_cut = 0.1; // in microseconds - 0.1 mu_s = 100 ns

///-------------------------------------
//------Light System Configuration------
///-------------------------------------
//--------------------------------------
///Config List:
///0 = Full Foils
///1 = Cath Foils
///2 = VUV only
const int config = 1; // cathode foils is the most likely candidate
//--------------------------------------
//--------------------------------------
//--------------------------------------
// These bools are used in the scintillation timing functions
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
TTree *data_tree_vuv = new TTree("data_tree_vuv", "data tree_vuv");
TTree *data_tree_vis = new TTree("data_tree_vis", "data tree_vis");

TTree *event_tree = new TTree("event_tree", "event tree");
double data_time;
double data_time_vuv;
double data_time_vis;

int data_pmt;
int data_pmt_vuv;
int data_pmt_vis;

int data_event;
int data_event_vuv;
int data_event_vis;

double data_x_pos;
double data_x_pos_vuv;
double data_x_pos_vis;

double data_y_pos;
double data_y_pos_vuv;
double data_y_pos_vis;

double data_z_pos;

int event_no;
int event_vox;
double event_x_pos;
double event_y_pos;
double event_z_pos;
double event_E;
//--------------------------------------
///-------------------------------------
//----LAr & Ar39 properties--------------------
///-------------------------------------
const double MassE = 0.510998910; // mass electron - MeV/c^2
const double Q_Ar = 0.565;//Q value of decay - Ar39

const double t_singlet = 0.000000006; //6ns
const double t_triplet = 0.0000015; //1.5 us
const double scint_time_window = 0.00001; //10 us

const int scint_yield_electron = 24000;//Scintillation yield of LAr at 500 V/cm
const double activity_Ar = 1.; //Ar39 roughly 1 Bq/k
///-------------------------------------
///-------------------------------------
//----SN properties---------------------
///-------------------------------------
const double Eav = 20.;//Average energy for SN spectrum
const double expected_sn = 2.8;//For poisson weighting
///-------------------------------------
///-------------------------------------
//----Radon properties---------------------
///-------------------------------------
const int scint_yield_alpha = 16800; // SY of alpha particles at 500 V/cm - from larsoft
const double activity_Rn = 0.000021; // Bq/kg
const double massAlpha = 3727.3794; // alpha particle mass - MeV/c^2
const double Q_Rn = 5.590; // deposited energy from a radon decay - Rn-222 --> Po-218
///-------------------------------------
///-------------------------------------
//----TPC and PMT properties---------------------
///-------------------------------------
const double quantum_efficiency = 0.2; //Expected value
const double mass = 112000.; //SBND 112ton LAr
const double time_window = 10.; //(0.0012 * 10.);//1.2 [ms] is the readout window
const double time_frames = time_window/0.0012;
///-------------------------------------
///-------------------------------------
//----Number of Events------------------
///-------------------------------------
// Fixed energy (electron like) events:
const int max_events_FE = 1;

// Ar-39 events:
//const int max_events_Ar = 1;
const int max_events_Ar = activity_Ar * mass/2 * time_window;//Half volume for 1 TPC
const int Ar_decays_per_sec = activity_Ar* mass/2; // decay rate in one TPC

// Radon events:
const int max_events_Rn = 1;
//const int max_events_Rn = activity_Rn * mass/2 * time_window;//Half volume for 1 (NOTE: for a small time window, this will probably return 0)
const double Rn_decays_per_sec = activity_Rn* mass/2; // decay rate in one TPC

// Supernova events:
const int max_events_SN = 10000;
//int max_events_SN = utility::poisson(expected_sn,gRandom->Uniform(1.),1.);




// Don't worry about the stuff below too much...
vector<vector<double>> myfile_data;
//--------------------------------------
//--------------------------------------
//---Of type LibraryAccess-------------
LibraryAccess lar_light;
//--------------------------------------
//--------------------------------------
//--Lists of variables for generating---
vector<double> energy_list;
vector<double> decay_time_list;
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
