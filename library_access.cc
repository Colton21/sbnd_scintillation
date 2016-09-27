#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TRandom3.h"
#include "TMath.h"

#include "library_access.h"
#include "utility_functions.h"

using namespace std;

LibraryAccess::LibraryAccess()
	: table_(std::vector<std::vector<float> >()),
	reflected_table_(std::vector<std::vector<float> >()),
	reflT_table_(std::vector<std::vector<float> >())
{

}



void LibraryAccess::LoadLibraryFromFile(std::string libraryfile, bool reflected, bool reflT0)
{
	cout << "Reading photon library from input file: " << libraryfile.c_str()<<endl;

	TFile *f = nullptr;
	TTree *tt = nullptr;

	try
	{
		f  =  TFile::Open(libraryfile.c_str());
		tt =  (TTree*)f->Get("PhotonLibraryData");

		if (!tt) {

			TKey *key = f->FindKeyAny("PhotonLibraryData");
			if (key)
				tt = (TTree*)key->ReadObj();
			else {
				cout << "PhotonLibraryData not found in file" <<libraryfile;
			}
		}
	}
	catch(...)
	{
		cout << "Error in ttree load, reading photon library: " << libraryfile.c_str()<<endl;
	}

	int voxel;
	int opChannel;
	float visibility;
	float reflVisibility;
	float reflT;
	int maxvoxel = tt->GetMaximum("Voxel")+1;
	int maxopChannel = tt->GetMaximum("OpChannel")+2;

	cout << "Photon lookup table size : " <<  maxvoxel << " voxels,  " << maxopChannel <<" channels " << endl;

	table_.resize(maxvoxel, std::vector<float>(maxopChannel, 0));
	reflected_table_.resize(maxvoxel, std::vector<float>(maxopChannel, 0));
	reflT_table_.resize(maxvoxel, std::vector<float>(maxopChannel, 0));


	tt->SetBranchAddress("Voxel",      &voxel);
	tt->SetBranchAddress("OpChannel",  &opChannel);
	tt->SetBranchAddress("Visibility", &visibility);
	if(reflected) {tt->SetBranchAddress("ReflVisibility", &reflVisibility); }
	if(reflT0) {tt->SetBranchAddress("ReflTfirst", &reflT); }

	size_t nentries = tt->GetEntries();

	for(size_t i=0; i!=nentries; ++i)
	{
		tt->GetEntry(i);
		if((voxel<0)||(voxel>= maxvoxel)||(opChannel<0)||(opChannel>= maxopChannel))
		{}
		else
		{
			table_.at(voxel).at(opChannel) = visibility;
			if(reflected) {reflected_table_.at(voxel).at(opChannel) = reflVisibility; }
			else{reflected_table_.at(voxel).at(opChannel) = 0; }
			if(reflT0) {reflT_table_.at(voxel).at(opChannel) = reflT; }
			else{reflT_table_.at(voxel).at(opChannel) = 0; }
		}
	}

	try
	{
		f->Close();
	}
	catch(...)
	{
		cout << "Error in closing file : " << libraryfile.c_str()<<endl;
	}
}

const float* LibraryAccess::GetReflT0(size_t voxel, int no_pmt)
{
	return &reflT_table_.at(voxel).at(no_pmt);
}

const float* LibraryAccess::GetReflCounts(size_t voxel, int no_pmt, bool reflected)
{
	if(reflected) {return &reflected_table_.at(voxel).at(no_pmt); }
	else{return 0; }
}

const float* LibraryAccess::GetCounts(size_t voxel, int no_pmt)
{
	return &table_.at(voxel).at(no_pmt);
}

const float* LibraryAccess::GetLibraryEntries(int voxID, bool reflected, int no_pmt)
{
	if(!reflected)
		return GetCounts(voxID, no_pmt);
	else
		return GetReflCounts(voxID, no_pmt, reflected);
}


vector<int> LibraryAccess::GetVoxelCoords(int id, double position[3])
{
	vector<int> returnvector;
	returnvector.resize(3);
	returnvector.at(0) =  id % gxSteps;
	returnvector.at(1) =  ((id - returnvector.at(0) ) / gxSteps) % gySteps;
	returnvector.at(2) =  ((id - returnvector.at(0) - (returnvector.at(1) * gxSteps)) / (gySteps * gxSteps)) % gzSteps;

	position[0] = gLowerCorner[0] + (returnvector.at(0) + 0.5)*(gUpperCorner[0] - gLowerCorner[0])/gxSteps;
	position[1] = gLowerCorner[1] + (returnvector.at(1) + 0.5)*(gUpperCorner[1] - gLowerCorner[1])/gySteps;
	position[2] = gLowerCorner[2] + (returnvector.at(2) + 0.5)*(gUpperCorner[2] - gLowerCorner[2])/gzSteps;

	return returnvector;

}

//This function takes is most of the information needed to calculate the number
//of photoelectrons on a given PMT
std::vector<double> LibraryAccess::PhotonLibraryAnalyzer(double _energy, const int _scint_yield, const double _quantum_efficiency, int _pmt_number, int _rand_voxel)
{
	//The number of photons created is determined by this formula:
	//Nphotons_created = Poisson < Scintillation Yield (24000/MeV) * dE/dX (MeV)>

	const int scint_yield = _scint_yield;
	const double quantum_efficiency = _quantum_efficiency;
	double energy = _energy;
  int i = _rand_voxel;
  vector<double> pmt_hits;

	int pre_Nphotons_created = scint_yield * energy;

  //Applying quantum efficiency right after calculation reduces total number
  //of photons working with initially
	int Nphotons_created = utility::poisson(pre_Nphotons_created, gRandom->Uniform(1.), energy) * quantum_efficiency;

  //Find position of event
	double position[3];
	vector<int> Coords = GetVoxelCoords(i, position);

  //Look up visibility parameter/timing by comparing the optical channel (PMT Number)
  //and the detector location (voxel, i)
	const float* Visibilities = GetLibraryEntries(i, false, _pmt_number);
	const float* ReflVisibilities = GetLibraryEntries(i, true, _pmt_number);
	const float* reflected_T0 = GetReflT0(i, _pmt_number);

	const float vis = *Visibilities;
	const float reflvis = *ReflVisibilities;

	//int ichan = _pmt_number;

  //Number of photoelectrons for a given PMT
	double hits_vuv = Nphotons_created * vis;
	double hits_vis = Nphotons_created * reflvis;
	double reflT0 = *reflected_T0;

	/// Cast the hits into int type
	int int_hits_vuv = hits_vuv;
	int int_hits_vis = hits_vis;

  //Push information back into a vector to readout in:
  //libraryanalyze_light_histo
	pmt_hits.push_back(int_hits_vuv);
	pmt_hits.push_back(int_hits_vis);
	pmt_hits.push_back(position[0]);//x position of voxel
	pmt_hits.push_back(position[1]);//y
	pmt_hits.push_back(position[2]);//z
	pmt_hits.push_back(reflT0);

	return pmt_hits;
}
