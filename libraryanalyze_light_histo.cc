#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "library_access.h"
#include "libraryanalyze_light_histo.h"


using namespace std;

int main()
{ 

  //////////////////////////////////////////////////////////////////////////////
  ////////////-----------------FUNCTIONS-----------------------------///////////
  //////////////////////////////////////////////////////////////////////////////

  TF1 *fSpectrum = new TF1("fSpectrum",utility::SpectrumFunction,0,Q_Ar,1);//-----Beta decay spectrum
  TF1 *flandau_sn = new TF1("flandau_sn",utility::fsn, 0, 50, 1);//--SN Nu spectrum
  fSpectrum->SetParameter(0, Q_Ar);
  flandau_sn->SetParameter(0, Eav);

  TRandom3 *fGauss = new TRandom3();
  
  TF1 *fScintillation_function = new TF1("Scintillation Timing", utility::Scintillation_function, 0, scint_time_window, 3);
  fScintillation_function->SetParameter(0, t_singlet);
  fScintillation_function->SetParameter(1, t_triplet);
  // the 3rd parameter is the type of particle 0 = electron, 1 = alpha
  // this needs to be set as it sets the ratio of the fast and slow components
  if(gen_argon == true || fixed_energy == true || supernova == true){ //i.e. an electron is the ionising particle
    fScintillation_function->FixParameter(2, 0); 
  }
  if(gen_radon == true){fScintillation_function->FixParameter(2, 1);} // an alpha particle is the ionising particle


  //////////////////////////////////////////////////////////////////////////////
  ////////////--------------INTRODUCE CONDITIONS---------------------///////////
  //////////////////////////////////////////////////////////////////////////////
  int max_events;
  int scint_yield;
  string particle;

  if(fixed_energy == true) {
    max_events = max_events_FE;
    scint_yield = scint_yield_electron;
    particle = "electron";
    cout << endl << "Generating " << max_events << " events, of fixed energy: " << fixedE << " MeV." << endl;
  }
  if(gen_argon == true) {
    max_events = max_events_Ar;
    scint_yield = scint_yield_electron;
    particle = "electron";
    cout << endl << "Generating " << max_events << ", Ar 39 decays in time window: " << time_window << " seconds." << endl;
    cout << "This is equal to " << time_frames << " PMT readout frames." << endl;

    cout << endl <<  "/////////////////////////////////////////////////////////////////////////////" << endl;
    cout << "***NOTE. In ONE TPC, We expect to see: " << Ar_decays_per_sec << " Ar 39 decays each second.***" << endl;
    cout << "//////////////////////////////////////////////////////////////////////////////" << endl << endl;
  }
  if(supernova == true){
    max_events = max_events_SN;
    scint_yield = scint_yield_electron;
    particle = "electron";
    cout << "\nGenerating " << max_events << ", supernova events.\n";
  }
  if(gen_radon == true){ //gen_radon == true
    max_events = max_events_Rn;
    scint_yield = scint_yield_alpha;
    particle = "alpha";
    cout << "\nGenerating " << max_events << ", Radon 222 decays in time window: " << time_window << " seconds." << endl;
    cout << "This is equal to " << time_frames << " PMT readout frames." << endl;

    cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
    cout << "***NOTE. In ONE TPC, we expect to see: " << Rn_decays_per_sec << " radon-222 decays each second.***" << endl;
    cout << "/////////////////////////////////////////////////////////////////////////////" << endl << endl;
  }


  cout << "With the following set up: ";
  if(config == 0) {cout << "Full Foils" << endl; }
  if(config == 1) {cout << "Cath Foils" << endl; }
  if(config == 2) {cout << "VUV only" << endl; }
  cout << endl;
  if (config == 0 || config ==1)
    {
      reflected = true;
      reflT = true;
    }
  if (config ==2)
    {
      reflected = false;
      reflT = false;
    }


  ////////////////////////////////////////////////////////////////////////////
  ////////////--------------CREATING THE TREES-----------------///////////////
  ////////////////////////////////////////////////////////////////////////////
  
  // IMPORTANT \\
  // data_tree will store entry for each and every photon from ALL events simulated
  // event_tree will store an entry for each event 
  // => data_tree is MUCH larger as each event will produce ~20,000 photons (but only a fraction of these are "detected")
  // It is the libraries which determine how many photons are detected by each PMT


  //data_tree was designed such that each photoelectron
  //would have their own entry in the Ttree with a time, which pmt, xyz, and
  //event number (data_event)
  //A current issue with the event number tracking is that it does not carry
  //over between several iterations of the code. For single uses this is fine,
  //but memory consuption and speed for running this code a single time
  //becomes a limiting factor, such that running many instances of this code
  //in parallel becomes much more useful. Therefore carrying over the event
  //number between iterations may be useful.

  // data tree for total light
  data_tree->Branch("data_time", &data_time, "data_time/D");
  data_tree->Branch("data_pmt", &data_pmt, "data_pmt/I");
  data_tree->Branch("data_x_pos", &data_x_pos, "data_x_pos/D");
  data_tree->Branch("data_y_pos", &data_y_pos, "data_y_pos/D");
  data_tree->Branch("data_z_pos", &data_z_pos, "data_z_pos/D"); 
  data_tree->Branch("data_event", &data_event, "data_event/I");

  // datatree for VUV light 
  data_tree_vuv->Branch("data_time_vuv", &data_time_vuv, "data_time_vuv/D");
  data_tree_vuv->Branch("data_pmt_vuv", &data_pmt_vuv, "data_pmt_vuv/I");
  data_tree_vuv->Branch("data_x_pos_vuv", &data_x_pos_vuv, "data_x_pos_vuv/D"); 
  data_tree_vuv->Branch("data_y_pos_vuv", &data_y_pos_vuv, "data_y_pos_vuv/D");
  data_tree_vuv->Branch("data_event_vuv", &data_event_vuv, "data_event_vuv/I");

  // datatree for visible light
  data_tree_vis->Branch("data_time_vis", &data_time_vis, "data_time_vis/D");
  data_tree_vis->Branch("data_pmt_vis", &data_pmt_vis, "data_pmt_vis/I");
  data_tree_vis->Branch("data_x_pos_vis", &data_x_pos_vis, "data_x_pos_vis/D");
  data_tree_vis->Branch("data_y_pos_vis", &data_y_pos_vis, "data_y_pos_vis/D");
  data_tree_vis->Branch("data_event_vis", &data_event_vis, "data_event_vis/I");


  //event_tree is designed to save one entry for every event (ar39 decay).
  //hopefully by matching event numbers across ar39 decay and several photoelectrons
  //the energy of the event, how many photoelectrons were detected, and its
  //detector location can be compared.
  event_tree->Branch("event_no", &event_no, "event_no/I");
  event_tree->Branch("event_vox", &event_vox, "event_vox/I");
  event_tree->Branch("event_x_pos", &event_x_pos, "event_x_pos/D");
  event_tree->Branch("event_y_pos", &event_y_pos, "event_y_pos/D");
  event_tree->Branch("event_z_pos", &event_z_pos, "event_z_pos/D");
  event_tree->Branch("event_E", &event_E, "event_E/D");

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////-----------------lOADING THE LIBRARY-----------------------------///////////
  ////////////////////////////////////////////////////////////////////////////////////////

  if(config == 0) {libraryfile = "Lib154PMTs8inch_FullFoilsTPB.root"; }
  if(config == 1) {libraryfile = "Lib154PMTs8inch_OnlyCathodeTPB.root"; }
  if(config == 2) {libraryfile = "Lib154PMTs8inch_NoCathodeNoFoils.root"; }
  lar_light.LoadLibraryFromFile(libraryfile, reflected, reflT);


  ////////////////////////////////////////////////////////////////////////////////////
  ////////////-------------READING OUT PMTS FROM .TXT FILE-----------------///////////
  ////////////////////////////////////////////////////////////////////////////////////

  //note: the text file simply has a list of the PMT positions, and they are
  //simply being filled into a vector.

  ifstream myfile;
  myfile.open("posPMTs_setup1.txt");
  if(myfile.is_open()) {cout << "File opened successfully" << endl; }
  while(!myfile.eof())
    {
      double num_pmt, x_pmt, y_pmt, z_pmt;
      if(myfile >> num_pmt >> x_pmt >> y_pmt >> z_pmt)
	{
	  vector<double> line_data({num_pmt, x_pmt, y_pmt, z_pmt});
	  myfile_data.push_back(line_data);
	}
      else{break; }
    }
  myfile.close();
  // End Reading out positions of PMT from txt file.



  //////////////////////////////////////////////////////////////////////////////
  ////////////-------------------MAIN CODE---------------------------///////////
  //////////////////////////////////////////////////////////////////////////////
  
  gRandom->SetSeed(0);  

  energy_list.reserve(max_events);
  decay_time_list.reserve(max_events);
  voxel_list.reserve(max_events);

  for(int events = 0; events < max_events; events++)
    {

      // determine the energy of the event
      double energy;
      if(fixed_energy == true) {energy = fixedE;} 
      if(gen_argon == true) {energy = fSpectrum->GetRandom();} // pull from the Ar beta spectrum (see utility_functions.cc)
      if(supernova == true) {energy = flandau_sn->GetRandom();} // Pull from the predicted SN spectrum (see utility_functions.cc)
      if(gen_radon == true) {energy = fGauss->Gaus(Q_Rn, 0.05);}// Gaus(av,sigma) - is a ROOT function, pulls from a Gaussian
     

      // determine the position of the event 
      double position[3];
      int rand_voxel;

      if(random_pos == true) { // get a random voxel
	rand_voxel = gRandom->Uniform(319999);
	lar_light.GetVoxelCoords(rand_voxel, position);
      }

      else if(fixed_xpos == true){
	double randomY = ((rand() % 400) - 200)+0.5; // random Y voxel
	double randomZ = (rand() % 500) + 0.5; // random Z voxel		
      
	position[0] = fixedX; position[1]= randomY; position[2] = randomZ; // fill the vector
	rand_voxel = lar_light.GetVoxelID(position); // get the ID of the voxel
      }
      else { // fixed_pos == true
	position[0] = fixedX; position[1]= fixedY; position[2] = fixedZ;
	rand_voxel = lar_light.GetVoxelID(position);
      }

      voxel_list.push_back(rand_voxel); // push back the position of the voxel onto an array
      energy_list.push_back(energy); // push back the energy of the event onto the array
  
      //double decay_time = time_window * gRandom->Uniform(1.); // selects a random time of the event to decay within the specified time window
      double decay_time = 0; // all decays happen at t = 0
      decay_time_list.push_back(decay_time);


      //////////////////////////////
      ///// Fill the event tree ////
      //////////////////////////////
      event_no = events; // event number is equal to the current loop iteration
      event_vox = rand_voxel;
      event_x_pos = position[0];
      event_y_pos = position[1];
      event_z_pos = position[2];
      event_E = energy_list.at(events); // pull from the energy list at the current iteration
      event_tree->Fill();

    } // end of event loop 


  ////////*********************************************************************************************************************************

  //vector<vector<double> > timecut100ns(max_events, vector<double>(60));


  //Loop over each PMT for each event
  for(int events = 0; events < max_events; events++) {
    //By printing the event number here I can track the progress of the generation.
    cout << "Event: " << events + 1 << endl;

    //Begin looping over the PMT array. SBND plans to implement 60 PMTs,
    //with the list of numbers in the .h file (realisticPMT_IDs).
    for(int pmt_loop = 0; pmt_loop < 60; pmt_loop++)
      {

	int num_pmt = realisticPMT_IDs[pmt_loop]; 

	double x_pmt = myfile_data.at(num_pmt).at(1);
	double y_pmt = myfile_data.at(num_pmt).at(2);
	double z_pmt = myfile_data.at(num_pmt).at(3);


	//Funciton defined in library_access.cc
	// - This function will determine how many photons hit the given PMT
	vector<double> pmt_hits = lar_light.PhotonLibraryAnalyzer(energy_list.at(events), scint_yield, quantum_efficiency, num_pmt, voxel_list.at(events));
	int num_VUV = pmt_hits.at(0);
	int num_VIS = pmt_hits.at(1);


	//If no photons from this event for this PMT, go to the next event.
	if(num_VUV+num_VIS == 0) {continue; } // forces the next iteration

	// distance to pmt = delta_x^2 + delta_y^2 + delta_z^2
	double distance_to_pmt =
	  std::sqrt((pmt_hits.at(2)-x_pmt)*(pmt_hits.at(2)-x_pmt) +
	 	    (pmt_hits.at(3)-y_pmt)*(pmt_hits.at(3)-y_pmt) +
	 	    (pmt_hits.at(4)-z_pmt)*(pmt_hits.at(4)-z_pmt));

	data_x_pos = pmt_hits.at(2);
	data_x_pos_vuv = pmt_hits.at(2);
	data_x_pos_vis = pmt_hits.at(2);

	data_y_pos = pmt_hits.at(3);
	data_y_pos_vuv = pmt_hits.at(3);
	data_y_pos_vis = pmt_hits.at(3);

	data_z_pos = pmt_hits.at(4);


	data_event = events;
	data_event_vuv = events;
	data_event_vis = events;

	
	///*************************** ///
	///***** TIMING OPERATIONS *** ///
	///*************************** ///
	// NOTE:  Distances in cm and times in ns, returns a vector of doubles for VUV arrival times OF EACH PHOTON
	
	///////////////////////
	// Get the VUV time
	//////////////////////

	vector<double> transport_time_vuv;
	if(num_VUV != 0) {transport_time_vuv = utility::GetVUVTime(distance_to_pmt, num_VUV);}


	//This statement is to prvent issues when the parameterisation is not well defined
	if(num_VUV != transport_time_vuv.size()) // is every VUV photon accounted for?
	  {
	    cout << "Param fail" << endl;
	    num_VUV = transport_time_vuv.size();
	  }


	double total_time_vuv; // total time for each VUV photon to hit a pmt from creation
	for(auto& x: transport_time_vuv)
	  {
	    //data_time has several parts:
	    // x is the output from the parameterisation -> the transport time (in NANOSECONDS) x * 0.001 converts from ns -> micros_s
	    // the scintillation function timing is also converted to microseconds
	    // the time window offset is when the decay occured, given already in microseconds

	    total_time_vuv = (x*0.001+(decay_time_list.at(events) + fScintillation_function->GetRandom())*1000000.); // in microseconds 

	    //////////////////////////100ns CUT//////////////////////////////////
	    if(total_time_vuv > time_cut && cut == true){ // 0.1 microseconds = 100 ns! 
	      continue; // continue to the next iteration without filling - don't bother filling it
	    }


	    data_time = total_time_vuv;
	    data_time_vuv = total_time_vuv;

	    data_pmt = num_pmt;
	    data_pmt_vuv = num_pmt;

	    data_tree_vuv->Fill();
	    data_tree->Fill();
	  }

	transport_time_vuv.clear();

	///////////////////////////////
	//calculating the visible light arrival times 
	///////////////////////////////

	vector<double> transport_time_vis;
	if(num_VIS != 0 && config == 1)
	  {
	    transport_time_vis = utility::GetVisibleTimeOnlyCathode(pmt_hits.at(5), num_VIS);
	    double total_time_vis;
	    for(auto &y : transport_time_vis)
	      {
		total_time_vis = (y*0.001+(decay_time_list.at(events) + fScintillation_function->GetRandom())*1000000.); // in microseconds

		data_time = total_time_vis; 
		data_time_vis = total_time_vis;

		if(total_time_vis > time_cut && cut == true){ // 0.1 microseconds = 100 ns! 
		  continue; // go onto the next interation - cut has been made
		}
	 
		data_pmt = num_pmt;
		data_pmt_vis = num_pmt;
		data_tree_vis->Fill();
		data_tree->Fill();
	      }
	    transport_time_vis.clear();
	  }
     

      }//end pmt loop



  }//end loop over events


  energy_list.clear();
  decay_time_list.clear();
  voxel_list.clear();

  //write the files
  data_file.Write();
  event_file.Write();

  return 0;

}//end main
