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
	if(gen_argon == true) {cout << "GEN Argon" << endl; }
	if(supernova == true) {cout << "GEN SN" << endl; }
	if(gen_argon == true) {cout << "Ar39 events: " << max_events << " in: " << time_window << " s" << endl; }
	if(supernova == true) {cout << "SN events: " << max_events_sn << endl; }

	if(config == 0) {cout << "Full Foils" << endl; }
	if(config == 1) {cout << "Cath Foils" << endl; }
	if(config == 2) {cout << "VUV only" << endl; }

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

	gRandom->SetSeed(0);

	TF1 *fSpectrum = new TF1("fSpectrum",utility::SpectrumFunction,0,Q,1);//-----Beta decay spectrum
	TF1 *flandau_sn = new TF1("flandau_sn",utility::fsn, 0, 50, 1);//--SN Nu spectrum
	fSpectrum->SetParameter(0, Q);
	flandau_sn->SetParameter(0, Eav);
	TF1 *fScintillation_function = new TF1("Ar Scintillation Timing", utility::Scintillation_function, 0, scint_time_window, 2);
	fScintillation_function->SetParameter(0, t_singlet);
	fScintillation_function->SetParameter(1, t_triplet);


	//*********** loading the library ****************
	if(config == 0) {libraryfile = "Lib154PMTs8inch_FullFoilsTPB.root"; }
  if(config == 1) {libraryfile = "Lib154PMTs8inch_OnlyCathodeTPB.root"; }
	if(config == 2) {libraryfile = "Lib154PMTs8inch_NoCathodeNoFoils.root"; }
	lar_light.LoadLibraryFromFile(libraryfile, reflected, reflT);

	// Reading out positions of PMTs from txt file *****
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
	// End Reading out positions of PMT from txt file
  //note: the text file simply has a list of the PMT positions, and they are
  //simply being filled into a vector.


  //data_tree preceeded event_tree and was designed such that each photoelectron
  //would have their own entry in the ttree with a time, which pmt, xyz, and
  //event number (data_event).
  //A current issue with the event number tracking is that it does not carry
  //over between several iterations of the code. For single uses this is fine,
  //but memory consuption and speed for running this code a single time
  //becomes a limiting factor, such that running many instances of this code
  //in parallel becomes much more useful. Therefore carrying over the event
  //number between iterations may be useful.
	data_tree->Branch("data_time", &data_time, "data_time/D");
	data_tree->Branch("data_pmt", &data_pmt, "data_pmt/I");
	data_tree->Branch("data_x_pos", &data_x_pos, "data_x_pos/D");
	data_tree->Branch("data_y_pos", &data_y_pos, "data_y_pos/D");
	data_tree->Branch("data_z_pos", &data_z_pos, "data_z_pos/D");

	data_tree->Branch("data_event", &data_event, "data_event/I");

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


//-----------------Ar39 Loop-----------------------------
	if(gen_argon == true)
	{
		energy_ar_list.reserve(max_events);
		scint_time_list.reserve(max_events);
		voxel_list.reserve(max_events);

		for(int events = 0; events < max_events; events++)
		{
			energy_ar_list.push_back(fSpectrum->GetRandom());
			double decay_time = time_window * gRandom->Uniform(1.);
			scint_time_list.push_back(decay_time);
			int rand_voxel = gRandom->Uniform(319999);
			voxel_list.push_back(rand_voxel);//select random voxel -> random position in detector

      //The goal with simulating the ar39 is to be close to the physical system
      //this means randomly creating many of its characteristics.
      //However, for looking specifically at the detector response in certain
      //areas, using a predefined set of energy and voxel can be useful.
      //Some preselected voxels are:
            // 161660 : (105, 7.5, 252.5)
						// 166575 (80, -177.5, 262.5)
						// 240393 (170, -152.5, 377.5)

			double position[3];
			lar_light.GetVoxelCoords(rand_voxel, position);
			event_no = events;
			event_vox = rand_voxel;
			event_x_pos = position[0];
			event_y_pos = position[1];
			event_z_pos = position[2];
			event_E = energy_ar_list.at(events);
			event_tree->Fill();

		}

    //Begin looping over the PMT array. SBND plans to implement 60 PMTs,
    //with the list of numbers in the .h file (realisticPMT_IDs).
		for(int pmt_loop = 0; pmt_loop < 60; pmt_loop++) {
			int num_pmt = realisticPMT_IDs[pmt_loop];

      //By printing the PMT number here I can track the progress of the generation.
			cout << "PMT: " << num_pmt << endl;

			double x_pmt = myfile_data.at(num_pmt).at(1);
			double y_pmt = myfile_data.at(num_pmt).at(2);
			double z_pmt = myfile_data.at(num_pmt).at(3);

			TVector3 optdet (x_pmt/100., y_pmt/100., z_pmt/100.);

      //Loop over all events for each PMT
			for(int events = 0; events < max_events; events++) {
        //Funciton defined in library_access.cc
				vector<double> pmt_hits = lar_light.PhotonLibraryAnalyzer(energy_ar_list.at(events), scint_yield, quantum_efficiency, num_pmt, voxel_list.at(events));
				int num_VUV = pmt_hits.at(0);
				int num_VIS = pmt_hits.at(1);

        //If no photoelectrons from this event for this PMT, go to the next event.
				if(num_VUV+num_VIS == 0) {continue; }

				double distance_to_pmt =
				        std::sqrt((pmt_hits.at(2)-x_pmt)*(pmt_hits.at(2)-x_pmt) +
				                  (pmt_hits.at(3)-y_pmt)*(pmt_hits.at(3)-y_pmt) +
				                  (pmt_hits.at(4)-z_pmt)*(pmt_hits.at(4)-z_pmt));

				data_x_pos = pmt_hits.at(2);
				data_y_pos = pmt_hits.at(3);
				data_z_pos = pmt_hits.at(4);
				data_event = events;

				///*************************** ///
				///***** TIMING OPERATIONS *** ///
				///*************************** ///

				vector<double> time_vuv;
				if(num_VUV != 0) { time_vuv = utility::GetVUVTime(distance_to_pmt, num_VUV); }

        //This statement is to prvent issues when the parameterisation is not
        //well defined.
				if(num_VUV != time_vuv.size())
				{
					cout << "Param fail" << endl;
					num_VUV = time_vuv.size();
				}

				for(auto& x: time_vuv)
				{
          //data_time has several parts:
          // x is the output from the parameterisation, and it's converted to microseconds
          // the scintillation function timing is also converted to microseconds
          // the time window offset is when the decay occured, given already in microseconds
					data_time = (x*0.001+(scint_time_list.at(events) + fScintillation_function->GetRandom())*1000000.);
					data_pmt = num_pmt;
					data_tree->Fill();
				}

				time_vuv.clear();

        //As seen in utility functions, the parameterisation for the case where
        //we have foils on the cathode and the case where we have foils on the
        //the whole TPC (not behind the PMTs), are different...

				//This is for foils only on cathode
				vector<double> time_vis;
				if(num_VIS != 0 && config == 1)
				{
					time_vis = utility::GetVisibleTimeOnlyCathode(pmt_hits.at(5), num_VIS);
					for(auto &y : time_vis)
					{
						data_time = (y*0.001+(scint_time_list.at(events) + fScintillation_function->GetRandom())*1000000.);
						data_pmt = num_pmt;
						data_tree->Fill();
					}
					time_vis.clear();
				}

        /// For full foil configuration only
				/// This function expects distances in m
				if(num_VIS != 0 && config == 0)
				{
					vector<double> time_vis;
					TVector3 scint_point (pmt_hits.at(2)/100., pmt_hits.at(3)/100., pmt_hits.at(4)/100.);
					double tmean = utility::TimingParamReflected(scint_point, optdet);

					double divide = 0.43*tmean + 7;
					double vuv_vgroup = 10.13; //cm/ns
					double t_direct = distance_to_pmt/vuv_vgroup;

					if(pmt_hits.at(5) < divide && !(pmt_hits.at(5) < divide && pmt_hits.at(5) > 20.))
					{
						time_vis = utility::GetVisibleTimeFullConfig1(pmt_hits.at(5), tmean, distance_to_pmt, num_VIS);
					}
					else{ time_vis = utility::GetVisibleTimeFullConfig2(pmt_hits.at(5), tmean, distance_to_pmt, num_VIS); }

					for(auto &y : time_vis)
					{
						data_time = (y*0.001+(scint_time_list.at(events) + fScintillation_function->GetRandom())*1000000.);
						data_pmt = num_pmt;
						data_tree->Fill();
					}
					time_vis.clear();
				}
				pmt_hits.clear();
			}//end looping events
		}//end loop over pmts

		data_tree->Fill();
	}//end if gen_argon

  /* Ion generation is another topic: this code was writen based on some
  tests for ions recombining on the cathode and scintillating. The numbers
  provided here are from a separate simulation which makes many assumptions!
  I would recommend in general not using this setting.
  */

	if(config == 0 && gen_ions == true)
	{
		///For full foils configuration!
		int total_extra_vis = 292416211 * time_window;
		for(int extra_vis = 0; extra_vis < total_extra_vis; extra_vis++)
		{
			double rand_time = gRandom->Uniform(time_window)*1000000.;
			int rand_num = gRandom->Uniform(60);
			int rand_pmt = realisticPMT_IDs[rand_num];

			data_time = rand_time;
			data_pmt = rand_pmt;
			data_tree->Fill();
		}
	}
	if(config == 1 && gen_ions == true)
	{
		///For cathode only foils!
		int total_extra_vis = 143210695 * time_window;
		for(int extra_vis = 0; extra_vis < total_extra_vis; extra_vis++)
		{
			double rand_time = gRandom->Uniform(time_window)*1000000.;
			int rand_num = gRandom->Uniform(60);
			int rand_pmt = realisticPMT_IDs[rand_num];
			data_time = rand_time;
			data_pmt = rand_pmt;
			data_tree->Fill();
		}
	}
	if(gen_ions == true)
	{
		///VUV - for all ion cases
		int total_extra_vuv = 80003251 * time_window;
		for(int extra_vuv = 0; extra_vuv < total_extra_vuv; extra_vuv++)
		{
			double rand_time = gRandom->Uniform(time_window)*1000000.;
			int rand_num = gRandom->Uniform(60);
			int rand_pmt = realisticPMT_IDs[rand_num];
			//total_timing.push_back(make_pair(rand_time, rand_pmt));
			data_time = rand_time;
			data_pmt = rand_pmt;
			data_tree->Fill();
		}
	}

	energy_ar_list.clear();
	scint_time_list.clear();
	voxel_list.clear();

  //Left-over plotting code to make a simple histogram to show the timing.
/*
   TCanvas *can_time2 = new TCanvas();
   can_time2->cd();
   h_single->GetYaxis()->SetTitle("Number of Photoelectrons / 60 PMTs");
   h_single->GetXaxis()->SetTitle("Time [#mu s]");
   h_single->Draw();
   can_time2->Print("pmt_time.pdf");
 */

	//------End Ar39 loop------------------
	//---------------------------------------------------
	//---------------------------------------------------
	//-----Supernova loop-------------------
	if(supernova == true)
	{

    //I want to treat each supernova event individually, as we expect few
    //events (~10) on a time-scale much larger than the longest scintillation
    //timing (~8 us or so), and likely won't overlap at all.
		for(int events_sn = 0; events_sn < max_events_sn; events_sn++)
		{
			double scint_time_list = time_window * gRandom->Uniform(1.);
			vector<double> pmt_hits_sn;
			double energy_sn = flandau_sn->GetRandom();
			int rand_voxel = gRandom->Uniform(319999);//select random voxel -> random position in detector

			//Loop over all the PMTs
			for(int loop_pmt = 0; loop_pmt < 60; loop_pmt++)
			{
				int num_pmt = realisticPMT_IDs[loop_pmt];

				pmt_hits_sn = lar_light.PhotonLibraryAnalyzer(energy_sn, scint_yield, quantum_efficiency, num_pmt, rand_voxel);

				double x_pmt = myfile_data.at(num_pmt).at(1);
				double y_pmt = myfile_data.at(num_pmt).at(2);
				double z_pmt = myfile_data.at(num_pmt).at(3);

				TVector3 optdet (x_pmt/100., y_pmt/100., z_pmt/100.);

				int num_VUV = pmt_hits_sn.at(0);
				int num_VIS = pmt_hits_sn.at(1);

				if(num_VUV+num_VIS == 0) {continue; }

				double distance_to_pmt =
				        std::sqrt((pmt_hits_sn.at(2)-x_pmt)*(pmt_hits_sn.at(2)-x_pmt) +
				                  (pmt_hits_sn.at(3)-y_pmt)*(pmt_hits_sn.at(3)-y_pmt) +
				                  (pmt_hits_sn.at(4)-z_pmt)*(pmt_hits_sn.at(4)-z_pmt));

				data_x_pos = pmt_hits_sn.at(2);
				data_y_pos = pmt_hits_sn.at(3);
				data_z_pos = pmt_hits_sn.at(4);


				///*************************** ///
				///***** TIMING OPERATIONS *** ///
				///*************************** ///

				vector<double> time_vuv;
				if(num_VUV != 0) { time_vuv = utility::GetVUVTime(distance_to_pmt, num_VUV); }

				if(num_VUV != time_vuv.size())
				{
					cout << "Param fail" << endl;
					num_VUV = time_vuv.size();
				}

				for(auto& x: time_vuv)
				{
					data_time = (x*0.001+(scint_time_list + fScintillation_function->GetRandom())*1000000.);
					data_pmt = num_pmt;
					data_event = events_sn;
					data_tree->Fill();
				}
				time_vuv.clear();


				// this is for foils only on cathode
				vector<double> time_vis;
				if(num_VIS != 0 && config == 1) {
					time_vis = utility::GetVisibleTimeOnlyCathode(pmt_hits_sn.at(5), num_VIS);
					for(auto &y : time_vis)
					{
						data_time = (y*0.001+(scint_time_list + fScintillation_function->GetRandom())*1000000.);
						data_pmt = num_pmt;
						data_event = events_sn;
						data_tree->Fill();
					}
					time_vis.clear();
				}

				/// For full foil configuration only
				/// This function expects distances in m
				if(num_VIS != 0 && config == 0)
				{
					vector<double> time_vis;
					TVector3 scint_point (pmt_hits_sn.at(2)/100., pmt_hits_sn.at(3)/100., pmt_hits_sn.at(4)/100.);
					double tmean = utility::TimingParamReflected(scint_point, optdet);


					double divide = 0.43*tmean + 7;
					double vuv_vgroup = 10.13;//cm/ns
					double t_direct = distance_to_pmt/vuv_vgroup;

					if(pmt_hits_sn.at(5) < divide && !(pmt_hits_sn.at(5) < divide && pmt_hits_sn.at(5) > 20.))
					{
						time_vis = utility::GetVisibleTimeFullConfig1(pmt_hits_sn.at(5), tmean, distance_to_pmt, num_VIS);
					}
					else{ time_vis = utility::GetVisibleTimeFullConfig2(pmt_hits_sn.at(5), tmean, distance_to_pmt, num_VIS); }

					/// Loop over number of photons generated
					for(auto &y : time_vis)
					{
						data_time = (y*0.001+(scint_time_list + fScintillation_function->GetRandom())*1000000.);
						data_pmt = num_pmt;
						data_event = events_sn;
						data_tree->Fill();
					}
					time_vis.clear();
				}

			}//End looping pmts



			if(config == 0 && gen_ions == true)
			{
				///For full foils configuration!
				int total_extra_vis = 292416211 * time_window;
				for(int extra_vis = 0; extra_vis < total_extra_vis; extra_vis++)
				{
					double rand_time = gRandom->Uniform(time_window)*1000000.;
					int rand_num = gRandom->Uniform(60);
					int rand_pmt = realisticPMT_IDs[rand_num];
					//total_timing.push_back(make_pair(rand_time, rand_pmt));
					data_time = rand_time;
					data_pmt = rand_pmt;
					data_tree->Fill();
				}
			}
			if(config == 1 && gen_ions == true)
			{
				///For cathode only foils!
				int total_extra_vis = 143210695 * time_window;
				for(int extra_vis = 0; extra_vis < total_extra_vis; extra_vis++)
				{
					double rand_time = gRandom->Uniform(time_window)*1000000.;
					int rand_num = gRandom->Uniform(60);
					int rand_pmt = realisticPMT_IDs[rand_num];
					//total_timing.push_back(make_pair(rand_time, rand_pmt));
					data_time = rand_time;
					data_pmt = rand_pmt;
					data_tree->Fill();
				}
			}
			if(gen_ions == true)
			{
				///For no foils! - VUV
				int total_extra_vuv = 80003251 * time_window;
				for(int extra_vuv = 0; extra_vuv < total_extra_vuv; extra_vuv++)
				{
					double rand_time = gRandom->Uniform(time_window)*1000000.;
					int rand_num = gRandom->Uniform(60);
					int rand_pmt = realisticPMT_IDs[rand_num];
					//total_timing.push_back(make_pair(rand_time, rand_pmt));
					data_time = rand_time;
					data_pmt = rand_pmt;
					data_tree->Fill();
				}
			}

		}//End loop events
	}//End if SN true
//-----------End Supernova loop--------------------------------------------
//--------------------------------------------------------------------------

	data_file.Write();
	event_file.Write();

	return 0;

}//end main
