//
// small program to read a Tree with physical events and write it
// into several (binned) files
// the files will be copied into folders named <binlow>.<binhigh>
// and must exist before executing this program!
// author P.Jasinski Promme@web.de jasinski@kph.uni-mainz.de
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <map>
#include <utility>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;

struct TParticle{
		int geantid;
		int charge;
		double px;
		double py;
		double pz;
		double E;
};

const int PIDGEANTTOCHARGE[16]={0,0,+1,-1,0,+1,-1,0,+1,-1,0,+1,-1,0,+1,-1};

// writes one event into stream
void WriteEventToStream(vector<TParticle>& particles, ofstream& stream);

double last_percent(-1.);
// draw a progress bar only when the length changes significantly
void DrawProgressBar(int len, double percent);

// input: root file name, tree name in the file, path to write to, nbins, binlow, binhigh [MeV], N particles
int main(int argc, char* argv[]){
	if (argc != 8){
		cout << " please check options" << endl;
		return 0;
	}
	int nbins = atoi(argv[4]);
	float binlow  = (float) atoi(argv[5]);
	float binhigh = (float) atoi(argv[6]);
	float binwidth = (binhigh-binlow)/nbins;
	int N_expected_events = atoi(argv[7]);
	string path = argv[3];
	TH1F checkhist("checkhist","invariant mass distribution", nbins, binlow, binhigh);
	// map of files that will be filled with streams to the text files	
	map<int, ofstream*> files; // ofstream must be stored as a pointer (try out what happens if not)
	// create the filestreams
	for (int i = (int) binlow; i < (int) binhigh; i+=(int)binwidth) {
		stringstream _filename;
		stringstream _filekey;
		_filekey << i << "." << i+(int)binwidth;
		_filename << path << _filekey.str() << "/" << _filekey.str() << ".evt";	
		cout << " opening file " << _filename.str() << endl;	
		// use the upper bound for the key		
		files[i+(int)binwidth] = new ofstream(_filename.str().c_str());
	}
	// open the root file
	TFile rootfile(argv[1]);
	// try to get the tree
	TTree* tree_events = (TTree*) rootfile.Get(argv[2]); 
	if (tree_events){
		cout << " found " << tree_events->GetEntries() << " events " << endl;
		// set variables to the leaves
		double	events_px;
		double	events_py;
		double	events_pz;
		double	events_E;
		int 	events_GeantPID;
		int		events_id;
		tree_events->SetBranchAddress("event_id", &events_id);//, "id/I");
		tree_events->SetBranchAddress("geantPID", &events_GeantPID);
		tree_events->SetBranchAddress("px", &events_px);
		tree_events->SetBranchAddress("py", &events_py);
		tree_events->SetBranchAddress("pz", &events_pz);
		tree_events->SetBranchAddress("E" , &events_E);
		// read the events and write the events to the file
		// assuming to have the correct order and number of events
		// but checking the first (beam) particle to stay the same
		int last_beampart_id(-1);
		int particle_counter(0);
		int events_id_cut(-1); // put a -1 if you want to filter the first or only event id in the tree
		// we need the datastructure to hold the particle information
		vector<TParticle> particles;
		int nentries = tree_events->GetEntries();
		for (int i = 0; i < nentries; i++){
			tree_events->GetEntry(i);
			if (events_id_cut < 0){ // take the first events_id to cut on
				events_id_cut = events_id;
				cout << " filtering events with events_id " << events_id << endl;
			} 
			DrawProgressBar(50, (i+1)/(float)nentries);
			if (events_id_cut != events_id){
				continue;
			}
			//cout << "here " << i << " particles " << particle_counter  << endl;
			// initilize the last_beampart_id if it appears to be the first entry found
			if (last_beampart_id < 0){
				last_beampart_id = events_GeantPID;
			}
			// count the number of particles of this event id
			particle_counter++;
			if (particle_counter > N_expected_events){
				// must be the first particle namely the beam particle again
				if (last_beampart_id != events_GeantPID){
					cout << " Error: structure does not recur! " << endl;
					break;
				}
				// reconstruct the invariant mass of the outgoing particle system
				// to get the bin entry
				double _px(0), _py(0), _pz(0), _E(0);
				for (int ipart = 1; ipart < N_expected_events; ipart++){
					_px += particles[ipart].px;
					_py += particles[ipart].py;
					_pz += particles[ipart].pz;
					_E  += particles[ipart].E;
				}
				int _M = (int) (sqrt(_E*_E - _px*_px - _py*_py - _pz*_pz) * 1000); // (MeV)
				// cut into the allowed range
				if ((binlow < _M) && (_M < binhigh)){
					// write all particles to the file
					WriteEventToStream(particles, *(*files.upper_bound(_M)).second);
					checkhist.Fill(_M);				
				}
				// clear the vector holding the particle information
				particles.clear();
				last_beampart_id = events_GeantPID;
				particle_counter = 1;
			}
			// store the new information available
			TParticle particle;
			particle.E = events_E;
			particle.geantid = events_GeantPID;
			particle.px 	 = events_px;
			particle.py 	 = events_py;
			particle.pz 	 = events_pz;
			if (events_GeantPID >= 16) {
				cout << " please check PIDGEANTTOCHARGE array " << endl;
				particle.charge = 0;
			} else {
				particle.charge = PIDGEANTTOCHARGE[events_GeantPID];
			}
			particles.push_back(particle);
		}
		// note: the last event sample is lost since I do not write it out
		// for a few thousand of them this might be negligible
		// close all files
	}
	// close all files
	for (int i = (int) binlow; i < (int) binhigh; i+=(int)binwidth) {		
		files[i+(int)binwidth]->close();
		delete files[i+(int)binwidth];
	}
	// Draw a crosscheck histogram
	TCanvas canvas("canvas", "canvas", 1000, 600);
	checkhist.Draw("P E Text");
	canvas.Print( (path+"checkhist.pdf").c_str());
	// finish the job
	rootfile.Close();
	return 0;
} 

void WriteEventToStream(vector<TParticle>& particles, ofstream& stream){
	//Final state particle 4-momenta in BNL format:
	//<Total number of particles in event (including beam)>
	//<geant3ID of beam><charge><px><py><pz><E>
	//<geant3ID of final state particle 1><charge><px><py><pz><E>
	//<geant3ID of final state particle 2><charge><px><py><pz><E>
	//<geant3ID of final state particle 3><charge><px><py><pz><E>
	// ...

	// set precision to a high value
	stream.precision(20);
	// writing into stream in BNL format
	stream << particles.size() << endl;
	for (unsigned int i = 0; i < particles.size(); i++){
		// use scientific output format to keep the maximum information
		stream << scientific << particles[i].geantid << " " << particles[i].charge << " "  << particles[i].px << " "  << particles[i].py << " "  << particles[i].pz << " "  << particles[i].E << endl;
	}
}

void DrawProgressBar(int len, double percent) {
  if ((int)(last_percent*100) == (int)(percent*100)) return;
  //cout << " drawing " << endl;
  cout << "\x1B[2K"; // Erase the entire current line.
  cout << "\x1B[0E"; // Move to the beginning of the current line.
  string progress;
  for (int i = 0; i < len; ++i) {
    if (i < static_cast<int>(len * percent)) {
      progress += "=";
    } else {
      progress += " ";
    }
  }
  cout << "[" << progress << "] " << (static_cast<int>(100 * percent)) << "%";
  flush(cout); // Required.
  last_percent = percent;
}

