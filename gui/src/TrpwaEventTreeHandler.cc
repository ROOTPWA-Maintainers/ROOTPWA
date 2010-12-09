/*
 * TrpwaEventTreeHandler.cc
 *
 *  Created on: Nov 11, 2010
 *      Author: Promme
 */

#include "TrpwaEventTreeHandler.h"
#include <dirent.h>
#include <cerrno>
#include <sys/types.h>
#include <sys/stat.h>

const int PIDGEANTTOCHARGE[16]={0,0,+1,-1,0,+1,-1,0,+1,-1,0,+1,-1,0,+1,-1};

TrpwaEventTreeHandler::TrpwaEventTreeHandler(){
	last_percent = -1.;
	_hist_data = NULL;
	_hist_data_mc = NULL;
	_hist_data_mc_acc = NULL;
	//_nentries_data.clear();
	//_nentries_data_mc.clear();
	//_nentires_data_mc_acc.clear();
}

bool TrpwaEventTreeHandler::Add_eventtreefile(string filename, string treename){
	bool result = true;
	if (!FileExists(filename)){
		cout << " Error in TrpwaEventTreeHandler::Add_eventtreefile(): File "<< filename <<" does not exist!" << endl;
		return false;
	}
	TFile* _file = new TFile(filename.c_str());
	if (_file->IsZombie()){
		cout << " Error in TrpwaEventTreeHandler::Add_eventtreefile(): File "<< filename <<" not accessible!" << endl;
		_file->Close();
		delete _file;
		return false;
	}
	if (!_file->Get(treename.c_str())){
		cout << " Error in TrpwaEventTreeHandler::Add_eventtreefile(): Tree " << treename << " not found in " << filename << endl;
		_file->Close();
		delete _file;
		return false;
	}
	_eventtreefiles.push_back(_file);
	_eventtreenames.push_back(treename);
	return result;
}



bool TrpwaEventTreeHandler::Add_eventtreefiles(vector<string> filenames, string treename){
	bool result = true;
	_eventtreefiles.clear();
	_eventtreenames.clear();
	for (vector<string>::iterator it = filenames.begin(); it != filenames.end(); it++){
		if (!Add_eventtreefile(*it, treename)){
			result = false;
			_eventtreefiles.clear();
			_eventtreenames.clear();
			return result;
		}
	}
	return result;
}



bool TrpwaEventTreeHandler::Write_Trees_to_BNL_events(bool overwrite){
	delete _hist_data;
	delete _hist_data_mc;
	delete _hist_data_mc_acc;
	bool result = true;
	cout << " Info: filtering events into bins " << endl;
	// map of files that will be filled with streams to the text files
	map<int, ofstream**> files; // ofstream must be stored as a pointer (try out what happens if not)
	// create the file streams and save upper and lower bounds
	// ensure in addition that the range is covered
	int binlow = -1;
	int binhigh= -1;
	int binhigh_previous = -1;
	float* xbins = new float[_bins.size()+1];
	int nbins(0);
	for (map<int, int>::iterator it = _bins.begin(); it != _bins.end(); it++){
		if (binlow == -1 || it->first < binlow ) binlow = it->first;
		if (binhigh== -1 || it->second> binhigh) binhigh= it->second;
		// check a closed range for binning
		if (binhigh_previous != -1){
			if (it->first != binhigh_previous){
				result = false;
				cout << " Error in TrpwaEventTreeHandler::Write_Trees_to_BNL_events(): Not all expected bins are given! " << endl;
				cout << "    Bin " << binhigh_previous << " to " << it->first << " is missing! " << endl;
			}
		}
		binhigh_previous = it->second;
		ofstream** _files = new ofstream*[3];
		_files[0] = NULL; _files[1] = NULL; _files[2] = NULL;
		for (int j = 0; j < 3; j++){ // real, MC, MC accepted
			stringstream _filename;
			stringstream _filekey;
			_filekey << it->first << "." << it->second;
			_filename << _dir << "/" << _bin_paths[it->first] << "/" << _filekey.str();
			switch (j) {
			case 0:
				break;// do nothing
			case 1:
				_filename << ".genbod";
				break;
			case 2:
				_filename << ".acc";
				break;
			}
			_filename << ".evt";
			cout << " opening file " << _filename.str() << endl;
			_files[j] = new ofstream(_filename.str().c_str());
		}
		// use the upper bound for the key
		files[it->second] = _files;
		// save the lower bound and in addition the last bin must have the upper bound
		xbins[nbins] = it->first;
		nbins++;
		xbins[nbins] = it->second; // I know that only the last bin will be not overwritten
	}
	// create some crosscheck hists
	TH1I* checkhist[5];
	checkhist[0] = new TH1I("checkhist0","invariant mass distribution of reconstructed events", nbins, xbins);
	checkhist[1] = new TH1I("checkhist1","invariant mass distribution of MC exact events", nbins, xbins);
	checkhist[2] = new TH1I("checkhist2","invariant mass distribution of MC reconstructed events", nbins, xbins);
	checkhist[3] = new TH1I("checkhist3","invariant mass distribution of MC truth reconstructed events", nbins, xbins);
	checkhist[4] = new TH1I("checkhist4","invariant mass distribution of MC false reconstructed events", nbins, xbins);
	checkhist[0]->SetDirectory(0); // to prevent root deleting this objects when closing files following
	checkhist[1]->SetDirectory(0);
	checkhist[2]->SetDirectory(0);
	checkhist[3]->SetDirectory(0);
	checkhist[4]->SetDirectory(0);

	for (unsigned int i = 0; i < _eventtreefiles.size(); i++){
		TFile& rootfile = *_eventtreefiles[i];
		TTree* tree_events = (TTree*) rootfile.Get(_eventtreenames[i].c_str());
		if (!tree_events){
			cout << "\n Sorry, tree " << _eventtreenames[i] << " not found! " << endl;
			result = false;
		}
		// in addition we do have 3 output streams for real data, MC generated and MC accepted events
		int openfile[3] = {0,0,0}; // variable to store whether the file contains real data, MC generated data or MC accepted events
		// search for entries in the tree
		cout << "\n scanning events in tree # " << i+1 << " of " << _eventtreefiles.size() << endl;
		openfile[0] = tree_events->GetEntries("isMC == 0");
		cout << " found " << openfile[0] << " reconstructed events " << endl;
		openfile[1] = tree_events->GetEntries("isMC > 0");
		cout << " found " << openfile[1] << " MC generated events " << endl;
		openfile[2] = tree_events->GetEntries("isMC == 2");
		cout << " found " << openfile[2] << " MC accepted events " << endl;
		int false_events = tree_events->GetEntries("isMC == -1");
		cout << " found " << false_events << " MC wrongly reconstructed events " << endl;

		// set variables to the leaves
		TClonesArray* events_lzvecs = NULL; // lz vectors of 0..X-1 particles out, X particle in, X+1 recoil particle
		vector <int>* events_charges = NULL;
		vector <int>* events_g3pids = NULL;
		TClonesArray* events_lzvecs_mc = NULL;
		int		events_id;
		int		events_isMC; // 0 not a MC event; 1 MC event; 2 accepted MC event; -1 MC but falsely reconstructed (since tree version 3)
		tree_events->SetBranchAddress("id"  , &events_id);
		tree_events->SetBranchAddress("p"   , &events_lzvecs);
		tree_events->SetBranchAddress("p_mc", &events_lzvecs_mc);
		tree_events->SetBranchAddress("q"   , &events_charges);
		tree_events->SetBranchAddress("pid" , &events_g3pids);
		tree_events->SetBranchAddress("isMC", &events_isMC);
		// read the events and write the events to the file
		int events_id_cut(-1); // put a -1 if you want to filter the first or only event id in the tree
		// we need the data structure to hold the particle information
		vector<TParticle> particles;
		vector<TParticle> particles_mc;
		int nentries = tree_events->GetEntries();
		cout << " running over " << nentries << " events " << endl;
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
			particles.clear();
			particles_mc.clear();
			TLorentzVector sumlz(0.,0.,0.,0.);
			TLorentzVector sumlz_mc(0.,0.,0.,0.);
			// loop through the particles
			// last two particles are the beam particle and the (neglected) recoil proton
			// the loop starts with -1 which is a special case to retrieve the
			// beam particle first that is not the first element in events_lzvecs
			// one could of course use also vector::insert but this would cost performance
			for (int i = -1; i < events_lzvecs->GetEntriesFast()-2; i++){
				int ipart = i;
				if (i < 0) ipart = events_lzvecs->GetEntriesFast()-2; // case recoil proton
				TLorentzVector* lzvec = (TLorentzVector*)(*events_lzvecs)[ipart];
				// store the new information available
				TParticle particle;
				particle.geantid = (*events_g3pids)[ipart];
				particle.charge  = (*events_charges)[ipart];
				particle.E       = lzvec->E();
				particle.px      = lzvec->X();
				particle.py      = lzvec->Y();
				particle.pz      = lzvec->Z();
				// reconstruct the invariant mass of the outgoing particle system
				// to get the bin entry
				if (i > -1) sumlz += *lzvec;
				particles.push_back(particle);
			}
			for (int i = -1; i < events_lzvecs_mc->GetEntriesFast()-2; i++){
				int ipart = i;
				if (i < 0) ipart = events_lzvecs_mc->GetEntriesFast()-2; // case recoil proton
				TLorentzVector* lzvec = (TLorentzVector*)(*events_lzvecs_mc)[ipart];
				// store the new information available
				TParticle particle;
				particle.geantid = (*events_g3pids)[ipart];
				particle.charge  = (*events_charges)[ipart];
				particle.E       = lzvec->E();
				particle.px      = lzvec->X();
				particle.py      = lzvec->Y();
				particle.pz      = lzvec->Z();
				// reconstruct the invariant mass of the outgoing particle system
				// to get the bin entry
				if (i > -1) sumlz_mc += *lzvec;
				particles_mc.push_back(particle);
			}
			int _M    = (int) floor(sumlz   .M() * 1000.); // (MeV)
			int _M_mc = (int) floor(sumlz_mc.M() * 1000.); // (MeV)
			// cut into the allowed range
			//if ((binlow <= _M) && (_M < binhigh)){
				// write all particles to the file
				if (events_isMC < -1 || events_isMC > 2){
					cout << " error: not defined MC status: " << events_isMC << endl;
					continue;
				}
				if (events_isMC == -1){ // falsely reconstructed events (mainly combinatorial background)
					checkhist[4]->Fill(_M);
				}
				if (events_isMC == 0){ // data events are filled as reconstructed
					checkhist[0]->Fill(_M);
					// reference on pointer on pointer... yepp, it is still working ;)
					if ((binlow <= _M) && (_M < binhigh)){
						ofstream* _pstream = (*files.upper_bound(_M)).second[0];
						WriteEventToStream(particles, *_pstream);
					}
				}
				if (events_isMC == 1 || events_isMC == 2){ // mc exact events are filled as generated
					checkhist[1]->Fill(_M_mc);
					if ((binlow <= _M_mc) && (_M_mc < binhigh)){
						ofstream* _pstream = (*files.upper_bound(_M_mc)).second[1];
						WriteEventToStream(particles_mc, *_pstream);
					}
				}
				if (events_isMC == 2){ // mc accepted events are filled as reconstructed
					checkhist[2]->Fill(_M);
					checkhist[3]->Fill(_M_mc); // check what the true distribution without moved bins was
					if ((binlow <= _M) && (_M < binhigh)){
						ofstream* _pstream = (*files.upper_bound(_M)).second[2];
						WriteEventToStream(particles, *_pstream);
					}
				}
				/*
				for (int i = 0; i < 3; i++){
					if ((i == events_isMC) || // not only the value it self is counted
						(i == 1 && events_isMC == 2)){ // but also accepted events are MC exact events
						// reference on pointer on pointer... yepp, it is still working ;)
						ofstream* _pstream = (*files.upper_bound(_M)).second[i];
						WriteEventToStream(particles, *_pstream);
						checkhist[i]->Fill(_M);
					}
				}*/
			//}
		}
		rootfile.Close();
	}

	// close all files
	cout << " deleting objects " << endl;
	for (map<int, ofstream**>::iterator it = files.begin(); it != files.end(); it++) {
		ofstream** _files = (it->second);
		for (int j = 0; j < 3; j++){
			//cout << i << " " << j << endl;
			ofstream* _file = _files[j];
			if (_file){
				_file->close();
				delete _file;
			}
		}
	}
	cout << " done " << endl;
	// Draw a crosscheck histogram
	TCanvas canvas("canvas", "canvas", 1000, 600);
	canvas.cd();
	for (int i=0; i<5; i++){
		checkhist[i]->DrawClone("P Text");
		string addstring("");
		if (i==0){
			_hist_data = (TH1I*)checkhist[i]->Clone("saved_hist_data");
		}
		if (i==1){
			addstring = ".genbod";
			_hist_data_mc = (TH1I*)checkhist[i]->Clone("saved_hist_data_mc");
		}
		if (i==2){
			addstring = ".acc";
			_hist_data_mc_acc = (TH1I*)checkhist[i]->Clone("saved_hist_data_mc_acc");
		}
		if (i < 3){
			canvas.Print( (_dir+"/checkhist"+addstring+".pdf").c_str());
		}
		else {
			if (i==3) canvas.Print((_dir+"/checkhist_mc_truth.pdf").c_str());
			if (i==4) canvas.Print((_dir+"/checkhist_mc_false_reconstructed.pdf").c_str());
		}
		delete checkhist[i];
	}
	return result;
}

int TrpwaEventTreeHandler::GetDir (string path,
		vector<string> &files, string filterext, bool rmext){
    DIR *dp;
    struct dirent *dirp;
    files.clear();
    if((dp  = opendir(path.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << path << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
    	string _filename = dirp->d_name;
    	if ((filterext != "dir") && (
    			(_filename.size() <= filterext.size()) ||
    			(filterext != "" && _filename.compare(_filename.size()-filterext.size(), filterext.size(), filterext) != 0))) continue;
    	if (rmext) _filename.erase(_filename.size()-filterext.size(), filterext.size());
    	// filter only directories if requested
    	if (filterext == "dir") {
    			if (!DirExists(path+"/"+_filename)) continue;
    			// filter paths to mother daughter directories out
    			if (_filename == "." || _filename == "..") continue;
    	};
    	//cout << _filename << endl;
    	files.push_back(_filename);
    }
    closedir(dp);
    return (signed) files.size();
}

bool TrpwaEventTreeHandler::Set_bin_path(string dirname){
	bool result = true;
	_bin_paths.clear();
	_bins.clear();
	if (!DirExists(dirname)){
		cout << " Error in TrpwaEventTreeHandler::Set_bin_path(): directory " << dirname << " does not exist!" << endl;
		return false;
	}
	_dir = dirname;
	// retrieve all folders in the given path
	vector<string> dirnames;
	GetDir(dirname, dirnames, "dir");
	for (vector<string>::iterator it = dirnames.begin(); it != dirnames.end(); it++){
		if (!DirExists(_dir+"/"+*it)){
			result = false;
			_bin_paths.clear();
			_bins.clear();
			return result;
		}
		int binlow;
		int binhigh;
		if (Path_to_bin(_dir+"/"+*it, binlow, binhigh)){
			_bin_paths[binlow]=*it;
			_bins[binlow]=binhigh;
		}
	}
	cout << " found " << _bin_paths.size() << " bin directories " << endl;
	return result;
}

void TrpwaEventTreeHandler::Reset(){
	for (vector<TFile*>::iterator it = _eventtreefiles.begin(); it != _eventtreefiles.end(); it++){
		(*it)->Clear();
		delete (*it);
	}
	_eventtreefiles.clear();
	_eventtreenames.clear();
	_bins.clear();
	_bin_paths.clear();
}



TrpwaEventTreeHandler::~TrpwaEventTreeHandler(){
	Reset();
}



bool TrpwaEventTreeHandler::Path_to_bin(string path, int & binlow, int & binhigh){
	bool result = true;
	binhigh = -1;
	binlow  = -1;
	if (!DirExists(path)){
		cout << " Error in TrpwaEventTreeHandler::Path_to_bin(): " << path << " does not exist! " << endl;
		return false;
	}
	// search for <binlow>.<binhigh> structure
	stringstream _binlow;
	stringstream _binhigh;
	// get rid of an eventual trailing Slash
	if (path[path.size()-1] == '/')
		path.erase(path.size()-1,1);
	// get rid of the path
	int slashpos = path.rfind('/');
	if (slashpos != (int) string::npos) {
		path.erase(0, slashpos + 1);
	}
	// find the '.'
	int pointpos = path.rfind('.');
	if (pointpos == (int) string::npos) {
		cout << " Error in TrpwaEventTreeHandler::Path_to_bin(): " << path << " is not a valid folder name! " << endl;
		return false;
	}
	// get the bin values
	_binlow.str(path.substr(0, pointpos));
	_binhigh.str(path.substr(pointpos+1, path.size()-pointpos-1));
	//cout << _binlow.str() << " " << _binlow.str().size() << " . " << _binhigh.str() << " " << _binhigh.str().size() << endl;
	_binlow >> binlow;
	_binhigh >> binhigh;
	//cout << binlow << " " << binhigh << endl;
	// check the result
	if (binlow < 0 || binhigh < 0 || binhigh <= binlow){
		binlow = -1;
		binhigh= -1;
		cout << " Error in TrpwaEventTreeHandler::Path_to_bin(): " << path << " is not a valid folder name! " << endl;
		return false;
	}
	return result;
}

// check whether a file exists
bool TrpwaEventTreeHandler::FileExists(string filename){
	  ifstream ifile(filename.c_str());
	  return ifile;
}

// check whether a directory exists
bool TrpwaEventTreeHandler::DirExists(string dirname){
	struct stat st;
	if(stat(dirname.c_str(),&st) == 0 && S_ISDIR(st.st_mode))
		return true;
	else
		return false;
}

void TrpwaEventTreeHandler::WriteEventToStream(vector<TParticle>& particles, ofstream& stream){
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

void TrpwaEventTreeHandler::DrawProgressBar(int len, double percent) {
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



