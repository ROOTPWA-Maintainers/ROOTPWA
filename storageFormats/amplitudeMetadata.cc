
#include "amplitudeMetadata.h"

#include <algorithm>

#include <boost/progress.hpp>

#include <TFile.h>
#include <TTree.h>

#include "amplitudeTreeLeaf.h"
#include "hashCalculator.h"
#include "reportingUtils.hpp"

using namespace rpwa;
using namespace std;

const std::string rpwa::amplitudeMetadata::amplitudeLeafName = "amplitude";

rpwa::amplitudeMetadata::amplitudeMetadata()
	: _contentHash(""),
	  _eventMetadata(),
	  _keyfileContent(""),
	  _rootpwaGitHash(""),
	  _objectBaseName(""),
	  _amplitudeTree(0) { }


rpwa::amplitudeMetadata::~amplitudeMetadata() { }


string rpwa::amplitudeMetadata::recalculateHash(const bool& printProgress) const
{
	amplitudeTreeLeaf* ampTreeLeaf = 0;
	hashCalculator hashor;
	if(not _amplitudeTree) {
		printWarn << "input tree not found in metadata." << endl;
		return "";
	}
	if(_amplitudeTree->SetBranchAddress(rpwa::amplitudeMetadata::amplitudeLeafName.c_str(), &ampTreeLeaf) < 0)
	{
		printWarn << "could not set address for branch '" << rpwa::amplitudeMetadata::amplitudeLeafName << "'." << endl;
		return "";
	}
	boost::progress_display* progressIndicator = printProgress ? new boost::progress_display(_amplitudeTree->GetEntries(), cout, "") : 0;
	for(long eventNumber = 0; eventNumber < _amplitudeTree->GetEntries(); ++eventNumber) {
		_amplitudeTree->GetEntry(eventNumber);
		if(progressIndicator) {
			++(*progressIndicator);
		}
		hashor.Update(ampTreeLeaf->amp());
	}
	return hashor.hash();
}


ostream& rpwa::amplitudeMetadata::print(ostream& out) const
{
	out << "amplitudeMetadata:" << endl
	    << "    contentHash ......... '" << _contentHash << "'"        << endl
	    << "    object base name .... '" << _objectBaseName << "'"     << endl
	    << "    rootpwa git hash .... '" << _rootpwaGitHash << "'"     << endl;
	if(_amplitudeTree) {
		out << "    amplitude entries ... "  << _amplitudeTree->GetEntries() << endl;
	}
	out << endl;
	out << "connected event metadata information:" << endl;
	for(unsigned int i = 0; i < _eventMetadata.size(); ++i) {
		out << _eventMetadata[i];
		out << endl;
	}
	return out;
}


const amplitudeMetadata* rpwa::amplitudeMetadata::readAmplitudeFile(TFile* inputFile, const string& objectBaseName, const bool& quiet)
{
	const pair<string, string> objectNames = amplitudeMetadata::getObjectNames(objectBaseName);
	amplitudeMetadata* amplitudeMeta = (amplitudeMetadata*)inputFile->Get(objectNames.second.c_str());
	if(not amplitudeMeta) {
		if(not quiet) {
			printWarn << "could not find amplitude metadata." << endl;
		}
		return 0;
	}
	amplitudeMeta->_amplitudeTree = (TTree*)inputFile->Get(objectNames.first.c_str());
	if(not amplitudeMeta->_amplitudeTree) {
		if(not quiet) {
			printWarn << "could not find amplitude tree." << endl;
		}
		return 0;
	}
	return amplitudeMeta;
}


pair<string, string> rpwa::amplitudeMetadata::getObjectNames(const string& objectBaseName)
{
	std::stringstream sstr;
	sstr << objectBaseName << ".amp";
	pair<string, string> retval = pair<string, string>();
	retval.first = sstr.str();
	sstr.str("");
	sstr << objectBaseName << ".meta";
	retval.second = sstr.str();
	return retval;
}


Int_t rpwa::amplitudeMetadata::Write(const char* name, Int_t option, Int_t bufsize) const
{
	Int_t retval = 0;
	if(_amplitudeTree) {
		retval = _amplitudeTree->Write();
	}
	return retval + TObject::Write(name, option, bufsize);
}

vector<vector<complex<double>>>
rpwa::loadAmplitudes(const vector<string>& ampFilenames,
                     const vector<string>& waveNames,
                     const string& eventFilename,
                     const multibinBoundariesType& otfBin,
                     long unsigned maxNmbEvents) {
	TFile* eventFile = TFile::Open(eventFilename.c_str());
	if (eventFile == nullptr or not eventFile->IsOpen()){
		printErr << "Cannot open event file '" << eventFilename << "'." << endl;
		throw;
	}

	const eventMetadata* eventMeta = eventMetadata::readEventFile(eventFile);
	if (eventMeta == nullptr) {
		printErr << "Cannot read event medatata from file '" << eventFilename << "'." << endl;
		throw;
	}

	additionalTreeVariables variables;
	if(otfBin.empty()) {
		printErr << "got event metadata but the binning map is emtpy." << endl;
		throw;
	}
	TTree* eventTree = eventMeta->eventTree();
	const unsigned long totNmbEvents = eventTree->GetEntries();
	if (totNmbEvents == 0) {
		printWarn << "event trees contain no amplitudes values." << endl;
		throw;
	}
	maxNmbEvents = (maxNmbEvents==0)? totNmbEvents: min(maxNmbEvents, totNmbEvents);

	if (not variables.setBranchAddresses(*eventMeta)) {
		printErr << "cannot set branch address to additional variables." << endl;
		throw;
	}

	// build list of event indices within the current multibin
	vector<long> eventIndicesInMultibin;
	for(unsigned long iEvent = 0; iEvent < totNmbEvents; ++iEvent){
		eventTree->GetEntry(iEvent);
		if (variables.inBoundaries(otfBin)) eventIndicesInMultibin.push_back(iEvent);
		if (eventIndicesInMultibin.size() >= maxNmbEvents) break;
	}
	printInfo << "load " << eventIndicesInMultibin.size() << " events of " << totNmbEvents << " events!" << endl;
	eventFile->Close();


	if (ampFilenames.empty()) {
		printErr << "did not receive any amplitude filenames." << endl;
		throw;
	}

	if (ampFilenames.size() != waveNames.size()){
		printErr << "Length of amplitude files (" << ampFilenames.size() << ") is different form number of wave names (" << waveNames.size() << ")." << endl;
		throw;
	}
	vector<vector<complex<double>>> amps;
	amps.reserve(ampFilenames.size());
	for(size_t waveIndex = 0; waveIndex < ampFilenames.size(); waveIndex++) {
		amplitudeTreeLeaf* ampTreeLeaf = nullptr;
		TFile* ampFile = TFile::Open(ampFilenames[waveIndex].c_str());
		if (ampFile == nullptr or not ampFile->IsOpen()){
			printErr << "Cannot open amplitude file '" << ampFilenames[waveIndex] << "'." << endl;
			throw;
		}
		const amplitudeMetadata* ampMeta = amplitudeMetadata::readAmplitudeFile(ampFile, waveNames[waveIndex]);

		if (totNmbEvents != static_cast<unsigned long>(ampMeta->amplitudeTree()->GetEntries())) {
			printErr << "amplitude trees do not all have the same entry count as the event tree." << endl;
			throw;
		}

		ampMeta->amplitudeTree()->SetBranchAddress(amplitudeMetadata::amplitudeLeafName.c_str(), &ampTreeLeaf);
		vector<complex<double>> ampsOfWave;
		ampsOfWave.reserve(eventIndicesInMultibin.size());
		for(const auto iEvent: eventIndicesInMultibin){
			ampMeta->amplitudeTree()->GetEntry(iEvent);

				if (ampTreeLeaf->nmbIncohSubAmps() != 1){
					printErr << "Cannot handle more than one subamplitude at the moment!" << endl;
					throw;
				}
				ampsOfWave.push_back(ampTreeLeaf->amp());
		}
		amps.push_back(ampsOfWave);
		ampFile->Close();
	}
	return amps;
}
