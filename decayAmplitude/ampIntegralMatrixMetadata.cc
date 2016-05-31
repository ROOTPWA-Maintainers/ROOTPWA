#include"ampIntegralMatrixMetadata.h"

#include"waveDescription.h"
#include"hashCalculator.h"

#include<TFile.h>


using namespace std;

const string rpwa::ampIntegralMatrixMetadata::objectNameInFile = "integralMetadata";

rpwa::ampIntegralMatrixMetadata::ampIntegralMatrixMetadata()
	: _contentHash(""),
	  _rootpwaGitHash(""),
	  _allZeroHash(""),
	  _ampIntegralMatrix(0),
	  _amplitudeHashes(),
	  _keyFileContents(),
	  _evtMetas() { }

rpwa::ampIntegralMatrixMetadata::~ampIntegralMatrixMetadata() { }


string rpwa::ampIntegralMatrixMetadata::recalculateHash() const {
	if (not _ampIntegralMatrix) {
		printWarn << "trying to calculate hash without an integral matrix." << endl;
		return "";
	}
	hashCalculator hashor;
	size_t nWaves = _ampIntegralMatrix->nmbWaves();
	for (size_t i = 0; i < nWaves; ++i) {
		for (size_t j = 0; j < nWaves; ++j) {
			hashor.Update(_ampIntegralMatrix->element(i,j));
		}
	}
	return hashor.hash();
}


ostream& rpwa::ampIntegralMatrixMetadata::print(ostream& out) const {
	out << "ampIntegralMatrixMetadata:"                                << endl
	    << "    contentHash ......... '" << _contentHash << "'"        << endl
	    << "    rootpwa git hash .... '" << _rootpwaGitHash << "'"     << endl
	    << "    binning map";
	if(_binningMap.empty()) {
		out << " ......... " << "<empty>" << endl;
	} else {
		out << ": " << endl;
		for(binningMapType::const_iterator it = _binningMap.begin(); it != _binningMap.end(); ++it) {
			out << "        variable '" << it->first << "' range " << it->second << endl;
		}
	}
	if(_ampIntegralMatrix) {
		out << "    nmbwaves ............ " << _ampIntegralMatrix->nmbWaves() << endl;
	}
	out << endl;
	out << "connected event metadata information:" << endl;
	for(unsigned int i = 0; i < _evtMetas.size(); ++i) {
		out << _evtMetas[i];
		out << endl;
	}
	return out;
}


const rpwa::ampIntegralMatrixMetadata*
rpwa::ampIntegralMatrixMetadata::readIntegralFile(TFile* inputFile,
                                                  const bool& quiet)
{
	rpwa::ampIntegralMatrixMetadata* integralMeta = (ampIntegralMatrixMetadata*)inputFile->Get(objectNameInFile.c_str());
	if(not integralMeta) {
		if(not quiet) {
			printWarn << "could not find amplitude metadata object '" << objectNameInFile << "'." << endl;
		}
		return 0;
	}
	integralMeta->_ampIntegralMatrix = (rpwa::ampIntegralMatrix*)inputFile->Get(rpwa::ampIntegralMatrix::integralObjectName.c_str());
	if(not integralMeta->_ampIntegralMatrix) {
		if(not quiet) {
			printWarn << "could not find ampIntegralMatrix '" << rpwa::ampIntegralMatrix::integralObjectName << "'." << endl;
		}
		return 0;
	}
	return integralMeta;
}


Int_t rpwa::ampIntegralMatrixMetadata::Write(const char* name, Int_t option, Int_t bufsize) const {
	Int_t retval = 0;
	if(_ampIntegralMatrix) {
		retval = _ampIntegralMatrix->Write(rpwa::ampIntegralMatrix::integralObjectName.c_str());
	}
	return retval + TObject::Write(name, option, bufsize);
}


bool rpwa::ampIntegralMatrixMetadata::mergeIntegralMatrix(const ampIntegralMatrixMetadata& second) {
	const ampIntegralMatrix* secondMatrix = second.getAmpIntegralMatrix();
	if (secondMatrix->nmbWaves() != _ampIntegralMatrix->nmbWaves()) {
		printErr << "mismatch in number of waves (" << secondMatrix->nmbWaves()
		         << "!=" << _ampIntegralMatrix->nmbWaves() << ")." << endl;
		return false;
	}
	if (second.rootpwaGitHash() != rootpwaGitHash()) {
		printErr << "mismatch in git hashes." << endl;
		return false;

	}
	{
		const vector<string> secondKeyFileContents = second.getKeyFileContents();
		for (size_t keyfiles_i = 0; keyfiles_i < secondKeyFileContents.size(); ++keyfiles_i) {
			if (not hasKeyFileContent(secondKeyFileContents[keyfiles_i])) {
				printErr << "mismatch in keyfiles." << endl;
				return false;
			}
		}
	}
	{
		const vector<string> secondAmplitudeHashes = second.getAmplitudeHashes();
		for (size_t hash_i = 0; hash_i < secondAmplitudeHashes.size(); ++hash_i) {
			if (not addAmplitudeHash(secondAmplitudeHashes[hash_i])) {
				printErr << "could not add amplitude hash." << endl;
				return false;
			}
		}
	}
	if (not mergeBinningMap(second.binningMap())) {
		printErr << "could not merge binning maps." << endl;
		return false;
	}
	{
		vector<rpwa::eventMetadata> secondMetas = second.evtMetas();
		for (size_t meta_i = 0; meta_i < secondMetas.size(); ++meta_i) {
			if (not addEventMetadata(secondMetas[meta_i])) {
			        printErr << "could not add event metadata." << endl;
			        return false;
			}
		}
	}
	*_ampIntegralMatrix += *secondMatrix;
	if (not setHash()) {
		printErr << "could not set hash." << endl;
		return false;
	}
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::setAmpIntegralMatrix(ampIntegralMatrix* matrix) {
	if (not matrix) {
		printErr << "got null-pointer." << endl;
		return false;
	}
	_ampIntegralMatrix = matrix;
	if (not setAllZeroHash()) {
		printErr << "could not set the hash for zero-amplitudes." << endl;
		return false;
	}
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::setHash() {
	if (not _ampIntegralMatrix) {
		printWarn << "integral matrix has not been set." << endl;
		return false;
	}
	_contentHash = recalculateHash();
	if (_contentHash == ""){
		printWarn << "hash calculation failed." << endl;
		return false;
	}
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::addKeyFileContent(const string& content) {
	if (hasKeyFileContent(content)) {
		printWarn << "cannot add keyfile content: already present." << endl;
		return false;
	}
	_keyFileContents.push_back(content);
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::mergeBinningMap(const map<string, pair<double, double> >& binningMap) {
	if (binningMap.size() != _binningMap.size()) {
		printErr << "numbers of binning variables do not match." << endl;
		return false;
	}
	map<string, pair<double, double> > newMap;
	typedef map<string, pair<double, double> >::const_iterator it_type;
	for(it_type iterator = binningMap.begin(); iterator != binningMap.end(); ++iterator) {
		const string binningVariable = iterator->first;
		if (_binningMap.find(binningVariable) == _binningMap.end()) {
			printErr << "variable '" << binningVariable << "' not in binning map." << endl;
			return false;
		}
		newMap[binningVariable] = pair<double, double>(std::min(iterator->second.first, _binningMap[binningVariable].first),
		                                               std::max(iterator->second.second, _binningMap[binningVariable].second));
	}
	_binningMap = newMap;
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::addAmplitudeHash(const string& hash) {
	if (hash == _allZeroHash) {
		return true;  // Do not add the allZerohash, but allow it
	}
	if (std::find(_amplitudeHashes.begin(), _amplitudeHashes.end(), hash) != _amplitudeHashes.end()) {
		printWarn << "hash '" << hash << "' already in amplitude hashes." << endl;
		return false;
	}
	_amplitudeHashes.push_back(hash);
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::setAllZeroHash() {
	if (not _ampIntegralMatrix) {
		_allZeroHash = "";
		printErr << "integral matrix has not been set." << endl;
		return false;
	}
	hashCalculator hashor;
	for (size_t i = 0; i < _ampIntegralMatrix->nmbEvents(); ++i) {
		hashor.Update(complex<double>(0.,0.));
	}
	_allZeroHash = hashor.hash();
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::addEventMetadata(const rpwa::eventMetadata& evtMeta) {
	if (std::find(_evtMetas.begin(), _evtMetas.end(), evtMeta) != _evtMetas.end()) {
		printWarn << "event metadata object has already been added." << endl;
		return false;
	}
	_evtMetas.push_back(evtMeta);
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::check() const {
	if (not _ampIntegralMatrix) {
		printWarn << "integral matrix not set." << endl;
		return false;
	}
	vector<string> waveNamesFromKeyFiles;
	bool retVal = true;
	for (size_t keyfile_i = 0; keyfile_i < _keyFileContents.size(); ++keyfile_i){
		vector<waveDescriptionPtr> waveDescriptions = rpwa::waveDescription::parseKeyFileContent(_keyFileContents[keyfile_i]);
		for (size_t wave_i = 0; wave_i < waveDescriptions.size(); ++wave_i) {
			isobarDecayTopologyPtr topo;
			if (not waveDescriptions[wave_i]->constructDecayTopology(topo)) {
				printWarn << "problems constructiong decayTopology" << endl;
				retVal = false;
				continue;
			}
			waveNamesFromKeyFiles.push_back(rpwa::waveDescription::waveNameFromTopology(*topo));
		}
	}
	for (size_t name_i = 0; name_i < waveNamesFromKeyFiles.size(); ++name_i) {
		if (not _ampIntegralMatrix->containsWave(waveNamesFromKeyFiles[name_i])) {
			printWarn << "_ampIntegralMatrix does not have wave '" << waveNamesFromKeyFiles[name_i] << "'."<< endl;
			retVal = false;
		}
	}
	size_t nmbWaves = _ampIntegralMatrix->nmbWaves();
	if (nmbWaves != waveNamesFromKeyFiles.size()) {
		printWarn << "mismatch between number of keyfiles (" << waveNamesFromKeyFiles.size() << ")"
		          << " and number of waves in matrix (" << nmbWaves << ")." << endl;
		retVal = false;
	}
	return retVal;
}


bool rpwa::ampIntegralMatrixMetadata::hasAmplitudeHash(const string& hash) const {
	return (std::find(_amplitudeHashes.begin(), _amplitudeHashes.end(), hash) != _amplitudeHashes.end());
}


bool rpwa::ampIntegralMatrixMetadata::hasKeyFileContent(const string& content) const {
	return (std::find(_keyFileContents.begin(), _keyFileContents.end(), content) != _keyFileContents.end());
}


bool rpwa::ampIntegralMatrixMetadata::writeToFile(TFile* outputFile) {
	if (not setHash()) {
		printErr << "could not setHash, not writing to file." << endl;
		return false;
	}
	if (not check()) {
		printErr << "metadata invalid, not writing to file." << endl;
		return false;
	}
	outputFile->cd();
	if (Write(objectNameInFile.c_str()) == 0) {
		printErr << "write failed." << endl;
		return false;
	}
	return true;
}
