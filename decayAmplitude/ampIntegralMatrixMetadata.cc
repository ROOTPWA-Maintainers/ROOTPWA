#include"ampIntegralMatrixMetadata.h"
#include"waveDescription.h"
#include"hashCalculator.h"
#include<TFile.h>

rpwa::ampIntegralMatrixMetadata::ampIntegralMatrixMetadata():
	_contentHash(""),
	_rootpwaGitHash(""),
	_objectBaseName(""),
	_allZeroHash(""),
	_ampIntegralMatrix(0),
	_amplitudeHashes(),
	_keyFileContents(),
	_evtMetas(){};

rpwa::ampIntegralMatrixMetadata::~ampIntegralMatrixMetadata(){};

std::string rpwa::ampIntegralMatrixMetadata::recalculateHash() const {
	if (!_ampIntegralMatrix) {
		printErr << "no ampIntegralMatrix given" << std::endl;
		return "";
	};
	hashCalculator heWhoHashes;
	size_t nWaves = _ampIntegralMatrix->nmbWaves();
	for (size_t i = 0; i < nWaves; ++i) {
		for (size_t j = 0; j < nWaves; ++j) {
			heWhoHashes.Update(_ampIntegralMatrix->element(i,j));
		};
	};	
	return heWhoHashes.hash();
};

std::ostream& rpwa::ampIntegralMatrixMetadata::print(std::ostream& out) const {
	out << "ampIntegralMatrixMetadata:" << std::endl
	    << "    contentHash ......... '" << _contentHash << "'"        << std::endl
	    << "    object base name .... '" << _objectBaseName << "'"     << std::endl
	    << "    rootpwa git hash .... '" << _rootpwaGitHash << "'"     << std::endl;
	if(_ampIntegralMatrix) {
		out << "    nmbwaves ... "  << _ampIntegralMatrix->nmbWaves() << std::endl;
	};
	out << std::endl;
	return out;
};


const rpwa::ampIntegralMatrixMetadata* rpwa::ampIntegralMatrixMetadata::readIntegralFile(TFile* inputFile, const std::string& objectBaseName, const bool& quiet) {
	const std::pair<std::string, std::string> objectNames = ampIntegralMatrixMetadata::getObjectNames(objectBaseName);
	rpwa::ampIntegralMatrixMetadata* integralMeta = (ampIntegralMatrixMetadata*)inputFile->Get(objectNames.second.c_str());
	if(not integralMeta) {
		if(not quiet) {
			printWarn << "could not find amplitude metadata." << std::endl;
		};
		return 0;
	};
	integralMeta->_ampIntegralMatrix = (rpwa::ampIntegralMatrix*)inputFile->Get(rpwa::ampIntegralMatrix::integralObjectName.c_str());
	if(not integralMeta->_ampIntegralMatrix) {
		if(not quiet) {
			printWarn << "could not find ampIntegralMatrix." << std::endl;
		};
		return 0;
	};
	return integralMeta;
};


std::pair<std::string, std::string> rpwa::ampIntegralMatrixMetadata::getObjectNames(const std::string& objectBaseName) {
	std::stringstream sstr;
	sstr << objectBaseName << ".amp";
	std::pair<std::string, std::string> retval = std::pair<std::string, std::string>();
	retval.first = sstr.str();
	sstr.str("");
	sstr << objectBaseName << ".meta";
	retval.second = sstr.str();
	return retval;
};


Int_t rpwa::ampIntegralMatrixMetadata::Write(const char* name, Int_t option, Int_t bufsize) const {
	Int_t retval = 0;
	if(_ampIntegralMatrix) {
		retval = _ampIntegralMatrix->Write(rpwa::ampIntegralMatrix::integralObjectName.c_str());
	};
	return retval + TObject::Write(name, option, bufsize);
};


bool rpwa::ampIntegralMatrixMetadata::mergeIntegralMatrix(const ampIntegralMatrixMetadata& second) {
	ampIntegralMatrix* secondMatrix = second.getAmpIntegralMatrix();
	if (secondMatrix->nmbWaves() != _ampIntegralMatrix->nmbWaves()) {
		printErr << "nmbWave() does not match" << std::endl;
		return false;
	};
	if (second.rootpwaGitHash() != rootpwaGitHash()) {
		printErr << "gitHashes differ" << std::endl;
		return false;

	};
	{
		std::vector<std::string> secondKeyFileContents = second.getKeyFileContents();
		for (size_t cont = 0; cont < secondKeyFileContents.size(); ++cont) {
			if (!hasKeyFileContent(secondKeyFileContents[cont])) {
				printErr << "unknown keyFileContent found" << std::endl;
				return false;
			};
		};
	};{
		std::vector<std::string> secondAmplitudeHashes = second.getAmplitudeHashes();
		for (size_t hash = 0; hash < secondAmplitudeHashes.size(); ++hash) {
			if (!addAmplitudeHash(secondAmplitudeHashes[hash])) {
				printErr << "hash already known. cant use same amplitudes twice" << std::endl;
				return false;
			};
		};
	};
	if (!mergeBinningMap(second.binningMap())) {
		printErr << "could not merge binningMaps" << std::endl;
		return false;
	};
	{
		std::vector<std::pair<rpwa::eventMetadata, std::vector<std::pair<size_t, size_t> > > > secondMetas = second.evtMetas();
		for (size_t meta = 0; meta < secondMetas.size(); ++meta) {
			for (size_t range = 0; range < secondMetas[meta].second.size(); ++range) {
				if (!addEventMetadata(secondMetas[meta].first, secondMetas[meta].second[range].first, secondMetas[meta].second[range].second)) {
					printErr << "could not add eventMetadata" << std::endl;
					return false;	
				};
			};
		};
	};
	*_ampIntegralMatrix+=*secondMatrix;
	if (!setHash()) {
		printErr << "could not set hash" << std::endl;
		return false;
	};
	return true;
};

bool rpwa::ampIntegralMatrixMetadata::setAmpIntegralMatrix(ampIntegralMatrix* matrix) {
	if (!matrix) {
		printErr << "no matrix given" << std::endl;
		return false;
	};
	_ampIntegralMatrix = matrix;	
	if (!setAllZeroHash()) {
		printErr << "could not set _allZeroHash" << std::endl;
		return false;
	};
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::setHash() {
	if (!_ampIntegralMatrix) {
		printErr << "no ampIntegralMatrix set" << std::endl;
		return false;
	};
	_contentHash = recalculateHash();
	if (_contentHash == ""){
		printErr << "_contentHash is empty" << std::endl;
		return false;
	};
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::setGitHash(const std::string &gitHash) {
	_rootpwaGitHash = gitHash;
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::setObjectBaseName(const std::string &baseName){
	_objectBaseName = baseName;
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::addKeyFileContent(const std::string &content) {
	for (size_t key = 0; key < _keyFileContents.size(); ++key) {
		if (_keyFileContents[key] == content) {
			printWarn << "keyFileContent already set" << std::endl;
			return false;
		};
	};
	_keyFileContents.push_back(content);
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::setBinningMap(const std::map<std::string, std::pair<double, double> > &binningMap) {
	_binningMap = binningMap;
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::mergeBinningMap(const std::map<std::string, std::pair<double, double> > &binningMapIn) {
	std::map<std::string, std::pair<double, double> > newMap;
	typedef std::map<std::string, std::pair<double, double> >::const_iterator it_type;
	for(it_type iterator = binningMapIn.begin(); iterator != binningMapIn.end(); ++iterator) {
		const std::string binningVariable = iterator->first;
		if (!_binningMap.count(binningVariable)) {
			printErr << binningVariable << " not in _binningMap" << std::endl;
			return false;
		};	
		newMap[binningVariable] = std::pair<double, double>(std::min(iterator->second.first, _binningMap[binningVariable].first), std::max(iterator->second.second, _binningMap[binningVariable].second));
	};
	if (newMap.size() != _binningMap.size()) {
		printErr << "sizes of _binningMap and newMap do not match. numbers of binning variables do not match" << std::endl;
		return false;
	};
	_binningMap = newMap;
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::addAmplitudeHash(const std::string &hash) {
	if (hash == _allZeroHash) {
		return true; // Do not add the allZerohash, but allow it
	}
	for (size_t i = 0; i < _amplitudeHashes.size(); ++i) {
		if (_amplitudeHashes[i] == hash) {
			printErr << "hash " << hash << " already in _amplitudeHashes" << std::endl;
			return false;
		};
	};
	_amplitudeHashes.push_back(hash);
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::setAllZeroHash() {
	if (!_ampIntegralMatrix) {
		_allZeroHash = "";
		printErr << "not integralMatrix set" << std::endl;
		return false;
	}
	size_t nEvts = _ampIntegralMatrix->nmbEvents();
	hashCalculator heWhoHashes;
	for (size_t i = 0; i < nEvts; ++i) {
		heWhoHashes.Update(std::complex<double>(0.,0.));
	}
	_allZeroHash = heWhoHashes.hash();
	return true;
}


bool rpwa::ampIntegralMatrixMetadata::addEventMetadata(const rpwa::eventMetadata& evtMeta, size_t eventMin, size_t eventMax) {
	if (eventMin >= eventMax) {
		printErr << "given event number borders are not sorted" << std::endl;
		return false;
	};
	std::string hash = evtMeta.contentHash();
	for (size_t meta = 0; meta < _evtMetas.size(); ++meta) {
		if (_evtMetas[meta].first.contentHash() == hash) {
			std::vector<std::pair<size_t, size_t> > newRanges;
			bool isIn = false;
			for (size_t set = 0; set < _evtMetas[meta].second.size(); ++set) {
				if (eventMin == _evtMetas[meta].second[set].first) {
					printErr << "two ranges start at the same event index" << std::endl;
					return false;
				};
				if (eventMin < _evtMetas[meta].second[set].first && !isIn) {
					isIn = true;
					newRanges.push_back(std::pair<size_t, size_t>(eventMin, eventMax));
				};
				newRanges.push_back(_evtMetas[meta].second[set]);
			};
			if (!isIn) {
				newRanges.push_back(std::pair<size_t, size_t>(eventMin, eventMax));
			};
			for (size_t  set = 1; set < newRanges.size(); ++set) {
				if (newRanges[set-1].second > newRanges[set].first) {
					printErr << "overlapping ranges found" << std::endl;
					return false;
				};
			};
			_evtMetas[meta].second = std::vector<std::pair<size_t, size_t> >(1, newRanges[0]);
			for (size_t set = 1; set < newRanges.size(); ++set) {
				if (newRanges[set].first == _evtMetas[meta].second[_evtMetas[meta].second.size()-1].second) {
					_evtMetas[meta].second[_evtMetas[meta].second.size()-1].second = newRanges[set].second;
				} else {
					_evtMetas[meta].second.push_back(newRanges[set]);
				};
			};
			return true;
		};
	};
	std::pair<eventMetadata, std::vector<std::pair<size_t, size_t> > > newEvtMeta(evtMeta, std::vector<std::pair<size_t, size_t> >(1, std::pair<size_t, size_t>(eventMin, eventMax)));
	return true;
};


bool rpwa::ampIntegralMatrixMetadata::check() const {
	if (!_ampIntegralMatrix) {
		printErr << "no _ampIntegralMatrix" << std::endl;
		return false;
	};
	std::vector<std::string> waveNamesFromKeyFiles;
	bool retVal = true;
	for (size_t key = 0; key < _keyFileContents.size(); ++key){ 
		std::vector<waveDescriptionPtr> waveDescriptions = rpwa::waveDescription::parseKeyFileContent(_keyFileContents[key]);
		for (size_t wave = 0; wave < waveDescriptions.size(); ++wave) {
			isobarDecayTopologyPtr topo;
			if (!waveDescriptions[wave]->constructDecayTopology(topo)) {
				printErr << "problems constructiong decayTopology" << std::endl;
				retVal = false;
				continue;
			};
			waveNamesFromKeyFiles.push_back(rpwa::waveDescription::waveNameFromTopology(*topo));
		};
	};
	for (size_t name = 0; name < waveNamesFromKeyFiles.size(); ++name) {
		if (!_ampIntegralMatrix->containsWave(waveNamesFromKeyFiles[name])) {
			printErr << "_ampIntegralMatrix does not have wave name: " << waveNamesFromKeyFiles[name] << std::endl;
			retVal = false;
		};
	};
	size_t nmbWaves = _ampIntegralMatrix->nmbWaves();
	if (nmbWaves != waveNamesFromKeyFiles.size()) {
		printErr << nmbWaves << " = nmbWaves != waveNamesFromKeyFiles.size() = " << waveNamesFromKeyFiles.size() << std::endl;
		retVal = false;
	};
	return retVal;
};



bool rpwa::ampIntegralMatrixMetadata::hasAmplitudeHash(const std::string& hash) const {
	for (size_t amp = 0; amp < _amplitudeHashes.size(); ++amp) { 
		if (_amplitudeHashes[amp] == hash) {
			return true;
		};
	};
	return false;
};


bool rpwa::ampIntegralMatrixMetadata::hasKeyFileContent(const std::string& content) const {
	for (size_t cont = 0; cont < _keyFileContents.size(); ++cont) {
		if (_keyFileContents[cont] == content) {
			return true;
		};
	};
	return false;
};


bool rpwa::ampIntegralMatrixMetadata::writeToFile(TFile* outputFile) {
	if (!setHash()) {
		printErr << "could not setHash. Abort writing..." << std::endl;
		return false;
	};
	if (!check()) {
		printErr << "check failed" << std::endl;
		return false;
	};	
	outputFile->cd();
	if (Write(getObjectNames(_objectBaseName).second.c_str()) == 0){
		printErr << "write failed" << std::endl;
		return false;
	};
	return true;
};
