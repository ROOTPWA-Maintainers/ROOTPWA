///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      container class for complex normalization integral matrices
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <fstream>

#include <boost/progress.hpp>

#include "TClass.h"

#include "reportingUtils.hpp"
#include "fileUtils.hpp"
#include "sumAccumulators.hpp"
#include "normalizationIntegral.h"


using namespace std;
using namespace boost;
using namespace boost::accumulators;
using namespace rpwa;


#if NORMALIZATIONINTEGRAL_ENABLED


ClassImp(normalizationIntegral);


bool normalizationIntegral::_debug = false;

    
normalizationIntegral::normalizationIntegral()
	: TObject   (),
	  _nmbWaves (0),
	  _nmbEvents(0)
{
	normalizationIntegral::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


normalizationIntegral::normalizationIntegral(const normalizationIntegral& integral)
{
	*this = integral;
}


normalizationIntegral::~normalizationIntegral()
{ }


normalizationIntegral&
normalizationIntegral::operator =(const normalizationIntegral& integral)
{
	if (this != &integral) {
		TObject::operator =(integral);
		_nmbWaves         = integral._nmbWaves;
		_waveNameIndexMap = integral._waveNameIndexMap;
		_nmbEvents        = integral._nmbEvents;
		_integrals        = integral._integrals;
	}
	return *this;
}


unsigned int
normalizationIntegral::waveIndex(const std::string& waveName) const
{
	waveNameIndexMapIterator entry = _waveNameIndexMap.find(waveName);
	if (entry == _waveNameIndexMap.end()) {
		printErr << "cannot find wave '" << waveName << "' in normalization integral" << endl;
		throw;
	}
	return entry->second;
}


const complex<double>&
normalizationIntegral::element(const unsigned int waveIndexI,
                               const unsigned int waveIndexJ) const
{
	if (waveIndexI >= _nmbWaves) {
		printErr << "wave index for i = " << waveIndexI << " out of range. "
		         << "number of waves = " << _nmbWaves << endl;
		throw;
	}
	if (waveIndexJ >= _nmbWaves) {
		printErr << "wave index for j = " << waveIndexJ << " out of range. "
		         << "number of waves = " << _nmbWaves << endl;
		throw;
	}
	return _integrals[waveIndexI][waveIndexJ];
	// el(iName, jName) / ((double)_nmbEvents);
}


const complex<double>&
normalizationIntegral::element(const std::string& waveNameI,
                               const std::string& waveNameJ) const
{
	return _integrals[waveIndex(waveNameI)][waveIndex(waveNameJ)];
	// el(iName, jName) / ((double)_nmbEvents);
}


bool
normalizationIntegral::integrate(const vector<string>& ampFileNames,
                                 const unsigned int    maxNmbEvents,
                                 const string&         weightFileName)
{
	// open amplitude files
	vector<ifstream*> ampFiles;
	streampos         ampFileLength = 0; 
	_nmbWaves = 0;  // counts waves
	for (vector<string>::const_iterator i = ampFileNames.begin(); i != ampFileNames.end(); ++i) {
		// get file path
		const string& ampFilePath = *i;
		// open amplitude file
		if (_debug)
			printInfo << "opening amplitude file '" << ampFilePath << "'" << endl;
		ifstream* ampFile = new ifstream();
		ampFile->open(ampFilePath.c_str());
		if(not *ampFile) {
			printWarn << "cannot open amplitude file '" << ampFilePath << "'. skipping." << endl;
			continue;
		}
		// check that all files have the same size
		const streampos fileSize = rpwa::fileSize(*ampFile);
		if (fileSize == 0) {
			printWarn << "amplitude file '" << ampFilePath << "' has zero size. skipping." << endl;
			continue;
    }
    if (ampFileLength == 0)
	    ampFileLength = fileSize;
    else if (fileSize != ampFileLength) {
			printWarn << "amplitude file '" << ampFilePath << "' has different size "
			          << "than previous file. skipping." << endl;
			continue;
    }
		ampFiles.push_back(ampFile);
		// get wave name and fill name -> index map
		const string waveName = fileNameFromPath(ampFilePath);
		_waveNameIndexMap[waveName] = _nmbWaves;
		++_nmbWaves;
  }
	if (ampFiles.size() > 0)
		printInfo << "calculating integral from " << ampFiles.size() << " amplitude file(s)" << endl;
	else {
		printWarn << "could not open any amplitude files. cannot calculate integral." << endl;
		return false;
	}

	// resize integral matrix
	_integrals.clear();
	_integrals.resize(_nmbWaves, vector<complex<double> >(_nmbWaves, 0));

	// open importance sampling weight file
	ifstream weightFile;
	bool     useWeight = false;
	if (weightFileName != "") {
		if (_debug)
			printInfo << "opening importance sampling weight file '" << weightFileName << "'" << endl;
		weightFile.open(weightFileName.c_str());
		if (not weightFile) { 
			printWarn << "cannot open importance sampling weight file '" << weightFileName << "' "
			          << "cannot calculate integral." << endl;
			return false;
		}
		useWeight = true;
	}

	// loop over events and calculate integral matrix
	_nmbEvents = 0;
	complex<double>* amps    = new complex<double>[_nmbWaves];
	bool             ampEof  = false;
	streampos        lastPos = ampFiles[0]->tellg();
	progress_display progressIndicator(ampFileLength, cout, "");
	accumulator_set<double,          stats<tag::sum(compensated)> > weightIntegral;
	accumulator_set<complex<double>, stats<tag::sum(compensated)> > integrals[_nmbWaves][_nmbWaves];
	while (not ampEof and ((maxNmbEvents) ? _nmbEvents < maxNmbEvents : true)) {
		// read importance sampling weight (no error handling yet!)
		double w = 1;
		if (useWeight)
			weightFile >> w;
		const double weight = 1 / w; // we have to de-weight the events!
		weightIntegral(weight);
		// read amplitude values for this event
		for (unsigned int waveIndex = 0; waveIndex < _nmbWaves; ++waveIndex) {
			ampFiles[waveIndex]->read((char*)&amps[waveIndex], sizeof(complex<double>));
			if ((ampEof = ampFiles[waveIndex]->eof()))
				break;
		}
		if (ampEof)
			break;
		++_nmbEvents;
		progressIndicator += ampFiles[0]->tellg() - lastPos;
		lastPos            = ampFiles[0]->tellg();
		// sum up integral matrix elements
		for (unsigned int waveIndexI = 0; waveIndexI < _nmbWaves; ++waveIndexI)
			for (unsigned int waveIndexJ = 0; waveIndexJ < _nmbWaves; ++waveIndexJ) {
				complex<double> val = amps[waveIndexI] * conj(amps[waveIndexJ]);
				if (useWeight)
					val *= weight;
				integrals[waveIndexI][waveIndexJ](val);
			}
	}
	// copy values from accumulators and (if necessary) renormalize to
	// integral of importance sampling weights
	const double weightNorm = sum(weightIntegral) / (double)_nmbEvents;
	for (unsigned int waveIndexI = 0; waveIndexI < _nmbWaves; ++waveIndexI)
		for (unsigned int waveIndexJ = 0; waveIndexJ < _nmbWaves; ++waveIndexJ) {
			_integrals[waveIndexI][waveIndexJ] = sum(integrals[waveIndexI][waveIndexJ]);
			if (useWeight)
				_integrals[waveIndexI][waveIndexJ] *= 1 / weightNorm;
		}


	for (unsigned int waveIndex = 0; waveIndex < _nmbWaves; ++waveIndex)
		delete ampFiles[waveIndex];
	ampFiles.clear();
	delete [] amps;
	return true;
}


bool
normalizationIntegral::writeAscii(ostream& out) const
{
	if (not out) {
		printWarn << "cannot write to output stream. cannot write integral." << endl;
		return false;
	}
	out << _nmbWaves << endl
	    << _nmbEvents << endl;
	// write integral matrix
  out << _nmbWaves << " " << _nmbWaves << endl;
  for (unsigned int waveIndexI = 0; waveIndexI < _nmbWaves; ++waveIndexI) {
	  for (unsigned int waveIndexJ = 0; waveIndexJ < _nmbWaves; ++waveIndexJ)
		  out << maxPrecisionDouble(_integrals[waveIndexI][waveIndexJ]) << "\t";
	  out << endl;
  }
  // write wave name -> index map
	out << _waveNameIndexMap.size() << endl;
	for (waveNameIndexMapIterator i = _waveNameIndexMap.begin(); i != _waveNameIndexMap.end(); ++i)
		out << i->first << " " << i->second << endl;
	return true;
}


bool
normalizationIntegral::readAscii(istream& in)
{
	if (not (in >> _nmbWaves >> _nmbEvents)) {
		printWarn << "could not read number of waves and events" << endl;
		return false;
	}
  // read matrix elements
	unsigned int nmbRows, nmbCols;
	if (not (in >> nmbRows >> nmbCols)) {
		printWarn << "could not read number of rows and columns" << endl;
		return false;
	}
	// resize integral matrix
	_integrals.clear();
	_integrals.resize(nmbRows, vector<complex<double> >(nmbCols, 0));
	for (unsigned int i = 0; i < nmbRows; ++i)
    for (unsigned int j = 0; j < nmbCols; ++j) {
	    if (not (in >> _integrals[i][j]) or in.eof()) {
		    printErr << "could not read integral values. stream seems trunkated." << endl;
        throw;
	    }
    }
  // read wave name -> index map
  unsigned int mapSize;
  in >> mapSize;
  while ((mapSize > 0) and in) {
	  string       waveName;
	  unsigned int waveIndex;
	  if (not (in >> waveName >> waveIndex)) {
		  printErr << "could not read wave name -> index map. stream seems trunkated." << endl;
		  return false;
	  }
    _waveNameIndexMap[waveName] = waveIndex;
    --mapSize;
  }
  return true;
}


bool
normalizationIntegral::writeAscii(const string& outFileName) const
{
	if (_debug)
		printInfo << "opening ASCII file '" << outFileName << "' for writing of integral matrix" << endl;
	std::ofstream outFile(outFileName.c_str());
	if (not outFile or not outFile.good()) {
		printWarn << "cannot open file '" << outFileName << "' for writing of integral matrix" << endl;
		return false;
	}
	const bool success = writeAscii(outFile);
	if (success)
		printInfo << "successfully wrote integral to ASCII file '" << outFileName << "'" << endl;
	else
		printWarn << "problems writing integral to ASCII file '" << outFileName << "'" << endl;
	return success;
}


bool
normalizationIntegral::readAscii(const string& inFileName)
{
	if (_debug)
		printInfo << "opening ASCII file '" << inFileName << "' for reading of integral matrix" << endl;
	std::fstream inFile(inFileName.c_str());
	if (not inFile or not inFile.good()) {
		printWarn << "cannot open file '" << inFileName << "' for reading of integral matrix" << endl;
		return false;
	}
	const bool success = readAscii(inFile);
	if (success)
		printInfo << "successfully read integral from ASCII file '" << inFileName << "'" << endl;
	else
		printWarn << "problems reading integral from ASCII file '" << inFileName << "'" << endl;
	return success;
}


#endif  // NORMALIZATIONINTEGRAL_ENABLED
