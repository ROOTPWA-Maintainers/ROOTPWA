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
	: TObject          (),
	  _weightFileName  (""),
	  _maxNmbEvents    (0),
	  _nmbWaves        (0),
	  _nmbEvents       (0)
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
		_weightFileName   = integral._weightFileName;
		_maxNmbEvents     = integral._maxNmbEvents;
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
normalizationIntegral::integrate()
{
	if (_nmbWaves < 1) {
		printWarn << "no waves to integrate over" << endl;
		return false;
	}

	// open importance sampling weight file
	ifstream weightFile;
	bool     useWeight = false;
	if (_weightFileName != "") {
		if (_debug)
			printInfo << "opening importance sampling weight file '" << _weightFileName << "'" << endl;
		weightFile.open(_weightFileName.c_str());
		if (not weightFile) { 
			printWarn << "cannot open importance sampling weight file '" << _weightFileName << "' "
			          << "cannot calculate integral." << endl;
			return false;
		}
		useWeight = true;
	}

	// open amplitude files
	bool      success       = true;
	ifstream* ampFiles      = new ifstream[_nmbWaves];
	streampos ampFileLength = 0;
	for (waveNameIndexMapIterator i = _waveNameIndexMap.begin(); i != _waveNameIndexMap.end(); ++i) {
		const string       ampFileName = i->first;
		const unsigned int waveIndex   = i->second;
		if (_debug)
			printInfo << "opening amplitude file '" << ampFileName << "'" << endl;
		ampFiles[waveIndex].open((ampFileName).c_str());
		if(not ampFiles[waveIndex]) {
			printWarn << "cannot open amplitude file '" << ampFileName << "'" << endl;
			success = false;
		}
		// check that all files have the same size
    ampFiles[waveIndex].seekg(0, ios::end);
    streampos fileSize = ampFiles[waveIndex].tellg();
    ampFiles[waveIndex].seekg(0, ios::beg);
    if (fileSize == 0) {
			printWarn << "amplitude file '" << ampFileName << "' has zero size" << endl;
			success = false;
    }
    if (ampFileLength == 0)
	    ampFileLength = fileSize;
    else if (fileSize != ampFileLength) {
			printWarn << "amplitude file '" << ampFileName << "' has different than previous file" << endl;
			success = false;
    }
	}
	if (not success) {
		printWarn << "problems opening amplitude file(s). cannot calculate integral." << endl;
		return false;
	}

	// loop over events and calculate integral matrix
	_nmbEvents = 0;
	complex<double>* amps    = new complex<double>[_nmbWaves];
	bool             ampEof  = false;
	streampos        lastPos = ampFiles[0].tellg();
	progress_display progressIndicator(ampFileLength, cout, "");
	accumulator_set<double,          stats<tag::sum(compensated)> > weightIntegral;
	accumulator_set<complex<double>, stats<tag::sum(compensated)> > integrals[_nmbWaves][_nmbWaves];
	while (not ampEof and ((_maxNmbEvents) ? _nmbEvents < _maxNmbEvents : true)) {
		// read importance sampling weight (no error handling yet!)
		double w = 1;
		if (useWeight)
			weightFile >> w;
		const double weight = 1 / w; // we have to de-weight the events!
		weightIntegral(weight);
		// read amplitude values for this event
		for (unsigned int waveIndex = 0; waveIndex < _nmbWaves; ++waveIndex) {
			ampFiles[waveIndex].read((char*)&amps[waveIndex], sizeof(complex<double>));
			if ((ampEof = ampFiles[waveIndex].eof()))
				break;
		}
		if (ampEof)
			break;
		++_nmbEvents;
		progressIndicator += ampFiles[0].tellg() - lastPos;
		lastPos            = ampFiles[0].tellg();
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

	delete [] ampFiles;
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
	if (not (in >> _nmbWaves >> _nmbEvents))
		return false;
  // read matrix elements
	unsigned int nmbRows, nmbCols;
	if (not (in >> nmbRows >> nmbCols))
		return false;
	// resize integral matrix
	_integrals.clear();
	_integrals.resize(nmbRows, vector<complex<double> >(nmbCols, 0));
	for (unsigned int i = 0; i < nmbRows; ++i)
    for (unsigned int j = 0; j < nmbCols; ++j) {
	    if (not (in >> _integrals[i][j]) or in.eof()) {
		    printErr << "could not read integral. stream seems trunkated." << endl;
        throw;
	    }
    }
  // read wave name -> index map
  unsigned int mapSize;
  in >> mapSize;
  while ((mapSize > 0) and in) {
	  string       waveName;
	  unsigned int waveIndex;
	  if (not (in >> waveName >> waveIndex))
		  return false;
    _waveNameIndexMap[waveName] = waveIndex;
    --mapSize;
  }
  return true;
}


#endif  // NORMALIZATIONINTEGRAL_ENABLED
