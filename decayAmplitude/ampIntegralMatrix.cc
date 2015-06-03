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
//
// Description:
//      container class for complex amplitude integral matrices
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <boost/progress.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "TClass.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TROOT.h"
#include "TBuffer.h"

#include "reportingUtils.hpp"
#include "fileUtils.hpp"
#include "sumAccumulators.hpp"
#include "amplitudeTreeLeaf.h"
#include "ampIntegralMatrix.h"
#include "amplitudeMetadata.h"


using namespace std;
using namespace boost;
using namespace boost::accumulators;
using namespace rpwa;


ClassImp(ampIntegralMatrix);


const std::string rpwa::ampIntegralMatrix::integralObjectName = "integral";
bool ampIntegralMatrix::_debug = false;


ampIntegralMatrix::ampIntegralMatrix()
	: TObject               (),
	  _nmbWaves             (0),
	  _waveNameToIndexMap   (),
	  _waveNames            (),
	  _nmbEvents            (0),
	  _integrals            (),
	  _intStorageShape      (),
	  _intStorageNmbElements(0),
	  _intStorageData       (0),
	  _waveDescriptions     ()
{
	//ampIntegralMatrix::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


ampIntegralMatrix::ampIntegralMatrix(const ampIntegralMatrix& integral)
	: TObject()
{
	*this = integral;
}


ampIntegralMatrix::~ampIntegralMatrix()
{ }


void
ampIntegralMatrix::clear()
{
	_nmbWaves  = 0;
	_nmbEvents = 0;
	_waveNameToIndexMap.clear();
	_waveNames.clear();
	_integrals.resize(extents[0][0]);
	_intStorageShape.clear();
	_intStorageNmbElements = 0;
	_intStorageData        = 0;
	_waveDescriptions.clear();
}


ampIntegralMatrix&
ampIntegralMatrix::operator =(const ampIntegralMatrix& integral)
{
	if (this != &integral) {
		TObject::operator   =(integral);
		_nmbWaves           = integral._nmbWaves;
		_waveNameToIndexMap = integral._waveNameToIndexMap;
		_waveNames          = integral._waveNames;
		_nmbEvents          = integral._nmbEvents;
		// multiarray's = operator does not seem to set shape correctly
		// setting manually prevents crash in glibc
		const sizeType* shape = integral._integrals.shape();
		_integrals.resize(extents[shape[0]][shape[1]]);
		_integrals        = integral._integrals;
		_waveDescriptions = integral._waveDescriptions;
	}
	return *this;
}


ampIntegralMatrix&
ampIntegralMatrix::operator +=(const ampIntegralMatrix& integral)
{
	if (not hasIdenticalWaveSet(integral)) {
		printErr << "cannot add " << *this << endl
		         << "and " << integral << endl
		         << "because the two integral matrices have different wave sets. Aborting..." << endl;
		throw;
	}
	for (unsigned int i = 0; i < _nmbWaves; ++i) {
		const string waveNameI = waveName(i);
		for (unsigned int j = 0; j < _nmbWaves; ++j)
			_integrals[i][j] += integral.matrix()[integral.waveIndex(waveNameI)][integral.waveIndex(waveName(j))];
	}
	_nmbEvents       += integral.nmbEvents();
	_waveDescriptions = integral._waveDescriptions;
	return *this;
}


ampIntegralMatrix&
ampIntegralMatrix::operator -=(const ampIntegralMatrix& integral)
{
	if (not hasIdenticalWaveSet(integral)) {
		printErr << "cannot subtract " << integral << endl
		         << "from " << *this << endl
		         << "because the two integral matrices have different wave sets. Aborting..." << endl;
		throw;
	}
	for (unsigned int i = 0; i < _nmbWaves; ++i) {
		const string waveNameI = waveName(i);
		for (unsigned int j = 0; j < _nmbWaves; ++j)
			_integrals[i][j] -= integral.matrix()[integral.waveIndex(waveNameI)][integral.waveIndex(waveName(j))];
	}
	_nmbEvents       -= integral.nmbEvents();
	_waveDescriptions = integral._waveDescriptions;
	return *this;
}


ampIntegralMatrix&
ampIntegralMatrix::operator *=(const double factor)
{
	for (unsigned int i = 0; i < _nmbWaves; ++i)
		for (unsigned int j = 0; j < _nmbWaves; ++j)
			_integrals[i][j] *= factor;
	_nmbEvents *= factor;
	return *this;
}


ampIntegralMatrix&
ampIntegralMatrix::operator /=(const double factor)
{
	for (unsigned int i = 0; i < _nmbWaves; ++i)
		for (unsigned int j = 0; j < _nmbWaves; ++j)
			_integrals[i][j] /= factor;
	_nmbEvents /= factor;
	return *this;
}


bool
ampIntegralMatrix::containsWave(const string& waveName) const
{
	waveNameToIndexMapIterator entry = _waveNameToIndexMap.find(waveName);
	if (entry == _waveNameToIndexMap.end())
		return false;
	return true;
}


unsigned int
ampIntegralMatrix::waveIndex(const string& waveName) const
{
	waveNameToIndexMapIterator entry = _waveNameToIndexMap.find(waveName);
	if (entry == _waveNameToIndexMap.end()) {
		printErr << "cannot find wave '" << waveName << "' in integral matrix. Aborting..." << endl;
		throw;
	}
	return entry->second;
}


const string&
ampIntegralMatrix::waveName(const unsigned int waveIndex) const
{
	if (waveIndex < _waveNames.size())
		return _waveNames[waveIndex];
	else {
		printErr << "wave index " << waveIndex << " is out of range. Aborting..." << endl;
		throw;
	}
}


const waveDescription*
ampIntegralMatrix::waveDesc(const unsigned int waveIndex) const
{
	if (waveIndex < _waveDescriptions.size())
		return &(_waveDescriptions[waveIndex]);
	return 0;
}


bool
ampIntegralMatrix::allWavesHaveDesc() const
{
	for (unsigned i = 0; i < _waveDescriptions.size(); ++i)
		if (not _waveDescriptions[i].keyFileParsed())
			return false;
	return true;
}


complex<double>
ampIntegralMatrix::element(const unsigned int waveIndexI,
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
	return _integrals[waveIndexI][waveIndexJ] / ((double)_nmbEvents);
}


bool
ampIntegralMatrix::integrate(const vector<TTree*>&           ampTrees,
                             const std::vector<std::string>& waveNames,
                             const unsigned long             maxNmbEvents,
                             const string&                   weightFileName)
{
	if (ampTrees.size() > 0)
		printInfo << "calculating integral for " << ampTrees.size() << " amplitude(s)" << endl;
	else {
		printWarn << "did not receive any amplitude trees. cannot calculate integral." << endl;
		return false;
	}
	if (waveNames.size() != ampTrees.size()) {
		printErr << "number of amplitude trees does not match the number of waves in integral matrix." << endl;
		return false;
	}
	for(unsigned int i=0; i < waveNames.size(); i++) {
		if (waveNames[i] == "") {
			printErr << "can not add empty wave name.";
			return false;
		}
		_nmbWaves += 1;
		_waveNames.push_back(waveNames[i]);
		_waveNameToIndexMap.insert( std::pair<std::string, unsigned int> (waveNames[i], _nmbWaves-1) );
	}
	const unsigned long nmbEvents = (unsigned long) ampTrees[0]->GetEntries();
	if (nmbEvents == 0) {
		printWarn << "amplitude trees contain no amplitudes values. cannot calculate integral." << endl;
		return false;
	}
	for (unsigned int i=1; i < ampTrees.size(); i++) {
		if (nmbEvents != (unsigned int) ampTrees[i]->GetEntries()) {
			printErr << "amplitude trees do not all have the same entry count." << endl;
			return false;
		}
	}
	if (maxNmbEvents == 0)
		_nmbEvents = nmbEvents;
	else
		_nmbEvents = min(nmbEvents, maxNmbEvents);
	// set amplitude tree leafs
	vector<amplitudeTreeLeaf*> ampTreeLeafs(_nmbWaves);
	for(unsigned int waveIndex = 0; waveIndex < _nmbWaves; waveIndex++) {
		ampTreeLeafs[waveIndex] = NULL;
		ampTrees[waveIndex]->SetBranchAddress(rpwa::amplitudeMetadata::amplitudeLeafName.c_str(), &ampTreeLeafs[waveIndex]);
	}

	// make sure that either all or none of the waves have description (needed?)
	if (not allWavesHaveDesc())
		_waveDescriptions.clear();

	// resize integral matrix
	//_integrals.clear();
	_integrals.resize(extents[_nmbWaves][_nmbWaves]);

	// open importance sampling weight file
	//!!! this should be provided as a friend tree for the amlitudes
	ifstream weightFile;
	bool     useWeight = false;
	if (weightFileName != "") {
		if (_debug)
			printDebug << "opening importance sampling weight file '" << weightFileName << "'" << endl;
		weightFile.open(weightFileName.c_str());
		if (not weightFile) {
			printWarn << "cannot open importance sampling weight file '" << weightFileName << "' "
			          << "cannot calculate integral." << endl;
			return false;
		}
		useWeight = true;
	}

	// loop over events and calculate integral matrix
	accumulator_set<double, stats<tag::sum(compensated)> > weightAcc;
	typedef accumulator_set<complex<double>, stats<tag::sum(compensated)> > complexAcc;
	vector<vector<complexAcc> > ampProdAcc(_nmbWaves, vector<complexAcc>(_nmbWaves));
	// process weight file and binary amplitude files
	vector<vector<complex<double> > > amps(_nmbWaves);
	progress_display progressIndicator(_nmbEvents, cout, "");
	bool             success = true;
	for (unsigned long iEvent = 0; iEvent < _nmbEvents; ++iEvent) {
		++progressIndicator;

		// sum up importance sampling weight
		double w = 1;
		if (useWeight)
			if (not(weightFile >> w)) {
				success = false;
				printWarn << "error reading weight file. stopping integration "
				          << "at event " << iEvent << " of total " << _nmbEvents << "." << endl;
				break;
			}
		const double weight = 1 / w; // we have to undo the weighting of the events!
		weightAcc(weight);

		// read amplitude values for this event from root trees
		for (unsigned int waveIndex = 0; waveIndex < _nmbWaves; ++waveIndex) {
			ampTrees[waveIndex]->GetEntry(iEvent);
			const unsigned int nmbSubAmps = ampTreeLeafs[waveIndex]->nmbIncohSubAmps();
			if (nmbSubAmps < 1) {
				printErr << "amplitude object for wave '" << _waveNames[waveIndex] << "' "
				         << "does not contain any amplitude values "
				         << "at event " << iEvent << " of total " << _nmbEvents << ". Aborting..." << endl;
				throw;
			}
			// get all incoherent subamps
			amps[waveIndex].resize(nmbSubAmps);
			for (unsigned int subAmpIndex = 0; subAmpIndex < nmbSubAmps; ++subAmpIndex)
				amps[waveIndex][subAmpIndex] = ampTreeLeafs[waveIndex]->incohSubAmp(subAmpIndex);
		}

		// sum up integral matrix elements
		for (unsigned int waveIndexI = 0; waveIndexI < _nmbWaves; ++waveIndexI)
			for (unsigned int waveIndexJ = 0; waveIndexJ < _nmbWaves; ++waveIndexJ) {
				// sum over incoherent subamps
				const unsigned int nmbSubAmps = amps[waveIndexI].size();
				if (nmbSubAmps != amps[waveIndexJ].size()) {
					printErr << "number of incoherent sub-amplitudes for wave '"
					         << _waveNames[waveIndexI] << "' = " << nmbSubAmps
					         << " differs from that of wave '" << _waveNames[waveIndexJ] << "' = "
					         << amps[waveIndexJ].size()
					         << " at event " << iEvent << " of total " << _nmbEvents << ". Aborting... "
					         << "be sure to use only .root amplitude files, "
					         << "if your channel has sub-amplitudes." << endl;
					throw;
				}
				complex<double> val = 0;
				for (unsigned int subAmpIndex = 0; subAmpIndex < nmbSubAmps; ++subAmpIndex)
					val += amps[waveIndexI][subAmpIndex] * conj(amps[waveIndexJ][subAmpIndex]);
				if (useWeight)
					val *= weight;
				ampProdAcc[waveIndexI][waveIndexJ](val);
			}
	}  // event loop

	// copy values from accumulators and (if necessary) renormalize to
	// integral of importance sampling weights
	const double weightNorm = sum(weightAcc) / (double)_nmbEvents;
	for (unsigned int waveIndexI = 0; waveIndexI < _nmbWaves; ++waveIndexI)
		for (unsigned int waveIndexJ = 0; waveIndexJ < _nmbWaves; ++waveIndexJ) {
			_integrals[waveIndexI][waveIndexJ] = sum(ampProdAcc[waveIndexI][waveIndexJ]);
			if (useWeight)
				_integrals[waveIndexI][waveIndexJ] *= 1 / weightNorm;
		}

	printSucc << "calculated integrals of " << _nmbWaves << " amplitude(s) "
	          << "for " << _nmbEvents << " events" << endl;
	return success;
}


void
ampIntegralMatrix::renormalize(const unsigned long nmbEventsRenorm)
{
	if (_debug)
		printDebug << "renormalizing integrals from " << _nmbEvents << " "
		           << "to " << nmbEventsRenorm << " events." << endl;
	*this *= (double)nmbEventsRenorm / _nmbEvents;
	_nmbEvents = nmbEventsRenorm;
}


bool
ampIntegralMatrix::writeAscii(ostream& out) const
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
	out << _waveNameToIndexMap.size() << endl;
	for (waveNameToIndexMapIterator i = _waveNameToIndexMap.begin();
	     i != _waveNameToIndexMap.end(); ++i)
		out << i->first << " " << i->second << endl;
	return true;
}


bool
ampIntegralMatrix::readAscii(istream& in)
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
	//_integrals.clear();
	_integrals.resize(extents[nmbRows][nmbCols]);
	for (unsigned int i = 0; i < nmbRows; ++i)
		for (unsigned int j = 0; j < nmbCols; ++j) {
			if (not (in >> _integrals[i][j]) or in.eof()) {
				printErr << "could not read integral values. stream seems truncated." << endl;
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
			printErr << "could not read wave name -> index map. stream seems truncated." << endl;
			return false;
		}
		_waveNameToIndexMap[waveName] = waveIndex;
		--mapSize;
	}
	_waveNames.resize(_nmbWaves);
	for (waveNameToIndexMapIterator i = _waveNameToIndexMap.begin();
	     i != _waveNameToIndexMap.end(); ++i)
		_waveNames[i->second] = i->first;
	return true;
}


bool
ampIntegralMatrix::writeAscii(const string& outFileName) const
{
	if (_debug)
		printDebug << "opening ASCII file '" << outFileName << "' for writing of integral matrix" << endl;
	ofstream outFile(outFileName.c_str());
	if (not outFile or not outFile.good()) {
		printWarn << "cannot open file '" << outFileName << "' for writing of integral matrix" << endl;
		return false;
	}
	const bool success = writeAscii(outFile);
	if (success)
		printSucc << "wrote integral to ASCII file '" << outFileName << "'" << endl;
	else
		printWarn << "problems writing integral to ASCII file '" << outFileName << "'" << endl;
	return success;
}


bool
ampIntegralMatrix::readAscii(const string& inFileName)
{
	if (_debug)
		printDebug << "opening ASCII file '" << inFileName << "' for reading of integral matrix" << endl;
	fstream inFile(inFileName.c_str());
	if (not inFile or not inFile.good()) {
		printWarn << "cannot open file '" << inFileName << "' for reading of integral matrix" << endl;
		return false;
	}
	const bool success = readAscii(inFile);
	if (success)
		printSucc << "read integral from ASCII file '" << inFileName << "'" << endl;
	else
		printWarn << "problems reading integral from ASCII file '" << inFileName << "'" << endl;
	return success;
}


ostream&
ampIntegralMatrix::print(ostream&   out,
                         const bool printIntegralValues) const
{
	out << "amplitude integral matrix:"  << endl
	    << "    number of waves .... " << _nmbWaves  << endl
	    << "    number of events ... " << _nmbEvents << endl
	    << "    wave name -> index map:" << endl;
	for (waveNameToIndexMapIterator i = _waveNameToIndexMap.begin();
	     i != _waveNameToIndexMap.end(); ++i)
		out << "        " << i->first << " -> " << i->second << " -> " << _waveNames[i->second] << endl;
	if (printIntegralValues) {
		out << "    integral matrix values:" << endl;
		for (unsigned int waveIndexI = 0; waveIndexI < _nmbWaves; ++waveIndexI)
			for (unsigned int waveIndexJ = 0; waveIndexJ < _nmbWaves; ++waveIndexJ)
				out << "        [" << setw(4) << waveIndexI << ", " << setw(4) << waveIndexJ << "] = "
				    << maxPrecisionDouble(_integrals[waveIndexI][waveIndexJ]) << endl;
	}
	return out;
}


void
ampIntegralMatrix::rebuildWaveNameToIndexMap()
{
	_waveNameToIndexMap.clear();
	for (unsigned int i = 0; i < _waveNames.size(); ++i)
		if (not containsWave(_waveNames[i]))
			_waveNameToIndexMap[_waveNames[i]] = i;
		else
			printWarn << "wave '" << _waveNames[i] << "' already exists in integral matrix. "
			          << "ignoring wave." << endl;
}


bool
ampIntegralMatrix::hasIdenticalWaveSet(const ampIntegralMatrix& integral) const
{
	if (nmbWaves() != integral.nmbWaves())
		return false;
	// check that all waves in this integral are also in other integral
	for (unsigned int i = 0; i < _nmbWaves; ++i)
		if (not integral.containsWave(waveName(i)))
			return false;
	// check that all waves in other integral are also in this integral
	for (unsigned int i = 0; i < integral.nmbWaves(); ++i)
		if (not containsWave(integral.waveName(i)))
			return false;
	return true;
}


void
ampIntegralMatrix::storeMultiArray()
{
	// set storage variables
	_intStorageShape.resize(_integrals.num_dimensions());
	for (sizeType i = 0; i < _integrals.num_dimensions(); ++i)
		_intStorageShape[i] = _integrals.shape()[i];
	_intStorageNmbElements = _integrals.num_elements();
	_intStorageData        = _integrals.data();
}


void
ampIntegralMatrix::readMultiArray()
{
	// rebuild multiarray from storage variables
	assert(_intStorageShape.size() == _integrals.num_dimensions());
	_integrals.resize(extents[_intStorageShape[0]][_intStorageShape[1]]);
	complex<double>* integralData = _integrals.data();
	for (unsigned int i = 0; i < _intStorageNmbElements; ++i)
		integralData[i] = _intStorageData[i];
}


// custom streamer that triggers copying of multiarray into C struct
// before writing and rebuilding of multiarrays after reading
//
// this might not work when ampIntegralMatrix is used as a TTree leaf
//
// a more elegant way would be to use read- and write-rule pragmas, but
// write rules are not yet supported
//
// see
//     http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=12213&p=54249&hilit=TStreamerInfo
//     http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=7334&p=30634&hilit=streamer+custom
void
ampIntegralMatrix::Streamer(TBuffer& R__b)
{
	if (R__b.IsReading()) {
		R__b.ReadClassBuffer(rpwa::ampIntegralMatrix::Class(), this);
		readMultiArray();
		rebuildWaveNameToIndexMap();
	} else {
		storeMultiArray();
		R__b.WriteClassBuffer(rpwa::ampIntegralMatrix::Class(), this);
	}
}
