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

#include "reportingUtils.hpp"
#include "fileUtils.hpp"
#include "sumAccumulators.hpp"
#include "amplitudeTreeLeaf.h"
#include "ampIntegralMatrix.h"


using namespace std;
using namespace boost;
using namespace boost::accumulators;
using namespace rpwa;


#ifdef USE_STD_COMPLEX_TREE_LEAFS
ClassImp(ampIntegralMatrix);
#endif


bool ampIntegralMatrix::_debug = false;

    
ampIntegralMatrix::ampIntegralMatrix()
	: TObject   (),
	  _nmbWaves (0),
	  _nmbEvents(0)
{
	ampIntegralMatrix::Class()->IgnoreTObjectStreamer();  // don't store TObject's fBits and fUniqueID
}


ampIntegralMatrix::ampIntegralMatrix(const ampIntegralMatrix& integral)
{
	*this = integral;
}


ampIntegralMatrix::~ampIntegralMatrix()
{ }


ampIntegralMatrix&
ampIntegralMatrix::operator =(const ampIntegralMatrix& integral)
{
	if (this != &integral) {
		TObject::operator     =(integral);
		_nmbWaves             = integral._nmbWaves;
		_waveNameWaveIndexMap = integral._waveNameWaveIndexMap;
		_waveIndexWaveNameMap = integral._waveIndexWaveNameMap;
		_nmbEvents            = integral._nmbEvents;
		_integrals            = integral._integrals;
	}
	return *this;
}


ampIntegralMatrix&
ampIntegralMatrix::operator +=(const ampIntegralMatrix& integral)
{
	if (not hasIdenticalWaveSet(integral)) {
		printErr << "cannot add " << *this << endl
		         << "and " << integral << endl
		         << "because the two integral matrices have different wave sets. aborting." << endl;
		throw;
	}
	for (unsigned int i = 0; i < nmbWaves(); ++i) {
		const string waveNameI = waveName(i);
		for (unsigned int j = 0; j < nmbWaves(); ++j)
			_integrals[i][j] += integral(waveIndex(waveNameI), waveIndex(waveName(j)));
	}
	_nmbEvents += integral.nmbEvents();
	return *this;
}


ampIntegralMatrix&
ampIntegralMatrix::operator -=(const ampIntegralMatrix& integral)
{
	if (not hasIdenticalWaveSet(integral)) {
		printErr << "cannot subtract " << integral << endl
		         << "from " << *this << endl
		         << "because the two integral matrices have different wave sets. aborting." << endl;
		throw;
	}
	for (unsigned int i = 0; i < nmbWaves(); ++i) {
		const string waveNameI = waveName(i);
		for (unsigned int j = 0; j < nmbWaves(); ++j)
			_integrals[i][j] -= integral(waveIndex(waveNameI), waveIndex(waveName(j)));
	}
	_nmbEvents -= integral.nmbEvents();
	return *this;
}


ampIntegralMatrix&
ampIntegralMatrix::operator *=(const double factor)
{
	for (unsigned int i = 0; i < nmbWaves(); ++i)
		for (unsigned int j = 0; j < nmbWaves(); ++j)
			_integrals[i][j] *= factor;
	_nmbEvents *= factor;
	return *this;
}


ampIntegralMatrix&
ampIntegralMatrix::operator /=(const double factor)
{
	for (unsigned int i = 0; i < nmbWaves(); ++i)
		for (unsigned int j = 0; j < nmbWaves(); ++j)
			_integrals[i][j] /= factor;
	_nmbEvents /= factor;
	return *this;
}


bool
ampIntegralMatrix::containsWave(const string& waveName) const
{
	waveNameWaveIndexMapIterator entry = _waveNameWaveIndexMap.find(waveName);
	if (entry == _waveNameWaveIndexMap.end())
		return false;
	return true;
}


unsigned int
ampIntegralMatrix::waveIndex(const string& waveName) const
{
	waveNameWaveIndexMapIterator entry = _waveNameWaveIndexMap.find(waveName);
	if (entry == _waveNameWaveIndexMap.end()) {
		printErr << "cannot find wave '" << waveName << "' in integral matrix. aborting." << endl;
		throw;
	}
	return entry->second;
}


const string&
ampIntegralMatrix::waveName(const unsigned int waveIndex) const
{
	if (waveIndex < _waveIndexWaveNameMap.size())
		return _waveIndexWaveNameMap[waveIndex];
	else {
		printErr << "wave index " << waveIndex << " is out of range. aborting." << endl;
		throw;
	}
}


complex<double>&
ampIntegralMatrix::operator ()(const unsigned int waveIndexI,
                               const unsigned int waveIndexJ)
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
}


const complex<double>&
ampIntegralMatrix::operator ()(const unsigned int waveIndexI,
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
}


bool
ampIntegralMatrix::integrate(const vector<string>& binAmpFileNames,
                             const vector<string>& rootAmpFileNames,
                             const unsigned long   maxNmbEvents,
                             const string&         weightFileName)
{
	// open amplitude files
	vector<ifstream*>          binAmpFiles;
	const unsigned long        nmbBinEvents  = openBinAmpFiles(binAmpFiles, binAmpFileNames, 0);
	const unsigned int         nmbBinWaves   = binAmpFiles.size();
	vector<TTree*>             ampTrees;
	vector<amplitudeTreeLeaf*> ampTreeLeafs;
	const unsigned long        nmbRootEvents = openRootAmpFiles(ampTrees, ampTreeLeafs,
	                                                            rootAmpFileNames, binAmpFiles.size());
	const unsigned int         nmbRootWaves  = ampTrees.size();
	_nmbWaves = nmbBinWaves + nmbRootWaves;
	if (_nmbWaves > 0)
		printInfo << "calculating integral for " << _nmbWaves << " amplitude(s)" << endl;
	else {
		printWarn << "could not open any amplitude files. cannot calculate integral." << endl;
		return false;
	}
	if ((nmbBinEvents > 0) and (nmbRootEvents > 0) and (nmbBinEvents != nmbRootEvents)) {
		printWarn << ".bin and .root amplitude files contain different number of amplitudes "
		          << "(" << nmbBinEvents << " vs. " << nmbRootEvents << "). "
		          << "cannot calculate integral." << endl;
		return false;
	}
	if ((nmbBinEvents == 0) and (nmbRootEvents == 0)) {
		printWarn << "amplitude files contain no amplitudes. cannot calculate integral." << endl;
		return false;
	}
	if (nmbBinEvents > 0)
		_nmbEvents = nmbBinEvents;
	else
		_nmbEvents = nmbRootEvents;
	_nmbEvents = (maxNmbEvents) ? min(_nmbEvents, maxNmbEvents) : _nmbEvents;

	// update wave index -> wave name map
  _waveIndexWaveNameMap.resize(_nmbWaves);
  for (waveNameWaveIndexMapIterator i = _waveNameWaveIndexMap.begin();
	     i != _waveNameWaveIndexMap.end(); ++i)
	  _waveIndexWaveNameMap[i->second] = i->first;

	// resize integral matrix
	_integrals.clear();
	_integrals.resize(_nmbWaves, vector<complex<double> >(_nmbWaves, 0));

	// open importance sampling weight file
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
	accumulator_set<double,          stats<tag::sum(compensated)> > weightIntegral;
	accumulator_set<complex<double>, stats<tag::sum(compensated)> > integrals[_nmbWaves][_nmbWaves];
	// process weight file and binary amplitude files
	complex<double>  amps[_nmbWaves];
  progress_display progressIndicator(_nmbEvents, cout, "");
  for (unsigned long iEvent = 0; iEvent < _nmbEvents; ++iEvent) {
	  ++progressIndicator;
	  
	  // read importance sampling weight
	  double w = 1;
	  if (useWeight)
		  if (not(weightFile >> w)) {
			  printWarn << "error reading weight file. stopping integration "
			            << "at event " << iEvent << " of total " << _nmbEvents << "." << endl;
			  break;
		  }
	  const double weight = 1 / w; // we have to de-weight the events!

	  // read amplitude values for this event from binary files
	  bool ampBinEof = false;
	  for (unsigned int waveIndex = 0; waveIndex < nmbBinWaves; ++waveIndex) {
		  binAmpFiles[waveIndex]->read((char*)&amps[waveIndex], sizeof(complex<double>));
		  if ((ampBinEof = binAmpFiles[waveIndex]->eof())) {
			  printWarn << "unexpected EOF while reading binary amplitude '"
			            << _waveIndexWaveNameMap[waveIndex] << "'. stopping integration "
			            << "at event " << iEvent << " of total " << _nmbEvents << "." << endl;
			  break;
		  }
	  }
	  if (ampBinEof)
		  break;

	  // read amplitude values for this event from root trees
	  for (unsigned int waveIndex = nmbBinWaves; waveIndex < _nmbWaves; ++waveIndex) {
		  const unsigned int index = waveIndex - nmbBinWaves;
		  ampTrees[index]->GetEntry(iEvent);
		  if (ampTreeLeafs[index]->nmbIncohSubAmps() < 1) {
			  printErr << "empty amplitude object for wave '" << _waveIndexWaveNameMap[waveIndex] << "' "
			           << "at event " << iEvent << " of total " << _nmbEvents << ". aborting." << endl;
			  throw;
		  }
		  //!!! at the moment only one amplitude per wave can be handled
		  amps[waveIndex] = ampTreeLeafs[index]->incohSubAmp(0);
	  }

	  // sum up integral matrix elements
	  weightIntegral(weight);
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

	printSucc << "calculated integrals of " << _nmbWaves << " amplitude(s) "
	          << "for " << _nmbEvents << " events" << endl;
  // cleanup
	for (unsigned int waveIndex = 0; waveIndex < _nmbWaves; ++waveIndex)
		delete binAmpFiles[waveIndex];
	binAmpFiles.clear();
	return true;
}


void
ampIntegralMatrix::renormalize(const unsigned long nmbEventsRenorm)
{
	if (_debug)
		printDebug << "renormalizing integrals from " << _nmbEvents << " "
		           << "to " << nmbEventsRenorm << " events." << endl;
  for (unsigned int waveIndexI = 0; waveIndexI < _nmbWaves; ++waveIndexI)
	  for (unsigned int waveIndexJ = 0; waveIndexJ < _nmbWaves; ++waveIndexJ)
		  _integrals[waveIndexI][waveIndexJ] *= (double)nmbEventsRenorm / _nmbEvents;
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
	out << _waveNameWaveIndexMap.size() << endl;
	for (waveNameWaveIndexMapIterator i = _waveNameWaveIndexMap.begin();
	     i != _waveNameWaveIndexMap.end(); ++i)
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
    _waveNameWaveIndexMap[waveName] = waveIndex;
    --mapSize;
  }
  _waveIndexWaveNameMap.resize(_nmbWaves);
  for (waveNameWaveIndexMapIterator i = _waveNameWaveIndexMap.begin();
	     i != _waveNameWaveIndexMap.end(); ++i)
	  _waveIndexWaveNameMap[i->second] = i->first;
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
	for (waveNameWaveIndexMapIterator i = _waveNameWaveIndexMap.begin();
	     i != _waveNameWaveIndexMap.end(); ++i)
		out << "        " << i->first << " -> " << i->second << " -> "
		    << _waveIndexWaveNameMap[i->second] << endl;
	if (printIntegralValues) {
		out << "    integral matrix values:" << endl;
		for (unsigned int waveIndexI = 0; waveIndexI < _nmbWaves; ++waveIndexI)
			for (unsigned int waveIndexJ = 0; waveIndexJ < _nmbWaves; ++waveIndexJ)
				out << "        [" << setw(4) << waveIndexI << ", " << setw(4) << waveIndexJ << "] = "
				    << maxPrecisionDouble(_integrals[waveIndexI][waveIndexJ]) << endl;
	}
	return out;
}


unsigned long
ampIntegralMatrix::openBinAmpFiles(vector<ifstream*>&    ampFiles,
                                   const vector<string>& ampFileNames,
                                   const unsigned int    waveIndexOffset)
{
	ampFiles.clear();
	streampos    ampFileSize = 0;
	unsigned int waveIndex   = 0;
	for (unsigned int i = 0; i < ampFileNames.size(); ++i) {

		// get file path
		const string& ampFilePath = ampFileNames[i];

		// open amplitude file
		if (_debug)
			printDebug << "opening binary amplitude file '" << ampFilePath << "'" << endl;
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
    if (ampFileSize == 0)
	    ampFileSize = fileSize;
    else if (fileSize != ampFileSize) {
			printWarn << "amplitude file '" << ampFilePath << "' has different size "
			          << "than previous file. skipping." << endl;
			continue;
    }
		ampFiles.push_back(ampFile);

		// get wave name and fill name -> index map
		const string waveName = fileNameFromPath(ampFilePath);
		_waveNameWaveIndexMap[waveName] = waveIndex + waveIndexOffset;
		++waveIndex;
  }
	return ampFileSize / sizeof(complex<double>);  // number of amplitudes
}


unsigned long
ampIntegralMatrix::openRootAmpFiles(vector<TTree*>&             ampTrees,
                                    vector<amplitudeTreeLeaf*>& ampTreeLeafs,
                                    const vector<string>&       ampFileNames,
                                    const unsigned int          waveIndexOffset,
                                    const string&               ampLeafName)
{
	ampTrees.clear    ();
	ampTreeLeafs.clear();
	unsigned long nmbAmps = 0;
#ifdef USE_STD_COMPLEX_TREE_LEAFS
	// force loading predefined std::complex dictionary
	// see http://root.cern.ch/phpBB3/viewtopic.php?f=5&t=9618&p=50164
	gROOT->ProcessLine("#include <complex>");
	unsigned int waveIndex = 0;
	for (unsigned int i = 0; i < ampFileNames.size(); ++i) {

		// get file path
		const string& ampFilePath = ampFileNames[i];

		// open amplitude file
		if (_debug)
			printDebug << "opening .root amplitude file '" << ampFilePath << "'" << endl;
		TFile* ampFile = TFile::Open(ampFilePath.c_str(), "READ");
		if (not ampFile or ampFile->IsZombie()) {
			printErr << "could open amplitude file '" << ampFilePath << "'. skipping." << endl;
			continue;
		}

		// find amplitude tree with name matching *.amp
		TTree*     ampTree = 0;
    TIterator* keys    = ampFile->GetListOfKeys()->MakeIterator();
    while (TKey* k = static_cast<TKey*>(keys->Next())) {
      if (!k) {
        printWarn << "NULL pointer to TKey in file '" << ampFilePath << "'. skipping." << endl;
        continue;
      }
      const string keyName      = k->GetName();
      unsigned int countMatches = 0;
      if (extensionFromPath(keyName) == "amp") {
	      ampFile->GetObject(keyName.c_str(), ampTree);
	      ++countMatches;
      }
      if (not ampTree) {
	      printWarn << "no key in file '" << ampFilePath << "' matches '*.amp'. skipping." << endl;
	      continue;
      }
      if (countMatches > 1) {
	      printWarn << "more than one key in file '" << ampFilePath << "' matches '*.amp'. "
	                << "skipping." << endl;
	      continue;
      }
    }

    // connect tree leaf
		amplitudeTreeLeaf* ampTreeLeaf = 0;
		ampTree->SetBranchAddress(ampLeafName.c_str(), &ampTreeLeaf);

		// check that all trees have the same number of entries
		const unsigned long nmbEntries = numeric_cast<unsigned long>(ampTree->GetEntriesFast());
		if (nmbEntries == 0) {
			printWarn << "amplitude tree '" << ampTree->GetName() << "' has zero entries. "
			          << "skipping." << endl;
			continue;
    }
    if (nmbAmps == 0)
	    nmbAmps = nmbEntries;
    else if (nmbEntries != nmbAmps) {
	    printWarn << "amplitude tree '" << ampTree->GetName() << "' has different number of entries "
			          << "than previous tree. skipping." << endl;
			continue;
    }
		ampTrees.push_back    (ampTree    );
		ampTreeLeafs.push_back(ampTreeLeaf);

		// get wave name and fill name -> index map
		string waveName = ampTree->GetName();
		waveName = waveName.substr(0, waveName.length() - 5);  // chop off ".amp"
		_waveNameWaveIndexMap[waveName] = waveIndex + waveIndexOffset;
		++waveIndex;
  }
#endif  // USE_STD_COMPLEX_TREE_LEAFS
	return nmbAmps;
}


bool
ampIntegralMatrix::hasIdenticalWaveSet(const ampIntegralMatrix& integral) const
{
	if (nmbWaves() != integral.nmbWaves())
		return false;
	// check that all waves in this integral are also in other integral
	for (unsigned int i = 0; i < nmbWaves(); ++i)
		if (not integral.containsWave(waveName(i)))
			return false;
	// check that all waves in other integral are also in this integral
	for (unsigned int i = 0; i < integral.nmbWaves(); ++i)
		if (not containsWave(integral.waveName(i)))
			return false;
	return true;
}
