#include "primeNumbers.h"

#include <algorithm>
#include <cmath>

#include <TFile.h>
#include <TTree.h>

#include <reportingUtils.hpp>

using namespace std;


rpwa::primeNumbers* rpwa::primeNumbers::_instance = 0;
const string rpwa::primeNumbers::TREE_NAME = "primeNumbers";
const string rpwa::primeNumbers::BRANCH_NAME = "primes";
const size_t rpwa::primeNumbers::_blockSize = 1e5;
const size_t rpwa::primeNumbers::_emergencyCacheSize = 10;


rpwa::primeNumbers::primeNumbers()
	: _cacheFile(0),
	  _cacheTree(0),
	  _treeEntry(0),
	  _primeNumberCache(_emergencyCacheSize, 0)
{
	// put some emergency primes here so that the static members in TFracNum.cc
	// can be initialized
	_primeNumberCache[0] = 2;
	size_t cacheSize = 1;
	entryType entry = 3;
	while(cacheSize < _emergencyCacheSize) {
		if(isPrimeNumber(entry)) {
			_primeNumberCache[cacheSize++] = entry;
		}
		entry += 2;
	}
}


rpwa::primeNumbers& rpwa::primeNumbers::instance()
{
	if(not _instance) {
		_instance = new rpwa::primeNumbers();
	}
	return *_instance;
}


const rpwa::primeNumbers::entryType& rpwa::primeNumbers::primeNumber(const size_t& index)
{
	const size_t initialSize = _primeNumberCache.size();
	if(index >= initialSize) {
		if(not _cacheTree) {
			printErr << "trying to get prime number from uninitialized object. Aborting..." << endl;
			throw;
		}
		size_t entriesToRead = ((size_t)(((double)index)/((double)_blockSize)) + 1) * _blockSize;
		if(entriesToRead >= (size_t)_cacheTree->GetEntries()) {
			if(index >= (size_t)_cacheTree->GetEntries()) {
				printErr << "prime number with index " << index << " is not in cache with size "
						 << _cacheTree->GetEntries() << ". Aborting..." << endl;
				throw;
			}
			entriesToRead = _cacheTree->GetEntries();
		}
		_primeNumberCache.resize(entriesToRead, 0);
		for(size_t i = initialSize; i < _primeNumberCache.size(); ++i) {
			if(_cacheTree->GetEntry(i) <= 0) {
				printErr << "could not get entry " << i << " from cache tree. Aborting..." << endl;
				throw;
			}
			_primeNumberCache[i] = _treeEntry;
		}
	}
	return _primeNumberCache[index];

}


bool rpwa::primeNumbers::readCacheFile(const std::string& fileName)
{
	if(_cacheFile or _cacheTree or _primeNumberCache.size() != _emergencyCacheSize) {
		printWarn << "we are already initialized." << endl;
		return false;
	}
	_cacheFile = TFile::Open(fileName.c_str(), "READ");
	if(not _cacheFile) {
		printWarn << "could not open prime number cache file '" << fileName << "'." << endl;
		return false;
	}
	_cacheFile->GetObject(TREE_NAME.c_str(), _cacheTree);
	if(not _cacheTree) {
		printWarn << "could not find tree '" << TREE_NAME << "' in file '" << fileName << "'." << endl;
		return false;
	}
	if(_cacheTree->GetEntries() <= 0) {
		printWarn << "the tree in file '" << fileName << "' is empty." << endl;
		return false;
	}
	int result = _cacheTree->SetBranchAddress(BRANCH_NAME.c_str(), const_cast<entryType*>(&_treeEntry));
	if(result) {
		printWarn << "could not set branch address for branch '" << BRANCH_NAME
		          << "' (" << result << ")"<< " in file '" << fileName << "'." << endl;
		return false;
	}
	size_t entriesToRead = min(_blockSize, (size_t)_cacheTree->GetEntries());
	_primeNumberCache.resize(entriesToRead, 0);
	for(size_t i = 0; i < _primeNumberCache.size(); ++i) {
		if(_cacheTree->GetEntry(i) <= 0) {
			printWarn << "could not get entry " << i << " from tree in file '" << fileName << "'." << endl;
			return false;
		}
		_primeNumberCache[i] = _treeEntry;
	}
	return true;
}


rpwa::primeNumbers::~primeNumbers()
{
	_cacheFile->Close();
	delete _cacheFile;
}


bool rpwa::primeNumbers::isPrimeNumber(const entryType& number)
{
	if(number == 2) {
		return true;
	}
	const entryType maxCandidate = (entryType)(sqrt((double)number));
	for(entryType i = 2; i <= max(maxCandidate, (entryType)2); ++i) {
		if(not (number % i)) {
			return false;
		}
	}
	return true;
}
