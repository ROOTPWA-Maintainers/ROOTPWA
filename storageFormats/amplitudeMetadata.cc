
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
	std::string name = objectBaseName;
	std::stringstream sstr;
	sstr << name << ".amp";
	pair<string, string> retval = pair<string, string>();
	retval.first = sstr.str();
	sstr.str("");
	sstr << name << ".meta";
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
