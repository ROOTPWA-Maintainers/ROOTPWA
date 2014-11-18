
#include "amplitudeFileWriter.h"

#include <TFile.h>

#include "amplitudeTreeLeaf.h"
#include "eventMetadata.h"
#include "reportingUtils.hpp"
#include "reportingUtilsEnvironment.h"

using namespace rpwa;
using namespace std;


rpwa::amplitudeFileWriter::amplitudeFileWriter()
	: _initialized(false),
	  _outputFile(0),
	  _metadata(),
	  _ampTreeLeaf(0),
	  _hashCalculator()
{

}


rpwa::amplitudeFileWriter::~amplitudeFileWriter()
{
	reset();
}


bool rpwa::amplitudeFileWriter::initialize(TFile&                       outputFile,
                                           const vector<eventMetadata*> eventMeta,
                                           const string&                keyfileContent,
                                           const string&                objectBaseName,
                                           const int&                   splitlevel,
                                           const int&                   buffsize)
{
	if(_initialized) {
		printWarn << "trying to initialized when already initialized." << endl;
		return false;
	}

	_outputFile = &outputFile;
	_outputFile->cd();

	vector<eventMetadata> eventMetaObjects;
	for(unsigned int i = 0; i < eventMeta.size(); ++i) {
		eventMetaObjects.push_back(*(eventMeta[i]));
	}
	_metadata.setEventMetadata(eventMetaObjects);
	_metadata.setKeyfileContent(keyfileContent);
	_metadata.setRootpwaGitHash(gitHash());
	_metadata.setObjectBaseName(objectBaseName);

	const string treeName = amplitudeMetadata::getObjectNames(objectBaseName).first;

	_metadata._amplitudeTree = new TTree(treeName.c_str(), treeName.c_str());
	_ampTreeLeaf = new rpwa::amplitudeTreeLeaf();
	_metadata._amplitudeTree->Branch(rpwa::amplitudeMetadata::amplitudeLeafName.c_str(), &_ampTreeLeaf, buffsize, splitlevel);

	_initialized = true;
	return _initialized;
}


void rpwa::amplitudeFileWriter::addAmplitude(const complex<double>& amplitude)
{
	if(not _initialized) {
		printWarn << "trying to add amplitude when not initialized." << endl;
		return;
	}
	_ampTreeLeaf->setAmp(amplitude);
	_hashCalculator.Update(amplitude);
	_metadata._amplitudeTree->Fill();
}


void rpwa::amplitudeFileWriter::addAmplitudes(const vector<complex<double> >& amplitudes)
{
	for(unsigned int i = 0; i < amplitudes.size(); ++i) {
		addAmplitude(amplitudes[i]);
	}
}


void rpwa::amplitudeFileWriter::reset()
{
	if(_ampTreeLeaf) {
		delete _ampTreeLeaf;
		_ampTreeLeaf = 0;
	}
	_outputFile = 0;
	_hashCalculator = hashCalculator();
	_initialized = false;
}


bool rpwa::amplitudeFileWriter::finalize()
{
	if(not _initialized) {
		printWarn << "trying to finalize when not initialized." << endl;
		return false;
	}
	_metadata.setContentHash(_hashCalculator.hash());
	_outputFile->cd();
	_metadata.Write(_metadata.getObjectNames().second.c_str());
	reset();
	return true;
}
