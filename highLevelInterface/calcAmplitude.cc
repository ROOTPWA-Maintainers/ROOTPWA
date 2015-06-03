
#include "calcAmplitude.h"

#include <boost/progress.hpp>

#include <TClonesArray.h>
#include <TTree.h>
#include <TTreePerfStats.h>

#include <reportingUtils.hpp>

using namespace std;
using namespace rpwa;


vector<complex<double> >
rpwa::hli::calcAmplitude(eventMetadata&            eventMeta,
                         const isobarAmplitudePtr& amplitude,
                         const long int            maxNmbEvents,
                         const bool                printProgress,
                         const string&             treePerfStatOutFileName,         // root file name for tree performance result
                         const long int            treeCacheSize)
{
	vector<complex<double> > retval;

	if(not amplitude) {
		printWarn << "null pointer to isobar decay amplitude. cannot process tree." << endl;
		return retval;
	}
	// initialize amplitude
	amplitude->init();
	const isobarDecayTopologyPtr& decayTopo = amplitude->decayTopology();

	TTree* tree = eventMeta.eventTree();
	if(not tree) {
		printErr << "event tree not found." << endl;
		return retval;
	}

	// create branch pointers and leaf variables
	TBranch*      prodKinMomentaBr  = 0;
	TBranch*      decayKinMomentaBr = 0;
	TClonesArray* prodKinMomenta    = 0;
	TClonesArray* decayKinMomenta   = 0;

	// connect leaf variables to tree branches
	tree->SetBranchAddress(eventMetadata::productionKinematicsMomentaBranchName.c_str(),  &prodKinMomenta,  &prodKinMomentaBr );
	tree->SetBranchAddress(eventMetadata::decayKinematicsMomentaBranchName.c_str(), &decayKinMomenta, &decayKinMomentaBr);
	tree->SetCacheSize(treeCacheSize);
	tree->AddBranchToCache(eventMetadata::productionKinematicsMomentaBranchName.c_str(),  true);
	tree->AddBranchToCache(eventMetadata::decayKinematicsMomentaBranchName.c_str(), true);
	tree->StopCacheLearningPhase();
	TTreePerfStats* treePerfStats = 0;
	if(treePerfStatOutFileName != "") {
		treePerfStats = new TTreePerfStats("ioPerf", tree);
	}

	// loop over events
	if(not decayTopo->initKinematicsData(eventMeta.productionKinematicsParticleNames(), eventMeta.decayKinematicsParticleNames())) {
		printWarn << "problems initializing input data. cannot read input data." << endl;
		return retval;
	}
	const long int    nmbEventsTree     = tree->GetEntries();
	const long int    nmbEvents         = ((maxNmbEvents > 0) ? min(maxNmbEvents, nmbEventsTree)
	                                       : nmbEventsTree);
	boost::progress_display* progressIndicator = (printProgress) ? new boost::progress_display(nmbEvents, cout, "") : 0;
	for (long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
		if(progressIndicator) {
			++(*progressIndicator);
		}

		if(tree->LoadTree(eventIndex) < 0) {
			break;
		}
		// read only required branches
		prodKinMomentaBr->GetEntry (eventIndex);
		decayKinMomentaBr->GetEntry(eventIndex);

		if(not prodKinMomenta or not decayKinMomenta) {
			printWarn << "at least one of the input data arrays is a null pointer: "
			          << "        production kinematics: " << "momenta = " << prodKinMomenta  << endl
			          << "        decay kinematics:      " << "momenta = " << decayKinMomenta << endl
			          << "skipping event." << endl;
			return vector<complex<double> >();
		}

		if(decayTopo->readKinematicsData(*prodKinMomenta, *decayKinMomenta)) {
			retval.push_back((*amplitude)());
		} else {
			printWarn << "problems reading event[" << eventIndex << "]" << endl;
			return vector<complex<double> >();
		}
	}

	if(printProgress) {
		tree->PrintCacheStats();
	}
	if(treePerfStats) {
		treePerfStats->SaveAs(treePerfStatOutFileName.c_str());
		delete treePerfStats;
	}
	return retval;
}
