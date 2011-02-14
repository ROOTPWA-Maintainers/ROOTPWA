//
// converts uDst tree to other formats
//


#include <string>
#include <vector>
#include <libgen.h>  // dirname()
#ifdef basename
#undef basename
#endif
#include <string.h>  // basename()

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TObject.h"

#include "reportingUtils.hpp"
#include "fileUtils.hpp"
#include "TFitResult.h"
#include "fitResult.h"


using namespace std;
using namespace rpwa;


#if TFITRESULT_ENABLED


void
convertTFitResultTree(const string&           inFileNamePattern = "./*.root",
                      const string&           outDirName        = "",  // if empty, files are modified in place
                      const unsigned long int maxNmbEvents      = -1,
                      const string&           inTreeName        = "pwa",
                      const string&           outTreeName       = "pwa",
                      const string&           inBranchName      = "fitResult",
                      const string&           outBranchName     = "fitResult_v2")
{ 

	// modify files in place, if outDirName is not specified
	string _outDirName = outDirName;
	if (outDirName == "") {
		char* s = strdup(inFileNamePattern.c_str());  // dummy required by dirname()
		_outDirName = dirname(s);
	}

	printInfo << "converting TFitResult trees in '" << inFileNamePattern << "' "
	          << "into fitResult trees." << endl;
	cout << "    writing to '" << _outDirName << "'" << endl;

	// get file list
	const vector<string> inFileNames = filesMatchingGlobPattern(inFileNamePattern);
	const unsigned int   nmbInFiles  = inFileNames.size();
	cout << "    found " << nmbInFiles << " files." << endl;
	if (inFileNames.size() == 0) {
		printErr << "no files do match pattern '" << inFileNamePattern << "'. exiting." << endl;
		return;
	}

	// process file list
	unsigned long int countEntriesWritten = 0;
	for (unsigned int fileIndex = 0; fileIndex < nmbInFiles; ++fileIndex) {
		const string outFileName = _outDirName + "/" + basename(inFileNames[fileIndex].c_str());
		const bool   updateMode  = (outFileName == inFileNames[fileIndex]) ? true : false;

		// open input file
		cout << "    reading from input file '" << inFileNames[fileIndex] << "' "
		     << "[" << fileIndex << "/" << nmbInFiles << "]" << endl;
		TFile* inFile = TFile::Open(inFileNames[fileIndex].c_str(), (updateMode) ? "UPDATE" : "READ");
		if (!inFile || inFile->IsZombie()) {
			printWarn << "cannot open file '" << inFileNames[fileIndex] << "'. skipping." << endl;
			continue;
		}

		// create output file, if necessary
		cout << "    writing to output file '" << outFileName << "'" << endl;
		TFile* outFile = 0;
		if (updateMode)
			outFile = inFile;
		else
			outFile = TFile::Open(outFileName.c_str(), "RECREATE");
		if (!outFile || outFile->IsZombie()) {
			printWarn << "cannot open file '" << outFileName << "'. skipping." << endl;
			continue;
		}

		// get input tree
		TTree* inTree = 0;
		inFile->GetObject(inTreeName.c_str(), inTree);
		if (!inTree) {
			printWarn << "cannot find tree '" << inTreeName << "'. skipping." << endl;
			continue;
		}
		TFitResult* inResult = 0;
		TBranch*    inBranch = 0;
		inTree->SetBranchAddress(inBranchName.c_str(), &inResult, &inBranch);
		// renaming old branch does not work for some reason; same name leafs are getting mixed up
		// if (updateMode && (outBranchName == inBranchName))
		//   inBranch->SetName((inBranchName + "_old").c_str());

		// create output tree
		outFile->cd();

		TTree* outTree = 0;
		if (updateMode)
			outTree = inTree;
		else
			//outTree = new TTree(outTreeName.c_str(), "converted from TFitResult");
			outTree = new TTree(outTreeName.c_str(), outTreeName.c_str());
		fitResult* outResult = 0;
		string _outBranchName = outBranchName;
		if (updateMode && (outBranchName == inBranchName))
			_outBranchName += "_v2";
		if (outTree->FindBranch(_outBranchName.c_str())) {
			printWarn << "branch '" << _outBranchName << "' already exists in output file. "
			          << "skipping." << endl;
			continue;
		}
		TBranch* outBranch = outTree->Branch(_outBranchName.c_str(), &outResult);

		// loop over events and fill output tree
		const unsigned long int nmbEvents = inTree->GetEntriesFast();
		for (unsigned long int eventIndex = 0; eventIndex < nmbEvents; ++eventIndex) {
			inTree->GetEntry(eventIndex);
			if (!inResult) {
				printWarn << "zero pointer for TFitResult at index " << eventIndex << ". skipping." << endl;
				continue;
			}
			// copy data
			outResult = new fitResult(*inResult);
			if (updateMode && (outTreeName == inTreeName))
				outBranch->Fill();
			else
				outTree->Fill();
			delete outResult;
			outResult = 0;
			++countEntriesWritten;
			if (countEntriesWritten >= maxNmbEvents)
				break;
		}

		outTree->Write("", TObject::kOverwrite);
		outFile->Close();
		delete outFile;
		if (inFile != outFile) {
			inFile->Close();
			delete inFile;
		}
    
		if (countEntriesWritten >= maxNmbEvents)
			break;
	}
	printInfo << "wrote " << countEntriesWritten << " entries to file." << endl;

}


#endif  // TFITRESULT_ENABLED
