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
//      basic test program for wave description class
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "TFile.h"
#include "TSystem.h"

#include "particleDataTable.h"
#include "waveDescription.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
	waveDescription::setDebug(true);

	printCompilerInfo();
	printGitHash();
	particleDataTable::readFile();

	if (1) {
		const string       keyFileName = "testWaveDescription.key";
		waveDescription    waveDesc;
		isobarAmplitudePtr amp;
		if (waveDesc.parseKeyFile(keyFileName) and waveDesc.constructAmplitude(amp)) {
			isobarDecayTopologyPtr topo = amp->decayTopology();
			printInfo << *amp;
			topo->writeGraphViz("testWaveDescription.dot");
			gSystem->Exec("dot -Tps -o testWaveDescription.ps testWaveDescription.dot");
			waveDesc.writeKeyFile("testWaveDescriptionWrite.key", *amp);  // test key file creation
			//waveDesc.writeKeyFile("testWaveDescriptionWrite.key", *(amp->decayTopology()));  // test key file creation
			// test file I/O of waveDescription
			const string waveName = waveDesc.waveNameFromTopology(*topo);
			{
				TFile* outFile = TFile::Open("testWaveDescription.root", "RECREATE");
				waveDesc.Write(waveName.c_str());
				outFile->Close();
			}
			amp.reset();
			cout << endl
			     << "--------------------------------------------------------------------------------"
			     << endl << endl;
			{
				TFile*           inFile    = TFile::Open("testWaveDescription.root", "READ");
				waveDescription* waveDesc2 = 0;
				inFile->GetObject(waveName.c_str(), waveDesc2);
				if (not waveDesc2)
					printErr << "cannot find wave description '" << waveName << "'" << endl;
				else
					printInfo << "key file:" << endl;
				waveDesc2->printKeyFileContents(cout);
				waveDesc2->constructAmplitude(amp);
				waveDesc2->writeKeyFile("testWaveDescriptionWrite2.key", *amp);  // test key file creation
				//waveDesc2->writeKeyFile("testWaveDescriptionWrite2.key", *(amp->decayTopology()));  // test key file creation
				printInfo << *amp;
				inFile->Close();
			}
		}
	}

}
