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
// $Rev:: 462                         $: revision of last commit
// $Author:: bgrube                   $: author of last commit
// $Date:: 2010-11-03 19:15:35 +0100 #$: date of last commit
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

#include "svnVersion.h"
#include "utilities.h"
#include "particleDataTable.h"
#include "waveDescription.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
	waveDescription::setDebug(true);

	printCompilerInfo();
	printSvnVersion();
	particleDataTable::readFile();

	if (1) {
		const string       keyFileName = "test.key";
		waveDescription    waveDesc;
		isobarAmplitudePtr amp;
		if (waveDesc.parseKeyFile(keyFileName) and waveDesc.constructAmplitude(amp)) {
			isobarDecayTopologyPtr topo = amp->decayTopology();
			printInfo << "key file:" << endl << waveDesc.keyFileContents();
			printInfo << *amp;
			topo->writeGraphViz("testAmplitude.dot");
			gSystem->Exec("dot -Tps -o testAmplitude.ps testAmplitude.dot");
			waveDesc.writeKeyFile("testWrite.key", *amp);  // test key file creation
			// test file I/O of waveDescription
			const string waveName = waveDesc.waveNameFromTopology(*topo);
			{
				TFile* outFile = TFile::Open("testWaveDescription.root", "RECREATE");
				waveDesc.Write(waveName.c_str());
				outFile->Close();
			}
			{
				TFile*           inFile    = TFile::Open("testWaveDescription.root", "READ");
				waveDescription* waveDesc2 = 0;
				inFile->GetObject(waveName.c_str(), waveDesc2);
				if (not waveDesc2)
					printErr << "cannot find wave description '" << waveName << "'" << endl;
				else
					printInfo << "key file:" << endl << waveDesc2->keyFileContents();
				inFile->Close();
			}
		}
	}

}
