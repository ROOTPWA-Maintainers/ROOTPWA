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
//      basic test program for vertex and decay topology
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <boost/assign/list_of.hpp>

#include "TSystem.h"

#include "particleDataTable.h"
#include "waveDescription.h"
#include "waveSetGenerator.h"


using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace rpwa;


int
main(int argc, char** argv)
{
	// switch on debug output
	//particleProperties::setDebug(true);
	//particleDataTable::setDebug(true);
	//particle::setDebug(true);
	//decayTopologyGraphType::setDebug(true);
	//isobarDecayVertex::setDebug(true);
	//decayTopology::setDebug(true);
	//isobarDecayTopology::setDebug(true);
	waveSetGenerator::setDebug(true);

	particleDataTable& pdt = particleDataTable::instance();
	pdt.readFile();

	// test selective comparision of particle properties
	if (0) {
		particleProperties prop1("bla",  2, -1, 0, +1, +1);
		particleProperties prop2("blub", 2, +1, 2, -1, -1);

		string opt="IGJPC";

		printDebug << "Comparison  result: " << (prop1 == prop2) << endl;
		printDebug << "Comparison with opt = " << opt << "  result: "
		           << (prop1 ==  pair<particleProperties, string>(prop2, opt)) << endl;

		vector<const particleProperties*> selection = pdt.entriesMatching(prop2, opt);
		printDebug << "Matching entries in pdt with prototype: " << endl
		           << prop2 << endl
		           << " with option " << opt << endl;

		for(unsigned int i = 0; i < selection.size(); ++i)
			cout << *selection[i] << endl;
	}


	// test particle comparison
	if (1) {
		particleProperties partProp1;
		partProp1.fillFromDataTable("pi+");
		particleProperties partProp2 = partProp1;

		particlePtr part1 = createParticle(partProp1);
		particlePtr part2 = createParticle(partProp2, 1, 2, +1);
		printDebug << "particle properties" << endl
		           << partProp1 << endl
		           << partProp2 << endl;
		printInfo << "comparing particle properties" << endl;
		partProp1 == partProp2;
		printDebug << "particles" << endl
		           << *part1 << endl
		           << *part2 << endl;
		printInfo << "comparing particles" << endl;
		*part1 == *part2;

		// const string       keyFileName = "../../keyfiles/key3pi/SET1_new/1-0-+0+rho770_11_pi-.key";
		const string       keyFileName = "testWaveDescription.key";
		waveDescription    waveDesc;
		isobarAmplitudePtr amp;
		if (waveDesc.parseKeyFile(keyFileName) and waveDesc.constructAmplitude(amp)) {
			isobarDecayTopologyPtr isoTopo = amp->decayTopology();
			decayTopologyPtr       topo    = static_pointer_cast<decayTopology>(isoTopo);
			// printInfo << *isoTopo << endl;
			// printInfo << *topo    << endl;
			isobarDecayVertexPtr isoVert = isoTopo->XIsobarDecayVertex();
			interactionVertexPtr vert    = topo->XDecayVertex();
			printInfo << "isobar decay vertex: " << *isoVert << endl;
			printInfo << "interaction vertex:  " << *vert    << endl;

			printInfo << "comparing isobar decay vertices" << endl;
			*isoVert == *isoVert;
			printInfo << "comparing interaction vertices" << endl;
			*vert == *vert;
		}

	}


	if (0) {
		waveSetGenerator waveSetGen;
		if (not waveSetGen.setWaveSetParameters("testWaveSetGenerator.key")) {
			cout << "could not initialize wave set generator. aborting." << endl;
			exit(1);
		}
		// const vector<string> isobarWhiteList = list_of("sigma")("rho(770)")("f2(1270)")
		// 	("pi(1300)")("a1(1260)")("a2(1320)")("pi2(1670)")
		// 	("f0(1500)")("f1(1285)")("rho(1450)")("rho(1700)")("rho3(1690)");
		// waveSetGen.setIsobarWhiteList(isobarWhiteList);
		cout << waveSetGen;
		waveSetGen.generateWaveSet();
		vector<isobarDecayTopology>& decays             = waveSetGen.waveSet();
		unsigned int                 consistentDecays   = 0;
		unsigned int                 inconsistentDecays = 0;
		for (unsigned int i = 0; i < decays.size(); ++i) {
			cout << decays[i];
			// decays[i].printPointers(cout);
			// for (decayTopologyGraphType::nodeIterator j = decays[i].nodes().first;
			// 	   j != decays[i].nodes().second; ++j)
			//   decays[i].vertex(*j)->printPointers(cout);
			//isobarDecayVertex::setDebug(true);
			bool isConsistent = decays[i].checkTopology() and decays[i].checkConsistency();
			//isobarDecayVertex::setDebug(false);
			if (isConsistent) {
				cout << "isobar decay topology is consistent" << endl;
				++consistentDecays;
			} else {
				cout << "isobar decay topology is NOT consistent" << endl;
				++inconsistentDecays;
			}
		}

		for (unsigned int i = 0; i < decays.size(); ++i)
			cout << setw(4) << i << ": " << waveDescription::waveNameFromTopology(decays[i]) << endl;

		gSystem->Exec("mkdir testWaveSetGenerator");
		waveSetGen.writeKeyFiles("testWaveSetGenerator");

		cout << "got " << inconsistentDecays << " inconsistent" << endl
		     << "and " << consistentDecays << " valid decays" << endl
		     << "out of " << decays.size() << " constructed decays" << endl;

		// decays.back().writeGraphViz("foo.dot");
		// gSystem->Exec("dot -Tps -o foo.ps foo.dot");
	}
}
