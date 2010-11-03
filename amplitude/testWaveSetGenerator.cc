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
#include "keyFileParser.h"
#include "waveSetGenerator.h"


using namespace std;
using namespace boost::assign;
using namespace rpwa;


int
main(int argc, char** argv)
{
	// switch on debug output
	//particleProperties::setDebug(true);
	//particleDataTable::setDebug(true);
	particle::setDebug(true);
	//decayTopologyGraphType::setDebug(true);
	//isobarDecayVertex::setDebug(true);
	//decayTopology::setDebug(true);
	//isobarDecayTopology::setDebug(true);
	waveSetGenerator::setDebug(true);

	particleDataTable& pdt = particleDataTable::instance();
	pdt.readFile();

	if (0) {
		particleProperties prop1("bla",  2, -1, 0, +1, +1);
		particleProperties prop2("blub", 2, +1, 2, -1, -1);
  
		string opt="IGJPC";
 
		cout << "Comparison  result: " << (prop1 == prop2) << endl;
		cout << "Comparison with opt = " << opt << "  result: "
		     << (prop1 ==  pair<particleProperties, string>(prop2, opt)) << endl;

		vector<const particleProperties*> selection = pdt.entriesMatching(prop2, opt);
		cout << "Matching entries in pdt with prototype: " << endl;
		cout << prop2 << endl;
		cout << " with option " << opt << endl;

		for(unsigned int i = 0; i < selection.size(); ++i)
			cout << *selection[i] << endl;
	}

	if (1) {
		waveSetGenerator waveSetGen;
		const vector<string> isobarWhiteList = list_of("sigma")("rho(770)")("f2(1270)")
			("pi(1300)")("a1(1260)")("a2(1320)")("pi2(1670)")
			("f0(1500)")("f1(1285)")("rho(1450)")("rho(1700)")("rho3(1690)");
		waveSetGen.setIsobarWhiteList(isobarWhiteList);
		cout << waveSetGen;
		waveSetGen.generateWaveSet("testWaveSetGenerator.key");
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
			cout << setw(4) << i << ": " << keyFileParser::keyFileNameFromTopology(decays[i]) << endl;
		cout << "got " << inconsistentDecays << " inconsistent" << endl
		     << "and " << consistentDecays << " valid decays" << endl
		     << "out of " << decays.size() << " constructed decays" << endl;

		waveSetGen.writeKeyFiles("testWaveSetGenerator");
		
		// decays.back().writeGraphViz("foo.dot");
		// gSystem->Exec("dot -Tps -o foo.ps foo.dot");
		// cout << parser.keyFileNameFromTopology(decays.back()) << endl;
	}
}
