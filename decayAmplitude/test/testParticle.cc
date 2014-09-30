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
//
// Description:
//      basic test program for particle and related classes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include "TVector3.h"
#include "TLorentzRotation.h"

#include "spinUtils.hpp"
#include "reportingUtilsRoot.hpp"
#include "conversionUtils.hpp"
#include "particleDataTable.h"
#include "particle.h"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{

	if (1) {
		// switch on debug output
		particleProperties::setDebug(true);
		particleDataTable::setDebug(true);
		particle::setDebug(true);
	}

	// test loading of particle data table
	particleDataTable& pdt = particleDataTable::instance();
	pdt.readFile("./particleDataTable.txt");
	pdt.readDecayModeFile("./testParticleDecays.txt");

	if (1)
		printInfo << "particle data table:" << endl
		          << pdt;

	// test filling of particle properties
	if (0) {
		particleProperties partProp;
		const string       partName = "pi+";
		partProp.fillFromDataTable(partName);
		printInfo << "particle properties for '" << partName << "':" << endl
		          << partProp << endl;
		pdt.addEntry(partProp);
		printInfo << "antiparticle properties for '" << partName << "':" << endl
		          << partProp.antiPartProperties() << endl;
	}

	// test construction of particles
	if (0) {
		TVector3 mom;
		mom = TVector3(1, 2, 3);
		const particle p1("pi+", true, 0,  0, 0, mom);
		mom = TVector3(2, 3, 4);
		const particle p2("pi-", true, 1, -1, 0, mom);
		particle p3 = p2;
		p3.setName("X+");
		p3.setSpinProj(+1);
		printInfo << "created particles: " << endl
		          << p1 << endl
		          << p2 << endl
		          << p3 << endl;
	}

	// checking charge name handling
	if (0) {
		for (int i = -2; i < 3; ++i) {
			stringstream c;
			c << "pi";
			if (abs(i) > 1)
				c << abs(i);
			c << sign(i);
			//const particle p(c.str(), i);
			int q;
			const string n = particle::chargeFromName(c.str(), q);
			printInfo << c.str() << ": charge = " << q << ", bare name = " << n << endl;
			const particle p(c.str());
			cout << "name = " << p.name() << endl
			     << endl;
		}
	}

	// check decay modes
	if (1) {
		cout << endl;
		// particleProperties::decayMode decay1 = particleProperties::decayMode();
		// printInfo << "decay mode: " << decay1 << endl;
		particleProperties partProp;

		// printDebug << "testing anti-particle decay modes" << endl;
		// partProp.fillFromDataTable("rho(770)-");
		// printInfo << "particle properties: " << partProp << endl;
		// particleProperties antiPartProp = partProp.antiPartProperties(true);
		// printInfo << "anti-particle properties: " << antiPartProp << endl;

		printDebug << "testing equality operator for decay modes" << endl;
		partProp.fillFromDataTable("sigma");
		printInfo << "particle properties: " << partProp << endl;
		for (unsigned int i = 0; i < partProp.nmbDecays(); ++i)
			for (unsigned int j = i; j < partProp.nmbDecays(); ++j) {
				cout << "    {decay mode[" << i << "]: " << partProp.decayModes()[i]
				     << " vs. decay mode[" << j << "]: " << partProp.decayModes()[j]
				     << "}: " << trueFalse(partProp.decayModes()[i] == partProp.decayModes()[j]) << endl;
			}
	}

	// checking spin-exotic
	if (0) {
		printInfo << "testing spin-exotic tag" << endl;
		for (particleDataTable::iterator i = pdt.begin(); i != pdt.end(); ++i) {
			const particleProperties& prop = i->second;
			const bool jpc  = jpcIsExotic(prop.J(), prop.P(), prop.C());
			const bool igjp = igjpIsExotic(prop.isospin(), prop.G(), prop.J(), prop.P());
			cout << prop.name() << ": " << yesNo(jpc) << " vs. " << yesNo(igjp)
			     << ((jpc != igjp) ? " <<<" : "") << endl;
		}
		for (int J = 0; J < 4; ++J)
			for (int P = -1; P <= 1; P += 2)
				for (int C = -1; C <= 1; C += 2) {
					const bool jpc   = jpcIsExotic(2 * J, P, C);
					const bool igjp1 = igjpIsExotic(0,  C, 2 * J, P);
					const bool igjp2 = igjpIsExotic(2, -C, 2 * J, P);
					cout << J << sign(P) << sign(C) << ": " << yesNo(jpc) << " vs. " << yesNo(igjp1)
					     << ", " << yesNo(igjp2) << endl;
				}
	}

	if (0) {

		{
			TLorentzVector p(0.5, 0.75, 1, 2);
			TVector3       n = TVector3(0, 0, 1).Cross(p.Vect());
			TRotation R1;
			R1.RotateZ(-n.Phi());
			R1.RotateY(piHalf - n.Theta());
			R1.RotateZ(piHalf);
			n *= R1;
			p *= R1;
			cout << "R1 -> " << n << "    " << p << endl;
			// rotate about yHfAxis so that daughter momentum is along z-axis
			TRotation R2;
			R2.RotateY(-signum(p.X()) * p.Theta());
			p *= R2;
			cout << "R2 -> " << p << endl;
			// boost into daughter RF
			TLorentzRotation L3(-p.BoostVector());
			cout << "L3 -> " << L3 * p << endl;

			R1.Transform(R2);
			TLorentzRotation L(R1);
			L.Boost(-p.BoostVector());
			p = TLorentzVector(0.5, 0.75, 1, 2);
			p *= L;
			cout << "L -> " << p << endl;
		}
	}

}
