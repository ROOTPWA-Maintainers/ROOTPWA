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

#include "Vec.h"
#include "lorentz.h"

#include "mathUtils.hpp"
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
	pdt.readFile();
	if (0)
		printInfo << "particle data table:" << endl
		          << pdt;

	// test filling of particle properties
	if (1) {
		particleProperties partProp;
		const string       partName = "pi";
		partProp.fillFromDataTable(partName);
		printInfo << "particle properties for '" << partName << "':" << endl
		          << partProp << endl;
		pdt.addEntry(partProp);
	}

	// test construction of particles
	if (1) {
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
	if (1) {
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

	if (0) {
		{
			fourVec  p(2, threeVec(0.5, 0.75, 1));
			threeVec n = threeVec(0, 0, 1) / p.V();
			cout << "before = " << n << "    " << p << endl;
			rotation         R;
			lorentzTransform L1;
			L1.set(R.set(n.phi(), n.theta() - M_PI / 2, -M_PI / 2));
			n *= R;
			p *= L1;
			cout << "L1 -> " << n << "    " << p << endl;
			lorentzTransform L2;
			L2.set(R.set(0, signof(p.x()) * p.theta(), 0));
			p *= L2;
			cout << "L2 -> " << p << endl;
			lorentzTransform L3;
			L3.set(p);
			p *= L3;
			cout << "L3 -> " << p << endl;

			matrix<double> X(4, 4);
			X = L3 * (L2 * L1);
			lorentzTransform L(X);
			p = fourVec(2, threeVec(0.5, 0.75, 1));
			p *= L;
			cout << "L -> " << p << endl;
		}

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
