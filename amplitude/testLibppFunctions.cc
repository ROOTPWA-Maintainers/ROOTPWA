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
//      basic test program for amplitude classes
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <vector>
#include <cstdlib>
#include <ctime>

#include "TRandom3.h"
#include "TStopwatch.h"

#include "svnVersion.h"
#include "utilities.h"
#include "factorial.hpp"
#include "dFunction.hpp"


using namespace std;
using namespace rpwa;


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printSvnVersion();

	if (0) {
		printInfo << "testing factorial" << endl;

		const unsigned int nmbIterations = 5;
		const unsigned int nmbData       = 10000000;
		//const unsigned int maxArg        = 13;  // 12! is maximum for unsigned int
		//const unsigned int maxArg        = 35;  // 34! is maximum for double
		const unsigned int maxArg        = 171;  // 170! is maximum for double

		vector<unsigned int> data  (nmbData);
		TRandom3             random(1234567890);
		for (unsigned int i = 0; i < nmbData; ++i)
			data[i] = (unsigned int)(random.Uniform() * maxArg);

		TStopwatch timer;
		timer.Reset();
		timer.Start();
		//vector<unsigned int> newFact(nmbData);
		//vector<float> newFact(nmbData);
		vector<double> newFact(nmbData);
		for (unsigned int it = 0; it < nmbIterations; ++it)
			for (unsigned int i = 0; i < nmbData; ++i)
				//newFact[i] = rpwa::factorial<unsigned int>::instance()(data[i]);
				//newFact[i] = rpwa::factorial<float>::instance()(data[i]);
		    newFact[i] = rpwa::factorial<double>::instance()(data[i]);
		timer.Stop();
		printInfo << "calculated mathUtil factorials for " << nmbIterations * nmbData << " values" << endl
		          << "    this consumed: ";
    timer.Print();

		timer.Reset();
		timer.Start();
		//vector<unsigned int> oldFact(nmbData);
		vector<double> oldFact(nmbData);
		for (unsigned int it = 0; it < nmbIterations; ++it)
			for (unsigned int i = 0; i < nmbData; ++i)
				//oldFact[i] = fact(data[i]);
				oldFact[i] = dfact(data[i]);
		timer.Stop();
		printInfo << "calculated libpp factorials for " << nmbIterations * nmbData << " values" << endl
		          << "    this consumed: ";
    timer.Print();

		bool success = true;
		for (unsigned int i = 0; i < nmbData; ++i)
			if (newFact[i] != oldFact[i]) {
				printWarn << "wrong result " << data[i] << "! is " << oldFact[i] << " not "
				          << newFact[i] << endl;
				success = false;
			}
		if (success)
			printInfo << "all factorials okay." << endl;
		else
			printInfo << "there were errors." << endl;
	}	


	if (0) {
		printInfo << "testing Wigner d-function" << endl;

		const unsigned int nmbAngles = 50000;
		// const unsigned int nmbAngles = 500;
		const int          maxJ      = 7;  // for larger values libpp implementation gives wrong results

		vector<double> angles(nmbAngles, 0);
		TRandom3       random(1234567890);
		for (unsigned int i = 0; i < nmbAngles; ++i)
			angles[i] = (random.Uniform(-piHalf, +piHalf));
		
    // determine size of data array
    unsigned int nmbVals = 0;
    for (int j = 0; j < maxJ; ++j)
	    for (int m = -j; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (unsigned int i = 0; i < angles.size(); ++i)
				    ++nmbVals;

    // compute mathUtils values
    //dFunction<double>::instance().setUseCache(false);
    TStopwatch timer;
    timer.Reset();
    timer.Start();
    unsigned int   valIndex = 0;
    vector<double> newVals(nmbVals, 0);
    for (int j = 0; j < maxJ; ++j)
	    for (int m = -j; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (unsigned int i = 0; i < angles.size(); ++i) {
				    newVals[valIndex] = dFunction<double>::instance()(2 * j, 2 * m, 2 * n, angles[i]);
				    ++valIndex;
			    }
    timer.Stop();
    printInfo << "calculated mathUtil d-Functions for " << newVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();

    // compute libpp values
    timer.Reset();
    timer.Start();
    valIndex = 0;
    vector<double> oldVals(nmbVals, 0);
    for (int j = 0; j < maxJ; ++j)
	    for (int m = -j; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (unsigned int i = 0; i < angles.size(); ++i) {
				    oldVals[valIndex] = d_jmn_b(2 * j, 2 * m, 2 * n, angles[i]);
				    ++valIndex;
			    }
    timer.Stop();
    printInfo << "calculated libpp d-Functions for " << oldVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();

    // check values
    double maxAbsDeviation = 0;
    double maxRelDeviation = 0;
    for (unsigned int i = 0; i < nmbVals; ++i) {
	    const double absDelta = oldVals[i] - newVals[i];
	    const double relDelta = absDelta / oldVals[i];
	    if (abs(absDelta) > maxAbsDeviation)
		    maxAbsDeviation = abs(absDelta);
	    if (abs(relDelta) > maxRelDeviation)
		    maxRelDeviation = abs(relDelta);
    }
    printInfo << "maximum absolute deviation is " << maxAbsDeviation << "; "
              << "maximum relative deviation is " << maxRelDeviation << endl;
	}
	

	if (0) {
		printInfo << "testing Wigner D-function" << endl;

		const unsigned int nmbAngles = 50000;
		const int          maxJ      = 7;  // for larger values libpp implementation gives wrong results

		vector<TVector3> angles(nmbAngles);
		TRandom3         random(1234567890);
		for (unsigned int i = 0; i < nmbAngles; ++i) {
			angles[i].SetX(random.Uniform(-pi,     +pi    ));
			angles[i].SetY(random.Uniform(-piHalf, +piHalf));
			angles[i].SetZ(random.Uniform(-pi,     +pi    ));
		}
		
    // determine size of data array
    unsigned int nmbVals = 0;
    for (int j = 0; j < maxJ; ++j)
	    for (int m = -j; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (unsigned int i = 0; i < angles.size(); ++i)
				    ++nmbVals;

    // compute mathUtils values
    TStopwatch timer;
    timer.Reset();
    timer.Start();
    unsigned int             valIndex = 0;
    vector<complex<double> > newVals(nmbVals, 0);
    for (int j = 0; j < maxJ; ++j)
	    for (int m = -j; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (unsigned int i = 0; i < angles.size(); ++i) {
				    const TVector3& a = angles[i];
				    // newVals[valIndex] = DFunction<complex<double> >
					  //                       (2 * j, 2 * m, 2 * n, a.X(), a.Y(), a.Z());
				    newVals[valIndex] = DFunctionConj<complex<double> >
					                        (2 * j, 2 * m, 2 * n, a.X(), a.Y(), a.Z());
				    ++valIndex;
			    }
    timer.Stop();
    printInfo << "calculated mathUtil D-Functions for " << newVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();

    // compute libpp values
    timer.Reset();
    timer.Start();
    valIndex = 0;
    vector<complex<double> > oldVals(nmbVals, 0);
    for (int j = 0; j < maxJ; ++j)
	    for (int m = -j; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (unsigned int i = 0; i < angles.size(); ++i) {
				    const TVector3& a = angles[i];
				    // oldVals[valIndex] = D(a.X(), a.Y(), a.Z(), 2 * j, 2 * m, 2 * n);
				    oldVals[valIndex] = conj(D(a.X(), a.Y(), a.Z(), 2 * j, 2 * m, 2 * n));
				    ++valIndex;
			    }
    timer.Stop();
    printInfo << "calculated libpp D-Functions for " << oldVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();

    // check values
    complex<double> maxDeviation = 0;
    for (unsigned int i = 0; i < nmbVals; ++i) {
	    const complex<double> delta = oldVals[i] - newVals[i];
	    if (abs(delta) > abs(maxDeviation))
		    maxDeviation = abs(delta);
    }
    printInfo << "maximum deviation is " << maxDeviation << endl;
	}


	if (1) {
		printInfo << "testing Wigner D-function in reflectivity basis" << endl;

		const unsigned int nmbAngles = 10000;
		const int          maxJ      = 7;  // for larger values libpp implementation gives wrong results

		vector<TVector3> angles(nmbAngles);
		TRandom3         random(1234567890);
		for (unsigned int i = 0; i < nmbAngles; ++i) {
			angles[i].SetX(random.Uniform(-pi,     +pi    ));
			angles[i].SetY(random.Uniform(-piHalf, +piHalf));
			angles[i].SetZ(random.Uniform(-pi,     +pi    ));
		}
		
    // determine size of data array
    unsigned int nmbVals = 0;
    for (int j = 0; j < maxJ; ++j)
	    for (int m = 0; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (int P = -1; P <= 1; P += 2)
				    for (int refl = -1; refl <= 1; refl += 2)
					    for (unsigned int i = 0; i < angles.size(); ++i)
						    ++nmbVals;

    // compute mathUtils values
    TStopwatch timer;
    timer.Reset();
    timer.Start();
    unsigned int             valIndex = 0;
    vector<complex<double> > newVals(nmbVals, 0);
    for (int j = 0; j < maxJ; ++j)
	    for (int m = 0; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (int P = -1; P <= 1; P += 2)
				    for (int refl = -1; refl <= 1; refl += 2)
					    for (unsigned int i = 0; i < angles.size(); ++i) {
						    const TVector3& a = angles[i];
						    // newVals[valIndex] = DFunctionRefl<complex<double> >
						    //                       (2 * j, 2 * m, 2 * n, P, refl, a.X(), a.Y(), a.Z());
						    newVals[valIndex] = DFunctionReflConj<complex<double> >
							                        (2 * j, 2 * m, 2 * n, P, refl, a.X(), a.Y(), a.Z());
						    ++valIndex;
					    }
    timer.Stop();
    printInfo << "calculated mathUtil D-Functions for " << newVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();

    // compute libpp values
    timer.Reset();
    timer.Start();
    valIndex = 0;
    vector<complex<double> > oldVals(nmbVals, 0);
    for (int j = 0; j < maxJ; ++j)
	    for (int m = 0; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (int P = -1; P <= 1; P += 2)
				    for (int refl = -1; refl <= 1; refl += 2)
					    for (unsigned int i = 0; i < angles.size(); ++i) {
						    const TVector3& a          = angles[i];
						    const int       _j         = 2 * j;
						    const int       _m         = 2 * m;
						    const int       _n         = 2 * n;
						    const double    preFactor  = (_m == 0 ? 0.5 : 1 / sqrt(2));
						    const double    reflFactor = (double)refl * (double)P * pow(-1, 0.5 * (_j - _m));
						    oldVals[valIndex] = preFactor * (                D(a.X(), a.Y(), a.Z(), _j,  _m, _n)
						                                      - reflFactor * D(a.X(), a.Y(), a.Z(), _j, -_m, _n));
						    oldVals[valIndex] = conj(oldVals[valIndex]);
						    ++valIndex;
					    }
    timer.Stop();
    printInfo << "calculated libpp D-Functions for " << oldVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();

    // check values
    complex<double> maxDeviation = 0;
    for (unsigned int i = 0; i < nmbVals; ++i) {
	    const complex<double> delta = oldVals[i] - newVals[i];
	    if (abs(delta) > abs(maxDeviation))
		    maxDeviation = abs(delta);
    }
    printInfo << "maximum deviation is " << maxDeviation << endl;
	}

}
