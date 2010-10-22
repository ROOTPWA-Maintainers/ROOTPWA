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
//      compare math functions to libpp implementations
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
#include <complex>

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include "TRandom3.h"
#include "TStopwatch.h"

#include "svnVersion.h"
#include "utilities.h"
#include "factorial.hpp"
#include "clebschGordanCoeff.hpp"
#include "dFunction.hpp"


using namespace std;
using namespace rpwa;
using namespace boost::math;


int
main(int argc, char** argv)
{
	printCompilerInfo();
	printSvnVersion  ();

	//////////////////////////////////////////////////////////////////////////////
	// factorial
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
				//newFact[i] = rpwa::factorial<unsigned int>(data[i]);
				//newFact[i] = rpwa::factorial<float>(data[i]);
		    newFact[i] = rpwa::factorial<double>(data[i]);
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


	//////////////////////////////////////////////////////////////////////////////
	// Clebsch-Gordan coefficients (j1 m1 j2 m2 | J M)
	if (1) {
		printInfo << "testing Clebsch-Gordan coefficients" << endl;

		const unsigned int nmbIterations = 100;
		const int          maxJ          = 4;

    // determine size of data array
		unsigned int nmbVals = maxJ * (2 * maxJ - 1);
		nmbVals = nmbVals * nmbVals * 2 * maxJ * (4 * maxJ - 1);

    // compute mathUtils values
    clebschGordanCoeffCached<double>::instance().setUseCache(true);
    TStopwatch timer;
    timer.Reset();
    timer.Start();
    vector<double> newVals(nmbVals, 0);
    for (unsigned int it = 0; it < nmbIterations; ++it) {
	    unsigned int valIndex = 0;
			for (int j1 = 0; j1 < maxJ; ++j1)
				for (int m1 = -j1; m1 <= j1; ++m1)
					for (int j2 = 0; j2 < maxJ; ++j2)
						for (int m2 = -j2; m2 <= j2; ++m2)
							for (int J = 0; J < 2 * maxJ; ++J)
								for (int M = -J; M <= J; ++M) {
									newVals[valIndex] =
										clebschGordanCoeff<double>(2 * j1, 2 * m1, 2 * j2, 2 * m2, 2 * J, 2 * M);
									++valIndex;
						    }
    }
    timer.Stop();
    printInfo << "calculated mathUtil Clebsch-Gordan coefficients for " << newVals.size()
              << " x " << nmbIterations << " calls " << endl << "    this consumed: ";
    timer.Print();
    printInfo << "size of cache is "
              << clebschGordanCoeffCached<double>::instance().cacheSize() / (1024. * 1024.)
              << " MBytes" << endl;

    // compute libpp values
    timer.Reset();
    timer.Start();
    vector<double> oldVals(nmbVals, 0);
    for (unsigned int it = 0; it < nmbIterations; ++it) {
	    unsigned int valIndex = 0;
	    for (int j1 = 0; j1 < maxJ; ++j1)
		    for (int m1 = -j1; m1 <= j1; ++m1)
			    for (int j2 = 0; j2 < maxJ; ++j2)
				    for (int m2 = -j2; m2 <= j2; ++m2)
					    for (int J = 0; J < 2 * maxJ; ++J)
						    for (int M = -J; M <= J; ++M) {
							    oldVals[valIndex] = clebsch(2 * j1, 2 * j2, 2 * J, 2 * m1, 2 * m2, 2 * M);
							    ++valIndex;
						    }
    }
    timer.Stop();
    printInfo << "calculated libpp Clebsch-Gordan coefficients for " << oldVals.size()
              << " x " << nmbIterations << " calls " << endl << "    this consumed: ";
    timer.Print();

    // check values
    double maxAbsDeviation = 0;
    double maxRelDeviation = 0;
    for (unsigned int i = 0; i < nmbVals; ++i) {
	    const double absDelta = oldVals[i] - newVals[i];
	    const double relDelta = absDelta / oldVals[i];
	    if (abs(absDelta) > maxAbsDeviation)
		    maxAbsDeviation = abs(absDelta);
	    if ((newVals[i] != 0) && (abs(relDelta) > maxRelDeviation))
		    maxRelDeviation = abs(relDelta);
    }
    printInfo << "maximum absolute deviation is " << maxAbsDeviation << "; "
              << "maximum relative deviation is " << maxRelDeviation << endl;
	}
	

	//////////////////////////////////////////////////////////////////////////////
	// Wigner d-function d^j_{m n}(theta)
	if (0) {
		printInfo << "testing Wigner d-function" << endl;

		const unsigned int nmbAngles = 50000;
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
    //dFunctionCached<double>::instance().setUseCache(false);
    TStopwatch timer;
    timer.Reset();
    timer.Start();
    unsigned int   valIndex = 0;
    vector<double> newVals(nmbVals, 0);
    for (int j = 0; j < maxJ; ++j)
	    for (int m = -j; m <= j; ++m)
		    for (int n = -j; n <= j; ++n)
			    for (unsigned int i = 0; i < angles.size(); ++i) {
				    newVals[valIndex] = dFunction(2 * j, 2 * m, 2 * n, angles[i]);
				    ++valIndex;
			    }
    timer.Stop();
    printInfo << "calculated mathUtil d-Functions for " << newVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();
    printInfo << "size of cache is "
              << dFunctionCached<double>::instance().cacheSize() / (1024. * 1024.) << " MBytes"
              << endl;

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
	

	//////////////////////////////////////////////////////////////////////////////
	// spherical harmonics Y_l^m(theta, phi)
	if (0) {
		printInfo << "testing spherical harmonics" << endl;

		const unsigned int nmbAngles = 150000;
		const int          maxL      = 10;

		vector<double> angles[2];
		angles[0].resize(nmbAngles, 0);
		angles[1].resize(nmbAngles, 0);
		TRandom3 random(1234567890);
		for (unsigned int i = 0; i < nmbAngles; ++i) {
			angles[0][i] = (random.Uniform(-piHalf, +piHalf));  // theta
			angles[1][i] = (random.Uniform(-pi,     +pi    ));  // phi
		}
		
    // determine size of data array
    unsigned int nmbVals = 0;
    for (int l = 0; l < maxL; ++l)
	    for (int m = -l; m <= l; ++m)
		    for (unsigned int i = 0; i < angles[0].size(); ++i)
			    ++nmbVals;

    // compute mathUtils values
    //dFunctionCached<double>::instance().setUseCache(false);
    TStopwatch timer;
    timer.Reset();
    timer.Start();
    unsigned int             valIndex = 0;
    vector<complex<double> > myVals(nmbVals, 0);
    for (int l = 0; l < maxL; ++l)
	    for (int m = -l; m <= l; ++m)
		    for (unsigned int i = 0; i < angles[0].size(); ++i) {
			    myVals[valIndex] = sphericalHarmonic<complex<double> >
				                       (2 * l, 2 * m, angles[0][i], angles[1][i]);
			    ++valIndex;
		    }
    timer.Stop();
    printInfo << "calculated mathUtil spherical harmonics for " << myVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();
    printInfo << "size of cache is "
              << dFunctionCached<double>::instance().cacheSize() / (1024. * 1024.) << " MBytes"
              << endl;

    // compute Boost values
    timer.Reset();
    timer.Start();
    valIndex = 0;
    vector<complex<double> > boostVals(nmbVals, 0);
    for (int l = 0; l < maxL; ++l)
	    for (int m = -l; m <= l; ++m)
		    for (unsigned int i = 0; i < angles[0].size(); ++i) {
			    boostVals[valIndex] = spherical_harmonic(l, m, angles[0][i], angles[1][i]);
			    ++valIndex;
		    }
    timer.Stop();
    printInfo << "calculated Boost spherical harmonics for " << boostVals.size() << " angles" << endl
		          << "    this consumed: ";
    timer.Print();

    // check values
    complex<double> maxDeviation = 0;
    for (unsigned int i = 0; i < nmbVals; ++i) {
	    const complex<double> delta = boostVals[i] - myVals[i];
	    if (abs(delta) > abs(maxDeviation))
		    maxDeviation = abs(delta);
    }
    printInfo << "maximum deviation is " << maxDeviation << endl;
	}


	//////////////////////////////////////////////////////////////////////////////
	// Wigner D-function D^j_{m n}(alpha, beta, gamma)
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
    printInfo << "size of cache is "
              << dFunctionCached<double>::instance().cacheSize() / (1024. * 1024.) << " MBytes"
              << endl;

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


	//////////////////////////////////////////////////////////////////////////////
	// Wigner D-function in reflectivity basis {^refl}D^j_{m n}(alpha, beta, gamma)
	if (0) {
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
    printInfo << "size of cache is "
              << dFunctionCached<double>::instance().cacheSize() / (1024. * 1024.) << " MBytes"
              << endl;

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
						    const double    preFactor  = (_m == 0 ? 0.5 : 1 / std::sqrt(2));
						    const double    reflFactor = (double)refl * (double)P * std::pow(-1, 0.5 * (_j - _m));
						    oldVals[valIndex] = preFactor * (               D(a.X(), a.Y(), a.Z(), _j,  _m, _n)
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
