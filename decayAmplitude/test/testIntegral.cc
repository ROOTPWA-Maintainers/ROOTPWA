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
//      basic test program for amplitude integral matrix
//
//
// Author List:
//      Boris Grube    TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <algorithm>

#include <boost/assign/list_of.hpp>

#include "TROOT.h"
#include "TFile.h"

#include "reportingUtils.hpp"
#include "ampIntegralMatrix.h"


using namespace std;
using namespace boost::assign;
using namespace rpwa;


int
main(int    argc,
     char** argv)
{
	// switch on debug output
	ampIntegralMatrix::setDebug(true);

	if (0) {
		// test integral calculation
		// get file list
		// vector<string> binAmpFileNames = filesMatchingGlobPattern("/data/compass/hadronData/massBins/2004/Q3PiData/r481.trunk/1260.1300/PSPAMPS/SYM/*.amp", true);
		// sort(binAmpFileNames.begin(), binAmpFileNames.end());
		// make sure order of waves is the same as used by int
		vector<string> binAmpFileNames = list_of
			("1-0-+0+f0980_00_pi-.amp"   )
			("1-0-+0+rho770_11_pi-.amp"  )
			("1-0-+0+sigma_00_pi-.amp"   )
			("1-1++0+f21270_12_pi-.amp"  )
			("1-1++0+rho770_01_pi-.amp"  )
			("1-1-+0-rho770_11_pi-.amp"  )
			("1-1++0+rho770_21_pi-.amp"  )
			("1-1++0+sigma_10_pi-.amp"   )
			("1-1++1+f21270_12_pi-.amp"  )
			("1-1++1-rho770_01_pi-.amp"  )
			("1-1++1+rho770_01_pi-.amp"  )
			("1-1-+1-rho770_11_pi-.amp"  )
			("1-1-+1+rho770_11_pi-.amp"  )
			("1-1++1+rho770_21_pi-.amp"  )
			("1-1++1+sigma_10_pi-.amp"   )
			("1-2-+0+f21270_02_pi-.amp"  )
			("1-2++0-f21270_12_pi-.amp"  )
			("1-2-+0+f21270_22_pi-.amp"  )
			("1-2-+0+rho770_11_pi-.amp"  )
			("1-2++0-rho770_21_pi-.amp"  )
			("1-2-+0+rho770_31_pi-.amp"  )
			("1-2-+0+sigma_20_pi-.amp"   )
			("1-2-+1-f21270_02_pi-.amp"  )
			("1-2-+1+f21270_02_pi-.amp"  )
			("1-2++1-f21270_12_pi-.amp"  )
			("1-2++1+f21270_12_pi-.amp"  )
			("1-2-+1+f21270_22_pi-.amp"  )
			("1-2-+1+rho770_11_pi-.amp"  )
			("1-2++1+rho770_21_pi-.amp"  )
			("1-2-+1+rho770_31_pi-.amp"  )
			("1-2-+1+sigma_20_pi-.amp"   )
			("1-3++0+f21270_12_pi-.amp"  )
			("1-3++0+rho31690_03_pi-.amp")
			("1-3++0+rho770_21_pi-.amp"  )
			("1-3++1+f21270_12_pi-.amp"  )
			("1-3++1+rho31690_03_pi-.amp")
			("1-3++1+rho770_21_pi-.amp"  )
			("1-4-+0+rho770_31_pi-.amp"  )
			("1-4++1+f21270_32_pi-.amp"  )
			("1-4-+1+rho770_31_pi-.amp"  )
			("1-4++1+rho770_41_pi-.amp"  );
		for (size_t i = 0; i < binAmpFileNames.size(); ++i)
			binAmpFileNames[i] = "/data/compass/hadronData/massBins/2004/Q3PiData/r481.trunk/1260.1300/PSPAMPS/SYM/" + binAmpFileNames[i];
		vector<string> rootAmpFileNames;
		ampIntegralMatrix integral;
		integral.integrate(binAmpFileNames, rootAmpFileNames);
		integral.writeAscii("testIntegral2.int");
	}


	if (1) {
		// test I/O and copying
		ampIntegralMatrix integral;
		// ascii I/O
		integral.readAscii("testIntegral.int");
		ampIntegralMatrix integral2(integral);
		integral2.writeAscii("testIntegral2.int");
		// root I/O
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 0, 0)
		// force loading predefined std::complex dictionary
		gROOT->ProcessLine("#include <complex>");
#endif
		{
			TFile* outFile = TFile::Open("testIntegral.root", "RECREATE");
			printInfo << "writing integral to 'testIntegral.root'" << endl;
			integral.Write("integral");
			outFile->Close();
		}
		{
			TFile*             inFile    = TFile::Open("testIntegral.root", "READ");
			ampIntegralMatrix* integral3 = 0;
			printInfo << "reading integral from 'testIntegral.root'" << endl;
			inFile->GetObject("integral", integral3);
			if (not integral3)
				printErr << "cannot find integral 'integral'" << endl;
			else {
				if (*integral3 == integral)
					printSucc << "ROOT file integral is identical" << endl;
				else
					printErr << "ROOT file intergral differs" << endl;
				integral3->writeAscii("testIntegral3.int");
			}
			inFile->Close();
		}
	}


	if (0) {
		// test arthmetic operations
		ampIntegralMatrix  integral, integrals[6];
		const unsigned int nmbInt = sizeof(integrals) / sizeof(integrals[0]);
		integral.readAscii("testIntegral.int");
		for (unsigned int i = 0; i < nmbInt; ++i)
			integrals[i].readAscii("testIntegral.int");
		integrals[0] = integral + integral;
		integrals[1] = integrals[0] - integral;
		integrals[2] = 2 * integral - integral;
		integrals[3] = integral * 0.1 + integral * 0.9;
		integrals[4] = integral / 2 + integral / 2;
		integrals[5] = (integral + 2 * integral) / 3;
		integrals[5].writeAscii("testIntegral1.int");
		for (unsigned int i = 1; i < nmbInt; ++i)
			if (integrals[i] == integral)
				printSucc << "integrals[" << i << "] is identical" << endl;
			else
				printErr << "integrals[" << i << "] differs" << endl;
	}


	if (0) {
		// test renormalization
		ampIntegralMatrix integral;
		integral.readAscii("testIntegral.int");
		const unsigned int nmbEvents = integral.nmbEvents();
		integral.renormalize(10000);
		integral.renormalize(nmbEvents);  // scale back to original
		integral.writeAscii("testIntegral2.int");
		cout << integral;
		//integral.print(cout, true);
	}


}
