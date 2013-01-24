///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
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
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
//
// Description:
//      example macro to create 5pi key file for gamp     
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//		Prometeusz Jasinski	 KPH
//
//
//-----------------------------------------------------------


#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <stdio.h>
#include <sstream>


#include "particleKey.h"
#include "waveKey.h"

#include "genKeyHelper.h"

using namespace std;
using namespace rpwa;

void genKey_5pi_crosscheck(const bool testKey = true, const string& dataFileName =
		"../keyfiles/key5pi/testdata_crosscheck.evt", // file with test data in .evt format
		const string& pdgTableFileName = "./pdgTable.txt") {

	const string thisFilePath = __FILE__;
	const string movetoFilePath = "${ROOTPWA}/keyfiles/keyKpipi/SETX";
	// define final state particles
	particleKey pi_minus_1("pi-");
	particleKey pi_plus_1("pi+");
	particleKey pi_minus_2("pi-");
	particleKey pi_plus_2("pi+");
	particleKey pi_minus_3("pi-");

	// define isobars: (name, daughter1, daughter2, L, S, mass dependence)

	// pi pi decay modes
	particleKey f0980_1   ("f0(980)"   , &pi_plus_1, &pi_minus_1, 0, 0); // 0++
	particleKey sigma_1   ("sigma"     , &pi_plus_1, &pi_minus_1, 0, 0, "amp_kach"); // 0++
	particleKey f01370_1  ("f0(1370)"  , &pi_plus_1, &pi_minus_1, 0, 0); // 0++
	particleKey f01500_1  ("f0(1500)"  , &pi_plus_1, &pi_minus_1, 0, 0); // 0++
	particleKey rho770_1  ("rho(770)"  , &pi_plus_1, &pi_minus_1, 1, 0); // 1--
	particleKey rho1450_1 ("rho(1450)" , &pi_plus_1, &pi_minus_1, 1, 0); // 1--
	particleKey f21270_1  ("f2(1270)"  , &pi_plus_1, &pi_minus_1, 2, 0); // 2++
	particleKey rho31690_1("rho3(1690)", &pi_plus_1, &pi_minus_1, 3, 0); // 3--

	particleKey f0980_2   ("f0(980)"   , &pi_plus_2, &pi_minus_2, 0, 0); // 0++
	particleKey sigma_2   ("sigma"     , &pi_plus_2, &pi_minus_2, 0, 0, "amp_kach"); // 0++
	particleKey f01370_2  ("f0(1370)"  , &pi_plus_2, &pi_minus_2, 0, 0); // 0++
	particleKey f01500_2  ("f0(1500)"  , &pi_plus_2, &pi_minus_2, 0, 0); // 0++
	particleKey rho770_2  ("rho(770)"  , &pi_plus_2, &pi_minus_2, 1, 0); // 1--
	particleKey rho1450_2 ("rho(1450)" , &pi_plus_2, &pi_minus_2, 1, 0); // 1--
	particleKey f21270_2  ("f2(1270)"  , &pi_plus_2, &pi_minus_2, 2, 0); // 2++
	particleKey rho31690_2("rho3(1690)", &pi_plus_2, &pi_minus_2, 3, 0); // 3--

	// pi pi pi decay modes
		// a1(1260) -> rho pi S-Wave
   	particleKey a11260_rhopiplus   ("a1(1269)", &rho770_1, &pi_plus_2 , 0, 1); // 1-(1++)
	particleKey a11260_rhopiminus  ("a1(1269)", &rho770_1, &pi_minus_2, 0, 1); // 1-(1++)
		// the same but the other way round
	particleKey a11260_piplusrho   ("a1(1269)", &pi_plus_2, &rho770_1 , 0, 1); // 1-(1++)
	particleKey a11260_piminusrho  ("a1(1269)", &pi_minus_2, &rho770_1, 0, 1); // 1-(1++)
		// a1(1260) -> sigma pi S-Wave
	particleKey a11260_sigmapiplus ("a1(1269)", &sigma_1 , &pi_plus_2 , 0, 1); // 1-(1++)
	particleKey a11260_sigmapiminus("a1(1269)", &sigma_1 , &pi_minus_2, 0, 1); // 1-(1++)
		// a2(1320) -> rho pi
	particleKey a21320_rhopiplus   ("a2(1320)", &rho770_1, &pi_plus_2 , 1, 1); // 1-(2++)
	particleKey a21320_rhopiminus  ("a2(1320)", &rho770_1, &pi_minus_2, 1, 1); // 1-(2++)
		// pi(1300) -> rho pi
	particleKey pi1300_rhopiplus   ("pi(1300)", &rho770_1, &pi_plus_2 , 1, 1); // 1-(0-+)
	particleKey pi1300_rhopiminus  ("pi(1300)", &rho770_1, &pi_minus_2, 1, 1); // 1-(0-+)
		// pi(1300) -> sigma pi
	particleKey pi1300_sigmapiplus   ("pi(1300)", &sigma_1, &pi_plus_2 , 0, 1); // 1-(0-+)
	particleKey pi1300_sigmapiminus  ("pi(1300)", &sigma_1, &pi_minus_2, 0, 1); // 1-(0-+)
		// pi(1800) -> rho pi
	particleKey pi1800_rhopiplus   ("pi(1300)", &rho770_1, &pi_plus_2 , 0, 1); // 1-(0-+)
	particleKey pi1800_rhopiminus  ("pi(1300)", &rho770_1, &pi_minus_2, 0, 1); // 1-(0-+)
		// pi(1800) -> sigma pi
	particleKey pi1800_sigmapiplus   ("pi(1800)", &sigma_1, &pi_plus_2 , 0, 1); // 1-(0-+)
	particleKey pi1800_sigmapiminus  ("pi(1800)", &sigma_1, &pi_minus_2, 0, 1); // 1-(0-+)
		// pi(1600) and pi(1670) were not picked by the evolution algorithm of s.neubert

	// pi pi pi pi decay modes
		// f0(1370) -> rho rho
	particleKey f01370_rhorho    ("f0(1370)", &rho770_1  , &rho770_2  , 0, 1); // 0+(0++)
		// f0(1370) -> sigma sigma
	particleKey f01370_sigmasigma("f0(1370)", &sigma_1, &sigma_2, 0, 1); // 0+(0++)
		// f01500 f01700 same and the decay into pi1300 pi was anyhow not observed

		// f2(1270) -> pi +/- a1(1260)->rho(770) pi -/+
	particleKey f21270_a11260_rhopi_piminus("f2(1270)", &a11260_rhopiplus , &pi_minus_2, 1, 1); // 0+(2++)
	particleKey f21270_a11260_rhopi_piplus ("f2(1270)", &a11260_rhopiminus, &pi_plus_2 , 1, 1); // 0+(2++)
		// the same but the other way round
	particleKey f21270_piminusa11260_pirho("f2(1270)", &pi_minus_2, &a11260_piplusrho , 1, 1); // 0+(2++)
	particleKey f21270_piplusa11260_pirho ("f2(1270)", &pi_plus_2 , &a11260_piminusrho, 1, 1); // 0+(2++)

if (1){ // reduced wave set
	cout << " creating the reduced 5pi wave set " << endl;
	{ 
		particleKey X("X", &f21270_a11260_rhopi_piminus, &pi_minus_3, 0, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ 
		particleKey X("X", &f21270_a11260_rhopi_piplus, &pi_minus_3, 0, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	// should give the same results as the previous
	{ 
		particleKey X("X", &pi_minus_3, &f21270_piminusa11260_pirho, 0, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ 
		particleKey X("X", &pi_minus_3, &f21270_piplusa11260_pirho, 0, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}

	{ 
		particleKey X("X", &rho770_2, &a11260_rhopiminus, 2, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 0,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	cout << " done " << endl;
	} // reduced wave set

	/*if (0){
	// let's create all possible combinations of decays according
	// to the rules of L,S coupling and parity

	//		pair<parity,   isobar>
	vector< pair<int,particleKey*> > isobar1[2]; // {K*, pi*}
	pair<int,particleKey*> isobar2[2]; 			 // {pi, K  }

	// fill all isobars
	isobar2[0] = make_pair(-1,&pi_minus);
	isobar2[1] = make_pair(-1,&K_minus);

	isobar1[0].push_back(make_pair(+1,&Kstar01430));
	isobar1[0].push_back(make_pair(-1,&Kstar892  ));
	isobar1[0].push_back(make_pair(+1,&Kstar21430));
	isobar1[0].push_back(make_pair(+1,&kappa     ));
	isobar1[0].push_back(make_pair(-1,&Kstar16800));
	isobar1[0].push_back(make_pair(-1,&Kstar17800));

	isobar1[1].push_back(make_pair(+1,&f0980   ));
	isobar1[1].push_back(make_pair(-1,&rho770  ));
	isobar1[1].push_back(make_pair(+1,&f21270  ));
	isobar1[1].push_back(make_pair(+1,&sigma   ));
	isobar1[1].push_back(make_pair(+1,&f01370  ));
	isobar1[1].push_back(make_pair(+1,&f01500  ));
	isobar1[1].push_back(make_pair(-1,&rho1450 ));
	isobar1[1].push_back(make_pair(-1,&rho31690));


	unsigned int Mmax    = 1;
	unsigned int lorbmax = 4;
	unsigned int Jmax    = 5;

	int wavecounter(0);

	for (unsigned int comb = 0; comb < 2; comb++){ // go through the to Combinations K* pi and pi* K
		for (unsigned int i = 0; i < isobar1[comb].size(); i++){ // loop over all K* isobars
			cout << endl;
			int parity_isobar1 = isobar1[comb][i].first;
			int parity_isobar2 = isobar2[comb].first;
			int J_isobar1 = isobar1[comb][i].second->L();
			int J_isobar2 = isobar2[comb].second->L();
			for (unsigned int spin = abs(J_isobar1-J_isobar2); spin <= (unsigned) abs(J_isobar1+J_isobar2); spin++){ // couple the spins
				for (unsigned int lorb = 0; lorb <= lorbmax; lorb++){ // allow only orbital angular momenta up to lorbmax
					// construct the parity
					int parity = parity_isobar1 * parity_isobar2 * pow(-1,lorb);
					// create the particle key
					particleKey X("X", isobar1[comb][i].second, isobar2[comb].second, lorb, spin);
					for (unsigned int J = abs(lorb-spin); (J <= (unsigned) abs(lorb+spin) && J <= Jmax); J++){ // couple now spin and orbital angular momentum
						for (unsigned int M = 0; (M <= J && M <= Mmax); M++){ // go through the J projections up to either the allowed value or Mmax
							for (int reflectivity = -1; reflectivity < 2; reflectivity+=2){
								// skip combinations that are forbidden in the reflectivity basis
								if (J==0 && parity == -1 && M == 0 && reflectivity == -1) continue;
								if (J==1 && parity == +1 && M == 0 && reflectivity == -1) continue;
								if (J==1 && parity == -1 && M == 0 && reflectivity == +1) continue;
								if (J==2 && parity == +1 && M == 0 && reflectivity == +1) continue;
								if (J==2 && parity == -1 && M == 0 && reflectivity == -1) continue;
								if (J==3 && parity == +1 && M == 0 && reflectivity == -1) continue;
								if (J==3 && parity == -1 && M == 0 && reflectivity == +1) continue;
								if (J==4 && parity == +1 && M == 0 && reflectivity == +1) continue;
								if (J==4 && parity == -1 && M == 0 && reflectivity == -1) continue;
								if (J==5 && parity == +1 && M == 0 && reflectivity == -1) continue;
								if (J==5 && parity == -1 && M == 0 && reflectivity == +1) continue;

								cout << "JPMe LS: " << J << " " << parity << " " << M << " " << reflectivity << " " << lorb << " " << spin << endl;
								waveKey wave(&X, J,parity,M,reflectivity);
								generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
								// move the .key file
								wavecounter++;
							}
						}
					}
				}
			}
		}
	}

	cout << " created " << wavecounter << " keys " << endl;

	}*/ // own wave set
}
