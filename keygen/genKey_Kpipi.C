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

void genKey_Kpipi(const bool testKey = true, const string& dataFileName =
		"../keyfiles/keyKpipi/testEventsKpipi.evt", // file with test data in .evt format
		const string& pdgTableFileName = "./pdgTable.txt") {

	const string thisFilePath = __FILE__;
	const string movetoFilePath = "${ROOTPWA}/keyfiles/keyKpipi/SETX";
	// define final state particles
	particleKey pi_minus("pi-");
	particleKey pi_plus("pi+");
	particleKey K_minus("K-");

	// define isobars: (name, daughter1, daughter2, L, S, mass dependence)

	// K pi decay modes
	// Daum et al. used a very controversial resonance, the kappa or K*(800)
	particleKey kappa     ("kappa"        , &K_minus, &pi_plus, 0, 0); //0+
	particleKey Kstar892  ("Kstar(892)0"  , &K_minus, &pi_plus, 1, 0); // 1-
	particleKey Kstar01430("Kstar0(1430)" , &K_minus, &pi_plus, 0, 0); // 0++
	particleKey Kstar21430("Kstar2(1430)0", &K_minus, &pi_plus, 2, 0); // 2++
	particleKey Kstar16800("Kstar(1680)"  , &K_minus, &pi_plus, 1, 0); // 1-
	particleKey Kstar17800("Kstar3(1780)" , &K_minus, &pi_plus, 3, 0); // 3-


	// pi pi decay modes
	particleKey f0980   ("f0(980)"   , &pi_plus, &pi_minus, 0, 0); // 0++
	particleKey sigma   ("sigma"     , &pi_plus, &pi_minus, 0, 0, "amp_kach"); // 0++
	particleKey f01370  ("f0(1370)"  , &pi_plus, &pi_minus, 0, 0); // 0++
	particleKey f01500  ("f0(1500)"  , &pi_plus, &pi_minus, 0, 0); // 0++
	particleKey rho770  ("rho(770)"  , &pi_plus, &pi_minus, 1, 0); // 1--
	particleKey rho1450 ("rho(1450)" , &pi_plus, &pi_minus, 1, 0); // 1--
	particleKey f21270  ("f2(1270)"  , &pi_plus, &pi_minus, 2, 0); // 2++
	particleKey rho31690("rho3(1690)", &pi_plus, &pi_minus, 3, 0); // 3--

	if (1){
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
								/*
								stringstream command;
								command << "mv " << wave.waveName(true) << " " << movetoFilePath << "/";
								cout << " executing " << command.str() << endl;
								system(command.str().c_str());
								command.str("");
								// remove the .C file
								command << "rm -f " << wave.waveName(true) << ".C" << endl;
								cout << " executing " << command.str() << endl;
								system(command.str().c_str());*/
								wavecounter++;
							}
						}
					}
				}
			}
		}
	}

	cout << " created " << wavecounter << " keys " << endl;

	} // own wave set


	if (0){
	// PWA wave set according to Daum et al.
	// ..JPLMpartiyexchange(isobar1 isobar2)
	{ // 0-P0+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 1, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 0,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+S0+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}

	{ // 0-S0+(eK)
		particleKey X("X", &sigma, &K_minus, 0, 0);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 0,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+P0+(eK)
		particleKey X("X", &sigma, &K_minus, 1, 0);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+P0+(kpi)
		particleKey X("X", &kappa, &pi_minus, 1, 0);
		//	    wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+S0+(pK)
		particleKey X("X", &rho770, &K_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+S1+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+S1+(pK)
		particleKey X("X", &rho770, &K_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 2+D1+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 2, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 2-P0+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 1, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}

	{ // 0-P0+(pK)
		particleKey X("X", &rho770, &K_minus, 1, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 0,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+D0+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 2, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+P1+(kpi)
		particleKey X("X", &kappa, &pi_minus, 1, 0);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}

	{ // 2+D1+(pK)
		particleKey X("X", &rho770, &K_minus, 2, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}

	{ // 1+D0+(pK)
		particleKey X("X", &rho770, &K_minus, 2, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 2-P0+(pK)
		particleKey X("X", &rho770, &K_minus, 1, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}

	{ // 2?S0+(K**pi)
		particleKey X("X", &Kstar21430, &pi_minus, 0, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1); // I'm not sure about the parity here
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}

	{ // 2-S0+(fK)
		particleKey X("X", &f21270, &K_minus, 0, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	// Flip waves, whatever this means
	/*
	{ // 1+S0+(pK)
		particleKey X("X", &rho770, &K_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}
	{ // 1+S1+(pK)
		particleKey X("X", &rho770, &K_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}

	{ // 1+P0+(eK)
		particleKey X("X", &f0980, &K_minus, 1, 0);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, false, pdgTableFileName);
	}*/
	}
}
