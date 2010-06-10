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
// File and Version Information:
// $Id$
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

#include "particleKey.h"
#include "waveKey.h"

#include "genKeyHelper.h"

using namespace std;
using namespace rpwa;

void genKey_Kpipi(const bool testKey = true, const string& dataFileName =
		"../keyfiles/keyKpipi/testEventsKpipi.evt", // file with test data in .evt format
		const string& pdgTableFileName = "./pdgTable.txt") {

	const string thisFilePath = __FILE__;
	// define final state particles
	particleKey pi_minus("pi-");
	particleKey pi_plus("pi+");
	particleKey K_minus("K-");

	// define isobars: (name, daughter1, daughter2, L, S, mass dependence)

	// K pi decay modes
	// Daum et al. used a very controversial resonance, the kappa or K*(800)
	//particleKey Kstar800("Kstar(800)",  &K_minus, &pi_plus, 0, 1); // not in the PDG table yet
	particleKey Kstar892("Kstar(892)0", &K_minus, &pi_plus, 0, 1); // 1-
	// probably I will have to treat both Kstar(1430) as a broad
	// Kpi-P wave since they are overlapping, but it is also very
	// likely that I'm talking bullshit
	particleKey Kstar01430("Kstar0(1430)" , &K_minus, &pi_plus, 0, 0); // 0++
	particleKey Kstar21430("Kstar2(1430)0", &K_minus, &pi_plus, 2, 0); // 2++

	// pi pi decay modes
	// f0(980) is also called as the epsilon in Daum et al. paper, I guess
	particleKey f0980 ("f0(980)" , &pi_plus, &pi_minus, 0, 0); // 0++
	particleKey rho770("rho(770)", &pi_plus, &pi_minus, 1, 0); // 1--
	particleKey f21270("f2(1270)", &pi_plus, &pi_minus, 2, 0); // 2++

	// PWA wave set according to Daum et al.
	// ..JPLMpartiyexchange(isobar1 isobar2)
	{ // 0-P0+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 1, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 0,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+S0+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}

	{ // 0-S0+(eK)
		particleKey X("X", &f0980, &K_minus, 0, 0);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 0,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+P0+(eK)
		particleKey X("X", &f0980, &K_minus, 1, 0);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+P0+(kpi)
//		particleKey X("X", &Kstar0800, &pi_minus, 1, 0);
		//      wave(&X, J, P, M, refl);
//		waveKey wave(&X, 1,+1, 0,+1);
//		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+S0+(pK)
		particleKey X("X", &rho770, &K_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+S1+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+S1+(pK)
		particleKey X("X", &rho770, &K_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 2+D1+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 2, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 2-P0+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 1, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}

	{ // 0?P0+(pK)
		particleKey X("X", &rho770, &K_minus, 1, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 0,-1, 0,+1); // I'm not sure about the parity
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+D0+(K*pi)
		particleKey X("X", &Kstar892, &pi_minus, 2, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+P1+(kpi)
//		particleKey X("X", &Kstar0800, &pi_minus, 1, 0);
		//      wave(&X, J, P, M, refl);
//		waveKey wave(&X, 1,+1, 1,+1);
//		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}

	{ // 2+D1+(pK)
		particleKey X("X", &rho770, &K_minus, 2, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}

	{ // 1+D0+(pK)
		particleKey X("X", &rho770, &K_minus, 2, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 2-P0+(pK)
		particleKey X("X", &rho770, &K_minus, 1, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}

	{ // 2?S0+(K**pi)
		particleKey X("X", &Kstar21430, &pi_minus, 0, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1); // I'm not sure about the parity here
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}

	{ // 2-S0+(fK)
		particleKey X("X", &f21270, &K_minus, 0, 2);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 2,-1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	// Flip waves, whatever this means
	/*
	{ // 1+S0+(pK)
		particleKey X("X", &rho770, &K_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}
	{ // 1+S1+(pK)
		particleKey X("X", &rho770, &K_minus, 0, 1);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 1,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}

	{ // 1+P0+(eK)
		particleKey X("X", &f0980, &K_minus, 1, 0);
		//      wave(&X, J, P, M, refl);
		waveKey wave(&X, 1,+1, 0,+1);
		generateKeyFile(wave, thisFilePath, testKey, dataFileName, pdgTableFileName);
	}*/
}
