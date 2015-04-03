///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2015 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      helper functions related to the fitResult class
//
//-------------------------------------------------------------------------


#include "partialWaveFitHelper.h"

#include "fitResult.h"

void
rpwa::partialWaveFitHelper::extractWaveList(const rpwa::fitResult& result,
                                            std::ostream& waveList)
{
	for(unsigned int i = 0; i < result.nmbWaves(); ++i) {
		const std::string waveName(result.waveName(i));
		if(waveName != "flat") {
			waveList << waveName;

			// test whether this wave is thresholded
			bool thresholded = true;
			for (unsigned int prodAmpIndex = 0; prodAmpIndex < result.nmbProdAmps(); ++prodAmpIndex)
				if (result.waveNameForProdAmp(prodAmpIndex) == waveName)
					thresholded &= (result.prodAmp(prodAmpIndex) == 0.);
			// set the threshold to something higher than the current mass bin
			if (thresholded)
				waveList << " " << 1.1 * result.massBinCenter();

			waveList << std::endl;
		}
	}
}
