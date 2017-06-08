///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2016 Sebastian Uhl (TUM)
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
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      implementation of the storage for the general input required by
//      the resonance fit
//
//-------------------------------------------------------------------------


#include "input.h"


rpwa::resonanceFit::input::bin::wave::wave(const std::string& waveName,
                                           const std::pair<double, double>& massLimits)
	: _waveName(waveName),
	  _massLimits(massLimits)
{
}


std::ostream&
rpwa::resonanceFit::input::bin::wave::print(std::ostream& out, const bool newLine) const
{
	out << "wave '" << _waveName << "'";

	if(_massLimits.first >= 0.0 and _massLimits.second >= 0.0) {
		out << " in mass range " << _massLimits.first << " to " << _massLimits.second << " GeV/c^2";
	} else if(_massLimits.first >= 0.0) {
		out << " for masses above " << _massLimits.first << " GeV/c^2";
	} else if(_massLimits.second >= 0.0) {
		out << " for masses below " << _massLimits.second << " GeV/c^2";
	}
	out << ".";

	if(newLine) {
		out << std::endl;
	}

	return out;
}


rpwa::resonanceFit::input::bin::bin(const std::string& fileName,
                                    const std::vector<rpwa::resonanceFit::input::bin::wave>& waves,
                                    const double tPrimeMean,
                                    const double rescaleErrors,
                                    const std::vector<std::string>& sysFileNames)
	: _fileName(fileName),
	  _waves(waves),
	  _tPrimeMean(tPrimeMean),
	  _rescaleErrors(rescaleErrors),
	  _sysFileNames(sysFileNames)
{
}


std::ostream&
rpwa::resonanceFit::input::bin::print(std::ostream& out, const bool newLine) const
{
	out << "fit result at '" << _fileName << "' with mean tPrime = " << _tPrimeMean;

	if(_rescaleErrors != 1.0) {
		out << ", errors are rescaled by " << _rescaleErrors;
	}

	if(_sysFileNames.size() > 0) {
		out << ", systematics bands are obtained from [ ";
		for(size_t idx = 0; idx < _sysFileNames.size(); ++idx) {
			out << ((idx != 0) ? ", " : "") << "'" << _sysFileNames[idx] << "'";
		}
		out << " ]";
	}

	out << ". ";

	const size_t nrWaves = this->nrWaves();
	out << nrWaves << " wave" << ((nrWaves != 1) ? "s" : "") << " used in this bin:";
	for(size_t idxWave = 0; idxWave < nrWaves; ++idxWave) {
		out << std::endl;
		out << getWave(idxWave);
	}

	if(newLine) {
		out << std::endl;
	}

	return out;
}


rpwa::resonanceFit::input::input(const std::vector<rpwa::resonanceFit::input::bin>& bins)
	: _bins(bins)
{
}


std::ostream&
rpwa::resonanceFit::input::print(std::ostream& out, const bool newLine) const
{
	const size_t nrBins = this->nrBins();

	out << "fit input for " << nrBins << " bin" << ((nrBins != 1) ? "s" : "") << ":";

	for(size_t idxBin = 0; idxBin < nrBins; ++idxBin) {
		out << std::endl;
		out << getBin(idxBin);
	}

	if(newLine) {
		out << std::endl;
	}

	return out;
}
