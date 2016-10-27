///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2015-2016 Sebastian Uhl (TUM)
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
//      implementation of final-state mass-dependence of resonance fit
//
//-------------------------------------------------------------------------


#include "fsmd.h"

#include <set>

#include <TFormula.h>

#include <reportingUtils.hpp>

#include "cache.h"
#include "parameters.h"
#include "resonanceFitHelper.h"


rpwa::resonanceFit::fsmd::fsmd(const size_t id,
                               const std::vector<size_t>& nrMassBins,
                               const boost::multi_array<double, 2>& massBinCenters,
                               const std::shared_ptr<TFormula>& function,
                               const boost::multi_array<rpwa::resonanceFit::parameter, 1>& parameters)
	: _id(id),
	  _sameFunctionForAllBins(true)
{
	// a single final-state mass-dependence for all bins is used
	// get dimensions from one array and make sure that all other arrays
	// have the same dimensions
	_nrBins = nrMassBins.size();
	if(_nrBins == 0) {
		printErr << "number of bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	const size_t maxMassBins = *(std::max_element(nrMassBins.begin(), nrMassBins.end()));
	if(maxMassBins == 0) {
		printErr << "maximal number of mass bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	checkSize(nrMassBins,
	          _nrBins, "number of bins is not correct for number of mass bins.");
	checkSize(massBinCenters,
	          _nrBins, "number of bins is not correct for centers of mass bins.",
	          maxMassBins, "maximal number of mass bins is not correct for centers of mass bins.");

	_functions.assign(_nrBins, function);

	_maxParameters = function->GetNpar();
	_nrParameters.assign(_nrBins, _maxParameters);
	_parametersIndex.assign(_nrBins, 0);

	_parameters.resize(boost::extents[_nrBins][_maxParameters]);
	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		adjustSizeAndSet(_parameters, idxBin, parameters);
	}

	for(size_t idxParameter = 0; idxParameter < _maxParameters; ++idxParameter) {
		if(_parameters[0][idxParameter].name() != _functions[0]->GetParName(idxParameter)) {
			printErr << "inconsistent naming of parameters for final-state mass-dependence, "
			         << "expected '" << _functions[0]->GetParName(idxParameter) << "', "
			         << "found '" << _parameters[0][idxParameter].name() << "' for parameter at index " << idxParameter << ". "
			         << "Aborting..." << std::endl;
			throw;
		}
	}

	// if all bins have the same mass binning, the value in each mass bin
	// is equal for all bins
	_binsEqualValues.assign(_nrBins, std::vector<size_t>());
	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		if(_binsEqualValues[idxBin].size() != 0) {
			continue;
		}

		std::vector<size_t> bins(1, idxBin);
		for(size_t jdxBin = idxBin+1; jdxBin < _nrBins; ++jdxBin) {
			if(nrMassBins[idxBin] != nrMassBins[jdxBin]) {
				continue;
			}

			bool sameBinning = true;
			for(size_t idxMass = 0; idxMass < nrMassBins[idxBin]; ++idxMass) {
				if(massBinCenters[idxBin][idxMass] != massBinCenters[jdxBin][idxMass]) {
					sameBinning = false;
					break;
				}
			}
			if(sameBinning) {
				bins.push_back(jdxBin);
			}
		}

		for(size_t jdxBin = 0; jdxBin < bins.size(); ++jdxBin) {
			if(_binsEqualValues[bins[jdxBin]].size() != 0) {
				printErr << "inconsistency when setting up bins used for caching." << std::endl;
				throw;
			}

			_binsEqualValues[bins[jdxBin]] = bins;
		}
	}
}


rpwa::resonanceFit::fsmd::fsmd(const size_t id,
                               const std::vector<size_t>& nrMassBins,
                               const boost::multi_array<double, 2>& massBinCenters,
                               const std::vector<std::shared_ptr<TFormula> >& functions,
                               const boost::multi_array<rpwa::resonanceFit::parameter, 2>& parameters)
	: _id(id),
	  _sameFunctionForAllBins(false),
	  _functions(functions),
	  _parameters(parameters)
{
	// a different final-state mass-dependence for each bin is used
	// get dimensions from one array and make sure that all other arrays
	// have the same dimensions
	_nrBins = nrMassBins.size();
	if(_nrBins == 0) {
		printErr << "number of bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	const size_t maxMassBins = *(std::max_element(nrMassBins.begin(), nrMassBins.end()));
	if(maxMassBins == 0) {
		printErr << "maximal number of mass bins is zero, cannot perform the fit. Aborting..." << std::endl;
		throw;
	}

	checkSize(nrMassBins,
	          _nrBins, "number of bins is not correct for number of mass bins.");
	checkSize(massBinCenters,
	          _nrBins, "number of bins is not correct for centers of mass bins.",
	          maxMassBins, "maximal number of mass bins is not correct for centers of mass bins.");

	checkSize(_functions,
	          _nrBins, "number of bins is not correct for functions.");

	_nrParameters.resize(_nrBins);
	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		_nrParameters[idxBin] = _functions[idxBin]->GetNpar();
	}
	_maxParameters = *(std::max_element(_nrParameters.begin(), _nrParameters.end()));

	checkSize(_nrParameters,
	          _nrBins, "number of bins is not correct for number of parameters.");

	_parametersIndex.assign(_nrBins, 0);
	for(size_t idxBin = 1; idxBin < _nrBins; ++idxBin) {
		_parametersIndex[idxBin] = _parametersIndex[idxBin-1] + _nrParameters[idxBin-1];
	}

	checkSize(_parametersIndex,
	          _nrBins, "number of bins is not correct for parameter indices.");

	checkSize(_parameters,
	          _nrBins, "number of bins is not correct for parameters.",
	          _maxParameters, "maximal number of parameters is not correct for parameters.");

	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		for(size_t idxParameter = 0; idxParameter < _nrParameters[idxBin]; ++idxParameter) {
			if(_parameters[idxBin][idxParameter].name() != _functions[idxBin]->GetParName(idxParameter)) {
				printErr << "inconsistent naming of parameters for final-state mass-dependence, "
				         << "expected '" << _functions[idxBin]->GetParName(idxParameter) << "', "
				         << "found '" << _parameters[idxBin][idxParameter].name() << "' "
				         << "for parameter at index " << idxParameter << " in bin " << idxBin << ". "
				         << "Aborting..." << std::endl;
				throw;
			}
		}
	}


	// each value is only valid for the bin it has been calculated for
	_binsEqualValues.assign(_nrBins, std::vector<size_t>(1, 0));
	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		_binsEqualValues[idxBin][0] = idxBin;
	}
}


rpwa::resonanceFit::fsmd::~fsmd()
{
}


size_t
rpwa::resonanceFit::fsmd::importParameters(const double* par,
                                           rpwa::resonanceFit::parameters& fitParameters,
                                           rpwa::resonanceFit::cache& cache)
{
	size_t sumNrParameters = 0;
	const size_t maxNrBins = _sameFunctionForAllBins ? 1 : _nrBins;
	for(size_t idxBin = 0; idxBin < maxNrBins; ++idxBin) {
		bool invalidateCache = false;
		for(size_t idxParameter = 0; idxParameter < _nrParameters[idxBin]; ++idxParameter) {
			const size_t parIndex = _parametersIndex[idxBin] + idxParameter;
			if(fitParameters.getParameter(_id, parIndex) != par[parIndex]) {
				fitParameters.setParameter(_id, parIndex, par[parIndex]);
				invalidateCache = true;
			}
		}
		sumNrParameters += _nrParameters[idxBin];

		if(invalidateCache) {
			if(_sameFunctionForAllBins or _binsEqualValues[idxBin].size() == _binsEqualValues.size()) {
				// the value is the same for all bins
				cache.setComponent(_id, std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
				cache.setProdAmp(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
			} else {
				for(std::vector<size_t>::const_iterator it = _binsEqualValues[idxBin].begin(); it != _binsEqualValues[idxBin].end(); ++it) {
					cache.setComponent(_id, *it, std::numeric_limits<size_t>::max(), 0.);
					cache.setProdAmp(std::numeric_limits<size_t>::max(), *it, std::numeric_limits<size_t>::max(), 0.);
				}
			}
		}
	}

	return sumNrParameters;
}


std::complex<double>
rpwa::resonanceFit::fsmd::val(const rpwa::resonanceFit::parameters& fitParameters,
                              rpwa::resonanceFit::cache& cache,
                              const size_t idxBin,
                              const double mass,
                              const size_t idxMass) const
{
	if(not _functions[idxBin]) {
		return 1.;
	}

	if(idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> fsmd = cache.getComponent(_id, idxBin, idxMass);
		if(fsmd != 0.) {
			return fsmd;
		}
	}

	const std::complex<double> fsmd = _functions[idxBin]->EvalPar(&mass, fitParameters.getParameters(_id)+_parametersIndex[idxBin]);

	if(idxMass != std::numeric_limits<size_t>::max()) {
		if(_binsEqualValues[idxBin].size() == _binsEqualValues.size()) {
			// the value is the same for all bins
			cache.setComponent(_id, std::numeric_limits<size_t>::max(), idxMass, fsmd);
		} else {
			for(std::vector<size_t>::const_iterator it = _binsEqualValues[idxBin].begin(); it != _binsEqualValues[idxBin].end(); ++it) {
				cache.setComponent(_id, *it, idxMass, fsmd);
			}
		}
	}

	return fsmd;
}


std::ostream&
rpwa::resonanceFit::fsmd::print(std::ostream& out) const
{
	out << "final-state mass-dependence (id " << _id << ")" << std::endl
	    << "    use equal values for each mass bin in group of bins: ";
	std::set<size_t> printed;
	for(size_t idxBin = 0; idxBin < _binsEqualValues.size(); ++idxBin) {
		if(printed.count(idxBin) > 0) {
			continue;
		}
		out << ((printed.size() > 0) ? ", " : "") << "{";
		for(size_t idx = 0; idx < _binsEqualValues[idxBin].size(); ++idx) {
			out << ((idx > 0) ? ", " : "") << _binsEqualValues[idxBin][idx];
			printed.insert(_binsEqualValues[idxBin][idx]);
		}
		out << "}";
	}
	out << std::endl;

	const size_t maxNrBins = _sameFunctionForAllBins ? 1 : _nrBins;
	for(size_t i = 0; i < maxNrBins; ++i) {
		printBin(i, out);
	}
	if(_sameFunctionForAllBins) {
		out << "this final-state mass-dependence is used for all bins." << std::endl;
	}

	return out;
}


std::ostream&
rpwa::resonanceFit::fsmd::printBin(const size_t idxBin,
                                   std::ostream& out) const
{
	out << "final-state mass-dependence for bin " << idxBin << std::endl;
	out << "formula: " << (_functions[idxBin] ? _functions[idxBin]->GetTitle() : "not set") << std::endl;

	for(size_t i = 0; i < _nrParameters[idxBin]; ++i) {
		out << "    [" << i << "] ";
		if(_parameters[idxBin][i].limitedLower() and _parameters[idxBin][i].limitedUpper()) {
			out << "limits: " << _parameters[idxBin][i].limitLower() << "-" << _parameters[idxBin][i].limitUpper() << " GeV/c^2";
		} else if(_parameters[idxBin][i].limitedLower()) {
			out << "lower limit: " << _parameters[idxBin][i].limitLower() << " GeV/c^2";
		} else if(_parameters[idxBin][i].limitedUpper()) {
			out << "upper limit: " << _parameters[idxBin][i].limitUpper() << " GeV/c^2";
		} else {
			out << "unlimited";
		}
		out << (_parameters[idxBin][i].fixed() ? " (FIXED) " : "") << std::endl;
	}

	return out;
}
