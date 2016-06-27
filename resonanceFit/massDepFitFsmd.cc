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


#include "massDepFitFsmd.h"

#include <map>

#include <boost/assign/std/vector.hpp>

#include <yaml-cpp/yaml.h>

#include <TFormula.h>

#include "conversionUtils.hpp"
#include "massDepFitCache.h"
#include "massDepFitParameters.h"
#include "reportingUtils.hpp"
#include "yamlCppUtils.hpp"


rpwa::massDepFit::fsmd::fsmd(const size_t id)
	: _id(id),
	  _nrBins(0),
	  _maxParameters(0),
	  _sameMassBinning(false)
{
}


rpwa::massDepFit::fsmd::~fsmd()
{
}


bool
rpwa::massDepFit::fsmd::init(const YAML::Node& configFsmd,
                             rpwa::massDepFit::parameters& fitParameters,
                             rpwa::massDepFit::parameters& fitParametersError,
                             const size_t nrBins,
                             const bool sameMassBinning,
                             const bool debug)
{
	if(debug) {
		printDebug << "start initializing final-state mass-dependence." << std::endl;
	}

	_sameMassBinning = sameMassBinning;

	if(not configFsmd.IsMap() and not configFsmd.IsSequence()) {
		printErr << "'finalStateMassDependence' is not a YAML map or sequence." << std::endl;
		return false;
	}

	if(configFsmd.IsSequence()) {
		if(nrBins != configFsmd.size()) {
			printErr << "incorrect number of final-state mass-dependencies provided, found " << configFsmd.size() << ", expected " << nrBins << "." << std::endl;
			return false;
		}

		_nrBins = nrBins;
		_functions.resize(_nrBins);
		_nrParameters.resize(_nrBins, 0);
		_parametersIndex.resize(_nrBins, 0);

		for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
			if(not initBin(configFsmd[idxBin], fitParameters, fitParametersError, idxBin, debug)) {
				printErr << "error while initializing final-state mass-dependence for bin " << idxBin << "." << std::endl;
				return false;
			}
		}
	}

	if(configFsmd.IsMap()) {
		_nrBins = 1;
		_functions.resize(nrBins);
		_nrParameters.resize(nrBins, 0);
		_parametersIndex.resize(nrBins, 0);

		if(not initBin(configFsmd, fitParameters, fitParametersError, 0, debug)) {
			printErr << "error while initializing final-state mass-dependence for bin 0." << std::endl;
			return false;
		}

		for(size_t idxBin = 1; idxBin < nrBins; ++idxBin) {
			_functions[idxBin] = _functions[0];
			_nrParameters[idxBin] = _nrParameters[0];
			_parametersIndex[idxBin] = _parametersIndex[0];
		}
	}

	if(debug) {
		print(printDebug);
		printDebug << "finished initializing final-state mass-dependence." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::fsmd::initBin(const YAML::Node& configFsmd,
                                rpwa::massDepFit::parameters& fitParameters,
                                rpwa::massDepFit::parameters& fitParametersError,
                                const size_t idxBin,
                                const bool debug)
{
	if(debug) {
		printDebug << "start initializing final-state mass-dependence for bin " << idxBin << "." << std::endl;
	}

	if(not configFsmd.IsMap()) {
		printErr << "'finalStateMassDependence' for bin " << idxBin << " is not a YAML map." << std::endl;
		return false;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("formula", YamlCppUtils::TypeString);
	if(not checkIfAllVariablesAreThere(configFsmd, mandatoryArguments)) {
		printErr << "'finalStateMassDependence' does not contain all required variables." << std::endl;
		return false;
	}

	const std::string formula = configFsmd["formula"].as<std::string>();

	_functions[idxBin] = std::make_shared<TFormula>("finalStateMassDependence", formula.c_str());

	_nrParameters[idxBin] = _functions[idxBin]->GetNpar();
	if(idxBin > 0) {
		_parametersIndex[idxBin] = _parametersIndex[idxBin-1] + _nrParameters[idxBin-1];
	}

	if(_nrParameters[idxBin] > _maxParameters) {
		_parametersFixed.resize(boost::extents[_nrBins][_nrParameters[idxBin]]);
		_parametersLimitLower.resize(boost::extents[_nrBins][_nrParameters[idxBin]]);
		_parametersLimitedLower.resize(boost::extents[_nrBins][_nrParameters[idxBin]]);
		_parametersLimitUpper.resize(boost::extents[_nrBins][_nrParameters[idxBin]]);
		_parametersLimitedUpper.resize(boost::extents[_nrBins][_nrParameters[idxBin]]);
		_parametersName.resize(boost::extents[_nrBins][_nrParameters[idxBin]]);
		_parametersStep.resize(boost::extents[_nrBins][_nrParameters[idxBin]]);

		for(size_t i = 0; i < _nrBins; ++i) {
			for(size_t j = _maxParameters; j < _nrParameters[idxBin]; ++j) {
				_parametersStep[i][j] = 0.0001;
			}
		}

		_maxParameters = _nrParameters[idxBin];
	}

	fitParameters.resize(_id+1, 0, _parametersIndex[idxBin]+_nrParameters[idxBin]+1, 0);
	fitParametersError.resize(_id+1, 0, _parametersIndex[idxBin]+_nrParameters[idxBin]+1, 0);

	for(size_t idxParameter = 0; idxParameter < _nrParameters[idxBin]; ++idxParameter) {
		const std::string parName = _functions[idxBin]->GetParName(idxParameter);
		_parametersName[idxBin][idxParameter] = parName;

		if(debug) {
			printDebug << "reading parameter '" << parName << "'." << std::endl;
		}

		const YAML::Node& configParameter = configFsmd[parName];
		if(not configParameter) {
			printErr << "final-state mass-dependence does not define parameter '" << parName << "'." << std::endl;
			return false;
		}

		std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
		boost::assign::insert(mandatoryArguments)
		                     ("val", YamlCppUtils::TypeFloat)
		                     ("fix", YamlCppUtils::TypeBoolean);
		if(not checkIfAllVariablesAreThere(configParameter, mandatoryArguments)) {
			printErr << "'" << parName << "' of final-state mass-dependence does not contain all required variables." << std::endl;
			return false;
		}

		const double parameter = configParameter["val"].as<double>();
		fitParameters.setParameter(_id, _parametersIndex[idxBin]+idxParameter, parameter);

		if(configParameter["error"]) {
			if(checkVariableType(configParameter["error"], YamlCppUtils::TypeFloat)) {
				const double error = configParameter["error"].as<double>();
				fitParametersError.setParameter(_id, _parametersIndex[idxBin]+idxParameter, error);
			} else {
				printErr << "variable 'error' for parameter '" << parName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				return false;
			}
		}

		const bool fixed = configParameter["fix"].as<bool>();
		_parametersFixed[idxBin][idxParameter] = fixed;

		if(configParameter["lower"]) {
			if(checkVariableType(configParameter["lower"], YamlCppUtils::TypeFloat)) {
				_parametersLimitedLower[idxBin][idxParameter] = true;
				_parametersLimitLower[idxBin][idxParameter] = configParameter["lower"].as<double>();
			} else {
				printErr << "variable 'lower' for parameter '" << parName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				return false;
			}
		} else {
			_parametersLimitedLower[idxBin][idxParameter] = false;
		}
		if(configParameter["upper"]) {
			if(checkVariableType(configParameter["upper"], YamlCppUtils::TypeFloat)) {
				_parametersLimitedUpper[idxBin][idxParameter] = true;
				_parametersLimitUpper[idxBin][idxParameter] = configParameter["upper"].as<double>();
			} else {
				printErr << "variable 'upper' for parameter '" << parName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				return false;
			}
		} else {
			_parametersLimitedUpper[idxBin][idxParameter] = false;
		}

		if(configParameter["step"]) {
			if(checkVariableType(configParameter["step"], YamlCppUtils::TypeFloat)) {
				_parametersStep[idxBin][idxParameter] = configParameter["step"].as<double>();
			} else {
				printErr << "variable 'step' for parameter '" << parName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				return false;
			}
		}
	}

	if(debug) {
		printDebug << "finished initializing final-state mass-dependence for bin " << idxBin << "." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::fsmd::write(YAML::Emitter& yamlOutput,
                              const rpwa::massDepFit::parameters& fitParameters,
                              const rpwa::massDepFit::parameters& fitParametersError,
                              const bool debug) const
{
	if(debug) {
		printDebug << "start writing final-state mass-dependence." << std::endl;
		print(printDebug);
	}

	if(_nrBins > 1) {
		yamlOutput << YAML::BeginSeq;
	}

	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
		if(not writeBin(yamlOutput, fitParameters, fitParametersError, idxBin, debug)) {
			printErr << "error while writing final-state mass-dependence for bin " << idxBin << "." << std::endl;
			return false;
		}
	}

	if(_nrBins > 1) {
		yamlOutput << YAML::EndSeq;
	}

	return true;
}


bool
rpwa::massDepFit::fsmd::writeBin(YAML::Emitter& yamlOutput,
                                 const rpwa::massDepFit::parameters& fitParameters,
                                 const rpwa::massDepFit::parameters& fitParametersError,
                                 const size_t idxBin,
                                 const bool debug) const
{
	if(debug) {
		printDebug << "start writing final-state mass-dependence for bin " << idxBin << "." << std::endl;
	}

	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "formula";
	yamlOutput << YAML::Value << _functions[idxBin]->GetTitle();

	for(size_t idxParameter = 0; idxParameter < _nrParameters[idxBin]; ++idxParameter) {
		yamlOutput << YAML::Key << _parametersName[idxBin][idxParameter];
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "val";
		yamlOutput << YAML::Value << fitParameters.getParameter(_id, _parametersIndex[idxBin]+idxParameter);

		yamlOutput << YAML::Key << "error";
		yamlOutput << YAML::Value << fitParametersError.getParameter(_id, _parametersIndex[idxBin]+idxParameter);

		if(_parametersLimitedLower[idxBin][idxParameter]) {
			yamlOutput << YAML::Key << "lower";
			yamlOutput << YAML::Value << _parametersLimitLower[idxBin][idxParameter];
		}

		if(_parametersLimitedUpper[idxBin][idxParameter]) {
			yamlOutput << YAML::Key << "upper";
			yamlOutput << YAML::Value << _parametersLimitUpper[idxBin][idxParameter];
		}

		yamlOutput << YAML::Key << "step";
		yamlOutput << YAML::Value << _parametersStep[idxBin][idxParameter];

		yamlOutput << YAML::Key << "fix";
		yamlOutput << YAML::Value << _parametersFixed[idxBin][idxParameter];

		yamlOutput << YAML::EndMap;
	}

	yamlOutput << YAML::EndMap;

	return true;
}


size_t
rpwa::massDepFit::fsmd::importParameters(const double* par,
                                         rpwa::massDepFit::parameters& fitParameters,
                                         rpwa::massDepFit::cache& cache)
{
	size_t sumNrParameters = 0;
	for(size_t idxBin = 0; idxBin < _nrBins; ++idxBin) {
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
			cache.setComponent(_id, ((_nrBins == 1) ? std::numeric_limits<size_t>::max() : idxBin), std::numeric_limits<size_t>::max(), 0.);
			cache.setProdAmp(std::numeric_limits<size_t>::max(), ((_nrBins == 1) ? std::numeric_limits<size_t>::max() : idxBin), std::numeric_limits<size_t>::max(), 0.);
		}
	}

	return sumNrParameters;
}


std::complex<double>
rpwa::massDepFit::fsmd::val(const rpwa::massDepFit::parameters& fitParameters,
                            rpwa::massDepFit::cache& cache,
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
		cache.setComponent(_id, ((_nrBins == 1 and _sameMassBinning) ? std::numeric_limits<size_t>::max() : idxBin), idxMass, fsmd);
	}

	return fsmd;
}


std::ostream&
rpwa::massDepFit::fsmd::print(std::ostream& out) const
{
	out << "final-state mass-dependence (id " << _id << ")" << std::endl
	    << "    all bins have the same mass binning: " << rpwa::yesNo(_sameMassBinning) << std::endl;

	for(size_t i = 0; i < _nrBins; ++i) {
		printBin(i, out);
	}

	return out;
}


std::ostream&
rpwa::massDepFit::fsmd::printBin(const size_t idxBin,
                                 std::ostream& out) const
{
	out << "final-state mass-dependence for bin " << idxBin << std::endl;
	out << "formula: " << ((_functions[idxBin] != NULL) ? _functions[idxBin]->GetTitle() : "not set") << std::endl;

	for(size_t i = 0; i < _nrParameters[idxBin]; ++i) {
		out << "    [" << i << "] ";
		if(_parametersLimitedLower[idxBin][i] and _parametersLimitedUpper[idxBin][i]) {
			out << "limits: " << _parametersLimitLower[idxBin][i] << "-" << _parametersLimitUpper[idxBin][i] << " GeV/c^2";
		} else if(_parametersLimitedLower[idxBin][i]) {
			out << "lower limit: " << _parametersLimitLower[idxBin][i] << " GeV/c^2";
		} else if(_parametersLimitedUpper[idxBin][i]) {
			out << "upper limit: " << _parametersLimitUpper[idxBin][i] << " GeV/c^2";
		} else {
			out << "unlimited";
		}
		out << (_parametersFixed[idxBin][i] ? " (FIXED) " : "") << std::endl;
	}

	return out;
}
