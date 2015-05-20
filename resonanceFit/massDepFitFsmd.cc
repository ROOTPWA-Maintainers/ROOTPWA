///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
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

#include "massDepFitCache.h"
#include "massDepFitParameters.h"
#include "reportingUtils.hpp"
#include "yamlCppUtils.hpp"


rpwa::massDepFit::fsmd::fsmd(const size_t id)
	: _id(id),
	  _function(NULL),
	  _nrParameters(0)
{
}


rpwa::massDepFit::fsmd::~fsmd()
{
	if(_function != NULL) {
		delete _function;
	}
}


bool
rpwa::massDepFit::fsmd::init(const YAML::Node& configFsmd,
                             rpwa::massDepFit::parameters& fitParameters,
                             rpwa::massDepFit::parameters& fitParametersError,
                             const bool debug)
{
	if(debug) {
		printDebug << "start initializing final-state mass-dependence." << std::endl;
	}

	std::map<std::string, YamlCppUtils::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("formula", YamlCppUtils::TypeString);
	if(not checkIfAllVariablesAreThere(configFsmd, mandatoryArguments)) {
		printErr << "'finalStateMassDependence' does not contain all required variables." << std::endl;
		return false;
	}

	const std::string formula = configFsmd["formula"].as<std::string>();

	if(_function != NULL) {
		delete _function;
	}
	_function = new TFormula("finalStateMassDependence", formula.c_str());

	_nrParameters = _function->GetNpar();

	_parametersFixed.resize(_nrParameters);
	_parametersLimitLower.resize(_nrParameters);
	_parametersLimitedLower.resize(_nrParameters);
	_parametersLimitUpper.resize(_nrParameters);
	_parametersLimitedUpper.resize(_nrParameters);
	_parametersName.resize(_nrParameters);
	_parametersStep.resize(_nrParameters, 0.0001);

	fitParameters.resize(_id+1, 0, _nrParameters, 0);
	fitParametersError.resize(_id+1, 0, _nrParameters, 0);

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		const std::string parName = _function->GetParName(idxParameter);
		_parametersName[idxParameter] = parName;

		if(debug) {
			printDebug << "reading parameter '" << parName << "'." << std::endl;
		}

		const YAML::Node& configParameter = configFsmd[parName];
		if (not configParameter) {
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
		fitParameters.setParameter(_id, idxParameter, parameter);

		if(configParameter["error"]) {
			if(checkVariableType(configParameter["error"], YamlCppUtils::TypeFloat)) {
				const double error = configParameter["error"].as<double>();
				fitParametersError.setParameter(_id, idxParameter, error);
			} else {
				printErr << "variable 'error' for parameter '" << parName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				return false;
			}
		}

		const bool fixed = configParameter["fix"].as<bool>();
		_parametersFixed[idxParameter] = fixed;

		if(configParameter["lower"]) {
			if(checkVariableType(configParameter["lower"], YamlCppUtils::TypeFloat)) {
				_parametersLimitedLower[idxParameter] = true;
				_parametersLimitLower[idxParameter] = configParameter["lower"].as<double>();
			} else {
				printErr << "variable 'lower' for parameter '" << parName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				return false;
			}
		} else {
			_parametersLimitedLower[idxParameter] = false;
		}
		if(configParameter["upper"]) {
			if(checkVariableType(configParameter["upper"], YamlCppUtils::TypeFloat)) {
				_parametersLimitedUpper[idxParameter] = true;
				_parametersLimitUpper[idxParameter] = configParameter["upper"].as<double>();
			} else {
				printErr << "variable 'upper' for parameter '" << parName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				return false;
			}
		} else {
			_parametersLimitedUpper[idxParameter] = false;
		}

		if(configParameter["step"]) {
			if(checkVariableType(configParameter["step"], YamlCppUtils::TypeFloat)) {
				_parametersStep[idxParameter] = configParameter["step"].as<double>();
			} else {
				printErr << "variable 'step' for parameter '" << parName << "' of final-state mass-dependence defined, but not a floating point number." << std::endl;
				return false;
			}
		}
	}

	if(debug) {
		print(printDebug);
		printDebug << "finished initializing final-state mass-dependence." << std::endl;
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

	yamlOutput << YAML::BeginMap;

	yamlOutput << YAML::Key << "formula";
	yamlOutput << YAML::Value << _function->GetTitle();

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		yamlOutput << YAML::Key << _parametersName[idxParameter];
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "val";
		yamlOutput << YAML::Value << fitParameters.getParameter(_id, idxParameter);

		yamlOutput << YAML::Key << "error";
		yamlOutput << YAML::Value << fitParametersError.getParameter(_id, idxParameter);

		if(_parametersLimitedLower[idxParameter]) {
			yamlOutput << YAML::Key << "lower";
			yamlOutput << YAML::Value << _parametersLimitLower[idxParameter];
		}

		if(_parametersLimitedUpper[idxParameter]) {
			yamlOutput << YAML::Key << "upper";
			yamlOutput << YAML::Value << _parametersLimitUpper[idxParameter];
		}

		yamlOutput << YAML::Key << "step";
		yamlOutput << YAML::Value << _parametersStep[idxParameter];

		yamlOutput << YAML::Key << "fix";
		yamlOutput << YAML::Value << _parametersFixed[idxParameter];

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
	bool invalidateCache = false;
	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		if (fitParameters.getParameter(_id, idxParameter) != par[idxParameter]) {
			fitParameters.setParameter(_id, idxParameter, par[idxParameter]);
			invalidateCache = true;
		}
	}

	if (invalidateCache) {
		cache.setComponent(_id, 0, std::numeric_limits<size_t>::max(), 0.);
		cache.setProdAmp(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), 0.);
	}

	return _nrParameters;
}


std::complex<double>
rpwa::massDepFit::fsmd::val(const rpwa::massDepFit::parameters& fitParameters,
                            rpwa::massDepFit::cache& cache,
                            const double mass,
                            const size_t idxMass) const
{
	if(not _function) {
		return 1.;
	}

	if(idxMass != std::numeric_limits<size_t>::max()) {
		const std::complex<double> fsmd = cache.getComponent(_id, 0, idxMass);
		if (fsmd != 0.) {
			return fsmd;
		}
	}

	const std::complex<double> fsmd = _function->EvalPar(&mass, fitParameters.getParameters(_id));

	if(idxMass != std::numeric_limits<size_t>::max()) {
		cache.setComponent(_id, 0, idxMass, fsmd);
	}

	return fsmd;
}


std::ostream&
rpwa::massDepFit::fsmd::print(std::ostream& out) const
{
	out << "final-state mass-dependence (id " << _id << ")" << std::endl;
	out << "formula: " << ((_function != NULL) ? _function->GetTitle() : "not set") << std::endl;

	for(size_t i=0; i<_nrParameters; ++i) {
		out << "    [" << i << "] ";
		if(_parametersLimitedLower[i] && _parametersLimitedUpper[i]) {
			out << "limits: " << _parametersLimitLower[i] << "-" << _parametersLimitUpper[i] << " GeV/c^2";
		} else if(_parametersLimitedLower[i]) {
			out << "lower limit: " << _parametersLimitLower[i] << " GeV/c^2";
		} else if(_parametersLimitedUpper[i]) {
			out << "upper limit: " << _parametersLimitUpper[i] << " GeV/c^2";
		} else {
			out << "unlimited";
		}
		out << (_parametersFixed[i] ? " (FIXED) " : "") << std::endl;
	}

	return out;
}
