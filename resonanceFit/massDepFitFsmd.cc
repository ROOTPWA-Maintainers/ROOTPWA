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

#include <TFormula.h>

#include "libConfigUtils.hpp"
#include "massDepFitCache.h"
#include "massDepFitParameters.h"
#include "reportingUtils.hpp"


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
rpwa::massDepFit::fsmd::init(const libconfig::Setting* configFsmd,
                             rpwa::massDepFit::parameters& fitParameters,
                             const std::vector<double>& massBinCenters,
                             const bool debug)
{
	if(debug) {
		printDebug << "starting initialization of final-state mass-dependence." << std::endl;
	}

	std::map<std::string, libconfig::Setting::Type> mandatoryArguments;
	boost::assign::insert(mandatoryArguments)
	                     ("formula", libconfig::Setting::TypeString)
	                     ("val", libconfig::Setting::TypeArray)
	                     ("lower", libconfig::Setting::TypeArray)
	                     ("upper", libconfig::Setting::TypeArray)
	                     ("error", libconfig::Setting::TypeArray)
	                     ("fix", libconfig::Setting::TypeArray);
	if(not checkIfAllVariablesAreThere(configFsmd, mandatoryArguments)) {
		printErr << "'finalStateMassDependence' section in configuration file contains errors." << std::endl;
		return false;
	}

	std::string formula;
	configFsmd->lookupValue("formula", formula);

	if(_function != NULL) {
		delete _function;
	}
	_function = new TFormula("finalStateMassDependence", formula.c_str());
	const int nrPar = _function->GetNpar();

	const libconfig::Setting& configFsmdValue = (*configFsmd)["val"];
	if(configFsmdValue.getLength() != nrPar) {
		printErr << "'val' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << std::endl;
		return false;
	}
	if(configFsmdValue.getLength() > 0 and configFsmdValue[0].getType() != libconfig::Setting::TypeFloat) {
		printErr << "'val' in 'finalStateMassDependence' has to be made up of numbers." << std::endl;
		return false;
	}

	const libconfig::Setting& configFsmdLower = (*configFsmd)["lower"];
	if(configFsmdLower.getLength() != nrPar) {
		printErr << "'lower' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << std::endl;
		return false;
	}
	if(configFsmdLower.getLength() > 0 and configFsmdLower[0].getType() != libconfig::Setting::TypeFloat) {
		printErr << "'lower' in 'finalStateMassDependence' has to be made up of numbers." << std::endl;
		return false;
	}

	const libconfig::Setting& configFsmdUpper = (*configFsmd)["upper"];
	if(configFsmdUpper.getLength() != nrPar) {
		printErr << "'upper' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << std::endl;
		return false;
	}
	if(configFsmdUpper.getLength() > 0 and configFsmdUpper[0].getType() != libconfig::Setting::TypeFloat) {
		printErr << "'upper' in 'finalStateMassDependence' has to be made up of numbers." << std::endl;
		return false;
	}

	const libconfig::Setting& configFsmdError = (*configFsmd)["error"];
	if(configFsmdError.getLength() != nrPar) {
		printErr << "'error' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << std::endl;
		return false;
	}
	if(configFsmdError.getLength() > 0 and configFsmdError[0].getType() != libconfig::Setting::TypeFloat) {
		printErr << "'error' in 'finalStateMassDependence' has to be made up of numbers." << std::endl;
		return false;
	}

	const libconfig::Setting& configFsmdFix = (*configFsmd)["fix"];
	if(configFsmdFix.getLength() != nrPar) {
		printErr << "'fix' in 'finalStateMassDependence' has to have a length of " << nrPar << "." << std::endl;
		return false;
	}
	if(configFsmdFix.getLength() > 0 and configFsmdFix[0].getType() != libconfig::Setting::TypeBoolean) {
		printErr << "'fix' in 'finalStateMassDependence' has to be made up of booleans." << std::endl;
		return false;
	}

	_nrParameters = _function->GetNpar();

	_parametersFixed.resize(_nrParameters);
	_parametersLimitLower.resize(_nrParameters);
	_parametersLimitedLower.resize(_nrParameters);
	_parametersLimitUpper.resize(_nrParameters);
	_parametersLimitedUpper.resize(_nrParameters);
	_parametersStep.resize(_nrParameters);

	fitParameters.resize(_id+1, 0, _nrParameters, 0);

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		fitParameters.setParameter(_id, idxParameter, configFsmdValue[idxParameter]);
		_parametersFixed[idxParameter] = configFsmdFix[idxParameter];
		_parametersLimitLower[idxParameter] = configFsmdLower[idxParameter];
		_parametersLimitedLower[idxParameter] = true;
		_parametersLimitUpper[idxParameter] = configFsmdUpper[idxParameter];
		_parametersLimitedUpper[idxParameter] = true;
		_parametersStep[idxParameter] = 0.0001;
	}

	if(debug) {
		print(printDebug);
		printDebug << "finished initialization of final-state mass-dependence." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::fsmd::update(const libconfig::Setting* configFsmd,
                               const rpwa::massDepFit::parameters& fitParameters,
                               const rpwa::massDepFit::parameters& fitParametersError,
                               const bool debug) const
{
	if(debug) {
		printDebug << "starting updating of final-state mass-dependence." << std::endl;
	}

	const libconfig::Setting& configFsmdValue = (*configFsmd)["val"];
	const libconfig::Setting& configFsmdError = (*configFsmd)["error"];

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		configFsmdValue[idxParameter] = fitParameters.getParameter(_id, idxParameter);
		configFsmdError[idxParameter] = fitParametersError.getParameter(_id, idxParameter);
	}

	if(debug) {
		printDebug << "finished updating of final-state mass-dependence." << std::endl;
	}

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
