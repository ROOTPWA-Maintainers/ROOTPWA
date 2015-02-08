//-----------------------------------------------------------
//
// Description:
//      Implementation of final-state mass-dependence
//      see massDepFitFsmd.h for details
//
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#include "massDepFitFsmd.h"

#include <map>

#include <boost/assign/std/vector.hpp>

#include <Math/Minimizer.h>
#include <TFormula.h>

#include "libConfigUtils.hpp"
#include "reportingUtils.hpp"


rpwa::massDepFit::fsmd::fsmd(const size_t id)
	: _id(id),
	  _function(NULL),
	  _functionFixed(false),
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

	_functionFixed = true;

	_nrParameters = _function->GetNpar();

	_parameters.resize(_nrParameters);
	_parametersFixed.resize(_nrParameters);
	_parametersLimitLower.resize(_nrParameters);
	_parametersLimitedLower.resize(_nrParameters);
	_parametersLimitUpper.resize(_nrParameters);
	_parametersLimitedUpper.resize(_nrParameters);
	_parametersStep.resize(_nrParameters);

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		_parameters[idxParameter] = configFsmdValue[idxParameter];
		_parametersFixed[idxParameter] = configFsmdFix[idxParameter];
		_parametersLimitLower[idxParameter] = configFsmdLower[idxParameter];
		_parametersLimitedLower[idxParameter] = true;
		_parametersLimitUpper[idxParameter] = configFsmdUpper[idxParameter];
		_parametersLimitedUpper[idxParameter] = true;
		_parametersStep[idxParameter] = 0.0001;

		if(not _parametersFixed[idxParameter]) {
			_functionFixed = false;
		}
	}

	// if this is final-state mass-dependence with no parameters at all, or
	// with all parameters fixed, pre-calculate the value for each mass bin
	_values.clear();
	if(_functionFixed) {
		const size_t nrMassBins = massBinCenters.size();

		_values.resize(nrMassBins);

		for(size_t idxMassBin=0; idxMassBin<nrMassBins; ++idxMassBin) {
			_values[idxMassBin] = val(massBinCenters[idxMassBin], std::numeric_limits<size_t>::max());
		}
	}

	if(debug) {
		printDebug << "finished initialization of final-state mass-dependence." << std::endl;
	}

	return true;
}


bool
rpwa::massDepFit::fsmd::update(const libconfig::Setting* configFsmd,
                               const ROOT::Math::Minimizer* minimizer,
                               const bool debug) const
{
	if(debug) {
		printDebug << "starting updating of final-state mass-dependence." << std::endl;
	}

	const libconfig::Setting& configFsmdValue = (*configFsmd)["val"];
	const libconfig::Setting& configFsmdError = (*configFsmd)["error"];

	for(size_t idxParameter=0; idxParameter<_nrParameters; ++idxParameter) {
		configFsmdValue[idxParameter] = _parameters[idxParameter];

		std::ostringstream sName;
		sName << "PSP__" << idxParameter;
		const int varIndex = minimizer->VariableIndex(sName.str());
		if(varIndex == -1) {
			printErr << "variable '" << sName.str() << "' used to extract the error for one parameter "
			         << "of the final-state mass-dependence not known to the minimizer." << std::endl;
			return false;
		}
		configFsmdError[idxParameter] = minimizer->Errors()[varIndex];
	}

	if(debug) {
		printDebug << "finished updating of final-state mass-dependence." << std::endl;
	}

	return true;
}


size_t
rpwa::massDepFit::fsmd::getParameters(double* par) const
{
	for(size_t idx=0; idx<_nrParameters; ++idx) {
		par[idx] = _parameters[idx];
	}

	return _nrParameters;
}


size_t
rpwa::massDepFit::fsmd::setParameters(const double* par)
{
	for(size_t idx=0; idx<_nrParameters; ++idx) {
		_parameters[idx] = par[idx];
	}

	return _nrParameters;
}


double
rpwa::massDepFit::fsmd::val(const double mass,
                            const size_t idxMass) const
{
	if(_functionFixed && idxMass != std::numeric_limits<size_t>::max()) {
		return _values[idxMass];
	}

	if(not _function) {
		return 1.;
	}

	return _function->EvalPar(&mass, &_parameters[0]);
}
