///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014-2016 Sebastian Uhl (TUM)
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
//      implementation of the wrapper around the minuit minimizer
//
//-------------------------------------------------------------------------


#include "minimizerRoot.h"

#include <boost/tokenizer.hpp>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <TMatrixT.h>

#include <reportingUtils.hpp>

#include "cache.h"
#include "components.h"
#include "fsmd.h"
#include "function.h"
#include "model.h"
#include "parameters.h"


rpwa::resonanceFit::minimizerRoot::functionAdaptor::functionAdaptor(const rpwa::resonanceFit::functionConstPtr& fitFunction)
	: _fitFunction(fitFunction)
{
}


rpwa::resonanceFit::minimizerRoot::functionAdaptor*
rpwa::resonanceFit::minimizerRoot::functionAdaptor::Clone() const
{
	return new rpwa::resonanceFit::minimizerRoot::functionAdaptor(*this);
}


unsigned int
rpwa::resonanceFit::minimizerRoot::functionAdaptor::NDim() const
{
	return _fitFunction->getNrParameters();
}


unsigned int
rpwa::resonanceFit::minimizerRoot::functionAdaptor::NPoint() const
{
	return _fitFunction->getNrDataPoints();
}


double
rpwa::resonanceFit::minimizerRoot::functionAdaptor::DoEval(const double* par) const
{
	return _fitFunction->chiSquare(par);
}


rpwa::resonanceFit::minimizerRoot::minimizerRoot(const rpwa::resonanceFit::modelConstPtr& fitModel,
                                                 const rpwa::resonanceFit::functionConstPtr& fitFunction,
                                                 const unsigned int maxNmbOfFunctionCalls,
                                                 const std::string minimizerType[],
                                                 const int minimizerStrategy,
                                                 const double minimizerTolerance,
                                                 const bool quiet)
	: _fitModel(fitModel),
	  _functionAdaptor(fitFunction),
	  _maxNmbOfIterations(20000),
	  _maxNmbOfFunctionCalls((maxNmbOfFunctionCalls > 0) ? maxNmbOfFunctionCalls : (5 * _maxNmbOfIterations * fitFunction->getNrParameters()))
{
	// setup minimizer
	printInfo << "creating and setting up minimizer '" << minimizerType[0] << "' "
	          << "using algorithm '" << minimizerType[1] << "'" << std::endl;
	_minimizer.reset(ROOT::Math::Factory::CreateMinimizer(minimizerType[0], minimizerType[1]));
	if(not _minimizer) {
		printErr << "could not create minimizer. exiting." << std::endl;
		throw;
	}
	_minimizer->SetFunction        (_functionAdaptor);
	_minimizer->SetStrategy        (minimizerStrategy);
	_minimizer->SetTolerance       (minimizerTolerance);
	_minimizer->SetPrintLevel      ((quiet) ? 0 : 3);
	_minimizer->SetMaxIterations   (_maxNmbOfIterations);
	_minimizer->SetMaxFunctionCalls(_maxNmbOfFunctionCalls);

	// special for Minuit2
	if(dynamic_cast<ROOT::Minuit2::Minuit2Minimizer*>(_minimizer.get())) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
		((ROOT::Minuit2::Minuit2Minimizer*)_minimizer.get())->SetStorageLevel(0);
#endif
	}
}


rpwa::resonanceFit::minimizerRoot::~minimizerRoot()
{
}


std::map<std::string, double>
rpwa::resonanceFit::minimizerRoot::minimize(std::vector<std::string>& freeParameters,
                                            rpwa::resonanceFit::parameters& fitParameters,
                                            rpwa::resonanceFit::parameters& fitParametersError,
                                            TMatrixT<double>& covarianceMatrix,
                                            rpwa::resonanceFit::cache& cache)
{
	// in case the freeParameters vector is empty, the default release
	// order is used
	if(freeParameters.size() == 0) {
		printWarn << "using default release order of parameters." << std::endl;
		freeParameters.push_back("coupling branching");
		freeParameters.push_back("coupling branching mass m0");
		freeParameters.push_back("*");
		freeParameters.push_back("hesse");
	}

	// keep list of parameters to free
	const size_t nrSteps = freeParameters.size();
	size_t removedSteps = 0;

	bool success = true;
	for(size_t step = 0; step < nrSteps; ++step) {
		// update number of allowed function calls
		_minimizer->SetMaxFunctionCalls(_maxNmbOfFunctionCalls);

		if(freeParameters[step-removedSteps] == "hesse") {
			if(step == 0) {
				printErr << "cannot calculate Hessian matrix without prior minimization." << std::endl;
				throw;
			}

			printInfo << "performing minimization step " << step << ": calculating Hessian matrix." << std::endl;
			success &= _minimizer->Hesse();

			if(not success) {
				printWarn << "calculation of Hessian matrix failed." << std::endl;
			} else {
				printInfo << "calculation of Hessian matrix successful." << std::endl;
			}
		} else {
			// set startvalues
			if(not initParameters(fitParameters, freeParameters[step-removedSteps])) {
				printErr << "error while setting start parameters for step " << step << "." << std::endl;
				return std::map<std::string, double>();
			}

			printInfo << "performing minimization step " << step << ": '" << freeParameters[step-removedSteps] << "' (" << _minimizer->NFree() << " free parameters)." << std::endl;
			success &= _minimizer->Minimize();

			if(not success) {
				printWarn << "minimization failed." << std::endl;
			} else {
				printInfo << "minimization successful." << std::endl;
			}
		}

		// copy current parameters from minimizer
		_fitModel->importParameters(_minimizer->Errors(), fitParametersError, cache);
		_fitModel->importParameters(_minimizer->X(), fitParameters, cache);

		// import covariance matrix after current step
		const unsigned int nmbPar = _functionAdaptor.NDim();
		covarianceMatrix.ResizeTo(nmbPar, nmbPar);
		_minimizer->GetCovMatrix(covarianceMatrix.GetMatrixArray());

		// remove finished release orders
		// - do not remove if the last finished step was the calculation
		//   of the Hessian
		while(step-removedSteps > 0 and freeParameters[step-removedSteps] != "hesse") {
			freeParameters.erase(freeParameters.begin());
			++removedSteps;
		}

		// number of maximal calls was exceeded or reached
		if(_minimizer->NCalls() >= _maxNmbOfFunctionCalls) {
			_maxNmbOfFunctionCalls = 0;

			// invalidate covariance matrix
			covarianceMatrix.ResizeTo(0, 0);

			break;
		}
		_maxNmbOfFunctionCalls -= _minimizer->NCalls();
	}

	printInfo << "minimizer status summary:" << std::endl
	          << "    total number of parameters .......................... " << _minimizer->NDim()             << std::endl
	          << "    number of free parameters ........................... " << _minimizer->NFree()            << std::endl
	          << "    maximum allowed number of iterations ................ " << _minimizer->MaxIterations()    << std::endl
	          << "    maximum allowed number of function calls ............ " << _minimizer->MaxFunctionCalls() << std::endl
	          << "    minimizer status .................................... " << _minimizer->Status()           << std::endl
	          << "    minimizer provides error and error matrix ........... " << _minimizer->ProvidesError()    << std::endl
	          << "    minimizer has performed detailed error validation ... " << _minimizer->IsValidError()     << std::endl
	          << "    estimated distance to minimum ....................... " << _minimizer->Edm()              << std::endl
	          << "    statistical scale used for error calculation ........ " << _minimizer->ErrorDef()         << std::endl
	          << "    minimizer strategy .................................. " << _minimizer->Strategy()         << std::endl
	          << "    absolute tolerance .................................. " << _minimizer->Tolerance()        << std::endl;

	// print results
	std::ostringstream output;
	const unsigned int nmbPar = _functionAdaptor.NDim();
	for(unsigned int i = 0; i<nmbPar; ++i) {
		output << "    parameter [" << std::setw(3) << i << "] "
		       << _minimizer->VariableName(i) << " "
		       << rpwa::maxPrecisionAlign(_minimizer->X()[i]) << " +- " << rpwa::maxPrecisionAlign(_minimizer->Errors()[i])
		       << std::endl;
	}
	printInfo << "minimization result:" << std::endl
	          << output.str();

	// store and print information about fit quality
	std::map<std::string, double> fitQuality;

	fitQuality["minStatus"] = _minimizer->Status() % 10;
	printInfo << "minimizer status = " << fitQuality["minStatus"] << std::endl;

	// if the last step was to calculate the Hessian matrix, and if that
	// finished, then store the status of this calculation
	if(freeParameters.back() == "hesse" and _maxNmbOfFunctionCalls != 0) {
		fitQuality["hesseStatus"] = _minimizer->Status() / 100;
		printInfo << "Hesse status = " << fitQuality["hesseStatus"] << std::endl;

		fitQuality["covMatrixStatus"] = _minimizer->CovMatrixStatus();
	}

	if(_maxNmbOfFunctionCalls == 0) {
		fitQuality["continue"] = 1.0;
		printWarn << "maximum allowed number of function calls exceeded." << std::endl;
	}

	fitQuality["chi2"] = _functionAdaptor(_minimizer->X());
	printInfo << "chi2 =" << rpwa::maxPrecisionAlign(fitQuality["chi2"]) << std::endl;

	fitQuality["ndf"] = _functionAdaptor.NPoint() - _minimizer->NFree();
	printInfo << "ndf = " << fitQuality["ndf"] << std::endl;

	fitQuality["redchi2"] = fitQuality["chi2"] / fitQuality["ndf"];
	printInfo << "chi2/ndf =" << rpwa::maxPrecisionAlign(fitQuality["redchi2"]) << std::endl;

	fitQuality["edm"] = _minimizer->Edm();

	return fitQuality;
}


bool
rpwa::resonanceFit::minimizerRoot::initParameters(const rpwa::resonanceFit::parameters& fitParameters,
                                                  const std::string& freeParameters) const
{
	// tokenize freeParameters string (default separators also include '*' and ',')
	boost::char_separator<char> separators(" \t\n");
	boost::tokenizer<boost::char_separator<char> > tokenizeFreeParameters(freeParameters, separators);

	// reset minimizer
	_minimizer->Clear();

	// changes status of variables (fixed/released)
	// * couplings have to be freed explicitely by adding 'coupling' to freeParameters
	// * branchings also have to be freed explicitely with the keyword 'branching'
	// * additional parameters can be freed with freeParameters
	// * fixed values from config remain fixed

	size_t parcount=0;
	// first add all couplings
	for(size_t idxComponent = 0; idxComponent < _fitModel->getNrComponents(); ++idxComponent) {
		const rpwa::resonanceFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
		for(size_t idxCoupling=0; idxCoupling<comp->getNrCouplings(); ++idxCoupling) {
			const rpwa::resonanceFit::component::channel& channel = comp->getChannelFromCouplingIdx(idxCoupling);
			const std::vector<size_t>& bins = channel.getBins();
			for(size_t i = 0; i < bins.size(); ++i) {
				const size_t idxBin = bins[i];
				std::ostringstream prefixBin;
				prefixBin << "coupling__bin"
				          << idxBin;
				std::ostringstream prefixName;
				prefixName << prefixBin.str()
				           << "__"
				           << comp->getName()
				           << "__";
				if(comp->getNrBranchings() > 1) {
					const std::string waveQN = channel.getWaveName().substr(0, channel.getWaveName().find("="));
					prefixName << waveQN;
				} else {
					prefixName << channel.getWaveName();
				}

				bool free = false;
				if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*") != tokenizeFreeParameters.end()
				   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "coupling") != tokenizeFreeParameters.end()
				   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), prefixBin.str()) != tokenizeFreeParameters.end()
				   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), prefixName.str()) != tokenizeFreeParameters.end()) {
					free = true;
				}
				bool fix = not free;

				const std::complex<double> parameter = fitParameters.getCoupling(idxComponent, idxCoupling, idxBin);

				if (fix) {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') fixed to " << parameter.real() << std::endl;
					_minimizer->SetFixedVariable(parcount,
					                             prefixName.str() + "__real",
					                             parameter.real());
					++parcount;

					if(not channel.isAnchor(idxBin)) {
						printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') fixed to " << parameter.imag() << std::endl;
						_minimizer->SetFixedVariable(parcount,
						                             prefixName.str() + "__imag",
						                             parameter.imag());
						++parcount;
					}
				} else {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') set to " << parameter.real() << std::endl;
					_minimizer->SetVariable(parcount,
					                        prefixName.str() + "__real",
					                        parameter.real(),
					                        0.1);
					++parcount;

					if(not channel.isAnchor(idxBin)) {
						printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') set to " << parameter.imag() << std::endl;
						_minimizer->SetVariable(parcount,
						                        prefixName.str() + "__imag",
						                        parameter.imag(),
						                        0.1);
						++parcount;
					}
				}
			}
		} // end loop over channels
	} // end loop over components

	// second eventually add all branchings
	for(size_t idxComponent = 0; idxComponent < _fitModel->getNrComponents(); ++idxComponent) {
		const rpwa::resonanceFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
		for(size_t idxBranching = 0; idxBranching < comp->getNrBranchings(); ++idxBranching) {
			// skip branchings that are always real and fixed to 1
			if(comp->isBranchingFixed(idxBranching)) {
				continue;
			}

			const rpwa::resonanceFit::component::channel& channel = comp->getChannelFromBranchingIdx(idxBranching);
			const std::string waveQN = channel.getWaveName().substr(0, channel.getWaveName().find("="));
			const std::string waveDecay = channel.getWaveName().substr(channel.getWaveName().find("=")+1);
			std::ostringstream prefixName;
			prefixName << "branching__"
			           << comp->getName()
			           << "__"
			           << waveQN
			           << "__"
			           << waveDecay;

			bool free = false;
			if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*") != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "branching") != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), prefixName.str()) != tokenizeFreeParameters.end()) {
				free = true;
			}
			bool fix = not free;

			const std::complex<double> parameter = fitParameters.getBranching(idxComponent, idxBranching);

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') fixed to " << parameter.real() << std::endl;
				_minimizer->SetFixedVariable(parcount,
				                             prefixName.str() + "__real",
				                             parameter.real());
				++parcount;

				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') fixed to " << parameter.imag() << std::endl;
				_minimizer->SetFixedVariable(parcount,
				                             prefixName.str() + "__imag",
				                             parameter.imag());
				++parcount;
			} else {
				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') set to " << parameter.real() << std::endl;
				_minimizer->SetVariable(parcount,
				                        prefixName.str() + "__real",
				                        parameter.real(),
				                        0.1);
				++parcount;

				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') set to " << parameter.imag() << std::endl;
				_minimizer->SetVariable(parcount,
				                        prefixName.str() + "__imag",
				                        parameter.imag(),
				                        0.1);
				++parcount;
			}
		} // end loop over channels
	} // end loop over components

	// third add parameters of the components, i.e. mass and width
	for(size_t idxComponent = 0; idxComponent < _fitModel->getNrComponents(); ++idxComponent) {
		const rpwa::resonanceFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
		for(size_t idxParameter=0; idxParameter<comp->getNrParameters(); ++idxParameter) {
			const rpwa::resonanceFit::parameter& parameter = comp->getParameter(idxParameter);
			const std::string name = comp->getName() + "__" + parameter.name();

			bool free = false;
			if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*") != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName()) != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), parameter.name()) != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), name) != tokenizeFreeParameters.end()) {
				free = true;
			}
			bool fix = not free;
			if(parameter.fixed()) {
				fix = true;
			}

			const double startValue = fitParameters.getParameter(idxComponent, idxParameter);

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << name << "') fixed to " << startValue << std::endl;
				_minimizer->SetFixedVariable(parcount,
				                             name,
				                             startValue);
			} else if(parameter.limitedLower() and parameter.limitedUpper()) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << startValue
				          << " (limited between " << parameter.limitLower()
				          << " and " << parameter.limitUpper() << ")" << std::endl;
				_minimizer->SetLimitedVariable(parcount,
				                               name,
				                               startValue,
				                               parameter.step(),
				                               parameter.limitLower(),
				                               parameter.limitUpper());
			} else if(parameter.limitedLower()) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << startValue
				          << " (limited larger than " << parameter.limitLower() << ")" << std::endl;
				_minimizer->SetLowerLimitedVariable(parcount,
				                                    name,
				                                    startValue,
				                                    parameter.step(),
				                                    parameter.limitLower());
			} else if(parameter.limitedUpper()) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << startValue
				          << " (limited smaller than " << parameter.limitUpper() << ")" << std::endl;
				_minimizer->SetUpperLimitedVariable(parcount,
				                                    name,
				                                    startValue,
				                                    parameter.step(),
				                                    parameter.limitUpper());
			} else {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << startValue << std::endl;
				_minimizer->SetVariable(parcount,
				                        name,
				                        startValue,
				                        parameter.step());
			}
			++parcount;
		}
	} // end loop over components

	// set parameters for final-state mass-dependence
	if(_fitModel->getFsmd()) {
		const rpwa::resonanceFit::fsmdConstPtr& fsmd = _fitModel->getFsmd();
		const size_t maxNrBins = fsmd->isSameFunctionForAllBins() ? 1 : fsmd->getNrBins();
		for(size_t idxBin = 0; idxBin < maxNrBins; ++idxBin) {
			for(size_t idxParameter = 0; idxParameter < fsmd->getNrParameters(idxBin); ++idxParameter) {
				const rpwa::resonanceFit::parameter& parameter = fsmd->getParameter(idxBin, idxParameter);
				std::ostringstream name;
				name << "fsmd__bin"
				     << idxBin
				     << "__"
				     << parameter.name();

				const bool fix = parameter.fixed();

				const double startValue = fitParameters.getParameter(_fitModel->getNrComponents(), fsmd->getParameterIndex(idxBin)+idxParameter);

				if(fix) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') fixed to " << startValue << std::endl;
					_minimizer->SetFixedVariable(parcount,
					                             name.str(),
					                             startValue);
				} else if(parameter.limitedLower() and parameter.limitedUpper()) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << startValue
					          << " (limited between " << parameter.limitLower()
					          << " and " << parameter.limitUpper() << ")" << std::endl;
					_minimizer->SetLimitedVariable(parcount,
					                               name.str(),
					                               startValue,
					                               parameter.step(),
					                               parameter.limitLower(),
					                               parameter.limitUpper());
				} else if(parameter.limitedLower()) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << startValue
					          << " (limited larger than " << parameter.limitLower() << ")" << std::endl;
					_minimizer->SetLowerLimitedVariable(parcount,
					                                    name.str(),
					                                    startValue,
					                                    parameter.step(),
					                                    parameter.limitLower());
				} else if(parameter.limitedUpper()) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << startValue
					          << " (limited smaller than " << parameter.limitUpper() << ")" << std::endl;
					_minimizer->SetUpperLimitedVariable(parcount,
					                                    name.str(),
					                                    startValue,
					                                    parameter.step(),
					                                    parameter.limitUpper());
				} else {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << startValue << std::endl;
					_minimizer->SetVariable(parcount,
					                        name.str(),
					                        startValue,
					                        parameter.step());
				}
				++parcount;
			}
		}
	}

	return true;
}
