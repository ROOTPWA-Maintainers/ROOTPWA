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


#include "massDepFitMinimizerRoot.h"

#include <boost/tokenizer.hpp>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Minuit2/Minuit2Minimizer.h>

#include "massDepFitCache.h"
#include "massDepFitComponents.h"
#include "massDepFitFsmd.h"
#include "massDepFitFunction.h"
#include "massDepFitModel.h"
#include "massDepFitParameters.h"
#include "reportingUtils.hpp"


rpwa::massDepFit::minimizerRoot::functionAdaptor::functionAdaptor(const rpwa::massDepFit::functionConstPtr& fitFunction)
	: _fitFunction(fitFunction)
{
}


rpwa::massDepFit::minimizerRoot::functionAdaptor*
rpwa::massDepFit::minimizerRoot::functionAdaptor::Clone() const
{
	return new rpwa::massDepFit::minimizerRoot::functionAdaptor(*this);
}


unsigned int
rpwa::massDepFit::minimizerRoot::functionAdaptor::NDim() const
{
	return _fitFunction->getNrParameters();
}


double
rpwa::massDepFit::minimizerRoot::functionAdaptor::DoEval(const double* par) const
{
	return _fitFunction->chiSquare(par);
}


rpwa::massDepFit::minimizerRoot::minimizerRoot(const rpwa::massDepFit::modelConstPtr& fitModel,
                                               const rpwa::massDepFit::functionConstPtr& fitFunction,
                                               const std::vector<std::string>& freeParameters,
                                               const unsigned int maxNmbOfFunctionCalls,
                                               const std::string minimizerType[],
                                               const int minimizerStrategy,
                                               const double minimizerTolerance,
                                               const bool quiet)
	: _fitModel(fitModel),
	  _functionAdaptor(fitFunction),
	  _freeParameters(freeParameters),
	  _maxNmbOfIterations(20000),
	  _maxNmbOfFunctionCalls((maxNmbOfFunctionCalls > 0) ? maxNmbOfFunctionCalls : (5 * _maxNmbOfIterations * fitFunction->getNrParameters())),
	  _runHesse(true)
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


rpwa::massDepFit::minimizerRoot::~minimizerRoot()
{
}


unsigned int
rpwa::massDepFit::minimizerRoot::getNrFreeParameters() const
{
	return _minimizer->NFree();
}


int
rpwa::massDepFit::minimizerRoot::minimize(rpwa::massDepFit::parameters& fitParameters,
                                          rpwa::massDepFit::parameters& fitParametersError,
                                          rpwa::massDepFit::cache& cache)
{
	// keep list of parameters to free
	const size_t nrSteps = _freeParameters.size();

	bool success = true;
	for(size_t step=0; step<nrSteps; ++step) {
		// set startvalues
		if(not initParameters(fitParameters, _freeParameters[step])) {
			printErr << "error while setting start parameters for step " << step << "." << std::endl;
			return false;
		}

		_minimizer->SetMaxFunctionCalls(_maxNmbOfFunctionCalls);

		printInfo << "performing minimization step " << step << ": '" << _freeParameters[step] << "' (" << _minimizer->NFree() << " free parameters)." << std::endl;
		success &= _minimizer->Minimize();

		if(not success) {
			printWarn << "minimization failed." << std::endl;
		} else {
			printInfo << "minimization successful." << std::endl;
		}

		// copy current parameters from minimizer
		_fitModel->importParameters(_minimizer->Errors(), fitParametersError, cache);
		_fitModel->importParameters(_minimizer->X(), fitParameters, cache);

		if(_minimizer->NCalls() >= _maxNmbOfFunctionCalls) {
			_maxNmbOfFunctionCalls = 0;
			break;
		}
		_maxNmbOfFunctionCalls -= _minimizer->NCalls();
	}

	if(_runHesse and _maxNmbOfFunctionCalls != 0) {
		printInfo << "calculating Hessian matrix." << std::endl;
		success &= _minimizer->Hesse();

		if(not success) {
			printWarn << "calculation of Hessian matrix failed." << std::endl;
		} else {
			printInfo << "calculation of Hessian matrix successful." << std::endl;
		}

		// copy current parameters from minimizer
		_fitModel->importParameters(_minimizer->Errors(), fitParametersError, cache);
		_fitModel->importParameters(_minimizer->X(), fitParameters, cache);
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

	return _minimizer->Status();
}


bool
rpwa::massDepFit::minimizerRoot::initParameters(const rpwa::massDepFit::parameters& fitParameters,
                                                const std::string& freeParameters) const
{
	// tokenize freeParameters string (default separators also include '*')
	boost::char_separator<char> separators(" ,\t\n");
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
		const rpwa::massDepFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
		for(size_t idxCoupling=0; idxCoupling<comp->getNrCouplings(); ++idxCoupling) {
			const rpwa::massDepFit::channel& channel = comp->getChannelFromCouplingIdx(idxCoupling);
			const std::vector<size_t>& bins = channel.getBins();
			for(size_t i = 0; i < bins.size(); ++i) {
				const size_t idxBin = bins[i];
				std::ostringstream prefixName;
				prefixName << "coupling__bin"
				           << idxBin
				           << "__"
				           << comp->getName()
				           << "__";
				if(_fitModel->useBranchings() and comp->getNrChannels() > 1) {
					const std::string waveQN = channel.getWaveName().substr(0, channel.getWaveName().find("="));
					prefixName << waveQN;
				} else {
					prefixName << channel.getWaveName();
				}

				bool free = false;
				if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*")!=tokenizeFreeParameters.end()
				   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "coupling")!=tokenizeFreeParameters.end() ) {
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
	if(_fitModel->useBranchings()) {
		for(size_t idxComponent = 0; idxComponent < _fitModel->getNrComponents(); ++idxComponent) {
			const rpwa::massDepFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
			for(size_t idxBranching = 0; idxBranching < comp->getNrBranchings(); ++idxBranching) {
				// skip branchings that are always real and fixed to 1
				if(comp->isBranchingFixed(idxBranching)) {
					continue;
				}

				const rpwa::massDepFit::channel& channel = comp->getChannelFromBranchingIdx(idxBranching);
				const std::string waveDecay = channel.getWaveName().substr(channel.getWaveName().find("=")+1);
				std::ostringstream prefixName;
				prefixName << "branching__"
				           << comp->getName()
				           << "__"
				           << waveDecay;

				bool free = false;
				if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*")!=tokenizeFreeParameters.end()
				   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "branching")!=tokenizeFreeParameters.end() ) {
					free = true;
				}
				bool fix = not free;

				const std::complex<double> parameter = fitParameters.getBranching(idxComponent, idxBranching);

				if (fix) {
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
	}

	// third add parameters of the components, i.e. mass and width
	for(size_t idxComponent = 0; idxComponent < _fitModel->getNrComponents(); ++idxComponent) {
		const rpwa::massDepFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
		for(size_t idxParameter=0; idxParameter<comp->getNrParameters(); ++idxParameter) {
			const std::string name = comp->getName() + "__" + comp->getParameterName(idxParameter);

			bool free = false;
			if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*")!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName())!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getParameterName(idxParameter))!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), name)!=tokenizeFreeParameters.end()) {
				free = true;
			}
			bool fix = not free;
			if(comp->getParameterFixed(idxParameter)) {
				fix = true;
			}

			const double parameter = fitParameters.getParameter(idxComponent, idxParameter);

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << name << "') fixed to " << parameter << std::endl;
				_minimizer->SetFixedVariable(parcount,
				                             name,
				                             parameter);
			} else if(comp->getParameterLimitedLower(idxParameter) && comp->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << parameter
				          << " (limited between " << comp->getParameterLimitLower(idxParameter)
				          << " and " << comp->getParameterLimitUpper(idxParameter) << ")" << std::endl;
				_minimizer->SetLimitedVariable(parcount,
				                               name,
				                               parameter,
				                               comp->getParameterStep(idxParameter),
				                               comp->getParameterLimitLower(idxParameter),
				                               comp->getParameterLimitUpper(idxParameter));
			} else if(comp->getParameterLimitedLower(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << parameter
				          << " (limited larger than " << comp->getParameterLimitLower(idxParameter) << ")" << std::endl;
				_minimizer->SetLowerLimitedVariable(parcount,
				                                    name,
				                                    parameter,
				                                    comp->getParameterStep(idxParameter),
				                                    comp->getParameterLimitLower(idxParameter));
			} else if(comp->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << parameter
				          << " (limited smaller than " << comp->getParameterLimitUpper(idxParameter) << ")" << std::endl;
				_minimizer->SetUpperLimitedVariable(parcount,
				                                    name,
				                                    parameter,
				                                    comp->getParameterStep(idxParameter),
				                                    comp->getParameterLimitUpper(idxParameter));
			} else {
				printInfo << "parameter " << parcount << " ('" << name << "') set to " << parameter << std::endl;
				_minimizer->SetVariable(parcount,
				                        name,
				                        parameter,
				                        comp->getParameterStep(idxParameter));
			}
			++parcount;
		}
	} // end loop over components

	// set parameters for final-state mass-dependence
	if(_fitModel->getFsmd()) {
		const rpwa::massDepFit::fsmdConstPtr& fsmd = _fitModel->getFsmd();
		for(size_t idxBin = 0; idxBin < fsmd->getNrBins(); ++idxBin) {
			for(size_t idxParameter = 0; idxParameter < fsmd->getNrParameters(idxBin); ++idxParameter) {
				std::ostringstream name;
				name << "fsmd__bin"
				     << idxBin
				     << "__"
				     << fsmd->getParameterName(idxBin, idxParameter);

				const bool fix = fsmd->getParameterFixed(idxBin, idxParameter);

				const double parameter = fitParameters.getParameter(_fitModel->getNrComponents(), fsmd->getParameterIndex(idxBin)+idxParameter);

				if(fix) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') fixed to " << parameter << std::endl;
					_minimizer->SetFixedVariable(parcount,
					                             name.str(),
					                             parameter);
				} else if(fsmd->getParameterLimitedLower(idxBin, idxParameter) and fsmd->getParameterLimitedUpper(idxBin, idxParameter)) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << parameter
					          << " (limited between " << fsmd->getParameterLimitLower(idxBin, idxParameter)
					          << " and " << fsmd->getParameterLimitUpper(idxBin, idxParameter) << ")" << std::endl;
					_minimizer->SetLimitedVariable(parcount,
					                               name.str(),
					                               parameter,
					                               fsmd->getParameterStep(idxBin, idxParameter),
					                               fsmd->getParameterLimitLower(idxBin, idxParameter),
					                               fsmd->getParameterLimitUpper(idxBin, idxParameter));
				} else if(fsmd->getParameterLimitedLower(idxBin, idxParameter)) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << parameter
					          << " (limited larger than " << fsmd->getParameterLimitLower(idxBin, idxParameter) << ")" << std::endl;
					_minimizer->SetLowerLimitedVariable(parcount,
					                                    name.str(),
					                                    parameter,
					                                    fsmd->getParameterStep(idxBin, idxParameter),
					                                    fsmd->getParameterLimitLower(idxBin, idxParameter));
				} else if(fsmd->getParameterLimitedUpper(idxBin, idxParameter)) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << parameter
					          << " (limited smaller than " << fsmd->getParameterLimitUpper(idxBin, idxParameter) << ")" << std::endl;
					_minimizer->SetUpperLimitedVariable(parcount,
					                                    name.str(),
					                                    parameter,
					                                    fsmd->getParameterStep(idxBin, idxParameter),
					                                    fsmd->getParameterLimitUpper(idxBin, idxParameter));
				} else {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << parameter << std::endl;
					_minimizer->SetVariable(parcount,
					                        name.str(),
					                        parameter,
					                        fsmd->getParameterStep(idxBin, idxParameter));
				}
				++parcount;
			}
		}
	}

	return true;
}
