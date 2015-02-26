///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014,2015 Sebastian Uhl (TUM)
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


const unsigned int rpwa::massDepFit::minimizerRoot::maxNmbOfIterations    = 20000;
const unsigned int rpwa::massDepFit::minimizerRoot::maxNmbOfFunctionCalls = 2000000;
const bool         rpwa::massDepFit::minimizerRoot::runHesse              = true;


rpwa::massDepFit::minimizerRoot::functionAdaptor::functionAdaptor(const rpwa::massDepFit::function& fitFunction)
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
	return _fitFunction.getNrParameters();
}


double
rpwa::massDepFit::minimizerRoot::functionAdaptor::DoEval(const double* par) const
{
	return _fitFunction.chiSquare(par);
}


rpwa::massDepFit::minimizerRoot::minimizerRoot(const rpwa::massDepFit::model& fitModel,
                                               const rpwa::massDepFit::function& fitFunction,
                                               const std::vector<std::string>& freeParameters,
                                               const std::string minimizerType[],
                                               const int minimizerStrategy,
                                               const double minimizerTolerance,
                                               const bool quiet)
	: _fitModel(fitModel),
	  _functionAdaptor(fitFunction),
	  _freeParameters(freeParameters)
{
	// setup minimizer
	printInfo << "creating and setting up minimizer '" << minimizerType[0] << "' "
	          << "using algorithm '" << minimizerType[1] << "'" << std::endl;
std::cout << _minimizer.get() << std::endl;
	_minimizer.reset(ROOT::Math::Factory::CreateMinimizer(minimizerType[0], minimizerType[1]));
std::cout << _minimizer.get() << std::endl;
	if(_minimizer.get() == NULL) {
		printErr << "could not create minimizer. exiting." << std::endl;
		throw;
	}
	_minimizer->SetFunction        (_functionAdaptor);
	_minimizer->SetStrategy        (minimizerStrategy);
	_minimizer->SetTolerance       (minimizerTolerance);
	_minimizer->SetPrintLevel      ((quiet) ? 0 : 3);
	_minimizer->SetMaxIterations   (maxNmbOfIterations);
	_minimizer->SetMaxFunctionCalls(maxNmbOfFunctionCalls);

	// special for Minuit2
	if(dynamic_cast<ROOT::Minuit2::Minuit2Minimizer*>(_minimizer.get())) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 34, 19)
		((ROOT::Minuit2::Minuit2Minimizer*)_minimizer.get())->SetStorageLevel(0);
#endif
	}
}

unsigned int
rpwa::massDepFit::minimizerRoot::NFree()
{
	return _minimizer->NFree();
}


bool
rpwa::massDepFit::minimizerRoot::minimize(rpwa::massDepFit::parameters& fitParameters,
                                          rpwa::massDepFit::parameters& fitParametersError,
                                          rpwa::massDepFit::cache& cache) const
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

		printInfo << "performing minimization step " << step << ": '" << _freeParameters[step] << "' (" << _minimizer->NFree() << " free parameters)." << std::endl;
		success &= _minimizer->Minimize();

		if(not success) {
			printWarn << "minimization failed." << std::endl;
		} else {
			printInfo << "minimization successful." << std::endl;
		}

		// copy current parameters from minimizer
		_fitModel.importParameters(_minimizer->Errors(), fitParametersError, cache);
		_fitModel.importParameters(_minimizer->X(), fitParameters, cache);
	}

	if (runHesse) {
		printInfo << "calculating Hessian matrix." << std::endl;
		success &= _minimizer->Hesse();

		if(not success) {
			printWarn << "calculation of Hessian matrix failed." << std::endl;
		} else {
			printInfo << "calculation of Hessian matrix successful." << std::endl;
		}

		// copy current parameters from minimizer
		_fitModel.importParameters(_minimizer->Errors(), fitParametersError, cache);
		_fitModel.importParameters(_minimizer->X(), fitParameters, cache);
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

	return success;
}


bool
rpwa::massDepFit::minimizerRoot::initParameters(const rpwa::massDepFit::parameters& fitParameters,
                                                const std::string& freeParameters) const
{
	// tokenize freeParameters string (default deparators also include '*')
	boost::char_separator<char> separators(" ,\t\n");
	boost::tokenizer<boost::char_separator<char> > tokenizeFreeParameters(freeParameters, separators);

	// reset minimizer
	_minimizer->Clear();

	// changes status of variables (fixed/released)
	// * couplings are always free
	// * additional parameters can be freed with freeParameters
	// * branchings also have to be freed explicitely
	// * fixed values from config remain fixed

	size_t parcount=0;
	// first add all couplings
	for(size_t idxComponent=0; idxComponent<_fitModel.getNrComponents(); ++idxComponent) {
		const rpwa::massDepFit::component* comp = _fitModel.getComponent(idxComponent);
		for(size_t idxChannel=0; idxChannel<comp->getNrChannels(); ++idxChannel) {
			// if branchings are used, not every channel has its own coupling
			if(idxChannel != comp->getChannelIdxCoupling(idxChannel)) {
				continue;
			}

			const rpwa::massDepFit::channel& channel = comp->getChannel(idxChannel);
			for(size_t idxBin=0; idxBin<channel.getNrBins(); ++idxBin) {
				std::ostringstream prefixName;
				prefixName << "coupling__bin"
				           << idxBin
				           << "__"
				           << comp->getName()
				           << "__";
				if(_fitModel.useBranchings() && comp->getNrChannels() > 1) {
					const std::string waveQN = channel.getWaveName().substr(0, 7);
					prefixName << waveQN;
				} else {
					prefixName << channel.getWaveName();
				}

				const std::complex<double> parameter = fitParameters.getCoupling(idxComponent, idxChannel, idxBin);

				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') set to " << parameter.real() << std::endl;
				_minimizer->SetVariable(parcount,
				                        prefixName.str() + "__real",
				                        parameter.real(),
				                        0.1);
				++parcount;

				if(not channel.isAnchor()) {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') set to " << parameter.imag() << std::endl;
					_minimizer->SetVariable(parcount,
					                        prefixName.str() + "__imag",
					                        parameter.imag(),
					                        0.1);
					++parcount;
				}
			}
		} // end loop over channels
	} // end loop over components

	// second eventually add all branchings
	if(_fitModel.useBranchings()) {
		for(size_t idxComponent=0; idxComponent<_fitModel.getNrComponents(); ++idxComponent) {
			const rpwa::massDepFit::component* comp = _fitModel.getComponent(idxComponent);
			// branching with idxChannel 0 is always real and fixed to 1
			for(size_t idxChannel=1; idxChannel<comp->getNrChannels(); ++idxChannel) {
				// if branchings are used, not every channel has its own coupling
				if(idxChannel != comp->getChannelIdxBranching(idxChannel)) {
					continue;
				}

				const rpwa::massDepFit::channel& channel = comp->getChannel(idxChannel);
				const std::string waveDecay = channel.getWaveName().substr(7);
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

				const std::complex<double> parameter = fitParameters.getBranching(idxComponent, idxChannel);

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
	for(size_t idxComponent=0; idxComponent<_fitModel.getNrComponents(); ++idxComponent) {
		const rpwa::massDepFit::component* comp = _fitModel.getComponent(idxComponent);
		for(size_t idxParameter=0; idxParameter<comp->getNrParameters(); ++idxParameter) {
			const std::string name = comp->getName() + "__" + comp->getParameterName(idxParameter);

			bool free = false;
			if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*")!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName())!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getParameterName(idxParameter))!=tokenizeFreeParameters.end()
			   || find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName() + "__" + comp->getParameterName(idxParameter))!=tokenizeFreeParameters.end()) {
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
	if(_fitModel.getFsmd() != NULL) {
		const rpwa::massDepFit::fsmd* fsmd = _fitModel.getFsmd();
		for(size_t idxParameter=0; idxParameter<fsmd->getNrParameters(); ++idxParameter) {
			std::ostringstream name;
			name << "PSP__" << idxParameter;

			const bool fix = fsmd->getParameterFixed(idxParameter);

			const double parameter = fitParameters.getParameter(_fitModel.getNrComponents(), idxParameter);

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') fixed to " << parameter << std::endl;
				_minimizer->SetFixedVariable(parcount,
				                             name.str(),
				                             parameter);
			} else if(fsmd->getParameterLimitedLower(idxParameter) && fsmd->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << parameter
				          << " (limited between " << fsmd->getParameterLimitLower(idxParameter)
				          << " and " << fsmd->getParameterLimitUpper(idxParameter) << ")" << std::endl;
				_minimizer->SetLimitedVariable(parcount,
				                               name.str(),
				                               parameter,
				                               fsmd->getParameterStep(idxParameter),
				                               fsmd->getParameterLimitLower(idxParameter),
				                               fsmd->getParameterLimitUpper(idxParameter));
			} else if(fsmd->getParameterLimitedLower(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << parameter
				          << " (limited larger than " << fsmd->getParameterLimitLower(idxParameter) << ")" << std::endl;
				_minimizer->SetLowerLimitedVariable(parcount,
				                                    name.str(),
				                                    parameter,
				                                    fsmd->getParameterStep(idxParameter),
				                                    fsmd->getParameterLimitLower(idxParameter));
			} else if(fsmd->getParameterLimitedUpper(idxParameter)) {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << parameter
				          << " (limited smaller than " << fsmd->getParameterLimitUpper(idxParameter) << ")" << std::endl;
				_minimizer->SetUpperLimitedVariable(parcount,
				                                    name.str(),
				                                    parameter,
				                                    fsmd->getParameterStep(idxParameter),
				                                    fsmd->getParameterLimitUpper(idxParameter));
			} else {
				printInfo << "parameter " << parcount << " ('" << name.str() << "') set to " << parameter << std::endl;
				_minimizer->SetVariable(parcount,
				                        name.str(),
				                        parameter,
				                        fsmd->getParameterStep(idxParameter));
			}
			++parcount;
		}
	}

	return true;
}
