#include "resonanceFit_py.h"

#include <boost/python.hpp>

#include <TDirectory.h>

#include <boostContainers_py.hpp>
#include <releaseGil_py.hpp>
#include <rootConverters_py.h>
#include <stlContainers_py.h>

#define RESONANCEFIT_FORWARD_HH_FROM_PYTHON
#include <resonanceFit/cache.h>
#include <resonanceFit/parameters.h>
#include <resonanceFit/resonanceFit.h>

namespace bp = boost::python;


namespace {


	bp::tuple
	readInputModel(const std::string& configFileName,
	               const double maxMassBinCenter,
	               const bool useBranchings)
	{
		rpwa::resonanceFit::inputConstPtr fitInput;
		rpwa::resonanceFit::modelConstPtr fitModel;
		rpwa::resonanceFit::parameters fitParameters;
		rpwa::resonanceFit::parameters fitParametersError;
		std::map<std::string, double> fitQuality;
		std::vector<std::string> freeParameters;

		{
			rpwa::py::releaseGil releaseGil;

			rpwa::resonanceFit::read(configFileName,
			                         maxMassBinCenter,
			                         fitInput,
			                         fitModel,
			                         fitParameters,
			                         fitParametersError,
			                         fitQuality,
			                         freeParameters,
			                         useBranchings);
		}

		bp::dict pyFitQuality;
		for(std::map<std::string, double>::const_iterator it = fitQuality.begin(); it != fitQuality.end(); ++it) {
			pyFitQuality[it->first] = it->second;
		}

		bp::list pyFreeParameters;
		for(std::vector<std::string>::const_iterator it = freeParameters.begin(); it != freeParameters.end(); ++it) {
			pyFreeParameters.append(*it);
		}

		return bp::make_tuple(std::const_pointer_cast<rpwa::resonanceFit::input>(fitInput),
		                      std::const_pointer_cast<rpwa::resonanceFit::model>(fitModel),
		                      fitParameters,
		                      fitParametersError,
		                      pyFitQuality,
		                      pyFreeParameters);
	}


	bp::tuple
	readInputDataModel(const std::string& configFileName,
	                   const bool useBranchings,
	                   const rpwa::resonanceFit::function::useCovarianceMatrix useCovariance,
	                   const std::string& valTreeName,
	                   const std::string& valBranchName)
	{
		rpwa::resonanceFit::inputConstPtr fitInput;
		rpwa::resonanceFit::dataConstPtr fitData;
		rpwa::resonanceFit::modelConstPtr fitModel;
		rpwa::resonanceFit::parameters fitParameters;
		rpwa::resonanceFit::parameters fitParametersError;
		std::map<std::string, double> fitQuality;
		std::vector<std::string> freeParameters;

		{
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 8, 0)
			rpwa::py::releaseGil releaseGil;
#endif

			rpwa::resonanceFit::read(configFileName,
			                         fitInput,
			                         fitData,
			                         fitModel,
			                         fitParameters,
			                         fitParametersError,
			                         fitQuality,
			                         freeParameters,
			                         useBranchings,
			                         useCovariance,
			                         valTreeName,
			                         valBranchName);
		}

		bp::dict pyFitQuality;
		for(std::map<std::string, double>::const_iterator it = fitQuality.begin(); it != fitQuality.end(); ++it) {
			pyFitQuality[it->first] = it->second;
		}

		bp::list pyFreeParameters;
		for(std::vector<std::string>::const_iterator it = freeParameters.begin(); it != freeParameters.end(); ++it) {
			pyFreeParameters.append(*it);
		}

		return bp::make_tuple(std::const_pointer_cast<rpwa::resonanceFit::input>(fitInput),
		                      std::const_pointer_cast<rpwa::resonanceFit::data>(fitData),
		                      std::const_pointer_cast<rpwa::resonanceFit::model>(fitModel),
		                      fitParameters,
		                      fitParametersError,
		                      pyFitQuality,
		                      pyFreeParameters);
	}


	void
	writeConfig(const std::string& configFileName,
	            const bp::object& pyFitInput,
	            const bp::object& pyFitModel,
	            const bp::object& pyFitParameters,
	            const bp::object& pyFitParametersError,
	            const bp::dict& pyFitQuality,
	            const bp::list& pyFreeParameters)
	{
		const rpwa::resonanceFit::inputConstPtr& fitInput = bp::extract<rpwa::resonanceFit::inputPtr>(pyFitInput)();
		const rpwa::resonanceFit::modelConstPtr& fitModel = bp::extract<rpwa::resonanceFit::modelPtr>(pyFitModel)();
		const rpwa::resonanceFit::parameters& fitParameters = bp::extract<rpwa::resonanceFit::parameters&>(pyFitParameters);
		const rpwa::resonanceFit::parameters& fitParametersError = bp::extract<rpwa::resonanceFit::parameters&>(pyFitParametersError);

		std::map<std::string, double> fitQuality;
		if(not rpwa::py::convertBPObjectToMap(pyFitQuality, fitQuality)) {
			PyErr_SetString(PyExc_TypeError, "Input for 'fitQuality' could not be converted to map from strings to doubles.");
			bp::throw_error_already_set();
		}

		std::vector<std::string> freeParameters;
		if(not rpwa::py::convertBPObjectToVector(pyFreeParameters, freeParameters)) {
			PyErr_SetString(PyExc_TypeError, "Input for 'freeParameters' could not be converted to vector of strings.");
			bp::throw_error_already_set();
		}

		{
			rpwa::py::releaseGil releaseGil;

			rpwa::resonanceFit::writeConfig(configFileName,
			                                fitInput,
			                                fitModel,
			                                fitParameters,
			                                fitParametersError,
			                                fitQuality,
			                                freeParameters);
		}
	}


	void
	createPlots(const bp::object& pyFitInput,
	            const bp::object& pyFitData,
	            const bp::object& pyFitModel,
	            const bp::object& pyFitParameters,
	            bp::object& pyCache,
	            PyObject* pyMainDirectory,
	            const bool rangePlotting,
	            const size_t extraBinning)
	{
		const rpwa::resonanceFit::inputConstPtr& fitInput = bp::extract<rpwa::resonanceFit::inputPtr>(pyFitInput)();
		const rpwa::resonanceFit::dataConstPtr& fitData = bp::extract<rpwa::resonanceFit::dataPtr>(pyFitData)();
		const rpwa::resonanceFit::modelConstPtr& fitModel = bp::extract<rpwa::resonanceFit::modelPtr>(pyFitModel)();
		const rpwa::resonanceFit::parameters& fitParameters = bp::extract<rpwa::resonanceFit::parameters&>(pyFitParameters);
		rpwa::resonanceFit::cache& cache = bp::extract<rpwa::resonanceFit::cache&>(pyCache);
		TDirectory* mainDirectory = rpwa::py::convertFromPy<TDirectory*>(pyMainDirectory);

		{
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 8, 0)
			rpwa::py::releaseGil releaseGil;
#endif

			rpwa::resonanceFit::createPlots(fitInput,
			                                fitData,
			                                fitModel,
			                                fitParameters,
			                                cache,
			                                mainDirectory,
			                                rangePlotting,
			                                extraBinning);
		}
	}


}


void rpwa::py::resonanceFit::exportResonanceFit() {

	// the classes and functions from the 'resonanceFit' namespace should
	// not be available from Python at the 'pyRootPwa.core.[...]' level,
	// but should also be included into their own module like
	// 'pyRootPwa.core.resonanceFit.[...]'. So below a submodule
	// 'resonanceFit' is added and the classes and functions are added at
	// that scope.

	bp::scope current;
	std::string submoduleName(bp::extract<std::string>(current.attr("__name__")));
	submoduleName.append(".resonanceFit");

	bp::object submodule(bp::borrowed(PyImport_AddModule(submoduleName.c_str())));
	current.attr("resonanceFit") = submodule;

	// switch the scope to the submodule
	bp::scope submoduleScope = submodule;

	bp::def(
		"readInputModel"
		, readInputModel
		, (bp::arg("configFileName"),
		   bp::arg("maxMassBinCenter"),
		   bp::arg("useBranchings"))
	);

	bp::def(
		"readInputDataModel"
		, readInputDataModel
		, (bp::arg("configFileName"),
		   bp::arg("useBranchings"),
		   bp::arg("useCovariance"),
		   bp::arg("valTreeName") = "pwa",
		   bp::arg("valBranchName") = "fitResult_v2")
	);

	bp::def(
		"writeConfig"
		, writeConfig
		, (bp::arg("configFileName"),
		   bp::arg("fitInput"),
		   bp::arg("fitModel"),
		   bp::arg("fitParameters"),
		   bp::arg("fitParametersError"),
		   bp::arg("fitQuality"),
		   bp::arg("freeParameters"))
	);

	bp::def(
		"createPlots"
		, createPlots
		, (bp::arg("fitInput"),
		   bp::arg("fitData"),
		   bp::arg("fitModel"),
		   bp::arg("fitParameters"),
		   bp::arg("cache"),
		   bp::arg("mainDirectory"),
		   bp::arg("rangePlotting") = false,
		   bp::arg("extraBinning") = 1)
	);

	bp::def(
		"debug"
		, rpwa::resonanceFit::debug
	);

	bp::def(
		"setDebug"
		, rpwa::resonanceFit::setDebug
		, (bp::arg("debug"))
	);

}
