#include "function_py.h"

#include <boost/python.hpp>

#include <stlContainers_py.h>

#define RESONANCEFIT_FORWARD_HH_FROM_PYTHON
#include <resonanceFit/cache.h>
#include <resonanceFit/function.h>
#include <resonanceFit/parameters.h>

namespace bp = boost::python;


namespace {


	double
	function_chiSquare_1(const rpwa::resonanceFit::function& self,
	                     const bp::list& pyPar)
	{
		std::vector<double> par;
		if(not rpwa::py::convertBPObjectToVector(pyPar, par)) {
			throw;
		}

		return self.chiSquare(par.data());
	}


	double
	function_chiSquare_2(const rpwa::resonanceFit::function& self,
	                     const rpwa::resonanceFit::parameters& fitParameters,
	                     rpwa::resonanceFit::cache& cache)
	{
		return self.chiSquare(fitParameters, cache);
	}


	double
	function_logLikelihood_1(const rpwa::resonanceFit::function& self,
	                         const bp::list& pyPar)
	{
		std::vector<double> par;
		if(not rpwa::py::convertBPObjectToVector(pyPar, par)) {
			throw;
		}

		return self.logLikelihood(par.data());
	}


	double
	function_logLikelihood_2(const rpwa::resonanceFit::function& self,
	                         const rpwa::resonanceFit::parameters& fitParameters,
	                         rpwa::resonanceFit::cache& cache)
	{
		return self.logLikelihood(fitParameters, cache);
	}


	double
	function_logPriorLikelihood_1(const rpwa::resonanceFit::function& self,
	                              const bp::list& pyPar)
	{
		std::vector<double> par;
		if(not rpwa::py::convertBPObjectToVector(pyPar, par)) {
			throw;
		}

		return self.logPriorLikelihood(par.data());
	}


	double
	function_logPriorLikelihood_2(const rpwa::resonanceFit::function& self,
	                              const rpwa::resonanceFit::parameters& fitParameters)
	{
		return self.logPriorLikelihood(fitParameters);
	}


}


void rpwa::py::resonanceFit::exportFunction() {

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

	bp::scope classScope = bp::class_<rpwa::resonanceFit::function, rpwa::resonanceFit::functionPtr>
		(
			"function"
			, bp::init<const rpwa::resonanceFit::dataPtr&, const rpwa::resonanceFit::modelPtr&, const bool>(bp::args("fitData",
			                                                                                                         "fitModel",
			                                                                                                         "useProductionAmplitudes"))
		)

		.def(
			"getNrParameters"
			, &rpwa::resonanceFit::function::getNrParameters
		)

		.def(
			"getNrDataPoints"
			, &rpwa::resonanceFit::function::getNrDataPoints
		)

		.def(
			"chiSquare"
			, &function_chiSquare_1
			, (bp::arg("par"))
		)

		.def(
			"chiSquare"
			, &function_chiSquare_2
			, (bp::arg("fitParameters"),
			   bp::arg("cache"))
		)

		.def(
			"logLikelihood"
			, &function_logLikelihood_1
			, (bp::arg("par"))
		)

		.def(
			"logLikelihood"
			, &function_logLikelihood_2
			, (bp::arg("fitParameters"),
			   bp::arg("cache"))
		)

		.def(
			"logPriorLikelihood"
			, &function_logPriorLikelihood_1
			, (bp::arg("par"))
		)

		.def(
			"logPriorLikelihood"
			, &function_logPriorLikelihood_2
			, (bp::arg("fitParameters"))
		)

		;

	bp::enum_<rpwa::resonanceFit::function::useCovarianceMatrix>("useCovarianceMatrix")
		.value("useDiagnalElementsOnly", rpwa::resonanceFit::function::useDiagnalElementsOnly)
		.value("useComplexDiagnalElementsOnly", rpwa::resonanceFit::function::useComplexDiagnalElementsOnly)
		.value("useFullCovarianceMatrix", rpwa::resonanceFit::function::useFullCovarianceMatrix)
		.value("useCovarianceMatrixDefault", rpwa::resonanceFit::function::useCovarianceMatrixDefault)
		.export_values();

}
