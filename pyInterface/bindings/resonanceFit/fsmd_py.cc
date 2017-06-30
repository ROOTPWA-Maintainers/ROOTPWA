#include "fsmd_py.h"

#include <boost/python.hpp>

#include <TFormula.h>
#include <TPython.h>

#include <boostContainers_py.hpp>
#include <rootConverters_py.h>
#include <stlContainers_py.h>

#define RESONANCEFIT_FORWARD_HH_FROM_PYTHON
#include <resonanceFit/cache.h>
#include <resonanceFit/forward.h>
#include <resonanceFit/fsmd.h>
#include <resonanceFit/parameters.h>

namespace bp = boost::python;


namespace {


	rpwa::resonanceFit::fsmdPtr
	fsmd_constructor_1(const size_t id,
	                   const bp::list& pyNrMassBins,
	                   const bp::object& pyMassBinCenters,
	                   const std::string& functionString,
	                   const bp::object& pyParameters)
	{
		std::vector<size_t> nrMassBins;
		if(not rpwa::py::convertBPObjectToVector(pyNrMassBins, nrMassBins)) {
			throw;
		}

		boost::multi_array<double, 2> massBinCenters;
		if(not rpwa::py::convertBPObjectToMultiArray(pyMassBinCenters, massBinCenters)) {
			throw;
		}

		std::shared_ptr<TFormula> function(new TFormula("finalStateMassDependence", functionString.c_str(), false));

		boost::multi_array<rpwa::resonanceFit::parameter, 1> parameters;
		if(not rpwa::py::convertBPObjectToMultiArray(pyParameters, parameters, false)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::fsmd>(id,
		                                                  nrMassBins,
		                                                  massBinCenters,
		                                                  function,
		                                                  parameters);
	}


	rpwa::resonanceFit::fsmdPtr
	fsmd_constructor_2(const size_t id,
	                   const bp::list& pyNrMassBins,
	                   const bp::object& pyMassBinCenters,
	                   const bp::list& pyFunctionsString,
	                   const bp::object& pyParameters)
	{
		std::vector<size_t> nrMassBins;
		if(not rpwa::py::convertBPObjectToVector(pyNrMassBins, nrMassBins)) {
			throw;
		}

		boost::multi_array<double, 2> massBinCenters;
		if(not rpwa::py::convertBPObjectToMultiArray(pyMassBinCenters, massBinCenters)) {
			throw;
		}

		std::vector<std::string> functionsString;
		if(not rpwa::py::convertBPObjectToVector(pyFunctionsString, functionsString)) {
			throw;
		}
		std::vector<std::shared_ptr<TFormula> > functions(functionsString.size());
		std::transform(functionsString.begin(), functionsString.end(), functions.begin(), [](const std::string& functionString){ return std::make_shared<TFormula>("finalStateMassDependence", functionString.c_str(), false); });

		boost::multi_array<rpwa::resonanceFit::parameter, 2> parameters;
		if(not rpwa::py::convertBPObjectToMultiArray(pyParameters, parameters, true)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::fsmd>(id,
		                                                  nrMassBins,
		                                                  massBinCenters,
		                                                  functions,
		                                                  parameters);
	}


	PyObject*
	fsmd_getFunction(const rpwa::resonanceFit::fsmd& self,
	                 const size_t idxBin)
	{
		TFormula* formula = self.getFunction(idxBin).get();
		return TPython::ObjectProxy_FromVoidPtr(formula, formula->ClassName(), false);
	}


	size_t
	fsmd_importParameters(const rpwa::resonanceFit::fsmd& self,
	                      const bp::list& pyPar,
	                      rpwa::resonanceFit::parameters& fitParameters,
	                      rpwa::resonanceFit::cache& cache)
	{
		std::vector<double> par;
		if(not rpwa::py::convertBPObjectToVector(pyPar, par)) {
			throw;
		}

		return self.importParameters(par.data(),
		                             fitParameters,
		                             cache);
	}


}


void rpwa::py::resonanceFit::exportFsmd() {

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

	bp::class_<rpwa::resonanceFit::fsmd, rpwa::resonanceFit::fsmdPtr>
		(
			"fsmd"
			, bp::no_init
		)

		.def(
			"__init__"
			, bp::make_constructor(fsmd_constructor_1,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("function"),
			                        bp::arg("parameters")))
		)

		.def(
			"__init__"
			, bp::make_constructor(fsmd_constructor_2,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("functions"),
			                        bp::arg("parameters")))
		)

		.def(
			bp::self_ns::str(bp::self)
		)

		.def(
			"getId"
			, &rpwa::resonanceFit::fsmd::getId
		)

		.def(
			"isSameFunctionForAllBins"
			, &rpwa::resonanceFit::fsmd::isSameFunctionForAllBins
		)

		.def(
			"getFunction"
			, &fsmd_getFunction
			, (bp::arg("idxBin"))
			, bp::with_custodian_and_ward_postcall<0, 1>()
		)

		.def(
			"getNrBins"
			, &rpwa::resonanceFit::fsmd::getNrBins
		)

		.def(
			"getNrParameters"
			, &rpwa::resonanceFit::fsmd::getNrParameters
			, (bp::arg("idxBin"))
		)

		.def(
			"getParameterIndex"
			, &rpwa::resonanceFit::fsmd::getParameterIndex
			, (bp::arg("idxBin"))
		)

		.def(
			"importParameters"
			, &fsmd_importParameters
			, (bp::arg("par"),
			   bp::arg("fitParameters"),
			   bp::arg("cache"))
		)

		.def(
			"getParameter"
			, &rpwa::resonanceFit::fsmd::getParameter
			, (bp::arg("idxBin"),
			   bp::arg("idxParameter"))
			, bp::return_internal_reference<1>()
		)

		.def(
			"val"
			, &rpwa::resonanceFit::fsmd::val
			, (bp::arg("fitParameters"),
			   bp::arg("cache"),
			   bp::arg("idxBin"),
			   bp::arg("mass"),
			   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
		)

		;

}
