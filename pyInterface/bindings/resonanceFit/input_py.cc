#include "input_py.h"

#include <boost/python.hpp>

#include <stlContainers_py.h>

#define RESONANCEFIT_FORWARD_HH_FROM_PYTHON
#include <resonanceFit/forward.h>
#include <resonanceFit/input.h>

namespace bp = boost::python;


namespace {


	rpwa::resonanceFit::inputPtr
	input_constructor(const bp::list& pyBins)
	{
		std::vector<rpwa::resonanceFit::input::bin> bins;
		if(not rpwa::py::convertBPObjectToVector(pyBins, bins)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::input>(bins);
	}


	std::shared_ptr<rpwa::resonanceFit::input::bin>
	bin_constructor(const std::string& fileName,
	                const bp::list& pyWaves,
	                const double tPrimeMean,
	                const double rescaleErrors,
	                const bp::list& pySysFileNames)
	{
		std::vector<rpwa::resonanceFit::input::bin::wave> waves;
		if(not rpwa::py::convertBPObjectToVector(pyWaves, waves)) {
			throw;
		}

		std::vector<std::string> sysFileNames;
		if(not rpwa::py::convertBPObjectToVector(pySysFileNames, sysFileNames)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::input::bin>(fileName,
		                                                        waves,
		                                                        tPrimeMean,
		                                                        rescaleErrors,
		                                                        sysFileNames);
	}


	bp::list
	bin_sysFileNames(const rpwa::resonanceFit::input::bin& self)
	{
		bp::list pySysFileNames;
		for(std::vector<std::string>::const_iterator it = self.sysFileNames().begin(); it != self.sysFileNames().end(); ++it) {
			pySysFileNames.append(*it);
		}

		return pySysFileNames;
	}


	std::shared_ptr<rpwa::resonanceFit::input::bin::wave>
	wave_constructor(const std::string& waveName,
	                 const bp::tuple& pyMassLimits)
	{
		std::pair<double, double> massLimits;
		if(not rpwa::py::convertBPObjectToPair(pyMassLimits, massLimits)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::input::bin::wave>(waveName,
		                                                              massLimits);
	}


	bp::tuple
	wave_massLimits(const rpwa::resonanceFit::input::bin::wave& self)
	{
		return bp::make_tuple(self.massLimits().first, self.massLimits().second);
	}


}


void rpwa::py::resonanceFit::exportInput() {

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

	bp::scope classScope = bp::class_<rpwa::resonanceFit::input, rpwa::resonanceFit::inputPtr>
		(
			"input"
			, bp::no_init
		)

		.def(
			"__init__"
			, bp::make_constructor(input_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("bins")))
		)

		.def(
			bp::self_ns::str(bp::self)
		)

		.def(
			"nrBins"
			, &rpwa::resonanceFit::input::nrBins
		)

		.def(
			"getBin"
			, &rpwa::resonanceFit::input::getBin
			, (bp::arg("idxBin"))
			, bp::return_internal_reference<>()
		)

		;

	bp::scope subClassScope = bp::class_<rpwa::resonanceFit::input::bin, std::shared_ptr<rpwa::resonanceFit::input::bin> >
		(
			"bin"
			, bp::no_init
		)

		.def(
			"__init__"
			, bp::make_constructor(bin_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("fileName"),
			                        bp::arg("waves"),
			                        bp::arg("tPrimeMean"),
			                        bp::arg("rescaleErrors") = 1.0,
			                        bp::arg("sysFileNames") = bp::list()))
		)

		.def(
			bp::self_ns::str(bp::self)
		)

		.def(
			"fileName"
			, &rpwa::resonanceFit::input::bin::fileName
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def(
			"nrWaves"
			, &rpwa::resonanceFit::input::bin::nrWaves
		)

		.def(
			"getWave"
			, &rpwa::resonanceFit::input::bin::getWave
			, (bp::arg("idxWave"))
			, bp::return_internal_reference<>()
		)

		.def(
			"tPrimeMean"
			, &rpwa::resonanceFit::input::bin::tPrimeMean
		)

		.def(
			"rescaleErrors"
			, &rpwa::resonanceFit::input::bin::rescaleErrors
		)

		.def(
			"sysFileNames"
			, &bin_sysFileNames
		)

		;

	bp::class_<rpwa::resonanceFit::input::bin::wave, std::shared_ptr<rpwa::resonanceFit::input::bin::wave> >
		(
			"wave"
			, bp::no_init
		)

		.def(
			"__init__"
			, bp::make_constructor(wave_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("waveName"),
			                        bp::arg("massLimits") = bp::make_tuple(-1.0, -1.0)))
		)

		.def(
			bp::self_ns::str(bp::self)
		)

		.def(
			"waveName"
			, &rpwa::resonanceFit::input::bin::wave::waveName
			, bp::return_value_policy<bp::copy_const_reference>()
		)

		.def(
			"massLimits"
			, &wave_massLimits
		)

		;

}
