#include "data_py.h"

#include <boost/python.hpp>

#include <boostContainers_py.hpp>
#include <rootConverters_py.h>
#include <stlContainers_py.h>

#define RESONANCEFIT_FORWARD_HH_FROM_PYTHON
#include <resonanceFit/data.h>

namespace bp = boost::python;


namespace {


	rpwa::resonanceFit::baseDataPtr
	baseData_constructor(const bp::list& pyNrWaves,
	                     const bp::object& pyWaveNames,
	                     const bp::list& pyNrMassBins,
	                     const bp::object& pyMassBinCenters,
	                     const bp::object& pyPhaseSpaceIntegrals)
	{
		std::vector<size_t> nrWaves;
		if(not rpwa::py::convertBPObjectToVector(pyNrWaves, nrWaves)) {
			throw;
		}

		boost::multi_array<std::string, 2> waveNames;
		if(not rpwa::py::convertBPObjectToMultiArray(pyWaveNames, waveNames)) {
			throw;
		}

		std::vector<size_t> nrMassBins;
		if(not rpwa::py::convertBPObjectToVector(pyNrMassBins, nrMassBins)) {
			throw;
		}

		boost::multi_array<double, 2> massBinCenters;
		if(not rpwa::py::convertBPObjectToMultiArray(pyMassBinCenters, massBinCenters)) {
			throw;
		}

		boost::multi_array<double, 3> phaseSpaceIntegrals;
		if(not rpwa::py::convertBPObjectToMultiArray(pyPhaseSpaceIntegrals, phaseSpaceIntegrals)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::baseData>(nrWaves,
		                                                      waveNames,
		                                                      nrMassBins,
		                                                      massBinCenters,
		                                                      phaseSpaceIntegrals);
	}


	bp::list
	baseData_nrWaves(const rpwa::resonanceFit::baseData& self)
	{
		const std::vector<size_t>& nrWaves = self.nrWaves();

		bp::list pyNrWaves;
		for(std::vector<size_t>::const_iterator it = nrWaves.begin(); it != nrWaves.end(); ++it) {
			pyNrWaves.append(*it);
		}

		return pyNrWaves;
	}


	bp::object
	baseData_waveNames(const rpwa::resonanceFit::baseData& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.waveNames());
	}


	bp::list
	baseData_nrMassBins(const rpwa::resonanceFit::baseData& self)
	{
		const std::vector<size_t>& nrMassBins = self.nrMassBins();

		bp::list pyNrMassBins;
		for(std::vector<size_t>::const_iterator it = nrMassBins.begin(); it != nrMassBins.end(); ++it) {
			pyNrMassBins.append(*it);
		}

		return pyNrMassBins;
	}


	bp::object
	baseData_massBinCenters(const rpwa::resonanceFit::baseData& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.massBinCenters());
	}


	bp::object
	baseData_phaseSpaceIntegrals(const rpwa::resonanceFit::baseData& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.phaseSpaceIntegrals());
	}


	rpwa::resonanceFit::dataPtr
	data_constructor(const bp::tuple& args,
	                 const bp::dict& kwargs)
	{
		std::vector<size_t> nrWaves;
		if(bp::len(args) > 0) {
			if(not rpwa::py::convertBPObjectToVector(args[1], nrWaves)) {
				throw;
			}
		} else if(kwargs.has_key("nrWaves")) {
			if(not rpwa::py::convertBPObjectToVector(kwargs["nrWaves"], nrWaves)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::string, 2> waveNames;
		if(bp::len(args) > 1) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[3])(), waveNames)) {
				throw;
			}
		} else if(kwargs.has_key("waveNames")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["waveNames"])(), waveNames)) {
				throw;
			}
		} else {
			throw;
		}

		std::vector<size_t> nrMassBins;
		if(bp::len(args) > 2) {
			if(not rpwa::py::convertBPObjectToVector(args[1], nrMassBins)) {
				throw;
			}
		} else if(kwargs.has_key("nrMassBins")) {
			if(not rpwa::py::convertBPObjectToVector(kwargs["nrMassBins"], nrMassBins)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<double, 2> massBinCenters;
		if(bp::len(args) > 3) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[2])(), massBinCenters)) {
				throw;
			}
		} else if(kwargs.has_key("massBinCenters")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["massBinCenters"])(), massBinCenters)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::pair<size_t, size_t>, 3> wavePairMassBinLimits;
		if(bp::len(args) > 4) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[4])(), wavePairMassBinLimits)) {
				throw;
			}
		} else if(kwargs.has_key("wavePairMassBinLimits")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["wavePairMassBinLimits"])(), wavePairMassBinLimits)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<double, 3> phaseSpaceIntegrals;
		if(bp::len(args) > 5) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[5])(), phaseSpaceIntegrals)) {
				throw;
			}
		} else if(kwargs.has_key("phaseSpaceIntegrals")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["phaseSpaceIntegrals"])(), phaseSpaceIntegrals)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::complex<double>, 3> productionAmplitudes;
		if(bp::len(args) > 6) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[6])(), productionAmplitudes)) {
				throw;
			}
		} else if(kwargs.has_key("productionAmplitudes")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["productionAmplitudes"])(), productionAmplitudes)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<TMatrixT<double>, 2> productionAmplitudesCovMatInv;
		if(bp::len(args) > 7) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[7])(), productionAmplitudesCovMatInv)) {
				throw;
			}
		} else if(kwargs.has_key("productionAmplitudesCovMatInv")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["productionAmplitudesCovMatInv"])(), productionAmplitudesCovMatInv)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::complex<double>, 4> spinDensityMatrixElements;
		if(bp::len(args) > 8) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[8])(), spinDensityMatrixElements)) {
				throw;
			}
		} else if(kwargs.has_key("spinDensityMatrixElements")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["spinDensityMatrixElements"])(), spinDensityMatrixElements)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<TMatrixT<double>, 2> spinDensityMatrixElementsCovMatInv;
		if(bp::len(args) > 9) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[9])(), spinDensityMatrixElementsCovMatInv)) {
				throw;
			}
		} else if(kwargs.has_key("spinDensityMatrixElementsCovMatInv")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["spinDensityMatrixElementsCovMatInv"])(), spinDensityMatrixElementsCovMatInv)) {
				throw;
			}
		} else {
			throw;
		}

		rpwa::resonanceFit::function::useCovarianceMatrix useCovariance;
		{
			if(bp::len(args) > 10) {
				useCovariance = bp::extract<rpwa::resonanceFit::function::useCovarianceMatrix>(args[10]);
			} else if(kwargs.has_key("useCovariance")) {
				useCovariance = bp::extract<rpwa::resonanceFit::function::useCovarianceMatrix>(kwargs["useCovariance"]);
			} else {
				throw;
			}
		}

		boost::multi_array<std::pair<double, double>, 3> plottingIntensities;
		if(bp::len(args) > 11) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[11])(), plottingIntensities)) {
				throw;
			}
		} else if(kwargs.has_key("plottingIntensities")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["plottingIntensities"])(), plottingIntensities)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::pair<double, double>, 4> plottingSpinDensityMatrixElementsReal;
		if(bp::len(args) > 12) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[12])(), plottingSpinDensityMatrixElementsReal)) {
				throw;
			}
		} else if(kwargs.has_key("plottingSpinDensityMatrixElementsReal")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["plottingSpinDensityMatrixElementsReal"])(), plottingSpinDensityMatrixElementsReal)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::pair<double, double>, 4> plottingSpinDensityMatrixElementsImag;
		if(bp::len(args) > 13) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[13])(), plottingSpinDensityMatrixElementsImag)) {
				throw;
			}
		} else if(kwargs.has_key("plottingSpinDensityMatrixElementsImag")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["plottingSpinDensityMatrixElementsImag"])(), plottingSpinDensityMatrixElementsImag)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::pair<double, double>, 4> plottingPhases;
		if(bp::len(args) > 14) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[14])(), plottingPhases)) {
				throw;
			}
		} else if(kwargs.has_key("plottingPhases")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["plottingPhases"])(), plottingPhases)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::pair<double, double>, 3> sysPlottingIntensities;
		if(bp::len(args) > 15) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[15])(), sysPlottingIntensities)) {
				throw;
			}
		} else if(kwargs.has_key("sysPlottingIntensities")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["sysPlottingIntensities"])(), sysPlottingIntensities)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::pair<double, double>, 4> sysPlottingSpinDensityMatrixElementsReal;
		if(bp::len(args) > 16) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[16])(), sysPlottingSpinDensityMatrixElementsReal)) {
				throw;
			}
		} else if(kwargs.has_key("sysPlottingSpinDensityMatrixElementsReal")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["sysPlottingSpinDensityMatrixElementsReal"])(), sysPlottingSpinDensityMatrixElementsReal)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::pair<double, double>, 4> sysPlottingSpinDensityMatrixElementsImag;
		if(bp::len(args) > 17) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[17])(), sysPlottingSpinDensityMatrixElementsImag)) {
				throw;
			}
		} else if(kwargs.has_key("sysPlottingSpinDensityMatrixElementsImag")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["sysPlottingSpinDensityMatrixElementsImag"])(), sysPlottingSpinDensityMatrixElementsImag)) {
				throw;
			}
		} else {
			throw;
		}

		boost::multi_array<std::pair<double, double>, 4> sysPlottingPhases;
		if(bp::len(args) > 18) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(args[18])(), sysPlottingPhases)) {
				throw;
			}
		} else if(kwargs.has_key("sysPlottingPhases")) {
			if(not rpwa::py::convertBPObjectToMultiArray(bp::extract<bp::object>(kwargs["sysPlottingPhases"])(), sysPlottingPhases)) {
				throw;
			}
		} else {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::data>(nrWaves,
		                                                  waveNames,
		                                                  nrMassBins,
		                                                  massBinCenters,
		                                                  wavePairMassBinLimits,
		                                                  phaseSpaceIntegrals,
		                                                  productionAmplitudes,
		                                                  productionAmplitudesCovMatInv,
		                                                  spinDensityMatrixElements,
		                                                  spinDensityMatrixElementsCovMatInv,
		                                                  useCovariance,
		                                                  plottingIntensities,
		                                                  plottingSpinDensityMatrixElementsReal,
		                                                  plottingSpinDensityMatrixElementsImag,
		                                                  plottingPhases,
		                                                  sysPlottingIntensities,
		                                                  sysPlottingSpinDensityMatrixElementsReal,
		                                                  sysPlottingSpinDensityMatrixElementsImag,
		                                                  sysPlottingPhases);
	}


	bp::object
	data_wavePairMassBinLimits(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.wavePairMassBinLimits());
	}


	bp::object
	data_productionAmplitudes(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.productionAmplitudes());
	}


	bp::object
	data_productionAmplitudesCovMatInv(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.productionAmplitudesCovMatInv());
	}


	bp::object
	data_spinDensityMatrixElements(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.spinDensityMatrixElements());
	}


	bp::object
	data_spinDensityMatrixElementsCovMatInv(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.spinDensityMatrixElementsCovMatInv());
	}


	bp::object
	data_spinDensityMatrixElementsCovMatInvArray(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.spinDensityMatrixElementsCovMatInvArray());
	}


	bp::object
	data_plottingIntensities(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.plottingIntensities());
	}


	bp::object
	data_plottingSpinDensityMatrixElementsReal(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.plottingSpinDensityMatrixElementsReal());
	}


	bp::object
	data_plottingSpinDensityMatrixElementsImag(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.plottingSpinDensityMatrixElementsImag());
	}


	bp::object
	data_plottingPhases(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.plottingPhases());
	}


	bp::object
	data_sysPlottingIntensities(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.sysPlottingIntensities());
	}


	bp::object
	data_sysPlottingSpinDensityMatrixElementsReal(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.sysPlottingSpinDensityMatrixElementsReal());
	}


	bp::object
	data_sysPlottingSpinDensityMatrixElementsImag(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.sysPlottingSpinDensityMatrixElementsImag());
	}


	bp::object
	data_sysPlottingPhases(const rpwa::resonanceFit::data& self)
	{
		return rpwa::py::convertMultiArrayToPy(self.sysPlottingPhases());
	}


}


void rpwa::py::resonanceFit::exportData() {

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

	bp::class_<rpwa::resonanceFit::baseData, rpwa::resonanceFit::baseDataPtr>
		(
			"baseData"
			, bp::no_init
		)

		.def(
			"__init__"
			, bp::make_constructor(baseData_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("nrWaves"),
			                        bp::arg("waveNames"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("phaseSpaceIntegrals")))
		)

		.def(
			"nrBins"
			, &rpwa::resonanceFit::baseData::nrBins
		)

		.def(
			"nrWaves"
			, &baseData_nrWaves
		)

		.def(
			"maxNrWaves"
			, &rpwa::resonanceFit::baseData::maxNrWaves
		)

		.def(
			"waveNames"
			, &baseData_waveNames
		)

		.def(
			"nrMassBins"
			, &baseData_nrMassBins
		)

		.def(
			"maxNrMassBins"
			, &rpwa::resonanceFit::baseData::maxNrMassBins
		)

		.def(
			"massBinCenters"
			, &baseData_massBinCenters
		)

		.def(
			"binsHaveEqualStructure"
			, &rpwa::resonanceFit::baseData::binsHaveEqualStructure
		)

		.def(
			"phaseSpaceIntegrals"
			, &baseData_phaseSpaceIntegrals
		)

		;

	bp::class_<rpwa::resonanceFit::data, rpwa::resonanceFit::dataPtr, bp::bases<rpwa::resonanceFit::baseData> >
		(
			"data"
			, bp::no_init
		)

		.def(
			"wavePairMassBinLimits"
			, &data_wavePairMassBinLimits
		)

		.def(
			"productionAmplitudes"
			, &data_productionAmplitudes
		)

		.def(
			"productionAmplitudesCovMatInv"
			, &data_productionAmplitudesCovMatInv
		)

		.def(
			"spinDensityMatrixElements"
			, &data_spinDensityMatrixElements
		)

		.def(
			"spinDensityMatrixElementsCovMatInv"
			, &data_spinDensityMatrixElementsCovMatInv
		)

		.def(
			"spinDensityMatrixElementsCovMatInvArray"
			, &data_spinDensityMatrixElementsCovMatInvArray
		)

		.def(
			"useCovariance"
			, &rpwa::resonanceFit::data::useCovariance
		)

		.def(
			"plottingIntensities"
			, &data_plottingIntensities
		)

		.def(
			"plottingSpinDensityMatrixElementsReal"
			, &data_plottingSpinDensityMatrixElementsReal
		)

		.def(
			"plottingSpinDensityMatrixElementsImag"
			, &data_plottingSpinDensityMatrixElementsImag
		)

		.def(
			"plottingPhases"
			, &data_plottingPhases
		)

		.def(
			"sysPlottingIntensities"
			, &data_sysPlottingIntensities
		)

		.def(
			"sysPlottingSpinDensityMatrixElementsReal"
			, &data_sysPlottingSpinDensityMatrixElementsReal
		)

		.def(
			"sysPlottingSpinDensityMatrixElementsImag"
			, &data_sysPlottingSpinDensityMatrixElementsImag
		)

		.def(
			"sysPlottingPhases"
			, &data_sysPlottingPhases
		)

		;

	def(
		"createDataObject"
		, bp::raw_function(data_constructor)
	);

}
