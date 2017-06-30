#include "components_py.h"

#include <boost/python.hpp>

#include <boostContainers_py.hpp>
#include <stlContainers_py.h>

#define RESONANCEFIT_FORWARD_HH_FROM_PYTHON
#include <resonanceFit/components.h>
#include <resonanceFit/forward.h>

namespace bp = boost::python;


namespace {


	std::shared_ptr<rpwa::resonanceFit::component::channel>
	component_channel_constructor(const std::string& waveName,
	                              const bp::list& pyWaveIndices,
	                              const bp::list& pyNrMassBins,
	                              const bp::object& pyMassBinCenters,
	                              const bp::object& pyPhaseSpaceIntegrals)
	{
		std::vector<size_t> waveIndices;
		if(not rpwa::py::convertBPObjectToVector(pyWaveIndices, waveIndices)) {
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

		return std::make_shared<rpwa::resonanceFit::component::channel>(waveName,
		                                                                waveIndices,
		                                                                nrMassBins,
		                                                                massBinCenters,
		                                                                phaseSpaceIntegrals);
	}


	bp::list
	component_channel_getWaveIndices(const rpwa::resonanceFit::component::channel& self)
	{
		const std::vector<size_t>& waveIndices = self.getWaveIndices();

		bp::list pyWaveIndices;
		for(std::vector<size_t>::const_iterator waveIndex = waveIndices.begin(); waveIndex != waveIndices.end(); ++waveIndex) {
			pyWaveIndices.append(*waveIndex);
		}

		return pyWaveIndices;
	}


	bp::list
	component_channel_getBins(const rpwa::resonanceFit::component::channel& self)
	{
		const std::vector<size_t>& bins = self.getBins();

		bp::list pyBins;
		for(std::vector<size_t>::const_iterator bin = bins.begin(); bin != bins.end(); ++bin) {
			pyBins.append(*bin);
		}

		return pyBins;
	}


	size_t
	component_importCouplings(const rpwa::resonanceFit::component& self,
	                          const bp::list& pyPar,
	                          rpwa::resonanceFit::parameters& fitParameters,
	                          rpwa::resonanceFit::cache& cache)
	{
		std::vector<double> par;
		if(not rpwa::py::convertBPObjectToVector(pyPar, par)) {
			throw;
		}

		return self.importCouplings(par.data(),
		                            fitParameters,
		                            cache);
	}


	size_t
	component_importBranchings(const rpwa::resonanceFit::component& self,
	                           const bp::list& pyPar,
	                           rpwa::resonanceFit::parameters& fitParameters,
	                           rpwa::resonanceFit::cache& cache)
	{
		std::vector<double> par;
		if(not rpwa::py::convertBPObjectToVector(pyPar, par)) {
			throw;
		}

		return self.importBranchings(par.data(),
		                             fitParameters,
		                             cache);
	}


	size_t
	component_importParameters(const rpwa::resonanceFit::component& self,
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


	std::complex<double>
	component_val(const rpwa::resonanceFit::component& self,
	              const rpwa::resonanceFit::parameters& fitParameters,
	              rpwa::resonanceFit::cache& cache,
	              const size_t idxBin,
	              const double mass,
	              const size_t idxMass)
	{
		return self.val(fitParameters,
		                cache,
		                idxBin,
		                mass,
		                idxMass);
	}


	bp::list
	fixedWidthBreitWigner_getDefaultParameters()
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = rpwa::resonanceFit::fixedWidthBreitWigner::getDefaultParameters();

		bp::list pyDefaultParameters;
		for(std::vector<rpwa::resonanceFit::parameter>::const_iterator it = defaultParameters.begin(); it != defaultParameters.end(); ++it) {
			pyDefaultParameters.append(*it);
		}

		return pyDefaultParameters;
	}


	std::shared_ptr<rpwa::resonanceFit::fixedWidthBreitWigner>
	fixedWidthBreitWigner_constructor(const size_t id,
	                                  const std::string& name,
	                                  const bp::list& pyParameters,
	                                  const bp::list& pyDecayChannels,
	                                  const bp::list& pyNrMassBins,
	                                  const bp::object& pyMassBinCenters,
	                                  const bool useBranchings)
	{
		std::vector<rpwa::resonanceFit::parameter> parameters;
		if(not rpwa::py::convertBPObjectToVector(pyParameters, parameters)) {
			throw;
		}

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		if(not rpwa::py::convertBPObjectToVector(pyDecayChannels, decayChannels)) {
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

		return std::make_shared<rpwa::resonanceFit::fixedWidthBreitWigner>(id,
		                                                                   name,
		                                                                   parameters,
		                                                                   decayChannels,
		                                                                   nrMassBins,
		                                                                   massBinCenters,
		                                                                   useBranchings);
	}


	bp::list
	dynamicWidthBreitWigner_getDefaultParameters()
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = rpwa::resonanceFit::dynamicWidthBreitWigner::getDefaultParameters();

		bp::list pyDefaultParameters;
		for(std::vector<rpwa::resonanceFit::parameter>::const_iterator it = defaultParameters.begin(); it != defaultParameters.end(); ++it) {
			pyDefaultParameters.append(*it);
		}

		return pyDefaultParameters;
	}


	std::shared_ptr<rpwa::resonanceFit::dynamicWidthBreitWigner>
	dynamicWidthBreitWigner_constructor(const size_t id,
	                                    const std::string& name,
	                                    const bp::list& pyParameters,
	                                    const bp::list& pyDecayChannels,
	                                    const bp::list& pyNrMassBins,
	                                    const bp::object& pyMassBinCenters,
	                                    const bool useBranchings,
	                                    const bp::list& pyBranchingRatio,
	                                    const bp::list& pyRelAngularMom,
	                                    const bp::list& pyMIsobar1,
	                                    const bp::list& pyMIsobar2)
	{
		std::vector<rpwa::resonanceFit::parameter> parameters;
		if(not rpwa::py::convertBPObjectToVector(pyParameters, parameters)) {
			throw;
		}

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		if(not rpwa::py::convertBPObjectToVector(pyDecayChannels, decayChannels)) {
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

		std::vector<double> branchingRatio;
		if(not rpwa::py::convertBPObjectToVector(pyBranchingRatio, branchingRatio)) {
			throw;
		}

		std::vector<int> relAngularMom;
		if(not rpwa::py::convertBPObjectToVector(pyRelAngularMom, relAngularMom)) {
			throw;
		}

		std::vector<double> mIsobar1;
		if(not rpwa::py::convertBPObjectToVector(pyMIsobar1, mIsobar1)) {
			throw;
		}

		std::vector<double> mIsobar2;
		if(not rpwa::py::convertBPObjectToVector(pyMIsobar2, mIsobar2)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::dynamicWidthBreitWigner>(id,
		                                                                     name,
		                                                                     parameters,
		                                                                     decayChannels,
		                                                                     nrMassBins,
		                                                                     massBinCenters,
		                                                                     useBranchings,
		                                                                     branchingRatio,
		                                                                     relAngularMom,
		                                                                     mIsobar1,
		                                                                     mIsobar2);
	}


	bp::list
	dynamicWidthBreitWigner_branchingRatio(const rpwa::resonanceFit::dynamicWidthBreitWigner& self)
	{
		const std::vector<double>& branchingRatio = self.branchingRatio();

		bp::list pyBranchingRatio;
		for(std::vector<double>::const_iterator it = branchingRatio.begin(); it != branchingRatio.end(); ++it) {
			pyBranchingRatio.append(*it);
		}

		return pyBranchingRatio;
	}


	bp::list
	dynamicWidthBreitWigner_mIsobar1(const rpwa::resonanceFit::dynamicWidthBreitWigner& self)
	{
		const std::vector<double>& mIsobar1 = self.mIsobar1();

		bp::list pyMIsobar1;
		for(std::vector<double>::const_iterator it = mIsobar1.begin(); it != mIsobar1.end(); ++it) {
			pyMIsobar1.append(*it);
		}

		return pyMIsobar1;
	}


	bp::list
	dynamicWidthBreitWigner_mIsobar2(const rpwa::resonanceFit::dynamicWidthBreitWigner& self)
	{
		const std::vector<double>& mIsobar2 = self.mIsobar2();

		bp::list pyMIsobar2;
		for(std::vector<double>::const_iterator it = mIsobar2.begin(); it != mIsobar2.end(); ++it) {
			pyMIsobar2.append(*it);
		}

		return pyMIsobar2;
	}


	bp::list
	dynamicWidthBreitWigner_relAngularMom(const rpwa::resonanceFit::dynamicWidthBreitWigner& self)
	{
		const std::vector<int>& relAngularMom = self.relAngularMom();

		bp::list pyRelAngularMom;
		for(std::vector<int>::const_iterator it = relAngularMom.begin(); it != relAngularMom.end(); ++it) {
			pyRelAngularMom.append(*it);
		}

		return pyRelAngularMom;
	}


	bp::list
	integralWidthBreitWigner_getDefaultParameters()
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = rpwa::resonanceFit::integralWidthBreitWigner::getDefaultParameters();

		bp::list pyDefaultParameters;
		for(std::vector<rpwa::resonanceFit::parameter>::const_iterator it = defaultParameters.begin(); it != defaultParameters.end(); ++it) {
			pyDefaultParameters.append(*it);
		}

		return pyDefaultParameters;
	}


	std::shared_ptr<rpwa::resonanceFit::integralWidthBreitWigner>
	integralWidthBreitWigner_constructor(const size_t id,
	                                     const std::string& name,
	                                     const bp::list& pyParameters,
	                                     const bp::list& pyDecayChannels,
	                                     const bp::list& pyNrMassBins,
	                                     const bp::object& pyMassBinCenters,
	                                     const bool useBranchings,
	                                     const bp::list& pyBranchingRatio,
	                                     const bp::list& pyMasses,
	                                     const bp::list& pyValues)
	{
		std::vector<rpwa::resonanceFit::parameter> parameters;
		if(not rpwa::py::convertBPObjectToVector(pyParameters, parameters)) {
			throw;
		}

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		if(not rpwa::py::convertBPObjectToVector(pyDecayChannels, decayChannels)) {
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

		std::vector<double> branchingRatio;
		if(not rpwa::py::convertBPObjectToVector(pyBranchingRatio, branchingRatio)) {
			throw;
		}

		std::vector<std::vector<double> > masses;
		for(ssize_t idxChannel = 0; idxChannel < bp::len(pyMasses); ++idxChannel) {
			std::vector<double> massesChannel;
			if(not rpwa::py::convertBPObjectToVector(pyMasses[idxChannel], massesChannel)) {
				throw;
			}
			masses.push_back(massesChannel);
		}

		std::vector<std::vector<double> > values;
		for(ssize_t idxChannel = 0; idxChannel < bp::len(pyValues); ++idxChannel) {
			std::vector<double> valuesChannel;
			if(not rpwa::py::convertBPObjectToVector(pyValues[idxChannel], valuesChannel)) {
				throw;
			}
			values.push_back(valuesChannel);
		}

		return std::make_shared<rpwa::resonanceFit::integralWidthBreitWigner>(id,
		                                                                      name,
		                                                                      parameters,
		                                                                      decayChannels,
		                                                                      nrMassBins,
		                                                                      massBinCenters,
		                                                                      useBranchings,
		                                                                      branchingRatio,
		                                                                      masses,
		                                                                      values);
	}


	bp::list
	integralWidthBreitWigner_branchingRatio(const rpwa::resonanceFit::integralWidthBreitWigner& self)
	{
		const std::vector<double>& branchingRatio = self.branchingRatio();

		bp::list pyBranchingRatio;
		for(std::vector<double>::const_iterator it = branchingRatio.begin(); it != branchingRatio.end(); ++it) {
			pyBranchingRatio.append(*it);
		}

		return pyBranchingRatio;
	}


	bp::list
	integralWidthBreitWigner_masses(const rpwa::resonanceFit::integralWidthBreitWigner& self)
	{
		const std::vector<std::vector<double> >& masses = self.masses();

		bp::list pyMasses;
		for(std::vector<std::vector<double> >::const_iterator itChannel = masses.begin(); itChannel != masses.end(); ++itChannel) {
			bp::list pyMassesChannel;
			for(std::vector<double>::const_iterator it = itChannel->begin(); it != itChannel->end(); ++it) {
				pyMassesChannel.append(*it);
			}
			pyMasses.append(pyMassesChannel);
		}

		return pyMasses;
	}


	bp::list
	integralWidthBreitWigner_values(const rpwa::resonanceFit::integralWidthBreitWigner& self)
	{
		const std::vector<std::vector<double> >& values = self.values();

		bp::list pyValues;
		for(std::vector<std::vector<double> >::const_iterator itChannel = values.begin(); itChannel != values.end(); ++itChannel) {
			bp::list pyValuesChannel;
			for(std::vector<double>::const_iterator it = itChannel->begin(); it != itChannel->end(); ++it) {
				pyValuesChannel.append(*it);
			}
			pyValues.append(pyValuesChannel);
		}

		return pyValues;
	}


	bp::list
	constantBackground_getDefaultParameters()
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = rpwa::resonanceFit::constantBackground::getDefaultParameters();

		bp::list pyDefaultParameters;
		for(std::vector<rpwa::resonanceFit::parameter>::const_iterator it = defaultParameters.begin(); it != defaultParameters.end(); ++it) {
			pyDefaultParameters.append(*it);
		}

		return pyDefaultParameters;
	}


	std::shared_ptr<rpwa::resonanceFit::constantBackground>
	constantBackground_constructor(const size_t id,
	                               const std::string& name,
	                               const bp::list& pyParameters,
	                               const bp::list& pyDecayChannels,
	                               const bp::list& pyNrMassBins,
	                               const bp::object& pyMassBinCenters,
	                               const bool useBranchings)
	{
		std::vector<rpwa::resonanceFit::parameter> parameters;
		if(not rpwa::py::convertBPObjectToVector(pyParameters, parameters)) {
			throw;
		}

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		if(not rpwa::py::convertBPObjectToVector(pyDecayChannels, decayChannels)) {
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

		return std::make_shared<rpwa::resonanceFit::constantBackground>(id,
		                                                                name,
		                                                                parameters,
		                                                                decayChannels,
		                                                                nrMassBins,
		                                                                massBinCenters,
		                                                                useBranchings);
	}


	bp::list
	exponentialBackground_getDefaultParameters()
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = rpwa::resonanceFit::exponentialBackground::getDefaultParameters();

		bp::list pyDefaultParameters;
		for(std::vector<rpwa::resonanceFit::parameter>::const_iterator it = defaultParameters.begin(); it != defaultParameters.end(); ++it) {
			pyDefaultParameters.append(*it);
		}

		return pyDefaultParameters;
	}


	std::shared_ptr<rpwa::resonanceFit::exponentialBackground>
	exponentialBackground_constructor(const size_t id,
	                                  const std::string& name,
	                                  const bp::list& pyParameters,
	                                  const bp::list& pyDecayChannels,
	                                  const bp::list& pyNrMassBins,
	                                  const bp::object& pyMassBinCenters,
	                                  const bool useBranchings,
	                                  const int relAngularMom,
	                                  const double mIsobar1,
	                                  const double mIsobar2,
	                                  const double exponent)
	{
		std::vector<rpwa::resonanceFit::parameter> parameters;
		if(not rpwa::py::convertBPObjectToVector(pyParameters, parameters)) {
			throw;
		}

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		if(not rpwa::py::convertBPObjectToVector(pyDecayChannels, decayChannels)) {
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

		return std::make_shared<rpwa::resonanceFit::exponentialBackground>(id,
		                                                                   name,
		                                                                   parameters,
		                                                                   decayChannels,
		                                                                   nrMassBins,
		                                                                   massBinCenters,
		                                                                   useBranchings,
		                                                                   relAngularMom,
		                                                                   mIsobar1,
		                                                                   mIsobar2,
		                                                                   exponent);
	}


	bp::list
	tPrimeDependentBackground_getDefaultParameters()
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = rpwa::resonanceFit::tPrimeDependentBackground::getDefaultParameters();

		bp::list pyDefaultParameters;
		for(std::vector<rpwa::resonanceFit::parameter>::const_iterator it = defaultParameters.begin(); it != defaultParameters.end(); ++it) {
			pyDefaultParameters.append(*it);
		}

		return pyDefaultParameters;
	}


	std::shared_ptr<rpwa::resonanceFit::tPrimeDependentBackground>
	tPrimeDependentBackground_constructor(const size_t id,
	                                      const std::string& name,
	                                      const bp::list& pyParameters,
	                                      const bp::list& pyDecayChannels,
	                                      const bp::list& pyNrMassBins,
	                                      const bp::object& pyMassBinCenters,
	                                      const bool useBranchings,
	                                      const bp::list& pyTPrimeMeans,
	                                      const int relAngularMom,
	                                      const double mIsobar1,
	                                      const double mIsobar2,
	                                      const double exponent)
	{
		std::vector<rpwa::resonanceFit::parameter> parameters;
		if(not rpwa::py::convertBPObjectToVector(pyParameters, parameters)) {
			throw;
		}

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		if(not rpwa::py::convertBPObjectToVector(pyDecayChannels, decayChannels)) {
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

		std::vector<double> tPrimeMeans;
		if(not rpwa::py::convertBPObjectToVector(pyTPrimeMeans, tPrimeMeans)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::tPrimeDependentBackground>(id,
		                                                                       name,
		                                                                       parameters,
		                                                                       decayChannels,
		                                                                       nrMassBins,
		                                                                       massBinCenters,
		                                                                       useBranchings,
		                                                                       tPrimeMeans,
		                                                                       relAngularMom,
		                                                                       mIsobar1,
		                                                                       mIsobar2,
		                                                                       exponent);
	}


	bp::list
	exponentialBackgroundIntegral_getDefaultParameters()
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = rpwa::resonanceFit::exponentialBackgroundIntegral::getDefaultParameters();

		bp::list pyDefaultParameters;
		for(std::vector<rpwa::resonanceFit::parameter>::const_iterator it = defaultParameters.begin(); it != defaultParameters.end(); ++it) {
			pyDefaultParameters.append(*it);
		}

		return pyDefaultParameters;
	}


	std::shared_ptr<rpwa::resonanceFit::exponentialBackgroundIntegral>
	exponentialBackgroundIntegral_constructor(const size_t id,
	                                          const std::string& name,
	                                          const bp::list& pyParameters,
	                                          const bp::list& pyDecayChannels,
	                                          const bp::list& pyNrMassBins,
	                                          const bp::object& pyMassBinCenters,
	                                          const bool useBranchings,
	                                          const bp::list& pyMasses,
	                                          const bp::list& pyValues,
	                                          const double exponent)
	{
		std::vector<rpwa::resonanceFit::parameter> parameters;
		if(not rpwa::py::convertBPObjectToVector(pyParameters, parameters)) {
			throw;
		}

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		if(not rpwa::py::convertBPObjectToVector(pyDecayChannels, decayChannels)) {
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

		std::vector<double> masses;
		if(not rpwa::py::convertBPObjectToVector(pyMasses, masses)) {
			throw;
		}

		std::vector<double> values;
		if(not rpwa::py::convertBPObjectToVector(pyValues, values)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::exponentialBackgroundIntegral>(id,
		                                                                           name,
		                                                                           parameters,
		                                                                           decayChannels,
		                                                                           nrMassBins,
		                                                                           massBinCenters,
		                                                                           useBranchings,
		                                                                           masses,
		                                                                           values,
		                                                                           exponent);
	}


	bp::list
	exponentialBackgroundIntegral_masses(const rpwa::resonanceFit::exponentialBackgroundIntegral& self)
	{
		const std::vector<double>& masses = self.masses();

		bp::list pyMasses;
		for(std::vector<double>::const_iterator it = masses.begin(); it != masses.end(); ++it) {
			pyMasses.append(*it);
		}

		return pyMasses;
	}


	bp::list
	exponentialBackgroundIntegral_values(const rpwa::resonanceFit::exponentialBackgroundIntegral& self)
	{
		const std::vector<double>& values = self.values();

		bp::list pyValues;
		for(std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it) {
			pyValues.append(*it);
		}

		return pyValues;
	}


	bp::list
	tPrimeDependentBackgroundIntegral_getDefaultParameters()
	{
		const std::vector<rpwa::resonanceFit::parameter> defaultParameters = rpwa::resonanceFit::tPrimeDependentBackgroundIntegral::getDefaultParameters();

		bp::list pyDefaultParameters;
		for(std::vector<rpwa::resonanceFit::parameter>::const_iterator it = defaultParameters.begin(); it != defaultParameters.end(); ++it) {
			pyDefaultParameters.append(*it);
		}

		return pyDefaultParameters;
	}


	std::shared_ptr<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>
	tPrimeDependentBackgroundIntegral_constructor(const size_t id,
	                                              const std::string& name,
	                                              const bp::list& pyParameters,
	                                              const bp::list& pyDecayChannels,
	                                              const bp::list& pyNrMassBins,
	                                              const bp::object& pyMassBinCenters,
	                                              const bool useBranchings,
	                                              const bp::list& pyTPrimeMeans,
	                                              const bp::list& pyMasses,
	                                              const bp::list& pyValues,
	                                              const double exponent)
	{
		std::vector<rpwa::resonanceFit::parameter> parameters;
		if(not rpwa::py::convertBPObjectToVector(pyParameters, parameters)) {
			throw;
		}

		std::vector<rpwa::resonanceFit::component::channel> decayChannels;
		if(not rpwa::py::convertBPObjectToVector(pyDecayChannels, decayChannels)) {
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

		std::vector<double> tPrimeMeans;
		if(not rpwa::py::convertBPObjectToVector(pyTPrimeMeans, tPrimeMeans)) {
			throw;
		}

		std::vector<double> masses;
		if(not rpwa::py::convertBPObjectToVector(pyMasses, masses)) {
			throw;
		}

		std::vector<double> values;
		if(not rpwa::py::convertBPObjectToVector(pyValues, values)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(id,
		                                                                               name,
		                                                                               parameters,
		                                                                               decayChannels,
		                                                                               nrMassBins,
		                                                                               massBinCenters,
		                                                                               useBranchings,
		                                                                               tPrimeMeans,
		                                                                               masses,
		                                                                               values,
		                                                                               exponent);
	}


	bp::list
	tPrimeDependentBackgroundIntegral_masses(const rpwa::resonanceFit::tPrimeDependentBackgroundIntegral& self)
	{
		const std::vector<double>& masses = self.masses();

		bp::list pyMasses;
		for(std::vector<double>::const_iterator it = masses.begin(); it != masses.end(); ++it) {
			pyMasses.append(*it);
		}

		return pyMasses;
	}


	bp::list
	tPrimeDependentBackgroundIntegral_values(const rpwa::resonanceFit::tPrimeDependentBackgroundIntegral& self)
	{
		const std::vector<double>& values = self.values();

		bp::list pyValues;
		for(std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it) {
			pyValues.append(*it);
		}

		return pyValues;
	}


}


void rpwa::py::resonanceFit::exportComponents() {

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

	{

		// scope for class 'component'
		// - 'channel' is a nested class in 'component'

		bp::scope classScope = bp::class_<rpwa::resonanceFit::component, rpwa::resonanceFit::componentPtr, boost::noncopyable>
			(
				"component"
				, bp::no_init
			)

			.def(
				bp::self_ns::str(bp::self)
			)

			.def(
				"getType"
				, &rpwa::resonanceFit::component::getType
			)

			.def(
				"getId"
				, &rpwa::resonanceFit::component::getId
			)

			.def(
				"getName"
				, &rpwa::resonanceFit::component::getName
				, bp::return_value_policy<bp::copy_const_reference>()
			)

			.def(
				"getNrChannels"
				, &rpwa::resonanceFit::component::getNrChannels
			)

			.def(
				"getChannel"
				, &rpwa::resonanceFit::component::getChannel
				, (bp::arg("idxDecayChannel"))
				, bp::return_internal_reference<1>()
			)

			.def(
				"getChannelFromCouplingIdx"
				, &rpwa::resonanceFit::component::getChannelFromCouplingIdx
				, (bp::arg("idxCoupling"))
				, bp::return_internal_reference<1>()
			)

			.def(
				"getChannelFromBranchingIdx"
				, &rpwa::resonanceFit::component::getChannelFromBranchingIdx
				, (bp::arg("idxBranching"))
				, bp::return_internal_reference<1>()
			)

			.def(
				"getTotalNrChannels"
				, &rpwa::resonanceFit::component::getTotalNrChannels
			)

			.def(
				"mapChannelToCoupling"
				, &rpwa::resonanceFit::component::mapChannelToCoupling
				, (bp::arg("idxDecayChannel"))
			)

			.def(
				"mapChannelToBranching"
				, &rpwa::resonanceFit::component::mapChannelToBranching
				, (bp::arg("idxDecayChannel"))
			)

			.def(
				"mapCouplingToMasterChannel"
				, &rpwa::resonanceFit::component::mapCouplingToMasterChannel
				, (bp::arg("idxCoupling"))
			)

			.def(
				"mapBranchingToMasterChannel"
				, &rpwa::resonanceFit::component::mapBranchingToMasterChannel
				, (bp::arg("idxBranching"))
			)

			.def(
				"getNrCouplings"
				, &rpwa::resonanceFit::component::getNrCouplings
			)

			.def(
				"importCouplings"
				, &component_importCouplings
				, (bp::arg("par"),
				   bp::arg("fitParameters"),
				   bp::arg("cache"))
			)

			.def(
				"getNrBranchings"
				, &rpwa::resonanceFit::component::getNrBranchings
			)

			.def(
				"importBranchings"
				, &component_importBranchings
				, (bp::arg("par"),
				   bp::arg("fitParameters"),
				   bp::arg("cache"))
			)

			.def(
				"getNrParameters"
				, &rpwa::resonanceFit::component::getNrParameters
			)

			.def(
				"importParameters"
				, &component_importParameters
				, (bp::arg("par"),
				   bp::arg("fitParameters"),
				   bp::arg("cache"))
			)

			.def(
				"isBranchingFixed"
				, &rpwa::resonanceFit::component::isBranchingFixed
				, (bp::arg("idxBranching"))
			)

			.def(
				"getParameter"
				, &rpwa::resonanceFit::component::getParameter
				, (bp::arg("idxParameter"))
				, bp::return_internal_reference<1>()
			)

			.def(
				"val"
				, &component_val
				, (bp::arg("fitParameters"),
				   bp::arg("cache"),
				   bp::arg("idxBin"),
				   bp::arg("mass"),
				   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
			)

			.def(
				"getCouplingPhaseSpace"
				, &rpwa::resonanceFit::component::getCouplingPhaseSpace
				, (bp::arg("fitParameters"),
				   bp::arg("cache"),
				   bp::arg("idxChannel"),
				   bp::arg("idxBin"),
				   bp::arg("mass"),
				   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
			)

			;

		bp::class_<rpwa::resonanceFit::component::channel>
			(
				"channel"
				, bp::no_init
			)

			.def(
				"__init__"
				, bp::make_constructor(component_channel_constructor,
				                       bp::default_call_policies(),
				                       (bp::arg("waveName"),
				                        bp::arg("waveIndices"),
				                        bp::arg("nrMassBins"),
				                        bp::arg("massBinCenters"),
				                        bp::arg("phaseSpaceIntegrals")))
			)

			.def(
				"getWaveName"
				, &rpwa::resonanceFit::component::channel::getWaveName
				, bp::return_value_policy<bp::copy_const_reference>()
			)

			.def(
				"getWaveIndices"
				, &component_channel_getWaveIndices
			)

			.def(
				"getBins"
				, &component_channel_getBins
			)

			.def(
				"getPhaseSpaceIntegral"
				, &rpwa::resonanceFit::component::channel::getPhaseSpaceIntegral
				, (bp::arg("idxBin"),
				   bp::arg("mass"),
				   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
			)

			;

	}

	bp::class_<rpwa::resonanceFit::fixedWidthBreitWigner, std::shared_ptr<rpwa::resonanceFit::fixedWidthBreitWigner>, bp::bases<rpwa::resonanceFit::component> >
		(
			"fixedWidthBreitWigner"
			, bp::no_init
		)

		.def(
			"getDefaultParameters"
			, &fixedWidthBreitWigner_getDefaultParameters
		)
		.staticmethod("getDefaultParameters")

		.def(
			"__init__"
			, bp::make_constructor(fixedWidthBreitWigner_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("name"),
			                        bp::arg("parameters"),
			                        bp::arg("decayChannels"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("useBranchings")))
		)

		;

	bp::class_<rpwa::resonanceFit::dynamicWidthBreitWigner, std::shared_ptr<rpwa::resonanceFit::dynamicWidthBreitWigner>, bp::bases<rpwa::resonanceFit::component> >
		(
			"dynamicWidthBreitWigner"
			, bp::no_init
		)

		.def(
			"getDefaultParameters"
			, &dynamicWidthBreitWigner_getDefaultParameters
		)
		.staticmethod("getDefaultParameters")

		.def(
			"__init__"
			, bp::make_constructor(dynamicWidthBreitWigner_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("name"),
			                        bp::arg("parameters"),
			                        bp::arg("decayChannels"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("useBranchings"),
			                        bp::arg("branchingRatio"),
			                        bp::arg("relAngularMom"),
			                        bp::arg("mIsobar1"),
			                        bp::arg("mIsobar2")))
		)

		.def(
			"branchingRatio"
			, &dynamicWidthBreitWigner_branchingRatio
		)

		.def(
			"relAngularMom"
			, &dynamicWidthBreitWigner_relAngularMom
		)

		.def(
			"mIsobar1"
			, &dynamicWidthBreitWigner_mIsobar1
		)

		.def(
			"mIsobar2"
			, &dynamicWidthBreitWigner_mIsobar2
		)

		;

	bp::class_<rpwa::resonanceFit::integralWidthBreitWigner, std::shared_ptr<rpwa::resonanceFit::integralWidthBreitWigner>, bp::bases<rpwa::resonanceFit::component> >
		(
			"integralWidthBreitWigner"
			, bp::no_init
		)

		.def(
			"getDefaultParameters"
			, &integralWidthBreitWigner_getDefaultParameters
		)
		.staticmethod("getDefaultParameters")

		.def(
			"__init__"
			, bp::make_constructor(integralWidthBreitWigner_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("name"),
			                        bp::arg("parameters"),
			                        bp::arg("decayChannels"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("useBranchings"),
			                        bp::arg("branchingRatio"),
			                        bp::arg("masses"),
			                        bp::arg("values")))
		)

		.def(
			"branchingRatio"
			, &integralWidthBreitWigner_branchingRatio
		)

		.def(
			"masses"
			, &integralWidthBreitWigner_masses
		)

		.def(
			"values"
			, &integralWidthBreitWigner_values
		)

		;

	bp::class_<rpwa::resonanceFit::constantBackground, std::shared_ptr<rpwa::resonanceFit::constantBackground>, bp::bases<rpwa::resonanceFit::component> >
		(
			"constantBackground"
			, bp::no_init
		)

		.def(
			"getDefaultParameters"
			, &constantBackground_getDefaultParameters
		)
		.staticmethod("getDefaultParameters")

		.def(
			"__init__"
			, bp::make_constructor(constantBackground_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("name"),
			                        bp::arg("parameters"),
			                        bp::arg("decayChannels"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("useBranchings")))
		)

		;

	bp::class_<rpwa::resonanceFit::exponentialBackground, std::shared_ptr<rpwa::resonanceFit::exponentialBackground>, bp::bases<rpwa::resonanceFit::component> >
		(
			"exponentialBackground"
			, bp::no_init
		)

		.def(
			"getDefaultParameters"
			, &exponentialBackground_getDefaultParameters
		)
		.staticmethod("getDefaultParameters")

		.def(
			"__init__"
			, bp::make_constructor(exponentialBackground_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("name"),
			                        bp::arg("parameters"),
			                        bp::arg("decayChannels"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("useBranchings"),
			                        bp::arg("relAngularMom"),
			                        bp::arg("mIsobar1"),
			                        bp::arg("mIsobar2"),
			                        bp::arg("exponent")))
		)

		.def(
			"relAngularMom"
			, &rpwa::resonanceFit::exponentialBackground::relAngularMom
		)

		.def(
			"mIsobar1"
			, &rpwa::resonanceFit::exponentialBackground::mIsobar1
		)

		.def(
			"mIsobar2"
			, &rpwa::resonanceFit::exponentialBackground::mIsobar2
		)

		.def(
			"exponent"
			, &rpwa::resonanceFit::exponentialBackground::exponent
		)

		;

	bp::class_<rpwa::resonanceFit::tPrimeDependentBackground, std::shared_ptr<rpwa::resonanceFit::tPrimeDependentBackground>, bp::bases<rpwa::resonanceFit::component> >
		(
			"tPrimeDependentBackground"
			, bp::no_init
		)

		.def(
			"getDefaultParameters"
			, &tPrimeDependentBackground_getDefaultParameters
		)
		.staticmethod("getDefaultParameters")

		.def(
			"__init__"
			, bp::make_constructor(tPrimeDependentBackground_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("name"),
			                        bp::arg("parameters"),
			                        bp::arg("decayChannels"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("useBranchings"),
			                        bp::arg("tPrimeMeans"),
			                        bp::arg("relAngularMom"),
			                        bp::arg("mIsobar1"),
			                        bp::arg("mIsobar2"),
			                        bp::arg("exponent")))
		)

		.def(
			"relAngularMom"
			, &rpwa::resonanceFit::tPrimeDependentBackground::relAngularMom
		)

		.def(
			"mIsobar1"
			, &rpwa::resonanceFit::tPrimeDependentBackground::mIsobar1
		)

		.def(
			"mIsobar2"
			, &rpwa::resonanceFit::tPrimeDependentBackground::mIsobar2
		)

		.def(
			"exponent"
			, &rpwa::resonanceFit::tPrimeDependentBackground::exponent
		)

		;

	bp::class_<rpwa::resonanceFit::exponentialBackgroundIntegral, std::shared_ptr<rpwa::resonanceFit::exponentialBackgroundIntegral>, bp::bases<rpwa::resonanceFit::component> >
		(
			"exponentialBackgroundIntegral"
			, bp::no_init
		)

		.def(
			"getDefaultParameters"
			, &exponentialBackgroundIntegral_getDefaultParameters
		)
		.staticmethod("getDefaultParameters")

		.def(
			"__init__"
			, bp::make_constructor(exponentialBackgroundIntegral_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("name"),
			                        bp::arg("parameters"),
			                        bp::arg("decayChannels"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("useBranchings"),
			                        bp::arg("masses"),
			                        bp::arg("values"),
			                        bp::arg("exponent")))
		)

		.def(
			"masses"
			, &exponentialBackgroundIntegral_masses
		)

		.def(
			"values"
			, &exponentialBackgroundIntegral_values
		)

		.def(
			"exponent"
			, &rpwa::resonanceFit::exponentialBackgroundIntegral::exponent
		)

		;

	bp::class_<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral, std::shared_ptr<rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>, bp::bases<rpwa::resonanceFit::component> >
		(
			"tPrimeDependentBackgroundIntegral"
			, bp::no_init
		)

		.def(
			"getDefaultParameters"
			, &tPrimeDependentBackgroundIntegral_getDefaultParameters
		)
		.staticmethod("getDefaultParameters")

		.def(
			"__init__"
			, bp::make_constructor(tPrimeDependentBackgroundIntegral_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("id"),
			                        bp::arg("name"),
			                        bp::arg("parameters"),
			                        bp::arg("decayChannels"),
			                        bp::arg("nrMassBins"),
			                        bp::arg("massBinCenters"),
			                        bp::arg("useBranchings"),
			                        bp::arg("tPrimeMeans"),
			                        bp::arg("masses"),
			                        bp::arg("values"),
			                        bp::arg("exponent")))
		)

		.def(
			"masses"
			, &tPrimeDependentBackgroundIntegral_masses
		)

		.def(
			"values"
			, &tPrimeDependentBackgroundIntegral_values
		)

		.def(
			"exponent"
			, &rpwa::resonanceFit::tPrimeDependentBackgroundIntegral::exponent
		)

		;

}
