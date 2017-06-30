#include "model_py.h"

#include <boost/python.hpp>

#include <stlContainers_py.h>

#define RESONANCEFIT_FORWARD_HH_FROM_PYTHON
#include <resonanceFit/components.h>
#include <resonanceFit/fsmd.h>
#include <resonanceFit/model.h>

namespace bp = boost::python;


namespace {


	rpwa::resonanceFit::modelPtr
	model_constructor(const rpwa::resonanceFit::inputPtr& fitInput,
	                  const bp::list& pyComp,
	                  const rpwa::resonanceFit::fsmdPtr& fsmd,
	                  const bp::list& pyAnchorWaveNames,
	                  const bp::list& pyAnchorComponentNames)
	{
		std::vector<rpwa::resonanceFit::componentPtr> componentsNC;
		if(not rpwa::py::convertBPObjectToVector(pyComp, componentsNC)) {
			throw;
		}

		// create a vector with pointers to const objects
		std::vector<rpwa::resonanceFit::componentConstPtr> components;
		for(std::vector<rpwa::resonanceFit::componentPtr>::const_iterator component = componentsNC.begin(); component != componentsNC.end(); ++component) {
			components.push_back(*component);
		}

		std::vector<std::string> anchorWaveNames;
		if(not rpwa::py::convertBPObjectToVector(pyAnchorWaveNames, anchorWaveNames)) {
			throw;
		}

		std::vector<std::string> anchorComponentNames;
		if(not rpwa::py::convertBPObjectToVector(pyAnchorComponentNames, anchorComponentNames)) {
			throw;
		}

		return std::make_shared<rpwa::resonanceFit::model>(fitInput,
		                                                   components,
		                                                   fsmd,
		                                                   anchorWaveNames,
		                                                   anchorComponentNames);
	}


	void
	model_importParameters(const rpwa::resonanceFit::model& self,
	                       const bp::list& pyPar,
	                       rpwa::resonanceFit::parameters& fitParameters,
	                       rpwa::resonanceFit::cache& cache)
	{
		std::vector<double> par;
		if(not rpwa::py::convertBPObjectToVector(pyPar, par)) {
			throw;
		}

		self.importParameters(par.data(),
		                      fitParameters,
		                      cache);
	}


	rpwa::resonanceFit::componentPtr
	model_getComponent(const rpwa::resonanceFit::modelPtr& self,
	                   const size_t idxComponent)
	{
		// cast away const-ness
		// - Python has no concept of const-ness
		// - boost::python does this for all other cases, but cannot
		//   handle smart pointers
		return std::const_pointer_cast<rpwa::resonanceFit::component>(self->getComponent(idxComponent));
	}


	rpwa::resonanceFit::fsmdPtr
	model_getFsmd(const rpwa::resonanceFit::modelPtr& self)
	{
		// cast away const-ness
		// - Python has no concept of const-ness
		// - boost::python does this for all other cases, but cannot
		//   handle smart pointers
		return std::const_pointer_cast<rpwa::resonanceFit::fsmd>(self->getFsmd());
	}


	bp::list
	model_anchorWaveNames(const rpwa::resonanceFit::modelPtr& self)
	{
		bp::list pyAnchorWaveNames;
		for(std::vector<std::string>::const_iterator it = self->anchorWaveNames().begin(); it != self->anchorWaveNames().end(); ++it) {
			pyAnchorWaveNames.append(*it);
		}

		return pyAnchorWaveNames;
	}


	bp::list
	model_anchorComponentNames(const rpwa::resonanceFit::modelPtr& self)
	{
		bp::list pyAnchorComponentNames;
		for(std::vector<std::string>::const_iterator it = self->anchorComponentNames().begin(); it != self->anchorComponentNames().end(); ++it) {
			pyAnchorComponentNames.append(*it);
		}

		return pyAnchorComponentNames;
	}


	bp::list
	model_getComponentChannel(const rpwa::resonanceFit::modelPtr& self,
	                          const size_t idxBin,
	                          const size_t idxWave)
	{
		const std::vector<std::pair<size_t, size_t> >& components = self->getComponentChannel(idxBin, idxWave);
		const size_t nrComponents = components.size();

		bp::list pyComponents;
		for(size_t idxComponents = 0; idxComponents < nrComponents; ++idxComponents) {
			const size_t idxComponent = components[idxComponents].first;
			const size_t idxChannel = components[idxComponents].second;

			pyComponents.append(bp::make_tuple(idxComponent, idxChannel));
		}

		return pyComponents;
	}


}


void rpwa::py::resonanceFit::exportModel() {

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

	bp::class_<rpwa::resonanceFit::model, rpwa::resonanceFit::modelPtr>
		(
			"model"
			, bp::no_init
		)

		.def(
			"__init__"
			, bp::make_constructor(model_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("fitInput"),
			                        bp::arg("comp"),
			                        bp::arg("fsmd"),
			                        bp::arg("anchorWaveNames"),
			                        bp::arg("anchorComponentNames")))
		)

		.def(
			bp::self_ns::str(bp::self)
		)

		.def(
			"isMappingEqualInAllBins"
			, &rpwa::resonanceFit::model::isMappingEqualInAllBins
		)

		.def(
			"getNrParameters"
			, &rpwa::resonanceFit::model::getNrComponents
		)

		.def(
			"importParameters"
			, &model_importParameters
			, (bp::arg("par"),
			   bp::arg("fitParameters"),
			   bp::arg("cache"))
		)

		.def(
			"getNrComponents"
			, &rpwa::resonanceFit::model::getNrComponents
		)

		.def(
			"getComponent"
			, &model_getComponent
			, (bp::arg("idxComponent"))
		)

		.def(
			"getFsmd"
			, &model_getFsmd
		)

		.def(
			"getMaxChannelsInComponent"
			, &rpwa::resonanceFit::model::getMaxChannelsInComponent
		)

		.def(
			"getMaxParametersInComponent"
			, &rpwa::resonanceFit::model::getMaxParametersInComponent
		)

		.def(
			"anchorWaveNames"
			, &model_anchorWaveNames
		)

		.def(
			"anchorComponentNames"
			, &model_anchorComponentNames
		)

		.def(
			"anchorWaveIndex"
			, &rpwa::resonanceFit::model::anchorWaveIndex
			, (bp::arg("idxBin"))
		)

		.def(
			"anchorComponentIndex"
			, &rpwa::resonanceFit::model::anchorComponentIndex
			, (bp::arg("idxBin"))
		)

		.def(
			"isAnchor"
			, &rpwa::resonanceFit::model::isAnchor
			, (bp::arg("idxBin"),
			   bp::arg("idxWave"),
			   bp::arg("idxComponent"))
		)

		.def(
			"getComponentChannel"
			, &model_getComponentChannel
			, (bp::arg("idxBin"),
			   bp::arg("idxWave"))
		)

		.def(
			"productionAmplitude"
			, &rpwa::resonanceFit::model::productionAmplitude
			, (bp::arg("fitParameters"),
			   bp::arg("cache"),
			   bp::arg("idxWave"),
			   bp::arg("idxBin"),
			   bp::arg("mass"),
			   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
		)

		.def(
			"intensity"
			, &rpwa::resonanceFit::model::intensity
			, (bp::arg("fitParameters"),
			   bp::arg("cache"),
			   bp::arg("idxWave"),
			   bp::arg("idxBin"),
			   bp::arg("mass"),
			   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
		)

		.def(
			"phaseAbsolute"
			, &rpwa::resonanceFit::model::phaseAbsolute
			, (bp::arg("fitParameters"),
			   bp::arg("cache"),
			   bp::arg("idxWave"),
			   bp::arg("idxBin"),
			   bp::arg("mass"),
			   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
		)

		.def(
			"spinDensityMatrix"
			, &rpwa::resonanceFit::model::spinDensityMatrix
			, (bp::arg("fitParameters"),
			   bp::arg("cache"),
			   bp::arg("idxWave"),
			   bp::arg("jdxWave"),
			   bp::arg("idxBin"),
			   bp::arg("mass"),
			   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
		)

		.def(
			"phase"
			, &rpwa::resonanceFit::model::phase
			, (bp::arg("fitParameters"),
			   bp::arg("cache"),
			   bp::arg("idxWave"),
			   bp::arg("jdxWave"),
			   bp::arg("idxBin"),
			   bp::arg("mass"),
			   bp::arg("idxMass") = std::numeric_limits<size_t>::max())
		)

		;

}
