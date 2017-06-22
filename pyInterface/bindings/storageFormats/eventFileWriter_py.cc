#include "eventFileWriter_py.h"

#include <boost/python.hpp>

#include <TFile.h>
#include <TVector3.h>

#include "eventFileWriter.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	bool eventFileWriter_initialize(rpwa::eventFileWriter& self,
	                                PyObject* pyOutputFile,
	                                std::string userString,
	                                rpwa::eventMetadata::eventsTypeEnum eventsType,
	                                bp::object pyProductionKinematicsParticleNames,
	                                bp::object pyDecayKinematicsParticleNames,
	                                bp::dict pyMultibinBoundaries,
	                                bp::object pyAdditionalVariableLabels,
	                                const int& splitlevel = 99,
	                                const int& buffsize = 256000)
	{
		TFile* outputFile = rpwa::py::convertFromPy<TFile*>(pyOutputFile);
		std::vector<std::string> productionKinematicsParticleNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyProductionKinematicsParticleNames, productionKinematicsParticleNames))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for initialStateParticleNames when executing rpwa::eventFileWriter::initialize()");
			bp::throw_error_already_set();
		}
		std::vector<std::string> decayKinematicsParticleNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyDecayKinematicsParticleNames, decayKinematicsParticleNames))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for finalStateParticleNames when executing rpwa::eventFileWriter::initialize()");
			bp::throw_error_already_set();
		}
		rpwa::multibinBoundariesType multibinBoundaries;
		{
			bp::list keys = pyMultibinBoundaries.keys();
			for(int i = 0; i < bp::len(keys); ++i) {
				std::pair<double, double> element;
				if(not rpwa::py::convertBPObjectToPair<double, double>(pyMultibinBoundaries[keys[i]], element)) {
					PyErr_SetString(PyExc_TypeError, "Got invalid pair for multibin boundaries when executing rpwa::eventFileWriter::initialize()");
					bp::throw_error_already_set();
				}
				bp::extract<std::string> getString(keys[i]);
				if(not getString.check()) {
					PyErr_SetString(PyExc_TypeError, "Got invalid key for multibin boundaries when executing rpwa::eventFileWriter::initialize()");
					bp::throw_error_already_set();
				}
				multibinBoundaries[getString()] = element;
			}
		}
		std::vector<std::string> additionalVariableLabels;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyAdditionalVariableLabels, additionalVariableLabels))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for additionalVariableLabels when executing rpwa::eventFileWriter::initialize()");
			bp::throw_error_already_set();
		}
		return self.initialize(*outputFile,
		                       userString,
		                       eventsType,
		                       productionKinematicsParticleNames,
		                       decayKinematicsParticleNames,
		                       multibinBoundaries,
		                       additionalVariableLabels,
		                       splitlevel,
		                       buffsize);
	}

	void eventFileWriter_addEvent(rpwa::eventFileWriter& self,
	                              PyObject* pyProductionKinematicsMomenta,
	                              PyObject* pyDecayKinematicsMomenta,
	                              bp::list  pyAdditionalVariablesToSave)
	{
		std::vector<double> additionalVariablesToSave;
		if(not rpwa::py::convertBPObjectToVector<double>(pyAdditionalVariablesToSave, additionalVariablesToSave)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for additionalVariablesToSave when executing rpwa::eventFileWriter::addEvent()");
			bp::throw_error_already_set();
		}

		TClonesArray* prodKinMomentaClonesArray = rpwa::py::convertFromPy<TClonesArray*>(pyProductionKinematicsMomenta);
		TClonesArray* decayKinMomentaClonesArray = rpwa::py::convertFromPy<TClonesArray*>(pyDecayKinematicsMomenta);
		if ((not prodKinMomentaClonesArray) or (not decayKinMomentaClonesArray)) {

			bp::extract<bp::list> getProdKinMomentaList(pyProductionKinematicsMomenta);
			if(not getProdKinMomentaList.check()) {
				printErr<<"Got invalid input for prodKinMomenta when executing rpwa::eventFileWriter::addEvent()"<<std::endl;
			}
			bp::list prodKinMomentaList = getProdKinMomentaList();

			std::vector<TVector3> productionKinematicsMomenta(len(prodKinMomentaList));
			for(unsigned int i = 0; i < len(prodKinMomentaList); ++i) {
				bp::object item = bp::extract<bp::object>(prodKinMomentaList[i]);
				productionKinematicsMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
			}

			bp::extract<bp::list> getDecayKinMomentaList(pyDecayKinematicsMomenta);
			if(not getDecayKinMomentaList.check()) {
				printErr<<"Got invalid input for decayKinMomenta when executing rpwa::eventFileWriter::addEvent()"<<std::endl;
			}
			bp::list decayKinMomentaList = getDecayKinMomentaList();

			std::vector<TVector3> decayKinematicsMomenta(len(decayKinMomentaList));
			for(unsigned int i = 0; i < len(decayKinMomentaList); ++i) {
				bp::object item = bp::extract<bp::object>(decayKinMomentaList[i]);
				decayKinematicsMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
			}

			self.addEvent(productionKinematicsMomenta, decayKinematicsMomenta, additionalVariablesToSave);
		} else {
			self.addEvent(*prodKinMomentaClonesArray, *decayKinMomentaClonesArray, additionalVariablesToSave);
		}
	}

}


void rpwa::py::exportEventFileWriter() {

	bp::class_<rpwa::eventFileWriter, boost::noncopyable>("eventFileWriter")
		.def(
			"initialize"
			, &eventFileWriter_initialize
			, (bp::arg("outputFile"),
			   bp::arg("userString"),
			   bp::arg("eventsType"),
			   bp::arg("productionKinematicsParticleNames"),
			   bp::arg("decayKinematicsParticleNames"),
			   bp::arg("multibinBoundaries"),
			   bp::arg("additionalVariableLabels"),
			   bp::arg("splitlevel")=99,
			   bp::arg("buffsize")=256000)
		)
		.def(
			"addEvent"
			, &eventFileWriter_addEvent
			, (bp::arg("productionKinematicsMomenta"),
			   bp::arg("decayKinematicsMomenta"),
			   bp::arg("additionalVariablesToSave")=bp::list())
		)
		.def("finalize", &rpwa::eventFileWriter::finalize)
		.def("reset", &rpwa::eventFileWriter::reset)
		.def(
			"initialized"
			, &rpwa::eventFileWriter::initialized
			, bp::return_value_policy<bp::copy_const_reference>()
		);

}
