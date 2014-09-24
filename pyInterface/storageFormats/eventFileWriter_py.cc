
#include "eventFileWriter_py.h"

#include <TFile.h>
#include <TVector3.h>

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	bool eventFileWriter_initialize(rpwa::eventFileWriter& self,
	                                PyObject* pyOutputFile,
	                                std::string userString,
	                                bp::object pyProductionKinematicsParticleNames,
	                                bp::object pyDecayKinematicsParticleNames,
	                                bp::dict pyBinningMap,
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
		std::map<std::string, std::pair<double, double> > binningMap;
		{
			bp::list keys = pyBinningMap.keys();
			for(int i = 0; i < bp::len(keys); ++i) {
				std::pair<double, double> element;
				if(not rpwa::py::convertBPObjectToPair<double, double>(pyBinningMap[keys[i]], element)) {
					PyErr_SetString(PyExc_TypeError, "Got invalid pair for binningMap when executing rpwa::eventFileWriter::initialize()");
					bp::throw_error_already_set();
				}
				bp::extract<std::string> getString(keys[i]);
				if(not getString.check()) {
					PyErr_SetString(PyExc_TypeError, "Got invalid key for binningMap when executing rpwa::eventFileWriter::initialize()");
					bp::throw_error_already_set();
				}
				binningMap[getString()] = element;
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
		                       productionKinematicsParticleNames,
		                       decayKinematicsParticleNames,
		                       binningMap,
		                       additionalVariableLabels,
		                       splitlevel,
		                       buffsize);
	}

	void eventFileWriter_addEvent(rpwa::eventFileWriter& self,
	                              bp::list pyProductionKinematicsMomenta,
	                              bp::list pyDecayKinematicsMomenta,
	                              bp::list pyAdditionalVariablesToSave)
	{
		std::vector<TVector3> productionKinematicsMomenta(len(pyProductionKinematicsMomenta));
		for(unsigned int i = 0; i < len(pyProductionKinematicsMomenta); ++i) {
			bp::object item = bp::extract<bp::object>(pyProductionKinematicsMomenta[i]);
			productionKinematicsMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}
		std::vector<TVector3> decayKinematicsMomenta(len(pyDecayKinematicsMomenta));
		for(unsigned int i = 0; i < len(pyDecayKinematicsMomenta); ++i) {
			bp::object item = bp::extract<bp::object>(pyDecayKinematicsMomenta[i]);
			decayKinematicsMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}
		std::vector<double> additionalVariablesToSave;
		if(not rpwa::py::convertBPObjectToVector<double>(pyAdditionalVariablesToSave, additionalVariablesToSave))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for additionalVariablesToSave when executing rpwa::eventFileWriter::addEvent()");
			bp::throw_error_already_set();
		}
		self.addEvent(productionKinematicsMomenta, decayKinematicsMomenta, additionalVariablesToSave);
	}

}


void rpwa::py::exportEventFileWriter() {

	bp::class_<rpwa::eventFileWriter, boost::noncopyable>("eventFileWriter")
		.def(
			"initialize"
			, &eventFileWriter_initialize
			, (bp::arg("outputFile"),
			   bp::arg("userString"),
			   bp::arg("productionKinematicsParticleNames"),
			   bp::arg("decayKinematicsParticleNames"),
			   bp::arg("binningMap"),
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
