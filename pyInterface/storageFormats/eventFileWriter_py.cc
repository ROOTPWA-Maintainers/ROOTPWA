
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
	                                bp::object pyInitialStateParticleNames,
	                                bp::object pyFinalStateParticleNames,
	                                bp::dict pyBinningMap,
	                                bp::object pyAdditionalVariableLabels,
	                                const std::string& eventTreeName = "rootPwaEvtTree",
	                                const std::string& initialStateMomentaBranchName = "prodKinMomenta",
	                                const std::string& finalStateMomentaBranchName   = "decayKinMomenta",
	                                const std::string& metadataName = "dataMetadata",
	                                const int& splitlevel = 99,
	                                const int& buffsize = 256000)
	{
		TFile* outputFile = rpwa::py::convertFromPy<TFile*>(pyOutputFile);
		std::vector<std::string> initialStateParticleNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyInitialStateParticleNames, initialStateParticleNames))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for initialStateParticleNames when executing rpwa::dataFileWriter::initialize()");
			bp::throw_error_already_set();
		}
		std::vector<std::string> finalStateParticleNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyFinalStateParticleNames, finalStateParticleNames))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for finalStateParticleNames when executing rpwa::dataFileWriter::initialize()");
			bp::throw_error_already_set();
		}
		std::map<std::string, std::pair<double, double> > binningMap;
		{
			bp::list keys = pyBinningMap.keys();
			for(int i = 0; i < bp::len(keys); ++i) {
				std::pair<double, double> element;
				if(not rpwa::py::convertBPObjectToPair<double, double>(pyBinningMap[keys[i]], element)) {
					PyErr_SetString(PyExc_TypeError, "Got invalid pair for binningMap when executing rpwa::dataFileWriter::initialize()");
					bp::throw_error_already_set();
				}
				bp::extract<std::string> getString(keys[i]);
				if(not getString.check()) {
					PyErr_SetString(PyExc_TypeError, "Got invalid key for binningMap when executing rpwa::dataFileWriter::initialize()");
					bp::throw_error_already_set();
				}
				binningMap[getString()] = element;
			}
		}
		std::vector<std::string> additionalVariableLabels;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyAdditionalVariableLabels, additionalVariableLabels))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for additionalVariableLabels when executing rpwa::dataFileWriter::initialize()");
			bp::throw_error_already_set();
		}
		return self.initialize(*outputFile,
		                       userString,
		                       initialStateParticleNames,
		                       finalStateParticleNames,
		                       binningMap,
		                       additionalVariableLabels,
                               eventTreeName,
                               initialStateMomentaBranchName,
                               finalStateMomentaBranchName,
                               metadataName,
                               splitlevel,
                               buffsize);
	}

	void eventFileWriter_addEvent(rpwa::eventFileWriter& self,
	                              bp::list pyInitialStateMomenta,
	                              bp::list pyFinalStateMomenta,
	                              bp::list pyAdditionalVariablesToSave)
	{
		std::vector<TVector3> initialStateMomenta(len(pyInitialStateMomenta));
		for(unsigned int i = 0; i < len(pyInitialStateMomenta); ++i) {
			bp::object item = bp::extract<bp::object>(pyInitialStateMomenta[i]);
			initialStateMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}
		std::vector<TVector3> finalStateMomenta(len(pyFinalStateMomenta));
		for(unsigned int i = 0; i < len(pyFinalStateMomenta); ++i) {
			bp::object item = bp::extract<bp::object>(pyFinalStateMomenta[i]);
			finalStateMomenta[i] = *rpwa::py::convertFromPy<TVector3*>(item.ptr());
		}
		std::vector<double> additionalVariablesToSave;
		if(not rpwa::py::convertBPObjectToVector<double>(pyAdditionalVariablesToSave, additionalVariablesToSave))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for additionalVariablesToSave when executing rpwa::dataFileWriter::addEvent()");
			bp::throw_error_already_set();
		}
		self.addEvent(initialStateMomenta, finalStateMomenta, additionalVariablesToSave);
	}

}


void rpwa::py::exportEventFileWriter() {

	bp::class_<rpwa::eventFileWriter>("dataFileWriter")
		.def(
			"initialize"
			, &eventFileWriter_initialize
			, (bp::arg("outputFile"),
			   bp::arg("userString"),
			   bp::arg("initialStateParticleNames"),
			   bp::arg("finalStateParticleNames"),
			   bp::arg("binningMap"),
			   bp::arg("additionalVariableLabels"),
			   bp::arg("eventTreeName")="rootPwaEvtTree",
			   bp::arg("initialStateMomentaBranchName")="prodKinMomenta",
			   bp::arg("finalStateMomentaBranchName")="decayKinMomenta",
			   bp::arg("metadataName")="dataMetadata",
			   bp::arg("splitlevel")=99,
			   bp::arg("buffsize")=256000)
		)
		.def(
			"addEvent"
			, &eventFileWriter_addEvent
			, (bp::arg("initialStateMomenta"),
			   bp::arg("finalStateMomenta"),
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
