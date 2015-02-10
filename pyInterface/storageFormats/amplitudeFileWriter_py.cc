
#include "amplitudeFileWriter_py.h"

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	bool amplitudeFileWriter_initialize(rpwa::amplitudeFileWriter& self,
	                                    PyObject*                  pyOutputFile,
                                        bp::object                 pyEventMetadata,
                                        const std::string&         keyfileContent,
                                        const std::string&         objectBasename,
                                        const int&                 splitlevel = 99,
                                        const int&                 buffsize = 256000)
	{
		TFile* outputFile = rpwa::py::convertFromPy<TFile*>(pyOutputFile);
		std::vector<rpwa::eventMetadata*> eventMeta;
		if(not rpwa::py::convertBPObjectToVector<rpwa::eventMetadata*>(pyEventMetadata, eventMeta))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for eventMetadata when executing rpwa::amplitudeFileWriter::initialize()");
			bp::throw_error_already_set();
		}
		return self.initialize(*outputFile, eventMeta, keyfileContent, objectBasename, splitlevel, buffsize);
	}

	void amplitudeFileWriter_addAmplitudes(rpwa::amplitudeFileWriter& self,
	                                       bp::list pyAmplitudes)
	{
		std::vector<std::complex<double> > amplitudes;
		if(not rpwa::py::convertBPObjectToVector<std::complex<double> >(pyAmplitudes, amplitudes))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for eventMetadata when executing rpwa::amplitudeFileWriter::addAmplitudes()");
			bp::throw_error_already_set();
		}
		self.addAmplitudes(amplitudes);
	}

}


void rpwa::py::exportAmplitudeFileWriter() {

	bp::class_<rpwa::amplitudeFileWriter, boost::noncopyable>("amplitudeFileWriter")
		.def(
			"initialize"
			, &amplitudeFileWriter_initialize
			, (bp::arg("outputFile"),
			   bp::arg("eventMetadata"),
			   bp::arg("keyfileContent"),
			   bp::arg("objectBaseName"),
			   bp::arg("splitlevel")=99,
			   bp::arg("buffsize")=256000)
		)

		.def("addAmplitude", &rpwa::amplitudeFileWriter::addAmplitude)
		.def("addAmplitudes", &amplitudeFileWriter_addAmplitudes)
		.def("reset", &rpwa::amplitudeFileWriter::reset)
		.def("finalize", &rpwa::amplitudeFileWriter::finalize)
		.def(
			"initialized"
			, &rpwa::amplitudeFileWriter::initialized
			, bp::return_value_policy<bp::copy_const_reference>()
		);

}
