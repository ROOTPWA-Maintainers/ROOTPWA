
#include "amplitudeMetadata_py.h"

#include "rootConverters_py.h"

namespace bp = boost::python;


namespace {

	bp::list amplitudeMetadata_eventMetadata(const rpwa::amplitudeMetadata& self)
	{
		return bp::list(self.eventMetadata());
	}

	const rpwa::amplitudeMetadata* amplitudeMetadata_readAmplitudeFile(PyObject* pyInputFile,
	                                                                   const std::string& objectBaseName,
	                                                                   const bool& quiet = false)
	{
		TFile* inputFile = rpwa::py::convertFromPy<TFile*>(pyInputFile);
		if(not inputFile) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for inputFile when executing rpwa::amplitudeMetadata::readAmplitudeFile()");
			bp::throw_error_already_set();
		}
		return rpwa::amplitudeMetadata::readAmplitudeFile(inputFile, objectBaseName, quiet);
	}

}


void rpwa::py::exportAmplitudeMetadata() {

	bp::class_<rpwa::amplitudeMetadata, boost::noncopyable>("amplitudeMetadata", bp::no_init)
		.def(bp::self_ns::str(bp::self))
		.def(
			"contentHash"
			, &rpwa::amplitudeMetadata::contentHash
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("eventMetadata", &amplitudeMetadata_eventMetadata)
		.def(
			"keyfileContent"
			, &rpwa::amplitudeMetadata::keyfileContent
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"rootpwaGitHash"
			, &rpwa::amplitudeMetadata::rootpwaGitHash
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"objectBaseName"
			, &rpwa::amplitudeMetadata::objectBaseName
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"recalculateHash"
			, &rpwa::amplitudeMetadata::recalculateHash
			, (bp::arg("printProgress")=false)
		)
		.def(
			"readAmplitudeFile"
			, &amplitudeMetadata_readAmplitudeFile
			, (bp::arg("inputFile"), bp::arg("objectBaseName"), bp::arg("quiet")=false)
			, bp::return_value_policy<bp::reference_existing_object>()
		)
		.staticmethod("readAmplitudeFile")
		.def_readonly("amplitudeLeafName", &rpwa::amplitudeMetadata::amplitudeLeafName)
		;

}
