#include "amplitudeMetadata_py.h"

#include <boost/python.hpp>

#include <TPython.h>
#include <TTree.h>

#include "amplitudeMetadata.h"
#include "rootConverters_py.h"

namespace bp = boost::python;


namespace {

	bp::list amplitudeMetadata_eventMetadata(const rpwa::amplitudeMetadata& self)
	{
		bp::list retval;
		const std::vector<rpwa::eventMetadata>& metas = self.eventMetadata();
		for(unsigned int i = 0; i < metas.size(); ++i) {
			retval.append(bp::object(boost::cref(metas[i])));
		}
		return retval;
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

	PyObject* amplitudeMetadata_amplitudeTree(rpwa::amplitudeMetadata& self)
	{
		TTree* tree = self.amplitudeTree();
		return TPython::ObjectProxy_FromVoidPtr(tree, tree->ClassName());
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
			, bp::return_value_policy<bp::manage_new_object, bp::with_custodian_and_ward_postcall<0, 1> >()
		)
		.staticmethod("readAmplitudeFile")
		.def("amplitudeTree", &amplitudeMetadata_amplitudeTree)
		.def_readonly("amplitudeLeafName", &rpwa::amplitudeMetadata::amplitudeLeafName)
		;

}
