#include "amplitudeMetadata_py.h"

#include <vector>
#include <complex>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <TPython.h>
#include <TTree.h>

#include "amplitudeMetadata.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"


namespace bp = boost::python;
namespace np = boost::python::numpy;


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


	np::ndarray
	amplitudeMetadata_loadAmplitudes(
	                                 bp::list& pyAmplitudesFilenames,
	                                 bp::list& pyWaveNames,
	                                 const std::string& eventFilename,
	                                 const bp::dict& pyOtfBin,
	                                 const long maxNmbEvents)
	                                 {
		const rpwa::multibinBoundariesType otfBin = rpwa::py::convertMultibinBoundariesFromPy(pyOtfBin);
		std::vector<std::string> amplitudeFilenames;
		if (not rpwa::py::convertBPObjectToVector<std::string>(pyAmplitudesFilenames, amplitudeFilenames)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for amplitudeFilenames when executing rpwa::amplitudeMetadata::loadAmplitudes()");
			bp::throw_error_already_set();
		}
		std::vector<std::string> waveNames;
		if (not rpwa::py::convertBPObjectToVector<std::string>(pyWaveNames, waveNames)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for waveNames when executing rpwa::amplitudeMetadata::loadAmplitudes()");
			bp::throw_error_already_set();
		}
		std::vector<std::vector<std::complex<double>>> amps = rpwa::loadAmplitudes(amplitudeFilenames, waveNames, eventFilename, otfBin, maxNmbEvents);
		np::ndarray pyIntMatrix = np::empty(bp::make_tuple(amps.size(), amps[0].size()), np::dtype::get_builtin<std::complex<double>>());
		for (unsigned int i = 0; i < amps.size(); ++i) {
			for (unsigned int j = 0; j < amps[i].size(); ++j) {
				pyIntMatrix[bp::make_tuple(i, j)] = amps[i][j];
			}
			amps[i].resize(0);
		}
		return pyIntMatrix;
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
	bp::def(
	       "loadAmplitudes",
	       &amplitudeMetadata_loadAmplitudes,
	      (bp::arg("amplitudeMetadata"),
	       bp::arg("waveNames"),
	       bp::arg("eventMeta"),
	       bp::arg("otfBin")=bp::dict(),
	       bp::arg("maxNmbEvents")=0));

}
