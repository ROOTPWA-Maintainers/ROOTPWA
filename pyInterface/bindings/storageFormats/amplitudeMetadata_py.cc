#include "amplitudeMetadata_py.h"

#include <vector>
#include <complex>
#include <boost/python.hpp>

#include <TPython.h>
#include <TTree.h>

#include "amplitudeMetadata.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

#if BOOST_VERSION < 106500
#include <boost/python/numeric.hpp>
#else
#include <boost/python/numpy.hpp>
namespace np = boost::python::numpy;
#endif

namespace bp = boost::python;


namespace {

#if BOOST_VERSION < 106500
	typedef boost::python::numeric::array nparray;
#else
	typedef boost::python::numpy::ndarray nparray;
#endif

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


	nparray
	amplitudeMetadata_loadAmplitudes(
	                                 const bp::object& pyAmplitudeMetadata,
	                                 const rpwa::eventMetadata* eventMeta,
	                                 const bp::dict& pyOtfBin,
	                                 const long maxNmbEvents)
	                                 {
		std::vector<const rpwa::amplitudeMetadata*> amplitudeMeta;
		if (not rpwa::py::convertBPObjectToVector<const rpwa::amplitudeMetadata*>(pyAmplitudeMetadata, amplitudeMeta)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for amplitudeMetadata when executing rpwa::ampIntegralMatrix::integrate()");
			bp::throw_error_already_set();
		}
		const rpwa::multibinBoundariesType otfBin = rpwa::py::convertMultibinBoundariesFromPy(pyOtfBin);
		std::vector<std::vector<std::complex<double>>> amps = rpwa::loadAmplitudes(amplitudeMeta, eventMeta, otfBin, maxNmbEvents);
#if BOOST_VERSION < 106500
		bp::list data;
		for (size_t i = 0; i < amps.size(); ++i) {
			bp::list row;
			for (size_t j = 0; j < amps[i].size(); ++j) {
				row.append(amps[i][j]);
			}
			data.append(row);
			amps[i].resize(0);
		}
		nparray pyIntMatrix(data);
#else
		nparray pyIntMatrix = np::empty(bp::make_tuple(amps.size(), amps[0].size()), np::dtype::get_builtin<std::complex<double>>());
		for (unsigned int i = 0; i < amps.size(); ++i) {
			for (unsigned int j = 0; j < amps[i].size(); ++j) {
				pyIntMatrix[bp::make_tuple(i, j)] = amps[i][j];
			}
			amps[i].resize(0);
		}
#endif
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
	       bp::arg("eventMeta")=bp::object(),
	       bp::arg("otfBin")=bp::dict(),
	       bp::arg("maxNmbEvents")=0));

}
