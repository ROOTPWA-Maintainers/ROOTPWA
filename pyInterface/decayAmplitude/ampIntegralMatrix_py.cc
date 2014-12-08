#include "ampIntegralMatrix_py.h"

#include <TDirectory.h>

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

/*
	const bp::list ampIntegralMatrix_waveDescriptions(const rpwa::ampIntegralMatrix& self) {
		return bp::list(self());
	}
*/
	const rpwa::waveDescription& ampIntegralMatrix_waveDesc1(const rpwa::ampIntegralMatrix& self,
	                                                         const unsigned int waveIndex)
	{
		return *(self.waveDesc(waveIndex));
	}

	const rpwa::waveDescription& ampIntegralMatrix_waveDesc2(const rpwa::ampIntegralMatrix& self,
	                                                         const std::string& waveName) {
		return *(self.waveDesc(waveName));
	}

	std::complex<double> ampIntegralMatrix_element1(const rpwa::ampIntegralMatrix& self,
	                                                const unsigned int waveIndexI,
	                                                const unsigned int waveIndexJ)
	{
		return self.element(waveIndexI, waveIndexJ);
	}

	std::complex<double> ampIntegralMatrix_element2(const rpwa::ampIntegralMatrix& self,
	                                                const std::string& waveNameI,
	                                                const std::string& waveNameJ)
	{
		return self.element(waveNameI, waveNameJ);
	}

	bool ampIntegralMatrix_integrate(rpwa::ampIntegralMatrix& self,
	                                 const bp::list&  ampTreeList,
	                                 const bp::object& waveNameList,
	                                 const unsigned long maxNmbEvents,
	                                 const std::string& weightFileName)
	{
		std::vector<TTree*> ampTreeVector(len(ampTreeList));
		for(unsigned int i = 0; i < len(ampTreeList); ++i) {
			bp::object item = bp::extract<bp::object>(ampTreeList[i]);
			ampTreeVector[i] = rpwa::py::convertFromPy<TTree*>(item.ptr());
		}
		std::vector<std::string> waveNamesVector;
		if(not rpwa::py::convertBPObjectToVector<std::string>(waveNameList, waveNamesVector))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input for waveNames when executing rpwa::ampIntegralMatrix::integrate()");
			bp::throw_error_already_set();
		}
		return self.integrate(ampTreeVector, waveNamesVector, maxNmbEvents, weightFileName);
	}

	bool ampIntegralMatrix_writeAscii(const rpwa::ampIntegralMatrix& self, const std::string& outFileName) {
		return self.writeAscii(outFileName);
	}

	bool ampIntegralMatrix_readAscii(rpwa::ampIntegralMatrix& self, const std::string& inFileName) {
		return self.readAscii(inFileName);
	}

	int ampIntegralMatrix_Write(const rpwa::ampIntegralMatrix& self, const char* name = 0) {
		return self.Write(name);
	}

}

void rpwa::py::exportAmpIntegralMatrix() {

	bp::class_<rpwa::ampIntegralMatrix>("ampIntegralMatrix")

		.def(bp::init<const rpwa::ampIntegralMatrix&>())

		.def(bp::self_ns::str(bp::self))
		.def(bp::self == bp::self)
		.def(bp::self != bp::self)
		.def(bp::self + bp::self)
		.def(bp::self - bp::self)
		.def(bp::self * double())
		.def(bp::self / double())
		.def(double() * bp::self)
		.def(bp::self += bp::self)
		.def(bp::self -= bp::self)
		.def(bp::self *= double())
		.def(bp::self /= double())

		.def("clear", &rpwa::ampIntegralMatrix::clear)
		.def("nmbWaves", &rpwa::ampIntegralMatrix::nmbWaves)
		.def("nmbEvents", &rpwa::ampIntegralMatrix::nmbEvents)
		.def("setNmbEvents", &rpwa::ampIntegralMatrix::setNmbEvents)
		.def("containsWave", &rpwa::ampIntegralMatrix::containsWave)
		.def("waveIndex", &rpwa::ampIntegralMatrix::waveIndex)
		.def(
			"waveName"
			, &rpwa::ampIntegralMatrix::waveName
			, bp::return_value_policy<bp::return_by_value>()
		)
//		Disabled because of missing == operator in rpwa::waveDescription
//		See also http://stackoverflow.com/questions/10680691/why-do-i-need-comparison-operators-in-boost-python-vector-indexing-suite
//		.def("waveDescriptions", &ampIntegralMatrix_waveDescriptions)
		.def(
			"waveDesc"
			, &ampIntegralMatrix_waveDesc1
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"waveDesc"
			, &ampIntegralMatrix_waveDesc2
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("allWavesHaveDesc", &rpwa::ampIntegralMatrix::allWavesHaveDesc)

//		Commenting this until it is decided how the boost::multi_array should be handled in python
//		.def("matrix", &rpwa::ampIntegralMatrix::matrix)

		.def("element", &ampIntegralMatrix_element1)
		.def("element", &ampIntegralMatrix_element2)

		.def("integrate"
		     , &ampIntegralMatrix_integrate
		     , (bp::arg("ampTreeList"),
		         bp::arg("waveNameList"),
		         bp::arg("maxNmbEvents")=0,
		         bp::arg("weightFileName")="")
		)

		.def("renormalize", &rpwa::ampIntegralMatrix::renormalize)
		.def("writeAscii", &ampIntegralMatrix_writeAscii)
		.def("readAscii", &ampIntegralMatrix_readAscii)

		.def("Write", &ampIntegralMatrix_Write, bp::arg("name")=0)
		.def("setBranchAddress", &rpwa::py::setBranchAddress<rpwa::ampIntegralMatrix*>)

		.add_static_property("debugAmpIntegralMatrix", &rpwa::ampIntegralMatrix::debug, &rpwa::ampIntegralMatrix::setDebug);

	bp::def(
		"getFromTDirectory"
		, &rpwa::py::getFromTDirectory<rpwa::ampIntegralMatrix>
		, bp::return_value_policy<bp::manage_new_object>()
	);

}
