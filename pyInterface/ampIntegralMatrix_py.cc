#include "ampIntegralMatrix_py.h"

namespace bp = boost::python;

namespace {

	struct ampIntegralMatrixWrapper : public rpwa::ampIntegralMatrix,
	                                         bp::wrapper<rpwa::ampIntegralMatrix>
	{

		ampIntegralMatrixWrapper()
			: rpwa::ampIntegralMatrix(),
			  bp::wrapper<rpwa::ampIntegralMatrix>() { };

		ampIntegralMatrixWrapper(const rpwa::ampIntegralMatrix& integral)
			: rpwa::ampIntegralMatrix(integral),
			  bp::wrapper<rpwa::ampIntegralMatrix>() { };
/*
		const bp::list waveDescriptions__() const {
			return bp::list(rpwa::ampIntegralMatrix());
		};
*/
		const rpwa::waveDescription* waveDesc__1(const unsigned int waveIndex) {
			return rpwa::ampIntegralMatrix::waveDesc(waveIndex);
		};

		const rpwa::waveDescription* waveDesc__2(const std::string& waveName) {
			return rpwa::ampIntegralMatrix::waveDesc(waveName);
		};

		std::complex<double> element__1(const unsigned int waveIndexI,
		                                const unsigned int waveIndexJ) const
		{
			return rpwa::ampIntegralMatrix::element(waveIndexI, waveIndexJ);
		};

		std::complex<double> element__2(const std::string& waveNameI,
		                                const std::string& waveNameJ) const
		{
			return rpwa::ampIntegralMatrix::element(waveNameI, waveNameJ);
		};

		bool writeAscii__(const std::string& outFileName) const {
			return rpwa::ampIntegralMatrix::writeAscii(outFileName);
		};

		bool readAscii__(const std::string& inFileName) {
			return rpwa::ampIntegralMatrix::readAscii(inFileName);
		};

		int Write__(std::string name) {
			return this->Write(name.c_str());
		};

	};

}

void rpwa::py::exportAmpIntegralMatrix() {

	bp::class_<ampIntegralMatrixWrapper>("ampIntegralMatrix")

		.def(bp::init<rpwa::ampIntegralMatrix&>())

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

		.def("clear", &ampIntegralMatrixWrapper::clear)
		.def("nmbWaves", &ampIntegralMatrixWrapper::nmbWaves)
		.def("nmbEvents", &ampIntegralMatrixWrapper::nmbEvents)
		.def("setNmbEvents", &ampIntegralMatrixWrapper::setNmbEvents)
		.def("containsWave", &ampIntegralMatrixWrapper::containsWave)
		.def("waveIndex", &ampIntegralMatrixWrapper::waveIndex)
		.def(
			"waveName"
			, &ampIntegralMatrixWrapper::waveName
			, bp::return_value_policy<bp::copy_const_reference>()
		)
//		Disabled because of missing == operator in rpwa::waveDescription
//		See also http://stackoverflow.com/questions/10680691/why-do-i-need-comparison-operators-in-boost-python-vector-indexing-suite
//		.def("waveDescriptions", &ampIntegralMatrixWrapper::waveDescriptions__)
		.def(
			"waveDesc"
			, &ampIntegralMatrixWrapper::waveDesc__1
			, bp::return_value_policy<bp::return_by_value>()
		)
		.def(
			"waveDesc"
			, &ampIntegralMatrixWrapper::waveDesc__2
			, bp::return_value_policy<bp::return_by_value>()
		)
		.def("allWavesHaveDesc", &ampIntegralMatrixWrapper::allWavesHaveDesc)

//		Commenting this until it is decided how the boost::multi_array should be handled in python
//		.def("matrix", &ampIntegralMatrixWrapper::matrix)

		.def("element", &ampIntegralMatrixWrapper::element__1)
		.def("element", &ampIntegralMatrixWrapper::element__2)
		.def("renormalize", &ampIntegralMatrixWrapper::renormalize)
		.def("writeAscii", &ampIntegralMatrixWrapper::writeAscii__)
		.def("readAscii", &ampIntegralMatrixWrapper::readAscii__)

		.def("Write", &ampIntegralMatrixWrapper::Write__)

		.add_static_property("debugAmpIntegralMatrix", &ampIntegralMatrixWrapper::debug, &ampIntegralMatrixWrapper::setDebug);

};
