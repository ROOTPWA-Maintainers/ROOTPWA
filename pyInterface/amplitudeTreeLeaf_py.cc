#include "amplitudeTreeLeaf_py.h"

namespace bp = boost::python;

namespace {

	struct amplitudeTreeLeafWrapper : public rpwa::amplitudeTreeLeaf,
	                                         bp::wrapper<rpwa::amplitudeTreeLeaf>
	{
	
		amplitudeTreeLeafWrapper()
			: rpwa::amplitudeTreeLeaf(),
			  bp::wrapper<rpwa::amplitudeTreeLeaf>() { };

		amplitudeTreeLeafWrapper(const rpwa::amplitudeTreeLeaf& ampTreeLeaf)
			: rpwa::amplitudeTreeLeaf(ampTreeLeaf),
			  bp::wrapper<rpwa::amplitudeTreeLeaf>() { };

		int Write__(std::string name) {
			return this->Write(name.c_str());
		};

		const std::complex<double>& incohSubAmp__1(const unsigned int index = 0) {
			return rpwa::amplitudeTreeLeaf::incohSubAmp(index);
		};

		const std::complex<double>& incohSubAmp__2(const std::string& subAmpLabel) {
			return rpwa::amplitudeTreeLeaf::incohSubAmp(subAmpLabel);
		};

	};

}

void rpwa::py::exportAmplitudeTreeLeaf() {

	bp::class_<amplitudeTreeLeafWrapper>("amplitudeTreeLeaf")

		.def(bp::self == bp::self)
		.def(bp::self += bp::self)
		.def(bp::self -= bp::self)

		.def(bp::self *= int())
		.def(bp::self *= double())
		.def(bp::self /= int())
		.def(bp::self /= double())

		.def(bp::self_ns::str(bp::self))

		.def("clear", &amplitudeTreeLeafWrapper::clear)

		.def("nmbIncohSubAmps", &amplitudeTreeLeafWrapper::nmbIncohSubAmps)

		.def("containsIncohSubAmp", &amplitudeTreeLeafWrapper::containsIncohSubAmp)
		.def("incohSubAmpIndex", &amplitudeTreeLeafWrapper::incohSubAmpIndex)
		.def(
			"incohSubAmpName"
			, &amplitudeTreeLeafWrapper::incohSubAmpName
			, bp::return_value_policy<bp::return_by_value>()
		)

		.def(
			"incohSubAmp"
			, &amplitudeTreeLeafWrapper::incohSubAmp__1
			, bp::return_value_policy<bp::return_by_value>()
		)

		.def(
			"incohSubAmp"
			, &amplitudeTreeLeafWrapper::incohSubAmp__2
			, bp::return_value_policy<bp::return_by_value>()
		)

		.def(
			"amp"
			, &amplitudeTreeLeafWrapper::amp
			, bp::return_value_policy<bp::return_by_value>()
		)

		.def("defineIncohSubAmps", &amplitudeTreeLeafWrapper::defineIncohSubAmps)
		.def("setIncohSubAmp", &amplitudeTreeLeafWrapper::setIncohSubAmp)
		.def("setAmp", &amplitudeTreeLeafWrapper::setAmp)

		.def("Write", &amplitudeTreeLeafWrapper::Write__)

		.add_static_property("debugAmplitudeTreeLeaf", &amplitudeTreeLeafWrapper::debug, &amplitudeTreeLeafWrapper::setDebug);

};

