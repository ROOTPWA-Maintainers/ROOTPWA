#include "amplitudeTreeLeaf_py.h"

namespace bp = boost::python;

namespace {

	int amplitudeTreeLeaf_Write(const rpwa::amplitudeTreeLeaf& self, std::string name) {
		return self.Write(name.c_str());
	};

	const std::complex<double>& amplitudeTreeLeaf_incohSubAmp1(const rpwa::amplitudeTreeLeaf& self, const unsigned int index = 0) {
		return self.incohSubAmp(index);
	};

	const std::complex<double>& amplitudeTreeLeaf_incohSubAmp2(const rpwa::amplitudeTreeLeaf& self, const std::string& subAmpLabel) {
		return self.incohSubAmp(subAmpLabel);
	};

}

void rpwa::py::exportAmplitudeTreeLeaf() {

	bp::class_<rpwa::amplitudeTreeLeaf>("amplitudeTreeLeaf")

		.def(bp::self == bp::self)
		.def(bp::self != bp::self)

		.def(bp::self + bp::self)
		.def(bp::self - bp::self)
		.def(bp::self += bp::self)
		.def(bp::self -= bp::self)

		.def(bp::self * int())
		.def(bp::self * double())
		.def(int() * bp::self)
		.def(double() * bp::self)
		.def(bp::self / double())
		.def(bp::self / int())

		.def(bp::self *= int())
		.def(bp::self *= double())
		.def(bp::self /= int())
		.def(bp::self /= double())

		.def(bp::self_ns::str(bp::self))

		.def("clear", &rpwa::amplitudeTreeLeaf::clear)

		.def("nmbIncohSubAmps", &rpwa::amplitudeTreeLeaf::nmbIncohSubAmps)

		.def("containsIncohSubAmp", &rpwa::amplitudeTreeLeaf::containsIncohSubAmp)
		.def("incohSubAmpIndex", &rpwa::amplitudeTreeLeaf::incohSubAmpIndex)
		.def(
			"incohSubAmpName"
			, &rpwa::amplitudeTreeLeaf::incohSubAmpName
			, bp::return_value_policy<bp::return_by_value>()
		)

		.def(
			"incohSubAmp"
			, &amplitudeTreeLeaf_incohSubAmp1
			, bp::return_value_policy<bp::return_by_value>()
		)

		.def(
			"incohSubAmp"
			, &amplitudeTreeLeaf_incohSubAmp2
			, bp::return_value_policy<bp::return_by_value>()
		)

		.def(
			"amp"
			, &rpwa::amplitudeTreeLeaf::amp
			, bp::return_value_policy<bp::return_by_value>()
		)

		.def("defineIncohSubAmps", &rpwa::amplitudeTreeLeaf::defineIncohSubAmps)
		.def("setIncohSubAmp", &rpwa::amplitudeTreeLeaf::setIncohSubAmp)
		.def("setAmp", &rpwa::amplitudeTreeLeaf::setAmp)

		.def("Write", &amplitudeTreeLeaf_Write)

		.add_static_property("debugAmplitudeTreeLeaf", &rpwa::amplitudeTreeLeaf::debug, &rpwa::amplitudeTreeLeaf::setDebug);

};

