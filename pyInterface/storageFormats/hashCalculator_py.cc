#include "hashCalculator_py.h"
#include "hashCalculator.h"

#include <boost/python.hpp>
#include <complex>

namespace bp = boost::python;

namespace {
	void hashCalculator_Update1(rpwa::hashCalculator& self,
	                       double                val) {
		self.Update(val);
	}
	void hashCalculator_Update2(rpwa::hashCalculator& self,
	                       std::complex<double>  val) {
		self.Update(val);
	}
}


void rpwa::py::exportHashCalculator() {

	bp::class_<rpwa::hashCalculator>("hashCalculator")
		.def(bp::init<const  rpwa::hashCalculator&>())
		.def("hash", &rpwa::hashCalculator::hash)
		.def("Update", &::hashCalculator_Update1, bp::arg("value"))
		.def("Update", &::hashCalculator_Update2, bp::arg("value"))
	;
};
