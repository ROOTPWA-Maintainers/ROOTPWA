
#include "randomNumberGenerator_py.h"

#include<TRandom3.h>

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	PyObject* randomNumberGenerator_getGenerator(rpwa::randomNumberGenerator& self) {
		return rpwa::py::convertToPy<TRandom3>(*(self.getGenerator()));
	}

}

void rpwa::py::exportRandomNumberGenerator() {

	bp::class_<rpwa::randomNumberGenerator, boost::noncopyable>("randomNumberGenerator", bp::no_init)

		.add_static_property(
			"instance"
			, bp::make_function( &rpwa::randomNumberGenerator::instance,  bp::return_value_policy<bp::reference_existing_object>() )
		)

		.def("getGenerator", &randomNumberGenerator_getGenerator)

		.def("seed", &rpwa::randomNumberGenerator::seed)
		.def("setSeed", &rpwa::randomNumberGenerator::setSeed)
		.def("rndm", &rpwa::randomNumberGenerator::rndm);

}
