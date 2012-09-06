#include "utils_py.h"

namespace bp = boost::python;

void rpwa::py::exportUtils() {

	bp::def("printSvnVersion", &rpwa::printSvnVersion);
	bp::def("printCompilerInfo", &rpwa::printCompilerInfo);
	bp::def("printLibraryInfo", &rpwa::printLibraryInfo);

};
