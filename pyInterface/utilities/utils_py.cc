#include "utils_py.h"

namespace bp = boost::python;

void rpwa::py::exportUtils() {

	bp::def("printGitHash", &rpwa::printGitHash);
	bp::def("printCompilerInfo", &rpwa::printCompilerInfo);
	bp::def("printLibraryInfo", &rpwa::printLibraryInfo);

}
