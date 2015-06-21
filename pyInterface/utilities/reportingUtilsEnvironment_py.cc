#include "reportingUtilsEnvironment_py.h"

namespace bp = boost::python;

void rpwa::py::exportReportingUtilsEnvironment() {

	bp::def("printGitHash", &rpwa::printGitHash);
	bp::def("printCompilerInfo", &rpwa::printCompilerInfo);
	bp::def("printLibraryInfo", &rpwa::printLibraryInfo);

}
