#include "reportingUtilsEnvironment_py.h"

#include <boost/python.hpp>

#include "reportingUtilsEnvironment.h"

namespace bp = boost::python;


void rpwa::py::exportReportingUtilsEnvironment() {

	bp::def("printGitHash", &rpwa::printGitHash);
	bp::def("printCompilerInfo", &rpwa::printCompilerInfo);
	bp::def("printLibraryInfo", &rpwa::printLibraryInfo);

	bp::def("gitHash", &rpwa::gitHash);

}
