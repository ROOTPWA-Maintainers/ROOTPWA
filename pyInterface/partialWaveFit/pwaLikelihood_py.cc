#include "pwaLikelihood_py.h"

#include <complex>

#include "rootConverters_py.h"
#include "stlContainers_py.h"


namespace bp = boost::python;

namespace {

	void
	pwaLikelihood_init(rpwa::pwaLikelihood<std::complex<double> >& self,
	     const unsigned int                          rank,
	     PyObject*                                   pyAmpFileList,
	     const double                                massBinCenter,
	     const std::string&                          waveListFileName,
	     const std::string&                          normIntFileName,
	     const std::string&                          accIntFileName,
	     const unsigned int                          numbAccEvents)
	{
		std::map<std::string, std::string> ampFilesMap;
		bp::dict pyDictAmpFilesDict = bp::extract<bp::dict>(pyAmpFileList);
		bp::list keys = pyDictAmpFilesDict.keys();
		for(int i = 0; i < bp::len(keys); ++i) {
			std::string waveName = bp::extract<std::string>(keys[i]);
			std::string fileName = bp::extract<std::string>(pyDictAmpFilesDict[keys[i]]);
			ampFilesMap.insert(std::pair<std::string, std::string>(waveName, fileName));
		}
		self.init(rank, ampFilesMap, massBinCenter, waveListFileName, normIntFileName, accIntFileName, numbAccEvents);
	}

	void
	pwaLikelihood_setQuiet(rpwa::pwaLikelihood<std::complex<double> >& self, const bool flag = true)
	{
		self.setQuiet(flag);
	}
}

void rpwa::py::exportPwaLikelihood() {

	bp::class_<rpwa::pwaLikelihood<std::complex<double> > >("pwaLikelihood")
		.def(
			"init"
			, ::pwaLikelihood_init
			, (bp::arg("rank"),
			   bp::arg("ampFileList"),
			   bp::arg("massBinCenter"),
			   bp::arg("waveListFileName"),
			   bp::arg("normIntFileName"),
			   bp::arg("accIntFileName"),
			   bp::arg("numbAccEvents") = 0)
			)
		.def("setQuiet"
		     ,::pwaLikelihood_setQuiet)
		     , (bp::arg("flag")=true);
}
