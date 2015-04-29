
#include "pwaFit_py.h"

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	rpwa::fitResultPtr pwaFit(PyObject*          ampTreesDict,
	                          PyObject*          normMatrix,
	                          PyObject*          accMatrix,
	                          const bp::object   pyWaveNames,
	                          const bp::object   pyWaveThresholds,
	                          const double       massBinMin = 0.,
	                          const double       massBinMax = 0.,
	                          const int          seed = 0,
	                          const std::string  startValFileName = "",
	                          const unsigned int accEventsOverride=0,
	                          const unsigned int rank = 1,
			                  const bool         runHesse=false,
	                          const bool         verbose = false)
	{
		std::map<std::string, TTree*> ampTreesMap;
		bp::dict pyDictAmpTreesDict = bp::extract<bp::dict>(ampTreesDict);
		bp::list keys = pyDictAmpTreesDict.keys();
		for(int i = 0; i < bp::len(keys); ++i) {
			bp::object curTree = pyDictAmpTreesDict[keys[i]];
			if(curTree) {
				std::string waveName = bp::extract<std::string>(keys[i]);
				TTree* treePtr = rpwa::py::convertFromPy<TTree*>(curTree.ptr());
				if(not treePtr) {
					PyErr_SetString(PyExc_TypeError, "Got invalid input for ampTreesDict when executing rpwa::pwaFit()");
					bp::throw_error_already_set();
				}
				ampTreesMap.insert(std::pair<std::string, TTree*>(waveName, treePtr));
			}
		}
		rpwa::ampIntegralMatrix* normMatrixPtr = rpwa::py::convertFromPy<rpwa::ampIntegralMatrix*>(normMatrix);
		rpwa::ampIntegralMatrix* accMatrixPtr = rpwa::py::convertFromPy<rpwa::ampIntegralMatrix*>(accMatrix);
		std::vector<std::string> vectorWaveNames;
		std::vector<double> vectorWaveThresholds;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyWaveNames, vectorWaveNames)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for waveNames when executing rpwa::pwaFit()");
			bp::throw_error_already_set();
		}
		if(not rpwa::py::convertBPObjectToVector<double>(pyWaveThresholds, vectorWaveThresholds)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for waveThresholds when executing rpwa::pwaFit()");
			bp::throw_error_already_set();
		}

		return rpwa::hli::pwaFit(ampTreesMap,
		                         *normMatrixPtr,
		                         *accMatrixPtr,
		                         vectorWaveNames,
		                         vectorWaveThresholds,
		                         massBinMin,
		                         massBinMax,
		                         seed,
		                         startValFileName,
		                         accEventsOverride,
		                         rank,
								 runHesse,
		                         verbose);
	}

}

void rpwa::py::exportPwaFit()
{

	bp::def(
		"pwaFit"
		, &::pwaFit
		, (bp::arg("ampTreesDict"),
		   bp::arg("normMatrix"),
		   bp::arg("accMatrix"),
		   bp::arg("waveNames"),
		   bp::arg("waveThresholds"),
		   bp::arg("massBinMin") = 0,
		   bp::arg("massBinMax") = 0,
		   bp::arg("seed") = 0,
		   bp::arg("startValFileName") = "",
		   bp::arg("accEventsOverride") = 0,
		   bp::arg("rank") = 1,
		   bp::arg("runHesse") = false,
		   bp::arg("verbose") = false)
	);

}
