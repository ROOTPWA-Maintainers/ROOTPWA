
#include "calcAmplitude_py.h"

namespace bp = boost::python;

namespace {

	bp::list calcAmplitude(rpwa::eventMetadata&            eventMeta,
	                       const rpwa::isobarAmplitudePtr& amplitude,
	                       const long int                  maxNmbEvents            = -1,
	                       const bool                      printProgress           = true,
	                       const std::string&              treePerfStatOutFileName = "",         // root file name for tree performance result
	                       const long int                  treeCacheSize           = 25000000)
	{
		return bp::list(rpwa::hli::calcAmplitude(eventMeta,
		                                         amplitude,
		                                         maxNmbEvents,
		                                         printProgress,
		                                         treePerfStatOutFileName,
		                                         treeCacheSize));
	}

}

void rpwa::py::exportCalcAmplitude()
{

	bp::def(
		"calcAmplitude"
		, &::calcAmplitude
		, (bp::arg("eventMeta"),
		   bp::arg("amplitude"),
		   bp::arg("maxNmbEvents")=-1,
		   bp::arg("printProgress")=true,
		   bp::arg("treePerfStatOutFileName")="",
		   bp::arg("treeCacheSize")=25000000)
	);

}
