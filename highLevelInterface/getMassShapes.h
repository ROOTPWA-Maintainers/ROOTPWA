#ifndef GETMASSSHAPES
#define GETMASSSHAPES
#include "isobarDecayTopology.h"
namespace rpwa {
	std::vector<std::complex<double> > getMassShapes(isobarDecayTopologyPtr &topo,
	                                                 const double            mass,
	                                                 const bool              useBarrierFactors);
}
#endif// GETMASSSHAPES
