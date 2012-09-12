#include "pythonAdministrator_py.h"

namespace bp = boost::python;

bool rpwa::py::pythonAdministrator::constructAmplitude(std::string keyFileName)
{
	_waveDescription.parseKeyFile(keyFileName);
	return _waveDescription.constructAmplitude(_amplitude);
};

bool rpwa::py::pythonAdministrator::initKinematicsData(PyObject* pyProdKinParticles,
                                                       PyObject* pyDecayKinParticles)
{
	TClonesArray* prodKinParticles = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinParticles);
	TClonesArray* decayKinParticles = rpwa::py::convertFromPy<TClonesArray*>(pyDecayKinParticles);
	if((prodKinParticles == NULL) || (decayKinParticles == NULL)) {
		printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::initKinematicsData()."<<std::endl;
		return false;
	}
	const isobarDecayTopologyPtr& topo = _amplitude->decayTopology();
	return topo->initKinematicsData(*prodKinParticles, *decayKinParticles);
};

bool rpwa::py::pythonAdministrator::readKinematicsData(PyObject* pyProdKinMomenta,
                                                       PyObject* pyDecayKinMomenta)
{
	TClonesArray* prodKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyProdKinMomenta);
	TClonesArray* decayKinMomenta = rpwa::py::convertFromPy<TClonesArray*>(pyDecayKinMomenta);
	if((prodKinMomenta == NULL) || (decayKinMomenta == NULL)) {
		printErr<<"Got invalid input when executing rpwa::diffractiveDissVertex::readKinematicsData()."<<std::endl;
		return false;
	}
	const isobarDecayTopologyPtr& topo = _amplitude->decayTopology();
	return topo->readKinematicsData(*prodKinMomenta, *decayKinMomenta);
};

void rpwa::py::exportPythonAdministrator() {

	bp::class_<rpwa::py::pythonAdministrator>("pythonAdministrator")
		.def("__call__", &rpwa::py::pythonAdministrator::operator())
		.def("constructAmplitude", &rpwa::py::pythonAdministrator::constructAmplitude)
		.def("initKinematicsData", &rpwa::py::pythonAdministrator::initKinematicsData)
		.def("readKinematicsData", &rpwa::py::pythonAdministrator::readKinematicsData)
		.def("reset", &rpwa::py::pythonAdministrator::reset)
		.def("setBranchAddress", &rpwa::py::pythonAdministrator::setBranchAddress)
		.staticmethod("setBranchAddress")
		.def(
			"branch"
			, &rpwa::py::pythonAdministrator::branch
			, (bp::arg("pyTree"),
			   bp::arg("pyAmplitudeTreeLeaf"),
			   bp::arg("name"),
			   bp::arg("bufsize")=32000,
			   bp::arg("splitlevel")=99)
		)
		.staticmethod("branch");

};
