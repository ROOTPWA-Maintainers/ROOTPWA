
#include "eventMetadata_py.h"

#include <TPython.h>
#include <TTree.h>

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	bp::dict eventMetadata_binningMap(const rpwa::eventMetadata& self)
	{
		bp::dict retval;
		std::map<std::string, std::pair<double, double> > binningMap = self.binningMap();
		for(std::map<std::string, std::pair<double, double> >::const_iterator it = binningMap.begin(); it != binningMap.end(); ++it)
		{
			bp::tuple val = bp::make_tuple<double, double>(it->second.first, it->second.second);
			retval[it->first] = val;
		}
		return retval;
	}

	bp::list eventMetadata_productionKinematicsParticleNames(const rpwa::eventMetadata& self)
	{
		return bp::list(self.productionKinematicsParticleNames());
	}

	bp::list eventMetadata_decayKinematicsParticleNames(const rpwa::eventMetadata& self)
	{
		return bp::list(self.decayKinematicsParticleNames());
	}

	bp::list eventMetadata_additionalSavedVariableLables(const rpwa::eventMetadata& self)
	{
		return bp::list(self.additionalSavedVariableLables());
	}

	PyObject* eventMetadata_eventTree(rpwa::eventMetadata& self)
	{
		TTree* tree = self.eventTree();
		return TPython::ObjectProxy_FromVoidPtr(tree, tree->ClassName());
	}

	const rpwa::eventMetadata* eventMetadata_readEventFile(PyObject* pyInputFile, const bool& quiet = false)
	{
		TFile* inputFile = rpwa::py::convertFromPy<TFile*>(pyInputFile);
		if(not inputFile) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for inputFile when executing rpwa::eventMetadata::readEventFile()");
			bp::throw_error_already_set();
		}
		return rpwa::eventMetadata::readEventFile(inputFile, quiet);
	}

}

void rpwa::py::exportEventMetadata() {

	bp::scope theScope = bp::class_<rpwa::eventMetadata, boost::noncopyable>("eventMetadata", bp::no_init)

		.def(bp::self_ns::str(bp::self))

		.def(
			"userString"
			, &rpwa::eventMetadata::userString
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"contentHash"
			, &rpwa::eventMetadata::contentHash
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("binningMap", &eventMetadata_binningMap)
		.def("productionKinematicsParticleNames", &eventMetadata_productionKinematicsParticleNames)
		.def("decayKinematicsParticleNames", &eventMetadata_decayKinematicsParticleNames)
		.def("additionalSavedVariableLables", &eventMetadata_additionalSavedVariableLables)
		.def("recalculateHash", &rpwa::eventMetadata::recalculateHash, (bp::arg("printProgress")=false))
		.def("eventTree", &eventMetadata_eventTree)
		.def(
			"readEventFile"
			, &eventMetadata_readEventFile
			, (bp::arg("inputFile"), bp::arg("quiet")=false)
			, bp::return_value_policy<bp::reference_existing_object>()
		)
		.staticmethod("readEventFile")
		.def_readonly("objectNameInFile", &rpwa::eventMetadata::objectNameInFile)
		.def_readonly("eventTreeName", &rpwa::eventMetadata::eventTreeName)
		.def_readonly("productionKinematicsMomentaBranchName", &rpwa::eventMetadata::productionKinematicsMomentaBranchName)
		.def_readonly("decayKinematicsMomentaBranchName", &rpwa::eventMetadata::decayKinematicsMomentaBranchName)

		;


	bp::enum_<rpwa::eventMetadata::eventsTypeEnum>("eventsTypeEnum")
		.value("OTHER", rpwa::eventMetadata::OTHER)
		.value("REAL", rpwa::eventMetadata::REAL)
		.value("GENERATED", rpwa::eventMetadata::GENERATED)
		.value("ACCEPTED", rpwa::eventMetadata::ACCEPTED)
		.export_values();

	theScope.attr("OTHER") = rpwa::eventMetadata::OTHER;
	theScope.attr("REAL") = rpwa::eventMetadata::REAL;
	theScope.attr("GENERATED") = rpwa::eventMetadata::GENERATED;
	theScope.attr("ACCEPTED") = rpwa::eventMetadata::ACCEPTED;

}
