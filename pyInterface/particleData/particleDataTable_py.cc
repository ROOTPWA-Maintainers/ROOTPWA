#include "particleDataTable_py.h"

#include<particleProperties.h>

#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	bool particleDataTable_read(bp::object pyLine) {
		std::string strLine = bp::extract<std::string>(pyLine);
		std::istringstream sstrLine(strLine, std::istringstream::in);
		return rpwa::particleDataTable::read(sstrLine);
	}

	bp::list particleDataTable_entriesMatching(const rpwa::particleProperties& prototype,
	                                           const std::string& sel,
	                                           const double minMass = 0.,
	                                           const double minMassWidthFactor = 0.,
	                                           const bp::object& pyWhiteList = bp::object(),
	                                           const bp::object& pyBlackList = bp::object(),
	                                           const bp::object& pyDecayProducts = bp::object(),
	                                           const bool& forceDecayCheck = true)
	{

		std::vector<std::string> whiteList;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyWhiteList, whiteList)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for whiteList when executing rpwa::particleDataTable::entriesMatching()");
			bp::throw_error_already_set();
		}
		std::vector<std::string> blackList;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyBlackList, blackList)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for blackList when executing rpwa::particleDataTable::entriesMatching()");
			bp::throw_error_already_set();
		}
		std::multiset<std::string> decayProducts;
		if(not rpwa::py::convertBPObjectToMultiSet<std::string>(pyDecayProducts, decayProducts)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for decayProducts when executing rpwa::particleDataTable::entriesMatching()");
			bp::throw_error_already_set();
		}

		std::vector<const rpwa::particleProperties*> retPtrVec = rpwa::particleDataTable::entriesMatching(prototype,
		                                                                                                  sel,
		                                                                                                  minMass,
		                                                                                                  minMassWidthFactor,
		                                                                                                  whiteList,
		                                                                                                  blackList,
		                                                                                                  decayProducts,
		                                                                                                  forceDecayCheck);
		std::vector<rpwa::particleProperties> retVec(retPtrVec.size());
		for(unsigned int i = 0; i < retVec.size(); ++i) {
			retVec[i] = *(retPtrVec[i]);
		}

		return bp::list(retVec);

	}
	bp::tuple particleDataTable_geantIdAndChargeFromParticleName(const std::string& name) {
			int id;
			int charge;
			rpwa::particleDataTable::geantIdAndChargeFromParticleName(name, id, charge);
			return bp::make_tuple(id, charge);
	}

}

void rpwa::py::exportParticleDataTable()
{

	bp::class_< rpwa::particleDataTable, boost::noncopyable >( "particleDataTable", bp::no_init )

		.def(bp::self_ns::str(bp::self))

		.add_static_property(
			"instance"
			, bp::make_function( &rpwa::particleDataTable::instance,  bp::return_value_policy<bp::reference_existing_object>() )
		)

		.def("isInTable",&rpwa::particleDataTable::isInTable)
		.staticmethod("isInTable")

		.def(
			"entry"
			, &rpwa::particleDataTable::entry
			, (bp::arg("partName"), bp::arg("warnIfNotExistent")=(bool const)(true))
			, bp::return_internal_reference<>()
		)
		.staticmethod( "entry" )

		.def("addEntry", &rpwa::particleDataTable::addEntry)
		.staticmethod("addEntry")

		.def(
			"entriesMatching"
			, &particleDataTable_entriesMatching
			, (bp::arg("prototype"),
			   bp::arg("sel"),
			   bp::arg("minMass")=0.,
			   bp::arg("minMassWidthFactor")=0.,
			   bp::arg("whiteList")=bp::list(),
			   bp::arg("blackList")=bp::list(),
			   bp::arg("decayProducts")=bp::list(),
			   bp::arg("forceDecayCheck")=true)
		)
		.staticmethod("entriesMatching")

		.def("nmbEntries", &rpwa::particleDataTable::nmbEntries)
		.staticmethod("nmbEntries")

		.def("__iter__", bp::iterator< rpwa::particleDataTable >())

		.def(
			"readFile"
			, &rpwa::particleDataTable::readFile
			, (bp::arg("fileName")="./particleDataTable.txt")
		)
		.staticmethod("readFile")
		.def("read", &particleDataTable_read)
		.staticmethod("read")
		.def("readDecayModeFile", &rpwa::particleDataTable::readDecayModeFile)
		.staticmethod("readDecayModeFile")

		.def("clear", &rpwa::particleDataTable::clear)
		.staticmethod("clear")

		.def("particleNameFromGeantId", &rpwa::particleDataTable::particleNameFromGeantId)
		.staticmethod("particleNameFromGeantId")

		.def("geantIdFromParticleName", &rpwa::particleDataTable::geantIdFromParticleName)
		.staticmethod("geantIdFromParticleName")

		.def("geantIdAndChargeFromParticleName"
		     , &particleDataTable_geantIdAndChargeFromParticleName
		     , (bp::arg("name")))
		.staticmethod("geantIdAndChargeFromParticleName")

		.add_static_property("debugParticleDataTable", &rpwa::particleDataTable::debug, &rpwa::particleDataTable::setDebug);

}
