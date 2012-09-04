#include "particleDataTable_py.h"

#include<particleProperties.h>

namespace bp = boost::python;

namespace {

	struct particleDataTableWrapper : public rpwa::particleDataTable,
	                                         bp::wrapper<rpwa::particleDataTable> {

		static bool read__(bp::object pyLine){
			std::string strLine = bp::extract<std::string>(pyLine);
			std::istringstream sstrLine(strLine, std::istringstream::in);
			return rpwa::particleDataTable::read(sstrLine);
		}

		static bp::list entriesMatching__(const rpwa::particleProperties& prototype,
		                                  const std::string& sel,
		                                  const double minMass = 0.,
		                                  const double minMassWidthFactor = 0.,
		                                  const bp::object& pyWhiteList = bp::object(),
		                                  const bp::object& pyBlackList = bp::object(),
		                                  const bp::object& pyDecayProducts = bp::object(),
		                                  const bool& forceDecayCheck = true)
		{
			bp::list pyListWhiteList = bp::extract<bp::list>(pyWhiteList);
			bp::list pyListBlackList = bp::extract<bp::list>(pyBlackList);
			bp::list pyListDecayProducts = bp::extract<bp::list>(pyDecayProducts);
			std::vector<std::string> whiteList(bp::len(pyListWhiteList), "");
			std::vector<std::string> blackList(bp::len(pyListBlackList), "");
			std::set<std::string> decayProducts;

			for(unsigned int i = 0; i < bp::len(pyListWhiteList); ++i) {
				whiteList[i] = bp::extract<std::string>(pyListWhiteList[i]);
			}
			for(unsigned int i = 0; i < bp::len(pyListBlackList); ++i) {
				blackList[i] = bp::extract<std::string>(pyListBlackList[i]);
			}
			for(unsigned int i = 0; i < bp::len(pyListDecayProducts); ++i) {
				decayProducts.insert(bp::extract<std::string>(pyListDecayProducts[i]));
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

			bp::object iter = bp::iterator<std::vector<rpwa::particleProperties> >()(retVec);
			return bp::list(iter);
		}

	};

}

void rpwa::py::exportParticleDataTable()
{

	bp::class_< particleDataTableWrapper, boost::noncopyable >( "particleDataTable", bp::no_init )    

		.def( bp::self_ns::str( bp::self ) )

		.add_static_property(
			"instance"
			, bp::make_function( &particleDataTableWrapper::instance,  bp::return_value_policy<bp::reference_existing_object>() )
		)

		.def("isInTable",&particleDataTableWrapper::isInTable)
		.staticmethod("isInTable") 

		.def( 
			"entry"
			, &particleDataTableWrapper::entry
			, (bp::arg("partName"), bp::arg("warnIfNotExistent")=(bool const)(true))
			, bp::return_value_policy<bp::reference_existing_object>()
		)
		.staticmethod( "entry" )    

		.def("addEntry", &particleDataTableWrapper::addEntry)
		.staticmethod("addEntry")    

		.def(
			"entriesMatching"
			, &particleDataTableWrapper::entriesMatching__
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

		.def("nmbEntries", &particleDataTableWrapper::nmbEntries) 
		.staticmethod("nmbEntries")

		.def("__iter__", bp::iterator< particleDataTableWrapper >())

		.def( 
			"readFile"
			, &particleDataTableWrapper::readFile
			, (bp::arg("fileName")="./particleDataTable.txt")
		)
		.staticmethod("readFile")
		.def("read", &particleDataTableWrapper::read__)
		.staticmethod("read")
		.def("readDecayFile", &particleDataTableWrapper::readDecayFile)
		.staticmethod("readDecayFile")

		.def("clear", &particleDataTableWrapper::clear)
		.staticmethod("clear")    

		.add_static_property("debugParticleDataTable", &particleDataTableWrapper::debug, &particleDataTableWrapper::setDebug);

}

