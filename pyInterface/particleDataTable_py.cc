#include<particleDataTable_py.h>

namespace bp = boost::python;

namespace {

	struct particleDataTableWrapper : public rpwa::particleDataTable {

		static bool read__(bp::object pyLine){
			std::string strLine = bp::extract<std::string>(pyLine);
			std::istringstream sstrLine(strLine, std::istringstream::in);
			return rpwa::particleDataTable::read(sstrLine);
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
				"entry", &particleDataTableWrapper::entry
				, ( bp::arg("partName"), bp::arg("warnIfNotExistent")=(bool const)(true) )
				, bp::return_value_policy<bp::reference_existing_object>()
		)
		.staticmethod( "entry" )    

		.def("addEntry", &particleDataTableWrapper::addEntry)
		.staticmethod("addEntry")    

//		.def("entriesMatching")

		.def("nmbEntries", &particleDataTableWrapper::nmbEntries) 
		.staticmethod("nmbEntries")

		.def("__iter__", bp::iterator< particleDataTableWrapper >())

		.def( 
				"readFile"
				, &particleDataTableWrapper::readFile
				, ( bp::arg("fileName")="./particleDataTable.txt" )
		)
		.staticmethod("readFile")
		.def("read", &particleDataTableWrapper::read__)
		.staticmethod("read")
		.def("readDecayFile", &particleDataTableWrapper::readDecayFile)
		.staticmethod("readDecayFile")

		.def("clear", &particleDataTableWrapper::clear)
		.staticmethod("clear")    

		.add_static_property("debug", &particleDataTableWrapper::debug, &particleDataTableWrapper::setDebug);

}

