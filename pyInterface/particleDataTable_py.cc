#include<particleDataTable_py.h>

namespace bp = boost::python;

void exportParticleDataTable()
{

	bp::class_< rpwa::particleDataTable, boost::noncopyable >( "particleDataTable", bp::no_init )    

		.def( bp::self_ns::str( bp::self ) )

		.add_static_property(
				"instance"
				, bp::make_function( &rpwa::particleDataTable::instance,  bp::return_value_policy<bp::reference_existing_object>() )
		)

		.def("isInTable",&rpwa::particleDataTable::isInTable)
		.staticmethod("isInTable") 

		.def( 
				"entry", &::rpwa::particleDataTable::entry
				, ( bp::arg("partName"), bp::arg("warnIfNotExistent")=(bool const)(true) )
				, bp::return_value_policy<bp::reference_existing_object>()
		)
		.staticmethod( "entry" )    

		.def("addEntry", &rpwa::particleDataTable::addEntry)
		.staticmethod("addEntry")    
		.def("nmbEntries", &rpwa::particleDataTable::nmbEntries) 
		.staticmethod("nmbEntries")

		.def("__iter__", bp::iterator< rpwa::particleDataTable >())

		.def( 
				"readFile"
				, &rpwa::particleDataTable::readFile
				, ( bp::arg("fileName")="./particleDataTable.txt" )
		)
		.staticmethod("readFile")
/*		.def("read", &rpwa::particleDataTable::read)
		.staticmethod("read")
*/		.def("readDecayFile", &rpwa::particleDataTable::readDecayFile)
		.staticmethod("readDecayFile")

		.def("clear", &rpwa::particleDataTable::clear)
		.staticmethod("clear")    

		.add_static_property("debug", &rpwa::particleDataTable::debug, &rpwa::particleDataTable::setDebug);

}

