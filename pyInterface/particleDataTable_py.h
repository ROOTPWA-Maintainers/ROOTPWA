#include "boost/python.hpp"

#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

#include "particleDataTable.h"

namespace bp = boost::python;

void exportParticleDataTable()
{
	{ //::std::vector< std::string >
		typedef bp::class_< std::vector< std::string > > vector_less__std_scope_string__greater__exposer_t;
		vector_less__std_scope_string__greater__exposer_t vector_less__std_scope_string__greater__exposer = vector_less__std_scope_string__greater__exposer_t( "vector_less__std_scope_string__greater_" );
		bp::scope vector_less__std_scope_string__greater__scope( vector_less__std_scope_string__greater__exposer );
		vector_less__std_scope_string__greater__exposer.def( bp::vector_indexing_suite< ::std::vector< std::string >, true >() );
	}

	bp::class_< rpwa::particleDataTable, boost::noncopyable >( "particleDataTable", bp::no_init )    
		.def( 
				"addEntry"
				, (bool (*)( ::rpwa::particleProperties const & ))( &::rpwa::particleDataTable::addEntry )
				, ( bp::arg("partProp") ) )    
		.def( 
				"begin"
				, (::std::_Rb_tree_const_iterator< std::pair< std::string const, rpwa::particleProperties > > (*)(  ))( &::rpwa::particleDataTable::begin ) )    
		.def( 
				"clear"
				, (void (*)(  ))( &::rpwa::particleDataTable::clear ) )    
		.def( 
				"debug"
				, (bool (*)(  ))( &::rpwa::particleDataTable::debug ) )    
		.def( 
				"dump"
				, (::std::ostream & (*)( ::std::ostream & ))( &::rpwa::particleDataTable::dump )
				, ( bp::arg("out") )
				, bp::return_internal_reference<1>() )
		.def( 
				"end"
				, (::std::_Rb_tree_const_iterator< std::pair< std::string const, rpwa::particleProperties > > (*)(  ))( &::rpwa::particleDataTable::end ) )    
/*		.def( 
				"entriesMatching"
				, (::std::vector< const rpwa::particleProperties* > (*)( ::rpwa::particleProperties const &,::std::string const &,double const,double const,::std::vector< std::string > const &,::std::vector< std::string > const &,::std::set< std::string > const &,bool const & ))( &::rpwa::particleDataTable::entriesMatching )
				, ( bp::arg("prototype"), bp::arg("sel"), bp::arg("minMass")=0, bp::arg("minMassWidthFactor")=0, bp::arg("whiteList")=std::vector<std::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >(), bp::arg("blackList")=std::vector<std::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >(), bp::arg("decayproducts")=std::set<std::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<std::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >(), bp::arg("forceDecayCheck")=true ) )    */
		.def( 
				"entry"
				, (::rpwa::particleProperties const * (*)( ::std::string const &,bool const ))( &::rpwa::particleDataTable::entry )
				, ( bp::arg("partName"), bp::arg("warnIfNotExistent")=(bool const)(true) )
				, bp::return_value_policy<bp::reference_existing_object>() )
		.def( 
				"instance"
				, (::rpwa::particleDataTable & (*)(  ))( &::rpwa::particleDataTable::instance )
				, bp::return_value_policy<bp::reference_existing_object>() )
		.def( 
				"isInTable"
				, (bool (*)( ::std::string const & ))( &::rpwa::particleDataTable::isInTable )
				, ( bp::arg("partName") ) )    
		.def( 
				"nmbEntries"
				, (unsigned int (*)(  ))( &::rpwa::particleDataTable::nmbEntries ) )    
		.def( 
				"print"
				, (::std::ostream & (*)( ::std::ostream & ))( &::rpwa::particleDataTable::print )
				, ( bp::arg("out") )
				, bp::return_internal_reference<1>() )
		.def( 
				"read"
				, (bool (*)( ::std::istream & ))( &::rpwa::particleDataTable::read )
				, ( bp::arg("in") ) )    
		.def( 
				"readDecayFile"
				, (bool (*)( ::std::string const & ))( &::rpwa::particleDataTable::readDecayFile )
				, ( bp::arg("fileName") ) )    
		.def( 
				"readFile"
				, (bool (*)( ::std::string const & ))( &::rpwa::particleDataTable::readFile )
				, ( bp::arg("fileName")="./particleDataTable.txt" ) )    
		.def( 
				"setDebug"
				, (void (*)( bool const ))( &::rpwa::particleDataTable::setDebug )
				, ( bp::arg("debug")=(bool const)(true) ) )    
		.staticmethod( "addEntry" )    
		.staticmethod( "begin" )    
		.staticmethod( "clear" )    
		.staticmethod( "debug" )    
		.staticmethod( "dump" )    
		.staticmethod( "end" )    
//		.staticmethod( "entriesMatching" )    
		.staticmethod( "entry" )    
		.staticmethod( "instance" )    
		.staticmethod( "isInTable" )    
		.staticmethod( "nmbEntries" )    
		.staticmethod( "print" )    
		.staticmethod( "read" )    
		.staticmethod( "readDecayFile" )    
		.staticmethod( "readFile" )    
		.staticmethod( "setDebug" )    
		.def( bp::self_ns::str( bp::self ) );

}

