#ifndef PARTICLEPROPERTIES_PY_H
#define PARTICLEPROPERTIES_PY_H

#include "boost/python.hpp"

#include "particleProperties.h"

namespace bp = boost::python;

struct particleProperties_wrapper : rpwa::particleProperties, bp::wrapper< rpwa::particleProperties > {

	particleProperties_wrapper( )
		: rpwa::particleProperties( )
		  , bp::wrapper< rpwa::particleProperties >(){
			  // null constructor

		  }

	particleProperties_wrapper(::rpwa::particleProperties const & partProp )
		: rpwa::particleProperties( boost::ref(partProp) )
		  , bp::wrapper< rpwa::particleProperties >(){
			  // copy constructor

		  }

	particleProperties_wrapper(::std::string const & partName, int const isospin, int const G, int const J, int const P, int const C )
		: rpwa::particleProperties( partName, isospin, G, J, P, C )
		  , bp::wrapper< rpwa::particleProperties >(){
			  // constructor

		  }

	virtual ::std::string qnSummary(  ) const  {
		if( bp::override func_qnSummary = this->get_override( "qnSummary" ) )
			return func_qnSummary(  );
		else
			return this->rpwa::particleProperties::qnSummary(  );
	}


	::std::string default_qnSummary(  ) const  {
		return rpwa::particleProperties::qnSummary( );
	}

};

void exportParticleProperties();

#endif
