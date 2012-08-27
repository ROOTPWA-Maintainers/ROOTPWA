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

void exportParticleProperties()
{
	{ //::rpwa::particleProperties
		typedef bp::class_< particleProperties_wrapper > particleProperties_exposer_t;
		particleProperties_exposer_t particleProperties_exposer = particleProperties_exposer_t( "particleProperties", bp::init< >() );
		bp::scope particleProperties_scope( particleProperties_exposer );
		particleProperties_exposer.def( bp::init< rpwa::particleProperties const & >(( bp::arg("partProp") )) );
		particleProperties_exposer.def( bp::init< std::string const &, int, int, int, int, int >(( bp::arg("partName"), bp::arg("isospin"), bp::arg("G"), bp::arg("J"), bp::arg("P"), bp::arg("C") )) );
		{ //::rpwa::particleProperties::C

			typedef int ( ::rpwa::particleProperties::*C_function_type )(  ) const;

			particleProperties_exposer.def( 
					"C"
					, C_function_type( &::rpwa::particleProperties::C ) );

		}
		{ //::rpwa::particleProperties::G

			typedef int ( ::rpwa::particleProperties::*G_function_type )(  ) const;

			particleProperties_exposer.def( 
					"G"
					, G_function_type( &::rpwa::particleProperties::G ) );

		}
		{ //::rpwa::particleProperties::J

			typedef int ( ::rpwa::particleProperties::*J_function_type )(  ) const;

			particleProperties_exposer.def( 
					"J"
					, J_function_type( &::rpwa::particleProperties::J ) );

		}
		{ //::rpwa::particleProperties::P

			typedef int ( ::rpwa::particleProperties::*P_function_type )(  ) const;

			particleProperties_exposer.def( 
					"P"
					, P_function_type( &::rpwa::particleProperties::P ) );

		}
/*		{ //::rpwa::particleProperties::addDecayMode

			typedef void ( ::rpwa::particleProperties::*addDecayMode_function_type )( ::std::set< std::string > const & ) ;

			particleProperties_exposer.def( 
					"addDecayMode"
					, addDecayMode_function_type( &::rpwa::particleProperties::addDecayMode )
					, ( bp::arg("daughters") ) );

		}*/
		{ //::rpwa::particleProperties::antiPartBareName

			typedef ::std::string ( ::rpwa::particleProperties::*antiPartBareName_function_type )(  ) const;

			particleProperties_exposer.def( 
					"antiPartBareName"
					, antiPartBareName_function_type( &::rpwa::particleProperties::antiPartBareName ) );

		}
		{ //::rpwa::particleProperties::antiPartName

			typedef ::std::string ( ::rpwa::particleProperties::*antiPartName_function_type )(  ) const;

			particleProperties_exposer.def( 
					"antiPartName"
					, antiPartName_function_type( &::rpwa::particleProperties::antiPartName ) );

		}
		{ //::rpwa::particleProperties::antiPartProperties

			typedef ::rpwa::particleProperties ( ::rpwa::particleProperties::*antiPartProperties_function_type )(  ) const;

			particleProperties_exposer.def( 
					"antiPartProperties"
					, antiPartProperties_function_type( &::rpwa::particleProperties::antiPartProperties ) );

		}
		{ //::rpwa::particleProperties::bareName

			typedef ::std::string ( ::rpwa::particleProperties::*bareName_function_type )(  ) const;

			particleProperties_exposer.def( 
					"bareName"
					, bareName_function_type( &::rpwa::particleProperties::bareName ) );

		}
		{ //::rpwa::particleProperties::bareNameLaTeX

			typedef ::std::string ( ::rpwa::particleProperties::*bareNameLaTeX_function_type )(  ) const;

			particleProperties_exposer.def( 
					"bareNameLaTeX"
					, bareNameLaTeX_function_type( &::rpwa::particleProperties::bareNameLaTeX ) );

		}
		{ //::rpwa::particleProperties::baryonNmb

			typedef int ( ::rpwa::particleProperties::*baryonNmb_function_type )(  ) const;

			particleProperties_exposer.def( 
					"baryonNmb"
					, baryonNmb_function_type( &::rpwa::particleProperties::baryonNmb ) );

		}
		{ //::rpwa::particleProperties::beauty

			typedef int ( ::rpwa::particleProperties::*beauty_function_type )(  ) const;

			particleProperties_exposer.def( 
					"beauty"
					, beauty_function_type( &::rpwa::particleProperties::beauty ) );

		}
		{ //::rpwa::particleProperties::charge

			typedef int ( ::rpwa::particleProperties::*charge_function_type )(  ) const;

			particleProperties_exposer.def( 
					"charge"
					, charge_function_type( &::rpwa::particleProperties::charge ) );

		}
		{ //::rpwa::particleProperties::chargeFromName

			typedef ::std::string ( *chargeFromName_function_type )( ::std::string const &,int & );

			particleProperties_exposer.def( 
					"chargeFromName"
					, chargeFromName_function_type( &::rpwa::particleProperties::chargeFromName )
					, ( bp::arg("name"), bp::arg("charge") ) );

		}
		{ //::rpwa::particleProperties::charm

			typedef int ( ::rpwa::particleProperties::*charm_function_type )(  ) const;

			particleProperties_exposer.def( 
					"charm"
					, charm_function_type( &::rpwa::particleProperties::charm ) );

		}
		{ //::rpwa::particleProperties::debug

			typedef bool ( *debug_function_type )(  );

			particleProperties_exposer.def( 
					"debug"
					, debug_function_type( &::rpwa::particleProperties::debug ) );

		}
		{ //::rpwa::particleProperties::dump

			typedef ::std::ostream & ( ::rpwa::particleProperties::*dump_function_type )( ::std::ostream & ) const;

			particleProperties_exposer.def( 
					"dump"
					, dump_function_type(&::rpwa::particleProperties::dump)
					, ( bp::arg("out") )
					, bp::return_internal_reference<1>() );

		}
		{ //::rpwa::particleProperties::fillFromDataTable

			typedef bool ( ::rpwa::particleProperties::*fillFromDataTable_function_type )( ::std::string const &,bool const ) ;

			particleProperties_exposer.def( 
					"fillFromDataTable"
					, fillFromDataTable_function_type( &::rpwa::particleProperties::fillFromDataTable )
					, ( bp::arg("name"), bp::arg("warnIfNotExistent")=(bool const)(true) ) );

		}
/*		{ //::rpwa::particleProperties::hasDecay

			typedef bool ( ::rpwa::particleProperties::*hasDecay_function_type )( ::std::set< std::string > const & ) const;

			particleProperties_exposer.def( 
					"hasDecay"
					, hasDecay_function_type( &::rpwa::particleProperties::hasDecay )
					, ( bp::arg("daughters") ) );

		}*/
		{ //::rpwa::particleProperties::isBaryon

			typedef bool ( ::rpwa::particleProperties::*isBaryon_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isBaryon"
					, isBaryon_function_type( &::rpwa::particleProperties::isBaryon ) );

		}
		{ //::rpwa::particleProperties::isItsOwnAntiPart

			typedef bool ( ::rpwa::particleProperties::*isItsOwnAntiPart_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isItsOwnAntiPart"
					, isItsOwnAntiPart_function_type( &::rpwa::particleProperties::isItsOwnAntiPart ) );

		}
		{ //::rpwa::particleProperties::isLepton

			typedef bool ( ::rpwa::particleProperties::*isLepton_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isLepton"
					, isLepton_function_type( &::rpwa::particleProperties::isLepton ) );

		}
		{ //::rpwa::particleProperties::isMeson

			typedef bool ( ::rpwa::particleProperties::*isMeson_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isMeson"
					, isMeson_function_type( &::rpwa::particleProperties::isMeson ) );

		}
		{ //::rpwa::particleProperties::isPhoton

			typedef bool ( ::rpwa::particleProperties::*isPhoton_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isPhoton"
					, isPhoton_function_type( &::rpwa::particleProperties::isPhoton ) );

		}
		{ //::rpwa::particleProperties::isSpinExotic

			typedef bool ( ::rpwa::particleProperties::*isSpinExotic_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isSpinExotic"
					, isSpinExotic_function_type( &::rpwa::particleProperties::isSpinExotic ) );

		}
		{ //::rpwa::particleProperties::isXParticle

			typedef bool ( ::rpwa::particleProperties::*isXParticle_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isXParticle"
					, isXParticle_function_type( &::rpwa::particleProperties::isXParticle ) );

		}
		{ //::rpwa::particleProperties::isospin

			typedef int ( ::rpwa::particleProperties::*isospin_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isospin"
					, isospin_function_type( &::rpwa::particleProperties::isospin ) );

		}
		{ //::rpwa::particleProperties::isospinProj

			typedef int ( ::rpwa::particleProperties::*isospinProj_function_type )(  ) const;

			particleProperties_exposer.def( 
					"isospinProj"
					, isospinProj_function_type( &::rpwa::particleProperties::isospinProj ) );

		}
		{ //::rpwa::particleProperties::mass

			typedef double ( ::rpwa::particleProperties::*mass_function_type )(  ) const;

			particleProperties_exposer.def( 
					"mass"
					, mass_function_type( &::rpwa::particleProperties::mass ) );

		}
		{ //::rpwa::particleProperties::mass2

			typedef double ( ::rpwa::particleProperties::*mass2_function_type )(  ) const;

			particleProperties_exposer.def( 
					"mass2"
					, mass2_function_type( &::rpwa::particleProperties::mass2 ) );

		}
		{ //::rpwa::particleProperties::nDecays

			typedef int ( ::rpwa::particleProperties::*nDecays_function_type )(  ) const;

			particleProperties_exposer.def( 
					"nDecays"
					, nDecays_function_type( &::rpwa::particleProperties::nDecays ) );

		}
		{ //::rpwa::particleProperties::name

			typedef ::std::string ( ::rpwa::particleProperties::*name_function_type )(  ) const;

			particleProperties_exposer.def( 
					"name"
					, name_function_type( &::rpwa::particleProperties::name ) );

		}
		{ //::rpwa::particleProperties::nameWithCharge

			typedef ::std::string ( *nameWithCharge_function_type )( ::std::string const &,int const );

			particleProperties_exposer.def( 
					"nameWithCharge"
					, nameWithCharge_function_type( &::rpwa::particleProperties::nameWithCharge )
					, ( bp::arg("bareName"), bp::arg("charge") ) );

		}
		{ //::rpwa::particleProperties::operator=

			typedef ::rpwa::particleProperties & ( ::rpwa::particleProperties::*assign_function_type )( ::rpwa::particleProperties const & ) ;

			particleProperties_exposer.def( 
					"assign"
					, assign_function_type( &::rpwa::particleProperties::operator= )
					, ( bp::arg("partProp") )
					, bp::return_self< >() );

		}
		{ //::rpwa::particleProperties::print

			typedef ::std::ostream & ( ::rpwa::particleProperties::*print_function_type )( ::std::ostream & ) const;

			particleProperties_exposer.def( 
					"print"
					, print_function_type(&::rpwa::particleProperties::print)
					, ( bp::arg("out") )
					, bp::return_internal_reference<1>() );

		}
		{ //::rpwa::particleProperties::qnSummary

			typedef ::std::string ( ::rpwa::particleProperties::*qnSummary_function_type )(  ) const;
			typedef ::std::string ( particleProperties_wrapper::*default_qnSummary_function_type )(  ) const;

			particleProperties_exposer.def( 
					"qnSummary"
					, qnSummary_function_type(&::rpwa::particleProperties::qnSummary)
					, default_qnSummary_function_type(&particleProperties_wrapper::default_qnSummary) );

		}
		{ //::rpwa::particleProperties::read

			typedef bool ( ::rpwa::particleProperties::*read_function_type )( ::std::istringstream & ) ;

			particleProperties_exposer.def( 
					"read"
					, read_function_type( &::rpwa::particleProperties::read )
					, ( bp::arg("line") ) );

		}
		{ //::rpwa::particleProperties::setAntiPartName

			typedef void ( ::rpwa::particleProperties::*setAntiPartName_function_type )( ::std::string const & ) ;

			particleProperties_exposer.def( 
					"setAntiPartName"
					, setAntiPartName_function_type( &::rpwa::particleProperties::setAntiPartName )
					, ( bp::arg("name") ) );

		}
		{ //::rpwa::particleProperties::setBaryonNmb

			typedef void ( ::rpwa::particleProperties::*setBaryonNmb_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setBaryonNmb"
					, setBaryonNmb_function_type( &::rpwa::particleProperties::setBaryonNmb )
					, ( bp::arg("baryonNmb") ) );

		}
		{ //::rpwa::particleProperties::setBeauty

			typedef void ( ::rpwa::particleProperties::*setBeauty_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setBeauty"
					, setBeauty_function_type( &::rpwa::particleProperties::setBeauty )
					, ( bp::arg("beauty") ) );

		}
		{ //::rpwa::particleProperties::setC

			typedef void ( ::rpwa::particleProperties::*setC_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setC"
					, setC_function_type( &::rpwa::particleProperties::setC )
					, ( bp::arg("C") ) );

		}
		{ //::rpwa::particleProperties::setCharge

			typedef void ( ::rpwa::particleProperties::*setCharge_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setCharge"
					, setCharge_function_type( &::rpwa::particleProperties::setCharge )
					, ( bp::arg("charge") ) );

		}
		{ //::rpwa::particleProperties::setCharm

			typedef void ( ::rpwa::particleProperties::*setCharm_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setCharm"
					, setCharm_function_type( &::rpwa::particleProperties::setCharm )
					, ( bp::arg("charm") ) );

		}
		{ //::rpwa::particleProperties::setDebug

			typedef void ( *setDebug_function_type )( bool const );

			particleProperties_exposer.def( 
					"setDebug"
					, setDebug_function_type( &::rpwa::particleProperties::setDebug )
					, ( bp::arg("debug")=(bool const)(true) ) );

		}
		{ //::rpwa::particleProperties::setG

			typedef void ( ::rpwa::particleProperties::*setG_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setG"
					, setG_function_type( &::rpwa::particleProperties::setG )
					, ( bp::arg("G") ) );

		}
		{ //::rpwa::particleProperties::setIGJPC

			typedef void ( ::rpwa::particleProperties::*setIGJPC_function_type )( int const,int const,int const,int const,int const ) ;

			particleProperties_exposer.def( 
					"setIGJPC"
					, setIGJPC_function_type( &::rpwa::particleProperties::setIGJPC )
					, ( bp::arg("isospin"), bp::arg("G"), bp::arg("J"), bp::arg("P"), bp::arg("C") ) );

		}
		{ //::rpwa::particleProperties::setIsospin

			typedef void ( ::rpwa::particleProperties::*setIsospin_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setIsospin"
					, setIsospin_function_type( &::rpwa::particleProperties::setIsospin )
					, ( bp::arg("isospin") ) );

		}
		{ //::rpwa::particleProperties::setJ

			typedef void ( ::rpwa::particleProperties::*setJ_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setJ"
					, setJ_function_type( &::rpwa::particleProperties::setJ )
					, ( bp::arg("J") ) );

		}
		{ //::rpwa::particleProperties::setMass

			typedef void ( ::rpwa::particleProperties::*setMass_function_type )( double const ) ;

			particleProperties_exposer.def( 
					"setMass"
					, setMass_function_type( &::rpwa::particleProperties::setMass )
					, ( bp::arg("mass") ) );

		}
		{ //::rpwa::particleProperties::setName

			typedef void ( ::rpwa::particleProperties::*setName_function_type )( ::std::string const & ) ;

			particleProperties_exposer.def( 
					"setName"
					, setName_function_type( &::rpwa::particleProperties::setName )
					, ( bp::arg("partName") ) );

		}
		{ //::rpwa::particleProperties::setP

			typedef void ( ::rpwa::particleProperties::*setP_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setP"
					, setP_function_type( &::rpwa::particleProperties::setP )
					, ( bp::arg("P") ) );

		}
		{ //::rpwa::particleProperties::setSCB

			typedef void ( ::rpwa::particleProperties::*setSCB_function_type )( int const,int const,int const ) ;

			particleProperties_exposer.def( 
					"setSCB"
					, setSCB_function_type( &::rpwa::particleProperties::setSCB )
					, ( bp::arg("strangeness"), bp::arg("charm"), bp::arg("beauty") ) );

		}
		{ //::rpwa::particleProperties::setStrangeness

			typedef void ( ::rpwa::particleProperties::*setStrangeness_function_type )( int const ) ;

			particleProperties_exposer.def( 
					"setStrangeness"
					, setStrangeness_function_type( &::rpwa::particleProperties::setStrangeness )
					, ( bp::arg("strangeness") ) );

		}
		{ //::rpwa::particleProperties::setWidth

			typedef void ( ::rpwa::particleProperties::*setWidth_function_type )( double const ) ;

			particleProperties_exposer.def( 
					"setWidth"
					, setWidth_function_type( &::rpwa::particleProperties::setWidth )
					, ( bp::arg("width") ) );

		}
		{ //::rpwa::particleProperties::strangeness

			typedef int ( ::rpwa::particleProperties::*strangeness_function_type )(  ) const;

			particleProperties_exposer.def( 
					"strangeness"
					, strangeness_function_type( &::rpwa::particleProperties::strangeness ) );

		}
		{ //::rpwa::particleProperties::stripChargeFromName

			typedef ::std::string ( *stripChargeFromName_function_type )( ::std::string const & );

			particleProperties_exposer.def( 
					"stripChargeFromName"
					, stripChargeFromName_function_type( &::rpwa::particleProperties::stripChargeFromName )
					, ( bp::arg("name") ) );

		}
		{ //::rpwa::particleProperties::width

			typedef double ( ::rpwa::particleProperties::*width_function_type )(  ) const;

			particleProperties_exposer.def( 
					"width"
					, width_function_type( &::rpwa::particleProperties::width ) );

		}
		particleProperties_exposer.staticmethod( "chargeFromName" );
		particleProperties_exposer.staticmethod( "debug" );
		particleProperties_exposer.staticmethod( "nameWithCharge" );
		particleProperties_exposer.staticmethod( "setDebug" );
		particleProperties_exposer.staticmethod( "stripChargeFromName" );
		particleProperties_exposer.def( bp::self != bp::other< std::pair< rpwa::particleProperties, std::string > >() );
		particleProperties_exposer.def( bp::self != bp::self );
		particleProperties_exposer.def( bp::self_ns::str( bp::self ) );
		particleProperties_exposer.def( bp::self == bp::other< std::pair< rpwa::particleProperties, std::string > >() );
		particleProperties_exposer.def( bp::self == bp::self );
	}
}
