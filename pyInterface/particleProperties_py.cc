#include<particleProperties_py.h>

namespace bp = boost::python;

void exportParticleProperties()
{

	bp::class_< rpwa::particleProperties >( "particleProperties" )

		.def(bp::init<rpwa::particleProperties>())
		.def(bp::init<std::string, int, int, int, int, int>())

		.def(bp::self == bp::self)
		.def(bp::self != bp::self)
//		.def(bp::self == bp::other< std::pair< rpwa::particleProperties, std::string > >())
//		.def(bp::self != bp::other< std::pair< rpwa::particleProperties, std::string > >())
		.def(bp::self_ns::str(bp::self))

		.add_property("name", &rpwa::particleProperties::name, &rpwa::particleProperties::setName)
		.add_property("bareName", &rpwa::particleProperties::bareName)
		.add_property("antiPartName", &rpwa::particleProperties::antiPartName, &rpwa::particleProperties::setAntiPartName)
		.add_property("antiPartBareName", &rpwa::particleProperties::antiPartBareName)
		.add_property("charge", &rpwa::particleProperties::charge, &rpwa::particleProperties::setCharge)
		.add_property("mass", &rpwa::particleProperties::mass, &rpwa::particleProperties::setMass)
		.add_property("mass2", &rpwa::particleProperties::mass2)
		.add_property("width", &rpwa::particleProperties::width, &rpwa::particleProperties::setWidth)
		.add_property("baryonNmb", &rpwa::particleProperties::baryonNmb, &rpwa::particleProperties::setBaryonNmb)
		.add_property("isospin", &rpwa::particleProperties::isospin, &rpwa::particleProperties::setIsospin)
		.add_property("isospinProj", &rpwa::particleProperties::isospinProj/*, &rpwa::particleProperties::setIsospinProj*/) //<- Uncomment as soon as branch isospin-sym is reintegrated
		.add_property("strangeness", &rpwa::particleProperties::strangeness, &rpwa::particleProperties::setStrangeness)
		.add_property("charm", &rpwa::particleProperties::charm, &rpwa::particleProperties::setCharm)
		.add_property("beauty", &rpwa::particleProperties::beauty, &rpwa::particleProperties::setBeauty)
		.add_property("G", &rpwa::particleProperties::G, &rpwa::particleProperties::setG)
		.add_property("J", &rpwa::particleProperties::J, &rpwa::particleProperties::setJ)
		.add_property("P", &rpwa::particleProperties::P, &rpwa::particleProperties::setP)
		.add_property("C", &rpwa::particleProperties::C, &rpwa::particleProperties::setC)

		.add_property("isXParticle", &rpwa::particleProperties::isXParticle)
		.add_property("isMeson", &rpwa::particleProperties::isMeson)
		.add_property("isBaryon", &rpwa::particleProperties::isBaryon)
		.add_property("isLepton", &rpwa::particleProperties::isLepton)
		.add_property("isPhoton", &rpwa::particleProperties::isPhoton)
		.add_property("isItsOwnAntiPart", &rpwa::particleProperties::isItsOwnAntiPart)
		.add_property("isSpinExotic", &rpwa::particleProperties::isSpinExotic)

		.def(
				"fillFromDataTable"
				, &rpwa::particleProperties::fillFromDataTable
				, (bp::arg("name"), bp::arg("warnIfNotExistent")=(bool const)(true))
		)

		.add_property("nDecays", &rpwa::particleProperties::nDecays)
//		.def("hasDecays", &rpwa::particleProperties::hasDecay)
//		.def("addDecayMode", &rpwa::particleProperties::addDecayMode)

		.def("setSCB", &rpwa::particleProperties::setSCB)
		.def("setIGJPC", &rpwa::particleProperties::setIGJPC)


		.add_property("antiPartProperties", &rpwa::particleProperties::antiPartProperties)

		.add_property("qnSummary", &rpwa::particleProperties::qnSummary)

		.add_property("bareNameLaTeX", &rpwa::particleProperties::bareNameLaTeX)

//		.def("read", &rpwa::particleProperties::read)

		.def("nameWithCharge", &rpwa::particleProperties::nameWithCharge)
		.staticmethod("nameWithCharge")
/*		.def("chargeFromName", &rpwa::particleProperties::chargeFromName)
		.staticmethod("chargeFromName")*/
		.def("stripChargeFromName", &rpwa::particleProperties::stripChargeFromName)
		.staticmethod("stripChargeFromName")

		.add_static_property("debug", &rpwa::particleProperties::debug, &rpwa::particleProperties::setDebug);

};

