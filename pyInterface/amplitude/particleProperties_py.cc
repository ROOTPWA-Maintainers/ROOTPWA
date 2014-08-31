#include "particleProperties_py.h"

#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	struct particlePropertiesWrapper : public rpwa::particleProperties,
	                                          bp::wrapper<rpwa::particleProperties>
	{
		particlePropertiesWrapper()
			: rpwa::particleProperties(),
			  bp::wrapper<rpwa::particleProperties>() { }

		particlePropertiesWrapper(const particleProperties& partProp)
			: rpwa::particleProperties(partProp),
			  bp::wrapper<rpwa::particleProperties>() { }

		particlePropertiesWrapper(const std::string& partName,
		                          const int isospin,
		                          const int G,
		                          const int J,
		                          const int P,
		                          const int C)
			: rpwa::particleProperties(partName, isospin, G, J, P, C),
			  bp::wrapper<rpwa::particleProperties>() { }

		std::string qnSummary() const {
			if(bp::override qnSummary = this->get_override("qnSummary")) {
				return qnSummary();
			}
			return rpwa::particleProperties::qnSummary();
		}

		std::string default_qnSummary() const {
			return rpwa::particleProperties::qnSummary();
		}

	};


	bool particleProperties_equal(const rpwa::particleProperties& self, const bp::object& rhsObj) {
		bp::extract<rpwa::particleProperties> get_partProp(rhsObj);
		if(get_partProp.check()) {
			return (self == get_partProp());
		}
		std::pair<rpwa::particleProperties, std::string> rhsPair;
		if(not rpwa::py::convertBPObjectToPair<rpwa::particleProperties, std::string>(rhsObj, rhsPair))
		{
			PyErr_SetString(PyExc_TypeError, "Got invalid input when executing rpwa::particleProperties::operator==()");
			bp::throw_error_already_set();
		}
		return (self == rhsPair);
	}

	bool particleProperties_nequal(const rpwa::particleProperties& self, const bp::object& rhsObj) {
		return not (self == rhsObj);
	}

	bool particleProperties_hasDecay(const rpwa::particleProperties& self, bp::object& pyDaughters) {
		std::multiset<std::string> daughters;
		if(not rpwa::py::convertBPObjectToMultiSet<std::string>(pyDaughters, daughters)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for daughters when executing rpwa::particleProperties::hasDecay()");
			bp::throw_error_already_set();
		}
		return self.hasDecay(daughters);
	}

	void particleProperties_addDecayMode(rpwa::particleProperties& self, bp::object& pyDaughters) {
		std::multiset<std::string> daughters;
		if(not rpwa::py::convertBPObjectToMultiSet<std::string>(pyDaughters, daughters)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for daughters when executing rpwa::particleProperties::addDecayMode()");
			bp::throw_error_already_set();
		}
		self.addDecayMode(daughters);
	}

	bool particleProperties_read(rpwa::particleProperties& self, bp::object& pyLine) {
		std::string strLine = bp::extract<std::string>(pyLine);
		std::istringstream sstrLine(strLine, std::istringstream::in);
		return self.read(sstrLine);
	}

	bp::tuple particleProperties_chargeFromName(const std::string& name) {
		int charge;
		std::string newName = rpwa::particleProperties::chargeFromName(name, charge);
		return bp::make_tuple(newName, charge);
	}

}


void rpwa::py::exportParticleProperties()
{

	bp::class_< particlePropertiesWrapper >( "particleProperties" )

		.def(bp::init<const rpwa::particleProperties&>())
		.def(bp::init<std::string, int, int, int, int, int>())

		.def("__eq__", &particleProperties_equal)
		.def("__neq__", &particleProperties_nequal)

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
		.add_property("isospinProj", &rpwa::particleProperties::isospinProj, &rpwa::particleProperties::setIsospinProj)
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

		.add_property("nmbDecays", &rpwa::particleProperties::nmbDecays)
		.def("hasDecay", &particleProperties_hasDecay)
		.def("addDecayMode", &particleProperties_addDecayMode)

		.add_property("isStable", &rpwa::particleProperties::isStable)

		.def("setSCB", &rpwa::particleProperties::setSCB)
		.def("setIGJPC", &rpwa::particleProperties::setIGJPC)


		.add_property("antiPartProperties", &rpwa::particleProperties::antiPartProperties)

		.def("qnSummary", &particlePropertiesWrapper::qnSummary, &particlePropertiesWrapper::default_qnSummary)
		.def("qnSummary", &rpwa::particleProperties::qnSummary)

		.add_property("bareNameLaTeX", &rpwa::particleProperties::bareNameLaTeX)

		.def("read", &particleProperties_read)

		.def("nameWithCharge", &rpwa::particleProperties::nameWithCharge)
		.staticmethod("nameWithCharge")

		.def("chargeFromName", &particleProperties_chargeFromName)
		.staticmethod("chargeFromName")

		.def("stripChargeFromName", &rpwa::particleProperties::stripChargeFromName)
		.staticmethod("stripChargeFromName")

		.add_static_property("debugParticleProperties", &rpwa::particleProperties::debug, &rpwa::particleProperties::setDebug);

}
