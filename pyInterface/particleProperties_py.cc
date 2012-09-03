#include "particleProperties_py.h"

#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	struct particlePropertiesWrapper : public rpwa::particleProperties,
	                                          bp::wrapper<rpwa::particleProperties>
	{
		particlePropertiesWrapper()
			: rpwa::particleProperties(),
			  bp::wrapper<rpwa::particleProperties>() { };

		particlePropertiesWrapper(const particleProperties& partProp)
			: rpwa::particleProperties(partProp),
			  bp::wrapper<rpwa::particleProperties>() { };

		particlePropertiesWrapper(const std::string& partName,
		                          const int isospin,
		                          const int G,
		                          const int J,
		                          const int P,
		                          const int C)
			: rpwa::particleProperties(partName, isospin, G, J, P, C),
			  bp::wrapper<rpwa::particleProperties>() { };

		bool equal__(const bp::object& rhsObj) {
			bp::extract<particleProperties> get_partProp(rhsObj);
			if(get_partProp.check()) {
				return (*(this) == get_partProp());
			}
			bp::tuple rhs = bp::extract<bp::tuple>(rhsObj);
			rpwa::particleProperties rhsProp = bp::extract<rpwa::particleProperties>(rhs[0]);
			std::string rhsString = bp::extract<std::string>(rhs[1]);
			std::pair<rpwa::particleProperties, std::string> rhsPair;
			rhsPair.first = rhsProp;
			rhsPair.second = rhsString;
			return (*(this) == rhsPair);
		}

		bool nequal__(const bp::object& rhsObj) {
			return not (*(this) == rhsObj);
		}

		bool hasDecay(bp::object& pyDaughters) const {
			bp::list pyDaughtersList = bp::extract<bp::list>(pyDaughters);
			std::set<std::string> daughters = rpwa::py::converBPObjectToStrSet(pyDaughters);
			return rpwa::particleProperties::hasDecay(daughters);
		}

		void addDecayMode(bp::object& pyDaughters) {
			bp::list pyDaughtersList = bp::extract<bp::list>(pyDaughters);
			std::set<std::string> daughters = rpwa::py::converBPObjectToStrSet(pyDaughters);
			rpwa::particleProperties::addDecayMode(daughters);
		}

		std::string qnSummary() const {
			if(bp::override qnSummary = this->get_override("qnSummary")) {
				return qnSummary();
			}
			return rpwa::particleProperties::qnSummary();
		}

		std::string default_qnSummary() const {
			return rpwa::particleProperties::qnSummary();
		}

		bool read__(bp::object& pyLine) {
			std::string strLine = bp::extract<std::string>(pyLine);
			std::istringstream sstrLine(strLine, std::istringstream::in);
			return rpwa::particleProperties::read(sstrLine);
		}

		static bp::tuple chargeFromName__(const std::string& name) {
			int charge;
			std::string new_name = rpwa::particleProperties::chargeFromName(name, charge);
			return bp::make_tuple(new_name, charge);
		}

	};

}

void rpwa::py::exportParticleProperties()
{

	bp::class_< particlePropertiesWrapper >( "particleProperties" )

		.def(bp::init<particlePropertiesWrapper&>())
		.def(bp::init<std::string, int, int, int, int, int>())

		.def("__eq__", &particlePropertiesWrapper::equal__)
		.def("__neq__", &particlePropertiesWrapper::nequal__)

		.def(bp::self_ns::str(bp::self))

		.add_property("name", &particlePropertiesWrapper::name, &particlePropertiesWrapper::setName)
		.add_property("bareName", &particlePropertiesWrapper::bareName)
		.add_property("antiPartName", &particlePropertiesWrapper::antiPartName, &particlePropertiesWrapper::setAntiPartName)
		.add_property("antiPartBareName", &particlePropertiesWrapper::antiPartBareName)
		.add_property("charge", &particlePropertiesWrapper::charge, &particlePropertiesWrapper::setCharge)
		.add_property("mass", &particlePropertiesWrapper::mass, &particlePropertiesWrapper::setMass)
		.add_property("mass2", &particlePropertiesWrapper::mass2)
		.add_property("width", &particlePropertiesWrapper::width, &particlePropertiesWrapper::setWidth)
		.add_property("baryonNmb", &particlePropertiesWrapper::baryonNmb, &particlePropertiesWrapper::setBaryonNmb)
		.add_property("isospin", &particlePropertiesWrapper::isospin, &particlePropertiesWrapper::setIsospin)
		.add_property("isospinProj", &particlePropertiesWrapper::isospinProj/*, &particlePropertiesWrapper::setIsospinProj*/) //<- Uncomment as soon as branch isospin-sym is reintegrated
		.add_property("strangeness", &particlePropertiesWrapper::strangeness, &particlePropertiesWrapper::setStrangeness)
		.add_property("charm", &particlePropertiesWrapper::charm, &particlePropertiesWrapper::setCharm)
		.add_property("beauty", &particlePropertiesWrapper::beauty, &particlePropertiesWrapper::setBeauty)
		.add_property("G", &particlePropertiesWrapper::G, &particlePropertiesWrapper::setG)
		.add_property("J", &particlePropertiesWrapper::J, &particlePropertiesWrapper::setJ)
		.add_property("P", &particlePropertiesWrapper::P, &particlePropertiesWrapper::setP)
		.add_property("C", &particlePropertiesWrapper::C, &particlePropertiesWrapper::setC)

		.add_property("isXParticle", &particlePropertiesWrapper::isXParticle)
		.add_property("isMeson", &particlePropertiesWrapper::isMeson)
		.add_property("isBaryon", &particlePropertiesWrapper::isBaryon)
		.add_property("isLepton", &particlePropertiesWrapper::isLepton)
		.add_property("isPhoton", &particlePropertiesWrapper::isPhoton)
		.add_property("isItsOwnAntiPart", &particlePropertiesWrapper::isItsOwnAntiPart)
		.add_property("isSpinExotic", &particlePropertiesWrapper::isSpinExotic)

		.def(
			"fillFromDataTable"
			, &particlePropertiesWrapper::fillFromDataTable
			, (bp::arg("name"), bp::arg("warnIfNotExistent")=(bool const)(true))
		)

		.add_property("nDecays", &particlePropertiesWrapper::nDecays)
		.def("hasDecay", &particlePropertiesWrapper::hasDecay)
		.def("addDecayMode", &particlePropertiesWrapper::addDecayMode)

		.def("setSCB", &particlePropertiesWrapper::setSCB)
		.def("setIGJPC", &particlePropertiesWrapper::setIGJPC)


		.add_property("antiPartProperties", &particlePropertiesWrapper::antiPartProperties)

		.def("qnSummary", &particlePropertiesWrapper::qnSummary, &particlePropertiesWrapper::default_qnSummary)

		.add_property("bareNameLaTeX", &particlePropertiesWrapper::bareNameLaTeX)

		.def("read", &particlePropertiesWrapper::read__)

		.def("nameWithCharge", &particlePropertiesWrapper::nameWithCharge)
		.staticmethod("nameWithCharge")

		.def("chargeFromName", &particlePropertiesWrapper::chargeFromName__)
		.staticmethod("chargeFromName")

		.def("stripChargeFromName", &particlePropertiesWrapper::stripChargeFromName)
		.staticmethod("stripChargeFromName")

		.add_static_property("debug", &particlePropertiesWrapper::debug, &particlePropertiesWrapper::setDebug);

};

