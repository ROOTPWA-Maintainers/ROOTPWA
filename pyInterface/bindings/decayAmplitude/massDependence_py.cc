#include "massDependence_py.h"

#include <boost/python.hpp>

#include "isobarDecayVertex.h"
#include "massDependence.h"

namespace bp = boost::python;


namespace {

	struct massDependenceWrapper : public rpwa::massDependence,
	                                      bp::wrapper<rpwa::massDependence>
	{

		massDependenceWrapper()
			: rpwa::massDependence(),
			  bp::wrapper<rpwa::massDependence>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			return this->get_override("amp")(v);
		}

		std::complex<double> operator()(const rpwa::isobarDecayVertex& v) {
			if(bp::override oper = this->get_override("__call__")) {
				return oper(v);
			}
			return rpwa::massDependence::operator()(v);
		}

		std::complex<double> default___call__(rpwa::isobarDecayVertex& v) {
			return rpwa::massDependence::operator()(v);
		}

		std::string name() const {
			return this->get_override("name")();
		}

	};

	struct flatMassDependenceWrapper : public rpwa::flatMassDependence,
	                                          bp::wrapper<rpwa::flatMassDependence>
	{

		flatMassDependenceWrapper()
			: rpwa::flatMassDependence(),
			  bp::wrapper<rpwa::flatMassDependence>() { }

		flatMassDependenceWrapper(const rpwa::flatMassDependence& dep)
			: rpwa::flatMassDependence(dep),
			  bp::wrapper<rpwa::flatMassDependence>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::flatMassDependence::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::flatMassDependence::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::flatMassDependence::name();
		}

		std::string default_name() const {
			return rpwa::flatMassDependence::name();
		}

	};

	struct relativisticBreitWignerWrapper : public rpwa::relativisticBreitWigner,
	                                               bp::wrapper<rpwa::relativisticBreitWigner>
	{

		relativisticBreitWignerWrapper()
			: rpwa::relativisticBreitWigner(),
			  bp::wrapper<rpwa::relativisticBreitWigner>() { }

		relativisticBreitWignerWrapper(const rpwa::relativisticBreitWigner& dep)
			: rpwa::relativisticBreitWigner(dep),
			  bp::wrapper<rpwa::relativisticBreitWigner>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::relativisticBreitWigner::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::relativisticBreitWigner::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::relativisticBreitWigner::name();
		}

		std::string default_name() const {
			return rpwa::relativisticBreitWigner::name();
		}

	};

	struct constWidthBreitWignerWrapper : public rpwa::constWidthBreitWigner,
		                                               bp::wrapper<rpwa::constWidthBreitWigner>
	{

		constWidthBreitWignerWrapper()
			: rpwa::constWidthBreitWigner(),
			  bp::wrapper<rpwa::constWidthBreitWigner>() { }

		constWidthBreitWignerWrapper(const rpwa::constWidthBreitWigner& dep)
			: rpwa::constWidthBreitWigner(dep),
			  bp::wrapper<rpwa::constWidthBreitWigner>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::constWidthBreitWigner::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::constWidthBreitWigner::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::constWidthBreitWigner::name();
		}

		std::string default_name() const {
			return rpwa::constWidthBreitWigner::name();
		}

	};

	struct rhoBreitWignerWrapper : public rpwa::rhoBreitWigner,
		                                  bp::wrapper<rpwa::rhoBreitWigner>
	{

		rhoBreitWignerWrapper()
			: rpwa::rhoBreitWigner(),
			  bp::wrapper<rpwa::rhoBreitWigner>() { }

		rhoBreitWignerWrapper(const rpwa::rhoBreitWigner& dep)
			: rpwa::rhoBreitWigner(dep),
			  bp::wrapper<rpwa::rhoBreitWigner>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::rhoBreitWigner::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::rhoBreitWigner::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::rhoBreitWigner::name();
		}

		std::string default_name() const {
			return rpwa::rhoBreitWigner::name();
		}

	};

	struct f0980BreitWignerWrapper : public rpwa::f0980BreitWigner,
		                                  bp::wrapper<rpwa::f0980BreitWigner>
	{

		f0980BreitWignerWrapper()
			: rpwa::f0980BreitWigner(),
			  bp::wrapper<rpwa::f0980BreitWigner>() { }

		f0980BreitWignerWrapper(const rpwa::f0980BreitWigner& dep)
			: rpwa::f0980BreitWigner(dep),
			  bp::wrapper<rpwa::f0980BreitWigner>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::f0980BreitWigner::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::f0980BreitWigner::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::f0980BreitWigner::name();
		}

		std::string default_name() const {
			return rpwa::f0980BreitWigner::name();
		}

	};

	struct piPiSWaveAuMorganPenningtonMWrapper : public rpwa::piPiSWaveAuMorganPenningtonM,
	                                                    bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonM>
	{

		piPiSWaveAuMorganPenningtonMWrapper()
			: rpwa::piPiSWaveAuMorganPenningtonM(),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonM>() { }

		piPiSWaveAuMorganPenningtonMWrapper(const rpwa::piPiSWaveAuMorganPenningtonM& dep)
			: rpwa::piPiSWaveAuMorganPenningtonM(dep),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonM>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::piPiSWaveAuMorganPenningtonM::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::piPiSWaveAuMorganPenningtonM::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::piPiSWaveAuMorganPenningtonM::name();
		}

		std::string default_name() const {
			return rpwa::piPiSWaveAuMorganPenningtonM::name();
		}

	};

	struct piPiSWaveAuMorganPenningtonVesWrapper : public rpwa::piPiSWaveAuMorganPenningtonVes,
	                                                      bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonVes>
	{

		piPiSWaveAuMorganPenningtonVesWrapper()
			: rpwa::piPiSWaveAuMorganPenningtonVes(),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonVes>() { }

		piPiSWaveAuMorganPenningtonVesWrapper(const rpwa::piPiSWaveAuMorganPenningtonVes& dep)
			: rpwa::piPiSWaveAuMorganPenningtonVes(dep),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonVes>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::piPiSWaveAuMorganPenningtonVes::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::piPiSWaveAuMorganPenningtonVes::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::piPiSWaveAuMorganPenningtonVes::name();
		}

		std::string default_name() const {
			return rpwa::piPiSWaveAuMorganPenningtonVes::name();
		}

	};

	struct piPiSWaveAuMorganPenningtonKachaevWrapper : public rpwa::piPiSWaveAuMorganPenningtonKachaev,
	                                                          bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonKachaev>
	{

		piPiSWaveAuMorganPenningtonKachaevWrapper()
			: rpwa::piPiSWaveAuMorganPenningtonKachaev(),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonKachaev>() { }

		piPiSWaveAuMorganPenningtonKachaevWrapper(const rpwa::piPiSWaveAuMorganPenningtonKachaev& dep)
			: rpwa::piPiSWaveAuMorganPenningtonKachaev(dep),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonKachaev>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::piPiSWaveAuMorganPenningtonKachaev::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::piPiSWaveAuMorganPenningtonKachaev::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::piPiSWaveAuMorganPenningtonKachaev::name();
		}

		std::string default_name() const {
			return rpwa::piPiSWaveAuMorganPenningtonKachaev::name();
		}

	};

	struct rhoPrimeMassDepWrapper : public rpwa::rhoPrimeMassDep,
	                                bp::wrapper<rpwa::rhoPrimeMassDep>
	{

		rhoPrimeMassDepWrapper()
			: rpwa::rhoPrimeMassDep(),
			  bp::wrapper<rpwa::rhoPrimeMassDep>() { }

		rhoPrimeMassDepWrapper(const rpwa::rhoPrimeMassDep& dep)
			: rpwa::rhoPrimeMassDep(dep),
			  bp::wrapper<rpwa::rhoPrimeMassDep>() { }

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::rhoPrimeMassDep::amp(v);
		}

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::rhoPrimeMassDep::amp(v);
		}

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::rhoPrimeMassDep::name();
		}

		std::string default_name() const {
			return rpwa::rhoPrimeMassDep::name();
		}

	};

}

void rpwa::py::exportMassDependence() {

	bp::class_<massDependenceWrapper, boost::noncopyable>("massDependence", bp::no_init)
		.def(bp::self_ns::str(bp::self))
		.def("amp", bp::pure_virtual(&rpwa::massDependence::amp))
		.def(
			"__call__"
			, &massDependenceWrapper::operator()
			, (::std::complex< double > ( massDependenceWrapper::* )( ::rpwa::isobarDecayVertex const & ) )(&massDependenceWrapper::default___call__)
		)
		.def("__call__", &rpwa::massDependence::operator())
		.def("name", bp::pure_virtual(&rpwa::massDependence::name))
		.add_static_property("debugMassDependence", &rpwa::massDependence::debug, &rpwa::massDependence::setDebug);

	bp::class_<flatMassDependenceWrapper, bp::bases<rpwa::massDependence> >("flatMassDependence")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &flatMassDependenceWrapper::amp, &flatMassDependenceWrapper::default_amp)
		.def("amp", &rpwa::flatMassDependence::amp)
		.def("name", &flatMassDependenceWrapper::name, &flatMassDependenceWrapper::default_name)
		.def("name", &rpwa::flatMassDependence::name);

	bp::class_<relativisticBreitWignerWrapper, bp::bases<rpwa::massDependence> >("relativisticBreitWigner")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &relativisticBreitWignerWrapper::amp, &relativisticBreitWignerWrapper::default_amp)
		.def("amp", &rpwa::relativisticBreitWigner::amp)
		.def("name", &relativisticBreitWignerWrapper::name, &relativisticBreitWignerWrapper::default_name)
		.def("name", &rpwa::relativisticBreitWigner::name);

	bp::class_<constWidthBreitWignerWrapper, bp::bases<rpwa::massDependence> >("constWidthBreitWigner")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &constWidthBreitWignerWrapper::amp, &constWidthBreitWignerWrapper::default_amp)
		.def("amp", &rpwa::constWidthBreitWigner::amp)
		.def("name", &constWidthBreitWignerWrapper::name, &constWidthBreitWignerWrapper::default_name)
		.def("name", &rpwa::constWidthBreitWigner::name);

	bp::class_<rhoBreitWignerWrapper, bp::bases<rpwa::massDependence> >("rhoBreitWigner")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &rhoBreitWignerWrapper::amp, &rhoBreitWignerWrapper::default_amp)
		.def("amp", &rpwa::rhoBreitWigner::amp)
		.def("name", &rhoBreitWignerWrapper::name, &rhoBreitWignerWrapper::default_name)
		.def("name", &rpwa::rhoBreitWigner::name);

	bp::class_<f0980BreitWignerWrapper, bp::bases<rpwa::massDependence> >("f0980BreitWigner")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &f0980BreitWignerWrapper::amp, &f0980BreitWignerWrapper::default_amp)
		.def("amp", &rpwa::f0980BreitWigner::amp)
		.def("name", &f0980BreitWignerWrapper::name, &f0980BreitWignerWrapper::default_name)
		.def("name", &rpwa::f0980BreitWigner::name);

	bp::class_<piPiSWaveAuMorganPenningtonMWrapper, bp::bases<rpwa::massDependence> >("piPiSWaveAuMorganPenningtonM")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &piPiSWaveAuMorganPenningtonMWrapper::amp, &piPiSWaveAuMorganPenningtonMWrapper::default_amp)
		.def("amp", &rpwa::piPiSWaveAuMorganPenningtonM::amp)
		.def("name", &piPiSWaveAuMorganPenningtonMWrapper::name, &piPiSWaveAuMorganPenningtonMWrapper::default_name)
		.def("name", &rpwa::piPiSWaveAuMorganPenningtonM::name);

	bp::class_<piPiSWaveAuMorganPenningtonVesWrapper, bp::bases<rpwa::massDependence> >("piPiSWaveAuMorganPenningtonVes")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &piPiSWaveAuMorganPenningtonVesWrapper::amp, &piPiSWaveAuMorganPenningtonVesWrapper::default_amp)
		.def("amp", &rpwa::piPiSWaveAuMorganPenningtonVes::amp)
		.def("name", &piPiSWaveAuMorganPenningtonVesWrapper::name, &piPiSWaveAuMorganPenningtonVesWrapper::default_name)
		.def("name", &rpwa::piPiSWaveAuMorganPenningtonVes::name);

	bp::class_<piPiSWaveAuMorganPenningtonKachaevWrapper, bp::bases<rpwa::massDependence> >("piPiSWaveAuMorganPenningtonKachaev")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &piPiSWaveAuMorganPenningtonKachaevWrapper::amp, &piPiSWaveAuMorganPenningtonKachaevWrapper::default_amp)
		.def("amp", &rpwa::piPiSWaveAuMorganPenningtonKachaev::amp)
		.def("name", &piPiSWaveAuMorganPenningtonKachaevWrapper::name, &piPiSWaveAuMorganPenningtonKachaevWrapper::default_name)
		.def("name", &rpwa::piPiSWaveAuMorganPenningtonKachaev::name);

	bp::class_<rhoPrimeMassDepWrapper, bp::bases<rpwa::massDependence> >("rhoPrimeMassDep")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &rhoPrimeMassDepWrapper::amp, &rhoPrimeMassDepWrapper::default_amp)
		.def("amp", &rpwa::rhoPrimeMassDep::amp)
		.def("name", &rhoPrimeMassDepWrapper::name, &rhoPrimeMassDepWrapper::default_name)
		.def("name", &rpwa::rhoPrimeMassDep::name);

	bp::register_ptr_to_python<rpwa::massDependencePtr>();
	bp::register_ptr_to_python<rpwa::flatMassDependencePtr>();
	bp::register_ptr_to_python<rpwa::relativisticBreitWignerPtr>();
	bp::register_ptr_to_python<rpwa::constWidthBreitWignerPtr>();
	bp::register_ptr_to_python<rpwa::rhoBreitWignerPtr>();
	bp::register_ptr_to_python<rpwa::f0980BreitWignerPtr>();
	bp::register_ptr_to_python<rpwa::piPiSWaveAuMorganPenningtonMPtr>();
	bp::register_ptr_to_python<rpwa::piPiSWaveAuMorganPenningtonVesPtr>();
	bp::register_ptr_to_python<rpwa::piPiSWaveAuMorganPenningtonKachaevPtr>();
	bp::register_ptr_to_python<rpwa::rhoPrimeMassDepPtr>();

}
