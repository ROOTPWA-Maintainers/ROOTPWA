#include "massDependence_py.h"

#include<isobarDecayVertex.h>

namespace bp = boost::python;

namespace {

	struct massDependenceWrapper : public rpwa::massDependence,
	                                      bp::wrapper<rpwa::massDependence>
	{

		massDependenceWrapper()
			: rpwa::massDependence(),
			  bp::wrapper<rpwa::massDependence>() { };

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			return this->get_override("amp")(v);
		};

		std::complex<double> operator()(const rpwa::isobarDecayVertex& v) {
			if(bp::override oper = this->get_override("__call__")) {
				return oper(v);
			}
			return rpwa::massDependence::operator()(v);
		};

		std::complex<double> default___call__(rpwa::isobarDecayVertex& v) {
			return rpwa::massDependence::operator()(v);
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::massDependence::name();
		};

		std::string default_name() const {
			return rpwa::massDependence::name();
		};

	};

	struct flatMassDependenceWrapper : public rpwa::flatMassDependence,
	                                          bp::wrapper<rpwa::flatMassDependence>
	{

		flatMassDependenceWrapper()
			: rpwa::flatMassDependence(),
			  bp::wrapper<rpwa::flatMassDependence>() { };

		flatMassDependenceWrapper(const flatMassDependence& dep)
			: rpwa::flatMassDependence(dep),
			  bp::wrapper<rpwa::flatMassDependence>() { };

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::flatMassDependence::amp(v);
		};

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::flatMassDependence::amp(v);
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::flatMassDependence::name();
		};

		std::string default_name() const {
			return rpwa::flatMassDependence::name();
		};

	};

	struct relativisticBreitWignerWrapper : public rpwa::relativisticBreitWigner,
	                                               bp::wrapper<rpwa::relativisticBreitWigner>
	{

		relativisticBreitWignerWrapper()
			: rpwa::relativisticBreitWigner(),
			  bp::wrapper<rpwa::relativisticBreitWigner>() { };

		relativisticBreitWignerWrapper(const relativisticBreitWigner& dep)
			: rpwa::relativisticBreitWigner(dep),
			  bp::wrapper<rpwa::relativisticBreitWigner>() { };

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::relativisticBreitWigner::amp(v);
		};

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::relativisticBreitWigner::amp(v);
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::relativisticBreitWigner::name();
		};

		std::string default_name() const {
			return rpwa::relativisticBreitWigner::name();
		};

	};

	struct piPiSWaveAuMorganPenningtonMWrapper : public rpwa::piPiSWaveAuMorganPenningtonM,
	                                                    bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonM>
	{

		piPiSWaveAuMorganPenningtonMWrapper()
			: rpwa::piPiSWaveAuMorganPenningtonM(),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonM>() { };

		piPiSWaveAuMorganPenningtonMWrapper(const piPiSWaveAuMorganPenningtonM& dep)
			: rpwa::piPiSWaveAuMorganPenningtonM(dep),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonM>() { };

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::piPiSWaveAuMorganPenningtonM::amp(v);
		};

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::piPiSWaveAuMorganPenningtonM::amp(v);
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::piPiSWaveAuMorganPenningtonM::name();
		};

		std::string default_name() const {
			return rpwa::piPiSWaveAuMorganPenningtonM::name();
		};

	};

	struct piPiSWaveAuMorganPenningtonVesWrapper : public rpwa::piPiSWaveAuMorganPenningtonVes,
	                                                      bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonVes>
	{

		piPiSWaveAuMorganPenningtonVesWrapper()
			: rpwa::piPiSWaveAuMorganPenningtonVes(),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonVes>() { };

		piPiSWaveAuMorganPenningtonVesWrapper(const piPiSWaveAuMorganPenningtonVes& dep)
			: rpwa::piPiSWaveAuMorganPenningtonVes(dep),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonVes>() { };

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::piPiSWaveAuMorganPenningtonVes::amp(v);
		};

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::piPiSWaveAuMorganPenningtonVes::amp(v);
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::piPiSWaveAuMorganPenningtonVes::name();
		};

		std::string default_name() const {
			return rpwa::piPiSWaveAuMorganPenningtonVes::name();
		};

	};

	struct piPiSWaveAuMorganPenningtonKachaevWrapper : public rpwa::piPiSWaveAuMorganPenningtonKachaev,
	                                                          bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonKachaev>
	{

		piPiSWaveAuMorganPenningtonKachaevWrapper()
			: rpwa::piPiSWaveAuMorganPenningtonKachaev(),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonKachaev>() { };

		piPiSWaveAuMorganPenningtonKachaevWrapper(const piPiSWaveAuMorganPenningtonKachaev& dep)
			: rpwa::piPiSWaveAuMorganPenningtonKachaev(dep),
			  bp::wrapper<rpwa::piPiSWaveAuMorganPenningtonKachaev>() { };

		std::complex<double> amp(const rpwa::isobarDecayVertex& v) {
			if(bp::override amp = this->get_override("amp")) {
				return amp(v);
			}
			return rpwa::piPiSWaveAuMorganPenningtonKachaev::amp(v);
		};

		std::complex<double> default_amp(const rpwa::isobarDecayVertex& v) {
			return rpwa::piPiSWaveAuMorganPenningtonKachaev::amp(v);
		};

		std::string name() const {
			if(bp::override name = this->get_override("name")) {
				return name();
			}
			return rpwa::piPiSWaveAuMorganPenningtonKachaev::name();
		};

		std::string default_name() const {
			return rpwa::piPiSWaveAuMorganPenningtonKachaev::name();
		};

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
		.def("name", &massDependenceWrapper::name, &massDependenceWrapper::default_name)
		.add_static_property("debugMassDependence", &massDependenceWrapper::debug, &massDependenceWrapper::setDebug);

	bp::class_<flatMassDependenceWrapper, bp::bases<rpwa::massDependence> >("flatMassDependence")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &flatMassDependenceWrapper::amp, &flatMassDependenceWrapper::default_amp)
		.def("name", &flatMassDependenceWrapper::name, &flatMassDependenceWrapper::default_name);

	bp::class_<relativisticBreitWignerWrapper, bp::bases<rpwa::massDependence> >("relativisticBreitWigner")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &relativisticBreitWignerWrapper::amp, &relativisticBreitWignerWrapper::default_amp)
		.def("name", &relativisticBreitWignerWrapper::name, &relativisticBreitWignerWrapper::default_name);

	bp::class_<piPiSWaveAuMorganPenningtonMWrapper, bp::bases<rpwa::massDependence> >("piPiSWaveAuMorganPenningtonM")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &piPiSWaveAuMorganPenningtonMWrapper::amp, &piPiSWaveAuMorganPenningtonMWrapper::default_amp)
		.def("name", &piPiSWaveAuMorganPenningtonMWrapper::name, &piPiSWaveAuMorganPenningtonMWrapper::default_name);

	bp::class_<piPiSWaveAuMorganPenningtonVesWrapper, bp::bases<rpwa::massDependence> >("piPiSWaveAuMorganPenningtonVes")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &piPiSWaveAuMorganPenningtonVesWrapper::amp, &piPiSWaveAuMorganPenningtonVesWrapper::default_amp)
		.def("name", &piPiSWaveAuMorganPenningtonVesWrapper::name, &piPiSWaveAuMorganPenningtonVesWrapper::default_name);

	bp::class_<piPiSWaveAuMorganPenningtonKachaevWrapper, bp::bases<rpwa::massDependence> >("piPiSWaveAuMorganPenningtonKachaev")
		.def(bp::self_ns::str(bp::self))
		.def("amp", &piPiSWaveAuMorganPenningtonKachaevWrapper::amp, &piPiSWaveAuMorganPenningtonKachaevWrapper::default_amp)
		.def("name", &piPiSWaveAuMorganPenningtonKachaevWrapper::name, &piPiSWaveAuMorganPenningtonKachaevWrapper::default_name);

	bp::register_ptr_to_python<rpwa::massDependencePtr>();
	bp::register_ptr_to_python<rpwa::flatMassDependencePtr>();
	bp::register_ptr_to_python<rpwa::relativisticBreitWignerPtr>();
	bp::register_ptr_to_python<rpwa::piPiSWaveAuMorganPenningtonMPtr>();
	bp::register_ptr_to_python<rpwa::piPiSWaveAuMorganPenningtonVesPtr>();
	bp::register_ptr_to_python<rpwa::piPiSWaveAuMorganPenningtonKachaevPtr>();

};

