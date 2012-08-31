#include<particle_py.h>

#include<TLorentzRotation.h>
#include<TPython.h>
#include<TVector3.h>

#include<rootConverters_py.h>

namespace bp = boost::python;

namespace {

	struct particleWrapper : public rpwa::particle,
	                                bp::wrapper<rpwa::particle>
	{

		particleWrapper()
			: rpwa::particle(),
			  bp::wrapper<rpwa::particle>() { };

		particleWrapper(const rpwa::particle& part)
			: rpwa::particle(part),
			  bp::wrapper<rpwa::particle>() { };

		particleWrapper(const rpwa::particleProperties& partProp,
		                const int                       index    = -1,
		                const int                       spinProj = 0,
		                const int                       refl     = 0,
		                const TVector3&                 momentum = TVector3())
			: rpwa::particle(partProp, index, spinProj, refl, momentum),
			  bp::wrapper<rpwa::particle>() { };

		particleWrapper(const std::string&        partName,
		                const bool                requirePartInTable = true,
		                const int                 index              = -1,
		                const int                 spinProj           = 0,
		                const int                 refl               = 0,
		                const TVector3&           momentum           = TVector3())
			: rpwa::particle(partName, requirePartInTable, index, spinProj, refl, momentum),
			  bp::wrapper<rpwa::particle>() { };

		particleWrapper(const std::string&        partName,
		                const int                 isospin,
		                const int                 G,
		                const int                 J,
		                const int                 P,
		                const int                 C,
		                const int                 spinProj,
		                const int                 refl  = 0,
		                const int                 index = -1)
			: rpwa::particle(partName, isospin, G, J, P, C, spinProj, refl, index),
			  bp::wrapper<rpwa::particle>() { };

		bool read__(bp::object& pyLine) {
			std::string strLine = bp::extract<std::string>(pyLine);
			std::istringstream sstrLine(strLine, std::istringstream::in);
			return rpwa::particleProperties::read(sstrLine);
		};

		PyObject* momentum() const {
			return rpwa::py::convertToPy<TVector3>(particle::momentum());
		}

		PyObject* transform__(PyObject* pyTrafo) {
			TVector3* trafoTV3 = rpwa::py::convertFromPy<TVector3*>(pyTrafo);
			if(trafoTV3 != NULL) {
				return rpwa::py::convertToPy<TLorentzVector>(rpwa::particle::transform(*trafoTV3));
			}
			TLorentzRotation* trafoTLR = rpwa::py::convertFromPy<TLorentzRotation*>(pyTrafo);
			if(trafoTLR != NULL) {
				return rpwa::py::convertToPy<TLorentzVector>(rpwa::particle::transform(*trafoTLR));
			}
			return NULL;
		};

		std::string qnSummary() const {
			if(bp::override qnSummary = this->get_override("qnSummary")) {
				return qnSummary();
			}
			return rpwa::particle::qnSummary();
		}

		std::string default_qnSummary() const {
			return rpwa::particle::qnSummary();
		}

		std::string label() const {
			if(bp::override label = this->get_override("label")) {
				return label();
			}
			return rpwa::particle::label();
		};

		std::string default_label() const {
			return rpwa::particle::label();
		};

	};

}

void rpwa::py::exportParticle() {

	bp::class_<particleWrapper, bp::bases<rpwa::particleProperties> >("particle")

		.def(bp::self_ns::str(bp::self))

		.add_property("spinProj", &particleWrapper::spinProj, &particleWrapper::setSpinProj)
		.def("momentum", &particleWrapper::momentum)

		.def(
	  		"transform"
  			, &particleWrapper::transform__
	  	)

		.def("qnSummary", &particleWrapper::qnSummary, &particleWrapper::default_qnSummary)

		.def("read", &particleWrapper::read__)

		.def("label", &particleWrapper::label, &particleWrapper::default_label);

};

