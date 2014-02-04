#include "particle_py.h"

#include<TLorentzRotation.h>
#include<TPython.h>
#include<TVector3.h>

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	struct particleWrapper : public rpwa::particle,
	                                bp::wrapper<rpwa::particle>
	{

		particleWrapper()
			: rpwa::particle(),
			  bp::wrapper<rpwa::particle>() { }

		particleWrapper(const rpwa::particle& part)
			: rpwa::particle(part),
			  bp::wrapper<rpwa::particle>() { }

		particleWrapper(const rpwa::particleProperties& partProp,
		                const int                       index    = -1,
		                const int                       spinProj = 0,
		                const int                       refl     = 0,
		                const bp::object&               momentum = bp::object())
			: rpwa::particle(partProp, index, spinProj, refl),
			  bp::wrapper<rpwa::particle>()
		{
			if(!(momentum.is_none())) {
				rpwa::particle::setMomentum(*(rpwa::py::convertFromPy<TVector3*>(momentum.ptr())));
			}
		}

		particleWrapper(const std::string&        partName,
		                const bool                requirePartInTable = true,
		                const int                 index              = -1,
		                const int                 spinProj           = 0,
		                const int                 refl               = 0,
		                const bp::object&         momentum           = bp::object())
			: rpwa::particle(partName, requirePartInTable, index, spinProj, refl),
			  bp::wrapper<rpwa::particle>()
		{
			if(!(momentum.is_none())) {
				rpwa::particle::setMomentum(*(rpwa::py::convertFromPy<TVector3*>(momentum.ptr())));
			}
		}

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
			  bp::wrapper<rpwa::particle>() { }

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
		}

		std::string default_label() const {
			return rpwa::particle::label();
		}

	};

	bool particle_equal(const rpwa::particle& self, const bp::object& rhsObj) {
		bp::extract<rpwa::particleProperties> get_partProp(rhsObj);
		if(get_partProp.check()) {
			return (self == get_partProp());
		}
		bp::tuple rhs = bp::extract<bp::tuple>(rhsObj);
		rpwa::particleProperties rhsProp = bp::extract<rpwa::particleProperties>(rhs[0]);
		std::string rhsString = bp::extract<std::string>(rhs[1]);
		std::pair<rpwa::particleProperties, std::string> rhsPair;
		rhsPair.first = rhsProp;
		rhsPair.second = rhsString;
		return (self == rhsPair);
	}

	bool particle_nequal(const rpwa::particle& self, const bp::object& rhsObj) {
		return not (self == rhsObj);
	}

	bool particle_read(rpwa::particle& self, bp::object& pyLine) {
		std::string strLine = bp::extract<std::string>(pyLine);
		std::istringstream sstrLine(strLine, std::istringstream::in);
		return self.read(sstrLine);
	}

	PyObject* particle_momentum(const rpwa::particle& self) {
		return rpwa::py::convertToPy<TVector3>(self.momentum());
	}

	void particle_setMomentum(rpwa::particle& self, PyObject* pyMom) {
		TVector3* newMom = rpwa::py::convertFromPy<TVector3*>(pyMom);
		if(newMom == NULL) {
			printErr<<"Got invalid input when executing rpwa::particle::setMomentum()."<<std::endl;
		} else {
			self.setMomentum(*newMom);
		}
	}

	PyObject* particle_lzVec(const rpwa::particle& self) {
		return rpwa::py::convertToPy<TLorentzVector>(self.lzVec());
	}

	void particle_setLzVec(rpwa::particle& self, PyObject* pyLzVec) {
		TLorentzVector* newLzVec = rpwa::py::convertFromPy<TLorentzVector*>(pyLzVec);
		if(newLzVec == NULL) {
			printErr<<"Got invalid input when executing rpwa::particle::setLzVec()."<<std::endl;
		} else {
			self.setLzVec(*newLzVec);
		}
	}

	PyObject* particle_transform(rpwa::particle& self, PyObject* pyTrafo) {
		TVector3* trafoTV3 = rpwa::py::convertFromPy<TVector3*>(pyTrafo);
		if(trafoTV3 != NULL) {
			return rpwa::py::convertToPy<TLorentzVector>(self.transform(*trafoTV3));
		}
		TLorentzRotation* trafoTLR = rpwa::py::convertFromPy<TLorentzRotation*>(pyTrafo);
		if(trafoTLR != NULL) {
			return rpwa::py::convertToPy<TLorentzVector>(self.transform(*trafoTLR));
		}
		printErr<<"Got invalid input when executing rpwa::particle::transform()."<<std::endl;
		return NULL;
	}

}

void rpwa::py::exportParticle() {

	bp::class_<particleWrapper, bp::bases<rpwa::particleProperties> >("particle")

		.def(bp::init<const rpwa::particle&>())
		.def(bp::init<const rpwa::particleProperties&, bp::optional<int, int, int, bp::object&> >())
		.def(bp::init<const std::string&, bp::optional<bool, int, int, int, const bp::object&> >())
		.def(bp::init<const std::string&, int, int, int, int, int, int, bp::optional<int, int> >())

		.def(bp::self_ns::str(bp::self))

		.def("__eq__", &particle_equal)
		.def("__neq__", &particle_nequal)

		.def("clone", &rpwa::particle::clone)

		.add_property("spinProj", &rpwa::particle::spinProj, &rpwa::particle::setSpinProj)
		.add_property("momentum", &particle_momentum, &particle_setMomentum)
		.add_property("lzVec", &particle_lzVec, &particle_setLzVec)
		.add_property("index", &rpwa::particle::index, &rpwa::particle::setIndex)
		.add_property("reflectivity", &rpwa::particle::reflectivity, &rpwa::particle::setReflectivity)

		.def("setProperties", &rpwa::particle::setProperties)

		.def("transform", &particle_transform)

		.def("qnSummary", &particleWrapper::qnSummary, &particleWrapper::default_qnSummary)
		.def("qnSummary", &rpwa::particle::qnSummary)

		.def("read", &particle_read)

		.def("label", &particleWrapper::label, &particleWrapper::default_label)
		.def("label", &rpwa::particle::label)

		.add_static_property("debugParticle", &rpwa::particle::debug, &rpwa::particle::setDebug);

	bp::register_ptr_to_python< rpwa::particlePtr >();

}
