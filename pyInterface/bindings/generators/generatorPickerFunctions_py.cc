#include "generatorPickerFunctions_py.h"

#include <boost/python.hpp>

#include "generatorPickerFunctions.h"

namespace bp = boost::python;


namespace {

	struct massAndTPrimePickerWrapper : public rpwa::massAndTPrimePicker,
	                                           bp::wrapper<rpwa::massAndTPrimePicker>
	{

		massAndTPrimePickerWrapper()
			: rpwa::massAndTPrimePicker(),
			  bp::wrapper<rpwa::massAndTPrimePicker>() { }

		void overrideMassRange(double lowerLimit, double upperLimit) {
			if(bp::override overrideMassRange = this->get_override("overrideMassRange")) {
				overrideMassRange(lowerLimit, upperLimit);
			} else {
				rpwa::massAndTPrimePicker::overrideMassRange(lowerLimit, upperLimit);
			}
		}

		void default_overrideMassRange(double lowerLimit, double upperLimit) {
			rpwa::massAndTPrimePicker::overrideMassRange(lowerLimit, upperLimit);
		}

		bp::tuple massRange__() {
			if(bp::override massRange = this->get_override("massRange")) {
				const std::pair<double, double>& retval = massRange();
				return bp::make_tuple(retval.first, retval.second);
			} else {
				const std::pair<double, double>& retval = rpwa::massAndTPrimePicker::massRange();
				return bp::make_tuple(retval.first, retval.second);
			}
		}

		bp::tuple default_massRange__() {
			const std::pair<double, double>& retval = rpwa::massAndTPrimePicker::massRange();
			return bp::make_tuple(retval.first, retval.second);
		}

		bp::tuple tPrimeRange__() {
			if(bp::override tPrimeRange = this->get_override("tPrimeRange")) {
				const std::pair<double, double>& retval = tPrimeRange();
				return bp::make_tuple(retval.first, retval.second);
			} else {
				const std::pair<double, double>& retval = rpwa::massAndTPrimePicker::tPrimeRange();
				return bp::make_tuple(retval.first, retval.second);
			}
		}

		bp::tuple default_tPrimeRange__() {
			const std::pair<double, double>& retval = rpwa::massAndTPrimePicker::tPrimeRange();
			return bp::make_tuple(retval.first, retval.second);
		}

	};

	bp::tuple massAndTPrimePicker_massRange(const rpwa::massAndTPrimePicker& self) {
		const std::pair<double, double>& retval = self.massRange();
		return bp::make_tuple(retval.first, retval.second);
	}

	bp::tuple massAndTPrimePicker_tPrimeRange(const rpwa::massAndTPrimePicker& self) {
		const std::pair<double, double>& retval = self.tPrimeRange();
		return bp::make_tuple(retval.first, retval.second);
	}

	struct uniformMassExponentialTPickerWrapper : public rpwa::uniformMassExponentialTPicker,
	                                                     bp::wrapper<rpwa::uniformMassExponentialTPicker>
	{

		uniformMassExponentialTPickerWrapper()
			: rpwa::uniformMassExponentialTPicker(),
			  bp::wrapper<rpwa::uniformMassExponentialTPicker>() { }

		uniformMassExponentialTPickerWrapper(const rpwa::uniformMassExponentialTPicker& picker)
			: rpwa::uniformMassExponentialTPicker(picker),
			  bp::wrapper<rpwa::uniformMassExponentialTPicker>() { }

		bool init(const libconfig::Setting& setting) {
			if(bp::override init = this->get_override("init")) {
				return init(setting);
			}
			return rpwa::uniformMassExponentialTPicker::init(setting);
		}

		bool default_init(const libconfig::Setting& setting) {
			return rpwa::uniformMassExponentialTPicker::init(setting);
		}

		bp::tuple call__() {
			bool success;
			double invariantMass;
			double tPrime;
			if(bp::override call = this->get_override("operator()")) {
				success = call(invariantMass, tPrime);
			} else {
				success = rpwa::uniformMassExponentialTPicker::operator()(invariantMass, tPrime);
			}
			return bp::make_tuple(success, invariantMass, tPrime);
		}

		bp::tuple default_call__() {
			bool success;
			double invariantMass;
			double tPrime;
			success = rpwa::uniformMassExponentialTPicker::operator()(invariantMass, tPrime);
			return bp::make_tuple(success, invariantMass, tPrime);
		}

	};

	bp::tuple uniformMassExponentialTPicker_call__(rpwa::uniformMassExponentialTPicker& self) {
		bool success;
		double invariantMass;
		double tPrime;
		success = self(invariantMass, tPrime);
		return bp::make_tuple(success, invariantMass, tPrime);
	}

	struct polynomialMassAndTPrimeSlopePickerWrapper : public rpwa::polynomialMassAndTPrimeSlopePicker,
	                                                          bp::wrapper<rpwa::polynomialMassAndTPrimeSlopePicker>
	{

		polynomialMassAndTPrimeSlopePickerWrapper()
			: rpwa::polynomialMassAndTPrimeSlopePicker(),
			  bp::wrapper<rpwa::polynomialMassAndTPrimeSlopePicker>() { }

		polynomialMassAndTPrimeSlopePickerWrapper(const rpwa::polynomialMassAndTPrimeSlopePicker& picker)
			: rpwa::polynomialMassAndTPrimeSlopePicker(picker),
			  bp::wrapper<rpwa::polynomialMassAndTPrimeSlopePicker>() { }

		bool init(const libconfig::Setting& setting) {
			if(bp::override init = this->get_override("init")) {
				return init(setting);
			}
			return rpwa::polynomialMassAndTPrimeSlopePicker::init(setting);
		}

		bool default_init(const libconfig::Setting& setting) {
			return rpwa::polynomialMassAndTPrimeSlopePicker::init(setting);
		}

		bp::tuple call__() {
			bool success;
			double invariantMass;
			double tPrime;
			if(bp::override call = this->get_override("operator()")) {
				success = call(invariantMass, tPrime);
			} else {
				success = rpwa::polynomialMassAndTPrimeSlopePicker::operator()(invariantMass, tPrime);
			}
			return bp::make_tuple(success, invariantMass, tPrime);
		}

		bp::tuple default_call__() {
			bool success;
			double invariantMass;
			double tPrime;
			success = rpwa::polynomialMassAndTPrimeSlopePicker::operator()(invariantMass, tPrime);
			return bp::make_tuple(success, invariantMass, tPrime);
		}

	};

	bp::tuple polynomialMassAndTPrimeSlopePicker_call__(rpwa::polynomialMassAndTPrimeSlopePicker& self) {
		bool success;
		double invariantMass;
		double tPrime;
		success = self(invariantMass, tPrime);
		return bp::make_tuple(success, invariantMass, tPrime);
	}

}

void rpwa::py::exportGeneratorPickerFunctions() {

	bp::class_<massAndTPrimePickerWrapper, boost::noncopyable>("massAndTPrimePicker", bp::no_init)
		.def("init", bp::pure_virtual(&rpwa::massAndTPrimePicker::init))
		.def(bp::self_ns::str(bp::self))
		.def("overrideMassRange"
			, &massAndTPrimePickerWrapper::overrideMassRange
			, &massAndTPrimePickerWrapper::default_overrideMassRange
		)
		.def("overrideMassRange", &rpwa::massAndTPrimePicker::overrideMassRange)
		.def("massRange", &massAndTPrimePickerWrapper::massRange__, &massAndTPrimePickerWrapper::default_massRange__)
		.def("massRange", &massAndTPrimePicker_massRange)
		.def("tPrimeRange", &massAndTPrimePickerWrapper::tPrimeRange__, &massAndTPrimePickerWrapper::default_tPrimeRange__)
		.def("tPrimeRange", &massAndTPrimePicker_tPrimeRange)
		.def("__call__", bp::pure_virtual(&rpwa::massAndTPrimePicker::operator()));

	bp::class_<uniformMassExponentialTPickerWrapper, bp::bases<rpwa::massAndTPrimePicker> >("uniformMassExponentialTPicker")
		.def(bp::init<const rpwa::uniformMassExponentialTPicker&>())
		.def("init", &uniformMassExponentialTPickerWrapper::init, &uniformMassExponentialTPickerWrapper::default_init)
		.def("init", &rpwa::uniformMassExponentialTPicker::init)
		.def(
			"__call__"
			, &uniformMassExponentialTPickerWrapper::call__
			, &uniformMassExponentialTPickerWrapper::default_call__
		)
		.def("__call__", &uniformMassExponentialTPicker_call__);

	bp::class_<polynomialMassAndTPrimeSlopePickerWrapper, bp::bases<rpwa::massAndTPrimePicker> >("polynomialMassAndTPrimeSlopePicker")
		.def(bp::init<const rpwa::polynomialMassAndTPrimeSlopePicker&>())
		.def("init", &polynomialMassAndTPrimeSlopePickerWrapper::init, &polynomialMassAndTPrimeSlopePickerWrapper::default_init)
		.def("init", &rpwa::polynomialMassAndTPrimeSlopePicker::init)
		.def(
			"__call__"
			, &polynomialMassAndTPrimeSlopePickerWrapper::call__
			, &polynomialMassAndTPrimeSlopePickerWrapper::default_call__
		)
		.def("__call__", &polynomialMassAndTPrimeSlopePicker_call__);

	bp::register_ptr_to_python<massAndTPrimePickerPtr >();
	bp::register_ptr_to_python<boost::shared_ptr<rpwa::uniformMassExponentialTPickerPtr> >();
	bp::register_ptr_to_python<boost::shared_ptr<rpwa::polynomialMassAndTPrimeSlopePickerPtr> >();

}
