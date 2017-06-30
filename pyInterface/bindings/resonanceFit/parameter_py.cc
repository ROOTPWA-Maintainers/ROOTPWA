#include "parameter_py.h"

#include <boost/python.hpp>

#include <resonanceFit/parameter.h>

namespace bp = boost::python;


namespace {


	std::shared_ptr<rpwa::resonanceFit::parameter>
	parameter_constructor(const std::string& name,
	                      const bp::object& startValue,
	                      const bp::object& startError,
	                      const bp::object& fixed,
	                      const bp::object& limitLower,
	                      const bp::object& limitUpper,
	                      const bp::object& step)
	{
		std::shared_ptr<rpwa::resonanceFit::parameter> parameter(new rpwa::resonanceFit::parameter);

		parameter->setName(name);

		if(not startValue.is_none()) {
			parameter->setStartValue(bp::extract<double>(startValue));
		}

		if(not startError.is_none()) {
			parameter->setStartError(bp::extract<double>(startError));
		}

		if(not fixed.is_none()) {
			parameter->setFixed(bp::extract<bool>(fixed));
		}

		if(not limitLower.is_none()) {
			parameter->setLimitLower(bp::extract<double>(limitLower));
			parameter->setLimitedLower(true);
		} else {
			parameter->setLimitedLower(false);
		}

		if(not limitUpper.is_none()) {
			parameter->setLimitUpper(bp::extract<double>(limitUpper));
			parameter->setLimitedUpper(true);
		} else {
			parameter->setLimitedUpper(false);
		}

		if(not step.is_none()) {
			parameter->setStep(bp::extract<double>(step));
		}

		return parameter;
	}


}


void rpwa::py::resonanceFit::exportParameter() {

	// the classes and functions from the 'resonanceFit' namespace should
	// not be available from Python at the 'pyRootPwa.core.[...]' level,
	// but should also be included into their own module like
	// 'pyRootPwa.core.resonanceFit.[...]'. So below a submodule
	// 'resonanceFit' is added and the classes and functions are added at
	// that scope.

	bp::scope current;
	std::string submoduleName(bp::extract<std::string>(current.attr("__name__")));
	submoduleName.append(".resonanceFit");

	bp::object submodule(bp::borrowed(PyImport_AddModule(submoduleName.c_str())));
	current.attr("resonanceFit") = submodule;

	// switch the scope to the submodule
	bp::scope submoduleScope = submodule;

	bp::class_<rpwa::resonanceFit::parameter, std::shared_ptr<rpwa::resonanceFit::parameter> >
		(
			"parameter"
			, bp::no_init
		)

		.def(
			"__init__"
			, bp::make_constructor(parameter_constructor,
			                       bp::default_call_policies(),
			                       (bp::arg("name"),
			                        bp::arg("startValue") = bp::object(),
			                        bp::arg("startError") = bp::object(),
			                        bp::arg("fixed") = bp::object(),
			                        bp::arg("limitLower") = bp::object(),
			                        bp::arg("limitUpper") = bp::object(),
			                        bp::arg("step") = bp::object()))
		)

		.def(
			bp::self_ns::str(bp::self)
		)

		.add_property(
			"name"
			, bp::make_function(&rpwa::resonanceFit::parameter::name,
			                    bp::return_value_policy<bp::copy_const_reference>())
			, &rpwa::resonanceFit::parameter::setName
		)

		.add_property(
			"startValue"
			, &rpwa::resonanceFit::parameter::startValue
			, &rpwa::resonanceFit::parameter::setStartValue
		)

		.add_property(
			"startError"
			, &rpwa::resonanceFit::parameter::startError
			, &rpwa::resonanceFit::parameter::setStartError
		)

		.add_property(
			"fixed"
			, &rpwa::resonanceFit::parameter::fixed
			, &rpwa::resonanceFit::parameter::setFixed
		)

		.add_property(
			"limitLower"
			, &rpwa::resonanceFit::parameter::limitLower
			, &rpwa::resonanceFit::parameter::setLimitLower
		)

		.add_property(
			"limitedLower"
			, &rpwa::resonanceFit::parameter::limitedLower
			, &rpwa::resonanceFit::parameter::setLimitedLower
		)

		.add_property(
			"limitUpper"
			, &rpwa::resonanceFit::parameter::limitUpper
			, &rpwa::resonanceFit::parameter::setLimitUpper
		)

		.add_property(
			"limitedUpper"
			, &rpwa::resonanceFit::parameter::limitedUpper
			, &rpwa::resonanceFit::parameter::setLimitedUpper
		)

		.add_property(
			"step"
			, &rpwa::resonanceFit::parameter::step
			, &rpwa::resonanceFit::parameter::setStep
		)

		;

}
