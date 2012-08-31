#include<rootConverters_py.h>

#include<TLorentzRotation.h>
#include<TPython.h>
#include<TVector3.h>

namespace bp = boost::python;

template<typename T>
PyObject* rpwa::py::convertToPy(const T& cxxObj) {
	T* newCxxObj = new T(cxxObj);
	return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName());
};

template<typename T>
T rpwa::py::convertFromPy(PyObject* pyObj) {
	TObject* TObj = (TObject*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
	T cxxObj = dynamic_cast<T>(TObj);
	return cxxObj;
};

void rpwa::py::exportRootConverters() {

	bp::def("__RootConverters_convertToPy_TVector3", &rpwa::py::convertToPy<TVector3>);
	bp::def(
		"__RootConverters_convertFromPy_TVector3", &rpwa::py::convertFromPy<TVector3*>
		, bp::return_internal_reference<1>()
	);

	bp::def("__RootConverters_convertToPy_TLorentzRotation", &rpwa::py::convertToPy<TLorentzRotation>);
	bp::def(
		"__RootConverters_convertFromPy_TLorentzRotation", &rpwa::py::convertFromPy<TLorentzRotation*>
		, bp::return_internal_reference<1>()
	);

	bp::def("__RootConverters_convertToPy_TLorentzVector", &rpwa::py::convertToPy<TLorentzVector>);
	bp::def(
		"__RootConverters_convertFromPy_TLorentzVector", &rpwa::py::convertFromPy<TLorentzVector*>
		, bp::return_internal_reference<1>()
	);

}
