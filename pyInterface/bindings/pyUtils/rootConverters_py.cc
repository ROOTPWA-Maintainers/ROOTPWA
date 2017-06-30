#include "rootConverters_py.h"

#include<TClonesArray.h>
#include<TDirectory.h>
#include<TF1.h>
#include<TFile.h>
#include<TLorentzRotation.h>
#include<TPython.h>
#include<TRandom3.h>
#include<TVector3.h>
#include<TTree.h>

#include<ampIntegralMatrix.h>
#include<amplitudeMetadata.h>
#include<amplitudeTreeLeaf.h>
#include<pwaLikelihood.h>
#include<fitResult.h>

namespace bp = boost::python;


template<typename T>
PyObject* rpwa::py::convertToPy(const T& cxxObj) {
	T* newCxxObj = new T(cxxObj);
	return TPython::ObjectProxy_FromVoidPtr(newCxxObj, newCxxObj->ClassName(), true);
}

template<typename T>
T rpwa::py::convertFromPy(PyObject* pyObj) {
	TObject* TObj = (TObject*)(TPython::ObjectProxy_AsVoidPtr(pyObj));
	T cxxObj = dynamic_cast<T>(TObj);
	return cxxObj;
}

template<typename T>
int rpwa::py::setBranchAddress(T objectPtr, PyObject* pyTree, const std::string& name)
{
		TTree* tree = rpwa::py::convertFromPy<TTree*>(pyTree);
		if(not tree) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for tree when executing rpwa::py::setBranchAddress()");
			bp::throw_error_already_set();
		}
		static std::map<T, T*> pointerMap;
		if(pointerMap.find(objectPtr) == pointerMap.end())
		{
			pointerMap[objectPtr] = new T(objectPtr);
		}
		return tree->SetBranchAddress(name.c_str(), pointerMap[objectPtr]);
}

// explicit template instantiation
template int rpwa::py::setBranchAddress<rpwa::fitResult*>(rpwa::fitResult* objectPtr, PyObject* pyTree, const std::string& name);
template int rpwa::py::setBranchAddress<rpwa::amplitudeTreeLeaf*>(rpwa::amplitudeTreeLeaf* objectPtr, PyObject* pyTree, const std::string& name);

template<typename T>
bool rpwa::py::branch(T objectPtr, PyObject* pyTree, const std::string& name, int bufsize, int splitlevel)
{
		TTree* tree = rpwa::py::convertFromPy<TTree*>(pyTree);
		if(not tree) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for tree when executing rpwa::py::branch()");
			bp::throw_error_already_set();
		}
		static std::map<T, T*> pointerMap;
		if(pointerMap.find(objectPtr) == pointerMap.end())
		{
			pointerMap[objectPtr] = new T(objectPtr);
		}
		TBranch* branch = tree->Branch(name.c_str(), pointerMap[objectPtr], bufsize, splitlevel);
		if (branch == 0) {
			return false;
		} else {
			return true;
		}
}

// explicit template instantiation
template bool rpwa::py::branch<rpwa::fitResult*>(rpwa::fitResult* objectPtr, PyObject* pyTree, const std::string& name, int bufsize, int splitlevel);
template bool rpwa::py::branch<rpwa::amplitudeTreeLeaf*>(rpwa::amplitudeTreeLeaf* objectPtr, PyObject* pyTree, const std::string& name, int bufsize, int splitlevel);

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

	bp::def("__RootConverters_convertToPy_TClonesArray", &rpwa::py::convertToPy<TClonesArray>);
	bp::def(
		"__RootConverters_convertFromPy_TClonesArray", &rpwa::py::convertFromPy<TClonesArray*>
		, bp::return_internal_reference<1>()
	);

	bp::def("__RootConverters_convertToPy_TRandom3", &rpwa::py::convertToPy<TRandom3>);
	bp::def(
		"__RootConverters_convertFromPy_TRandom3", &rpwa::py::convertFromPy<TRandom3*>
		, bp::return_internal_reference<1>()
	);

	bp::def(
		"__RootConverters_convertFromPy_TTree", &rpwa::py::convertFromPy<TTree*>
		, bp::return_internal_reference<1>()
	);

	bp::def("__RootConverters_convertToPy_TMatrixD", &rpwa::py::convertToPy<TMatrixT<double> >);
	bp::def(
		"__RootConverters_convertFromPy_TMatrixD", &rpwa::py::convertFromPy<TMatrixT<double>*>
		, bp::return_internal_reference<1>()
	);

	bp::def("__RootConverters_convertToPy_TVectorD", &rpwa::py::convertToPy<TVectorT<double> >);
	bp::def(
		"__RootConverters_convertFromPy_TVectorD", &rpwa::py::convertFromPy<TVectorT<double>*>
		, bp::return_internal_reference<1>()
	);

	bp::def(
		"__RootConverters_convertFromPy_TFile", &rpwa::py::convertFromPy<TFile*>
		, bp::return_internal_reference<1>()
	);

	bp::def(
		"__RootConverters_convertFromPy_TDirectory", &rpwa::py::convertFromPy<TDirectory*>
		, bp::return_internal_reference<1>()
	);

	bp::def(
		"__RootConverters_convertFromPy_TF1", &rpwa::py::convertFromPy<TF1*>
		, bp::return_internal_reference<1>()
	);

	bp::def(
		"__RootConverters_convertFromPy_rpwaAmpIntegralMatrix", &rpwa::py::convertFromPy<rpwa::ampIntegralMatrix*>
		, bp::return_internal_reference<1>()
	);

	bp::def(
		"__RootConverters_convertFromPy_eventMetadata", &rpwa::py::convertFromPy<rpwa::eventMetadata*>
		, bp::return_internal_reference<1>()
	);

	bp::def(
		"__RootConverters_convertFromPy_rpwaAmplitudeMetadata", &rpwa::py::convertFromPy<rpwa::amplitudeMetadata*>
		, bp::return_internal_reference<1>()
	);

	bp::def(
		"__RootConverters_convertFromPy_rpwaPwaLikelihood", &rpwa::py::convertFromPy<rpwa::pwaLikelihood<std::complex<double> >* >
		, bp::return_internal_reference<1>()
	);

}
