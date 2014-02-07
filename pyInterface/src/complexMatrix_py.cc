#include "complexMatrix_py.h"

#include "rootConverters_py.h"

namespace bp = boost::python;

namespace {

	struct complexMatrixWrapper : rpwa::complexMatrix,
	                              bp::wrapper<rpwa::complexMatrix>
	{

		complexMatrixWrapper()
			: rpwa::complexMatrix(),
			  bp::wrapper<rpwa::complexMatrix> () { }

		complexMatrixWrapper(const rpwa::complexMatrix& cMatrix)
			: rpwa::complexMatrix(cMatrix),
			  bp::wrapper<rpwa::complexMatrix> () { }

		complexMatrixWrapper(PyObject* pyRe, PyObject* pyIm)
			: rpwa::complexMatrix(),
			  bp::wrapper<rpwa::complexMatrix> ()
		{
			TMatrixD* re = rpwa::py::convertFromPy<TMatrixD*>(pyRe);
			TMatrixD* im = rpwa::py::convertFromPy<TMatrixD*>(pyIm);
			if((not re) or (not im)) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for im/re when executing rpwa::complexMatrix()");
				bp::throw_error_already_set();
			}
			(*this) = rpwa::complexMatrix(*re, *im);
		}

		complexMatrixWrapper(const int i, const int j)
			: rpwa::complexMatrix(i, j),
			  bp::wrapper<rpwa::complexMatrix> () { }

		void Print__(const char* option = "") {
			if(bp::override Print = this->get_override("Print")) {
				Print(option);
			} else {
				rpwa::complexMatrix::Print(option);
			}
		}

		void default_Print__(const char* option = "") {
			rpwa::complexMatrix::Print(option);
		}

	};

	void complexMatrix_set(rpwa::complexMatrix& self,
	                       const int i,
	                       const int j,
	                       const std::complex<double>& c)
	{
		self.set(i, j, c);
	}

	std::complex<double> complexMatrix_get(const rpwa::complexMatrix& self, const int i, const int j)
	{
		TComplex retval = self.get(i, j);
		return std::complex<double>(retval.Re(), retval.Im());
	}

	std::complex<double> complexMatrix__call__(const rpwa::complexMatrix& self, const int i, const int j)
	{
		TComplex retval = self(i, j);
		return std::complex<double>(retval.Re(), retval.Im());
	}

}

void rpwa::py::exportComplexMatrix() {

	bp::class_<complexMatrixWrapper>("complexMatrix")
		.def(bp::init<const rpwa::complexMatrix&>())
		.def(bp::init<PyObject*, PyObject*>())
		.def(bp::init<const int, const int>())
		.def(bp::self_ns::str(bp::self))
		.def(bp::self * bp::self)
		.def(bp::self - bp::self)
		.def(bp::self + bp::self)
		.def("resizeTo", &rpwa::complexMatrix::resizeTo)
		.def("set", &complexMatrix_set)
		.def("get", &complexMatrix_get)
		.def("__call__", complexMatrix__call__)
		.def("nRows", &rpwa::complexMatrix::nRows)
		.def("nCols", &rpwa::complexMatrix::nCols)
		.def("Print", &complexMatrixWrapper::Print__, &complexMatrixWrapper::default_Print__)
		.def("Print", &rpwa::complexMatrix::Print)
		.def("t", &rpwa::complexMatrix::t)
		.def("dagger", &rpwa::complexMatrix::dagger);

}
