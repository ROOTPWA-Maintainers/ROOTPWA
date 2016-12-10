#ifndef BOOSTCONTAINERS_PY_H
#define BOOSTCONTAINERS_PY_H


#include <boost/multi_array.hpp>
#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>

#include <TMatrixT.h>
#include <TPython.h>

// set up numpy for usage in this file
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL RPWA_PyArray_API
#include <numpy/arrayobject.h>

#include <reportingUtils.hpp>
#include <rootConverters_py.h>
#include <stlContainers_py.h>


namespace rpwa {

	namespace py {

		template<typename T1>
		bool convertBPTupleToTuple(const boost::python::tuple& pyTuple, boost::tuples::tuple<T1>& tuple, const unsigned int index = 0) {
			if (boost::python::len(pyTuple) != 1) {
				printWarn << "tuple has wrong length, expected 1, got " << boost::python::len(pyTuple) << "." << std::endl;
				return false;
			}

			boost::python::extract<T1> element(pyTuple[0]);
			if (not element.check()) {
				printWarn << "cannot convert tuple element at index " << index << " to desired type." << std::endl;
				return false;
			}
			tuple.head = element();

			return true;
		}


		template<typename T1, typename T2>
		bool convertBPTupleToTuple(const boost::python::tuple& pyTuple, boost::tuples::tuple<T1, T2>& tuple, const unsigned int index = 0) {
			if (boost::python::len(pyTuple) != 2) {
				printWarn << "tuple has wrong length, expected 2, got " << boost::python::len(pyTuple) << "." << std::endl;
				return false;
			}

			boost::python::extract<T1> element(pyTuple[0]);
			if (not element.check()) {
				printWarn << "cannot convert tuple element at index " << index << " to desired type." << std::endl;
				return false;
			}
			tuple.head = element();

			boost::tuples::tuple<T2> tail;
			if (not convertBPTupleToTuple<T2>(boost::python::tuple(pyTuple.slice(1, boost::python::len(pyTuple))), tail, index+1)) {
				return false;
			}
			tuple.tail = tail;

			return true;
		}


		template<typename T1, typename T2, typename T3>
		bool convertBPTupleToTuple(const boost::python::tuple& pyTuple, boost::tuples::tuple<T1, T2, T3>& tuple, const unsigned int index = 0) {
			if (boost::python::len(pyTuple) != 3) {
				printWarn << "tuple has wrong length, expected 3, got " << boost::python::len(pyTuple) << "." << std::endl;
				return false;
			}

			boost::python::extract<T1> element(pyTuple[0]);
			if (not element.check()) {
				printWarn << "cannot convert tuple element at index " << index << " to desired type." << std::endl;
				return false;
			}
			tuple.head = element();

			boost::tuples::tuple<T2, T3> tail;
			if (not convertBPTupleToTuple<T2, T3>(boost::python::tuple(pyTuple.slice(1, boost::python::len(pyTuple))), tail, index+1)) {
				return false;
			}
			tuple.tail = tail;

			return true;
		}

		template<typename T, typename npyT, const int typeenum, const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<T, dim>& array) {
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_FromObject(pyArray.ptr(), typeenum, dim, dim);
			if(tempArray == NULL) {
				return false;
			}

			array.resize(std::vector<npy_intp>(PyArray_SHAPE(tempArray), PyArray_SHAPE(tempArray) + dim));
			array.assign((npyT*) PyArray_DATA(tempArray), (npyT*) PyArray_DATA(tempArray) + PyArray_SIZE(tempArray));

			return true;
		}

		template<const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<double, dim>& array) {
			return rpwa::py::convertBPObjectToMultiArray<double, npy_double, NPY_DOUBLE, dim>(pyArray, array);
		}

		template<typename T, typename npyT, const int typeenum, const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<std::complex<T>, dim>& array) {
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_FromObject(pyArray.ptr(), typeenum, dim, dim);
			if(tempArray == NULL) {
				return false;
			}

			array.resize(std::vector<npy_intp>(PyArray_SHAPE(tempArray), PyArray_SHAPE(tempArray) + dim));
			for(npy_intp idx = 0; idx < PyArray_SIZE(tempArray); ++idx) {
				const npyT element = ((npyT*) PyArray_DATA(tempArray))[idx];
				array.data()[idx] = std::complex<T>(element.real, element.imag);
			}

			return true;
		}

		template<const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<std::complex<double>, dim>& array) {
			return rpwa::py::convertBPObjectToMultiArray<double, npy_cdouble, NPY_CDOUBLE, dim>(pyArray, array);
		}

		template<typename T, typename npyT, const int typeenum, const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<std::pair<T, T>, dim>& array) {
			// there are tow possibilities the pairs are stored:
			// 1. in a dim+1 dimensional array where the last two
			//    dimensions correspond to one pair
			// 2. in a dim dimensitonal array where each element
			//    is a Python tuple

			PyArrayObject* tempArray = (PyArrayObject*) PyArray_FromObject(pyArray.ptr(), typeenum, dim + 1, dim + 1);
			if(tempArray != NULL) {
				// first case
				if(PyArray_SHAPE(tempArray)[dim - 1] != 2) {
					return false;
				}

				array.resize(std::vector<npy_intp>(PyArray_SHAPE(tempArray), PyArray_SHAPE(tempArray) + dim));
				for(npy_intp idx = 0; idx < PyArray_SIZE(tempArray); idx += 2) {
					const T first = ((npyT*) PyArray_DATA(tempArray))[idx];
					const T second = ((npyT*) PyArray_DATA(tempArray))[idx + 1];

					array.data()[idx / 2] = std::make_pair(first, second);
				}

				return true;
			} else {
				PyErr_Clear();
			}

			tempArray = (PyArrayObject*) PyArray_FromObject(pyArray.ptr(), NPY_OBJECT, dim, dim);
			if(tempArray != NULL) {
				// second case
				array.resize(std::vector<npy_intp>(PyArray_SHAPE(tempArray), PyArray_SHAPE(tempArray) + dim));
				for(npy_intp idx = 0; idx < PyArray_SIZE(tempArray); ++idx) {
					boost::python::object element(boost::python::handle<>(boost::python::borrowed(((PyObject**) PyArray_DATA(tempArray))[idx])));

					if(not rpwa::py::convertBPObjectToPair(element, array.data()[idx])) {
						return false;
					}
				}

				return true;
			} else {
				PyErr_Clear();
			}

			return false;
		}

		template<const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<std::pair<double, double>, dim>& array) {
			return convertBPObjectToMultiArray<double, npy_double, NPY_DOUBLE, dim>(pyArray, array);
		}

		template<const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<std::pair<size_t, size_t>, dim>& array) {
			return convertBPObjectToMultiArray<size_t, npy_ulong, NPY_ULONG, dim>(pyArray, array);
		}

		template<const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<std::string, dim>& array) {
			// there are two possibilities the strings are stored:
			// 1. as arrays of strings with fixed size padded with null bytes
			// 2. as objects

			PyArrayObject* tempArray = (PyArrayObject*) PyArray_FromObject(pyArray.ptr(), NPY_STRING, dim, dim);
			if(tempArray != NULL) {
				// first case
				array.resize(std::vector<npy_intp>(PyArray_SHAPE(tempArray), PyArray_SHAPE(tempArray) + dim));
				for(npy_intp idx = 0; idx < PyArray_SIZE(tempArray); ++idx) {
					std::string element((char*) PyArray_DATA(tempArray) + idx * PyArray_ITEMSIZE(tempArray), PyArray_ITEMSIZE(tempArray));
					// all elements in the Python array have the
					// same size and seem to be padded with null
					// bytes, remove those
					while(element.back() == '\0') {
						element.pop_back();
					}

					array.data()[idx] = element;
				}

				return true;
			} else {
				PyErr_Clear();
			}

			tempArray = (PyArrayObject*) PyArray_FromObject(pyArray.ptr(), NPY_OBJECT, dim, dim);
			if(tempArray != NULL) {
				// second case
				array.resize(std::vector<npy_intp>(PyArray_SHAPE(tempArray), PyArray_SHAPE(tempArray) + dim));
				for(npy_intp idx = 0; idx < PyArray_SIZE(tempArray); ++idx) {
					boost::python::object element(boost::python::handle<>(boost::python::borrowed(((PyObject**) PyArray_DATA(tempArray))[idx])));
					array.data()[idx] = boost::python::extract<std::string>(element)();
				}

				return true;
			} else {
				PyErr_Clear();
			}

			return false;
		}

		template<typename T, const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<TMatrixT<T>, dim>& array) {
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_FromObject(pyArray.ptr(), NPY_OBJECT, dim, dim);
			if(tempArray == NULL) {
				return false;
			}

			array.resize(std::vector<npy_intp>(PyArray_SHAPE(tempArray), PyArray_SHAPE(tempArray) + dim));
			for(npy_intp idx = 0; idx < PyArray_SIZE(tempArray); ++idx) {
				TMatrixT<T>* element = rpwa::py::convertFromPy<TMatrixT<T>*>(((PyObject**) PyArray_DATA(tempArray))[idx]);
				array.data()[idx].ResizeTo(*element);
				array.data()[idx] = *element;
			}

			return true;
		}

		template<typename T, const size_t dim>
		bool convertBPObjectToMultiArray(const boost::python::object& pyArray, boost::multi_array<T, dim>& array, const bool ignoreNone = false) {
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_FromObject(pyArray.ptr(), NPY_OBJECT, dim, dim);
			if(tempArray == NULL) {
				return false;
			}

			array.resize(std::vector<npy_intp>(PyArray_SHAPE(tempArray), PyArray_SHAPE(tempArray) + dim));
			for(npy_intp idx = 0; idx < PyArray_SIZE(tempArray); ++idx) {
				boost::python::object element(boost::python::handle<>(boost::python::borrowed(((PyObject**) PyArray_DATA(tempArray))[idx])));

				if(ignoreNone and element.is_none()) {
					continue;
				}

				array.data()[idx] = boost::python::extract<T>(element);
			}

			return true;
		}

		template<typename T, typename npyT, const int typeenum, const size_t dim>
		boost::python::object
		convertMultiArrayToPy(const boost::multi_array<T, dim>& array) {
			std::vector<npy_intp> shape(array.shape(), array.shape() + dim);
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_ZEROS(dim, shape.data(), typeenum, 0);

			for(size_t idx = 0; idx < array.num_elements(); ++idx) {
				((npyT*) PyArray_DATA(tempArray))[idx] = array.data()[idx];
			}

			return boost::python::object(boost::python::handle<>(boost::python::borrowed((PyObject*) tempArray)));
		}

		template<const size_t dim>
		boost::python::object
		convertMultiArrayToPy(const boost::multi_array<double, dim>& array) {
			return convertMultiArrayToPy<double, npy_double, NPY_DOUBLE, dim>(array);
		}

		template<typename T, typename npyT, const int typeenum, const size_t dim>
		boost::python::object
		convertMultiArrayToPy(const boost::multi_array<std::complex<T>, dim>& array) {
			std::vector<npy_intp> shape(array.shape(), array.shape() + dim);
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_ZEROS(dim, shape.data(), typeenum, 0);

			for(size_t idx = 0; idx < array.num_elements(); ++idx) {
				((npyT*) PyArray_DATA(tempArray))[idx].real = array.data()[idx].real();
				((npyT*) PyArray_DATA(tempArray))[idx].imag = array.data()[idx].imag();
			}

			return boost::python::object(boost::python::handle<>(boost::python::borrowed((PyObject*) tempArray)));
		}

		template<const size_t dim>
		boost::python::object
		convertMultiArrayToPy(const boost::multi_array<std::complex<double>, dim>& array) {
			return convertMultiArrayToPy<double, npy_cdouble, NPY_CDOUBLE, dim>(array);
		}

		template<typename T, const size_t dim>
		boost::python::object
		convertMultiArrayToPy(const boost::multi_array<std::pair<T, T>, dim>& array) {
			std::vector<npy_intp> shape(array.shape(), array.shape() + dim);
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_ZEROS(dim, shape.data(), NPY_OBJECT, 0);

			for(size_t idx = 0; idx < array.num_elements(); ++idx) {
				boost::python::tuple element = boost::python::make_tuple(array.data()[idx].first,
				                                                         array.data()[idx].second);

				((PyObject**) PyArray_DATA(tempArray))[idx] = element.ptr();
				PyArray_Item_INCREF((char*) ((PyObject**) PyArray_DATA(tempArray) + idx), PyArray_DESCR(tempArray));
			}

			return boost::python::object(boost::python::handle<>(boost::python::borrowed((PyObject*) tempArray)));
		}

		template<typename T, const size_t dim>
		boost::python::object
		convertMultiArrayToPy(const boost::multi_array<TMatrixT<T>, dim>& array) {
			std::vector<npy_intp> shape(array.shape(), array.shape() + dim);
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_ZEROS(dim, shape.data(), NPY_OBJECT, 0);

			for(size_t idx = 0; idx < array.num_elements(); ++idx) {
				TMatrixT<T>* element = new TMatrixT<T>(array.data()[idx]);
				((PyObject**) PyArray_DATA(tempArray))[idx] = TPython::ObjectProxy_FromVoidPtr(element, element->ClassName(), true);
			}

			return boost::python::object(boost::python::handle<>(boost::python::borrowed((PyObject*) tempArray)));
		}

		template<typename T, const size_t dim>
		boost::python::object
		convertMultiArrayToPy(const boost::multi_array<T, dim>& array) {
			std::vector<npy_intp> shape(array.shape(), array.shape() + dim);
			PyArrayObject* tempArray = (PyArrayObject*) PyArray_ZEROS(dim, shape.data(), NPY_OBJECT, 0);

			for(size_t idx = 0; idx < array.num_elements(); ++idx) {
				boost::python::object element(array.data()[idx]);
				((PyObject**) PyArray_DATA(tempArray))[idx] = element.ptr();
				PyArray_Item_INCREF((char*) ((PyObject**) PyArray_DATA(tempArray) + idx), PyArray_DESCR(tempArray));
			}

			return boost::python::object(boost::python::handle<>(boost::python::borrowed((PyObject*)tempArray)));
		}

	}

}

#endif // BOOSTCONTAINERS_PY_H
