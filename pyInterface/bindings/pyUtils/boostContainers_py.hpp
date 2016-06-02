#ifndef BOOSTCONTAINERS_PY_H
#define BOOSTCONTAINERS_PY_H


#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>

#include <reportingUtils.hpp>


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

	}

}

#endif // BOOSTCONTAINERS_PY_H
