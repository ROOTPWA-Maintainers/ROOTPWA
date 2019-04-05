#ifndef STLCONTAINERS_PY_H
#define STLCONTAINERS_PY_H

#include <map>
#include <set>
#include <utility>
#include <vector>

#include <boost/python.hpp>

#include "reportingUtils.hpp"
#include "multibinTypes.h"


namespace rpwa {

	namespace py {

		void exportStlContainers();

		template<typename T>
		bool convertBPObjectToSet(const boost::python::object& pyList, std::set<T>& set) {
			boost::python::extract<boost::python::list> getList(pyList);
			if(not getList.check()) {
				printWarn<<"cannot convert boost::python::object to list."<<std::endl;
				return false;
			}
			boost::python::list pyListList = getList();
			for(int i = 0; i < boost::python::len(pyListList); ++i) {
				boost::python::extract<T> getListItem(pyListList[i]);
				if(not getListItem.check()) {
					printWarn<<"cannot convert list item."<<std::endl;
					return false;
				}
				set.insert(getListItem());
			}
			return true;
		}

		template<typename T>
		bool convertBPObjectToMultiSet(const boost::python::object& pyList, std::multiset<T>& multiset) {
			boost::python::extract<boost::python::list> getList(pyList);
			if(not getList.check()) {
				printWarn<<"cannot convert boost::python::object to list."<<std::endl;
				return false;
			}
			boost::python::list pyListList = getList();
			for(int i = 0; i < boost::python::len(pyListList); ++i) {
				boost::python::extract<T> getListItem(pyListList[i]);
				if(not getListItem.check()) {
					printWarn<<"cannot convert list item."<<std::endl;
					return false;
				}
				multiset.insert(getListItem());
			}
			return true;
		}

		template<typename T>
		bool convertBPObjectToVector(const boost::python::object& pyList, std::vector<T>& vector)
		{
			boost::python::extract<boost::python::list> getList(pyList);
			if(not getList.check()) {
				printWarn<<"cannot convert boost::python::object to list."<<std::endl;
				return false;
			}
			boost::python::list pyListList = getList();
			vector.resize(boost::python::len(pyListList));
			if(boost::python::len(pyListList) != 0) {
				for(unsigned int i = 0; i < boost::python::len(pyListList); ++i) {
					boost::python::extract<T> getListItem(pyListList[i]);
					if(not getListItem.check()) {
						printWarn<<"cannot convert list item "<<i<<"."<<std::endl;
						return false;
					}
					vector[i] = getListItem();
				}
			}
			return true;
		}

		template<typename T, typename U>
		bool convertBPObjectToPair(const boost::python::object& pyTupleObject, std::pair<T, U>& pair) {
			boost::python::extract<boost::python::tuple> getTuple(pyTupleObject);
			if(not getTuple.check()) {
				printWarn<<"cannot convert boost::python::object to tuple."<<std::endl;
				return false;
			}
			boost::python::tuple pyTuple = getTuple();
			if(boost::python::len(pyTuple) != 2) {
				printWarn<<"tuple has wrong length, expected 2, got "<<boost::python::len(pyTuple)<<"."<<std::endl;
				return false;
			}
			boost::python::extract<T> getFirstElement(pyTuple[0]);
			if(not getFirstElement.check()) {
				printWarn<<"cannot convert first tuple element to desired type."<<std::endl;
				return false;
			}
			pair.first = getFirstElement();
			boost::python::extract<U> getSecondElement(pyTuple[1]);
			if(not getSecondElement.check()) {
				printWarn<<"cannot convert second tuple element to desired type."<<std::endl;
				return false;
			}
			pair.second = getSecondElement();
			return true;
		}

		template<typename T, typename U>
		bool convertBPObjectToMap(const boost::python::object& pyDict, std::map<T,U>& map)
		{
			boost::python::extract<boost::python::dict> getDict(pyDict);
			if(not getDict.check()) {
				printWarn<<"cannot convert boost::python::object to dict."<<std::endl;
				return false;
			}
			boost::python::dict pyDictDict = getDict();
			map.clear();
			boost::python::list pyKeys = pyDictDict.keys();
			for (int i = 0; i < boost::python::len(pyKeys); ++i) {
				boost::python::extract<T> getKey(pyKeys[i]);
				if (not getKey.check()) {
					PyErr_SetString(PyExc_TypeError, "Got invalid key");
					boost::python::throw_error_already_set();
				}
				boost::python::extract<U> getValue(pyDictDict[pyKeys[i]]);
				if (not getValue.check()) {
					PyErr_SetString(PyExc_TypeError, "Got invalid value");
					boost::python::throw_error_already_set();
				}
				map[getKey()] = getValue();
			}
			return true;
		}

		rpwa::multibinBoundariesType convertMultibinBoundariesFromPy(const boost::python::dict& pyMultibinBoundaries);
		boost::python::dict convertMultibinBoundariesToPy(const rpwa::multibinBoundariesType& multibinBoundaries);
		/**
		 * @return a Python Multibin object with the multibinBoundaries as boundaries
		 */
		boost::python::object convertMultibinBoundariesToPyMultibin(const rpwa::multibinBoundariesType& multibinBoundaries);
	}
}

#endif
