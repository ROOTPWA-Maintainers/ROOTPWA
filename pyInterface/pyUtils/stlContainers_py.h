#ifndef STLCONTAINERS_PY_H
#define STLCONTAINERS_PY_H

#include "boost/python.hpp"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include<map>
#include<set>
#include<utility>
#include<vector>

#include<particle.h>
#include<isobarDecayVertex.h>
#include<interactionVertex.h>

namespace rpwa {
	namespace py {

		void exportStlContainers();

		template<typename T>
		bool convertBPObjectToSet(const boost::python::object& pyList, std::set<T>& set) {
			boost::python::extract<boost::python::list> getList(pyList);
			if(not getList.check()) {
				printWarn<<"Cannot convert boost::python::object to list."<<std::endl;
				return false;
			}
			boost::python::list pyListList = getList();
			for(int i = 0; i < boost::python::len(pyListList); ++i) {
				boost::python::extract<T> getListItem(pyListList[i]);
				if(not getListItem.check()) {
					printWarn<<"Cannot convert list item."<<std::endl;
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
				printWarn<<"Cannot convert boost::python::object to list."<<std::endl;
				return false;
			}
			boost::python::list pyListList = getList();
			for(int i = 0; i < boost::python::len(pyListList); ++i) {
				boost::python::extract<T> getListItem(pyListList[i]);
				if(not getListItem.check()) {
					printWarn<<"Cannot convert list item."<<std::endl;
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
				printWarn<<"Cannot convert boost::python::object to list."<<std::endl;
				return false;
			}
			boost::python::list pyListList = getList();
			vector.resize(boost::python::len(pyListList));
			if(boost::python::len(pyListList) != 0) {
				for(unsigned int i = 0; i < boost::python::len(pyListList); ++i) {
					boost::python::extract<T> getListItem(pyListList[i]);
					if(not getListItem.check()) {
						printWarn<<"Cannot convert list item."<<std::endl;
						return false;
					}
					vector[i] = getListItem();
				}
			}
			return true;
		}

	}

}

#endif
