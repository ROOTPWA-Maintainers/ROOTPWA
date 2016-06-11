#ifndef ROOTCONVERTERS_PY_H
#define ROOTCONVERTERS_PY_H

#include <boost/python.hpp>

#include <string>


namespace rpwa {

	namespace py {

		template<typename T>
		PyObject* convertToPy(const T& cxxObj);

		template<typename T>
		T convertFromPy(PyObject* pyObj);

		template<typename T>
		int setBranchAddress(T objectPtr, PyObject* pyTree, const std::string& name);

		template<typename T>
		bool branch(T objectPtr,
		            PyObject* pyTree,
		            const std::string& name,
		            int bufsize = 32000,
		            int splitlevel = 99);

		void exportRootConverters();

	}

}

#endif
