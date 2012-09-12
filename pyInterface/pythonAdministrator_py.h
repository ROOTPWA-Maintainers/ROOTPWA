#ifndef PYTHONADMINISTRATOR_PY_H
#define PYTHONADMINISTRATOR_PY_H

#include<boost/python.hpp>

#include<TTree.h>

#include<amplitudeTreeLeaf.h>
#include<waveDescription.h>

#include "rootConverters_py.h"

namespace rpwa {

	namespace py {

		struct pythonAdministrator {

			rpwa::waveDescription _waveDescription;
			rpwa::isobarAmplitudePtr _amplitude;

			pythonAdministrator()
				: _waveDescription(rpwa::waveDescription()) { };

			bool constructAmplitude(std::string keyFileName);

			bool initKinematicsData(PyObject* pyProdKinParticles,
			                        PyObject* pyDecayKinParticles);

			bool readKinematicsData(PyObject* pyProdKinMomenta,
			                        PyObject* pyDecayKinMomenta);

			std::complex<double> operator()() const {
				return _amplitude->amplitude();
			}

			void reset() {
				_waveDescription.clear();
				_amplitude = rpwa::isobarAmplitudePtr();
			};

			static int setBranchAddress(PyObject* pyTree, PyObject* pyAmplitudeTreeLeaf, std::string branchName) {
				rpwa::amplitudeTreeLeaf& amplitudeTreeLeaf = boost::python::extract<rpwa::amplitudeTreeLeaf&>(pyAmplitudeTreeLeaf);
				TTree* tree = rpwa::py::convertFromPy<TTree*>(pyTree);
				return tree->SetBranchAddress(branchName.c_str(), &amplitudeTreeLeaf);
			};

			static void branch(PyObject* pyTree, PyObject* pyAmplitudeTreeLeaf, std::string name, int bufsize = 32000, int splitlevel = 99) {
				rpwa::amplitudeTreeLeaf& amplitudeTreeLeaf = boost::python::extract<rpwa::amplitudeTreeLeaf&>(pyAmplitudeTreeLeaf);
				TTree* tree = rpwa::py::convertFromPy<TTree*>(pyTree);
				tree->Branch(name.c_str(), &amplitudeTreeLeaf, bufsize, splitlevel);
			}

		};

		inline std::ostream& operator << (std::ostream& out,
		                                  const rpwa::py::pythonAdministrator& pyAdmin)
		{
			return pyAdmin._amplitude->print(out);
		};

		void exportPythonAdministrator();

	}

}

#endif

