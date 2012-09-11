#ifndef PYTHONADMINISTRATOR_PY_H
#define PYTHONADMINISTRATOR_PY_H

#include<boost/python.hpp>

#include<waveDescription.h>

#include "rootConverters_py.h"
#include<TClonesArray.h>

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

		};

		void exportPythonAdministrator();

	}

}

#endif

