#include "pwaLikelihood_py.h"

#include <boost/python.hpp>

#include "amplitudeMetadata.h"
#include "boostContainers_py.hpp"
#include "pwaLikelihood.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

#include <complex>

namespace bp = boost::python;


namespace {

	bool
	pwaLikelihood_init(rpwa::pwaLikelihood<std::complex<double> >& self,
	                   const bp::list&                             pyWaveDescriptionThresholds,
	                   const unsigned int                          rank,
	                   const double                                massBinCenter)
	{
		std::vector<boost::python::tuple> vectorPyTuples;
		if(not rpwa::py::convertBPObjectToVector<boost::python::tuple>(pyWaveDescriptionThresholds, vectorPyTuples)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for waveDescriptionThresholds when executing rpwa::pwaFit()");
			bp::throw_error_already_set();
		}

		std::vector<rpwa::pwaLikelihood<std::complex<double> >::waveDescThresType> vectorWaveDescriptionThresholds;
		for(size_t i = 0; i < vectorPyTuples.size(); ++i) {
			boost::tuples::tuple<std::string, rpwa::waveDescription, double> waveDescriptionThresholds;
			if(not rpwa::py::convertBPTupleToTuple<std::string, rpwa::waveDescription, double>(vectorPyTuples[i], waveDescriptionThresholds)) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for waveDescriptionThresholds when executing rpwa::pwaFit()");
				bp::throw_error_already_set();
			}
			vectorWaveDescriptionThresholds.push_back(waveDescriptionThresholds);
		}

		return self.init(vectorWaveDescriptionThresholds, rank, massBinCenter);
	}


	bool
	pwaLikelihood_addNormIntegral(rpwa::pwaLikelihood<std::complex<double> >& self,
	                              PyObject*                                   pyNormMatrix)
	{
		rpwa::ampIntegralMatrix* normMatrix = rpwa::py::convertFromPy<rpwa::ampIntegralMatrix* >(pyNormMatrix);
		if(not normMatrix) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for normMatrix when executing rpwa::pwaLikelihood::addNormIntegral()");
			bp::throw_error_already_set();
		}
		return self.addNormIntegral(*normMatrix);
	}


	bool
	pwaLikelihood_addAmplitude(rpwa::pwaLikelihood<std::complex<double> >& self,
	                           bp::list                                    pyMetas)
	{
		std::vector<const rpwa::amplitudeMetadata*> metas;
		if (not rpwa::py::convertBPObjectToVector<const rpwa::amplitudeMetadata*>(pyMetas, metas)){
			PyErr_SetString(PyExc_TypeError, "could not extract vector of amplitude metadata");
			bp::throw_error_already_set();
		}
		return self.addAmplitude(metas);
	}


	bool
	pwaLikelihood_addAccIntegral(rpwa::pwaLikelihood<std::complex<double> >& self,
	                             PyObject*                                   pyAccMatrix,
	                             const unsigned int                          accEventsOverride)
	{
		rpwa::ampIntegralMatrix* accMatrix = rpwa::py::convertFromPy<rpwa::ampIntegralMatrix* >(pyAccMatrix);
		if(not accMatrix) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for accMatrix when executing rpwa::pwaLikelihood::addAccIntegral()");
			bp::throw_error_already_set();
		}
		return self.addAccIntegral(*accMatrix, accEventsOverride);
	}


	bp::list
	pwaLikelihood_Gradient(rpwa::pwaLikelihood<std::complex<double> >& self,
	                       const bp::list&                             pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		std::vector<double> gradient(nmbPar, 0.);
		self.Gradient(par.data(), gradient.data());
		return bp::list(gradient);
	}


	bp::tuple
	pwaLikelihood_FdF(rpwa::pwaLikelihood<std::complex<double> >& self,
	                  const bp::list&                             pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		std::vector<double> gradient(nmbPar, 0);
		double funcVal = 0.;
		self.FdF(par.data(), funcVal, gradient.data());
		return bp::make_tuple(funcVal, bp::list(gradient));
	}


	double
	pwaLikelihood_DoEval(rpwa::pwaLikelihood<std::complex<double> >& self,
	                     const bp::list&                             pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return self.DoEval(par.data());
	}


	double
	pwaLikelihood_DoDerivative(rpwa::pwaLikelihood<std::complex<double> >& self,
	                           const bp::list&                             pyPar,
	                           const unsigned int                          derivIndex)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return self.DoDerivative(par.data(), derivIndex);
	}


	PyObject*
	pwaLikelihood_Hessian(rpwa::pwaLikelihood<std::complex<double> >& self,
	                      const bp::list&                             pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return rpwa::py::convertToPy<TMatrixT<double> >(self.Hessian(par.data()));
	}


	bp::list
	pwaLikelihood_HessianEigenVectors(rpwa::pwaLikelihood<std::complex<double> >& self,
	                                  PyObject*                                   pyHessian)
	{
		TMatrixT<double>* hessian = rpwa::py::convertFromPy<TMatrixT<double>*>(pyHessian);
		if(not hessian) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for hessian when executing rpwa::pwaLikelihood::HessianEigenVectors()");
			bp::throw_error_already_set();
		}
		std::vector<std::pair<TVectorT<double>, double> > eigenVectors = self.HessianEigenVectors(*hessian);
		bp::list pyEigenVectors;
		for(std::vector<std::pair<TVectorT<double>, double> >::const_iterator it = eigenVectors.begin(); it != eigenVectors.end(); ++it) {
			pyEigenVectors.append(bp::make_tuple(boost::python::handle<>(rpwa::py::convertToPy(it->first)), it->second));
		}
		return pyEigenVectors;
	}


	PyObject*
	pwaLikelihood_CovarianceMatrixFromPar(rpwa::pwaLikelihood<std::complex<double> >& self,
	                                      const bp::list&                             pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return rpwa::py::convertToPy<TMatrixT<double> >(self.CovarianceMatrix(par.data()));
	}


	PyObject*
	pwaLikelihood_CovarianceMatrixFromMatrix(rpwa::pwaLikelihood<std::complex<double> >& self,
	                                         PyObject*                                   pyHessian)
	{
		TMatrixT<double>* hessian = rpwa::py::convertFromPy<TMatrixT<double>*>(pyHessian);
		if(not hessian) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for hessian when executing rpwa::pwaLikelihood::CovarianceMatrix()");
			bp::throw_error_already_set();
		}
		return rpwa::py::convertToPy<TMatrixT<double> >(self.CovarianceMatrix(*hessian));
	}


	bp::list
	pwaLikelihood_CorrectParamSigns(rpwa::pwaLikelihood<std::complex<double> >& self,
	                                const bp::list&                             pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return bp::list(self.CorrectParamSigns(par.data()));
	}


	bool
	pwaLikelihood_setOnTheFlyBinning(rpwa::pwaLikelihood<std::complex<double> >& self,
	                                 bp::dict                                    pyBinningMap,
	                                 bp::list                                    pyEvtMetas)
	{
		std::map<std::string, std::pair<double, double> > binningMap;
		const bp::list keys = pyBinningMap.keys();
		for(unsigned int i = 0; i < bp::len(keys); i++){
			std::string binningVar = bp::extract<std::string>(keys[i]);
			double lowerBound      = bp::extract<double>(pyBinningMap[binningVar][0]);
			double upperBound      = bp::extract<double>(pyBinningMap[binningVar][1]);
			binningMap.insert(std::pair<std::string, std::pair<double, double> >(binningVar, std::pair<double, double>(lowerBound, upperBound)));
		}
		std::vector<const rpwa::eventMetadata*> evtMetas;
		if (not rpwa::py::convertBPObjectToVector<const rpwa::eventMetadata*>(pyEvtMetas, evtMetas)){
			PyErr_SetString(PyExc_TypeError, "could not extract event metadatas");
			bp::throw_error_already_set();
		}
		return self.setOnTheFlyBinning(binningMap, evtMetas);
	}

}


void rpwa::py::exportPwaLikelihood() {

	bp::scope theScope = bp::class_<rpwa::pwaLikelihood<std::complex<double> > >("pwaLikelihood")
		.def(
			"init"
			, ::pwaLikelihood_init
			, (bp::arg("waveDescriptionThresholds"),
			   bp::arg("rank") = 1,
			   bp::arg("massBinCenter") = 0.)
		)
		.def("addNormIntegral", ::pwaLikelihood_addNormIntegral)
		.def("addNormIntegral", &rpwa::pwaLikelihood<std::complex<double> >::addNormIntegral)
		.def(
			"addAccIntegral"
			, ::pwaLikelihood_addAccIntegral
			, (bp::arg("accMatrix"),
			   bp::arg("accEventsOverride") = 0)
		)
		.def(
			"addAccIntegral"
			, &rpwa::pwaLikelihood<std::complex<double> >::addAccIntegral
			, (bp::arg("accMatrix"),
			   bp::arg("accEventsOverride") = 0)
		)
		.def("addAmplitude", ::pwaLikelihood_addAmplitude)
		.def("setOnTheFlyBinning", ::pwaLikelihood_setOnTheFlyBinning)
		.def("finishInit", &rpwa::pwaLikelihood<std::complex<double> >::finishInit)
		.def("Gradient", ::pwaLikelihood_Gradient)
		.def("FdF", ::pwaLikelihood_FdF)
		.def("DoEval", ::pwaLikelihood_DoEval)
		.def("DoDerivative", ::pwaLikelihood_DoDerivative)
		.def("Hessian", ::pwaLikelihood_Hessian)
		.def("HessianEigenVectors", ::pwaLikelihood_HessianEigenVectors)
		.def("CovarianceMatrix", ::pwaLikelihood_CovarianceMatrixFromMatrix)
		.def("CovarianceMatrix", ::pwaLikelihood_CovarianceMatrixFromPar)
		.def("CorrectParamSigns", ::pwaLikelihood_CorrectParamSigns)
		.def("nmbEvents", &rpwa::pwaLikelihood<std::complex<double> >::nmbEvents)
		.def("rank", &rpwa::pwaLikelihood<std::complex<double> >::rank)
		.def(
			"nmbWaves"
			, &rpwa::pwaLikelihood<std::complex<double> >::nmbWaves
			, (bp::arg("reflectivity") = 0)
		)
		.def("nmbPars", &rpwa::pwaLikelihood<std::complex<double> >::nmbPars)
		.def("nmbParsFixed", &rpwa::pwaLikelihood<std::complex<double> >::nmbParsFixed)
		.def(
			"parameter"
			, &rpwa::pwaLikelihood<std::complex<double> >::parameter
			, bp::return_internal_reference<>()
		)
		.def(
			"parameters"
			, &rpwa::pwaLikelihood<std::complex<double> >::parameters
			, bp::return_internal_reference<>()
		)
		.def("useNormalizedAmps", &rpwa::pwaLikelihood<std::complex<double> >::useNormalizedAmps)
		.def("setPriorType", &rpwa::pwaLikelihood<std::complex<double> >::setPriorType)
		.def("priorType", &rpwa::pwaLikelihood<std::complex<double> >::priorType)
		.def("setCauchyWidth", &rpwa::pwaLikelihood<std::complex<double> >::setCauchyWidth)
		.def("cauchyWidth", &rpwa::pwaLikelihood<std::complex<double> >::cauchyWidth)
		.def(
			"setQuiet"
			, &rpwa::pwaLikelihood<std::complex<double> >::setQuiet
			, (bp::arg("flag") = true)
		)
		.staticmethod("setQuiet")
		.def(bp::self_ns::str(bp::self))
		;


	bp::class_<rpwa::pwaLikelihood<std::complex<double> >::fitParameter>("fitParameter", bp::no_init)
		.def(
			"waveName"
			, &rpwa::pwaLikelihood<std::complex<double> >::fitParameter::waveName
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("rank", &rpwa::pwaLikelihood<std::complex<double> >::fitParameter::rank)
		.def("threshold", &rpwa::pwaLikelihood<std::complex<double> >::fitParameter::threshold)
		.def("fixed", &rpwa::pwaLikelihood<std::complex<double> >::fitParameter::fixed)
		.def("realPart", &rpwa::pwaLikelihood<std::complex<double> >::fitParameter::realPart)
		.def("parName", &rpwa::pwaLikelihood<std::complex<double> >::fitParameter::parName)
		;


	bp::enum_<rpwa::pwaLikelihood<std::complex<double> >::priorEnum>("priorEnum")
		.value("FLAT", rpwa::pwaLikelihood<std::complex<double> >::FLAT)
		.value("HALF_CAUCHY", rpwa::pwaLikelihood<std::complex<double> >::HALF_CAUCHY)
		.export_values();

}
