#include "pwaLikelihood_py.h"

#include <complex>

#include "rootConverters_py.h"
#include "stlContainers_py.h"


namespace bp = boost::python;

namespace {

	void
	pwaLikelihood_init(rpwa::pwaLikelihood<std::complex<double> >& self,
	                   const unsigned int                          rank,
	                   PyObject*                                   pyAmpFileList,
	                   const double                                massBinCenter,
	                   const std::string&                          waveListFileName,
	                   const std::string&                          normIntFileName,
	                   const std::string&                          accIntFileName,
	                   const unsigned int                          numbAccEvents)
	{
		std::map<std::string, std::string> ampFilesMap;
		bp::dict pyDictAmpFilesDict = bp::extract<bp::dict>(pyAmpFileList);
		bp::list keys = pyDictAmpFilesDict.keys();
		for(int i = 0; i < bp::len(keys); ++i) {
			std::string waveName = bp::extract<std::string>(keys[i]);
			std::string fileName = bp::extract<std::string>(pyDictAmpFilesDict[keys[i]]);
			ampFilesMap.insert(std::pair<std::string, std::string>(waveName, fileName));
		}
		self.init(rank, ampFilesMap, massBinCenter, waveListFileName, normIntFileName, accIntFileName, numbAccEvents);
	}

	bp::list
	pwaLikelihood_Gradient(rpwa::pwaLikelihood<std::complex<double> >& self, bp::list& pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		printErr << nmbPar << std::endl;
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		std::vector<double> gradient(nmbPar, 0.);
		self.Gradient(par.data(), gradient.data());
		return bp::list(gradient);
	}

	bp::tuple
	pwaLikelihood_FdF(rpwa::pwaLikelihood<std::complex<double> >& self, bp::list& pyPar)
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
	pwaLikelihood_DoEval(rpwa::pwaLikelihood<std::complex<double> >& self, bp::list& pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return self.DoEval(par.data());
	}

	double
	pwaLikelihood_DoDerivative(rpwa::pwaLikelihood<std::complex<double> >& self, bp::list& pyPar, unsigned int derivIndex)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return self.DoDerivative(par.data(), derivIndex);
	}

	PyObject*
	pwaLikelihood_Hessian(rpwa::pwaLikelihood<std::complex<double> >& self, bp::list& pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return rpwa::py::convertToPy<TMatrixT<double> >(self.Hessian(par.data()));
	}

	PyObject*
	pwaLikelihood_CovarianceMatrixFromPar(rpwa::pwaLikelihood<std::complex<double> >& self, bp::list& pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return rpwa::py::convertToPy<TMatrixT<double> >(self.CovarianceMatrix(par.data()));
	}

	PyObject*
	pwaLikelihood_CovarianceMatrixFromMatrix(rpwa::pwaLikelihood<std::complex<double> >& self, PyObject* pyCovMatrix)
	{
		TMatrixT<double>* covMatrix = rpwa::py::convertFromPy<TMatrixT<double>* >(pyCovMatrix);
		if(not covMatrix) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for covMatrix when executing rpwa::pwaLikelihood::CovarianceMatrix()");
			bp::throw_error_already_set();
		}
		return rpwa::py::convertToPy<TMatrixT<double> >(self.CovarianceMatrix(*covMatrix));
	}

	bp::list
	pwaLikelihood_CorrectParamSigns(rpwa::pwaLikelihood<std::complex<double> >& self, bp::list& pyPar)
	{
		const unsigned int nmbPar = bp::len(pyPar);
		std::vector<double> par(nmbPar, 0);
		for(unsigned int i = 0; i < nmbPar; ++i) {
			par[i] = bp::extract<double>(pyPar[i]);
		}
		return bp::list(self.CorrectParamSigns(par.data()));
	}
}

void rpwa::py::exportPwaLikelihood() {

	bp::class_<rpwa::pwaLikelihood<std::complex<double> > >("pwaLikelihood")
		.def("init",
		     ::pwaLikelihood_init,
		     (bp::arg("rank"),
		      bp::arg("ampFileList"),
		      bp::arg("massBinCenter"),
		      bp::arg("waveListFileName"),
		      bp::arg("normIntFileName"),
		      bp::arg("accIntFileName"),
		      bp::arg("numbAccEvents") = 0)
		     )
		.def("Gradient", ::pwaLikelihood_Gradient)
		.def("FdF", ::pwaLikelihood_FdF)
		.def("DoEval", ::pwaLikelihood_DoEval)
		.def("DoDerivative", ::pwaLikelihood_DoDerivative)
		.def("Hessian", ::pwaLikelihood_Hessian)
		.def("CovarianceMatrix", ::pwaLikelihood_CovarianceMatrixFromPar)
		.def("CovarianceMatrix", ::pwaLikelihood_CovarianceMatrixFromMatrix)
		.def("CorrectParamSigns", ::pwaLikelihood_CorrectParamSigns)
		.def("nmbEvents", &rpwa::pwaLikelihood<std::complex<double> >::nmbEvents)
		.def("rank", &rpwa::pwaLikelihood<std::complex<double> >::rank)
		.def("nmbWaves", &rpwa::pwaLikelihood<std::complex<double> >::nmbWaves,
		     (bp::arg("reflectivity")=0))
		.def("nmbPars", &rpwa::pwaLikelihood<std::complex<double> >::nmbPars)
		.def("nmbParsFixed", &rpwa::pwaLikelihood<std::complex<double> >::nmbParsFixed)
		.def("parName", &rpwa::pwaLikelihood<std::complex<double> >::parName)
		.def("parThreshold", &rpwa::pwaLikelihood<std::complex<double> >::parThreshold)
		.def("parFixed", &rpwa::pwaLikelihood<std::complex<double> >::parFixed)
		.def("useNormalizedAmps", &rpwa::pwaLikelihood<std::complex<double> >::useNormalizedAmps)
		.def("setPriorType", &rpwa::pwaLikelihood<std::complex<double> >::setPriorType)
		.def("priorType", &rpwa::pwaLikelihood<std::complex<double> >::priorType)
		.def(
			"setQuiet"
			, &rpwa::pwaLikelihood<std::complex<double> >::setQuiet
			, (bp::arg("flag")=true)
		)
		.staticmethod("setQuiet")
		.def(bp::self_ns::str(bp::self));
	bp::enum_<rpwa::pwaLikelihood<std::complex<double> >::priorEnum>("priorEnum")
		.value("FLAT", rpwa::pwaLikelihood<std::complex<double> >::FLAT)
		.value("HALF_CAUCHY", rpwa::pwaLikelihood<std::complex<double> >::HALF_CAUCHY)
		.export_values();
}
