
#include "fitResult_py.h"

#include <TTree.h>

#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;

namespace {

	void fitResult_fill(rpwa::fitResult& self,
	                    const unsigned int                        nmbEvents,
	                    const unsigned int                        normNmbEvents,
	                    const double                              massBinCenter,
	                    const double                              logLikelihood,
	                    const int                                 rank,
	                    const bp::object&                         pyProdAmps,
	                    const bp::object&                         pyProdAmpNames,
	                    PyObject*                                 pyFitParCovMatrix,
	                    const bp::object&                         pyFitParCovMatrixIndices,
	                    const rpwa::complexMatrix&                normIntegral,
	                    const rpwa::complexMatrix&                acceptedNormIntegral,
	                    const bp::object&                         pyPhaseSpaceIntegral,
	                    const bool                                converged,
	                    const bool                                hasHessian)
	{
		std::vector<std::complex<double> > prodAmps;
		if(not rpwa::py::convertBPObjectToVector<std::complex<double> >(pyProdAmps, prodAmps)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for prodAmps when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		std::vector<std::string> prodAmpNames;
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyProdAmpNames, prodAmpNames)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for prodAmpNames when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		TMatrixT<double>* fitParCovMatrix = rpwa::py::convertFromPy<TMatrixT<double>* >(pyFitParCovMatrix);
		if(not fitParCovMatrix) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for fitParCovMatrix when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		bp::list pyListFitParCovMatrixIndices = bp::extract<bp::list>(pyFitParCovMatrixIndices);
		std::vector<std::pair<int, int> > fitParCovMatrixIndices(bp::len(pyListFitParCovMatrixIndices));
		for(int i = 0; i < bp::len(pyListFitParCovMatrixIndices); ++i) {
			if(not rpwa::py::convertBPObjectToPair<int, int>(pyListFitParCovMatrixIndices[i], fitParCovMatrixIndices[i]))
			{
				std::stringstream strStr;
				strStr<<"Could not convert element "<<i<<" when executing rpwa::fitResult::fill()";
				PyErr_SetString(PyExc_TypeError, strStr.str().c_str());
				bp::throw_error_already_set();
			}
		}
		std::vector<double> phaseSpaceIntegral;
		if(not rpwa::py::convertBPObjectToVector<double>(pyPhaseSpaceIntegral, phaseSpaceIntegral)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for phaseSpaceIntegral when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		self.fill(nmbEvents, normNmbEvents, massBinCenter, logLikelihood, rank, prodAmps, prodAmpNames, *fitParCovMatrix,
		          fitParCovMatrixIndices, normIntegral, acceptedNormIntegral, phaseSpaceIntegral, converged, hasHessian);
	}

	bp::list fitResult_evidenceComponents(const rpwa::fitResult& self) {
		return bp::list(self.evidenceComponents());
	}

	std::string fitResult_waveName(const rpwa::fitResult& self, const unsigned int waveIndex)
	{
		return std::string((const char*)self.waveName(waveIndex));
	}

	std::string fitResult_waveNameEsc(const rpwa::fitResult& self, const unsigned int waveIndex)
	{
		return std::string((const char*)self.waveNameEsc(waveIndex));
	}

	std::string fitResult_prodAmpName(const rpwa::fitResult& self, const unsigned int prodAmpIndex)
	{
		return std::string((const char*)self.prodAmpName(prodAmpIndex));
	}

	std::string fitResult_prodAmpNameEsc(const rpwa::fitResult& self, const unsigned int prodAmpIndex)
	{
		return std::string((const char*)self.prodAmpNameEsc(prodAmpIndex));
	}

	std::string fitResult_waveNameForProdAmp(const rpwa::fitResult& self, const unsigned int prodAmpIndex)
	{
		return std::string((const char*)self.waveNameForProdAmp(prodAmpIndex));
	}

	PyObject* fitResult_prodAmpCov_1(const rpwa::fitResult& self, const unsigned int prodAmpIndex)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.prodAmpCov(prodAmpIndex));
	}

	PyObject* fitResult_prodAmpCov_2(const rpwa::fitResult& self, const std::vector<unsigned int>& prodAmpIndices)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.prodAmpCov(prodAmpIndices));
	}

	PyObject* fitResult_prodAmpCov_3(const rpwa::fitResult& self, const bp::object& pyProdAmpIndexPairs)
	{
		bp::list pyListProdAmpIndexPairs = bp::extract<bp::list>(pyProdAmpIndexPairs);
		std::vector<std::pair<unsigned int, unsigned int> > prodAmpIndexPairs(bp::len(pyListProdAmpIndexPairs));
		for(int i = 0; i < bp::len(pyListProdAmpIndexPairs); ++i) {
			bp::tuple tuple = bp::extract<bp::tuple>(pyListProdAmpIndexPairs[i]);
			prodAmpIndexPairs[i].first  = bp::extract<unsigned int>(tuple[0]);
			prodAmpIndexPairs[i].second = bp::extract<unsigned int>(tuple[1]);
		}
		return rpwa::py::convertToPy<TMatrixT<double> >(self.prodAmpCov(prodAmpIndexPairs));
	}

	PyObject* fitResult_prodAmpCov_4(const rpwa::fitResult& self,
	                                 const std::vector<unsigned int>& prodAmpIndicesA,
	                                 const std::vector<unsigned int>& prodAmpIndicesB)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.prodAmpCov(prodAmpIndicesA, prodAmpIndicesB));
	}

	double fitResult_phaseSpaceIntegral_1(const rpwa::fitResult& self, const unsigned int waveIndex)
	{
		return self.phaseSpaceIntegral(waveIndex);
	}

	double fitResult_phaseSpaceIntegral_2(const rpwa::fitResult& self, const std::string& waveName)
	{
		return self.phaseSpaceIntegral(waveName);
	}

	PyObject* fitResult_spinDensityMatrixElemCov(const rpwa::fitResult& self,
	                                             const unsigned int waveIndexA,
	                                             const unsigned int waveIndexB)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.spinDensityMatrixElemCov(waveIndexA, waveIndexB));
	}

	double fitResult_intensity_1(const rpwa::fitResult& self, const unsigned int waveIndex) {
		return self.intensity(waveIndex);
	}

	double fitResult_intensityErr_1(const rpwa::fitResult& self, const unsigned int waveIndex) {
			return self.intensityErr(waveIndex);
	}

	double fitResult_intensity_2(const rpwa::fitResult& self, const char* waveNamePattern) {
		return self.intensity(waveNamePattern);
	}

	double fitResult_intensityErr_2(const rpwa::fitResult& self, const char* waveNamePattern) {
			return self.intensityErr(waveNamePattern);
	}

	double fitResult_intensity_3(const rpwa::fitResult& self) {
		return self.intensity();
	}

	double fitResult_intensityErr_3(const rpwa::fitResult& self) {
			return self.intensityErr();
	}

	double fitResult_phase_1(const rpwa::fitResult& self,
	                         const unsigned int waveIndexA,
	                         const unsigned int waveIndexB)
	{
		return self.phase(waveIndexA, waveIndexB);
	}

	double fitResult_phaseErr_1(const rpwa::fitResult& self,
	                            const unsigned int waveIndexA,
	                            const unsigned int waveIndexB)
	{
		return self.phaseErr(waveIndexA, waveIndexB);
	}

	double fitResult_phase_2(const rpwa::fitResult& self,
	                         const std::string waveNameA,
	                         const std::string waveNameB)
	{
		return self.phase(waveNameA, waveNameB);
	}

	double fitResult_phaseErr_2(const rpwa::fitResult& self,
	                            const std::string waveNameA,
	                            const std::string waveNameB)
	{
		return self.phaseErr(waveNameA, waveNameB);
	}

	bp::list fitResult_prodAmps(const rpwa::fitResult& self)
	{
		bp::list retval;
		const std::vector<TComplex>& prodAmps = self.prodAmps();
		for(unsigned int i = 0; i < prodAmps.size(); ++i) {
			retval.append(std::complex<double>(prodAmps[i].Re(), prodAmps[i].Im()));
		}
		return retval;
	}

	bp::list fitResult_prodAmpNames(const rpwa::fitResult& self)
	{
		return bp::list(self.prodAmpNames());
	}

	bp::list fitResult_waveNames(const rpwa::fitResult& self)
	{
		return bp::list(self.waveNames());
	}

	PyObject* fitResult_fitParCovMatrix(const rpwa::fitResult self)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.fitParCovMatrix());
	}

	bp::list fitResult_fitParCovIndices(const rpwa::fitResult self)
	{
		bp::list retval;
		const std::vector<std::pair<Int_t, Int_t> >& fitParCovIndices = self.fitParCovIndices();
		for(unsigned int i = 0; i < fitParCovIndices.size(); ++i) {
			const std::pair<Int_t, Int_t>& item = fitParCovIndices[i];
			retval.append(bp::make_tuple(item.first, item.second));
		}
		return retval;
	}

	bp::dict fitResult_normIntIndexMap(const rpwa::fitResult self)
	{
		bp::dict retval;
		const std::map<Int_t, Int_t>& normIntIndexMap = self.normIntIndexMap();
		for(std::map<Int_t, Int_t>::const_iterator it = normIntIndexMap.begin(); it != normIntIndexMap.end(); ++it)
		{
			retval[it->first] = it->second;
		}
		return retval;
	}

	std::string fitResult_printProdAmpNames(const rpwa::fitResult self)
	{
		std::stringstream sstr;
		self.printProdAmpNames(sstr);
		return sstr.str();
	}

	std::string fitResult_printWaveNames(const rpwa::fitResult self)
	{
		std::stringstream sstr;
		self.printWaveNames(sstr);
		return sstr.str();
	}

	std::string fitResult_printProdAmps(const rpwa::fitResult self)
	{
		std::stringstream sstr;
		self.printProdAmps(sstr);
		return sstr.str();
	}

	std::string fitResult_printWaves(const rpwa::fitResult self)
	{
		std::stringstream sstr;
		self.printWaves(sstr);
		return sstr.str();
	}

	std::string fitResult_printAmpsGenPW(const rpwa::fitResult self)
	{
		std::stringstream sstr;
		self.printAmpsGenPW(sstr);
		return sstr.str();
	}

	int fitResult_Write(const rpwa::fitResult& self, const char* name = 0)
	{
		return self.Write(name);
	}

}

void rpwa::py::exportFitResult() {

	bp::class_<rpwa::fitResult>("fitResult")
		.def(bp::init<const rpwa::fitResult&>())
		.def(bp::self_ns::str(bp::self))
		.def("reset", &rpwa::fitResult::reset)
		.def("fill", &fitResult_fill)
		.def("massBinCenter", &rpwa::fitResult::massBinCenter)
		.def("logLikelihood", &rpwa::fitResult::logLikelihood)
		.def("evidence", &rpwa::fitResult::evidence)
		.def("evidenceComponents", &fitResult_evidenceComponents)
		.def("converged", &rpwa::fitResult::converged)
		.def("hasHessian", &rpwa::fitResult::hasHessian)
		.def("rank", &rpwa::fitResult::rank)
		.def("nmbEvents", &rpwa::fitResult::nmbEvents)
		.def("normNmbEvents", &rpwa::fitResult::normNmbEvents)
		.def("nmbWaves", &rpwa::fitResult::nmbWaves)
		.def("nmbProdAmps", &rpwa::fitResult::nmbProdAmps)
		.def("waveName", &fitResult_waveName)
		.def("waveNameEsc", &fitResult_waveNameEsc)
		.def("prodAmpName", &fitResult_prodAmpName)
		.def("prodAmpNameEsc", &fitResult_prodAmpNameEsc)
		.def("waveNameForProdAmp", &fitResult_waveNameForProdAmp)
		.def("waveIndex", &rpwa::fitResult::waveIndex)
		.def("prodAmpIndex", &rpwa::fitResult::prodAmpIndex)
		.def("fitParameter", &rpwa::fitResult::fitParameter)
		.def("fitParameterCov", &rpwa::fitResult::fitParameterCov)
		.def("prodAmp", &fitResult::prodAmp)
		.def("prodAmpCov", &fitResult_prodAmpCov_3)
		.def("prodAmpCov", &fitResult_prodAmpCov_4)
		.def("prodAmpCov", &fitResult_prodAmpCov_2)
		.def("prodAmpCov", &fitResult_prodAmpCov_1)
		.def("covMatrixValid", &rpwa::fitResult::covMatrixValid)
		.def("normIntegral", &rpwa::fitResult::normIntegral)
		.def("acceptedNormIntegral", &rpwa::fitResult::acceptedNormIntegral)
		.def("phaseSpaceIntegral", &fitResult_phaseSpaceIntegral_1)
		.def("phaseSpaceIntegral", &fitResult_phaseSpaceIntegral_2)
		.def("spinDensityMatrixElem", &rpwa::fitResult::spinDensityMatrixElem)
		.def("spinDensityMatrixElemCov", &fitResult_spinDensityMatrixElemCov)
		.def("intensity", &fitResult_intensity_1)
		.def("intensityErr", &fitResult_intensityErr_1)
		.def("intensity", &fitResult_intensity_2)
		.def("intensityErr", &fitResult_intensityErr_2)
		.def("intensity", &fitResult_intensity_3)
		.def("intensityErr", &fitResult_intensityErr_3)
		.def("phase", &fitResult_phase_1)
		.def("phaseErr", &fitResult_phaseErr_1)
		.def("phase", &fitResult_phase_2)
		.def("phaseErr", &fitResult_phaseErr_2)
		.def("coherence", &rpwa::fitResult::coherence)
		.def("coherenceErr", &rpwa::fitResult::coherenceErr)
		.def("overlap", &rpwa::fitResult::overlap)
		.def("overlapErr", &rpwa::fitResult::overlapErr)
		.def("prodAmps", &fitResult_prodAmps)
		.def("prodAmpNames", &fitResult_prodAmpNames)
		.def("waveNames", &fitResult_waveNames)
		.def("fitParCovMatrix", &fitResult_fitParCovMatrix)
		.def("fitParCovIndices", &fitResult_fitParCovIndices)
		.def(
			"normIntegralMatrix"
			, &rpwa::fitResult::normIntegralMatrix
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"acceptedNormIntegralMatrix"
			, &rpwa::fitResult::acceptedNormIntegralMatrix
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("normIntIndexMap", &fitResult_normIntIndexMap)
		.def("printProdAmpNames", &fitResult_printProdAmpNames)
		.def("printWaveNames", &fitResult_printWaveNames)
		.def("printProdAmps", &fitResult_printProdAmps)
		.def("printWaves", &fitResult_printWaves)
		.def("printAmpsGenPW", &fitResult_printAmpsGenPW)

		.def("Write", &fitResult_Write, bp::arg("name")=0)
		.def("setBranchAddress", &rpwa::py::setBranchAddress<rpwa::fitResult*>)
		.def("branch", &rpwa::py::branch<rpwa::fitResult*>);
	bp::register_ptr_to_python<rpwa::fitResultPtr>();
}
