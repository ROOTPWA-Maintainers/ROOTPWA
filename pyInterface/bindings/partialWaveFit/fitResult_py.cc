#include "fitResult_py.h"

#define BOOST_PYTHON_MAX_ARITY 20
#include <boost/python.hpp>
#include <vector>
#include <list>
#include <map>

#include <TTree.h>

#include "fitResult.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	void fitResult_fill_1(rpwa::fitResult& self,
	                      const unsigned int                        nmbEvents,
	                      const unsigned int                        normNmbEvents,
	                      const bp::dict&                           pyMultibinBoundaries,
	                      const double                              logLikelihood,
	                      const int                                 rank,
	                      const bp::object&                         pyProdAmps,
	                      const bp::object&                         pyProdAmpNames,
	                      PyObject*                                 pyFitParCovMatrix,
	                      const bp::object&                         pyFitParCovMatrixIndices,
	                      const bp::object&                         pyNormIntegral,
	                      const bp::object&                         pyAcceptedNormIntegral,
	                      const bp::object&                         pyPhaseSpaceIntegral,
	                      const bool                                converged,
	                      const bool                                hasHessian,
	                      const bp::object&                         pyDatasetRatios,
	                      const bp::object&                         pyDatasetRatiosCovMatrixIndices,
	                      const bp::object&                         pyAdditionalFitParameters,
	                      const bp::object&                         pyAdditionalFitParametersCovMatrixIndices )
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
		TMatrixT<double>* fitParCovMatrix = nullptr;
		if (pyFitParCovMatrix != Py_None){
			fitParCovMatrix = rpwa::py::convertFromPy<TMatrixT<double>* >(pyFitParCovMatrix);
			if(not fitParCovMatrix) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for fitParCovMatrix when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
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
		std::vector<double>* phaseSpaceIntegralPtr = nullptr;
		if (not pyPhaseSpaceIntegral.is_none()){
			if(not rpwa::py::convertBPObjectToVector<double>(pyPhaseSpaceIntegral, phaseSpaceIntegral)) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for phaseSpaceIntegral when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
			phaseSpaceIntegralPtr = &phaseSpaceIntegral;
		}
		rpwa::complexMatrix* normIntegral = nullptr;
		if (not pyNormIntegral.is_none()){
			normIntegral = bp::extract<rpwa::complexMatrix*>(pyNormIntegral);
			if (normIntegral == nullptr){
				PyErr_SetString(PyExc_TypeError, "Got invalid input for normIntegral when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
		}
		rpwa::complexMatrix* acceptedNormIntegral = nullptr;
		if (not pyAcceptedNormIntegral.is_none()){
			acceptedNormIntegral = bp::extract<rpwa::complexMatrix*>(pyAcceptedNormIntegral);
			if (acceptedNormIntegral == nullptr){
				PyErr_SetString(PyExc_TypeError, "Got invalid input for acceptedNormIntegral when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
		}
		const rpwa::multibinBoundariesType multibinBoundaries = rpwa::py::convertMultibinBoundariesFromPy(pyMultibinBoundaries);

		std::unique_ptr<std::map<std::string, double>> datasetRatios = nullptr;
		if (not pyDatasetRatios.is_none()){
			datasetRatios.reset(new std::map<std::string, double>);
			rpwa::py::convertBPObjectToMap(pyDatasetRatios, *datasetRatios);
		}
		std::unique_ptr<std::map<std::string, int>> datasetRatiosCovMatrixIndices = nullptr;
		if (not pyDatasetRatiosCovMatrixIndices.is_none()){
			datasetRatiosCovMatrixIndices.reset(new std::map<std::string, int>);
			rpwa::py::convertBPObjectToMap(pyDatasetRatiosCovMatrixIndices, *datasetRatiosCovMatrixIndices);
		}
		std::unique_ptr<std::map<std::string, double>> additionalFitParameters = nullptr;
		if (not pyAdditionalFitParameters.is_none()){
			additionalFitParameters.reset(new std::map<std::string, double>);
			rpwa::py::convertBPObjectToMap(pyAdditionalFitParameters, *additionalFitParameters);
		}
		std::unique_ptr<std::map<std::string, int>> additionalFitParametersCovMatrixIndices = nullptr;
		if (not pyAdditionalFitParametersCovMatrixIndices.is_none()){
			additionalFitParametersCovMatrixIndices.reset(new std::map<std::string, int>);
			rpwa::py::convertBPObjectToMap(pyAdditionalFitParametersCovMatrixIndices, *additionalFitParametersCovMatrixIndices);
		}
		self.fill(nmbEvents, normNmbEvents, multibinBoundaries, logLikelihood, rank, prodAmps, prodAmpNames, fitParCovMatrix,
		          fitParCovMatrixIndices, normIntegral, acceptedNormIntegral, phaseSpaceIntegralPtr, converged, hasHessian,
		          datasetRatios.get(), datasetRatiosCovMatrixIndices.get(), additionalFitParameters.get(), additionalFitParametersCovMatrixIndices.get());
	}

	void fitResult_fill_2(rpwa::fitResult& self, const rpwa::fitResult& result, const bool fillCovMatrix, const bool fillIntegralMatrices) {
		self.fill(result, fillCovMatrix, fillIntegralMatrices);
	}

	void fitResult_fill_3(rpwa::fitResult& self, const rpwa::fitResult& result,
	                      PyObject* pyFitParCovMatrix,
	                      const bp::object& pyFitParCovMatrixIndices,
	                      const bp::object& pyNormIntegral,
	                      const bp::object& pyAcceptedNormIntegral,
	                      const bp::object& pyPhaseSpaceIntegral)
	                      {
		TMatrixT<double>* fitParCovMatrix = nullptr;
		if (pyFitParCovMatrix != nullptr) {
			fitParCovMatrix = rpwa::py::convertFromPy<TMatrixT<double>*>(pyFitParCovMatrix);
			if (not fitParCovMatrix) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for fitParCovMatrix when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
		}

		std::vector < std::pair<int, int> > fitParCovMatrixIndices;
		std::vector < std::pair<int, int> > *fitParCovMatrixIndicesPtr = nullptr;
		if (not pyFitParCovMatrixIndices.is_none()) { // not NONE
			fitParCovMatrixIndicesPtr = &fitParCovMatrixIndices;
			bp::extract<std::vector<std::pair<int, int>>> value(pyFitParCovMatrixIndices);
			if (not value.check()) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for fitParCovMatrixIndices when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
			fitParCovMatrixIndices = value;
		}

		rpwa::complexMatrix normIntegral;
		rpwa::complexMatrix* normIntegralPtr = nullptr;
		if (not pyNormIntegral.is_none()) { // not NONE
			normIntegralPtr = &normIntegral;
			bp::extract<rpwa::complexMatrix> value(pyNormIntegral);
			if (not value.check()) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for normIntegral when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
			normIntegral = value;
		}

		rpwa::complexMatrix acceptedNormIntegral;
		rpwa::complexMatrix* acceptedNormIntegralPtr = nullptr;
		if (not pyAcceptedNormIntegral.is_none()) { // not NONE
			acceptedNormIntegralPtr = &acceptedNormIntegral;
			bp::extract<rpwa::complexMatrix> value(pyAcceptedNormIntegral);
			if (not value.check()) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for acceptedNormIntegral when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
			acceptedNormIntegral = value;
		}

		std::vector<double> phaseSpaceIntegral;
		std::vector<double>* phaseSpaceIntegralPtr = nullptr;
		if (not pyPhaseSpaceIntegral.is_none()) { // not NONE
			phaseSpaceIntegralPtr = &phaseSpaceIntegral;
			if (not rpwa::py::convertBPObjectToVector<double>(pyPhaseSpaceIntegral, phaseSpaceIntegral)) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for phaseSpaceIntegral when executing rpwa::fitResult::fill()");
				bp::throw_error_already_set();
			}
		}

		self.fill(result, fitParCovMatrix, fitParCovMatrixIndicesPtr, normIntegralPtr, acceptedNormIntegralPtr, phaseSpaceIntegralPtr);
	}

	bp::dict fitResult_multibinCenter(const rpwa::fitResult& self) {
		const rpwa::multibinCenterType multibinCenter = self.multibinCenter();

		bp::dict pyMultibinCenter;
		for(rpwa::multibinCenterType::const_iterator it = multibinCenter.begin(); it != multibinCenter.end(); ++it){
			pyMultibinCenter[it->first] = it->second;
		}

		return pyMultibinCenter;
	}

	bp::dict fitResult_multibinBoundaries(const rpwa::fitResult& self) {
		return rpwa::py::convertMultibinBoundariesToPy(self.multibinBoundaries());
	}

	bp::list fitResult_evidenceComponents(const rpwa::fitResult& self) {
		return bp::list(self.evidenceComponents());
	}

	bp::list fitResult_waveIndicesMatchingPattern(const rpwa::fitResult& self,
	                                              const std::string&     pattern) {
		return bp::list(self.waveIndicesMatchingPattern(pattern));
	}

	PyObject* fitResult_prodAmpCov_1(const rpwa::fitResult& self, const unsigned int prodAmpIndex)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.prodAmpCov(prodAmpIndex));
	}

	PyObject* fitResult_prodAmpCov_2(const rpwa::fitResult& self, const std::vector<unsigned int>& prodAmpIndices)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.prodAmpCov(prodAmpIndices));
	}

	double fitResult_phaseSpaceIntegral_1(const rpwa::fitResult& self, const unsigned int waveIndex)
	{
		return self.phaseSpaceIntegral(waveIndex);
	}

	double fitResult_phaseSpaceIntegral_2(const rpwa::fitResult& self, const std::string& waveName)
	{
		return self.phaseSpaceIntegral(waveName);
	}

	std::complex<double> fitResult_spinDensityMatrixElem_1(const rpwa::fitResult& self,
	                                                       const unsigned int waveIndexA,
	                                                       const unsigned int waveIndexB)
	{
		return self.spinDensityMatrixElem(waveIndexA, waveIndexB);
	}

	PyObject* fitResult_spinDensityMatrixElemCov_1(const rpwa::fitResult& self,
	                                               const unsigned int waveIndexA,
	                                               const unsigned int waveIndexB)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.spinDensityMatrixElemCov(waveIndexA, waveIndexB));
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

	double fitResult_coherence_1(const rpwa::fitResult& self,
	                             const unsigned int waveIndexA,
	                             const unsigned int waveIndexB)
	{
		return self.coherence(waveIndexA, waveIndexB);
	}

	double fitResult_coherenceErr_1(const rpwa::fitResult& self,
	                                const unsigned int waveIndexA,
	                                const unsigned int waveIndexB)
	{
		return self.coherenceErr(waveIndexA, waveIndexB);
	}

	double fitResult_overlap_1(const rpwa::fitResult& self,
	                           const unsigned int waveIndexA,
	                           const unsigned int waveIndexB)
	{
		return self.overlap(waveIndexA, waveIndexB);
	}

	double fitResult_overlapErr_1(const rpwa::fitResult& self,
	                              const unsigned int waveIndexA,
	                              const unsigned int waveIndexB)
	{
		return self.overlapErr(waveIndexA, waveIndexB);
	}

	std::complex<double> fitResult_spinDensityMatrixElem_2(const rpwa::fitResult& self,
	                                                       const std::string& waveNameA,
	                                                       const std::string& waveNameB)
	{
		return self.spinDensityMatrixElem(waveNameA, waveNameB);
	}

	PyObject* fitResult_spinDensityMatrixElemCov_2(const rpwa::fitResult& self,
	                                               const std::string& waveNameA,
	                                               const std::string& waveNameB)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.spinDensityMatrixElemCov(waveNameA, waveNameB));
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

	double fitResult_coherence_2(const rpwa::fitResult& self,
	                             const std::string waveNameA,
	                             const std::string waveNameB)
	{
		return self.coherence(waveNameA, waveNameB);
	}

	double fitResult_coherenceErr_2(const rpwa::fitResult& self,
	                                const std::string waveNameA,
	                                const std::string waveNameB)
	{
		return self.coherenceErr(waveNameA, waveNameB);
	}

	double fitResult_overlap_2(const rpwa::fitResult& self,
	                           const std::string waveNameA,
	                           const std::string waveNameB)
	{
		return self.overlap(waveNameA, waveNameB);
	}

	double fitResult_overlapErr_2(const rpwa::fitResult& self,
	                              const std::string waveNameA,
	                              const std::string waveNameB)
	{
		return self.overlapErr(waveNameA, waveNameB);
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

	bp::list fitResult_phaseSpaceIntegralVector(const rpwa::fitResult& self)
	{
		return bp::list(self.phaseSpaceIntegralVector());
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

	int fitResult_Write(const rpwa::fitResult& self, const char* name = 0)
	{
		return self.Write(name);
	}


	bp::dict
	fitResult_getFitResultsFromFilesInMutibins(bp::list& fileNamesPy, const std::string& treeName,
	                                           const std::string& branchName,
	                                           const bool onlyBestInMultibin, const bool stripMatricesFromFurtherAttempts,
	                                           const bool onlyConvergedResults) {
		bp::dict pyFitResultsInMultibins;
		std::vector<std::string> fileNames;
		if( rpwa::py::convertBPObjectToVector(fileNamesPy, fileNames)){

			std::map<rpwa::multibinBoundariesType, std::list<rpwa::fitResult> > fitResultsInMultibins =
					rpwa::getFitResultsFromFilesInMultibins(fileNames,
															treeName,
															branchName,
															onlyBestInMultibin,
															stripMatricesFromFurtherAttempts,
															onlyConvergedResults);


			for (auto& elem : fitResultsInMultibins) {
				const rpwa::multibinBoundariesType& multibinBoundaries = elem.first;
				std::list<rpwa::fitResult>& results = elem.second;
				bp::list pyResults;
				while(not results.empty()){
					pyResults.append(results.front());
					results.pop_front();
				}

				bp::object pyMultiBin = rpwa::py::convertMultibinBoundariesToPyMultibin(multibinBoundaries);
				pyFitResultsInMultibins[pyMultiBin] = pyResults;
			}
		}
		return pyFitResultsInMultibins;
	}

}

void rpwa::py::exportFitResult() {

	bp::def("escapeRegExpSpecialChar", &rpwa::escapeRegExpSpecialChar);
	bp::def("unescapeRegExpSpecialChar", &rpwa::unescapeRegExpSpecialChar);

	bp::class_<rpwa::fitResult>("fitResult")
		.def(bp::init<const rpwa::fitResult&>())
		.def(bp::self_ns::str(bp::self))
		.def("reset", &rpwa::fitResult::reset)
		.def("fill", &fitResult_fill_1, (bp::arg("nmbEvents"),
		                                 bp::arg("normNmbEvents"),
		                                 bp::arg("pyMultibinBoundaries"),
		                                 bp::arg("logLikelihood"),
		                                 bp::arg("rank"),
		                                 bp::arg("pyProdAmps"),
		                                 bp::arg("pyProdAmpNames"),
		                                 bp::arg("pyFitParCovMatrix"),
		                                 bp::arg("pyFitParCovMatrixIndices"),
		                                 bp::arg("pyNormIntegral"),
		                                 bp::arg("pyAcceptedNormIntegral"),
		                                 bp::arg("pyPhaseSpaceIntegral"),
		                                 bp::arg("converged"),
		                                 bp::arg("hasHessian"),
		                                 bp::arg("pyDatasetRatios") = bp::object(),
		                                 bp::arg("pyDatasetRatiosCovMatrixIndices") = bp::object(),
		                                 bp::arg("pyAdditionalFitParameters") = bp::object(),
		                                 bp::arg("pyAdditionalFitParametersCovMatrixIndices") = bp::object() ))
		.def("fill", &fitResult_fill_2, (bp::arg("fillCovMatrix") = true, bp::arg("fillIntegralMatrices") = true))
		.def("fill", &fitResult_fill_3)
		.def("multibinBoundaries", &::fitResult_multibinBoundaries)
		.def("multibinCenter", &::fitResult_multibinCenter)
		.def("nmbEvents", &rpwa::fitResult::nmbEvents)
		.def("normNmbEvents", &rpwa::fitResult::normNmbEvents)
		.def("massBinCenter", &rpwa::fitResult::massBinCenter)
		.def("logLikelihood", &rpwa::fitResult::logLikelihood)
		.def("evidence", &rpwa::fitResult::evidence)
		.def("evidenceComponents", &fitResult_evidenceComponents)
		.def("rank", &rpwa::fitResult::rank)
		.def("covMatrixValid", &rpwa::fitResult::covMatrixValid)
		.def("converged", &rpwa::fitResult::converged)
		.def("hasHessian", &rpwa::fitResult::hasHessian)
		.def("nmbWaves", &rpwa::fitResult::nmbWaves)
		.def("nmbProdAmps", &rpwa::fitResult::nmbProdAmps)
		.def(
			"waveName"
			, &rpwa::fitResult::waveName
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("waveNameEsc", &rpwa::fitResult::waveNameEsc)
		.def(
			"prodAmpName"
			, &rpwa::fitResult::prodAmpName
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("prodAmpNameEsc", &rpwa::fitResult::prodAmpNameEsc)
		.def("waveNameForProdAmp", &rpwa::fitResult::waveNameForProdAmp)
		.def("rankOfProdAmp", &rpwa::fitResult::rankOfProdAmp)
		.def("waveIndex", &rpwa::fitResult::waveIndex)
		.def("prodAmpIndex", &rpwa::fitResult::prodAmpIndex)
		.def("waveIndicesMatchingPattern", &fitResult_waveIndicesMatchingPattern)
		.def("fitParameter", &rpwa::fitResult::fitParameter)
		.def("datasetRatio", &rpwa::fitResult::datasetRatio)
		.def("datasetLabels", &rpwa::fitResult::datasetLabels)
		.def("prodAmp", &fitResult::prodAmp)
		.def("prodAmpCov", &fitResult_prodAmpCov_2)
		.def("prodAmpCov", &fitResult_prodAmpCov_1)
		.def("normIntegral", &rpwa::fitResult::normIntegral)
		.def("acceptedNormIntegral", &rpwa::fitResult::acceptedNormIntegral)
		.def("phaseSpaceIntegral", &fitResult_phaseSpaceIntegral_1)
		.def("phaseSpaceIntegral", &fitResult_phaseSpaceIntegral_2)

		.def("spinDensityMatrixElem", &fitResult_spinDensityMatrixElem_1)
		.def("spinDensityMatrixElemCov", &fitResult_spinDensityMatrixElemCov_1)
		.def("phase", &fitResult_phase_1)
		.def("phaseErr", &fitResult_phaseErr_1)
		.def("coherence", &fitResult_coherence_1)
		.def("coherenceErr", &fitResult_coherenceErr_1)
		.def("overlap", &fitResult_overlap_1)
		.def("overlapErr", &fitResult_overlapErr_1)

		.def("spinDensityMatrixElem", &fitResult_spinDensityMatrixElem_2)
		.def("spinDensityMatrixElemCov", &fitResult_spinDensityMatrixElemCov_2)
		.def("phase", &fitResult_phase_2)
		.def("phaseErr", &fitResult_phaseErr_2)
		.def("coherence", &fitResult_coherence_2)
		.def("coherenceErr", &fitResult_coherenceErr_2)
		.def("overlap", &fitResult_overlap_2)
		.def("overlapErr", &fitResult_overlapErr_2)

		.def("intensity", &fitResult_intensity_1)
		.def("intensityErr", &fitResult_intensityErr_1)
		.def("intensity", &fitResult_intensity_2)
		.def("intensityErr", &fitResult_intensityErr_2)
		.def("intensity", &fitResult_intensity_3)
		.def("intensityErr", &fitResult_intensityErr_3)

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
		.def("phaseSpaceIntegralVector", &fitResult_phaseSpaceIntegralVector)
		.def("normIntIndexMap", &fitResult_normIntIndexMap)
		.def("printProdAmps", &fitResult_printProdAmps)
		.def("printWaves", &fitResult_printWaves)

		.def("Write", &fitResult_Write, bp::arg("name")=0)
		.def("setBranchAddress", &rpwa::py::setBranchAddress<rpwa::fitResult*>)
		.def(
			"branch"
			, &rpwa::py::branch<rpwa::fitResult*>
			, (bp::arg("fitResult"),
			   bp::arg("tree"),
			   bp::arg("name"),
			   bp::arg("bufsize")=32000,
			   bp::arg("splitlevel")=99)
		)
		;

	bp::register_ptr_to_python<rpwa::fitResultPtr>();

	bp::def("getFitResultsFromFilesInMultibins", &fitResult_getFitResultsFromFilesInMutibins);

}
