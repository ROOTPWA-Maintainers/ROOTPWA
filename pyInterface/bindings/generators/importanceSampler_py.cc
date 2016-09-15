#include"importanceSampler_py.h"

#include<TFile.h>

#include"importanceSampler.h"
#include"rootConverters_py.h"
#include"stlContainers_py.h"

namespace bp = boost::python;


namespace {


	bool
	importanceSampler_initializeFileWriter(rpwa::importanceSampler& self,
	                                       PyObject*                outFilePy,
	                                       const std::string&       userString,
	                                       const bool               storeMassAndTPrime,
	                                       const std::string&       massVariableName,
	                                       const std::string&       tPrimeVariableName)
	{
		TFile* outputFile = rpwa::py::convertFromPy<TFile*>(outFilePy);
		if(not outputFile) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for outputFile when executing rpwa::importanceSampler::initializeFileWriter(...)");
			bp::throw_error_already_set();
		}

		return self.initializeFileWriter(outputFile, userString, storeMassAndTPrime, massVariableName, tPrimeVariableName);
	}


	void
	importanceSampler_setMassPrior(rpwa::importanceSampler& self,
	                               bp::object&              pyPrior)
	{
		TF1* prior = 0;
		if(not pyPrior.is_none()) {
			prior = rpwa::py::convertFromPy<TF1*>(pyPrior.ptr());
			if(not prior) {
				PyErr_SetString(PyExc_TypeError, "Got invalid input for prior when executing rpwa::importanceSampler::setMassPrior(...)");
				bp::throw_error_already_set();
			}
		}

		self.setMassPrior(prior);
	}


	void
	importanceSampler_printFuncInfo(rpwa::importanceSampler& self)
	{
		self.printFuncInfo(std::cout);
	}


	double
	importanceSampler_LogLikelihood(rpwa::importanceSampler& self,
	                                const bp::list&          pyParameters)
	{
		std::vector<double> parameters;
		if(not rpwa::py::convertBPObjectToVector<double>(pyParameters, parameters)) {
			PyErr_SetString(PyExc_TypeError, "invalid parameters gotten");
			bp::throw_error_already_set();
		}

		return self.LogLikelihood(parameters);
	}


	double
	importanceSampler_LogAPrioriProbability(rpwa::importanceSampler& self,
	                                        const bp::list&          pyParameters)
	{
		std::vector<double> parameters;
		if(not rpwa::py::convertBPObjectToVector<double>(pyParameters, parameters)) {
			PyErr_SetString(PyExc_TypeError, "invalid parameters gotten");
			bp::throw_error_already_set();
		}

		return self.LogAPrioriProbability(parameters);
	}


	void
	importanceSampler_SetPrecision(rpwa::importanceSampler& self,
	                               const std::string&       value = "medium")
	{
		if (value == "medium") {
			self.SetPrecision(BCEngineMCMC::kMedium);
		} else {
			PyErr_SetString(PyExc_TypeError, "unknown precision value");
			bp::throw_error_already_set();
		}
	}


	int
	importanceSampler_MarginalizeAll(rpwa::importanceSampler& self)
	{
		return self.MarginalizeAll(BCIntegrate::kMargMetropolis);
	}


	bp::list
	importanceSampler_FindMode(rpwa::importanceSampler& self)
	{
		return bp::list(self.FindMode());
	}


	void
	importanceSampler_WriteMarkovChain(rpwa::importanceSampler& self,
	                                   const std::string&       fileName,
	                                   const std::string&       option,
	                                   const bool               flag_run,
	                                   const bool               flag_prerun)
	{
		self.WriteMarkovChain(fileName, option, flag_run, flag_prerun);
	}


	std::string
	importanceSampler_GetSafeName(rpwa::importanceSampler& self)
	{
		std::string retVal = self.GetSafeName();
		return retVal;
	}


	void
	importanceSampler_SetInitialPositions(rpwa::importanceSampler& self,
	                                      const bp::list&          positionPy)
	{
		std::vector<double> position;
		if (not rpwa::py::convertBPObjectToVector<double>(positionPy, position)) {
			PyErr_SetString(PyExc_TypeError, "invalid startign points gotten");
			bp::throw_error_already_set();
		}

		self.SetInitialPositions(position);
	}


}


void rpwa::py::exportImportanceSampler() {

	bp::class_<rpwa::importanceSampler>("importanceSampler", bp::init<rpwa::modelIntensityPtr,
	                                                                  rpwa::beamAndVertexGeneratorPtr,
	                                                                  rpwa::massAndTPrimePickerPtr,
	                                                                  const rpwa::Beam&,
	                                                                  const rpwa::Target&,
	                                                                  const rpwa::FinalState&>())

		.def(
			"initializeFileWriter"
			, &::importanceSampler_initializeFileWriter
			, (bp::arg("outFile"),
			   bp::arg("userString") = "importanceSampledEvents",
			   bp::arg("storeMassAndTPrime") = true,
			   bp::arg("massVariableName") = "mass",
			   bp::arg("tPrimeVariableName") = "tPrime")
		)
		.def("finalizeFileWriter", &rpwa::importanceSampler::finalizeFileWriter)

		.def(
			"setPhaseSpaceOnly"
			, &rpwa::importanceSampler::setPhaseSpaceOnly
			, (bp::arg("inputValue") = true)
		)
		.def(
			"setMassPrior"
			, &importanceSampler_setMassPrior
			, bp::with_custodian_and_ward<1,2>()
			, (bp::arg("prior") = boost::python::object())
		)

		.def("nCalls", &rpwa::importanceSampler::nCalls)
		.def("resetFuncInfo", &rpwa::importanceSampler::resetFuncInfo)
		.def("printFuncInfo", &importanceSampler_printFuncInfo)

		// From here BAT stuff. Therefore names start with capital letters
		.def(
			"LogLikelihood"
			, &::importanceSampler_LogLikelihood
			, (bp::arg("parameters"))
		)
		.def(
			"LogAPrioriProbability"
			, &::importanceSampler_LogAPrioriProbability
			, (bp::arg("parameters"))
		)
		.def("FindMode", &::importanceSampler_FindMode)
		.def(
			"SetPrecision"
			, &::importanceSampler_SetPrecision
			, (bp::arg("value") = "medium")
		)
		.def(
			"PrintAllMarginalized"
			, &rpwa::importanceSampler::PrintAllMarginalized
			, (bp::arg("fileName"),
			   bp::arg("hdiv") = 1,
			   bp::arg("vdiv") = 1)
		)
		.def("PrintSummary", &rpwa::importanceSampler::PrintSummary)
		.def("MarginalizeAll", &::importanceSampler_MarginalizeAll)
		.def(
			"WriteMarkovChain"
			, &::importanceSampler_WriteMarkovChain
			, (bp::arg("fileName"),
			   bp::arg("option") = "CREATE",
			   bp::arg("flag_run") = true,
			   bp::arg("flag_prerun") = true)
		)
		.def("GetSafeName", &::importanceSampler_GetSafeName)
		.def(
			"SetNIterationsRun"
			, &rpwa::importanceSampler::SetNIterationsRun
			, (bp::arg("nIterations"))
		)
		.def(
			"SetNChains"
			, &rpwa::importanceSampler::SetNChains
			, (bp::arg("nChains"))
		)
		.def(
			"PrintCorrelationMatrix"
			, &rpwa::importanceSampler::PrintCorrelationMatrix
			, (bp::arg("fileName"))
		)
		.def(
			"SetRandomSeed"
			, &rpwa::importanceSampler::SetRandomSeed
			, (bp::arg("seed"))
		)
		.def(
			"SetInitialPositions"
			, &::importanceSampler_SetInitialPositions
			, (bp::arg("position"))
		)
		.def(
			"SetNLag"
			, &rpwa::importanceSampler::SetNLag
			, (bp::arg("nLag"))
		)

	;

	bp::register_ptr_to_python<rpwa::importanceSamplerPtr>();

}
