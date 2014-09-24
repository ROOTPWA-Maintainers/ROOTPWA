#include "waveDescription_py.h"

#include "amplitudeMetadata.h"

namespace bp = boost::python;

namespace {

	std::string waveDescription_printKeyFileContent(const rpwa::waveDescription& self) {
		std::stringstream sstr;
		self.printKeyFileContent(sstr);
		return sstr.str();
	}

	bp::tuple waveDescription_constructDecayTopology(const rpwa::waveDescription& self, bool fromTemplate = false) {
		rpwa::isobarDecayTopologyPtr topo;
		bool result = self.constructDecayTopology(topo, fromTemplate);
		return bp::make_tuple(result, topo);
	}

	bp::tuple waveDescription_constructAmplitude1(const rpwa::waveDescription& self) {
		rpwa::isobarAmplitudePtr amplitude;
		bool result = self.constructAmplitude(amplitude);
		return bp::make_tuple(result, amplitude);
	}

	bp::tuple waveDescription_constructAmplitude2(const rpwa::waveDescription& self, rpwa::isobarDecayTopologyPtr& topo) {
		rpwa::isobarAmplitudePtr amplitude;
		bool result = self.constructAmplitude(amplitude, topo);
		return bp::make_tuple(result, amplitude);
	}

	bool waveDescription_writeKeyFile(const std::string& keyFileName,
	                                  const bp::object   pyTopoOrAmp,
							          const bool         writeProdVert = true)
	{
		bp::extract<rpwa::isobarDecayTopology> get_iDT(pyTopoOrAmp);
		if(get_iDT.check()) {
			rpwa::isobarDecayTopology topoOrAmpiDT = get_iDT();
			return rpwa::waveDescription::writeKeyFile(keyFileName, topoOrAmpiDT, writeProdVert);
		}
		rpwa::isobarAmplitude* topoOrAmp = bp::extract<rpwa::isobarAmplitude*>(pyTopoOrAmp);
		return rpwa::waveDescription::writeKeyFile(keyFileName, *topoOrAmp, writeProdVert);
	}

	int waveDescription_Write(const rpwa::waveDescription& self, const char* name = 0) {
		return self.Write(name);
	}

}

void rpwa::py::exportWaveDescription() {

	bp::class_<rpwa::waveDescription>("waveDescription")

		.def(bp::init<const rpwa::amplitudeMetadata*>())
		.def("parseKeyFile", &rpwa::waveDescription::parseKeyFile)
		.def("parseKeyFileContent", &rpwa::waveDescription::parseKeyFileContent)
		.def("keyFileParsed", &rpwa::waveDescription::keyFileParsed)

		.def("keyFileContent", &rpwa::waveDescription::keyFileContent)
		.def("printKeyFileContent", &waveDescription_printKeyFileContent)

		.def(
			"constructDecayTopology"
			, &waveDescription_constructDecayTopology
			, bp::arg("fromTemplate")=false
		)

		.def("constructAmplitude", &waveDescription_constructAmplitude1)
		.def("constructAmplitude", &waveDescription_constructAmplitude2)

		.def(
			"writeKeyFile"
			, &waveDescription_writeKeyFile
			, (bp::arg("keyFileName"),
			   bp::arg("topoOrAmp"),
			   bp::arg("writeProdVert")=false)
		)
		.staticmethod("writeKeyFile")

		.def(
			"waveNameFromTopology"
			, &rpwa::waveDescription::waveNameFromTopology
			, (bp::arg("topo"),
			   bp::arg("newConvention")=false,
			   bp::arg("currentVertex")=rpwa::isobarDecayVertexPtr())
		)
		.staticmethod("waveNameFromTopology")

		.def(
			"waveLaTeXFromTopology"
			, &rpwa::waveDescription::waveLaTeXFromTopology
			, (bp::arg("topo"),
			   bp::arg("currentVertex")=rpwa::isobarDecayVertexPtr())
		)
		.staticmethod("waveLaTeXFromTopology")

		.def("Write", &waveDescription_Write, bp::arg("name")=0)

		.add_static_property("debugWaveDescription", &rpwa::waveDescription::debug, &rpwa::waveDescription::setDebug);

}
