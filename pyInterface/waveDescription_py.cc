#include "waveDescription_py.h"

namespace bp = boost::python;

namespace {

	struct waveDescriptionWrapper : public rpwa::waveDescription,
	                                       bp::wrapper<rpwa::waveDescription>
	{

		waveDescriptionWrapper()
			: rpwa::waveDescription(),
			  bp::wrapper<rpwa::waveDescription>() { };

		waveDescriptionWrapper(const rpwa::waveDescription& waveDesc)
			: rpwa::waveDescription(waveDesc),
			  bp::wrapper<rpwa::waveDescription>() { };

		std::string printKeyFileContents__() const {
			std::stringstream sstr;
			printKeyFileContents(sstr);
			return sstr.str();
		};

		bool constructAmplitude__1(rpwa::isobarAmplitudePtr& amplitude) const {
			return rpwa::waveDescription::constructAmplitude(amplitude);
		};

		bool constructAmplitude__2(rpwa::isobarAmplitudePtr& amplitude, rpwa::isobarDecayTopologyPtr& topo) const {
			return rpwa::waveDescription::constructAmplitude(amplitude, topo);
		};

		static bool writeKeyFile__(const std::string& keyFileName,
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
		};

	};

}

void rpwa::py::exportWaveDescription() {

	bp::class_<waveDescriptionWrapper>("waveDescription")

		.def("parseKeyFile", &waveDescriptionWrapper::parseKeyFile)
		.def("keyFileParsed", &waveDescriptionWrapper::keyFileParsed)

		.def("keyFileContents", &waveDescriptionWrapper::keyFileContents)
		.def("printKeyFileContents", &waveDescriptionWrapper::printKeyFileContents__)

		.def(
			"constructDecayTopology"
			, &waveDescriptionWrapper::constructDecayTopology
			, (bp::arg("topo"),
			   bp::arg("fromTemplate")=false)
		)

		.def("constructAmplitude", &waveDescriptionWrapper::constructAmplitude__1)
		.def("constructAmplitude", &waveDescriptionWrapper::constructAmplitude__2)

		.def(
			"writeKeyFile"
			, &waveDescriptionWrapper::writeKeyFile__
			, (bp::arg("keyFileName"),
			   bp::arg("topoOrAmp"),
			   bp::arg("writeProdVert")=false)
		)

		.def(
			"waveNameFromTopology"
			, &waveDescriptionWrapper::waveNameFromTopology
			, (bp::arg("topo"),
			   bp::arg("newConvention")=false,
			   bp::arg("currentVertex")=isobarDecayVertexPtr())
		)

		.def(
			"waveLaTeXFromTopology"
			, &waveDescriptionWrapper::waveLaTeXFromTopology
			, (bp::arg("topo"),
			   bp::arg("currentVertex")=isobarDecayVertexPtr())
		)

		.add_static_property("debugWaveDescription", &waveDescriptionWrapper::debug, &waveDescriptionWrapper::setDebug);

};

