///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014-2016 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      implementation of the master class of the resonance fit
//
//-------------------------------------------------------------------------


#include "resonanceFit.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include <TFormula.h>

#include <reportingUtils.hpp>

#include "components.h"
#include "forward.h"
#include "fsmd.h"
#include "input.h"
#include "model.h"
#include "parameters.h"


namespace {


	void
	writeFitQuality(YAML::Emitter& yamlOutput,
	                const std::map<std::string, double>& fitQuality)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing 'fitquality'." << std::endl;
		}

		yamlOutput << YAML::Key << "fitquality";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;
		for(std::map<std::string, double>::const_iterator it = fitQuality.begin(); it != fitQuality.end(); ++it) {
			yamlOutput << YAML::Key << it->first;
			yamlOutput << YAML::Value << it->second;
		}
		yamlOutput << YAML::EndMap;
	}


	void
	writeFreeParameters(YAML::Emitter& yamlOutput,
	                    const std::vector<std::string>& freeParameters)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing 'freeparameters'." << std::endl;
		}

		if(freeParameters.size() > 0) {
			yamlOutput << YAML::Key << "freeparameters";
			yamlOutput << YAML::Value << freeParameters;
		}
	}


	void
	writeInput(YAML::Emitter& yamlOutput,
	           const rpwa::resonanceFit::inputConstPtr& fitInput)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing 'input'." << std::endl;
		}

		yamlOutput << YAML::Key << "input";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginSeq;

		for(size_t idxBin = 0; idxBin < fitInput->nrBins(); ++idxBin) {
			const rpwa::resonanceFit::input::bin& bin = fitInput->getBin(idxBin);

			if(rpwa::resonanceFit::debug()) {
				printDebug << "writing bin " << idxBin << " to 'input'." << std::endl;
			}

			yamlOutput << YAML::BeginMap;

			yamlOutput << YAML::Key << "name";
			yamlOutput << YAML::Value << bin.fileName();

			yamlOutput << YAML::Key << "tPrimeMean";
			yamlOutput << YAML::Value << bin.tPrimeMean();

			if(bin.rescaleErrors() != 1.) {
				yamlOutput << YAML::Key << "rescaleErrors";
				yamlOutput << YAML::Value << bin.rescaleErrors();
			}

			if(bin.sysFileNames().size() > 0) {
				yamlOutput << YAML::Key << "systematics";
				yamlOutput << YAML::Value << bin.sysFileNames();
			}

			// check whether the same waves and ranges have already been written
			size_t idxAnchor = 0;
			for( ; idxAnchor < idxBin; ++idxAnchor) {
				// continue if number of waves does not match
				if(bin.nrWaves() != fitInput->getBin(idxAnchor).nrWaves())
					continue;

				// compare the individual waves
				bool match = true;
				for(size_t idxWave = 0; idxWave < bin.nrWaves(); ++idxWave) {
					if(bin.getWave(idxWave).waveName() != fitInput->getBin(idxAnchor).getWave(idxWave).waveName()) {
						match = false;
						break;
					}
					if(bin.getWave(idxWave).massLimits() != fitInput->getBin(idxAnchor).getWave(idxWave).massLimits()) {
						match = false;
						break;
					}
				}

				// all wave names and mass ranges match
				if(match)
					break;
			}

			yamlOutput << YAML::Key << "waves";
			yamlOutput << YAML::Value;

			std::ostringstream anchorName;
			anchorName << "waveset" << idxAnchor;

			if(idxAnchor == idxBin) {
				// this combination of waves and mass ranges was not yet encountered
				yamlOutput << YAML::Anchor(anchorName.str());
				yamlOutput << YAML::BeginSeq;
				for(std::vector<rpwa::resonanceFit::input::bin::wave>::const_iterator wave = bin.waves().begin(); wave != bin.waves().end(); ++wave) {
					yamlOutput << YAML::BeginMap;

					yamlOutput << YAML::Key << "name";
					yamlOutput << YAML::Value << wave->waveName();

					if(wave->massLimits().first >= 0) {
						yamlOutput << YAML::Key << "massLower";
						yamlOutput << YAML::Value << wave->massLimits().first;
					}
					if(wave->massLimits().second >= 0) {
						yamlOutput << YAML::Key << "massUpper";
						yamlOutput << YAML::Value << wave->massLimits().second;
					}

					yamlOutput << YAML::EndMap;
				}
				yamlOutput << YAML::EndSeq;
			} else {
				yamlOutput << YAML::Alias(anchorName.str());
			}

			yamlOutput << YAML::EndMap;
		}

		yamlOutput << YAML::EndSeq;
	}


	void
	writeModelAnchors(YAML::Emitter& yamlOutput,
	                  const std::vector<std::string>& anchorWaveNames,
	                  const std::vector<std::string>& anchorComponentNames)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing 'anchorwave'." << std::endl;
		}

		yamlOutput << YAML::Key << "anchorwave";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginSeq;

		for(size_t idxBin = 0; idxBin < anchorWaveNames.size(); ++idxBin) {
			// check whether the same wave and component names have already been written
			size_t idxAnchor = 0;
			for( ; idxAnchor < idxBin; ++idxAnchor) {
				if(anchorWaveNames[idxBin] != anchorWaveNames[idxAnchor]) {
					continue;
				}
				if(anchorComponentNames[idxBin] != anchorComponentNames[idxAnchor]) {
					continue;
				}

				// anchor wave and component names match
				break;
			}

			std::ostringstream anchorName;
			anchorName << "anchor" << idxAnchor;

			if(idxAnchor == idxBin) {
				yamlOutput << YAML::Anchor(anchorName.str());

				yamlOutput << YAML::Flow;
				yamlOutput << YAML::BeginSeq;

				yamlOutput << YAML::Value << anchorWaveNames[idxBin];
				yamlOutput << YAML::Value << anchorComponentNames[idxBin];

				yamlOutput << YAML::EndSeq;
				yamlOutput << YAML::Block;
			} else {
				yamlOutput << YAML::Alias(anchorName.str());
			}
		}

		yamlOutput << YAML::EndSeq;
	}


	void
	writeModelParameter(YAML::Emitter& yamlOutput,
	                    const rpwa::resonanceFit::parameter& parameter,
	                    const double value,
	                    const double error)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing parameter '" << parameter.name() << "'." << std::endl;
		}

		yamlOutput << YAML::Key << parameter.name();
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "val";
		yamlOutput << YAML::Value << value;

		yamlOutput << YAML::Key << "error";
		yamlOutput << YAML::Value << error;

		if(parameter.limitedLower()) {
			yamlOutput << YAML::Key << "lower";
			yamlOutput << YAML::Value << parameter.limitLower();
		}

		if(parameter.limitedUpper()) {
			yamlOutput << YAML::Key << "upper";
			yamlOutput << YAML::Value << parameter.limitUpper();
		}

		yamlOutput << YAML::Key << "step";
		yamlOutput << YAML::Value << parameter.step();

		yamlOutput << YAML::Key << "fix";
		yamlOutput << YAML::Value << parameter.fixed();

		yamlOutput << YAML::EndMap;
	}


	template<typename T>
	void
	writeModelComponentDecayChannelBranchingRatio(YAML::Emitter& yamlOutput,
	                                              const std::shared_ptr<const T>& component,
	                                              const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "branchingRatio";
		yamlOutput << YAML::Value << component->branchingRatio()[idxDecayChannel];
	}


	template<typename T>
	void
	writeModelComponentDecayChannelExponent(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "exponent";
		yamlOutput << YAML::Value << component->exponent();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelIntegral(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "integral";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::Flow;
		yamlOutput << YAML::BeginSeq;

		const std::vector<double> masses =  component->masses();
		const std::vector<double> values =  component->values();
		for(size_t idx = 0; idx < masses.size(); ++idx) {
			yamlOutput << YAML::BeginSeq;
			yamlOutput << masses[idx];
			yamlOutput << values[idx];
			yamlOutput << YAML::EndSeq;
		}

		yamlOutput << YAML::EndSeq;
		yamlOutput << YAML::Block;
	}


	template<typename T>
	void
	writeModelComponentDecayChannelIntegral(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component,
	                                        const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "integral";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::Flow;
		yamlOutput << YAML::BeginSeq;

		const std::vector<double> masses =  component->masses()[idxDecayChannel];
		const std::vector<double> values =  component->values()[idxDecayChannel];
		for(size_t idx = 0; idx < masses.size(); ++idx) {
			yamlOutput << YAML::BeginSeq;
			yamlOutput << masses[idx];
			yamlOutput << values[idx];
			yamlOutput << YAML::EndSeq;
		}

		yamlOutput << YAML::EndSeq;
		yamlOutput << YAML::Block;
	}


	template<typename T>
	void
	writeModelComponentDecayChannelMIsobar1(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "mIsobar1";
		yamlOutput << YAML::Value << component->mIsobar1();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelMIsobar1(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component,
	                                        const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "mIsobar1";
		yamlOutput << YAML::Value << component->mIsobar1()[idxDecayChannel];
	}


	template<typename T>
	void
	writeModelComponentDecayChannelMIsobar2(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "mIsobar2";
		yamlOutput << YAML::Value << component->mIsobar2();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelMIsobar2(YAML::Emitter& yamlOutput,
	                                        const std::shared_ptr<const T>& component,
	                                        const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "mIsobar2";
		yamlOutput << YAML::Value << component->mIsobar2()[idxDecayChannel];
	}


	template<typename T>
	void
	writeModelComponentDecayChannelRelAngularMom(YAML::Emitter& yamlOutput,
	                                             const std::shared_ptr<const T>& component)
	{
		yamlOutput << YAML::Key << "relAngularMom";
		yamlOutput << YAML::Value << component->relAngularMom();
	}


	template<typename T>
	void
	writeModelComponentDecayChannelRelAngularMom(YAML::Emitter& yamlOutput,
	                                             const std::shared_ptr<const T>& component,
	                                             const size_t idxDecayChannel)
	{
		yamlOutput << YAML::Key << "relAngularMom";
		yamlOutput << YAML::Value << component->relAngularMom()[idxDecayChannel];
	}


	void
	writeModelComponent(YAML::Emitter& yamlOutput,
	                    const rpwa::resonanceFit::componentConstPtr& component,
	                    const rpwa::resonanceFit::parameters& fitParameters,
	                    const rpwa::resonanceFit::parameters& fitParametersError)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing component and its parameters." << std::endl;
		}

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "name";
		yamlOutput << YAML::Value << component->getName();

		yamlOutput << YAML::Key << "type";
		yamlOutput << YAML::Value << component->getType();

		for(size_t idxParameter = 0; idxParameter < component->getNrParameters(); ++idxParameter) {
			writeModelParameter(yamlOutput,
			                    component->getParameter(idxParameter),
			                    fitParameters.getParameter(component->getId(), idxParameter),
			                    fitParametersError.getParameter(component->getId(), idxParameter));
		}

		yamlOutput << YAML::Key << "decaychannels";
		yamlOutput << YAML::Value;
		yamlOutput << YAML::BeginSeq;

		for(size_t idxDecayChannel = 0; idxDecayChannel < component->getNrChannels(); ++idxDecayChannel) {
			yamlOutput << YAML::BeginMap;

			yamlOutput << YAML::Key << "amp";
			yamlOutput << YAML::Value << component->getChannel(idxDecayChannel).getWaveName();

			if(component->mapCouplingToMasterChannel(component->mapChannelToCoupling(idxDecayChannel)) == idxDecayChannel) {
				yamlOutput << YAML::Key << "couplings";
				yamlOutput << YAML::Value;
				yamlOutput << YAML::BeginSeq;

				const std::vector<size_t>& bins = component->getChannel(idxDecayChannel).getBins();
				for(size_t i = 0; i < bins.size(); ++i) {
					const size_t idxBin = bins[i];
					yamlOutput << YAML::Flow;
					yamlOutput << YAML::BeginSeq;

					yamlOutput << fitParameters.getCoupling(component->getId(), component->mapChannelToCoupling(idxDecayChannel), idxBin).real();
					yamlOutput << fitParameters.getCoupling(component->getId(), component->mapChannelToCoupling(idxDecayChannel), idxBin).imag();

					yamlOutput << YAML::EndSeq;
					yamlOutput << YAML::Block;
				}

				yamlOutput << YAML::EndSeq;
			}

			if(component->getNrBranchings() > 1) {
				if(component->mapBranchingToMasterChannel(component->mapChannelToBranching(idxDecayChannel)) == idxDecayChannel) {
					yamlOutput << YAML::Key << "branching";
					yamlOutput << YAML::Value;

					yamlOutput << YAML::Flow;
					yamlOutput << YAML::BeginSeq;

					yamlOutput << fitParameters.getBranching(component->getId(), component->mapChannelToBranching(idxDecayChannel)).real();
					yamlOutput << fitParameters.getBranching(component->getId(), component->mapChannelToBranching(idxDecayChannel)).imag();

					yamlOutput << YAML::EndSeq;
					yamlOutput << YAML::Block;
				}
			}

			if(std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component)) {
				writeModelComponentDecayChannelMIsobar1(yamlOutput,
				                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
				                                        idxDecayChannel);
				writeModelComponentDecayChannelMIsobar2(yamlOutput,
				                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
				                                        idxDecayChannel);
				writeModelComponentDecayChannelRelAngularMom(yamlOutput,
				                                             std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
				                                             idxDecayChannel);
				writeModelComponentDecayChannelBranchingRatio(yamlOutput,
				                                              std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
				                                              idxDecayChannel);
			}
			if(std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component)) {
				writeModelComponentDecayChannelIntegral(yamlOutput,
				                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component),
				                                        idxDecayChannel);
				writeModelComponentDecayChannelBranchingRatio(yamlOutput,
				                                              std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component),
				                                              idxDecayChannel);
			}

			yamlOutput << YAML::EndMap;
		}

		yamlOutput << YAML::EndSeq;

		if(component->getTotalNrChannels() > component->getNrChannels()) {
			yamlOutput << YAML::Key << "extradecaychannels";
			yamlOutput << YAML::Value;
			yamlOutput << YAML::BeginSeq;

			for(size_t idxDecayChannel = component->getNrChannels(); idxDecayChannel < component->getTotalNrChannels(); ++idxDecayChannel) {
				yamlOutput << YAML::BeginMap;

				if(std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component)) {
					writeModelComponentDecayChannelMIsobar1(yamlOutput,
					                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
					                                        idxDecayChannel);
					writeModelComponentDecayChannelMIsobar2(yamlOutput,
					                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
					                                        idxDecayChannel);
					writeModelComponentDecayChannelRelAngularMom(yamlOutput,
					                                             std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
					                                             idxDecayChannel);
					writeModelComponentDecayChannelBranchingRatio(yamlOutput,
					                                              std::dynamic_pointer_cast<const rpwa::resonanceFit::dynamicWidthBreitWigner>(component),
					                                              idxDecayChannel);
				}
				if(std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component)) {
					writeModelComponentDecayChannelIntegral(yamlOutput,
					                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component),
					                                        idxDecayChannel);
					writeModelComponentDecayChannelBranchingRatio(yamlOutput,
					                                              std::dynamic_pointer_cast<const rpwa::resonanceFit::integralWidthBreitWigner>(component),
					                                              idxDecayChannel);
				}

				yamlOutput << YAML::EndMap;
			}

			yamlOutput << YAML::EndSeq;
		}

		if(std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component)) {
			writeModelComponentDecayChannelMIsobar1(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component));
			writeModelComponentDecayChannelMIsobar2(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component));
			writeModelComponentDecayChannelRelAngularMom(yamlOutput,
			                                             std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component));
			writeModelComponentDecayChannelExponent(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackground>(component));
		}
		if(std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component)) {
			writeModelComponentDecayChannelMIsobar1(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component));
			writeModelComponentDecayChannelMIsobar2(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component));
			writeModelComponentDecayChannelRelAngularMom(yamlOutput,
			                                             std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component));
			writeModelComponentDecayChannelExponent(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackground>(component));
		}
		if(std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackgroundIntegral>(component)) {
			writeModelComponentDecayChannelIntegral(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackgroundIntegral>(component));
			writeModelComponentDecayChannelExponent(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::exponentialBackgroundIntegral>(component));
		}
		if(std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(component)) {
			writeModelComponentDecayChannelIntegral(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(component));
			writeModelComponentDecayChannelExponent(yamlOutput,
			                                        std::dynamic_pointer_cast<const rpwa::resonanceFit::tPrimeDependentBackgroundIntegral>(component));
		}

		yamlOutput << YAML::EndMap;
	}


	void
	writeModelComponents(YAML::Emitter& yamlOutput,
	                     const std::vector<rpwa::resonanceFit::componentConstPtr>& components,
	                     const rpwa::resonanceFit::parameters& fitParameters,
	                     const rpwa::resonanceFit::parameters& fitParametersError)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing 'components'." << std::endl;
		}

		yamlOutput << YAML::Key << "components";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginSeq;

		for(size_t idxComponent = 0; idxComponent < components.size(); ++idxComponent) {
			writeModelComponent(yamlOutput,
			                    components[idxComponent],
			                    fitParameters,
			                    fitParametersError);
		}

		yamlOutput << YAML::EndSeq;
	}


	void
	writeModelFsmdBin(YAML::Emitter& yamlOutput,
	                  const rpwa::resonanceFit::fsmdConstPtr& fsmd,
	                  const size_t idxBin,
	                  const rpwa::resonanceFit::parameters& fitParameters,
	                  const rpwa::resonanceFit::parameters& fitParametersError)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing 'finalStateMassDependence' for an individual bin." << std::endl;
		}

		yamlOutput << YAML::BeginMap;

		yamlOutput << YAML::Key << "formula";
		yamlOutput << YAML::Value << fsmd->getFunction(idxBin)->GetTitle();

		for(size_t idxParameter = 0; idxParameter < fsmd->getNrParameters(idxBin); ++idxParameter) {
			writeModelParameter(yamlOutput,
			                    fsmd->getParameter(idxBin, idxParameter),
			                    fitParameters.getParameter(fsmd->getId(), fsmd->getParameterIndex(idxBin)+idxParameter),
			                    fitParametersError.getParameter(fsmd->getId(), fsmd->getParameterIndex(idxBin)+idxParameter));
		}

		yamlOutput << YAML::EndMap;
	}


	void
	writeModelFsmd(YAML::Emitter& yamlOutput,
	               const rpwa::resonanceFit::fsmdConstPtr& fsmd,
	               const rpwa::resonanceFit::parameters& fitParameters,
	               const rpwa::resonanceFit::parameters& fitParametersError)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing 'finalStateMassDependence'." << std::endl;
		}

		yamlOutput << YAML::Key << "finalStateMassDependence";
		yamlOutput << YAML::Value;

		if(not fsmd->isSameFunctionForAllBins()) {
			yamlOutput << YAML::BeginSeq;
		}

		const size_t maxNrBins = fsmd->isSameFunctionForAllBins() ? 1 : fsmd->getNrBins();
		for(size_t idxBin = 0; idxBin < maxNrBins; ++idxBin) {
			writeModelFsmdBin(yamlOutput,
			                  fsmd,
			                  idxBin,
			                  fitParameters,
			                  fitParametersError);
		}

		if(not fsmd->isSameFunctionForAllBins()) {
			yamlOutput << YAML::EndSeq;
		}
	}


	void
	writeModel(YAML::Emitter& yamlOutput,
	           const rpwa::resonanceFit::modelConstPtr& fitModel,
	           const rpwa::resonanceFit::parameters& fitParameters,
	           const rpwa::resonanceFit::parameters& fitParametersError)
	{
		yamlOutput << YAML::Key << "model";
		yamlOutput << YAML::Value;

		yamlOutput << YAML::BeginMap;

		writeModelAnchors(yamlOutput,
		                  fitModel->anchorWaveNames(),
		                  fitModel->anchorComponentNames());

		std::vector<rpwa::resonanceFit::componentConstPtr> components;
		for(size_t idxComponent = 0; idxComponent < fitModel->getNrComponents(); ++idxComponent) {
			components.push_back(fitModel->getComponent(idxComponent));
		}
		writeModelComponents(yamlOutput,
		                     components,
		                     fitParameters,
		                     fitParametersError);

		writeModelFsmd(yamlOutput,
		               fitModel->getFsmd(),
		               fitParameters,
		               fitParametersError);

		yamlOutput << YAML::EndMap;
	}


	void
	writeConfig(YAML::Emitter& yamlOutput,
	            const rpwa::resonanceFit::inputConstPtr& fitInput,
	            const rpwa::resonanceFit::modelConstPtr& fitModel,
	            const rpwa::resonanceFit::parameters& fitParameters,
	            const rpwa::resonanceFit::parameters& fitParametersError,
	            const std::map<std::string, double>& fitQuality,
	            const std::vector<std::string>& freeParameters)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "writing configuration file." << std::endl;
		}

		yamlOutput << YAML::BeginMap;

		writeFitQuality(yamlOutput, fitQuality);
		writeFreeParameters(yamlOutput, freeParameters);
		writeInput(yamlOutput, fitInput);
		writeModel(yamlOutput, fitModel, fitParameters, fitParametersError);

		yamlOutput << YAML::EndMap;
	}


}


void
rpwa::resonanceFit::writeConfig(const std::string& configFileName,
                                const rpwa::resonanceFit::inputConstPtr& fitInput,
                                const rpwa::resonanceFit::modelConstPtr& fitModel,
                                const rpwa::resonanceFit::parameters& fitParameters,
                                const rpwa::resonanceFit::parameters& fitParametersError,
                                const std::map<std::string, double>& fitQuality,
                                const std::vector<std::string>& freeParameters)
{
	std::ofstream configFile(configFileName);

	YAML::Emitter yamlOutput(configFile);
	::writeConfig(yamlOutput,
	              fitInput,
	              fitModel,
	              fitParameters,
	              fitParametersError,
	              fitQuality,
	              freeParameters);

	// newline at end-of-file
	configFile << std::endl;

	configFile.close();
}
