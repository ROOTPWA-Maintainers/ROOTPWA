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

#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TMultiGraph.h>

#include <reportingUtils.hpp>

#include "components.h"
#include "data.h"
#include "fsmd.h"
#include "information.h"
#include "model.h"


namespace {


	void
	createPlotsWave(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                const rpwa::resonanceFit::dataConstPtr& fitData,
	                const rpwa::resonanceFit::modelConstPtr& fitModel,
	                const rpwa::resonanceFit::parameters& fitParameters,
	                rpwa::resonanceFit::cache& cache,
	                TDirectory* outDirectory,
	                const bool rangePlotting,
	                const size_t extraBinning,
	                const size_t idxWave,
	                const size_t idxBin)
	{
		const rpwa::resonanceFit::information::wave& wave = fitInformation->getWave(idxWave);
		if(rpwa::resonanceFit::debug()) {
			printDebug << "start creating plots for wave '" << wave.waveName() << "' in bin " << idxBin << "." << std::endl;
		}

		TMultiGraph graphs;
		graphs.SetName(wave.waveName().c_str());
		graphs.SetTitle(wave.waveName().c_str());

		TGraphErrors* systematics = NULL;
		if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
			systematics = new TGraphErrors;
			systematics->SetName((wave.waveName() + "__sys").c_str());
			systematics->SetTitle((wave.waveName() + "__sys").c_str());
			systematics->SetLineColor(kAzure-9);
			systematics->SetFillColor(kAzure-9);
			graphs.Add(systematics, "2");
		}

		TGraphErrors* data = new TGraphErrors;
		data->SetName((wave.waveName() + "__data").c_str());
		data->SetTitle((wave.waveName() + "__data").c_str());
		graphs.Add(data, "P");

		TGraph* fit = new TGraph;
		fit->SetName((wave.waveName() + "__fit").c_str());
		fit->SetTitle((wave.waveName() + "__fit").c_str());
		fit->SetLineColor(kRed);
		fit->SetLineWidth(2);
		fit->SetMarkerColor(kRed);
		graphs.Add(fit, "L");

		TGraph* phaseSpace = new TGraph;
		phaseSpace->SetName((wave.waveName() + "__ps").c_str());
		phaseSpace->SetTitle((wave.waveName() + "__ps").c_str());
		graphs.Add(phaseSpace, "L");

		const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel->getComponentChannel(idxBin, idxWave);
		std::vector<TGraph*> components;
		for(size_t idxComponents = 0; idxComponents < compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			TGraph* component = new TGraph;
			component->SetName((wave.waveName() + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());
			component->SetTitle((wave.waveName() + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());

			Color_t color = kBlue;
			if(fitModel->getComponent(idxComponent)->getName().find("bkg") != std::string::npos) {
				color = kMagenta;
			}
			component->SetLineColor(color);
			component->SetMarkerColor(color);

			graphs.Add(component, "L");
			components.push_back(component);
		}

		// plot data
		double maxIE = -std::numeric_limits<double>::max();
		for(size_t point = 0; point <= (fitData->nrMassBins()[idxBin]-1); ++point) {
			const size_t idxMass = point;
			const double mass = fitData->massBinCenters()[idxBin][idxMass];
			const double halfBin = 0.5 * (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1);

			data->SetPoint(point, mass, fitData->plottingIntensities()[idxBin][idxMass][idxWave].first);
			data->SetPointError(point, halfBin, fitData->plottingIntensities()[idxBin][idxMass][idxWave].second);
			maxIE = std::max(maxIE, fitData->plottingIntensities()[idxBin][idxMass][idxWave].first+fitData->plottingIntensities()[idxBin][idxMass][idxWave].second);

			if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
				const double minSI = fitData->sysPlottingIntensities()[idxBin][idxMass][idxWave].first;
				const double maxSI = fitData->sysPlottingIntensities()[idxBin][idxMass][idxWave].second;
				systematics->SetPoint(point, mass, (maxSI+minSI)/2.);
				systematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);
				maxIE = std::max(maxIE, maxSI);
			}
		}

		// plot fit, either over full or limited mass range
		const size_t firstPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].first) : 0;
		const size_t lastPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
		for(size_t point = firstPoint; point <= lastPoint; ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			const double intensity = fitModel->intensity(fitParameters, cache, idxWave, idxBin, mass, idxMass);
			fit->SetPoint(point-firstPoint, mass, intensity);
			maxIE = std::max(maxIE, intensity);

			for(size_t idxComponents = 0; idxComponents < compChannel.size(); ++idxComponents) {
				const size_t idxComponent = compChannel[idxComponents].first;
				const size_t idxChannel = compChannel[idxComponents].second;

				std::complex<double> prodAmp = fitModel->getComponent(idxComponent)->val(fitParameters, cache, idxBin, mass, idxMass);
				prodAmp *= fitModel->getComponent(idxComponent)->getCouplingPhaseSpace(fitParameters, cache, idxChannel, idxBin, mass, idxMass);
				if(fitModel->getFsmd()) {
					prodAmp *= fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass);
				}

				components[idxComponents]->SetPoint(point-firstPoint, mass, norm(prodAmp));
				maxIE = std::max(maxIE, norm(prodAmp));
			}
		}

		boost::multi_array<double, 2>::const_array_view<1>::type viewM = fitData->massBinCenters()[boost::indices[idxBin][boost::multi_array<double, 2>::index_range(0, fitData->nrMassBins()[idxBin])]];
		boost::multi_array<double, 3>::const_array_view<1>::type viewInt = fitData->phaseSpaceIntegrals()[boost::indices[idxBin][boost::multi_array<double, 3>::index_range(0, fitData->nrMassBins()[idxBin])][idxWave]];
		ROOT::Math::Interpolator phaseSpaceInterpolator(std::vector<double>(viewM.begin(), viewM.end()), std::vector<double>(viewInt.begin(), viewInt.end()), ROOT::Math::Interpolation::kLINEAR);

		// plot phase-space
		double maxP = -std::numeric_limits<double>::max();
		for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			double ps = pow((idxMass != std::numeric_limits<size_t>::max()) ? fitData->phaseSpaceIntegrals()[idxBin][idxMass][idxWave] : phaseSpaceInterpolator.Eval(mass), 2);
			if(fitModel->getFsmd()) {
				ps *= std::norm(fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass));
			}
			phaseSpace->SetPoint(point, mass, ps);
			maxP = std::max(maxP, ps);
		}

		// scale phase-space graph to half-height of intensity graphs
		for(Int_t idx = 0; idx < phaseSpace->GetN(); ++idx) {
			double x, y;
			phaseSpace->GetPoint(idx, x, y);
			phaseSpace->SetPoint(idx, x, y * 0.5 * maxIE/maxP);
		}

		outDirectory->cd();
		graphs.Write();
	}


	void
	createPlotsWaveSum(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                   const rpwa::resonanceFit::dataConstPtr& fitData,
	                   const rpwa::resonanceFit::modelConstPtr& fitModel,
	                   const rpwa::resonanceFit::parameters& fitParameters,
	                   rpwa::resonanceFit::cache& cache,
	                   TDirectory* outDirectory,
	                   const bool rangePlotting,
	                   const size_t extraBinning,
	                   const size_t idxWave)
	{
		const rpwa::resonanceFit::information::wave& wave = fitInformation->getWave(idxWave);
		if(rpwa::resonanceFit::debug()) {
			printDebug << "start creating plots for wave '" << wave.waveName() << "' for sum over all bins." << std::endl;
		}

		// all mass binnings must be the same to be able to create the sum plots
		if(not fitData->hasSameMassBinning() or not fitModel->isMappingEqualInAllBins()) {
			printErr << "cannot create plots for wave '" << wave.waveName() << "' for sum over all bins if the bins used different mass binnings." << std::endl;
			throw;
		}
		const size_t idxBin = 0;

		TMultiGraph graphs;
		graphs.SetName(wave.waveName().c_str());
		graphs.SetTitle(wave.waveName().c_str());

		TGraphErrors* data = new TGraphErrors;
		data->SetName((wave.waveName() + "__data").c_str());
		data->SetTitle((wave.waveName() + "__data").c_str());
		graphs.Add(data, "P");

		TGraph* fit = new TGraph;
		fit->SetName((wave.waveName() + "__fit").c_str());
		fit->SetTitle((wave.waveName() + "__fit").c_str());
		fit->SetLineColor(kRed);
		fit->SetLineWidth(2);
		fit->SetMarkerColor(kRed);
		graphs.Add(fit, "L");

		const std::vector<std::pair<size_t, size_t> >& compChannel = fitModel->getComponentChannel(idxBin, idxWave);
		std::vector<TGraph*> components;
		for(size_t idxComponents = 0; idxComponents < compChannel.size(); ++idxComponents) {
			const size_t idxComponent = compChannel[idxComponents].first;
			TGraph* component = new TGraph;
			component->SetName((wave.waveName() + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());
			component->SetTitle((wave.waveName() + "__" + fitModel->getComponent(idxComponent)->getName()).c_str());

			Color_t color = kBlue;
			if(fitModel->getComponent(idxComponent)->getName().find("bkg") != std::string::npos) {
				color = kMagenta;
			}
			component->SetLineColor(color);
			component->SetMarkerColor(color);

			graphs.Add(component, "L");
			components.push_back(component);
		}

		// plot data
		for(size_t point = 0; point <= (fitData->nrMassBins()[idxBin]-1); ++point) {
			const size_t idxMass = point;
			const double mass = fitData->massBinCenters()[idxBin][idxMass];
			const double halfBin = 0.5 * (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1);

			double sum = 0.;
			double error2 = 0.;
			for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
				sum += fitData->plottingIntensities()[idxBin][idxMass][idxWave].first;
				error2 += std::pow(fitData->plottingIntensities()[idxBin][idxMass][idxWave].second, 2);
			}
			data->SetPoint(point, mass, sum);
			data->SetPointError(point, halfBin, sqrt(error2));
		}

		// plot fit, either over full or limited mass range
		const size_t firstPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].first) : 0;
		const size_t lastPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][idxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
		for(size_t point = firstPoint; point <= lastPoint; ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			double sum = 0.;
			for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
				sum += fitModel->intensity(fitParameters, cache, idxWave, idxBin, mass, idxMass);
			}
			fit->SetPoint(point-firstPoint, mass, sum);

			for(size_t idxComponents = 0; idxComponents < compChannel.size(); ++idxComponents) {
				const size_t idxComponent = compChannel[idxComponents].first;
				const size_t idxChannel = compChannel[idxComponents].second;

				double sum = 0.;
				for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
					std::complex<double> prodAmp = fitModel->getComponent(idxComponent)->val(fitParameters, cache, idxBin, mass, idxMass);
					prodAmp *= fitModel->getComponent(idxComponent)->getCouplingPhaseSpace(fitParameters, cache, idxChannel, idxBin, mass, idxMass);
					if(fitModel->getFsmd()) {
						prodAmp *= fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass);
					}
					sum += norm(prodAmp);
				}
				components[idxComponents]->SetPoint(point-firstPoint, mass, sum);
			}
		}

		outDirectory->cd();
		graphs.Write();
	}


	void
	createPlotsWavePair(const rpwa::resonanceFit::informationConstPtr& fitInformation,
	                    const rpwa::resonanceFit::dataConstPtr& fitData,
	                    const rpwa::resonanceFit::modelConstPtr& fitModel,
	                    const rpwa::resonanceFit::parameters& fitParameters,
	                    rpwa::resonanceFit::cache& cache,
	                    TDirectory* outDirectory,
	                    const bool rangePlotting,
	                    const size_t extraBinning,
	                    const size_t idxWave,
	                    const size_t jdxWave,
	                    const size_t idxBin)
	{
		const rpwa::resonanceFit::information::wave& waveI = fitInformation->getWave(idxWave);
		const rpwa::resonanceFit::information::wave& waveJ = fitInformation->getWave(jdxWave);
		if(rpwa::resonanceFit::debug()) {
			printDebug << "start creating plots for wave pair '" << waveI.waveName() << "' and '" << waveJ.waveName() << "' in bin " << idxBin << "." << std::endl;
		}

		const std::string realName = waveI.waveName() + "__" + waveJ.waveName() + "__real";
		const std::string imagName = waveI.waveName() + "__" + waveJ.waveName() + "__imag";
		const std::string phaseName = waveI.waveName() + "__" + waveJ.waveName() + "__phase";

		TMultiGraph real;
		real.SetName(realName.c_str());
		real.SetTitle(realName.c_str());

		TMultiGraph imag;
		imag.SetName(imagName.c_str());
		imag.SetTitle(imagName.c_str());

		TMultiGraph phase;
		phase.SetName(phaseName.c_str());
		phase.SetTitle(phaseName.c_str());

		TGraphErrors* realSystematics = NULL;
		TGraphErrors* imagSystematics = NULL;
		TGraphErrors* phaseSystematics = NULL;
		if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
			realSystematics = new TGraphErrors;
			realSystematics->SetName((realName + "__sys").c_str());
			realSystematics->SetTitle((realName + "__sys").c_str());
			realSystematics->SetLineColor(kAzure-9);
			realSystematics->SetFillColor(kAzure-9);
			real.Add(realSystematics, "2");

			imagSystematics = new TGraphErrors;
			imagSystematics->SetName((imagName + "__sys").c_str());
			imagSystematics->SetTitle((imagName + "__sys").c_str());
			imagSystematics->SetLineColor(kAzure-9);
			imagSystematics->SetFillColor(kAzure-9);
			imag.Add(imagSystematics, "2");

			phaseSystematics = new TGraphErrors;
			phaseSystematics->SetName((phaseName + "__sys").c_str());
			phaseSystematics->SetTitle((phaseName + "__sys").c_str());
			phaseSystematics->SetLineColor(kAzure-9);
			phaseSystematics->SetFillColor(kAzure-9);
			phase.Add(phaseSystematics, "2");
		}

		TGraphErrors* realData = new TGraphErrors;
		realData->SetName((realName + "__data").c_str());
		realData->SetTitle((realName + "__data").c_str());
		real.Add(realData, "P");

		TGraphErrors* imagData = new TGraphErrors;
		imagData->SetName((imagName + "__data").c_str());
		imagData->SetTitle((imagName + "__data").c_str());
		imag.Add(imagData, "P");

		TGraphErrors* phaseData = new TGraphErrors;
		phaseData->SetName((phaseName + "__data").c_str());
		phaseData->SetTitle((phaseName + "__data").c_str());
		phase.Add(phaseData, "P");

		TGraph* realFit = new TGraph;
		realFit->SetName((realName + "__fit").c_str());
		realFit->SetTitle((realName + "__fit").c_str());
		realFit->SetLineColor(kRed);
		realFit->SetLineWidth(2);
		realFit->SetMarkerColor(kRed);
		real.Add(realFit, "L");

		TGraph* imagFit = new TGraph;
		imagFit->SetName((imagName + "__fit").c_str());
		imagFit->SetTitle((imagName + "__fit").c_str());
		imagFit->SetLineColor(kRed);
		imagFit->SetLineWidth(2);
		imagFit->SetMarkerColor(kRed);
		imag.Add(imagFit, "L");

		TGraph* phaseFit = new TGraph;
		phaseFit->SetName((phaseName + "__fit").c_str());
		phaseFit->SetTitle((phaseName + "__fit").c_str());
		phaseFit->SetLineColor(kRed);
		phaseFit->SetLineWidth(2);
		phaseFit->SetMarkerColor(kRed);
		phase.Add(phaseFit, "L");

		// plot data
		for(size_t point = 0; point <= (fitData->nrMassBins()[idxBin]-1); ++point) {
			const size_t idxMass = point;
			const double mass = fitData->massBinCenters()[idxBin][idxMass];
			const double halfBin = 0.5 * (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1);

			realData->SetPoint(point, mass, fitData->plottingSpinDensityMatrixElementsReal()[idxBin][idxMass][idxWave][jdxWave].first);
			realData->SetPointError(point, halfBin, fitData->plottingSpinDensityMatrixElementsReal()[idxBin][idxMass][idxWave][jdxWave].second);

			imagData->SetPoint(point, mass, fitData->plottingSpinDensityMatrixElementsImag()[idxBin][idxMass][idxWave][jdxWave].first);
			imagData->SetPointError(point, halfBin, fitData->plottingSpinDensityMatrixElementsImag()[idxBin][idxMass][idxWave][jdxWave].second);

			phaseData->SetPoint(point, mass, fitData->plottingPhases()[idxBin][idxMass][idxWave][jdxWave].first);
			phaseData->SetPointError(point, halfBin, fitData->plottingPhases()[idxBin][idxMass][idxWave][jdxWave].second);

			if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
				const double minSR = fitData->sysPlottingSpinDensityMatrixElementsReal()[idxBin][idxMass][idxWave][jdxWave].first;
				const double maxSR = fitData->sysPlottingSpinDensityMatrixElementsReal()[idxBin][idxMass][idxWave][jdxWave].second;
				realSystematics->SetPoint(point, mass, (maxSR+minSR)/2.);
				realSystematics->SetPointError(point, halfBin, (maxSR-minSR)/2.);

				const double minSI = fitData->sysPlottingSpinDensityMatrixElementsImag()[idxBin][idxMass][idxWave][jdxWave].first;
				const double maxSI = fitData->sysPlottingSpinDensityMatrixElementsImag()[idxBin][idxMass][idxWave][jdxWave].second;
				imagSystematics->SetPoint(point, mass, (maxSI+minSI)/2.);
				imagSystematics->SetPointError(point, halfBin, (maxSI-minSI)/2.);

				const double minSP = fitData->sysPlottingPhases()[idxBin][idxMass][idxWave][jdxWave].first;
				const double maxSP = fitData->sysPlottingPhases()[idxBin][idxMass][idxWave][jdxWave].second;
				phaseSystematics->SetPoint(point, mass, (maxSP+minSP)/2.);
				phaseSystematics->SetPointError(point, halfBin, (maxSP-minSP)/2.);
			}
		}

		// plot fit, either over full or limited mass range
		const size_t firstPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][jdxWave].first) : 0;
		const size_t lastPoint = rangePlotting ? (extraBinning*fitData->wavePairMassBinLimits()[idxBin][idxWave][jdxWave].second) : (extraBinning*(fitData->nrMassBins()[idxBin]-1));
		for(size_t point = firstPoint; point <= lastPoint; ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			const std::complex<double> element = fitModel->spinDensityMatrix(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass);
			realFit->SetPoint(point-firstPoint, mass, element.real());
			imagFit->SetPoint(point-firstPoint, mass, element.imag());
		}

		// keep track of phase over full mass range
		TGraph phaseFitAll;
		for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			const double phase = fitModel->phase(fitParameters, cache, idxWave, jdxWave, idxBin, mass, idxMass) * TMath::RadToDeg();

			if(point != 0) {
				int bestOffs = 0;
				double bestDiff = std::numeric_limits<double>::max();

				double x;
				double prev;
				phaseFitAll.GetPoint(point-1, x, prev);
				for(int offs = -5; offs < 6; ++offs) {
					if(std::abs(phase + offs*360. - prev) < bestDiff) {
						bestDiff = std::abs(phase + offs*360. - prev);
						bestOffs = offs;
					}
				}

				phaseFitAll.SetPoint(point, mass, phase + bestOffs*360.);
			} else {
				phaseFitAll.SetPoint(point, mass, phase);
			}
		}

		// rectify phase graphs
		for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			double x;
			double valueFit;
			phaseFitAll.GetPoint(point, x, valueFit);

			if(idxMass != std::numeric_limits<size_t>::max()) {
				int bestOffs = 0;
				double bestDiff = std::numeric_limits<double>::max();

				double data;
				phaseData->GetPoint(idxMass, x, data);
				for(int offs = -5; offs < 6; ++offs) {
					if(std::abs(data + offs*360. - valueFit) < bestDiff) {
						bestDiff = std::abs(data + offs*360. - valueFit);
						bestOffs = offs;
					}
				}

				phaseData->SetPoint(idxMass, x, data + bestOffs*360.);
				if(fitInformation->getBin(idxBin).sysFileNames().size() > 0) {
					phaseSystematics->GetPoint(idxMass, x, data);
					phaseSystematics->SetPoint(idxMass, x, data + bestOffs*360.);
				}
			}

			// check that this mass bin should be taken into account for this
			// combination of waves
			if(point < firstPoint or point > lastPoint) {
				continue;
			}

			phaseFit->SetPoint(point-firstPoint, mass, valueFit);
		}

		outDirectory->cd();
		real.Write();
		imag.Write();
		phase.Write();
	}


	void
	createPlotsFsmd(const rpwa::resonanceFit::dataConstPtr& fitData,
	                const rpwa::resonanceFit::modelConstPtr& fitModel,
	                const rpwa::resonanceFit::parameters& fitParameters,
	                rpwa::resonanceFit::cache& cache,
	                TDirectory* outDirectory,
	                const size_t extraBinning,
	                const size_t idxBin)
	{
		if(rpwa::resonanceFit::debug()) {
			printDebug << "start creating plots for final-state mass-dependence." << std::endl;
		}

		TGraph graph;
		graph.SetName("finalStateMassDependence");
		graph.SetTitle("finalStateMassDependence");

		for(size_t point = 0; point <= (extraBinning*(fitData->nrMassBins()[idxBin]-1)); ++point) {
			const size_t idxMass = (point%extraBinning == 0) ? (point/extraBinning) : std::numeric_limits<size_t>::max();
			const double massStep = (fitData->massBinCenters()[idxBin][fitData->nrMassBins()[idxBin]-1] - fitData->massBinCenters()[idxBin][0]) / (fitData->nrMassBins()[idxBin] - 1) / extraBinning;
			const double mass = (idxMass != std::numeric_limits<size_t>::max()) ? fitData->massBinCenters()[idxBin][idxMass] : (fitData->massBinCenters()[idxBin][point/extraBinning] + (point%extraBinning) * massStep);

			graph.SetPoint(point, mass, std::norm(fitModel->getFsmd()->val(fitParameters, cache, idxBin, mass, idxMass)));
		}

		outDirectory->cd();
		graph.Write();
	}


}


void
rpwa::resonanceFit::createPlots(const rpwa::resonanceFit::informationConstPtr& fitInformation,
                                const rpwa::resonanceFit::dataConstPtr& fitData,
                                const rpwa::resonanceFit::modelConstPtr& fitModel,
                                const rpwa::resonanceFit::parameters& fitParameters,
                                rpwa::resonanceFit::cache& cache,
                                TFile* outFile,
                                const bool rangePlotting,
                                const size_t extraBinning)
{
	if(rpwa::resonanceFit::debug()) {
		printDebug << "start creating plots." << std::endl;
	}

	const bool sameMassBinning = fitData->hasSameMassBinning();
	for(size_t idxBin = 0; idxBin < fitInformation->nrBins(); ++idxBin) {
		TDirectory* outDirectory = NULL;
		if(fitInformation->nrBins() == 1) {
			outDirectory = outFile;
		} else {
			std::ostringstream name;
			name << "bin" << idxBin;
			outDirectory = outFile->mkdir(name.str().c_str());
		}

		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			createPlotsWave(fitInformation,
			                fitData,
			                fitModel,
			                fitParameters,
			                cache,
			                outDirectory,
			                rangePlotting,
			                extraBinning,
			                idxWave,
			                idxBin);
		}

		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			for(size_t jdxWave = idxWave+1; jdxWave < fitInformation->nrWaves(); ++jdxWave) {
				createPlotsWavePair(fitInformation,
				                    fitData,
				                    fitModel,
				                    fitParameters,
				                    cache,
				                    outDirectory,
				                    rangePlotting,
				                    extraBinning,
				                    idxWave,
				                    jdxWave,
				                    idxBin);
			}
		}

		if(fitModel->getFsmd() and (fitModel->getFsmd()->getNrBins() != 1 or not sameMassBinning)) {
			createPlotsFsmd(fitData,
			                fitModel,
			                fitParameters,
			                cache,
			                outDirectory,
			                extraBinning,
			                idxBin);
		}
	}

	if(fitInformation->nrBins() != 1 and sameMassBinning and fitModel->isMappingEqualInAllBins()) {
		for(size_t idxWave = 0; idxWave < fitInformation->nrWaves(); ++idxWave) {
			createPlotsWaveSum(fitInformation,
			                   fitData,
			                   fitModel,
			                   fitParameters,
			                   cache,
			                   outFile,
			                   rangePlotting,
			                   extraBinning,
			                   idxWave);
		}
	}

	if(fitModel->getFsmd() and (fitModel->getFsmd()->getNrBins() == 1 and sameMassBinning)) {
		createPlotsFsmd(fitData,
		                fitModel,
		                fitParameters,
		                cache,
		                outFile,
		                extraBinning,
		                0);
	}

	if(rpwa::resonanceFit::debug()) {
		printDebug << "finished creating plots." << std::endl;
	}
}
