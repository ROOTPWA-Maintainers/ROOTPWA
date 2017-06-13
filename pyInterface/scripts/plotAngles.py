#!/usr/bin/env python

import argparse
import itertools
import os
import sys

import numpy

import pyRootPwa
import pyRootPwa.utils
ROOT = pyRootPwa.ROOT


def combineSyms(boseSyms, isoSyms):
	newSyms = []
	for isoSym in isoSyms:
		for boseSym in boseSyms:
			if len(isoSym["fsPartPermMap"]) != len(boseSym["fsPartPermMap"]):
				pyRootPwa.utils.printErr("final state particle mismatch when combining symmetrizations. Aborting...")
				sys.exit(1)
			newPerm = { "fsPartPermMap": [], "factor": boseSym["factor"]*isoSym["factor"] }
			for i in xrange(len(isoSym["fsPartPermMap"])):
				newPerm["fsPartPermMap"].append(boseSym["fsPartPermMap"][isoSym["fsPartPermMap"][i]])
			newSyms.append(newPerm)
	return newSyms


def getPermutations(topology):
	boseSyms = topology.getBoseSymmetrization()
	isoSyms = topology.getIsospinSymmetrization()
	if isoSyms:
		boseSyms = combineSyms(boseSyms, isoSyms)
	keys = [tuple(x["fsPartPermMap"]) for x in boseSyms]
	permutations = {}
	for key in keys:
		permutations[key] = [[None, None] for x in range(topology.nmbDecayVertices())]
	vertex_i = 0
	for vertex in topology.decayVertices():
		daughter1 = vertex.daughter1()
		daughter2 = vertex.daughter2()
		fsParts1 = []
		fsParts2 = []
		if topology.isFsParticle(daughter1):
			fsParts1.append(topology.fsParticlesIndex(daughter1))
		else:
			fsParts1 = topology.getFsPartIndicesConnectedToVertex(topology.toVertex(daughter1))
		if topology.isFsParticle(daughter2):
			fsParts2.append(topology.fsParticlesIndex(daughter2))
		else:
			fsParts2 = topology.getFsPartIndicesConnectedToVertex(topology.toVertex(daughter2))
		fsPartsParent = topology.getFsPartIndicesConnectedToVertex(vertex)
		keys = permutations.keys()
		for i in range(len(keys)):
			permutation = keys[i]
			# this is the part for the masses
			massGroup = sorted([permutation[x] for x in fsPartsParent])
			# this is the part for the angles
			group1 = sorted([permutation[x] for x in fsParts1])
			group2 = sorted([permutation[x] for x in fsParts2])
			# and now check all other permutations
			massAlreadyThere = False
			anglesAlreadyThere = False
			for j in range(i):
				testPermutation = keys[j]
				# this is the part for the masses
				testMassGroup = sorted([testPermutation[x] for x in fsPartsParent])
				if massGroup == testMassGroup:
					massAlreadyThere = True
				# this is the part for the angles
				testGroup1 = sorted([testPermutation[x] for x in fsParts1])
				testGroup2 = sorted([testPermutation[x] for x in fsParts2])
				if group1 == testGroup1 and group2 == testGroup2:
					anglesAlreadyThere = True
				if massAlreadyThere and anglesAlreadyThere:
					break
			permutations[permutation][vertex_i] = [not massAlreadyThere, not anglesAlreadyThere]
		vertex_i += 1
	return permutations


def getPermutationsForIsobarCombinations(topology, isobarCombinations):
	boseSyms = topology.getBoseSymmetrization()
	isoSyms = topology.getIsospinSymmetrization()
	if isoSyms:
		boseSyms = combineSyms(boseSyms, isoSyms)
	keys = [tuple(x["fsPartPermMap"]) for x in boseSyms]
	permutations = {}
	for key in keys:
		permutations[key] = [None for x in range(len(isobarCombinations))]
	for isobarCombination_i in range(len(isobarCombinations)):
		isobar1 = topology.isobarDecayVertices()[isobarCombinations[isobarCombination_i][0]]
		isobar2 = topology.isobarDecayVertices()[isobarCombinations[isobarCombination_i][1]]
		fsParts1 = topology.getFsPartIndicesConnectedToVertex(isobar1)
		fsParts2 = topology.getFsPartIndicesConnectedToVertex(isobar2)
		for i in range(len(keys)):
			permutation = keys[i]
			massGroup1 = sorted([permutation[x] for x in fsParts1])
			massGroup2 = sorted([permutation[x] for x in fsParts2])
			alreadyThere = False
			for j in range(i):
				testPermutation = keys[j]
				testMassGroup1 = sorted([testPermutation[x] for x in fsParts1])
				testMassGroup2 = sorted([testPermutation[x] for x in fsParts2])
				if massGroup1 == testMassGroup1 and massGroup2 == testMassGroup2:
					alreadyThere = True
					break
			permutations[permutation][isobarCombination_i] = not alreadyThere
	return permutations


if __name__ == "__main__":

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="No help text available."
	                                )
	parser.add_argument("inputFile", type=str, metavar="inputFile", help="inputFile")
	parser.add_argument("outputFile", metavar="output-file", help="path to output file")
	parser.add_argument("templateFile", metavar="template-file", help="path to template file")
	parser.add_argument("-b", action="append", metavar="massBin(s)", default=[], dest="massBins", help="mass bins to be calculated (default: all)")
	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-g", "--n-bins", type=int, metavar="n-bins", default=100, dest="nHistogramBins", help="number of bins for the histograms (default: 100)")
	parser.add_argument("-g2d", "--n-2D-bins", type=int, metavar="n-bins", default=50, dest="nHistogram2DBins", help="number of bins for the 2D histograms (default: 50)")
	parser.add_argument("-i", "--min-mass", type=float, metavar="min-mass", default=0.0, dest="massMin", help="minimum mass in GeV/c^{2} for the histograms (default: 0)")
	parser.add_argument("-a", "--max-mass", type=float, metavar="max-mass", default=10.0, dest="massMax", help="maximum mass in GeV/c^{2} for the histograms (default: 10)")
	parser.add_argument("-t", "--type", type=str, metavar="type", default="data", dest="type", help='type of input, can be "data", "gen" or "acc" (default: "data")')
	parser.add_argument("-m", "--mass-hist", type=str, metavar="xMassHist", dest="xMassHistogramPath",
	                    help='histogram of the X mass, specified as pathToRootFile.root/pathOfHistInRootFile')
	parser.add_argument("--disable-bose-symmetrization", action="store_true", dest="disableBoseSymmetrization", help="do not consider Bose-symmetric permutations")
	parser.add_argument("--weight-file", type=str, metavar="weightFile", dest="weightFile", help="weightFile")
	arguments = parser.parse_args()
	if not arguments.massBins:
		arguments.massBins.append("all")

	if arguments.type not in ["data", "acc", "gen"]:
		pyRootPwa.utils.printErr("invalid type '" + arguments.type + "' found, must be either 'data', 'gen' or 'acc'. Aborting...")
		sys.exit(5)
	inputTypeIndex = {"data": 0, "gen": 1, "acc": 2}

	xMassHistogram = None
	massSliceHistogram = None
	if arguments.xMassHistogramPath != None:
		fullPath = arguments.xMassHistogramPath
		if fullPath.find('.root/') < 0:
			pyRootPwa.utils.printWarn("Path to X mass histogram '" + fullPath +
			                          "' not valid. Not using X mass histogram.")
		else:
			(rootFilePath, histogramPath) = fullPath.rsplit('.root/', 1)
			rootFilePath += '.root'
			if not os.path.isfile(rootFilePath):
				pyRootPwa.utils.printWarn("X mass histogram root file '" + rootFilePath +
				                          "' not found. Not using X mass histogram.")
			else:
				rootFile = pyRootPwa.ROOT.TFile.Open(rootFilePath, "READ")
				xMassHistogram = rootFile.Get(histogramPath)
				if not xMassHistogram or xMassHistogram.ClassName() != "TH1D":
					pyRootPwa.utils.printWarn("X mass histogram '" + histogramPath +
					                          "' not found in root file '" + rootFilePath +
					                          "'. Not using X mass histogram.")
					xMassHistogram = None
	if xMassHistogram is not None:
		massSliceHistogram = xMassHistogram.Clone("mass_slice")
		massSliceHistogram.Reset()
		pyRootPwa.utils.printInfo("X mass histogram found. Switching root to batch mode.")
		pyRootPwa.ROOT.gROOT.SetBatch(True)

	pyRootPwa.config = pyRootPwa.rootPwaConfig()
	if not pyRootPwa.config.initialize(arguments.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + arguments.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(pyRootPwa.config.pdgFileName)

	waveDescs = pyRootPwa.core.waveDescription.parseKeyFile(arguments.templateFile)

	if len(waveDescs) != 1:
		pyRootPwa.utils.printErr("could not read wave description. Aborting...")
		sys.exit(5)

	(result, topology) = waveDescs[0].constructDecayTopology(True)

	if not result:
		pyRootPwa.utils.printErr("could not construct topology. Aborting...")
		sys.exit(5)

	isobarCombinations = list(itertools.combinations(range(topology.nmbDecayVertices()), 2))
	for i in range(len(isobarCombinations)):
		vertex1 = topology.isobarDecayVertices()[isobarCombinations[i][0]]
		vertex2 = topology.isobarDecayVertices()[isobarCombinations[i][1]]
		fsParts1 = topology.getFsPartIndicesConnectedToVertex(vertex1)
		fsParts2 = topology.getFsPartIndicesConnectedToVertex(vertex2)
		if len(fsParts2) > len(fsParts1):
			isobarCombinations[i] = (isobarCombinations[i][1], isobarCombinations[i][0])
	isobarPermutations = getPermutationsForIsobarCombinations(topology, isobarCombinations)

	if arguments.disableBoseSymmetrization:
		permutations = {}
		permutation = [x for x in range(topology.nmbFsParticles())]
		permutations[tuple(permutation)] = [(True, True) for i in range(topology.nmbDecayVertices())]
	else:
		permutations = getPermutations(topology)

	outputFile = pyRootPwa.ROOT.TFile.Open(arguments.outputFile, "RECREATE")

	inputFileRanges = {}
	hists = {}
	isobarCombinationHistograms = {}
	massBinWidth = '%.2f' % (1000.0*((arguments.massMax - arguments.massMin)/arguments.nHistogramBins))
	cosBinWidth = '%.2f' % (2000.0 / arguments.nHistogramBins)
	thetaBinWidth = '%.2f' % (pyRootPwa.ROOT.TMath.TwoPi() / arguments.nHistogramBins)
	for binRange in arguments.massBins:
		rangeName = ""
		eventFile = ROOT.TFile.Open(arguments.inputFile)
		evtMeta = pyRootPwa.core.eventMetadata.readEventFile(eventFile)
		if 'mass' not in evtMeta.binningMap():
			pyRootPwa.utils.printErr("'mass' is not in the binning map. Aborting...")
			sys.exit(1)
		rangeName = str(int(round((evtMeta.binningMap()['mass'][0] + evtMeta.binningMap()['mass'][1])/2.)))
		inputFileRanges[rangeName] = [ arguments.inputFile ]
		eventFile.Close()
		if not inputFileRanges[rangeName]:
			pyRootPwa.utils.printErr("no input files found for mass bins " + str(rangeName) + ".")
		outputFile.mkdir(rangeName)
		outputFile.cd(rangeName)
		hists[rangeName] = []
		for i in range(topology.nmbDecayVertices()):
			parent = topology.isobarDecayVertices()[i].parent()
			daughter1 = topology.isobarDecayVertices()[i].daughter1()
			daughter2 = topology.isobarDecayVertices()[i].daughter2()
			label = parent.name + " -> [" + daughter1.name + " " + daughter2.name + "]"
			name = parent.name + "_" + daughter1.name + "_" + daughter2.name
			nEntriesMass = 0
			nEntriesAngles = 0
			for permutationKey in permutations:
				tuples = permutations[permutationKey][i]
				if tuples[0]:
					nEntriesMass += 1
				if tuples[1]:
					nEntriesAngles += 1
			massYAxisTitle = ""
			if nEntriesMass == 1:
				massYAxisTitle = "events / " + massBinWidth + " MeV"
			else:
				massYAxisTitle = str(nEntriesMass) + " entries per event / " + massBinWidth + " MeV"
			cosYAxisTitle = ""
			thetaYAxisTitle = ""
			if nEntriesAngles == 1:
				cosYAxisTitle = "events / " + cosBinWidth + " mrad"
				thetaYAxisTitle = "events"
			else:
				cosYAxisTitle = str(nEntriesAngles) + " entries per event / " + cosBinWidth + " mrad"
				thetaYAxisTitle = str(nEntriesAngles) + " entries per event"

			massTitle = "m(" + parent.name + ");m(" + parent.name + ") [GeV/c^{2}];" + massYAxisTitle
			phiTitle = "#phi(" + label + ");#phi(" + label + ") [rad];" + cosYAxisTitle
			thetaTitle = "cos #theta(" + label + ");cos #theta(" + label + ");" + thetaYAxisTitle
			phiVsMassTitle = "#phi(" + label + ") vs. m(" + parent.name + ");m(" + parent.name + ") [GeV/c^{2}];#phi(" + label + ") [rad];" + thetaYAxisTitle
			thetaVsMassTitle = "cos #theta(" + label + ") vs. m(" + parent.name + ");m(" + parent.name + ") [GeV/c^{2}];cos #theta(" + label + ");" + thetaYAxisTitle
			phiVsThetaTitle = "#phi vs. cos(#theta) (" + label + ");cos #theta(" + label + ");#phi(" + label + ") [rad];" + thetaYAxisTitle

			hists[rangeName].append([pyRootPwa.ROOT.TH1D("m_" + parent.name, massTitle, arguments.nHistogramBins, arguments.massMin, arguments.massMax),
									 pyRootPwa.ROOT.TH1D(name + "_phi", phiTitle, arguments.nHistogramBins, -pyRootPwa.ROOT.TMath.Pi(), pyRootPwa.ROOT.TMath.Pi()),
			                         pyRootPwa.ROOT.TH1D(name + "_cosTheta", thetaTitle, arguments.nHistogramBins, -1, 1),
			                         pyRootPwa.ROOT.TH2D(name + "_phi_vs_m_" + parent.name, phiVsMassTitle, arguments.nHistogram2DBins, arguments.massMin,
			                                             arguments.massMax, arguments.nHistogram2DBins, -pyRootPwa.ROOT.TMath.Pi(), pyRootPwa.ROOT.TMath.Pi()),
			                         pyRootPwa.ROOT.TH2D(name + "_cosTheta_vs_m_" + parent.name, thetaVsMassTitle, arguments.nHistogram2DBins, arguments.massMin,
			                                             arguments.massMax, arguments.nHistogram2DBins, -1, 1),
			                         pyRootPwa.ROOT.TH2D(name + "_phi_vs_cosTheta", phiVsThetaTitle, arguments.nHistogram2DBins, -1, 1, arguments.nHistogram2DBins,
			                                             -pyRootPwa.ROOT.TMath.Pi(), pyRootPwa.ROOT.TMath.Pi())])
			for hist in hists[rangeName][-1]:
				hist.SetMinimum(0)
		isobarCombinationHistograms[rangeName] = {}
		for isobarCombination_i in range(len(isobarCombinations)):
			isobarCombination = isobarCombinations[isobarCombination_i]
			index1 = isobarCombination[0]
			index2 = isobarCombination[1]
			isobar1Name = topology.isobarDecayVertices()[index1].parent().name
			isobar2Name = topology.isobarDecayVertices()[index2].parent().name
			nEntries = 0
			for permutationKey in isobarPermutations:
				if isobarPermutations[permutationKey][isobarCombination_i]:
					nEntries += 1
			if nEntries == 1:
				zAxisTitle = "events"
			else:
				zAxisTitle = str(nEntries) + " entries per event"
			histName = "m(" + isobar2Name + ") vs. m(" + isobar1Name + ");m(" + isobar1Name + ") [GeV/c^{2}];m(" + isobar2Name + ") [GeV/c^{2}];" + zAxisTitle
			isobarCombinationHistograms[rangeName][isobarCombination] = pyRootPwa.ROOT.TH2D("m_" + isobar2Name + "_vs_m_" + isobar1Name, histName, arguments.nHistogram2DBins,
			                                                                                arguments.massMin, arguments.massMax, arguments.nHistogram2DBins,
			                                                                                arguments.massMin, arguments.massMax)

	assert inputFileRanges.keys() == hists.keys()

	for rangeName in inputFileRanges:

		pyRootPwa.utils.printInfo("Processing bin range " + rangeName)
		outputFile.cd(rangeName)

		for dataFileName in inputFileRanges[rangeName]:

			# Do all the initialization
			pyRootPwa.utils.printInfo("Opening input file " + dataFileName)
			dataFile = pyRootPwa.ROOT.TFile.Open(dataFileName)
			evtMeta = pyRootPwa.core.eventMetadata.readEventFile(dataFile)
			dataTree = evtMeta.eventTree()

			if arguments.weightFile is not None:
				dataTree.AddFriend(pyRootPwa.config.weightTreeName, arguments.weightFile)

			topology.initKinematicsData(evtMeta.productionKinematicsParticleNames(), evtMeta.decayKinematicsParticleNames())
			nEvents = dataTree.GetEntries()
			prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			dataTree.SetBranchAddress(pyRootPwa.core.eventMetadata.productionKinematicsMomentaBranchName, prodKinMomenta)
			dataTree.SetBranchAddress(pyRootPwa.core.eventMetadata.decayKinematicsMomentaBranchName, decayKinMomenta)

			# Handle weighted MC
			weight = numpy.array(1, dtype = float)
			if arguments.weightFile is not None:
				dataTree.SetBranchAddress("weight", weight)
			else:
				weight = 1.

			progressbar = pyRootPwa.utils.progressBar(0, nEvents-1, sys.stdout)
			progressbar.start()

			# Loop over Events
			for i in range(nEvents):
				dataTree.GetEntry(i)

				# Read input data
				topology.readKinematicsData(prodKinMomenta, decayKinMomenta)

				for permutationKey in permutations:

					permutation = list(permutationKey)

					topology.revertMomenta(permutation)

					# Transform all particles
					topology.calcIsobarLzVec()
					beamLv = topology.productionVertex().referenceLzVec()
					XLv = topology.XParticle().lzVec
					gjTrans = pyRootPwa.core.isobarAmplitude.gjTransform(beamLv, XLv)
					for vertex in topology.isobarDecayVertices():
						vertex.transformOutParticles(gjTrans)
					for vertex in topology.isobarDecayVertices()[1:]:
						hfTrans = pyRootPwa.core.isobarHelicityAmplitude.hfTransform(vertex.parent().lzVec)
						for subVertex in topology.subDecay(vertex).decayVertices():
							subVertex.transformOutParticles(hfTrans)

					# Fill the histograms
					hist_i = 0

					masses = []
					for vertex in topology.isobarDecayVertices():
						parent = vertex.parent()
						mass = parent.lzVec.M()
						masses.append(mass)
						if permutations[permutationKey][hist_i][0]:
							hists[rangeName][hist_i][0].Fill(mass, weight)
							if hist_i == 0 and xMassHistogram is not None:
								massSliceHistogram.Fill(mass, weight)
						if permutations[permutationKey][hist_i][1]:
							daughter = vertex.daughter1()
							phi = daughter.lzVec.Phi()
							theta = daughter.lzVec.Theta()
							cosTheta = pyRootPwa.ROOT.TMath.Cos(theta)
							hists[rangeName][hist_i][1].Fill(phi, weight)
							hists[rangeName][hist_i][2].Fill(cosTheta, weight)
							hists[rangeName][hist_i][3].Fill(mass, phi, weight)
							hists[rangeName][hist_i][4].Fill(mass, cosTheta, weight)
							hists[rangeName][hist_i][5].Fill(cosTheta, phi, weight)
						hist_i += 1

					for isobarCombination in isobarCombinations:
						if not isobarPermutations[permutationKey]:
							continue
						index1 = isobarCombination[0]
						index2 = isobarCombination[1]
						isobarCombinationHistograms[rangeName][isobarCombination].Fill(masses[index1], masses[index2], weight)

				progressbar.update(i)

			dataFile.Close()

		if xMassHistogram is not None:
			xMassCan = pyRootPwa.ROOT.TCanvas("X_mass", "X mass", 1000, 700)
			xMassCan.cd()
			xMassHistogram.Draw()
			massSliceHistogram.Draw("same")
			massSliceHistogram.SetFillColor(pyRootPwa.ROOT.kOrange)
			outputFile.cd(rangeName)
			xMassCan.Write()
			xMassCan.Close()
			allHistograms = []
			for i in range(len(hists[rangeName])):
				for j in range(len(hists[rangeName][i])):
					allHistograms.append(hists[rangeName][i][j])
			for isobarCombination in isobarCombinations:
				allHistograms.append(isobarCombinationHistograms[rangeName][isobarCombination])
			for hist in allHistograms:
				canvas = pyRootPwa.ROOT.TCanvas(hist.GetName() + "_canvas", hist.GetTitle(), 1000, 1400)
				canvas.Divide(1, 2)
				canvas.cd(1)
				xMassHistogram.Draw()
				massSliceHistogram.Draw("same")
				massSliceHistogram.SetFillColor(pyRootPwa.ROOT.kOrange)
				canvas.cd(2)
				hist.Draw()
				if hist.ClassName() == "TH2D":
					hist.SetDrawOption("colz")
				else:
					hist.SetFillColor(pyRootPwa.ROOT.kOrange)
				canvas.Write()
				canvas.Close()
			massSliceHistogram.Reset()

	outputFile.Write()
	outputFile.Close()

	pyRootPwa.utils.printPrintingSummary(pyRootPwa.utils.printingCounter)
