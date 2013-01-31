#!/usr/bin/env python

import argparse
import glob
import itertools
import sys

import pyRootPwa
import pyRootPwa.utils


def getPermutations(topology):
	boseSyms = topology.getBoseSymmetrization()
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
	parser.add_argument("outputFile", metavar="output-file", help="path to output file")
	parser.add_argument("templateFile", metavar="template-file", help="path to template file")
	parser.add_argument("-b", action="append", metavar="massBin(s)", default=[], dest="massBins", help="mass bins to be calculated (default: all)")
	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-g", "--n-bins", type=int, metavar="n-bins", default=100, dest="nHistogramBins", help="number of bins for the histograms (default: 100)")
	parser.add_argument("-g2d", "--n-2D-bins", type=int, metavar="n-bins", default=50, dest="nHistogram2DBins", help="number of bins for the 2D histograms (default: 50)")
	parser.add_argument("-i", "--min-mass", type=float, metavar="min-mass", default=0.0, dest="massMin", help="minimum mass in GeV/c^{2} for the histograms (default: 0)")
	parser.add_argument("-a", "--max-mass", type=float, metavar="max-mass", default=10.0, dest="massMax", help="maximum mass in GeV/c^{2} for the histograms (default: 10)")
	parser.add_argument("-t", "--type", type=str, metavar="type", default="data", dest="type", help='type of input, can be "data", "gen" or "acc" (default: "data")')
	parser.add_argument("--disable-bose-symmetrization", action="store_true", dest="disableBoseSymmetrization", help="do not consider Bose-symmetric permutations")
	arguments = parser.parse_args()
	if len(arguments.massBins) == 0:
		arguments.massBins.append("all")

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	printingCounter = 5 * [0]
	pyRootPwa.utils.printErr = pyRootPwa.utils.printErrClass(printingCounter)
	pyRootPwa.utils.printWarn = pyRootPwa.utils.printWarnClass(printingCounter)
	pyRootPwa.utils.printSucc = pyRootPwa.utils.printSuccClass(printingCounter)
	pyRootPwa.utils.printInfo = pyRootPwa.utils.printInfoClass(printingCounter)
	pyRootPwa.utils.printDebug = pyRootPwa.utils.printDebugClass(printingCounter)

	if arguments.type not in ["data", "acc", "gen"]:
		pyRootPwa.utils.printErr("invalid type '" + arguments.type + "' found, must be either 'data', 'gen' or 'acc'. Aborting...")
		sys.exit(5)
	inputTypeIndex = {"data": 0, "gen": 1, "acc": 2}

	pyRootPwa.config = pyRootPwa.rootPwaConfig(arguments.configFileName)
	pyRootPwa.core.particleDataTable.readFile(pyRootPwa.config.pdgFileName)

	waveDesc = pyRootPwa.core.waveDescription()
	waveDesc.parseKeyFile(arguments.templateFile)
	(result, topology) = waveDesc.constructDecayTopology(True)

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
	massBinWidth = '%.2f' % (1000.0*((arguments.massMax - arguments.massMin)/arguments.nHistogramBins))
	cosBinWidth = '%.2f' % (2000.0 / arguments.nHistogramBins)
	thetaBinWidth = '%.2f' % (pyRootPwa.ROOT.TMath.TwoPi() / arguments.nHistogramBins)
	for binRange in arguments.massBins:
		allMassBins = sorted(glob.glob(pyRootPwa.config.dataDirectory + '/' + pyRootPwa.config.massBinDirectoryNamePattern))
		massBins = pyRootPwa.utils.parseMassBinArgs(allMassBins, binRange)
		rangeName = massBins[0].rsplit('/', 1)[-1] + '_' + massBins[-1].rsplit('/', 1)[-1]
		inputFileRanges[rangeName] = pyRootPwa.utils.getListOfInputFiles(massBins)[inputTypeIndex[arguments.type]]
		if not inputFileRanges[rangeName]:
			pyRootPwa.utils.printErr("No input files found for mass bins " + str(massBins) + ".")
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
			for permutationKey in permutations.keys():
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
			                         pyRootPwa.ROOT.TH2D(name + "_phi_vs_m_" + parent.name, phiVsMassTitle, arguments.nHistogram2DBins, arguments.massMin, arguments.massMax, arguments.nHistogram2DBins, -pyRootPwa.ROOT.TMath.Pi(), pyRootPwa.ROOT.TMath.Pi()),
			                         pyRootPwa.ROOT.TH2D(name + "_cosTheta_vs_m_" + parent.name, thetaVsMassTitle, arguments.nHistogram2DBins, arguments.massMin, arguments.massMax, arguments.nHistogram2DBins, -1, 1),
			                         pyRootPwa.ROOT.TH2D(name + "_phi_vs_cosTheta", phiVsThetaTitle, arguments.nHistogram2DBins, -1, 1, arguments.nHistogram2DBins, -pyRootPwa.ROOT.TMath.Pi(), pyRootPwa.ROOT.TMath.Pi())])
			for hist in hists[rangeName][-1]:
				hist.SetMinimum(0)
		isobarCombinationHistograms = {}
		for isobarCombination_i in range(len(isobarCombinations)):
			isobarCombination = isobarCombinations[isobarCombination_i]
			index1 = isobarCombination[0]
			index2 = isobarCombination[1]
			isobar1Name = topology.isobarDecayVertices()[index1].parent().name
			isobar2Name = topology.isobarDecayVertices()[index2].parent().name
			nEntries = 0
			for permutationKey in isobarPermutations.keys():
				if isobarPermutations[permutationKey][isobarCombination_i]:
					nEntries += 1
			if nEntries == 1:
				zAxisTitle = "events"
			else:
				zAxisTitle = str(nEntries) + " entries per event"
			histName = "m(" + isobar2Name + ") vs. m(" + isobar1Name + ");m(" + isobar1Name + ") [GeV/c^{2}];m(" + isobar2Name + ") [GeV/c^{2}];" + zAxisTitle
			isobarCombinationHistograms[isobarCombination] = pyRootPwa.ROOT.TH2D("m_" + isobar2Name + "_vs_m_" + isobar1Name, histName, arguments.nHistogram2DBins, arguments.massMin, arguments.massMax, arguments.nHistogram2DBins, arguments.massMin, arguments.massMax)

	assert(inputFileRanges.keys() == hists.keys())

	for rangeName in inputFileRanges.keys():

		pyRootPwa.utils.printInfo("Processing bin range " + rangeName)
		outputFile.cd(rangeName)

		for dataFileName in inputFileRanges[rangeName]:

			# Do all the initialization
			pyRootPwa.utils.printInfo("Opening input file " + dataFileName)
			dataFile = pyRootPwa.ROOT.TFile.Open(dataFileName)
			dataTree = dataFile.Get(pyRootPwa.config.inTreeName)

			prodKinParticles = dataFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
			decayKinParticles = dataFile.Get(pyRootPwa.config.decayKinPartNamesObjName)
			topology.initKinematicsData(prodKinParticles, decayKinParticles)
			nEvents = dataTree.GetEntries()
			prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
			dataTree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, prodKinMomenta)
			dataTree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, decayKinMomenta)

			progressbar = pyRootPwa.utils.progressBar(0, nEvents-1, sys.stdout)
			progressbar.start()

			# Loop over Events
			for i in range(nEvents):
				dataTree.GetEntry(i)

				# Read input data
				topology.readKinematicsData(prodKinMomenta, decayKinMomenta)

				for permutationKey in permutations.keys():

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
						if(permutations[permutationKey][hist_i][0]):
							hists[rangeName][hist_i][0].Fill(mass)
						if(permutations[permutationKey][hist_i][1]):
							daughter = vertex.daughter1()
							phi = daughter.lzVec.Phi()
							theta = daughter.lzVec.Theta()
							cosTheta = pyRootPwa.ROOT.TMath.Cos(theta)
							hists[rangeName][hist_i][1].Fill(phi)
							hists[rangeName][hist_i][2].Fill(cosTheta)
							hists[rangeName][hist_i][3].Fill(mass, phi)
							hists[rangeName][hist_i][4].Fill(mass, cosTheta)
							hists[rangeName][hist_i][5].Fill(cosTheta, phi)
						hist_i += 1

					for isobarCombination_i in range(len(isobarCombinations)):
						isobarCombination = isobarCombinations[isobarCombination_i]
						if not isobarPermutations[permutationKey]:
							continue
						index1 = isobarCombination[0]
						index2 = isobarCombination[1]
						isobarCombinationHistograms[isobarCombination].Fill(masses[index1], masses[index2])

				progressbar.update(i)

			dataFile.Close()

	outputFile.Write()
	outputFile.Close()

	pyRootPwa.utils.printPrintingSummary(printingCounter)
