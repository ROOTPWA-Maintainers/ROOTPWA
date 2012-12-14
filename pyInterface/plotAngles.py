#!/usr/bin/python2.7

import argparse
import sys
import time

import pyRootPwa
import pyRootPwa.utils

def randomizeTClonesArray(array):
	entries = array.GetEntries()
	if entries == 0:
		return array
	retval = pyRootPwa.ROOT.TClonesArray("TVector3")
	map = []
	for i in range(entries):
		while True:
			num = int(round(pyRootPwa.ROOT.gRandom.Uniform(-0.5, entries-0.5)))
			if num not in map:
				map.append(num)
				break
	for i in range(entries):
		retval[i] = array[map[i]]
	return retval

if __name__ == "__main__":

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="No help text available."
	                                )
	parser.add_argument("dataFile", metavar="data-file", help="path to data file")
	parser.add_argument("outputFile", metavar="output-file", help="path to output file")
	parser.add_argument("templateFile", metavar="template-file", help="path to template file")
	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	arguments = parser.parse_args()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	printingCounter = 5 * [0]
	pyRootPwa.utils.printErr = pyRootPwa.utils.printErrClass(printingCounter)
	pyRootPwa.utils.printWarn = pyRootPwa.utils.printWarnClass(printingCounter)
	pyRootPwa.utils.printSucc = pyRootPwa.utils.printSuccClass(printingCounter)
	pyRootPwa.utils.printInfo = pyRootPwa.utils.printInfoClass(printingCounter)
	pyRootPwa.utils.printDebug = pyRootPwa.utils.printDebugClass(printingCounter)

	pyRootPwa.config = pyRootPwa.rootPwaConfig(arguments.configFileName)
	pyRootPwa.core.particleDataTable.readFile(pyRootPwa.config.pdgFileName)

	waveDesc = pyRootPwa.core.waveDescription()
	waveDesc.parseKeyFile(arguments.templateFile)
	(result, topology) = waveDesc.constructDecayTopology(True)

	if not result:
		pyRootPwa.utils.printErr("could not construct topology. Aborting...")
		sys.exit(5)

	outputFile = pyRootPwa.ROOT.TFile.Open(arguments.outputFile, "RECREATE")

	hists = []
	for i in range(topology.nmbDecayVertices()):
		hists.append([pyRootPwa.ROOT.TH1D("phi" + str(i), "phi" + str(i), 100, -3.142, 3.142),
		              pyRootPwa.ROOT.TH1D("theta" + str(i), "theta" + str(i), 100, -1, 1),
		              pyRootPwa.ROOT.TH1D("m" + str(i), "m" + str(i), 100, 0, 5)])

	dataFile = pyRootPwa.ROOT.TFile.Open(arguments.dataFile)
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

	for i in range(nEvents):
		dataTree.GetEntry(i)

		# Transform all particles
		topology.readKinematicsData(prodKinMomenta, randomizeTClonesArray(decayKinMomenta))
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

		hist_i = 0
		for vertex in topology.isobarDecayVertices():
			daughter = vertex.daughter1()
			parent = vertex.parent()
			mass = parent.lzVec.M()
			phi = daughter.lzVec.Phi()
			theta = daughter.lzVec.Theta()
			hists[hist_i][0].Fill(phi)
			hists[hist_i][1].Fill(pyRootPwa.ROOT.TMath.Cos(theta))
			hists[hist_i][2].Fill(mass)
			hist_i += 1

		progressbar.update(i)

	outputFile.Write()
	outputFile.Close()

	pyRootPwa.utils.printPrintingSummary(printingCounter)
