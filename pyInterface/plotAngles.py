#!/usr/bin/python2.7

import argparse
import sys
import time

import pyRootPwa
import pyRootPwa.utils

if __name__ == "__main__":

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="No help text available."
	                                )
	parser.add_argument("dataFile", metavar="data-file", help="path to data file")
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
		topology.readKinematicsData(prodKinMomenta, decayKinMomenta)
		progressbar.update(i)

	pyRootPwa.utils.printPrintingSummary(printingCounter)
