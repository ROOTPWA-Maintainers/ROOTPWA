#!/usr/bin/env python

import argparse
import os.path

import pyRootPwa
ROOT = pyRootPwa.ROOT

if __name__ == "__main__":

	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	parser = argparse.ArgumentParser(description="convert event file")
	parser.add_argument("inputFile", type=str, metavar="inputFile", help="input file in ROOTPWA format without meta data")
	parser.add_argument("outputFile", type=str, metavar="outputFile", help="input file in ROOTPWA format with meta data")
	parser.add_argument("-l", type=str, metavar="string", dest="userString", help="label which is saved to the metadata (default: output file name)")
	parser.add_argument("-v", action="store_true", dest="debug", help="verbose; print debug output (default: false)")
	args = parser.parse_args()

	prodKinPartNamesObjName = "prodKinParticles"
	prodKinMomentaLeafName = "prodKinMomenta"
	decayKinPartNamesObjName = "decayKinParticles"
	decayKinMomentaLeafName = "decayKinMomenta"
	inTreeName = "rootPwaEvtTree"

	userString = "" if not args.userString else args.userString

	inputFileName = os.path.abspath(args.inputFile)
	inputFile = ROOT.TFile(inputFileName, "READ")

	outputFileName = os.path.abspath(args.outputFile)
	outputFile = ROOT.TFile(outputFileName, "NEW")

	initialStateParticleNames = []
	for name in inputFile.Get(prodKinPartNamesObjName):
		initialStateParticleNames.append(str(name))
	finalStateParticleNames = []
	for name in inputFile.Get(decayKinPartNamesObjName):
		finalStateParticleNames.append(str(name))

	fileWriter = pyRootPwa.core.eventFileWriter()
	fileWriter.initialize(outputFile,
	                      userString,
	                      initialStateParticleNames,
	                      finalStateParticleNames,
	                      {},
	                      [],
	                      inTreeName,
	                      prodKinMomentaLeafName,
	                      decayKinMomentaLeafName)

	inputTree = inputFile.Get(inTreeName)
	for event in inputTree:
		initialStateMomenta = []
		for particle in event.__getattr__(prodKinMomentaLeafName):
			initialStateMomenta.append(ROOT.TVector3(particle))
		finalStateMomenta = []
		for particle in event.__getattr__(decayKinMomentaLeafName):
			finalStateMomenta.append(ROOT.TVector3(particle))
		fileWriter.addEvent(initialStateMomenta, finalStateMomenta)

	fileWriter.finalize()
