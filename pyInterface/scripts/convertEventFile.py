#!/usr/bin/env python

import argparse
import os.path
import sys

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

if __name__ == "__main__":

	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	parser = argparse.ArgumentParser(description="convert event file")
	parser.add_argument("inputFile", type=str, metavar="inputFile", help="input file in ROOTPWA format without meta data")
	parser.add_argument("outputFile", type=str, metavar="outputFile", help="input file in ROOTPWA format with meta data")
	parser.add_argument("-b", "--binning", action='append', help="declare current bin in the form 'binningVariable;lowerBound;upperBound' (e.g. 'mass;1.0;1.1'). " +
	                                                             "You can use the argument multiple times for multiple binning variables")
	parser.add_argument("-a", "--variable", action='append', help="add an additional variable to import from the old tree. " +
	                                                              "You can use the argument multiple times for multiple binning variables")
	parser.add_argument("-l", type=str, metavar="string", dest="auxString", help="auxiliary string stored in metadata (default: output file name)")
	parser.add_argument("-t", type=str, metavar="eventsType", dest="eventsType", help="type of data (can be 'real', 'generated' or 'accepted', default: 'other')")
	parser.add_argument("-v", action="store_true", dest="debug", help="verbose; print debug output (default: false)")
	parser.add_argument("--datasetLabel", type=str, metavar="<label>", dest="datasetLabel", default="not-defined", help="Data-set label (default: 'not-defined')")
	parser.add_argument("--datasetDescription", type=str, metavar="<description>", dest="datasetDescription", default="", help="Data-set description (default: empty)")
	args = parser.parse_args()

	prodKinPartNamesObjName = "prodKinParticles"
	prodKinMomentaLeafName = "prodKinMomenta"
	decayKinPartNamesObjName = "decayKinParticles"
	decayKinMomentaLeafName = "decayKinMomenta"
	inTreeName = "rootPwaEvtTree"

	auxString = "" if not args.auxString else args.auxString

	eventTypeTranslation = { "other": pyRootPwa.core.eventMetadata.OTHER,
	                         "real": pyRootPwa.core.eventMetadata.REAL,
	                         "generated": pyRootPwa.core.eventMetadata.GENERATED,
	                         "accepted": pyRootPwa.core.eventMetadata.ACCEPTED,
	                       }
	eventsTypeString = str.lower(args.eventsType)
	if args.eventsType not in eventTypeTranslation.keys():
		pyRootPwa.utils.printErr("type '" + args.eventsType + "' is invalid as an event data type.")
		sys.exit(1)

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

	multibinBoundaries = pyRootPwa.utils.multibinBoundariesFromArgList(args.binning)
	if not multibinBoundaries:
		pyRootPwa.utils.printWarn("received no valid binning map argument")

	additionalVars = [] if not args.variable else args.variable

	fileWriter = pyRootPwa.core.eventFileWriter()
	fileWriter.initialize(outputFile,
	                      auxString,
	                      eventTypeTranslation[eventsTypeString],
	                      initialStateParticleNames,
	                      finalStateParticleNames,
	                      args.datasetLabel,
	                      args.datasetDescription,
	                      multibinBoundaries,
	                      additionalVars)

	inputTree = inputFile.Get(inTreeName)
	progressBar = pyRootPwa.utils.progressBar(0, inputTree.GetEntries(), sys.stdout)
	progressBar.start()
	i = 0
	for event in inputTree:
		initialStateMomenta = []
		for particle in event.__getattr__(prodKinMomentaLeafName):
			initialStateMomenta.append(ROOT.TVector3(particle))
		finalStateMomenta = []
		for particle in event.__getattr__(decayKinMomentaLeafName):
			finalStateMomenta.append(ROOT.TVector3(particle))
		variables = []
		for variable in additionalVars:
			variables.append(event.__getattr__(variable))
		fileWriter.addEvent(initialStateMomenta, finalStateMomenta, variables)
		progressBar.update(i)
		i += 1

	fileWriter.finalize()
