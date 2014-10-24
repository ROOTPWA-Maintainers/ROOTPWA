#!/usr/bin/env python

import argparse
import multiprocessing
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="Deweights ROOTPWA .root file."
	                                 )
	parser.add_argument("inputFilename", type=str, metavar="<inputFile>", help="input root file")
	parser.add_argument("-o", "--outputFile", type=str, metavar="<outputFile>", default="out.root", help="deweighted output root file (if not specified, generated automatically)")
	parser.add_argument("-f", "--weightFactor", type=float, default=1., metavar="#", dest="weightFactor", help="weight factor (default = 1")
	parser.add_argument("-s", "--seed", type=int, metavar="#", dest="seed", default=2451083, help="random number seed (default = 2451083")

	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo

	inputFile = pyRootPwa.ROOT.TFile.Open(args.inputFilename, "READ")
	if not inputFile:
		printErr("error opening input file. Aborting...")
		sys.exit(1)

	metaData = pyRootPwa.core.eventMetadata.readEventFile(inputFile)
	if metaData == 0:
		printErr("error reading metaData. Input file is not a RootPWA root file.")
	inputTree = metaData.eventTree()
	maxWeight = 0.

	for event in inputTree:
		if event.weight>maxWeight:
			maxWeight=event.weight

	printInfo("maxWeight: " + str(maxWeight))
	maxWeight *= args.weightFactor
	printInfo("adjusted maxWeight to " + str(maxWeight))

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFile, "NEW")
	if not outputFile:
		printErr("output file already exists. Aborting...")
		sys.exit(1)

	fileWriter = pyRootPwa.core.eventFileWriter()
	if fileWriter.initialize(
				                                outputFile,
				                                metaData.userString(),
				                                pyRootPwa.core.eventMetadata.GENERATED,
				                                metaData.productionKinematicsParticleNames(),
				                                metaData.decayKinematicsParticleNames(),
# TODO: FILL THESE
				                                {},
				                                []
				                                ):
		pyRootPwa.ROOT.gRandom.SetSeed(args.seed)
		acceptedEntries = 0
		overallEntries = 0

		for event in inputTree:
			normWeight = event.weight / maxWeight
			cut = pyRootPwa.ROOT.gRandom.Rndm()
			if normWeight > cut:
				fileWriter.addEvent(event.prodKinMomenta, event.decayKinMomenta)
				acceptedEntries += 1
			overallEntries += 1

		printSucc("efficiency = %.2f" % (acceptedEntries*100./overallEntries) + "%")
		printSucc("# accepted events = " + str(acceptedEntries))

		fileWriter.finalize()
	else:
		printErr("could not initialize fileWriter. Aborting...")
	inputFile.Close()
	outputFile.Close()
