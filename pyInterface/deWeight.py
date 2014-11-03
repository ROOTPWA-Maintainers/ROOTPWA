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
	parser.add_argument("inputFileName", type=str, metavar="<inputFile>", help="input RootPwa file")
	parser.add_argument("outputFileName", type=str, metavar="<outputFile>", help="deweighted output RootPwa file")
	parser.add_argument("-f", "--weightFactor", type=float, default=1., metavar="#", help="weight factor (default = 1)")
	parser.add_argument("-s", "--seed", type=int, metavar="#", dest="seed", default=2451083, help="random number seed (default = 2451083)")

	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo

	inputFile = pyRootPwa.ROOT.TFile.Open(args.inputFileName, "READ")
	if not inputFile:
		printErr("error opening input file '" + args.inputFileName + "'. Aborting...")
		sys.exit(1)

	metaData = pyRootPwa.core.eventMetadata.readEventFile(inputFile)
	if metaData == 0:
		printErr("error reading metaData. Input file is not a RootPWA root file.")
	inputTree = metaData.eventTree()
	maxWeight = 0.

	for event in inputTree:
		if event.weight > maxWeight:
			maxWeight = event.weight

	printInfo("maxWeight: " + str(maxWeight))
	maxWeight *= args.weightFactor
	printInfo("adjusted maxWeight to " + str(maxWeight))

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	if not outputFile:
		printErr("output file '" + args.outputFileName + "' already exists. Aborting...")
		sys.exit(1)

	fileWriter = pyRootPwa.core.eventFileWriter()
	if not fileWriter.initialize(outputFile,
	                             metaData.userString(),
	                             pyRootPwa.core.eventMetadata.GENERATED,
	                             metaData.productionKinematicsParticleNames(),
	                             metaData.decayKinematicsParticleNames(),
	                             metaData.binningMap(),
	                             ["weight"]):
		printErr("could not initialize fileWriter. Aborting...")
		inputFile.Close()
		sys.exit(1)
	pyRootPwa.ROOT.gRandom.SetSeed(args.seed)
	acceptedEntries = 0
	overallEntries = 0

	for event in inputTree:
		normWeight = event.weight / maxWeight
		cut = pyRootPwa.ROOT.gRandom.Rndm()
		if normWeight > cut:
			fileWriter.addEvent(event.prodKinMomenta, event.decayKinMomenta, [event.weight])
			acceptedEntries += 1
		overallEntries += 1
	fileWriter.finalize()
	inputFile.Close()

	printSucc("efficiency = %.2f" % (acceptedEntries*100./overallEntries) + "%")
	printSucc("# accepted events = " + str(acceptedEntries))
	printSucc("successfully deweighted input file '" + args.inputFileName + "' and saved as '" + args.outputFileName + "'.")
