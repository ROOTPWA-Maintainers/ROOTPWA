#!/usr/bin/env python

import argparse
import sys

import numpy

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
	printWarn = pyRootPwa.utils.printWarn
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

	additionalVariableNames = metaData.additionalTreeVariableNames()
	additionalVariables = [None] * len(additionalVariableNames)
	weightIndex = additionalVariableNames.index("weight")
	for i, additionalVariableLabel in enumerate(additionalVariableNames):
		additionalVariables[i] = numpy.array(1, dtype = float)
		inputTree.SetBranchAddress(additionalVariableLabel, additionalVariables[i])

	if metaData.hasAuxValue("weightNorm"): # stored during generation of pseudoData
		printWarn("The input event file is already deweighted!")

	if metaData.hasAuxValue("maxWeight"): # stored during generation of pseudoData
		maxWeight = metaData.auxValue("maxWeight")
	else: # determine maxWeight yourself
		maxWeight = 0.
		for i in xrange(inputTree.GetEntries()):
			inputTree.GetEntry(i)
			if float(additionalVariables[weightIndex]) > maxWeight:
				maxWeight = float(additionalVariables[weightIndex])

	printInfo("maxWeight: " + str(maxWeight))
	weightNorm = maxWeight * args.weightFactor
	printInfo("adjusted maxWeight to " + str(weightNorm))

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	if not outputFile:
		printErr("output file '" + args.outputFileName + "' already exists. Aborting...")
		sys.exit(1)

	fileWriter = pyRootPwa.core.eventFileWriter()
	if not fileWriter.initialize(outputFile,
	                             metaData.auxString(),
	                             pyRootPwa.core.eventMetadata.REAL,
	                             metaData.productionKinematicsParticleNames(),
	                             metaData.decayKinematicsParticleNames(),
	                             metaData.multibinBoundaries(),
	                             [ additionalVariableLabel for i, additionalVariableLabel in enumerate(additionalVariableNames)
	                                                           if i != weightIndex ]):
		printErr("could not initialize fileWriter. Aborting...")
		inputFile.Close()
		sys.exit(1)
	pyRootPwa.ROOT.gRandom.SetSeed(args.seed)
	acceptedEntries = 0
	overallEntries = 0

	prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
	inputTree.SetBranchAddress(pyRootPwa.core.eventMetadata.productionKinematicsMomentaBranchName, prodKinMomenta)
	inputTree.SetBranchAddress(pyRootPwa.core.eventMetadata.decayKinematicsMomentaBranchName, decayKinMomenta)

	for eventIndex in xrange(inputTree.GetEntries()):
		inputTree.GetEntry(eventIndex)
		normWeight = additionalVariables[weightIndex] / weightNorm
		cut = pyRootPwa.ROOT.gRandom.Rndm()
		if normWeight > cut:
			fileWriter.addEvent(prodKinMomenta, decayKinMomenta,
			                    [float(variable) for i, variable in enumerate(additionalVariables) if i != weightIndex])
			acceptedEntries += 1
		overallEntries += 1
	fileWriter.getEventMetadata().setAuxValue("maxWeight", maxWeight)
	fileWriter.getEventMetadata().setAuxValue("weightNorm", weightNorm)
	fileWriter.finalize()
	inputFile.Close()

	printSucc("efficiency = %.2f" % (acceptedEntries*100./overallEntries) + "%")
	printSucc("# accepted events = " + str(acceptedEntries))
	printSucc("successfully deweighted input file '" + args.inputFileName + "' and saved as '" + args.outputFileName + "'.")
