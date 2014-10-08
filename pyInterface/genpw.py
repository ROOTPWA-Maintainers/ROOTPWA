#!/usr/bin/env python

import argparse
import multiprocessing
import os.path
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="generate phase space Monte Carlo events"
	                                )

	parser.add_argument("reactionFile", type=str, metavar="reactionFile", help="reaction config file")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=100, help="(max) number of events to generate (default: 100)")
	parser.add_argument("-a", type=int, metavar="#", dest="maxAttempts", default=0, help="(max) number of attempts to do (default: infinity)")
	parser.add_argument("-o", type=str, metavar="<outputFile>", dest="outputFileName", default="", help="ROOT output file (if not specified, generated automatically)")
	parser.add_argument("-p", type=str, metavar="<particleDataTable>", dest="particleDataTableFileName", default="./particleDataTable.txt", help="path to particle data table file (default: ./particleDataTable.txt)")
	parser.add_argument("-c", action="store_true", dest="comgeantOutput", help="if present, a comgeant eventfile (.fort.26) is written with same naming as the root file")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=123456, help="random number generator seed (default: 123456)")
	parser.add_argument("-M", type=float, metavar="#", dest="massLowerBinBoundary", help="lower boundary of mass range in MeV (overwrites values from reaction file)")
	parser.add_argument("-B", type=float, metavar="#", dest="massBinWidth", help="width of mass bin in MeV")
	parser.add_argument("-u", "--userString", type=str, metavar="#", dest="userString", help="metadata user string", default="")
	parser.add_argument("--beamfile", type=str, metavar="<beamFile>", dest="beamFileName", help="path to beam file (overrides values from config file)")
	parser.add_argument("--noRandomBeam", action="store_true", dest="noRandomBeam", help="read the events from the beamfile sequentially")
	parser.add_argument("--randomBlockBeam", action="store_true", dest="randomBlockBeam", help="like --noRandomBeam but with random starting position")

	args = parser.parse_args()

	# print some info
	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	if args.maxAttempts and args.maxAttempts < args.nEvents:
		printWarn("Maximum attempts is smaller than the number of events. Setting it to infinity.")
		args.maxAttempts = 0

	overrideMass = (args.massLowerBinBoundary is not None) or (args.massBinWidth is not None)

	if overrideMass and not ((args.massLowerBinBoundary is not None) and (args.massBinWidth is not None)):
		printErr("'-M' and '-B' can only be set together. Aborting...")
		sys.exit(2)

	printInfo("Setting random seed to " + str(args.seed))
	pyRootPwa.core.randomNumberGenerator.instance.setSeed(args.seed)

	pyRootPwa.core.particleDataTable.instance.readFile(args.particleDataTableFileName)

	generatorManager = pyRootPwa.core.generatorManager()
	if args.beamFileName is not None:
		generatorManager.overrideBeamFile(args.beamFileName)

	if not generatorManager.readReactionFile(args.reactionFile):
		printErr("could not read reaction file. Aborting...")
		sys.exit(1)

	if overrideMass:
		generatorManager.overrideMassRange(args.massLowerBinBoundary / 1000., (args.massLowerBinBoundary + args.massBinWidth) / 1000.)
	if args.noRandomBeam:
		generatorManager.readBeamfileSequentially()
	if args.randomBlockBeam:
		generatorManager.readBeamfileSequentially()
		generatorManager.randomizeBeamfileStartingPosition()

	if not generatorManager.initializeGenerator():
		printErr("could not initialize generator. Aborting...")
		sys.exit(1)

	if args.outputFileName == "":
		if overrideMass:
			args.outputFileName = str(args.massLowerBinBoundary) + "." + str(args.massLowerBinBoundary + args.massBinWidth) + ".phaseSpace.root"
		else:
			index = 1
			filename = "1.phaseSpace.root"
			while os.path.exists(filename):
				index += 1
				filename = str(index) + ".phaseSpace.root"
			args.outputFileName = filename

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	if not outputFile:
		printErr("could not open output file. Aborting...")
		sys.exit(1)
	printInfo("opened output root file: " + args.outputFileName)

	outputComgeantFile = None
	if args.comgeantOutput:
		outputComgeantFile = open(args.outputFileName.replace(".evt", ".fort.26"), 'w')
		printInfo("opened output comgeant file: " + args.outputFileName)

	try:
		print(generatorManager)
		progressBar = pyRootPwa.utils.progressBar(0, args.nEvents, sys.stdout)
		progressBar.start()
		attempts = 0
		eventsGenerated = 0

		fileWriter = pyRootPwa.core.eventFileWriter()
		generator = generatorManager.getGenerator()
		beam = generator.getGeneratedBeam()
		finalState = generator.getGeneratedFinalState()

		prodKinNames = [ beam.name ]
		decayKinNames = [ particle.name for particle in finalState ]
		success = fileWriter.initialize(
		                                outputFile,
		                                args.userString,
		                                pyRootPwa.core.eventMetadata.GENERATED,
		                                prodKinNames,
		                                decayKinNames,
# TODO: FILL THESE
		                                {},
		                                []
		                                )
		if not success:
			printErr('could not initialize file writer. Aborting...')
			sys.exit(1)

		while eventsGenerated < args.nEvents:
			attempts += generatorManager.event()
			if args.maxAttempts and attempts > args.maxAttempts:
				printWarn("reached maximum attempts. Aborting...")
				break
			beam = generator.getGeneratedBeam()
			finalState = generator.getGeneratedFinalState()
			prodKin = [ beam.lzVec.Vect() ]
			decayKin = [ particle.lzVec.Vect() for particle in finalState ]
			fileWriter.addEvent(prodKin, decayKin)

			if args.comgeantOutput:
				recoil = generator.getGeneratedRecoil()
				vertex = generator.getGeneratedVertex()
				outputComgeantFile.write(generator.convertEventToComgeant(beam, recoil, vertex, finalState, False))

			eventsGenerated += 1
			progressBar.update(eventsGenerated)
	finally:
		fileWriter.finalize()
		if outputComgeantFile is not None:
			outputComgeantFile.close()

	printSucc("generated " + str(eventsGenerated) + " events.")
	printInfo("attempts: " + str(attempts))
	printInfo("efficiency: {0:.2%}".format(float(eventsGenerated) / float(attempts)))

	pyRootPwa.utils.printPrintingSummary(pyRootPwa.utils.printingCounter)
