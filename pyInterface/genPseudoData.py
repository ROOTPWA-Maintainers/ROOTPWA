#!/usr/bin/env python

import argparse
import math
import sys

import pyRootPwa
import pyRootPwa.core

def norm(cNumber):
	return cNumber.real*cNumber.real + cNumber.imag*cNumber.imag


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="generate phase space Monte Carlo events"
	                                )

	parser.add_argument("reactionFile", type=str, metavar="reactionFile", help="reaction config file")
	parser.add_argument("fitResult", type=str, metavar="fitResult", help="fitResult to get the production amplitudes")
	parser.add_argument("integralFile", type=str, metavar="integralFile", help="integral file")
	parser.add_argument("outputFile", type=str, metavar="outputFile", help="output root file")
	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName",
	                    help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=100, help="(max) number of events to generate (default: 100)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random number generator seed (default: 0)")
	parser.add_argument("-M", type=float, metavar="#", dest="massLowerBinBoundary",
	                    help="lower boundary of mass range in MeV (overwrites values from reaction file)")
	parser.add_argument("-B", type=float, metavar="#", dest="massBinWidth", help="width of mass bin in MeV")
	parser.add_argument("-u", "--userString", type=str, metavar="#", dest="userString", help="metadata user string", default="")

	args = parser.parse_args()

	# print some info
	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.instance.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	pyRootPwa.core.integralTableContainer.setDirectory(config.phaseSpaceIntegralDirectory)
	pyRootPwa.core.integralTableContainer.setUpperMassBound(config.phaseSpaceUpperMassBound)

	# read integral matrix from ROOT file
	integralFile = pyRootPwa.ROOT.TFile.Open(args.integralFile)
	integral = pyRootPwa.core.ampIntegralMatrix.getFromTDirectory(integralFile, pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	integralFile.Close()
	nmbNormEvents = integral.nmbEvents()

	overrideMass = (args.massLowerBinBoundary is not None) or (args.massBinWidth is not None)

	if overrideMass and not ((args.massLowerBinBoundary is not None) and (args.massBinWidth is not None)):
		printErr("'-M' and '-B' can only be set together. Aborting...")
		sys.exit(2)

	printInfo("Setting random seed to " + str(args.seed))
	pyRootPwa.core.randomNumberGenerator.instance.setSeed(args.seed)

	generatorManager = pyRootPwa.core.generatorManager()

	if not generatorManager.readReactionFile(args.reactionFile):
		printErr("could not read reaction file. Aborting...")
		sys.exit(1)

	if overrideMass:
		generatorManager.overrideMassRange(args.massLowerBinBoundary / 1000., (args.massLowerBinBoundary + args.massBinWidth) / 1000.)

	if not generatorManager.initializeGenerator():
		printErr("could not initialize generator. Aborting...")
		sys.exit(1)

	massRange = generatorManager.getGenerator().getTPrimeAndMassPicker().massRange()
	# unit of mass is GeV in generator, and MeV in the fit result
	massBinCenter = 1000. * (massRange[0] + massRange[1]) / 2.
	fitResult = pyRootPwa.utils.getBestFitResultFromFile(fitResultFileName = args.fitResult,
	                                                     massBinCenter = massBinCenter,
	                                                     fitResultTreeName = config.fitResultTreeName,
	                                                     fitResultBranchName = config.fitResultBranchName)
	if not fitResult:
		printErr("could not find fit result in file '" + args.fitResult +
		         "' for mass bin " + str(massBinCenter) + ". Aborting...")
		sys.exit(1)

	waveDescriptions = []
	amplitudes = []
	reflectivities = []
	prodAmps = []
	waveNames = []
	waveNames = fitResult.waveNames()

	if fitResult.nmbProdAmps() != len(waveNames):
		pyRootPwa.utils.printErr('Number of production amplitudes (' + str(fitResult.nmbProdAmps()) +
		                         ') not equal to number of wave names (' + str(len(waveNames)) + '). Aborting...')
		sys.exit(1)

	if 'flat' in waveNames:
		waveNames.remove('flat')  # ignore flat wave

	for waveName in waveNames:
		waveDescription = fileManager.getWaveDescription(waveName)
		reflectivities.append(pyRootPwa.core.partialWaveFitHelper.getReflectivity(waveName))
		waveIndex = fitResult.waveIndex(waveName)

		if not waveName == fitResult.prodAmpName(waveIndex)[3:]:
			printErr("mismatch between waveName '" + waveName + "' and prodAmpName '" + fitResult.prodAmpName(waveIndex)[3:] + "'. Aborting...")
		prodAmps.append(fitResult.prodAmp(waveIndex))
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			printErr('could not construct amplitude for wave "' + waveName + '".')
			sys.exit(1)
		amplitude.init()
		printInfo(amplitude)
		waveDescriptions.append(waveDescription)
		amplitudes.append(amplitude)

	printSucc("read and constructed amplitudes for " + str(len(waveDescriptions)) + " keyfiles.")

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFile, "NEW")
	if not outputFile:
		printErr("could not open output file. Aborting...")
		sys.exit(1)

	fileWriter = pyRootPwa.core.eventFileWriter()

	printInfo("opened output root file: " + args.outputFile)
	try:
		printInfo(generatorManager)
		progressBar = pyRootPwa.utils.progressBar(0, args.nEvents, sys.stdout)
		progressBar.start()
		attempts = 0
		decayKin = None
		prodKin = None
		first = True
		weight = 0.

		eventsGenerated = 0
		for eventsGenerated in range(args.nEvents):

			attempts += generatorManager.event()
			generator = generatorManager.getGenerator()
			beam = generator.getGeneratedBeam()
			finalState = generator.getGeneratedFinalState()

			if first:
				prodKinNames = [ beam.name ]
				decayKinNames = [ particle.name for particle in finalState]
				success = fileWriter.initialize(
				                                outputFile,
				                                args.userString,
				                                pyRootPwa.core.eventMetadata.REAL,
				                                prodKinNames,
				                                decayKinNames,
# TODO: FILL THESE
				                                { "mass": (1000. * massRange[0], 1000. * massRange[1]) },
				                                ["weight"]
				                                )
				if not success:
					printErr('could not initialize file writer. Aborting...')
					sys.exit(1)

				for amplitude in amplitudes:
					topo = amplitude.decayTopology()
					if not topo.initKinematicsData(prodKinNames, decayKinNames):
						printErr('could not initialize kinematics Data. Aborting...')
						sys.exit(1)
				first = False

			prodKin = [ beam.lzVec.Vect() ]
			decayKin = [ particle.lzVec.Vect() for particle in finalState ]

			posReflAmpSum = 0
			negReflAmpSum = 0

			for i, _ in enumerate(amplitudes):
				amplitude = amplitudes[i]
				topo = amplitude.decayTopology()
				if not topo.readKinematicsData(prodKin, decayKin):
					progressBar.cancel()
					printErr('could not read kinematics data. Aborting...')
					sys.exit(1)

				amp = (amplitude() * prodAmps[i]) / math.sqrt(integral.element(waveNames[i], waveNames[i]).real * nmbNormEvents)
				if reflectivities[i] > 0:
					posReflAmpSum += amp
				else:
					negReflAmpSum += amp

			weight = norm(posReflAmpSum) + norm(negReflAmpSum)
			fileWriter.addEvent(prodKin, decayKin, [weight])

			progressBar.update(eventsGenerated)
	except:
		raise
	finally:
		fileWriter.finalize()

	eventsGenerated += 1
	printSucc("generated " + str(eventsGenerated) + " events.")
	printInfo("attempts: " + str(attempts))
	printInfo("efficiency: " + str(100. * (float(eventsGenerated) / float(attempts))) + "%")

	pyRootPwa.utils.printPrintingSummary(pyRootPwa.utils.printingCounter)
