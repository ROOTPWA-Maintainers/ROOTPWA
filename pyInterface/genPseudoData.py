#!/usr/bin/env python

import argparse
import math
import multiprocessing
import numpy
import os.path
import sys

import pyRootPwa

def norm(c):
	return (c.real*c.real + c.imag*c.imag)

def getBestFitResult(massBinCenter, fitResultTree):
	fitResult = pyRootPwa.core.fitResult().getAsRootObject()
	fitResultTree.SetBranchAddress(config.fitResultBranchName, fitResult)
	massBinCenterBest = 0.
	bestIndex = 0
	bestLikeli = 0.
	for i in range(fitResultTree.GetEntries()):
		fitResultTree.GetEntry(i)
		mass = fitResult.massBinCenter()
		logLike = fitResult.logLikelihood()
		if i == 0 or abs(massBinCenter - mass) < abs(massBinCenter - massBinCenterBest):
			bestIndex = i
			bestLikeli = logLike
			massBinCenterBest = mass
		elif abs(massBinCenter - mass) == abs(massBinCenter - massBinCenterBest):
			if logLike < bestLike:
				bestIndex = i
				massBinCenterBest = mass
				logLikeBest = logLike
	fitResultTree.GetEntry(bestIndex)
	return pyRootPwa.core.fitResult(fitResult)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="generate phase space Monte Carlo events"
	                                )

	parser.add_argument("reactionFile", type=str, metavar="reactionFile", help="reaction config file")
	parser.add_argument("fitResult", type=str, metavar="fitResult", help="fitResult to get the production amplitudes")
	parser.add_argument("integralFile", type=str, metavar="integralFile", help="integral file")
	parser.add_argument("outputFile", type=str, metavar="outputFile", help="output root file")
	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=100, help="(max) number of events to generate (default: 100)")
	parser.add_argument("-p", type=str, metavar="<particleDataTable>", dest="particleDataTableFileName", default="./particleDataTable.txt", help="path to particle data table file (default: ./particleDataTable.txt)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=123456, help="random number generator seed (default: 123456)")
	parser.add_argument("-M", type=float, metavar="#", dest="massLowerBinBoundary", help="lower boundary of mass range in MeV (overwrites values from reaction file)")
	parser.add_argument("-B", type=float, metavar="#", dest="massBinWidth", help="width of mass bin in MeV")
	parser.add_argument("-k", "--keyfiles", type=str, metavar="keyfiles", dest="keyfiles", nargs="*", help="keyfiles to calculate amplitude for (overrides settings from the config file)")

	args = parser.parse_args()

	# print some info
	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	#initialize the printing functors
	printingCounter = multiprocessing.Array('i', [0]*5)
	pyRootPwa.utils.printErr = pyRootPwa.utils.printErrClass(printingCounter)
	pyRootPwa.utils.printWarn = pyRootPwa.utils.printWarnClass(printingCounter)
	pyRootPwa.utils.printSucc = pyRootPwa.utils.printSuccClass(printingCounter)
	pyRootPwa.utils.printInfo = pyRootPwa.utils.printInfoClass(printingCounter)
	pyRootPwa.utils.printDebug = pyRootPwa.utils.printDebugClass(printingCounter)

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	pyRootPwa.core.particleDataTable.instance.readFile(args.particleDataTableFileName)

	config = pyRootPwa.rootPwaConfig(args.configFileName)

	integral = pyRootPwa.core.ampIntegralMatrix()
	if not integral.readAscii(args.integralFile):
		printErr("Cannot read normalization integral from file '" + args.integralFile + "'. Aborting...")
		sys.exit(1)
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

	keyfiles = []
	if args.keyfiles is not None:
		keyfiles = args.keyfiles
		for keyfile in keyfiles:
			if not (os.path.isfile(keyfile) and keyfile.endswith(".key")):
				pyRootPwa.utils.printErr("Keyfile '" + keyfile + "' not valid. Aborting...")
				sys.exit(1)
	else:
		keyfiles = pyRootPwa.utils.getListOfKeyfiles(pyRootPwa.config.keyfilePattern)
		if len(keyfiles) == 0:
			pyRootPwa.utils.printErr("No keyfiles found with valid file extension. Aborting...")
			sys.exit(1)

	fitResultFile = pyRootPwa.ROOT.TFile.Open(args.fitResult, "READ")
	if not fitResultFile:
		printErr("Could not open fit result file. Aborting...")
		sys.exit(1)
	fitResultTree = fitResultFile.Get(config.fitResultTreeName)
	if not fitResultTree:
		printErr("Could not find fit result tree '" + config.fitResultTreeName +
		         "' in file '" + args.fitResult + "'. Aborting...")
		sys.exit(1)
	massBinCenterBest = 0
	massRange = generatorManager.getGenerator().getTPrimeAndMassPicker().massRange()
	massBinCenter = (massRange[0] + massRange[1]) / 2.
	fitResult = getBestFitResult(massBinCenter, fitResultTree)

	waveDescriptions = []
	amplitudes = []
	reflectivities = []
	prodAmps = []
	waveNames = []
	waveNames = fitResult.waveNames()

	if(fitResult.nmbProdAmps() != len(waveNames)):
		pyRootPwa.utils.printErr('Number of production amplitudes (' + str(fitResult.nmbProdAmps()) +
		                         ') not equal to number of wave names (' + str(len(waveNames)) + '). Aborting...')
		sys.exit(1)

	waveNames.remove('flat')  # ignore flat wave

	for waveName in waveNames:
		keyfile = 'keyfiles/' + waveName.replace(".amp", ".key")
		if not os.path.isfile(keyfile):
			pyRootPwa.utils.printErr('Keyfile "' + keyfile + '" does not exist. Aborting...')
			sys.exit(1)
		reflectivities.append(1 if waveName[6] == '+' else -1)
		waveIndex = fitResult.waveIndex(waveName)

		if not waveName == fitResult.prodAmpName(waveIndex)[3:]:
			printErr("Mismatch between waveName '" + waveName + "' and prodAmpName '" + fitResult.prodAmpName(waveIndex)[3:] + "'. Aborting...")
		prodAmps.append(fitResult.prodAmp(waveIndex))
		waveDescription = pyRootPwa.core.waveDescription()
		waveDescription.parseKeyFile(keyfile)
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			pyRootPwa.utils.printErr('Could not construct amplitude for keyfile "' + keyfile + '".')
			sys.exit(1)
		amplitude.init()
		print(amplitude)
		waveDescriptions.append(waveDescription)
		amplitudes.append(amplitude)

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFile, "NEW")
	if not outputFile:
		printErr("could not open output file. Aborting...")
		sys.exit(1)

	printInfo("opened output root file: " + args.outputFile)
	try:
		print(generatorManager)
		progressBar = pyRootPwa.utils.progressBar(0, args.nEvents-1, sys.stdout)
		progressBar.start()
		attempts = 0
		decayKin = None
		prodKin = None
		first = True
		outTree = None
		weight = numpy.zeros(1, dtype = float)

		for eventsGenerated in range(args.nEvents):

				attempts += generatorManager.event()
				generator = generatorManager.getGenerator()
				beam = generator.getGeneratedBeam()
				finalState = generator.getGeneratedFinalState()

				if first:
					nDecayParticles = len(finalState)
					decayKinNames = pyRootPwa.ROOT.TClonesArray("TObjString", len(finalState))
					for i in range(len(finalState)):
						decayKinNames[i] = pyRootPwa.ROOT.TObjString(finalState[i].name)
					prodKinNames = pyRootPwa.ROOT.TClonesArray("TObjString", 1)
					prodKinNames[0] = pyRootPwa.ROOT.TObjString(beam.name)
					for amplitude in amplitudes:
						topo = amplitude.decayTopology()
						if not topo.initKinematicsData(prodKinNames, decayKinNames):
							pyRootPwa.utils.printErr('Could not initialize kinematics Data. Aborting...')
							sys.exit(1)
					decayKin = pyRootPwa.ROOT.TClonesArray("TVector3", len(finalState))
					prodKin = pyRootPwa.ROOT.TClonesArray("TVector3", 1)
					outTree = pyRootPwa.ROOT.TTree(config.inTreeName, config.inTreeName)
					outTree.Branch(config.prodKinMomentaLeafName, "TClonesArray", prodKin, 256000, 99)
					outTree.Branch(config.decayKinMomentaLeafName, "TClonesArray", decayKin, 256000, 99)
					outTree.Branch("weight", weight, "weight/D")
					prodKinNames.Write(config.prodKinPartNamesObjName, pyRootPwa.ROOT.TObject.kSingleKey)
					decayKinNames.Write(config.decayKinPartNamesObjName, pyRootPwa.ROOT.TObject.kSingleKey)
					first = False

				for i in range(len(finalState)):
					decayKin[i] = finalState[i].lzVec.Vect()
				prodKin[0] = beam.lzVec.Vect()

				posReflAmpSum = 0
				negReflAmpSum = 0

				for i in range(len(amplitudes)):
					amplitude = amplitudes[i]
					topo = amplitude.decayTopology()
					if not topo.readKinematicsData(prodKin, decayKin):
						progressBar.cancel()
						pyRootPwa.utils.printErr('Could not read kinematics data. Aborting...')
						sys.exit(1)

					amp = (amplitude() * prodAmps[i]) / math.sqrt(integral.element(waveNames[i], waveNames[i]).real * nmbNormEvents)
					if reflectivities[i] > 0:
						posReflAmpSum += amp
					else:
						negReflAmpSum += amp

				weight[0] = norm(posReflAmpSum) + norm(negReflAmpSum)
				outTree.Fill()

				progressBar.update(eventsGenerated)
	except:
		raise
	else:
		outTree.Write()
	finally:
		outputFile.Close()

	eventsGenerated += 1
	printSucc("generated " + str(eventsGenerated) + " events.")
	printInfo("attempts: " + str(attempts))
	printInfo("efficiency: " + str(100. * (float(eventsGenerated) / float(attempts))) + "%")

	pyRootPwa.utils.printPrintingSummary(printingCounter)
