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
	bestIndex = 0
	bestMass = 0.
	bestLikeli = 0.
	for i in range(fitResultTree.GetEntries()):
		fitResultTree.GetEntry(i)
		mass = fitResult.massBinCenter()
		likeli = fitResult.logLikelihood()
		if i == 0 or abs(massBinCenter - mass) < abs(massBinCenter - bestMass):
			bestIndex = i
			bestMass = mass
			bestLikeli = likeli
		elif abs(massBinCenter - mass) == abs(massBinCenter - bestMass):
			if likeli < bestLikeli:
				bestIndex = i
				bestMass = mass
				bestLikeli = likeli
	fitResultTree.GetEntry(bestIndex)
	return pyRootPwa.core.fitResult(fitResult)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="generate phase space Monte Carlo events"
	                                )

	parser.add_argument("reactionFile", type=str, metavar="reactionFile", help="reaction config file")
	parser.add_argument("fitResult", type=str, metavar="fitResult", help="fitResult to get the production amplitudes")
	parser.add_argument("integralFile", type=str, metavar="integralFile", help="integral file")
	parser.add_argument("keyfileDirectory", type=str, metavar="keyfileDir", help="directory with all keyfiles which appear in the fitResult")
	parser.add_argument("outputFile", type=str, metavar="outputFile", help="output root file")
	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=100, help="(max) number of events to generate (default: 100)")
	parser.add_argument("-p", type=str, metavar="<particleDataTable>", dest="particleDataTableFileName", default="./particleDataTable.txt", help="path to particle data table file (default: ./particleDataTable.txt)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=123456, help="random number generator seed (default: 123456)")
	parser.add_argument("-M", type=float, metavar="#", dest="massLowerBinBoundary", help="lower boundary of mass range in MeV (overwrites values from reaction file)")
	parser.add_argument("-B", type=float, metavar="#", dest="massBinWidth", help="width of mass bin in MeV")

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

	pyRootPwa.core.particleDataTable.instance.readFile(args.particleDataTableFileName)

	config = pyRootPwa.rootPwaConfig(args.configFileName)

	if config.outputFileFormat == "root":
		# read integral matrix from ROOT file
		integralFile = pyRootPwa.ROOT.TFile.Open(args.integralFile)
		integral = pyRootPwa.core.ampIntegralMatrix(integralFile.Get("integral"))
		integralFile.Close()
	else:
		# read integral matrix from ASCII file
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

	keyfileDirectory = os.path.abspath(args.keyfileDirectory)
	if not os.path.isdir(keyfileDirectory):
		printErr("keyfile directory '" + keyfileDirectory + "' does not exist. Aborting...")
		sys.exit(1)

	fitResultFile = pyRootPwa.ROOT.TFile.Open(args.fitResult, "READ")
	if not fitResultFile:
		printErr("could not open fit result file. Aborting...")
		sys.exit(1)
	fitResultTree = fitResultFile.Get(config.fitResultTreeName)
	if not fitResultTree:
		printErr("could not find fit result tree '" + config.fitResultTreeName +
		         "' in file '" + args.fitResult + "'. Aborting...")
		sys.exit(1)
	massRange = generatorManager.getGenerator().getTPrimeAndMassPicker().massRange()
	# unit of mass is GeV in generator, and MeV in the fit result
	massBinCenter = 1000. * (massRange[0] + massRange[1]) / 2.
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

	if 'flat' in waveNames:
		waveNames.remove('flat')  # ignore flat wave

	for waveName in waveNames:
		if config.outputFileFormat == "root":
			keyfile = keyfileDirectory + "/" + waveName + ".key"
		else:
			keyfile = keyfileDirectory + "/" + waveName.replace(".amp", ".key")
		if not os.path.isfile(keyfile):
			printErr('keyfile "' + keyfile + '" does not exist. Aborting...')
			sys.exit(1)
		reflectivities.append(1 if waveName[6] == '+' else -1)
		waveIndex = fitResult.waveIndex(waveName)

		if not waveName == fitResult.prodAmpName(waveIndex)[3:]:
			printErr("mismatch between waveName '" + waveName + "' and prodAmpName '" + fitResult.prodAmpName(waveIndex)[3:] + "'. Aborting...")
		prodAmps.append(fitResult.prodAmp(waveIndex))
		waveDescription = pyRootPwa.core.waveDescription()
		waveDescription.parseKeyFile(keyfile)
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			printErr('could not construct amplitude for keyfile "' + keyfile + '".')
			sys.exit(1)
		amplitude.init()
		print(amplitude)
		waveDescriptions.append(waveDescription)
		amplitudes.append(amplitude)

	printSucc("read and constructed amplitudes for " + str(len(waveDescriptions)) + " keyfiles.")

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
							printErr('could not initialize kinematics Data. Aborting...')
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
						printErr('could not read kinematics data. Aborting...')
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

	pyRootPwa.utils.printPrintingSummary(pyRootPwa.utils.printingCounter)
