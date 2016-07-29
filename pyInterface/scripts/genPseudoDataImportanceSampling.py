#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="generate phase-space Monte Carlo events"
	                                )

	parser.add_argument("reactionFile", type=str, metavar="reactionFile", help="reaction config file")
	parser.add_argument("fitResult", type=str, metavar="fitResult", help="fitResult to get the production amplitudes")
	parser.add_argument("outputFile", type=str, metavar="outputFile", help="output root file")
	parser.add_argument("-c", type=str, metavar="config-file", default="rootpwa.config", dest="configFileName",
	                    help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-i", "--integralFile", type=str, metavar="integralFile", help="integral file")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=100, help="(max) number of events to generate (default: 100)")
	parser.add_argument("-C", type=int, metavar="#", dest="nChains", default=5, help="number od MCMC chains to be run")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random number generator seed (default: 0)")
	parser.add_argument("-l", type=int, metavar="#", dest="lag", default=1, help="take just every nth generated event into the sample to avoid correlations")
	parser.add_argument("-M", type=float, metavar="#", dest="massLowerBinBoundary",
	                    help="lower boundary of mass range in MeV (!) (overwrites values from reaction file)")
	parser.add_argument("-B", type=float, metavar="#", dest="massBinWidth", help="width of mass bin in MeV (!)")
	parser.add_argument("-u", "--userString", type=str, metavar="#", dest="userString", help="metadata user string", default="importanceSampledEvents")
	parser.add_argument("--massTPrimeVariableNames", type=str, dest="massTPrimeVariableNames", help="Name of the mass and t' variable (default: %(default)s)",
	                    default="mass,tPrime")
	parser.add_argument("--noStoreMassTPrime", action="store_true", dest="noStoreMassTPrime", help="Do not store mass and t' variable of each event.")
	parser.add_argument("-p", "--phaseSpaceOnly", help="do phase space only (default: false)", action="store_true")
	parser.add_argument("--beamfile", type=str, metavar="<beamFile>", dest="beamFileName", help="path to beam file (overrides values from config file)")
	parser.add_argument("--noRandomBeam", action="store_true", dest="noRandomBeam", help="read the events from the beamfile sequentially")
	parser.add_argument("--randomBlockBeam", action="store_true", dest="randomBlockBeam", help="like --noRandomBeam but with random starting position")

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

	overrideMass = (args.massLowerBinBoundary is not None) or (args.massBinWidth is not None)

	if overrideMass and not ((args.massLowerBinBoundary is not None) and (args.massBinWidth is not None)):
		printErr("'-M' and '-B' can only be set together. Aborting...")
		sys.exit(2)

	printInfo("Setting random seed to " + str(args.seed))
	pyRootPwa.core.randomNumberGenerator.instance.setSeed(args.seed)

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

	massRange = generatorManager.getGenerator().getTPrimeAndMassPicker().massRange()
	massBinCenter = (massRange[0] + massRange[1]) / 2.
	fitResult = pyRootPwa.utils.getBestFitResultFromFile(fitResultFileName = args.fitResult,
	                                                     massBinCenter = massBinCenter,
	                                                     fitResultTreeName = config.fitResultTreeName,
	                                                     fitResultBranchName = config.fitResultBranchName)
	if not fitResult:
		printErr("could not find fit result in file '" + args.fitResult +
		         "' for mass bin " + str(massBinCenter) + ". Aborting...")
		sys.exit(1)

	waveNames = fitResult.waveNames()
	if 'flat' in waveNames:
		waveNames.remove('flat')  # ignore flat wave

	model = pyRootPwa.core.modelIntensity(fitResult)
	for waveName in waveNames:
		waveDescription = fileManager.getWaveDescription(waveName)
		(result, amplitude) = waveDescription.constructAmplitude()
		if not result:
			printErr('could not construct amplitude for wave "' + waveName + '".')
			sys.exit(1)
		if not model.addDecayAmplitude(amplitude):
			printErr('could not add amplitude for wave "' + waveName + '".')

	# overwrite integral matrix from fit result with one read from a file
	if args.integralFile:
		integralFile = pyRootPwa.ROOT.TFile.Open(args.integralFile, "READ")
		integralMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(integralFile)
		integral = integralMeta.getAmpIntegralMatrix()
		model.loadPhaseSpaceIntegral(integral)

	# do not let BAT create histograms in the eventfile, otherwise this script
	# will exit with a segmentation violation due to ROOT ownership
	pyRootPwa.ROOT.TH1.AddDirectory(False)

	modelSampler = generatorManager.getImportanceSampler(model)
	if args.phaseSpaceOnly:
		modelSampler.setPhaseSpaceOnly()

	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFile, "NEW")
	if not outputFile:
		printErr("could not open output file. Aborting...")
		sys.exit(1)
	if len(args.massTPrimeVariableNames.split(',')) != 2:
		printErr("Option --massTPrimeVariableNames has wrong format '" + args.massTPrimeVariableNames + "'. Aborting...")
		sys.exit(1)
	success = modelSampler.initializeFileWriter(outputFile,
	                                            args.userString,
	                                            not args.noStoreMassTPrime,
	                                            args.massTPrimeVariableNames.split(',')[0],
	                                            args.massTPrimeVariableNames.split(',')[1])
	if not success:
		printErr('could not initialize file writer. Aborting...')
		sys.exit(1)

	modelSampler.SetNChains(args.nChains)
	modelSampler.SetNIterationsRun(args.nEvents/args.nChains*args.lag)
	modelSampler.SetNLag(args.lag)
	modelSampler.SetRandomSeed(args.seed)
	modelSampler.MarginalizeAll()

	modelSampler.finalizeFileWriter()

	modelSampler.PrintAllMarginalized(args.outputFile.replace(".root",".pdf").replace(".ROOT",".pdf"),2,4)
	modelSampler.PrintCorrelationMatrix(args.outputFile.replace(".root","_coma.pdf").replace(".ROOT","_coma.pdf"))

	modelSampler.printFuncInfo()

	realEfficiency = float(args.nEvents)/float(modelSampler.nCalls())
	printInfo("generated "+str(args.nEvents)+" with "+str(args.nChains)+" chains, needed "+str(modelSampler.nCalls()) + " actual calls => efficiency:"+str(realEfficiency))
	if realEfficiency < .5 and realEfficiency > .3:
		printSucc("efficiency in acceptable range")
	elif realEfficiency <= .3:
		printWarn("low efficiency encountered. Try decreasing the lag to speed up the process")
	elif realEfficiency >= .5 and realEfficiency <= 1.:
		printWarn("high efficiency encountered. Try increasing the lag to avoid correlations")
	elif realEfficiency > 1.:
		printErr("efficiency > 1 encountered. Data are highly correlated. Increase the lag.")
