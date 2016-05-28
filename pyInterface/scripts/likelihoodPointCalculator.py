#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="pwa fit executable"
	                                )

	parser.add_argument("inputFileName", type=str, metavar="inputFile", help="path to fit result")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="integralBin", default=0, help="integral bin id of fit (default: 0)")
	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
	parser.add_argument("-P", "--cauchyPriorWidth", type=float, metavar ="WIDTH", default=0.5, help="width of half-Cauchy prior (default: 0.5)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0,
	                    help="number of input events to normalize acceptance to (default: use number of events from normalization integral file)")
	parser.add_argument("--noAcceptance", help="do not take acceptance into account (default: false)", action="store_true")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	args = parser.parse_args()

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	if args.integralBin < 0:
		pyRootPwa.utils.printErr("bin < 0 (" + str(args.integralBin) + "). Aborting...")
		sys.exit(1)
	elif args.integralBin >= len(fileManager.binList):
		pyRootPwa.utils.printErr("bin out of range (" + str(args.integralBin) + ">=" + str(len(fileManager.binList)) + "). Aborting...")
		sys.exit(1)
	multiBin = fileManager.binList[args.integralBin]
	eventAndAmpFileDict = fileManager.getEventAndAmplitudeFilePathsInBin(multiBin, pyRootPwa.core.eventMetadata.REAL)
	if not eventAndAmpFileDict:
		pyRootPwa.utils.printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)

	psIntegralPath  = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.GENERATED)
	accIntegralPath = psIntegralPath
	if not args.noAcceptance:
		accIntegralPath = fileManager.getIntegralFilePath(multiBin, pyRootPwa.core.eventMetadata.ACCEPTED)
	elif args.accEventsOverride != 0:
		# for a fit without acceptance corrections the number of events
		# the acceptance matrix is normalized to needs to be equal to
		# the number of events in the normalization matrix
		intFile = pyRootPwa.ROOT.TFile.Open(psIntegralPath, "READ")
		intMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(intFile)
		intMatrix = intMeta.getAmpIntegralMatrix()
		if args.accEventsOverride != intMatrix.nmbEvents():
			pyRootPwa.utils.printErr("incorrect number of events for normalization of integral matrix for "
			                         "a fit without acceptance (got: {:d}, expected: {:d}). Aborting...".format(args.accEventsOverride, intMatrix.nmbEvents()))
			sys.exit(1)

	result = pyRootPwa.utils.getFitResultFromFile(fitResultFileName = args.inputFileName,
	                                              fitResultTreeName = config.fitResultTreeName,
	                                              fitResultBranchName = config.fitResultBranchName)
	if not result:
		pyRootPwa.utils.printErr("could not get fit result from file '" + args.inputFileName + "'. Aborting...")
		sys.exit(1)

	waveDescThres = pyRootPwa.utils.getWaveDescThresFromFitResult(result, fileManager.getWaveDescriptions())
	if not waveDescThres:
		pyRootPwa.utils.printErr("error while getting wave names, descriptions and thresholds. Aborting...")
		sys.exit(1)

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = result.massBinCenter(),
	                                      eventAndAmpFileDict = eventAndAmpFileDict,
	                                      normIntegralFileName = psIntegralPath,
	                                      accIntegralFileName = accIntegralPath,
	                                      multiBin = multiBin,
	                                      accEventsOverride = args.accEventsOverride,
	                                      cauchy = args.cauchyPriors,
	                                      cauchyWidth = args.cauchyPriorWidth,
	                                      rank = result.rank(),
	                                      verbose = args.verbose)
	if not likelihood:
		pyRootPwa.utils.printErr("error while initializing likelihood. Aborting...")
		sys.exit(1)

	pars = []
	for parameter in likelihood.parameters():
		pars.append(result.fitParameter(parameter.parName()))

	pyRootPwa.utils.printSucc("likelihood is {: .15e}.".format(likelihood.DoEval(pars)))
