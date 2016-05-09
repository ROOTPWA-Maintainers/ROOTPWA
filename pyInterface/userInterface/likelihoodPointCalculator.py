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
	parser.add_argument("-b", type=int, metavar="#", dest="binID", default=0, help="bin ID of fit (default: 0)")
	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
	parser.add_argument("-P", "--cauchyPriorWidth", type=float, metavar ="WIDTH", default=0.5, help="width of half-Cauchy prior (default: 0.5)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0,
	                    help="number of input events to normalize acceptance to (default: use number of events from normalization integral file)")
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

	ampFileList = fileManager.getAmplitudeFilePaths(args.binID, pyRootPwa.core.eventMetadata.REAL)
	if not ampFileList:
		pyRootPwa.utils.printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)
	binningMap = fileManager.getBinFromID(args.binID)

	psIntegralPath  = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.GENERATED)
	accIntegralPath = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.ACCEPTED)

	result = pyRootPwa.utils.getFitResultFromFile(fitResultFileName = args.inputFileName,
	                                              fitResultTreeName = config.fitResultTreeName,
	                                              fitResultBranchName = config.fitResultBranchName)
	if not result:
		pyRootPwa.utils.printErr("could not get fit result from file '" + args.inputFileName + "'. Aborting...")
		sys.exit(1)

	waveDescThres = pyRootPwa.utils.getWaveDescThresFromFitResult(result, fileManager.getKeyFiles())
	if not waveDescThres:
		pyRootPwa.utils.printErr("error while getting wave names, descriptions and thresholds. Aborting...")
		sys.exit(1)

	likelihood = pyRootPwa.initLikelihood(waveDescThres = waveDescThres,
	                                      massBinCenter = result.massBinCenter(),
	                                      ampFileList = ampFileList,
	                                      normIntegralFileName = psIntegralPath,
	                                      accIntegralFileName = accIntegralPath,
	                                      accEventsOverride = args.accEventsOverride,
	                                      cauchy = args.cauchyPriors,
	                                      cauchyWidth = args.cauchyPriorWidth,
	                                      rank = result.rank(),
	                                      verbose = args.verbose)

	pars = []
	for i in range(likelihood.nmbPars()):
		parName = likelihood.parName(i)
		pars.append(result.fitParameter(parName))

	pyRootPwa.utils.printSucc("likelihood is {: .15e}.".format(likelihood.DoEval(pars)))
