#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="pwa NLopt fit executable"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="integralBin", default=0, help="integral bin id of fit (default: 0)")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random seed (default: 0)")
	parser.add_argument("-N", type=int, metavar="#", dest="nAttempts", default=1, help="number of fit attempts to perform")
	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
	parser.add_argument("-P", "--cauchyPriorWidth", type=float, metavar ="WIDTH", default=0.5, help="width of half-Cauchy prior (default: 0.5)")
	parser.add_argument("-w", type=str, metavar="path", dest="waveListFileName", default="", help="path to wavelist file (default: none)")
	parser.add_argument("-S", type=str, metavar="path", dest="startValFileName", default="", help="path to start value fit result file (default: none)")
	parser.add_argument("-r", type=int, metavar="#", dest="rank", default=1, help="rank of spin density matrix (default: 1)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0,
	                    help="number of input events to normalize acceptance to (default: use number of events from normalization integral file)")
	parser.add_argument("--noAcceptance", help="do not take acceptance into account (default: false)", action="store_true")
	parser.add_argument("-H", "--checkHessian", help="check analytical Hessian eigenvalues (default: false)", action="store_true")
	parser.add_argument("-z", "--saveSpace", help="save space by not saving integral and covariance matrices (default: false)", action="store_true")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		printErr("loading the file manager failed. Aborting...")
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
		printErr("could not retrieve valid amplitude file list. Aborting...")
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
			printErr("incorrect number of events for normalization of integral matrix for "
			         "a fit without acceptance (got: {:d}, expected: {:d}). Aborting...".format(args.accEventsOverride, intMatrix.nmbEvents()))
			sys.exit(1)

	fitResults = pyRootPwa.pwaNloptFit(
	                                   eventAndAmpFileDict = eventAndAmpFileDict,
	                                   normIntegralFileName = psIntegralPath,
	                                   accIntegralFileName = accIntegralPath,
	                                   multiBin = multiBin,
	                                   waveListFileName = args.waveListFileName,
	                                   waveDescriptions = fileManager.getWaveDescriptions(),
	                                   seed = args.seed,
	                                   cauchy = args.cauchyPriors,
	                                   cauchyWidth = args.cauchyPriorWidth,
	                                   startValFileName = args.startValFileName,
	                                   accEventsOverride = args.accEventsOverride,
	                                   checkHessian = args.checkHessian,
	                                   saveSpace = args.saveSpace,
	                                   rank = args.rank,
	                                   verbose = args.verbose,
	                                   attempts = args.nAttempts
	                                  )
	if not fitResults:
		printErr("didn't get valid fit result(s). Aborting...")
		sys.exit(1)
	printInfo("writing result(s) to '" + args.outputFileName + "'")
	valTreeName   = "pwa"
	valBranchName = "fitResult_v2"
	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "UPDATE")
	if (not outputFile) or outputFile.IsZombie():
		printErr("cannot open output file '" + args.outputFileName + "'. Aborting...")
		sys.exit(1)
	fitResult = pyRootPwa.core.fitResult()
	tree = outputFile.Get(valTreeName)
	if not tree:
		printInfo("file '" + args.outputFileName + "' is empty. "
		        + "creating new tree '" + valTreeName + "' for PWA result.")
		tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
		if not fitResult.branch(tree, valBranchName):
			printErr("failed to create new branch '" + valBranchName + "' in file '" + args.outputFileName + "'.")
			sys.exit(1)
	else:
		fitResult.setBranchAddress(tree, valBranchName)
	for result in fitResults:
		fitResult.fill(result)
		tree.Fill()
	nmbBytes = tree.Write()
	outputFile.Close()
	if nmbBytes == 0:
		printErr("problems writing fit result to TKey 'fitResult' "
		       + "in file '" + args.outputFileName + "'")
		sys.exit(1)
	else:
		printSucc("wrote fit result to TKey 'fitResult' "
		        + "in file '" + args.outputFileName + "'")
