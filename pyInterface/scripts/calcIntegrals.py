#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="calculate integral matrices"
	                                )

	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="integralBin", default=-1, dest="integralBin", help="bin to be calculated (default: all)")
	parser.add_argument("-e", type=str, metavar="eventsType", default="all", dest="eventsType", help="events type to be calculated ('generated' or 'accepted', default: both)")
	parser.add_argument("-d", metavar="dataset", default=-1, dest="dataset", help="data-set id or data-set label (default: all)")
	parser.add_argument("-w", type=str, metavar="path", dest="weightsFileName", default="", help="path to MC weight file for de-weighting (default: none)")
	args = parser.parse_args()

	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	if args.integralBin >= 0:
		if args.integralBin >= len(fileManager.binList):
			pyRootPwa.utils.printErr("bin out of range (" + str(args.integralBin) + ">=" + str(len(fileManager.binList)) + "). Aborting...")
			sys.exit(1)
		binList = [ fileManager.binList[args.integralBin] ]
	else:
		binList = fileManager.binList

	eventsTypes = []
	if args.eventsType == "generated":
		eventsTypes = [ pyRootPwa.core.eventMetadata.GENERATED ]
	elif args.eventsType == "accepted":
		eventsTypes = [ pyRootPwa.core.eventMetadata.ACCEPTED ]
	elif args.eventsType == "all":
		eventsTypes = [ pyRootPwa.core.eventMetadata.GENERATED,
		                pyRootPwa.core.eventMetadata.ACCEPTED ]
	else:
		pyRootPwa.utils.printErr("Invalid events type given ('" + args.eventsType + "'). Aborting...")
		sys.exit(1)

	datasets = []
	if args.dataset == -1:
		datasets = fileManager.datasetLabels
	else:
		datasets = [args.dataset]

	for multiBin in binList:
		for eventsType in eventsTypes:
			for dataset in datasets:
				outputFileName = fileManager.getIntegralFilePath(multiBin, eventsType, dataset)
				eventAndAmpFileDict = fileManager.getEventAndAmplitudeFilePathsInBin(multiBin, eventsType, dataset)
				if not eventAndAmpFileDict:
					printErr("could not retrieve valid amplitude file list. Aborting...")
					sys.exit(1)
				printInfo("calculating integral matrix from " + str(len(eventAndAmpFileDict)) + " amplitude files:")
				if not pyRootPwa.calcIntegrals(outputFileName, eventAndAmpFileDict, multiBin, args.weightsFileName):
					printErr("integral calculation failed. Aborting...")
					sys.exit(1)
				printSucc("wrote integral to TKey '" + pyRootPwa.core.ampIntegralMatrix.integralObjectName + "' "
							+ "in file '" + outputFileName + "'")
