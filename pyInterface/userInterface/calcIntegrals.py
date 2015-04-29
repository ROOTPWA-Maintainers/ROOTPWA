#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="generate phase space Monte Carlo events"
	                                )

	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-i", type=str, metavar="integralName", dest="integralName", default="integral", help="integral TKey name (default: 'integral')")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=0, help="maximum number of events to process (default: all)")
# 	parser.add_argument("-r", type=int, metavar="#", dest="nEventsRenorm", default=0, help="number of events to renormalize to (default: no renormalization)")
	parser.add_argument("-w", type=str, metavar="path", dest="weightsFileName", default="", help="path to MC weight file for de-weighting (default: none)")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	args = parser.parse_args()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

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

	for binID in fileManager.getBinIDList():
		for eventsType in [ pyRootPwa.core.eventMetadata.GENERATED, pyRootPwa.core.eventMetadata.ACCEPTED ]:
			ampFileList = fileManager.getAmpFilePaths(binID, eventsType)
			if not ampFileList:
				printErr("could not retrieve valid amplitude file list. Aborting...")
				sys.exit(1)
			printInfo("calculating integral matrix from " + str(len(ampFileList)) + " amplitude files:")
			integral = pyRootPwa.calcIntegrals(ampFileList, args.nEvents, args.weightsFileName)

			outputFileName = fileManager.getIntegralFilePath(binID, eventsType)
			outputFile = pyRootPwa.ROOT.TFile.Open(outputFileName, "RECREATE")
			if not outputFile:
				printErr("cannot open output file '" + outputFileName + "'. Aborting...")
				sys.exit(1)
			nmbBytes = integral.Write(args.integralName)
			outputFile.Close()
			if nmbBytes == 0:
				printErr("problems writing integral to TKey '" + args.integralName + "' "
				         + "in file '" + outputFileName + "'")
				continue
			else:
				printSucc("wrote integral to TKey '" + args.integralName + "' "
				          + "in file '" + outputFileName + "'")
