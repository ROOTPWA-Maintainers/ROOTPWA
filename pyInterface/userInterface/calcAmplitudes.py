#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

def readWaveList(waveListFileName):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	                        + "'" + waveListFileName + "'.")
	try:
		with open(waveListFileName, 'r') as waveListFile:
			waveNamesFromWavelist = []
			waveThresholds = []
			lineNmb = 0
			for line in waveListFile:
				if (line[0] == '#'):  # comments start with #
					continue
				line = line.replace('\n', '')
				lineArray = line.split(" ")
				if(len(lineArray) >= 1 and len(lineArray) <= 2):
					waveName = lineArray[0]
					if(len(lineArray) == 1):
						threshold = 0
					else:
						threshold = lineArray[1]
					waveNamesFromWavelist.append(waveName)
					waveThresholds.append(float(threshold))
				else:
					pyRootPwa.utils.printWarn("cannot parse line '" + line + "' in wave list file "
					                        + "'" + waveListFileName + "'.")
				lineNmb += 1
	except IOError:
		pyRootPwa.utils.printError("cannot open file '" + waveListFileName + "'. Aborting...")
		sys.exit(1)
	pyRootPwa.utils.printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return (waveNamesFromWavelist, waveThresholds)

if __name__ == "__main__":

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="calculates decay amplitudes "
	                                             "for given wave for events in "
	                                             "input data files and "
	                                             "writes amplitudes to file"
	                                )

	parser.add_argument("-c", type=str, metavar="configFileName", default="rootpwa.config", dest="configFileName", help="path to config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", default=-1, dest="maxNmbEvents",  help="maximum number of events to read (default: all)")
	parser.add_argument("-b", type=int, metavar="massBin", default=-1, dest="massBin", help="mass bin to be calculated (default: all)")
	parser.add_argument("-e", type=str, metavar="eventsType", default="all", dest="eventsType", help="events type to be calculated ('real', 'generated' or 'accepted', default: all)")
	parser.add_argument("-f", "--no-progress-bar", action="store_true", dest="noProgressBar", help="disable progress bars (decreases computing time)")
	parser.add_argument("-k", "--keyfileIndex", type=int, metavar="#", default=-1, help="keyfile index to calculate amplitude for (overrides settings from the config file, index from 0 to number of keyfiles - 1)")
	parser.add_argument("-w", type=str, metavar="wavelistFileName", default="", dest="wavelistFileName", help="path to wavelist file (default: none)")
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

	if not args.wavelistFileName == "" and not args.keyfileIndex == -1:
		pyRootPwa.utils.printErr("Setting both options -k and -w is conflicting. Aborting...")
		sys.exit(1)

	waveList = []
	if not args.wavelistFileName == "":
		(waveList, waveThresholds) = readWaveList(args.wavelistFileName)
	if not args.keyfileIndex == -1:
		allWaveNames = fileManager.getWaveNameList()
		if not args.keyfileIndex < len(allWaveNames):
			pyRootPwa.utils.printErr("keyfileIndex from command line argument out of range. Maximum value is " + str(len(allWaveNames)-1) + ". Aborting...")
			sys.exit(1)
		waveList = [allWaveNames[args.keyfileIndex]]
		pyRootPwa.utils.printInfo("using keyfile index " + str(args.keyfileIndex) + " resulting in the wave name '" + waveList[0] + "'.")
	if len(waveList) == 0:
		waveList = fileManager.getWaveNameList()

	binIDList = fileManager.getBinIDList()
	if (not args.massBin == -1):
		binIDList = [args.massBin]

	eventsTypes = []
	if (args.eventsType == "real"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.REAL ]
	elif (args.eventsType == "generated"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.GENERATED ]
	elif (args.eventsType == "accepted"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.ACCEPTED ]
	elif (args.eventsType == "all"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.REAL,
		                pyRootPwa.core.eventMetadata.GENERATED,
		                pyRootPwa.core.eventMetadata.ACCEPTED ]
	else:
		pyRootPwa.utils.printErr("Invalid events type given ('" + args.eventsType + "'). Aborting...")
		sys.exit(1)

	for binID in binIDList:
		for waveName in waveList:
			for eventsType in eventsTypes:
				dataFile = fileManager.getDataFile(binID, eventsType)
				if not dataFile:
					continue
				if not pyRootPwa.calcAmplitude(dataFile.dataFileName, fileManager.getKeyFile(waveName)[0], fileManager.getAmplitudeFilePath(binID, waveName, eventsType), args.maxNmbEvents, not args.noProgressBar, fileManager.getKeyFile(waveName)[1]):
					pyRootPwa.utils.printWarn("could not calculate amplitude.")
