#!/usr/bin/python2.7

import argparse
import ConfigParser
import glob
import os
import sys

import pyRootPwa
import pyRootPwa.configuration
import pyRootPwa.utils

if __name__ == "__main__":

	# print some info
	pyRootPwa.printCompilerInfo()
	pyRootPwa.printLibraryInfo()
	pyRootPwa.printSvnVersion()

	# parse command line arguments
	parser = argparse.ArgumentParser(
	                                 description="calculates decay amplitudes "
	                                             "for given wave for events in "
	                                             "input data files and "
	                                             "writes amplitudes to file",
	                                 formatter_class=argparse.RawTextHelpFormatter
	                                )

	parser.add_argument("-c", type=str, metavar="file", default="rootpwa.config", dest="configFileName", help="path ot config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", default=-1, dest="maxNmbEvents",  help="maximum number of events to read (default: all)")
	parser.add_argument("-mb", type=str, metavar="massBin(s)", default="all", dest="massBins", help="mass bins to be calculated (default: all)")
	parser.add_argument("-v", action="store_true", dest="debug", help="verbose; print debug output (default: false)")

	arguments = parser.parse_args()

	# config file
	pyRootPwa.config = pyRootPwa.configuration.rootPwaConfig(arguments.configFileName)

	# get the massBins as specified on the command line
	allMassBins = sorted(glob.glob(pyRootPwa.config.dataDirectory + '/' + pyRootPwa.config.massBinDirectoryNamePattern))
	massBins = pyRootPwa.utils.parseMassBinArgs(allMassBins, arguments.massBins)
	if not massBins:
		pyRootPwa.utils.printErr("Mass bin command line option was invalid. Aborting...")
		sys.exit(1)

	(inputDataFiles, inputPSFiles, inputAccPSFiles) = pyRootPwa.utils.getListOfInputFiles(massBins)

	sys.exit(1)



	for filename in inputFiles:
		if filename.rfind(".") < 0:
			pyRootPwa.utils.printWarn('file extension missing in filename "' + filename + '". Skipping...')
			continue
		fileExt = filename[filename.rfind(".")+1:]
		if (fileExt != "root") and (fileExt != "evt"):
			pyRootPwa.utils.printWarn('input file "' + filename + '" is neither a .root nor a .evt file. Skipping...')
			continue
		inputFiles.append(filename)

	if len(inputFiles) <= 0:
		pyRootPwa.utils.printErr("no valid input files found. Aborting...")
		sys.exit(1)


	# get list of keyfiles
	keyfiles = []
	if os.path.isdir(keyfilePattern):
		keyfiles = glob.glob(keyfilePattern + "/*.key")
	elif os.path.isfile(keyfilePattern) and keyfilePattern.find(".key") > 0:
		keyfiles.append(keyfilePattern)
	else:
		globbedKeyfiles = glob.glob(keyfilePattern)
		for keyfile in globbedKeyfiles:
			if os.path.isfile(keyfile) and keyfile.find(".key") > 0:
				keyfiles.append(keyfile)
			else:
				pyRootPwa.utils.printWarn("Keyfile " + keyfile + " is not valid. Skipping...")
	if len(keyfiles) == 0:
		pyRootPwa.utils.printErr("No keyfiles found with valid file extension. Aborting...")
		sys.exit(1)

	#get list of input files
	inputFiles = []




	# check output file extension
	if arguments.ampFileName.rfind(".") < 0:
		pyRootPwa.utils.printErr('amplitude file "' + arguments.ampFileName + '" has no file extension. aborting.')
		sys.exit(1)
	fileExt = arguments.ampFileName[arguments.ampFileName.rfind(".")+1:]
	writeRootFile = False
	if fileExt == "root":
		writeRootFile = True
	elif fileExt != "amp":
		pyRootPwa.utils.printErr('amplitude file "' + arguments.ampFileName + '" is neither a .root nor a .amp file. Aborting...')
		sys.exit(1)

	# initialize the particleDataTable
	pyRootPwa.particleDataTable.readFile(pdgFileName)

	# open .root and .evt input files
	dataTrees = []
	for filename in inputFiles:
		if(filename[filename.rfind(".")+1:] == "root"):
			pyRootPwa.printErr("Root files not yet supported. Aborting...")
			sys.exit(1)
		else:
			dataTrees.append(pyRootPwa.utils.getTreeFromEvtFile(filename))






#	# Parse the keyfile and get the amplitude
#	waveDesc = pyRootPwa.waveDescription()
#	if not waveDesc.parseKeyFile(arguments.keyFileName):
#		pyRootPwa.utils.printErr('problems reading key file "' + arguments.keyFileName + '". aborting')
#		sys.exit(1)
#
#	(waveDescConstructionSuccess, amplitude) = waveDesc.constructAmplitude()
#	if not waveDescConstructionSuccess:
#		pyRootPwa.utils.printErr('problems constructing decay topology from key file. aborting.')
#		sys.exit(1)
#
#	printString = 'creating amplitude file "' + arguments.ampFileName + '"'
#	if writeRootFile:
#		pyRootPwa.utils.printInfo(printString)
#		pyRootPwa.ROOT.TFile.Open(arguments.ampFileName, "RECREATE")
#		amplitude.decayTopology()
#		waveName = waveDesc.waveNameFromTopology(amplitude.decayTopology(), arguments.newKeyFileNameConvention, None)
#		waveDesc.Write(waveName)
#	else:
#		if arguments.asciiOutput:
#			printString += "; ASCII mode"
#		else:
#			printString += "; binary mode"
#
#		pyRootPwa.utils.printInfo(printString)
#		pass
#
