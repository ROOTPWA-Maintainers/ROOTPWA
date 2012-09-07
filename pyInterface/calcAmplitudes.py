#!/usr/bin/python2.7

import argparse
import ConfigParser
import glob
import os
import sys

import pyRootPwa
import pyRootPwa.utils

if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="calculates decay amplitudes "
	                                             "for given wave for events in "
	                                             "input data files and "
	                                             "writes amplitudes to file",
	                                 formatter_class=argparse.RawTextHelpFormatter
	                                )

	parser.add_argument("-c", type=str, metavar="file", default="rootpwa.config", dest="configFileName", help="path ot config file (default: ./rootpwa.config)")
	parser.add_argument("-n", type=int, metavar="#", default=-1, dest="maxNmbEvents",  help="maximum number of events to read (default: all)")
	parser.add_argument("-o", type=str, metavar="file", default="./out.root", dest="ampFileName", help="path to amplitude file (.amp or .root format; default: ./out.root)")
	parser.add_argument("-v", action="store_true", dest="debug", help="verbose; print debug output (default: false)")
	parser.add_argument("input_files", nargs="+", help="input data file(s) (.evt or .root format)")

	arguments = parser.parse_args()

	# config file
	pyRootPwa.config = ConfigParser.ConfigParser()
	try:
		with open(arguments.configFileName, 'r') as configFile:
			pyRootPwa.config.readfp(configFile)
	except IOError:
		pyRootPwa.utils.printErr("config file could not be opened. Aborting...")
		sys.exit(1)
	except ConfigParser.Error:
		pyRootPwa.utils.printErr("config file could not be parsed. Aborting...")
		sys.exit(1)

	try:
		pdgFileName              = os.path.expanduser(os.path.expandvars(pyRootPwa.config.get('general', 'particleDataTable')))
		keyfilePattern           = os.path.expanduser(os.path.expandvars(pyRootPwa.config.get('amplitudes', 'keyfiles')))
		prodKinPartNamesObjName  = pyRootPwa.config.get('amplitudes', 'prodKinPartNamesObjName')
		prodKinMomentaLeafName   = pyRootPwa.config.get('amplitudes', 'prodKinMomentaLeafName')
		decayKinPartNamesObjName = pyRootPwa.config.get('amplitudes', 'decayKinPartNamesObjName')
		decayKinMomentaLeafName  = pyRootPwa.config.get('amplitudes', 'decayKinMomentaLeafName')
	except ConfigParser.Error:
		pyRootPwa.utils.printErr("a required entry was missing from the config file. Aborting...")
		sys.exit(1)

	keyfiles = []
	if os.path.isdir(keyfilePattern):
		keyfiles = glob.glob(keyfilePattern + "/*.key")
	elif os.path.isfile(keyfilePattern) and keyfilePattern.find(".key") > 0:
		keyfiles.append(keyfilePattern)
	else:
		pyRootPwa.utils.printErr("No keyfiles found with valid file extension. Aborting...")

	# print some info
	pyRootPwa.printCompilerInfo()
	pyRootPwa.printLibraryInfo()
	pyRootPwa.printSvnVersion()

	inputFiles = []

	# check input file extensions
	for filename in arguments.input_files:
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

#	# open .root and .evt input files
#	inTrees = []
#	inChain = pyRootPwa.ROOT.TChain(arguments.inTreeName)
#
#	for filename in inputFiles:
#		if(filename[filename.rfind(".")+1:] == "root"):
#			pyRootPwa.utils.printInfo('opening ROOT input file "' + filename + '".')
#			if(inChain.Add(filename) < 1):
#				pyRootPwa.utils.printWarn('no events in ROOT input file "' + filename + '".')
#		else:
#			with open(filename, "r") is infile:
#				prodNames = pyRootPwa.ROOT.TClonesArray("TObjString")
#				decayNames = pyRootPwa.ROOT.TClonesArray("TObjString")
#
#
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
