#!/usr/bin/python2.7

import argparse
import ConfigParser
import glob
import os
import sys

import pyRootPwa
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
		pdgFileName                 = os.path.expanduser(os.path.expandvars(pyRootPwa.config.get('general', 'particleDataTable')))

		if pyRootPwa.config.has_option('general', 'dataDirectory'):
			dataDirectory = os.path.expanduser(os.path.expandvars(pyRootPwa.config.get('general', 'dataDirectory')))
			if dataDirectory == "":
				dataDirectory = os.path.abspath(os.path.dirname((arguments.configFileName)))
		else:
			dataDirectory = os.path.abspath(os.path.dirname((arguments.configFileName)))

		if not pyRootPwa.config.has_option('general', 'dataFileExtensionQualifier'):
			dataFileExtensionQualifier = ""
		else:
			dataFileExtensionQualifier = pyRootPwa.config.get('general', 'dataFileExtensionQualifier')

		phaseSpaceEventFileExtenisonQualifier = pyRootPwa.config.get('general', 'phaseSpaceEventFileExtenisonQualifier')
		accCorrPSEventFileExtensionQualifier  = pyRootPwa.config.get('general', 'accCorrPSEventFileExtensionQualifier')
		massBinDirectoryNamePattern           = pyRootPwa.config.get('general', 'massBinDirectoryNamePattern')
		keyfilePattern                        = os.path.expanduser(os.path.expandvars(pyRootPwa.config.get('amplitudes', 'keyfiles')))
		prodKinPartNamesObjName               = pyRootPwa.config.get('amplitudes', 'prodKinPartNamesObjName')
		prodKinMomentaLeafName                = pyRootPwa.config.get('amplitudes', 'prodKinMomentaLeafName')
		decayKinPartNamesObjName              = pyRootPwa.config.get('amplitudes', 'decayKinPartNamesObjName')
		decayKinMomentaLeafName               = pyRootPwa.config.get('amplitudes', 'decayKinMomentaLeafName')
	except ConfigParser.Error:
		pyRootPwa.utils.printErr("a required entry was missing from the config file. Aborting...")
		sys.exit(1)

	if not os.path.isdir(dataDirectory):
		pyRootPwa.utils.printErr("Data directory invalid. Aborting...")
		sys.exit(1)

	if phaseSpaceEventFileExtenisonQualifier == "":
		pyRootPwa.utils.printWarn("File extension qualifier for the phase space events is empty, no phase space events will be calculated...")
	if accCorrPSEventFileExtensionQualifier == "":
		pyRootPwa.utils.printWarn("File extension qualifier for the acceptance corrected phase space events is empty, no acc. cor. phase space events will be calculated...")
	# get the massBins as specified on the command line
	allMassBins = sorted(glob.glob(dataDirectory + '/' + massBinDirectoryNamePattern))
	massBins = []
	massBinIndizesProblem = False
	if arguments.massBins == "all":
		massBins = allMassBins
	elif arguments.massBins.find("-") > 0 or arguments.massBins.find(",") > 0:
		rawMassBinIndizes = arguments.massBins.split(",")
		massBinIndizes = []
		for massBinIndex in rawMassBinIndizes:
			if massBinIndex.find("-") > 0:
				(lb, tmp, ub) = massBinIndex.partition("-")
				try:
					lb = int(lb)
					ub = int(ub)
				except ValueError:
					massBinIndizesProblem = True
					break
				for i in range(lb, ub+1):
					massBinIndizes.append(i)
			else:
				try:
					mbi = int(massBinIndex)
				except ValueError:
					massBinIndizesProblem = True
					break
				massBinIndizes.append(mbi)
		for index in massBinIndizes:
			try:
				massBins.append(allMassBins[index-1])
			except IndexError:
				pyRootPwa.utils.printErr("Mass bin command line option out of range. Aborting...")
				sys.exit(1)
		print(massBinIndizes)
	else:
		try:
			mbi = int(arguments.massBins)
			massBins.append(allMassBins[mbi-1])
		except ValueError:
			massBinIndizesProblem = True
		except IndexError:
			pyRootPwa.utils.printErr("Mass bin command line option out of range. Aborting...")
			sys.exit(1)
	if massBinIndizesProblem:
		pyRootPwa.utils.printErr("Mass bin command line option was invalid. Aborting...")
		sys.exit(1)

	# collect all the input files
	inputDataFiles = []
	inputPSFiles = []
	inputAccPSFiles = []
	for massBin in massBins:
		inputFile = massBin + "/" + massBin.rsplit('/', 1)[-1]
		if dataFileExtensionQualifier != "":
			inputFile += "." + dataFileExtensionQualifier
		if os.path.isfile(inputFile + ".root"):
			inputFile += ".root"
		elif os.path.isfile(inputFile + ".evt"):
			inputFile += ".evt"
		else:
			pyRootPwa.utils.printErr('Mass bin "' + massBin + '" does not contain data input file "' + inputFile + '{.root/.evt}". Aborting...')
			sys.exit(1)
		inputDataFiles.append(inputFile)
		if phaseSpaceEventFileExtenisonQualifier != "":
			inputFile = massBin + "/" + massBin.rsplit('/', 1)[-1] + "." + phaseSpaceEventFileExtenisonQualifier
			if os.path.isfile(inputFile + ".root"):
				inputFile += ".root"
				inputPSFiles.append(inputFile)
			elif os.path.isfile(inputFile + ".evt"):
				inputFile += ".evt"
				inputPSFiles.append(inputFile)
			else:
				pass
#				pyRootPwa.utils.printWarn('Mass bin "' + massBin + '" does not contain phase space input file "' + inputFile + '{.root/.evt}".')
		if accCorrPSEventFileExtensionQualifier != "":
			inputFile = massBin + "/" + massBin.rsplit('/', 1)[-1] + "." + accCorrPSEventFileExtensionQualifier
			if os.path.isfile(inputFile + ".root"):
				inputFile += ".root"
				inputAccPSFiles.append(inputFile)
			elif os.path.isfile(inputFile + ".evt"):
				inputFile += ".evt"
				inputAccPSFiles.append(inputFile)
			else:
				pass
#				pyRootPwa.utils.printWarn('Mass bin "' + massBin + '" does not data contain acc. cor. phase space input file "' + inputFile + '{.root/.evt}".')




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
