#!/usr/bin/python2.7

import argparse
import ConfigParser
import glob
import os
import sys

import pyRootPwa
import pyRootPwa.amplitude
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

	# get the lists of input files
	(inputDataFiles, inputPSFiles, inputAccPSFiles) = pyRootPwa.utils.getListOfInputFiles(massBins)
	if not inputDataFiles:
		pyRootPwa.utils.printErr('No input data files found. Aborting...')
		sys.exit(1)
	if len(inputDataFiles) != len(massBins):
		pyRootPwa.utils.printErr('Not all data input files are present (should have ' + str(len(massBins)) + ', found ' + str(len(inputDataFiles)) + '). Aborting...')
		sys.exit(1)
	if len(inputPSFiles) != len(massBins):
		pyRootPwa.utils.printWarn('Not all phase space input files are present (should have ' + str(len(massBins)) + ', found ' + str(len(inputPSFiles)) + ').')
	if len(inputAccPSFiles) != len(massBins):
		pyRootPwa.utils.printWarn('Not all acceptance corrected phase space input files are present (should have ' + str(len(massBins)) + ', found ' + str(len(inputAccPSFiles)) + ').')


	# get list of keyfiles
	keyfiles = pyRootPwa.utils.getListOfKeyfiles(pyRootPwa.config.keyfilePattern)
	if len(keyfiles) == 0:
		pyRootPwa.utils.printErr("No keyfiles found with valid file extension. Aborting...")
		sys.exit(1)

	# initialize the particleDataTable
	pyRootPwa.particleDataTable.readFile(pyRootPwa.config.pdgFileName)

	for inputFile in (inputDataFiles + inputPSFiles + inputAccPSFiles):

		# check and if necessary create the output direcotries
		outDir = inputFile.rsplit('/', 1)[0] + '/'
		if inputFile in inputDataFiles:
			outDir += pyRootPwa.config.dataAmplitudeDirectoryName
		elif inputFile in inputPSFiles:
			outDir += pyRootPwa.config.phaseSpaceAmpDirectoryName
		elif inputFile in inputAccPSFiles:
			outDir += pyRootPwa.config.accCorrPSAmpDirectoryName
		if not os.path.exists(outDir):
			os.mkdir(outDir)
		else:
			if not os.path.isdir(outDir):
				pyRootPwa.utils.printErr('Output data directory does not appear to be a directory. Aborting...')
				sys.exit(1)

		for keyfile in keyfiles:

			if pyRootPwa.config.outputFileFormat == 'root':
				outFileExtension = '.root'
			else:
				outFileExtension = '.amp'

			outFileName = outDir + '/' + keyfile.rsplit('/',1)[-1].replace('.key', outFileExtension)

			if os.path.exists(outFileName):
				pyRootPwa.utils.printInfo('Output file "' + outFileName + '" already present. Skipping...')
				continue

			pyRootPwa.utils.printInfo('Opening input file "' + inputFile + '".')

			success = False
			if outFileExtension == '.amp':
				with open(outFileName, 'w') as outputFile:
					inFile = pyRootPwa.ROOT.TFile.Open(inputFile)
					success = pyRootPwa.amplitude.calcAmplitudes(inFile, keyfile, outputFile)
			else:
				outFile = pyRootPwa.ROOT.TFile.Open(outFileName, 'RECREATE')
				inFile = pyRootPwa.ROOT.TFile.Open(inputFile)
				success = pyRootPwa.amplitude.calcAmplitudes(inFile, keyfile, outFile)

			if success:
				pyRootPwa.utils.printSucc('Created amplitude file "' + outFileName + '".')
			else:
				pyRootPwa.utils.printErr('Amplitude calculation failed for input file "' + inputFile + '" and keyfile "' + keyfile + '".')
				if os.path.exists(outFileName):
					os.remove(outFileName)
			print
