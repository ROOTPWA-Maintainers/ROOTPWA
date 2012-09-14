#!/usr/bin/python2.7

import argparse
import ConfigParser
import glob
import multiprocessing
import os
import Queue
import sys
import traceback

import pyRootPwa
import pyRootPwa.amplitude
import pyRootPwa.configuration
import pyRootPwa.infile
import pyRootPwa.utils

if __name__ == "__main__":

	# print some info
	pyRootPwa.printCompilerInfo()
	pyRootPwa.printLibraryInfo()
	pyRootPwa.printSvnVersion()

	# making the output nice for multi-threading
	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	printingCounter = multiprocessing.Array('i', [0]*5)

	pyRootPwa.utils.printErr = pyRootPwa.utils.printErrClass(printingCounter)
	pyRootPwa.utils.printWarn = pyRootPwa.utils.printWarnClass(printingCounter)
	pyRootPwa.utils.printSucc = pyRootPwa.utils.printSuccClass(printingCounter)
	pyRootPwa.utils.printInfo = pyRootPwa.utils.printInfoClass(printingCounter)
	pyRootPwa.utils.printDebug = pyRootPwa.utils.printDebugClass(printingCounter)

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
	parser.add_argument("-j", type=int, metavar=("jobs"), default=1, dest="nJobs", help="EXPERIMENTAL: number of jobs used (default: 1)")
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

	inputDataDict = {}

	for inputFile in (inputDataFiles + inputPSFiles + inputAccPSFiles):

		try:
			inFile = pyRootPwa.infile.inputFile(inputFile)
		except pyRootPwa.exception.pyRootPwaException as exc:
			pyRootPwa.utils.printErr('Could not open input file "' + inputFile + '": ' + str(exc) + '. Skipping...')
			continue
		except KeyboardInterrupt:
			pyRootPwa.utils.printInfo('Recieved keyboard interrupt. Aborting...')
			sys.exit(1)

		inputDataDict[inFile] = []

		for keyfile in keyfiles:

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
					pyRootPwa.utils.printErr('Output data directory "' + outDir + '" does not appear to be a directory. Skipping...')
					continue

			if pyRootPwa.config.outputFileFormat == 'root':
				outFileExtension = '.root'
			else:
				outFileExtension = '.amp'

			outFileName = outDir + '/' + keyfile.rsplit('/',1)[-1].replace('.key', outFileExtension)

			if os.path.exists(outFileName):
				pyRootPwa.utils.printInfo('Output file "' + outFileName + '" already present. Skipping...')
				continue

			inputDataDict[inFile].append((keyfile, outFileName))

	processQueue = multiprocessing.JoinableQueue()

	# Fill the queue, reorder so different input files come after each other
	while inputDataDict:
		for inFl in inputDataDict.keys():
			val = inputDataDict[inFl]
			if not val:
				del inputDataDict[inFl]
				continue
			(keyfile, outFileName) = val.pop(0)
			processQueue.put((inFl, keyfile, outFileName))

	jobs = []
	for i in range(arguments.nJobs):
		jobs.append(pyRootPwa.amplitude.AmplitudeCalculator(processQueue))
	for job in jobs:
		job.daemon = True
		job.start()

	try:
		processQueue.join()
		processQueue.close()
	except KeyboardInterrupt:
		pyRootPwa.utils.printInfo('Recieved keyboard interrupt. Aborting...')
		while True:
			try:
				processQueue.get(True, 2)
			except IOError:
				pass
			except Queue.Empty:
				break
	except:
		pyRootPwa.utils.printErr('Unknown exception during amplitude calculation. Aborting...')
		traceback.print_exc()
		while True:
			try:
				processQueue.get(True, 2)
			except IOError:
				pass
			except Queue.Empty:
				break

	pyRootPwa.utils.printPrintingSummary(printingCounter)

