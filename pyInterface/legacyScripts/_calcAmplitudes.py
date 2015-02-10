
import glob
import multiprocessing
import os
import Queue
import sys
import traceback

import pyRootPwa.amplitude
import pyRootPwa.core
import pyRootPwa.utils

def calcAmplitudes(configFileName, massBins, **arguments):

	# print some info
	pyRootPwa.core.printCompilerInfo()
	pyRootPwa.core.printLibraryInfo()
	pyRootPwa.core.printGitHash()

	# making the output nice for multi-threading
	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()

	# Check and get the arguments
	nJobs = 1
	if 'nJobs' in arguments:
		nJobs = int(arguments['nJobs'])
		if nJobs < 1:
			nJobs = 1
	progressBar = True
	if 'progressBar' in arguments:
		progressBar = bool(arguments['progressBar'])
	proFile = ''
	if 'proFile' in arguments:
		if nJobs != 1:
			pyRootPwa.utils.printErr('Profiler cannot be run with nJobs > 1. Aborting...')
			return
		proFile = str(arguments['proFile'])
	maxNmbEvents = -1
	if 'maxNmbEvents' in arguments:
		maxNmbEvents = int(arguments['maxNmbEvents'])

	# config file
	pyRootPwa.config = pyRootPwa.rootPwaConfig(configFileName)

	# get the massBins as specified on the command line
	allMassBins = sorted(glob.glob(pyRootPwa.config.dataDirectory + '/' + pyRootPwa.config.massBinDirectoryNamePattern))
	massBins = pyRootPwa.utils.parseMassBinArgs(allMassBins, massBins)
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
	if "keyfiles" in arguments:
		keyfiles = arguments["keyfiles"]
		for keyfile in keyfiles:
			if not (os.path.isfile(keyfile) and keyfile.endswith(".key")):
				pyRootPwa.utils.printErr("Keyfile '" + keyfile + "' not valid. Aborting...")
				sys.exit(1)
	else:
		keyfiles = pyRootPwa.utils.getListOfKeyfiles(pyRootPwa.config.keyfilePattern)
		if len(keyfiles) == 0:
			pyRootPwa.utils.printErr("No keyfiles found with valid file extension. Aborting...")
			sys.exit(1)

	# initialize the particleDataTable
	pyRootPwa.core.particleDataTable.readFile(pyRootPwa.config.pdgFileName)

	inputDataDict = {}

	for inputFile in (inputDataFiles + inputPSFiles + inputAccPSFiles):

		try:
			inFile = pyRootPwa.amplitude.inputFile(inputFile)
		except pyRootPwa.rootPwaException as exc:
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
	silence = (nJobs != 1)
	for i in range(nJobs):
		jobs.append(pyRootPwa.amplitude.amplitudeCalculator(processQueue, silence, progressBar, maxNmbEvents))
	for job in jobs:
		job.daemon = True
		if proFile != '':
			job.run(True)
		else:
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

	pyRootPwa.utils.printPrintingSummary(pyRootPwa.utils.printingCounter)
