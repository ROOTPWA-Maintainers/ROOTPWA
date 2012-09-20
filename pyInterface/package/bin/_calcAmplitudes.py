
import ConfigParser
import glob
import os
import sys
import traceback

import pyRootPwa
import pyRootPwa.amplitude
import pyRootPwa.configuration
import pyRootPwa.infile
import pyRootPwa.utils

def calcAmplitudes(configFileName, massBins):

	# print some info
	pyRootPwa.printCompilerInfo()
	pyRootPwa.printLibraryInfo()
	pyRootPwa.printSvnVersion()

	# config file
	pyRootPwa.config = pyRootPwa.configuration.rootPwaConfig(configFileName)

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
				pyRootPwa.utils.printErr('Output data directory "' + outDir + '" does not appear to be a directory. Skipping...')
				continue

		try:
			inFile = pyRootPwa.infile.inputFile(inputFile)
		except pyRootPwa.exception.pyRootPwaException as exc:
			pyRootPwa.utils.printErr('Could not open input file "' + inputFile + '": ' + str(exc) + '. Skipping...')
			continue
		except KeyboardInterrupt:
			pyRootPwa.utils.printInfo('Recieved keyboard interrupt. Aborting...')
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
			try:
				if outFileExtension == '.amp':
					with open(outFileName, 'w') as outFile:
						success = pyRootPwa.amplitude.calcAmplitudes(inFile, keyfile, outFile)
				else:
					outFile = pyRootPwa.ROOT.TFile.Open(outFileName, 'RECREATE')
					success = pyRootPwa.amplitude.calcAmplitudes(inFile, keyfile, outFile)
			except KeyboardInterrupt:
				pyRootPwa.utils.printInfo('Recieved keyboard interrupt. Aborting...')
				if os.path.exists(outFileName):
					os.remove(outFileName)
				sys.exit(1)
			except:
				pyRootPwa.utils.printErr('Unknown exception during amplitude calculation. Aborting...')
				traceback.print_exc()
				if os.path.exists(outFileName):
					os.remove(outFileName)
				sys.exit(1)

			if success:
				pyRootPwa.utils.printSucc('Created amplitude file "' + outFileName + '".')
			else:
				pyRootPwa.utils.printErr('Amplitude calculation failed for input file "' + inputFile + '" and keyfile "' + keyfile + '".')
				if os.path.exists(outFileName):
					os.remove(outFileName)

			print

