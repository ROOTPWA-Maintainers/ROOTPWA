
import os as _os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def calcAmplitude(inputFileName,
                  keyFileName,
                  outputFileName,
                  maxNumberOfEvents = -1,
                  printProgress = True):

	printDebug = pyRootPwa.utils.printDebug
	printInfo = pyRootPwa.utils.printInfo
	printSucc = pyRootPwa.utils.printSucc
	printWarn = pyRootPwa.utils.printWarn
	printErr = pyRootPwa.utils.printErr

	printInfo("Calculating amplitude for key file '" + keyFileName +
	                          "' with input file '" + inputFileName + "', output file '" + outputFileName + ".")

	if 'ROOTPWA' not in _os.environ:
		printWarn("$ROOTPWA not set.")
		return False

	inputFile = ROOT.TFile.Open(inputFileName, "READ")
	if not inputFile:
		printWarn("could not open input file '" + inputFileName + "'.")
		return False
	eventMeta = pyRootPwa.core.eventMetadata.readEventFile(inputFile, True)
	if not eventMeta:
		printWarn("could not read metadata from input file '" + inputFileName + "'.")
		return False
	nEvents = eventMeta.eventTree().GetEntries()
	if maxNumberOfEvents > 0:
		nEvents = min(nEvents, maxNumberOfEvents)
	outputFile = ROOT.TFile.Open(outputFileName, "NEW")
	if not outputFile:
		printWarn("could not open output file '" + outputFileName + "'.")
		return False

	if not _os.path.isfile(keyFileName):
		printWarn("key file '" + keyFileName + "' not found.")
		outputFile.Close()
		return False
	waveDescription = pyRootPwa.core.waveDescription()
	waveDescription.parseKeyFile(keyFileName)
	(result, amplitude) = waveDescription.constructAmplitude()
	if not result:
		printWarn("could not construct amplitude from keyfile '" + keyFileName + "'.")
		outputFile.Close()
		return False

	ampFileWriter = pyRootPwa.core.amplitudeFileWriter()
	objectBaseName = waveDescription.waveNameFromTopology(amplitude.decayTopology(), True)
	if not ampFileWriter.initialize(outputFile, [eventMeta], waveDescription.keyFileContent(), objectBaseName):
		printWarn("could not initialize amplitudeFileWriter.")
		outputFile.Close()
		return False
	amplitudes = pyRootPwa.core.calcAmplitude(eventMeta, amplitude, nEvents, printProgress)
	if not amplitudes:
		printWarn("could not calculate amplitudes.")
		outputFile.Close()
		return False
	if nEvents != len(amplitudes):
		printWarn("number of events (" + str(nEvents) +
		          ") does not match with number of amplitudes (" + str(len(amplitudes)) + ").")
		return False
	ampFileWriter.addAmplitudes(amplitudes)
	if not ampFileWriter.finalize():
		printWarn("could not finalize amplitudeFileWriter.")
		outputFile.Close()
		return False

	outputFile.Close()
	printSucc("successfully calculated amplitude for " + str(nEvents) + " events.")
	return True
