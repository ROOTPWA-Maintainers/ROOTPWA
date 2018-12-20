import collections
import glob
import os
import itertools
import cPickle as pickle

import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT


def saveFileManager(fileManagerObject, path):
	if not os.path.isfile(path):
		try:
			pickle.dump(fileManagerObject, open(path, "wb"))
		except:
			pyRootPwa.utils.printErr("error saving file manager.")
			return False
		return True

	pyRootPwa.utils.printErr("cannot open file manager file '" + path + "'. File already exists.")
	return False


def loadFileManager(path):
	try:
		fileManagerObject = pickle.load(open(path, "rb"))
	except IOError:
		pyRootPwa.utils.printErr("error loading file manager. File manager not found at '" + path + "'")
		return None
	except:
		pyRootPwa.utils.printErr("error loading file manager.")
		raise
	if not fileManagerObject.areKeyFilesSynced():
		pyRootPwa.utils.printErr("key files are not the same as in file manager.")
		return None
	if not fileManagerObject.areEventFilesSynced():
		pyRootPwa.utils.printErr("event files are not the same as in file manager.")
		return None
	return fileManagerObject


class InputFile(object):

	def __init__(self, dataFileName, multibinBoundaries, eventsType, additionalVariables):
		self.dataFileName = dataFileName
		self.multibinBoundaries = multibinBoundaries
		self.eventsType = eventsType
		self.additionalVariables = additionalVariables


	def __str__(self):
		retval = ""
		retval += "'" + self.dataFileName + "': ("
		retval += str(self.multibinBoundaries) + ", "
		retval += str(self.eventsType) + ", "
		retval += str(self.additionalVariables) + ")"
		return retval


class EventsType(object):
# pylint: disable=C0103
	OTHER     = 0
	REAL      = 1
	GENERATED = 2
	ACCEPTED  = 3
# pylint: enable=C0103


class fileManager(object):


	def __init__(self):
		self.dataDirectory      = ""
		self.keyDirectory       = ""
		self.amplitudeDirectory = ""
		self.integralDirectory  = ""

		self.dataFiles              = collections.OrderedDict()
		self.keyFiles               = collections.OrderedDict()
		self._specialAmplitudeFiles = collections.OrderedDict()
		self._mergedAmplitudes      = collections.defaultdict(dict) # {<events-type>: {<tag>: [<wave-name> ...]}}
		self.integralFiles          = collections.OrderedDict()
		self.binList = []


	def initialize(self, configObject):
		self.dataDirectory      = configObject.dataDirectory
		self.keyDirectory       = configObject.keyDirectory
		self.amplitudeDirectory = configObject.ampDirectory
		self.integralDirectory  = configObject.intDirectory
		pyRootPwa.utils.printInfo("data file dir read from config file: '" + self.dataDirectory + "'.")
		pyRootPwa.utils.printInfo("key file dir read from config file: '" + self.keyDirectory + "'.")
		pyRootPwa.utils.printInfo("amplitude file dir read from config file: '" + self.amplitudeDirectory + "'.")
		pyRootPwa.utils.printInfo("integral file dir read from config file: '" + self.integralDirectory + "'.")

		self.binList = sorted(configObject.integralBinning)
		if not self.binList:
			pyRootPwa.utils.printWarn("no bins found, falling back to file binning")
		pyRootPwa.utils.printInfo("created a list with " + str(len(self.binList)) + " bins from config file.")
		self.dataFiles = self._openDataFiles()
		if not self.dataFiles:
			pyRootPwa.utils.printErr("no event files found.")
			return False
		if not self.binList:
			for _, inputFiles in self.dataFiles.iteritems():
				for inputFile in inputFiles:
					pyRootPwa.utils.printInfo("checking bin '" + str(inputFile.multibinBoundaries) + "'.")
					try:
						latestBin = pyRootPwa.utils.multiBin(inputFile.multibinBoundaries)
					except (TypeError, ValueError):
						pyRootPwa.utils.printErr("no binning given in config file and no multibin boundaries" +
						                         " found in data file '" + str(inputFile.dataFileName) + "'.")
						return False
					alreadyPresent = False
					for multiBin in self.binList:
						if latestBin == multiBin:
							alreadyPresent = True
							continue
						elif latestBin.overlap(multiBin, False):
							pyRootPwa.utils.printWarn("overlap found of bin '" + str(latestBin) + "' and bin '" + str(multiBin) + "'.")
					if not alreadyPresent:
						pyRootPwa.utils.printInfo("adding bin '" + str(inputFile.multibinBoundaries) + "'.")
						self.binList.append(latestBin)
		else:
			for _, inputFiles in self.dataFiles.iteritems():
				for inputFile in inputFiles:
					for variableName in self.binList[0].boundaries.keys():
						if variableName not in inputFile.additionalVariables:
							pyRootPwa.utils.printErr("variable '" + str(variableName) + "' required by binning, but not " +
							                         "in additional variables of event file '" + inputFile.dataFileName + "' (found " + str(inputFile.additionalVariables) + ").")
							return False
		self.keyFiles = self._openKeyFiles()
		if not self.keyFiles:
			pyRootPwa.utils.printErr("error loading keyfiles.")
			return False

		amplitudeFiles = self._getAmplitudeFilePaths(self.dataFiles.keys(), self.getWaveNameList())
		# make sure directories exist (create if necessary) and check if they are empty
		if not os.path.exists(self.amplitudeDirectory):
			os.makedirs(self.amplitudeDirectory)
		if not os.listdir(self.amplitudeDirectory) == []:
			pyRootPwa.utils.printWarn("directory '" + self.amplitudeDirectory + "' is not empty.")
		pyRootPwa.utils.printInfo("number of amplitude files: " + str(len(amplitudeFiles)))
		self.integralFiles = self._getIntegralFilePaths()
		pyRootPwa.utils.printInfo("number of integral files: " + str(len(self.integralFiles)))
		return True


	def getEventAndAmplitudePairPathsForWave(self, eventsType, waveName):
		eventsType = fileManager.pyEventsType(eventsType)
		retval = []
		if eventsType not in self.dataFiles:
			pyRootPwa.utils.printWarn("events type '" + str(eventsType) + "' not found.")
			return []
		if waveName not in self.keyFiles:
			pyRootPwa.utils.printWarn("key file for wave name '" + waveName + "' not found.")
			return []
		waveNameIndex = self.keyFiles.keys().index(waveName)
		eventFiles = self.dataFiles[eventsType]
		amplitudeFiles = self._getAmplitudeFilePaths([eventsType], [waveName])
		for eventFileId, eventFile in enumerate(eventFiles):
			amplitudeFile = amplitudeFiles[ (eventsType, eventFileId, waveNameIndex) ]
			retval.append( (eventFile.dataFileName, amplitudeFile) )
		return retval


	def getKeyFile(self, waveName):
		return self.keyFiles[waveName]


	def getWaveDescription(self, waveName):
		if waveName not in self.keyFiles:
			pyRootPwa.utils.printErr("wave name '" + str(waveName) + "' cannot be produced by any of the keyfiles.")
			return None
		keyFileName = self.keyFiles[waveName]
		waveDescriptions = pyRootPwa.core.waveDescription.parseKeyFile(keyFileName)
		for waveDescriptionID, waveDescription in enumerate(waveDescriptions):
			(success, amplitude) = waveDescription.constructAmplitude()
			if not success:
				pyRootPwa.utils.printErr("could not construct decay topology for wave description at index " + str(waveDescriptionID) + " of key file '" + keyFileName + "'.")
				return None
			constructedWaveName = waveDescription.waveNameFromTopology(amplitude.decayTopology())
			if waveName == constructedWaveName:
				return waveDescription
		pyRootPwa.utils.printErr("none of the constructed topologies matched the given wave name.")
		return None


	def getWaveDescriptions(self):
		retval = collections.OrderedDict()
		for waveName in self.keyFiles:
			retval[waveName] = self.getWaveDescription(waveName)
		pyRootPwa.utils.printSucc("constructed all {0} topologies for the known wave names.".format(len(retval)))
		return retval


	def getEventAndAmplitudeFilePathsInBin(self, multiBin, eventsType):
		# returns { "dataFileName" : { "waveName" : "amplitudeFileName" } }
		eventsType = fileManager.pyEventsType(eventsType)
		eventFileIds = self._getEventFileIdsForIntegralBin(multiBin, eventsType)
		if not eventFileIds:
			pyRootPwa.utils.printWarn("no matching event files found in bin '" + str(multiBin) + "' and events type '" + str(eventsType) + "'.")
			return collections.OrderedDict()
		retval = collections.OrderedDict()
		amplitudeFiles = self._getAmplitudeFilePaths([eventsType], self.getWaveNameList(), {eventsType: eventFileIds})
		for eventFileId in eventFileIds:
			eventFileName = self.dataFiles[eventsType][eventFileId].dataFileName
			retval[eventFileName] = collections.OrderedDict()
			for waveName_i, waveName in enumerate(self.getWaveNameList()):
				retval[eventFileName][waveName] = amplitudeFiles[(eventsType, eventFileId, waveName_i)]
		return retval


	def getIntegralFilePath(self, multiBin, eventsType):
		eventsType = fileManager.pyEventsType(eventsType)
		return self.integralFiles[(self.binList.index(multiBin), eventsType)]


	def _getEventFileIdsForIntegralBin(self, multiBin, eventsType):
		eventsType = fileManager.pyEventsType(eventsType)
		if eventsType not in self.dataFiles:
			pyRootPwa.utils.printWarn("events type '" + str(eventsType) + "' not in data files.")
			return []
		eventFileIds = []
		for eventFileId, inputFile in enumerate(self.dataFiles[eventsType]):
			found = True
			for variableName in inputFile.multibinBoundaries:
				if variableName in multiBin.boundaries:
					if (inputFile.multibinBoundaries[variableName][1] < multiBin.boundaries[variableName][0]) or \
					   (inputFile.multibinBoundaries[variableName][0] > multiBin.boundaries[variableName][1]):
						found = False
						break
			if found:
				eventFileIds.append(eventFileId)
		return eventFileIds


	def getWaveNameList(self):
		return self.keyFiles.keys()

	def getWaveIndex(self, waveName):
		return self.getWaveNameList().index(waveName)


	def areEventFilesSynced(self):
		return set(self._getEventFilePaths()) == set(glob.glob(self.dataDirectory + "/*.root"))


	def areKeyFilesSynced(self):
		return set(self._getKeyFilePaths()) == set(glob.glob(self.keyDirectory + "/*.key"))


	def areFilesSynced(self):
		return self.areEventFilesSynced() and self.areKeyFilesSynced()


	def _nmbAmplitudeFiles(self):
		nmbEventFiles = 0
		for key in self.dataFiles:
			nmbEventFiles += len(self.dataFiles[key])
		return nmbEventFiles * len(self.keyFiles)


	def _getAmplitudeFilePaths(self, eventsTypes, waveNames, eventFileIds = None):
		'''
		Build ordered dictionary of AmplitudeFilePaths for the given eventsTypes and the given waves
		@param eventFileIds: Dictionary of {<eventsType>: [<eventFieldID>]. If None, all event-field IDs of the given eventsTypes.
		'''
		amplitudeFiles = collections.OrderedDict()
		if not self.dataFiles:
			pyRootPwa.utils.printWarn("cannot create amplitude file path collection without data files.")
			return collections.OrderedDict()
		waves = [(self.getWaveIndex(w), w) for w in waveNames]

		for eventsType in eventsTypes:
			if eventFileIds is None:
				inputFileIndices = range(len(self.dataFiles[eventsType]))
			else:
				if eventsType not in eventFileIds:
					pyRootPwa.utils.printErr("Events type '{0}' not in dictionary of event-field IDs!".format(eventsType))
					raise Exception()
				inputFileIndices = eventFileIds[eventsType]

			if eventsType in self._mergedAmplitudes:
				wavesWithMergedAmplitudes = list(itertools.chain(*self._mergedAmplitudes[eventsType].values()))
			else:
				wavesWithMergedAmplitudes = []

			for inputFile_i in inputFileIndices:
				for waveName_i, waveName in waves:
					key = (eventsType, inputFile_i, waveName_i)
					if key not in self._specialAmplitudeFiles:
						if waveName not in wavesWithMergedAmplitudes:
							amplitudeFileName = "{waveName}_eventFileId-{eventsId}_{eventsType}.root".format(waveName=waveName,
							                                                                                 eventsId=str(inputFile_i),
							                                                                                 eventsType=str(eventsType))
						else:
							tags = [tag for tag, waveNames in self._mergedAmplitudes[eventsType].iteritems() if waveName in waveNames]
							if len(tags) != 1:
								pyRootPwa.utils.printErr("Cannot find unique merge tag for wave '{0}'. Aborting...".format(waveName))
								raise Exception()
							amplitudeFileName = "{tag}_eventFileId-{eventsId}_{eventsType}.root".format(tag=tags[0],
							                                                                            eventsId=str(inputFile_i),
							                                                                            eventsType=str(eventsType))
						amplitudeFilePath = os.path.join(self.amplitudeDirectory, amplitudeFileName)
					else:
						amplitudeFilePath = self._specialAmplitudeFiles[key]
					amplitudeFiles[key] = amplitudeFilePath

		return amplitudeFiles


	def _getIntegralFilePaths(self):
		integralFiles = collections.OrderedDict()
		for binID, _ in enumerate(self.binList):
			for eventsType in [EventsType.GENERATED, EventsType.ACCEPTED]:
				integralFiles[(binID, eventsType)] = self.integralDirectory + "/integral_binID-" + str(binID) + "_" + str(eventsType) + ".root"
		return integralFiles


	def _openDataFiles(self):
		dataFileNames = sorted(glob.glob(self.dataDirectory + "/*.root"))
		inputFiles = {}
		for dataFileName in dataFileNames:
			dataFile = ROOT.TFile.Open(dataFileName, "READ")
			if not dataFile:
				pyRootPwa.utils.printErr("could not open event file '" + dataFileName + "'.")
				return collections.OrderedDict()
			eventMeta = pyRootPwa.core.eventMetadata.readEventFile(dataFile)
			if not eventMeta:
				pyRootPwa.utils.printErr("could not find metadata in event file '" + dataFileName + "'.")
				return  collections.OrderedDict()
			inputFile = InputFile(dataFileName,
			                      eventMeta.multibinBoundaries(),
			                      fileManager.pyEventsType(eventMeta.eventsType()),
			                      eventMeta.additionalTreeVariableNames())
			if inputFile.eventsType not in inputFiles:
				inputFiles[inputFile.eventsType] = []
			inputFiles[inputFile.eventsType].append(inputFile)
			dataFile.Close()
		retval =  collections.OrderedDict()
		for eventsType in sorted(inputFiles):
			retval[eventsType] = inputFiles[eventsType]
		return inputFiles


	def _openKeyFiles(self):
		keyFileNames = glob.glob(self.keyDirectory + "/*.key")
		keyFiles = {}

		for keyFileName in keyFileNames:
			waveDescriptions = pyRootPwa.core.waveDescription.parseKeyFile(keyFileName)
			if not waveDescriptions:
				pyRootPwa.utils.printErr("could not read wave description from key file '" + keyFileName + "'.")
				return collections.OrderedDict()
			for waveDescriptionID, waveDescription in enumerate(waveDescriptions):
				(success, amplitude) = waveDescription.constructAmplitude()
				if not success:
					pyRootPwa.utils.printErr("could not construct decay topology for wave description at index " + str(waveDescriptionID) + " of key file '" + keyFileName + "'.")
					return collections.OrderedDict()
				waveName = waveDescription.waveNameFromTopology(amplitude.decayTopology())
				if waveName in keyFiles.keys():
					pyRootPwa.utils.printErr("duplicate wave name ('" + waveName + "' from files '" + keyFiles[waveName] +
					                         "' and '" + keyFileName + "' (index " + str(waveDescriptionID) + ").")
					return collections.OrderedDict()
				keyFiles[waveName] = keyFileName
		retval = collections.OrderedDict()
		for waveName in sorted(keyFiles):
			retval[waveName] = keyFiles[waveName]
		return retval


	def _getEventFilePaths(self):
		retval = []
		for _, inputFiles in self.dataFiles.iteritems():
			for inputFile in inputFiles:
				retval.append(inputFile.dataFileName)
		return retval


	def _getKeyFilePaths(self):
		retval = []
		for _, keyfilePath in self.keyFiles.iteritems():
			if keyfilePath not in retval:
				retval.append(keyfilePath)
		return retval


	def __repr__(self):
		retStr = "keyfiles:\n"
		for waveName in self.keyFiles.keys():
			retStr += waveName + " >> " + self.keyFiles[waveName] + "\n"
		retStr += "\nDataFiles:\n"
		for eventsType in self.dataFiles:
			for dataFile in self.dataFiles[eventsType]:
				retStr += ("eventsType [" + str(eventsType) + "], bin [" + str(dataFile.multibinBoundaries) +
				           "] >> " + dataFile.dataFileName + "\n")
		retStr += "\nAmpFiles:\n"
		for eventsType in self.dataFiles:
			for eventFileId, _ in enumerate(self.dataFiles[eventsType]):
				amplitudeFiles = self._getAmplitudeFilePaths([eventsType], self.getWaveNameList(), {eventsType: [eventFileId]})
				for waveName_i, waveName in enumerate(self.keyFiles.keys()):
					retStr += ("eventsType [" + str(eventsType) + "], eventFileId [" + str(eventFileId) +
					           "], waveName [" + waveName + "] >> " + amplitudeFiles[(eventsType, eventFileId, waveName_i)] + "\n")
		retStr += "\nIntFiles:\n"
		for eventsType in [EventsType.GENERATED, EventsType.ACCEPTED]:
			for binID, _ in enumerate(self.binList):
				retStr += "eventsType [" + str(eventsType) + "], binID [" + str(binID) + "] >> " + self.integralFiles[(binID, eventsType)] + "\n"
		return retStr


	@staticmethod
	def pyEventsType(eventsType):
		if not isinstance(eventsType, pyRootPwa.core.eventMetadata.eventsTypeEnum):
			return eventsType
		if eventsType == pyRootPwa.core.eventMetadata.OTHER:
			return EventsType.OTHER
		elif eventsType == pyRootPwa.core.eventMetadata.REAL:
			return EventsType.REAL
		elif eventsType == pyRootPwa.core.eventMetadata.GENERATED:
			return EventsType.GENERATED
		elif eventsType == pyRootPwa.core.eventMetadata.ACCEPTED:
			return EventsType.ACCEPTED
		return -1
