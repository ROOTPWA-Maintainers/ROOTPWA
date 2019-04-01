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

	pyRootPwa.utils.printErr("cannot save file manager file '" + path + "'. File already exists.")
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

	def __init__(self, dataFileName, datasetLabel, multibinBoundaries, eventsType, additionalVariables):
		self.dataFileName = dataFileName
		self.datasetLabel = datasetLabel
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
	'''
	@brief: This class handles all information about file paths, file-name conventions, defined waves, and other administrative information

	EventFiles: Event files contain trees with kinematic variables, e.g. momenta, and binning variables, e.g. mass and t'.
	            Events of different events types (real data, generated, accepted) are stored in individual event files.
	            Event files can contain events only from one kinematic bin. The binning of the event files must not be mixed up
	            with the binning of the final fit (see `Bins` below).
	            During the initialization of the file manager, an unique id, the eventFileId, is assigned to each event file.

	KeyFiles: KeyFiles represent key files in which partial waves, e.g. amplitudes are defined

	Amplitudes: Amplitudes represent amplitude files in which precalculated amplitudes for the input events are stored.
	            By default, the amplitude file names follow the naming scheme '<waveName>_eventFileId-<eventFileId>_<eventsType>.root'.
	            In addition, special amplitude file names can be defined for each amplitude, event-field id,
	            and amplitude name individually (see member variable `_specialAmplitudeFiles`)
	            In addition, one amplitude file can contain more amplitudes of the same event file. The member
	            variable `_mergedAmplitudes` can contain lists of wave names for each events type. The amplitudes
	            of all waves in this list are merged in one amplitude file per event file, for all event files.
	            The tag identifies the merged amplitude file name containing the amplitudes.

	Integrals: Integrals represent integral files in which the precalculated integral matrices for the generated and
	           accepted Monte Carlo events are stored.
	           One integral per bin and events type (generated or accepted)

	Bins: `Bins` are the bins in which the fit is performed. Also the integral matrices are calculated in these bins.
	      The event files have their own binning, which is independent from the binning represented by `bins`.

	'''

	_deprecatedMembers = {'dataDirectory': 'eventDirectory', 'dataFiles': 'eventFiles'}

	def __init__(self):
		self.eventDirectory     = ""
		self.keyDirectory       = ""
		self.amplitudeDirectory = ""
		self.integralDirectory  = ""

		self.eventFiles             = collections.OrderedDict() # {<events-type>: [<event-filepath>]}
		self.keyFiles               = collections.OrderedDict() # {<wavename>: <key-filepath>}
		self._specialAmplitudeFiles = collections.OrderedDict() # {(<events-type>, <eventfile-id>, <wavename>): <amplitude-filepath>}
		self._mergedAmplitudes      = collections.defaultdict(dict) # {<events-type>: {<tag>: [<wave-name> ...]}}
		self.integralFiles          = collections.OrderedDict() # {(<dataset-id>, <bin-id>, <events-type>): <integral-file path>}
		self.binList = []
		self.datasetLabels = []

	def __getattribute__(self, name, *args, **kwargs):
		if name in fileManager._deprecatedMembers:
			pyRootPwa.utils.printWarn("'{0}' member is deprecated. Use '{1}'!".format(name, fileManager._deprecatedMembers[name]))
			return object.__getattribute__(self, fileManager._deprecatedMembers[name])
		return object.__getattribute__(self, name, *args, **kwargs)

	def __setattr__(self, name, *args, **kwargs):
		if name in fileManager._deprecatedMembers:
			pyRootPwa.utils.printWarn("'{0}' member is deprecated. Use '{1}'!".format(name, fileManager._deprecatedMembers[name]))
			return object.__setattr__(self, fileManager._deprecatedMembers[name], *args, **kwargs)
		return object.__setattr__(self, name, *args, **kwargs)

	def initialize(self, configObject):
		self.eventDirectory      = configObject.eventDirectory
		self.keyDirectory       = configObject.keyDirectory
		self.amplitudeDirectory = configObject.ampDirectory
		self.integralDirectory  = configObject.intDirectory
		pyRootPwa.utils.printInfo("event file dir read from config file: '" + self.eventDirectory + "'.")
		pyRootPwa.utils.printInfo("key file dir read from config file: '" + self.keyDirectory + "'.")
		pyRootPwa.utils.printInfo("amplitude file dir read from config file: '" + self.amplitudeDirectory + "'.")
		pyRootPwa.utils.printInfo("integral file dir read from config file: '" + self.integralDirectory + "'.")

		self.binList = sorted(configObject.integralBinning)
		if not self.binList:
			pyRootPwa.utils.printWarn("no bins found, falling back to file binning")
		pyRootPwa.utils.printInfo("created a list with " + str(len(self.binList)) + " bins from config file.")
		self.eventFiles = self._openDataFiles()
		if not self.eventFiles:
			pyRootPwa.utils.printErr("no event files found.")
			return False
		if not self.datasetLabels:
			self.datasetLabels = set()
			for _, eventFiles in self.eventFiles.iteritems():
				for eventFile in eventFiles:
					self.datasetLabels.add(eventFile.datasetLabel)
			self.datasetLabels = sorted(list(self.datasetLabels))
		if not self.binList:
			for _, eventFiles in self.eventFiles.iteritems():
				for eventFile in eventFiles:
					pyRootPwa.utils.printInfo("checking bin '" + str(eventFile.multibinBoundaries) + "'.")
					try:
						latestBin = pyRootPwa.utils.multiBin(eventFile.multibinBoundaries)
					except (TypeError, ValueError):
						pyRootPwa.utils.printErr("no binning given in config file and no multibin boundaries" +
						                         " found in event file '" + str(eventFile.dataFileName) + "'.")
						return False
					alreadyPresent = False
					for multiBin in self.binList:
						if latestBin == multiBin:
							alreadyPresent = True
							continue
						elif latestBin.overlap(multiBin, False):
							pyRootPwa.utils.printWarn("overlap found of bin '" + str(latestBin) + "' and bin '" + str(multiBin) + "'.")
					if not alreadyPresent:
						pyRootPwa.utils.printInfo("adding bin '" + str(eventFile.multibinBoundaries) + "'.")
						self.binList.append(latestBin)
		else:
			for _, eventFiles in self.eventFiles.iteritems():
				for eventFile in eventFiles:
					for variableName in self.binList[0].boundaries.keys():
						if variableName not in eventFile.additionalVariables:
							pyRootPwa.utils.printErr("variable '" + str(variableName) + "' required by binning, but not " +
							                         "in additional variables of event file '" + eventFile.dataFileName + "' (found " + str(eventFile.additionalVariables) + ").")
							return False
		self.keyFiles = self._openKeyFiles()
		if not self.keyFiles:
			pyRootPwa.utils.printErr("error loading keyfiles.")
			return False

		amplitudeFiles = self._getAmplitudeFilePaths(self.eventFiles.keys(), self.getWaveNameList())
		# check if amplitude directory is empty
		if not os.listdir(self.amplitudeDirectory) == []:
			pyRootPwa.utils.printWarn("directory '" + self.amplitudeDirectory + "' is not empty.")
		pyRootPwa.utils.printInfo("number of amplitude files: " + str(len(amplitudeFiles)))
		self.integralFiles = self._getIntegralFilePaths()
		pyRootPwa.utils.printInfo("number of integral files: " + str(len(self.integralFiles)))
		return True

	def _datasetToDatasetID(self, dataset):
		'''
		@param dataset: data-set ID or data-set label
		@return datasetID
		'''
		if isinstance(dataset, int):
			if dataset >=0 and dataset < len(self.datasetLabels):
				dataset = dataset
			else:
				raise pyRootPwa.utils.exception("Data-set id '{0}' out of range!".format(dataset))
		elif isinstance(dataset, str):
			if dataset in self.datasetLabels:
				dataset = self.datasetLabels.index(dataset)
			else:
				raise pyRootPwa.utils.exception("Data-set label '{0}' not in list of data-set labels '{1}'!".format(dataset, self.datasetLabels))
		elif dataset is None: # for backwards compatibility
			if len(self.datasetLabels) == 1:
				dataset = 0
			else:
				raise pyRootPwa.utils.exception("More the one data-set available. Data-set needs to be defined!")
		else:
			raise pyRootPwa.utils.exception("Cannot handle data-set identifier type '{0}'!".format(type(dataset)))
		return dataset

	def getEventAndAmplitudePairPathsForWave(self, eventsType, waveName):
		'''
		@param eventsType: Return EventAndAmplitudePairPaths for this events type
		@return: List of (eventfilepath, amplitudefilepath) for the given wave, for all all event-field IDs
		'''
		eventsType = fileManager.pyEventsType(eventsType)
		retVal = []
		if eventsType not in self.eventFiles:
			pyRootPwa.utils.printWarn("events type '" + str(eventsType) + "' not found.")
			return []
		if waveName not in self.keyFiles:
			pyRootPwa.utils.printWarn("key file for wave name '" + waveName + "' not found.")
			return []
		waveNameIndex = self.keyFiles.keys().index(waveName)
		eventFiles = self.eventFiles[eventsType]
		amplitudeFiles = self._getAmplitudeFilePaths([eventsType], [waveName])
		for eventFileId, eventFile in enumerate(eventFiles):
			amplitudeFile = amplitudeFiles[ (eventsType, eventFileId, waveNameIndex) ]
			retVal.append( (eventFile.dataFileName, amplitudeFile) )
		return retVal


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


	def getEventAndAmplitudeFilePathsInBin(self, multiBin, eventsType, dataset = None):
		'''
		@return: { "eventFileName" : { "waveName" : "amplitudeFileName" } }
		'''
		eventsType = fileManager.pyEventsType(eventsType)
		datasetID = self._datasetToDatasetID(dataset)
		eventFileIds = self._getEventFileIdsForIntegralBin(multiBin, eventsType, datasetID)
		if not eventFileIds:
			pyRootPwa.utils.printWarn("no matching event files found in bin '" + str(multiBin) + "' and events type '" + str(eventsType) + "'.")
			return collections.OrderedDict()
		retval = collections.OrderedDict()
		amplitudeFiles = self._getAmplitudeFilePaths([eventsType], self.getWaveNameList(), {eventsType: eventFileIds})
		for eventFileId in eventFileIds:
			eventFileName = self.eventFiles[eventsType][eventFileId].dataFileName
			retval[eventFileName] = collections.OrderedDict()
			for waveName_i, waveName in enumerate(self.getWaveNameList()):
				retval[eventFileName][waveName] = amplitudeFiles[(eventsType, eventFileId, waveName_i)]
		return retval


	def getIntegralFilePath(self, multiBin, eventsType, dataset = None):
		datasetID = self._datasetToDatasetID(dataset)
		eventsType = fileManager.pyEventsType(eventsType)
		return self.integralFiles[(datasetID, self.binList.index(multiBin), eventsType)]


	def openIntegralFile(self, multiBin, eventsType, dataset = None):
		'''
		@return: integralMatrix, integralMetadata, integralFile
		'''
		integralFile = ROOT.TFile.Open(self.getIntegralFilePath(multiBin, eventsType, dataset), "READ")
		integralMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(integralFile)
		integralMatrix = integralMeta.getAmpIntegralMatrix()
		return integralMatrix, integralMeta, integralFile


	def getIntegralMatrix(self, multiBin, eventsType, dataset = None):
		'''
		@return: integralMatrix
		'''
		integralMatrix, _, integralFile = self.openIntegralFile(multiBin, eventsType, dataset)
		integralFile.Close()
		return integralMatrix


	def _getEventFileIdsForIntegralBin(self, multiBin, eventsType, datasetID):
		'''
		@return: List of all event-field IDs of the given events type that lie within or overlay with the given multibin.
		'''
		eventsType = fileManager.pyEventsType(eventsType)
		if eventsType not in self.eventFiles:
			pyRootPwa.utils.printWarn("events type '" + str(eventsType) + "' not in event files.")
			return []
		eventFileIds = []
		for eventFileId, eventFile in enumerate(self.eventFiles[eventsType]):
			overlaps = True
			if eventFile.datasetLabel != self.datasetLabels[datasetID]:
				overlaps = False
			for variableName in eventFile.multibinBoundaries:
				if variableName in multiBin.boundaries:
					# if event-file bin is completely below or above the multibin
					if (eventFile.multibinBoundaries[variableName][1] < multiBin.boundaries[variableName][0]) or \
					   (eventFile.multibinBoundaries[variableName][0] > multiBin.boundaries[variableName][1]):
						overlaps = False
						break
			if overlaps:
				eventFileIds.append(eventFileId)
		return eventFileIds


	def getWaveNameList(self):
		return self.keyFiles.keys()

	def getWaveIndex(self, waveName):
		return self.getWaveNameList().index(waveName)


	def areEventFilesSynced(self):
		return set(self._getEventFilePaths()) == set(glob.glob(self.eventDirectory + "/*.root"))


	def areKeyFilesSynced(self):
		return set(self._getKeyFilePaths()) == set(glob.glob(self.keyDirectory + "/*.key"))


	def areFilesSynced(self):
		return self.areEventFilesSynced() and self.areKeyFilesSynced()


	def _nmbAmplitudeFiles(self):
		nmbEventFiles = 0
		for key in self.eventFiles:
			nmbEventFiles += len(self.eventFiles[key])
		return nmbEventFiles * len(self.keyFiles)


	def _getAmplitudeFilePaths(self, eventsTypes, waveNames, eventFileIds = None):
		'''
		Build ordered dictionary of AmplitudeFilePaths for the given eventsTypes and the given waves
		@param eventFileIds: Dictionary of {<eventsType>: [<eventFieldID>]. If None, all event-field IDs of the given eventsTypes.
		'''
		amplitudeFiles = collections.OrderedDict()
		if not self.eventFiles:
			pyRootPwa.utils.printWarn("cannot create amplitude file path collection without event files.")
			return collections.OrderedDict()
		waves = [(self.getWaveIndex(w), w) for w in waveNames]

		for eventsType in eventsTypes:
			if eventFileIds is None:
				outputEventFileIds = range(len(self.eventFiles[eventsType]))
			else:
				if eventsType not in eventFileIds:
					pyRootPwa.utils.printErr("Events type '{0}' not in dictionary of event-field IDs!".format(eventsType))
					raise Exception()
				outputEventFileIds = eventFileIds[eventsType]

			if eventsType in self._mergedAmplitudes:
				wavesWithMergedAmplitudes = list(itertools.chain(*self._mergedAmplitudes[eventsType].values()))
			else:
				wavesWithMergedAmplitudes = []

			for eventFileId in outputEventFileIds:
				for waveName_i, waveName in waves:
					key = (eventsType, eventFileId, waveName_i)
					if key not in self._specialAmplitudeFiles:
						if waveName not in wavesWithMergedAmplitudes:
							amplitudeFileName = "{waveName}_eventFileId-{eventsId}_{eventsType}.root".format(waveName=waveName,
							                                                                                 eventsId=str(eventFileId),
							                                                                                 eventsType=str(eventsType))
						else:
							tags = [tag for tag, waveNames in self._mergedAmplitudes[eventsType].iteritems() if waveName in waveNames]
							if len(tags) != 1:
								pyRootPwa.utils.printErr("Cannot find unique merge tag for wave '{0}'. Aborting...".format(waveName))
								raise Exception()
							amplitudeFileName = "{tag}_eventFileId-{eventsId}_{eventsType}.root".format(tag=tags[0],
							                                                                            eventsId=str(eventFileId),
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
				for datasetID, _ in enumerate(self.datasetLabels):
					filename = "integral_dataset-"+str(datasetID)+"_binID-" + str(binID) + "_" + str(eventsType) + ".root"
					integralFiles[(datasetID, binID, eventsType)] = os.path.join(self.integralDirectory, filename)
		return integralFiles


	def _openDataFiles(self):
		dataFileNames = sorted(glob.glob(self.eventDirectory + "/*.root"))
		eventFiles = {}
		for dataFileName in dataFileNames:
			eventFile = ROOT.TFile.Open(dataFileName, "READ")
			if not eventFile:
				pyRootPwa.utils.printErr("could not open event file '" + dataFileName + "'.")
				return collections.OrderedDict()
			eventMeta = pyRootPwa.core.eventMetadata.readEventFile(eventFile)
			if not eventMeta:
				pyRootPwa.utils.printErr("could not find metadata in event file '" + dataFileName + "'.")
				return  collections.OrderedDict()
			inputFile = InputFile(dataFileName,
			                      eventMeta.datasetLabel(),
			                      eventMeta.multibinBoundaries(),
			                      fileManager.pyEventsType(eventMeta.eventsType()),
			                      eventMeta.additionalTreeVariableNames())
			if inputFile.eventsType not in eventFiles:
				eventFiles[inputFile.eventsType] = []
			eventFiles[inputFile.eventsType].append(inputFile)
			eventFile.Close()
		retval =  collections.OrderedDict()
		for eventsType in sorted(eventFiles):
			retval[eventsType] = eventFiles[eventsType]
		return eventFiles


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
		for _, eventFiles in self.eventFiles.iteritems():
			for eventFile in eventFiles:
				retval.append(eventFile.dataFileName)
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
		for eventsType in self.eventFiles:
			for eventFile in self.eventFiles[eventsType]:
				retStr += ("eventsType [" + str(eventsType) + "], bin [" + str(eventFile.multibinBoundaries) +
				           "] >> " + eventFile.dataFileName + "\n")
		retStr += "\nAmpFiles:\n"
		for eventsType in self.eventFiles:
			for eventFileId, _ in enumerate(self.eventFiles[eventsType]):
				amplitudeFiles = self._getAmplitudeFilePaths([eventsType], self.getWaveNameList(), {eventsType: [eventFileId]})
				for waveName_i, waveName in enumerate(self.keyFiles.keys()):
					retStr += ("eventsType [" + str(eventsType) + "], eventFileId [" + str(eventFileId) +
					           "], waveName [" + waveName + "] >> " + amplitudeFiles[(eventsType, eventFileId, waveName_i)] + "\n")
		retStr += "\nIntFiles:\n"
		for datasetID, _ in enumerate(self.datasetLabels):
			for eventsType in [EventsType.GENERATED, EventsType.ACCEPTED]:
				for binID, _ in enumerate(self.binList):
					retStr += "datasetID [" + str(datasetID) + "], eventsType [" + str(eventsType) + "], binID [" + str(binID) + "] >> "
					retStr += self.integralFiles[(datasetID, binID, eventsType)] + "\n"
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
