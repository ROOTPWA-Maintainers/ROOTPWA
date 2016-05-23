import collections
import glob
import os
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
	else:
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

	def __init__(self, dataFileName, binningMap, eventsType, additionalVariables):
		self.dataFileName = dataFileName
		self.binningMap = binningMap
		self.eventsType = eventsType
		self.additionalVariables = additionalVariables


	def __str__(self):
		retval = ""
		retval += "'" + self.dataFileName + "': ("
		retval += str(self.binningMap) + ", "
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
		self.limitFilesInDir = -1

		self.dataFiles       = collections.OrderedDict()
		self.keyFiles        = collections.OrderedDict()
		self.amplitudeFiles  = collections.OrderedDict()
		self.integralFiles   = collections.OrderedDict()
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
					self.binList.append(inputFile.binningMap)
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

		if int(configObject.limitFilesPerDir) <= 0:
			self.limitFilesInDir = -1
			pyRootPwa.utils.printInfo("limit for files per directory set to infinite")
		else:
			self.limitFilesInDir = int(configObject.limitFilesPerDir)
			pyRootPwa.utils.printInfo("limit for files per directory set to " + str(self.limitFilesInDir))
		self.amplitudeFiles = self._getAmplitudeFilePaths()
		pyRootPwa.utils.printInfo("number of amplitude files: " + str(len(self.amplitudeFiles)))
		self.integralFiles = self._getIntegralFilePaths()
		pyRootPwa.utils.printInfo("number of integral files: " + str(len(self.integralFiles)))
		return True


	def checkBinning(self):
		raise Exception("not implemented. " + repr(self))


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
		for eventFileId, eventFile in enumerate(eventFiles):
			amplitudeFile = self.amplitudeFiles[ (eventsType, eventFileId, waveNameIndex) ]
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
				pyRootPwa.utils.printErr("could not construct decay topology for wave descrption at index " + str(waveDescriptionID) + " of key file '" + keyFileName + "'.")
				return collections.OrderedDict()
			constructedWaveName = waveDescription.waveNameFromTopology(amplitude.decayTopology())
			if waveName == constructedWaveName:
				return waveDescription
		pyRootPwa.utils.printErr("none of the constructed topologies matched the given wave name.")
		return None


	def getWaveDescriptions(self):
		retval = collections.OrderedDict()
		for waveName in self.keyFiles:
			retval[waveName] = self.getWaveDescription(waveName)
		return retval


#	def getAmplitudeFilePaths(self, eventFileId, eventsType):
#		eventsType = fileManager.pyEventsType(eventsType)
#		retval = {}
#		for waveName_i, waveName in enumerate(self.keyFiles.keys()):
#			if (eventsType, eventFileId, waveName_i) in self.amplitudeFiles:
#				retval[waveName] = self.amplitudeDirectory + "/" + self.amplitudeFiles[(eventsType, eventFileId, waveName_i)]
#		return ampFileList


#	def getAmplitudeFilePath(self, eventFileId, waveName, eventsType):
#		eventsType = fileManager.pyEventsType(eventsType)
#		if (eventsType, eventFileId, waveName) not in self.amplitudeFiles:
#			pyRootPwa.utils.printWarn("no matching amplitude found.")
#			return ""
#		return self.amplitudeDirectory + "/" + self.amplitudeFiles[(eventsType, eventFileId, self.keyFiles.keys().index(waveName))]


#	def getAllAmplitudeFilePaths(self):
#		retval = []
#		for _, amplitudeFilePath in self.amplitudeFiles.iteritems():
#			retval.append(amplitudeFilePath)
#		return retval


	def getEventAndAmplitudeFilePathsInBin(self, multiBin, eventsType):
		# returns { "dataFileName" : { "waveName" : "amplitudeFileName" } }
		eventsType = fileManager.pyEventsType(eventsType)
		eventFileIds = self._getEventFileIdsForIntegralBin(multiBin, eventsType)
		if not eventFileIds:
			pyRootPwa.utils.printWarn("no matching event files found in bin '" + str(multiBin) + "' and events type '" + str(eventsType) + "'.")
			return collections.OrderedDict()
		retval = collections.OrderedDict()
		for eventFileId in eventFileIds:
			eventFileName = self.dataFiles[eventsType][eventFileId].dataFileName
			retval[eventFileName] = collections.OrderedDict()
			for waveName_i, waveName in enumerate(self.keyFiles.keys()):
				retval[eventFileName][waveName] = self.amplitudeFiles[(eventsType, eventFileId, waveName_i)]
		return retval


	def getIntegralFilePath(self, multiBin, eventsType):
		eventsType = fileManager.pyEventsType(eventsType)
		return self.integralFiles[(self.binList.index(multiBin), eventsType)]


#	def getBinIndex(self, binningInfo, checkConsistency = True):
#		# binInformation = { "variableName": value }
#		if not self.binList:
#			pyRootPwa.utils.printWarn("bin list is empty.")
#			return -1
#		foundIndex = -1
#		for binIndex, multiBin in enumerate(self.binList):
#			if multiBin.inBin(binningInfo):
#				if not checkConsistency:
#					return binIndex
#				if foundIndex > 0:
#					pyRootPwa.utils.printWarn("point " + str(binningInfo) + " lies in at least two bins: " + str(binIndex) + " and " + str(foundIndex))
#					return -1
#				foundIndex = binIndex
#		pyRootPwa.utils.printWarn("point " + str(binningInfo) + " does not lay in any bin.")
#		return -1
#
#
#	def getBin(self, binningInfo, checkConsistency = True):
#		binIndex = self.getBinIndex(binningInfo, checkConsistency)
#		if binIndex < 0:
#			return None
#		return self.binList[self.getBinIndex(binningInfo, checkConsistency)]


	def _getEventFileIdsForIntegralBin(self, multiBin, eventsType):
		eventsType = fileManager.pyEventsType(eventsType)
		if eventsType not in self.dataFiles:
			pyRootPwa.utils.printWarn("events type '" + str(eventsType) + "' not in data files.")
			return []
		eventFileIds = []
		for eventFileId, inputFile in enumerate(self.dataFiles[eventsType]):
			found = True
			for variableName in inputFile.binningMap:
				if variableName in multiBin.boundaries:
					if (inputFile.binningMap[variableName][1] < multiBin.boundaries[variableName][0]) or \
					   (inputFile.binningMap[variableName][0] > multiBin.boundaries[variableName][1]):
						found = False
						break
			if found:
				eventFileIds.append(eventFileId)
		return eventFileIds


#	def getEventFilePathsForForIntegralBin(self, multiBin, eventsType = None):
#		eventFileIds = self.getEventFileIds(multiBin, eventsType)
#		retval = collections.OrderedDict()
#		for evTyp, evFileIds in eventFileIds.iteritems():
#			retval[evTyp] = []
#			for evFileId in evFileIds:
#				retval[evTyp].append(self.dataFiles[evTyp][evFileId])
#		return retval


	def getWaveNameList(self):
		return self.keyFiles.keys()


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


	def _getAmplitudeFilePaths(self):
		amplitudeFiles = collections.OrderedDict()
		if not self.dataFiles:
			pyRootPwa.utils.printWarn("cannot create amplitude file path collection without data files.")
			return collections.OrderedDict()
		needSubdirs = self.limitFilesInDir != -1 and self._nmbAmplitudeFiles() > self.limitFilesInDir
		fileCount = 0
		dirSet = []
		for eventsType in self.dataFiles:
			for inputFile_i, _ in enumerate(self.dataFiles[eventsType]):
				for waveName_i, waveName in enumerate(self.keyFiles.keys()):
					amplitudeFileName = waveName + "_eventFileId-" + str(inputFile_i) + "_" + str(eventsType) + ".root"
					if needSubdirs:
						subDir = fileCount/self.limitFilesInDir
						fullDir = self.amplitudeDirectory + "/" + str(subDir)
						if not fullDir in dirSet:
							dirSet.append(fullDir)
						amplitudeFileName = str(subDir) + "/" + amplitudeFileName
						fileCount += 1
					else:
						amplitudeFileName = self.amplitudeDirectory + "/" + amplitudeFileName
					amplitudeFiles[(eventsType, inputFile_i, waveName_i)] = amplitudeFileName

		# make sure directories exist (create if neccessary) and check if they are empty
		if needSubdirs:
			for subDir in dirSet:
				if not os.path.isdir(subDir):
					os.mkdir(subDir)
					pyRootPwa.utils.printInfo("created folder for amplitude files: '" + subDir + "'.")
				else:
					if not os.listdir(subDir) == []:
						pyRootPwa.utils.printWarn("directory '" + subDir + "' is not empty.")
		else:
			if not os.listdir(self.amplitudeDirectory) == []:
				pyRootPwa.utils.printWarn("directory '" + self.amplitudeDirectory + "' is not empty.")

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
			                      eventMeta.binningMap(),
			                      fileManager.pyEventsType(eventMeta.eventsType()),
			                      eventMeta.additionalSavedVariableLables())
			if inputFile.eventsType not in inputFiles:
				inputFiles[inputFile.eventsType] = []
			inputFiles[inputFile.eventsType].append(inputFile)
			dataFile.Close()
		retval =  collections.OrderedDict()
		for eventsType in sorted(inputFiles.keys()):
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
					pyRootPwa.utils.printErr("could not construct decay topology for wave descrption at index " + str(waveDescriptionID) + " of key file '" + keyFileName + "'.")
					return collections.OrderedDict()
				waveName = waveDescription.waveNameFromTopology(amplitude.decayTopology())
				if waveName in keyFiles.keys():
					pyRootPwa.utils.printErr("duplicate wave name ('" + waveName +"' from files '" +
					                         keyFiles[waveName][0] + "' (index " + str(keyFiles[waveName][1]) +
					                         ") and '" + keyFileName + "' (index " + str(waveDescriptionID) + ").")
					return collections.OrderedDict()
				keyFiles[waveName] = keyFileName
		retval = collections.OrderedDict()
		for waveName in sorted(keyFiles.keys()):
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
				retStr += ("eventsType [" + str(eventsType) + "], bin [" + str(dataFile.binningMap) +
				           "] >> " + dataFile.dataFileName + "\n")
		retStr += "\nAmpFiles:\n"
		for eventsType in self.dataFiles:
			for eventFileId, _ in enumerate(self.dataFiles[eventsType]):
				for waveName_i, waveName in enumerate(self.keyFiles.keys()):
					retStr += ("eventsType [" + str(eventsType) + "], eventFileId [" + str(eventFileId) +
					           "], waveName [" + waveName + "] >> " + self.amplitudeFiles[(eventsType, eventFileId, waveName_i)] + "\n")
		retStr += "\nIntFiles:\n"
		for eventsType in [EventsType.GENERATED, EventsType.ACCEPTED]:
			for binID, _ in enumerate(self.binList):
				retStr += "eventsType [" + str(eventsType) + "], binID [" + str(binID) + "] >> " + self.integralFiles[(binID, eventsType)] + "\n"
		return retStr


	def isFilePerDirLimitReached(self):
		return len(self.binList) * len(self.keyFiles) * len(self.dataFiles) > self.limitFilesInDir and not self.limitFilesInDir == -1


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
		else:
			return -1
