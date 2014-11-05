
import glob
import os

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT


class fileManager:

	class InputFile:
		dataFileName = ""
		binningMap = {}
		eventsType = None
		def __init__(self, dataFileName, binningMap, eventsType):
			self.dataFileName = dataFileName
			self.binningMap = binningMap
			self.eventsType = eventsType

	class KeyFile:
		keyFileName = ""
		waveName = ""
		def __init__(self, keyFileName, waveName):
			self.keyFileName = keyFileName
			self.waveName = waveName

	dataDirectory      = ""
	keyDirectory       = ""
	amplitudeDirectory = ""

	limitFilesInDir = -1

	dataFiles       = {}
	keyFiles        = []
	amplitudeFiles  = {}
	globalAxes = {}
	binList = []

	def initialize(self, configFileName):
		config = pyRootPwa.rootPwaConfig()
		if not config.initialize(configFileName):
			pyRootPwa.utils.printErr("loading config file '" + configFileName + "' failed.")
			return False
		self.dataDirectory      = config.dataDirectory
		self.keyDirectory       = config.keyDirectory
		self.amplitudeDirectory = config.ampDirectory
		pyRootPwa.utils.printInfo("data file dir read from config file: '" + self.dataDirectory + "'.")
		pyRootPwa.utils.printInfo("key file dir read from config file: '" + self.keyDirectory + "'.")
		pyRootPwa.utils.printInfo("amplitude file dir read from config file: '" + self.amplitudeDirectory + "'.")

		self.dataFiles = self._openDataFiles()
		allAxes = []
		for eventsType in self.dataFiles:
			allAxes.append(self._getBinningAxes(self.dataFiles[eventsType]))
		self.globalAxes = self._combineAxes(allAxes)
		self.binList = self._createBinIDs()

		particleDataTableFileName = os.environ['ROOTPWA'] + "/particleData/particleDataTable.txt"
		if not pyRootPwa.core.particleDataTable.readFile(particleDataTableFileName):
			pyRootPwa.utils.printErr("error reading particle data table file '" + particleDataTableFileName + "'.")
		self.keyFiles = self._openKeyFiles()

		if int(config.limitFilesPerDir) == 0:
			self.limitFilesInDir = len(self.keyFiles) * len(self.dataFiles)
			pyRootPwa.utils.printInfo("limit for files per directory set to " + str(self.limitFilesInDir))
		elif int(config.limitFilesPerDir) == -1:
			self.limitFilesInDir = -1
			pyRootPwa.utils.printInfo("limit for files per directory set to infinite")
		else:
			self.limitFilesInDir = int(config.limitFilesPerDir)
			pyRootPwa.utils.printInfo("limit for files per directory set to " + str(self.limitFilesInDir))
		return True


	def getMissingBins(self):
		missingBins = []
		for binID in self.getBinIDList():
			bin = self.binList[binID]
			for eventsType in self.dataFiles:
				found = False
				for file in self.dataFiles[eventsType]:
					if file.binningMap == bin:
						found = True
						break
				if not found:
					missingBins.append(binID)
		return missingBins


	def getDataFilePaths(self):
		return fileManager.convertDataFilesToPaths(self.dataFiles)


	def getKeyFilePaths(self):
		return fileManager.convertKeyFilesToPaths(self.keyFiles)


	def getAmpFilePaths(self, eventsType):
		ampFileList = []
		for key in sorted(self.amplitudeFiles):
			if (key[2] == eventsType):
				ampFileList.append(self.amplitudeDirectory + "/" + self.amplitudeFiles[key])
		return ampFileList


	def getDataFile(self, binID, eventsType):
		foundFiles = []
		for file in self.dataFiles[eventsType]:
			found = True
			if not file.binningMap.keys() == self.binList[binID].keys():
				continue
			for binningVariable in file.binningMap:
				if not self.binList[binID][binningVariable] == file.binningMap[binningVariable]:
					found = False
					break
			if found: return file
		pyRootPwa.utils.printWarn("no data file found for binID = " + str(binID) + "and eventsType = '" + str(eventsType) + "'.")
		return False


	def getKeyFile(self, keyFileID):
		return self.keyFiles[keyFileID]


	def getAmplitudeFilePaths(self):
		pass


	def getAmplitudeFilePath(self, binID, keyFileID, eventsType):
		binInformation = self.getBinFromID(binID)
		keyFile = self.keyFiles[keyFileID]
		if self.isFilePerDirLimitReached():
			return self.amplitudeDirectory + "/" + self.getAmplitudeFileSubdir(binID, keyFileID) + "/" + keyFile.waveName + "_binID-" + str(binID) + "_" + str(eventsType) + ".root"
		else:
			return self.amplitudeDirectory + "/" + keyFile.waveName + "_binID-" + str(binID) + "_" + str(eventsType) + ".root"


	def getAmplitudeFileSubdir(self, binID, keyFileID):
		subDirNumber = (binID * len(self.keyFiles) + keyFileID) * len(self.dataFiles) / self.limitFilesInDir
		return str(subDirNumber)


	def getBinID(self, binInformation):
		foundBins = []
		for binID in range(len(self.binList)):
			found = True
			bin = self.binList[binID]
			for binningVariable in binInformation:
				lower = bin[binningVariable][0]
				upper = bin[binningVariable][1]
				given = binInformation[binningVariable]
				if lower > given or upper < given:
					found = False
					break
			if found: foundBins.append(binID)
		return foundBins


	def getBinFromID(self, id):
		if id not in self.getBinIDList():
			pyRootPwa.utils.printErr("id not found: " + str(id) + ".")
			raise Exception("do this properly")
		return self.binList[id]


	def getBinIDList(self):
		return range(len(self.binList))


	def getKeyFileIDList(self):
		return range(len(self.keyFiles))


	def _getBinningAxes(self, fileList):
		if not fileList:
			pyRootPwa.utils.printErr("got no input files")
			raise Exception("do this properly")
		binAxes = {}
		for binningVariable in fileList[0].binningMap:
			binAxes[binningVariable] = []
		for inputFile in fileList:
			if inputFile.binningMap.keys() != binAxes.keys():
				pyRootPwa.utils.printErr("data file '" + inputFile.dataFileName +
				                         "' seems to have different binning variables.")
				return {}
			for key in inputFile.binningMap:
				if inputFile.binningMap[key] not in binAxes[key]:
					binAxes[key].append(inputFile.binningMap[key])
		for binningVariable in binAxes:
			binAxes[binningVariable] = sorted(binAxes[binningVariable], key=lambda t: t[0])
			bins = binAxes[binningVariable]
			if not bins:
				pyRootPwa.util.printErr("no bins found for variable '" + binningVariable + "'.")
				return {}
			bound = bins[0]
			for bin in bins[1:]:
				if bound[1] < bin[0]:
					pyRootPwa.utils.printWarn("gap in bin structure found for binned variable '" +
					                          binningVariable + "' between bin '" + str(bound) +
					                          "' and bin '" + str(bin) + "'.")
				bound = bin
		return binAxes


	def _combineAxes(self, axes):
		if not axes:
			pyRootPwa.utils.printErr("got no axes.")
			raise Exception("do this properly")
		binningVariables = axes[0].keys()
		for axis in axes:
			if not axis.keys() == binningVariables:
				pyRootPwa.utils.printWarn("found different binning variables in different axes.")
		globalAxes = {}
		for binningVariable in binningVariables:
			binAxes = [ axis[binningVariable] for axis in axes if binningVariable in axis.keys() ]
			globalBinList = []
			for axis in binAxes:
				for bin in axis:
					if bin not in globalBinList:
						globalBinList.append(bin)
			globalBinList = sorted(globalBinList, key=lambda t: t[0])
			if not globalBinList:
				pyRootPwa.util.printErr("global bin list empty for '" + binningVariable + "'.")
				return {}
			bound = globalBinList[0]
			for bin in globalBinList[1:]:
				if bound[1] > bin[0]:
					pyRootPwa.utils.printErr("overlap in bin structure found for binned variable '" +
					                         binningVariable + "' between bin '" + str(bound) +
					                         "' and bin '" + str(bin) + "'.")
					return {}
				bound = bin
			globalAxes[binningVariable] = globalBinList
		pyRootPwa.utils.printSucc("combined bin axes: " + str(globalAxes.keys()))
		return globalAxes


	def _openKeyFiles(self):
		keyFileNames = glob.glob(self.keyDirectory + "/*.key")
		keyFiles = []

		for keyFileID in range(len(keyFileNames)):
			keyFileName = keyFileNames[keyFileID]
			waveDescription = pyRootPwa.core.waveDescription()
			waveDescription.parseKeyFile(keyFileName)
			(success, amplitude) = waveDescription.constructAmplitude()
			if not success:
				pyRootPwa.utils.printErr("could not construct decay topology for key file '" + keyFileName + "'.")
				return []
			waveName = waveDescription.waveNameFromTopology(amplitude.decayTopology(), True)
			if waveName in keyFiles:
				pyRootPwa.utils.printErr("duplicate wave name ('" + waveName + "' from files '" + keyFiles[waveName] + "' and '" + keyFileName + "'.")
				return []
			keyFile = fileManager.KeyFile(keyFileName, waveName)
			keyFiles.append(keyFile)
		return keyFiles


	def _openDataFiles(self):
		dataFileNames = glob.glob(self.dataDirectory + "/*.root")
		inputFiles = {}

		for dataFileName in dataFileNames:
			dataFile = ROOT.TFile.Open(dataFileName, "READ")
			if not dataFile:
				pyRootPwa.utils.printErr("could not open event file '" + dataFileName + "'.")
				return {}
			eventMeta = pyRootPwa.core.eventMetadata.readEventFile(dataFile)
			if not eventMeta:
				pyRootPwa.utils.printErr("could not find metadata in event file '" + dataFileName + "'.")
				return {}
			inputFile = fileManager.InputFile(dataFileName, eventMeta.binningMap(), eventMeta.eventsType())
			if inputFile.eventsType not in inputFiles:
				inputFiles[inputFile.eventsType] = []
			inputFiles[inputFile.eventsType].append(inputFile)
			dataFile.Close()

		if not inputFiles:
			pyRootPwa.utils.printErr("no binning maps found.")
			return {}
		for type in inputFiles:
			if not inputFiles[type]:
				pyRootPwa.utils.printErr("no binning maps found for type '" + type + "'.")
				return {}

		return inputFiles


	def _createBinIDs(self):
		binList = self._iterateBins(0, {}, [])
		return binList


	def _iterateBins(self, dim, currentBin, result):
		keys = self.globalAxes.keys()
		if dim < len(self.globalAxes):
			for bin in self.globalAxes[keys[dim]]:
				newCurrentBin = currentBin.copy()
				newCurrentBin[keys[dim]] = bin
				result = self._iterateBins(dim+1, newCurrentBin, result)
		elif dim == len(self.globalAxes):
			result.append(currentBin)
			return result
		else: raise Exception("do this properly")
		return result


	def isFilePerDirLimitReached(self):
		return len(self.binList) * len(self.keyFiles) * len(self.dataFiles) > self.limitFilesInDir or not self.limitFilesInDir == -1


	def areDataFilesSynced(self):
		return self.getDataFilePaths() == fileManager.convertDataFilesToPaths(self._openDataFiles())


	def areKeyFilesSynced(self):
		return self.getKeyFilePaths() == fileManager.convertKeyFilesToPaths(self._openKeyFiles())


	def areAmpFilesSynced(self):
		return False


	def areFilesSynced(self):
		return self.areDataFilesSynced() and self.areDataFilesSynced() and self.areAmpFilesSynced()


	@staticmethod
	def convertDataFilesToPaths(dataFilesDict):
		allDataFiles = []
		for eventsType in dataFilesDict:
			for dataFile in dataFilesDict[eventsType]:
				allDataFiles.append(dataFile.dataFileName)
		return allDataFiles

	@staticmethod
	def convertKeyFilesToPaths(keyFilesList):
		allKeyFiles = []
		for keyFile in keyFilesList:
			allKeyFiles.append(keyFile.keyFileName)
		return allKeyFiles
