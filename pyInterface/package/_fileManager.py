
import glob
import os.path

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT


class fileManager:

	class InputFile:
		dataFileName = ""
		binningMap = {}
		eventsType = None

	dataDirectory = ""
	amplitudeDirectory = ""

	dataFiles = {}
	globalAxes = {}
	binList = []

	def __init__(self, dataDirectory):
		self.dataDirectory = dataDirectory
		self.dataFiles = self._openInputFiles()
		allAxes = []
		for eventsType in self.dataFiles:
			allAxes.append(self._getBinningAxes(self.dataFiles[eventsType]))
		self.globalAxes = self._combineAxes(allAxes)
		self.binList = self._createBinIDs()

	def getMissingBins(self):
		missingBins = {}
		for bin in self.binList:
			for eventsType in self.dataFiles:
				found = False
				for file in self.dataFiles[eventsType]:
					if file.binningMap == bin:
						found = True
						break
				if not found:
					if not eventsType in missingBins: missingBins[eventsType]=[]
					missingBins[eventsType].append(bin.copy())
		return missingBins

	def getDataFilePaths(self):
		allDataFiles =[]
		for eventsType in self.dataFiles:
			for dataFile in self.dataFiles[eventsType]:
				allDataFiles.append(dataFile.dataFileName)
		return allDataFiles

	def getDataFilePath(self, binInformation, eventsType):
		foundFiles = []
		for file in self.dataFiles[eventsType]:
			found = True
			for binningVariable in binInformation:
				lower = file.binningMap[binningVariable][0]
				upper = file.binningMap[binningVariable][1]
				givenValue = binInformation[binningVariable]
				if givenValue < lower or givenValue > upper:
					found = False
					break
				elif givenValue == lower or givenValue == upper:
					pyRootPwa.utils.printErr("hit boundary of the " + binningVariable + "-bin.")
					return []
			if found: foundFiles.append(file.dataFileName)
		return foundFiles


	def getAmplitudeFilePaths(self):
		pass


	def getAmplitudeFilePath(self, binInformation, waveName, eventsType):
		pass


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
		if (not id < len(self.binList)) or (id < 0):
			pyRootPwa.utils.printErr("id not found: " + str(id))
			raise Exception("do this properly")
		return self.binList[id]

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


	def _openInputFiles(self):
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
			inputFile = self.InputFile()
			inputFile.dataFileName = dataFileName
			inputFile.binningMap = eventMeta.binningMap()
			inputFile.eventsType = eventMeta.eventsType()
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
		if not dim >= len(self.globalAxes):
			for bin in self.globalAxes[keys[dim]]:
				newCurrentBin = currentBin.copy()
				newCurrentBin[keys[dim]] = bin
				result = self._iterateBins(dim+1, newCurrentBin, result)
		else:
			result.append(currentBin)
			return result
		return result
