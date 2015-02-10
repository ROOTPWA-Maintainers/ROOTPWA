
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

	def __init__(self):
		pass


	def getDataFilePaths(self):
		pass


	def getDataFilePath(self, binInformation, eventsType):
		pass


	def getAmplitudeFilePaths(self):
		pass


	def getAmplitudeFilePath(self, binInformation, waveName, eventsType):
		pass


	def _getBinningAxes(self, inputFiles):
		if not inputFiles:
			pyRootPwa.utils.printErr("got no input files")
			raise Exception("do this properly")
		binAxes = {}
		for binningVariable in inputFiles[0].binningMap:
			binAxes[binningVariable] = []
		for inputFile in inputFiles:
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
