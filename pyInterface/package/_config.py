
import ast
import ConfigParser
import os

import pyRootPwa.utils

class rootPwaConfig(object):

	config = None
	configFileName = ""

	# general section
	pdgFileName                            = ""
	fileManagerPath                        = ""
	dataDirectory                          = ""
	keyDirectory                           = ""
	ampDirectory                           = ""
	intDirectory                           = ""
	limitFilesPerDir                       = -1

	integralBinning                        = []

	# amplitude section
	phaseSpaceIntegralDirectory            = ""
	phaseSpaceUpperMassBound               = 0.

	# fit section
	fitResultTreeName                      = ""
	fitResultBranchName                    = ""

	# other
	phaseSpaceWeightFileExtensionQualifier = ""
	accCorrPSWeightFileExtensionQualifier  = ""
	phaseSpaceAmpDirectoryName             = ""
	accCorrPSAmpDirectoryName              = ""
	weightTreeName                         = ""

	def getPathFromConfig(self, category, varName, defaultValue):
		if self.config.has_option(category, varName):
			path = os.path.expanduser(os.path.expandvars(self.config.get(category, varName)))
			if not os.path.isabs(path):
				path = os.path.abspath(os.path.dirname(os.path.abspath(self.configFileName)) + "/" + path)
		else:
			path = defaultValue
			pyRootPwa.utils.printWarn("can not find '" + varName + "' in category '" + category + "'. Using default value '" + defaultValue + "'.")
		return path


	def initialize(self, configFileName):
		self.config = ConfigParser.ConfigParser()
		self.configFileName = configFileName
		configDir = os.path.dirname(os.path.abspath(self.configFileName))

		try:
			with open(configFileName, 'r') as configFile:
				self.config.readfp(configFile)
		except IOError:
			pyRootPwa.utils.printErr("config file '" + configFileName + "' could not be opened.")
			return False
		except ConfigParser.Error:
			pyRootPwa.utils.printErr("config file '" + configFileName + "' could not be parsed.")
			return False

		try:
			self.pdgFileName     = self.getPathFromConfig("general", "particleDataTable", os.path.expandvars("$ROOTPWA/particleData/particleDataTable.txt"))
			self.fileManagerPath = self.getPathFromConfig("general", "fileManagerPath"  , configDir + "/fileManager.pkl")
			self.dataDirectory   = self.getPathFromConfig("general", "dataFileDirectory", configDir + "/data")
			self.keyDirectory    = self.getPathFromConfig("general", "keyFileDirectory" , configDir + "/keyfiles")
			self.ampDirectory    = self.getPathFromConfig("general", "ampFileDirectory" , configDir + "/amps")
			self.intDirectory    = self.getPathFromConfig("general", "intFileDirectory" , configDir + "/ints")

			if self.config.has_option('general', 'limitFilesPerDir'):
				self.limitFilesPerDir                   = int(self.config.get('general', 'limitFilesPerDir'))
			else:
				self.limitFilesPerDir                   = -1

			self.integralBinning = _readBinning(self.config.get("general", "integralBinning"))
			if not self.integralBinning:
				pyRootPwa.utils.printErr("could not read integral binning string.")
				return False

			self.phaseSpaceIntegralDirectory = self.config.get("amplitude", "phaseSpaceIntegralDirectory")
			self.phaseSpaceUpperMassBound    = float(self.config.get("amplitude", "phaseSpaceUpperMassBound"))

			self.fitResultTreeName = self.config.get('fit', 'treeName')
			self.fitResultBranchName = self.config.get('fit', 'fitResultBranch')

			self.phaseSpaceWeightFileExtensionQualifier = self.config.get('other', 'phaseSpaceWeightFileExtensionQualifier')
			self.accCorrPSWeightFileExtensionQualifier  = self.config.get('other', 'accCorrPSWeightFileExtensionQualifier')
			self.phaseSpaceAmpDirectoryName             = self.config.get('other', 'phaseSpaceAmpDirectoryName')
			self.accCorrPSAmpDirectoryName              = self.config.get('other', 'accCorrPSAmpDirectoryName')
			self.weightTreeName                         = self.config.get('other', 'weightTreeName')

		except ValueError as exc:
			pyRootPwa.utils.printErr("a variable had the wrong type ('" + str(exc) + "').")
			return False
		except ConfigParser.Error as exc:
			pyRootPwa.utils.printErr("a required entry was missing from the config file ('" + str(exc) + "').")
			return False

		if not os.path.isdir(self.dataDirectory):
			os.mkdir(self.dataDirectory)
			pyRootPwa.utils.printInfo("created data directory '" + self.dataDirectory + "'.")
		if not os.path.isdir(self.keyDirectory):
			os.mkdir(self.keyDirectory)
			pyRootPwa.utils.printInfo("created key directory '" + self.keyDirectory + "'.")
		if not os.path.isdir(self.ampDirectory):
			os.mkdir(self.ampDirectory)
			pyRootPwa.utils.printInfo("created amplitude directory '" + self.ampDirectory + "'.")
		if not os.path.isdir(self.intDirectory):
			os.mkdir(self.intDirectory)
			pyRootPwa.utils.printInfo("created integral directory '" + self.intDirectory + "'.")
		return True


def _readBinning(binningString):
	# Three cases:
	#
	# 1. list of bins [ { "name": (lowBound, highBound)} ]
	# 2. grid coordinates for axes { "name": [ binBoundaries ] }
	# 3. edges of rectangular grid { "name": (lowBound, highBound, nBins) }
	#
	# Case 3. should be converted to case 2. and case 2. to case 1.
	#
	# returns list of bins as in 1. or empty list in case of errors.
	try:
		inputVal = ast.literal_eval(binningString)
	except (SyntaxError, ValueError):
		pyRootPwa.utils.printWarn("error when converting binning string to python.")
		return []
	retval = []
	if isinstance(inputVal, list): # case 1.
		for binningItem in inputVal:
			try:
				retval.append(pyRootPwa.utils.multiBin(binningItem))
			except (TypeError, ValueError):
				pyRootPwa.utils.printWarn("could not convert entry in binning list to multiBin.")
				return []
	elif isinstance(inputVal, dict): # case 2. + 3.
		try:
			inputVal = _binningInputHandlingForCaseTwoAndThree(inputVal)
		except (ValueError, TypeError):
			return []
		keys = inputVal.keys()
		indices = [1] * len(keys)
		limits = [ len(inputVal[key]) for key in keys ]
		while indices[len(indices)-1] != limits[len(limits)-1]:
			boundaries = {}
			for i, key in enumerate(keys):
				boundaries[key] = (inputVal[key][indices[i]-1], inputVal[key][indices[i]])
			retval.append(pyRootPwa.utils.multiBin(boundaries))
			indices[0] += 1
			for metaIndex in xrange(0, len(indices)-1):
				if indices[metaIndex] == limits[metaIndex]:
					indices[metaIndex] = 1
					indices[metaIndex+1] += 1
	else:
		pyRootPwa.utils.printWarn("binning string is neither a list nor a dict.")
		return []
	return retval


def _binningInputHandlingForCaseTwoAndThree(inputVal):
	# check the input dictionary for binning cases 2. and 3. (see above)
	# return a dictionary which conforms to case 2..
	for key in inputVal:
		if not isinstance(key, str):
			errMsg = "key in binning definition is not of type 'str'."
			pyRootPwa.utils.printWarn(errMsg)
			raise TypeError(errMsg)
		axesPartition = inputVal[key]
		if not (isinstance(axesPartition, list) or isinstance(axesPartition, tuple)):
			errMsg = "axes partition for variable '" + key + "' is not of type 'tuple' or 'list'."
			pyRootPwa.utils.printWarn(errMsg)
			raise TypeError(errMsg)
		for i in axesPartition:
			if not (isinstance(i, float) or isinstance(i, int)):
				errMsg = "axes parition for variable '" + key + "' is not a number."
				pyRootPwa.utils.printWarn(errMsg)
				raise TypeError(errMsg)
		if isinstance(axesPartition, tuple):
			if len(axesPartition) != 3:
				errMsg = "axes partition for variable '" + key + "' does not have length 3."
				pyRootPwa.utils.printWarn(errMsg)
				raise ValueError(errMsg)
			expandedBoundaryList = []
			for i in xrange(axesPartition[2]+1):
				expandedBoundaryList.append(float(axesPartition[0]) + (float(axesPartition[1] - axesPartition[0]) / float(axesPartition[2])) * i)
			inputVal[key] = expandedBoundaryList
	return inputVal


def _testReadBinning():
	testCases = []
	testString = "["
	for i in xrange(4):
		for j in xrange(4):
			testString += "{ \"mass\": (" + str(1.0+0.1*i) + ", " + str(1.0+0.1*(i+1)) + "), \"tPrime\": (" + str(1.0+0.1*j) + ", " + str(1.0+0.1*(j+1)) + ") }, "
	testString = testString[:-2]
	testString += "]"
	testCases.append(testString)
	testCases.append("{ \"mass\": [ 1.0, 1.5, 1.9, 1.95, 2.0 ], \"tPrime\": [ 1.0, 1.05, 1.1, 2.0 ] }")
	testCases.append("{ \"mass\": ( 1.0, 1.4, 4 ), \"tPrime\": ( 2.0, 2.4, 4 ) }")
	testCases.append("{ \"mass\": ( 1.0, 1.4, 4 ), \"tPrime\": [ 1.0, 1.05, 1.1, 2.0 ] }")
	testCases.append("{ \"mass\": [ 1.0, 1.5, 1.9, 1.95, 2.0 ], \"tPrime\": ( 2.0, 2.4, 4 ) }")
	for testCase in testCases:
		binning = _readBinning(testCase)
		pyRootPwa.utils.printDebug(testCase + "\n")
		for mBin in binning:
			pyRootPwa.utils.printDebug(str(mBin))
		pyRootPwa.utils.printDebug("##################################################")
