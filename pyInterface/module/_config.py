import math
import ast
import ConfigParser
import os
import errno

import pyRootPwa.utils

class rootPwaConfig(object):

	config = None
	configFileName = ""

	# general section
	pdgFileName                            = ""
	fileManagerPath                        = ""
	eventDirectory                         = ""
	keyDirectory                           = ""
	ampDirectory                           = ""
	intDirectory                           = ""

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

	_deprecatedMembers = {'dataDirectory': 'eventDirectory'}

	def __getattribute__(self, name, *args, **kwargs):
		if name in rootPwaConfig._deprecatedMembers:
			pyRootPwa.utils.printWarn("'{0}' member is deprecated. Use '{1}'!".format(name, rootPwaConfig._deprecatedMembers[name]))
			return object.__getattribute__(self, rootPwaConfig._deprecatedMembers[name])
		return object.__getattribute__(self, name, *args, **kwargs)

	def __setattr__(self, name, *args, **kwargs):
		if name in rootPwaConfig._deprecatedMembers:
			pyRootPwa.utils.printWarn("'{0}' member is deprecated. Use '{1}'!".format(name, rootPwaConfig._deprecatedMembers[name]))
			return object.__setattr__(self, rootPwaConfig._deprecatedMembers[name], *args, **kwargs)
		return object.__setattr__(self, name, *args, **kwargs)

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
			self.eventDirectory   = self.getPathFromConfig("general", "dataFileDirectory", configDir + "/data")
			self.keyDirectory    = self.getPathFromConfig("general", "keyFileDirectory" , configDir + "/keyfiles")
			self.ampDirectory    = self.getPathFromConfig("general", "ampFileDirectory" , configDir + "/amps")
			self.intDirectory    = self.getPathFromConfig("general", "intFileDirectory" , configDir + "/ints")

			if self.config.has_option('general', 'limitFilesPerDir'):
				pyRootPwa.utils.printWarn("'limitFilesPerDir' options is no longer supported! Ignoring it.")

			self.integralBinning = _readBinning(self.config.get("general", "integralBinning"))
			if not self.integralBinning:
				pyRootPwa.utils.printWarn("could not read integral binning string.")

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

		_createDirectory(self.eventDirectory, 'data')
		_createDirectory(self.keyDirectory, 'key')
		_createDirectory(self.ampDirectory, 'amplitude')
		_createDirectory(self.intDirectory, 'integral')
		return True


def _createDirectory(directoryPath, label):
	if not os.path.isdir(directoryPath):
		pyRootPwa.utils.printInfo("creating " + label + " directory '" + directoryPath + "'.")
		try:
			os.makedirs(directoryPath)
		except OSError as exception:  # ignore exception if the the directory was created meanwhile
			if exception.errno == errno.EEXIST and not os.path.isdir(directoryPath):  # this path was created meanwhile, but it is no directory
				pyRootPwa.utils.printErr("Path '{0}' exists but is no directory!".format(directoryPath))
				raise exception
			elif exception.errno != errno.EEXIST:
				raise exception


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
	binningString = binningString.rstrip('"')
	binningString = binningString.lstrip('"')
	binningString.replace("\\\"", "\"")
	try:
		inputVal = ast.literal_eval(binningString)
	except (SyntaxError, ValueError):
		pyRootPwa.utils.printWarn("error when converting binning string to python.")
		return []
	if not inputVal:
		return []
	retVal = []
	if isinstance(inputVal, list): # case 1.
		for binningItem in inputVal:
			try:
				retVal.append(pyRootPwa.utils.multiBin(binningItem))
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
			retVal.append(pyRootPwa.utils.multiBin(boundaries))
			indices[0] += 1
			for metaIndex in xrange(0, len(indices)-1):
				if indices[metaIndex] == limits[metaIndex]:
					indices[metaIndex] = 1
					indices[metaIndex+1] += 1
	else:
		pyRootPwa.utils.printWarn("binning string is neither a list nor a dict.")
		return []
	if retVal:
		for multiBin in retVal:
			if not multiBin.sameBinningVariables(retVal[0]):
				pyRootPwa.utils.printWarn("not all bins have the same binning variables.")
				return []
	return retVal


def _roundToNSignificantDigits(value, significantDigits):
	oomValue = int( math.floor(math.log10(abs(value))) )
	nDigits = significantDigits - 1 - oomValue
	return round(value, nDigits)


def _binningInputHandlingForCaseTwoAndThree(inputVal):
	# check the input dictionary for binning cases 2. and 3. (see above)
	# return a dictionary which conforms to case 2..
	for key in inputVal:
		if not isinstance(key, str):
			failMsg = "key in binning definition is not of type 'str'."
			pyRootPwa.utils.printWarn(failMsg)
			raise TypeError(failMsg)
		axesPartition = inputVal[key]
		if not isinstance(axesPartition, (list,tuple)):
			failMsg = "axes partition for variable '" + key + "' is not of type 'tuple' or 'list'."
			pyRootPwa.utils.printWarn(failMsg)
			raise TypeError(failMsg)
		for i in axesPartition:
			if not isinstance(i, (float,int)):
				failMsg = "axes parition for variable '" + key + "' is not a number."
				pyRootPwa.utils.printWarn(failMsg)
				raise TypeError(failMsg)
		if isinstance(axesPartition, tuple):
			if len(axesPartition) != 3:
				failMsg = "axes partition for variable '" + key + "' does not have length 3."
				pyRootPwa.utils.printWarn(failMsg)
				raise ValueError(failMsg)
			expandedBoundaryList = []
			xMin = float(axesPartition[0])
			xMax = float(axesPartition[1])
			nBins = int(axesPartition[2])
			for i in xrange(nBins+1):
				# round to 14 significant digits to get rid of numeric artifacts of the bin borders
				expandedBoundaryList.append(_roundToNSignificantDigits(xMin + (xMax-xMin)/nBins * i, 14))
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
