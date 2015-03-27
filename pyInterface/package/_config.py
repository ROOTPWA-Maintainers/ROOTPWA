
import ConfigParser
import os

import pyRootPwa.utils

class rootPwaConfig:

	config = None
	configFileName = ""

	# general section
	pdgFileName                            = ""
	fileManagerPath                        = ""
	dataDirectory                          = ""
	keyDirectory                           = ""
	ampDirectory                           = ""
	limitFilesPerDir                       = -1

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

			if self.config.has_option('general', 'limitFilesPerDir'):
				self.limitFilesPerDir                   = self.config.get('general', 'limitFilesPerDir')
			else:
				self.limitFilesPerDir                   = -1

			self.fitResultTreeName = self.config.get('fit', 'treeName')
			self.fitResultBranchName = self.config.get('fit', 'fitResultBranch')

			self.phaseSpaceWeightFileExtensionQualifier = self.config.get('other', 'phaseSpaceWeightFileExtensionQualifier')
			self.accCorrPSWeightFileExtensionQualifier  = self.config.get('other', 'accCorrPSWeightFileExtensionQualifier')
			self.phaseSpaceAmpDirectoryName             = self.config.get('other', 'phaseSpaceAmpDirectoryName')
			self.accCorrPSAmpDirectoryName              = self.config.get('other', 'accCorrPSAmpDirectoryName')
			self.weightTreeName                         = self.config.get('other', 'weightTreeName')

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
		return True
