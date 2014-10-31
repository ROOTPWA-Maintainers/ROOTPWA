
import ConfigParser
import os
import sys

import pyRootPwa.utils

class rootPwaConfig:

	config = None
	configFileName = ""

	# general section
	pdgFileName                            = ""
	massBinDirectoryNamePattern            = ""
	dataDirectory                          = ""
	phaseSpaceWeightFileExtensionQualifier = ""
	accCorrPSWeightFileExtensionQualifier  = ""
	fileManagerPath                        = ""

	# amplitude section
	keyDirectory                           = ""
	ampDirectory                           = ""
	dataAmplitudeDirectoryName             = ""
	phaseSpaceAmpDirectoryName             = ""
	accCorrPSAmpDirectoryName              = ""
	weightTreeName                         = ""

	# fit section
	fitResultTreeName                      = ""
	fitResultBranchName                    = ""


	def __init__(self, configFileName):
		self.config = ConfigParser.ConfigParser()
		self.configFileName = configFileName

		try:
			with open(configFileName, 'r') as configFile:
				self.config.readfp(configFile)
		except IOError:
			pyRootPwa.utils.printErr("config file '" + configFileName + "' could not be opened. Aborting...")
			sys.exit(1)
		except ConfigParser.Error:
			pyRootPwa.utils.printErr("config file '" + configFileName + "' could not be parsed. Aborting...")
			sys.exit(1)

		try:
			self.pdgFileName = os.path.expanduser(os.path.expandvars(self.config.get('general', 'particleDataTable')))

			self.dataDirectory                          = os.path.expanduser(os.path.expandvars(self.config.get('general', 'dataFileDirectory')))
			self.phaseSpaceWeightFileExtensionQualifier = self.config.get('general', 'phaseSpaceWeightFileExtensionQualifier')
			self.accCorrPSWeightFileExtensionQualifier  = self.config.get('general', 'accCorrPSWeightFileExtensionQualifier')
			self.massBinDirectoryNamePattern            = self.config.get('general', 'massBinDirectoryNamePattern')
			if self.config.has_option('general', 'fileManagerPath'):
				self.fileManagerPath                   = os.path.expanduser(os.path.expandvars(self.config.get('general', 'fileManagerPath')))
			else:
				self.fileManagerPath                   = os.path.expanduser("./fileManager.p")

			rawKeyfilePattern                           = os.path.expanduser(os.path.expandvars(self.config.get('amplitudes', 'keyfiles')))
			self.keyDirectory                           = os.path.expanduser(os.path.expandvars(self.config.get('amplitudes', 'keyFileDirectory')))
			self.ampDirectory                           = os.path.expanduser(os.path.expandvars(self.config.get('amplitudes', 'ampFileDirectory')))
			if self.config.has_option('amplitudes', 'limitFilesPerDir'):
				self.limitFilesPerDir                   = self.config.get('amplitudes', 'limitFilesPerDir')
			else:
				self.limitFilesPerDir                   = 0

			self.dataAmplitudeDirectoryName             = self.config.get('amplitudes', 'dataAmplitudeDirectoryName')
			self.phaseSpaceAmpDirectoryName             = self.config.get('amplitudes', 'phaseSpaceAmpDirectoryName')
			self.accCorrPSAmpDirectoryName              = self.config.get('amplitudes', 'accCorrPSAmpDirectoryName')

			self.weightTreeName                         = self.config.get('amplitudes', 'weightTreeName')

			rawKeyfilePattern = rawKeyfilePattern.replace('\n', '').replace(' ', '')
			if ';' in rawKeyfilePattern:
				self.keyfilePattern = rawKeyfilePattern.split(';')
			elif ',' in rawKeyfilePattern:
				self.keyfilePattern = rawKeyfilePattern.split(',')
			else:
				self.keyfilePattern = [rawKeyfilePattern]

			self.fitResultTreeName = self.config.get('fit', 'treeName')
			self.fitResultBranchName = self.config.get('fit', 'fitResultBranch')

		except ConfigParser.Error as exc:
			pyRootPwa.utils.printErr("a required entry was missing from the config file ('" + str(exc) + "'). Aborting...")
			sys.exit(1)

		if not os.path.isdir(self.dataDirectory):
			pyRootPwa.utils.printErr("Data directory invalid. Aborting...")
			sys.exit(1)
