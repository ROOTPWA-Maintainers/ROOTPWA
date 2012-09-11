
import ConfigParser
import os
import sys

import pyRootPwa.utils

class rootPwaConfig:

	config = None
	configFileName = ""

	pdgFileName = ""
	dataDirectory = ""
	dataFileExtensionQualifier = ""

	phaseSpaceEventFileExtenisonQualifier = ""
	accCorrPSEventFileExtensionQualifier  = ""
	massBinDirectoryNamePattern           = ""
	keyfilePattern                        = ""
	prodKinPartNamesObjName               = ""
	prodKinMomentaLeafName                = ""
	decayKinPartNamesObjName              = ""
	decayKinMomentaLeafName               = ""


	def __init__(self, configFileName):
		self.config = ConfigParser.ConfigParser()
		self.configFileName = configFileName

		try:
			with open(configFileName, 'r') as configFile:
				self.config.readfp(configFile)
		except IOError:
			pyRootPwa.utils.printErr("config file could not be opened. Aborting...")
			sys.exit(1)
		except ConfigParser.Error:
			pyRootPwa.utils.printErr("config file could not be parsed. Aborting...")
			sys.exit(1)

		try:
			self.pdgFileName                 = os.path.expanduser(os.path.expandvars(self.config.get('general', 'particleDataTable')))

			if self.config.has_option('general', 'dataDirectory'):
				self.dataDirectory = os.path.expanduser(os.path.expandvars(self.config.get('general', 'dataDirectory')))
				if self.dataDirectory == "":
					self.dataDirectory = os.path.abspath(os.path.dirname(configFileName))
			else:
				self.dataDirectory = os.path.abspath(os.path.dirname(configFileName))

			if not self.config.has_option('general', 'dataFileExtensionQualifier'):
				self.dataFileExtensionQualifier = ""
			else:
				self.dataFileExtensionQualifier = self.config.get('general', 'dataFileExtensionQualifier')

			self.phaseSpaceEventFileExtenisonQualifier = self.config.get('general', 'phaseSpaceEventFileExtenisonQualifier')
			self.accCorrPSEventFileExtensionQualifier  = self.config.get('general', 'accCorrPSEventFileExtensionQualifier')
			self.massBinDirectoryNamePattern           = self.config.get('general', 'massBinDirectoryNamePattern')
			self.keyfilePattern                        = os.path.expanduser(os.path.expandvars(self.config.get('amplitudes', 'keyfiles')))
			self.prodKinPartNamesObjName               = self.config.get('amplitudes', 'prodKinPartNamesObjName')
			self.prodKinMomentaLeafName                = self.config.get('amplitudes', 'prodKinMomentaLeafName')
			self.decayKinPartNamesObjName              = self.config.get('amplitudes', 'decayKinPartNamesObjName')
			self.decayKinMomentaLeafName               = self.config.get('amplitudes', 'decayKinMomentaLeafName')
		except ConfigParser.Error:
			pyRootPwa.utils.printErr("a required entry was missing from the config file. Aborting...")
			sys.exit(1)

		if not os.path.isdir(self.dataDirectory):
			pyRootPwa.utils.printErr("Data directory invalid. Aborting...")
			sys.exit(1)

		if self.phaseSpaceEventFileExtenisonQualifier == "":
			pyRootPwa.utils.printWarn("File extension qualifier for the phase space events is empty, no phase space events will be calculated...")
		if self.accCorrPSEventFileExtensionQualifier == "":
			pyRootPwa.utils.printWarn("File extension qualifier for the acceptance corrected phase space events is empty, no acc. cor. phase space events will be calculated...")

