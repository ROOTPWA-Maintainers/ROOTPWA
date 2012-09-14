
import collections
import multiprocessing

import pyRootPwa
import pyRootPwa.exception
import pyRootPwa.utils

class inputFile():

	inFileName = ""

	tree = None
	prodKinParticles = None
	decayKinParticles = None

	_readRootFile = None
	_first = None
	_entriesInTree = 0

	def __init__(self, inFileName):
		self.inFileName = inFileName
		if pyRootPwa.config is None:
			raise pyRootPwa.exception.pyRootPwaException("pyRootPwa configuration not initialized")
		if inFileName.endswith('.root'):
			self._readRootFile = True
			self._inFile = pyRootPwa.ROOT.TFile.Open(inFileName)
			self.prodKinParticles = self._inFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
			self.decayKinParticles = self._inFile.Get(pyRootPwa.config.decayKinPartNamesObjName)
			self.tree = self._inFile.Get(pyRootPwa.config.inTreeName)
			self._entriesInTree = int(self.tree.GetEntries())
		elif inFileName.endswith('.evt'):
			self._readRootFile = False
			self._first = True
			(self.prodKinParticles, self.decayKinParticles, self.tree) = pyRootPwa.utils.getTreeFromEvtFile(inFileName, inFileName)
			self._entriesInTree = int(self.tree.GetEntries())
		else:
			raise pyRootPwa.exception.pyRootPwaException("Unknown file extension")

		self._prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
		self._decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
		self.tree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, self._prodKinMomenta)
		self.tree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, self._decayKinMomenta)

	def __len__(self):
		return self._entriesInTree

	def __getitem__(self, inslice):
		retval = []
		for index in range(inslice.start, inslice.stop):
			if self._readRootFile and not self._inFile.IsOpen():
				self._inFile = pyRootPwa.ROOT.TFile.Open(self.inFileName)
				self.tree = self._inFile.Get(pyRootPwa.config.inTreeName)
				self.tree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, self._prodKinMomenta)
				self.tree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, self._decayKinMomenta)
			elif self._first:
				self.tree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, self._prodKinMomenta)
				self.tree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, self._decayKinMomenta)
				self._first = False
			self.tree.GetEntry(index)
			retval.append((self._prodKinMomenta.Clone(), self._decayKinMomenta.Clone()))
		return retval

	def __str__(self):
		retval = ""
		retval += "pyRootPwa.inputFile\n"
		retval += "-------------------\n\n"
		retval += "inFileName = " + self.inFileName + "\n"
		retval += str(self.tree) + "\n"
		retval += str(self.prodKinParticles) + "\n"
		retval += str(self.decayKinParticles) + "\n\n"
		return retval
