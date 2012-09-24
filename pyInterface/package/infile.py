
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
	_entriesInTree = 0
	_inFile = None

	_inWith = False

	def __init__(self, inFileName):
		self.inFileName = inFileName
		if pyRootPwa.config is None:
			raise pyRootPwa.exception.pyRootPwaException("pyRootPwa configuration not initialized")
		if inFileName.endswith('.root'):
			self._readRootFile = True
		elif inFileName.endswith('.evt'):
			self._readRootFile = False
			(self.prodKinParticles, self.decayKinParticles, self.tree) = pyRootPwa.utils.getTreeFromEvtFile(inFileName, inFileName)
			self._entriesInTree = int(self.tree.GetEntries())
		else:
			raise pyRootPwa.exception.pyRootPwaException("Unknown file extension")
		self._prodKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")
		self._decayKinMomenta = pyRootPwa.ROOT.TClonesArray("TVector3")

	def __len__(self):
		return self._entriesInTree

	def __enter__(self):
		if self._readRootFile:
			self._inFile = pyRootPwa.ROOT.TFile.Open(self.inFileName)
			self.prodKinParticles = self._inFile.Get(pyRootPwa.config.prodKinPartNamesObjName)
			self.decayKinParticles = self._inFile.Get(pyRootPwa.config.decayKinPartNamesObjName)
			self.tree = self._inFile.Get(pyRootPwa.config.inTreeName)
			self._entriesInTree = int(self.tree.GetEntries())
			self.tree = self._inFile.Get(pyRootPwa.config.inTreeName)
			self.tree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, self._prodKinMomenta)
			self.tree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, self._decayKinMomenta)
		else:
			self.tree.SetBranchAddress(pyRootPwa.config.prodKinMomentaLeafName, self._prodKinMomenta)
			self.tree.SetBranchAddress(pyRootPwa.config.decayKinMomentaLeafName, self._decayKinMomenta)
		self._inWith = True
		return self._inFile

	def __exit__(self, *args):
		if self._readRootFile:
			self._inFile.Close()
		self._inWith = False

	def __getitem__(self, index):
		if not self._inWith:
			raise pyRootPwa.exception.pyRootPwaException("Not in with statement")
		self.tree.GetEntry(index)
		return (self._prodKinMomenta, self._decayKinMomenta)

	def __str__(self):
		retval = ""
		retval += "pyRootPwa.inputFile\n"
		retval += "-------------------\n\n"
		retval += "inFileName = " + self.inFileName + "\n"
		retval += str(self.tree) + "\n"
		retval += str(self.prodKinParticles) + "\n"
		retval += str(self.decayKinParticles) + "\n\n"
		return retval
